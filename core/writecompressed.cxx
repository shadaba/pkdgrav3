/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 */
#include "writecompressed.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "pkd.h"
#include "snapcomp.h"

static_assert(std::is_standard_layout_v<ServiceWriteCompressed::input>);
static_assert(std::is_standard_layout_v<ServiceWriteCompressed::output>);

namespace {
// One FileWriter per file, shared by the threads of this process. Created by the first thread that asks,
// closed (index, trailer, fsync) by the last one to release it.
struct Shared {
    std::unique_ptr<snapcomp::FileWriter> w;
    int users = 0;
};
std::mutex registry_mutex;
std::map<std::string, Shared> registry;

snapcomp::FileWriter* acquire(const std::string& path, const snapcomp::Header& h, int users, bool write_through) {
    std::lock_guard<std::mutex> lock(registry_mutex);
    auto& e = registry[path];
    if (!e.w) {
        e.w = std::make_unique<snapcomp::FileWriter>(path, h);
        e.w->set_write_through(write_through);
        e.users = users;
    }
    return e.w.get();
}

// MemAvailable from /proc/meminfo in GB (Linux); a huge value elsewhere.
double mem_available_gb() {
    double gb = 1e30;
    if (FILE* f = std::fopen("/proc/meminfo", "r")) {
        char line[256];
        unsigned long long kb;
        while (std::fgets(line, sizeof line, f))
            if (std::sscanf(line, "MemAvailable: %llu kB", &kb) == 1) {
                gb = double(kb) / (1024.0 * 1024.0);
                break;
            }
        std::fclose(f);
    }
    return gb;
}

// Decided once per snapshot and process, by whichever thread gets here first.
std::mutex decide_mutex;
std::string decided_for;
bool decided_through = false;
double decided_avail = 1e30;
bool write_through_for(const char* prefix, int mode, double bytes) {
    std::lock_guard<std::mutex> lock(decide_mutex);
    if (decided_for != prefix) {
        decided_for = prefix;
        decided_avail = mem_available_gb();
        decided_through = mode == 1 || (mode < 0 && decided_avail * 1e9 * 1.0737 < 2.0 * bytes);
    }
    return decided_through;
}

// Returns true if this call closed the file (i.e. this thread was the last user).
bool release(const std::string& path) {
    std::unique_ptr<snapcomp::FileWriter> last;
    {
        std::lock_guard<std::mutex> lock(registry_mutex);
        auto it = registry.find(path);
        if (it == registry.end()) return false;
        if (--it->second.users == 0) {
            last = std::move(it->second.w);
            registry.erase(it);
        }
    }
    if (!last) return false;
    last->close(); // outside the lock: fsync can take a while
    return true;
}

// "<prefix>.done": all files of this process are complete. The staging drain (tools/staging/pkd_drain.py)
// copies <prefix>.* when it sees it, whatever the log verbosity. Written by the thread that closed the last
// file; every thread handles the files in the same order, so all earlier files are closed by then.
void write_marker(const char* prefix, unsigned nfiles, int process) {
    std::string path = std::string(prefix) + ".done";
    if (FILE* f = std::fopen(path.c_str(), "w")) {
        std::fprintf(f, "{\"files\":%u,\"process\":%d}\n", nfiles, process);
        std::fclose(f);
    }
    else throw std::runtime_error("cannot write " + path);
}

std::string file_name(const char* prefix, double percent, int process) {
    char suffix[64];
    std::snprintf(suffix, sizeof suffix, ".s%02d.%d.pkc", int(std::lround(percent)), process);
    return std::string(prefix) + suffix;
}
} // namespace

int ServiceWriteCompressed::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto& in = *static_cast<input*>(vin);
    auto& out = *static_cast<output*>(vout);
    out = output();
    auto pkd = pst->plcl->pkd;
    auto mdl = pst->mdl;
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    try {
        snapcomp::Header h;
        h.n_total = in.nTotal;
        h.lattice_n = in.nLattice;
        h.ncell = in.nCell;
        h.sub = in.nSub;
        h.vel_step = in.dVelStep;
        h.vel_to_kms = in.dVelToKms;
        h.box_size = in.dBoxMpch;
        h.redshift = in.dRedshift;
        h.time = in.dTime;
        h.omega_m = in.dOmega0;
        h.omega_l = in.dLambda;
        h.h = in.dHubble;
        h.n_subsamples = in.nPercents;
        for (unsigned k = 0; k < in.nPercents; ++k) h.percents[k] = in.dPercents[k];
        h.subsample_seed = in.iSeed;
        h.writer = uint32_t(mdl->Proc());
        h.n_writers = uint32_t(mdl->Procs());
        h.block_size = in.nBlock;
        h.zstd_level = in.iLevel;
        h.step = in.iStep;
        h.info.assign(in.achInfo, strnlen(in.achInfo, sizeof in.achInfo));
        h.pot_step = in.dPotStep;
        h.pot_to_kms2 = in.dPotToKms2;
        const bool with_pot = h.has_potential();
        if (with_pot && !pkd->particles.present(PKD_FIELD::oPotential))
            throw std::runtime_error("dCompPotStep > 0 needs the potential in memory (bMemPotential)");

        snapcomp::SubsampleTable table(in.dPercents, int(in.nPercents), in.iSeed);
        snapcomp::BlockEncoder enc(h);
        const size_t B = in.nBlock;
        std::vector<double> bp;
        std::vector<float> bv;
        std::vector<uint64_t> bi;
        std::vector<float> bf;
        bp.reserve(3 * B);
        bv.reserve(3 * B);
        bi.reserve(B);
        if (with_pot) bf.reserve(B);
        const double dvFac = in.dvFac;
        const unsigned thread = unsigned(mdl->Core());

        const bool through = write_through_for(in.achPrefix, in.iWriteThrough, in.dBytesPerProcess);
        out.nWriteThrough = through ? 1 : 0;
        out.nThreads = 1;
        {
            std::lock_guard<std::mutex> lock(decide_mutex);
            out.dMemAvailGB = decided_avail;
        }
        for (unsigned k = 0; k < in.nPercents; ++k) {
            snapcomp::Header hk = h;
            hk.subsample = k;
            const std::string path = file_name(in.achPrefix, in.dPercents[k], mdl->Proc());
            snapcomp::FileWriter* w = acquire(path, hk, mdl->Cores(), through);
            auto flush = [&] {
                if (bi.empty()) return;
                auto e0 = clock::now();
                const auto& blk = enc.encode(bi.size(), bp.data(), bv.data(), bi.data(), with_pot ? bf.data() : nullptr);
                auto e1 = clock::now();
                w->append(blk.data(), blk.size(), uint32_t(bi.size()), thread);
                out.dEncode += std::chrono::duration<double>(e1 - e0).count();
                out.dWrite += std::chrono::duration<double>(clock::now() - e1).count();
                out.nBytes += blk.size();
                out.nPotBytes += enc.stats.pot;
                out.nWritten += bi.size();
                bp.clear(), bv.clear(), bi.clear(), bf.clear();
            };
            // One pass per subsample keeps the scratch to a single block.
            for (auto& p : pkd->particles) {
                if (p.species() != FIO_SPECIES_DARK) {
                    if (k == 0) ++out.nSkipped;
                    continue;
                }
                uint64_t id;
                if (p.have_particle_id()) id = p.ParticleID();
                else if (!pkd->bNoParticleOrder) id = p.order();
                else throw std::runtime_error("compressed output needs particle IDs (bMemParticleID) or particle order");
                if (table(id) != int(k)) continue;
                const auto r = p.position();
                const auto v = p.velocity();
                for (int d = 0; d < 3; ++d) {
                    double x = (r[d] - in.dCenter[d]) / in.dPeriod[d];
                    bp.push_back(x - std::floor(x + 0.5)); // periodic wrap into [-0.5, 0.5)
                    bv.push_back(float(v[d] * dvFac));
                }
                bi.push_back(id);
                if (with_pot) bf.push_back(p.potential());
                if (bi.size() == B) flush();
            }
            flush();
            if (release(path) && k + 1 == in.nPercents) write_marker(in.achPrefix, in.nPercents, mdl->Proc());
        }
        out.nScratch = enc.scratch_bytes() + bp.capacity() * sizeof(double) + bv.capacity() * sizeof(float) + bi.capacity() * sizeof(uint64_t) +
                       bf.capacity() * sizeof(float);
        out.nScratchMax = out.nScratch;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ERROR: compressed output on thread %d: %s\n", mdl->Self(), e.what());
        std::fflush(stderr);
        std::abort();
    }
    out.dMaxThread = std::chrono::duration<double>(clock::now() - t0).count();
    return sizeof(output);
}
