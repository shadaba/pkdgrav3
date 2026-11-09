// pkc_to_gadget: compressed snapshot (<prefix>.s*.pkc) -> Gadget-2 binary files (SnapFormat 1) for
// Rockstar (FILE_FORMAT = "GADGET2"). Stopgap until Rockstar reads .pkc natively.
//
// usage: pkc_to_gadget [--threads T] [--max-per-file M] [--percents 1,3] prefix out_prefix
//   writes out_prefix.0 .. out_prefix.<K-1> and prints the Rockstar settings to use.
//
// Gadget conventions used (as Rockstar's io_gadget.c expects):
//   positions  float32, Mpc/h, in [0, L)            -> GADGET_LENGTH_CONVERSION = 1
//   velocities float32, v_peculiar / sqrt(a), km/s  (Rockstar multiplies by sqrt(a))
//   IDs        uint64                               -> GADGET_ID_BYTES = 8
//   mass       MassTable[1] in 1e10 Msun/h           -> GADGET_MASS_CONVERSION = 1e10
// The particle mass is rho_mean * L^3 / (particles written), so a subsample keeps the mean density.
// Each Fortran record must stay below 2^31 bytes, hence at most M particles per file (default 2^24).
#include "snapcomp.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace snapcomp;

#pragma pack(push, 1)
struct GadgetHeader { // 256 bytes
    uint32_t npart[6];
    double mass[6];
    double time, redshift;
    int32_t flag_sfr, flag_feedback;
    uint32_t npartTotal[6];
    int32_t flag_cooling, num_files;
    double BoxSize, Omega0, OmegaLambda, HubbleParam;
    int32_t flag_stellarage, flag_metals;
    uint32_t npartTotalHighWord[6];
    int32_t flag_entropy_instead_u;
    char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4 - 4];
};
#pragma pack(pop)
static_assert(sizeof(GadgetHeader) == 256, "Gadget header must be 256 bytes");

static std::vector<std::string> files_of(const std::string& prefix) {
    auto slash = prefix.rfind('/');
    std::string dir = slash == std::string::npos ? "." : prefix.substr(0, slash);
    std::string base = slash == std::string::npos ? prefix : prefix.substr(slash + 1);
    std::vector<std::string> out;
    DIR* d = opendir(dir.c_str());
    if (!d) throw std::runtime_error("cannot open directory " + dir);
    while (auto* e = readdir(d)) {
        std::string n = e->d_name;
        if (n.size() > base.size() + 6 && n.compare(0, base.size() + 2, base + ".s") == 0 && n.compare(n.size() - 4, 4, ".pkc") == 0)
            out.push_back(dir + "/" + n);
    }
    closedir(d);
    std::sort(out.begin(), out.end());
    if (out.empty()) throw std::runtime_error("no files " + prefix + ".s*.pkc");
    return out;
}

static void record(FILE* f, const void* data, uint64_t bytes) {
    if (bytes >= (1ull << 31)) throw std::runtime_error("Gadget record too large; lower --max-per-file");
    uint32_t b = uint32_t(bytes);
    if (fwrite(&b, 4, 1, f) != 1 || (bytes && fwrite(data, 1, bytes, f) != bytes) || fwrite(&b, 4, 1, f) != 1) throw std::runtime_error("write failed");
}

int main(int argc, char** argv) {
    int T = int(std::thread::hardware_concurrency());
    uint64_t maxper = 1ull << 24;
    std::vector<double> want;
    std::string prefix, out;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--threads" && i + 1 < argc) T = std::stoi(argv[++i]);
        else if (a == "--max-per-file" && i + 1 < argc) maxper = std::stoull(argv[++i]);
        else if (a == "--percents" && i + 1 < argc) {
            std::string s = argv[++i];
            for (size_t p = 0; p < s.size();) {
                size_t q = s.find(',', p);
                want.push_back(std::stod(s.substr(p, q - p)));
                p = q == std::string::npos ? s.size() : q + 1;
            }
        }
        else if (prefix.empty()) prefix = a;
        else out = a;
    }
    if (prefix.empty() || out.empty()) {
        std::fprintf(stderr, "usage: pkc_to_gadget [--threads T] [--max-per-file M] [--percents 1,3] prefix out_prefix\n");
        return 2;
    }
    T = std::max(1, T);

    // selected files and their blocks
    std::vector<std::unique_ptr<FileReader>> readers;
    for (auto& p : files_of(prefix)) {
        auto r = std::make_unique<FileReader>(p);
        const Header& h = r->header();
        double pc = h.percents[h.subsample];
        bool keep = want.empty();
        for (double w : want) keep = keep || std::fabs(w - pc) < 1e-9;
        if (keep) readers.push_back(std::move(r));
    }
    if (readers.empty()) throw std::runtime_error("no file matches --percents");
    const Header& h0 = readers[0]->header();
    struct Job {
        const FileReader* r;
        size_t b;
        uint64_t n;
    };
    std::vector<std::vector<Job>> groups(1);
    uint64_t in_group = 0, total = 0;
    for (auto& r : readers)
        for (size_t b = 0; b < r->n_blocks(); ++b) {
            uint64_t n = r->block_particles(b);
            if (in_group + n > maxper && in_group) groups.emplace_back(), in_group = 0;
            groups.back().push_back({r.get(), b, n});
            in_group += n, total += n;
        }

    const double a = 1.0 / (1.0 + h0.redshift), L = h0.box_size, Q = double(h0.q());
    const double rho_crit = 27.7536627; // 1e10 Msun/h per (Mpc/h)^3
    const double mpart = rho_crit * h0.omega_m * L * L * L / double(total);
    const float vfac = float(h0.vel_step / std::sqrt(a));

    std::atomic<size_t> next{0};
    std::vector<std::thread> th;
    std::vector<std::string> err(T);
    for (int t = 0; t < T; ++t)
        th.emplace_back([&, t] {
            try {
                Block blk;
                std::vector<float> pos, vel;
                std::vector<uint64_t> ids;
                for (size_t g; (g = next.fetch_add(1)) < groups.size();) {
                    uint64_t n = 0;
                    for (auto& j : groups[g]) n += j.n;
                    pos.resize(3 * n), vel.resize(3 * n), ids.resize(n);
                    uint64_t k = 0;
                    for (auto& j : groups[g]) {
                        j.r->read_block(j.b, blk);
                        for (size_t i = 0; i < blk.size(); ++i, ++k) {
                            for (int d = 0; d < 3; ++d) {
                                pos[3 * k + d] = float((double(blk.X[3 * i + d]) + 0.5) / Q * L);
                                vel[3 * k + d] = float(blk.V[3 * i + d]) * vfac;
                            }
                            ids[k] = blk.id[i];
                        }
                    }
                    GadgetHeader gh;
                    std::memset(&gh, 0, sizeof gh);
                    gh.npart[1] = uint32_t(n);
                    gh.mass[1] = mpart;
                    gh.time = a, gh.redshift = h0.redshift;
                    gh.npartTotal[1] = uint32_t(total & 0xffffffffu);
                    gh.npartTotalHighWord[1] = uint32_t(total >> 32);
                    gh.num_files = int32_t(groups.size());
                    gh.BoxSize = L, gh.Omega0 = h0.omega_m, gh.OmegaLambda = h0.omega_l, gh.HubbleParam = h0.h;
                    std::string path = out + "." + std::to_string(g);
                    FILE* f = std::fopen(path.c_str(), "wb");
                    if (!f) throw std::runtime_error("cannot create " + path);
                    record(f, &gh, sizeof gh);
                    record(f, pos.data(), 12 * n);
                    record(f, vel.data(), 12 * n);
                    record(f, ids.data(), 8 * n);
                    if (std::fclose(f) != 0) throw std::runtime_error("close failed " + path);
                }
            } catch (std::exception& e) {
                err[t] = e.what();
            }
        });
    for (auto& x : th) x.join();
    for (auto& e : err)
        if (!e.empty()) {
            std::fprintf(stderr, "pkc_to_gadget: %s\n", e.c_str());
            return 1;
        }
    std::printf("wrote %zu Gadget-2 files %s.0 .. %s.%zu: %llu particles, z = %.4f, L = %g Mpc/h, m_p = %.6e Msun/h\n", groups.size(), out.c_str(), out.c_str(),
                groups.size() - 1, (unsigned long long)total, h0.redshift, L, mpart * 1e10);
    std::printf("Rockstar: FILE_FORMAT = \"GADGET2\"  NUM_BLOCKS = %zu  FILENAME = \"<name>.<block>\"  GADGET_ID_BYTES = 8\n"
                "          GADGET_LENGTH_CONVERSION = 1  GADGET_MASS_CONVERSION = 1e+10  GADGET_HALO_PARTICLE_TYPE = 1\n",
                groups.size());
    return 0;
}
