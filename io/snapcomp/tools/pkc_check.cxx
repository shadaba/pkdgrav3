// pkc_check: verify a compressed snapshot (all <prefix>.s*.pkc files) without loading it into memory.
//   - every ID 0..N^3-1 present exactly once (1-bit map, N^3/8 bytes), none out of range
//   - every particle in the subsample its ID hashes to
//   - distance of each particle from its Lagrangian lattice point (lattice IDs only), in lattice cells
//   - with --ref <prefix2>: persistence. The smallest subsample of both snapshots must hold the same IDs;
//     reports how far each of those particles moved between the two snapshots.
// usage: pkc_check [--threads T] [--ref other_prefix] prefix
#include "snapcomp.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace snapcomp;

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

struct Job {
    const FileReader* r;
    size_t block;
};

// Run fn(block, thread) over all blocks of the given readers with T threads.
template <typename F> static void parallel_blocks(const std::vector<std::unique_ptr<FileReader>>& readers, int T, F fn) {
    std::vector<Job> jobs;
    for (auto& r : readers)
        for (size_t b = 0; b < r->n_blocks(); ++b) jobs.push_back({r.get(), b});
    std::atomic<size_t> next{0};
    std::vector<std::thread> th;
    for (int t = 0; t < T; ++t)
        th.emplace_back([&, t] {
            Block blk;
            for (size_t j; (j = next.fetch_add(1)) < jobs.size();) {
                jobs[j].r->read_block(jobs[j].block, blk);
                fn(*jobs[j].r, blk, t);
            }
        });
    for (auto& x : th) x.join();
}

int main(int argc, char** argv) {
    int T = int(std::thread::hardware_concurrency());
    std::string ref, prefix;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--threads" && i + 1 < argc) T = std::stoi(argv[++i]);
        else if (a == "--ref" && i + 1 < argc) ref = argv[++i];
        else prefix = a;
    }
    if (prefix.empty()) {
        std::fprintf(stderr, "usage: pkc_check [--threads T] [--ref other_prefix] prefix\n");
        return 2;
    }
    T = std::max(1, T);
    auto paths = files_of(prefix);
    std::vector<std::unique_ptr<FileReader>> readers;
    for (auto& p : paths) readers.emplace_back(new FileReader(p));
    const Header& h0 = readers[0]->header();
    const uint64_t N = h0.lattice_n, ntot = h0.n_total, Q = h0.q();
    if (N == 0 || N * N * N != ntot) throw std::runtime_error("pkc_check needs lattice IDs (lattice_n^3 == n_total)");
    SubsampleTable table(h0.percents, int(h0.n_subsamples), h0.subsample_seed);
    std::printf("%s: %zu files, N=%llu^3, z=%.4f, writers=%u, subsamples=%u, grid %ux%u (dx %.4f kpc/h), dv %g km/s, %s\n", prefix.c_str(), paths.size(),
                (unsigned long long)N, h0.redshift, h0.n_writers, h0.n_subsamples, h0.ncell, h0.sub, h0.box_size * 1e3 / double(Q), h0.vel_step,
                h0.has_potential() ? ("potential dphi " + std::to_string(h0.pot_step) + " (km/s)^2").c_str() : "no potential");

    // ---- pass 1: every ID once, membership, lattice offsets
    std::vector<std::atomic<uint64_t>> bits((ntot + 63) / 64);
    for (auto& b : bits) b.store(0, std::memory_order_relaxed);
    std::atomic<uint64_t> n_read{0}, dup{0}, out_of_range{0}, wrong_sub{0};
    const int NB = 8; // histogram of max-dimension lattice offset: <0.5,<1,<2,<4,<8,<16,<64,>=64 cells
    const double edges[NB - 1] = {0.5, 1, 2, 4, 8, 16, 64};
    std::vector<std::vector<uint64_t>> hist(T, std::vector<uint64_t>(NB, 0));
    std::vector<double> sum2(T, 0), maxoff(T, 0);
    parallel_blocks(readers, T, [&](const FileReader& r, const Block& b, int t) {
        const int k = int(r.header().subsample);
        n_read += b.size();
        for (size_t j = 0; j < b.size(); ++j) {
            uint64_t id = b.id[j];
            if (id >= ntot) {
                ++out_of_range;
                continue;
            }
            uint64_t m = 1ull << (id & 63);
            if (bits[id >> 6].fetch_or(m, std::memory_order_relaxed) & m) ++dup;
            if (table(id) != k) ++wrong_sub;
            uint64_t q[3] = {id % N, (id / N) % N, id / (N * N)};
            double mx = 0;
            for (int d = 0; d < 3; ++d) {
                double o = (double(b.X[3 * j + d]) + 0.5) / double(Q) * double(N) - (double(q[d]) + 0.5); // lattice cells
                o -= double(N) * std::floor(o / double(N) + 0.5);
                sum2[t] += o * o;
                mx = std::max(mx, std::fabs(o));
            }
            maxoff[t] = std::max(maxoff[t], mx);
            int bin = 0;
            while (bin < NB - 1 && mx >= edges[bin]) ++bin;
            hist[t][bin]++;
        }
    });
    uint64_t present = 0;
    for (auto& b : bits) present += uint64_t(__builtin_popcountll(b.load()));
    std::vector<uint64_t> H(NB, 0);
    double s2 = 0, mo = 0;
    for (int t = 0; t < T; ++t) {
        for (int k = 0; k < NB; ++k) H[k] += hist[t][k];
        s2 += sum2[t], mo = std::max(mo, maxoff[t]);
    }
    const uint64_t missing = ntot - present;
    std::printf("  IDs: %llu read, %llu distinct of %llu; duplicates %llu, missing %llu, out of range %llu; wrong subsample %llu\n",
                (unsigned long long)n_read.load(), (unsigned long long)present, (unsigned long long)ntot, (unsigned long long)dup.load(),
                (unsigned long long)missing, (unsigned long long)out_of_range.load(), (unsigned long long)wrong_sub.load());
    std::printf("  distance from own lattice point (largest of x,y,z, in lattice cells): rms %.3f (3-D), max %.2f;", std::sqrt(s2 / double(present)), mo);
    const char* lab[NB] = {"<0.5", "<1", "<2", "<4", "<8", "<16", "<64", ">=64"};
    for (int k = 0; k < NB; ++k) std::printf(" %s:%.4f%%", lab[k], 100.0 * double(H[k]) / double(present));
    std::printf("\n  (random IDs would give rms ~%.0f cells)\n", double(N) / 2);
    bool ok = dup == 0 && missing == 0 && out_of_range == 0 && wrong_sub == 0 && n_read == ntot;

    // ---- pass 2: persistence against --ref, smallest subsample
    if (!ref.empty()) {
        auto rpaths = files_of(ref);
        std::vector<std::unique_ptr<FileReader>> rr, cur;
        for (auto& p : rpaths) {
            auto r = std::make_unique<FileReader>(p);
            if (r->header().subsample == 0) rr.push_back(std::move(r));
        }
        for (auto& p : paths) {
            auto r = std::make_unique<FileReader>(p);
            if (r->header().subsample == 0) cur.push_back(std::move(r));
        }
        struct Rec {
            uint64_t id;
            int64_t X[3];
        };
        std::vector<Rec> a;
        std::mutex mtx;
        parallel_blocks(rr, T, [&](const FileReader&, const Block& b, int) {
            std::vector<Rec> loc(b.size());
            for (size_t j = 0; j < b.size(); ++j) loc[j] = {b.id[j], {b.X[3 * j], b.X[3 * j + 1], b.X[3 * j + 2]}};
            std::lock_guard<std::mutex> l(mtx);
            a.insert(a.end(), loc.begin(), loc.end());
        });
        std::sort(a.begin(), a.end(), [](const Rec& x, const Rec& y) { return x.id < y.id; });
        std::atomic<uint64_t> found{0}, notfound{0}, ncur{0};
        std::vector<std::vector<float>> disp(T);
        const double L = h0.box_size;
        parallel_blocks(cur, T, [&](const FileReader&, const Block& b, int t) {
            ncur += b.size();
            for (size_t j = 0; j < b.size(); ++j) {
                auto it = std::lower_bound(a.begin(), a.end(), b.id[j], [](const Rec& x, uint64_t v) { return x.id < v; });
                if (it == a.end() || it->id != b.id[j]) {
                    ++notfound;
                    continue;
                }
                ++found;
                double r2 = 0;
                for (int d = 0; d < 3; ++d) {
                    double dx = double(b.X[3 * j + d] - it->X[d]) / double(Q);
                    dx -= std::floor(dx + 0.5);
                    r2 += dx * dx;
                }
                disp[t].push_back(float(std::sqrt(r2) * L));
            }
        });
        std::vector<float> all;
        for (auto& d : disp) all.insert(all.end(), d.begin(), d.end());
        std::sort(all.begin(), all.end());
        auto pct = [&](double p) { return all.empty() ? 0.0 : double(all[size_t(p * double(all.size() - 1))]); };
        bool same = notfound == 0 && ncur == a.size();
        std::printf("  persistence vs %s (smallest subsample, %.0f%%): %zu IDs there, %llu here, %llu found, %llu not found -> same set: %s\n", ref.c_str(),
                    h0.percents[0], a.size(), (unsigned long long)ncur.load(), (unsigned long long)found.load(), (unsigned long long)notfound.load(),
                    same ? "yes" : "NO");
        std::printf("  same-ID displacement (Mpc/h): median %.3f, 90%% %.3f, 99%% %.3f, max %.3f; random pairs would be ~%.1f\n", pct(0.5), pct(0.9), pct(0.99),
                    pct(1.0), 0.48 * L);
        ok = ok && same;
    }
    std::printf("  %s\n", ok ? "CHECK OK" : "CHECK FAILED");
    return ok ? 0 : 1;
}
