// pkc_convert: pkdgrav3 HDF5 snapshot -> snapcomp files (one per subsample), written the way pkdgrav3
// will write them: T threads, each owning a contiguous slice of particles (as in pkdgrav3 memory),
// each looping over subsamples and appending blocks to the shared per-subsample file.
// With --verify, every file is read back and compared with the input.
//
// usage: pkc_convert [--dx kpc/h] [--dv km/s] [--percents 1,3,5,9,10,12,15,20,25 | --percents 100]
//                    [--seed S] [--block B] [--level L] [--threads T] [--verify] in.hdf5 out_prefix
// Output files: <out_prefix>.s<NN>.0.pkc (NN = subsample percent).
#include "snapcomp.h"

#include <hdf5.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <vector>

using namespace snapcomp;

static double attr(hid_t g, const char* name) {
    hid_t a = H5Aopen(g, name, H5P_DEFAULT);
    if (a < 0) throw std::runtime_error(std::string("missing attribute ") + name);
    double v = 0;
    hid_t t = H5Aget_type(a);
    if (H5Tget_class(t) == H5T_FLOAT || H5Tget_class(t) == H5T_INTEGER) H5Aread(a, H5T_NATIVE_DOUBLE, &v);
    H5Tclose(t);
    H5Aclose(a);
    return v;
}

template <typename T> static std::vector<T> dataset(hid_t f, const char* name, hid_t type, size_t per) {
    hid_t d = H5Dopen(f, name, H5P_DEFAULT);
    if (d < 0) throw std::runtime_error(std::string("missing dataset ") + name);
    hid_t s = H5Dget_space(d);
    hsize_t dims[2] = {0, 1};
    H5Sget_simple_extent_dims(s, dims, nullptr);
    std::vector<T> v(size_t(dims[0]) * per);
    if (H5Dread(d, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0) throw std::runtime_error(std::string("read failed ") + name);
    H5Sclose(s);
    H5Dclose(d);
    return v;
}

static size_t file_size(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0 ? size_t(st.st_size) : 0;
}

int main(int argc, char** argv) {
    double dx_kpc = 1.95, dv = 2.44;
    std::vector<double> pc = {1, 3, 5, 9, 10, 12, 15, 20, 25};
    uint64_t seed = 0;
    uint32_t block = 65536;
    int level = 3, T = int(std::thread::hardware_concurrency());
    bool verify = false;
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + a);
            return argv[++i];
        };
        if (a == "--dx") dx_kpc = std::stod(next());
        else if (a == "--dv") dv = std::stod(next());
        else if (a == "--seed") seed = std::stoull(next());
        else if (a == "--block") block = uint32_t(std::stoul(next()));
        else if (a == "--level") level = std::stoi(next());
        else if (a == "--threads") T = std::stoi(next());
        else if (a == "--verify") verify = true;
        else if (a == "--percents") {
            pc.clear();
            std::string s = next();
            for (size_t p = 0; p < s.size();) {
                size_t q = s.find(',', p);
                pc.push_back(std::stod(s.substr(p, q == std::string::npos ? std::string::npos : q - p)));
                p = q == std::string::npos ? s.size() : q + 1;
            }
        }
        else args.push_back(a);
    }
    if (args.size() != 2) {
        std::fprintf(stderr, "usage: pkc_convert [--dx kpc/h] [--dv km/s] [--percents a,b,..] [--seed S] [--block B] [--level L] [--threads T] [--verify] "
                             "in.hdf5 out_prefix\n");
        return 2;
    }
    T = std::max(1, T);

    // ---- read the snapshot
    auto t0 = std::chrono::steady_clock::now();
    hid_t f = H5Fopen(args[0].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0) throw std::runtime_error("cannot open " + args[0]);
    hid_t g = H5Gopen(f, "/Header", H5P_DEFAULT);
    const double L = attr(g, "BoxSize"), z = attr(g, "Redshift"), time = attr(g, "Time");
    const double om = attr(g, "Omega0"), ol = attr(g, "OmegaLambda"), hh = attr(g, "HubbleParam");
    H5Gclose(g);
    auto xf = dataset<float>(f, "/PartType1/Coordinates", H5T_NATIVE_FLOAT, 3);
    auto vf = dataset<float>(f, "/PartType1/Velocities", H5T_NATIVE_FLOAT, 3);
    auto ids = dataset<uint64_t>(f, "/PartType1/ParticleIDs", H5T_NATIVE_UINT64, 1);
    H5Fclose(f);
    const size_t n = ids.size();
    std::vector<double> xd(xf.begin(), xf.end());
    xf.clear();
    xf.shrink_to_fit();
    const double a = 1.0 / (1.0 + z);
    const double vel_to_kms = a * L * 100.0 / std::sqrt(8.0 * M_PI / 3.0); // pkdgrav3 dx/dt -> peculiar km/s
    uint64_t N = uint64_t(std::llround(std::cbrt(double(n))));
    bool lattice = N * N * N == n;
    for (size_t i = 0; i < n && lattice; ++i) lattice = ids[i] < n;
    double tread = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // ---- header
    Header h;
    h.n_total = n;
    h.lattice_n = lattice ? N : 0;
    choose_grid(L, dx_kpc / 1000.0, n, h.ncell, h.sub);
    h.vel_step = dv;
    h.vel_to_kms = vel_to_kms;
    h.box_size = L;
    h.redshift = z;
    h.time = time;
    h.omega_m = om;
    h.omega_l = ol;
    h.h = hh;
    h.n_subsamples = uint32_t(pc.size());
    for (size_t k = 0; k < pc.size(); ++k) h.percents[k] = pc[k];
    h.subsample_seed = seed;
    h.block_size = block;
    h.zstd_level = level;
    h.info = "{\"source\":\"" + args[0] + "\",\"tool\":\"pkc_convert\"}";
    SubsampleTable table(pc.data(), int(pc.size()), seed);
    std::printf("%s: %zu particles, L=%g Mpc/h, z=%.3f, lattice N=%llu; ncell=%u sub=%u -> dx=%.4f kpc/h, dv=%g km/s, %d threads, block %u\n", args[0].c_str(),
                n, L, z, (unsigned long long)h.lattice_n, h.ncell, h.sub, L * 1000.0 / double(h.q()), dv, T, block);

    // ---- write: T threads, contiguous slices, subsample loop, shared writers
    std::vector<std::string> paths;
    std::vector<std::unique_ptr<FileWriter>> w;
    for (size_t k = 0; k < pc.size(); ++k) {
        Header hk = h;
        hk.subsample = uint32_t(k);
        char name[64];
        std::snprintf(name, sizeof name, ".s%02d.0.pkc", int(std::lround(pc[k])));
        paths.push_back(args[1] + name);
        w.emplace_back(new FileWriter(paths.back(), hk));
    }
    t0 = std::chrono::steady_clock::now();
    std::vector<double> enc_sec(T, 0);
    std::vector<std::thread> th;
    for (int t = 0; t < T; ++t)
        th.emplace_back([&, t] {
            size_t lo = n * size_t(t) / size_t(T), hi = n * size_t(t + 1) / size_t(T);
            BlockEncoder enc(h);
            std::vector<double> bp;
            std::vector<float> bv;
            std::vector<uint64_t> bi;
            bp.reserve(3 * size_t(block));
            bv.reserve(3 * size_t(block));
            bi.reserve(block);
            for (size_t k = 0; k < pc.size(); ++k) {
                auto flush = [&] {
                    if (bi.empty()) return;
                    auto e0 = std::chrono::steady_clock::now();
                    auto& out = enc.encode(bi.size(), bp.data(), bv.data(), bi.data());
                    enc_sec[t] += std::chrono::duration<double>(std::chrono::steady_clock::now() - e0).count();
                    w[k]->append(out.data(), out.size(), uint32_t(bi.size()), uint32_t(t));
                    bp.clear(), bv.clear(), bi.clear();
                };
                for (size_t i = lo; i < hi; ++i) {
                    if (table(ids[i]) != int(k)) continue;
                    for (int d = 0; d < 3; ++d) bp.push_back(xd[3 * i + d]), bv.push_back(vf[3 * i + d]);
                    bi.push_back(ids[i]);
                    if (bi.size() == block) flush();
                }
                flush();
            }
        });
    for (auto& x : th) x.join();
    for (auto& x : w) x->close();
    double twrite = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double tenc = 0;
    for (double s : enc_sec) tenc += s;
    size_t total = 0;
    std::printf("  %-8s %12s %12s %10s\n", "percent", "particles", "bytes", "B/particle");
    for (size_t k = 0; k < pc.size(); ++k) {
        FileReader r(paths[k]);
        size_t b = file_size(paths[k]);
        total += b;
        std::printf("  %-8g %12llu %12zu %10.3f\n", pc[k], (unsigned long long)r.header().n_particles, b, double(b) / double(r.header().n_particles));
    }
    std::printf("  total %zu bytes = %.3f B/particle (raw float32+uint64: 32)\n", total, double(total) / double(n));
    std::printf("  time: read %.2f s, write %.2f s wall with %d threads (%.1f M particles/s); encode %.2f thread-s = %.2f M particles/s per thread\n", tread,
                twrite, T, double(n) / twrite / 1e6, tenc, double(n) / tenc / 1e6);

    if (!verify) return 0;
    // ---- verify: every particle once, quantized values match, error bounds
    t0 = std::chrono::steady_clock::now();
    std::vector<size_t> where;
    std::vector<std::pair<uint64_t, size_t>> sorted;
    if (lattice) {
        where.assign(n, size_t(-1));
        for (size_t i = 0; i < n; ++i) where[ids[i]] = i;
    }
    else {
        sorted.resize(n);
        for (size_t i = 0; i < n; ++i) sorted[i] = {ids[i], i};
        std::sort(sorted.begin(), sorted.end());
    }
    std::vector<uint8_t> seen(n, 0);
    const double Q = double(h.q());
    double emax = 0, vmax = 0, e2 = 0, v2 = 0;
    size_t bad = 0, count = 0;
    Block b;
    std::vector<double> x, v;
    for (size_t k = 0; k < pc.size(); ++k) {
        FileReader r(paths[k]);
        for (size_t bi = 0; bi < r.n_blocks(); ++bi) {
            r.read_block(bi, b);
            b.positions(r.header(), x);
            b.velocities(r.header(), v);
            for (size_t j = 0; j < b.size(); ++j, ++count) {
                size_t i;
                if (lattice) i = where[b.id[j]];
                else {
                    auto it = std::lower_bound(sorted.begin(), sorted.end(), std::make_pair(b.id[j], size_t(0)));
                    i = (it != sorted.end() && it->first == b.id[j]) ? it->second : size_t(-1);
                }
                if (i == size_t(-1) || seen[i]++ || table(b.id[j]) != int(k)) {
                    ++bad;
                    continue;
                }
                for (int d = 0; d < 3; ++d) {
                    double dx = x[3 * j + d] - xd[3 * i + d];
                    dx -= std::floor(dx + 0.5);
                    double dvk = v[3 * j + d] - double(vf[3 * i + d]) * vel_to_kms;
                    emax = std::max(emax, std::fabs(dx)), vmax = std::max(vmax, std::fabs(dvk));
                    e2 += dx * dx, v2 += dvk * dvk;
                }
            }
        }
    }
    for (size_t i = 0; i < n; ++i) bad += seen[i] != 1;
    double tver = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("  verify: %zu particles read, %zu problems; pos max %.4f rms %.4f kpc/h (step %.4f); vel max %.4f rms %.4f km/s (step %g); %.2f s\n", count,
                bad, emax * L * 1000, std::sqrt(e2 / (3.0 * n)) * L * 1000, L * 1000 / Q, vmax, std::sqrt(v2 / (3.0 * n)), dv, tver);
    bool ok = bad == 0 && count == n && emax <= 0.5 / Q * (1 + 1e-9) && vmax <= 0.5 * dv * (1 + 1e-6) + 1e-6;
    std::printf("  %s\n", ok ? "VERIFY OK" : "VERIFY FAILED");
    return ok ? 0 : 1;
}
