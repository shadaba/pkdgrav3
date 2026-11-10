// pkc_to_gadget: compressed snapshot (<prefix>.s*.pkc) -> Gadget-2 binary files (SnapFormat 1) for
// Rockstar (FILE_FORMAT = "GADGET2"). Used to validate Rockstar's native PKC reader (sahyadri-sandbox).
//
// usage: pkc_to_gadget [--threads T] [--max-per-file M | --blocks K] [--percents 1,3]
//                      [--from-hdf5 snap.hdf5] prefix out_prefix
//   writes out_prefix.0 .. out_prefix.<K-1> and out_prefix.rockstar.cfg (Rockstar settings incl. W0/WA).
//   --blocks K      K files; file k holds exactly the particles, in the same order, that Rockstar's native
//                   PKC reader gives block k with NUM_BLOCKS = K (files ordered by subsample and writer,
//                   blocks in index order; pkc block i goes to file floor(start_i * K / total)).
//   --from-hdf5 F   particles from the pkdgrav3 HDF5 snapshot F (unquantized, in file order) instead of the
//                   .pkc blocks; box, redshift and cosmology still come from the .pkc header of prefix (the
//                   same output step), so both routes give Rockstar identical settings. Needs HDF5.
//   Cosmology (format version 3): OmegaLambda = omega_l + omega_de (w0-wa runs store their dark energy in
//   omega_de); W0/WA go to the .rockstar.cfg file because a Gadget header cannot hold them.
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
#ifdef SNAPCOMP_HDF5
    #include <hdf5.h>
#endif

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

#ifdef SNAPCOMP_HDF5
static double h5attr(hid_t g, const char* name) {
    hid_t a = H5Aopen(g, name, H5P_DEFAULT);
    if (a < 0) throw std::runtime_error(std::string("HDF5: missing attribute ") + name);
    double v = 0;
    H5Aread(a, H5T_NATIVE_DOUBLE, &v);
    H5Aclose(a);
    return v;
}
template <typename T> static std::vector<T> h5data(hid_t f, const char* name, hid_t type, size_t per) {
    hid_t d = H5Dopen(f, name, H5P_DEFAULT);
    if (d < 0) throw std::runtime_error(std::string("HDF5: missing dataset ") + name);
    hid_t sp = H5Dget_space(d);
    hsize_t dims[2] = {0, 1};
    H5Sget_simple_extent_dims(sp, dims, nullptr);
    std::vector<T> v(size_t(dims[0]) * per);
    if (H5Dread(d, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0) throw std::runtime_error(std::string("HDF5: read failed ") + name);
    H5Sclose(sp);
    H5Dclose(d);
    return v;
}
#endif

int main(int argc, char** argv) {
    try {
        int T = int(std::thread::hardware_concurrency());
        uint64_t maxper = 1ull << 24;
        int64_t K = 0;
        std::vector<double> want;
        std::string prefix, out, h5;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--threads" && i + 1 < argc) T = std::stoi(argv[++i]);
            else if (a == "--max-per-file" && i + 1 < argc) maxper = std::stoull(argv[++i]);
            else if (a == "--blocks" && i + 1 < argc) K = std::stoll(argv[++i]);
            else if (a == "--from-hdf5" && i + 1 < argc) h5 = argv[++i];
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
            std::fprintf(stderr,
                         "usage: pkc_to_gadget [--threads T] [--max-per-file M | --blocks K] [--percents 1,3] [--from-hdf5 snap.hdf5] prefix out_prefix\n");
            return 2;
        }
        T = std::max(1, T);
#ifndef SNAPCOMP_HDF5
        if (!h5.empty()) throw std::runtime_error("--from-hdf5: this build has no HDF5");
#endif
        if (!h5.empty() && !want.empty()) throw std::runtime_error("--from-hdf5 reads every particle; do not combine with --percents");

        // selected files, in the native reader's order: by subsample, then writer
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
        std::sort(readers.begin(), readers.end(), [](const std::unique_ptr<FileReader>& x, const std::unique_ptr<FileReader>& y) {
            const Header &a = x->header(), &b = y->header();
            return a.subsample != b.subsample ? a.subsample < b.subsample : a.writer < b.writer;
        });
        const Header& h0 = readers[0]->header();

        // cosmology for Rockstar
        if (!h0.has_cosmology)
            std::fprintf(stderr, "pkc_to_gadget: warning: format version %s without a cosmology block: OmegaLambda = omega_l, W0/WA unknown (-1, 0 written)\n",
                         h0.has_potential() ? "2" : "1");
        if (h0.omega_l > 0 && h0.omega_de > 0 && (h0.w0 != -1.0 || h0.wa != 0.0))
            throw std::runtime_error("run has both a cosmological constant and w0-wa dark energy; Rockstar has one dark-energy component");
        const double OmegaLambda = h0.omega_l + h0.omega_de;
        const double W0 = h0.omega_de > 0 ? h0.w0 : -1.0, WA = h0.omega_de > 0 ? h0.wa : 0.0;

        struct Job {
            const FileReader* r;
            size_t b;
            uint64_t n;
        };
        std::vector<Job> jobs;
        uint64_t total = 0;
        for (auto& r : readers)
            for (size_t b = 0; b < r->n_blocks(); ++b) {
                jobs.push_back({r.get(), b, r->block_particles(b)});
                total += r->block_particles(b);
            }

        // particles from HDF5 (route "exact snapshot"): checked against the .pkc header of the same step
        std::vector<float> hx, hv;
        std::vector<uint64_t> hid;
#ifdef SNAPCOMP_HDF5
        if (!h5.empty()) {
            hid_t f = H5Fopen(h5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (f < 0) throw std::runtime_error("cannot open " + h5);
            hid_t g = H5Gopen(f, "/Header", H5P_DEFAULT);
            const double L5 = h5attr(g, "BoxSize"), z5 = h5attr(g, "Redshift");
            H5Gclose(g);
            auto x = h5data<double>(f, "/PartType1/Coordinates", H5T_NATIVE_DOUBLE, 3);
            auto v = h5data<double>(f, "/PartType1/Velocities", H5T_NATIVE_DOUBLE, 3);
            hid = h5data<uint64_t>(f, "/PartType1/ParticleIDs", H5T_NATIVE_UINT64, 1);
            H5Fclose(f);
            if (hid.size() != h0.n_total || L5 != h0.box_size || std::fabs(z5 - h0.redshift) > 1e-6 * (1 + z5))
                throw std::runtime_error("HDF5 file and .pkc header disagree (particle count, box or redshift): not the same output step?");
            const double a = 1.0 / (1.0 + h0.redshift), L = h0.box_size;
            const double vel_to_kms = a * L * 100.0 / std::sqrt(8.0 * M_PI / 3.0); // pkdgrav3 dx/dt -> peculiar km/s (as pkc_convert)
            const double vg = vel_to_kms / std::sqrt(a);                           // Gadget velocity unit
            hx.resize(x.size()), hv.resize(v.size());
            for (size_t i = 0; i < x.size(); ++i) {
                double u = x[i] - std::floor(x[i] + 0.5); // box units in [-0.5, 0.5)
                hx[i] = float((u + 0.5) * L);
                if (hx[i] >= float(L)) hx[i] -= float(L);
                hv[i] = float(v[i] * vg);
            }
            total = hid.size();
            jobs.clear();
        }
#endif

        // groups = output files
        std::vector<std::vector<Job>> groups;
        std::vector<std::pair<uint64_t, uint64_t>> ranges; // HDF5 mode: [begin, end) per file
        if (!h5.empty()) {
            const int64_t k = K > 0 ? K : int64_t((total + maxper - 1) / maxper);
            for (int64_t g = 0; g < k; ++g) ranges.push_back({total * uint64_t(g) / uint64_t(k), total * uint64_t(g + 1) / uint64_t(k)});
        }
        else if (K > 0) {
            if (uint64_t(K) > jobs.size()) throw std::runtime_error("--blocks larger than the number of pkc blocks");
            groups.resize(size_t(K));
            uint64_t start = 0;
            for (auto& j : jobs) {
                groups[size_t(start * uint64_t(K) / total)].push_back(j);
                start += j.n;
            }
        }
        else {
            groups.resize(1);
            uint64_t in_group = 0;
            for (auto& j : jobs) {
                if (in_group + j.n > maxper && in_group) groups.emplace_back(), in_group = 0;
                groups.back().push_back(j);
                in_group += j.n;
            }
        }
        const size_t nfiles = h5.empty() ? groups.size() : ranges.size();

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
                    for (size_t g; (g = next.fetch_add(1)) < nfiles;) {
                        uint64_t n = 0;
                        if (h5.empty()) {
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
                        }
                        else {
                            const uint64_t b0 = ranges[g].first;
                            n = ranges[g].second - b0;
                            pos.assign(hx.begin() + 3 * b0, hx.begin() + 3 * (b0 + n));
                            vel.assign(hv.begin() + 3 * b0, hv.begin() + 3 * (b0 + n));
                            ids.assign(hid.begin() + b0, hid.begin() + b0 + n);
                        }
                        GadgetHeader gh;
                        std::memset(&gh, 0, sizeof gh);
                        gh.npart[1] = uint32_t(n);
                        gh.mass[1] = mpart;
                        gh.time = a, gh.redshift = h0.redshift;
                        gh.npartTotal[1] = uint32_t(total & 0xffffffffu);
                        gh.npartTotalHighWord[1] = uint32_t(total >> 32);
                        gh.num_files = int32_t(nfiles);
                        gh.BoxSize = L, gh.Omega0 = h0.omega_m, gh.OmegaLambda = OmegaLambda, gh.HubbleParam = h0.h;
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
            if (!e.empty()) throw std::runtime_error(e);

        // Rockstar settings, incl. what a Gadget header cannot carry
        const std::string cfg = out + ".rockstar.cfg";
        FILE* c = std::fopen(cfg.c_str(), "w");
        if (!c) throw std::runtime_error("cannot create " + cfg);
        std::fprintf(c,
                     "# written by pkc_to_gadget from %s%s%s; append to the Rockstar config\n"
                     "FILE_FORMAT = \"GADGET2\"\nNUM_BLOCKS = %zu\nGADGET_ID_BYTES = 8\nGADGET_LENGTH_CONVERSION = 1\n"
                     "GADGET_MASS_CONVERSION = 1e+10\nGADGET_HALO_PARTICLE_TYPE = 1\nW0 = %.17g\nWA = %.17g\n"
                     "# from the Gadget header (Rockstar reads them there): Om = %.17g, Ol = %.17g, h0 = %.17g\n",
                     prefix.c_str(), h5.empty() ? "" : " and ", h5.c_str(), nfiles, W0, WA, h0.omega_m, OmegaLambda, h0.h);
        if (std::fclose(c) != 0) throw std::runtime_error("close failed " + cfg);

        std::printf("wrote %zu Gadget-2 files %s.0 .. %s.%zu: %llu particles%s, z = %.4f, L = %g Mpc/h, m_p = %.6e Msun/h\n", nfiles, out.c_str(), out.c_str(),
                    nfiles - 1, (unsigned long long)total, h5.empty() ? "" : " (from HDF5)", h0.redshift, L, mpart * 1e10);
        std::printf("Om = %g, Ol = %g (omega_l %g + omega_de %g), h = %g, W0 = %g, WA = %g\nRockstar settings: %s\n", h0.omega_m, OmegaLambda, h0.omega_l,
                    h0.omega_de, h0.h, W0, WA, cfg.c_str());
        return 0;
    } catch (std::exception& e) {
        std::fprintf(stderr, "pkc_to_gadget: %s\n", e.what());
        return 1;
    }
}
