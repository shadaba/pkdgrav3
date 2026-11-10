// Unit tests for snapcomp. usage: snapcomp_test [scratch_dir]
#include "snapcomp.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

using namespace snapcomp;

static int failures = 0;
#define CHECK(c)                                                                                                                                               \
    do {                                                                                                                                                       \
        if (!(c)) {                                                                                                                                            \
            std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c);                                                                                  \
            ++failures;                                                                                                                                        \
        }                                                                                                                                                      \
    } while (0)

static Header make_header(uint32_t ncell, uint32_t sub, double dv, uint64_t lattice) {
    Header h;
    h.ncell = ncell;
    h.sub = sub;
    h.vel_step = dv;
    h.vel_to_kms = 1.0;
    h.lattice_n = lattice;
    h.box_size = 100;
    h.zstd_level = 3;
    return h;
}

struct Parts {
    std::vector<double> pos;
    std::vector<float> vel;
    std::vector<uint64_t> id;
};

// Lattice particles displaced by a random field: a crude stand-in for a simulation.
static Parts lattice_parts(uint64_t N, double disp, double vrms, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> g(0, 1);
    Parts p;
    for (uint64_t iz = 0; iz < N; ++iz)
        for (uint64_t iy = 0; iy < N; ++iy)
            for (uint64_t ix = 0; ix < N; ++ix) {
                uint64_t q[3] = {ix, iy, iz};
                for (int d = 0; d < 3; ++d) {
                    double x = (q[d] + 0.5) / N - 0.5 + disp * g(rng) / N;
                    x -= std::floor(x + 0.5); // wrap into [-0.5, 0.5)
                    p.pos.push_back(x);
                    p.vel.push_back(float(vrms * g(rng)));
                }
                p.id.push_back((iz * N + iy) * N + ix);
            }
    // scramble order, as in pkdgrav3 memory
    std::vector<size_t> perm(p.id.size());
    for (size_t i = 0; i < perm.size(); ++i) perm[i] = i;
    std::shuffle(perm.begin(), perm.end(), rng);
    Parts s;
    for (size_t i : perm) {
        for (int d = 0; d < 3; ++d) {
            s.pos.push_back(p.pos[3 * i + d]);
            s.vel.push_back(p.vel[3 * i + d]);
        }
        s.id.push_back(p.id[i]);
    }
    return s;
}

// Encode in blocks, decode, compare with independent quantization; check error bounds.
static void roundtrip(const Header& h, const Parts& p, size_t block, const char* name, bool print = true) {
    BlockEncoder enc(h);
    const double Q = double(h.q());
    size_t n = p.id.size(), bytes = 0;
    double emax = 0, vmax = 0;
    Block b;
    std::vector<double> x, v;
    int coders[3] = {0, 0, 0};
    for (size_t s = 0; s < n; s += block) {
        size_t m = std::min(block, n - s);
        auto& out = enc.encode(m, &p.pos[3 * s], &p.vel[3 * s], &p.id[s]);
        bytes += out.size();
        decode_block(h, out.data(), out.size(), b);
        coders[b.id_coder]++;
        CHECK(b.size() == m);
        // every input particle must appear exactly once: match by ID
        std::vector<long> where(m, -1);
        std::vector<std::pair<uint64_t, size_t>> in(m), dec(m);
        for (size_t i = 0; i < m; ++i) in[i] = {p.id[s + i], i}, dec[i] = {b.id[i], i};
        std::sort(in.begin(), in.end());
        std::sort(dec.begin(), dec.end());
        b.positions(h, x);
        b.velocities(h, v);
        for (size_t i = 0; i < m; ++i) {
            CHECK(in[i].first == dec[i].first);
            size_t a = s + in[i].second, k = dec[i].second;
            for (int d = 0; d < 3; ++d) {
                int64_t X = int64_t(std::floor((p.pos[3 * a + d] + 0.5) * Q)) % int64_t(h.q());
                if (X < 0) X += int64_t(h.q());
                CHECK(b.X[3 * k + d] == X);
                double dx = x[3 * k + d] - p.pos[3 * a + d];
                dx -= std::floor(dx + 0.5);
                emax = std::max(emax, std::fabs(dx));
                vmax = std::max(vmax, std::fabs(v[3 * k + d] - p.vel[3 * a + d]));
            }
        }
    }
    CHECK(emax <= 0.5 / Q * (1 + 1e-9));
    CHECK(vmax <= 0.5 * h.vel_step * (1 + 1e-6) + 1e-6);
    if (print)
        std::printf("  %-34s n=%8zu  %6.2f B/particle  max err %.3g step(pos) %.3g step(vel)  coders raw/delta/lagr %d/%d/%d\n", name, n, double(bytes) / n,
                    emax * Q, vmax / h.vel_step, coders[0], coders[1], coders[2]);
}

static void test_subsamples() {
    struct G {
        uint64_t id, hash;
        int k0, k7;
    };
    // golden values from tools/snapcomp_proto/snapcomp.py (numpy)
    const G g[] = {{0ull, 0xe220a8397b1dcdafull, 8, 8},
                   {1ull, 0x910a2dec89025cc1ull, 7, 6},
                   {2ull, 0x975835de1c9756ceull, 7, 7},
                   {3ull, 0x1d0b14e4db018fedull, 3, 8},
                   {12345ull, 0x22118258a9d111a0ull, 3, 6},
                   {4294967295ull, 0x73b13ba2aff181c0ull, 6, 7},
                   {1099511627783ull, 0xb94083aa31091a4eull, 7, 6},
                   {8589934591ull, 0xb0caed0a7bf9c2c8ull, 7, 7},
                   {9223372036854775807ull, 0x2a67d7552e039ea7ull, 3, 6},
                   {18446744073709551615ull, 0xe4d971771b652c20ull, 8, 4}};
    const double pc[] = {1, 3, 5, 9, 10, 12, 15, 20, 25};
    SubsampleTable t0(pc, 9, 0), t7(pc, 9, 7);
    for (auto& x : g) {
        CHECK(splitmix64(x.id, 0) == x.hash);
        CHECK(t0(x.id) == x.k0);
        CHECK(t7(x.id) == x.k7);
    }
    const long want[] = {10257, 30096, 49943, 89668, 100182, 119969, 150152, 200382, 249351}; // numpy, ids 0..1e6-1
    long cnt[9] = {};
    for (uint64_t i = 0; i < 1000000; ++i) cnt[t0(i)]++;
    for (int k = 0; k < 9; ++k) CHECK(cnt[k] == want[k]);
    bool threw = false;
    try {
        const double bad[] = {50, 40};
        SubsampleTable t(bad, 2, 0);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw);
    std::printf("  subsample hash and membership match numpy (10 golden IDs, 1e6 counts)\n");
}

static void test_edge_cases() {
    Header h = make_header(16, 256, 1.0, 16);
    // single particle at the box edges
    for (double x : {-0.5, -0.5 + 1e-15, 0.0, 0.5 - 1e-12, 0.5, 0.75, -0.75}) {
        Parts p;
        p.pos = {x, -x, 0.25};
        p.vel = {1e4f, -1e4f, 0.f};
        p.id = {4095};
        roundtrip(h, p, 1, "edge", false);
    }
    // all particles in one cell, IDs not on the lattice, huge values
    Parts p;
    std::mt19937_64 rng(1);
    for (int i = 0; i < 5000; ++i) {
        for (int d = 0; d < 3; ++d) p.pos.push_back(0.01 + 1e-3 * double(rng() % 1000) / 1000), p.vel.push_back(float(int64_t(rng() % 200001) - 100000));
        p.id.push_back(rng());
    }
    roundtrip(h, p, 5000, "one cell, random 64-bit IDs");
    // lattice_n set but IDs outside the lattice -> must fall back
    Header h2 = make_header(8, 64, 2.0, 8);
    roundtrip(h2, p, 1000, "lattice_n set, IDs outside lattice");
    // tiny grid, sub = 1
    Header h3 = make_header(3, 1, 0.5, 0);
    Parts r = lattice_parts(6, 0.3, 50, 9);
    roundtrip(h3, r, 7, "ncell=3 sub=1, odd block size");
    // corrupted block must throw, not crash
    BlockEncoder enc(h);
    Parts q = lattice_parts(16, 0.5, 100, 3);
    auto out = enc.encode(q.id.size(), q.pos.data(), q.vel.data(), q.id.data());
    std::vector<uint8_t> bad(out.begin(), out.end());
    int thrown = 0;
    for (size_t k : {size_t(4), size_t(20), bad.size() / 2, bad.size() - 3}) {
        auto c = bad;
        c[k] ^= 0x5a;
        Block b;
        try {
            decode_block(h, c.data(), c.size(), b);
        } catch (const std::exception&) {
            ++thrown;
        }
    }
    CHECK(thrown == 4);
    Block b;
    bool threw = false;
    try {
        decode_block(h, out.data(), out.size() / 2, b);
    } catch (const std::exception&) {
        threw = true;
    }
    CHECK(threw);
    std::printf("  edge cases: box edges, one cell, non-lattice IDs, sub=1, corruption detected (%d/4 + truncation)\n", thrown);
}

static void test_header() {
    Header h = make_header(400, 256, 2.44, 2048);
    h.n_total = 2048ull * 2048 * 2048;
    h.redshift = 1.5;
    h.subsample = 3;
    h.n_subsamples = 9;
    const double pc[] = {1, 3, 5, 9, 10, 12, 15, 20, 25};
    for (int i = 0; i < 9; ++i) h.percents[i] = pc[i];
    h.subsample_seed = 42;
    h.info = "{\"code\":\"pkdgrav3\"}";
    h.step = 100;
    std::vector<uint8_t> b(HEADER_SIZE);
    h.pack(b.data());
    Header u = Header::unpack(b.data());
    CHECK(u.ncell == 400 && u.sub == 256 && u.vel_step == 2.44 && u.lattice_n == 2048 && u.n_total == h.n_total);
    CHECK(u.redshift == 1.5 && u.subsample == 3 && u.n_subsamples == 9 && u.percents[8] == 25 && u.subsample_seed == 42);
    CHECK(u.info == h.info && u.step == 100 && u.block_size == 65536 && u.zstd_level == 3);
    b[0] = 'X';
    bool threw = false;
    try {
        Header::unpack(b.data());
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw);
    uint32_t nc, sub;
    choose_grid(200.0, 0.00195, 2048ull * 2048 * 2048, nc, sub);
    CHECK(nc == 2048 && sub == 50);
    choose_grid(12.5, 0.00195, 128ull * 128 * 128, nc, sub);
    CHECK(nc == 128 && sub == 50);
    choose_grid(1000.0, 10.0, 1000, nc, sub); // cells of the mean spacing (100), 10 steps each
    CHECK(nc == 10 && sub == 10);
    choose_grid(1000.0, 400.0, 1000, nc, sub); // step much larger than a mean-spacing cell -> coarser cells
    CHECK(nc == 3 && sub == 1);
    std::printf("  header pack/unpack, bad magic, choose_grid\n");
}

// Many threads append to one file; the reader must return every particle exactly once.
// Optional potential: version 2 header, exact quantization, error <= step/2, smaller than raw floats;
// files without it stay version 1 and the encoder refuses a potential the header does not announce.
static void test_potential(const std::string& dir) {
    Header h = make_header(64, 50, 2.44, 64);
    Parts p = lattice_parts(64, 2.0, 300, 7);
    const size_t n = p.id.size();
    std::vector<float> pot(n); // smooth field + small noise, in a writer unit of 10 (km/s)^2
    std::mt19937_64 rng(3);
    std::normal_distribution<double> g(0, 1);
    for (size_t i = 0; i < n; ++i) {
        const double* x = &p.pos[3 * i];
        pot[i] = float(-3e4 * (1 + std::cos(6.283 * x[0]) * std::cos(6.283 * x[1]) * std::cos(6.283 * x[2])) + 50 * g(rng));
    }
    std::vector<uint8_t> hb(HEADER_SIZE);
    h.pack(hb.data());
    CHECK(hb[8] == 1 && !Header::unpack(hb.data()).has_potential());
    bool threw = false;
    try {
        BlockEncoder e0(h);
        e0.encode(10, p.pos.data(), p.vel.data(), p.id.data(), pot.data());
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw);

    h.pot_step = 100;
    h.pot_to_kms2 = 10;
    h.info = std::string(5000, 'x'); // must not run into the potential fields
    h.pack(hb.data());
    Header u = Header::unpack(hb.data());
    CHECK(hb[8] == 2 && u.pot_step == 100 && u.pot_to_kms2 == 10 && u.info.size() == 4080 - 288 - 1);
    h.info.clear();

    const std::string path = dir + "/snapcomp_pot.pkc";
    BlockEncoder enc(h);
    size_t bytes = 0, pbytes = 0;
    {
        FileWriter w(path, h);
        for (size_t s = 0; s < n; s += 65536) {
            size_t m = std::min<size_t>(65536, n - s);
            auto& out = enc.encode(m, &p.pos[3 * s], &p.vel[3 * s], &p.id[s], &pot[s]);
            w.append(out.data(), out.size(), uint32_t(m));
            bytes += out.size();
            pbytes += enc.stats.pot;
        }
        w.close();
    }
    FileReader r(path);
    CHECK(r.header().has_potential() && r.header().pot_step == 100);
    Block b;
    std::vector<double> ph;
    double emax = 0;
    bool exact = true;
    size_t seen = 0;
    for (size_t i = 0; i < r.n_blocks(); ++i) {
        r.read_block(i, b);
        b.potentials(r.header(), ph);
        CHECK(ph.size() == b.size());
        seen += b.size();
    }
    // the input is shuffled, so map ID -> input index
    std::vector<size_t> at(n);
    for (size_t i = 0; i < n; ++i) at[p.id[i]] = i;
    for (size_t i = 0; i < r.n_blocks(); ++i) {
        r.read_block(i, b);
        b.potentials(r.header(), ph);
        for (size_t k = 0; k < b.size(); ++k) {
            const double truth = double(pot[at[b.id[k]]]) * 10;
            exact = exact && b.P[k] == std::llrint(truth / 100);
            emax = std::max(emax, std::fabs(ph[k] - truth));
        }
    }
    CHECK(seen == n && exact && emax <= 50 * (1 + 1e-9));
    CHECK(double(pbytes) / n < 2.0);
    unlink(path.c_str());
    std::printf("  potential: version 2 header, %zu particles, max err %.1f (km/s)^2 for step 100, %.2f B/particle for the potential, %.2f in all\n", n, emax,
                double(pbytes) / n, double(bytes) / n);
}

static void test_file_threads(const std::string& dir) {
    Header h = make_header(32, 64, 1.0, 32);
    h.n_total = 32 * 32 * 32;
    Parts p = lattice_parts(32, 1.0, 200, 11);
    const std::string path = dir + "/snapcomp_test.pkc";
    const int T = 8;
    const size_t block = 1000;
    {
        FileWriter w(path, h);
        std::vector<std::thread> th;
        for (int t = 0; t < T; ++t)
            th.emplace_back([&, t] {
                BlockEncoder enc(h);
                for (size_t s = size_t(t) * block; s < p.id.size(); s += T * block) {
                    size_t m = std::min(block, p.id.size() - s);
                    auto& out = enc.encode(m, &p.pos[3 * s], &p.vel[3 * s], &p.id[s]);
                    w.append(out.data(), out.size(), uint32_t(m), uint32_t(t));
                }
            });
        for (auto& x : th) x.join();
        w.close();
    }
    FileReader r(path);
    CHECK(r.header().n_particles == p.id.size());
    std::vector<int> seen(p.id.size(), 0);
    Block b;
    for (size_t i = 0; i < r.n_blocks(); ++i) {
        r.read_block(i, b);
        for (size_t k = 0; k < b.size(); ++k) seen[b.id[k]]++;
    }
    bool once = true;
    for (int s : seen) once = once && s == 1;
    CHECK(once);
    // truncated file (no trailer) must be rejected
    if (truncate(path.c_str(), off_t(HEADER_SIZE + 100)) == 0) {
        bool threw = false;
        try {
            FileReader bad(path);
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(threw);
    }
    unlink(path.c_str());
    std::printf("  file: %d threads, %zu blocks, every particle read back once; truncated file rejected\n", T, r.n_blocks());
}

static void test_speed() {
    Header h = make_header(64, 50, 2.44, 64);
    Parts p = lattice_parts(64, 2.0, 300, 5); // 262144 particles
    BlockEncoder enc(h);
    const size_t block = 65536;
    auto t0 = std::chrono::steady_clock::now();
    size_t bytes = 0;
    for (int rep = 0; rep < 4; ++rep)
        for (size_t s = 0; s < p.id.size(); s += block) bytes += enc.encode(block, &p.pos[3 * s], &p.vel[3 * s], &p.id[s]).size();
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double mps = 4.0 * p.id.size() / sec / 1e6;
    Block b;
    auto& out = enc.encode(block, p.pos.data(), p.vel.data(), p.id.data());
    t0 = std::chrono::steady_clock::now();
    for (int rep = 0; rep < 16; ++rep) decode_block(h, out.data(), out.size(), b);
    double dsec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("  speed (one thread, block 65536, zstd-3): encode %.1f M particles/s, decode %.1f M particles/s, %.2f B/particle\n", mps,
                16.0 * block / dsec / 1e6, double(bytes) / (4.0 * p.id.size()));
}

static void test_cosmology_header(const std::string& dir) {
    // version 3: binary cosmology block at 4000, with and without the potential; info ends before 4000
    Header h = make_header(32, 25, 2.44, 32);
    h.has_cosmology = true;
    h.omega_b = 0.049, h.omega_rad = 9e-5, h.omega_de = 0.69, h.w0 = -0.9, h.wa = 0.1;
    h.sigma8 = 0.81, h.n_s = 0.965, h.running = 0.001, h.normalization = 0, h.pivot = 0.05;
    h.info = std::string(5000, 'y');
    std::vector<uint8_t> hb(HEADER_SIZE);
    h.pack(hb.data());
    Header u = Header::unpack(hb.data());
    CHECK(hb[8] == 3 && u.has_cosmology && !u.has_potential() && u.info.size() == 4000 - 288 - 1);
    CHECK(u.omega_b == 0.049 && u.omega_rad == 9e-5 && u.omega_de == 0.69 && u.w0 == -0.9 && u.wa == 0.1);
    CHECK(u.sigma8 == 0.81 && u.n_s == 0.965 && u.running == 0.001 && u.normalization == 0 && u.pivot == 0.05);
    h.pot_step = 50, h.pot_to_kms2 = 1;
    h.info = "{}";
    h.pack(hb.data());
    u = Header::unpack(hb.data());
    CHECK(hb[8] == 3 && u.has_potential() && u.pot_step == 50 && u.w0 == -0.9 && u.info == "{}");
    // a version-3 file round trip (no potential)
    h.pot_step = 0;
    Parts p = lattice_parts(32, 1.0, 200, 11);
    const std::string path = dir + "/snapcomp_v3.pkc";
    {
        BlockEncoder enc(h);
        FileWriter w(path, h);
        auto& out = enc.encode(p.id.size(), p.pos.data(), p.vel.data(), p.id.data());
        w.append(out.data(), out.size(), uint32_t(p.id.size()));
        w.close();
    }
    FileReader r(path);
    Block b;
    r.read_block(0, b);
    CHECK(r.header().has_cosmology && r.header().wa == 0.1 && b.size() == p.id.size() && b.P.empty());
    std::remove(path.c_str());
    std::printf("  cosmology header (version 3): ok\n");
}

int main(int argc, char** argv) {
    std::string dir = argc > 1 ? argv[1] : "/tmp";
    std::printf("snapcomp tests\n");
    test_subsamples();
    test_header();
    test_edge_cases();
    Header h = make_header(64, 50, 2.44, 64);
    roundtrip(h, lattice_parts(64, 0.3, 100, 1), 65536, "lattice, small displacement");
    roundtrip(h, lattice_parts(64, 4.0, 500, 2), 65536, "lattice, large displacement");
    roundtrip(h, lattice_parts(64, 4.0, 500, 2), 1000, "lattice, block 1000");
    test_file_threads(dir);
    test_potential(dir);
    test_cosmology_header(dir);
    test_speed();
    std::printf(failures ? "FAILED: %d checks\n" : "all tests passed\n", failures);
    return failures ? 1 : 0;
}
