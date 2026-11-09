/*  snapcomp implementation. See snapcomp.h and docs/snapcomp_format.md.
 *  This file is part of pkdgrav3 (GPLv3).
 */
#include "snapcomp.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <sys/stat.h>
#include <system_error>
#include <unistd.h>
#include <zstd.h>

namespace snapcomp {

namespace {

const char HEADER_MAGIC[8] = {'P', 'K', 'C', 'O', 'M', 'P', 0, 1};
const char TRAILER_MAGIC[8] = {'P', 'K', 'C', 'I', 'D', 'X', 0, 1};
const char BLOCK_MAGIC[4] = {'P', 'K', 'B', '1'};
constexpr size_t BLOCK_HEADER = 16, STREAM_ENTRY = 12, INDEX_ENTRY = 32, TRAILER = 32, INFO_OFFSET = 288;

enum StreamKind : uint8_t { K_CELLKEY = 1, K_COUNT = 2, K_OFFSET = 3, K_VEL = 4, K_ID_RAW = 5, K_ID_DELTA = 6, K_ID_LAGR = 7 };

template <typename T> void put(uint8_t* b, size_t off, T v) {
    std::memcpy(b + off, &v, sizeof(T));
}
template <typename T> T get(const uint8_t* b, size_t off) {
    T v;
    std::memcpy(&v, b + off, sizeof(T));
    return v;
}

inline uint64_t zigzag(int64_t s) {
    return (uint64_t(s) << 1) ^ uint64_t(s >> 63);
}
inline int64_t unzigzag(uint64_t u) {
    return int64_t(u >> 1) ^ -int64_t(u & 1);
}

inline int64_t floor_mod(int64_t a, int64_t m) {
    int64_t r = a % m;
    return r < 0 ? r + m : r;
}

void check(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(std::string("snapcomp: ") + what);
}

int width_for(uint64_t maxv) {
    return maxv < (1ull << 8) ? 1 : maxv < (1ull << 16) ? 2 : maxv < (1ull << 32) ? 4 : 8;
}

} // namespace

// ---------------------------------------------------------------------------------------------
// Header
// ---------------------------------------------------------------------------------------------
void Header::pack(uint8_t* b) const {
    std::memset(b, 0, HEADER_SIZE);
    std::memcpy(b, HEADER_MAGIC, 8);
    put<uint32_t>(b, 8, FORMAT_VERSION);
    put<uint32_t>(b, 12, HEADER_SIZE);
    put<uint64_t>(b, 16, n_particles);
    put<uint64_t>(b, 24, n_blocks);
    put<uint64_t>(b, 32, lattice_n);
    put<uint64_t>(b, 40, n_total);
    put<uint32_t>(b, 48, ncell);
    put<uint32_t>(b, 52, sub);
    put<double>(b, 56, vel_step);
    put<double>(b, 64, vel_to_kms);
    put<double>(b, 72, box_size);
    put<double>(b, 80, redshift);
    put<double>(b, 88, time);
    put<double>(b, 96, omega_m);
    put<double>(b, 104, omega_l);
    put<double>(b, 112, h);
    put<uint32_t>(b, 120, subsample);
    put<uint32_t>(b, 124, n_subsamples);
    for (int i = 0; i < MAX_SUBSAMPLES; ++i) put<double>(b, 128 + 8 * i, percents[i]);
    put<uint64_t>(b, 256, subsample_seed);
    put<uint32_t>(b, 264, writer);
    put<uint32_t>(b, 268, n_writers);
    put<uint32_t>(b, 272, block_size);
    put<int32_t>(b, 276, zstd_level);
    put<uint64_t>(b, 280, step);
    size_t n = std::min(info.size(), size_t(HEADER_SIZE - INFO_OFFSET - 1));
    std::memcpy(b + INFO_OFFSET, info.data(), n);
}

Header Header::unpack(const uint8_t* b) {
    check(std::memcmp(b, HEADER_MAGIC, 8) == 0, "not a snapcomp file (bad header magic)");
    check(get<uint32_t>(b, 8) == FORMAT_VERSION, "unsupported format version");
    check(get<uint32_t>(b, 12) == HEADER_SIZE, "unexpected header size");
    Header h;
    h.n_particles = get<uint64_t>(b, 16);
    h.n_blocks = get<uint64_t>(b, 24);
    h.lattice_n = get<uint64_t>(b, 32);
    h.n_total = get<uint64_t>(b, 40);
    h.ncell = get<uint32_t>(b, 48);
    h.sub = get<uint32_t>(b, 52);
    h.vel_step = get<double>(b, 56);
    h.vel_to_kms = get<double>(b, 64);
    h.box_size = get<double>(b, 72);
    h.redshift = get<double>(b, 80);
    h.time = get<double>(b, 88);
    h.omega_m = get<double>(b, 96);
    h.omega_l = get<double>(b, 104);
    h.h = get<double>(b, 112);
    h.subsample = get<uint32_t>(b, 120);
    h.n_subsamples = get<uint32_t>(b, 124);
    for (int i = 0; i < MAX_SUBSAMPLES; ++i) h.percents[i] = get<double>(b, 128 + 8 * i);
    h.subsample_seed = get<uint64_t>(b, 256);
    h.writer = get<uint32_t>(b, 264);
    h.n_writers = get<uint32_t>(b, 268);
    h.block_size = get<uint32_t>(b, 272);
    h.zstd_level = get<int32_t>(b, 276);
    h.step = get<uint64_t>(b, 280);
    const char* s = reinterpret_cast<const char*>(b + INFO_OFFSET);
    h.info.assign(s, strnlen(s, HEADER_SIZE - INFO_OFFSET));
    check(h.ncell > 0 && h.sub > 0 && h.vel_step > 0, "invalid grid or velocity step in header");
    return h;
}

void choose_grid(double box_size, double pos_step, uint64_t n_total, uint32_t& ncell, uint32_t& sub) {
    check(box_size > 0 && pos_step > 0, "box size and position step must be positive");
    double c = std::cbrt(double(n_total));
    ncell = uint32_t(std::max(1.0, std::round(c)));
    ncell = std::min<uint32_t>(ncell, 1u << 20); // keys < 2^60
    // never make a cell smaller than the step
    double s = std::round(box_size / (double(ncell) * pos_step));
    if (s < 1) {
        ncell = uint32_t(std::max(1.0, std::round(box_size / pos_step)));
        s = 1;
    }
    sub = uint32_t(std::min(s, 4294967295.0 / ncell));
}

// ---------------------------------------------------------------------------------------------
// Subsamples
// ---------------------------------------------------------------------------------------------
uint64_t splitmix64(uint64_t x, uint64_t seed) {
    uint64_t z = x + 0x9E3779B97F4A7C15ull * (seed + 1);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

SubsampleTable::SubsampleTable(const double* percents, int n_, uint64_t seed_) : n(n_), seed(seed_) {
    check(n >= 1 && n <= MAX_SUBSAMPLES, "number of subsamples must be 1..16");
    double c = 0;
    for (int i = 0; i < n; ++i) {
        check(percents[i] > 0, "subsample percents must be positive");
        c += percents[i];
        edges[i] = c / 100.0; // same arithmetic as numpy cumsum(...)/100
    }
    check(std::fabs(c - 100.0) < 1e-9, "subsample percents must sum to 100");
}

int SubsampleTable::operator()(uint64_t id) const {
    double u = double(splitmix64(id, seed) >> 11) / 9007199254740992.0; // 2^53
    int k = 0;
    while (k < n && edges[k] <= u) ++k; // searchsorted(side='right')
    return k < n ? k : n - 1;           // u < 1 <= edges[n-1] up to rounding
}

// ---------------------------------------------------------------------------------------------
// Encoder
// ---------------------------------------------------------------------------------------------
struct BlockEncoder::Impl {
    Header h;
    ZSTD_CCtx* cctx = nullptr;
    struct Rec {
        uint64_t key, id;
        uint32_t i;
    };
    std::vector<Rec> rec;
    std::vector<int64_t> X, V;   // 3n, unsorted
    std::vector<uint64_t> tmp;   // stream values
    std::vector<uint8_t> planes; // byte planes of one stream
    struct Stream {
        uint8_t width, kind;
        uint32_t count;
        std::vector<uint8_t> data; // compressed
    };
    std::vector<Stream> streams, idbest, idtry;
    std::vector<uint8_t> out;

    void make_stream(std::vector<Stream>& dst, uint8_t kind, const std::vector<uint64_t>& v) {
        uint64_t m = 0;
        for (auto x : v) m = std::max(m, x);
        int w = width_for(m);
        size_t n = v.size();
        planes.resize(n * w);
        for (int b = 0; b < w; ++b) {
            uint8_t* p = planes.data() + size_t(b) * n;
            const int sh = 8 * b;
            for (size_t i = 0; i < n; ++i) p[i] = uint8_t(v[i] >> sh);
        }
        Stream s;
        s.width = uint8_t(w);
        s.kind = kind;
        s.count = uint32_t(n);
        s.data.resize(ZSTD_compressBound(planes.size()));
        size_t c = ZSTD_compress2(cctx, s.data.data(), s.data.size(), planes.data(), planes.size());
        check(!ZSTD_isError(c), ZSTD_isError(c) ? ZSTD_getErrorName(c) : "");
        s.data.resize(c);
        dst.push_back(std::move(s));
    }
};

BlockEncoder::BlockEncoder(const Header& h) : p(new Impl) {
    check(h.ncell > 0 && h.sub > 0 && h.vel_step > 0, "encoder needs ncell, sub and vel_step");
    p->h = h;
    p->cctx = ZSTD_createCCtx();
    check(p->cctx, "ZSTD_createCCtx failed");
    ZSTD_CCtx_setParameter(p->cctx, ZSTD_c_compressionLevel, h.zstd_level);
    ZSTD_CCtx_setParameter(p->cctx, ZSTD_c_checksumFlag, 1);
}

size_t BlockEncoder::scratch_bytes() const {
    const Impl& s = *p;
    size_t b = s.rec.capacity() * sizeof(Impl::Rec) + (s.X.capacity() + s.V.capacity()) * sizeof(int64_t) + s.tmp.capacity() * sizeof(uint64_t) +
               s.planes.capacity() + s.out.capacity();
    for (const auto* v : {&s.streams, &s.idbest, &s.idtry})
        for (const auto& x : *v) b += x.data.capacity();
    return b + ZSTD_sizeof_CCtx(s.cctx);
}

BlockEncoder::~BlockEncoder() {
    ZSTD_freeCCtx(p->cctx);
    delete p;
}

const std::vector<uint8_t>& BlockEncoder::encode(size_t n, const double* pos, const float* vel, const uint64_t* ids) {
    check(n >= 1 && n <= 0xffffffffu, "block must hold 1 .. 2^32-1 particles");
    Impl& s = *p;
    const Header& h = s.h;
    const int64_t Q = int64_t(h.q()), sub = h.sub, nc = h.ncell;
    const double vfac = h.vel_to_kms / h.vel_step;
    s.X.resize(3 * n);
    s.V.resize(3 * n);
    s.rec.resize(n);
    for (size_t i = 0; i < n; ++i) {
        int64_t c[3];
        for (int d = 0; d < 3; ++d) {
            int64_t X = floor_mod(int64_t(std::floor((pos[3 * i + d] + 0.5) * double(Q))), Q);
            s.X[3 * i + d] = X;
            c[d] = X / sub;
            s.V[3 * i + d] = std::llrint(double(vel[3 * i + d]) * vfac);
        }
        s.rec[i] = {uint64_t((c[2] * nc + c[1]) * nc + c[0]), ids[i], uint32_t(i)};
    }
    std::sort(s.rec.begin(), s.rec.end(), [](const Impl::Rec& a, const Impl::Rec& b) { return a.key < b.key || (a.key == b.key && a.id < b.id); });

    s.streams.clear();
    auto& v = s.tmp;
    // cells: keys (first value, then differences) and counts
    v.clear();
    std::vector<uint64_t> counts;
    uint64_t prev = 0;
    for (size_t j = 0; j < n;) {
        size_t k = j;
        while (k < n && s.rec[k].key == s.rec[j].key) ++k;
        v.push_back(s.rec[j].key - prev);
        prev = s.rec[j].key;
        counts.push_back(k - j);
        j = k;
    }
    const uint32_t n_cells = uint32_t(v.size());
    s.make_stream(s.streams, K_CELLKEY, v);
    s.make_stream(s.streams, K_COUNT, counts);
    // offsets in the cell
    for (int d = 0; d < 3; ++d) {
        v.resize(n);
        for (size_t j = 0; j < n; ++j) v[j] = uint64_t(s.X[3 * s.rec[j].i + d] % sub);
        s.make_stream(s.streams, K_OFFSET, v);
    }
    // velocities
    for (int d = 0; d < 3; ++d) {
        v.resize(n);
        for (size_t j = 0; j < n; ++j) v[j] = zigzag(s.V[3 * s.rec[j].i + d]);
        s.make_stream(s.streams, K_VEL, v);
    }
    // IDs: try every applicable coder, keep the smallest
    auto total = [](const std::vector<Impl::Stream>& ss) {
        size_t t = 0;
        for (auto& x : ss) t += x.data.size() + STREAM_ENTRY;
        return t;
    };
    int best = ID_RAW;
    s.idbest.clear();
    v.resize(n);
    for (size_t j = 0; j < n; ++j) v[j] = s.rec[j].id;
    s.make_stream(s.idbest, K_ID_RAW, v);

    s.idtry.clear();
    {
        uint64_t prev_first = 0;
        for (size_t j = 0; j < n; ++j) {
            bool first = (j == 0 || s.rec[j].key != s.rec[j - 1].key);
            uint64_t ref = first ? prev_first : s.rec[j - 1].id;
            v[j] = zigzag(int64_t(s.rec[j].id - ref)); // wrapping difference
            if (first) prev_first = s.rec[j].id;
        }
        s.make_stream(s.idtry, K_ID_DELTA, v);
        if (total(s.idtry) < total(s.idbest)) {
            std::swap(s.idbest, s.idtry);
            best = ID_DELTA_IN_CELL;
        }
    }
    const uint64_t N = h.lattice_n;
    bool lattice = N > 0 && N < (1ull << 21);
    if (lattice) {
        const uint64_t N3 = N * N * N;
        for (size_t j = 0; j < n && lattice; ++j) lattice = s.rec[j].id < N3;
    }
    if (lattice) {
        s.idtry.clear();
        for (int d = 0; d < 3; ++d) {
            for (size_t j = 0; j < n; ++j) {
                uint64_t id = s.rec[j].id;
                int64_t q = int64_t(d == 0 ? id % N : d == 1 ? (id / N) % N : id / (N * N));
                int64_t l = int64_t((unsigned __int128)s.X[3 * s.rec[j].i + d] * N / uint64_t(Q));
                int64_t r = floor_mod(q - l + int64_t(N / 2), int64_t(N)) - int64_t(N / 2);
                v[j] = zigzag(r);
            }
            s.make_stream(s.idtry, K_ID_LAGR, v);
        }
        if (total(s.idtry) < total(s.idbest)) {
            std::swap(s.idbest, s.idtry);
            best = ID_LAGRANGIAN;
        }
    }
    for (auto& x : s.idbest) s.streams.push_back(std::move(x));

    // assemble
    size_t size = BLOCK_HEADER + STREAM_ENTRY * s.streams.size();
    for (auto& x : s.streams) size += x.data.size();
    s.out.assign(size, 0);
    uint8_t* b = s.out.data();
    std::memcpy(b, BLOCK_MAGIC, 4);
    put<uint32_t>(b, 4, uint32_t(n));
    put<uint32_t>(b, 8, n_cells);
    put<uint8_t>(b, 12, uint8_t(s.streams.size()));
    put<uint8_t>(b, 13, uint8_t(best));
    size_t off = BLOCK_HEADER, data = BLOCK_HEADER + STREAM_ENTRY * s.streams.size();
    stats = Stats{};
    stats.id_coder = best;
    for (auto& x : s.streams) {
        put<uint8_t>(b, off, x.width);
        put<uint8_t>(b, off + 1, x.kind);
        put<uint32_t>(b, off + 4, x.count);
        put<uint32_t>(b, off + 8, uint32_t(x.data.size()));
        off += STREAM_ENTRY;
        std::memcpy(b + data, x.data.data(), x.data.size());
        data += x.data.size();
        size_t t = x.data.size() + STREAM_ENTRY;
        (x.kind <= K_OFFSET ? stats.pos : x.kind == K_VEL ? stats.vel : stats.id) += t;
        stats.raw_bytes += size_t(x.count) * x.width;
    }
    stats.pos += BLOCK_HEADER;
    stats.total = size;
    return s.out;
}

// ---------------------------------------------------------------------------------------------
// Decoder
// ---------------------------------------------------------------------------------------------
namespace {
struct DCtx {
    ZSTD_DCtx* d = ZSTD_createDCtx();
    ~DCtx() { ZSTD_freeDCtx(d); }
};

void read_stream(const uint8_t* src, size_t csize, int w, uint32_t count, std::vector<uint64_t>& v, std::vector<uint8_t>& planes) {
    static thread_local DCtx dctx;
    check(w == 1 || w == 2 || w == 4 || w == 8, "bad stream width");
    size_t raw = size_t(count) * w;
    planes.resize(raw);
    size_t r = ZSTD_decompressDCtx(dctx.d, planes.data(), raw, src, csize);
    check(!ZSTD_isError(r) && r == raw, "stream decompression failed or size mismatch");
    v.assign(count, 0);
    for (int b = 0; b < w; ++b) {
        const uint8_t* p = planes.data() + size_t(b) * count;
        const int sh = 8 * b;
        for (uint32_t i = 0; i < count; ++i) v[i] |= uint64_t(p[i]) << sh;
    }
}
} // namespace

void decode_block(const Header& h, const uint8_t* b, size_t size, Block& out) {
    check(size >= BLOCK_HEADER && std::memcmp(b, BLOCK_MAGIC, 4) == 0, "bad block magic");
    const uint32_t n = get<uint32_t>(b, 4), n_cells = get<uint32_t>(b, 8);
    const int ns = get<uint8_t>(b, 12), coder = get<uint8_t>(b, 13);
    const int expect = 8 + (coder == ID_LAGRANGIAN ? 3 : 1);
    check(ns == expect && coder <= ID_LAGRANGIAN, "unexpected stream count or ID coder");
    check(size >= BLOCK_HEADER + STREAM_ENTRY * size_t(ns), "block truncated (stream table)");
    std::vector<std::vector<uint64_t>> st(ns);
    std::vector<uint8_t> planes;
    size_t off = BLOCK_HEADER, data = BLOCK_HEADER + STREAM_ENTRY * ns;
    for (int k = 0; k < ns; ++k) {
        int w = get<uint8_t>(b, off);
        uint32_t count = get<uint32_t>(b, off + 4), cs = get<uint32_t>(b, off + 8);
        check(count == (k < 2 ? n_cells : n), "stream count mismatch");
        check(data + cs <= size, "block truncated (stream data)");
        read_stream(b + data, cs, w, count, st[k], planes);
        off += STREAM_ENTRY;
        data += cs;
    }
    const int64_t Q = int64_t(h.q()), sub = h.sub, nc = h.ncell;
    out.id_coder = coder;
    out.X.resize(3 * size_t(n));
    out.V.resize(3 * size_t(n));
    out.id.resize(n);
    // cells -> positions
    uint64_t key = 0;
    size_t j = 0;
    for (uint32_t c = 0; c < n_cells; ++c) {
        key += st[0][c];
        int64_t cell[3] = {int64_t(key % nc), int64_t((key / nc) % nc), int64_t(key / (uint64_t(nc) * nc))};
        check(j + st[1][c] <= n, "cell counts exceed block size");
        for (uint64_t m = 0; m < st[1][c]; ++m, ++j)
            for (int d = 0; d < 3; ++d) {
                int64_t o = int64_t(st[2 + d][j]);
                check(o < sub, "offset out of range");
                out.X[3 * j + d] = cell[d] * sub + o;
            }
    }
    check(j == n, "cell counts do not add up to block size");
    for (size_t i = 0; i < n; ++i)
        for (int d = 0; d < 3; ++d) out.V[3 * i + d] = unzigzag(st[5 + d][i]);
    if (coder == ID_RAW) {
        for (size_t i = 0; i < n; ++i) out.id[i] = st[8][i];
    }
    else if (coder == ID_DELTA_IN_CELL) {
        uint64_t prev_first = 0;
        size_t i = 0;
        for (uint32_t c = 0; c < n_cells; ++c)
            for (uint64_t m = 0; m < st[1][c]; ++m, ++i) {
                uint64_t ref = (m == 0) ? prev_first : out.id[i - 1];
                out.id[i] = ref + uint64_t(unzigzag(st[8][i]));
                if (m == 0) prev_first = out.id[i];
            }
    }
    else {
        const uint64_t N = h.lattice_n;
        check(N > 0, "Lagrangian ID coder needs lattice_n");
        for (size_t i = 0; i < n; ++i) {
            int64_t q[3];
            for (int d = 0; d < 3; ++d) {
                int64_t l = int64_t((unsigned __int128)out.X[3 * i + d] * N / uint64_t(Q));
                q[d] = floor_mod(unzigzag(st[8 + d][i]) + l, int64_t(N));
            }
            out.id[i] = (uint64_t(q[2]) * N + uint64_t(q[1])) * N + uint64_t(q[0]);
        }
    }
}

void Block::positions(const Header& h, std::vector<double>& x) const {
    const double Q = double(h.q());
    x.resize(X.size());
    for (size_t i = 0; i < X.size(); ++i) x[i] = (double(X[i]) + 0.5) / Q - 0.5;
}

void Block::velocities(const Header& h, std::vector<double>& v) const {
    v.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) v[i] = double(V[i]) * h.vel_step;
}

// ---------------------------------------------------------------------------------------------
// Files
// ---------------------------------------------------------------------------------------------
FileWriter::FileWriter(const std::string& path, const Header& h) : path_(path), h_(h) {
    fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
    if (fd_ < 0) throw std::system_error(errno, std::generic_category(), "snapcomp: cannot create " + path);
    std::vector<uint8_t> b(HEADER_SIZE);
    h_.n_particles = h_.n_blocks = 0;
    h_.pack(b.data());
    pwrite_all(b.data(), b.size(), 0);
}

FileWriter::~FileWriter() {
    if (fd_ >= 0) {
        try {
            close();
        } catch (...) {
        }
    }
}

void FileWriter::pwrite_all(const void* data, size_t size, uint64_t offset) {
    const uint8_t* p = static_cast<const uint8_t*>(data);
    while (size) {
        ssize_t w = ::pwrite(fd_, p, size, off_t(offset));
        if (w < 0) {
            if (errno == EINTR) continue;
            throw std::system_error(errno, std::generic_category(), "snapcomp: write to " + path_);
        }
        p += w;
        size -= size_t(w);
        offset += uint64_t(w);
    }
}

void FileWriter::append(const void* data, size_t size, uint32_t n_particles, uint32_t thread) {
    check(fd_ >= 0, "append after close");
    uint64_t off = next_.fetch_add(size);
    pwrite_all(data, size, off);
    std::lock_guard<std::mutex> lock(mtx_);
    index_.push_back({off, size, n_particles, thread});
}

void FileWriter::close() {
    if (fd_ < 0) return;
    std::sort(index_.begin(), index_.end(), [](const Entry& a, const Entry& b) { return a.offset < b.offset; });
    uint64_t np = 0;
    std::vector<uint8_t> b(index_.size() * INDEX_ENTRY + TRAILER, 0);
    for (size_t i = 0; i < index_.size(); ++i) {
        put<uint64_t>(b.data(), i * INDEX_ENTRY, index_[i].offset);
        put<uint64_t>(b.data(), i * INDEX_ENTRY + 8, index_[i].size);
        put<uint32_t>(b.data(), i * INDEX_ENTRY + 16, index_[i].n);
        put<uint32_t>(b.data(), i * INDEX_ENTRY + 20, index_[i].thread);
        np += index_[i].n;
    }
    const uint64_t idx = next_.load();
    uint8_t* t = b.data() + index_.size() * INDEX_ENTRY;
    std::memcpy(t, TRAILER_MAGIC, 8);
    put<uint64_t>(t, 8, idx);
    put<uint64_t>(t, 16, index_.size());
    put<uint64_t>(t, 24, np);
    pwrite_all(b.data(), b.size(), idx);
    h_.n_particles = np;
    h_.n_blocks = index_.size();
    std::vector<uint8_t> hb(HEADER_SIZE);
    h_.pack(hb.data());
    pwrite_all(hb.data(), hb.size(), 0);
    if (::fsync(fd_) != 0) throw std::system_error(errno, std::generic_category(), "snapcomp: fsync " + path_);
    int fd = fd_;
    fd_ = -1;
    if (::close(fd) != 0) throw std::system_error(errno, std::generic_category(), "snapcomp: close " + path_);
}

namespace {
void pread_all(int fd, void* data, size_t size, uint64_t offset) {
    uint8_t* p = static_cast<uint8_t*>(data);
    while (size) {
        ssize_t r = ::pread(fd, p, size, off_t(offset));
        if (r < 0 && errno == EINTR) continue;
        check(r > 0, "short read");
        p += r;
        size -= size_t(r);
        offset += uint64_t(r);
    }
}
} // namespace

FileReader::FileReader(const std::string& path) {
    fd_ = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd_ < 0) throw std::system_error(errno, std::generic_category(), "snapcomp: cannot open " + path);
    try {
        struct stat st;
        check(::fstat(fd_, &st) == 0, "fstat failed");
        uint64_t size = uint64_t(st.st_size);
        check(size >= HEADER_SIZE + TRAILER, "file too small (incomplete?)");
        uint8_t t[TRAILER];
        pread_all(fd_, t, TRAILER, size - TRAILER);
        check(std::memcmp(t, TRAILER_MAGIC, 8) == 0, "no trailer: file incomplete (writer did not close)");
        uint64_t idx = get<uint64_t>(t, 8), nb = get<uint64_t>(t, 16);
        check(idx + nb * INDEX_ENTRY + TRAILER == size, "index does not match file size");
        std::vector<uint8_t> hb(HEADER_SIZE);
        pread_all(fd_, hb.data(), HEADER_SIZE, 0);
        h_ = Header::unpack(hb.data());
        check(h_.n_blocks == nb, "header and trailer disagree on block count");
        std::vector<uint8_t> b(nb * INDEX_ENTRY);
        if (nb) pread_all(fd_, b.data(), b.size(), idx);
        index_.resize(nb);
        for (uint64_t i = 0; i < nb; ++i) {
            index_[i] = {get<uint64_t>(b.data(), i * INDEX_ENTRY), get<uint64_t>(b.data(), i * INDEX_ENTRY + 8), get<uint32_t>(b.data(), i * INDEX_ENTRY + 16),
                         get<uint32_t>(b.data(), i * INDEX_ENTRY + 20)};
            check(index_[i].offset >= HEADER_SIZE && index_[i].offset + index_[i].size <= idx, "index entry out of range");
        }
    } catch (...) {
        ::close(fd_);
        throw;
    }
}

FileReader::~FileReader() {
    if (fd_ >= 0) ::close(fd_);
}

void FileReader::read_block(size_t i, Block& out) const {
    check(i < index_.size(), "block index out of range");
    std::vector<uint8_t> b(index_[i].size);
    pread_all(fd_, b.data(), b.size(), index_[i].offset);
    decode_block(h_, b.data(), b.size(), out);
    check(out.size() == index_[i].n, "block particle count differs from index");
}

} // namespace snapcomp
