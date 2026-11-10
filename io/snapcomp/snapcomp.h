/*  snapcomp: compressed, parallel snapshot output for N-body particles (ID, position, velocity and,
 *  optionally, gravitational potential).
 *
 *  Self-contained: depends only on the C++17 standard library, POSIX and zstd; no pkdgrav3 headers.
 *  File format: docs/snapcomp_format.md.
 *
 *  Typical use, one writer (file) per subsample shared by all threads of a process:
 *
 *      snapcomp::Header h = ...;                    // grid, steps, cosmology, subsample, ...
 *      snapcomp::FileWriter w(path, h);             // creates the file
 *      // in each thread:
 *      snapcomp::BlockEncoder enc(h);               // per-thread scratch, about 10 MB at 65536 particles
 *      auto& bytes = enc.encode(n, pos, vel, ids);  // one block, pos in box units [-0.5, 0.5)
 *      w.append(bytes.data(), bytes.size(), n);     // thread-safe, positional write
 *      // once all threads are done:
 *      w.close();                                   // index, trailer, header, fsync
 *
 *  This file is part of pkdgrav3 (GPLv3).
 */
#ifndef SNAPCOMP_H
#define SNAPCOMP_H

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

namespace snapcomp {

constexpr uint32_t FORMAT_VERSION = 2; // written only when the potential is stored; otherwise 1 (same bytes as before)
constexpr uint32_t HEADER_SIZE = 4096;
constexpr int MAX_SUBSAMPLES = 16;

enum IdCoder : uint8_t { ID_RAW = 0, ID_DELTA_IN_CELL = 1, ID_LAGRANGIAN = 2 };

/// File header (see docs/snapcomp_format.md for the byte layout).
struct Header {
    uint64_t n_particles = 0; // set by FileWriter::close
    uint64_t n_blocks = 0;    // set by FileWriter::close
    uint64_t lattice_n = 0;   // N if IDs are (iz*N+iy)*N+ix, else 0
    uint64_t n_total = 0;     // particles in the whole snapshot
    uint32_t ncell = 0;       // cell grid per dimension
    uint32_t sub = 0;         // position steps per cell
    double vel_step = 0;      // km/s
    double vel_to_kms = 1;    // writer velocity unit -> km/s
    double box_size = 0;      // Mpc/h
    double redshift = 0, time = 0, omega_m = 0, omega_l = 0, h = 0;
    uint32_t subsample = 0, n_subsamples = 1;
    double percents[MAX_SUBSAMPLES] = {100};
    uint64_t subsample_seed = 0;
    uint32_t writer = 0, n_writers = 1;
    uint32_t block_size = 65536;
    int32_t zstd_level = 3;
    uint64_t step = 0;
    double pot_step = 0;    // (km/s)^2; 0 = no potential stored (version 1 file)
    double pot_to_kms2 = 1; // writer potential unit -> peculiar potential in (km/s)^2
    std::string info;       // free text, at most 4080-289 bytes

    bool has_potential() const { return pot_step > 0; }

    uint64_t q() const { return uint64_t(ncell) * sub; } // position steps per box length

    void pack(uint8_t* buf) const;            // HEADER_SIZE bytes
    static Header unpack(const uint8_t* buf); // throws std::runtime_error on bad magic/version
};

/// Choose ncell and sub for a requested position step (same length unit as box_size).
/// Cells of about the mean particle spacing (ncell ~ cbrt(n_total)), as tested in Phase 1.
void choose_grid(double box_size, double pos_step, uint64_t n_total, uint32_t& ncell, uint32_t& sub);

/// Subsample index of a particle (0 .. n-1) from its ID; identical to the Python reference.
uint64_t splitmix64(uint64_t x, uint64_t seed);
struct SubsampleTable {
    double edges[MAX_SUBSAMPLES] = {};
    int n = 1;
    uint64_t seed = 0;
    SubsampleTable() = default;
    SubsampleTable(const double* percents, int n, uint64_t seed); // percents must sum to 100
    int operator()(uint64_t id) const;
};

/// Per-thread block encoder. Not thread-safe; one per thread. Reuses its buffers.
class BlockEncoder {
  public:
    explicit BlockEncoder(const Header& h);
    ~BlockEncoder();
    BlockEncoder(const BlockEncoder&) = delete;
    BlockEncoder& operator=(const BlockEncoder&) = delete;

    /// Encode n >= 1 particles. pos: 3n doubles in box units (periodic, wrapped into [-0.5, 0.5));
    /// vel: 3n floats in the writer's unit (converted with vel_to_kms); ids: n values;
    /// pot: n floats in the writer's unit (converted with pot_to_kms2), required iff h.has_potential().
    /// Returns the block bytes, valid until the next call.
    const std::vector<uint8_t>& encode(size_t n, const double* pos, const float* vel, const uint64_t* ids, const float* pot = nullptr);

    /// Heap bytes held by this encoder (capacity of its buffers), for memory accounting.
    size_t scratch_bytes() const;

    /// Bytes per particle of the last block, per stream group (pos, vel, id) and the ID coder used.
    struct Stats {
        size_t pos = 0, vel = 0, id = 0, pot = 0, total = 0;
        int id_coder = 0;
        size_t raw_bytes = 0; // uncompressed stream bytes given to zstd
    } stats;

  private:
    struct Impl;
    Impl* p;
};

/// Decoded block: particles in block order (by cell, then ID).
struct Block {
    std::vector<int64_t> X; // 3n quantized positions in [0, ncell*sub)
    std::vector<int64_t> V; // 3n quantized velocities (units of vel_step)
    std::vector<uint64_t> id;
    std::vector<int64_t> P; // n quantized potentials (units of pot_step); empty if not stored
    int id_coder = 0;
    size_t size() const { return id.size(); }
    void positions(const Header& h, std::vector<double>& x) const;  // box units [-0.5, 0.5)
    void velocities(const Header& h, std::vector<double>& v) const; // km/s
    void potentials(const Header& h, std::vector<double>& p) const; // peculiar potential, (km/s)^2
};
void decode_block(const Header& h, const uint8_t* data, size_t size, Block& out); // throws on corruption

/// One output file shared by many threads. append() is thread-safe and lock-free for the data;
/// blocks land at increasing offsets reserved with an atomic counter (pwrite, no seeking).
class FileWriter {
  public:
    FileWriter(const std::string& path, const Header& h); // creates/truncates; throws on error
    ~FileWriter();                                        // closes if still open (no exception)
    void append(const void* data, size_t size, uint32_t n_particles, uint32_t thread = 0);
    void close(); // index + trailer + final header + fsync
    uint64_t bytes_written() const { return next_.load(); }
    /// Write-through (Linux): each append waits until its block is on disk and drops it from the page
    /// cache, so a snapshot needs no free RAM. Use when the node's memory is full (else it only adds
    /// the disk time to the stall). Off by default.
    void set_write_through(bool on) { write_through_ = on; }

  private:
    struct Entry {
        uint64_t offset, size;
        uint32_t n, thread;
    };
    std::string path_;
    Header h_;
    int fd_ = -1;
    std::atomic<uint64_t> next_{HEADER_SIZE};
    std::mutex mtx_;
    std::vector<Entry> index_;
    bool write_through_ = false;
    void pwrite_all(const void* data, size_t size, uint64_t offset);
};

class FileReader {
  public:
    explicit FileReader(const std::string& path); // throws if the file is incomplete or invalid
    ~FileReader();
    const Header& header() const { return h_; }
    size_t n_blocks() const { return index_.size(); }
    uint64_t block_particles(size_t i) const { return index_[i].n; }
    void read_block(size_t i, Block& out) const;

  private:
    struct Entry {
        uint64_t offset, size;
        uint32_t n, thread;
    };
    int fd_ = -1;
    Header h_;
    std::vector<Entry> index_;
};

} // namespace snapcomp
#endif
