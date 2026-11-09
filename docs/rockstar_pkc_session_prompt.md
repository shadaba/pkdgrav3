# Prompt for a new session: Rockstar reads compressed pkdgrav3 snapshots natively

Paste everything below the line into a new Claude Code session started in the sahyadri-sandbox clone
(https://github.com/a-paranjape/sahyadri-sandbox, branch master).

---

## Context: where pkdgrav3 output stands (2026-11-10)

We run Sahyadri-style dark-matter simulations with pkdgrav3 3.5.1 on pawna (4 nodes x 192 cores, 755 GB
each, slow NFS /storage at ~30 MB/s). Writing HDF5 snapshots took ~2.5 h per 2048^3 snapshot. We added a
compressed, parallel snapshot format to pkdgrav3:

- Code: https://github.com/shadaba/pkdgrav3, branch `feature/compressed-output`. Main docs there:
  `docs/compressed_io_plan.md` (plan and all measurements), `docs/snapcomp_format.md` (byte-level format
  spec, version 1), `docs/running_on_pawna.md` (setup, parameters, staging/drain),
  `docs/memory_persistent_ids.md`.
- Enabled with `bCompressedOutput = True` in the pkdgrav3 parameter file. Every thread compresses its own
  particles and writes in parallel. A 512^3 snapshot takes 1.4-1.9 s instead of 250-280 s for HDF5;
  2048^3 is expected at about 40-60 s.
- Size: about 5-8 bytes per particle with subsamples (HDF5 with float32 + IDs: 32 B; Sahyadri's own
  compressed format: ~9.9 B).
- Status: verified at 128^3 and 512^3 (every ID once, errors within half a quantization step, IDs
  persistent across snapshots). A 2048^3 test on one node is running. Python reader and a Gadget-2
  converter exist (below).

### The format in brief (read `docs/snapcomp_format.md` in the pkdgrav3 repo for the exact layout)

- One snapshot = several files: `<name>.<step>.s<pct>.<writer>.pkc`, e.g. `snap.00100.s01.0.pkc`. One file
  per subsample per writer (writer = MPI process; on pawna one per node). Plus a marker
  `snap.00100.done` when all files of a writer are complete. Every file stands alone.
- Subsamples (Sahyadri convention): 1, 3, 5, 9, 10, 12, 15, 20, 25 percent (sum 100). A particle's
  subsample is a hash of its ID (splitmix64, seed in the header), so it is the same particle set in every
  snapshot. Halo finding needs all subsamples together; analyses can read only the 1% file.
- File: 4096-byte header, then compressed blocks (<= 65536 particles each, written by many threads in
  any order), then a block index (32 B per block) and a 32-byte trailer. A file without a trailer is
  incomplete and must be rejected.
- Header holds everything a halo finder needs: n_total, box_size (Mpc/h), redshift, omega_m, omega_l, h,
  lattice_n (N, IDs are Lagrangian lattice indices `(iz*N+iy)*N+ix`), ncell, sub, vel_step, subsample
  index and percents, writer, block_size, step. Particle mass is not stored: m_p =
  2.77536627e11 * Omega_m * L^3 / n_total Msun/h.
- Block: particles sorted by (cell, ID). Streams (byte-plane shuffled, one zstd frame each, with
  checksums): occupied cell keys, particles per cell, position offset within cell (x, y, z), quantized
  velocities (zigzag), IDs (one of three coders, chosen per block).
- Decoded values: integer position X in [0, ncell*sub) per axis; position in Mpc/h
  = (X + 0.5) / (ncell*sub) * box_size, in [0, L). Velocity = V * vel_step, peculiar km/s (not Gadget's
  v/sqrt(a)). IDs uint64.
- Accuracy (Sahyadri default, set in the parameter file): position step 1.95 kpc/h, velocity step
  2.44 km/s (max error is half a step).
- IDs: persistent. pkdgrav3 runs in ordered mode (`bMemUnordered = False`), so the ID is iOrder =
  Lagrangian lattice index, at no memory cost. Needed for consistent-trees.

### Tools that already exist (pkdgrav3 repo)

- C++17 library `io/snapcomp/` (zstd + POSIX, no pkdgrav3 headers, standalone CMake build:
  `cmake -S io/snapcomp -B build`). `snapcomp::FileReader(path)`: `header()`, `n_blocks()`,
  `block_particles(i)`, `read_block(i, Block&)` (const and thread-safe: several threads can decode
  blocks of one file at once). `Block` has `X` (3n int64), `V` (3n int64), `id` (n uint64).
- `pkc_check`: verifies a snapshot without loading it (IDs once, subsamples, persistence with `--ref`).
- `pkc_to_gadget`: .pkc -> Gadget-2 files for Rockstar (`FILE_FORMAT = "GADGET2"`, `GADGET_ID_BYTES = 8`,
  v/sqrt(a)). Works, but needs 32 B/particle of extra disk (275 GB at 2048^3) and a conversion pass.
  This is the reference to compare the native reader against.
- Python reader `tools/pkcreader/pkcreader.py` (numpy + zstandard): `PkcSnapshot(prefix).read()` /
  `read_block("pos")` in SnapshotReader style.

### How Rockstar is run today

- sahyadri-sandbox `code/Rockstar/gfcstanford-rockstar-36ce9eea36ee/` holds only the files Aseem
  patched (see `scripts/rockstar/readme_aseem.txt`: AREPO/HDF5 reader for Gadget-4, TIRPC guards, curvature
  term); the rest is upstream Rockstar at bitbucket.org/gfcstanford/rockstar commit 36ce9eea36ee.
- On pawna, pkdgrav3 runs are processed from TIPSY output (`FILE_FORMAT = "TIPSY"`, IDs from file order),
  one Rockstar call per snapshot (`scripts/run_rockstar_batch.sh` in a colleague's clone; config like
  `config/halos/rockstar_n1024_tipsy.cfg`). Production directories of other users are read-only.

## Goal of this session

Plan (do not implement yet) a native reader so Rockstar takes `FILE_FORMAT = "PKC"` and reads the
compressed snapshot directly, with no conversion step and no extra disk. Write the plan as a doc in this
repo and stop for review. The plan should cover:

1. Where the code goes: a new `io/io_pkc.c` (+ header) in our Rockstar tree, hooked into
   `io/meta_io.c` and `config.template.h` with the fewest edits; Makefile target. Choose between
   (a) a small C99 decoder written from the format spec with the zstd C API, and (b) an `extern "C"`
   wrapper around the C++ snapcomp library linked with g++. Say which and why (Rockstar is C99; keeping
   the spec as the single source of truth argues for (a); test against (b) / `pkc_to_gadget`).
2. How Rockstar's blocks map onto .pkc files: `NUM_BLOCKS` files per snapshot vs. our files per
   subsample and writer of very different sizes (1% to 25%). Probably: Rockstar block k reads every
   (file, pkc-block) pair with index = k mod NUM_BLOCKS, so all readers get similar loads and any
   NUM_BLOCKS works. Check this against how `PARALLEL_IO` readers and `load_particles()` are used.
3. Header -> Rockstar globals: SCALE_NOW = 1/(1+z), BOX_SIZE, Om, Ol, h0, TOTAL_PARTICLES,
   PARTICLE_MASS, AVG_PARTICLE_SPACING; units (Mpc/h and peculiar km/s, i.e. no sqrt(a)); IDs 64-bit.
   Reject incomplete files (no trailer) and snapshots whose files do not add up to n_total.
4. Options: a subsample list (e.g. `PKC_SUBSAMPLES`, default all) for quick tests on 1-10%;
   filename pattern (`FILENAME = "snap.<snap>"`, files found by `<name>.s*.pkc`).
5. Memory and speed: decode one pkc block at a time into Rockstar's particle array; threads if useful.
   Estimate at 2048^3 on one pawna node.
6. Validation: on a 512^3 run with both outputs (bCompAlsoStandard = True), compare Rockstar catalogs
   from (i) the exact snapshot (TIPSY or HDF5), (ii) pkc_to_gadget, (iii) native PKC. (ii) and (iii)
   must agree exactly; (i) vs (iii) gives the effect of quantization (halo counts, mass function,
   matched halo masses, positions, velocities, Vmax; consistent-trees on a few snapshots).
7. Risks and open questions.

Rules: work on a branch of this repo (or a fork), never push to master without asking. Do not run
simulations larger than 512^3 without approval. Other users' files on pawna are read-only; jobs go
through PBS. Building third-party source downloaded from the internet may need explicit approval in
Claude Code; ask before doing it. State how every measurement was made.
