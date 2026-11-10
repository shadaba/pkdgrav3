# snapcomp file format, versions 1 and 2

Version 2 = version 1 plus the gravitational potential (two header fields, one stream per block). A
file without the potential is written as version 1, byte for byte as before. Readers accept both.

Compressed dark-matter snapshot written in parallel by pkdgrav3 threads (Option A of
[compressed_io_plan.md](compressed_io_plan.md)). Reference implementation: `io/snapcomp/`
(C++, no pkdgrav3 headers); independent Python reader: `tools/snapcomp_proto/pkc_reader.py`.
All integers are little-endian. Sizes are in bytes.

## Snapshot = set of files

One snapshot is written as one file per (subsample, writer). A writer is a group of threads that
share a file (an MPI rank on pawna: one writer for the whole node). The pkdgrav3 integration names them
`<name>.<step>.s<pct>.<writer>.pkc`, e.g. `snap.00100.s01.0.pkc`. Each file is self-contained: it can be
read, copied or deleted on its own. Subsamples are disjoint; together they hold every particle once.

## File layout

```
offset 0      header      4096 bytes (fixed)
offset 4096   blocks      any number, in any order, written concurrently by many threads
              index       nblocks x 32 bytes
end - 32      trailer     32 bytes
```

### Header (4096 bytes)

| offset | type | name | meaning |
|---|---|---|---|
| 0 | char[8] | magic | `PKCOMP\0\1` |
| 8 | u32 | version | 1, or 2 if the potential is stored |
| 12 | u32 | header_size | 4096 |
| 16 | u64 | n_particles | particles in this file |
| 24 | u64 | n_blocks | blocks in this file |
| 32 | u64 | lattice_n | N if IDs are Lagrangian indices `(iz*N+iy)*N+ix`, else 0 |
| 40 | u64 | n_total | particles in the whole snapshot (all subsamples, all writers) |
| 48 | u32 | ncell | cell grid per dimension |
| 52 | u32 | sub | position steps per cell |
| 56 | f64 | vel_step | velocity quantum (km/s) |
| 64 | f64 | vel_to_kms | factor from the writer's velocity unit to km/s |
| 72 | f64 | box_size | comoving box (Mpc/h) |
| 80 | f64 | redshift | |
| 88 | f64 | time | pkdgrav3 code time |
| 96 | f64 | omega_m | |
| 104 | f64 | omega_l | |
| 112 | f64 | h | |
| 120 | u32 | subsample | index of this subsample (0-based) |
| 124 | u32 | n_subsamples | |
| 128 | f64[16] | percents | percent of particles per subsample (first n_subsamples used) |
| 256 | u64 | subsample_seed | seed of the ID hash |
| 264 | u32 | writer | writer index |
| 268 | u32 | n_writers | |
| 272 | u32 | block_size | maximum particles per block |
| 276 | i32 | zstd_level | |
| 280 | u64 | step | simulation step |
| 288 | char[...] | info | free text (NUL-terminated), ends before 4080 in version 2 (at most 3791 bytes). pkdgrav3 writes JSON (about 1.7-2.6 KB) with `code`, the steps, `name`, and the groups `sim` (name, set, comment from achSimName/achSimSet/achSimComment, code version and git describe, UTC write time, host, MPI processes, threads, gpu), `box` (L in Mpc/h and Mpc, particle count, nGrid, particle mass in Msun/h and Msun, softening in comoving kpc/h and its source, softening limits), `ic` (initial redshift, LPT order, seed, fixed amplitude/phase, input/transfer/CLASS files), `cosmology` (h, H0, Omega_m, Omega_lambda, Omega_b, Omega_rad, Omega_DE, Omega_K, w0, wa, omega_b, omega_cdm, sigma8/normalization, n_s, running, pivot, CLASS flag and files; and at this snapshot: redshift, scale factor, E(z), H(z), Omega_m(z), Omega_DE(z), growth D(z)/D(0), growth rate f, age in Gyr, Hubble time, comoving distance in Mpc/h), `numerics` (current theta, replicas, Ewald and its cut, eta, eps_acc_step, new_kdk) and `steps` (this step, total steps, steps per segment, segment end redshifts, output interval). The binary omega_m/omega_l/h stay for simple readers; for w0-wa dark energy omega_l is 0. |
| 4080 | f64 | pot_step | version 2 only: potential quantum, (km/s)^2 (> 0) |
| 4088 | f64 | pot_to_kms2 | version 2 only: factor from the writer's potential unit to (km/s)^2 |

Position step = box / (ncell * sub); velocity step = vel_step km/s; potential step = pot_step (km/s)^2.

### Block

```
block header   16 bytes
stream table   n_streams x 12 bytes
stream data    n_streams zstd frames, back to back
```

Block header:

| offset | type | name |
|---|---|---|
| 0 | char[4] | `PKB1` |
| 4 | u32 | n (particles in the block, >= 1) |
| 8 | u32 | n_cells (occupied cells) |
| 12 | u8 | n_streams |
| 13 | u8 | id_coder (0 raw, 1 delta_in_cell, 2 lagrangian residual) |
| 14 | u16 | reserved (0) |

Stream table entry: `u8 width` (1, 2, 4 or 8), `u8 kind`, `u16 reserved`, `u32 count`,
`u32 compressed_size`. A stream holds `count` unsigned integers of `width` bytes, stored as byte planes
(all lowest bytes, then all second bytes, ...), compressed as one zstd frame (with content checksum).

Signed values are zigzag-coded: `u = (s << 1) ^ (s >> 63)`.

Streams in order (`kind` value in brackets):

| # | kind | count | content |
|---|---|---|---|
| 0 | 1 | n_cells | occupied cell keys, first value then differences (keys ascending) |
| 1 | 2 | n_cells | particles per occupied cell |
| 2-4 | 3 | n | position offset in its cell, x, y, z (0 .. sub-1) |
| 5-7 | 4 | n | velocity x, y, z, zigzag of the quantized value |
| 8 | 5 | n | id_coder 0: IDs as they are |
| 8 | 6 | n | id_coder 1: zigzag of ID differences (see below) |
| 8-10 | 7 | n | id_coder 2: zigzag of Lagrangian residual x, y, z |
| last | 8 | n | version 2 only: potential, zigzag of the difference to the previous particle in block order (0 before the first) |

Particles in a block are sorted by cell key, then by ID. Decoding, for particle j in cell c:

- cell key `k = (cz*ncell + cy)*ncell + cx`; integer position `X = c*sub + offset` (per dimension, in
  `[0, ncell*sub)`); position in box units `x = (X + 0.5)/(ncell*sub) - 0.5`; in Mpc/h multiply by box_size.
- velocity `v = V * vel_step` km/s (peculiar).
- potential (version 2): `P_j = P_{j-1} + unzigzag(u_j)` with `P_{-1} = 0`; `phi = P * pot_step` in (km/s)^2.
  This is the peculiar gravitational potential in physical units (pkdgrav3: code potential * U^2 / a, with
  U the code velocity unit in km/s); its zero point is pkdgrav3's (periodic, mean density removed).
- id_coder 1: the stream holds `d_j = id_j - id_{j-1}` inside a cell and, for the first particle of a
  cell, `id_first(c) - id_first(previous cell)` (previous first ID = 0 for the first cell).
- id_coder 2 (needs lattice_n = N > 0): `l = floor(X*N/(ncell*sub))` per dimension (Eulerian position in
  lattice units), lattice coordinate `q = (r + l) mod N`, `id = (q_z*N + q_y)*N + q_x`.

The writer picks the smallest ID coding per block.

### Index entry (32 bytes)

`u64 offset`, `u64 size` (whole block), `u32 n`, `u32 writer_thread`, `u64 reserved`.

### Trailer (32 bytes)

`char[8] "PKCIDX\0\1"`, `u64 index_offset`, `u64 n_blocks`, `u64 n_particles`.

A file without a valid trailer is incomplete (writer crashed); its blocks can still be recovered by
scanning for `PKB1` from offset 4096.

## Subsample membership

`subsample(id) = number of cumulative edges <= u`, where `u = (splitmix64(id, seed) >> 11) / 2^53`,
edges = cumulative sum of percents / 100 (as doubles, in order), and

```
splitmix64(x, seed): z = x + 0x9E3779B97F4A7C15 * (seed + 1)  (mod 2^64)
                     z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
                     z = (z ^ (z >> 27)) * 0x94D049BB133111EB
                     return z ^ (z >> 31)
```

Depends only on the ID, so a particle is in the same subsample in every snapshot.
