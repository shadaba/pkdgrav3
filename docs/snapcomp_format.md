# snapcomp file format, versions 1 to 3

| version | adds | written when |
|---|---|---|
| 1 | - | the potential is not stored and the writer gives no binary cosmology (standalone tools such as `pkc_convert`) |
| 2 | the gravitational potential: `pot_step`, `pot_to_kms2` in the header, one stream per block | potential stored, no binary cosmology |
| 3 | a binary cosmology block (10 doubles at 4000); the potential is optional (present if `pot_step > 0`) | **pkdgrav3 from 2026-11-11 on** (every snapshot) |

Readers accept all three. The info text ends before 4000 (version 3), 4080 (version 2) or 4096 (version 1).

Compressed dark-matter snapshot written in parallel by pkdgrav3 threads (Option A of
[compressed_io_plan.md](compressed_io_plan.md)). Reference implementation: `io/snapcomp/`
(C++, no pkdgrav3 headers); independent Python reader: `tools/pkcreader/pkcreader.py`.
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

| offset | type | name | unit | meaning |
|---|---|---|---|---|
| 0 | char[8] | magic | - | `PKCOMP\0\1` |
| 8 | u32 | version | - | 1, or 2 if the potential is stored |
| 12 | u32 | header_size | bytes | 4096 |
| 16 | u64 | n_particles | - | particles in this file |
| 24 | u64 | n_blocks | - | blocks in this file |
| 32 | u64 | lattice_n | - | N if IDs are Lagrangian indices `(iz*N+iy)*N+ix`, else 0 |
| 40 | u64 | n_total | - | particles in the whole snapshot (all subsamples, all writers) |
| 48 | u32 | ncell | - | cell grid per dimension |
| 52 | u32 | sub | - | position steps per cell |
| 56 | f64 | vel_step | km/s | velocity quantum (max error vel_step/2) |
| 64 | f64 | vel_to_kms | km/s per writer unit | factor from the writer's velocity unit to km/s (informational) |
| 72 | f64 | box_size | Mpc/h, comoving | box side L |
| 80 | f64 | redshift | - | z of the snapshot |
| 88 | f64 | time | pkdgrav3 code time | for physical time use `cosmology.age_gyr` in the info JSON |
| 96 | f64 | omega_m | - | Omega_m today |
| 104 | f64 | omega_l | - | Omega_Lambda today (0 for w0-wa runs, see info JSON `Omega_DE`) |
| 112 | f64 | h | - | H0 / (100 km/s/Mpc) |
| 120 | u32 | subsample | - | index of this subsample (0-based) |
| 124 | u32 | n_subsamples | - | |
| 128 | f64[16] | percents | percent | share of particles per subsample (first n_subsamples used) |
| 256 | u64 | subsample_seed | - | seed of the ID hash |
| 264 | u32 | writer | - | writer index (MPI process) |
| 268 | u32 | n_writers | - | |
| 272 | u32 | block_size | particles | maximum particles per block |
| 276 | i32 | zstd_level | - | |
| 280 | u64 | step | - | simulation step (from the snapshot name `<name>.<step>`) |
| 288 | char[...] | info | - | NUL-terminated text, ends before 4000 (version 3, at most 3711 characters), 4080 (version 2) or 4096 (version 1); pkdgrav3 writes JSON, see "Info JSON" below |
| 4000 | f64 | omega_b | - | version 3: Omega_b today |
| 4008 | f64 | omega_rad | - | version 3: Omega_rad today (photons + massless neutrinos as pkdgrav3 uses them; 0 if not used) |
| 4016 | f64 | omega_de | - | version 3: Omega of the dynamical dark energy component (0 if only a cosmological constant, see omega_l) |
| 4024 | f64 | w0 | - | version 3: w(a) = w0 + wa (1 - a) of omega_de |
| 4032 | f64 | wa | - | version 3 |
| 4040 | f64 | sigma8 | - | version 3: power spectrum normalisation (0 if `normalization` is used) |
| 4048 | f64 | n_s | - | version 3: primordial spectral index |
| 4056 | f64 | running | - | version 3: running of n_s |
| 4064 | f64 | normalization | - | version 3: pkdgrav3 dNormalization (alternative to sigma8; 0 if sigma8 is used) |
| 4072 | f64 | pivot | 1/Mpc | version 3: pivot scale of n_s and running |
| 4080 | f64 | pot_step | (km/s)^2 | versions 2, 3: potential quantum (max error pot_step/2); > 0 in version 2; 0 in version 3 = no potential |
| 4088 | f64 | pot_to_kms2 | (km/s)^2 per writer unit | versions 2, 3 (if pot_step > 0): factor from the writer's potential unit to (km/s)^2 (informational) |

The background expansion of any run follows from the binary header alone (version 3):
E(a)^2 = omega_m a^-3 + omega_rad a^-4 + omega_l + omega_de a^(-3(1+w0+wa)) exp(-3 wa (1-a)) + omega_k a^-2, with
omega_k = 1 - omega_m - omega_l - omega_de - omega_rad.

Position step = box_size / (ncell * sub) Mpc/h; velocity step = vel_step km/s; potential step = pot_step (km/s)^2.

### Particle data: quantities and units

| quantity | stored as | unit after decoding | accuracy |
|---|---|---|---|
| position | cell key + offset in cell (integers) | comoving Mpc/h, in [0, L) (box units [-0.5, 0.5) before multiplying by box_size and shifting) | max error box_size/(ncell*sub)/2 per coordinate (0.98 kpc/h for the default 1.953 kpc/h step) |
| velocity | integer multiple of vel_step | km/s, peculiar velocity v = a dx/dt (not Gadget's v/sqrt(a)) | max error vel_step/2 per component |
| ID | u64 (lossless) | - | exact; for pkdgrav3 ICs the Lagrangian index `(iz*N+iy)*N+ix` |
| potential (pot_step > 0) | integer multiple of pot_step | (km/s)^2, peculiar gravitational potential in physical units (pkdgrav3 code potential * U^2 / a, U = code velocity unit in km/s); zero point: periodic, mean density removed; includes pkdgrav3's constant softened self term -(105/64) m/eps | max error pot_step/2 |
| mass | not stored | Msun/h: all particles have `box.particle_mass_msunh` = rho_crit Omega_m L^3 / n_total (rho_crit = 2.77536627e11 h^2 Msun/Mpc^3) | exact |

### Info JSON (pkdgrav3, snapshots written from 2026-11-10 on)

Top level:

| key | unit | meaning |
|---|---|---|
| `code` | - | "pkdgrav3 <version>" |
| `pos_step_kpch_requested` | kpc/h, comoving | dCompPosStep as requested (the actual step is box_size/(ncell*sub)) |
| `vel_step_kms` | km/s | dCompVelStep |
| `pot_step_kms2` | (km/s)^2 | potential step of this snapshot, 0 = not stored |
| `name` | - | snapshot path as written |

`sim`: provenance

| key | unit | meaning |
|---|---|---|
| `name`, `set`, `comment` | - | parameters achSimName (default: run directory name), achSimSet, achSimComment |
| `code`, `code_git` | - | pkdgrav3 version; `git describe` of the source at configure time (`-dirty` = local changes) |
| `written_utc` | ISO 8601, UTC | when the snapshot was written |
| `host` | - | host of the master process |
| `mpi_processes`, `threads` | - | MPI processes; worker threads in total (pkdgrav3 Threads()) |
| `gpu` | - | true if built with CUDA and bGPU is on |

`box`

| key | unit | meaning |
|---|---|---|
| `L_mpch`, `L_mpc` | Mpc/h, Mpc (comoving) | box side |
| `n_particles` | - | particles in the snapshot |
| `n_grid` | - | IC grid per dimension (0 if ICs were read from a file) |
| `particle_mass_msunh`, `particle_mass_msun` | Msun/h, Msun | particle mass |
| `softening_kpch` | kpc/h, comoving | gravitational softening length (pkdgrav3 kernel parameter eps) |
| `softening_source` | - | "dSoft", "default 1/(50 nGrid)" or "per particle (input file)" (then softening_kpch is not meaningful) |
| `softening_max` | - | dSoftMax (code units; a multiplier if `softening_max_is_multiplier`) |
| `max_physical_softening` | code units | dMaxPhysicalSoft (0 = none) |

`ic`

| key | unit | meaning |
|---|---|---|
| `z_initial` | - | dRedFrom |
| `lpt_order` | - | iLPT |
| `seed` | - | iSeed |
| `fixed_amplitude`, `phase_pi` | -, units of pi | bFixedAmpIC, dFixedAmpPhasePI |
| `input_file`, `tf_file`, `class`, `class_file` | - | achInFile (ICs read from a file), achTfFile, bClass, achClassFilename |

`cosmology`: parameters (today)

| key | unit | meaning |
|---|---|---|
| `h`, `H0` | -, km/s/Mpc | Hubble parameter |
| `Omega_m`, `Omega_lambda`, `Omega_b`, `Omega_rad`, `Omega_DE`, `Omega_K` | - | density parameters; Omega_K = 1 - Omega_m - Omega_lambda - Omega_DE - Omega_rad |
| `w0`, `wa` | - | dark energy w(a) = w0 + wa (1 - a), for the Omega_DE component |
| `omega_b`, `omega_cdm` | - | Omega_b h^2, (Omega_m - Omega_b) h^2 |
| `sigma8`, `normalization` | -, - | power spectrum normalisation (one of them is used) |
| `n_s`, `running`, `pivot` | -, -, 1/Mpc | primordial spectrum |
| `class`, `class_file`, `tf_file` | - | background/transfer source |

`cosmology`: values at this snapshot (from pkdgrav3's own background, so w0-wa, radiation and CLASS runs are covered)

| key | unit | meaning |
|---|---|---|
| `redshift`, `scale_factor` | - | z, a |
| `E_z`, `H_z_kms_mpc` | -, km/s/Mpc | H(z)/H0, H(z) |
| `Omega_m_z`, `Omega_DE_z` | - | Omega_m a^-3/E^2; (Omega_lambda + DE density)/E^2 |
| `growth_D` | - | linear growth factor D(z)/D(z = 0) |
| `growth_f` | - | growth rate dlnD/dlna |
| `age_gyr`, `hubble_time_gyr` | Gyr | cosmic age at a; 1/H0 |
| `comoving_distance_mpch` | Mpc/h | line-of-sight comoving distance from z = 0 |

`numerics`

| key | unit | meaning |
|---|---|---|
| `theta` | - | tree opening angle at this time (dTheta, may switch with redshift) |
| `replicas` | - | periodic replicas in the tree walk (nReplicas at this time) |
| `ewald`, `ewald_cut` | -, box lengths | Ewald correction on/off and dEwCut |
| `eta`, `eps_acc_step`, `new_kdk` | - | time-step criterion dEta, bEpsAccStep, bNewKDK |

`steps`

| key | unit | meaning |
|---|---|---|
| `step`, `total_steps` | - | step of this snapshot; sum of nSteps |
| `steps_per_segment` | - | nSteps (one number if all segments are equal, else a list) |
| `segment_end_redshifts` | - | dRedTo list (6 significant digits) |
| `output_interval` | steps | iOutInterval |

Files written before these keys existed have only `code`, the steps and `name`; readers should treat
every key as optional (`tools/pkcreader` returns empty dicts).

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
| last | 8 | n | only if pot_step > 0 (versions 2, 3): potential, zigzag of the difference to the previous particle in block order (0 before the first) |

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
