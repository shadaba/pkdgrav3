# pkdgrav3: compressed, parallel snapshot output (plan and Claude Code brief)

Status: draft v1, 2026-11-06. Open items at the end.

## 1. Goal
Write 100 snapshots of 2048^3 particles (ID, position, velocity) from pkdgrav3 in a compressed format, in parallel, with negligible extra RAM. Target: cut the simulation stall per snapshot by at least 20x versus today (about 30 MB/s, roughly 2.5 h per 275 GB snapshot) and cut bytes on disk by at least 3x versus raw float32 + uint64 ID.

## 2. Environment and constraints
- pawna: single node, `-np 1`, 192 pthreads, NFS output. Multi-rank MPI fails there (SIGABRT in MPI cache layer).
- pegasus: 4 nodes x 32 threads (MPI).
- Current output: HDF5, gzip level 1, `bParaWrite` does not work, potential written only when requested.
- IDs are required in all snapshots (Rockstar + consistent-trees).
- The simulation already uses most of RAM: extra memory from the writer must stay under about 1 percent of the run's RSS.
- Must work for both `-np 1` pthreads and multi-node MPI. Each MDL thread owns its own PKD, so per-thread compress-and-write works for both.
- pkdgrav3 is GPLv3; anything linked into it is GPLv3.

## 3. Findings from the code (wullm/pkdgrav3, commit 7a5041e)
- `MSR::Write` (master.cxx) calls `Reorder()` (global sort by ID) before every write, then sets `iLastRungRT = iLastRungDD = -1`. The next step therefore redoes a full domain decomposition. A cell-sorted format does not need either.
- `pstWrite` (pst.c): with `bParaWrite=0` one thread loops over all other threads (`PST_SENDPARTICLES`) and writes every particle itself through `writeParticle` -> `fioWriteDark` (pkd.c, io/fio.c). One thread does all conversion, checksumming and I/O.
- HDF5 writer (io/fio.c): per-particle double-to-float conversion, 65536-row chunks, Fletcher32, extendable datasets, no MPI-IO. `AllNodeWrite` sets `FIO_FLAG_POTENTIAL | FIO_FLAG_DENSITY` unconditionally, so zero-filled potential/density columns look like they are written.
- The tree I read has NO gzip/deflate filter on the HDF5 datasets. The user's build has gzip level 1, so it is a patched or different tree. If gzip -1 runs in the single writer thread it alone may cap throughput near 30-60 MB/s.
- Only the HDF5 read-open is mutex-protected (`fio.c` around line 3567). Concurrent HDF5 writes from pthreads are unsafe unless HDF5 is built thread-safe (`H5_HAVE_THREADSAFE`). This is a likely reason `bParaWrite` fails.
- Existing Gadget compression (sahyadri-sandbox `scripts/post-process/Nbody_Compression.py`, `compression12bit.py`): positions sorted by global cell (Ngrid), 12-bit cell counts with overflow list, uint8 offset within cell (3 B); velocities 12-bit (4.5 B); IDs 32/64-bit, uncompressed; gzip applied by h5py. With 8 B IDs this is about 16 B per particle (2x vs 32 B).

## 4. Design
Default (Option A):
1. Own block format, one file per writer (Gadget-style multi-file), plain `pwrite` of large buffers, small header + block-index sidecar. No HDF5 or MPI-IO in the hot path.
2. Each MDL thread compresses and writes its own particles. No gathering to one writer, no `Reorder()`.
3. Blocks of about 2^18 particles in local tree order (spatially coherent). Within a block: sort by cell, then encode as in the Gadget scheme (cell counts, uint8 offsets, 12-bit velocities), so the existing decoder logic carries over. Per-thread scratch about 10 MB; a global permutation array would cost about 34 GB, so do not use one.
4. Entropy stage per block: zstd or blosc2 on the quantized streams (replaces gzip).
5. IDs (largest remaining cost; required): per block, if IDs are Lagrangian grid indices split into (qx,qy,qz) and store per-block minimum plus bit-packed offsets. Expected 2-4 B per particle, to be measured. Fallback: 33-bit packing.
6. File layout for NFS: about 16-64 writer files, writes of 4-16 MB, one fsync per file at the end, no small writes.
7. Optional staging: write to node-local scratch (if available), then copy to NFS in the background so the simulation does not wait on NFS.

Option B to test offline (decision experiment, not a commitment): Lagrangian (ID) order with implicit IDs and lattice-predictor coding of positions/velocities (Lorenzo-style residuals + zstd). Compare bytes per particle against Option A at several redshifts.

## 5. Phases and exit criteria
0. Measure (no code changes to pkdgrav3 beyond parameters).
   - `nfsstat -m` (version, rsize/wsize, nconnect, sync/async). fio on the NFS directory: bs=4M, numjobs=1,4,16,64, report MB/s. This sets the realistic ceiling: speedup is capped at ceiling / 30 MB/s.
   - Time single-thread `gzip -1` on a float32 dump to test whether gzip alone explains 30 MB/s.
   - 512^3 baseline run: write to NFS and to `/dev/shm`, compare, to separate software cost from storage cost.
   - `h5cc -showconfig | grep -i thread`; find out which patch adds gzip.
   Exit: written table of ceilings and the identified bottleneck. (Results: section 10.)
1. Offline prototype (Python/numba) on real 512^3 Gadget snapshots at 3 redshifts: bytes per particle per stream (positions, velocities, IDs), max and rms error, Option A vs Option B, ID coding variants. Verify the decoder against `SnapshotReader`.
   Exit: chosen scheme, bytes per particle, error no worse than the current pipeline.
2. C/C++ codec library with no pkdgrav3 dependency (quantizers, bit-packers, ID coder, block index) plus round-trip unit tests.
3. Integration in the pkdgrav3 fork: new parameter (for example `bCompressedOutput`), new leaf-thread service called from `MSR::Write`, bypassing FIO and `Reorder()`. Hook points: master.cxx (`Write`/`AllNodeWrite`), pst.c (`pstWrite` branch), pkd.c (new `pkdWriteCompressed`). Keep the diff to existing files small.
   Exit: 512^3 and 1024^3 round trip on pawna (`-np 1`) and pegasus (MPI).
4. Benchmark 2048^3: GB/s, stall per snapshot, peak RSS delta, bytes per snapshot.
   Exit: at least 20x less stall or NFS ceiling reached; RSS delta under 1 percent.
5. Reader: Python class in sahyadri-sandbox (next to `SnapshotReader`), converter to the format Rockstar reads, optional background drain/staging.

## 6. Repository layout
- Fork wullm/pkdgrav3 (or whichever tree is actually run) on a feature branch. All new code in new files/directory (for example `io/snapcomp/`); minimal edits to existing files.
- Keep the codec free of pkdgrav3 headers so it can be split into its own repo later.
- Python reader and format spec in sahyadri-sandbox.

## 7. Open items
- Which pkdgrav3 tree and build carries gzip level 1 in HDF5?
- How exactly does `bParaWrite` fail (crash, hang, corrupt file)?
- Are IDs Lagrangian grid indices (id = (i*N+j)*N+k)?
- Does pawna have node-local scratch disk?
- Accuracy target: keep current uint8 cell offset + 12-bit velocity and current Ngrid?

## 8. Notes after rebasing onto pkdgrav3 3.5.1 (2026-11-06)
Section 3 was read from wullm 7a5041e (3.1.2). The branch is now based on upstream v3.5.1 (dcce745f), the release production runs; section 3 must be re-checked against it.
- Production binary: /user/swarajit/apps_user/pkdgrav3_openmpi/bin/pkdgrav3, built from /user/swarajit/src/pkdgrav3-3.5.1, linked to conda env pkdgrav_no_mpi (HDF5 2.1.0, threadsafe ON, parallel OFF) and /apps/openmpi-5.0.1; run with `mpirun -np 1`.
- gzip source (open item 1): 3.5.1 hardcodes shuffle + H5Pset_deflate(level 9) + Fletcher32 (io/fio.c:2535); the production job LD_PRELOADs libh5shim.so with PKD_H5_DEFLATE=1 to force level 1.
- HDF5 chunks are (65536, 1): every column is chunked and compressed separately. Density and Potential columns are written even with bDoPotOutput=0.
- Mac, 256^3, 12 threads, local SSD, pkdgrav3 "Output file ... Wallclock" line: HDF5 (gzip 9) 66 s/snapshot, 326 MB; tipsy 2.9 s, 604 MB. Software cost dominates on fast storage.
- Production 2048^3 snapshot on disk: 166 GB (gzip 1), not 275 GB.
- Restart "second pass": with `from PKDGRAV import restore` in the restart file, pkdgrav3 does not see the PKDGRAV module in the script's namespace, so after `restore()` returns it generates ICs and runs the whole simulation again (msrpython.cxx `is_PKDGRAV_imported`; main.cxx). `import PKDGRAV; PKDGRAV.restore(...)` avoids it (tested on the Mac, 32^3). This is what the production guard scripts detect and kill.
- A checkpoint stores `achIoPath`; a restart reuses it unless `restore(..., achIoPath=...)` overrides it (tested).

## 9. Storage on pawna and node-local staging (2026-11-06)
Probe job on node4 (inside PBS) plus login-node probes:
- `/storage`: NFS 4.2 from 10.20.255.85, rsize/wsize 1 MiB, hard, TCP, no nconnect; 421 TB, 216 TB free. Production output goes here.
- Home `/user`: ext4 on the login node (164 TB, 36 TB free), NFS 4.2 to compute nodes (10.20.255.80). Same protocol, smaller server, less space: not a better target.
- Node-local: each compute node has a 1.5 TB SSD array (MegaRAID MR9540-8i) as `/`, about 1.4 TB free; `/tmp`, `/var/tmp` and PBS's `$TMPDIR=/var/tmp/pbs.<jobid>` are writable. `/dev/shm` is 378 GB of RAM (too valuable at 2048^3).
- Write ceilings, measured with parallel `dd` (bs=4M, incompressible data from a /dev/shm file, `conv=fdatasync`, one file per stream, 64 GB in total per row, wall time from start to the last fdatasync; fio is not installed). Node4, 2026-11-06:
  - node-local `/var/tmp`: 427 / 797 / 1458 / 1035 MB/s for 1 / 4 / 16 / 64 streams (job 2536379).
  - NFS `/storage`: 28 MB/s (1 stream), 25 MB/s (4 streams) (job 2536378; a pawna build was downloading to /storage at the same time). More streams do not help: the total is capped at about 25-30 MB/s, which by itself explains production's 30 MB/s.
  - Links are 10 GbE (1.25 GB/s), so the cap is the NFS server side or other users' load. Ask the admins: sync export? disks? nconnect?
  - NFS `/storage`, 8 GB per row (jobs 2536391): 28 MB/s (16 streams), 33 MB/s (64 streams).
  - Home `/user` (NFS from the login node), 8 GB per row (jobs 2536390/2536391): 270 / 325 / 379-395 / 426-488 MB/s for 1 / 4 / 16 / 64 streams. 8 GB may partly sit in the login node's cache if the export is async, so treat this as an upper bound; still about 10x `/storage`.
- Consequence: a single 166 GB snapshot needs about 1.5-1.8 h on NFS whatever we do in pkdgrav3. Staging removes it from the simulation's critical path; compression shortens the background copy. The 20x stall target can only be met with staging (or a faster NFS).

Staging design (decided: separate drain process first; in-process copy only if measurements show interference):
- Write side needs no code change. pkdgrav3 3.5.1 sends snapshots to `achIoPath` while logs and checkpoints follow `achOutPath`/`achCheckpointPath`. Set `achIoPath = "<stage>/current/{name}.{step:05d}{type}"` with `<stage>` on local disk; checkpoints stay on NFS so restarts never depend on a node's disk.
- `tools/staging/pkd_drain.py` (stdlib Python) runs next to pkdgrav3 in the same PBS job. A snapshot counts as finished when the log shows "Writing <path> ..." followed by "Output file has been successfully written", or when `<path>.done` exists (our compressed writer will create one). It copies every piece to NFS as `.part` in 16 MiB chunks with an optional MB/s limit, fsyncs, checks the size, renames into place, fsyncs the directory, then deletes the local copy.
- Disk guard: `<stage>/current` is a symlink to `<stage>/local`; if local free space drops below `--min-free-gb` it is atomically switched to the NFS destination (later snapshots go straight to NFS) and switched back when space recovers.
- End of job: the drain exits when pkdgrav3 has ended and everything finished is copied (exit 2 if anything is left). The job script must leave time for this: set pkdgrav3's `iWallRunTime` below the PBS walltime by at least one snapshot's copy time.
- Run the drain with `nice -n 19 ionice -c 3`; it uses about one core and NFS bandwidth while 191 cores compute.
- Works with the production binary today and with the compressed format later (multi-file pieces + `.done`).
- Tested on the Mac (64^3, 4 HDF5 snapshots): drain after the run and during it; all copies valid; rate limit and space switch triggered as designed. `--recover` drains leftovers after a crash.
- Usage, job script and examples: [running_on_pawna.md](running_on_pawna.md) (`tools/staging/run_staged.pbs`).
- Risks: node failure loses snapshots not yet drained (checkpoints are safe); the copy competes with the run's MPI/NFS traffic (to be measured at 512^3 on pawna); PBS kill at walltime can interrupt a copy (the `.part` file is left, the local original is kept).
- In-process alternative (later, only if needed): a low-priority pthread in pkdgrav3 that copies between steps; it knows exactly when writes finish and can pause during communication-heavy phases, but needs code in pkdgrav3 and shares the process's failure modes.

## 10. Phase 0 results (2026-11-07)
512^3 (134M particles), our 3.5.1 build, 64 threads, `mpirun -np 1`, node3/node4. Time = pkdgrav3's
"Output file has been successfully written, Wallclock" per snapshot, which is the simulation stall.
HDF5 = stock 3.5.1 (shuffle + gzip 9 + Fletcher32), 2.59 GB per snapshot; Tipsy = raw, 4.83 GB.

| Snapshot written to | HDF5 gzip 9 | Tipsy |
|---|---|---|
| NFS `/storage` | 370 s | 183-209 s (23-26 MB/s) |
| `/dev/shm` (RAM) | 277-287 s | 18.2 s (265 MB/s) |
| node-local SSD + drain (job 2536385) | 279-290 s | not run |

- Staging works at 512^3: stall equals the RAM case; the drain copied each 2.41 GiB snapshot in 92-101 s (25-27 MB/s, the NFS ceiling) and finished before the next snapshot; all 4 drained, local disk cleaned, peak RSS unchanged (24.7 GB).
- Two bottlenecks, both large:
  1. **Software**: the single writer thread. Even raw Tipsy into RAM is 265 MB/s; HDF5 with gzip 9 is 15x slower (about 19 MB/s of raw data). At 2048^3 (64x) that is about 20 min (Tipsy) and 5 h (gzip 9) of stall from software alone.
  2. **NFS `/storage`**: 25-33 MB/s in total for 1-64 streams; Tipsy to NFS runs exactly at that rate. Links are 10 GbE, so the limit is the server or its load.
- Storage ceilings: node-local SSD 0.4-1.5 GB/s (best at 16 streams); home 0.27-0.49 GB/s (upper bound); `/storage` 0.03 GB/s.
- Implication: per-thread parallel compression (Phases 1-3) removes bottleneck 1; staging removes bottleneck 2 from the stall; compression shortens the background copy. Home is a possible faster destination (36 TB free; shared) to discuss.
- Not yet measured: production gzip level 1 at 512^3 (`H5LEVEL=1`), standalone `gzip -1` speed on a float32 dump.

## 11. Phase 1 results (2026-11-08, offline prototype `tools/snapcomp_proto/`)
Data: our own pkdgrav3 3.5.1 runs on the Mac, 128^3, production cosmology and accuracy settings,
2LPT ICs at z=49, snapshots at z=3, 1, 0, in two boxes: 200 Mpc/h (production box, realistic large-scale
flows) and 12.5 Mpc/h (production mass resolution, realistic halo interiors).

Findings about the data:
- IDs are Lagrangian lattice indices, id = (iz*N + iy)*N + ix (rms displacement 2.1 cells at z=3 with this
  mapping vs 7.4 for others). Open item closed.
- Files are in pkdgrav3 memory order, not ID order (no Reorder() in this path).
- HDF5 `Time` is code time; a = 1/(1+z). Velocity conversion verified against Zel'dovich (slope 1.02).

Quantizers (identical for A and B, so errors are identical): positions with step L/(256 N) (the uint8 offset
in a mean-spacing cell of the Gadget scheme), velocities with step dv = 1 km/s (max error 0.5 km/s).

Bytes per particle, zstd level 3, all round trips exact:

| Box, z | A: blocks of 2^18 in file order | B: ID order, Lorenzo, implicit IDs | B without prediction |
|---|---|---|---|
| 200, z=3 | 7.65 (pos 3.35, vel 3.64, id 0.66) | 5.95 | 7.24 |
| 200, z=1 | 8.16 (3.34, 3.89, 0.93) | 7.45 | 7.84 |
| 200, z=0 | 8.38 (3.29, 4.01, 1.08) | 8.00 | 8.25 |
| 12.5, z=3 | 7.43 (3.30, 3.02, 1.10) | 6.44 | 7.32 |
| 12.5, z=1 | 7.81 (3.22, 3.27, 1.32) | 7.46 | 8.01 |
| 12.5, z=0 | 7.91 (3.17, 3.37, 1.38) | 8.01 | 8.37 |

References: raw float32 + uint64 = 32 B; stock HDF5 (gzip 9) 19.8-20.0 B; lossless float + zstd 20.4-20.6 B;
production 2048^3 gzip-1 file: 166 GB / 2048^3 = 19.3 B.

Sensitivity (Option A, z=0): each halving of dv adds 0.35-0.4 B; each doubling of position resolution adds
0.4-0.5 B. dv = 0.5 / 1 / 2 / 4 km/s -> 8.72 / 8.38 / 8.03 / 7.61 B (200 Mpc/h box); sub = 128 / 256 / 512 /
1024 -> 8.01 / 8.38 / 8.92 / 9.30 B.

ID coding (Option A): sorting each cell's particles by ID and storing ID differences within a cell
(`delta_in_cell`) wins at z <= 1 (0.9-1.4 B); Lagrangian-minus-Eulerian cell residuals win at z=3 in the big
box (0.66 B); raw uint64 + zstd costs 1.0-1.8 B.

Encode speed: zstd-3 runs at 500-1300 MB/s of quantized input per core, i.e. tens of millions of particles
per second per thread; compression will not be the bottleneck in situ.

Recommendation: **Option A**. It meets the 3x target (3.8-4.3x vs raw, about 2.4x smaller than today's
gzip-1 HDF5) and fits the per-thread design. Option B is better only at high redshift (by up to 1.7 B at
z=3, nothing at z=0) and needs particles in global ID order, i.e. the gather/Reorder() we want to remove.
Expected 2048^3 snapshot: about 8 B x 8.6e9 = 69 GB (vs 166 GB today); 41 min to drain at 28 MB/s.

To verify before Phase 2: (1) ID cost at larger N (displacements in cell units grow with N at fixed box;
test 256^3 on the Mac and 512^3 on pawna); (2) the accuracy target: is dv = 1 km/s and L/(256 N) at least as
good as the current Gadget 12-bit pipeline? (open item); (3) block size (2^16-2^20) and zstd level.

## 12. Subsamples and user-chosen accuracy (requirements for Phases 2-3, 2026-11-08)
- Accuracy is chosen by the user in the run's parameter file (position step as steps per mean
  inter-particle spacing, velocity step in km/s, block size, zstd level), with defaults
  sub = 256, dv = 1 km/s. The values are stored in each file header so readers need nothing else.
- Subsamples as in Sahyadri: disjoint percentages that sum to 100 (default 1, 3, 5, 9, 10, 12, 15, 20,
  25), one file set per subsample, so users download only what they need and the large ones can be
  deleted. Sahyadri's `create_subsamples()` (sahyadri-sandbox `scripts/post-process/Nbody_Compression.py`)
  shuffles the array index with `np.random.seed(0)`; that needs a global permutation (about 69 GB at
  2048^3) and cannot be done in situ. In situ, each particle's subsample is decided from a hash of its
  ID (splitmix64 with a seed), so every thread decides alone and membership is the same in every
  snapshot. Memberships differ from existing Sahyadri files.
- Cost (prototype `run_subsamples.py`, Option A, 128^3, z=0): the 9 subsample files together take 10%
  more than one file (9.19 vs 8.38 B/particle in the 200 Mpc/h box; 8.69 vs 7.91 at 12.5 Mpc/h). The
  1% file costs about 10 B/particle, the 25% file about 9. All subsample round trips exact.
- Sahyadri for reference (arXiv 2601.07924): 400^3 grid, 8-bit offsets (dx = L/102400 = 1.95 kpc/h),
  12-bit counts with overflow, velocities 12 bit over +-5000 km/s (dv = 2.44 km/s), IDs uncompressed.
  Our default (L/(256*2048) = 0.38 kpc/h, 1 km/s) is 5x finer in position and 2.4x finer in velocity.
- Option A at exactly Sahyadri's accuracy (dx = 1.95 kpc/h, i.e. 50 steps per mean spacing at production
  mass resolution; dv = 2.44 km/s), 12.5 Mpc/h box, zstd-3, IDs lossless:
  z=3 6.03, z=1 6.41, z=0 6.51 B/particle in one file (pos 2.31, vel 2.82, id 1.38 at z=0);
  the 9 subsample files together 7.28 B/particle at z=0 (+12%).
  Sahyadri's files (paper section 2.1.2: 56 GB without IDs, 85 GB with, for 2048^3) are about 6.5 + 3.4 =
  9.9 B/particle. Same accuracy, all subsamples: about 26% smaller, with lossless IDs and no +-5000 km/s
  velocity clip. Projected 2048^3 snapshot with subsamples: about 63 GB.
- Parameters (Phase 3) should take the position step as a length (kpc/h) and the velocity step in km/s,
  so presets are easy: "sahyadri" = 1.95 kpc/h, 2.44 km/s.
- Note for the paper: section 2.1.2 says "ten subsamples" but lists nine percentages (they sum to 100).

### Phase 1 exit (2026-11-09)
Block size (Option A, Sahyadri accuracy, z=0, 12.5 / 200 Mpc/h): 4096: 6.61 / 7.28; 16384: 6.55 / 7.21;
65536: 6.59 / 7.19; 262144: 6.51 / 7.07; 1048576: 6.50 / 7.03 B/particle. Default 65536 (+1-2% vs 2^18,
4x less scratch memory, about 10 MB per thread, 2 GB for 192 threads).

Decisions (user, 2026-11-08): subsample membership from an ID hash; default accuracy = Sahyadri
(1.95 kpc/h, 2.44 km/s), overridable in the parameter file.

Exit criteria:
- Chosen scheme: Option A (per-thread blocks, cell sort, 8-bit-style offsets, quantized velocities,
  per-block ID coding, zstd). Met.
- Bytes per particle: 6.0-6.5 at Sahyadri accuracy, 7.3 with the 9 subsample files (z=0), vs about 9.9
  for the current Sahyadri files and 19.3 for production gzip-1 HDF5. Met.
- Error no worse than the current pipeline: identical quantization steps to Sahyadri (max 0.98 kpc/h,
  1.22 km/s), lossless IDs, no velocity clipping. Met.
- Not done: decoder check against `SnapshotReader` (our format differs; the Python reader is Phase 5),
  and the ID cost at N > 128 (displacements in cell units grow with N). Carried into Phase 2 tests:
  the C++ codec will be run on a 512^3 pawna snapshot.

## 13. Phase 2: codec library (2026-11-09)
`io/snapcomp/` (C++17, zstd, POSIX; no pkdgrav3 headers), format in [snapcomp_format.md](snapcomp_format.md).
- `BlockEncoder` (one per thread, reusable buffers): quantize, sort by (cell, ID), streams for cell keys,
  counts, offsets, velocities and IDs (raw / delta in cell / Lagrangian residual, smallest kept), zstd with
  checksums. `decode_block` checks every size and checksum and throws on corruption.
- `FileWriter`: one file shared by all threads of a process; `append()` reserves the offset with an atomic
  counter and writes with `pwrite` (no lock on the data path); `close()` writes the block index, trailer
  and final header, then fsync. `FileReader` rejects files without a trailer (writer crashed).
- `SubsampleTable`: Sahyadri percentages by splitmix64 hash of the ID; identical to the Python version.
- `choose_grid`: cells of the mean particle spacing, steps per cell from the requested position step.
- Standalone build: `cmake -S io/snapcomp -B build`; uses system zstd or fetches zstd 1.5.6.

Tests (Mac, `snapcomp_test`, all pass): subsample hash vs numpy golden values (10 IDs + counts over 10^6
IDs); header round trip; exact round trips with error <= half a step for lattice particles, random 64-bit
IDs, all particles in one cell, box edges, sub = 1, odd block sizes; 4/4 byte corruptions and a truncated
block detected; 8 threads appending to one file, every particle read back once; truncated file rejected.

Real snapshots (Mac, `pkc_convert --verify`, 128^3, Sahyadri accuracy, 8 threads, 9 subsample files):
all particles back exactly once, max errors 0.98 kpc/h and 1.22 km/s (half a step) for all 6 snapshots.
12.5 Mpc/h box: 6.95 / 7.32 / 7.40 B/particle at z = 3 / 1 / 0 (Python prototype with 2^18 blocks over the
whole sample: 7.28 at z = 0; the difference is per-thread partial blocks for the small subsamples).
Single file: 6.599 B/particle, identical to the prototype (6.59).
Independent Python reader (`tools/snapcomp_proto/pkc_reader.py`, from the spec only) decodes the C++ files
exactly: IDs once, quantized positions identical, errors half a step, membership equal to the hash.

Speed (Mac M4 Pro, one thread, block 65536, zstd-3): encode 8-10 M particles/s, decode 44-46 M particles/s.
8 threads: 56 M particles/s for one file (including all I/O). At 2048^3 on 192 threads: about 5 s of
compression per snapshot.
