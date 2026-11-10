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

Pawna (job 2536419, node3, 2x EPYC 9654, gcc 11, zstd 1.5.6 fetched): all unit tests pass; one thread
encodes 5.7 M particles/s and decodes 25 M/s. A 512^3 run at production mass resolution (50 Mpc/h, z = 49
-> 0, 64 threads, 3 h 12 min, 26 GB peak RSS) gave HDF5 snapshots at z = 3, 1, 0, converted with 64
threads to local SSD and verified (every particle once, max errors 0.98 kpc/h and 1.22 km/s):

| 512^3, B/particle | z = 3 | z = 1 | z = 0 |
|---|---|---|---|
| 9 subsample files, Sahyadri accuracy | 7.34 | 7.77 | 7.95 |
| single file, Sahyadri accuracy | 6.35 | 6.89 | 7.14 |
| single file, 0.38 kpc/h and 1 km/s | 7.73 | 8.27 | 8.55 |

Streams at z = 0 (9 files): pos 2.57, vel 3.17, id 2.22 B/particle (128^3 Mac: 2.52, 2.91, 1.94). IDs cost
0.28 B more at 4x N per side, velocities 0.25 B more (stronger large-scale flows in the larger box).
Expected at 2048^3 (another 4x): about 8.5-9 B/particle with subsamples, i.e. 73-77 GB per snapshot
(Sahyadri format: 85 GB; production HDF5: 166 GB).
Throughput: 64 threads compressed, split into 9 subsamples and wrote 134M particles in 1.6-2.2 s
(61-84 M particles/s, 4.1-4.6 M particles/s per thread). Scaled to 2048^3 on 192 threads: about 45 s
per snapshot, against 2.5 h today.

### Phase 2 exit (2026-11-09)
Codec library with no pkdgrav3 dependency (quantizers, bit/byte packing, three ID coders, block index,
subsamples), file writer for many threads, reader, format spec, unit tests on Mac and pawna, independent
Python reader, and verified round trips on real 128^3 and 512^3 pkdgrav3 snapshots. Done.
Kept for Phase 3/5 work: 512^3 z = 0 subsample files in /storage/shadab/pkdgrav3_io/bench/p2_pkc (1 GB).

## 14. Phase 3: integration in pkdgrav3 (2026-11-09)
`bCompressedOutput = True` makes every snapshot (not checkpoints) go through `MSR::WriteCompressed`
(`msr_writecompressed.cxx`) and `ServiceWriteCompressed` (`core/writecompressed.{h,cxx}`): all leaf
threads run at once, each compresses its own particles (one pass per subsample, one block of scratch) and
appends to one file per subsample shared by the threads of its process; the last thread closes each file and,
after the last one, writes `<name>.done`. No gathering, no `Reorder()`, no HDF5. File names follow
`BuildIoName`, so `achIoPath` and staging work unchanged. Parameters: `dCompPosStep` (1.95 kpc/h),
`dCompVelStep` (2.44 km/s), `lstCompSubsamples` (Sahyadri), `iCompSubsampleSeed`, `iCompBlockSize`,
`iCompLevel`, and `bCompAlsoStandard` (also write HDF5/Tipsy from the same state, for validation).
Edits to existing files: one branch in `MSR::Write`, one declaration in `master.h`, one service registration
in `main.cxx`, CMake (library + two sources; reuses cpptrace's zstd), and the new parameter entries.

Tests:
- Mac, 64^3 (6.25 Mpc/h), 8 threads, `bCompAlsoStandard`: snapshots at z = 3, 1, 0 checked against the HDF5
  of the same step with `compare_pkc_hdf5.py` (spec-only reader): every ID once, correct subsample, errors
  half a step. 9 files: 6.84 / 7.12 / 7.20 B/particle; one file: 5.90 / 6.21 / 6.29.
- Mac, MPI, 2 processes x 4 threads: 18 files per snapshot (2 writers), all checks pass.
- Mac, staged to a local directory with the drain: all files drained via the `.done` marker (verbose off)
  and via log lines (verbose on). Bug found on pawna and fixed: with both triggers the drain left the marker
  behind and exited 2 (data was complete).
- pawna, job 2536426, 512^3 (50 Mpc/h), `-np 1`, 64 threads, snapshots at z = 20, 10, 5:

| | z = 20 | z = 10 | z = 5 |
|---|---|---|---|
| A: compressed write (stall), local SSD | 1.44 s | 1.90 s | 1.35 s |
| A: HDF5 (stock gzip 9), same run, same state, local SSD | 250.6 s | 256.6 s | 281.9 s |
| A: compare with that HDF5 | OK | OK | OK |
| B: compressed, staged (stall) | 2.71 s | 1.90 s | (similar) |
| B: drain to NFS while the run continues | 54 s | 28 s | 23 s |
| size (9 subsample files) | 5.89 B | 6.43 B | 6.98 B |

  130-200x less stall than HDF5 on the same disk; peak RSS 26.2-26.7 GB (Phase 2 run without the writer:
  26.0 GB). Drained snapshots read back complete (134,217,728 unique IDs each).
- IDs are persistent: same ID set at z = 3 and z = 0, same-ID displacement median 1.4 Mpc/h (12.5 box) /
  5.4 Mpc/h (200 box) vs 6.0 / 98.5 for random pairs; stored losslessly; subsample membership depends on the
  ID only, so a particle stays in its subsample in every snapshot (consistent-trees works per subsample).

Exit criterion status: 512^3 on pawna `-np 1` done; MPI tested on the Mac (2 processes). Not done: 1024^3
(needs approval, above the 512^3 limit) and pegasus (no access from this setup).

## 15. Phase 4: full-size tests with persistent IDs at no memory cost (2026-11-09)
Setup: a fresh root on pawna (`/storage/shadab/pkdgrav3_io/fresh`) built from GitHub with
`tools/pawna/setup_pawna.sh` exactly as in docs/running_on_pawna.md section 2a (setup and smoke test passed).
Runs: `bench/phase4/p4.pbs` with `p4.par.in`: production mass resolution (L = 200 N/2048 Mpc/h), production
cosmology/accuracy, 2LPT ICs from pkdgrav3 at z = 49, 6 steps, snapshots at z = 45 / 25 / 12, **ordered mode
(`bMemUnordered = False`, no `bMemParticleID`: ID = iOrder)**, compressed output (Sahyadri accuracy, 9
subsamples). Part A (N <= 512): compressed + HDF5 from the same state, exact comparison. Part B: staged on local
SSD, drained to NFS, memory sampled every 15 s, `pkc_check` on every drained snapshot (+ persistence vs z = 45).

| | 128^3 (32 thr, job 2536429) | 512^3 (64 thr, job 2536430) |
|---|---|---|
| A: compressed vs HDF5, same state | OK (3/3) | OK (3/3), max err 0.9773 kpc/h incl. float32 of HDF5 |
| A: stall compressed / HDF5 (gzip 9) | 0.04-0.55 s / 4.0-4.4 s | 1.56-1.69 s / 216-230 s |
| B: stall (staged) | 0.04 s | 1.33-1.66 s |
| B: size (9 files), z = 45 / 25 / 12 | 4.88 / 5.30 / 5.84 B | 5.31 / 5.72 / 6.28 B |
| B: drain to NFS | 7-8 s | 19-24 s (31-35 MB/s) |
| peak RSS (time -v) | 5.8 GB | 25.0 GB |
| writer scratch | 6.2 MB/thread | 14.0 MB/thread (0.88 GB total) |
| IDs once, right subsample (all 3) | yes, 2,097,152 | yes, 134,217,728 |
| z = 45: max distance from own lattice point | 1.81 cells | 4.09 cells (random IDs: ~N/2) |
| 1% subsample same IDs at z = 45, 25, 12 | yes (21,250) | yes (1,342,802) |
| same-ID displacement z = 45 -> 12, max | 0.50 Mpc/h | 1.11 Mpc/h (random pairs ~24) |

No swap was added by any run (the node already had 15-159 GB of swap in use by other jobs).
### 2048^3, first attempt (job 2536431, node2 to itself, 192 threads, 2026-11-10)
Measured with the job's 15 s sampler (`ps -o rss` of pkdgrav3, `free -m` for the node) and pkdgrav3's own
timers in the run log.
- Peak pkdgrav3 RSS **746.0 GiB** (ordered mode, persistent IDs = iOrder). The run without persistent IDs
  measured 742 GB earlier and the run with bMemParticleID 870 GB: persistent IDs now cost about nothing.
- The node has 755 GiB, so the run is at the edge: about 1.5-2 GB free from the first gravity step on;
  swap grew to 25-27 GB and stayed flat. All 192 cores stayed busy (PBS cpupercent 18933), i.e. no thrashing.
- Timings: IC (2LPT) 350 s; first gravity 603 s; later gravity 328-604 s; domain decomposition 27-51 s;
  tree 14-39 s.
- At the first snapshot (z = 45) swap jumped by 7 GB within 15 s while the compressed files were written
  to the local SSD, and the job's swap guard (stop if swap grows by more than 30 GB) killed pkdgrav3 in the
  middle of the write. No snapshot was completed, so no ID checks yet.
- Cause: the written data sat in the page cache, and on a full node the kernel swapped out simulation
  memory to keep it. Fix: the writer now pushes each block to disk and drops it from the page cache
  (`sync_file_range` + `posix_fadvise(DONTNEED)`, Linux only), so writing a snapshot needs no extra RAM.

### 2048^3, second attempt (job 2536432, node4 to itself, 2026-11-10), with the page-cache fix
Same job, swap guard raised to 60 GB of growth. Measured as above; the ID check ran as a separate job
(`bench/phase4/check.pbs`, 32 threads, reading the drained files from NFS).
- RSS 744 GiB from the first gravity step on, swap 25 GB (as in the first attempt).
- Snapshot z = 45 (step 1): **stall 50.0 s** (production HDF5 ~2.5 h); 45,624 MB = **5.31 B/particle**
  for 8.59e9 particles with 9 subsample files; writer scratch 2.14 GB over 192 threads (11.5 MB per
  thread, 0.29% of the RSS). Encode 2258 s and append (write + flush to SSD) 3690 s summed over threads,
  i.e. the stall is set by the local SSD at about 0.9 GB/s. No swap jump during the write.
- Drain to NFS: 42.5 GiB in 2003 s (22 MB/s). Swap grew from 25 to 41 GB while it ran: the drain's
  reads and NFS writes filled the page cache. Fixed afterwards in `pkd_drain.py` (fadvise DONTNEED on
  both files, fdatasync every 256 MiB); the running job still used the old drain.
- `pkc_check` z = 45: 8,589,934,592 IDs read, all distinct, none missing or out of range, all in the
  right subsample; distance from own lattice point rms 2.52 cells, max 6.91 (random IDs: ~1024).
  4 min 50 s, 1.3 GB RSS.
- Side effect found with the 512^3 run of section 16 (straight to NFS, plenty of RAM): unconditional
  write-through made the stall the NFS time (30-35 s instead of 1.4-1.9 s). Now a parameter,
  `iCompWriteThrough` (-1 auto, 0, 1); auto uses it only when MemAvailable < 2x the process's snapshot
  (10 B/particle), i.e. on for 2048^3 on one node, off for 512^3.
  Regression check (job 2536435, commit 7b2fb67f, `bench/phase4/p4.pbs` with N=512, 64 threads, same
  as the earlier 512^3 Phase 4 run): the run log prints "write-through: 0 of 63 threads" (node
  MemAvailable 189-235 GiB) at all three snapshots; stalls 1.43-1.52 s straight to NFS (A) and staged
  (B); HDF5 219-233 s; 5.31/5.72/6.28 B/particle; max error within half a step; `pkc_check` all
  134,217,728 IDs once, 1% set identical in all three snapshots; RSS 23.9 GiB; drain 25-29 MB/s.
- Whole job (`p4.pbs`, /usr/bin/time -v and the 15 s sampler; pkc_check by the job after the run, 192
  threads, reading the drained files on NFS): pkdgrav3 rc=0, 12 h 50 min wall for 6 steps (z = 49 -> 12).
  **Peak RSS 746.7 GiB** (782,991,992 kB; sampler 746.5 GB); node swap 42.7 GB at the peak (25 GB of it
  before the first snapshot, the rest from the old drain, see above).

  | snapshot | z | stall | size | B/particle | encode / append (thread-s) | drain to NFS |
  |---|---|---|---|---|---|---|
  | snap.00001 | 45 | 50.0 s | 45.6 GB | 5.311 | 2258 / 3690 | 42.5 GiB, 2003 s, 22 MB/s |
  | snap.00003 | 25 | 54.2 s | 49.3 GB | 5.734 | 2249 / 3908 | 45.9 GiB, 2077 s, 23 MB/s |
  | snap.00006 | 12 | 58.2 s | 54.1 GB | 6.302 | 2222 / 4580 | 50.4 GiB, 2517 s, 21 MB/s |

  Writer scratch 2.14-2.15 GB over 192 threads (0.29% of RSS) at every snapshot.
- IDs: in all three snapshots 8,589,934,592 IDs read, all distinct, none missing or out of range, all in
  the right subsample. Distance from own lattice point rms 2.52 / 4.47 / 8.95 cells (random: ~1024).
  Persistence of the 1% subsample against z = 45: 85,912,182 IDs there and here, all found, same set; same-ID
  displacement median 0.17 Mpc/h (z = 25) and 0.56 Mpc/h (z = 12), max 1.9 (random pairs: ~96).
- **Phase 4 exit met:** stall 50-58 s vs ~1.6-2.5 h for production HDF5 (> 100x), 5.3-6.3 B/particle vs
  19.3 (production gzip-1 HDF5) and 32 (raw float32 + uint64 ID), writer RAM 0.29%, persistent IDs at no
  memory cost (746.7 GiB vs 742 GB production without IDs and 870 GB with bMemParticleID).

## 16. Phase 5: readers (2026-11-10)
- Python: `tools/pkcreader/pkcreader.py` (numpy + zstandard, from the spec only, one file).
  `PkcSnapshot(prefix)`: header values (redshift, scale, Lbox, npart, mpart, Om, Olam, hubble, percents),
  `read(fields, subsamples, units, dtype)` decoding blocks on threads into preallocated arrays,
  `read_block()` in SnapshotReader style ((3, n) for pos/vel), `iter_blocks()` for streaming. Files are
  memory-mapped. The prototype's `pkc_reader.py` is now a shim to it.
- Rockstar, stopgap: `io/snapcomp/tools/pkc_to_gadget` (C++, threads) writes Gadget-2 files
  (float32 Mpc/h in [0, L), v_pec/sqrt(a), 8-byte IDs, mass in 1e10 Msun/h from Omega_m L^3 / n, at most
  2^24 particles per file so every record stays below 2^31 bytes) and prints the Rockstar settings.
- Tests (Mac, `tools/pkcreader/test_pkcreader.py`, 64^3 snapshot at z = 3 with HDF5 written from the same
  state): every ID once and found in HDF5; position error 0.9766 kpc/h and velocity error 1.22 km/s
  (half a step); read / read_block / iter_blocks give the same particles in the same order; 1% subsample
  is exactly the hashed set; Gadget files (1 and 6 files): totals, header, mass, IDs, positions and
  v*sqrt(a) equal to the reader's. All pass.
- pawna, 512^3 in 200 Mpc/h to z = 0 (job 2536434, `bench/phase5/rs512.pbs`, 64 threads, 1 h 46 min,
  RSS 25.0 GiB), compressed + HDF5 of the same states at z = 3, 2, 1, 0.5, 0, straight to NFS:
  8.25 / 8.47 / 8.71 / 8.82 / 8.89 B/particle; HDF5 write 384-443 s; compressed stall 30-35 s (with
  the unconditional write-through of that build, see section 15). `pkc_check`: all 134,217,728 IDs once
  and in the right subsample at every redshift; the 1% set (1,342,802 IDs) is the same at all five;
  displacement since z = 3 at z = 0: median 5.4, max 20.2 Mpc/h (random pairs ~96).
  `test_pkcreader.py` at z = 0: all checks pass (position error 0.9794 kpc/h = half step + float32
  rounding of the HDF5 file at L = 200; velocity 1.22 km/s); `pkc_to_gadget`: 9 files in 1.8 s
  (64 threads, 4.3 GB RSS), read back equal. Snapshots kept in `fresh/bench/rs512/` for Rockstar.
- Not done here: running Rockstar on converted files. Building Rockstar from the downloaded upstream
  source on the Mac was refused by the Claude Code permission check; needs the user's approval or the
  pawna build. Native Rockstar reading is planned in a separate session:
  [rockstar_pkc_session_prompt.md](rockstar_pkc_session_prompt.md).
- Sahyadri-sandbox: the reader is meant to sit next to `scripts/post-process/readers.py`; not copied
  there yet (that repo belongs to a-paranjape; needs a branch and the owner's agreement).

## 17. Potential in the compressed format: status and cost estimate (2026-11-10)
- Status: not written. Format v1 holds positions, velocities and IDs only. HDF5 snapshots always carry
  `Potential` (io/fio.c, FIO_FLAG_POTENTIAL); `bDoPotOutput` adds a separate `.pot` array file (text).
  Production has one, `snap.02020.pot`, 111.8 GB (13 B/particle).
- Computing it costs nothing extra: the gravity kernels always accumulate `fPot` (it feeds the energy
  diagnostics, gravity/grav2.cxx), and with gravity on `bMemPotential` defaults to 1 (simulate.cxx), so
  every run already stores it, 4 B of the 48-byte particle. The 2048^3 log shows `POTENTIAL` in the
  memory model, i.e. it is already inside the 742-746 GiB. Snapshots are written at top-level step
  boundaries, where gravity has just been computed for all rungs, so the stored value is current for
  every particle (from reading the code; to be checked when implemented).
- Bytes (Mac, `tools/snapcomp_proto/potcomp.py`, 128^3 HDF5 snapshots of section 11): potential in
  peculiar (km/s)^2, spatial order (Morton key on 1024^3, blocks of 2^18), quantised with step dphi,
  delta, zigzag, byte shuffle, zstd 3:

  | dphi (km/s)^2 | L=12.5, z=3 | L=12.5, z=0 | L=200, z=0 |
  |---|---|---|---|
  | 1    | 1.62 | 1.66 | 2.36 |
  | 10   | 1.26 | 1.29 | 1.88 |
  | 100  | 0.78 | 0.83 | 1.50 |
  | 1000 | 0.43 | 0.46 | 1.09 |

  (B/particle; float32 + zstd: 3.5.) The velocity step (2.44 km/s) already gives a kinetic-energy error
  of about +-370 (km/s)^2 at 300 km/s, so dphi = 100 (max error 50) is ample: about 0.8-1.5 B/particle,
  10-20% on top of 5.3-8.9 B/particle. Python encode 5-9 ns/particle (zstd-bound) against 263 ns/particle
  for the current encoder (2048^3 z=45: 2258 s summed over threads for 8.59e9 particles).
- Estimate for 2048^3, 200 Mpc/h, 100 snapshots to z = 0 (single node, 192 threads, local SSD + drain):
  RAM +0 (already stored; writer scratch +~1 MB/thread). Encode +~70 core-s per snapshot. Stall is
  SSD-bound, so it grows with bytes: +~5-10 s per snapshot, ~15 min wall for the run, ~50 core-h.
  Disk +~0.9 TB (vs ~7 TB for positions/velocities/IDs). Drain +~6 min per snapshot in the background.
  Production 2048^3 took ~17.4 days of running on one node (file timestamps of
  /storage/swarajit/pkdgrav3_sims/pkdgrav3_n2048_sahyadri, gaps between restarts removed), ~80,000
  core-h at 192 threads: the potential adds ~0.06% (~0.1-0.15% of compute alone).

### Implemented (2026-11-10, commit 9f07c9a6)
- `dCompPotStep` (flag `compdphi`, (km/s)^2, default 0 = off). Format version 2 (docs/snapcomp_format.md):
  header `pot_step`, `pot_to_kms2` at 4080/4088, one stream per block (kind 8, zigzag delta in block
  order). Files without the potential are version 1, byte for byte as before. Readers: C++ `Block::P` /
  `potentials()`, Python `PkcSnapshot.read(fields=("pot",))`, `has_potential`, `pot_step`.
- Unit conversion `phi_code * U^2 / a` (U = code velocity unit in km/s) checked on the section-11 snapshots:
  <v^2> of the 1000 particles around the deepest well over the well depth (ring at 2-3x their radius) is
  2.00 / 2.04 / 2.54 at z = 3 / 1 / 0 with 1/a; with no factor or with a it drifts by 4x and 16x.
- Mac test (128^3, 12.5 Mpc/h = production mass resolution, 50 steps z = 49 -> 0, 12 threads, 13.3 min,
  compressed + HDF5 of the same states, dCompPotStep = 100): `test_pkcreader.py` against the HDF5
  `Potential` at z = 3, 1, 0: error <= 50.00 (km/s)^2 (half step), same order as the other fields,
  header unit = KmPerSecUnit^2/a; all other checks unchanged. Bytes from the files' stream tables:

  | z | all fields | potential | on top of the rest |
  |---|---|---|---|
  | 3 | 8.05 B/particle | 1.15 | 16.7% |
  | 1 | 8.45 | 1.17 | 16.1% |
  | 0 | 8.57 | 1.20 | 16.2% |

  (Higher than the 0.8 of the prototype estimate: the 9 subsample files thin out neighbours.)
- pawna 512^3 run with checkpoints and a restart (section 19) uses dCompPotStep = 100.

### Is pkdgrav3's potential right? (2026-11-10)
- Independent check: periodic FFT Poisson solve (CIC on 256^3, deconvolved, code units: lap phi =
  4 pi (rho - Omega_m)) on the section-11 HDF5 snapshots, interpolated to the particles, against pkdgrav3's
  `Potential` (same as the compressed one within half a step), both in (km/s)^2:

  | snapshot | slope, all | corr, all | slope, density < 2x mean | corr | offset (km/s)^2 |
  |---|---|---|---|---|---|
  | 200 Mpc/h, z = 3 | 1.001 | 0.9998 | 1.000 | 1.0000 | -2.86e5 |
  | 200 Mpc/h, z = 0 | 1.017 | 0.9952 | 1.000 | 0.9973 | -7.15e4 |
  | 12.5 Mpc/h, z = 3 | 1.001 | 0.9979 | 1.000 | 0.9908 | -1.17e3 |
  | 12.5 Mpc/h, z = 0 | 1.005 | 0.9987 | 0.994 | 0.9653 | -6.9e2 |

  Shape and normalisation agree (the 12.5 Mpc/h z = 0 grid does not resolve that box's structure).
- The constant offset is each particle's softened interaction with itself (gravity/pp.h at r = 0:
  -(105/64) m/eps = -1.64 m/eps; the fit gives -1.57 m/eps in the 200 Mpc/h box at both redshifts). It is
  the same for every particle of a snapshot, so it changes no force and no energy difference, but it
  is not zero: production 2048^3 (same m/eps in code units as 128^3 at 12.5 Mpc/h) has -293/a (km/s)^2
  (-293 at z = 0, -1170 at z = 3). Subtract it to compare with codes that leave it out. Stored as
  pkdgrav3 computes it (same as its HDF5).
- `bDoPotOutput` does not change the calculation. It only turns on the potential field (on anyway with
  gravity, `bMemPotential` defaults to 1) and, after every snapshot, does `Reorder()` and writes
  `<name>.<step>.pot` (an array file through the single writer; 13 B/particle in production). Leave it off
  when the compressed files hold the potential. The PM part of pkdgrav3 (`nGridLin`, linear species only)
  adds forces but not potential; production does not use it.

### Potential only at some redshifts (2026-11-10)
- `dCompPotMaxRedshift` (flag `compphiz`, default 1e30 = every snapshot): with `dCompPotStep` > 0, the
  potential is stored only in snapshots at z <= this. Mac 64^3 test (z = 3, 1, 0 with the cut at 1.5):
  z = 3 written as version 1 without the potential (log: "no potential (z > dCompPotMaxRedshift = 1.5)"),
  z = 1 and 0 as version 2 with it (1.07, 1.09 B/particle); pkc_check and the Python reader report it
  per snapshot.
- How much it saves depends on the output times. Snapshots evenly spaced in time from z = 49 (production:
  one every 20 of 2020 equal steps, if the steps are equal in time) put 38% of them at z <= 0.5, 58% at
  z <= 1, 77% at z <= 2: a cut at z = 1 saves ~40% of the potential's 16%, i.e. +~9% bytes instead of +16%.

## 18. Run-time estimate: 2048^3, 200 Mpc/h, z = 49 -> 0, 100 snapshots, one pawna node (2026-11-10)
From the production run's file timestamps (`ls --time-style=+%s` of
/storage/swarajit/pkdgrav3_sims/pkdgrav3_n2048_sahyadri; its logs and snap.timing are not readable):
- Production: 2020 steps, a snapshot every 20 steps, a checkpoint every 200. 17.9 days from the start of
  the first log (Oct 18 21:49) to the last snapshot, including ~0.5 d of gaps between restarts.
- Snapshot-to-snapshot interval (20 steps + one HDF5 write): median 2.98 h, 2.85 h after step 1000;
  4.8-5.4 h at high z (smaller opening angle). Intervals with a checkpoint: 6.4-8.7 h, i.e. ~4.5 h more.
- A 173 GB HDF5 snapshot at the measured NFS 25-30 MB/s is 1.6-1.9 h, so compute is ~1.0-1.25 h per 20
  steps late in the run (~3-4 min per step).
- Split of the 17.9 d: snapshot writes ~7.4 d (101 x ~1.75 h), checkpoints ~1.9 d, gaps ~0.5 d,
  compute ~8 d.
Estimate with compressed output staged on the local SSD + background drain:
- Compute ~8 d (unchanged) + snapshot stall 100 x 1-1.5 min ~2-3 h (+~0.25 h with the potential).
- Checkpoints written to NFS as today: +1.9 d -> **~10 days**. Checkpoints staged and drained like the
  snapshots (~400 GB each, minutes on the SSD; the 4-5 h drain fits in the ~10-12 h between checkpoints,
  if the SSD has room for one checkpoint + one snapshot): **~8-8.5 days**. ~37,000-46,000 core-h at 192.
- Drain: ~75-85 GB per late snapshot (8.5-9.5 B/particle with subsamples and potential) at 25 MB/s is
  50-60 min, against ~60-75 min of compute between snapshots: it keeps up, with little margin.
- Uncertain: the production HDF5 write time is inferred, not measured; the thread count of the
  production run is not known to me (192 assumed); our 2048^3 test runs with ~40 GB in swap.

## 20. Estimate: 2580^3 and 3250^3 on one node (2026-11-10)
2580^3 = 1.72e10 particles (2.00x 2048^3), 3250^3 = 3.43e10 (4.00x).
- Memory, measured at 2048^3 (job 2536432, ordered mode, no ID field, `mem.2536432.log` sampled every 15 s):
  pkdgrav3 allocates 732 GB within 30 s of start (particle store with `dExtraStore`, tree, buffers) and peaks at
  801.6 GB RSS during the first gravity step, with 25-43 GB more in swap: about 98 B/particle needed in total.
  Memory scales with N:

  | | particles | RAM needed (98 B/p) | with `bMemIntegerPosition` (~78 B/p) | RAM to buy (+15%) |
  |---|---|---|---|---|
  | 2048^3 | 8.59e9 | 0.84 TB | 0.67 TB | 1 TB |
  | 2580^3 | 1.72e10 | 1.68 TB | 1.34 TB | 2 TB (1.5 TB with integer positions, tight) |
  | 3250^3 | 3.43e10 | 3.36 TB | 2.68 TB | 4 TB (3 TB with integer positions, tight) |

  The integer-position saving (19.6 B/particle) was measured at 256^3 on the Mac (memory_persistent_ids.md)
  and checked against the compressed writer; not yet in a full-length run.
- Run time (one node, 192 threads, scaled from the ~8 d compute of 2048^3 in section 18; RAM adds no speed):
  same mass resolution (box 252 / 317 Mpc/h): cost ~ N log N, 2.06x and 4.2x: **~17 d** and **~34 d**.
  Same 200 Mpc/h box (finer mass): more particles in deep timestep rungs, guessed +20-50% more:
  **~20-25 d** and **~40-55 d**. Uncertain by at least +-25%.
- Output (100 snapshots, 8.5 B/particle late, +16% with the potential at every snapshot): ~146 / ~290 GB per
  snapshot, 14-17 / 29-34 TB per run. Drain at 25 MB/s: 1.6 / 3.2 h per snapshot, against ~4 / ~8 h between
  snapshots: keeps up.
- Checkpoints (48 B/particle): 0.82 / 1.65 TB. 3250^3 does not fit the 1.4 TB local SSD for staging: write
  checkpoints straight to home (0.27-0.49 GB/s: 1-1.7 h each) or add NVMe.
