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
   Exit: written table of ceilings and the identified bottleneck.
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
