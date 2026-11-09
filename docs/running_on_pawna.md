# Running pkdgrav3 on pawna: build, I/O options and snapshot staging

How to build our pkdgrav3 (3.5.1 + branch `feature/compressed-output`) on pawna and run it with the
different snapshot I/O options, including staging snapshots on the node-local SSD and draining them
to NFS during the run. Everything here was tested; measurements are in
[compressed_io_plan.md](compressed_io_plan.md) sections 8-9.

## 1. Where things are

| What | Path |
|---|---|
| Source (git clone of the fork) | `/storage/shadab/pkdgrav3_io/pkdgrav3` |
| Build | `/storage/shadab/pkdgrav3_io/build/pkdgrav3` |
| Python venv (Cython 3, dill, tomli; system NumPy) | `/storage/shadab/pkdgrav3_io/venv` |
| Runtime environment (source it in every job) | `bench/phase0/env_pawna.sh` |
| Staged-run job script | `tools/staging/run_staged.pbs` |
| Drain tool | `tools/staging/pkd_drain.py` |
| HDF5 gzip-level shim | `tools/h5shim/` |
| Benchmark jobs and logs | `bench/phase0/`, `/storage/shadab/pkdgrav3_io/bench/logs/` |

Storage seen from a compute node:

| Path | What it is | Free | Write speed (Phase 0, 4 MB blocks, flushed) |
|---|---|---|---|
| `/storage` | NFS 4.2 from 10.20.255.85 | 216 TB | 28 / 25 / 28 / 33 MB/s with 1 / 4 / 16 / 64 streams (total is capped) |
| `/user` (home) | NFS 4.2 from the login node | 36 TB | 270 / 325 / 395 / 488 MB/s (8 GB tests; upper bound) |
| `/var/tmp`, `/tmp` | node-local SSD array (ext4, also `/`) | 1.4 TB per node | 427 / 797 / 1458 / 1035 MB/s with 1 / 4 / 16 / 64 streams |
| `/dev/shm` | RAM | 378 GB | memory speed; needed by the simulation at 2048^3 |

The links are 10 GbE, so the NFS ceiling comes from the server side or other users' load, not the
network. Questions for the admins: is `/storage` exported `sync`; what disks back it; is `nconnect`
possible.

## 2. Setting up from scratch (from GitHub)

Everything lives in one root directory of your choice (called `$ROOT` below; for example
`/storage/$USER/pkdgrav3_io`). Nothing from an existing build is needed.

### 2a. One command

```bash
ssh pawna
git clone -b feature/compressed-output https://github.com/shadaba/pkdgrav3.git /tmp/$USER-pkd-setup
bash /tmp/$USER-pkd-setup/tools/pawna/setup_pawna.sh /storage/$USER/pkdgrav3_io      # 20-40 min
rm -rf /tmp/$USER-pkd-setup
qsub -v PKD_IO_ROOT=/storage/$USER/pkdgrav3_io /storage/$USER/pkdgrav3_io/pkdgrav3/bench/phase0/smoke32.pbs
cat /storage/$USER/pkdgrav3_io/bench/logs/p0_smoke32.<job>.log                       # expect rc=0
```

### 2b. The same, step by step

```bash
ROOT=/storage/$USER/pkdgrav3_io
mkdir -p $ROOT/bench/logs && cd $ROOT
# 1. source (the fork, branch with compressed output; never push to upstream)
git clone -b feature/compressed-output https://github.com/shadaba/pkdgrav3.git pkdgrav3
(cd pkdgrav3 && git submodule update --init)          # EOSlib, ANEOSmaterial
# 2. Python: system 3.10 with the system NumPy 1.21 (its headers are the ones the compiler finds)
python3 -m venv --system-site-packages venv
venv/bin/python -m pip install "Cython>=3.0" "dill>=0.3.6" tomli "h5py<3.12" zstandard
# 3. build pkdgrav3 -> $ROOT/build/pkdgrav3, snapcomp tools -> $ROOT/build_snapcomp/, gzip shim
PKD_IO_ROOT=$ROOT bash -l pkdgrav3/bench/phase0/compile_pawna.sh
```

What `compile_pawna.sh` does (see the script): `source bench/phase0/env_pawna.sh` (module
`openmpi-5.0.1`, unset the Intel MPI paths, `PYTHONPATH` to the venv), then
`cmake -S pkdgrav3 -B build -DCMAKE_BUILD_TYPE=Release -DPython3_EXECUTABLE=$ROOT/venv/bin/python`,
`cmake --build build -j 16`, the same for `pkdgrav3/io/snapcomp` into `build_snapcomp`, and `make` in
`tools/h5shim`. pkdgrav3 downloads cpptrace, libdwarf and zstd during the first cmake; snapcomp downloads
zstd for its tools (pawna has no zstd headers). Dependencies used from the system: GCC 11, CMake 3.22,
HDF5 1.10.7, GSL 2.7, FFTW 3, Boost 1.74, hwloc.

Updating later: `cd $ROOT/pkdgrav3 && git pull --ff-only && PKD_IO_ROOT=$ROOT bash -l bench/phase0/compile_pawna.sh`
(`CLEAN=1` in front for a full rebuild of pkdgrav3).

Job scripts in this repository take the root from `PKD_IO_ROOT` (default `/storage/shadab/pkdgrav3_io`):
`qsub -v PKD_IO_ROOT=$ROOT,... script.pbs`. In your own job scripts: `export PKD_IO_ROOT=$ROOT;
source $ROOT/pkdgrav3/bench/phase0/env_pawna.sh`, then
`mpirun --bind-to none -np 1 -x PYTHONPATH $ROOT/build/pkdgrav3 -sz $NCPUS run.par`.

The login node is for builds and quick probes only; runs go through PBS.

Pitfalls that cost time (all handled by the scripts):
- The login profile puts Intel oneAPI MPI into `CPATH` and `LIBRARY_PATH`. gcc searches `CPATH`
  before CMake's `-isystem`, so Intel's `mpi.h` silently replaces Open MPI's and the binary crashes
  at start-up in `MPI_Comm_size` (fault address `0x440000e0` = MPICH's `MPI_COMM_WORLD`).
  `compile_pawna.sh` and `env_pawna.sh` unset these.
- `/usr/include/python3.10/numpy` points at the system NumPy 1.21 headers and is searched first, so
  the venv uses the system NumPy too (a venv NumPy 2.x gives Cython compile errors).
- pkdgrav3 embeds the system Python; `env_pawna.sh` sets `PYTHONPATH` to the venv so it finds `dill`.
- `.gitignore` ignores `build*`, so scripts must not be named `build_*`.
- PBS `-o` copy-back from compute nodes loses output here; every job script writes its own log with
  `exec > logfile`.

### 2c. Persistent IDs without extra memory (use this for 2048^3)

`bMemParticleID = True` grows every particle record from 48 to 56 bytes (+76 to +128 GB at 2048^3, which no
longer fits a 755 GB node). Use ordered mode instead; the particle header then carries `iOrder`, a persistent
ID, at no extra cost:

```python
bMemUnordered  = False   # particles keep iOrder (42-43 bits) in their 8-byte header
bMemParticleID = False   # no separate 8-byte ID field
```

- `iOrder` is set from the initial conditions and moves with the particle; for pkdgrav3-generated ICs it is
  the Lagrangian lattice index `(iz*N + iy)*N + ix`, exactly what `ParticleID` would hold.
- Compressed output (`bCompressedOutput`) writes `iOrder` as the ID automatically and skips any re-sorting.
- Standard HDF5/Tipsy output in ordered mode has no ID column: pkdgrav3 sorts particles by `iOrder`
  (`Reorder()`) before writing, so the row number is the ID.
- ICs read from a file: `iOrder` is the particle's row in the IC file. If the file's IDs are its row numbers
  (0- or 1-based) this is the same thing; otherwise keep the IC file's ID column as the lookup table.
- Measurements and the record layouts: [memory_persistent_ids.md](memory_persistent_ids.md).

## 3. I/O parameters that matter

| Parameter | Default | Meaning |
|---|---|---|
| `bHDF5` (`+hdf5` / `-hdf5`) | false | HDF5 (shuffle + gzip 9 + Fletcher32, 65536x1 chunks) or Tipsy (raw) |
| `achOutPath` | `{name}.{step:05d}{type}` | logs, and the default for everything else |
| `achIoPath` | = `achOutPath` | **snapshots only**: point this at local disk to stage |
| `achCheckpointPath` | = `achOutPath` | checkpoints (keep on NFS) |
| `bParaWrite` / `nParaWrite` | false / 0 | parallel writing; does not work on pawna (open item) |
| `bDoPotOutput` | 0 | potential is only computed on request, but Density and Potential columns are written anyway |
| `iWallRunTime` | 0 | **minutes** of wall time before pkdgrav3 checkpoints and stops |
| `iCheckInterval` | 0 | steps between checkpoints |

The gzip level is hardcoded to 9 in 3.5.1. To change it without recompiling, preload
`tools/h5shim/libh5_deflate_shim.so` and set `PKD_H5_DEFLATE` (0 = no gzip, 1 = what production
uses). Build it once:

```bash
cd /storage/shadab/pkdgrav3_io/pkdgrav3/tools/h5shim && make
```

## 4. Running

All runs go through `tools/staging/run_staged.pbs`, configured with `qsub -v`. It writes in `RUNDIR`:
`staged.<job>.log` (summary), `run.<job>.log` (pkdgrav3 output), `time.<job>.log` (peak memory),
`drain.<job>.log` (copy events).

Set the size with `-l select=1:ncpus=N:mem=M` and the wall time with `-l walltime=` (defaults 64
cores, 96 GB, 2 h). Keep `iWallRunTime` (minutes) below the PBS wall time by at least one snapshot's
drain time so the last snapshot reaches NFS.

Example run directory (fresh start):

```bash
mkdir -p /storage/shadab/runs/test512 && cd /storage/shadab/runs/test512
cp /storage/shadab/pkdgrav3_io/pkdgrav3/bench/phase0/baseline512.par test512.par
cp /storage/shadab/pkdgrav3_io/pkdgrav3/examples/euclid_z0_transfer_combined.dat .
S=/storage/shadab/pkdgrav3_io/pkdgrav3/tools/staging/run_staged.pbs
```

### 4a. Staged on local SSD, drained to NFS (recommended)

```bash
qsub -v RUNDIR=$PWD,PAR=test512.par $S
```
Snapshots go to `/var/tmp/shadab_stage_test512/current/` and end up in `$PWD/snapshots/`. The
simulation waits only for the local write; the copy runs at low priority (`nice 19`, `ionice idle`).

### 4b. Straight to NFS (what production does)

```bash
qsub -v RUNDIR=$PWD,PAR=test512.par,STAGING=0 $S
```

### 4c. HDF5 compression level

```bash
qsub -v RUNDIR=$PWD,PAR=test512.par,H5LEVEL=1 $S   # production setting
qsub -v RUNDIR=$PWD,PAR=test512.par,H5LEVEL=0 $S   # HDF5 without gzip
```
For Tipsy instead of HDF5 set `bHDF5 = False` in the parameter file.

### 4d. Restart from a checkpoint

```bash
qsub -v RUNDIR=$PWD,RESTART=snap.00002.chk $S
```
The script writes
```python
import PKDGRAV
PKDGRAV.restore("snap.00002.chk", achIoPath="/var/tmp/shadab_stage_test512/current/{name}.{step:05d}{type}")
```
Two things learned the hard way:
- Use `import PKDGRAV` and `PKDGRAV.restore(...)`, **not** `from PKDGRAV import restore`. pkdgrav3
  only treats the script as being in control if the `PKDGRAV` module object is in its namespace;
  otherwise, after `restore()` returns, it generates the ICs and runs the whole simulation again,
  overwriting snapshots (tested on the Mac: steps 3-4 after the restart, then IC + steps 1-4 again).
  This is the "second pass" the production guard scripts kill.
- A checkpoint stores `achIoPath`. Without the override, a restart writes snapshots to the old
  path. The staging directory is therefore named after the run, not the job, so it stays the same.

### 4e. Running the production binary with staging

```bash
cat > /storage/shadab/runs/env_production.sh <<'EOF'
module purge; module load openmpi-5.0.1
source /user/swarajit/anaconda3/etc/profile.d/conda.sh; conda activate pkdgrav_no_mpi
EOF
qsub -v RUNDIR=$PWD,PAR=test512.par,BIN=/user/swarajit/apps_user/pkdgrav3_openmpi/bin/pkdgrav3,ENVSH=/storage/shadab/runs/env_production.sh,H5LEVEL=1 $S
```
(Read-only use of the production install; nothing there is modified.)

### 4f. Compressed parallel output (bCompressedOutput)

Add to the parameter file (all but the first line are optional; values shown are the defaults):

```python
bCompressedOutput  = True
dCompPosStep       = 1.95            # kpc/h, max position error is half of it
dCompVelStep       = 2.44            # km/s,  max velocity error is half of it
lstCompSubsamples  = [1, 3, 5, 9, 10, 12, 15, 20, 25]   # percent, sum 100; [100] = one file
iCompSubsampleSeed = 0
iCompBlockSize     = 65536           # particles per block (scratch about 150 B/particle per thread)
iCompLevel         = 3               # zstd level
```

Every thread compresses its own particles and writes in parallel; nothing is gathered. Each snapshot
becomes `<name>.s<percent>.<process>.pkc` (one file per subsample and MPI process) and a
`<name>.done` marker per process once its files are complete. Checkpoints are unchanged.
It works with staging (4a) as is: the drain copies `<name>.*` when it sees the marker.
Dark matter only (other species are skipped with a warning). `bCompAlsoStandard = True` also writes
the normal HDF5/Tipsy snapshot from the same state, for validation; do not combine it with staging.

Reading (Python, numpy + zstandard): see `tools/snapcomp_proto/pkc_reader.py`:

```python
import sys; sys.path.insert(0, "/storage/shadab/pkdgrav3_io/pkdgrav3/tools/snapcomp_proto")
from pkc_reader import read_snapshot, PKCFile
pos, vel, ids = read_snapshot("snapshots/snap.00100", percents=[1])   # Mpc/h, km/s, uint64
```

### 4g. Options summary

| Variable | Default | Meaning |
|---|---|---|
| `RUNDIR` | required | run directory on NFS (logs, checkpoints) |
| `PAR` / `RESTART` | one required | fresh parameter file, or checkpoint to restart from |
| `DEST` | `$RUNDIR/snapshots` | final snapshot directory |
| `STAGING` | 1 | 1 = local SSD + drain, 0 = straight to `DEST` |
| `STAGE` | `/var/tmp/$USER_stage_<rundir name>` | local staging directory |
| `BW` | 0 | drain rate limit, MB/s |
| `MINFREE` | 200 | GB; below this, new snapshots go straight to NFS |
| `H5LEVEL` | unset | HDF5 gzip level via the shim |
| `BIN`, `ENVSH`, `NTH` | our build, `env_pawna.sh`, PBS ncpus | executable, environment, threads |

## 5. The drain tool

`pkd_drain.py` (Python standard library only) runs next to pkdgrav3 in the same job:
- **finished snapshot**: the log shows `Writing <path> ...` followed by
  `Output file has been successfully written`, or a `<path>.done` file exists (for the future
  compressed writer);
- **copy**: every piece (`<name>`, `<name>.*`) to `DEST/<file>.part` in 16 MiB chunks, fsync, size
  check, rename, fsync the directory, then delete the local copy;
- **disk guard**: `<stage>/current` is a symlink to `<stage>/local`; below `--min-free-gb` it is
  switched to `DEST` so new snapshots go straight to NFS, and switched back when space recovers;
- **end**: exits when pkdgrav3 has exited and everything is copied (exit 0), or exits 2 and lists
  what is left. `run_staged.pbs` deletes the staging directory only on exit 0.

Recovering snapshots left on a node (job killed, drain error): submit a job **on the same node**:

```bash
qsub -l select=1:host=node4:ncpus=1 -l walltime=04:00:00 -- /usr/bin/python3 \
  /storage/shadab/pkdgrav3_io/pkdgrav3/tools/staging/pkd_drain.py --recover \
  --stage /var/tmp/shadab_stage_test512 --dest /storage/shadab/runs/test512/snapshots \
  --report /storage/shadab/runs/test512/recover.log
```
A snapshot that was interrupted mid-write is copied too; check it before use.

Watch a running job:

```bash
qstat -u shadab
tail -f /storage/shadab/runs/test512/drain.<job>.log
grep "successfully written" /storage/shadab/runs/test512/run.<job>.log   # stall per snapshot
```

## 6. What staging buys at the current NFS speed

Measured at 512^3 (HDF5 gzip 9, 2.4 GiB per snapshot): 370 s stall writing straight to NFS, 280 s
staged (the same as writing to RAM); each snapshot was drained in about 95 s during the run.
The remaining 280 s is the single-threaded HDF5/gzip writer, which Phases 1-3 replace.


With NFS at about 25-30 MB/s in total, a 166 GB production snapshot needs about 1.5-1.8 h to copy.
The simulation no longer waits for it, but the drain must finish before the next snapshot fills the
local disk (1.4 TB holds about 8 such snapshots; the disk guard sends the overflow straight to NFS).
A 3-5x smaller compressed snapshot (Phases 1-3) cuts the copy to about 20-35 min.
