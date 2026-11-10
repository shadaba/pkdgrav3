# Running pkdgrav3 on pegasus (IUCAA)

Pegasus setup is kept separate from pawna's ([running_on_pawna.md](running_on_pawna.md)). Its files are
`tools/pegasus/` and `bench/pegasus/`. The feature guide is [compressed_output_guide.md](compressed_output_guide.md).

## Cluster facts (probed 2026-11-10 from login node ln4)
- **Compute nodes** (2 sockets × 16 cores each):

  | Nodes | CPU | RAM | Queue |
  |---|---|---|---|
  | hpc001-060 | Skylake | 384 GB | |
  | hpc061-080 | Ice Lake | 512 GB | |
  | hpc081-100 | Sapphire Rapids | 512 GB | debug queue uses hpc081-091 |

  Login ln3 and ln4 are Ice Lake (Xeon Gold 6326). InfiniBand is mlx5 (UCX transports rc/dc/ud_mlx5).
- **Software:** no system MPI, CMake, HDF5 or FFTW; only `gcc/11.2.0`, `boost/1.77.0` and `python3` modules,
  plus conda in `~/miniconda3`. The build therefore takes everything from one conda environment (conda-forge).
- **Storage:**
  - Home: 150 GB quota, nearly full.
  - `/data`: Lustre, 942 TB free.
  - Work directory: `/data/project/hpc2502016/shadab/pkdgrav3_io`.
- **Debug queue:** at most 2 running and 2 queued jobs per user, 128 cores and 4 h.

## Set-up (one command, on a login node)

```bash
git clone -b main https://github.com/shadaba/pkdgrav3.git /tmp/$USER-pkd
bash /tmp/$USER-pkd/tools/pegasus/setup_pegasus.sh /data/project/hpc2502016/shadab/pkdgrav3_io
rm -rf /tmp/$USER-pkd
```

`setup_pegasus.sh` does the following:
- Creates `<root>/env`, a conda environment of about 3 GB. Its package cache goes in `<root>/conda_pkgs`,
  not in home. It contains:
  - Python 3.11, numpy < 2, Cython, dill, tomli, h5py, zstandard;
  - the conda C/C++ compilers, make, cmake, git;
  - Open MPI 5 (with UCX), FFTW with MPI and threads, HDF5 (serial), GSL, Boost, libhwloc, zstd, zlib.
- Clones the repository into `<root>/pkdgrav3`.
- Runs `tools/pegasus/compile_pegasus.sh`, which builds `<root>/build/pkdgrav3` and
  `<root>/build_snapcomp/{pkc_check,pkc_to_gadget,snapcomp_test}`.

The build is for `-march=skylake-avx512`, so one binary should run on every node type; it has run only on the
Sapphire Rapids debug nodes so far. `ARCH=auto` builds for the login node's CPU instead (Ice Lake).

The first build took about 1 h 15 min on Lustre with `-j 8`. To update:

```bash
cd $ROOT/pkdgrav3 && git pull --ff-only && PKD_IO_ROOT=$ROOT bash -l tools/pegasus/compile_pegasus.sh
```

Pitfalls met (all handled by the scripts):
- conda-forge calls hwloc `libhwloc`.
- libdwarf (fetched by pkdgrav3) needs `zlib.h`, and the conda compilers do not see the system headers.
- Home is nearly full, hence `CONDA_PKGS_DIRS` in the work directory.

## Job scripts
Start a job script with:

```bash
export PKD_IO_ROOT=/data/project/hpc2502016/shadab/pkdgrav3_io
source $PKD_IO_ROOT/pkdgrav3/tools/pegasus/env_pegasus.sh
```

Then launch:
- **One node:** `mpirun -np 1 --bind-to none $PKD_IO_ROOT/build/pkdgrav3 -sz 32 run.par`, or
  `-np 2 --map-by ppr:1:socket ... -sz 16` for one process per socket.
- **Several nodes:**

  ```bash
  mpirun -np <2*nodes> --hostfile $PBS_NODEFILE --map-by ppr:2:node --bind-to none --mca pml ucx \
         $mpirun_pegasus_args $PKD_IO_ROOT/build/pkdgrav3 -sz 16 run.par
  ```

  Ask PBS for `select=<nodes>:ncpus=32:mpiprocs=2` and `place=scatter`. `$mpirun_pegasus_args` (set by
  `env_pegasus.sh`) passes `--prefix` and the environment to the remote nodes, which are started over ssh.

The UCX library prints `cuDeviceGetCount ... initialization error` at start-up on CPU nodes (it probes for
CUDA). This is harmless.

## Tests (2026-11-10)

| Job | Script | Nodes | Run | Result |
|---|---|---|---|---|
| 37408 | `bench/pegasus/smoke.pbs` | 1 (hpc081, Sapphire Rapids) | 64^3, L = 6.25 Mpc/h, z = 49 -> 12 in 6 steps, compressed + potential + HDF5 | 1 process × 32 threads: 18 s, 6.08/6.69/7.17 B/particle at z = 45/25/12 (potential 1.03), compressed stall 0.02 s. 2 processes × 16 threads: 12 s, 18 files (2 writers). `pkc_check` OK on every snapshot (all IDs once, right subsamples, same 1% set). `test_pkcreader` against HDF5: ALL OK. |
| 37409 | `bench/pegasus/mpi2.pbs` | 2 (hpc081 + hpc082) | 128^3, L = 12.5 Mpc/h, same steps, 4 processes × 16 threads, `--mca pml ucx` | 78.5 s wall, 3.1 GB peak RSS per process. 6.06/6.48/7.02 B/particle (potential 1.16). 36 files per snapshot (4 writers). `pkc_check` OK on all three snapshots. `test_pkcreader`: ALL OK. |

These builds include the Ewald potential fix (commit 2219e678; compressed_io_plan.md section 19).

## Write to /data directly: local staging is not needed (job 37410, `bench/pegasus/io512.pbs`)
Setup: 4 debug nodes (hpc081-084), 8 MPI processes × 16 threads. 512^3, L = 50 Mpc/h (production mass resolution),
z = 49 -> 12 in 6 steps, compressed snapshots with the potential at z = 45, 25 and 12.

**Disks.** The compute nodes have a node-local ext4 `/local` with 514 GB free; `/tmp` and `/var/tmp` have about
93 GB. `dd` measurements (8 GiB, O_DIRECT):

| Target | Writers | Speed |
|---|---|---|
| `/local` | 1 | 1.3 GB/s |
| `/data` (Lustre) | 1 | 1.0 GB/s |
| `/data` | 32 (8 per node, 2 GiB each) | 64 GiB in 1.7 s, about 37 GB/s |

The 37 GB/s is an upper bound: the Lustre servers may cache, and the data was zeros.

**Snapshot stall** (time the simulation waits), the same three snapshots in each run:

| Run | Where | Stall per snapshot | Run wall time |
|---|---|---|---|
| A | `/data`, write-through auto (off: plenty of memory) | 0.69 / 0.59 / 0.62 s | 33 min 21 s (includes HDF5) |
| A, HDF5 of the same states (serial writer) | `/data` | 306 / 283 / 258 s | |
| C | `/data`, write-through forced (every block on disk before the next) | 0.61 / 0.60 / 0.63 s | 18 min 53 s |
| B | `/local` | 1.19 / 0.80 / 0.84 s | 18 min 30 s |
| B, then copying the files to `/data` | all 4 nodes at once | 2.7 GiB in 2.9 s (944 MiB/s, 1 copy per node) | |

Snapshots are 880 / 941 / 1016 MB (6.55 / 7.01 / 7.57 B/particle, of which the potential is 1.26-1.30).
`pkc_check` passes on all three runs: all 134,217,728 IDs once, 72 files from 8 writers.

Conclusion: on pegasus, write snapshots and checkpoints straight to `/data`. Writing through to Lustre is as fast
as the node's own disk, so staging plus drain adds nothing; `run_staged.pbs` is for pawna's slow NFS. HDF5
remains about 280 s per 512^3 snapshot because it is the serial writer, not because of the disk. At 2048^3
(about 55-80 GB per snapshot over 4+ nodes) the stall is expected to be seconds (not yet measured: jobs
37412/37413 below ended before their snapshot).

## Large N on 4 debug nodes: memory and step cost (jobs 37412, 37413, `bench/pegasus/big.pbs`)
Setup: 4 Sapphire Rapids debug nodes (hpc081-084, 128 cores), 8 MPI processes x 16 threads (120 worker
threads), `place=scatter:excl`. Production accuracy and stepping (sahyadri2048.par: classic theta and replica
switches, Ewald, dEta 0.3, bEpsAccStep, 20 equal-time steps per output from z = 49 to 11.98), 2LPT ICs
generated at z = 49, ordered-mode IDs, softening 8.138e-6 x 2048/N. Memory: `free -m` on every node every
15 s, and pkdgrav3's own "resident size" / "free memory" report. Times: pkdgrav3's timers in the run log
and `snap.log` (wall seconds per step).

| | 2048^3, 200 Mpc/h (37412) | 2500^3, 400 Mpc/h (37413) |
|---|---|---|
| Particles | 8.59e9 | 1.56e10 |
| Peak used per node (rank-0 node / others) | 220 / 196-199 GB | 387 / 364-366 GB |
| Lowest MemAvailable per node | 289 GB | 125 GB |
| pkdgrav3 resident size per process | 94.0 GB | 175.8 GB |
| Total, per particle | ~94 B | ~95 B |
| Swap | < 0.3 GB (unchanged from start) | < 0.3 GB |
| IC generation (2LPT) | 81 s | 148 s |
| First domain decomposition / tree | 113 s / 16 s | 335 s / 32 s |
| Gravity, all particles active (rung 0-2) | 1011-1113 s (1870-2050 Gflop/s) | 1529-1536 s (2380-2410 Gflop/s) |
| Gravity, rung 3 (42-53% active) | 862-965 s | 1103 s |
| Domain + tree per substep | ~31 s | ~57 s |
| One production-size step at z ~ 45 | **8416 s (2.34 h)**, measured (step 1) | ~11,000 s (3.05 h), from the substeps |
| Ended | walltime at step 1.6 of 4 (z = 36.6), no snapshot | cancelled after 3 gravity calls (data enough) |

Both fit in 4 nodes; memory is ~95 B/particle on pegasus (pawna measured 98 with swap in use). At
z ~ 45 nearly every particle is on rung 2 or 3, so each production step costs 8 gravity calls; this is set by
the timestep criterion, not by our test. Per particle, a step of 2500^3/400 costs 0.72x that of 2048^3/200
(coarser mass resolution).

Compared with pawna (same accuracy, 2048^3, job 2536432, 192 threads): first gravity 603 s there vs
1111 s here, i.e. a Sapphire Rapids core does ~0.82x the work of a pawna core. pawna ran z = 49 -> 12
in 12 h 50 min (production's first 20 steps took ~13 h of compute), so the 2.34 h measured here for one of
those steps on 128 cores is consistent (the first steps are the costliest: theta 0.4 above z = 20).

Consequence: the debug queue (4 h, 128 cores) cannot run a production model; 2048^3 on 4 nodes would need
~15 days (~90 restarts of 4 h). A production queue is needed. The full-run estimate for 3200^3 is in
compressed_io_plan.md section 22.

## Not done yet
- Strong scaling beyond 4 nodes (the debug queue allows 128 cores): needed before choosing between 8, 16 and
  32 nodes for 3200^3.
- The compressed-write stall at 2048^3 and above on pegasus.
- Timing on the Skylake and Ice Lake nodes, and the GPU nodes (a CUDA build would need `cuda-nvcc` 12.x in the
  environment; the P100 is Pascal).
