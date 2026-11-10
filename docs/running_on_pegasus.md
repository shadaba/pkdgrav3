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

## Not done yet
- Staging with `run_staged.pbs` across several nodes (one drain per node; checkpoint pieces spread over nodes).
  It is probably unnecessary on Lustre: test writing straight to `/data` first.
- Timing at 512^3 and above, the Skylake and Ice Lake nodes, and the GPU nodes (a CUDA build would need
  `cuda-nvcc` 12.x in the environment; the P100 is Pascal).
