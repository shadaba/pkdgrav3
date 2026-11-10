# Compressed parallel output for pkdgrav3: what was added and how to use it

Branch `feature/compressed-output` of https://github.com/shadaba/pkdgrav3 (pkdgrav3 3.5.1 plus the changes
below). This page is the user guide. For details see:
- set-up and staging on pawna: [running_on_pawna.md](running_on_pawna.md);
- file format: [snapcomp_format.md](snapcomp_format.md);
- design, measurements and their methods: [compressed_io_plan.md](compressed_io_plan.md) (sections in brackets below).

## 1. Features

| Feature | What it does | Where |
|---|---|---|
| Compressed parallel snapshots | `bCompressedOutput = True`. Every thread compresses its own particles and appends blocks to shared files at the same time. Nothing is gathered and nothing is sorted. Positions and velocities are quantised to a step you choose (the maximum error is half the step). IDs are lossless. Data is zstd compressed. 5.3-8.9 B/particle, against 32-40 for HDF5. Checkpoints are unchanged. | `io/snapcomp/` (codec, no pkdgrav3 headers), `core/writecompressed.*`, `msr_writecompressed.cxx` [13, 14] |
| Subsamples | `lstCompSubsamples`, for example `[1,3,5,9,10,12,15,20,25]` %. There is one file per subsample. A particle's subsample depends only on its ID (a hash with seed `iCompSubsampleSeed`), so it is the same in every snapshot. Reading the 1% file reads 1% of the bytes. | [12] |
| Persistent IDs at no memory cost | Ordered mode (`bMemUnordered = False`, `bMemParticleID = False`). The compressed writer stores `iOrder` as the ID. For pkdgrav3 ICs this equals the Lagrangian lattice index, the same as `ParticleID`. It saves the +128 GB that `bMemParticleID` costs at 2048^3. | [memory_persistent_ids.md](memory_persistent_ids.md), [15] |
| Write-through when memory is short | `iCompWriteThrough = -1` (auto). Blocks bypass the page cache only when the node's free memory is less than twice the snapshot size. This avoids swapping out the simulation at 2048^3. A 512^3 run keeps its 1.5 s stall. | [15] |
| Gravitational potential (optional) | `dCompPotStep > 0` stores the peculiar potential in (km/s)^2 with that step (format version 2). It costs about 1.2 B/particle (+16%) at a step of 100. `dCompPotMaxRedshift` stores it only at low redshift. Checked against pkdgrav3's HDF5 `Potential` and an independent FFT Poisson solve. | [17] |
| Snapshot staging and drain | `tools/staging/run_staged.pbs` with `pkd_drain.py`. Snapshots are written to the node's local SSD (about 1 GB/s) and copied to NFS in the background. The disk guard switches back to writing straight to NFS when the SSD is full. | [9], running_on_pawna.md section 4-5 |
| Checkpoint staging | `CHKSTAGE=1`. Checkpoints are also written to the local SSD and copied by a second worker to their own destination, which may be another file server (for example home). The restart file is copied last. Only the newest `KEEPCHK` complete checkpoints are kept. | [19], `tools/staging/test_drain.py` |
| Readers and tools | `tools/pkcreader/pkcreader.py` (Python: numpy + zstandard). `pkc_check` (IDs, subsamples, persistence between snapshots). `pkc_to_gadget` (Gadget-2 files for Rockstar). `pkc_convert` (HDF5 to .pkc). `bench/phase5/cmp_runs.py` (compares two runs particle by particle). | [16] |
| pawna set-up | `tools/pawna/setup_pawna.sh`: clone, venv and build in a new directory with one command. | running_on_pawna.md section 2 |

Not covered: other species (gas and stars are skipped with a warning); more than one MPI process per node
was tested only on the Mac; Rockstar reading `.pkc` natively (planned, see
[rockstar_pkc_session_prompt.md](rockstar_pkc_session_prompt.md)).

> **Potential fix (2026-11-10):** pkdgrav3 3.5.1 adds part of the Ewald potential with the wrong sign (forces are
> fine). Fixed in commit 2219e678 of this branch; potentials from older builds, HDF5 included, are wrong for
> some particles. See compressed_io_plan.md section 19.

## 2. Measured results (how: see the sections cited)

| Run | Bytes per particle | Stall per snapshot | Notes |
|---|---|---|---|
| 2048^3, 200 Mpc/h, z = 45 / 25 / 12, one pawna node, 192 threads, local SSD + drain [15] | 5.31 / 5.73 / 6.30 | 50 / 54 / 58 s | 45.6-54.1 GB per snapshot; peak RSS 747 GiB; all 8.59e9 IDs once; the 1% set is the same in every snapshot; drain to NFS 21-23 MB/s |
| 512^3, 200 Mpc/h, z = 3 to 0, 64 threads, straight to NFS [15, 16] | 8.25-8.89 | 1.5 s (auto write-through off) | HDF5 of the same state: 384-443 s |
| 128^3, 12.5 Mpc/h, Mac, potential step 100 [17] | 8.05-8.57 in total, of which 1.15-1.20 is the potential | - | potential error at most 50 (km/s)^2 = half a step |

Defaults `dCompPosStep = 1.95` kpc/h and `dCompVelStep = 2.44` km/s give maximum errors of 0.98 kpc/h and
1.22 km/s.

## 3. How to run

### 3.1 Parameter-file settings

Add to a normal pkdgrav3 parameter file. Only the first line is required; the other values shown are the defaults:

```python
bCompressedOutput   = True
dCompPosStep        = 1.95       # kpc/h; max error 0.98
dCompVelStep        = 2.44       # km/s;  max error 1.22
lstCompSubsamples   = [1, 3, 5, 9, 10, 12, 15, 20, 25]   # percent, sum 100; [100] = one file per process
iCompSubsampleSeed  = 0
iCompBlockSize      = 65536      # particles per block
iCompLevel          = 3          # zstd level
iCompWriteThrough   = -1         # -1 auto, 0 off, 1 on
dCompPotStep        = 0.0        # (km/s)^2; 0 = no potential
dCompPotMaxRedshift = 1.0e30     # potential only in snapshots at z <= this
bCompAlsoStandard   = False      # also write the HDF5/Tipsy snapshot (validation only)
```

For persistent IDs without extra memory (recommended, essential at 2048^3):

```python
bMemUnordered  = False
bMemParticleID = False
```

Leave `bDoPotOutput` off. It does not change the potential, but it adds a slow single-writer `.pot` file
after every snapshot [17].

Each snapshot `<name>.<step>` becomes `<name>.<step>.s<percent>.<process>.pkc` (one file per subsample and
MPI process) and `<name>.<step>.done` (one per process, written when its files are complete).

### 3.2 Variants

Examples. The comments say what changes.

```python
# a) Defaults: subsamples, no potential
bCompressedOutput = True

# b) One file per process, coarser positions (max error 2 kpc/h): fewer bytes
bCompressedOutput = True
lstCompSubsamples = [100]
dCompPosStep      = 4.0

# c) Potential in every snapshot, max error 50 (km/s)^2 (about +16% bytes)
bCompressedOutput = True
dCompPotStep      = 100.0

# d) Potential only from z = 1 on (about +9% bytes over the run if snapshots are evenly spaced in time)
bCompressedOutput   = True
dCompPotStep        = 100.0
dCompPotMaxRedshift = 1.0

# e) Validation: compressed and HDF5 from the same state (do not combine with staging)
bCompressedOutput = True
bCompAlsoStandard = True
bHDF5             = True
```

The run log shows what each snapshot got, for example:
`Writing snap.00040 compressed ..., potential dphi 100 (km/s)^2` or
`... no potential (z > dCompPotMaxRedshift = 1)`, then
`Compressed N particles into X MB (Y B/particle) ...` and `of which potential: ...`.

### 3.3 Laptop or workstation (no PBS)

```bash
cd /path/to/run
cp $SRC/examples/euclid_z0_transfer_combined.dat .
# test.par: a 64^3-256^3 parameter file with the lines of 3.1
mpirun -np 1 /path/to/build/pkdgrav3 -sz 8 test.par        # or: /path/to/build/pkdgrav3 test.par
python3 -c 'import sys; sys.path.insert(0, "'$SRC'/tools/pkcreader"); from pkcreader import PkcSnapshot; \
  s = PkcSnapshot("snap.00010"); print(s.redshift, s.npart, s.has_potential)'
```

### 3.4 pawna: one-time set-up

```bash
git clone -b feature/compressed-output https://github.com/shadaba/pkdgrav3.git /tmp/$USER-pkd-setup
bash /tmp/$USER-pkd-setup/tools/pawna/setup_pawna.sh /storage/$USER/pkdgrav3_io      # 20-40 min
rm -rf /tmp/$USER-pkd-setup
```
After this, `$ROOT=/storage/$USER/pkdgrav3_io` holds `pkdgrav3/` (source), `build/pkdgrav3`, `build_snapcomp/`
(`pkc_check`, `pkc_to_gadget`, `pkc_convert`) and `venv/`. To update: `cd $ROOT/pkdgrav3 && git pull --ff-only && PKD_IO_ROOT=$ROOT bash -l bench/phase0/compile_pawna.sh`.

### 3.5 pawna: runs with `run_staged.pbs`

All settings are passed with `qsub -v`. Set the size with `-l`.

```bash
ROOT=/storage/$USER/pkdgrav3_io; S=$ROOT/pkdgrav3/tools/staging/run_staged.pbs
mkdir -p /storage/$USER/runs/test512 && cd /storage/$USER/runs/test512
cp $ROOT/pkdgrav3/examples/euclid_z0_transfer_combined.dat .     # and write test512.par
```

| Goal | Command |
|---|---|
| Snapshots staged on local SSD and drained to `./snapshots` (recommended) | `qsub -v PKD_IO_ROOT=$ROOT,RUNDIR=$PWD,PAR=test512.par $S` |
| Snapshots straight to NFS | `qsub -v ...,STAGING=0 $S` |
| Snapshots elsewhere | `qsub -v ...,DEST=/storage/$USER/snaps/test512 $S` |
| Checkpoints staged too, copied to home, newest 2 kept | `qsub -v ...,CHKSTAGE=1,CHKDEST=/user/$USER/chk/test512,KEEPCHK=2 $S` |
| Keep all checkpoints | `...,KEEPCHK=0` |
| Restart (snapshots staged) | `qsub -v PKD_IO_ROOT=$ROOT,RUNDIR=$PWD,RESTART=snap.00100.chk $S` |
| Restart from a staged checkpoint | `qsub -v PKD_IO_ROOT=$ROOT,RUNDIR=$PWD,CHKSTAGE=1,CHKDEST=/user/$USER/chk/test512,RESTART=/user/$USER/chk/test512/snap.00100.chk $S` |
| Limit drain bandwidth / change the disk guard | `...,BW=200,MINFREE=300` (MB/s, GB) |
| Bigger node and longer wall time | `qsub -l select=1:ncpus=192:mem=750gb -l walltime=240:00:00 -v ... $S` |

`RUNDIR` holds the logs: `staged.<job>.log` (summary: stalls, drain events), `run.<job>.log` (pkdgrav3),
`time.<job>.log` (peak memory) and `drain.<job>.log`. The staging directory `/var/tmp/$USER_stage_<rundir name>`
is named after the run, so it is the same at every restart. It is deleted only when everything has been
copied.

For a restart, pass the same `CHKSTAGE`/`CHKDEST` as the first job. The script overrides the paths stored
in the checkpoint (`achIoPath`, `achCheckpointPath`), so new snapshots and checkpoints are staged again.

### 3.6 Stopping a run cleanly

- `touch $RUNDIR/STOP`: at the next step, pkdgrav3 writes a checkpoint and a snapshot and exits. It deletes the file.
- `iWallRunTime = <minutes>` in the parameter file does the same automatically. Set it below the PBS wall
  time by enough to copy the last checkpoint and snapshot (2048^3: about 2 h, because the checkpoint is
  412 GB to home at 0.3-0.5 GB/s and the last snapshot about 1 h to NFS).

### 3.7 Example: 2048^3, 200 Mpc/h, z = 49 to 0, 100 snapshots, one pawna node

Needs the user's approval before submitting (runs above 512^3). Parameter file `n2048.par`:

```python
from accuracy import classic_theta_switch, classic_replicas_switch
achOutName      = "snap"
achTfFile       = "euclid_z0_transfer_combined.dat"
nGrid           = 2048
dBoxSize        = 200.0
iLPT            = 2
iSeed           = 314159265
dRedFrom        = 49
dRedTo          = 0.0
nSteps          = 2000
iOutInterval    = 20          # 100 snapshots
iCheckInterval  = 200         # about every 10-12 h; 412 GB each
iWallRunTime    = 14280       # minutes: 240 h wall time minus 2 h for the final copies
h               = 0.6736
dOmega0         = 0.3137721
dLambda         = 0.6862279
dSigma8         = 0.81
dSpectral       = 0.965
bComove         = True
bPeriodic       = True
bEwald          = True
bEpsAccStep     = True
bNewKDK         = True
dEta            = 0.3
dTheta          = classic_theta_switch()
nReplicas       = classic_replicas_switch()
bMemUnordered   = False       # persistent IDs (iOrder) without the +128 GB of bMemParticleID
bMemParticleID  = False
bDoDensity      = False
bCompressedOutput   = True
dCompPotStep        = 100.0   # optional
dCompPotMaxRedshift = 1.0     # optional: potential only from z = 1 on
```

```bash
qsub -l select=1:ncpus=192:mem=750gb -l walltime=240:00:00 \
     -v PKD_IO_ROOT=$ROOT,RUNDIR=$PWD,PAR=n2048.par,CHKSTAGE=1,CHKDEST=/user/$USER/chk/n2048,KEEPCHK=2 $S
# if the wall time ends first: restart from the newest checkpoint
qsub -l select=1:ncpus=192:mem=750gb -l walltime=240:00:00 \
     -v PKD_IO_ROOT=$ROOT,RUNDIR=$PWD,CHKSTAGE=1,CHKDEST=/user/$USER/chk/n2048,RESTART=$(ls /user/$USER/chk/n2048/*.chk | sort | tail -1) $S
```

Expected (estimate in [18]): about 8-8.5 days with staged checkpoints; about 10 days with checkpoints written
straight to NFS. Disk: about 7 TB of snapshots (+0.4-0.9 TB with the potential). Local SSD needed: one
checkpoint and a few snapshots (about 0.6 TB of the 1.4 TB free). Ask for 11-12 days of PBS wall time to
finish in one job; otherwise use the restart line above.

### 3.8 After the run

```bash
B=$ROOT/build_snapcomp
$B/pkc_check --threads 64 snapshots/snap.00100                          # IDs once, subsamples, potential step
$B/pkc_check --threads 64 --ref snapshots/snap.00020 snapshots/snap.00100   # same 1% set, displacements
$B/pkc_to_gadget --threads 64 [--percents 1] snapshots/snap.00100 /var/tmp/$USER/gadget/snap.00100   # for Rockstar; 32 B/particle
$B/pkc_convert --threads 64 old_snapshot.hdf5 converted/snap.00100      # existing HDF5 snapshot to .pkc
$ROOT/venv/bin/python $ROOT/pkdgrav3/bench/phase5/cmp_runs.py runA/snap.00040 runB/snap.00040   # two runs, particle by particle
```

Python:

```python
import sys; sys.path.insert(0, "/storage/shadab/pkdgrav3_io/pkdgrav3/tools/pkcreader")
from pkcreader import PkcSnapshot
s = PkcSnapshot("snapshots/snap.00100")      # s.redshift, s.Lbox, s.npart, s.mpart, s.has_potential, s.pot_step
d = s.read(subsamples=[1])                   # pos (n,3) Mpc/h in [0, L), vel (n,3) km/s peculiar, ids
d = s.read(fields=("pos", "pot", "ids"))     # pot: peculiar potential in (km/s)^2 (version-2 files)
for pos, vel, ids in s.iter_blocks():        # streaming
    ...
```

The stored potential includes pkdgrav3's constant softened self term, -(105/64) m/eps. At 2048^3 in
200 Mpc/h this is -293/a (km/s)^2, the same for every particle. Subtract it before comparing with codes that
leave it out [17].

### 3.9 Recovering after a crash

If the job ended before the drain finished, `staged.<job>.log` prints the command to run. Run it as a job on the same node:

```bash
qsub -l select=1:host=node4:ncpus=1 -l walltime=04:00:00 -- /usr/bin/python3 $ROOT/pkdgrav3/tools/staging/pkd_drain.py \
  --recover --stage /var/tmp/${USER}_stage_n2048 --dest $PWD/snapshots --chk-dest /user/$USER/chk/n2048 --report $PWD/recover.log
```

A checkpoint is complete only if its restart file (`snap.NNNNN.chk`, without a suffix) is at the destination. The
drain copies that file last.

## 4. Tests

| Test | Command | Result |
|---|---|---|
| Codec unit tests (round trip, errors, potential v1/v2) | `build_snapcomp/snapcomp_test` | all pass |
| Drain with a fake pkdgrav3 (snapshots + checkpoints, keep-2, interrupted checkpoint, recover) | `python3 tools/staging/test_drain.py` | ALL OK (Mac and pawna) |
| Reader against HDF5 of the same state | `tools/pkcreader/test_pkcreader.py snap.NNNNN [gadget_prefix]` | all pass [16, 17] |
| Checkpoint/restart on pawna, 512^3 (staged checkpoints to home, STOP, restart, comparison with an uninterrupted run) | `bench/phase5/chk512.pbs` | [19] |
