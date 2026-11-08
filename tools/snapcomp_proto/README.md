# Phase 1 offline prototype

`snapcomp.py` (Option A, Option B, quantizers, decoders) and `run_phase1.py` (driver). Python with
numpy, h5py, zstandard; no pkdgrav3 code. Every encoder has a decoder and the driver checks that
each round trip reproduces the quantized data exactly (`--nocheck` skips this).

```bash
python run_phase1.py --dv 1.0 --sub 256 --json out.json snap.00040 snap.00065 snap.00090
```
- `--dv`: velocity step in km/s (peculiar velocity; max error dv/2).
- `--sub`: position steps per mean inter-particle spacing; step = L/(sub*N) (max error half a step).
- `--block`: particles per block for Option A (default 2^18), `--level`: zstd level (default 3).

Snapshot conventions verified on pkdgrav3 3.5.1 output:
- positions in box units [-0.5, 0.5); IDs are Lagrangian: id = (iz*N + iy)*N + ix, cell-centred;
- velocities are dx/dt in box units per code time; peculiar km/s = a * v * L[Mpc/h] * 100 / sqrt(8 pi/3)
  (checked against the Zel'dovich relation at z=3: slope 1.02);
- the HDF5 header attribute `Time` is code time, not the expansion factor; use a = 1/(1+Redshift).
