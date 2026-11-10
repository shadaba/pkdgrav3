# pkcreader: Python reader for compressed pkdgrav3 snapshots

`pkcreader.py` reads the `.pkc` files written with `bCompressedOutput = True`
(format: [docs/snapcomp_format.md](../../docs/snapcomp_format.md)). It needs only numpy and
`zstandard` (`pip install zstandard`), is written from the format spec alone, and is one file so it
can be copied next to Sahyadri's `scripts/post-process/readers.py`.

```python
from pkcreader import PkcSnapshot
s = PkcSnapshot("snapshots/snap.00100")   # every <prefix>.s*.pkc file of the snapshot
s.redshift, s.scale, s.Lbox, s.npart, s.mpart, s.Om, s.Olam, s.hubble, s.percents
d = s.read()                              # dict: pos (n,3) Mpc/h in [0, L), vel (n,3) km/s, ids uint64
d = s.read(fields=("pos",), subsamples=[1, 3], dtype=np.float64)
pos = s.read_block("pos", subsamples=[1]) # like SnapshotReader.read_block: (3, n); "vel" (3, n); "ids" (n,)
for pos, vel, ids in s.iter_blocks():     # streaming, <= 65536 particles at a time
    ...
```

- Cosmology: `s.Om`, `s.Olam`, `s.hubble` from the binary header (version 3 files also have Omega_b, Omega_rad,
  Omega_DE, w0, wa, sigma8, n_s, running, normalization, pivot in `s.header["cosmology"]`); the full set (w0, wa, Omega_DE, Omega_rad,
  Omega_b, sigma8, n_s, ...) and values at the snapshot (E(z), H(z), Omega_m(z), growth D and f, age, comoving
  distance) in `s.cosmology`. Everything else pkdgrav3 records is in `s.info`: `sim`, `box`, `ic`, `numerics`,
  `steps` (see docs/snapcomp_format.md). Snapshots written before 2026-11-10 give empty dicts. For w0-wa runs
  `s.Olam` is 0.
- Velocities are peculiar, in km/s (Gadget files store v/sqrt(a); this reader does not).
- `mpart` = rho_crit * Omega_m * L^3 / N in Msun/h.
- Subsamples are given by percent. A particle's subsample depends only on its ID, so it is the same in
  every snapshot. The order of particles is the same on every call.
- Decoding runs on `threads` threads (default: up to 32). Memory: the output arrays only
  (`read()` with float32: 32 B/particle; one 1% subsample of 2048^3 is 2.7 GB).
- `PkcSnapshot` refuses a snapshot whose files do not add up to the total (missing or incomplete files).

`test_pkcreader.py snap.00006 [gadget_prefix]` checks the reader (and `pkc_to_gadget` output) against
the HDF5 snapshot of the same step, written with `bCompAlsoStandard = True`.
