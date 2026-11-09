# Persistent particle IDs without the extra memory (2026-11-09)

Problem: a 2048^3 run without `bMemParticleID` needs about 742 GB (fits one 755 GB pawna node); with it,
about 870 GB (does not fit; swapping makes it very slow).

## Why the ID costs more than 8 bytes
pkdgrav3 3.5.1 packs each particle into one record whose size is rounded up to its largest alignment
(8 bytes here). For a dark-matter run with `bNewKDK = True`:

| layout | fields (bytes) | record |
|---|---|---|
| `bMemUnordered = True`, no ID | header 4 + velocity 12 + position 24 + new-rung 1 = 41 | 48 |
| `bMemUnordered = True`, `bMemParticleID = True` (production) | 41 + ID 8 = 49 | **56** |
| `bMemUnordered = False` (ordered), no ID | header 8 (**contains iOrder** and the new rung) + 24 + 12 = 44 | **48** |
| ordered + `bMemIntegerPosition = True` | 8 + 12 + 12 = 32 | **32** |

The ID makes every record 8 bytes larger, and the extra particle store (`dExtraStore`, default 10%) and the
domain-decomposition buffers scale with the record size.

## Measured (Mac, 256^3 = 16.8M particles, ICs + 2 steps, 8 threads, peak RSS from /usr/bin/time)

| layout | peak RSS | per particle vs A | scaled to 2048^3 |
|---|---|---|---|
| A unordered, no ID | 2117.6 MB | - | (742 GB measured by the user) |
| B unordered + ID (production) | 2258.6 MB | +8.8 B | +76 GB |
| C ordered, no ID | 2120.8 MB | +0.2 B | +2 GB |
| E ordered + integer positions | 1804.5 MB | -19.6 B | -168 GB |

(At 2048^3 the user saw +128 GB for the ID, more than the +76 GB scaled here; buffers that scale with the
record size are larger in a big run. Ordered mode keeps the record size, so it avoids all of it.)

## Ordered mode gives a persistent ID for free
In ordered mode every particle header holds `iOrder` (42/43 bits, up to 4.4e12 particles), set from the
initial conditions and carried with the particle for the whole run. Verified on the Mac (64^3, 8 threads,
compressed output + HDF5 from the same state, z = 3 and z = 0):
- compressed IDs (`iOrder`) are exactly 0..N^3-1 and each particle's compressed position equals the HDF5 row
  with that index (the standard writer sorts by iOrder before writing in ordered mode) to half a step;
- the HDF5 rows are close to their lattice points (rms 6 cells at z = 3 vs about 32 for a random order), so
  for pkdgrav3-generated ICs `iOrder = (iz*N + iy)*N + ix`, the same as `ParticleID`;
- the same holds with integer positions.
The compressed writer already uses `iOrder` when `bMemParticleID` is off, and skips the `Reorder()` that the
standard writer does in ordered mode before every snapshot.

Caveat for ICs read from a file (production reads `snapshot_ics_000_v4corrected.hdf5`): `iOrder` is the
particle's position in the IC file, not the ID stored in it. If the file's IDs are its row numbers (0- or
1-based), the mapping is trivial; otherwise the IC file's ID column is the lookup table, needed once.

## How to use option 1 (ordered mode)

In the parameter file:

```python
bMemUnordered  = False
bMemParticleID = False
```

Nothing else changes; the run does not need a rebuild. To check a run: compressed output with
`pkc_check` (built in `build_snapcomp/`) on an early snapshot shows every particle within about a
cell of the lattice point its ID names, and `--ref <earlier snapshot>` shows the same IDs in the same
subsample with physical displacements; standard output: row i of the HDF5 file is the particle with
`iOrder = i`. Results of these checks on pawna (128^3, 512^3, 2048^3) are in
[compressed_io_plan.md](compressed_io_plan.md), Phase 4.

Further saving, not yet validated for physics: `bMemIntegerPosition = True` (32-bit integer positions,
resolution box/2^31, about 0.1 pc/h for 200 Mpc/h) takes the record to 32 bytes; measured -19.6 B/particle
in total (-168 GB at 2048^3). Compare a 512^3 run with and without it (power spectrum, halo mass function)
before using it in production.
