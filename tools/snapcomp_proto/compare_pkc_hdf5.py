#!/usr/bin/env python3
"""Check a compressed snapshot written by pkdgrav3 (bCompressedOutput + bCompAlsoStandard) against the
HDF5 snapshot of the same step: every ID exactly once, every file's members in the right subsample,
position error <= half a step (+ float32 rounding of the HDF5 file), velocity error <= half a step.

usage: compare_pkc_hdf5.py snap.00010            (reads snap.00010 and snap.00010.s*.pkc)
"""
import glob
import sys

import numpy as np

import snapcomp as sc
from pkc_reader import PKCFile


def main(prefix):
    s = sc.load(prefix)
    n = len(s["ids"])
    files = sorted(glob.glob(prefix + ".s*.pkc"))
    assert files, "no .pkc files"
    xs, vs, ids, ok_member = [], [], [], True
    nwriters = set()
    for p in files:
        f = PKCFile(p)
        h = f.header
        x, v, i = f.read(units="box")
        ok_member &= bool(np.all(sc.subsample_of(i, h["percents"], h["subsample_seed"]) == h["subsample"]))
        nwriters.add(h["n_writers"])
        xs.append(x), vs.append(v), ids.append(i)
    x, v, i = np.concatenate(xs), np.concatenate(vs), np.concatenate(ids)
    order = np.argsort(s["ids"])
    sid = s["ids"][order]
    pos = np.searchsorted(sid, i)
    pos = np.minimum(pos, n - 1)
    found = sid[pos] == i
    j = order[pos]
    once = len(i) == n and found.all() and len(np.unique(i)) == n
    Q = h["q"]
    dx = (x - s["x"][j].astype(np.float64) + 0.5) % 1.0 - 0.5
    dvk = v - s["v"][j].astype(np.float64) * h["vel_to_kms"]
    tol_x = 0.5 / Q + 4e-8  # float32 rounding of HDF5 positions near |x| = 0.5
    ok = once and ok_member and np.abs(dx).max() <= tol_x and np.abs(dvk).max() <= 0.5 * h["vel_step"] * (1 + 1e-6) + 1e-4
    L = h["box_size"]
    nbytes = sum(len(open(p, "rb").read()) for p in files)
    print(f"{prefix}: z={h['redshift']:.3f} {len(files)} files ({sorted(nwriters)} writer), {len(i)} of {n} particles; IDs once: {once}; "
          f"membership: {ok_member}")
    print(f"  pos max err {np.abs(dx).max()*L*1e3:.4f} kpc/h (half step {L*1e3/Q/2:.4f}); vel max err {np.abs(dvk).max():.4f} km/s "
          f"(half step {h['vel_step']/2}); {nbytes/n:.3f} B/particle vs HDF5 file -> {'OK' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    sys.exit(0 if all([main(p) for p in sys.argv[1:]]) else 1)
