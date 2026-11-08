#!/usr/bin/env python3
"""Option A applied to each Sahyadri subsample separately: bytes per particle per subsample file.
usage: run_subsamples.py [--dv 1] [--sub 256] snap ..."""
import argparse

import snapcomp as sc

p = argparse.ArgumentParser()
p.add_argument("snaps", nargs="+")
p.add_argument("--dv", type=float, default=1.0)
p.add_argument("--sub", type=int, default=256)
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
for path in a.snaps:
    s = sc.load(path)
    X, V, Q = sc.quantize(s, a.dv, a.sub)
    k = sc.subsample_of(s["ids"], seed=a.seed)
    full = sc.option_A(s, X, V, Q, a.sub, check=False)["total"]
    print(f"== {path} z={s['z']:.2f}: one file {full:.2f} B/particle")
    tot = 0.0
    for i, f in enumerate(sc.SAHYADRI_FRACTIONS):
        m = k == i
        t = dict(s, ids=s["ids"][m])
        r = sc.option_A(t, X[m], V[m], Q, a.sub, check=True)
        tot += r["total"] * m.sum()
        print(f"  {f:3d}%: {m.sum():9d} particles ({100*m.mean():5.2f}%)  pos {r['pos']:.2f} vel {r['vel']:.2f}"
              f" id {r['id']:.2f}  total {r['total']:.2f}")
    print(f"  all subsample files: {tot/len(k):.2f} B/particle ({100*(tot/len(k)/full-1):+.0f}% vs one file)")
