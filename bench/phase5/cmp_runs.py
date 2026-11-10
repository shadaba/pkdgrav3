#!/usr/bin/env python3
"""Compare the same snapshot from two runs (e.g. uninterrupted vs restarted from a checkpoint), particle by
particle (matched by ID), in quantization steps. Identical physics gives 0 everywhere.

usage: cmp_runs.py prefixA prefixB      (compressed snapshots: <prefix>.s*.pkc)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools", "pkcreader"))
from pkcreader import PkcSnapshot  # noqa: E402


def load(prefix):
    s = PkcSnapshot(prefix)
    f = ("pos", "vel", "ids") + (("pot",) if s.has_potential else ())
    d = s.read(fields=f, units="box", dtype=np.float64)
    o = np.argsort(d["ids"])
    return s, {k: v[o] for k, v in d.items()}


def main(a, b):
    sa, A = load(a)
    sb, B = load(b)
    h = sa.header
    ok = np.array_equal(A["ids"], B["ids"])
    print(f"{os.path.basename(a)}: z = {sa.redshift:.6f} vs {sb.redshift:.6f}; N = {len(A['ids'])}; same IDs: {ok}")
    if not ok:
        return 1
    dx = A["pos"] - B["pos"]
    dx -= np.round(dx)
    nx = np.rint(np.abs(dx) * h["q"]).max(1)          # position steps
    nv = np.rint(np.abs(A["vel"] - B["vel"]) / h["vel_step"]).max(1)
    rows = [("position", nx), ("velocity", nv)]
    if "pot" in A and "pot" in B:
        rows.append(("potential", np.rint(np.abs(A["pot"] - B["pot"]) / sa.pot_step)))
    same = True
    for name, n in rows:
        frac = np.mean(n > 0)
        same &= frac == 0
        print(f"  {name:9s}: differ for {100*frac:.4f}% of particles; max {int(n.max())} steps; "
              f"> 1 step: {100*np.mean(n > 1):.4f}%")
    print("  IDENTICAL" if same else "  DIFFERENT")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
