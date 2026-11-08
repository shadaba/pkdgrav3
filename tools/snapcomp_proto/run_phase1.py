#!/usr/bin/env python3
"""Phase 1: compare Option A and Option B on pkdgrav3 HDF5 snapshots.

usage: run_phase1.py [--dv 1.0] [--sub 256] [--block 262144] [--level 3] [--json out.json] snap ...
Prints one table per snapshot: bytes per particle per stream, errors, encode speed.
"""
import argparse
import json
import sys
import time

import numpy as np
import zstandard

import snapcomp as sc


def lossless_ref(s, level):
    """float32 positions/velocities + uint64 IDs, byte-shuffled, zstd: what lossless gets."""
    c = zstandard.ZstdCompressor(level=level)
    tot = 0
    for a in (s["x"], s["v"], s["ids"]):
        b = np.ascontiguousarray(a)
        w = b.dtype.itemsize
        tot += len(c.compress(b.view(np.uint8).reshape(-1, w).T.copy().tobytes()))
    return tot / len(s["ids"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("snaps", nargs="+")
    p.add_argument("--dv", type=float, default=1.0, help="velocity quantum, km/s")
    p.add_argument("--sub", type=int, default=256, help="position steps per mean inter-particle spacing")
    p.add_argument("--block", type=int, default=1 << 18)
    p.add_argument("--level", type=int, default=3, help="zstd level")
    p.add_argument("--nocheck", action="store_true")
    p.add_argument("--json")
    a = p.parse_args()
    results = []
    for path in a.snaps:
        t0 = time.time()
        s = sc.load(path)
        X, V, Q = sc.quantize(s, a.dv, a.sub)
        err = sc.errors(s, X, V, Q, a.dv)
        n = len(s["ids"])
        ref = dict(raw=32.0, hdf5_gzip9=sum(s["hdf5_bytes"].values()) / n, lossless_zstd=lossless_ref(s, a.level))
        A = sc.option_A(s, X, V, Q, a.sub, a.block, a.level, check=not a.nocheck)
        B = sc.option_B(s, X, V, Q, a.sub, a.level, check=not a.nocheck)
        B0 = sc.option_B_plain(s, X, V, Q, a.sub, a.level)
        vrms = float(np.sqrt((s["v"].astype(np.float64) ** 2).mean()) * s["vkms"])
        r = dict(snap=path, N=s["N"], L=s["L"], z=s["z"], dv=a.dv, sub=a.sub, level=a.level, block=a.block,
                 vrms_kms=vrms, err=err, ref=ref, A=A, B=B, B_plain=B0, seconds=time.time() - t0)
        results.append(r)
        print(f"\n== {path}  N={s['N']}^3  L={s['L']} Mpc/h  z={s['z']:.2f}  "
              f"dx={s['L']*1000/Q:.3f} kpc/h  dv={a.dv} km/s  (3-D vel rms {vrms:.0f} km/s)  zstd-{a.level}")
        print(f"   errors: pos max {err['pos_max_kpch']:.3f} rms {err['pos_rms_kpch']:.3f} kpc/h;"
              f"  vel max {err['vel_max_kms']:.3f} rms {err['vel_rms_kms']:.3f} km/s")
        print(f"   refs B/particle: raw {ref['raw']:.1f}  hdf5(gzip9) {ref['hdf5_gzip9']:.2f}"
              f"  lossless float+zstd {ref['lossless_zstd']:.2f}")
        print(f"   {'scheme':<22}{'pos':>7}{'vel':>7}{'id':>7}{'total':>8}{'x vs raw':>10}")
        for name, R in (("A (id: " + A["id_best"] + ")", A), ("B Lorenzo", B), ("B no prediction", B0)):
            print(f"   {name:<22}{R['pos']:7.2f}{R['vel']:7.2f}{R['id']:7.2f}{R['total']:8.2f}{32/R['total']:9.1f}x")
        print("   A id variants: " + "  ".join(f"{v} {A['id_'+v]:.2f}" for v in sc.ID_VARIANTS)
              + f"   | zstd encode {A['enc_MBps']:.0f} MB/s (A), {B['enc_MBps']:.0f} MB/s (B), one core")
        sys.stdout.flush()
    if a.json:
        with open(a.json, "w") as f:
            json.dump(results, f, indent=1, default=float)


if __name__ == "__main__":
    main()
