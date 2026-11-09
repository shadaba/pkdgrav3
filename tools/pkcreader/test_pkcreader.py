#!/usr/bin/env python3
"""Check PkcSnapshot (and optionally pkc_to_gadget output) against the HDF5 snapshot of the same step.

usage: test_pkcreader.py snap.00006 [gadget_prefix]
  snap.00006        pkdgrav3 HDF5 file; snap.00006.s*.pkc next to it (bCompressedOutput + bCompAlsoStandard)
  gadget_prefix     files <gadget_prefix>.0 .. written by pkc_to_gadget from the same .pkc files

Checks: header values; every ID once; read() / read_block() / iter_blocks() agree and are in the same
order; position and velocity errors within half a step of the HDF5 values; subsample selection; and for
the Gadget files: header, particle mass, v/sqrt(a) convention, same particles as the reader.
"""
import glob
import os
import struct
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
sys.path.insert(0, os.path.join(here, "..", "snapcomp_proto"))
import snapcomp as sc  # noqa: E402  (HDF5 loader of the prototype)
from pkcreader import RHO_CRIT, PkcSnapshot  # noqa: E402


def check(cond, what):
    print(("  ok    " if cond else "  FAIL  ") + what)
    return bool(cond)


def read_gadget(prefix):
    """Minimal Gadget-2 (SnapFormat 1) reader, following Rockstar's io_gadget.c."""
    files = sorted(glob.glob(prefix + ".*"), key=lambda p: int(p.rsplit(".", 1)[1]))
    hs, P, V, I = [], [], [], []
    for p in files:
        with open(p, "rb") as f:
            def rec():
                (n,) = struct.unpack("<I", f.read(4))
                b = f.read(n)
                assert struct.unpack("<I", f.read(4))[0] == n, "record markers differ"
                return b
            hb = rec()
            h = dict(npart=struct.unpack_from("<6I", hb, 0), mass=struct.unpack_from("<6d", hb, 24),
                     time=struct.unpack_from("<d", hb, 72)[0], redshift=struct.unpack_from("<d", hb, 80)[0],
                     npartTotal=struct.unpack_from("<6I", hb, 96), num_files=struct.unpack_from("<i", hb, 124)[0],
                     box=struct.unpack_from("<4d", hb, 128), high=struct.unpack_from("<6I", hb, 168))
            n = h["npart"][1]
            P.append(np.frombuffer(rec(), "<f4").reshape(n, 3))
            V.append(np.frombuffer(rec(), "<f4").reshape(n, 3))
            ib = rec()
            assert len(ib) == 8 * n, "IDs must be 8 bytes"
            I.append(np.frombuffer(ib, "<u8"))
            hs.append(h)
    return hs, np.concatenate(P), np.concatenate(V), np.concatenate(I)


def main(h5, gadget=None):
    ok = True
    s = PkcSnapshot(h5)
    ref = sc.load(h5)
    n = len(ref["ids"])
    print(f"{h5}: {len(s.files)} files, N = {s.npart}, z = {s.redshift:.4f}, L = {s.Lbox}, m_p = {s.mpart:.4e} Msun/h")
    ok &= check(s.npart == n, f"particle count matches HDF5 ({n})")
    ok &= check(abs(s.mpart - RHO_CRIT * s.Om * s.Lbox**3 / n) < 1e-6 * s.mpart, "particle mass from Omega_m, L, N")

    d = s.read(dtype=np.float64)
    ids = d["ids"]
    ok &= check(len(np.unique(ids)) == n, "every ID exactly once")
    order = np.argsort(ref["ids"])
    j = order[np.searchsorted(ref["ids"][order], ids)]
    ok &= check(np.array_equal(ref["ids"][j], ids), "all IDs found in HDF5")
    h = s.header
    Q, L = h["q"], s.Lbox
    xref = (ref["x"][j].astype(np.float64) + 0.5) * L  # HDF5 box units [-0.5, 0.5) -> Mpc/h [0, L)
    dx = (d["pos"] - xref + 0.5 * L) % L - 0.5 * L
    dv = d["vel"] - ref["v"][j].astype(np.float64) * h["vel_to_kms"]
    ok &= check(np.abs(dx).max() <= 0.5 * L / Q * (1 + 1e-6) + 4e-8 * L,
                f"position error {np.abs(dx).max()*1e3:.4f} kpc/h <= half step {0.5*L/Q*1e3:.4f}")
    ok &= check(np.abs(dv).max() <= 0.5 * h["vel_step"] * (1 + 1e-6) + 1e-4,
                f"velocity error {np.abs(dv).max():.4f} km/s <= half step {0.5*h['vel_step']}")
    ok &= check(d["pos"].min() >= 0 and d["pos"].max() < L, "positions in [0, L)")

    p32 = s.read_block("pos")
    i2 = s.read_block("ids")
    ok &= check(p32.shape == (3, n) and np.array_equal(i2, ids), "read_block: (3, n) and same order as read()")
    ok &= check(np.allclose(p32.T, d["pos"], atol=1e-5 * L, rtol=0), "read_block positions = read() positions")
    it = np.concatenate([b[2] for b in s.iter_blocks()])
    ok &= check(np.array_equal(it, ids), "iter_blocks: same particles, same order")

    small = s.percents[0]
    sub = s.read(fields=("ids",), subsamples=[small])["ids"]
    want = sc.subsample_of(ids, h["percents"], h["subsample_seed"]) == 0
    ok &= check(np.array_equal(np.sort(sub), np.sort(ids[want])), f"subsample {small}%: {len(sub)} particles, the right ones")
    ok &= check(s.count([small]) == len(sub), "count() matches")

    if gadget:
        hs, P, V, I = read_gadget(gadget)
        tot = (hs[0]["high"][1] << 32) + hs[0]["npartTotal"][1]
        a = 1 / (1 + s.redshift)
        ok &= check(tot == n and sum(x["npart"][1] for x in hs) == n and hs[0]["num_files"] == len(hs),
                    f"Gadget: {len(hs)} files, {n} particles, totals and num_files consistent")
        ok &= check(abs(hs[0]["mass"][1] * 1e10 - s.mpart) < 1e-6 * s.mpart and abs(hs[0]["time"] - a) < 1e-12
                    and hs[0]["box"][0] == L and abs(hs[0]["box"][1] - s.Om) < 1e-12,
                    "Gadget header: mass (1e10 Msun/h), a, box, Omega_m")
        g = np.argsort(I)
        k = np.argsort(ids)
        ok &= check(np.array_equal(I[g], ids[k]), "Gadget: same IDs as the reader")
        ok &= check(np.abs(P[g] - d["pos"][k]).max() <= 2e-7 * L, "Gadget positions = reader positions (float32)")
        ok &= check(np.abs(V[g] * np.sqrt(a) - d["vel"][k]).max() <= 1e-5 * (1 + np.abs(d["vel"]).max()),
                    "Gadget velocities = peculiar km/s / sqrt(a)")
    print("ALL OK" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(*sys.argv[1:3]) else 1)
