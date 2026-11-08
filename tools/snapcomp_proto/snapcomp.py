"""Phase 1 offline prototype: compression schemes for pkdgrav3 dark-matter snapshots.

Reads a pkdgrav3 HDF5 snapshot (positions in box units [-0.5, 0.5), velocities dx/dt in box units per
internal time, uint64 IDs = Lagrangian lattice index id = (iz*N + iy)*N + ix) and compares:

  Option A  blocks of BLOCK particles in file (= pkdgrav3 memory) order. Inside a block particles
            are re-ordered by (cell, id); positions = occupied cells (delta keys + counts) and 8-bit
            offsets per dimension; velocities quantized; IDs coded per block (several variants);
            each stream zstd-compressed.
  Option B  global Lagrangian (ID) order, IDs implicit; displacement and velocity fields predicted
            with a 3-D Lorenzo predictor from lattice neighbours; residuals zstd-compressed.

Both options use the same quantizers, so they have identical errors; only bytes differ:
  positions: integer X = floor(x * Q), Q = 256 * N per dimension (uint8 offset in a cell of the
             N^3 mean-spacing grid, as in the existing Gadget scheme); max error 1/(2Q) box after
             decoding to the cell centre... we decode to (X + 0.5)/Q.
  velocities: V = round(v_kms / dv) with dv in km/s (peculiar velocity).

Every encoder has a decoder; tests assert the round trip reproduces the quantized data exactly.
No pkdgrav3 code is used.
"""
import time

import h5py
import numpy as np
import zstandard

HUBBLE0 = np.sqrt(8 * np.pi / 3)  # pkdgrav3 internal H0


# ---------------------------------------------------------------------------------------------
# snapshot reading and quantization
# ---------------------------------------------------------------------------------------------
def load(path):
    with h5py.File(path, "r") as f:
        h = dict(f["Header"].attrs)
        g = f["PartType1"]
        x = g["Coordinates"][:]
        v = g["Velocities"][:]
        ids = g["ParticleIDs"][:]
        hdf5_bytes = {k: g[k].id.get_storage_size() for k in ("Coordinates", "Velocities", "ParticleIDs")}
    n = int(round(len(ids) ** (1 / 3)))
    assert n**3 == len(ids), "expect a full N^3 lattice"
    a = 1.0 / (1.0 + float(h["Redshift"]))  # header "Time" is pkdgrav3 code time, not a
    L = float(h["BoxSize"])  # Mpc/h
    vkms = a * L * 100.0 / HUBBLE0  # box units / internal time  ->  peculiar km/s
    return dict(x=x, v=v, ids=ids, N=n, L=L, a=a, z=float(h["Redshift"]), vkms=vkms, hdf5_bytes=hdf5_bytes)


def quantize(s, dv_kms=1.0, sub=256):
    """Integer positions (per dim, Q = sub*N) and velocities (units of dv_kms)."""
    Q = sub * s["N"]
    X = np.floor((s["x"].astype(np.float64) + 0.5) * Q).astype(np.int64) % Q
    V = np.rint(s["v"].astype(np.float64) * s["vkms"] / dv_kms).astype(np.int64)
    return X, V, Q


def dequantize(X, V, Q, s, dv_kms):
    x = (X + 0.5) / Q - 0.5
    v = V * dv_kms / s["vkms"]
    return x, v


def errors(s, X, V, Q, dv_kms):
    x, v = dequantize(X, V, Q, s, dv_kms)
    dx = (x - s["x"] + 0.5) % 1.0 - 0.5
    dvk = (v - s["v"]) * s["vkms"]
    Lkpc = s["L"] * 1000.0
    return dict(pos_max_kpch=np.abs(dx).max() * Lkpc, pos_rms_kpch=np.sqrt((dx**2).mean()) * Lkpc,
                vel_max_kms=np.abs(dvk).max(), vel_rms_kms=np.sqrt((dvk**2).mean()))


# ---------------------------------------------------------------------------------------------
# stream helpers
# ---------------------------------------------------------------------------------------------
def zigzag(a):
    a = a.astype(np.int64)
    return ((a << 1) ^ (a >> 63)).astype(np.uint64)


def unzigzag(u):
    u = u.astype(np.uint64)
    return ((u >> np.uint64(1)).astype(np.int64)) ^ (-(u & np.uint64(1)).astype(np.int64))


def to_bytes(u):
    """Smallest unsigned width that holds u, byte-shuffled (byte planes) for the entropy coder."""
    u = np.ascontiguousarray(u, dtype=np.uint64)
    m = int(u.max()) if u.size else 0
    w = 1 if m < 1 << 8 else 2 if m < 1 << 16 else 4 if m < 1 << 32 else 8
    t = u.astype({1: np.uint8, 2: np.uint16, 4: np.uint32, 8: np.uint64}[w])
    return w, t.view(np.uint8).reshape(-1, w).T.copy().tobytes()


def from_bytes(b, w, n):
    t = np.frombuffer(b, np.uint8).reshape(w, n).T.copy()
    return t.view({1: np.uint8, 2: np.uint16, 4: np.uint32, 8: np.uint64}[w]).reshape(n).astype(np.uint64)


class Z:
    """zstd with timing; one compressor per level."""

    def __init__(self, level):
        self.c = zstandard.ZstdCompressor(level=level)
        self.d = zstandard.ZstdDecompressor()
        self.t = 0.0
        self.raw = 0

    def enc(self, u):
        w, b = to_bytes(u)
        t0 = time.perf_counter()
        c = self.c.compress(b)
        self.t += time.perf_counter() - t0
        self.raw += len(b)
        return (w, len(u), c)

    def dec(self, item):
        w, n, c = item
        return from_bytes(self.d.decompress(c), w, n)


def nbytes(items):
    return sum(len(c) + 3 for (_, _, c) in items)  # + small per-stream header


# ---------------------------------------------------------------------------------------------
# Option A: blocks in file order
# ---------------------------------------------------------------------------------------------
ID_VARIANTS = ("raw64", "delta_in_cell", "lagr_resid")


def encode_block_A(X, V, ids, N, Q, sub, z):
    """One block. Returns dict stream_name -> list of zstd items."""
    cell = X // sub  # (n,3) cell indices on the N grid
    key = (cell[:, 2] * N + cell[:, 1]) * N + cell[:, 0]
    order = np.lexsort((ids, key))  # by cell, then id inside a cell
    key, X, V, ids, cell = key[order], X[order], V[order], ids[order], cell[order]
    ukey, start, cnt = np.unique(key, return_index=True, return_counts=True)
    out = {}
    # positions: occupied-cell keys (delta) + counts + per-particle 8-bit offsets
    out["pos"] = [z.enc(np.diff(ukey, prepend=0).astype(np.uint64)), z.enc(cnt.astype(np.uint64))]
    off = X - cell * sub
    out["pos"] += [z.enc(off[:, d].astype(np.uint64)) for d in range(3)]
    # velocities: per dimension, zigzag of quantized values (no prediction; order is by cell)
    out["vel"] = [z.enc(zigzag(V[:, d])) for d in range(3)]
    # IDs, three variants (the best one is used for the total)
    out["id_raw64"] = [z.enc(ids.astype(np.uint64))]
    d = np.diff(ids.astype(np.int64), prepend=0)
    d[start] = ids[start].astype(np.int64) - np.concatenate(([0], ids[start[:-1]].astype(np.int64)))
    out["id_delta_in_cell"] = [z.enc(zigzag(d))]
    q = np.stack([ids % N, (ids // N) % N, ids // (N * N)], 1).astype(np.int64)  # (ix,iy,iz)
    r = (q - cell + N // 2) % N - N // 2  # Lagrangian lattice minus Eulerian cell, periodic
    out["id_lagr_resid"] = [z.enc(zigzag(r[:, k])) for k in range(3)]
    return out, order


def decode_block_A(out, N, sub, z, variant):
    ukey = np.cumsum(z.dec(out["pos"][0]).astype(np.int64))
    cnt = z.dec(out["pos"][1]).astype(np.int64)
    key = np.repeat(ukey, cnt)
    cell = np.stack([key % N, (key // N) % N, key // (N * N)], 1)
    off = np.stack([z.dec(out["pos"][2 + d]).astype(np.int64) for d in range(3)], 1)
    X = cell * sub + off
    V = np.stack([unzigzag(z.dec(out["vel"][d])) for d in range(3)], 1)
    if variant == "raw64":
        ids = z.dec(out["id_raw64"][0])
    elif variant == "delta_in_cell":
        d = unzigzag(z.dec(out["id_delta_in_cell"][0]))
        start = np.concatenate(([0], np.cumsum(cnt)[:-1]))
        first = np.cumsum(d[start])  # ids of the first particle in each cell
        dd = d.copy()
        dd[start] = 0
        within = np.cumsum(dd) - np.repeat(np.cumsum(dd)[start], cnt)
        ids = (np.repeat(first, cnt) + within).astype(np.uint64)
    else:
        r = np.stack([unzigzag(z.dec(out["id_lagr_resid"][k])) for k in range(3)], 1)
        q = (r + cell) % N
        ids = ((q[:, 2] * N + q[:, 1]) * N + q[:, 0]).astype(np.uint64)
    return X, V, ids


def option_A(s, X, V, Q, sub, block=1 << 18, level=3, check=True):
    z = Z(level)
    N = s["N"]
    ids = s["ids"]
    tot = {k: 0 for k in ("pos", "vel", *("id_" + v for v in ID_VARIANTS))}
    for b0 in range(0, len(ids), block):
        sl = slice(b0, b0 + block)
        out, order = encode_block_A(X[sl], V[sl], ids[sl].astype(np.int64), N, Q, sub, z)
        for k in out:
            tot[k] += nbytes(out[k])
        if check:
            for var in ID_VARIANTS:
                Xd, Vd, idd = decode_block_A(out, N, sub, z, var)
                assert np.array_equal(Xd, X[sl][order]) and np.array_equal(Vd, V[sl][order]), "A pos/vel"
                assert np.array_equal(idd, ids[sl][order]), f"A ids {var}"
    n = len(ids)
    best = min(ID_VARIANTS, key=lambda v: tot["id_" + v])
    res = {k: tot[k] / n for k in tot}
    res.update(id=res["id_" + best], id_best=best, total=(tot["pos"] + tot["vel"] + tot["id_" + best]) / n,
               enc_MBps=z.raw / 2**20 / max(z.t, 1e-9))
    return res


# ---------------------------------------------------------------------------------------------
# Option B: global Lagrangian order, implicit IDs, Lorenzo prediction
# ---------------------------------------------------------------------------------------------
def lorenzo_residual(F):
    """F: integer field (N,N,N) indexed [iz,iy,ix], non-periodic 3-D Lorenzo predictor."""
    P = np.zeros((F.shape[0] + 1, F.shape[1] + 1, F.shape[2] + 1), np.int64)
    P[1:, 1:, 1:] = F
    pred = (P[:-1, 1:, 1:] + P[1:, :-1, 1:] + P[1:, 1:, :-1] - P[:-1, :-1, 1:] - P[:-1, 1:, :-1]
            - P[1:, :-1, :-1] + P[:-1, :-1, :-1])
    return F - pred


def lorenzo_inverse(R):
    # The Lorenzo residual is the 3-D mixed finite difference; its inverse is a 3-D prefix sum.
    return np.cumsum(np.cumsum(np.cumsum(R, 0), 1), 2)


def option_B(s, X, V, Q, sub, level=3, check=True):
    z = Z(level)
    N = s["N"]
    order = np.argsort(s["ids"], kind="stable")
    assert np.array_equal(s["ids"][order], np.arange(N**3)), "IDs must be 0..N^3-1"
    Xo, Vo = X[order], V[order]
    i = np.arange(N**3)
    lat = np.stack([i % N, (i // N) % N, i // (N * N)], 1) * sub + sub // 2  # lattice in Q units
    D = (Xo - lat + Q // 2) % Q - Q // 2  # displacement, periodic, in Q units
    out = {"pos": [], "vel": []}
    for d in range(3):
        out["pos"].append(z.enc(zigzag(lorenzo_residual(D[:, d].reshape(N, N, N)).ravel())))
        out["vel"].append(z.enc(zigzag(lorenzo_residual(Vo[:, d].reshape(N, N, N)).ravel())))
    if check:
        for d in range(3):
            Dd = lorenzo_inverse(unzigzag(z.dec(out["pos"][d])).reshape(N, N, N)).ravel()
            assert np.array_equal((Dd + lat[:, d]) % Q, Xo[:, d]), "B pos"
            Vd = lorenzo_inverse(unzigzag(z.dec(out["vel"][d])).reshape(N, N, N)).ravel()
            assert np.array_equal(Vd, Vo[:, d]), "B vel"
    n = N**3
    res = {k: nbytes(out[k]) / n for k in out}
    res.update(id=0.0, total=res["pos"] + res["vel"], enc_MBps=z.raw / 2**20 / max(z.t, 1e-9))
    return res


def option_B_plain(s, X, V, Q, sub, level=3):
    """B without prediction (displacements and velocities as they are) to show what Lorenzo buys."""
    z = Z(level)
    N = s["N"]
    order = np.argsort(s["ids"], kind="stable")
    Xo, Vo = X[order], V[order]
    i = np.arange(N**3)
    lat = np.stack([i % N, (i // N) % N, i // (N * N)], 1) * sub + sub // 2
    D = (Xo - lat + Q // 2) % Q - Q // 2
    pos = sum(nbytes([z.enc(zigzag(D[:, d]))]) for d in range(3)) / N**3
    vel = sum(nbytes([z.enc(zigzag(Vo[:, d]))]) for d in range(3)) / N**3
    return dict(pos=pos, vel=vel, id=0.0, total=pos + vel)
