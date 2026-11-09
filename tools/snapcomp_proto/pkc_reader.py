"""Reader for snapcomp (.pkc) files, written from docs/snapcomp_format.md only (numpy + zstandard).

    from pkc_reader import PKCFile, read_snapshot
    f = PKCFile("snap.00100.s01.0.pkc")
    f.header["box_size"], f.header["redshift"]
    pos, vel, ids = f.read()                 # Mpc/h in [0, L), peculiar km/s, uint64
    pos, vel, ids = read_snapshot("snap.00100", percents=[1, 3])   # several files of one snapshot

Positions are returned in Mpc/h in [0, box_size) (cell-centre convention of the format); velocities in
km/s. Use `units="box"` for box units in [-0.5, 0.5), as pkdgrav3 uses internally.
"""
import glob
import struct

import numpy as np
import zstandard

HEADER_SIZE = 4096
_HEADER_MAGIC = b"PKCOMP\x00\x01"
_TRAILER_MAGIC = b"PKCIDX\x00\x01"
_BLOCK_MAGIC = b"PKB1"
ID_RAW, ID_DELTA_IN_CELL, ID_LAGRANGIAN = 0, 1, 2

_FIELDS = [  # (offset, struct code, name)
    (16, "Q", "n_particles"), (24, "Q", "n_blocks"), (32, "Q", "lattice_n"), (40, "Q", "n_total"),
    (48, "I", "ncell"), (52, "I", "sub"), (56, "d", "vel_step"), (64, "d", "vel_to_kms"),
    (72, "d", "box_size"), (80, "d", "redshift"), (88, "d", "time"), (96, "d", "omega_m"),
    (104, "d", "omega_l"), (112, "d", "h"), (120, "I", "subsample"), (124, "I", "n_subsamples"),
    (256, "Q", "subsample_seed"), (264, "I", "writer"), (268, "I", "n_writers"), (272, "I", "block_size"),
    (276, "i", "zstd_level"), (280, "Q", "step"),
]


def parse_header(b):
    if b[:8] != _HEADER_MAGIC:
        raise ValueError("not a snapcomp file")
    version, hsize = struct.unpack_from("<II", b, 8)
    if version != 1 or hsize != HEADER_SIZE:
        raise ValueError(f"unsupported version {version} / header size {hsize}")
    h = {name: struct.unpack_from("<" + code, b, off)[0] for off, code, name in _FIELDS}
    h["percents"] = list(struct.unpack_from("<16d", b, 128))[: h["n_subsamples"]]
    h["info"] = b[288:].split(b"\x00", 1)[0].decode(errors="replace")
    h["q"] = h["ncell"] * h["sub"]
    return h


def _unzigzag(u):
    u = u.astype(np.uint64)
    return (u >> np.uint64(1)).astype(np.int64) ^ -(u & np.uint64(1)).astype(np.int64)


_DCTX = zstandard.ZstdDecompressor()


def decode_block(h, b):
    """Decode one block: integer positions X (n,3), velocities V (n,3) in vel_step units, IDs (n)."""
    if b[:4] != _BLOCK_MAGIC:
        raise ValueError("bad block magic")
    n, n_cells, ns, coder = struct.unpack_from("<IIBB", b, 4)
    streams, off, data = [], 16, 16 + 12 * ns
    for _ in range(ns):
        w, kind, _, count, cs = struct.unpack_from("<BBHII", b, off)
        raw = _DCTX.decompress(b[data:data + cs], max_output_size=count * w)
        if len(raw) != count * w:
            raise ValueError("stream size mismatch")
        planes = np.frombuffer(raw, np.uint8).reshape(w, count)
        v = np.zeros(count, np.uint64)
        for k in range(w):
            v |= planes[k].astype(np.uint64) << np.uint64(8 * k)
        streams.append(v)
        off += 12
        data += cs
    nc, sub, Q = h["ncell"], h["sub"], h["q"]
    keys = np.cumsum(streams[0].astype(np.int64))
    counts = streams[1].astype(np.int64)
    key = np.repeat(keys, counts)
    cell = np.stack([key % nc, (key // nc) % nc, key // (nc * nc)], 1)
    X = cell * sub + np.stack([streams[2 + d].astype(np.int64) for d in range(3)], 1)
    V = np.stack([_unzigzag(streams[5 + d]) for d in range(3)], 1)
    if coder == ID_RAW:
        ids = streams[8]
    elif coder == ID_DELTA_IN_CELL:
        d = _unzigzag(streams[8]).astype(np.uint64)  # wrapping differences
        start = np.concatenate(([0], np.cumsum(counts)[:-1]))
        first = np.cumsum(d[start], dtype=np.uint64)
        within = d.copy()
        within[start] = 0
        cw = np.cumsum(within, dtype=np.uint64)
        ids = np.repeat(first, counts) + (cw - np.repeat(cw[start], counts))
    elif coder == ID_LAGRANGIAN:
        N = h["lattice_n"]
        l = (X * N) // Q  # X < 2^32 and N < 2^21 fit in int64
        q = (np.stack([_unzigzag(streams[8 + d]) for d in range(3)], 1) + l) % N
        ids = ((q[:, 2] * N + q[:, 1]) * N + q[:, 0]).astype(np.uint64)
    else:
        raise ValueError(f"unknown ID coder {coder}")
    return X, V, ids.astype(np.uint64)


class PKCFile:
    def __init__(self, path):
        self.path = path
        with open(path, "rb") as f:
            data = f.read()
        if len(data) < HEADER_SIZE + 32 or data[-32:-24] != _TRAILER_MAGIC:
            raise ValueError(f"{path}: no trailer (incomplete file)")
        idx, nb, npart = struct.unpack_from("<QQQ", data, len(data) - 24)
        self.header = parse_header(data[:HEADER_SIZE])
        if self.header["n_blocks"] != nb or self.header["n_particles"] != npart:
            raise ValueError(f"{path}: header and trailer disagree")
        ent = np.frombuffer(data, np.uint8, nb * 32, idx).reshape(nb, 32)
        self.index = [struct.unpack_from("<QQI", ent[i].tobytes(), 0) for i in range(nb)]  # (offset, size, n)
        self._data = data

    def blocks(self):
        for off, size, n in self.index:
            X, V, ids = decode_block(self.header, self._data[off:off + size])
            if len(ids) != n:
                raise ValueError("block count mismatch")
            yield X, V, ids

    def read(self, units="mpc"):
        h = self.header
        parts = list(self.blocks())
        if not parts:
            return np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0, np.uint64)
        X = np.concatenate([p[0] for p in parts])
        V = np.concatenate([p[1] for p in parts])
        ids = np.concatenate([p[2] for p in parts])
        x = (X + 0.5) / h["q"]  # [0, 1)
        x = x - 0.5 if units == "box" else x * h["box_size"]
        return x, V * h["vel_step"], ids


def read_snapshot(prefix, percents=None, units="mpc"):
    """Read all files of one snapshot (or only the listed subsample percents) and concatenate."""
    files = sorted(glob.glob(prefix + ".s*.pkc"))
    out = []
    for p in files:
        f = PKCFile(p)
        pc = f.header["percents"][f.header["subsample"]]
        if percents is None or any(abs(pc - q) < 1e-9 for q in percents):
            out.append(f.read(units))
    if not out:
        raise FileNotFoundError(prefix)
    return tuple(np.concatenate([o[k] for o in out]) for k in range(3))
