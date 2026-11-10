"""How many bytes per particle would the potential cost in the compressed format?

Emulates the codec's layout: particles in spatial order (Morton key of a 1024^3 grid, standing in for
pkdgrav3's tree order + the codec's (cell, id) sort), blocks of 2^18, inside a block the potential is
quantised with step dphi, delta-coded along the order, zigzagged, byte-shuffled, zstd level 3.
Potential converted to physical peculiar (km/s)^2: phi_code * KmPerSecUnit^2 / a.
"""
import sys, time
import h5py, numpy as np, zstandard as zstd


def morton(ix, iy, iz, bits=10):
    k = np.zeros(ix.shape, np.uint64)
    for b in range(bits):
        for d, c in enumerate((ix, iy, iz)):
            k |= ((c >> np.uint64(b)) & np.uint64(1)) << np.uint64(3 * b + d)
    return k


def encode(q, level=3):
    d = np.diff(q, prepend=np.int64(0))
    z = ((d << 1) ^ (d >> 63)).astype(np.uint64)
    w = 1 if z.max() < 256 else 2 if z.max() < 65536 else 4
    planes = z.astype({1: np.uint8, 2: np.uint16, 4: np.uint32}[w]).view(np.uint8).reshape(-1, w).T.copy()
    return len(zstd.ZstdCompressor(level=level).compress(planes.tobytes())) + 1


for path in sys.argv[1:]:
    with h5py.File(path, "r") as f:
        h = dict(f["Header"].attrs)
        U = float(f["Units"].attrs["KmPerSecUnit"])
        x = f["PartType1/Coordinates"][:].astype(np.float64)
        phi = f["PartType1/Potential"][:].astype(np.float64)
    z = float(h["Redshift"]); a = 1 / (1 + z); n = len(phi)
    phi *= U * U / a
    g = np.floor((x + 0.5) * 1024).astype(np.uint64) % np.uint64(1024)
    order = np.argsort(morton(g[:, 0], g[:, 1], g[:, 2]), kind="stable")
    p = phi[order]
    print(f"{path}: N={n}, z={z:.2f}, L={h['BoxSize']:.1f}; phi (km/s)^2 min {p.min():.3g} median {np.median(p):.3g} max {p.max():.3g}")
    raw = len(zstd.ZstdCompressor(level=3).compress(phi.astype(np.float32).tobytes()))
    print(f"  float32 + zstd, ID order: {raw / n:.2f} B/particle")
    for dphi in (1.0, 10.0, 100.0, 1000.0):
        t = time.perf_counter(); nb = 0
        for s in range(0, n, 1 << 18):
            nb += encode(np.rint(p[s:s + (1 << 18)] / dphi).astype(np.int64))
        dt = time.perf_counter() - t
        # error in the escape speed sqrt(2|phi|) at the deepest point, for scale
        vesc = np.sqrt(2 * -p.min())
        print(f"  dphi {dphi:6.0f} (km/s)^2: {nb / n:.2f} B/particle  (max error {dphi/2:.0f} (km/s)^2; "
              f"v_esc {vesc:.0f} km/s at the deepest point, error {dphi/2/vesc:.3f} km/s); python encode {dt / n * 1e9:.0f} ns/particle")
