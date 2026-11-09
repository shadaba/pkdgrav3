#!/usr/bin/env python3
"""Compressed bytes per particle by stream group (pos / vel / id) and ID-coder use, from .pkc block
tables only (no decompression). Standard library only. usage: pkc_stats.py file.pkc ..."""
import struct
import sys

KIND = {1: "pos", 2: "pos", 3: "pos", 4: "vel", 5: "id", 6: "id", 7: "id"}
tot = {"pos": 0, "vel": 0, "id": 0, "hdr": 0}
coders = [0, 0, 0]
npart = 0
for path in sys.argv[1:]:
    with open(path, "rb") as f:
        f.seek(-32, 2)
        t = f.read(32)
        assert t[:8] == b"PKCIDX\x00\x01", path
        idx, nb, n = struct.unpack_from("<QQQ", t, 8)
        f.seek(idx)
        index = f.read(nb * 32)
        npart += n
        for i in range(nb):
            off, size = struct.unpack_from("<QQ", index, 32 * i)
            f.seek(off)
            head = f.read(16)
            ns, coder = struct.unpack_from("<BB", head, 12)
            coders[coder] += 1
            table = f.read(12 * ns)
            tot["hdr"] += 16
            for k in range(ns):
                kind = table[12 * k + 1]
                cs = struct.unpack_from("<I", table, 12 * k + 8)[0]
                tot[KIND[kind]] += cs + 12
print(f"{len(sys.argv)-1} files, {npart} particles: " + "  ".join(f"{k} {v/npart:.3f}" for k, v in tot.items())
      + f"  total {sum(tot.values())/npart:.3f} B/particle;  ID coder blocks raw/delta/lagr {coders[0]}/{coders[1]}/{coders[2]}")
