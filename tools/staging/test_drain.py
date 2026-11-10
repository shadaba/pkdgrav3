#!/usr/bin/env python3
"""Test pkd_drain.py with a fake pkdgrav3 (no simulation): snapshots and checkpoints of the same steps,
checkpoints to a separate destination with --keep-chk, an interrupted checkpoint, and --recover.

usage: test_drain.py [scratch_dir]      (Python standard library only; runs in a few seconds)
"""
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import time

here = os.path.dirname(os.path.abspath(__file__))
DRAIN = os.path.join(here, "pkd_drain.py")
fails = 0


def check(cond, what):
    global fails
    print(("  ok    " if cond else "  FAIL  ") + what)
    fails += not cond


def blob(path, size, seed):
    with open(path, "wb") as f:
        f.write(hashlib.sha256(str(seed).encode()).digest() * (size // 32 + 1))


def digest(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


FAKE = r'''
import os, sys, time
from hashlib import sha256
stage, log, interrupted = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
def blob(path, size, seed):
    with open(path, "wb") as f:
        f.write(sha256(str(seed).encode()).digest() * (size // 32 + 1))
with open(log, "w", buffering=1) as L:
    for step in (10, 20, 30, 40):
        # HDF5-style snapshot via the log, then a compressed one with a .done marker
        p = f"{stage}/current/snap.{step:05d}"
        L.write(f"Writing {p} in HDF5 format serially ...\n"); blob(p, 300000 + step, step)
        L.write("Output file has been successfully written, Wallclock: 0.1 secs.\n")
        q = f"{stage}/current/c.{step:05d}"
        for s in (1, 99):
            blob(f"{q}.s{s:02d}.0.pkc", 50000 + s, step * 100 + s)
        open(q + ".done", "w").write("{}\n")
        # checkpoint: data pieces, .pkl, restart file last (as MSR::Checkpoint / writeParameters)
        c = f"{stage}/chk/snap.{step:05d}.chk"
        L.write(f"Writing checkpoint for Step: {step} Time:0.1 Redshift:1\n")
        for i in range(3):
            blob(f"{c}.{i}", 200000 + i, step * 10 + i)
        blob(c + ".pkl", 3000, step)
        open(c, "w").write("import PKDGRAV as msr\nmsr.restore(__file__)\n")
        L.write("Checkpoint has been successfully written, Wallclock: 0.2 secs.\n")
        time.sleep(0.5)
    if interrupted:  # killed in the middle of a checkpoint: no .pkl, no restart file
        c = f"{stage}/chk/snap.00050.chk"
        L.write("Writing checkpoint for Step: 50 Time:0.1 Redshift:1\n")
        blob(f"{c}.0", 200000, 500)
'''


def run(root, interrupted):
    stage, dest, chk = (os.path.join(root, d) for d in ("stage", "dest", "chkdest"))
    log = os.path.join(root, "run.log")
    rep = os.path.join(root, "drain.log")
    subprocess.check_call([sys.executable, DRAIN, "--stage", stage, "--dest", dest, "--chk-dest", chk, "--init-only"])
    fake = subprocess.Popen([sys.executable, "-c", FAKE, stage, log, "1" if interrupted else "0"])
    d = subprocess.Popen([sys.executable, DRAIN, "--stage", stage, "--dest", dest, "--chk-dest", chk, "--keep-chk", "2",
                          "--log", log, "--pid", str(fake.pid), "--poll", "0.2", "--report", rep, "--min-free-gb", "0"])
    fake.wait()
    rc = d.wait(timeout=60)
    return stage, dest, chk, rc, rep


def main(scratch):
    root = tempfile.mkdtemp(prefix="drain_test_", dir=scratch)
    try:
        print("normal run: 4 snapshots (HDF5-style + compressed) and 4 checkpoints, --keep-chk 2")
        os.makedirs(os.path.join(root, "a"))
        stage, dest, chk, rc, rep = run(os.path.join(root, "a"), interrupted=False)
        check(rc == 0, f"drain exit code 0 (got {rc})")
        snaps = sorted(os.listdir(dest))
        check(all(".chk" not in n for n in snaps), "no checkpoint file in the snapshot destination")
        check(sum(n.startswith("snap.") for n in snaps) == 4 and sum(n.endswith(".pkc") for n in snaps) == 8,
              f"all snapshot files in dest ({len(snaps)} files)")
        for step in (10, 40):
            want = os.path.join(root, "w")
            blob(want, 300000 + step, step)
            check(digest(os.path.join(dest, f"snap.{step:05d}")) == digest(want), f"snap.{step:05d} identical")
        ck = sorted(os.listdir(chk))
        check(ck == sorted([f"snap.{s:05d}.chk{x}" for s in (30, 40) for x in ("", ".0", ".1", ".2", ".pkl")]),
              "checkpoint destination holds exactly checkpoints 30 and 40 (5 files each); 10 and 20 deleted")
        blob(os.path.join(root, "w"), 200001, 401)
        check(digest(os.path.join(chk, "snap.00040.chk.1")) == digest(os.path.join(root, "w")), "checkpoint piece identical")
        check(not os.listdir(os.path.join(stage, "local")) and not os.listdir(os.path.join(stage, "local_chk")),
              "local staging directories empty")
        text = open(rep).read()
        order = [l.split()[2] for l in text.splitlines() if " copied " in l and "snap.00040.chk" in l]
        check(order and os.path.basename(order[-1]) == "snap.00040.chk", "restart file copied last")
        check(text.count("deleted old checkpoint") == 2, "two old checkpoints deleted (log)")

        print("interrupted checkpoint at the end, then --recover")
        os.makedirs(os.path.join(root, "b"))
        stage, dest, chk, rc, rep = run(os.path.join(root, "b"), interrupted=True)
        check(rc == 2, f"drain exit code 2: something left on local disk (got {rc})")
        left = os.listdir(os.path.join(stage, "local_chk"))
        check(left == ["snap.00050.chk.0"], f"only the interrupted piece left locally ({left})")
        check(not os.path.exists(os.path.join(chk, "snap.00050.chk")), "no restart file for the interrupted checkpoint")
        r = subprocess.call([sys.executable, DRAIN, "--stage", stage, "--dest", dest, "--chk-dest", chk, "--recover",
                             "--report", rep, "--poll", "0.2"])
        check(r == 0 and os.path.exists(os.path.join(chk, "snap.00050.chk.0")), "--recover copied the piece")
        check(sorted(n for n in os.listdir(chk) if n.endswith(".chk")) == ["snap.00030.chk", "snap.00040.chk"],
              "recovered piece did not count as a complete checkpoint")
    finally:
        shutil.rmtree(root)
    print("ALL OK" if not fails else f"FAILED: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
