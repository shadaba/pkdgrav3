#!/usr/bin/env python3
"""Drain pkdgrav3 snapshots from node-local staging disk to shared storage (NFS).

pkdgrav3 (3.5.1+) writes snapshots to  achIoPath = "<stage>/current/{name}.{step:05d}{type}"
while logs and checkpoints stay in the run directory on NFS. This process runs next to
pkdgrav3 in the same PBS job and:

  1. detects finished snapshots, either from the pkdgrav3 log ("Writing <path> in ..."
     followed by "Output file has been successfully written") or from a "<path>.done"
     marker file (for writers that create one);
  2. copies each one (the file and any "<path>.<suffix>" pieces) to <dest> as "<file>.part",
     throttled to --bwlimit MB/s, fsyncs it, checks the size, renames it into place,
     fsyncs the directory, and only then deletes the local copy;
  3. protects the node's shared local disk: <stage>/current is a symlink to
     <stage>/local; when local free space drops below --min-free-gb it is switched to <dest>,
     so later snapshots go straight to NFS (slower, but the node never fills up). It is switched
     back once space recovers. pkdgrav3 resolves the path afresh on every write.
  4. exits once the watched pkdgrav3 process (--pid) has ended and everything finished has
     been drained. Exit code 0 = all drained; 2 = something is left on local disk.

Python standard library only (pawna has Python 3.10). One line per event in --report.
"""
import argparse
import os
import re
import shutil
import signal
import sys
import time

WRITING = re.compile(r"Writing (\S+) in \S+ format")
WRITTEN = re.compile(r"Output file has been successfully written")
CHUNK = 16 << 20  # 16 MiB copies: large sequential NFS writes


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S")


class Drainer:
    def __init__(self, a):
        self.a = a
        self.local = os.path.join(a.stage, "local")
        self.current = os.path.join(a.stage, "current")
        self.pending = []      # finished snapshots (paths as pkdgrav3 wrote them), oldest first
        self.seen = set()
        self.log_pos = 0
        self.last_writing = None
        self.stop = False
        self.rep = open(a.report, "a", buffering=1) if a.report else sys.stdout
        os.makedirs(self.local, exist_ok=True)
        os.makedirs(a.dest, exist_ok=True)
        if not os.path.islink(self.current):
            self.point_current(self.local)

    def say(self, *msg):
        print(now(), *msg, file=self.rep, flush=True)

    # ---- staging target switch -------------------------------------------------------
    def point_current(self, target):
        tmp = self.current + ".new"
        if os.path.lexists(tmp):
            os.unlink(tmp)
        os.symlink(os.path.abspath(target), tmp)
        os.replace(tmp, self.current)  # atomic swap
        self.say("current ->", target)

    def check_space(self):
        free_gb = shutil.disk_usage(self.local).free / 2**30
        to_nfs = os.path.realpath(self.current) == os.path.realpath(self.a.dest)
        if not to_nfs and free_gb < self.a.min_free_gb:
            self.say(f"local free {free_gb:.0f} GB < {self.a.min_free_gb} GB: new snapshots go to NFS")
            self.point_current(self.a.dest)
        elif to_nfs and free_gb > 2 * self.a.min_free_gb:
            self.point_current(self.local)

    # ---- detect finished snapshots ----------------------------------------------------
    def scan_log(self):
        if not self.a.log or not os.path.exists(self.a.log):
            return
        with open(self.a.log, "r", errors="replace") as f:
            f.seek(self.log_pos)
            data = f.read()
            # only consume complete lines
            cut = data.rfind("\n") + 1
            self.log_pos += len(data[:cut].encode(errors="replace"))
        for line in data[:cut].splitlines():
            m = WRITING.search(line)
            if m:
                self.last_writing = m.group(1)
            elif WRITTEN.search(line) and self.last_writing:
                self.add(self.last_writing)
                self.last_writing = None

    def scan_markers(self):
        for d in (self.local,):
            for n in sorted(os.listdir(d)):
                if n.endswith(".done"):
                    self.add(n[:-5], marker=os.path.join(d, n))

    def add(self, path, marker=None):
        path = os.path.basename(path)  # located later in local/ or dest/ (symlink may have flipped)
        if path in self.seen:
            # Seen before (e.g. from the log); a marker arriving later must not be left behind.
            if marker:
                for k, (p, m) in enumerate(self.pending):
                    if p == path:
                        self.pending[k] = (p, marker)
                        break
                else:  # already drained
                    if os.path.exists(marker):
                        os.unlink(marker)
            return
        self.seen.add(path)
        self.pending.append((path, marker))
        self.say("finished", path)

    # ---- copy ---------------------------------------------------------------------------
    def pieces(self, path):
        d, base = os.path.split(path)
        if not os.path.isdir(d):
            return []
        return sorted(os.path.join(d, n) for n in os.listdir(d)
                      if (n == base or n.startswith(base + ".")) and not n.endswith((".part", ".done")))

    def copy_one(self, src, dst):
        part = dst + ".part"
        size = os.path.getsize(src)
        t0 = time.time()
        done = 0
        with open(src, "rb") as fi, open(part, "wb") as fo:
            while True:
                buf = fi.read(CHUNK)
                if not buf:
                    break
                fo.write(buf)
                done += len(buf)
                if self.a.bwlimit > 0:  # simple rate limit
                    ahead = done / (self.a.bwlimit * 2**20) - (time.time() - t0)
                    if ahead > 0:
                        time.sleep(ahead)
            fo.flush()
            os.fsync(fo.fileno())
        if os.path.getsize(part) != size:
            raise IOError(f"size mismatch {part}: {os.path.getsize(part)} != {size}")
        os.replace(part, dst)
        dfd = os.open(os.path.dirname(dst), os.O_RDONLY)
        try:
            os.fsync(dfd)
        finally:
            os.close(dfd)
        dt = time.time() - t0
        self.say(f"copied {src} -> {dst} {size/2**30:.2f} GiB {dt:.1f} s {size/2**20/max(dt,1e-9):.0f} MB/s")
        return size, dt

    def drain_one(self, name, marker):
        files = self.pieces(os.path.join(self.local, name))
        if not files:
            if self.pieces(os.path.join(self.a.dest, name)):
                self.say("already on NFS", name)  # written while current -> dest
            else:
                self.say("WARNING nothing found for", name)
            return True
        total = 0
        t0 = time.time()
        for src in files:
            dst = os.path.join(self.a.dest, os.path.basename(src))
            size, _ = self.copy_one(src, dst)
            total += size
        for src in files:      # delete only after every piece is safely on NFS
            os.unlink(src)
        if marker and os.path.exists(marker):
            os.unlink(marker)
        dt = time.time() - t0
        self.say(f"drained {name} {len(files)} file(s) {total/2**30:.2f} GiB "
                 f"{dt:.1f} s {total/2**20/max(dt,1e-9):.0f} MB/s")
        return True

    # ---- main loop ---------------------------------------------------------------------
    def alive(self):
        if not self.a.pid:
            return not self.stop
        try:
            os.kill(self.a.pid, 0)
            return not self.stop
        except ProcessLookupError:
            return False
        except PermissionError:
            return True

    def run(self):
        self.say(f"start stage={self.a.stage} dest={self.a.dest} log={self.a.log} pid={self.a.pid} "
                 f"bwlimit={self.a.bwlimit} MB/s min_free={self.a.min_free_gb} GB")
        while True:
            running = self.alive()
            self.scan_log()
            self.scan_markers()
            self.check_space()
            while self.pending:
                path, marker = self.pending[0]
                try:
                    self.drain_one(path, marker)
                    self.pending.pop(0)
                except OSError as e:
                    self.say("ERROR", path, e, "- retry in 60 s")
                    time.sleep(60)
                    break
            if not running:
                break
            time.sleep(self.a.poll)
        self.scan_log()
        self.scan_markers()
        for path, marker in list(self.pending):
            self.drain_one(path, marker)
        left = [os.path.join(self.local, n) for n in os.listdir(self.local)]
        if left:
            self.say("LEFT ON LOCAL DISK (incomplete or unrecognised):", " ".join(left))
            return 2
        self.say("all drained")
        return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", required=True, help="node-local staging root (creates local/ and current)")
    p.add_argument("--dest", required=True, help="final directory on shared storage")
    p.add_argument("--log", help="pkdgrav3 stdout log to watch for finished writes")
    p.add_argument("--pid", type=int, help="pkdgrav3 process to wait for")
    p.add_argument("--bwlimit", type=float, default=0, help="copy rate limit in MB/s (0 = none)")
    p.add_argument("--min-free-gb", type=float, default=200, help="switch new writes to NFS below this")
    p.add_argument("--poll", type=float, default=10, help="seconds between checks")
    p.add_argument("--report", help="append events here (default stdout)")
    p.add_argument("--init-only", action="store_true", help="create stage dirs and symlink, then exit")
    p.add_argument("--recover", action="store_true",
                   help="after pkdgrav3 has ended: drain every file left in <stage>/local, then exit "
                        "(a snapshot interrupted mid-write is copied too; check it)")
    a = p.parse_args()
    d = Drainer(a)
    if a.init_only:
        return 0
    if a.recover:
        for n in sorted(os.listdir(d.local)):
            if not n.endswith((".part", ".done")):
                d.add(n)
        d.stop = True  # no pid to wait for: drain once and exit
    signal.signal(signal.SIGTERM, lambda *_: setattr(d, "stop", True))
    return d.run()


if __name__ == "__main__":
    sys.exit(main())
