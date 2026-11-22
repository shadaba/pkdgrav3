#!/usr/bin/env python3
"""Drain pkdgrav3 snapshots (and checkpoints) from node-local staging disk to shared storage (NFS).

pkdgrav3 (3.5.1+) writes snapshots to  achIoPath = "<stage>/current/{name}.{step:05d}{type}"
and, optionally, checkpoints to achCheckpointPath = "<stage>/chk/{name}.{step:05d}{type}"; logs stay in
the run directory on NFS. This process runs next to pkdgrav3 in the same PBS job and:

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

Checkpoints (<stage>/chk -> <stage>/local_chk, its own directory so snapshot and checkpoint files of the
same step never mix) are found from the log ("Checkpoint has been successfully written": the restart
file <name>.<step>.chk is written last) and copied to --chk-dest (default --dest) by a second worker,
in parallel with the snapshots, so they can go to a different file server. The restart file is copied
last, so a checkpoint whose restart file is on the destination is complete. With --keep-chk K, older
checkpoints in --chk-dest are deleted once K newer complete ones are there (restart file first).

Python standard library only (pawna has Python 3.10). One line per event in --report.
"""
import argparse
import os
import re
import select
import shutil
import signal
import subprocess
import sys
import threading
import time

WRITING = re.compile(r"Writing (\S+) in \S+ format")
WRITTEN = re.compile(r"Output file has been successfully written")
CHK_WRITTEN = re.compile(r"Checkpoint has been successfully written")
CHK_NAME = re.compile(r"^(?P<name>.+)\.(?P<step>\d+)\.chk$")  # restart file of a checkpoint
CHUNK = 16 << 20  # 16 MiB copies: large sequential NFS writes

# One file copy, run as a child process by Drainer.copy_one: argv src part chunk bwlimit_MBps. Prints the bytes copied
# after every chunk (the parent's progress watch). Writes with O_DIRECT where the file system takes it (NFS: no page
# cache, so a node with little free memory does not throttle or stall the writes in its dirty-page writeback; the
# nodes of the 2048^3 runs had 2-8 GB free), else buffered with fdatasync + dropping the written pages; the source
# pages are dropped after reading (fadvise), so the copy does not push the simulation's memory into swap.
COPY_CHILD = r"""
import fcntl, mmap, os, sys, time
src, part, chunk, bw = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
fl = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
direct = hasattr(os, "O_DIRECT")
try:
    fo = os.open(part, fl | (os.O_DIRECT if direct else 0), 0o600)
except OSError:
    direct = False
    fo = os.open(part, fl, 0o600)
_t = os.environ.get("PKD_DRAIN_TEST_STALL")      # tests: the first attempt of a matching file hangs
if _t and part.endswith(".part") and _t in os.path.basename(part):
    time.sleep(3600)
buf = mmap.mmap(-1, chunk)
mv = memoryview(buf)
done = flushed = 0
t0 = time.time()
with open(src, "rb", buffering=0) as fi:
    while True:
        k = fi.readinto(mv)
        if not k:
            break
        w = 0
        while w < k:
            try:
                w += os.write(fo, mv[w:k])
            except OSError:
                if not direct:
                    raise
                # e.g. a local file system that wants aligned O_DIRECT sizes (the last piece): go on buffered
                fcntl.fcntl(fo, fcntl.F_SETFL, fcntl.fcntl(fo, fcntl.F_GETFL) & ~os.O_DIRECT)
                direct = False
        if hasattr(os, "posix_fadvise"):
            os.posix_fadvise(fi.fileno(), done, k, os.POSIX_FADV_DONTNEED)
        done += k
        if not direct and done - flushed >= (256 << 20):
            os.fdatasync(fo)
            if hasattr(os, "posix_fadvise"):
                os.posix_fadvise(fo, flushed, done - flushed, os.POSIX_FADV_DONTNEED)
            flushed = done
        print(done, flush=True)
        if bw > 0:
            ahead = done / (bw * 2**20) - (time.time() - t0)
            if ahead > 0:
                time.sleep(ahead)
os.fsync(fo)
os.close(fo)
print(done, flush=True)
"""


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S")


class Drainer:
    def __init__(self, a):
        self.a = a
        a.chk_dest = a.chk_dest or a.dest
        self.local = os.path.join(a.stage, "local")
        self.current = os.path.join(a.stage, "current")
        self.local_chk = os.path.join(a.stage, "local_chk")
        self.chk = os.path.join(a.stage, "chk")
        # per kind: local dir, symlink, destination
        self.kinds = {"snap": (self.local, self.current, a.dest), "chk": (self.local_chk, self.chk, a.chk_dest)}
        self.pending = {"snap": [], "chk": []}  # finished items (names), oldest first
        self.finishing = False  # pkdgrav3 has ended: workers drain what is queued, then return
        self.lock = threading.Lock()
        self.seen = set()
        self.log_pos = 0
        self.last_writing = None
        self.stop = False
        self.attempts = {}  # dst -> copy attempts so far (each writes a new temporary name)
        self.abandoned = []  # copy processes given up on (stalled)
        self.rep = open(a.report, "a", buffering=1) if a.report else sys.stdout
        for local, link, dest in self.kinds.values():
            os.makedirs(local, exist_ok=True)
            os.makedirs(dest, exist_ok=True)
            if not os.path.islink(link):
                self.point(link, local)

    def say(self, *msg):
        with self.lock:
            print(now(), *msg, file=self.rep, flush=True)

    # ---- staging target switch -------------------------------------------------------
    def point(self, link, target):
        tmp = link + ".new"
        if os.path.lexists(tmp):
            os.unlink(tmp)
        os.symlink(os.path.abspath(target), tmp)
        os.replace(tmp, link)  # atomic swap
        self.say(os.path.basename(link), "->", target)

    def check_space(self):
        free_gb = shutil.disk_usage(self.local).free / 2**30
        for kind, (local, link, dest) in self.kinds.items():
            to_nfs = os.path.realpath(link) == os.path.realpath(dest)
            if not to_nfs and free_gb < self.a.min_free_gb:
                self.say(f"local free {free_gb:.0f} GB < {self.a.min_free_gb} GB: new {kind} files go to {dest}")
                self.point(link, dest)
            elif to_nfs and free_gb > 2 * self.a.min_free_gb:
                self.point(link, local)

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
            elif CHK_WRITTEN.search(line):
                self.scan_checkpoints()

    def scan_checkpoints(self):
        """Complete checkpoints in local_chk: their restart file (written last) exists."""
        for n in sorted(os.listdir(self.local_chk)):
            if CHK_NAME.match(n):
                self.add(n, kind="chk")

    def scan_markers(self):
        for d in (self.local,):
            for n in sorted(os.listdir(d)):
                if n.endswith(".done"):
                    self.add(n[:-5], marker=os.path.join(d, n))

    def add(self, path, marker=None, kind="snap"):
        path = os.path.basename(path)  # located later in local/ or dest/ (symlink may have flipped)
        with self.lock:
            if (kind, path) in self.seen:
                # Seen before (e.g. from the log); a marker arriving later must not be left behind.
                if marker:
                    q = self.pending[kind]
                    for k, (p, m) in enumerate(q):
                        if p == path:
                            q[k] = (p, marker)
                            break
                    else:  # already drained
                        if os.path.exists(marker):
                            os.unlink(marker)
                return
            self.seen.add((kind, path))
            self.pending[kind].append((path, marker))
        self.say("finished", kind, path)

    # ---- copy ---------------------------------------------------------------------------
    def pieces(self, path):
        d, base = os.path.split(path)
        if not os.path.isdir(d):
            return []
        names = [n for n in os.listdir(d) if (n == base or n.startswith(base + ".")) and not n.endswith(".done")
                 and not re.search(r"\.part\d*$", n)]
        # the base file last: for a checkpoint it is the restart file, whose presence means "complete"
        return [os.path.join(d, n) for n in sorted(names, key=lambda n: (n == base, n))]

    def copy_one(self, src, dst):
        """One copy attempt of src -> dst in a child process (COPY_CHILD), watched through its progress pipe. A copy
        with no progress for --stall-s is abandoned (its process may be stuck in the kernel on NFS and cannot be
        killed) and fails with OSError; the worker retries, and every attempt writes a new temporary name
        (<dst>.part, .part2, ...): reopening a file whose writes hang blocks too (pawna, 20-22 Nov 2026)."""
        n = self.attempts.get(dst, 0) + 1
        self.attempts[dst] = n
        part = dst + (".part" if n == 1 else ".part%d" % n)
        size = os.path.getsize(src)
        t0 = time.time()
        p = subprocess.Popen([sys.executable, "-c", COPY_CHILD, src, part, str(CHUNK), str(self.a.bwlimit)],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL)
        done, t_last, buf = 0, time.time(), b""
        while True:
            r, _, _ = select.select([p.stdout], [], [], 5.0)
            if r:
                chunk = os.read(p.stdout.fileno(), 4096)
                if not chunk:          # child closed its output: finished (or died)
                    break
                buf += chunk
                *lines, buf = buf.split(b"\n")
                if lines:
                    done = int(lines[-1] or done)
                    t_last = time.time()
            elif time.time() - t_last > self.a.stall_s:
                try:
                    p.kill()
                except OSError:
                    pass
                self.abandoned.append(p)
                self.say(f"STALLED copy {src} -> {part}: no progress for {self.a.stall_s:.0f} s after "
                         f"{done/2**30:.2f} of {size/2**30:.2f} GiB; abandoned (pid {p.pid}), next attempt as a new file")
                raise IOError(f"copy stalled: {part}")
        try:
            rc = p.wait(timeout=max(60.0, self.a.stall_s))
        except subprocess.TimeoutExpired:
            self.abandoned.append(p)
            self.say(f"STALLED copy {src} -> {part}: output closed but the process does not exit; abandoned")
            raise IOError(f"copy stalled at exit: {part}")
        if rc != 0:
            raise IOError(f"copy {src} -> {part} failed (rc {rc}): {p.stderr.read().decode(errors='replace')[-300:]}")
        if os.path.getsize(part) != size:
            raise IOError(f"size mismatch {part}: {os.path.getsize(part)} != {size}")
        os.replace(part, dst)
        dfd = os.open(os.path.dirname(dst), os.O_RDONLY)
        try:
            os.fsync(dfd)
        finally:
            os.close(dfd)
        self.attempts.pop(dst, None)
        dt = time.time() - t0
        self.say(f"copied {src} -> {dst} {size/2**30:.2f} GiB {dt:.1f} s {size/2**20/max(dt,1e-9):.0f} MB/s"
                 + (f" (attempt {n})" if n > 1 else ""))
        return size, dt

    def drain_one(self, name, marker, kind="snap"):
        local, _, dest = self.kinds[kind]
        files = self.pieces(os.path.join(local, name))
        if not files:
            if self.pieces(os.path.join(dest, name)):
                self.say("already on NFS", name)  # written while the link pointed to dest
            else:
                self.say("WARNING nothing found for", name)
            if kind == "chk":
                self.prune_checkpoints()
            return True
        total = 0
        t0 = time.time()
        for src in files:
            dst = os.path.join(dest, os.path.basename(src))
            size, _ = self.copy_one(src, dst)
            total += size
        for src in files:      # delete only after every piece is safely on NFS
            os.unlink(src)
        if marker and os.path.exists(marker):
            os.unlink(marker)
        dt = time.time() - t0
        self.say(f"drained {name} {len(files)} file(s) {total/2**30:.2f} GiB "
                 f"{dt:.1f} s {total/2**20/max(dt,1e-9):.0f} MB/s")
        if kind == "chk":
            self.prune_checkpoints()
        return True

    def prune_checkpoints(self):
        """Keep the newest --keep-chk complete checkpoints in chk_dest; delete older ones (restart file first)."""
        K = self.a.keep_chk
        if K <= 0:
            return
        d = self.a.chk_dest
        by_name = {}
        for n in os.listdir(d):
            m = CHK_NAME.match(n)
            if m and os.path.isfile(os.path.join(d, n)):
                by_name.setdefault(m.group("name"), []).append((int(m.group("step")), n))
        for lst in by_name.values():
            lst.sort()
            if len(lst) <= K:
                continue
            for _, base in lst[:-K]:
                gone = [os.path.join(d, base)] + [p for p in self.pieces(os.path.join(d, base)) if os.path.basename(p) != base]
                total = 0
                for p in gone:
                    total += os.path.getsize(p)
                    os.unlink(p)
                self.say(f"deleted old checkpoint {os.path.join(d, base)} ({len(gone)} files, {total/2**30:.2f} GiB); "
                         f"kept {', '.join(n for _, n in lst[-K:])}")

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

    def worker(self, kind):
        """Drain one kind (snapshots or checkpoints) in order; they may go to different file servers."""
        failures = 0
        while True:
            with self.lock:
                item = self.pending[kind][0] if self.pending[kind] else None
                if item is None and self.finishing:
                    return
            if item is None:
                time.sleep(1)
                continue
            path, marker = item
            try:
                self.drain_one(path, marker, kind)
                with self.lock:
                    self.pending[kind].pop(0)
                failures = 0
            except OSError as e:
                failures += 1
                if self.finishing and failures >= self.a.retries:
                    self.say("ERROR", kind, path, e, f"- giving up after {failures} tries; left on local disk")
                    with self.lock:
                        self.pending[kind].pop(0)
                    failures = 0
                    continue
                self.say("ERROR", kind, path, e, f"- retry in {self.a.retry_wait:g} s")
                time.sleep(self.a.retry_wait)

    def run(self):
        self.say(f"start stage={self.a.stage} dest={self.a.dest} chk_dest={self.a.chk_dest} keep_chk={self.a.keep_chk} "
                 f"log={self.a.log} pid={self.a.pid} bwlimit={self.a.bwlimit} MB/s min_free={self.a.min_free_gb} GB")
        workers = [threading.Thread(target=self.worker, args=(k,), daemon=True) for k in self.kinds]
        for w in workers:
            w.start()
        while True:
            running = self.alive()
            self.scan_log()
            self.scan_markers()
            self.check_space()
            if not running:
                break
            time.sleep(self.a.poll)
        self.scan_log()
        self.scan_markers()
        self.scan_checkpoints()
        with self.lock:
            self.finishing = True
        for w in workers:
            w.join()
        left = [os.path.join(d, n) for d in (self.local, self.local_chk) for n in os.listdir(d)]
        if left:
            self.say("LEFT ON LOCAL DISK (incomplete or unrecognised):", " ".join(left))
            return 2
        self.say("all drained")
        return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", required=True, help="node-local staging root (creates local/ and current)")
    p.add_argument("--dest", required=True, help="final directory on shared storage")
    p.add_argument("--chk-dest", help="final directory for checkpoints (default: --dest); may be another file server")
    p.add_argument("--keep-chk", type=int, default=0,
                   help="keep only the newest K complete checkpoints in --chk-dest (0 = keep all)")
    p.add_argument("--log", help="pkdgrav3 stdout log to watch for finished writes")
    p.add_argument("--pid", type=int, help="pkdgrav3 process to wait for")
    p.add_argument("--bwlimit", type=float, default=0, help="copy rate limit in MB/s (0 = none)")
    p.add_argument("--min-free-gb", type=float, default=200, help="switch new writes to NFS below this")
    p.add_argument("--poll", type=float, default=10, help="seconds between checks")
    p.add_argument("--retries", type=int, default=10, help="after pkdgrav3 has ended: copy attempts per item")
    p.add_argument("--retry-wait", type=float, default=60, help="seconds between copy attempts")
    p.add_argument("--stall-s", type=float, default=900,
                   help="a copy with no progress for this long is abandoned and tried again as a new file")
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
            if not n.endswith(".done") and not re.search(r"\.part\d*$", n):
                d.add(n)
        # complete checkpoints as one unit each; anything else (an interrupted checkpoint) file by file
        names = sorted(os.listdir(d.local_chk))
        bases = [n for n in names if CHK_NAME.match(n)]
        for n in names:
            if re.search(r"\.part\d*$", n) or (n not in bases and any(n.startswith(b + ".") for b in bases)):
                continue
            d.add(n, kind="chk")
        d.stop = True  # no pid to wait for: drain once and exit
    signal.signal(signal.SIGTERM, lambda *_: setattr(d, "stop", True))
    return d.run()


if __name__ == "__main__":
    sys.exit(main())
