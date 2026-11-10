#!/usr/bin/env python3
"""Build a manifest covering only the seeds that are still missing.

    python3 hpc/requeue.py              # look, change nothing
    python3 hpc/requeue.py --write      # write hpc/jobs/manifest.requeue.tsv
    bash hpc/submit.sh --requeue        # submit that manifest

Why not just rerun submit.sh. Resubmitting a full array relaunches every task,
and although a finished seed exits in under a second on the skip-if-exists test,
5,000 of those still queue, start and bill. This writes job lists holding only
the gaps, so the arrays are the size of the work left.

It also does three things the first manifest got wrong:

  * --time is scaled by how badly each class actually did, read from sacct, not
    by one multiplier across the board. A class that lost a few tasks off the
    tail of its distribution gets 2x; a class where the limit cut through the
    middle gets 4x, because there the cost model is wrong rather than tight.
    A TIMEOUT loses the seed, while an over-request costs nothing: billing is on
    time used, not time asked for.
  * --mem is 2 GB, which is what ARCSIM asked for after measuring a 207 MiB
    peak across the first pass. The original requests were 4, 8 and 16 G by
    matrix size, and those requests, not cores, are what stopped tasks from
    starting: an epyc node has 503 GB, so 16 G a task caps it at 31 tasks
    however many cores sit idle.
  * the %N concurrency caps are shared out of one global budget (CAP_TOTAL,
    default 128) instead of being set per array, since %N is per array and the
    first pass could have asked for 96+64+96 at once.
"""
import os
import subprocess
import sys
from math import ceil

CAP_TOTAL = int(os.environ.get("CAP_TOTAL", 128))
MIN_WALL = 3600                      # seconds; nothing is asked for less
SINCE = os.environ.get("SINCE", "2026-10-29")

# Wall-time multiplier as a function of the class's observed timeout rate. Below
# the threshold the model is roughly right and only the tail overran; above it the
# limit fell inside the bulk of the distribution and the model itself is wrong.
#
# TAIL_MULT is 1.25 rather than 2.0 because ARCSIM traced the first pass's
# scattered timeouts --- all 59 outside the refinement _lg arrays --- to node-153,
# whose management controller had a bad temperature setting and was throttling the
# CPUs about 8x. 104 of the 135 tasks that landed there were killed. The node is
# fixed and back in service, so those classes' original limits were adequate and
# need margin, not doubling. Only the refinement _lg timeouts were real: 274 of
# their 319 kills were on healthy nodes.
TAIL_MULT, BULK_MULT, BULK_AT = 1.25, 4.0, 0.20

# TIME_MULT overrides the whole thing, for when sacct is unavailable or you want
# one flat factor.
TIME_MULT = os.environ.get("TIME_MULT")
TIME_MULT = float(TIME_MULT) if TIME_MULT else None

MANIFEST = "hpc/jobs/manifest.tsv"
OUT_MANIFEST = "hpc/jobs/manifest.requeue.tsv"

# Memory. Two independent measurements, and the larger one sets the request:
#
#   getrusage on the binary alone, i7-12800H, all eight programs: 21 bytes per
#   matrix nonzero plus ~8 MB, so 101 MB on s3dkq4m2 (4,820,891 nnz), 55 MB on
#   nasasrb, 28 MB on bcsstk37.
#
#   Slurm MaxRSS on the first pass, reported by ARCSIM: 207 MiB peak, refinement
#   arrays included. Higher because the cgroup is charged for the whole job step
#   -- bash, sed and the binary -- not the binary on its own.
#
# ARCSIM asked for 2 GB a task on that basis, which is what MEM_FLOOR_MB sets.
# Scaling the smaller matrices below it would buy nothing: an epyc node has
# 503 GB, so 128 tasks at 2 GB is half of one node's memory, and the first pass
# showed the 16 G requests were what kept tasks from starting.
MEM_PER_NNZ = 21
MEM_FIXED_MB = 8
MEM_SAFETY = 6
MEM_FLOOR_MB = int(os.environ.get("MEM_FLOOR_MB", 2048))


def nnz_of(matrix):
    """Nonzero count from the Matrix Market banner, so no table can drift."""
    path = "data/matrices/%s.mtx" % matrix
    with open(path) as f:
        for line in f:
            if line.startswith("%"):
                continue
            return int(line.split()[2])
    raise ValueError("no size line in " + path)


def hhmmss(seconds):
    s = int(ceil(max(seconds, MIN_WALL) / 600.0) * 600)
    return "%02d:%02d:00" % (s // 3600, (s % 3600) // 60)


def secs(t):
    """[dd-]hh:mm:ss -> seconds. sacct uses the day form past 24 hours."""
    d = 0
    if "-" in t:
        d, t = t.split("-", 1)
    try:
        h, m, s = (int(x) for x in t.split(":"))
    except ValueError:
        return 0
    return int(d) * 86400 + h * 3600 + m * 60 + s


def timeout_rates():
    """Per-class timeout fraction and longest completed task, from sacct.

    Returns {class: (fraction_timed_out, max_completed_seconds)}. An empty dict
    when sacct is missing or says nothing, in which case every class falls back to
    TAIL_MULT.
    """
    try:
        out = subprocess.run(
            ["sacct", "-u", os.environ.get("USER", ""), "-S", SINCE, "-X", "-P",
             "-n", "-o", "JobName,State,Elapsed"],
            capture_output=True, text=True, timeout=120).stdout
    except (OSError, subprocess.SubprocessError):
        return {}

    tally = {}
    for line in out.splitlines():
        f = line.split("|")
        if len(f) < 3 or not f[0].startswith("pcg_"):
            continue
        cls, state = f[0][4:], f[1].split()[0]
        t, done, longest = tally.setdefault(cls, [0, 0, 0])
        tally[cls][0] = t + (state == "TIMEOUT")
        if state == "COMPLETED":
            tally[cls][1] = done + 1
            tally[cls][2] = max(longest, secs(f[2]))
    return {c: (v[0] / float(v[0] + v[1]) if v[0] + v[1] else 0.0, v[2])
            for c, v in tally.items()}


def logname(matrix, seed, extra):
    if extra:
        return "%s_seed%s_it%s.log" % (matrix, seed, extra)
    return "%s_seed%s.log" % (matrix, seed)


write = "--write" in sys.argv
rows, nnz_cache = [], {}
rates = {} if TIME_MULT else timeout_rates()

with open(MANIFEST) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        cls, binary, outdir, joblist, ntasks, time, mem, cap = \
            line.rstrip("\n").split("\t")

        missing, peak_nnz = [], 0
        with open(joblist) as jl:
            for jline in jl:
                parts = jline.split()
                if not parts:
                    continue
                matrix, seed = parts[0], parts[1]
                extra = parts[2] if len(parts) > 2 else ""
                out = os.path.join(outdir, logname(matrix, seed, extra))
                if os.path.exists(out) and os.path.getsize(out) > 0:
                    continue
                missing.append(jline.rstrip("\n"))
                if matrix not in nnz_cache:
                    nnz_cache[matrix] = nnz_of(matrix)
                peak_nnz = max(peak_nnz, nnz_cache[matrix])

        if not missing:
            continue

        mb = (MEM_PER_NNZ * peak_nnz / 1e6 + MEM_FIXED_MB) * MEM_SAFETY
        newmem = "%dM" % max(MEM_FLOOR_MB, int(ceil(mb / 256.0) * 256))

        # The new limit is built on the old one, not on the longest task that
        # survived: in a class the limit cut through, the survivors all stopped
        # just short of it, so their elapsed times say nothing about how long the
        # killed ones needed.
        if TIME_MULT:
            mult, why = TIME_MULT, "TIME_MULT"
        else:
            frac, _ = rates.get(cls, (0.0, 0))
            if frac >= BULK_AT:
                mult, why = BULK_MULT, "%.0f%% killed" % (100 * frac)
            else:
                mult, why = TAIL_MULT, "%.0f%% killed" % (100 * frac)

        rows.append(dict(cls=cls, binary=binary, outdir=outdir,
                         joblist="hpc/jobs/rq_%s.txt" % cls,
                         lines=missing, ntasks=len(missing),
                         time=hhmmss(secs(time) * mult), why=why,
                         oldtime=time, mem=newmem, oldmem=mem))

if not rows:
    print("nothing missing; every seed in the manifest has a non-empty log.")
    sys.exit(0)

# Share the concurrency budget to minimize the finish time of the slowest array.
#
# Splitting it in proportion to task count starves a small array holding long
# tasks: 9 sweep_xl seeds at a 60 h limit would draw 3 of 128 slots and take three
# waves, 180 h, while nothing else needed them. An array's finish is
# ceil(ntasks/cap) waves of its wall limit, so each spare slot goes to whichever
# array is currently projected to finish last. That is the quantity that matters
# -- the run is done when the last array is.
total = sum(r["ntasks"] for r in rows)
for r in rows:
    r["wall"] = secs(r["time"])


def caps_for(target):
    """Slots each array needs to finish by `target`, or None if impossible.

    An array of n tasks at cap c runs in ceil(n/c) waves of its wall limit, so
    the waves it can afford is floor(target/wall) and the cap it needs to fit n
    tasks into that many waves is ceil(n/waves).
    """
    need = []
    for r in rows:
        waves = int(target // r["wall"])
        if waves < 1:
            return None
        need.append(min(r["ntasks"], ceil(r["ntasks"] / float(waves))))
    return need


# Pick the earliest finish the budget can buy. Handing a slot to the worst array
# one at a time stalls, because a slot only helps when it drops the wave count by
# a whole step: 87 tasks at 22 slots and at 28 both take 4 waves. Candidate
# targets are the finish times actually reachable, w waves of some array's limit.
targets = sorted({w * r["wall"] for r in rows
                  for w in range(1, r["ntasks"] + 1)})
chosen = None
for t in targets:
    need = caps_for(t)
    if need is not None and sum(need) <= CAP_TOTAL:
        chosen = need
        break
if chosen is None:                   # budget cannot reach any target; one each
    chosen = [1] * len(rows)
for r, c in zip(rows, chosen):
    r["cap"] = c

# Spend whatever the chosen target left over, longest-finishing array first. It
# cannot improve the worst case -- that is what the target fixed -- but it pulls
# the other arrays in, which frees the cluster sooner.
spare = CAP_TOTAL - sum(r["cap"] for r in rows)
while spare > 0:
    cand = [r for r in rows if r["cap"] < r["ntasks"]
            and ceil(r["ntasks"] / float(r["cap"] + 1)) < ceil(r["ntasks"] / float(r["cap"]))]
    if not cand:
        break
    r = max(cand, key=lambda r: ceil(r["ntasks"] / float(r["cap"])) * r["wall"])
    step = ceil(r["ntasks"] / float(ceil(r["ntasks"] / float(r["cap"])) - 1)) - r["cap"]
    if step > spare:
        break
    r["cap"] += step
    spare -= step

print("%-18s %6s %22s %16s %5s  %s"
      % ("class", "left", "--time", "mem", "cap", "why"))
for r in rows:
    print("%-18s %6d %10s -> %-9s %6s -> %-6s %5d  %s"
          % (r["cls"], r["ntasks"], r["oldtime"], r["time"],
             r["oldmem"], r["mem"], r["cap"], r["why"]))
print("\n%d array(s), %d tasks left, %d concurrent at most"
      % (len(rows), total, sum(r["cap"] for r in rows)))

# Worst-case finish: each array's tasks run in ceil(ntasks/cap) waves of its
# wall limit. Tasks that finish early shorten this, and arrays overlap, so it is
# an upper bound rather than a forecast.
worst = max(ceil(r["ntasks"] / float(r["cap"])) * secs(r["time"]) for r in rows)
print("longest single array, worst case: %.1f h" % (worst / 3600.0))

if not write:
    print("\nNothing written. Re-run with --write, then: bash hpc/submit.sh --requeue")
    sys.exit(0)

stale = 0
for r in rows:
    with open(r["joblist"], "w") as f:
        f.write("\n".join(r["lines"]) + "\n")
    # A killed task leaves <log>.tmp behind. The program truncates it on reopen,
    # so it is harmless, but removing it keeps progress.sh's partial count honest.
    for line in r["lines"]:
        parts = line.split()
        tmp = os.path.join(r["outdir"], logname(
            parts[0], parts[1], parts[2] if len(parts) > 2 else "") + ".tmp")
        if os.path.exists(tmp):
            os.remove(tmp)
            stale += 1
if stale:
    print("removed %d stale .tmp file(s) for seeds being redone" % stale)
with open(OUT_MANIFEST, "w") as f:
    f.write("#class\tbinary\toutdir\tjoblist\tntasks\ttime\tmem\tcap\n")
    for r in rows:
        f.write("%(cls)s\t%(binary)s\t%(outdir)s\t%(joblist)s\t%(ntasks)d\t"
                "%(time)s\t%(mem)s\t%(cap)d\n" % r)
print("\nwrote %s and %d job list(s)" % (OUT_MANIFEST, len(rows)))
print("Submit with: bash hpc/submit.sh --requeue")
