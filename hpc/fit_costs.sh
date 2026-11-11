#!/bin/bash
# Per-matrix, per-program cost in seconds, from what the run actually did.
#
#   bash hpc/fit_costs.sh                 # all classes, since Oct 29
#   bash hpc/fit_costs.sh 2026-11-02      # one pass only
#
# Prints a median and a maximum per (class, matrix), plus a REFINE_SECONDS block
# ready to paste into make_joblists.py. Run it on the login node after a pass
# finishes; it only reads sacct and the job lists.
#
# Why this is not a one-liner: sacct knows the array task index, not the matrix.
# The mapping is the job list -- task 17 ran line 17 -- so each class's elapsed
# times have to be joined against the list it was submitted with. Both manifests
# are consulted, since a resubmitted class used hpc/jobs/rq_<class>.txt while
# the first pass used hpc/jobs/<class>.txt.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
source "$REPO/hpc/config.sh"

SINCE="${1:-2026-10-29}"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

sacct -u "$USER" -S "$SINCE" -X -P -n -o JobID,JobName,State,Elapsed \
    | awk -F'|' '$2 ~ /^pcg_/ && $3 == "COMPLETED"' > "$TMP/done.psv"

if [ ! -s "$TMP/done.psv" ]; then
    echo "no COMPLETED pcg_* tasks since $SINCE"
    exit 1
fi

python3 - "$TMP/done.psv" <<'PY'
import os, sys, glob
from collections import defaultdict


def secs(t):
    d = 0
    if "-" in t:
        d, t = t.split("-", 1)
    try:
        h, m, s = (int(x) for x in t.split(":"))
    except ValueError:
        return None
    return int(d) * 86400 + h * 3600 + m * 60 + s


# class -> [matrix per 1-based task index]. A resubmitted class has both lists;
# prefer rq_, since that is what the later array ran.
lists = {}
for path in sorted(glob.glob("hpc/jobs/*.txt")):
    name = os.path.basename(path)[:-4]
    cls = name[3:] if name.startswith("rq_") else name
    if name.startswith("rq_") or cls not in lists:
        lists[cls] = [l.split()[0] for l in open(path) if l.split()]

# A class that was resubmitted has two arrays with two different job lists, and
# a task index means different matrices in each. Group by array id so each is
# resolved against the list its own array used.
by_array = defaultdict(list)
for line in open(sys.argv[1]):
    jobid, name, _, elapsed = line.rstrip("\n").split("|")
    cls = name[4:]
    if "_" not in jobid:
        continue
    arr, idx = jobid.split("_", 1)
    if not idx.isdigit():
        continue
    e = secs(elapsed)
    if e is not None:
        by_array[(cls, arr)].append((int(idx), e))

# The first pass's arrays span the full job list; a requeue array is shorter. Use
# length to decide which list an array's indices address.
cost = defaultdict(list)
for (cls, arr), rows in by_array.items():
    full = lists.get(cls)
    rq = lists.get(cls) if not os.path.exists("hpc/jobs/rq_%s.txt" % cls) \
        else [l.split()[0] for l in open("hpc/jobs/rq_%s.txt" % cls) if l.split()]
    orig = [l.split()[0] for l in open("hpc/jobs/%s.txt" % cls)
            if l.split()] if os.path.exists("hpc/jobs/%s.txt" % cls) else None
    hi = max(i for i, _ in rows)
    table = None
    for cand in (orig, rq):
        if cand and hi <= len(cand):
            table = cand
            break
    if table is None:
        table = full
    if not table:
        continue
    for i, e in rows:
        if i <= len(table):
            cost[(cls, table[i - 1])].append(e)

print("%-18s %-11s %6s %10s %10s %10s" %
      ("class", "matrix", "seeds", "median", "max", "max/med"))
for (cls, m), v in sorted(cost.items()):
    v.sort()
    med, mx = v[len(v) // 2], v[-1]
    print("%-18s %-11s %6d %9ds %9ds %9.2f"
          % (cls, m, len(v), med, mx, mx / float(med) if med else 0))

# One number per matrix for the refinement programs, taken as the largest median
# across the four of them: they agree within a few percent, and the wall limit
# has to cover the slowest.
ref = defaultdict(list)
for (cls, m), v in cost.items():
    if cls.startswith(("refine11", "ref12quire", "ref12plain", "refexactb")):
        v2 = sorted(v)
        ref[m].append(v2[len(v2) // 2])
if ref:
    print("\n# Paste into make_joblists.py. Median seconds per seed on one core,")
    print("# worst of the four refinement programs, measured.")
    print("REFINE_SECONDS = {")
    for m in sorted(ref, key=lambda m: -max(ref[m])):
        print('    "%s": %6d,' % (m, max(ref[m])))
    print("}")
PY
