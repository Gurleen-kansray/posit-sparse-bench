#!/bin/bash
# What sacct and the result files say about a run: which tasks died, how close the
# survivors came to their limit, and whether any output is malformed.
#
#   bash hpc/audit.sh                  # everything since Oct 29
#   bash hpc/audit.sh 2026-11-02       # since a given date
#   MEM=1 bash hpc/audit.sh            # also peak MaxRSS per array (slow)
#
# Run it WITHOUT redirecting stderr, so the section markers reach your terminal
# and you can see where it is:
#
#   bash hpc/audit.sh > audit.txt      # markers on screen, report in the file
#
# Every external command here is wrapped in a timeout. Nothing in this script can
# hang: a section that exceeds its budget says so and the script moves on.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
source "$REPO/hpc/config.sh"

SINCE="${1:-2026-10-29}"
SACCT_T="${SACCT_T:-120}"        # seconds allowed per sacct call
SCAN_T="${SCAN_T:-180}"          # seconds allowed per filesystem scan

RAW=$(mktemp)
trap 'rm -f "$RAW"' EXIT

mark() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$1" >&2; }

# -X is what keeps this fast: one row per job rather than one per step. Asking for
# steps across 5,000 array tasks is tens of thousands of rows and minutes of
# slurmdbd time, which is why MaxRSS is behind MEM=1.
mark "querying sacct (one row per task)"
if ! timeout "$SACCT_T" sacct -u "$USER" -S "$SINCE" -X -P -n \
        -o JobID,JobName,State,Elapsed,Timelimit,ReqMem,ExitCode > "$RAW"; then
    echo "sacct did not return within ${SACCT_T}s. Raise it with SACCT_T=300,"
    echo "or narrow the window with: bash hpc/audit.sh <YYYY-MM-DD>"
    exit 1
fi
mark "got $(wc -l < "$RAW") rows"

echo "=== state by array ==="
awk -F'|' '$2 ~ /^pcg_/ {
        split($3, s, " "); key = $2 "|" s[1]; c[key]++
    }
    END { for (k in c) { split(k, p, "|"); printf "%-20s %-12s %5d\n", p[1], p[2], c[k] } }' \
    "$RAW" | sort

n=$(awk -F'|' '$3 ~ /TIMEOUT/' "$RAW" | wc -l)
echo
echo "=== tasks that hit the wall limit ($n total, first 40) ==="
awk -F'|' '$3 ~ /TIMEOUT/ {printf "%-22s %-20s %10s / %s\n", $1, $2, $4, $5}' "$RAW" \
    | head -40

echo
echo "=== tasks that failed for another reason ==="
awk -F'|' '$3 ~ /FAILED|CANCELLED|NODE_FAIL|OUT_OF_ME/ {
        printf "%-22s %-20s %-14s exit=%s\n", $1, $2, $3, $7}' "$RAW" | head -40

echo
echo "=== which nodes the kills landed on ==="
# A node holding far more than its share is a hardware fault to report, not a
# cost model to refit. On the first pass node-153 was throttling its CPUs 8x and
# accounted for 104 of 135 tasks that ran there.
if timeout "$SACCT_T" sacct -u "$USER" -S "$SINCE" -X -P -n \
        -o JobName,State,NodeList > "$RAW.nodes" 2>/dev/null; then
    awk -F'|' '$2 ~ /TIMEOUT/ {print $3}' "$RAW.nodes" \
        | sort | uniq -c | sort -rn | head -10
    rm -f "$RAW.nodes"
else
    echo "(sacct timed out; skipped)"
fi

echo
echo "=== closest to the wall limit, among completed ==="
awk -F'|' '
    function sec(t) {
        d = 0
        if (t ~ /-/) { split(t, a, "-"); d = a[1]; t = a[2] }
        split(t, b, ":")
        return d*86400 + b[1]*3600 + b[2]*60 + b[3]
    }
    $3 == "COMPLETED" && $5 ~ /:/ {
        e = sec($4); l = sec($5)
        if (l > 0) printf "%.3f\t%-20s %10s / %-10s\n", e/l, $2, $4, $5
    }' "$RAW" | sort -rn | head -15 | cut -f2-

if [ -n "${MEM:-}" ]; then
    echo
    echo "=== memory: requested vs peak actually used ==="
    mark "querying sacct for MaxRSS -- this needs every job STEP, so it is slow"
    if timeout "$SACCT_T" sacct -u "$USER" -S "$SINCE" -P -n \
            -o JobID,JobName,MaxRSS,ReqMem > "$RAW.mem" 2>/dev/null; then
        awk -F'|' '
            function kb(v) {
                x = v + 0
                if (v ~ /K$/) return x
                if (v ~ /M$/) return x * 1024
                if (v ~ /G$/) return x * 1024 * 1024
                return x / 1024
            }
            function reqmb(v) {
                x = v + 0
                if (v ~ /[Gg]/) return x * 1024
                return x
            }
            $2 != "batch" && $2 ~ /^pcg_/ {
                name[$1] = substr($2, 5)
                if ($4 != "") req[substr($2, 5)] = reqmb($4)
                next
            }
            $2 == "batch" {
                id = $1; sub(/\.batch$/, "", id)
                a = name[id]
                if (a == "" || $3 == "") next
                k = kb($3)
                if (k > peak[a]) peak[a] = k
                n[a]++
            }
            END {
                printf "%-20s %6s %12s %12s %8s\n", \
                       "array", "steps", "peak MaxRSS", "requested", "used"
                for (a in peak)
                    printf "%-20s %6d %9.1f MB %9.0f MB %7.2f%%\n", \
                           a, n[a], peak[a]/1024, req[a],
                           (req[a] > 0 ? 100 * (peak[a]/1024) / req[a] : 0)
            }' "$RAW.mem" | sort
        rm -f "$RAW.mem"
    else
        echo "(sacct timed out after ${SACCT_T}s; raise with SACCT_T=600)"
    fi
fi

DIRS="results/seed_sweep_50 results/ladder_50 results/qmv_50/qmv
      results/qmv_50/plain results/refine11_50 results/refine12_quire_50
      results/refine12_plain_50 results/refine_exactb_50"

echo
echo "=== results on disk vs .tmp left behind ==="
mark "counting result files"
for d in $DIRS; do
    [ -d "$d" ] || continue
    ok=$(timeout "$SCAN_T" find "$d" -name '*.log' -size +0 | wc -l)
    empty=$(timeout "$SCAN_T" find "$d" -name '*.log' -empty | wc -l)
    tmp=$(timeout "$SCAN_T" find "$d" -name '*.log.tmp' | wc -l)
    printf '%-28s %5d complete  %3d EMPTY  %3d .tmp\n' "$d" "$ok" "$empty" "$tmp"
done

echo
echo "=== log integrity: header width and body length ==="
# A log with a short header, or with no data rows, is worse than a missing one:
# the aggregation scripts index by column position.
#
# This reads three lines per file and stops, in one process. Reading the files
# whole is about 10 million lines over a shared filesystem, and an awk or shell
# per file is ~10,000 process launches; both took long enough to look like a
# hang. A timeout on `xargs` does not help either -- it kills xargs and leaves
# the awk it spawned reading, so the pipeline never ends and Ctrl-C cannot
# reach it.
mark "scanning log headers (3 lines per file)"
timeout -k 5 "$SCAN_T" python3 -I - $DIRS <<'PY' || \
    echo "(scan exceeded ${SCAN_T}s; skipped)"
import os, sys
from collections import Counter

for d in sys.argv[1:]:
    if not os.path.isdir(d):
        continue
    cols, short, n = Counter(), 0, 0
    for root, _, files in os.walk(d):
        for fn in files:
            if not fn.endswith(".log"):
                continue
            p = os.path.join(root, fn)
            try:
                if os.path.getsize(p) == 0:
                    continue
                with open(p, errors="replace") as f:
                    head = [next(f, None) for _ in range(3)]
            except OSError:
                continue
            n += 1
            if head[1] is not None:
                cols[len(head[1].split())] += 1
            if head[2] is None:
                short += 1
    width = " ".join("%dx%d" % (c, w) for w, c in sorted(cols.items()))
    print("%-28s header cols: %-22s %d of %d with no data rows"
          % (d, width or "none", short, n))
PY

mark "done"
