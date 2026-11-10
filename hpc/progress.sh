#!/bin/bash
# Completion per program, counted against the job lists rather than against a
# hardcoded 50, so it stays correct when the manifest changes.
#
#   bash hpc/progress.sh            # one line per array class
#   bash hpc/progress.sh -v         # also per matrix, for incomplete classes
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
source "$REPO/hpc/config.sh"
MANIFEST=hpc/jobs/manifest.tsv
[ -s "$MANIFEST" ] || { echo "no $MANIFEST; run python3 hpc/make_joblists.py"; exit 1; }
VERBOSE="${1:-}"

# A job list line is "matrix seed [maxiter]". The runner names the log
# <matrix>_seed<seed>.log, or <matrix>_seed<seed>_it<maxiter>.log when a
# non-default budget is given, so the expected filename is derivable here.
logname() {
    if [ -n "${3:-}" ]; then echo "${1}_seed${2}_it${3}.log"
    else echo "${1}_seed${2}.log"; fi
}

total_done=0; total_want=0
printf '%-18s %11s  %s\n' class done bar
while IFS=$'\t' read -r cls binary outdir joblist ntasks time mem cap; do
    case "$cls" in \#*|"") continue ;; esac
    have=0
    while read -r m s extra; do
        [ -n "${m:-}" ] || continue
        [ -s "$outdir/$(logname "$m" "$s" "${extra:-}")" ] && have=$((have+1))
    done < "$joblist"
    total_done=$((total_done+have)); total_want=$((total_want+ntasks))
    pct=$(( ntasks > 0 ? have*100/ntasks : 0 ))
    filled=$((pct/4)); empty=$((25-filled))
    bar=""
    [ "$filled" -gt 0 ] && bar=$(printf '#%.0s' $(seq 1 "$filled"))
    [ "$empty" -gt 0 ] && bar="$bar$(printf '.%.0s' $(seq 1 "$empty"))"
    printf '%-18s %5d/%-5d [%s] %3d%%\n' "$cls" "$have" "$ntasks" "$bar" "$pct"

    if [ "$VERBOSE" = "-v" ] && [ "$have" -lt "$ntasks" ]; then
        awk '{print $1}' "$joblist" | sort -u | while read -r m; do
            want=$(awk -v m="$m" '$1==m' "$joblist" | wc -l)
            got=$(ls "$outdir/${m}_seed"*.log 2>/dev/null | wc -l)
            [ "$got" -lt "$want" ] && printf '    %-12s %d/%d\n' "$m" "$got" "$want"
        done
    fi
done < "$MANIFEST"

printf '\n%d/%d logs (%d%%)\n' "$total_done" "$total_want" \
    $(( total_want > 0 ? total_done*100/total_want : 0 ))

# A .tmp means a task died mid-write; that seed will be redone on resubmission.
tmp=$(awk -F'\t' '!/^#/ {print $3}' "$MANIFEST" | sort -u \
      | xargs -I{} sh -c 'ls {}/*.tmp 2>/dev/null' | wc -l)
[ "$tmp" -gt 0 ] && echo "$tmp partial .tmp file(s) -- those seeds are incomplete"

# Scanning every task log for FAIL costs one read per file, and after two passes
# $SLURM_OUT holds thousands of them on a shared filesystem that the running
# tasks are also writing to. Only the recent ones can be new, so bound the scan;
# -a forces the full sweep.
nout=$(find "$SLURM_OUT" -name '*.out' 2>/dev/null | wc -l)
if [ "${1:-}" = "-a" ] || [ "$nout" -le 600 ]; then
    scan=$(find "$SLURM_OUT" -name '*.out' 2>/dev/null)
else
    scan=$(find "$SLURM_OUT" -name '*.out' -mmin -1440 2>/dev/null)
    echo "(scanning the last day's task logs only, of $nout; -a for all)"
fi
if [ -n "$scan" ]; then
    hits=$(printf '%s\n' "$scan" | timeout 60 xargs -r grep -h "^FAIL" 2>/dev/null \
           | sort | uniq -c | sort -rn)
    if [ -n "$hits" ]; then
        echo "task log(s) report FAIL:"
        printf '%s\n' "$hits" | head
    fi
fi

echo
echo "Gaps are filled by resubmitting: completed logs are skipped."
echo "  bash hpc/submit.sh              # everything still missing"
echo "  bash hpc/submit.sh sweep        # just one program's classes"
