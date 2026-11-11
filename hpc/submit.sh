#!/bin/bash
# Submit the experiment set as one SLURM array per (program, size class), reading
# hpc/jobs/manifest.tsv. Regenerate that with:
#     python3 hpc/make_joblists.py
#
#   bash hpc/submit.sh                 # everything, longest arrays first
#   bash hpc/submit.sh sweep ladder    # only classes whose name starts with these
#   DRYRUN=1 bash hpc/submit.sh        # print the sbatch commands, submit nothing
#
# The manifest is ordered longest-task-first. A task cannot be split, so the longest
# one floors the finish time; starting those first lets the cheap tasks backfill
# around them. Per-class --time requests are what make that backfill possible.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
source "$REPO/hpc/config.sh"

MANIFEST=hpc/jobs/manifest.tsv
GEN="python3 hpc/make_joblists.py"
ARGS=()
for a in "$@"; do
    if [ "$a" = "--requeue" ]; then
        MANIFEST=hpc/jobs/manifest.requeue.tsv
        GEN="python3 hpc/requeue.py --write"
    else
        ARGS+=("$a")
    fi
done

if [ ! -s "$MANIFEST" ]; then
    echo "ERROR: $MANIFEST is missing. Generate it with:"
    echo "    $GEN"
    exit 1
fi

# The account and partition come from hpc/config.sh, which honors anything
# already exported. Most clusters reject jobs without an account.
ACCOUNT="$HPC_ACCOUNT"
if [ -z "$ACCOUNT" ]; then
    echo "ERROR: no account set. The scheduler rejects jobs without one."
    echo "  Find yours with:   sacctmgr show assoc where user=\$USER format=account"
    echo "  Then add to ~/.bashrc:"
    echo "      export SBATCH_ACCOUNT=<account>"
    echo "  and re-login, or: source ~/.bashrc"
    exit 1
fi

WANT=("${ARGS[@]+"${ARGS[@]}"}")

# One array per class is fine; two sets of arrays over the same job lists are not,
# because both would write the same <log>.tmp. Refuse rather than corrupt a seed.
running=$(squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -c '^pcg_' || true)
if [ "${running:-0}" -gt 0 ] && [ -z "${FORCE:-}" ]; then
    echo "ERROR: $running pcg_* job(s) are still queued or running."
    echo "  Two arrays over the same job list share one .tmp per seed and will"
    echo "  clobber each other. Wait for them, or scancel, or FORCE=1 to override."
    squeue -u "$USER" -o '%.18i %.20j %.8T %.10M %R' | head
    exit 1
fi

# Slurm will not create these. A missing output directory fails the job before it
# runs, with JobLaunchFailure and no log to explain it.
mkdir -p "$SLURM_OUT"
awk -F'\t' '!/^#/ {print $3}' "$MANIFEST" | sort -u | while read -r d; do
    mkdir -p "$d"
done

submitted=0
skipped=0
while IFS=$'\t' read -r cls binary outdir joblist ntasks time mem cap; do
    case "$cls" in \#*|"") continue ;; esac

    if [ ${#WANT[@]} -gt 0 ]; then
        match=0
        for w in "${WANT[@]}"; do
            case "$cls" in "$w"*) match=1 ;; esac
        done
        [ $match -eq 1 ] || { skipped=$((skipped+1)); continue; }
    fi

    if [ ! -x "$binary" ]; then
        echo "SKIP $cls: $binary is not built. Run: bash hpc/build.sh"
        skipped=$((skipped+1))
        continue
    fi

    cmd=(sbatch
         --account="$ACCOUNT"
         --partition="$PARTITION"
         --job-name="pcg_${cls}"
         --time="$time"
         --mem="$mem"
         --array="1-${ntasks}%${cap}"
         --output="$SLURM_OUT/%x_%A_%a.out"
         --error="$SLURM_OUT/%x_%A_%a.out"
         --export=ALL,JOBLIST="$joblist",BIN="$binary",OUTDIR="$outdir"
         hpc/sweep.sbatch)

    printf '%-18s %5s tasks  %10s  %5s  cap %-3s -> %s\n' \
        "$cls" "$ntasks" "$time" "$mem" "$cap" "$outdir"
    if [ -n "${DRYRUN:-}" ]; then
        printf '    '; printf '%q ' "${cmd[@]}"; echo
    else
        "${cmd[@]}"
    fi
    submitted=$((submitted+1))
done < "$MANIFEST"

echo
echo "submitted $submitted array(s), skipped $skipped"
cat <<'EOF'

Monitor:
  squeue -u $USER
  bash hpc/progress.sh
Cancel one program's arrays:
  scancel --name=pcg_sweep_xl
EOF
