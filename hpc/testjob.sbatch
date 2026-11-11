#!/bin/bash
# Self-contained test job, for ARCSIM staff to check resource requests before the
# full run. Same shape as the production arrays but 4 tasks instead of 5,249 and
# about two seconds of compute each.
#
#   mkdir -p hpc/slurm results/testjob
#   sbatch --account=$SBATCH_ACCOUNT hpc/testjob.sbatch
#   cat hpc/slurm/pcgtest_*.out
#
# WHAT THE PRODUCTION RUN LOOKS LIKE
#   27 job arrays, 5,249 tasks total, ~15,400 core-hours.
#   Every task: 1 node, 1 task, 1 core, no MPI, no threads, no GPU.
#   Each task is one independent process reading one matrix file and writing one
#   log file. Nothing is shared between tasks and nothing communicates.
#   Per-task --mem is 4G, 8G or 16G by matrix size; --time is 00:30:00 to
#   53:50:00 by per-seed cost. Arrays are throttled with %N (9 to 96 concurrent).
#   Longest single task: 54 h (one seed of a 4.8M-nonzero matrix at 2000 CG
#   iterations in software-emulated posit arithmetic).
#
# The programs do not checkpoint, so a task killed at its wall limit loses that
# seed's work. Completed logs are skipped on resubmission, so a preempted or
# requeued task costs only the seed in flight.
#SBATCH --job-name=pcgtest
#SBATCH --partition=epyc
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=hpc/slurm/%x_%A_%a.out
#SBATCH --error=hpc/slurm/%x_%A_%a.out
#SBATCH --requeue

set -uo pipefail
REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO"

# Four cheap seeds of the smallest matrix, one per array task.
MATRIX=bcsstk03
SEED=$SLURM_ARRAY_TASK_ID
OUTDIR=results/testjob
BIN=bin/cg_compare_seeded

echo "=== task $SLURM_ARRAY_TASK_ID on $(hostname) at $(date -Is) ==="
scontrol show config 2>/dev/null | grep -i ClusterName || true
lscpu | grep -E "^Model name|^CPU\(s\):" || true

if [ ! -x "$BIN" ]; then
    echo "ERROR: $BIN not built. Run: bash hpc/build.sh"
    exit 1
fi

mkdir -p "$OUTDIR"
OUT="$OUTDIR/${MATRIX}_seed${SEED}.log"
if [ -s "$OUT" ]; then
    echo "skip $MATRIX seed $SEED (already done)"
    exit 0
fi

t0=$SECONDS
"$BIN" "data/matrices/${MATRIX}.mtx" "${OUT}.tmp" "$SEED" > /dev/null
rc=$?
if [ $rc -eq 0 ] && [ -s "${OUT}.tmp" ]; then
    mv "${OUT}.tmp" "$OUT"
    echo "done  $MATRIX seed $SEED in $((SECONDS-t0))s, $(wc -l < "$OUT") log lines"
    head -1 "$OUT"
else
    rm -f "${OUT}.tmp"
    echo "FAIL  $MATRIX seed $SEED rc=$rc"
    exit $rc
fi
