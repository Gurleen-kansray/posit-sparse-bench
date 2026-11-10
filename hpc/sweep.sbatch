#!/bin/bash
# SLURM array job for the 50-seed x 13-matrix posit32 CG sweep.
#
# Submitted by hpc/submit.sh, which sets --array, --time and JOBLIST per
# size class. Do not sbatch this directly unless you pass those yourself.
#
# One task = one (matrix, seed) pair = one independent process. No
# communication, no shared state, so the array scales to whatever the
# scheduler gives it.
#
# NOTE: --account is supplied by submit.sh from the account in hpc/config.sh. Most clusters require
# it for usage tracking and billing, and a job without it is rejected.
#
# NOTE: Slurm does NOT create the output directory. If hpc/slurm/ is missing
# the job fails before it runs, with JobLaunchFailure and no explanation in the
# log (there is no log). submit.sh mkdir -p's it; if you sbatch this file by
# hand, create it first.
#SBATCH --job-name=pcgsweep
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=hpc/slurm/%x_%A_%a.out
#SBATCH --error=hpc/slurm/%x_%A_%a.out
#SBATCH --requeue

set -uo pipefail

REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO"

JOBLIST="${JOBLIST:?set JOBLIST=hpc/joblist_<class>.txt}"
BIN="${BIN:-bin/cg_compare_seeded}"
OUTDIR="${OUTDIR:-results/seed_sweep_50}"

# Job lists are "matrix seed" or "matrix seed maxiter". The third field is
# optional: absent means the program's default budget of 2000.
read -r MATRIX SEED MAXITER < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
if [ -z "${MATRIX:-}" ]; then
    echo "no job at line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

mkdir -p "$OUTDIR"
# A non-default budget goes in the filename, so a 4000-iteration run never
# overwrites the 2000-iteration one it is being compared against.
if [ -n "${MAXITER:-}" ]; then
    OUT="$OUTDIR/${MATRIX}_seed${SEED}_it${MAXITER}.log"
else
    OUT="$OUTDIR/${MATRIX}_seed${SEED}.log"
fi

# Idempotent: a completed log is never recomputed, so a requeued or
# preempted task costs only the seed in flight. Write to .tmp and move, so
# a kill mid-write cannot leave a half-log that looks finished.
if [ -s "$OUT" ]; then
    echo "skip  $MATRIX seed $SEED (exists)"
    exit 0
fi

echo "start $MATRIX seed $SEED maxiter=${MAXITER:-default}  task=$SLURM_ARRAY_TASK_ID  host=$(hostname)  $(date -Is)"
t0=$SECONDS
"$REPO/$BIN" "data/matrices/${MATRIX}.mtx" "${OUT}.tmp" "$SEED" ${MAXITER:+$MAXITER} > /dev/null
rc=$?

# Three separate failures, three nonzero exits. A program that returns 0 while
# writing nothing must not be recorded as a success: `exit $rc` would exit 0 and
# the task would show COMPLETED in sacct with no result behind it.
if [ $rc -ne 0 ]; then
    rm -f "${OUT}.tmp"
    echo "FAIL  $MATRIX seed $SEED  rc=$rc  $((SECONDS-t0))s"
    exit $rc
fi
if [ ! -s "${OUT}.tmp" ]; then
    rm -f "${OUT}.tmp"
    echo "FAIL  $MATRIX seed $SEED  empty or missing result  $((SECONDS-t0))s"
    exit 90
fi
if ! mv "${OUT}.tmp" "$OUT"; then
    echo "FAIL  $MATRIX seed $SEED  mv to $OUT failed  $((SECONDS-t0))s"
    exit 91
fi
echo "done  $MATRIX seed $SEED  $((SECONDS-t0))s"
