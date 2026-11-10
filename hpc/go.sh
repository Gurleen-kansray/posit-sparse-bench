#!/bin/bash
# One command to run everything.
#
#   bash hpc/go.sh
#
#  If you want the main results sooner, submit in two waves:
#
#  bash hpc/submit.sh small medium large      # tables: ~1 day
#  bash hpc/submit.sh maxn budget2k           # the extras, later
#
#
# Builds the binaries against the vendored Universal headers, checks the quire
# is the standard 512 bits, then submits all five job arrays. Prints what to
# type to check progress. Safe to re-run: finished seeds are skipped, so
# running it again after a failure or a wall-time kill just fills the gaps.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
source "$REPO/hpc/config.sh"

echo "=============================================================="
echo " posit32 CG sweep  --  50 seeds x 13 matrices, all experiments"
echo "=============================================================="
echo

# ---- 1. build -------------------------------------------------------------
if [ ! -x bin/cg_compare_seeded ]; then
    echo "[1/3] building"
    if ! bash hpc/build.sh; then
        echo
        echo "BUILD FAILED. Nothing was submitted."
        echo "Most likely cause: no GCC with -std=c++20. Try:"
        echo "    module avail gcc"
        echo "then edit the 'module load' line in hpc/build.sh."
        exit 1
    fi
else
    echo "[1/3] binaries already built (delete bin/ to force a rebuild)"
fi
echo

# ---- 2. submit ------------------------------------------------------------
echo "[2/3] submitting job arrays"
if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found -- are you on a login node?"
    exit 1
fi
# Regenerate the job lists so the manifest matches the current cost model and
# wall-time ceiling, then submit every array. The manifest is already ordered
# longest-task-first: a task cannot be split, so the longest one floors the finish
# time, and starting those first lets the cheap tasks backfill around them.
python3 hpc/make_joblists.py || exit 1
echo
bash hpc/submit.sh || exit 1
echo

# ---- 3. what to do next ---------------------------------------------------
cat <<'EOF'
[3/3] submitted. Nothing else to do until it finishes.

Check progress (any time, from the login node):

    cd <repo> && bash hpc/progress.sh

That prints a bar per matrix out of 50. Expect roughly 24-36 hours for
everything, depending on how much of the cluster you get.

When every matrix reads 50/50, download from your laptop:

    rsync -avPz <you>@<host>:<repo>/results/seed_sweep_50/ \
        ~/Documents/research/pcg_gurleen/posit-sparse-bench/results/seed_sweep_50/

If some seeds are missing because tasks hit the wall-time limit, do NOT run this
script again: it rebuilds the full manifest and resubmits all 27 arrays. Use the
requeue path instead, which submits only the gaps, with a longer --time and the
measured memory request:

    bash hpc/audit.sh                # why they died; MaxRSS vs ReqMem
    python3 hpc/requeue.py           # dry run
    python3 hpc/requeue.py --write
    bash hpc/submit.sh --requeue
EOF
