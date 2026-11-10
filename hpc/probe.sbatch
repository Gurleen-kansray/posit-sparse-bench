#!/bin/bash
# Reports what you actually get from a partition: which cluster, which node,
# which CPU, how many cores, which compilers. Run this before committing 509
# tasks, because ARCSIM's epyc and haswell partitions exist on both Penobscot
# and Katahdin and the Getting Started guide does not say how to choose.
#
#   mkdir -p hpc/slurm
#   sbatch --account=$SBATCH_ACCOUNT --partition=epyc hpc/probe.sbatch
#   # then, when it finishes:
#   cat hpc/slurm/probe_*.out
#
# Try it per partition to compare: epyc, epyc-hm, skylake, haswell.
#SBATCH --job-name=probe
#SBATCH --output=hpc/slurm/%x_%j.out
#SBATCH --error=hpc/slurm/%x_%j.out
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

echo "=== cluster and node ==="
scontrol show config 2>/dev/null | grep -i "ClusterName" || echo "ClusterName: unavailable"
echo "hostname        : $(hostname)"
echo "partition       : ${SLURM_JOB_PARTITION:-?}"
echo "nodelist        : ${SLURM_JOB_NODELIST:-?}"

echo
echo "=== cpu ==="
lscpu | grep -E "^Model name|^CPU\(s\):|^Core\(s\) per socket|^Socket|^Thread|^CPU max MHz|^NUMA node\(s\)"

echo
echo "=== memory ==="
free -g | head -2

echo
echo "=== compilers available ==="
module avail gcc 2>&1 | head -12
echo "--- currently loaded ---"
module list 2>&1 | head -5
echo "--- default g++ ---"
g++ --version 2>/dev/null | head -1 || echo "no g++ in PATH without a module"

echo
echo "=== julia ==="
module avail julia 2>&1 | head -6

echo
echo "=== this partition's limits ==="
scontrol show partition "${SLURM_JOB_PARTITION}" 2>/dev/null \
    | tr ' ' '\n' | grep -E "MaxTime|DefaultTime|MaxNodes|TotalCPUs|TotalNodes" \
    || echo "unavailable"

echo
echo "=== single-core speed check (lower is faster) ==="
# Rough integer/FP benchmark, comparable across partitions. The sweep is
# single-threaded and ALU-bound, so this predicts per-seed cost better than
# clock speed does. Reference: 3.0 s on an i7-12800H (the workstation the cost
# model was calibrated on), so a result of 5.4 s means roughly 1.8x the per-seed
# times quoted in README.md -- which is the factor those estimates already assume.
python3 -c "
import time
t0 = time.perf_counter()
s = 0.0
for i in range(1, 60000000):
    s += (i % 7) * 1.0000001
print('  elapsed: %.1f s   (i7-12800H reference: 3.0 s)' % (time.perf_counter() - t0))
" 2>/dev/null || echo "  (no python3; skipped)"
