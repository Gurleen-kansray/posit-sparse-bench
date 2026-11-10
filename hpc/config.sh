#!/bin/bash
# Site configuration for the SLURM sweep harness. Every other script in hpc/
# sources this file, so this is the only place to edit when moving to a new
# cluster. Nothing here is specific to one site except the defaults.
#
#   source hpc/config.sh
#
# Anything already exported in your environment wins, so a one-off override is
#   PARTITION=haswell bash hpc/submit.sh
# and a permanent one belongs in ~/.bashrc, not here.

# ---------------------------------------------------------------- directories
# Where this harness lives, relative to the repo root. Change it only if you
# rename the directory, and see the note on #SBATCH below before you do.
: "${HPC_DIR:=hpc}"

# Slurm's own per-task .out files. Slurm does NOT create this directory, and a
# job whose output directory is missing fails before it runs, with
# JobLaunchFailure and no log to explain it. submit.sh mkdir -p's it.
: "${SLURM_OUT:=$HPC_DIR/slurm}"

# ------------------------------------------------------------------ scheduler
# Partition to submit to. Measured throughput on ARCSIM Penobscot, relative to
# an i7-12800H laptop core: epyc 2.00x, epyc-milan 1.97x, epyc-genoa 1.70x,
# skylake 2.60x, haswell 3.23x. We use epyc despite haswell being faster per
# core: haswell nodes have 62 GB, which fits three 16 GB tasks, against epyc's
# 503 GB and 31 tasks. Memory, not clock, set the throughput here.
: "${PARTITION:=epyc}"

# Account for usage tracking. Most clusters reject jobs without one:
#   "Invalid account or account/partition combination specified"
# Find yours with:  sacctmgr show assoc where user=$USER format=account
# Slurm reads SBATCH_ACCOUNT natively, so exporting that in ~/.bashrc covers
# every sbatch you ever run, not just this harness. HPC_ACCOUNT and the older
# ARCSIM_ACCOUNT are honored first, in that order.
: "${HPC_ACCOUNT:=${ARCSIM_ACCOUNT:-${SBATCH_ACCOUNT:-}}}"

# Compiler module. `module avail gcc` lists what a site offers; any GCC with
# -std=c++20 works, which means 11 or newer. Empty skips the module load.
: "${GCC_MODULE:=gcc/12.2.0}"

# ----------------------------------------------------------------- cost model
# These four are read by make_joblists.py and requeue.py, which take them from
# the environment, so they must be exported rather than merely set.

# Per-core speed of the cluster relative to the machine the cost constants in
# make_joblists.py were fitted on. Recalibrate after any run with
#   bash hpc/fit_costs.sh
# which refits from the measured per-matrix times in results/costs.txt. Do not
# scale timings taken from a contended parallel laptop run: measured against one
# idle cluster core, a 12-way-contended laptop run is SLOWER, not faster, and
# scaling the wrong way projected 23 h for jobs that took 5.
: "${CLUSTER:=2.0}"

# Billing rate in dollars per core-hour, used only for the cost estimate the
# manifest prints. ARCSIM quoted $0.0051/core-hour at the Internal/Federal rate.
: "${RATE:=0.008}"

# Wall-clock ceiling in hours for any single task. A task that would exceed it
# is dropped from the manifest and reported, rather than submitted to be killed.
# Check the real limit with:  sinfo -o '%P %l'
: "${MAX_WALL_H:=168}"

# Seeds per matrix. One seed is one draw of the solution vector x*.
: "${SEEDS:=50}"

export HPC_DIR SLURM_OUT PARTITION HPC_ACCOUNT GCC_MODULE
export CLUSTER RATE MAX_WALL_H SEEDS

# ------------------------------------------------------------ #SBATCH caveat
# sbatch reads #SBATCH lines as literal text before any shell runs, so a
# variable there does NOT expand: `#SBATCH --output=$HPC_DIR/slurm/%x.out`
# creates a directory literally named '$HPC_DIR'. The .sbatch files therefore
# carry hardcoded defaults, and submit.sh overrides them on the command line
# from the values above, which is what takes effect for the production arrays.
# If you sbatch testjob.sbatch or probe.sbatch by hand, the literals in those
# two files apply, so edit them there or pass --partition and --output yourself.
