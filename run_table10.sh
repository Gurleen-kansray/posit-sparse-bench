#!/bin/bash
# usage: run_table10.sh <binary> <outdir> <nseeds>
cd /home/acer/posit-sparse-bench
for m in nasa4704 bodyy4 sts4098 bcsstk14 nos2 mhd4800b bcsstk03; do for s in $(seq 1 $3); do echo "$m $s"; done; done \
 | xargs -P 4 -L1 ./table10_worker.sh $1 $2
echo "$2 n=$3 finished $(date)"
