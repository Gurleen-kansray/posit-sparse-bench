#!/bin/bash
# usage: table10_worker.sh <binary> <outdir> <matrix> <seed>
cd /home/acer/posit-sparse-bench
bin=$1; out=$2; m=$3; s=$4
mkdir -p results/table10_q512/$out
f=results/table10_q512/$out/${m}_seed$s.log
[ -s "$f" ] && exit 0
./$bin data/matrices/$m.mtx "$f.tmp" $s >/dev/null 2>&1 && mv "$f.tmp" "$f"
