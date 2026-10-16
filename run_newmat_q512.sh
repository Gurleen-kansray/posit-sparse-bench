#!/bin/bash
cd /home/acer/posit-sparse-bench
m=$1
mkdir -p results/newmat_q512
for s in $(seq 1 25); do
  f=results/newmat_q512/${m}_seed$s.log
  [ -s "$f" ] && continue
  ./cg_compare_seeded_q512 data/matrices/$m.mtx "$f.tmp" $s >/dev/null 2>&1 && mv "$f.tmp" "$f"
done
echo "$m done $(date)"
