#!/bin/bash
set -e
mkdir -p results/ladder_logs/seed_sweep
for mat in bcsstk03 bcsstk14 bcsstk36 bcsstk37 bcsstk38 nasasrb mhd4800b s3dkt3m2 s3dkq4m2 sts4098 nasa4704 nos2 bodyy4; do
  echo "=== $mat ==="
  for seed in $(seq 1 50); do
    ./cg_compare_seeded data/matrices/${mat}.mtx results/ladder_logs/seed_sweep/${mat}_seed${seed}.log ${seed}
  done
done
python3 aggregate_seeded_results.py
