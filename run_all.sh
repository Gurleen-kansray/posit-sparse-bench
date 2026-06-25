#!/bin/bash
# Reproduce all CoNGA benchmark results
# Usage: ./run_all.sh (inside Docker container with -v /mnt/d/posits-work:/mnt/d/posits-work)

set -e
echo "=== posit-sparse-bench: CoNGA reproducibility run ==="
echo "Universal v3.80, quire<N,ES,2>, 300 CG iterations each"
echo ""

echo "--- bcsstk03 ladder ---"
./bcsstk03_ladder

echo "--- bcsstk14 ladder ---"
./bcsstk14_ladder

echo "--- bcsstk38 ladder ---"
./bcsstk38_ladder

echo "--- mhd4800b ladder ---"
./mhd4800b_ladder

echo ""
echo "=== Done. Check /mnt/d/posits-work/ for logs ==="
