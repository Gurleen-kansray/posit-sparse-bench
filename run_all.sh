#!/bin/bash
# Reproduce all CoNGA benchmark results
# Usage inside container: bash run_all.sh
set -e
echo "=== posit-sparse-bench: CoNGA reproducibility run ==="
echo "Universal v3.80, quire<N,ES,2>, 300 CG iterations each"
echo ""

for matrix in bcsstk03 bcsstk14 bcsstk38 mhd4800b; do
    echo "--- ${matrix} ladder ---"
    ./${matrix}_ladder
    echo "Done ${matrix}"
done

echo ""
echo "=== All done. Logs in results/ ==="
