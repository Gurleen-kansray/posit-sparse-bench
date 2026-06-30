#!/bin/bash
set -e
echo "=== posit-sparse-bench: CoNGA reproducibility run ==="
echo "Universal v3.80, quire<N,ES,2>, 300 CG iterations each"
echo ""
for matrix in bcsstk03 bcsstk14 bcsstk36 bcsstk37 bcsstk38 mhd4800b nasasrb nasa4704 nos2 bodyy4 s3dkt3m2 s3dkq4m2 sts4098; do
    echo "--- ${matrix} ladder ---"
    ./generic_ladder data/matrices/${matrix}.mtx results/${matrix}_ladder.log
    echo "Done ${matrix}"
done
echo ""
echo "=== All done. Logs in results/ ==="
