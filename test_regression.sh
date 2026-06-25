#!/bin/bash
# Regression tests
set -e
echo "=== Regression Tests ==="

make INCLUDES=-I/mnt/d/posits-work/universal/include bcsstk03_ladder 2>/dev/null

# Test 1: bcsstk03 iter=0 exact value check
./bcsstk03_ladder
GOT=$(grep "iter=0" results/bcsstk03_ladder.log | grep -o "p32q=[^ ]*" | cut -d= -f2)
EXPECTED="2.1359165657e-05"
if [ "$GOT" = "$EXPECTED" ]; then
    echo "PASS: bcsstk03 p32q exact match"
else
    echo "FAIL: bcsstk03 p32q got $GOT expected $EXPECTED"
    exit 1
fi

# Test 2: check log has 300 iterations
COUNT=$(grep -c "^iter=" results/bcsstk03_ladder.log)
if [ "$COUNT" -eq 300 ]; then
    echo "PASS: bcsstk03 has 300 iterations"
else
    echo "FAIL: bcsstk03 has $COUNT iterations (expected 300)"
    exit 1
fi

# Test 3: p64q should match double exactly at iter=0
P64Q=$(grep "iter=0" results/bcsstk03_ladder.log | grep -o "p64q=[^ ]*" | cut -d= -f2)
PAP_D=$(grep "iter=0" results/bcsstk03_ladder.log | grep -o "pAp_d=[^ ]*" | cut -d= -f2)
if [ "$P64Q" = "$PAP_D" ]; then
    echo "PASS: p64q matches double64 at iter=0"
else
    echo "FAIL: p64q=$P64Q != pAp_d=$PAP_D"
    exit 1
fi

echo "=== All tests passed ==="
