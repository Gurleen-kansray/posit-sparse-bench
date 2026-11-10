#!/bin/bash
# Build the sweep binaries. Run once from the repo root:
#   bash hpc/build.sh
#
# Universal is header-only, so it needs a clone and an -I, nothing more. The
# only real requirement is a compiler with -std=c++20; GCC 11+ is enough.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
source "$REPO/hpc/config.sh"

# GCC_MODULE comes from hpc/config.sh; `module avail gcc` lists the choices.
if [ -n "${GCC_MODULE:-}" ]; then
    module load "$GCC_MODULE" 2>/dev/null || module load gcc 2>/dev/null || true
fi
echo "compiler: $(g++ --version | head -1)"

# Universal headers.
#
# DO NOT clone upstream Universal here. These programs use the pre-generalized
# quire API, `quire<nbits, es, capacity>`. Current upstream main replaced it with
# `quire<NumberType, capacity, LimbType>`, so a fresh clone fails with
#   "type/value mismatch at argument 1 ... expected a type, got 32"
# on every `quire<32,2,31>`. The headers that built the reported results are the
# older release installed at /usr/local/include/universal on the workstation;
# ship those instead, which also removes "which Universal version" as a variable.
#
# Create the tarball on the workstation:
#   tar czf universal-headers.tgz -C /usr/local/include universal
# Copy it next to this script, then:
#   rsync -avP universal-headers.tgz <you>@<host>:posit-sparse-bench/hpc/
UNIVERSAL_INC="${UNIVERSAL_INC:-$REPO/external/universal-headers}"
if [ ! -d "$UNIVERSAL_INC/universal" ]; then
    if [ -f "$REPO/hpc/universal-headers.tgz" ]; then
        echo "unpacking hpc/universal-headers.tgz"
        mkdir -p "$UNIVERSAL_INC"
        tar xzf "$REPO/hpc/universal-headers.tgz" -C "$UNIVERSAL_INC"
    elif [ -d /usr/local/include/universal ]; then
        echo "using system headers at /usr/local/include"
        UNIVERSAL_INC=/usr/local/include
    else
        echo "ERROR: no Universal headers found."
        echo "  expected $REPO/hpc/universal-headers.tgz"
        echo "  build it on the workstation with:"
        echo "    tar czf universal-headers.tgz -C /usr/local/include universal"
        exit 1
    fi
fi
echo "universal include: $UNIVERSAL_INC"

# Fail loudly if the API is the generalized one: better than 200 lines of
# template errors per program.
if ! printf '#include <universal/number/posit/posit.hpp>\n#include <universal/number/posit/fdp.hpp>\nint main(){ sw::universal::quire<32,2,31> q; (void)q; }\n' \
     | g++ -std=c++20 -x c++ -fsyntax-only -I"$UNIVERSAL_INC" - 2>/dev/null; then
    echo "ERROR: these headers do not support quire<32,2,31>."
    echo "  They are probably current upstream main, which uses"
    echo "  quire<posit<32,2>,31> instead. Use the older release headers."
    exit 1
fi

mkdir -p bin
for p in cg_compare_seeded cg_compare_seeded_q512_qmv cg_compare_maxn \
         generic_ladder_seeded ladder_quiregain cg_refinement_seeded_v2 \
         cg_refinement_seeded_v3_q512 cg_refinement_exactb_q512 \
         cg_refinement_noquire bnorm_tool; do
    if [ -f "src/$p.cpp" ]; then
        printf '  %-34s' "$p"
        if g++ -std=c++20 -O2 -march=native -I"$UNIVERSAL_INC" \
               -o "bin/$p" "src/$p.cpp" 2>"bin/$p.buildlog"; then
            echo "ok"
        else
            echo "FAILED (see bin/$p.buildlog)"
        fi
    fi
done

# Table 10's plain-matvec arm is the same source with -DPLAIN_MV: the matvec policy
# is a compile-time switch, not a runtime flag, so it needs its own binary.
printf '  %-34s' "qmv_plain (-DPLAIN_MV)"
if g++ -std=c++20 -O2 -march=native -DPLAIN_MV -I"$UNIVERSAL_INC" \
       -o bin/qmv_plain src/cg_compare_seeded_q512_qmv.cpp 2>bin/qmv_plain.buildlog; then
    echo "ok"
else
    echo "FAILED (see bin/qmv_plain.buildlog)"
fi

# Sanity check: the quire must be the standard 512 bits, not 511 or 513.
# quire_size<nbits,es,capacity>() returns range + capacity, i.e. the quire
# WITHOUT the sign bit, which Universal manages separately. The standard's
# width is 16*nbits, so for posit32: 480 + 31 + 1 = 512. Note quire_traits
# does not exist in these headers -- that is the newer generalized API.
cat > /tmp/qcheck_$$.cpp <<'EOF'
#include <cstdio>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;
int main(){
    int qbits = quire_size<32,2,31>();   // range + capacity
    int total = qbits + 1;               // + sign
    printf("posit32 quire: qbits=%d + sign = %d bits  %s\n", qbits, total,
           total == 512 ? "== 16n OK" : "!= 16n MISMATCH");
    return total == 512 ? 0 : 1;
}
EOF
if g++ -std=c++20 -O1 -I"$UNIVERSAL_INC" -o /tmp/qcheck_$$ /tmp/qcheck_$$.cpp 2>/tmp/qcheck_$$.err; then
    if /tmp/qcheck_$$; then
        rm -f /tmp/qcheck_$$ /tmp/qcheck_$$.cpp /tmp/qcheck_$$.err
    else
        echo "STOP: the quire is not 512 bits. Do not run the sweep."
        rm -f /tmp/qcheck_$$ /tmp/qcheck_$$.cpp /tmp/qcheck_$$.err
        exit 1
    fi
else
    echo "WARNING: quire width check did not compile:"
    head -3 /tmp/qcheck_$$.err
    rm -f /tmp/qcheck_$$.cpp /tmp/qcheck_$$.err
fi

echo
echo "built into $REPO/bin -- now: bash hpc/submit.sh"
