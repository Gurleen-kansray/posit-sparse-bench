#include <stddef.h>
#include <cstdint>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

extern "C" int dotX_posit(double *C, double *A, double *B, int n)
{
    quire<32,2,2> q = 0;
    for (int i = 0; i < n; ++i) {
        posit<32,2> pa(A[i]);
        posit<32,2> pb(B[i]);
        q += quire_mul(pa, pb);
    }
    posit<32,2> presult;
    convert(q.to_value(), presult);
    *C = double(presult);
    return 1;
}
