// posit_precision_curve.cpp
// Computes relative precision (ULP/value) of posit<32,2> vs float32 across a magnitude sweep.
// Answers: at what magnitude does posit32's local precision fall below float32's fixed 2^-23?
#include <cstdio>
#include <cmath>
#include <universal/number/posit/posit.hpp>
using namespace sw::universal;
using P32 = posit<32,2>;

int main(){
    printf("magnitude      posit32_relprec   float32_relprec   posit32_better\n");
    // sweep magnitudes from 1e-20 to 1e20, log-spaced
    for(double exp10=-20; exp10<=20; exp10+=0.5){
        double mag = pow(10.0, exp10);
        P32 pv(mag);
        P32 pv_next = pv;
        // increment to next representable posit (ULP step)
        pv_next++;
        double posit_ulp = double(pv_next) - double(pv);
        double posit_relprec = fabs(posit_ulp / double(pv));

        float fv = (float)mag;
        float fv_next = std::nextafter(fv, INFINITY);
        double float_ulp = (double)fv_next - (double)fv;
        double float_relprec = fabs(float_ulp / (double)fv);

        printf("%.3e   %.6e     %.6e     %s\n",
            mag, posit_relprec, float_relprec,
            (posit_relprec < float_relprec) ? "posit32" : "float32");
    }
    return 0;
}
