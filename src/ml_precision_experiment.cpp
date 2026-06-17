// ml_precision_experiment.cpp
//
// Tests Gustafson's HPEC claim ("exact dot products obviate the need for
// mixed precision") at ML scale, not FEM scale.
//
// Setup mirrors a real neural net layer: weights w ~ N(0, 1/L) (He/Xavier-style
// init scaled so output variance stays O(1) regardless of layer width L),
// activations x ~ N(0,1) clipped to [-3,3] (post-LayerNorm/BatchNorm range).
// This is the standard "mixed precision training" pattern: fp16 storage,
// fp32 accumulate, because long dot-product reduction in fp16 alone drifts
// due to rounding/cancellation even though every individual value is small
// and well-scaled.
//
// Compared, for dot-product length L and many random trials:
//   ref          - double64 accumulate (ground truth)
//   fp16_naive   - fp16 storage, fp16 accumulate (the thing nobody does, because it's bad)
//   fp32_mixed   - fp16 storage, fp32 accumulate (the actual industry-standard mixed-precision pattern)
//   posit16_naive- posit16 storage, posit16 accumulate
//   posit16_quire- posit16 storage, EXACT quire accumulate, round once at the end
//
// Question: does posit16_quire match or beat fp32_mixed, using only a 16-bit
// accumulator instead of a 32-bit one? If yes, that's direct evidence for the
// "obviates mixed precision" claim in the domain it's actually about.

#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

struct Stats { double max_rel; double mean_rel; };

Stats summarize(const std::vector<double>& rel_errs) {
    double mx = 0, sum = 0;
    for (double e : rel_errs) { mx = std::max(mx, e); sum += e; }
    return { mx, sum / rel_errs.size() };
}

int main() {
    std::vector<int> lengths = {256, 512, 1024, 4096};
    const int TRIALS = 300;
    std::mt19937 rng(42);

    FILE* log = fopen("results/logs/ml_precision_experiment.log", "w");

    auto banner = [&](){
        printf("%-8s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s\n",
               "L", "fp16_max", "fp16_mean", "fp32mix_max", "fp32mix_mean",
               "p16naive_max", "p16naive_mean", "p16quire_max", "p16quire_mean");
        fprintf(log, "%-8s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s\n",
               "L", "fp16_max", "fp16_mean", "fp32mix_max", "fp32mix_mean",
               "p16naive_max", "p16naive_mean", "p16quire_max", "p16quire_mean");
    };
    banner();

    for (int L : lengths) {
        std::normal_distribution<double> wdist(0.0, 1.0 / std::sqrt((double)L));
        std::normal_distribution<double> xdist(0.0, 1.0);

        std::vector<double> e_fp16, e_fp32mix, e_p16naive, e_p16quire;
        int skipped = 0;

        for (int t = 0; t < TRIALS; t++) {
            std::vector<double> w(L), x(L);
            for (int i = 0; i < L; i++) {
                w[i] = wdist(rng);
                double xv = xdist(rng);
                xv = std::clamp(xv, -3.0, 3.0);
                x[i] = xv;
            }

            // Reference: double64
            double ref = 0, scale = 0;
            for (int i = 0; i < L; i++) { ref += w[i] * x[i]; scale += std::fabs(w[i] * x[i]); }

            // Skip trials where catastrophic cancellation makes relative error
            // meaningless for ALL methods (same filter rationale as the FEM work:
            // near-zero denominator inflates relative error without reflecting
            // precision quality)
            if (std::fabs(ref) < 1e-3 * scale) { skipped++; continue; }

            // fp16 naive: storage + accumulate both in _Float16
            {
                _Float16 acc = 0;
                for (int i = 0; i < L; i++) {
                    _Float16 wi = (_Float16)w[i];
                    _Float16 xi = (_Float16)x[i];
                    acc = acc + (_Float16)(wi * xi);
                }
                double got = (double)acc;
                e_fp16.push_back(std::fabs(got - ref) / std::fabs(ref));
            }

            // fp32 mixed: storage in fp16, accumulate in fp32 (industry standard)
            {
                float acc = 0;
                for (int i = 0; i < L; i++) {
                    float wi = (float)(_Float16)w[i];
                    float xi = (float)(_Float16)x[i];
                    acc += wi * xi;
                }
                double got = (double)acc;
                e_fp32mix.push_back(std::fabs(got - ref) / std::fabs(ref));
            }

            // posit16 naive: storage + accumulate both in posit16
            {
                posit<16,1> acc(0);
                for (int i = 0; i < L; i++) {
                    posit<16,1> wi(w[i]), xi(x[i]);
                    posit<16,1> prod = wi * xi;
                    acc = acc + prod;
                }
                double got = (double)acc;
                e_p16naive.push_back(std::fabs(got - ref) / std::fabs(ref));
            }

            // posit16 + quire: storage in posit16, EXACT accumulate, round once
            {
                quire<16,1,2> q = 0;
                for (int i = 0; i < L; i++) {
                    posit<16,1> wi(w[i]), xi(x[i]);
                    q += quire_mul(wi, xi);
                }
                posit<16,1> result;
                convert(q.to_value(), result);
                double got = (double)result;
                e_p16quire.push_back(std::fabs(got - ref) / std::fabs(ref));
            }
        }

        auto s1 = summarize(e_fp16);
        auto s2 = summarize(e_fp32mix);
        auto s3 = summarize(e_p16naive);
        auto s4 = summarize(e_p16quire);

        printf("%-8d %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e  (skipped=%d/%d)\n",
               L, s1.max_rel, s1.mean_rel, s2.max_rel, s2.mean_rel,
               s3.max_rel, s3.mean_rel, s4.max_rel, s4.mean_rel, skipped, TRIALS);
        fprintf(log, "%-8d %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e  (skipped=%d/%d)\n",
               L, s1.max_rel, s1.mean_rel, s2.max_rel, s2.mean_rel,
               s3.max_rel, s3.mean_rel, s4.max_rel, s4.mean_rel, skipped, TRIALS);
        fflush(stdout); fflush(log);
    }

    fclose(log);
    printf("\nLog saved: results/logs/ml_precision_experiment.log\n");
    return 0;
}
