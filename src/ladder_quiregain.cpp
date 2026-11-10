// Quire gain for the two CG inner products (Table 4).
//
//   ladder_quiregain <matrix.mtx> <logfile> <seed> [maxiter]
//
// A double-driven ladder. The CG recurrence runs entirely in double; at every
// iteration the four vectors are rounded to posit32 and each inner product is
// formed twice, once by a naive sequential loop and once in the quire, both
// compared against the double-precision value. This isolates the accumulation:
// the two policies see the SAME vectors, so the ratio measures rounding and not
// trajectory divergence.
//
// That distinction is why Table 4 cannot come from cg_compare_seeded. There the
// quire and naive arms are separate CG runs whose iterates diverge, so their
// inner products differ for reasons unrelated to how the sum was accumulated.
//
// Replaces generic_ladder_seeded.cpp for this purpose, which differs from the
// definition in three ways: it seeds the residual from uniform_real(0.5, 1.5)
// rather than from b = A x_true, it never records a double reference for r'z,
// and it runs a fixed 300 iterations instead of truncating at convergence.
//
// Log format: line 1 metadata, line 2 header, line 3+ one row per iteration,
// then conv_iter and window trailers. Relative errors are |posit - double| /
// |double|, so a LARGER naive error over the quire error is a gain above 1.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

// Capacity 31: with es=2 the quire's value range is 2^es(4N-8) = 16N-32 bits, so
// a 31-bit carry guard plus the sign bit gives the standard's 16N, 512 bits for
// posit32.
template<unsigned N, unsigned ES> using Q = quire<N, ES, 31>;

struct MTX { int n; std::vector<int> row, col; std::vector<double> val; };

MTX read_mtx(const char* f) {
    MTX m; FILE* fp = fopen(f, "r");
    if (!fp) { printf("ERROR: cannot open %s\n", f); exit(1); }
    char buf[256];
    while (fgets(buf, 256, fp) && buf[0] == '%');
    int rows, cols, nnz; sscanf(buf, "%d %d %d", &rows, &cols, &nnz);
    m.n = rows; int r, c; double v;
    while (fscanf(fp, "%d %d %lf", &r, &c, &v) == 3) {
        m.row.push_back(r - 1); m.col.push_back(c - 1); m.val.push_back(v);
        if (r != c) { m.row.push_back(c - 1); m.col.push_back(r - 1); m.val.push_back(v); }
    }
    fclose(fp); return m;
}

void matvec(const MTX& A, const double* x, double* y) {
    for (int i = 0; i < A.n; i++) y[i] = 0;
    for (size_t k = 0; k < A.row.size(); k++) y[A.row[k]] += A.val[k] * x[A.col[k]];
}

double dot_d(const double* a, const double* b, int n) {
    double s = 0; for (int i = 0; i < n; i++) s += a[i] * b[i]; return s;
}

// Both dot products of a and b in posit<N,ES>, from identical rounded inputs.
template<unsigned N, unsigned ES>
void dot_both(const std::vector<double>& a, const std::vector<double>& b, int n,
              double& quire_val, double& naive_val) {
    Q<N, ES> q = 0;
    for (int i = 0; i < n; i++) { posit<N,ES> x(a[i]), y(b[i]); q += quire_mul(x, y); }
    posit<N,ES> qr; convert(q.to_value(), qr); quire_val = double(qr);
    posit<N,ES> s = 0;
    for (int i = 0; i < n; i++) { posit<N,ES> x(a[i]), y(b[i]); s = s + x * y; }
    naive_val = double(s);
}

double relerr(double v, double ref) {
    if (ref == 0.0) return (v == 0.0) ? 0.0 : INFINITY;
    return fabs(v - ref) / fabs(ref);
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        printf("usage: ladder_quiregain <matrix.mtx> <logfile> <seed> [maxiter]\n");
        return 1;
    }
    const char* mtx_path = argv[1];
    const char* log_path = argv[2];
    unsigned seed = (unsigned)atoi(argv[3]);
    const int MAXITER = (argc > 4 && atoi(argv[4]) > 0) ? atoi(argv[4]) : 300;
    const double RTOL = 1e-10;

    MTX A = read_mtx(mtx_path);
    const int n = A.n;

    std::vector<double> diagA(n, 1.0);
    for (size_t k = 0; k < A.row.size(); k++)
        if (A.row[k] == A.col[k]) diagA[A.row[k]] = A.val[k];

    // Same construction as cg_compare_seeded: a known x_true, then b = A x_true,
    // so a given seed means the same right-hand side across the whole suite.
    std::mt19937 rng(seed);
    std::normal_distribution<double> ndist(0.0, 1.0);
    std::vector<double> x_true(n);
    for (int i = 0; i < n; i++) x_true[i] = ndist(rng);

    std::vector<double> b(n), x(n, 0.0), r(n), z(n), p(n), Ap(n);
    matvec(A, x_true.data(), b.data());
    for (int i = 0; i < n; i++) r[i] = b[i];
    for (int i = 0; i < n; i++) z[i] = r[i] / diagA[i];
    for (int i = 0; i < n; i++) p[i] = z[i];
    double rz_d = dot_d(r.data(), z.data(), n);
    const double bn = sqrt(dot_d(b.data(), b.data(), n));

    FILE* log = fopen(log_path, "w");
    if (!log) { printf("ERROR: cannot write %s\n", log_path); return 1; }
    fprintf(log, "matrix=%s n=%d seed=%u maxiter=%d rtol=%.1e\n",
            mtx_path, n, seed, MAXITER, RTOL);
    fprintf(log, "iter pAp_rel_n pAp_rel_q rz_rel_n rz_rel_q resid_rel\n");

    int conv = -1;
    for (int iter = 0; iter < MAXITER; iter++) {
        matvec(A, p.data(), Ap.data());
        const double pAp_d = dot_d(p.data(), Ap.data(), n);

        double pq, pn, rq, rn;
        dot_both<32,2>(p, Ap, n, pq, pn);
        dot_both<32,2>(r, z,  n, rq, rn);

        const double rr = sqrt(dot_d(r.data(), r.data(), n)) / bn;
        fprintf(log, "%d %.10e %.10e %.10e %.10e %.10e\n", iter,
                relerr(pn, pAp_d), relerr(pq, pAp_d),
                relerr(rn, rz_d),  relerr(rq, rz_d), rr);

        const double alpha = rz_d / pAp_d;
        for (int i = 0; i < n; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
        for (int i = 0; i < n; i++) z[i] = r[i] / diagA[i];
        const double rz2 = dot_d(r.data(), z.data(), n);
        const double beta = rz2 / rz_d;
        rz_d = rz2;
        for (int i = 0; i < n; i++) p[i] = z[i] + beta * p[i];

        if (conv < 0 && sqrt(dot_d(r.data(), r.data(), n)) / bn < RTOL) conv = iter;
    }

    // The statistics window is the pre-convergence stretch: iterations 0 through
    // conv inclusive, or all of them if the reference never converged. Rows past
    // the window stay in the log so the choice can be revisited without rerunning.
    const int window = (conv < 0) ? MAXITER : conv + 1;
    fprintf(log, "conv_iter=%d\n", conv);
    fprintf(log, "window=%d\n", window);
    fclose(log);
    printf("Done. %s seed=%u conv=%d window=%d -> %s\n",
           mtx_path, seed, conv, window, log_path);
    return 0;
}
