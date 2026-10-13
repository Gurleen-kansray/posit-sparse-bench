// Iterative refinement using quire-exact residual (John Gustafson's suggestion, Aug 13 email)
// Properly integrated version -- matches log conventions of cg_compare_seeded / cg_compare_cond_probe
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <numeric>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

using P32 = posit<32,2>;

struct MTX { int n; std::vector<int> row, col; std::vector<double> val; };

MTX read_mtx(const char* path){
    MTX A; FILE* fp = fopen(path, "r");
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int r,c,nnz; sscanf(buf,"%d %d %d",&r,&c,&nnz);
    A.n = r;
    int rr, cc; double v;
    while(fscanf(fp,"%d %d %lf",&rr,&cc,&v)==3){
        A.row.push_back(rr-1); A.col.push_back(cc-1); A.val.push_back(v);
        if(rr != cc){ A.row.push_back(cc-1); A.col.push_back(rr-1); A.val.push_back(v); }
    }
    fclose(fp);
    return A;
}

void matvec_d(const MTX& A, const std::vector<double>& x, std::vector<double>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]] += A.val[k]*x[A.col[k]];
}
void matvec_p32n(const MTX& A, const std::vector<P32>& x, std::vector<P32>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]=y[A.row[k]]+P32(A.val[k])*x[A.col[k]];
}
void exact_residual_quire(const MTX& A, const std::vector<P32>& x0,
                           const std::vector<P32>& b, std::vector<double>& d_out){
    int n = A.n;
    std::vector<quire<32,2,2>> q(n);
    for(int i=0;i<n;i++) q[i]=0;
    for(int k=0;k<(int)A.row.size();k++)
        q[A.row[k]] += quire_mul(P32(A.val[k]), x0[A.col[k]]);
    for(int i=0;i<n;i++){
        double Ax0_i = double(q[i].to_value());
        d_out[i] = double(b[i]) - Ax0_i;
    }
}
std::vector<double> cg_solve_double(const MTX& A, const std::vector<double>& b,
                                     int maxiter=2000, double tol=1e-10){
    int n = A.n;
    std::vector<double> x(n,0.0), r(b), diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++) if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];
    double bnorm = std::sqrt(std::inner_product(b.begin(),b.end(),b.begin(),0.0));
    if (bnorm < 1e-300) bnorm = 1.0;
    std::vector<double> z(n), p(n), Ap(n);
    for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
    p = z;
    double rz_old = 0; for(int i=0;i<n;i++) rz_old += r[i]*z[i];
    for(int iter=0; iter<maxiter; iter++){
        matvec_d(A, p, Ap);
        double pAp=0; for(int i=0;i<n;i++) pAp += p[i]*Ap[i];
        if (std::fabs(pAp) < 1e-300) break;
        double alpha = rz_old/pAp;
        for(int i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }
        double rnorm=0; for(int i=0;i<n;i++) rnorm+=r[i]*r[i];
        if (std::sqrt(rnorm)/bnorm < tol) break;
        for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
        double rz_new=0; for(int i=0;i<n;i++) rz_new += r[i]*z[i];
        double beta = rz_new/rz_old;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
        rz_old = rz_new;
    }
    return x;
}
std::vector<P32> cg_solve_posit(const MTX& A, const std::vector<P32>& b,
                                 int maxiter=2000, double tol=1e-6){
    int n = A.n;
    std::vector<P32> x(n,P32(0)), r(b), diagA(n,P32(1));
    for(int k=0;k<(int)A.row.size();k++) if(A.row[k]==A.col[k]) diagA[A.row[k]]=P32(A.val[k]);
    double bnorm_p=0; for(int i=0;i<n;i++) bnorm_p += double(b[i])*double(b[i]); bnorm_p = std::sqrt(bnorm_p);
    if (bnorm_p < 1e-300) bnorm_p = 1.0;
    std::vector<P32> z(n), p(n), Ap(n);
    for(int i=0;i<n;i++) z[i]=P32(1.0/double(diagA[i]))*r[i];
    p = z;
    P32 rz_old=P32(0); for(int i=0;i<n;i++) rz_old = rz_old + r[i]*z[i];
    for(int iter=0; iter<maxiter; iter++){
        if (iter % 50 == 0) { fprintf(stderr, "[cg_solve_posit] n=%d iter=%d/%d\n", n, iter, maxiter); fflush(stderr); }
        matvec_p32n(A, p, Ap);
        P32 pAp=P32(0); for(int i=0;i<n;i++) pAp = pAp + p[i]*Ap[i];
        if (double(pAp) == 0.0) break;
        P32 alpha = rz_old/pAp;
        for(int i=0;i<n;i++){ x[i]=x[i]+alpha*p[i]; r[i]=r[i]-alpha*Ap[i]; }
        double rnorm=0; for(int i=0;i<n;i++) rnorm += double(r[i])*double(r[i]);
        if (std::sqrt(rnorm)/bnorm_p < tol) break;
        for(int i=0;i<n;i++) z[i]=P32(1.0/double(diagA[i]))*r[i];
        P32 rz_new=P32(0); for(int i=0;i<n;i++) rz_new = rz_new + r[i]*z[i];
        P32 beta = rz_new/rz_old;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
        rz_old = rz_new;
    }
    return x;
}

int main(int argc, char** argv){
    if (argc < 4){ printf("usage: %s matrix.mtx logfile seed\n", argv[0]); return 1; }
    MTX A = read_mtx(argv[1]);
    const char* logfile = argv[2];
    int seed = atoi(argv[3]);
    int n = A.n;
    std::mt19937 rng(seed);
    std::normal_distribution<double> ndist(0.0,1.0);
    std::vector<double> x_true(n);
    for(int i=0;i<n;i++) x_true[i] = ndist(rng);
    std::vector<double> b_d(n);
    matvec_d(A, x_true, b_d);
    double x_true_norm = std::sqrt(std::inner_product(x_true.begin(),x_true.end(),x_true.begin(),0.0));
    std::vector<P32> b_p(n), x0_p(n,P32(0)), diagA_p(n,P32(1));
    for(int i=0;i<n;i++) b_p[i]=P32(b_d[i]);
    for(int k=0;k<(int)A.row.size();k++) if(A.row[k]==A.col[k]) diagA_p[A.row[k]]=P32(A.val[k]);
    std::vector<P32> r_p(b_p), z_p(n), p_p(n), Ap_p(n);
    for(int i=0;i<n;i++) z_p[i]=P32(1.0/double(diagA_p[i]))*r_p[i];
    p_p = z_p;
    P32 rz_old=P32(0); for(int i=0;i<n;i++) rz_old = rz_old + r_p[i]*z_p[i];
    int maxiter=2000;
    double bnorm_main=0; for(int i=0;i<n;i++) bnorm_main += b_d[i]*b_d[i]; bnorm_main = std::sqrt(bnorm_main);
    if (bnorm_main < 1e-300) bnorm_main = 1.0;
    for(int iter=0; iter<maxiter; iter++){
        if (iter % 50 == 0) { fprintf(stderr, "[main initial posit solve] seed=%d n=%d iter=%d/%d\n", seed, n, iter, maxiter); fflush(stderr); }
        matvec_p32n(A, p_p, Ap_p);
        P32 pAp=P32(0); for(int i=0;i<n;i++) pAp = pAp + p_p[i]*Ap_p[i];
        if (double(pAp) == 0.0) break;
        P32 alpha = rz_old/pAp;
        for(int i=0;i<n;i++){ x0_p[i]=x0_p[i]+alpha*p_p[i]; r_p[i]=r_p[i]-alpha*Ap_p[i]; }
        double rnorm=0; for(int i=0;i<n;i++) rnorm += double(r_p[i])*double(r_p[i]);
        if (std::sqrt(rnorm)/bnorm_main < 1e-6) break;
        for(int i=0;i<n;i++) z_p[i]=P32(1.0/double(diagA_p[i]))*r_p[i];
        P32 rz_new=P32(0); for(int i=0;i<n;i++) rz_new = rz_new + r_p[i]*z_p[i];
        P32 beta = rz_new/rz_old;
        for(int i=0;i<n;i++) p_p[i]=z_p[i]+beta*p_p[i];
        rz_old = rz_new;
    }
    double se0=0; for(int i=0;i<n;i++){ double e=double(x0_p[i])-x_true[i]; se0+=e*e; }
    double solerr_before = std::sqrt(se0)/x_true_norm;

    std::vector<double> d(n);
    exact_residual_quire(A, x0_p, b_p, d);

    std::vector<double> y_d = cg_solve_double(A, d);
    std::vector<double> x_refined_d(n);
    for(int i=0;i<n;i++) x_refined_d[i] = double(x0_p[i]) + y_d[i];
    double se1=0; for(int i=0;i<n;i++){ double e=x_refined_d[i]-x_true[i]; se1+=e*e; }
    double solerr_after_double = std::sqrt(se1)/x_true_norm;

    std::vector<P32> d_p(n);
    for(int i=0;i<n;i++) d_p[i]=P32(d[i]);
    std::vector<P32> y_p = cg_solve_posit(A, d_p);
    std::vector<double> x_refined_p(n);
    for(int i=0;i<n;i++) x_refined_p[i] = double(x0_p[i]) + double(y_p[i]);
    double se2=0; for(int i=0;i<n;i++){ double e=x_refined_p[i]-x_true[i]; se2+=e*e; }
    double solerr_after_posit = std::sqrt(se2)/x_true_norm;

    FILE* out = fopen(logfile, "w");
    fprintf(out, "matrix=%s n=%d seed=%d\n", argv[1], n, seed);
    fprintf(out, "solerr_before solerr_after_double solerr_after_posit improvement_double improvement_posit\n");
    fprintf(out, "%.10e %.10e %.10e %.4f %.4f\n",
            solerr_before, solerr_after_double, solerr_after_posit,
            solerr_before/solerr_after_double, solerr_before/solerr_after_posit);
    fclose(out);

    printf("matrix=%s n=%d seed=%d\n", argv[1], n, seed);
    printf("solerr_before_refinement       = %.10e\n", solerr_before);
    printf("solerr_after_refinement_double = %.10e\n", solerr_after_double);
    printf("solerr_after_refinement_posit  = %.10e\n", solerr_after_posit);
    printf("improvement_factor_double      = %.4f\n", solerr_before/solerr_after_double);
    printf("improvement_factor_posit       = %.4f\n", solerr_before/solerr_after_posit);
    return 0;
}
