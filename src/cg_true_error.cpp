#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;
using P32 = posit<32,2>;

struct MTX { int n; std::vector<int> row,col; std::vector<double> val; std::vector<double> diag; };
MTX read_mtx(const char* f){
    MTX m; FILE* fp=fopen(f,"r");
    if(!fp){ printf("ERROR: cannot open %s\n",f); exit(1); }
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int rows,cols,nnz; sscanf(buf,"%d %d %d",&rows,&cols,&nnz);
    m.n=rows; m.diag.assign(rows,0.0);
    int r,c; double v;
    while(true){
        int nread = fscanf(fp,"%d %d %lf",&r,&c,&v);
        if(nread==2) v=1.0; // pattern format: no value column, implicit 1.0
        else if(nread!=3) break;
        m.row.push_back(r-1); m.col.push_back(c-1); m.val.push_back(v);
        if(r==c) m.diag[r-1]=v;
        if(r!=c){m.row.push_back(c-1); m.col.push_back(r-1); m.val.push_back(v);}
    } fclose(fp); return m;
}
void matvec_d(const MTX& A, const std::vector<double>& x, std::vector<double>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=A.val[k]*x[A.col[k]];
}
double dot_d(const std::vector<double>& a, const std::vector<double>& b, int n){
    double s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
double norm2(const std::vector<double>& a, int n){
    double s=0; for(int i=0;i<n;i++) s+=a[i]*a[i]; return sqrt(s);
}

// Jacobi-preconditioned CG in double precision, tracking TRUE error against
// a known x_true (not just residual). Per James Quinlan's suggestion:
// generate random x_true, set b = A*x_true, start x_0 = 0.
int main(int argc, char** argv){
    if(argc<3){ printf("usage: %s matrix.mtx seed [maxiter]\n", argv[0]); return 1; }
    MTX A = read_mtx(argv[1]);
    unsigned seed = atoi(argv[2]);
    int maxiter = argc>3 ? atoi(argv[3]) : 2000;
    int n = A.n;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.5,1.5);

    // ground-truth x_true, random
    std::vector<double> x_true(n);
    for(int i=0;i<n;i++) x_true[i]=dist(rng);
    double x_true_norm = norm2(x_true, n);

    // b = A * x_true
    std::vector<double> b(n);
    matvec_d(A, x_true, b);
    double b_norm = norm2(b, n);

    // Jacobi-preconditioned CG, x_0 = 0
    std::vector<double> x(n,0.0), r(n), z(n), p(n), Ap(n), diff(n);
    r = b; // r_0 = b - A*x_0 = b since x_0=0

    for(int i=0;i<n;i++) z[i] = r[i]/A.diag[i];
    p = z;
    double rz_old = dot_d(r,z,n);

    FILE* log = fopen(("results/true_error_logs/" + std::string(argv[1]).substr(std::string(argv[1]).find_last_of('/')+1) + "_seed" + std::to_string(seed) + ".log").c_str(), "w");
    fprintf(log, "iter rel_error rel_residual\n");

    for(int iter=0; iter<maxiter; iter++){
        matvec_d(A, p, Ap);
        double pAp = dot_d(p, Ap, n);
        double alpha = rz_old / pAp;
        for(int i=0;i<n;i++) x[i]+=alpha*p[i];
        for(int i=0;i<n;i++) r[i]-=alpha*Ap[i];

        for(int i=0;i<n;i++) diff[i]=x[i]-x_true[i];
        double rel_error = norm2(diff,n)/x_true_norm;
        double rel_residual = norm2(r,n)/b_norm;

        fprintf(log, "%d %.10e %.10e\n", iter, rel_error, rel_residual);

        if(rel_residual < 1e-10) break;

        for(int i=0;i<n;i++) z[i]=r[i]/A.diag[i];
        double rz_new = dot_d(r,z,n);
        double beta = rz_new/rz_old;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
        rz_old = rz_new;
    }
    fclose(log);
    printf("done: %s seed=%u\n", argv[1], seed);
    return 0;
}
