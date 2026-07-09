// hybrid_probe.cpp
// Runs THREE otherwise-identical double64 CG solvers, differing ONLY in how pAp
// (the denominator of alpha) is rounded each iteration:
//   (A) pure double64 pAp (ground truth control)
//   (B) pAp rounded through float32 precision (cast pAp to float, back to double)
//   (C) pAp rounded through posit32-naive precision (cast pAp to P32, back to double)
// Everything else (vectors, matvec, all other ops) stays double64 throughout.
// If (C)'s p-vector diverges from (A) faster than (B)'s does, that ISOLATES pAp
// rounding as the causal mechanism, rather than a coincidental correlation.
#include <cstdio>
#include <cmath>
#include <vector>
#include <universal/number/posit/posit.hpp>
using namespace sw::universal;
using P32 = posit<32,2>;

struct MTX { int n; std::vector<int> row,col; std::vector<double> val; };

MTX read_mtx(const char* f){
    MTX m; FILE* fp=fopen(f,"r");
    if(!fp){ printf("ERROR: cannot open %s\n",f); exit(1); }
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int rows,cols,nnz; sscanf(buf,"%d %d %d",&rows,&cols,&nnz);
    m.n=rows; int r,c; double v;
    while(fscanf(fp,"%d %d %lf",&r,&c,&v)==3){
        m.row.push_back(r-1); m.col.push_back(c-1); m.val.push_back(v);
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
double relL2(const std::vector<double>& a, const std::vector<double>& b, int n){
    double num=0, den=0;
    for(int i=0;i<n;i++){ double d=a[i]-b[i]; num+=d*d; den+=a[i]*a[i]; }
    return (den>0) ? sqrt(num/den) : sqrt(num);
}

double round_via_float32(double x){ return (double)(float)x; }
double round_via_posit32(double x){ P32 p(x); return double(p); }

void run_variant(const MTX& A, const std::vector<double>& diagA, int n, int maxiter,
                  double (*rounder)(double), std::vector<std::vector<double>>& p_history){
    std::vector<double> x(n,0),r(n,1),p(n),Ap(n),z(n);
    for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
    for(int i=0;i<n;i++) p[i]=z[i];
    double rz=dot_d(r,z,n);
    for(int iter=0; iter<maxiter; iter++){
        matvec_d(A,p,Ap);
        double pAp_exact = dot_d(p,Ap,n);
        double pAp = rounder ? rounder(pAp_exact) : pAp_exact;
        double alpha = rz/pAp;
        for(int i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }
        for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
        double rz2=dot_d(r,z,n);
        double beta=rz2/rz; rz=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
        p_history[iter]=p;
    }
}

int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: hybrid_probe <matrix.mtx> <logfile> [maxiter=30]\n"); return 1; }
    int maxiter = argc>3 ? atoi(argv[3]) : 30;
    MTX A=read_mtx(argv[1]);
    int n=A.n;
    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    std::vector<std::vector<double>> p_ctrl(maxiter, std::vector<double>(n));
    std::vector<std::vector<double>> p_f32(maxiter, std::vector<double>(n));
    std::vector<std::vector<double>> p_p32(maxiter, std::vector<double>(n));

    run_variant(A, diagA, n, maxiter, nullptr, p_ctrl);
    run_variant(A, diagA, n, maxiter, round_via_float32, p_f32);
    run_variant(A, diagA, n, maxiter, round_via_posit32, p_p32);

    FILE* log=fopen(argv[2],"w");
    fprintf(log,"iter  relL2_f32pAp_vs_ctrl   relL2_p32pAp_vs_ctrl   ratio(p32/f32)\n");
    for(int iter=0; iter<maxiter; iter++){
        double rf = relL2(p_f32[iter], p_ctrl[iter], n);
        double rp = relL2(p_p32[iter], p_ctrl[iter], n);
        double ratio = (rf>0) ? rp/rf : 0;
        fprintf(log,"%4d  %.6e            %.6e            %.3f\n", iter, rf, rp, ratio);
    }
    fclose(log);
    printf("Done. Log: %s\n",argv[2]);
    printf("If p32 column grows faster than f32 column, pAp rounding alone explains the divergence.\n");
}
