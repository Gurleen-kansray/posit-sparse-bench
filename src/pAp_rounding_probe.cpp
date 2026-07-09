// pAp_rounding_probe.cpp
// Runs the pure double64 CG (ground truth), and at each iteration compares the
// SINGLE-VALUE rounding error of pAp when cast through posit32 vs float32.
// This isolates precisely how much rounding error each format introduces on
// this one scalar, at its actual observed magnitude each iteration.
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

int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: pAp_rounding_probe <matrix.mtx> <logfile> [maxiter=30]\n"); return 1; }
    int maxiter = argc>3 ? atoi(argv[3]) : 30;
    MTX A=read_mtx(argv[1]);
    int n=A.n;
    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    FILE* log=fopen(argv[2],"w");
    fprintf(log,"iter  pAp_true         relerr_float32     relerr_posit32     ratio(p32/f32)\n");

    std::vector<double> x(n,0),r(n,1),p(n),Ap(n),z(n);
    for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
    for(int i=0;i<n;i++) p[i]=z[i];
    double rz=dot_d(r,z,n);

    for(int iter=0; iter<maxiter; iter++){
        matvec_d(A,p,Ap);
        double pAp_true = dot_d(p,Ap,n);

        double pAp_f32 = (double)(float)pAp_true;
        double pAp_p32 = double(P32(pAp_true));

        double relerr_f32 = fabs(pAp_f32 - pAp_true) / fabs(pAp_true);
        double relerr_p32 = fabs(pAp_p32 - pAp_true) / fabs(pAp_true);
        double ratio = (relerr_f32>0) ? relerr_p32/relerr_f32 : 0;

        fprintf(log,"%4d  %.6e  %.6e      %.6e      %.3f\n",
            iter, pAp_true, relerr_f32, relerr_p32, ratio);

        double alpha = rz/pAp_true;
        for(int i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }
        for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
        double rz2=dot_d(r,z,n);
        double beta=rz2/rz; rz=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
    }
    fclose(log);
    printf("Done. Log: %s\n",argv[2]);
}
