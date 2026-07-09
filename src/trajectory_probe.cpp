// trajectory_probe.cpp
// Runs double64, float32, and posit32-naive CG in lockstep and logs, per iteration:
//   - alpha, beta for each variant
//   - relative L2 difference of p vector: float32 vs double64, posit32 vs double64, float32 vs posit32
// This directly shows WHEN and HOW SHARPLY the float32 and posit32-naive trajectories
// diverge from each other and from ground truth, rather than assuming compounding.
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
void matvec_f(const MTX& A, const std::vector<float>& x, std::vector<float>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=(float)A.val[k]*x[A.col[k]];
}
void matvec_p(const MTX& A, const std::vector<P32>& x, std::vector<P32>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]=y[A.row[k]]+P32(A.val[k])*x[A.col[k]];
}
double dot_d(const std::vector<double>& a, const std::vector<double>& b, int n){
    double s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
float dot_f(const std::vector<float>& a, const std::vector<float>& b, int n){
    float s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
double dot_p_naive(const std::vector<P32>& a, const std::vector<P32>& b, int n){
    P32 s=0; for(int i=0;i<n;i++) s=s+a[i]*b[i]; return double(s);
}

double relL2(const std::vector<double>& a, const std::vector<double>& b, int n){
    double num=0, den=0;
    for(int i=0;i<n;i++){ double d=a[i]-b[i]; num+=d*d; den+=a[i]*a[i]; }
    return (den>0) ? sqrt(num/den) : sqrt(num);
}

int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: trajectory_probe <matrix.mtx> <logfile> [maxiter=30]\n"); return 1; }
    int maxiter = argc>3 ? atoi(argv[3]) : 30;
    MTX A=read_mtx(argv[1]);
    int n=A.n;
    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    FILE* log=fopen(argv[2],"w");
    fprintf(log,"iter  alpha_d          alpha_f          alpha_p          relL2_p_f_vs_d   relL2_p_p_vs_d   relL2_p_f_vs_p\n");

    std::vector<double> xd(n,0),rd(n,1),pd(n),Apd(n),zd(n);
    for(int i=0;i<n;i++) zd[i]=rd[i]/diagA[i];
    for(int i=0;i<n;i++) pd[i]=zd[i];
    double rzd=dot_d(rd,zd,n);

    std::vector<float> xf(n,0),rf(n,1),pf(n),Apf(n),zf(n);
    for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i]);
    for(int i=0;i<n;i++) pf[i]=zf[i];
    float rzf=dot_f(rf,zf,n);

    std::vector<P32> xn(n,0),rn(n,1),pn(n),Apn(n),zn(n);
    for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i]);
    for(int i=0;i<n;i++) pn[i]=zn[i];
    double rzn=dot_p_naive(rn,zn,n);

    for(int iter=0; iter<maxiter; iter++){
        matvec_d(A,pd,Apd);
        double pApd=dot_d(pd,Apd,n);
        double alphad=rzd/pApd;
        for(int i=0;i<n;i++){xd[i]+=alphad*pd[i]; rd[i]-=alphad*Apd[i];}
        for(int i=0;i<n;i++) zd[i]=rd[i]/diagA[i];
        double rzd2=dot_d(rd,zd,n);
        double betad=rzd2/rzd; rzd=rzd2;
        for(int i=0;i<n;i++) pd[i]=zd[i]+betad*pd[i];

        matvec_f(A,pf,Apf);
        float pApf=dot_f(pf,Apf,n);
        float alphaf=rzf/pApf;
        for(int i=0;i<n;i++){xf[i]+=alphaf*pf[i]; rf[i]-=alphaf*Apf[i];}
        for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i])*rf[i];
        float rzf2=dot_f(rf,zf,n);
        float betaf=rzf2/rzf; rzf=rzf2;
        for(int i=0;i<n;i++) pf[i]=zf[i]+betaf*pf[i];

        matvec_p(A,pn,Apn);
        double pApn=dot_p_naive(pn,Apn,n);
        P32 alphan(rzn/pApn);
        for(int i=0;i<n;i++){xn[i]=xn[i]+alphan*pn[i]; rn[i]=rn[i]-alphan*Apn[i];}
        for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i])*rn[i];
        double rzn2=dot_p_naive(rn,zn,n);
        P32 betan(rzn2/rzn); rzn=rzn2;
        for(int i=0;i<n;i++) pn[i]=zn[i]+betan*pn[i];

        std::vector<double> pd_d(pd);
        std::vector<double> pf_d(n), pn_d(n);
        for(int i=0;i<n;i++){ pf_d[i]=(double)pf[i]; pn_d[i]=double(pn[i]); }

        double r_f_vs_d = relL2(pf_d, pd_d, n);
        double r_p_vs_d = relL2(pn_d, pd_d, n);
        double r_f_vs_p = relL2(pf_d, pn_d, n);

        fprintf(log,"%4d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
            iter,alphad,(double)alphaf,double(alphan), r_f_vs_d, r_p_vs_d, r_f_vs_p);
    }
    fclose(log);
    printf("Done. Log: %s\n",argv[2]);
}
