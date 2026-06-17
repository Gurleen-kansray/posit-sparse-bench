#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
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

// matvec for each precision
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

// dot products
double dot_d(const std::vector<double>& a, const std::vector<double>& b, int n){
    double s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
float dot_f(const std::vector<float>& a, const std::vector<float>& b, int n){
    float s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
double dot_p(const std::vector<P32>& a, const std::vector<P32>& b, int n){
    quire<32,2,2> q=0;
    for(int i=0;i<n;i++) q+=quire_mul(a[i],b[i]);
    P32 r; convert(q.to_value(),r); return double(r);
}

int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: cg_compare <matrix.mtx> <logfile>\n"); return 1; }
    MTX A=read_mtx(argv[1]);
    int n=A.n;

    // print matrix info
    double vmin=1e300,vmax=0;
    for(double v:A.val){double av=fabs(v);if(av>0&&av<vmin)vmin=av;if(av>vmax)vmax=av;}
    printf("Matrix: %s  n=%d  nnz=%zu  val_ratio=%.3e\n",argv[1],n,A.row.size(),vmax/vmin);

    // diagonal preconditioner
    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    FILE* log=fopen(argv[2],"w");
    fprintf(log,"iter  res_double64      res_float32      res_posit32q\n");

    // --- double64 CG ---
    std::vector<double> xd(n,0),rd(n,1),pd(n),Apd(n),zd(n);
    for(int i=0;i<n;i++) zd[i]=rd[i]/diagA[i];
    for(int i=0;i<n;i++) pd[i]=zd[i];
    double rzd=dot_d(rd,zd,n);

    // --- float32 CG ---
    std::vector<float> xf(n,0),rf(n,1),pf(n),Apf(n),zf(n);
    for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i]);
    for(int i=0;i<n;i++) pf[i]=zf[i];
    float rzf=dot_f(rf,zf,n);

    // --- posit32+quire CG ---
    std::vector<P32> xp(n,0),rp(n,1),pp(n),App(n),zp(n);
    for(int i=0;i<n;i++) zp[i]=P32(1.0/diagA[i]);
    for(int i=0;i<n;i++) pp[i]=zp[i];
    double rzp=dot_p(rp,zp,n);

    int maxiter=500;
    for(int iter=0;iter<maxiter;iter++){
        // double step
        matvec_d(A,pd,Apd);
        double pApd=dot_d(pd,Apd,n);
        double alphad=rzd/pApd;
        for(int i=0;i<n;i++){xd[i]+=alphad*pd[i]; rd[i]-=alphad*Apd[i];}
        for(int i=0;i<n;i++) zd[i]=rd[i]/diagA[i];
        double rzd2=dot_d(rd,zd,n);
        double betad=rzd2/rzd; rzd=rzd2;
        for(int i=0;i<n;i++) pd[i]=zd[i]+betad*pd[i];
        double resd=sqrt(dot_d(rd,rd,n));

        // float step
        matvec_f(A,pf,Apf);
        float pApf=dot_f(pf,Apf,n);
        float alphaf=rzf/pApf;
        for(int i=0;i<n;i++){xf[i]+=alphaf*pf[i]; rf[i]-=alphaf*Apf[i];}
        for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i])*rf[i];
        float rzf2=dot_f(rf,zf,n);
        float betaf=rzf2/rzf; rzf=rzf2;
        for(int i=0;i<n;i++) pf[i]=zf[i]+betaf*pf[i];
        double resf=sqrt((double)dot_f(rf,rf,n));

        // posit step
        matvec_p(A,pp,App);
        double pApp=dot_p(pp,App,n);
        P32 alphap(rzp/pApp);
        for(int i=0;i<n;i++){xp[i]=xp[i]+alphap*pp[i]; rp[i]=rp[i]-alphap*App[i];}
        for(int i=0;i<n;i++) zp[i]=P32(1.0/diagA[i])*rp[i];
        double rzp2=dot_p(rp,zp,n);
        P32 betap(rzp2/rzp); rzp=rzp2;
        for(int i=0;i<n;i++) pp[i]=zp[i]+betap*pp[i];
        std::vector<P32> rp_tmp(rp);
        double resp=sqrt(dot_p(rp_tmp,rp_tmp,n));

        fprintf(log,"%4d  %.6e  %.6e  %.6e\n",iter,resd,resf,resp);

        // stop if all converged
        if(resd<1e-10 && resf<1e-10 && resp<1e-10) break;
    }
    fclose(log);
    printf("Done. Log: %s\n",argv[2]);
}
