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
float dot_f_fma(const std::vector<float>& a, const std::vector<float>& b, int n){
    float s=0; for(int i=0;i<n;i++) s=std::fma(a[i],b[i],s); return s;
}
double dot_p(const std::vector<P32>& a, const std::vector<P32>& b, int n){
    quire<32,2,2> q=0;
    for(int i=0;i<n;i++) q+=quire_mul(a[i],b[i]);
    P32 r; convert(q.to_value(),r); return double(r);
}
double dot_range(const std::vector<P32>& a, const std::vector<P32>& b, int n){
    double vmin=1e300, vmax=0;
    for(int i=0;i<n;i++){
        double t = fabs(double(a[i])*double(b[i]));
        if(t>0 && t<vmin) vmin=t;
        if(t>vmax) vmax=t;
    }
    return (vmin>0) ? vmax/vmin : 0;
}
double dot_p_naive(const std::vector<P32>& a, const std::vector<P32>& b, int n){
    P32 s=0;
    for(int i=0;i<n;i++) s=s+a[i]*b[i];
    return double(s);
}

double sat_fraction(const std::vector<P32>& v, int n){
    static P32 maxp( SpecificValue::maxpos );
    static P32 minp( SpecificValue::minpos );
    int count = 0;
    for(int i=0;i<n;i++){
        P32 av = (v[i] < P32(0)) ? -v[i] : v[i];
        if(av == maxp || av == minp) count++;
    }
    return (double)count / n;
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
    fprintf(log,"iter  res_double64      res_float32      res_posit32q      res_posit32n      res_float32fma      pAp_dynrange      pAp_mag      sat_frac_pq      sat_frac_pn      pAp_double        pAp_naive\n");

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

    // --- float32 FMA CG ---
    std::vector<float> xg(n,0),rg(n,1),pg(n),Apg(n),zg(n);
    for(int i=0;i<n;i++) zg[i]=(float)(1.0/diagA[i]);
    for(int i=0;i<n;i++) pg[i]=zg[i];
    float rzg=dot_f_fma(rg,zg,n);

    // --- posit32+quire CG ---
    std::vector<P32> xp(n,0),rp(n,1),pp(n),App(n),zp(n);
    for(int i=0;i<n;i++) zp[i]=P32(1.0/diagA[i]);
    for(int i=0;i<n;i++) pp[i]=zp[i];
    double rzp=dot_p(rp,zp,n);

    // --- posit32 naive (no quire) CG ---
    std::vector<P32> xn(n,0),rn(n,1),pn(n),Apn(n),zn(n);
    for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i]);
    for(int i=0;i<n;i++) pn[i]=zn[i];
    double rzn=dot_p_naive(rn,zn,n);

    int maxiter=2000;
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

        // float32 FMA step
        matvec_f(A,pg,Apg);
        float pApg=dot_f_fma(pg,Apg,n);
        float alphag=rzg/pApg;
        for(int i=0;i<n;i++){xg[i]+=alphag*pg[i]; rg[i]-=alphag*Apg[i];}
        for(int i=0;i<n;i++) zg[i]=(float)(1.0/diagA[i])*rg[i];
        float rzg2=dot_f_fma(rg,zg,n);
        float betag=rzg2/rzg; rzg=rzg2;
        for(int i=0;i<n;i++) pg[i]=zg[i]+betag*pg[i];
        double resg=sqrt((double)dot_f_fma(rg,rg,n));

        // posit+quire step
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
        double pAp_range = dot_range(pp,App,n);
        double pAp_mag = fabs(pApp);
        double satp = sat_fraction(pp, n);

        // posit naive (no quire) step
        matvec_p(A,pn,Apn);
        double pApn=dot_p_naive(pn,Apn,n);
        P32 alphan(rzn/pApn);
        for(int i=0;i<n;i++){xn[i]=xn[i]+alphan*pn[i]; rn[i]=rn[i]-alphan*Apn[i];}
        for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i])*rn[i];
        double rzn2=dot_p_naive(rn,zn,n);
        P32 betan(rzn2/rzn); rzn=rzn2;
        for(int i=0;i<n;i++) pn[i]=zn[i]+betan*pn[i];
        std::vector<P32> rn_tmp(rn);
        double resn=sqrt(dot_p_naive(rn_tmp,rn_tmp,n));
        double satn = sat_fraction(pn, n);

        fprintf(log,"%4d  %.10e  %.10e  %.10e  %.10e  %.10e  %.6e  %.6e  %.6e  %.6e  %.10e  %.10e\n",iter,resd,resf,resp,resn,resg,pAp_range,pAp_mag,satp,satn,pApd,pApn);

        // stop if all converged
        if(resd<1e-10 && resf<1e-10 && resp<1e-10 && resn<1e-10 && resg<1e-10) break;
    }
    fclose(log);
    printf("Done. Log: %s\n",argv[2]);
}
