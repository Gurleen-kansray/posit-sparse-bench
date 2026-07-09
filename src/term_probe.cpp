// term_probe.cpp
// Dumps individual p[i]*Ap[i] term magnitudes at specified iterations
// for float32 and posit32-naive CG, to compare against the posit32/float32
// precision crossover zone (~3.16e-5 to ~1e5, from posit_precision_curve).
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
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
void matvec_f(const MTX& A, const std::vector<float>& x, std::vector<float>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=(float)A.val[k]*x[A.col[k]];
}
void matvec_p(const MTX& A, const std::vector<P32>& x, std::vector<P32>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]=y[A.row[k]]+P32(A.val[k])*x[A.col[k]];
}
float dot_f(const std::vector<float>& a, const std::vector<float>& b, int n){
    float s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
double dot_p_naive(const std::vector<P32>& a, const std::vector<P32>& b, int n){
    P32 s=0; for(int i=0;i<n;i++) s=s+a[i]*b[i]; return double(s);
}

int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: term_probe <matrix.mtx> <logfile> [iter_start=15] [iter_end=25]\n"); return 1; }
    int iter_start = argc>3 ? atoi(argv[3]) : 15;
    int iter_end   = argc>4 ? atoi(argv[4]) : 25;
    MTX A=read_mtx(argv[1]);
    int n=A.n;
    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    FILE* log=fopen(argv[2],"w");
    fprintf(log,"iter  min_term_f  max_term_f  min_term_p  max_term_p  pct_terms_outside_posit_zone  pAp_f  pAp_p\n");

    std::vector<float> xf(n,0),rf(n,1),pf(n),Apf(n),zf(n);
    for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i]);
    for(int i=0;i<n;i++) pf[i]=zf[i];
    float rzf=dot_f(rf,zf,n);

    std::vector<P32> xn(n,0),rn(n,1),pn(n),Apn(n),zn(n);
    for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i]);
    for(int i=0;i<n;i++) pn[i]=zn[i];
    double rzn=dot_p_naive(rn,zn,n);

    const double ZONE_LO=3.16e-5, ZONE_HI=1e5;

    for(int iter=0; iter<=iter_end; iter++){
        // float step
        matvec_f(A,pf,Apf);
        float pApf=dot_f(pf,Apf,n);
        float alphaf=rzf/pApf;
        for(int i=0;i<n;i++){xf[i]+=alphaf*pf[i]; rf[i]-=alphaf*Apf[i];}
        for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i])*rf[i];
        float rzf2=dot_f(rf,zf,n);
        float betaf=rzf2/rzf; rzf=rzf2;
        for(int i=0;i<n;i++) pf[i]=zf[i]+betaf*pf[i];

        // posit naive step
        matvec_p(A,pn,Apn);
        double pApn=dot_p_naive(pn,Apn,n);
        P32 alphan(rzn/pApn);
        for(int i=0;i<n;i++){xn[i]=xn[i]+alphan*pn[i]; rn[i]=rn[i]-alphan*Apn[i];}
        for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i])*rn[i];
        double rzn2=dot_p_naive(rn,zn,n);
        P32 betan(rzn2/rzn); rzn=rzn2;
        for(int i=0;i<n;i++) pn[i]=zn[i]+betan*pn[i];

        if(iter>=iter_start && iter<=iter_end){
            double minf=1e300,maxf=0,minp=1e300,maxp=0;
            int outside=0;
            for(int i=0;i<n;i++){
                double tf=fabs((double)pf[i]*(double)Apf[i]);
                double tp=fabs(double(pn[i])*double(Apn[i]));
                if(tf>0){ minf=std::min(minf,tf); maxf=std::max(maxf,tf); }
                if(tp>0){ minp=std::min(minp,tp); maxp=std::max(maxp,tp);
                    if(tp<ZONE_LO || tp>ZONE_HI) outside++;
                }
            }
            double pct_outside = 100.0*outside/n;
            fprintf(log,"%4d  %.6e  %.6e  %.6e  %.6e  %.2f  %.6e  %.6e\n",
                iter,minf,maxf,minp,maxp,pct_outside,(double)pApf,pApn);
        }
    }
    fclose(log);
    printf("Done. Log: %s\n",argv[2]);
}
