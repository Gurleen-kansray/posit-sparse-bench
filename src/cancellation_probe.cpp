// cancellation_probe.cpp
// For target iterations, dumps the accumulation trajectory of the pAp dot product
// (partial sum after each term) for float32 and posit32-naive, plus running max |term|,
// to directly show catastrophic cancellation: large terms of both signs summing to a
// much smaller result, and how rounding error grows relative to the shrinking partial sum.
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
double dot_exact(const std::vector<double>& a, const std::vector<double>& b, int n){
    double s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}

int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: cancellation_probe <matrix.mtx> <logfile> <target_iter>\n"); return 1; }
    int target_iter = atoi(argv[3]);
    MTX A=read_mtx(argv[1]);
    int n=A.n;
    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    std::vector<float> xf(n,0),rf(n,1),pf(n),Apf(n),zf(n);
    for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i]);
    for(int i=0;i<n;i++) pf[i]=zf[i];
    float rzf=dot_f(rf,zf,n);

    std::vector<P32> xn(n,0),rn(n,1),pn(n),Apn(n),zn(n);
    for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i]);
    for(int i=0;i<n;i++) pn[i]=zn[i];
    double rzn=dot_p_naive(rn,zn,n);

    for(int iter=0; iter<=target_iter; iter++){
        matvec_f(A,pf,Apf);
        float pApf;
        if(iter==target_iter){
            std::vector<double> terms(n);
            for(int i=0;i<n;i++) terms[i]=(double)pf[i]*(double)Apf[i];
            double posmax=0, negmax=0;
            for(int i=0;i<n;i++){ if(terms[i]>posmax) posmax=terms[i]; if(terms[i]<negmax) negmax=terms[i]; }
            double exact_sum = dot_exact(std::vector<double>(pf.begin(),pf.end()),
                                          std::vector<double>(Apf.begin(),Apf.end()), n);
            float sf=0; double max_partial_f=0;
            for(int i=0;i<n;i++){ sf+=pf[i]*Apf[i]; if(fabs((double)sf)>max_partial_f) max_partial_f=fabs((double)sf); }
            printf("=== iter %d (float32) ===\n", iter);
            printf("largest positive term: %.6e\n", posmax);
            printf("largest negative term: %.6e\n", negmax);
            printf("exact (double) sum:     %.6e\n", exact_sum);
            printf("float32 accumulated:    %.6e\n", (double)sf);
            printf("max |partial sum| during accumulation: %.6e\n", max_partial_f);
            printf("cancellation ratio (max|partial|/|final|): %.3e\n", max_partial_f/fabs(exact_sum));
        }
        pApf=dot_f(pf,Apf,n);
        float alphaf=rzf/pApf;
        for(int i=0;i<n;i++){xf[i]+=alphaf*pf[i]; rf[i]-=alphaf*Apf[i];}
        for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i])*rf[i];
        float rzf2=dot_f(rf,zf,n);
        float betaf=rzf2/rzf; rzf=rzf2;
        for(int i=0;i<n;i++) pf[i]=zf[i]+betaf*pf[i];

        matvec_p(A,pn,Apn);
        double pApn;
        if(iter==target_iter){
            std::vector<double> terms(n);
            for(int i=0;i<n;i++) terms[i]=double(pn[i])*double(Apn[i]);
            double posmax=0, negmax=0;
            for(int i=0;i<n;i++){ if(terms[i]>posmax) posmax=terms[i]; if(terms[i]<negmax) negmax=terms[i]; }
            std::vector<double> pnd(n), Apnd(n);
            for(int i=0;i<n;i++){ pnd[i]=double(pn[i]); Apnd[i]=double(Apn[i]); }
            double exact_sum = dot_exact(pnd, Apnd, n);
            P32 sp=0; double max_partial_p=0;
            for(int i=0;i<n;i++){ sp=sp+pn[i]*Apn[i]; if(fabs(double(sp))>max_partial_p) max_partial_p=fabs(double(sp)); }
            printf("=== iter %d (posit32 naive) ===\n", iter);
            printf("largest positive term: %.6e\n", posmax);
            printf("largest negative term: %.6e\n", negmax);
            printf("exact (double) sum:     %.6e\n", exact_sum);
            printf("posit32 accumulated:    %.6e\n", double(sp));
            printf("max |partial sum| during accumulation: %.6e\n", max_partial_p);
            printf("cancellation ratio (max|partial|/|final|): %.3e\n", max_partial_p/fabs(exact_sum));
        }
        pApn=dot_p_naive(pn,Apn,n);
        P32 alphan(rzn/pApn);
        for(int i=0;i<n;i++){xn[i]=xn[i]+alphan*pn[i]; rn[i]=rn[i]-alphan*Apn[i];}
        for(int i=0;i<n;i++) zn[i]=P32(1.0/diagA[i])*rn[i];
        double rzn2=dot_p_naive(rn,zn,n);
        P32 betan(rzn2/rzn); rzn=rzn2;
        for(int i=0;i<n;i++) pn[i]=zn[i]+betan*pn[i];
    }
    printf("Done.\n");
}
