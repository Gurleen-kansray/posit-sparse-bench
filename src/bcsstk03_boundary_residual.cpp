#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;
struct MTX { int n; std::vector<int> row,col; std::vector<double> val; };
MTX read_mtx(const char* f){
    MTX m; FILE* fp=fopen(f,"r");
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int rows,cols,nnz; sscanf(buf,"%d %d %d",&rows,&cols,&nnz);
    m.n=rows; int r,c; double v;
    while(fscanf(fp,"%d %d %lf",&r,&c,&v)==3){
        m.row.push_back(r-1); m.col.push_back(c-1); m.val.push_back(v);
        if(r!=c){m.row.push_back(c-1); m.col.push_back(r-1); m.val.push_back(v);}
    } fclose(fp); return m;
}
void matvec(const MTX& A,const std::vector<double>& vs,const double* x,double* y,int n){
    for(int i=0;i<n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=vs[k]*x[A.col[k]];
}
int main(){
    MTX A=read_mtx("data/matrices/bcsstk03.mtx");
    int n=A.n;
    std::vector<double> diagA(n,0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];
    double diag_max=0;
    for(int i=0;i<n;i++) diag_max=std::max(diag_max,diagA[i]);

    // Test these three scales — one clean converge, one ambiguous, one clear DIV
    std::vector<double> scales={1e8, 2e8, 5e8};

    for(double S : scales){
        std::vector<double> vs=A.val;
        std::vector<double> ds=diagA;
        for(int k=0;k<(int)A.row.size();k++){
            if(A.row[k]==A.col[k] && A.val[k]>=diag_max*0.1){
                vs[k]*=S; ds[A.row[k]]=vs[k];
            }
        }
        // Run naive only, print every 10th residual
        printf("\n=== scale=%.0e ===  (naive posit32)\n",S);
        printf("iter  residual\n");
        std::vector<double> x(n,0),r(n,1),p(n),Ap(n),z(n);
        for(int i=0;i<n;i++) z[i]=ds[i]>0?r[i]/ds[i]:r[i];
        for(int i=0;i<n;i++) p[i]=z[i];
        double rz=0; for(int i=0;i<n;i++) rz+=r[i]*z[i];
        for(int iter=0;iter<500;iter++){
            matvec(A,vs,p.data(),Ap.data(),n);
            posit<32,2> pn=0;
            for(int i=0;i<n;i++){posit<32,2> pa(p[i]),pb(Ap[i]); pn=pn+pa*pb;}
            double pAp=double(pn);
            if(!std::isfinite(pAp)||pAp<=0){printf("%4d  BREAKDOWN\n",iter); break;}
            double alpha=rz/pAp;
            for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
            for(int i=0;i<n;i++) z[i]=ds[i]>0?r[i]/ds[i]:r[i];
            double rz2=0; for(int i=0;i<n;i++) rz2+=r[i]*z[i];
            double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
            if(iter%20==0||sqrt(rr)<1e-6) printf("%4d  %.3e\n",iter,sqrt(rr));
            if(sqrt(rr)<1e-6){printf("CONVERGED at %d\n",iter); break;}
            double beta=rz2/rz; rz=rz2;
            for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
        }
    }
}
