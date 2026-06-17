
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
void matvec(const MTX& A,const double* x,double* y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=A.val[k]*x[A.col[k]];
}
int main(){
    MTX A=read_mtx("/mnt/d/posits-work/bcsstk14/bcsstk14.mtx");
    int n=A.n;
    std::vector<double> x(n,0),r(n,1),p(n,1),Ap(n);
    double rr_d=0; for(int i=0;i<n;i++) rr_d+=r[i]*r[i];
    FILE* log=fopen("/mnt/d/posits-work/bcsstk14_dotx64.log","w");
    FILE* hist=fopen("/mnt/d/posits-work/bcsstk14_vals.log","w");
    for(auto v:A.val) fprintf(hist,"%.10e\n",v);
    fclose(hist);
    for(int iter=0;iter<50;iter++){
        matvec(A,p.data(),Ap.data());
        double pAp_d=0; for(int i=0;i<n;i++) pAp_d+=p[i]*Ap[i];
        quire<64,2,2> q=0;
        for(int i=0;i<n;i++){posit<64,2> pa(p[i]),pb(Ap[i]); q+=quire_mul(pa,pb);}
        posit<64,2> pr; convert(q.to_value(),pr); double pAp_p=double(pr);
        fprintf(log,"iter=%d pAp_d=%.10e pAp_p=%.10e diff=%.3e\n",iter,pAp_d,pAp_p,fabs(pAp_d-pAp_p));
        double alpha=rr_d/pAp_d;
        for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
        double rr2=0; for(int i=0;i<n;i++) rr2+=r[i]*r[i];
        double beta=rr2/rr_d; rr_d=rr2;
        for(int i=0;i<n;i++) p[i]=r[i]+beta*p[i];
        fprintf(log,"  residual=%.3e\n",sqrt(rr_d));
    }
    fclose(log); printf("Done\n");
}
