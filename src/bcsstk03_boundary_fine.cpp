#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
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

void matvec(const MTX& A, const std::vector<double>& vs,
            const double* x, double* y, int n){
    for(int i=0;i<n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=vs[k]*x[A.col[k]];
}

struct CGResult { int iters; double final_residual; bool converged; };

CGResult run_cg(const MTX& A, const std::vector<double>& vs,
                const std::vector<double>& diag, int method,
                int max_iter, double tol){
    int n=A.n;
    std::vector<double> x(n,0),r(n,1),p(n),Ap(n),z(n);
    for(int i=0;i<n;i++) z[i]=diag[i]>0?r[i]/diag[i]:r[i];
    for(int i=0;i<n;i++) p[i]=z[i];
    double rz=0; for(int i=0;i<n;i++) rz+=r[i]*z[i];
    for(int iter=0;iter<max_iter;iter++){
        matvec(A,vs,p.data(),Ap.data(),n);
        double pAp=0;
        if(method==0){
            for(int i=0;i<n;i++) pAp+=p[i]*Ap[i];
        } else if(method==1){
            quire<32,2,2> q=0;
            for(int i=0;i<n;i++){
                posit<32,2> pa(p[i]),pb(Ap[i]); q+=quire_mul(pa,pb);
            }
            posit<32,2> pr; convert(q.to_value(),pr); pAp=double(pr);
        } else {
            posit<32,2> pn=0;
            for(int i=0;i<n;i++){
                posit<32,2> pa(p[i]),pb(Ap[i]); pn=pn+pa*pb;
            }
            pAp=double(pn);
        }
        if(!std::isfinite(pAp)||pAp<=0) return {iter,1e99,false};
        double alpha=rz/pAp;
        for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
        for(int i=0;i<n;i++) z[i]=diag[i]>0?r[i]/diag[i]:r[i];
        double rz2=0; for(int i=0;i<n;i++) rz2+=r[i]*z[i];
        double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
        if(sqrt(rr)<tol){ return {iter,sqrt(rr),true}; }
        double beta=rz2/rz; rz=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
    }
    double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
    return {max_iter,sqrt(rr),false};
}

int main(){
    MTX A=read_mtx("data/matrices/bcsstk03.mtx");
    int n=A.n;
    std::vector<double> diagA(n,0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];
    double diag_max=0;
    for(int i=0;i<n;i++) diag_max=std::max(diag_max,diagA[i]);

    // Fine sweep: 1e8 to 1e10 in half-decade steps
    // This is where naive diverges — we want to pinpoint exactly where
    std::vector<double> scales = {
        1e8, 2e8, 5e8,
        1e9, 2e9, 5e9,
        1e10
    };

    printf("Fine-grained divergence boundary sweep (bcsstk03)\n");
    printf("%-12s %-14s %-10s %-10s %-10s\n",
           "scale","diag_range","dbl","quire","naive");
    printf("%-12s %-14s %-10s %-10s %-10s\n",
           "------","----------","---","-----","-----");

    FILE* log=fopen("results/logs/bcsstk03_boundary_fine.log","w");
    fprintf(log,"scale diag_range dbl_iters quire_iters naive_iters dbl_res quire_res naive_res\n");

    for(double S : scales){
        std::vector<double> vs=A.val;
        std::vector<double> ds=diagA;
        double threshold=diag_max*0.1;
        for(int k=0;k<(int)A.row.size();k++){
            if(A.row[k]==A.col[k] && A.val[k]>=threshold){
                vs[k]*=S; ds[A.row[k]]=vs[k];
            }
        }
        double dmin=1e300,dmax=0;
        for(int i=0;i<n;i++){
            if(ds[i]>0){dmin=std::min(dmin,ds[i]); dmax=std::max(dmax,ds[i]);}
        }
        double range=dmax/dmin;
        auto rd=run_cg(A,vs,ds,0,500,1e-6);
        auto rq=run_cg(A,vs,ds,1,500,1e-6);
        auto rn=run_cg(A,vs,ds,2,500,1e-6);
        auto fi=[](CGResult& r)->std::string{
            return r.converged?std::to_string(r.iters):"DIV";
        };
        printf("%-12.1e %-14.2e %-10s %-10s %-10s\n",
               S,range,fi(rd).c_str(),fi(rq).c_str(),fi(rn).c_str());
        fprintf(log,"%.1e %.2e %s %s %s %.3e %.3e %.3e\n",
               S,range,fi(rd).c_str(),fi(rq).c_str(),fi(rn).c_str(),
               rd.final_residual,rq.final_residual,rn.final_residual);
        fflush(stdout); fflush(log);
    }
    fclose(log);
    printf("\nLog saved: results/logs/bcsstk03_boundary_fine.log\n");
}
