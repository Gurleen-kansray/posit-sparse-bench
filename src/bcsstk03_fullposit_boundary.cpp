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

// double matvec
void matvec_double(const MTX& A, const std::vector<double>& val_scaled,
                   const double* x, double* y, int n){
    for(int i=0;i<n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++)
        y[A.row[k]] += val_scaled[k] * x[A.col[k]];
}

// posit32 matvec - converts entries and vector to posit32, accumulates in posit32
void matvec_posit32(const MTX& A, const std::vector<double>& val_scaled,
                    const double* x, double* y, int n){
    for(int i=0;i<n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++){
        posit<32,2> pa(val_scaled[k]);
        posit<32,2> pb(x[A.col[k]]);
        posit<32,2> py(y[A.row[k]]);
        py = py + pa*pb;
        y[A.row[k]] = double(py);
    }
}

struct CGResult { int iters; double final_residual; bool converged; };

// method: 0=double/double, 1=double-mv/quire-dot, 2=double-mv/naive-dot
//         3=posit32-mv/quire-dot, 4=posit32-mv/naive-dot
CGResult run_cg(const MTX& A, const std::vector<double>& val_scaled,
                const std::vector<double>& diagScaled, int method,
                int max_iter, double tol){
    int n = A.n;
    std::vector<double> x(n,0), r(n,1), p(n), Ap(n), z(n);
    for(int i=0;i<n;i++) z[i] = diagScaled[i]>0 ? r[i]/diagScaled[i] : r[i];
    for(int i=0;i<n;i++) p[i] = z[i];
    double rz = 0; for(int i=0;i<n;i++) rz += r[i]*z[i];

    bool use_posit_mv = (method==3 || method==4);

    for(int iter=0; iter<max_iter; iter++){
        if(use_posit_mv)
            matvec_posit32(A, val_scaled, p.data(), Ap.data(), n);
        else
            matvec_double(A, val_scaled, p.data(), Ap.data(), n);

        double pAp = 0;
        if(method==0 || method==1 || method==3){
            // quire dot (or pure double)
            if(method==0){
                for(int i=0;i<n;i++) pAp += p[i]*Ap[i];
            } else {
                quire<32,2,2> q=0;
                for(int i=0;i<n;i++){
                    posit<32,2> pa(p[i]), pb(Ap[i]);
                    q += quire_mul(pa,pb);
                }
                posit<32,2> pr; convert(q.to_value(),pr); pAp=double(pr);
            }
        } else {
            // naive dot (method 2 or 4)
            posit<32,2> s=0;
            for(int i=0;i<n;i++){
                posit<32,2> pa(p[i]), pb(Ap[i]);
                s = s + pa*pb;
            }
            pAp = double(s);
        }

        if(!std::isfinite(pAp) || pAp <= 0){
            double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
            return {iter, sqrt(rr), false};
        }

        double alpha = rz/pAp;
        for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
        for(int i=0;i<n;i++) z[i] = diagScaled[i]>0 ? r[i]/diagScaled[i] : r[i];
        double rz2=0; for(int i=0;i<n;i++) rz2+=r[i]*z[i];

        double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
        double res=sqrt(rr);

        if(!std::isfinite(res) || res > 1e30)
            return {iter, res, false};
        if(res < tol)
            return {iter+1, res, true};

        double beta = rz2/rz; rz=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
    }
    double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
    return {max_iter, sqrt(rr), false};
}

int main(){
    MTX A = read_mtx("/mnt/d/posits-work/bcsstk03/bcsstk03.mtx");
    int n = A.n;

    std::vector<double> scales = {1.0, 1e2, 1e4, 1e6, 1e8, 1e9, 1e10, 1e11, 1e12};

    std::vector<double> diagA(n,0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];
    double diag_max=0;
    for(int i=0;i<n;i++) diag_max=std::max(diag_max,diagA[i]);

    FILE* log = fopen("/mnt/d/posits-work/bcsstk03_fullposit_boundary.log","w");

    const char* hdr = "%-10s %-12s %-10s %-14s %-14s %-14s %-14s\n";
    const char* row = "%-10.0e %-12.2e %-10s %-14s %-14s %-14s %-14s\n";
    printf(hdr,  "scale","diag_range","dbl/dbl","dbl-mv/quire","dbl-mv/naive","p32-mv/quire","p32-mv/naive");
    fprintf(log, hdr, "scale","diag_range","dbl/dbl","dbl-mv/quire","dbl-mv/naive","p32-mv/quire","p32-mv/naive");

    const int MAX_ITER=600;
    const double TOL=1e-6;

    for(double S : scales){
        std::vector<double> val_scaled = A.val;
        std::vector<double> diag_scaled = diagA;
        double threshold = diag_max * 0.1;
        for(int k=0;k<(int)A.row.size();k++){
            if(A.row[k]==A.col[k] && A.val[k] >= threshold){
                val_scaled[k] *= S;
                diag_scaled[A.row[k]] = val_scaled[k];
            }
        }
        double dmin=1e300, dmax=0;
        for(int i=0;i<n;i++){
            if(diag_scaled[i]>0){
                dmin=std::min(dmin,diag_scaled[i]);
                dmax=std::max(dmax,diag_scaled[i]);
            }
        }
        double actual_range = dmax/dmin;

        auto r0 = run_cg(A,val_scaled,diag_scaled,0,MAX_ITER,TOL);
        auto r1 = run_cg(A,val_scaled,diag_scaled,1,MAX_ITER,TOL);
        auto r2 = run_cg(A,val_scaled,diag_scaled,2,MAX_ITER,TOL);
        auto r3 = run_cg(A,val_scaled,diag_scaled,3,MAX_ITER,TOL);
        auto r4 = run_cg(A,val_scaled,diag_scaled,4,MAX_ITER,TOL);

        auto fmt = [](CGResult& r) -> std::string {
            return r.converged ? std::to_string(r.iters) : "DIV";
        };

        printf(row, S, actual_range,
               fmt(r0).c_str(), fmt(r1).c_str(), fmt(r2).c_str(),
               fmt(r3).c_str(), fmt(r4).c_str());
        fprintf(log, row, S, actual_range,
               fmt(r0).c_str(), fmt(r1).c_str(), fmt(r2).c_str(),
               fmt(r3).c_str(), fmt(r4).c_str());
        fflush(stdout); fflush(log);
    }

    fclose(log);
    printf("Done.\n");
    return 0;
}
