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

// p32-mv/naive, prints residual every iteration
void run_diagnostic(const MTX& A, const std::vector<double>& val_scaled,
                    const std::vector<double>& diagScaled,
                    double S, FILE* log){
    int n = A.n;
    std::vector<double> x(n,0), r(n,1), p(n), Ap(n), z(n);
    for(int i=0;i<n;i++) z[i] = diagScaled[i]>0 ? r[i]/diagScaled[i] : r[i];
    for(int i=0;i<n;i++) p[i] = z[i];
    double rz = 0; for(int i=0;i<n;i++) rz += r[i]*z[i];

    fprintf(log, "\n=== scale=%.0e (p32-mv/naive) ===\n", S);
    printf(      "\n=== scale=%.0e (p32-mv/naive) ===\n", S);
    fprintf(log, "%-6s %-14s %-10s\n", "iter", "residual", "status");
    printf(      "%-6s %-14s %-10s\n", "iter", "residual", "status");

    const double TOL = 1e-6;
    const int MAX_ITER = 600;

    for(int iter=0; iter<MAX_ITER; iter++){
        matvec_posit32(A, val_scaled, p.data(), Ap.data(), n);

        // naive dot
        posit<32,2> s=0;
        for(int i=0;i<n;i++){
            posit<32,2> pa(p[i]), pb(Ap[i]);
            s = s + pa*pb;
        }
        double pAp = double(s);

        const char* status = "";
        if(!std::isfinite(pAp) || pAp <= 0){
            double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
            fprintf(log, "%-6d %-14.6e %-10s  <- pAp=%.3e BREAKDOWN\n",
                    iter, sqrt(rr), "DIV", pAp);
            printf(      "%-6d %-14.6e %-10s  <- pAp=%.3e BREAKDOWN\n",
                    iter, sqrt(rr), "DIV", pAp);
            return;
        }

        double alpha = rz/pAp;
        for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
        for(int i=0;i<n;i++) z[i] = diagScaled[i]>0 ? r[i]/diagScaled[i] : r[i];
        double rz2=0; for(int i=0;i<n;i++) rz2+=r[i]*z[i];

        double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
        double res = sqrt(rr);

        // print every iter but condense once stable
        bool print_this = (iter < 20) || (iter % 50 == 0) ||
                          (res < TOL*100) || !std::isfinite(res) || res > 1e20;

        if(print_this){
            if(!std::isfinite(res) || res > 1e20) status = "<- DIVERGING";
            else if(res < TOL) status = "<- CONVERGED";
            fprintf(log, "%-6d %-14.6e %s\n", iter, res, status);
            printf(      "%-6d %-14.6e %s\n", iter, res, status);
        }

        if(!std::isfinite(res) || res > 1e30) return;
        if(res < TOL){
            fprintf(log, "  Converged at iter %d, residual %.3e\n", iter+1, res);
            printf(      "  Converged at iter %d, residual %.3e\n", iter+1, res);
            return;
        }

        double beta = rz2/rz; rz=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
    }
    double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
    fprintf(log, "  Did not converge in %d iters, final residual %.3e\n", MAX_ITER, sqrt(rr));
    printf(      "  Did not converge in %d iters, final residual %.3e\n", MAX_ITER, sqrt(rr));
}

int main(){
    MTX A = read_mtx("/mnt/d/posits-work/bcsstk03/bcsstk03.mtx");
    int n = A.n;

    std::vector<double> diagA(n,0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];
    double diag_max=0;
    for(int i=0;i<n;i++) diag_max=std::max(diag_max,diagA[i]);

    FILE* log = fopen("/mnt/d/posits-work/bcsstk03_anomaly_diag.log","w");

    // Only the three scales around the anomaly
    for(double S : {1e8, 1e9, 1e10}){
        std::vector<double> val_scaled = A.val;
        std::vector<double> diag_scaled = diagA;
        double threshold = diag_max * 0.1;
        for(int k=0;k<(int)A.row.size();k++){
            if(A.row[k]==A.col[k] && A.val[k] >= threshold){
                val_scaled[k] *= S;
                diag_scaled[A.row[k]] = val_scaled[k];
            }
        }
        run_diagnostic(A, val_scaled, diag_scaled, S, log);
    }

    fclose(log);
    printf("\nDone. Full log: /mnt/d/posits-work/bcsstk03_anomaly_diag.log\n");
    return 0;
}
