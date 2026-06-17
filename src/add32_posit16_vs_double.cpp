#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

struct MTX { 
    int rows, cols; 
    std::vector<int> row, col; 
    std::vector<double> val; 
};

MTX read_mtx(const char* f, bool symmetric_expand){
    MTX m; FILE* fp=fopen(f,"r");
    if(!fp){ printf("Cannot open %s\n",f); exit(1); }
    char buf[256];
    bool is_symmetric = false;
    while(fgets(buf,256,fp)){
        if(buf[0]=='%'){
            if(strstr(buf,"symmetric")) is_symmetric=true;
            continue;
        }
        break;
    }
    int rows,cols,nnz; sscanf(buf,"%d %d %d",&rows,&cols,&nnz);
    m.rows=rows; m.cols=cols;
    int r,c; double v;
    while(fscanf(fp,"%d %d %lf",&r,&c,&v)==3){
        m.row.push_back(r-1); m.col.push_back(c-1); m.val.push_back(v);
        if(symmetric_expand && is_symmetric && r!=c){
            m.row.push_back(c-1); m.col.push_back(r-1); m.val.push_back(v);
        }
    }
    fclose(fp); return m;
}

std::vector<double> read_rhs(const char* f, int n){
    std::vector<double> b(n, 1.0);
    FILE* fp = fopen(f,"r");
    if(!fp){ printf("No rhs file, using b=ones\n"); return b; }
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int rows,cols; sscanf(buf,"%d %d",&rows,&cols);
    for(int i=0;i<rows && i<n;i++){
        int r; double v;
        if(fscanf(fp,"%d %lf",&r,&v)==2) b[r-1]=v;
        else { double vv; if(fscanf(fp,"%lf",&vv)==1) b[i]=vv; }
    }
    fclose(fp);
    printf("Loaded rhs from file, n=%d\n", n);
    return b;
}

void matvec(const MTX& A, const double* x, double* y){
    for(int i=0;i<A.rows;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++)
        y[A.row[k]] += A.val[k]*x[A.col[k]];
}

struct CGResult { int iters; double final_residual; bool converged; double b_norm; };

// method: 0=double, 1=posit32+quire, 2=posit16+quire, 3=posit16+naive
CGResult run_cg(const MTX& A, const std::vector<double>& b,
                const std::vector<double>& diagA,
                int method, int max_iter, double tol_rel,
                FILE* trace, const char* label){
    int n = A.rows;
    std::vector<double> x(n,0), r(n), p(n), Ap(n), z(n);
    
    // r = b - A*x = b (since x=0)
    for(int i=0;i<n;i++) r[i]=b[i];
    
    double b_norm=0; for(int i=0;i<n;i++) b_norm+=b[i]*b[i]; b_norm=sqrt(b_norm);
    double tol = tol_rel * b_norm;
    
    // diagonal preconditioner
    for(int i=0;i<n;i++) z[i] = diagA[i]>0 ? r[i]/diagA[i] : r[i];
    for(int i=0;i<n;i++) p[i] = z[i];
    double rz=0; for(int i=0;i<n;i++) rz+=r[i]*z[i];

    if(trace) fprintf(trace,"iter,residual_rel\n");

    for(int iter=0; iter<max_iter; iter++){
        matvec(A, p.data(), Ap.data());

        double pAp=0;
        if(method==0){
            for(int i=0;i<n;i++) pAp+=p[i]*Ap[i];
        } else if(method==1){
            quire<32,2,2> q=0;
            for(int i=0;i<n;i++){
                posit<32,2> pa(p[i]),pb(Ap[i]);
                q+=quire_mul(pa,pb);
            }
            posit<32,2> pr; convert(q.to_value(),pr); pAp=double(pr);
        } else if(method==2){
            quire<16,1,2> q=0;
            for(int i=0;i<n;i++){
                posit<16,1> pa(p[i]),pb(Ap[i]);
                q+=quire_mul(pa,pb);
            }
            posit<16,1> pr; convert(q.to_value(),pr); pAp=double(pr);
        } else {
            posit<16,1> s=0;
            for(int i=0;i<n;i++){
                posit<16,1> pa(p[i]),pb(Ap[i]);
                s=s+pa*pb;
            }
            pAp=double(s);
        }

        if(!std::isfinite(pAp)||pAp<=0){
            double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
            if(trace) fprintf(trace,"%d,%.6e\n",iter,sqrt(rr)/b_norm);
            return {iter, sqrt(rr), false, b_norm};
        }

        double alpha=rz/pAp;
        for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
        for(int i=0;i<n;i++) z[i]=diagA[i]>0?r[i]/diagA[i]:r[i];
        double rz2=0; for(int i=0;i<n;i++) rz2+=r[i]*z[i];

        double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
        double res=sqrt(rr);
        double res_rel=res/b_norm;

        if(trace) fprintf(trace,"%d,%.6e\n",iter,res_rel);

        if(!std::isfinite(res)||res>1e30) return {iter,res,false,b_norm};
        if(res_rel<tol_rel) return {iter+1,res,true,b_norm};

        double beta=rz2/rz; rz=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
    }
    double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
    return {max_iter,sqrt(rr),false,b_norm};
}

int main(){
    MTX A = read_mtx("/mnt/d/posits-work/add32/add32.mtx", false);
    int n = A.rows;
    printf("Matrix: add32, n=%d, nnz=%d\n", n, (int)A.row.size());

    std::vector<double> b = read_rhs("/mnt/d/posits-work/add32/add32_b.mtx", n);

    std::vector<double> diagA(n,0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    // print matrix stats
    double dmin=1e300,dmax=0,vmin=1e300,vmax=0;
    for(int i=0;i<n;i++) if(diagA[i]>0){dmin=std::min(dmin,diagA[i]);dmax=std::max(dmax,diagA[i]);}
    for(auto v:A.val) if(v!=0){vmin=std::min(vmin,fabs(v));vmax=std::max(vmax,fabs(v));}
    printf("Diag range: %.3e to %.3e (ratio %.3e)\n",dmin,dmax,dmax/dmin);
    printf("Val range:  %.3e to %.3e (ratio %.3e)\n",vmin,vmax,vmax/vmin);

    const int MAX_ITER=2000;
    const double TOL=1e-8;

    struct Method { int id; const char* label; const char* tracefile; };
    std::vector<Method> methods = {
        {0, "double",       "/mnt/d/posits-work/add32_trace_double.csv"},
        {1, "posit32+quire","/mnt/d/posits-work/add32_trace_p32q.csv"},
        {2, "posit16+quire","/mnt/d/posits-work/add32_trace_p16q.csv"},
        {3, "posit16+naive","/mnt/d/posits-work/add32_trace_p16n.csv"},
    };

    printf("\n%-16s %-10s %-14s %-10s\n","method","iters","final_res_rel","converged");
    FILE* summary = fopen("/mnt/d/posits-work/add32_posit16_summary.log","w");
    fprintf(summary,"method,iters,final_res_rel,converged\n");

    for(auto& m : methods){
        FILE* trace = fopen(m.tracefile,"w");
        auto res = run_cg(A, b, diagA, m.id, MAX_ITER, TOL, trace, m.label);
        fclose(trace);

        printf("%-16s %-10d %-14.3e %-10s\n",
               m.label, res.iters, res.final_residual/res.b_norm,
               res.converged?"YES":"NO");
        fprintf(summary,"%s,%d,%.6e,%s\n",
                m.label, res.iters, res.final_residual/res.b_norm,
                res.converged?"YES":"NO");
        fflush(stdout);
    }
    fclose(summary);
    printf("\nTrace CSVs written for plotting.\n");
    return 0;
}
