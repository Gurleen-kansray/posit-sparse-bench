// verify_ilu.cpp
// Independent verification, separate from the CG loop:
// 1. Checks ILU(0) factorization sanity: does L*U match A on A's sparsity pattern?
// 2. Runs CG (Jacobi and ILU0) and independently recomputes ||b - Ax|| from
//    scratch at the end, rather than trusting the loop's self-reported residual.
//
// Usage: ./verify_ilu <matrix.mtx> [maxiter]
// Build: g++ -O2 -std=c++17 -o verify_ilu verify_ilu.cpp

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

struct CSR {
    int n;
    std::vector<int> row_ptr, col_idx;
    std::vector<double> val;
};

CSR read_mtx_csr(const char* fname){
    FILE* fp = fopen(fname, "r");
    if(!fp){ fprintf(stderr,"ERROR: cannot open %s\n", fname); exit(1); }
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int rows, cols, nnz_declared;
    sscanf(buf, "%d %d %d", &rows, &cols, &nnz_declared);
    std::vector<std::map<int,double>> adj(rows);
    int r,c; double v;
    while(fscanf(fp,"%d %d %lf",&r,&c,&v)==3){
        r--; c--;
        adj[r][c] += v;
        if(r!=c) adj[c][r] += v;
    }
    fclose(fp);
    CSR A; A.n = rows;
    A.row_ptr.assign(rows+1,0);
    for(int i=0;i<rows;i++) A.row_ptr[i+1] = A.row_ptr[i] + (int)adj[i].size();
    A.col_idx.resize(A.row_ptr[rows]);
    A.val.resize(A.row_ptr[rows]);
    int k=0;
    for(int i=0;i<rows;i++) for(auto& kv: adj[i]){ A.col_idx[k]=kv.first; A.val[k]=kv.second; k++; }
    return A;
}

void matvec(const CSR& A, const std::vector<double>& x, std::vector<double>& y){
    for(int i=0;i<A.n;i++){
        double s=0;
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++) s += A.val[k]*x[A.col_idx[k]];
        y[i]=s;
    }
}
double dot(const std::vector<double>& a, const std::vector<double>& b){
    double s=0; for(size_t i=0;i<a.size();i++) s+=a[i]*b[i]; return s;
}

struct ILU0 { CSR L, U; };

ILU0 compute_ilu0(const CSR& A){
    int n = A.n;
    std::vector<double> LU_val = A.val;
    auto find_idx = [&](int row, int col)->int{
        for(int k=A.row_ptr[row]; k<A.row_ptr[row+1]; k++) if(A.col_idx[k]==col) return k;
        return -1;
    };
    for(int i=1;i<n;i++){
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
            int col_k = A.col_idx[k];
            if(col_k >= i) continue;
            int piv_idx = find_idx(col_k, col_k);
            if(piv_idx<0 || LU_val[piv_idx]==0.0) continue;
            double factor = LU_val[k] / LU_val[piv_idx];
            LU_val[k] = factor;
            for(int j=A.row_ptr[i]; j<A.row_ptr[i+1]; j++){
                int col_j = A.col_idx[j];
                if(col_j <= col_k) continue;
                int uk = find_idx(col_k, col_j);
                if(uk>=0) LU_val[j] -= factor * LU_val[uk];
            }
        }
    }
    ILU0 result; result.L.n=n; result.U.n=n;
    result.L.row_ptr.assign(n+1,0); result.U.row_ptr.assign(n+1,0);
    for(int i=0;i<n;i++) for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
        if(A.col_idx[k]<i) result.L.row_ptr[i+1]++; else result.U.row_ptr[i+1]++;
    }
    for(int i=0;i<n;i++){ result.L.row_ptr[i+1]+=result.L.row_ptr[i]; result.U.row_ptr[i+1]+=result.U.row_ptr[i]; }
    result.L.col_idx.resize(result.L.row_ptr[n]); result.L.val.resize(result.L.row_ptr[n]);
    result.U.col_idx.resize(result.U.row_ptr[n]); result.U.val.resize(result.U.row_ptr[n]);
    std::vector<int> lpos(n), upos(n);
    for(int i=0;i<n;i++){ lpos[i]=result.L.row_ptr[i]; upos[i]=result.U.row_ptr[i]; }
    for(int i=0;i<n;i++) for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
        int col=A.col_idx[k];
        if(col<i){ result.L.col_idx[lpos[i]]=col; result.L.val[lpos[i]]=LU_val[k]; lpos[i]++; }
        else { result.U.col_idx[upos[i]]=col; result.U.val[upos[i]]=LU_val[k]; upos[i]++; }
    }
    return result;
}

void ilu_solve(const ILU0& ilu, const std::vector<double>& r, std::vector<double>& z){
    int n = ilu.L.n;
    std::vector<double> w(n);
    for(int i=0;i<n;i++){
        double s=r[i];
        for(int k=ilu.L.row_ptr[i]; k<ilu.L.row_ptr[i+1]; k++) s -= ilu.L.val[k]*w[ilu.L.col_idx[k]];
        w[i]=s;
    }
    for(int i=n-1;i>=0;i--){
        double s=w[i]; double diag=1.0;
        for(int k=ilu.U.row_ptr[i]; k<ilu.U.row_ptr[i+1]; k++){
            if(ilu.U.col_idx[k]==i) diag=ilu.U.val[k];
            else if(ilu.U.col_idx[k]>i) s -= ilu.U.val[k]*z[ilu.U.col_idx[k]];
        }
        z[i]=s/diag;
    }
}

int main(int argc, char* argv[]){
    if(argc<2){ printf("usage: verify_ilu <matrix.mtx> [maxiter]\n"); return 1; }
    int maxiter = (argc>=3) ? atoi(argv[2]) : 2000;
    CSR A = read_mtx_csr(argv[1]);
    int n = A.n;
    printf("=== Matrix: %s  n=%d  nnz=%zu ===\n", argv[1], n, A.val.size());

    ILU0 ilu = compute_ilu0(A);

    // CHECK 1: L*U reconstruction vs A, on A's own sparsity pattern only.
    // This is exactly what ILU(0) promises: exact match on the existing
    // pattern (up to floating point), error only in dropped fill-in positions.
    printf("\n--- CHECK 1: ILU(0) pattern-match reconstruction ---\n");
    double max_pattern_err = 0, sum_pattern_err = 0;
    int checked = 0;
    for(int i=0;i<n;i++){
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
            int j = A.col_idx[k];
            double a_ij = A.val[k];
            double lu_ij = 0;
            for(int lk=ilu.L.row_ptr[i]; lk<ilu.L.row_ptr[i+1]; lk++){
                int m = ilu.L.col_idx[lk];
                double l_im = ilu.L.val[lk];
                for(int uk=ilu.U.row_ptr[m]; uk<ilu.U.row_ptr[m+1]; uk++){
                    if(ilu.U.col_idx[uk]==j){ lu_ij += l_im * ilu.U.val[uk]; break; }
                }
            }
            for(int uk=ilu.U.row_ptr[i]; uk<ilu.U.row_ptr[i+1]; uk++){
                if(ilu.U.col_idx[uk]==j){ lu_ij += ilu.U.val[uk]; break; }
            }
            double err = fabs(lu_ij - a_ij);
            max_pattern_err = std::max(max_pattern_err, err);
            sum_pattern_err += err;
            checked++;
        }
    }
    printf("Checked %d pattern entries. Max |LU - A| = %.6e, Mean |LU - A| = %.6e\n",
           checked, max_pattern_err, sum_pattern_err/checked);
    printf("(Should be ~machine epsilon if ILU(0) is implemented correctly)\n");

    // CHECK 2: independently recompute true residual ||b - Ax|| after CG,
    // not trusting the loop's self-tracked residual.
    printf("\n--- CHECK 2: independent residual recomputation ---\n");
    std::vector<double> diagA(n);
    for(int i=0;i<n;i++) for(int k=A.row_ptr[i];k<A.row_ptr[i+1];k++) if(A.col_idx[k]==i) diagA[i]=A.val[k];

    auto run_and_verify = [&](bool use_ilu, const char* label){
        std::vector<double> x(n,0), r(n,1), p(n), Ap(n), z(n), b(n,1);
        if(use_ilu) ilu_solve(ilu,r,z); else for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
        p=z; double rz=dot(r,z);
        int final_iter=0;
        for(int iter=0; iter<maxiter; iter++){
            matvec(A,p,Ap);
            double pAp=dot(p,Ap);
            double alpha=rz/pAp;
            for(int i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }
            if(use_ilu) ilu_solve(ilu,r,z); else for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
            double rz2=dot(r,z);
            double beta=rz2/rz; rz=rz2;
            for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
            final_iter=iter;
            double loop_res = sqrt(dot(r,r));
            if(loop_res<1e-10) break;
        }
        std::vector<double> Ax(n);
        matvec(A,x,Ax);
        std::vector<double> true_r(n);
        for(int i=0;i<n;i++) true_r[i] = b[i]-Ax[i];
        double true_res = sqrt(dot(true_r,true_r));
        printf("%s: stopped at iter %d, independently recomputed ||b-Ax|| = %.6e\n",
               label, final_iter, true_res);
    };

    run_and_verify(false, "jacobi");
    run_and_verify(true, "ilu0");

    return 0;
}
