// cg_compare_ilu_double.cpp
// Double64-only CG with ILU(0) preconditioning, to test whether a stronger
// preconditioner fixes the oscillation/non-convergence seen under Jacobi
// on bcsstk36, bcsstk37, nasasrb. If this works, extend to posit32 variants.
//
// Usage: ./cg_compare_ilu_double <matrix.mtx> <logfile> [maxiter]
//
// Build: g++ -O2 -std=c++17 -o cg_compare_ilu_double cg_compare_ilu_double.cpp

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>

struct CSR {
    int n;
    std::vector<int> row_ptr;   // size n+1
    std::vector<int> col_idx;   // size nnz
    std::vector<double> val;    // size nnz
};

// Read MatrixMarket triplet format (symmetric, mirrored), build CSR directly.
CSR read_mtx_csr(const char* fname){
    FILE* fp = fopen(fname, "r");
    if(!fp){ fprintf(stderr,"ERROR: cannot open %s\n", fname); exit(1); }
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int rows, cols, nnz_declared;
    sscanf(buf, "%d %d %d", &rows, &cols, &nnz_declared);

    // Use a map<(row,col), val> to naturally handle symmetric mirroring
    // and duplicate-entry accumulation (MatrixMarket allows repeated indices).
    std::vector<std::map<int,double>> adj(rows);

    int r,c; double v;
    while(fscanf(fp,"%d %d %lf",&r,&c,&v)==3){
        r--; c--; // 1-indexed -> 0-indexed
        adj[r][c] += v;
        if(r!=c) adj[c][r] += v;
    }
    fclose(fp);

    CSR A;
    A.n = rows;
    A.row_ptr.resize(rows+1, 0);
    for(int i=0;i<rows;i++) A.row_ptr[i+1] = A.row_ptr[i] + (int)adj[i].size();
    A.col_idx.resize(A.row_ptr[rows]);
    A.val.resize(A.row_ptr[rows]);
    int k=0;
    for(int i=0;i<rows;i++){
        for(auto& kv : adj[i]){
            A.col_idx[k] = kv.first;
            A.val[k] = kv.second;
            k++;
        }
    }
    return A;
}

void matvec(const CSR& A, const std::vector<double>& x, std::vector<double>& y){
    for(int i=0;i<A.n;i++){
        double s=0;
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++)
            s += A.val[k]*x[A.col_idx[k]];
        y[i]=s;
    }
}

double dot(const std::vector<double>& a, const std::vector<double>& b){
    double s=0; for(size_t i=0;i<a.size();i++) s+=a[i]*b[i]; return s;
}

// ---- ILU(0): same sparsity pattern as A, standard incomplete factorization ----
// Produces L (unit lower triangular, stored without the implicit 1s) and
// U (upper triangular), both in CSR sharing A's sparsity pattern.
// Reference: Saad, "Iterative Methods for Sparse Linear Systems", Alg 10.4.
struct ILU0 {
    CSR L; // strictly lower, unit diagonal implicit
    CSR U; // upper, including diagonal
};

ILU0 compute_ilu0(const CSR& A){
    int n = A.n;
    // Work on a dense-per-row sparse copy of A's values, factorize in place.
    // Use a full copy of A's CSR values as the working factorization array.
    std::vector<double> LU_val = A.val; // in-place factorization, same pattern as A

    // For fast column lookup within a row's sparsity pattern:
    auto find_idx = [&](int row, int col)->int{
        for(int k=A.row_ptr[row]; k<A.row_ptr[row+1]; k++)
            if(A.col_idx[k]==col) return k;
        return -1;
    };

    for(int i=1;i<n;i++){
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
            int col_k = A.col_idx[k];
            if(col_k >= i) continue; // only process entries left of diagonal
            int piv_idx = find_idx(col_k, col_k); // diagonal of row col_k
            if(piv_idx<0 || LU_val[piv_idx]==0.0) continue; // skip if pattern missing/zero pivot
            double factor = LU_val[k] / LU_val[piv_idx];
            LU_val[k] = factor;
            for(int j=A.row_ptr[i]; j<A.row_ptr[i+1]; j++){
                int col_j = A.col_idx[j];
                if(col_j <= col_k) continue;
                int uk = find_idx(col_k, col_j);
                if(uk>=0) LU_val[j] -= factor * LU_val[uk];
                // if uk<0, entry (col_k,col_j) not in pattern -> drop (this IS the "0" in ILU(0))
            }
        }
    }

    ILU0 result;
    result.L.n = n; result.U.n = n;
    result.L.row_ptr.assign(n+1,0);
    result.U.row_ptr.assign(n+1,0);
    for(int i=0;i<n;i++){
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
            if(A.col_idx[k] < i) result.L.row_ptr[i+1]++;
            else result.U.row_ptr[i+1]++;
        }
    }
    for(int i=0;i<n;i++){ result.L.row_ptr[i+1]+=result.L.row_ptr[i]; result.U.row_ptr[i+1]+=result.U.row_ptr[i]; }
    result.L.col_idx.resize(result.L.row_ptr[n]);
    result.L.val.resize(result.L.row_ptr[n]);
    result.U.col_idx.resize(result.U.row_ptr[n]);
    result.U.val.resize(result.U.row_ptr[n]);
    std::vector<int> lpos(n), upos(n);
    for(int i=0;i<n;i++){ lpos[i]=result.L.row_ptr[i]; upos[i]=result.U.row_ptr[i]; }
    for(int i=0;i<n;i++){
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
            int col = A.col_idx[k];
            if(col < i){
                result.L.col_idx[lpos[i]] = col;
                result.L.val[lpos[i]] = LU_val[k];
                lpos[i]++;
            } else {
                result.U.col_idx[upos[i]] = col;
                result.U.val[upos[i]] = LU_val[k];
                upos[i]++;
            }
        }
    }
    return result;
}

// Solve (LU) z = r  via forward substitution (Lw=r) then back substitution (Uz=w)
// L has implicit unit diagonal.
void ilu_solve(const ILU0& ilu, const std::vector<double>& r, std::vector<double>& z){
    int n = ilu.L.n;
    std::vector<double> w(n);
    // forward: L w = r
    for(int i=0;i<n;i++){
        double s = r[i];
        for(int k=ilu.L.row_ptr[i]; k<ilu.L.row_ptr[i+1]; k++)
            s -= ilu.L.val[k]*w[ilu.L.col_idx[k]];
        w[i]=s; // unit diagonal, no division
    }
    // backward: U z = w
    for(int i=n-1;i>=0;i--){
        double s = w[i];
        double diag=1.0;
        for(int k=ilu.U.row_ptr[i]; k<ilu.U.row_ptr[i+1]; k++){
            if(ilu.U.col_idx[k]==i) diag = ilu.U.val[k];
            else if(ilu.U.col_idx[k] > i) s -= ilu.U.val[k]*z[ilu.U.col_idx[k]];
        }
        z[i] = s/diag;
    }
}

int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: cg_compare_ilu_double <matrix.mtx> <logfile> [maxiter]\n"); return 1; }
    int maxiter = (argc>=4) ? atoi(argv[3]) : 2000;

    CSR A = read_mtx_csr(argv[1]);
    int n = A.n;
    printf("Matrix: %s  n=%d  nnz=%zu\n", argv[1], n, A.val.size());

    printf("Computing ILU(0) factorization...\n");
    ILU0 ilu = compute_ilu0(A);
    printf("ILU(0) done.\n");

    // --- Jacobi-preconditioned double64 CG (baseline, for direct comparison) ---
    std::vector<double> diagA(n);
    for(int i=0;i<n;i++)
        for(int k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++)
            if(A.col_idx[k]==i) diagA[i]=A.val[k];

    auto run_cg = [&](bool use_ilu, FILE* log, const char* label){
        std::vector<double> x(n,0), r(n,1), p(n), Ap(n), z(n);
        if(use_ilu) ilu_solve(ilu, r, z);
        else for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
        p = z;
        double rz = dot(r,z);
        for(int iter=0; iter<maxiter; iter++){
            matvec(A,p,Ap);
            double pAp = dot(p,Ap);
            double alpha = rz/pAp;
            for(int i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }
            if(use_ilu) ilu_solve(ilu, r, z);
            else for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
            double rz2 = dot(r,z);
            double beta = rz2/rz; rz=rz2;
            for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
            double res = sqrt(dot(r,r));
            fprintf(log, "%s %4d %.10e\n", label, iter, res);
            if(res < 1e-10){ printf("%s converged at iter %d\n", label, iter); break; }
            if(iter==maxiter-1) printf("%s did NOT converge within %d iters (final res=%.3e)\n", label, maxiter, res);
        }
    };

    FILE* log = fopen(argv[2], "w");
    fprintf(log, "variant iter residual\n");
    printf("Running Jacobi baseline...\n");
    run_cg(false, log, "jacobi");
    printf("Running ILU(0)...\n");
    run_cg(true, log, "ilu0");
    fclose(log);

    printf("Done. Log: %s\n", argv[2]);
    return 0;
}
// (verification helper — not part of the main flow, just documenting the check)
