// static_conditioning.cpp  (v4 - RCM reordering + true CMSW via sparse Cholesky)
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <map>
#include <queue>
#include <algorithm>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

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

template <size_t nbits, size_t es>
MTX quantize_matrix(const MTX& A){
    MTX Aq = A;
    for(size_t k=0;k<Aq.val.size();k++){
        posit<nbits,es> p(Aq.val[k]);
        Aq.val[k] = double(p);
    }
    return Aq;
}

void matvec_d(const MTX& A, const std::vector<double>& x, std::vector<double>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=A.val[k]*x[A.col[k]];
}

double power_iteration_lambda_max(const MTX& A, int max_iters=500, double tol=1e-10){
    int n = A.n;
    std::vector<double> v(n), Av(n);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0,1.0);
    for(int i=0;i<n;i++) v[i]=dist(rng);
    auto normalize=[&](std::vector<double>& x){
        double norm=0; for(double xi:x) norm+=xi*xi; norm=std::sqrt(norm);
        if(norm>0) for(double& xi:x) xi/=norm;
    };
    normalize(v);
    double lambda=0, lambda_old=0;
    for(int it=0; it<max_iters; it++){
        matvec_d(A, v, Av);
        double rayleigh=0, vv=0;
        for(int i=0;i<n;i++){ rayleigh += v[i]*Av[i]; vv += v[i]*v[i]; }
        lambda = rayleigh/vv;
        normalize(Av);
        v = Av;
        if(std::fabs(lambda-lambda_old) < tol*std::fabs(lambda)) break;
        lambda_old = lambda;
    }
    return lambda;
}

std::vector<std::vector<int>> build_adjacency(const MTX& A){
    std::vector<std::vector<int>> adj(A.n);
    for(size_t k=0;k<A.row.size();k++){
        int i=A.row[k], j=A.col[k];
        if(i!=j) adj[i].push_back(j);
    }
    for(auto& nbrs : adj){
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
    return adj;
}

std::vector<int> rcm_order(const MTX& A){
    int n = A.n;
    auto adj = build_adjacency(A);
    std::vector<int> degree(n);
    for(int i=0;i<n;i++) degree[i] = (int)adj[i].size();

    std::vector<bool> visited(n, false);
    std::vector<int> order; order.reserve(n);

    for(int start=0; start<n; start++){
        if(visited[start]) continue;
        int root = start;
        for(int i=start;i<n;i++){
            if(!visited[i] && degree[i] < degree[root]) root = i;
        }
        std::queue<int> q;
        visited[root] = true;
        q.push(root);
        while(!q.empty()){
            int u = q.front(); q.pop();
            order.push_back(u);
            std::vector<int> nbrs;
            for(int v : adj[u]) if(!visited[v]) nbrs.push_back(v);
            std::sort(nbrs.begin(), nbrs.end(), [&](int a,int b){ return degree[a]<degree[b]; });
            for(int v : nbrs){
                if(!visited[v]){ visited[v]=true; q.push(v); }
            }
        }
    }
    std::reverse(order.begin(), order.end());
    return order;
}

MTX permute_matrix(const MTX& A, const std::vector<int>& perm){
    int n = A.n;
    std::vector<int> inv(n);
    for(int i=0;i<n;i++) inv[perm[i]] = i;
    MTX Ap; Ap.n = n;
    Ap.row.reserve(A.row.size());
    Ap.col.reserve(A.col.size());
    Ap.val = A.val;
    for(size_t k=0;k<A.row.size();k++){
        Ap.row.push_back(inv[A.row[k]]);
        Ap.col.push_back(inv[A.col[k]]);
    }
    return Ap;
}

struct SparseChol {
    int n;
    std::vector<std::map<int,double>> L;
    bool ok = true;
};

SparseChol sparse_cholesky(const MTX& A){
    int n = A.n;
    SparseChol C; C.n = n; C.L.assign(n, {});

    std::vector<std::map<int,double>> Amap(n);
    for(size_t k=0;k<A.row.size();k++){
        int i=A.row[k], j=A.col[k];
        if(j<=i) Amap[i][j] = A.val[k];
    }

    for(int j=0;j<n;j++){
        double sum = 0.0;
        for(auto& kv : C.L[j]) sum += kv.second * kv.second;
        double diag = 0.0;
        auto it = Amap[j].find(j);
        if(it != Amap[j].end()) diag = it->second;
        double d = diag - sum;
        if(d <= 0.0 || !std::isfinite(d)){
            C.ok = false;
            return C;
        }
        double Ljj = std::sqrt(d);
        C.L[j][j] = Ljj;

        for(int i=j+1;i<n;i++){
            auto ait = Amap[i].find(j);
            if(ait == Amap[i].end() && C.L[j].size()<=1) continue;
            double aij = (ait!=Amap[i].end()) ? ait->second : 0.0;

            double dot = 0.0;
            const auto& Li = C.L[i];
            const auto& Lj = C.L[j];
            if(!Li.empty()){
                for(auto& kv : Lj){
                    if(kv.first>=j) continue;
                    auto lit = Li.find(kv.first);
                    if(lit != Li.end()) dot += lit->second * kv.second;
                }
            }
            double val = (aij - dot) / Ljj;
            if(val != 0.0) C.L[i][j] = val;
        }
    }
    return C;
}

void chol_forward_greedy(const SparseChol& C, std::vector<double>& z, std::vector<double>& e_out){
    int n = C.n;
    z.assign(n, 0.0);
    e_out.assign(n, 0.0);
    for(int i=0;i<n;i++){
        double sum = 0.0;
        for(auto& kv : C.L[i]){
            if(kv.first < i) sum += kv.second * z[kv.first];
        }
        double Lii = C.L[i].at(i);
        double zp = ( 1.0 - sum) / Lii;
        double zm = (-1.0 - sum) / Lii;
        if(std::fabs(zp) >= std::fabs(zm)){
            z[i] = zp; e_out[i] = 1.0;
        } else {
            z[i] = zm; e_out[i] = -1.0;
        }
    }
}

void chol_back_solve(const SparseChol& C, const std::vector<double>& z, std::vector<double>& y){
    int n = C.n;
    std::vector<double> x = z;
    y.assign(n, 0.0);
    for(int i=n-1;i>=0;i--){
        double Lii = C.L[i].at(i);
        double val = x[i] / Lii;
        y[i] = val;
        for(auto& kv : C.L[i]){
            if(kv.first < i) x[kv.first] -= kv.second * val;
        }
    }
}

double cmsw_condest_true(const MTX& A_reordered, bool& success){
    success = true;
    SparseChol C = sparse_cholesky(A_reordered);
    if(!C.ok){ success = false; return -1.0; }

    int n = A_reordered.n;
    std::vector<double> colsum(n,0.0);
    for(size_t k=0;k<A_reordered.row.size();k++) colsum[A_reordered.col[k]] += std::fabs(A_reordered.val[k]);
    double norm1_A = 0.0;
    for(double cs : colsum) norm1_A = std::max(norm1_A, cs);

    std::vector<double> z, e_chosen;
    chol_forward_greedy(C, z, e_chosen);
    std::vector<double> y;
    chol_back_solve(C, z, y);

    double norm1_y=0, norm1_z=0;
    for(double zi : z) norm1_z += std::fabs(zi);
    for(double yi : y) norm1_y += std::fabs(yi);

    if(norm1_z <= 0.0 || !std::isfinite(norm1_y)){ success = false; return -1.0; }

    double est_norm_Ainv = norm1_y / norm1_z;
    return est_norm_Ainv * norm1_A;
}

template <size_t nbits, size_t es>
void saturation_and_nnz(const std::vector<double>& vals, long& sat, long& nnz_q, long& total){
    posit<nbits,es> maxp, minp;
    maxpos(maxp); minpos(minp);
    double maxv = double(maxp), minv = double(minp);
    sat = 0; nnz_q = 0; total = (long)vals.size();
    for(double v : vals){
        if(v==0.0) continue;
        posit<nbits,es> p(v);
        double pv = double(p);
        if(pv != 0.0) nnz_q++;
        if(std::fabs(pv)==maxv && std::fabs(v)>maxv) sat++;
        else if(std::fabs(pv)==minv && std::fabs(v)<minv) sat++;
    }
}

template <size_t nbits, size_t es>
void run_bitwidth(const MTX& A, const MTX& A_rcm, const std::string& name, int bits, std::ofstream& csv){
    MTX Aq = quantize_matrix<nbits,es>(A);
    MTX Aq_rcm = quantize_matrix<nbits,es>(A_rcm);
    double lambda_max = power_iteration_lambda_max(Aq);

    bool chol_ok;
    double condest = cmsw_condest_true(Aq_rcm, chol_ok);
    std::string condest_str;
    if(chol_ok){
        char buf[64]; snprintf(buf,64,"%.6e", condest); condest_str = buf;
    } else {
        condest_str = "CHOL_FAIL";
    }

    long sat, nnz_q, total;
    saturation_and_nnz<nbits,es>(A.val, sat, nnz_q, total);

    printf("[static_conditioning] %s bits=%d lambda_max=%.6e condest=%s sat_frac=%.4f\n",
           name.c_str(), bits, lambda_max, condest_str.c_str(), total?double(sat)/total:0.0);

    csv << name << "," << bits << "," << lambda_max << "," << condest_str << ","
        << nnz_q << "," << sat << "," << total << "," << (total?double(sat)/total:0) << "\n";
}

int main(int argc, char* argv[]){
    if(argc<2){ printf("usage: static_conditioning <matrix.mtx>\n"); return 1; }
    std::string path = argv[1];
    std::string name = path.substr(path.find_last_of("/\\")+1);
    name = name.substr(0, name.find_last_of('.'));

    MTX A = read_mtx(argv[1]);
    printf("[static_conditioning] %s: n=%d nnz=%zu\n", name.c_str(), A.n, A.val.size());

    printf("[static_conditioning] computing RCM ordering...\n");
    std::vector<int> perm = rcm_order(A);
    MTX A_rcm = permute_matrix(A, perm);
    printf("[static_conditioning] RCM ordering done.\n");

    bool write_header = false;
    { std::ifstream test("results/csv/static_conditioning.csv"); write_header = !test.good(); }
    std::ofstream csv("results/csv/static_conditioning.csv", std::ios::app);
    if(write_header) csv << "matrix,bitwidth,lambda_max,condest_cmsw,nnz_static,saturated_count,total_entries,sat_fraction\n";

    run_bitwidth<8,0>(A, A_rcm, name, 8, csv);
    run_bitwidth<16,1>(A, A_rcm, name, 16, csv);
    run_bitwidth<32,2>(A, A_rcm, name, 32, csv);
    run_bitwidth<64,2>(A, A_rcm, name, 64, csv);

    csv.close();
    printf("[static_conditioning] appended results for %s to results/csv/static_conditioning.csv\n", name.c_str());
    return 0;
}
