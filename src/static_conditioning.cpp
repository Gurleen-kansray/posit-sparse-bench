// static_conditioning.cpp
// Prof. Quinlan's Part A: static quantization/conditioning sweep, per matrix per bit-width T
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <random>
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

void jacobi_pcg(const MTX& A, const std::vector<double>& b, std::vector<double>& x,
                 int max_iters=2000, double tol=1e-10){
    int n = A.n;
    std::vector<double> diag(n,0.0);
    for(int k=0;k<(int)A.row.size();k++) if(A.row[k]==A.col[k]) diag[A.row[k]] += A.val[k];
    x.assign(n,0.0);
    std::vector<double> r=b, z(n), p, Ap(n);
    for(int i=0;i<n;i++) z[i] = r[i]/diag[i];
    p = z;
    double rz_old=0; for(int i=0;i<n;i++) rz_old += r[i]*z[i];
    double bnorm=0; for(double bi:b) bnorm += bi*bi; bnorm = std::sqrt(bnorm);
    if(bnorm==0) bnorm=1.0;

    for(int it=0; it<max_iters; it++){
        matvec_d(A, p, Ap);
        double pAp=0; for(int i=0;i<n;i++) pAp += p[i]*Ap[i];
        if(pAp==0) break;
        double alpha = rz_old/pAp;
        for(int i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }
        double rnorm=0; for(double ri:r) rnorm+=ri*ri; rnorm=std::sqrt(rnorm);
        if(rnorm/bnorm < tol) break;
        for(int i=0;i<n;i++) z[i]=r[i]/diag[i];
        double rz_new=0; for(int i=0;i<n;i++) rz_new += r[i]*z[i];
        double beta = rz_new/rz_old;
        for(int i=0;i<n;i++) p[i] = z[i] + beta*p[i];
        rz_old = rz_new;
    }
}

double cmsw_condest(const MTX& A){
    int n = A.n;
    std::vector<double> colsum(n,0.0);
    for(int k=0;k<(int)A.row.size();k++) colsum[A.col[k]] += std::fabs(A.val[k]);
    double norm1_A = 0.0;
    for(double cs : colsum) norm1_A = std::max(norm1_A, cs);

    std::vector<double> e(n);
    for(int i=0;i<n;i++) e[i] = (i%2==0) ? 1.0 : -1.0;

    std::vector<double> y;
    jacobi_pcg(A, e, y);
    std::vector<double> z;
    jacobi_pcg(A, y, z);

    double norm1_y=0, norm1_z=0;
    for(double yi:y) norm1_y += std::fabs(yi);
    for(double zi:z) norm1_z += std::fabs(zi);

    double est_norm_Ainv = (norm1_y>0) ? (norm1_z/norm1_y) : 0.0;
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

int main(int argc, char* argv[]){
    if(argc<2){ printf("usage: static_conditioning <matrix.mtx>\n"); return 1; }
    std::string path = argv[1];
    std::string name = path.substr(path.find_last_of("/\\")+1);
    name = name.substr(0, name.find_last_of('.'));

    MTX A = read_mtx(argv[1]);
    printf("[static_conditioning] %s: n=%d nnz=%zu\n", name.c_str(), A.n, A.val.size());

    double lambda_max = power_iteration_lambda_max(A);
    double condest = cmsw_condest(A);
    printf("[static_conditioning] lambda_max=%.6e condest=%.6e\n", lambda_max, condest);

    bool write_header = false;
    { std::ifstream test("results/csv/static_conditioning.csv"); write_header = !test.good(); }
    std::ofstream csv("results/csv/static_conditioning.csv", std::ios::app);
    if(write_header) csv << "matrix,bitwidth,lambda_max,condest_cmsw,nnz_static,saturated_count,total_entries,sat_fraction\n";

    long sat, nnz_q, total;

    saturation_and_nnz<8,0>(A.val, sat, nnz_q, total);
    csv << name << ",8," << lambda_max << "," << condest << "," << nnz_q << "," << sat << "," << total << "," << (total?double(sat)/total:0) << "\n";

    saturation_and_nnz<16,1>(A.val, sat, nnz_q, total);
    csv << name << ",16," << lambda_max << "," << condest << "," << nnz_q << "," << sat << "," << total << "," << (total?double(sat)/total:0) << "\n";

    saturation_and_nnz<32,2>(A.val, sat, nnz_q, total);
    csv << name << ",32," << lambda_max << "," << condest << "," << nnz_q << "," << sat << "," << total << "," << (total?double(sat)/total:0) << "\n";

    saturation_and_nnz<64,2>(A.val, sat, nnz_q, total);
    csv << name << ",64," << lambda_max << "," << condest << "," << nnz_q << "," << sat << "," << total << "," << (total?double(sat)/total:0) << "\n";

    csv.close();
    printf("[static_conditioning] appended results for %s to results/csv/static_conditioning.csv\n", name.c_str());
    return 0;
}
