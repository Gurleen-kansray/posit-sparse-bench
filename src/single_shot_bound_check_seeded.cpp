#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;
using P32 = posit<32,2>;

struct MTX { int n; std::vector<int> row,col; std::vector<double> val; };
MTX read_mtx(const char* f){
    MTX m; FILE* fp=fopen(f,"r");
    if(!fp){ printf("ERROR: cannot open %s\n",f); exit(1); }
    char buf[256];
    while(fgets(buf,256,fp) && buf[0]=='%');
    int rows,cols,nnz; sscanf(buf,"%d %d %d",&rows,&cols,&nnz);
    m.n=rows; int r,c; double v;
    while(true){
        int nread = fscanf(fp,"%d %d %lf",&r,&c,&v);
        if(nread==2) v=1.0; // pattern format: no value column, implicit 1.0
        else if(nread!=3) break;
        m.row.push_back(r-1); m.col.push_back(c-1); m.val.push_back(v);
        if(r!=c){m.row.push_back(c-1); m.col.push_back(r-1); m.val.push_back(v);}
    } fclose(fp); return m;
}
void matvec_d(const MTX& A, const std::vector<double>& x, std::vector<double>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=A.val[k]*x[A.col[k]];
}
double dot_d(const std::vector<double>& a, const std::vector<double>& b, int n){
    double s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
double dot_p_quire(const std::vector<P32>& a, const std::vector<P32>& b, int n){
    quire<32,2,2> q=0;
    for(int i=0;i<n;i++) q+=quire_mul(a[i],b[i]);
    P32 r; convert(q.to_value(),r); return double(r);
}
int main(int argc, char** argv){
    if(argc<3){ printf("usage: %s matrix.mtx seed\n", argv[0]); return 1; }
    MTX A = read_mtx(argv[1]);
    int n = A.n;
    unsigned seed = atoi(argv[2]);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.5,1.5);
    std::vector<double> p(n);
    for(int i=0;i<n;i++) p[i]=dist(rng);
    std::vector<double> Ap(n);
    matvec_d(A, p, Ap);
    double s_exact = dot_d(p, Ap, n);
    std::vector<P32> p32(n), Ap32(n);
    for(int i=0;i<n;i++){ p32[i]=P32(p[i]); Ap32[i]=P32(Ap[i]); }
    double s_quire = dot_p_quire(p32, Ap32, n);
    double rel_err = fabs(s_quire - s_exact) / fabs(s_exact);
    double u = 3.725290298e-09;
    printf("matrix=%s seed=%u rel_err=%.6e ratio_to_u=%.3f %s\n",
        argv[1], seed, rel_err, rel_err/u, (rel_err<=u ? "PASS" : "FAIL"));
    return 0;
}
