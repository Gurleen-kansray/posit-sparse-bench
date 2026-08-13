#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

using P8  = posit<8,2>;
using P16 = posit<16,2>;
using P32 = posit<32,2>;
using P64 = posit<64,2>;

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
void matvec_f(const MTX& A, const std::vector<float>& x, std::vector<float>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=(float)A.val[k]*x[A.col[k]];
}
template<size_t N>
void matvec_p(const MTX& A, const std::vector<posit<N,2>>& x, std::vector<posit<N,2>>& y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]=y[A.row[k]]+posit<N,2>(A.val[k])*x[A.col[k]];
}

double dot_d(const std::vector<double>& a, const std::vector<double>& b, int n){
    double s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
float dot_f(const std::vector<float>& a, const std::vector<float>& b, int n){
    float s=0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}
float dot_f_fma(const std::vector<float>& a, const std::vector<float>& b, int n){
    float s=0; for(int i=0;i<n;i++) s=std::fma(a[i],b[i],s); return s;
}
// N=8/16/32/64, ES=2 uniformly (per ratified standard)
template<size_t N>
double dot_p_quire(const std::vector<posit<N,2>>& a, const std::vector<posit<N,2>>& b, int n){
    quire<N,2,2> q=0;
    for(int i=0;i<n;i++) q+=quire_mul(a[i],b[i]);
    posit<N,2> r; convert(q.to_value(),r); return double(r);
}
template<size_t N>
double dot_p_naive(const std::vector<posit<N,2>>& a, const std::vector<posit<N,2>>& b, int n){
    posit<N,2> s=0;
    for(int i=0;i<n;i++) s=s+a[i]*b[i];
    return double(s);
}

int main(int argc, char* argv[]){
    if(argc<4){ printf("usage: cg_compare_seeded <matrix.mtx> <logfile> <seed>\n"); return 1; }
    const char* mtx_path = argv[1];
    const char* log_path = argv[2];
    unsigned int seed = (unsigned int)atoi(argv[3]);

    MTX A = read_mtx(mtx_path);
    int n = A.n;

    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    // (1) Known ground truth x_true, then b = A*x_true -> real solution error is computable
    std::mt19937 rng(seed);
    std::normal_distribution<double> ndist(0.0, 1.0);
    std::vector<double> x_true(n);
    for(int i=0;i<n;i++) x_true[i] = ndist(rng);
    std::vector<double> b(n);
    matvec_d(A, x_true, b);
    double x_true_norm = sqrt(dot_d(x_true, x_true, n));

    FILE* log=fopen(log_path,"w");
    fprintf(log,"matrix=%s n=%d seed=%u\n", mtx_path, n, seed);
    fprintf(log,"iter res_d res_f res_ffma res_p8q res_p8n res_p16q res_p16n res_p32q res_p32n res_p64q res_p64n "
                "solerr_d solerr_f solerr_ffma solerr_p32q solerr_p32n "
                "pAp_d pAp_p8q pAp_p8n pAp_p16q pAp_p16n pAp_p32q pAp_p32n pAp_p64q pAp_p64n abs_dot_d cond_xy_d "
                "rz_d rz_p32q rz_p32n "
                "guard_p8q guard_p8n guard_p16q guard_p16n guard_p32q guard_p32n guard_p64q guard_p64n\n");

    // double
    std::vector<double> xd(n,0),rd(b),pd(n),Apd(n),zd(n);
    for(int i=0;i<n;i++) zd[i]=rd[i]/diagA[i];
    for(int i=0;i<n;i++) pd[i]=zd[i];
    double rzd=dot_d(rd,zd,n);
    double bnorm = sqrt(dot_d(b,b,n));
    fprintf(log,"bnorm=%.10e\n", bnorm);

    // float32 naive
    std::vector<float> xf(n,0),rf(n),pf(n),Apf(n),zf(n);
    for(int i=0;i<n;i++) rf[i]=(float)b[i];
    for(int i=0;i<n;i++) zf[i]=rf[i]/(float)diagA[i];
    for(int i=0;i<n;i++) pf[i]=zf[i];
    float rzf=dot_f(rf,zf,n);

    // float32 fma
    std::vector<float> xg(n,0),rg(rf),pg(n),Apg(n),zg(n);
    for(int i=0;i<n;i++) zg[i]=zf[i];
    for(int i=0;i<n;i++) pg[i]=zg[i];
    float rzg=dot_f_fma(rg,zg,n);

#define SETUP_POSIT(W) \
    std::vector<P##W> xp##W##q(n,0), rp##W##q(n), pp##W##q(n), App##W##q(n), zp##W##q(n); \
    for(int i=0;i<n;i++) rp##W##q[i]=P##W(b[i]); \
    for(int i=0;i<n;i++) zp##W##q[i]=P##W(1.0/diagA[i])*rp##W##q[i]; \
    for(int i=0;i<n;i++) pp##W##q[i]=zp##W##q[i]; \
    double rzp##W##q = dot_p_quire<W>(rp##W##q,zp##W##q,n); \
    std::vector<P##W> xp##W##n(n,0), rp##W##n(rp##W##q), pp##W##n(n), App##W##n(n), zp##W##n(n); \
    for(int i=0;i<n;i++) zp##W##n[i]=zp##W##q[i]; \
    for(int i=0;i<n;i++) pp##W##n[i]=zp##W##n[i]; \
    double rzp##W##n = dot_p_naive<W>(rp##W##n,zp##W##n,n);

    SETUP_POSIT(8)
    SETUP_POSIT(16)
    SETUP_POSIT(32)
    SETUP_POSIT(64)

    int maxiter = 2000;
    int conv_iter = -1;
    int conv_iter_p32q = -1;
    int conv_iter_p32n = -1;
    int conv_iter_p32q_loose = -1;
    int conv_iter_p32n_loose = -1;

    for(int iter=0; iter<maxiter; iter++){
        double rzd_used = rzd;
        double rzp32q_used = rzp32q;
        double rzp32n_used = rzp32n;
        // double
        matvec_d(A,pd,Apd);
        double pApd = dot_d(pd,Apd,n);
        double abs_dot_d = 0; for(int i=0;i<n;i++) abs_dot_d += fabs(pd[i])*fabs(Apd[i]);
        double cond_xy_d = fabs(pApd) > 0 ? abs_dot_d / fabs(pApd) : INFINITY;
        double alphad = rzd/pApd;
        for(int i=0;i<n;i++){ xd[i]+=alphad*pd[i]; rd[i]-=alphad*Apd[i]; }
        for(int i=0;i<n;i++) zd[i]=rd[i]/diagA[i];
        double rzd2=dot_d(rd,zd,n);
        double betad=rzd2/rzd; rzd=rzd2;
        for(int i=0;i<n;i++) pd[i]=zd[i]+betad*pd[i];
        double resd=sqrt(dot_d(rd,rd,n));
        std::vector<double> ed(n); for(int i=0;i<n;i++) ed[i]=xd[i]-x_true[i];
        double solerr_d = sqrt(dot_d(ed,ed,n)) / x_true_norm;

        // float32 naive
        matvec_f(A,pf,Apf);
        float pApf = dot_f(pf,Apf,n);
        float alphaf = rzf/pApf;
        for(int i=0;i<n;i++){ xf[i]+=alphaf*pf[i]; rf[i]-=alphaf*Apf[i]; }
        for(int i=0;i<n;i++) zf[i]=(float)(1.0/diagA[i])*rf[i];
        float rzf2=dot_f(rf,zf,n);
        float betaf=rzf2/rzf; rzf=rzf2;
        for(int i=0;i<n;i++) pf[i]=zf[i]+betaf*pf[i];
        double resf=sqrt((double)dot_f(rf,rf,n));
        double se=0; for(int i=0;i<n;i++){ double e=(double)xf[i]-x_true[i]; se+=e*e; }
        double solerr_f = sqrt(se)/x_true_norm;

        // float32 fma
        matvec_f(A,pg,Apg);
        float pApg = dot_f_fma(pg,Apg,n);
        float alphag = rzg/pApg;
        for(int i=0;i<n;i++){ xg[i]+=alphag*pg[i]; rg[i]-=alphag*Apg[i]; }
        for(int i=0;i<n;i++) zg[i]=(float)(1.0/diagA[i])*rg[i];
        float rzg2=dot_f_fma(rg,zg,n);
        float betag=rzg2/rzg; rzg=rzg2;
        for(int i=0;i<n;i++) pg[i]=zg[i]+betag*pg[i];
        double resffma=sqrt((double)dot_f_fma(rg,rg,n));
        double se_ffma=0; for(int i=0;i<n;i++){ double e=(double)xg[i]-x_true[i]; se_ffma+=e*e; }
        double solerr_ffma = sqrt(se_ffma)/x_true_norm;

#define STEP_POSIT(W) \
        matvec_p<W>(A,pp##W##q,App##W##q); \
        double pApp##W##q = dot_p_quire<W>(pp##W##q,App##W##q,n); \
        int guard_##W##q = (std::isnan(pApp##W##q) || pApp##W##q <= 0.0) ? 1 : 0; \
        if(guard_##W##q){ pApp##W##q = 1.0; } \
        P##W alphap##W##q(rzp##W##q/pApp##W##q); \
        for(int i=0;i<n;i++){ xp##W##q[i]=xp##W##q[i]+alphap##W##q*pp##W##q[i]; rp##W##q[i]=rp##W##q[i]-alphap##W##q*App##W##q[i]; } \
        for(int i=0;i<n;i++) zp##W##q[i]=P##W(1.0/diagA[i])*rp##W##q[i]; \
        double rzp##W##q2 = dot_p_quire<W>(rp##W##q,zp##W##q,n); \
        P##W betap##W##q(rzp##W##q2/rzp##W##q); rzp##W##q=rzp##W##q2; \
        for(int i=0;i<n;i++) pp##W##q[i]=zp##W##q[i]+betap##W##q*pp##W##q[i]; \
        double resp##W##q = sqrt(dot_p_quire<W>(rp##W##q,rp##W##q,n)); \
        \
        matvec_p<W>(A,pp##W##n,App##W##n); \
        double pApp##W##n = dot_p_naive<W>(pp##W##n,App##W##n,n); \
        int guard_##W##n = (std::isnan(pApp##W##n) || pApp##W##n <= 0.0) ? 1 : 0; \
        if(guard_##W##n){ pApp##W##n = 1.0; } \
        P##W alphap##W##n(rzp##W##n/pApp##W##n); \
        for(int i=0;i<n;i++){ xp##W##n[i]=xp##W##n[i]+alphap##W##n*pp##W##n[i]; rp##W##n[i]=rp##W##n[i]-alphap##W##n*App##W##n[i]; } \
        for(int i=0;i<n;i++) zp##W##n[i]=P##W(1.0/diagA[i])*rp##W##n[i]; \
        double rzp##W##n2 = dot_p_naive<W>(rp##W##n,zp##W##n,n); \
        P##W betap##W##n(rzp##W##n2/rzp##W##n); rzp##W##n=rzp##W##n2; \
        for(int i=0;i<n;i++) pp##W##n[i]=zp##W##n[i]+betap##W##n*pp##W##n[i]; \
        double resp##W##n = sqrt(dot_p_naive<W>(rp##W##n,rp##W##n,n));

        STEP_POSIT(8)
        STEP_POSIT(16)
        STEP_POSIT(32)
        STEP_POSIT(64)

        double se32q=0, se32n=0;
        for(int i=0;i<n;i++){ double e=double(xp32q[i])-x_true[i]; se32q+=e*e; }
        for(int i=0;i<n;i++){ double e=double(xp32n[i])-x_true[i]; se32n+=e*e; }
        double solerr_p32q = sqrt(se32q)/x_true_norm;
        double solerr_p32n = sqrt(se32n)/x_true_norm;

        fprintf(log,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e "
                    "%.10e %.10e %.10e %.10e %.10e "
                    "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e "
                    "%.10e %.10e %.10e "
                    "%d %d %d %d %d %d %d %d\n",
                iter, resd, resf, resffma, resp8q, resp8n, resp16q, resp16n, resp32q, resp32n, resp64q, resp64n,
                solerr_d, solerr_f, solerr_ffma, solerr_p32q, solerr_p32n,
                pApd, pApp8q, pApp8n, pApp16q, pApp16n, pApp32q, pApp32n, pApp64q, pApp64n, abs_dot_d, cond_xy_d,
                rzd_used, rzp32q_used, rzp32n_used,
                guard_8q, guard_8n, guard_16q, guard_16n, guard_32q, guard_32n, guard_64q, guard_64n);

        if(conv_iter < 0 && (resd/bnorm) < 1e-10) conv_iter = iter;
        if(conv_iter_p32q < 0 && (resp32q/bnorm) < 1e-10) conv_iter_p32q = iter;
        if(conv_iter_p32n < 0 && (resp32n/bnorm) < 1e-10) conv_iter_p32n = iter;
        if(conv_iter_p32q_loose < 0 && (resp32q/bnorm) < 1e-6) conv_iter_p32q_loose = iter;
        if(conv_iter_p32n_loose < 0 && (resp32n/bnorm) < 1e-6) conv_iter_p32n_loose = iter;
        if(conv_iter >= 0 && iter > conv_iter + 20) break;
    }
    fprintf(log,"conv_iter_double=%d\n", conv_iter);
    fprintf(log,"conv_iter_p32q=%d\n", conv_iter_p32q);
    fprintf(log,"conv_iter_p32n=%d\n", conv_iter_p32n);
    fprintf(log,"conv_iter_p32q_loose_1e-6=%d\n", conv_iter_p32q_loose);
    fprintf(log,"conv_iter_p32n_loose_1e-6=%d\n", conv_iter_p32n_loose);
    fclose(log);
    printf("Done. Log: %s\n", log_path);
}
