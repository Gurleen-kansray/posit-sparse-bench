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

void matvec(const MTX& A,const double* x,double* y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(int k=0;k<(int)A.row.size();k++) y[A.row[k]]+=A.val[k]*x[A.col[k]];
}

template<size_t N, size_t ES>
void compute_pAp(const std::vector<double>& p, const std::vector<double>& Ap, int n,
                 double& pAp_quire, double& pAp_naive){
    quire<N,ES,2> q=0;
    for(int i=0;i<n;i++){ posit<N,ES> pa(p[i]),pb(Ap[i]); q+=quire_mul(pa,pb); }
    posit<N,ES> pr; convert(q.to_value(),pr); pAp_quire=double(pr);
    posit<N,ES> pnaive=0;
    for(int i=0;i<n;i++){ posit<N,ES> pa(p[i]),pb(Ap[i]); pnaive=pnaive+pa*pb; }
    pAp_naive=double(pnaive);
}

template<size_t N, size_t ES>
void compute_rz(const std::vector<double>& r, const std::vector<double>& z, int n,
                double& rz_quire, double& rz_naive){
    quire<N,ES,2> q=0;
    for(int i=0;i<n;i++){ posit<N,ES> pa(r[i]),pb(z[i]); q+=quire_mul(pa,pb); }
    posit<N,ES> pr; convert(q.to_value(),pr); rz_quire=double(pr);
    posit<N,ES> pnaive=0;
    for(int i=0;i<n;i++){ posit<N,ES> pa(r[i]),pb(z[i]); pnaive=pnaive+pa*pb; }
    rz_naive=double(pnaive);
}
int main(int argc, char* argv[]){
    if(argc<3){ printf("usage: generic_ladder <matrix.mtx> <logfile>\n"); return 1; }
    const char* mtx_path = argv[1];
    const char* log_path = argv[2];

    MTX A=read_mtx(mtx_path);
    int n=A.n;

    double vmin=1e300, vmax=0;
    for(double v: A.val){ double av=fabs(v); if(av>0&&av<vmin) vmin=av; if(av>vmax) vmax=av; }
    printf("Loaded %s: n=%d nnz=%zu\n", mtx_path, n, A.row.size());
    printf("Value range: min_abs=%.3e max_abs=%.3e ratio=%.3e\n", vmin, vmax, vmax/vmin);

    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    std::vector<double> x(n,0),r(n,1),p(n),Ap(n),z(n);
    for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
    for(int i=0;i<n;i++) p[i]=z[i];
    double rz_d=0; for(int i=0;i<n;i++) rz_d+=r[i]*z[i];

    FILE* log=fopen(log_path,"w");
    fprintf(log,"matrix=%s n=%d val_range=%.3e/%.3e\n", mtx_path, n, vmin, vmax);

    for(int iter=0;iter<300;iter++){
        matvec(A,p.data(),Ap.data());
        double pAp_d=0; for(int i=0;i<n;i++) pAp_d+=p[i]*Ap[i];
        double p8q,p8n,p16q,p16n,p32q,p32n,p64q,p64n;
        double p8q_e2,p8n_e2,p16q_e2,p16n_e2;
        compute_pAp<8,0> (p,Ap,n,p8q,   p8n);
        compute_pAp<8,2> (p,Ap,n,p8q_e2, p8n_e2);
        compute_pAp<16,1>(p,Ap,n,p16q,  p16n);
        compute_pAp<16,2>(p,Ap,n,p16q_e2,p16n_e2);
        compute_pAp<32,2>(p,Ap,n,p32q,  p32n);
        compute_pAp<64,2>(p,Ap,n,p64q,  p64n);
        double rz32q, rz32n;
        compute_rz<32,2>(r,z,n,rz32q,rz32n);
        double alpha_d    = rz_d/pAp_d;
        double alpha_p32q = (p32q != 0.0) ? rz_d/p32q : 0.0;
        double alpha_p32n = (p32n != 0.0) ? rz_d/p32n : 0.0;
        double alpha_full_q = (p32q != 0.0) ? rz32q/p32q : 0.0;
        double alpha_full_n = (p32n != 0.0) ? rz32n/p32n : 0.0;
        fprintf(log,"iter=%d pAp_d=%.10e p8q=%.10e p8n=%.10e p8q_e2=%.10e p8n_e2=%.10e p16q=%.10e p16n=%.10e p16q_e2=%.10e p16n_e2=%.10e p32q=%.10e p32n=%.10e p64q=%.10e p64n=%.10e alpha_d=%.10e alpha_p32q=%.10e alpha_p32n=%.10e rz32q=%.10e rz32n=%.10e alpha_full_q=%.10e alpha_full_n=%.10e\n",
                iter,pAp_d,p8q,p8n,p8q_e2,p8n_e2,p16q,p16n,p16q_e2,p16n_e2,p32q,p32n,p64q,p64n,alpha_d,alpha_p32q,alpha_p32n,rz32q,rz32n,alpha_full_q,alpha_full_n);
        double alpha=alpha_d;
        for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
        for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
        double rz2=0; for(int i=0;i<n;i++) rz2+=r[i]*z[i];
        double beta=rz2/rz_d; rz_d=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
        double true_rr=0; for(int i=0;i<n;i++) true_rr+=r[i]*r[i];
        fprintf(log,"  residual=%.3e\n",sqrt(true_rr));
    }
    fclose(log);
    printf("Done. Log written to %s\n", log_path);
}
