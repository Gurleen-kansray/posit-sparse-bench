#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
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

void matvec(const MTX& A,const double* x,double* y){
    for(int i=0;i<A.n;i++) y[i]=0;
    for(size_t k=0;k<A.row.size();k++) y[A.row[k]]+=A.val[k]*x[A.col[k]];
}

template<size_t N, size_t ES>
void compute_pAp(const std::vector<double>& p,const std::vector<double>& Ap,int n,
                 double& pAp_quire,double& pAp_naive){
    quire<N,ES,2> q=0;
    for(int i=0;i<n;i++){ posit<N,ES> pa(p[i]),pb(Ap[i]); q+=quire_mul(pa,pb); }
    posit<N,ES> pr; convert(q.to_value(),pr); pAp_quire=double(pr);
    posit<N,ES> pnaive=0;
    for(int i=0;i<n;i++){ posit<N,ES> pa(p[i]),pb(Ap[i]); pnaive=pnaive+pa*pb; }
    pAp_naive=double(pnaive);
}

void matrix_stats(const MTX& A, double& dyn_range, double& diag_ratio){
    // dynamic range = max/min of ALL absolute values (full matrix)
    double vmax=0,vmin=1e300;
    // diagonal ratio = max_diag/min_diag
    double dmax=0,dmin=1e300;
    for(size_t k=0;k<A.row.size();k++){
        double av=fabs(A.val[k]);
        if(av>0){ vmax=std::max(vmax,av); vmin=std::min(vmin,av); }
        if(A.row[k]==A.col[k]){
            dmax=std::max(dmax,av); dmin=std::min(dmin,av);
        }
    }
    // Use diagonal-only dynamic range — physically meaningful for FEM
    dyn_range=(dmin>0)?dmax/dmin:0;
    diag_ratio=(dmin>0)?dmax/dmin:0;
}

int main(int argc,char* argv[]){
    if(argc<2){ printf("Usage: %s matrix.mtx [iters]\n",argv[0]); return 1; }
    int max_iters=(argc>=3)?atoi(argv[2]):300;
    std::string path(argv[1]);
    std::string name=path.substr(path.find_last_of("/\\")+1);
    if(name.size()>4) name=name.substr(0,name.size()-4);

    MTX A=read_mtx(argv[1]);
    int n=A.n;
    double dyn_range,diag_ratio;
    matrix_stats(A,dyn_range,diag_ratio);

    printf("Matrix:%s n=%d log10_dynrange=%.1f diag_ratio=%.3e\n",
           name.c_str(),n,log10(dyn_range),diag_ratio);

    std::vector<double> diagA(n,1.0);
    for(int k=0;k<(int)A.row.size();k++)
        if(A.row[k]==A.col[k]) diagA[A.row[k]]=A.val[k];

    std::vector<double> x(n,0),r(n),p(n),Ap(n),z(n);
    for(int i=0;i<n;i++) r[i]=1.0;
    for(int i=0;i<n;i++) p[i]=z[i]=r[i]/diagA[i];
    double rz_d=0; for(int i=0;i<n;i++) rz_d+=r[i]*z[i];

    double mp8q=0,mp8n=0,mp16q=0,mp16n=0,mp32q=0,mp32n=0,mp64q=0,mp64n=0,pmax=0;

    for(int iter=0;iter<max_iters;iter++){
        matvec(A,p.data(),Ap.data());
        double pAp_d=0; for(int i=0;i<n;i++) pAp_d+=p[i]*Ap[i];
        pmax=std::max(pmax,fabs(pAp_d));
        if(fabs(pAp_d)>pmax*1e-6){
            double q8q,q8n,q16q,q16n,q32q,q32n,q64q,q64n;
            compute_pAp<8,0> (p,Ap,n,q8q, q8n);
            compute_pAp<16,1>(p,Ap,n,q16q,q16n);
            compute_pAp<32,2>(p,Ap,n,q32q,q32n);
            compute_pAp<64,2>(p,Ap,n,q64q,q64n);
            auto rel=[&](double e){ return fabs(pAp_d)>0?fabs(e-pAp_d)/fabs(pAp_d):0; };
            mp8q=std::max(mp8q,rel(q8q)); mp8n=std::max(mp8n,rel(q8n));
            mp16q=std::max(mp16q,rel(q16q)); mp16n=std::max(mp16n,rel(q16n));
            mp32q=std::max(mp32q,rel(q32q)); mp32n=std::max(mp32n,rel(q32n));
            mp64q=std::max(mp64q,rel(q64q)); mp64n=std::max(mp64n,rel(q64n));
        }
        double alpha=rz_d/pAp_d;
        for(int i=0;i<n;i++){x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];}
        for(int i=0;i<n;i++) z[i]=r[i]/diagA[i];
        double rz2=0; for(int i=0;i<n;i++) rz2+=r[i]*z[i];
        double rr=0; for(int i=0;i<n;i++) rr+=r[i]*r[i];
        if(sqrt(rr)<1e-10){ printf("Converged iter=%d\n",iter); break; }
        double beta=rz2/rz_d; rz_d=rz2;
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];
    }

    auto V=[](double e)->const char*{
        if(e>1e4)  return "OVERFLOW";
        if(e>1e-2) return "FAIL";
        if(e>1e-4) return "MARGINAL";
        return "PASS";
    };
    const char* mv="posit64";
    if(strcmp(V(mp8q),"PASS")==0)  mv="posit8";
    else if(strcmp(V(mp16q),"PASS")==0) mv="posit16";
    else if(strcmp(V(mp32q),"PASS")==0) mv="posit32";

    printf("posit8  quire=%.3e(%s) naive=%.3e\n",mp8q,V(mp8q),mp8n);
    printf("posit16 quire=%.3e(%s) naive=%.3e\n",mp16q,V(mp16q),mp16n);
    printf("posit32 quire=%.3e(%s) naive=%.3e\n",mp32q,V(mp32q),mp32n);
    printf("posit64 quire=%.3e(%s) naive=%.3e\n",mp64q,V(mp64q),mp64n);
    printf("MIN_VIABLE=%s\n",mv);
    printf("CSV:%s,%d,%.1f,%.3e,%.3e,%.3e,%.3e,%s\n",
           name.c_str(),n,log10(dyn_range),mp8q,mp16q,mp32q,mp64q,mv);
}
