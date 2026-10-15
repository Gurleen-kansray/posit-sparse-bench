import re, glob, numpy as np
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = {"double":6,"p32q":9,"p32n":10}
E = {k:[] for k in S}; q=[]; n=[]
for f in sorted(glob.glob("nasa4704_seed*.log")):
    if not re.search(r"seed\d+\.log$", f): continue
    D = np.loadtxt(f, skiprows=2, max_rows=4704, usecols=range(11))
    for k,c in S.items(): E[k].append(D[:,c])
    c = dict(l.strip().split("=") for l in open(f) if l.startswith("conv_iter_"))
    q.append(int(c["conv_iter_p32q_loose_1e-6"])); n.append(int(c["conv_iter_p32n_loose_1e-6"]))
q=np.array(q); n=np.array(n); d=n-q
print("Wilcoxon naive vs quire iters:", stats.wilcoxon(n,q))
print("sign test (quire fewer):", stats.binomtest(int((d>0).sum()), int((d!=0).sum()), 0.5).pvalue)
print("diff as pct of naive: median {:.2f}".format(float(np.median(100*d/n))))
print("seeds with |diff|>100:", [(i+1,int(q[i]),int(n[i])) for i in range(20) if abs(d[i])>100])
for k in S:
    M=np.array(E[k]); print(k,"median err at iters 500/1000/2000/3000/4703:",["{:.2e}".format(float(np.median(M[:,i]))) for i in (500,1000,2000,3000,4703)])
it=np.arange(4704); plt.figure(figsize=(6,4))
for k,l in (("double","double"),("p32q","posit32+quire"),("p32n","posit32 naive")):
    M=np.array(E[k]); plt.plot(it,np.median(M,0),label=l); plt.fill_between(it,M.min(0),M.max(0),alpha=0.15)
plt.yscale("log"); plt.axvline(2000,ls="--",c="k",lw=0.8); plt.xlabel("iteration"); plt.ylabel("relative solution error")
plt.title("nasa4704, 20 seeds (median, min-max)"); plt.legend(); plt.tight_layout(); plt.savefig("nasa4704_maxn_curves.png",dpi=200)
