import re, glob, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = {"double":6, "p32q":9, "p32n":10}          # solerr columns
rows = []
for f in sorted(glob.glob("nasa4704_seed*.log")):
    m = re.search(r"seed(\d+)\.log$", f)
    if not m: continue                           # skips *_check copies
    s = int(m.group(1))
    D = np.loadtxt(f, skiprows=2, max_rows=4704, usecols=range(11))
    it = D[:,0].astype(int)
    conv = {}
    for l in open(f):
        if l.startswith("conv_iter_"):
            k,v = l.strip().split("="); conv[k]=int(v)
    r = {"seed":s, "last_iter":int(it[-1])}
    for name,c in V.items():
        e = D[:,c]
        r[f"{name}_min"] = e.min(); r[f"{name}_min_iter"] = int(it[e.argmin()])
        r[f"{name}_e2000"] = e[it==2000][0]; r[f"{name}_elast"] = e[-1]
    r["q_1e6"] = conv["conv_iter_p32q_loose_1e-6"]
    r["n_1e6"] = conv["conv_iter_p32n_loose_1e-6"]
    rows.append(r)

keys = list(rows[0].keys())
with open("maxn_summary.csv","w") as o:
    o.write(",".join(keys)+"\n")
    for r in rows: o.write(",".join(str(r[k]) for k in keys)+"\n")

A = lambda k: np.array([r[k] for r in rows], float)
print("seeds:", len(rows))
for name in V:
    print(f"{name}: min err median {np.median(A(name+'_min')):.3e} (min at iter median {np.median(A(name+'_min_iter')):.0f}); "
          f"err@2000 median {np.median(A(name+'_e2000')):.3e}; err@{rows[0]['last_iter']} median {np.median(A(name+'_elast')):.3e}")
print("double: min at last iter in", int(sum(r['double_min_iter']==r['last_iter'] for r in rows)), "of", len(rows))
print("double err ratio e2000/elast median:", np.median(A('double_e2000')/A('double_elast')))
q,n = A("q_1e6"), A("n_1e6")
print("iters to 1e-6: quire median", np.median(q), "naive median", np.median(n))
print("quire fewer / equal / more:", int((q<n).sum()), int((q==n).sum()), int((q>n).sum()))
print("paired diff (n-q) median", np.median(n-q), "mean", (n-q).mean())
print("seeds with 1e-6 iter > 2000: quire", int((q>2000).sum()), "naive", int((n>2000).sum()))
print("solerr ratio naive/quire at last iter: median", np.median(A('p32n_elast')/A('p32q_elast')),
      "| at min:", np.median(A('p32n_min')/A('p32q_min')))

fig,ax = plt.subplots(1,2,figsize=(10,4))
ax[0].boxplot([A('double_elast'),A('p32q_elast'),A('p32n_elast')], labels=["double","posit32+quire","posit32 naive"])
ax[0].set_yscale("log"); ax[0].set_ylabel("rel. solution error at iter n-1"); ax[0].set_title("nasa4704, 20 seeds")
ax[1].boxplot([q,n], labels=["quire","naive"]); ax[1].set_ylabel("iterations to 1e-6 residual")
plt.tight_layout(); plt.savefig("nasa4704_maxn_box.png", dpi=200)
