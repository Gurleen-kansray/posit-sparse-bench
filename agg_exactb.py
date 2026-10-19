import glob, os, statistics as st
d = "results/exactb_q512"
mats = sorted({os.path.basename(f).rsplit("_seed",1)[0] for f in glob.glob(d+"/*_seed*.log")})
print("matrix,seeds,metric,mean,sd,median")
for m in mats:
    cols = {"improvement_posit":[], "improvement_double":[], "after_double/after_posit":[]}
    for f in sorted(glob.glob(f"{d}/{m}_seed*.log")):
        L = open(f).read().splitlines()
        h, v = L[1].split(), L[2].split()
        g = lambda k: float(v[h.index(k)])
        cols["improvement_posit"].append(g("improvement_posit"))
        cols["improvement_double"].append(g("improvement_double"))
        cols["after_double/after_posit"].append(g("solerr_after_double")/g("solerr_after_posit"))
    for k, r in cols.items():
        print(f"{m},{len(r)},{k},{st.mean(r):.4f},{st.stdev(r):.4f},{st.median(r):.4f}")
