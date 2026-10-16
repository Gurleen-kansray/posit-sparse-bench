import sys, glob, re, statistics as st, os
mode = sys.argv[1] if len(sys.argv) > 1 else "conv"   # conv | last
d = "results/newmat_q512"
mats = sorted({os.path.basename(f).rsplit("_seed",1)[0] for f in glob.glob(d+"/*_seed*.log") if not f.endswith(".tmp")})
print("matrix,seeds,mean,sd,median")
for m in mats:
    r = []
    for f in sorted(glob.glob(f"{d}/{m}_seed*.log")):
        if f.endswith(".tmp"): continue
        L = open(f).read().splitlines()
        hdr = L[1].split()
        rows = [l.split() for l in L if l[:1].isdigit()]
        ci = int(next(l for l in L if l.startswith("conv_iter_double=")).split("=")[1])
        row = rows[-1] if mode == "last" else next((x for x in rows if int(x[0]) == ci), rows[-1])
        q, n = float(row[hdr.index("solerr_p32q")]), float(row[hdr.index("solerr_p32n")])
        r.append(q / n)
    print(f"{m},{len(r)},{st.mean(r):.4f},{st.stdev(r):.4f},{st.median(r):.4f}")
