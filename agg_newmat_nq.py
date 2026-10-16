import glob, os, csv, statistics as st
d = "results/newmat_q512"
mats = sorted({os.path.basename(f).rsplit("_seed",1)[0] for f in glob.glob(d+"/*_seed*.log")})
out = [["matrix","seeds","median_naive_over_quire_last","mean","sd","quire_better_seeds","median_at_double_conv"]]
for m in mats:
    last, conv = [], []
    for f in sorted(glob.glob(f"{d}/{m}_seed*.log")):
        L = open(f).read().splitlines(); h = L[1].split()
        rows = [l.split() for l in L if l[:1].isdigit()]
        ci = int(next(l for l in L if l.startswith("conv_iter_double=")).split("=")[1])
        iq, inn = h.index("solerr_p32q"), h.index("solerr_p32n")
        r = rows[-1]; last.append(float(r[inn]) / float(r[iq]))
        c = next((x for x in rows if int(x[0]) == ci), rows[-1]); conv.append(float(c[inn]) / float(c[iq]))
    out.append([m, len(last), round(st.median(last),4), round(st.mean(last),4), round(st.stdev(last),4), sum(x > 1 for x in last), round(st.median(conv),4)])
csv.writer(open("results/csv/newmat_naive_over_quire.csv","w")).writerows(out)
for r in out: print(*r, sep="\t")
