import glob, re, os, numpy as np
from collections import defaultdict
R = defaultdict(list); tot = defaultdict(int)
for f in glob.glob("results/ladder_logs/seed_sweep/*_seed*.log"):
    m = re.match(r"(.+)_seed(\d+)\.log$", os.path.basename(f))
    if not m: continue
    mat = m.group(1); tot[mat] += 1
    c = {}
    for l in open(f, errors="ignore"):
        if l.startswith("conv_iter_"):
            k, v = l.strip().split("=")
            try: c[k] = int(v)
            except: pass
    q = c.get("conv_iter_p32q_loose_1e-6"); n = c.get("conv_iter_p32n_loose_1e-6")
    if q is None or n is None: continue
    R[mat].append((q, n))
out = ["matrix,seeds_total,seeds_with_counter,both_reached,quire_median,naive_median,quire_fewer,equal,quire_more,median_pct_diff"]
print("%-10s %5s %5s %5s %8s %8s  q<n =  q>n" % ("matrix","tot","ctr","both","q_med","n_med"))
for mat in sorted(tot):
    A = [(q, n) for q, n in R[mat] if q > 0 and n > 0]
    if not A:
        print("%-10s %5d %5d %5d   (none reached 1e-6)" % (mat, tot[mat], len(R[mat]), 0))
        out.append("%s,%d,%d,0,,,,,," % (mat, tot[mat], len(R[mat]))); continue
    q = np.array([a[0] for a in A], float); n = np.array([a[1] for a in A], float)
    row = (mat, tot[mat], len(R[mat]), len(A), np.median(q), np.median(n), int((q<n).sum()), int((q==n).sum()), int((q>n).sum()), np.median(100*(n-q)/n))
    print("%-10s %5d %5d %5d %8.1f %8.1f  %d / %d / %d  (%.2f%%)" % (row[:6] + row[6:9] + (row[9],)))
    out.append("%s,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%.2f" % row)
open("results/csv/iter_counts_all_matrices.csv", "w").write("\n".join(out) + "\n")
