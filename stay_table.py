import glob, re, os, subprocess, csv, numpy as np
from collections import defaultdict
def stay(r, t):
    above = np.where(~(r < t))[0]
    if len(above) == 0: return 0
    return int(above[-1]) + 1 if above[-1] + 1 < len(r) else -1
R = defaultdict(list)
for f in sorted(glob.glob("results/ladder_logs/seed_sweep/*_seed*.log")):
    m = re.match(r"(.+)_seed(\d+)\.log$", os.path.basename(f))
    if not m: continue
    mat, s = m.group(1), int(m.group(2))
    L = open(f, errors="ignore").read().splitlines()
    if len(L) < 3: continue
    h = L[1].split(); iq, inn = h.index("res_p32q"), h.index("res_p32n")
    rows = [l.split() for l in L[2:] if l[:1].isdigit()]
    b = float(subprocess.check_output(["./bnorm_tool", f"data/matrices/{mat}.mtx", "/tmp/x.log", str(s)]))
    R[mat].append((np.array([float(r[iq]) for r in rows]) / b, np.array([float(r[inn]) for r in rows]) / b))
out = [["tol","matrix","seeds","both_stay","never_stay","quire_fewer","equal","quire_more","median_paired_diff_iters","median_paired_diff_pct_of_naive"]]
for t in (1e-6, 1e-5):
    for mat in sorted(R):
        Q = np.array([stay(a, t) for a, _ in R[mat]]); N = np.array([stay(b, t) for _, b in R[mat]])
        ok = (Q >= 0) & (N >= 0)
        if not ok.any(): out.append([t, mat, len(Q), 0, int((~ok).sum())] + [""]*5); continue
        d = (N - Q)[ok]
        out.append([t, mat, len(Q), int(ok.sum()), int((~ok).sum()), int((d>0).sum()), int((d==0).sum()), int((d<0).sum()), float(np.median(d)), round(float(np.median(100*d/np.maximum(N[ok],1))), 2)])
csv.writer(open("results/csv/iter_stay_below.csv", "w")).writerows(out)
for r in out: print(*r, sep="\t")
