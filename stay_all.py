import glob, re, os, subprocess, numpy as np
from collections import defaultdict
from scipy import stats
def stay(r, t):
    above = np.where(~(r < t))[0]
    if len(above) == 0: return 0
    return int(above[-1]) + 1 if above[-1] + 1 < len(r) else -1
R = defaultdict(list)
for f in sorted(glob.glob("results/ladder_logs/seed_sweep/*_seed*.log")):
    m = re.match(r"(.+)_seed(\d+)\.log$", os.path.basename(f))
    if not m: continue
    mat, s = m.group(1), int(m.group(2))
    lines = open(f, errors="ignore").read().splitlines()
    if len(lines) < 3: continue
    h = lines[1].split(); iq, inn = h.index("res_p32q"), h.index("res_p32n")
    rows = [l.split() for l in lines[2:] if l[:1].isdigit()]
    q = np.array([float(r[iq]) for r in rows]); n = np.array([float(r[inn]) for r in rows])
    try: b = float(subprocess.check_output(["./bnorm_tool", f"data/matrices/{mat}.mtx", "/tmp/x.log", str(s)]))
    except Exception: continue
    R[mat].append((s, len(rows), q / b, n / b))
for t in (1e-6, 1e-5):
    print(f"\ntol {t:g}\n%-10s %4s %5s %6s %7s %7s  q<n/eq/q>n  p_sign" % ("matrix","seeds","both","never","q_med","n_med"))
    for mat in sorted(R):
        Q = np.array([stay(r[2], t) for r in R[mat]]); N = np.array([stay(r[3], t) for r in R[mat]])
        ok = (Q >= 0) & (N >= 0); never = int((~ok).sum())
        if not ok.any(): print("%-10s %4d %5d %6d" % (mat, len(Q), 0, never)); continue
        d = (N - Q)[ok]; k = int((d > 0).sum()); m_ = int((d != 0).sum())
        p = stats.binomtest(k, m_, 0.5).pvalue if m_ else 1.0
        print("%-10s %4d %5d %6d %7.0f %7.0f  %d/%d/%d  %.3f" % (mat, len(Q), int(ok.sum()), never, np.median(Q[ok]), np.median(N[ok]), k, int((d==0).sum()), int((d<0).sum()), p))
