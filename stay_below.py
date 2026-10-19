import subprocess, numpy as np, glob, re
def stay(r, t):
    above = np.where(r >= t)[0]
    if len(above) == 0: return 0
    return int(above[-1]) + 1 if above[-1] + 1 < len(r) else -1   # -1 = never stays below
for t in (1e-6, 1e-5):
    Q, N = [], []
    for s in range(1, 21):
        f = f"results/maxn_all/nasa4704_seed{s}.log"
        b = float(subprocess.check_output(["./bnorm_tool", "data/matrices/nasa4704.mtx", "/tmp/x.log", str(s)]))
        D = np.loadtxt(f, skiprows=2, max_rows=4704, usecols=range(11))
        Q.append(stay(D[:, 4] / b, t)); N.append(stay(D[:, 5] / b, t))
    Q = np.array(Q); N = np.array(N)
    ok = (Q > 0) & (N > 0)
    print(f"tol {t:g}: quire never-stays {int((Q<0).sum())}, naive never-stays {int((N<0).sum())}; both stay: {int(ok.sum())}")
    if ok.any():
        d = (N - Q)[ok]
        print(f"  median quire {np.median(Q[ok]):.0f}, naive {np.median(N[ok]):.0f}; quire fewer/equal/more: {int((d>0).sum())}/{int((d==0).sum())}/{int((d<0).sum())}; median diff {np.median(d):.0f}")
    print("  quire:", Q.tolist()); print("  naive:", N.tolist())
