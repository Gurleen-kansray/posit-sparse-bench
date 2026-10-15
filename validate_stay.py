import glob, re, os, subprocess, numpy as np
def stay(r, t):
    above = np.where(~(r < t))[0]
    if len(above) == 0: return 0
    return int(above[-1]) + 1 if above[-1] + 1 < len(r) else -1
def stay2(r, t):
    r = np.where(np.isnan(r), np.inf, r)
    suf = np.maximum.accumulate(r[::-1])[::-1]
    ok = np.where(suf < t)[0]
    return int(ok[0]) if len(ok) else -1
S = [(np.array([2,2,.5,2,.5,.5]),1,4), (np.array([.5,2]),1,-1), (np.array([.5,.4]),1,0), (np.array([np.nan,.5]),1,1)]
for r,t,e in S: assert stay(r,t)==e==stay2(r,t), (r,t,e,stay(r,t),stay2(r,t))
print("synthetic cases: OK")
def load(f):
    L = open(f, errors="ignore").read().splitlines()
    h = L[1].split(); iq, inn = h.index("res_p32q"), h.index("res_p32n")
    rows = [l.split() for l in L[2:] if l[:1].isdigit()]
    return np.array([float(x[iq]) for x in rows]), np.array([float(x[inn]) for x in rows])
files = []
for m, ss in (("bcsstk03", range(1,11)), ("sts4098", range(1,6)), ("s3dkt3m2", range(1,4))):
    files += [(m, s, f"results/ladder_logs/seed_sweep/{m}_seed{s}.log") for s in ss]
files += [("nasa4704", s, f"results/maxn_all/nasa4704_seed{s}.log") for s in range(1,21)]
bad = tot = 0
for m, s, f in files:
    if not os.path.exists(f): print("missing", f); continue
    b = float(subprocess.check_output(["./bnorm_tool", f"data/matrices/{m}.mtx", "/tmp/x.log", str(s)]))
    for r in load(f):
        for t in (1e-6, 1e-5):
            a, c = stay(r/b, t), stay2(r/b, t); tot += 1
            if a != c: bad += 1; print("MISMATCH", f, t, a, c)
print("real-log comparisons:", tot, "mismatches:", bad)
