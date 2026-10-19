import glob, re, os, numpy as np
from collections import defaultdict
SKIP = ("corrupted", "excluded", "BACKUP", "cond_probe", "rz_sanity", "maxn")
D = defaultdict(dict)   # (mat,seed) -> {src: (q,n)}
for f in glob.glob("results/**/*.log", recursive=True) + glob.glob("results/**/*.out", recursive=True):
    if any(s in f for s in SKIP): continue
    m = re.match(r"([A-Za-z0-9]+)_seed(\d+)", os.path.basename(f))
    if not m: continue
    q = n = None
    for l in open(f, errors="ignore"):
        if l.startswith("conv_iter_p32q_loose_1e-6="): q = int(l.split("=")[1])
        elif l.startswith("conv_iter_p32n_loose_1e-6="): n = int(l.split("=")[1])
    if q is not None and n is not None:
        D[(m.group(1), int(m.group(2)))][os.path.dirname(f).split("/")[-1]] = (q, n)
conf = [(k, v) for k, v in D.items() if len(set(v.values())) > 1]
print("pairs:", len(D), "| conflicting duplicates:", len(conf))
for k, v in conf[:10]: print("  CONFLICT", k, v)
tot = {"bcsstk03":50,"bcsstk14":30,"bcsstk36":30,"bcsstk37":30,"bcsstk38":30,"bodyy4":30,"mhd4800b":50,"nasa4704":30,"nasasrb":30,"nos2":40,"s3dkq4m2":25,"s3dkt3m2":25,"sts4098":30}
by = defaultdict(list)
for (mat, s), v in D.items(): by[mat].append(list(v.values())[0])
print("%-10s %4s %4s %5s %7s %7s  q<n/eq/q>n" % ("matrix","tot","got","both","q_med","n_med"))
for mat in sorted(tot):
    A = [(q, n) for q, n in by.get(mat, []) if q > 0 and n > 0]
    if not A: print("%-10s %4d %4d %5d" % (mat, tot[mat], len(by.get(mat, [])), 0)); continue
    q = np.array([a[0] for a in A], float); n = np.array([a[1] for a in A], float)
    print("%-10s %4d %4d %5d %7.1f %7.1f  %d/%d/%d" % (mat, tot[mat], len(by[mat]), len(A), np.median(q), np.median(n), (q<n).sum(), (q==n).sum(), (q>n).sum()))
