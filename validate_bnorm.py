import glob, re, os, subprocess, numpy as np
from collections import defaultdict
SKIP = ("corrupted","excluded","BACKUP","cond_probe","rz_sanity","maxn")
C = {}
for f in glob.glob("results/**/*.log", recursive=True):
    if any(s in f for s in SKIP): continue
    m = re.match(r"([A-Za-z0-9]+)_seed(\d+)", os.path.basename(f))
    if not m: continue
    q = n = None
    for l in open(f, errors="ignore"):
        if l.startswith("conv_iter_p32q_loose_1e-6="): q = int(l.split("=")[1])
        elif l.startswith("conv_iter_p32n_loose_1e-6="): n = int(l.split("=")[1])
    if q is not None and n is not None: C[(m.group(1), int(m.group(2)))] = (q, n)
tally = defaultdict(lambda: [0, 0]); bad = []
for f in sorted(glob.glob("results/ladder_logs/seed_sweep/*_seed*.log")):
    m = re.match(r"(.+)_seed(\d+)\.log$", os.path.basename(f))
    if not m or (m.group(1), int(m.group(2))) not in C: continue
    mat, s = m.group(1), int(m.group(2))
    L = open(f, errors="ignore").read().splitlines()
    h = L[1].split(); iq, inn = h.index("res_p32q"), h.index("res_p32n")
    rows = [l.split() for l in L[2:] if l[:1].isdigit()]
    b = float(subprocess.check_output(["./bnorm_tool", f"data/matrices/{mat}.mtx", "/tmp/x.log", str(s)]))
    def first(col):
        for r in rows:
            if float(r[col]) / b < 1e-6: return int(r[0])
        return -1
    got = (first(iq), first(inn)); ok = got == C[(mat, s)]
    tally[mat][0] += ok; tally[mat][1] += 1
    if not ok: bad.append((mat, s, got, C[(mat, s)]))
for mat in sorted(tally): print("%-10s match %d / %d" % (mat, *tally[mat]))
print("mismatches:", len(bad)); [print(" ", x) for x in bad[:10]]
