import re, glob, os, csv
import statistics as stats

LOG_DIR = "results/ladder_logs/seed_dist"
OUT_CSV = "results/csv/divergence_distribution.csv"
pattern = re.compile(r"iter=(\d+) pAp_d=([\d.eE+-]+).*?p32q=([\d.eE+-]+) p32n=([\d.eE+-]+)")

def divergence_iter_for_log(path):
    with open(path) as f:
        lines = f.readlines()
    iters = []
    for line in lines:
        m = pattern.search(line)
        if not m: continue
        it, pAp_d, p32q, p32n = int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))
        if abs(pAp_d) > 1e-30:
            err_q, err_n = abs(p32q-pAp_d)/abs(pAp_d), abs(p32n-pAp_d)/abs(pAp_d)
        else:
            err_q, err_n = abs(p32q-pAp_d), abs(p32n-pAp_d)
        iters.append((it, err_q, err_n))
    if not iters: return None
    for i in range(len(iters)):
        it, eq, en = iters[i]
        if en > 0 and eq > 0 and en/eq > 10:
            window = iters[i:i+5]
            if all(w[2]/w[1] > 10 if w[1] > 0 else False for w in window):
                return it
    return None

by_matrix = {}
for path in glob.glob(os.path.join(LOG_DIR, "*_seed*.log")):
    base = os.path.basename(path)
    m = re.match(r"(.+)_seed(\d+)\.log", base)
    if not m: continue
    matrix, seed = m.group(1), int(m.group(2))
    by_matrix.setdefault(matrix, []).append((seed, divergence_iter_for_log(path)))

rows = []
for matrix, trials in sorted(by_matrix.items()):
    vals = [d for _, d in trials if d is not None]
    n_total = len(trials)
    row = {"matrix": matrix, "n_trials": n_total, "n_diverged": len(vals),
           "frac_diverged": len(vals)/n_total if n_total else 0.0,
           "mean_div_iter": stats.mean(vals) if vals else "",
           "median_div_iter": stats.median(vals) if vals else "",
           "stdev_div_iter": stats.pstdev(vals) if len(vals) > 1 else "",
           "min_div_iter": min(vals) if vals else "", "max_div_iter": max(vals) if vals else ""}
    rows.append(row)
    print(f"{matrix:<12} n_trials={n_total:<5} diverged={len(vals):<5} mean={row['mean_div_iter']!s:<8} stdev={row['stdev_div_iter']!s:<8} range=[{row['min_div_iter']},{row['max_div_iter']}]")

if rows:
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    print(f"\nWrote {OUT_CSV}")
