import csv
from fractions import Fraction
import statistics as stats

# Re-read the original pairs to compute an exact (rational) reference dot product
def read_pairs(path):
    pairs = []
    with open(path) as f:
        total = int(f.readline())
        for _ in range(total):
            tag, n = f.readline().split()
            n = int(n)
            x = [float(v) for v in f.readline().split()]
            w = [float(v) for v in f.readline().split()]
            pairs.append((tag, x, w))
    return pairs

pairs = read_pairs("dotpairs.txt")

exact = []
for tag, x, w in pairs:
    acc = Fraction(0)
    for xi, wi in zip(x, w):
        acc += Fraction(xi) * Fraction(wi)   # exact: float64 values are exact dyadic rationals
    exact.append(float(acc))  # cast only for reporting; comparison done in relative-error space below

rows = []
with open("harness_results.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)

assert len(rows) == len(exact)

methods = ["float32", "double64", "posit32", "posit32_quire"]
errs = {m: {"L1": [], "L2": [], "ALL": []} for m in methods}

for row, ref in zip(rows, exact):
    tag = row["layer"]
    for m in methods:
        val = float(row[m])
        rel = abs(val - ref) / max(abs(ref), 1e-300)
        errs[m][tag].append(rel)
        errs[m]["ALL"].append(rel)

def summarize(vals):
    if not vals:
        return (0, 0, 0)
    return (stats.mean(vals), stats.median(vals), max(vals))

print(f"{'method':16s} {'scope':5s} {'mean_rel_err':>14s} {'median_rel_err':>16s} {'max_rel_err':>14s}")
for m in methods:
    for scope in ["L1", "L2", "ALL"]:
        mean_e, med_e, max_e = summarize(errs[m][scope])
        print(f"{m:16s} {scope:5s} {mean_e:14.3e} {med_e:16.3e} {max_e:14.3e}")
    print()

# Cross-check vs prior synthetic-Gaussian ML result and FEM result framing:
# report ratio of posit32_quire error to double64 error (does quire match double for these real, bounded values?)
q_all = errs["posit32_quire"]["ALL"]
d_all = errs["double64"]["ALL"]
p_all = errs["posit32"]["ALL"]
fl_all = errs["float32"]["ALL"]

print("posit32+quire mean err / double64 mean err :", stats.mean(q_all) / stats.mean(d_all))
print("posit32 (no quire) mean err / quire mean err:", stats.mean(p_all) / stats.mean(q_all))
print("float32 mean err / posit32+quire mean err   :", stats.mean(fl_all) / stats.mean(q_all))
