import glob, re
from collections import defaultdict
import statistics as stats

files = sorted(glob.glob("results/refinement_v2/*_seed*.log"))
data = defaultdict(list)

for f in files:
    with open(f) as fh:
        lines = fh.read().strip().split("\n")
    m = re.match(r"matrix=\S+/(\w+)\.mtx", lines[0])
    matname = m.group(1)
    vals = lines[2].split()
    solerr_before, solerr_after_double, solerr_after_posit, imp_double, imp_posit = map(float, vals)
    data[matname].append((solerr_before, solerr_after_double, solerr_after_posit, imp_double, imp_posit))

print(f"{'matrix':<12}{'n_seeds':<9}{'imp_posit_mean':<16}{'imp_posit_std':<15}{'n_worse(<1.0)':<15}{'imp_double_mean':<16}")
for mat, rows in data.items():
    imp_posit = [r[4] for r in rows]
    imp_double = [r[3] for r in rows]
    n_worse = sum(1 for x in imp_posit if x < 1.0)
    mean_p = stats.mean(imp_posit)
    std_p = stats.stdev(imp_posit) if len(imp_posit) > 1 else 0.0
    mean_d = stats.mean(imp_double)
    print(f"{mat:<12}{len(rows):<9}{mean_p:<16.4f}{std_p:<15.4f}{n_worse:<15}{mean_d:<16.4f}")
