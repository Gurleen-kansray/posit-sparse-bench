import glob, re
from collections import defaultdict
import statistics as stats
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.edgecolor'] = '#333333'
mpl.rcParams['axes.linewidth'] = 0.8

files = sorted(glob.glob("results/refinement_v2/*_seed*.log"))
data = defaultdict(list)

for f in files:
    with open(f) as fh:
        lines = fh.read().strip().split("\n")
    m = re.match(r"matrix=\S+/(\w+)\.mtx", lines[0])
    matname = m.group(1)
    vals = lines[2].split()
    solerr_before, solerr_after_double, solerr_after_posit, imp_double, imp_posit = map(float, vals)
    data[matname].append(imp_posit)

rows = []
for mat, imps in data.items():
    n = len(imps)
    mean = stats.mean(imps)
    std = stats.stdev(imps) if n > 1 else 0.0
    rows.append((mat, n, mean, std))

rows.sort(key=lambda r: r[2])

with open("results/refinement_v2/refinement_summary_FINAL.csv", "w") as out:
    out.write("matrix,n_seeds,imp_posit_mean,imp_posit_std\n")
    for mat, n, mean, std in rows:
        out.write(f"{mat},{n},{mean:.4f},{std:.4f}\n")

labels = [f"{mat} (n={n})" for mat, n, _, _ in rows]
means = [r[2] for r in rows]
stds = [r[3] for r in rows]

fig, ax = plt.subplots(figsize=(9, 6.5))
bar_color = "#3d5a80"
err_color = "#1b1b1b"

y_pos = range(len(labels))
ax.barh(y_pos, means, xerr=stds, color=bar_color, edgecolor="none",
        error_kw=dict(ecolor=err_color, elinewidth=1.1, capsize=3, capthick=1.1),
        height=0.62, zorder=3)

ax.set_yticks(y_pos)
ax.set_yticklabels(labels)
ax.set_xscale("log")
ax.set_xlim(0.8, max(m + s for m, s in zip(means, stds)) * 6)
ax.set_xlabel("Solution error improvement factor (log scale)", fontsize=11)
ax.set_title("Iterative refinement improvement factor\n(mean ± std, posit32 correction solve)",
             fontsize=12, pad=12)

ax.grid(axis="x", which="major", linestyle="-", linewidth=0.4, color="#cccccc", zorder=0)
ax.grid(axis="x", which="minor", linestyle=":", linewidth=0.3, color="#e5e5e5", zorder=0)
ax.set_axisbelow(True)

for i, (mean, std) in enumerate(zip(means, stds)):
    label_x = (mean + std) * 1.25
    ax.text(label_x, i, f"{mean:.2f}\u00d7", va="center", fontsize=9.5, color="#1b1b1b")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("fig1_refinement_hbar.png", dpi=300, bbox_inches="tight")
print("Wrote fig1_refinement_hbar.png (v3, clean labels)")
