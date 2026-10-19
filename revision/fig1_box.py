import glob, re
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = []
for f in sorted(glob.glob("results/refinement_v2/*_seed*.log")):
    L = open(f).read().strip().split("\n")
    m = re.match(r"matrix=\S+/(\w+)\.mtx", L[0]).group(1)
    seed = int(re.search(r"seed(\d+)\.log", f).group(1))
    v = L[2].split()
    rows.append((m, seed, float(v[4])))
df = pd.DataFrame(rows, columns=["matrix", "seed", "improv"])
df.to_csv("revision/tables/refinement_v2_perseed.csv", index=False)

s = df.groupby("matrix").improv.agg(n="count", median="median", mean="mean",
                                    worse=lambda x: (x < 1).sum())
print(s.sort_values("median", ascending=False).round(3).to_string())
print("total runs:", s.n.sum(), " improved:", int((df.improv >= 1).sum()))

order = s.sort_values("median").index.tolist()
fig, ax = plt.subplots(figsize=(8, 5))
ax.boxplot([df.loc[df.matrix == m, "improv"] for m in order], vert=False,
           tick_labels=[f"{m} (n={s.loc[m,'n']})" for m in order], showfliers=True)
ax.set_xscale("log"); ax.axvline(1, ls="--", lw=0.8)
ax.set_xlabel("solution error improvement factor (before / after), log scale")
ax.set_title("Iterative refinement, posit32 correction solve")
plt.tight_layout(); plt.savefig("revision/figs/fig1_refinement_box_all13.png", dpi=300)
