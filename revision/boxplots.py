import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
# usage: python boxplots.py file.csv matrix_col value_col out.png "title" [log]
f, mcol, vcol, out, title = sys.argv[1:6]
df = pd.read_csv(f)
names = list(df[mcol].unique())
data = [df.loc[df[mcol]==m, vcol].dropna() for m in names]
fig, ax = plt.subplots(figsize=(9,4))
ax.boxplot(data, labels=names, showfliers=True)
ax.axhline(1.0, ls="--", lw=0.8)
if len(sys.argv) > 6: ax.set_yscale("log")
ax.set_ylabel(vcol); ax.set_title(title)
plt.xticks(rotation=45, ha="right"); plt.tight_layout()
plt.savefig(out, dpi=300)
