import glob, os, re, sys
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def parse(path):
    L = open(path).read().splitlines()
    if len(L) < 3: return None
    h = L[1].split()
    if not h or h[0] != "iter": return None
    rows = []
    for ln in L[2:]:
        v = ln.split()
        if len(v) == len(h) and re.match(r"^-?\d+$", v[0]):
            rows.append(v)
    if not rows: return None
    df = pd.DataFrame(rows, columns=h).apply(pd.to_numeric, errors="coerce")
    return df

def ratios(dirs, at=("last",)):
    out = []
    for d in dirs:
        for f in glob.glob(os.path.join(d, "*_seed*.log")):
            m, s = os.path.basename(f)[:-4].rsplit("_seed", 1)
            df = parse(f)
            if df is None or "solerr_p32q" not in df: continue
            for a in at:
                r = df.iloc[-1] if a == "last" else (df[df.iter == a].iloc[0] if (df.iter == a).any() else None)
                if r is None: continue
                if r.solerr_p32q > 0 and r.solerr_p32n > 0:
                    out.append((m, int(s), a, r.solerr_p32n / r.solerr_p32q))
    return pd.DataFrame(out, columns=["matrix", "seed", "iter", "ratio"])

def box(df, col, out, title, ylabel, log=True, order=None):
    names = order or list(df.matrix.unique())
    data = [df.loc[df.matrix == m, col].dropna() for m in names]
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.boxplot(data, tick_labels=names, showfliers=True)
    ax.axhline(1.0, ls="--", lw=0.8)
    if log: ax.set_yscale("log")
    ax.set_ylabel(ylabel); ax.set_title(title)
    plt.xticks(rotation=45, ha="right"); plt.tight_layout()
    plt.savefig(out, dpi=300); plt.close()

# Fig A: 13 main matrices, naive/quire solution-error ratio at last row
main = ratios(["results/ladder_logs/seed_sweep"])
main.to_csv("revision/tables/main13_ratio_perseed.csv", index=False)
box(main, "ratio", "revision/figs/main13_ratio.png",
    "Solution error ratio naive/quire (>1: quire better)", "naive / quire")

# Fig B: 7 extra matrices
nm = ratios(["results/newmat_q512"])
nm.to_csv("revision/tables/newmat_ratio_perseed.csv", index=False)
box(nm, "ratio", "revision/figs/newmat_ratio.png",
    "Extra matrices: naive/quire error ratio", "naive / quire")

# Fig C: refinement, no-quire vs quire
r = pd.read_csv("revision/tables/refinement_baseline_per_seed.csv")
ms = list(r.matrix.unique())
fig, ax = plt.subplots(figsize=(9, 4))
pos = np.arange(len(ms))
ax.boxplot([r.loc[r.matrix == m, "noquire_improv_posit"] for m in ms], positions=pos - 0.18, widths=0.3,
           patch_artist=True, boxprops=dict(facecolor="lightgray"), tick_labels=[""] * len(ms))
ax.boxplot([r.loc[r.matrix == m, "q512_improv_posit"] for m in ms], positions=pos + 0.18, widths=0.3,
           patch_artist=True, boxprops=dict(facecolor="lightblue"), tick_labels=[""] * len(ms))
ax.set_xticks(pos); ax.set_xticklabels(ms, rotation=45, ha="right")
ax.set_yscale("log"); ax.axhline(1, ls="--", lw=0.8)
ax.set_ylabel("improvement factor (before/after)")
ax.set_title("Iterative refinement: no-quire (gray) vs quire (blue)")
plt.tight_layout(); plt.savefig("revision/figs/refinement_box.png", dpi=300); plt.close()

# Table: sensitivity of naive/quire ratio to CG iteration (Review 3, Table 3)
its = [1, 10, 100, 500, 1000, "last"]
sens = ratios(["results/ladder_logs/seed_sweep"], at=its)
tab = sens.groupby(["matrix", "iter"]).ratio.median().unstack("iter")
tab = tab[[str(i) if str(i) in tab.columns else i for i in its if (str(i) in tab.columns or i in tab.columns)]]
tab.round(4).to_csv("revision/tables/iter_sensitivity_median_ratio.csv")
print(tab.round(4).to_string())

# Fig D: exact-b robustness (quire arm only)
e = pd.read_csv("revision/tables/exactb_perseed.csv")
box(e, "improvement_posit", "revision/figs/exactb_box.png",
    "Refinement improvement with exact-b (quire)", "solerr before / after")
