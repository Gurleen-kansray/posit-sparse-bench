import pandas as pd
from scipy.stats import spearmanr

condest = pd.read_csv("results/csv/condest_fallback_results_full.csv")
divergence = pd.read_csv("results/csv/divergence_summary.csv")

df = condest[condest["condest_cmsw"] != "FAIL"].copy()
df["condest_cmsw"] = df["condest_cmsw"].astype(float)

method_by_matrix = df.groupby("matrix")["condest_method"].nunique()
pure_matrices = method_by_matrix[method_by_matrix == 1].index
df_pure = df[df["matrix"].isin(pure_matrices)]
invpower_matrices = df_pure[df_pure["condest_method"] == "INVPOWER"]["matrix"].unique()

merged = df.merge(divergence, on="matrix")
merged = merged[pd.to_numeric(merged["divergence_iter"], errors="coerce").notna()]
merged["divergence_iter"] = merged["divergence_iter"].astype(float)

print("Matrices classified as method-pure INVPOWER:", sorted(invpower_matrices))
print()

for bits in [8, 16, 32, 64]:
    sub = merged[(merged["bitwidth"] == bits) & (merged["matrix"].isin(invpower_matrices))]
    sub = sub[sub["matrix"] != "bcsstk14"]
    print(f"--- bits={bits}, invpower-only, bcsstk14 excluded ---")
    print(sub[["matrix", "condest_cmsw", "divergence_iter"]].sort_values("condest_cmsw").to_string(index=False))
    if len(sub) >= 3:
        rho, p = spearmanr(sub["condest_cmsw"], sub["divergence_iter"])
        print(f"n={len(sub)} rho={rho:.3f} p={p:.3f}")
    else:
        print(f"n={len(sub)} — too few points for a meaningful correlation")
    print()
