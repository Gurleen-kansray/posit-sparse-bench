import pandas as pd
from scipy.stats import spearmanr

condest = pd.read_csv("results/csv/condest_fallback_results_full.csv")
divergence = pd.read_csv("results/csv/divergence_summary.csv")

df = condest[condest["condest_cmsw"] != "FAIL"].copy()
df["condest_cmsw"] = df["condest_cmsw"].astype(float)

# looser filter: just require condest_method == INVPOWER at bitwidth==16
# (don't require purity across all bitwidths)
invpower16_matrices = df[(df["bitwidth"] == 16) & (df["condest_method"] == "INVPOWER")]["matrix"].unique()

merged = df.merge(divergence, on="matrix")
merged = merged[pd.to_numeric(merged["divergence_iter"], errors="coerce").notna()]
merged["divergence_iter"] = merged["divergence_iter"].astype(float)

print("Matrices with condest_method==INVPOWER at bits=16:", sorted(invpower16_matrices))
print()

sub = merged[(merged["bitwidth"] == 16) & (merged["matrix"].isin(invpower16_matrices))]
sub_excl = sub[sub["matrix"] != "bcsstk14"]

print("--- bits=16, INVPOWER-at-16 only, bcsstk14 INCLUDED ---")
print(sub[["matrix", "condest_cmsw", "divergence_iter"]].sort_values("condest_cmsw").to_string(index=False))
if len(sub) >= 3:
    rho, p = spearmanr(sub["condest_cmsw"], sub["divergence_iter"])
    print(f"n={len(sub)} rho={rho:.3f} p={p:.3f}")
print()

print("--- bits=16, INVPOWER-at-16 only, bcsstk14 EXCLUDED ---")
print(sub_excl[["matrix", "condest_cmsw", "divergence_iter"]].sort_values("condest_cmsw").to_string(index=False))
if len(sub_excl) >= 3:
    rho, p = spearmanr(sub_excl["condest_cmsw"], sub_excl["divergence_iter"])
    print(f"n={len(sub_excl)} rho={rho:.3f} p={p:.3f}")
else:
    print(f"n={len(sub_excl)} — too few points for a meaningful correlation")
