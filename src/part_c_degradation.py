import pandas as pd
from scipy.stats import spearmanr

condest = pd.read_csv("results/csv/condest_fallback_results_full.csv")
divergence = pd.read_csv("results/csv/divergence_summary.csv")

# pivot to get condest at each bitwidth per matrix
pivot = condest[condest["condest_cmsw"] != "FAIL"].copy()
pivot["condest_cmsw"] = pivot["condest_cmsw"].astype(float)
wide = pivot.pivot(index="matrix", columns="bitwidth", values="condest_cmsw")

# degradation ratio: how much worse is condest at low precision vs 64-bit
wide["ratio_8_64"] = wide[8] / wide[64]
wide["ratio_16_64"] = wide[16] / wide[64]

merged = wide.merge(divergence, on="matrix")
merged = merged[pd.to_numeric(merged["divergence_iter"], errors="coerce").notna()]
merged["divergence_iter"] = merged["divergence_iter"].astype(float)

print(merged[["matrix", "ratio_8_64", "ratio_16_64", "divergence_iter"]].sort_values("ratio_8_64").to_string(index=False))
print()

for col in ["ratio_8_64", "ratio_16_64"]:
    sub = merged.dropna(subset=[col])
    rho, p = spearmanr(sub[col], sub["divergence_iter"])
    print(f"{col}: n={len(sub)} rho={rho:.3f} p={p:.3f}")
