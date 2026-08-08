import pandas as pd
from scipy.stats import spearmanr

condest = pd.read_csv("results/csv/condest_fallback_results_full.csv")
divergence = pd.read_csv("results/csv/divergence_summary.csv")

merged = condest.merge(divergence, on="matrix")
merged = merged[merged["condest_cmsw"] != "FAIL"]
merged["condest_cmsw"] = merged["condest_cmsw"].astype(float)

merged = merged[pd.to_numeric(merged["divergence_iter"], errors="coerce").notna()]
merged["divergence_iter"] = merged["divergence_iter"].astype(float)

print(f"usable rows per bitwidth:")
print(merged.groupby("bitwidth").size())
print()

for bits in [8, 16, 32, 64]:
    sub = merged[merged["bitwidth"] == bits]
    if len(sub) < 3:
        print(f"bits={bits}: insufficient data (n={len(sub)})")
        continue
    rho_c, p_c = spearmanr(sub["condest_cmsw"], sub["divergence_iter"])
    rho_l, p_l = spearmanr(sub["lambda_max"], sub["divergence_iter"])
    print(f"bits={bits}: n={len(sub)} | condest rho={rho_c:.3f} p={p_c:.3f} | lambda_max rho={rho_l:.3f} p={p_l:.3f}")

print()
print("=== Raw data per bitwidth (sanity check) ===")
for bits in [8, 16]:
    sub = merged[merged["bitwidth"] == bits][["matrix", "condest_cmsw", "lambda_max", "divergence_iter"]]
    sub = sub.sort_values("condest_cmsw")
    print(f"\n--- bits={bits} ---")
    print(sub.to_string(index=False))

print()
print("=== Group comparison: divergence vs no-divergence ===")
for bits in [8, 16, 32, 64]:
    sub = merged[merged["bitwidth"] == bits]
    zero = sub[sub["divergence_iter"] == 0]["condest_cmsw"]
    nonzero = sub[sub["divergence_iter"] > 0]["condest_cmsw"]
    print(f"bits={bits}: divergence=0 group (n={len(zero)}) condest range: "
          f"{zero.min():.2e} - {zero.max():.2e}" if len(zero) else f"bits={bits}: no zero-group")
    print(f"          divergence>0 group (n={len(nonzero)}) condest range: "
          f"{nonzero.min():.2e} - {nonzero.max():.2e}" if len(nonzero) else "")
