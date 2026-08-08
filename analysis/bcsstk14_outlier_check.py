import csv
from pathlib import Path

condest_rows = list(csv.DictReader(open("results/csv/condest_fallback_results_full.csv")))
diverge_rows = {r["matrix"]: r for r in csv.DictReader(open("results/csv/divergence_summary.csv"))}

print(f"{'matrix':<12}{'bitwidth':<9}{'method':<10}{'sat_frac':<12}{'condest':<14}{'div_iter':<10}")
for row in condest_rows:
    m = row["matrix"]
    div = diverge_rows.get(m, {}).get("divergence_iter", "n/a")
    print(f"{m:<12}{row['bitwidth']:<9}{row['condest_method']:<10}"
          f"{float(row['sat_fraction']):<12.4f}{row['condest_cmsw']:<14}{div:<10}")

print("\n=== Rows using LU fallback ===")
for r in [r for r in condest_rows if r["condest_method"] == "LU"]:
    print(f"  {r['matrix']} @ {r['bitwidth']}-bit: condest={r['condest_cmsw']}, sat_fraction={float(r['sat_fraction']):.4f}")

print("\n=== 8-bit saturation, sorted descending ===")
for r in sorted([r for r in condest_rows if r["bitwidth"]=="8"], key=lambda r: float(r["sat_fraction"]), reverse=True):
    print(f"  {r['matrix']:<12} sat_fraction={float(r['sat_fraction']):.4f}  method={r['condest_method']}")
