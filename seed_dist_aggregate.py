import re, glob, os, csv
import statistics as st

results = {}

for path in sorted(glob.glob("results/ladder_logs/seed_dist/*_seed*.log")):
    fname = os.path.basename(path)
    m = re.match(r"(.+)_seed(\d+)\.log", fname)
    if not m:
        continue
    matrix, seed = m.group(1), int(m.group(2))

    with open(path) as f:
        lines = f.readlines()

    iters = []
    conv_iter = None

    for line in lines:
        if line.startswith("iter="):
            d = {}
            for tok in line.split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    d[k] = v
            if "iter" not in d or "p32q" not in d or "p32n" not in d:
                continue
            it = int(d["iter"])
            pAp_d = float(d.get("pAp_d", d.get("p64q", "nan")))
            p32q_v = float(d["p32q"])
            p32n_v = float(d["p32n"])
            if abs(pAp_d) > 1e-30:
                err_q = abs(p32q_v - pAp_d) / abs(pAp_d)
                err_n = abs(p32n_v - pAp_d) / abs(pAp_d)
            else:
                err_q = abs(p32q_v - pAp_d)
                err_n = abs(p32n_v - pAp_d)
            iters.append((it, err_q, err_n))
        elif line.startswith("  residual="):
            resid = float(line.split("=")[1])
            if resid < 1e-10 and conv_iter is None:
                conv_iter = iters[-1][0] if iters else None

    div_iter = None
    for i in range(len(iters)):
        it, eq, en = iters[i]
        if en > 0 and eq > 0 and en / eq > 10:
            window = iters[i:i+5]
            if len(window) == 5 and all(
                (w[2] / w[1] > 10 if w[1] > 0 else False) for w in window
            ):
                div_iter = it
                break

    results.setdefault(matrix, []).append({
        "seed": seed,
        "div_iter": div_iter,
        "conv_iter": conv_iter,
        "n_iters": len(iters),
    })

out_rows = []
for matrix in sorted(results.keys()):
    runs = results[matrix]
    n_seeds = len(runs)
    div_vals = [r["div_iter"] for r in runs if r["div_iter"] is not None]
    conv_vals = [r["conv_iter"] for r in runs if r["conv_iter"] is not None]
    seed1 = next((r["div_iter"] for r in runs if r["seed"] == 1), None)

    row = {
        "matrix": matrix,
        "n_seeds": n_seeds,
        "n_div_detected": len(div_vals),
        "div_mean": round(st.mean(div_vals), 1) if div_vals else "N/A",
        "div_std": round(st.pstdev(div_vals), 1) if len(div_vals) > 1 else ("0.0" if div_vals else "N/A"),
        "div_min": min(div_vals) if div_vals else "N/A",
        "div_max": max(div_vals) if div_vals else "N/A",
        "seed1_div": seed1 if seed1 is not None else "None",
        "pct_div_detected": round(100 * len(div_vals) / n_seeds, 1) if n_seeds else "N/A",
        "conv_mean": round(st.mean(conv_vals), 1) if conv_vals else "N/A",
        "conv_std": round(st.pstdev(conv_vals), 1) if len(conv_vals) > 1 else ("0.0" if conv_vals else "N/A"),
        "pct_converged": round(100 * len(conv_vals) / n_seeds, 1) if n_seeds else "N/A",
    }
    out_rows.append(row)

os.makedirs("results/csv", exist_ok=True)
out_path = "results/csv/seed_distribution_summary.csv"
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=out_rows[0].keys())
    writer.writeheader()
    writer.writerows(out_rows)

print(f"Wrote {out_path}")
print()
header = f"{'matrix':<12}{'n_seeds':<9}{'div_mean':<10}{'div_std':<9}{'div_min':<9}{'div_max':<9}{'pct_div':<9}{'conv_mean':<11}{'pct_conv':<9}"
print(header)
for r in out_rows:
    line = f"{r['matrix']:<12}{r['n_seeds']:<9}{str(r['div_mean']):<10}{str(r['div_std']):<9}{str(r['div_min']):<9}{str(r['div_max']):<9}{str(r['pct_div_detected']):<9}{str(r['conv_mean']):<11}{str(r['pct_converged']):<9}"
    print(line)
