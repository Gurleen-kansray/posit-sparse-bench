import re, glob, os, csv

logs = glob.glob("results/*_ladder.log")
out_rows = []

for log_path in logs:
    matrix = os.path.basename(log_path).replace("_ladder.log", "")
    with open(log_path) as f:
        lines = f.readlines()

    iters = []
    for line in lines:
        m = re.search(
            r"iter=(\d+) pAp_d=([\d.eE+-]+).*?p32q=([\d.eE+-]+) p32n=([\d.eE+-]+)",
            line)
        if m:
            it, pAp_d, p32q, p32n = int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))
            if abs(pAp_d) > 1e-30:
                err_q = abs(p32q - pAp_d) / abs(pAp_d)
                err_n = abs(p32n - pAp_d) / abs(pAp_d)
            else:
                err_q = abs(p32q - pAp_d)
                err_n = abs(p32n - pAp_d)
            iters.append((it, err_q, err_n))

    if not iters:
        continue

    div_iter = None
    for i in range(len(iters)):
        it, eq, en = iters[i]
        if en > 0 and eq > 0 and en / eq > 10:
            window = iters[i:i+5]
            if all(w[2] / w[1] > 10 if w[1] > 0 else False for w in window):
                div_iter = it
                break

    final_eq = iters[-1][1]
    final_en = iters[-1][2]
    max_eq = max(e[1] for e in iters)
    max_en = max(e[2] for e in iters)

    out_rows.append({
        "matrix": matrix,
        "n_iters_logged": len(iters),
        "divergence_iter": div_iter if div_iter is not None else "no_divergence_detected",
        "final_err_quire": final_eq,
        "final_err_naive": final_en,
        "max_err_quire": max_eq,
        "max_err_naive": max_en,
        "gain_final": (final_en / final_eq) if final_eq > 0 else None,
    })

out_rows.sort(key=lambda r: r["matrix"])

with open("results/csv/divergence_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=out_rows[0].keys())
    writer.writeheader()
    writer.writerows(out_rows)

for r in out_rows:
    print(r)
