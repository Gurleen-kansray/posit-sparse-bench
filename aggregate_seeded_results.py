import re, glob, os, csv
import statistics as st
from collections import defaultdict

LOG_DIR = "results/ladder_logs/seed_sweep"
OUT_DIR = "results/csv"
os.makedirs(OUT_DIR, exist_ok=True)

# widths whose columns we track for max/median pre-convergence error and gain
WIDTHS = ["p8", "p16", "p32", "p64"]

def parse_log(path):
    rows = []
    conv_iter = None
    with open(path) as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith("matrix="):
                continue
            if line.startswith("conv_iter_double="):
                conv_iter = int(line.split("=")[1])
                continue
            if line.startswith("iter "):
                header = line.split()
                continue
            if header is None:
                continue
            parts = line.split()
            if len(parts) != len(header):
                continue
            d = dict(zip(header, parts))
            rows.append(d)
    return rows, conv_iter

def to_float(s):
    try:
        v = float(s)
        if v != v:  # NaN
            return None
        return v
    except:
        return None

results = defaultdict(list)  # matrix -> list of per-seed summary dicts

for path in sorted(glob.glob(f"{LOG_DIR}/*_seed*.log")):
    fname = os.path.basename(path)
    m = re.match(r"(.+)_seed(\d+)\.log", fname)
    if not m:
        continue
    matrix, seed = m.group(1), int(m.group(2))
    rows, conv_iter = parse_log(path)
    if conv_iter is None:
        continue

    # restrict to pre-convergence iterations only (item 6)
    # NOTE: conv_iter_double == -1 means double precision never converged
    # within maxiter. Previously this made `iter < conv_iter` exclude EVERY
    # row (since no iter is < -1), silently dropping the whole seed. When
    # double doesn't converge, fall back to using all logged rows instead.
    if conv_iter == -1:
        pre_rows = rows
    else:
        pre_rows = [r for r in rows if int(r["iter"]) < conv_iter]
    if not pre_rows:
        continue

    summary = {"seed": seed, "conv_iter": conv_iter}

    for w in WIDTHS:
        q_key, n_key = f"res_{w}q", f"res_{w}n"
        pAp_q_key, pAp_n_key = f"pAp_{w}q", f"pAp_{w}n"
        guard_q_key, guard_n_key = f"guard_{w}q", f"guard_{w}n"

        q_vals = [to_float(r.get(q_key)) for r in pre_rows if to_float(r.get(q_key)) is not None]
        n_vals = [to_float(r.get(n_key)) for r in pre_rows if to_float(r.get(n_key)) is not None]

        summary[f"{w}q_max_res"] = max(q_vals) if q_vals else None
        summary[f"{w}q_med_res"] = st.median(q_vals) if q_vals else None
        summary[f"{w}n_max_res"] = max(n_vals) if n_vals else None
        summary[f"{w}n_med_res"] = st.median(n_vals) if n_vals else None

        # gain = naive residual / quire residual, per-iteration then take max/median of ratio
        gains = []
        for r in pre_rows:
            qv, nv = to_float(r.get(q_key)), to_float(r.get(n_key))
            if qv and nv and qv > 0:
                gains.append(nv / qv)
        summary[f"{w}_gain_max"] = max(gains) if gains else None
        summary[f"{w}_gain_med"] = st.median(gains) if gains else None

        # guard fire counts (item 4)
        guard_q_fires = sum(1 for r in pre_rows if r.get(guard_q_key) == "1")
        guard_n_fires = sum(1 for r in pre_rows if r.get(guard_n_key) == "1")
        summary[f"{w}q_guard_fires"] = guard_q_fires
        summary[f"{w}n_guard_fires"] = guard_n_fires

    # pAp accuracy gain vs. double-precision reference (this is the actual
    # p^T A p accuracy comparison — previously this table used res_p32q/res_p32n,
    # which is the CG residual norm, not pAp, and had no double reference at all)
    pAp_gains = []
    for r in pre_rows:
        pd_ = to_float(r.get("pAp_d"))
        pq_ = to_float(r.get("pAp_p32q"))
        pn_ = to_float(r.get("pAp_p32n"))
        if pd_ is not None and pq_ is not None and pn_ is not None and pd_ != 0:
            relq = abs(pq_ - pd_) / abs(pd_)
            reln = abs(pn_ - pd_) / abs(pd_)
            if relq > 0:
                pAp_gains.append(reln / relq)
    summary["pAp_gain_max"] = max(pAp_gains) if pAp_gains else None
    summary["pAp_gain_med"] = st.median(pAp_gains) if pAp_gains else None

    # solution error (item 1) — value at the FINAL iteration of the run,
    # not the last pre-convergence row (bug fix: use full `rows`, not `pre_rows`)
    for key in ["solerr_d", "solerr_f", "solerr_p32q", "solerr_p32n"]:
        vals = [to_float(r.get(key)) for r in rows if to_float(r.get(key)) is not None]
        summary[f"{key}_final"] = vals[-1] if vals else None

    # solution-error gain (naive/quire), same convention as pAp gain —
    # this is the direct test of James's hypothesis: does the pAp gain transfer?
    sn = summary.get("solerr_p32n_final")
    sq = summary.get("solerr_p32q_final")
    summary["solerr_gain"] = (sn / sq) if (sn is not None and sq is not None and sq > 0) else None

    results[matrix].append(summary)

# Aggregate mean+-std across seeds per matrix (item 5)
out_rows = []
for matrix in sorted(results.keys()):
    runs = results[matrix]
    n_seeds = len(runs)
    row = {"matrix": matrix, "n_seeds": n_seeds}

    for w in WIDTHS:
        for metric in [f"{w}_gain_max", f"{w}_gain_med", f"{w}q_max_res", f"{w}n_max_res"]:
            vals = [r[metric] for r in runs if r.get(metric) is not None]
            row[f"{metric}_mean"] = round(st.mean(vals), 4) if vals else "N/A"
            row[f"{metric}_std"] = round(st.pstdev(vals), 4) if len(vals) > 1 else ("0.0" if vals else "N/A")

        gf_q = sum(r.get(f"{w}q_guard_fires", 0) or 0 for r in runs)
        gf_n = sum(r.get(f"{w}n_guard_fires", 0) or 0 for r in runs)
        row[f"{w}q_total_guard_fires"] = gf_q
        row[f"{w}n_total_guard_fires"] = gf_n

    for metric in ["pAp_gain_max", "pAp_gain_med"]:
        vals = [r[metric] for r in runs if r.get(metric) is not None]
        row[f"{metric}_mean"] = round(st.mean(vals), 4) if vals else "N/A"
        row[f"{metric}_std"] = round(st.pstdev(vals), 4) if len(vals) > 1 else ("0.0" if vals else "N/A")

    # solution-error gain — reported side-by-side with p32 pAp gain so the
    # "does it transfer to the solution" comparison James asked for is explicit
    solerr_vals = [r["solerr_gain"] for r in runs if r.get("solerr_gain") is not None]
    row["solerr_gain_mean"] = round(st.mean(solerr_vals), 4) if solerr_vals else "N/A"
    row["solerr_gain_std"] = round(st.pstdev(solerr_vals), 4) if len(solerr_vals) > 1 else ("0.0" if solerr_vals else "N/A")
    row["solerr_gain_median"] = round(st.median(solerr_vals), 4) if solerr_vals else "N/A"

    out_rows.append(row)

out_path = f"{OUT_DIR}/seeded_aggregate_summary.csv"
if out_rows:
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_rows[0].keys())
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote {out_path}")
    for r in out_rows:
        print(r["matrix"],
              "pAp_gain_mean=", r.get("pAp_gain_max_mean"), "+-", r.get("pAp_gain_max_std"),
              "| solerr_gain_mean=", r.get("solerr_gain_mean"), "+-", r.get("solerr_gain_std"))
else:
    print("No results parsed — check log directory/format")
