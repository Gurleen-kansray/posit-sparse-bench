#!/usr/bin/env python3
"""
Verify cond(x,y) worst-case Wilkinson conditioning results across all matrices.
Log format: line1 = 'matrix=... n=... seed=...', line2 = header, then space-delimited rows.
"""
import os, re, sys, glob, csv
import pandas as pd

LOG_DIR = "results/cond_probe"
OUT_CSV = "results/csv/cond_probe_summary.csv"
STATIC_COND_CSV = "results/csv/static_conditioning.csv"  # for matrix cond(A), if present
TOL = 1e-10

# arithmetic -> (residual column, solerr column)
ARMS = {
    "d":     ("res_d", "solerr_d"),
    "f":     ("res_f", "solerr_f"),
    "ffma":  ("res_ffma", "solerr_ffma"),
    "p32q":  ("res_p32q", "solerr_p32q"),
    "p32n":  ("res_p32n", "solerr_p32n"),
}

def load_static_cond():
    if not os.path.exists(STATIC_COND_CSV):
        return {}
    df = pd.read_csv(STATIC_COND_CSV)
    # try common column names
    name_col = next((c for c in df.columns if "matrix" in c.lower() or "name" in c.lower()), None)
    cond_col = next((c for c in df.columns if "cond" in c.lower()), None)
    if not name_col or not cond_col:
        return {}
    return dict(zip(df[name_col], df[cond_col]))

def parse_log(path):
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()
    if len(lines) < 3:
        return None, None
    meta = lines[0].strip()
    header = lines[1].split()
    data_lines = [l for l in lines[2:] if l.strip()]
    rows = [l.split() for l in data_lines]
    df = pd.DataFrame(rows, columns=header).apply(pd.to_numeric, errors="coerce")
    seed = int(re.search(r"seed=(\d+)", meta).group(1))
    return df, seed

def conv_iter(df, res_col, tol=TOL):
    ok = df[df[res_col] < tol]
    return int(ok.iloc[0]["iter"]) if not ok.empty else -1

def main():
    static_cond = load_static_cond()
    matrices = sorted(set(
        os.path.basename(f).split("_seed")[0]
        for f in glob.glob(f"{LOG_DIR}/*_seed*.log")
    ))
    if not matrices:
        print(f"No logs found in {LOG_DIR}/"); sys.exit(1)

    print(f"Found {len(matrices)} matrices: {matrices}\n")
    rows_out = []

    for mat in matrices:
        logs = sorted(glob.glob(f"{LOG_DIR}/{mat}_seed*.log"))
        seed_nums = sorted(set(int(re.search(r"_seed(\d+)\.log$", f).group(1)) for f in logs))
        missing = [s for s in range(1, 51) if s not in seed_nums]

        print(f"--- {mat} ---")
        print(f"  Logs: {len(logs)} | Seeds present: {len(seed_nums)}/50")
        if missing:
            print(f"  MISSING seeds: {missing}")
        else:
            print(f"  All 50 seeds present.")

        best_row, best_val, best_seed, best_df = None, -1, None, None
        for f in logs:
            df, seed = parse_log(f)
            if df is None or "cond_xy_d" not in df.columns:
                print(f"    (skipping empty/corrupt log: {f})")
                continue
            idx = df["cond_xy_d"].idxmax()
            val = df.loc[idx, "cond_xy_d"]
            if val > best_val:
                best_val, best_row, best_seed, best_df = val, df.loc[idx], seed, df

        if best_row is None:
            print(f"  WARNING: no cond_xy_d parsed for {mat}\n")
            continue

        worst_iter = int(best_row["iter"])
        print(f"  Worst cond(x,y): {best_val:.4g} (seed {best_seed}, iter {worst_iter})")
        cond_A = static_cond.get(mat, None)
        print(f"  cond(A): {cond_A if cond_A else 'not found in ' + STATIC_COND_CSV}")

        out = {
            "matrix": mat,
            "cond_A": cond_A,
            "worst_cond_xy": best_val,
            "worst_seed": best_seed,
            "worst_iter": worst_iter,
        }
        for arm, (res_col, solerr_col) in ARMS.items():
            ci = conv_iter(best_df, res_col)
            se = best_row.get(solerr_col, None)
            print(f"    {arm:6s} conv_iter(tol=1e-10)={ci:6d}  solerr@worst_iter={se}")
            out[f"conv_iter_{arm}"] = ci
            out[f"solerr_{arm}"] = se
        rows_out.append(out)
        print()

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    if rows_out:
        fieldnames = sorted(set(k for r in rows_out for k in r.keys()))
        with open(OUT_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows_out)
        print(f"Summary written to {OUT_CSV}")
    else:
        print("Nothing parsed.")

if __name__ == "__main__":
    main()
