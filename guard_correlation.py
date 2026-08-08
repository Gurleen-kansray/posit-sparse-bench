#!/usr/bin/env python3
"""
Correlates posit16 guard-trip fraction (from guard_rate_analysis.py logic)
against matrix size (n) and condition number (condest_cmsw, F64raw row)
using Spearman rank correlation.
"""
import sys, glob, os, re, csv
from scipy.stats import spearmanr
import math

GUARD_COLS = ["guard_p16q","guard_p16n"]

def parse_file(path):
    with open(path) as f:
        lines = f.readlines()
    if len(lines) < 3:
        return None
    m = re.search(r'n=(\d+)', lines[0])
    n = int(m.group(1)) if m else None
    header = lines[1].split()
    if header[0] != 'iter':
        return None
    col_idx = {c: header.index(c) for c in GUARD_COLS if c in header}
    if not col_idx:
        return None
    guard_any = {c: 0 for c in col_idx}
    for line in lines[2:]:
        line = line.strip()
        if not line or not re.match(r'^-?\d+\s', line):
            continue
        vals = line.split()
        if len(vals) != len(header):
            continue
        for c, idx in col_idx.items():
            try:
                v = int(float(vals[idx]))
            except (ValueError, IndexError):
                continue
            if v == 1:
                guard_any[c] = 1
    return n, guard_any

def load_condest(csv_path):
    condest = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['bitwidth'] == 'F64raw':
                try:
                    condest[row['matrix']] = float(row['condest_cmsw'])
                except ValueError:
                    pass
    return condest

def main(seed_dir, csv_path):
    condest = load_condest(csv_path)
    matrices = sorted(set(
        os.path.basename(p).rsplit('_seed', 1)[0]
        for p in glob.glob(os.path.join(seed_dir, '*_seed*.log'))
    ))

    rows = []
    for m in matrices:
        files = sorted(glob.glob(os.path.join(seed_dir, f'{m}_seed*.log')))
        n_val = None
        trips = {c: [] for c in GUARD_COLS}
        for fpath in files:
            parsed = parse_file(fpath)
            if parsed is None:
                continue
            n_val, guard_any = parsed
            for c in GUARD_COLS:
                trips[c].append(guard_any.get(c, 0))
        if n_val is None or len(files) < 2:
            continue  # skip under-seeded matrices (s3dk*, nasasrb) - not enough data for a fraction
        frac16q = sum(trips['guard_p16q']) / len(trips['guard_p16q'])
        frac16n = sum(trips['guard_p16n']) / len(trips['guard_p16n'])
        c = condest.get(m)
        rows.append((m, n_val, c, frac16q, frac16n, len(files)))

    print(f"{'matrix':12s} {'n':>8s} {'condest':>14s} {'16q_frac':>9s} {'16n_frac':>9s} {'n_seeds':>8s}")
    for r in rows:
        cstr = f"{r[2]:.3e}" if r[2] else "N/A"
        print(f"{r[0]:12s} {r[1]:>8d} {cstr:>14s} {r[3]:>9.2f} {r[4]:>9.2f} {r[5]:>8d}")

    ns = [r[1] for r in rows]
    conds = [r[2] for r in rows if r[2] is not None]
    conds_matched = [(r[1], r[2], r[3], r[4]) for r in rows if r[2] is not None]

    frac16q_all = [r[3] for r in rows]
    frac16n_all = [r[4] for r in rows]

    print("\n--- Spearman: guard-trip fraction vs n (all matrices with >=2 seeds) ---")
    rho_n_q, p_n_q = spearmanr(ns, frac16q_all)
    rho_n_n, p_n_n = spearmanr(ns, frac16n_all)
    print(f"16q vs n: rho={rho_n_q:.4f} p={p_n_q:.4f}")
    print(f"16n vs n: rho={rho_n_n:.4f} p={p_n_n:.4f}")

    if len(conds_matched) >= 3:
        ns_m = [x[0] for x in conds_matched]
        conds_m = [x[1] for x in conds_matched]
        f16q_m = [x[2] for x in conds_matched]
        f16n_m = [x[3] for x in conds_matched]
        print(f"\n--- Spearman: guard-trip fraction vs condest_cmsw (n={len(conds_matched)} matrices) ---")
        rho_c_q, p_c_q = spearmanr(conds_m, f16q_m)
        rho_c_n, p_c_n = spearmanr(conds_m, f16n_m)
        print(f"16q vs condest: rho={rho_c_q:.4f} p={p_c_q:.4f}")
        print(f"16n vs condest: rho={rho_c_n:.4f} p={p_c_n:.4f}")
    else:
        print("\nNot enough matrices with condest data for correlation.")

if __name__ == '__main__':
    seed_dir = sys.argv[1] if len(sys.argv) > 1 else 'results/ladder_logs/seed_sweep'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'results/csv/static_conditioning.csv'
    main(seed_dir, csv_path)
