#!/usr/bin/env python3
"""
Computes, per matrix, from results/ladder_logs/seed_sweep/<matrix>_seed*.log:
- solution-error gain ratio (solerr_p32n / solerr_p32q) at last logged iteration per seed
- pAp accuracy gain ratio (relerr_p32n / relerr_p32q vs pAp_d) at last logged iteration per seed
- mean/std/median across seeds, plus n_seeds actually used (non-empty, parseable)

Usage: python3 compute_gain_ratios.py <path_to_seed_sweep_dir>
"""
import sys, glob, os, statistics, re

def parse_log(path):
    """Return list of dict rows (last row only needed, but parse all for safety)."""
    with open(path) as f:
        lines = f.readlines()
    if len(lines) < 3:
        return None
    header = lines[1].split()
    if header[0] != 'iter':
        return None
    for line in reversed(lines[2:]):
        line = line.strip()
        if not line or not re.match(r'^-?\d+\s', line):
            continue
        vals = line.split()
        if len(vals) != len(header):
            continue
        return dict(zip(header, vals))
    return None

def safe_float(x):
    try:
        v = float(x)
        if v != v:  # NaN check
            return None
        return v
    except (ValueError, TypeError):
        return None

def main(seed_dir):
    matrices = sorted(set(
        os.path.basename(p).rsplit('_seed', 1)[0]
        for p in glob.glob(os.path.join(seed_dir, '*_seed*.log'))
    ))

    print(f"{'matrix':12s} {'n_seeds':8s} {'n_valid':8s} {'solerr_ratio_mean':18s} {'solerr_ratio_std':17s} {'solerr_ratio_median':19s} {'pAp_gain_mean':14s} {'pAp_gain_std':13s}")
    for m in matrices:
        files = sorted(glob.glob(os.path.join(seed_dir, f'{m}_seed*.log')))
        solerr_ratios = []
        pap_gains = []
        for fpath in files:
            row = parse_log(fpath)
            if row is None:
                continue
            solq = safe_float(row.get('solerr_p32q'))
            soln = safe_float(row.get('solerr_p32n'))
            if solq and soln and solq > 0:
                solerr_ratios.append(soln / solq)

            pd = safe_float(row.get('pAp_d'))
            pq = safe_float(row.get('pAp_p32q'))
            pn = safe_float(row.get('pAp_p32n'))
            if pd and pq and pn:
                relq = abs(pq - pd) / abs(pd) if pd != 0 else None
                reln = abs(pn - pd) / abs(pd) if pd != 0 else None
                if relq is not None and reln is not None and relq > 0:
                    pap_gains.append(reln / relq)

        n_valid = len(solerr_ratios)
        if n_valid == 0:
            print(f"{m:12s} {len(files):<8d} {0:<8d} {'--':18s} {'--':17s} {'--':19s} {'--':14s} {'--':13s}")
            continue

        mean_s = statistics.mean(solerr_ratios)
        std_s = statistics.stdev(solerr_ratios) if n_valid > 1 else 0.0
        med_s = statistics.median(solerr_ratios)
        mean_p = statistics.mean(pap_gains) if pap_gains else float('nan')
        std_p = statistics.stdev(pap_gains) if len(pap_gains) > 1 else 0.0

        print(f"{m:12s} {len(files):<8d} {n_valid:<8d} {mean_s:<18.4f} {std_s:<17.4f} {med_s:<19.4f} {mean_p:<14.4f} {std_p:<13.4f}")

if __name__ == '__main__':
    seed_dir = sys.argv[1] if len(sys.argv) > 1 else 'results/ladder_logs/seed_sweep'
    main(seed_dir)
