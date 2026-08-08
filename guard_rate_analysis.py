#!/usr/bin/env python3
"""
Per matrix, across all seeds in results/ladder_logs/seed_sweep/<matrix>_seed*.log:
- fraction of seeds where each guard (8q,8n,16q,16n,32q,32n,64q,64n) fires at least once
- mean number of iterations with guard=1 per seed (severity, not just onset)
- matrix size n (from header) for correlation with instability

Usage: python3 guard_rate_analysis.py <path_to_seed_sweep_dir>
"""
import sys, glob, os, re, statistics

GUARD_COLS = ["guard_p8q","guard_p8n","guard_p16q","guard_p16n",
              "guard_p32q","guard_p32n","guard_p64q","guard_p64n"]

def parse_file(path):
    with open(path) as f:
        lines = f.readlines()
    if len(lines) < 3:
        return None
    header_line = lines[0]
    m = re.search(r'n=(\d+)', header_line)
    n = int(m.group(1)) if m else None
    header = lines[1].split()
    if header[0] != 'iter':
        return None
    col_idx = {c: header.index(c) for c in GUARD_COLS if c in header}
    if not col_idx:
        return None
    guard_any = {c: 0 for c in col_idx}
    guard_count = {c: 0 for c in col_idx}
    n_iters = 0
    for line in lines[2:]:
        line = line.strip()
        if not line or not re.match(r'^-?\d+\s', line):
            continue
        vals = line.split()
        if len(vals) != len(header):
            continue
        n_iters += 1
        for c, idx in col_idx.items():
            try:
                v = int(float(vals[idx]))
            except (ValueError, IndexError):
                continue
            if v == 1:
                guard_count[c] += 1
                guard_any[c] = 1
    return n, guard_any, guard_count, n_iters

def main(seed_dir):
    matrices = sorted(set(
        os.path.basename(p).rsplit('_seed', 1)[0]
        for p in glob.glob(os.path.join(seed_dir, '*_seed*.log'))
    ))

    print(f"{'matrix':12s} {'n':8s} {'n_seeds':8s} " + " ".join(f"{c.replace('guard_p',''):>8s}" for c in GUARD_COLS))
    results = []
    for m in matrices:
        files = sorted(glob.glob(os.path.join(seed_dir, f'{m}_seed*.log')))
        n_val = None
        trip_frac = {c: [] for c in GUARD_COLS}
        for fpath in files:
            parsed = parse_file(fpath)
            if parsed is None:
                continue
            n_val, guard_any, guard_count, n_iters = parsed
            for c in GUARD_COLS:
                if c in guard_any:
                    trip_frac[c].append(guard_any[c])

        n_seeds = len(files)
        row = [m, str(n_val or '?'), str(n_seeds)]
        for c in GUARD_COLS:
            vals = trip_frac[c]
            frac = (sum(vals) / len(vals)) if vals else 0.0
            row.append(f"{frac:8.2f}")
        print(f"{row[0]:12s} {row[1]:8s} {row[2]:8s} " + " ".join(row[3:]))
        results.append((m, n_val, n_seeds, trip_frac))

    return results

if __name__ == '__main__':
    seed_dir = sys.argv[1] if len(sys.argv) > 1 else 'results/ladder_logs/seed_sweep'
    main(seed_dir)
