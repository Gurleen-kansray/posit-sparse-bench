import sys, glob, os, re

def parse_log(path):
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
        if v != v:
            return None
        return v
    except (ValueError, TypeError):
        return None

def main(seed_dir):
    matrices = sorted(set(
        os.path.basename(p).rsplit('_seed', 1)[0]
        for p in glob.glob(os.path.join(seed_dir, '*_seed*.log'))
    ))
    results = []
    for m in matrices:
        files = sorted(glob.glob(os.path.join(seed_dir, f'{m}_seed*.log')))
        for fpath in files:
            row = parse_log(fpath)
            if row is None:
                continue
            pd = safe_float(row.get('pAp_d'))
            pq = safe_float(row.get('pAp_p32q'))
            pn = safe_float(row.get('pAp_p32n'))
            if pd and pq and pn:
                relq = abs(pq - pd) / abs(pd) if pd != 0 else None
                reln = abs(pn - pd) / abs(pd) if pd != 0 else None
                if relq is not None and reln is not None and relq > 0:
                    gain = reln / relq
                    results.append((gain, m, os.path.basename(fpath)))
    results.sort(reverse=True)
    print("Top 15 pAp gains (matrix, file, gain):")
    for gain, m, f in results[:15]:
        print(f"{gain:16.4f}  {m:12s}  {f}")
    print(f"\nTrue max pAp gain across ALL matrices/seeds: {results[0][0]:.4f}  (matrix={results[0][1]}, file={results[0][2]})")
    print(f"Total valid seed-level gain values: {len(results)}")

if __name__ == '__main__':
    seed_dir = sys.argv[1] if len(sys.argv) > 1 else 'results/ladder_logs/seed_sweep'
    main(seed_dir)
