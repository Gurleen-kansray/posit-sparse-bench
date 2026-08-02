import re, glob

u = 3.725290298e-09
floor = 1e-15  # tighter floor

for fname in sorted(glob.glob("results/ladder_logs/*_ladder.log")):
    fails = 0
    total = 0
    filtered = 0
    max_ratio = 0
    with open(fname) as f:
        for line in f:
            m32 = re.search(r'p32q=([\d.eE+-]+)', line)
            m64 = re.search(r'p64q=([\d.eE+-]+)', line)
            if m32 and m64:
                p32q = float(m32.group(1))
                p64q = float(m64.group(1))
                if abs(p64q) < floor:
                    filtered += 1
                    continue
                rel_err = abs(p32q - p64q) / abs(p64q)
                total += 1
                ratio = rel_err / u
                max_ratio = max(max_ratio, ratio)
                if rel_err > u:
                    fails += 1
    print(f"{fname}: {fails}/{total} FAIL (filtered {filtered}), max ratio_to_u={max_ratio:.2f}")
