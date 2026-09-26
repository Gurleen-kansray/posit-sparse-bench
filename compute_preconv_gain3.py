#!/usr/bin/env python3
import glob, os, statistics

LOG_DIR = "results/ladder_logs/seed_sweep"
MATRICES = ["bcsstk03","bcsstk14","bcsstk36","bcsstk37","bcsstk38","bodyy4",
            "mhd4800b","nasa4704","nasasrb","nos2","s3dkq4m2","s3dkt3m2","sts4098"]
TOL = 1e-6

def parse_log(path):
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()
    header = None
    rows = []
    for line in lines:
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "iter":
            header = parts
            continue
        if header is None or len(parts) != len(header):
            continue
        rows.append(dict(zip(header, parts)))
    return rows

def preconv_gain_for_seed(rows):
    if not rows:
        return None
    res0 = float(rows[0]["res_d"])
    if res0 == 0:
        return None
    conv_iter = None
    for r in rows:
        rel = float(r["res_d"]) / res0
        if rel < TOL:
            conv_iter = int(r["iter"])
            break
    if conv_iter is None:
        return None

    gains = []
    for r in rows:
        it = int(r["iter"])
        if it >= conv_iter:
            continue
        try:
            pap_d = float(r["pAp_d"]); pap_q = float(r["pAp_p32q"]); pap_n = float(r["pAp_p32n"])
        except (KeyError, ValueError):
            continue
        if any(v != v for v in (pap_d, pap_q, pap_n)):
            continue
        if pap_d == 0:
            continue
        err_q = abs(pap_q - pap_d) / abs(pap_d)
        err_n = abs(pap_n - pap_d) / abs(pap_d)
        if err_q == 0:
            continue
        gains.append(err_n / err_q)
    if not gains:
        return None
    return statistics.mean(gains), conv_iter

def main():
    print(f"{'matrix':<12}{'n_total':>9}{'n_used':>8}{'preconv_gain_mean':>20}{'preconv_gain_std':>18}")
    for matrix in MATRICES:
        files = sorted(glob.glob(os.path.join(LOG_DIR, f"{matrix}_seed*.log")))
        seed_gains = []
        for fpath in files:
            rows = parse_log(fpath)
            result = preconv_gain_for_seed(rows)
            if result is not None:
                g, _ = result
                seed_gains.append(g)
        if not seed_gains:
            print(f"{matrix:<12}{len(files):>9}{0:>8}{'NO DATA':>20}{'':>18}")
            continue
        mean_g = statistics.mean(seed_gains)
        std_g = statistics.stdev(seed_gains) if len(seed_gains) > 1 else 0.0
        print(f"{matrix:<12}{len(files):>9}{len(seed_gains):>8}{mean_g:>20.4f}{std_g:>18.4f}")

if __name__ == "__main__":
    main()
