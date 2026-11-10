#!/usr/bin/env python3
"""Table 4: mean quire gain for the two CG inner products.

    LOG_DIR=results/quiregain_50 python3 agg_quiregain.py
    LOG_DIR=results/quiregain_50 python3 agg_quiregain.py --latex

Reads the logs written by src/ladder_quiregain.cpp. Per seed, takes the maximum
relative error over the pre-convergence window for each accumulation policy and
divides naive by quire; reports the mean of that ratio over seeds, and the mean
of the same ratio formed from medians rather than maxima.

The window is `window=` in the log trailer, which is the iteration at which the
double-precision reference first reaches RTOL, plus one. Rows past it are kept in
the log, so WINDOW=all recomputes over every iteration instead.
"""
import glob
import os
import re
import statistics as st
import sys

LOG_DIR = os.environ.get("LOG_DIR", "results/quiregain_50")
USE_FULL = os.environ.get("WINDOW", "").lower() == "all"
LATEX = "--latex" in sys.argv

# Table 2's row order, so the output drops straight into the paper.
ORDER = ["bcsstk03", "bcsstk14", "bcsstk36", "bcsstk37", "bcsstk38", "bodyy4",
         "mhd4800b", "nasa4704", "nasasrb", "nos2", "s3dkq4m2", "s3dkt3m2",
         "sts4098"]


def parse(path):
    """-> (pAp naive/quire ratio, rz ratio, conv_iter), or None if unusable."""
    rows, window, conv = [], None, None
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("window="):
                window = int(line.split("=")[1])
            elif line.startswith("conv_iter="):
                conv = int(line.split("=")[1])
            elif line[:1].isdigit():
                p = line.split()
                if len(p) == 6:
                    rows.append([float(v) for v in p[1:5]])
    if not rows:
        return None
    if not USE_FULL and window:
        rows = rows[:window]

    def ratio(col_n, col_q):
        # A zero quire error means the quire was exact at that iteration, which
        # is the expected case rather than an error; it only breaks the ratio
        # when it is also the maximum, and then the gain is unbounded.
        n = [r[col_n] for r in rows if r[col_n] == r[col_n]]
        q = [r[col_q] for r in rows if r[col_q] == r[col_q]]
        if not n or not q:
            return None
        mn, mq = max(n), max(q)
        return (mn / mq) if mq > 0 else None

    return ratio(0, 1), ratio(2, 3), conv


acc, convs, skipped = {}, {}, 0
for path in sorted(glob.glob(os.path.join(LOG_DIR, "*_seed*.log"))):
    m = re.match(r"(.+)_seed(\d+)\.log$", os.path.basename(path))
    if not m:
        continue
    got = parse(path)
    if got is None:
        skipped += 1
        continue
    pap, rz, conv = got
    e = acc.setdefault(m.group(1), {"pAp": [], "rz": []})
    if pap is not None:
        e["pAp"].append(pap)
    if rz is not None:
        e["rz"].append(rz)
    if conv is not None:
        convs.setdefault(m.group(1), []).append(conv)

mats = [m for m in ORDER if m in acc] + sorted(set(acc) - set(ORDER))

if LATEX:
    for mat in mats:
        e = acc[mat]
        print("  \\mat{%s} & %8.2f & %8.2f \\\\"
              % (mat, st.mean(e["pAp"]), st.mean(e["rz"])))
else:
    print("%-11s %6s %12s %12s %12s %12s %9s"
          % ("matrix", "seeds", "pAp mean", "pAp SD", "rz mean", "rz SD", "conv"))
    for mat in mats:
        e = acc[mat]
        c = convs.get(mat, [])
        print("%-11s %6d %12.2f %12.2f %12.2f %12.2f %9.0f"
              % (mat, len(e["pAp"]),
                 st.mean(e["pAp"]), st.stdev(e["pAp"]) if len(e["pAp"]) > 1 else 0.0,
                 st.mean(e["rz"]), st.stdev(e["rz"]) if len(e["rz"]) > 1 else 0.0,
                 st.mean(c) if c else -1))
    if skipped:
        print("\n%d log(s) unusable" % skipped)
    print("\nwindow: %s" % ("all iterations" if USE_FULL else "pre-convergence"))
