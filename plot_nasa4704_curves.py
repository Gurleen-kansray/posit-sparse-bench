#!/usr/bin/env python3
"""Fig. 2: solution error against CG iteration for nasa4704, run to n iterations.

    LOG_DIR=results/maxn_50 python3 plot_nasa4704_curves.py

Reads <matrix>_seed*_it<n>.log written by cg_compare_seeded with maxiter = n.
Plots, per arithmetic, the median over seeds with an interquartile band and a
min-to-max band, and prints the figures the surrounding text quotes.
"""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOG_DIR = os.environ.get("LOG_DIR", "results/maxn_50")
MATRIX = os.environ.get("MATRIX", "nasa4704")
OUT = os.environ.get("OUT_PNG", "fig_nasa4704_curves.png")

matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.size"] = 11
matplotlib.rcParams["axes.edgecolor"] = "#333333"
matplotlib.rcParams["axes.linewidth"] = 0.8

ARMS = [("solerr_d", "double", "#1b1b1b"),
        ("solerr_p32q", "posit32 + quire", "#3d5a80"),
        ("solerr_p32n", "posit32 naive", "#bb6b2a")]

files = sorted(glob.glob(os.path.join(LOG_DIR, "%s_seed*.log" % MATRIX)))
if not files:
    sys.exit("no logs matching %s/%s_seed*.log" % (LOG_DIR, MATRIX))

cols, series, nit = None, {k: [] for k, _, _ in ARMS}, None
for f in files:
    with open(f, errors="replace") as fh:
        lines = fh.read().splitlines()
    hi = next(i for i, l in enumerate(lines) if l.startswith("iter "))
    header = lines[hi].split()
    if cols is None:
        cols = {k: header.index(k) for k, _, _ in ARMS}
    rows = [l.split() for l in lines[hi + 1:] if l[:1].isdigit()]
    rows = [r for r in rows if len(r) == len(header)]
    if nit is None:
        nit = len(rows)
    if len(rows) != nit:          # a truncated seed would misalign the bands
        print("skipping %s: %d rows, expected %d" % (f, len(rows), nit))
        continue
    arr = np.array([[float(r[cols[k]]) for k, _, _ in ARMS] for r in rows])
    for j, (k, _, _) in enumerate(ARMS):
        series[k].append(arr[:, j])

nseeds = len(series[ARMS[0][0]])
it = np.arange(nit)

fig, ax = plt.subplots(figsize=(9, 5.2))
for k, label, colour in ARMS:
    M = np.vstack(series[k])
    med = np.median(M, axis=0)
    lo, hi = np.percentile(M, 25, axis=0), np.percentile(M, 75, axis=0)
    ax.fill_between(it, M.min(axis=0), M.max(axis=0), color=colour, alpha=0.10,
                    linewidth=0)
    ax.fill_between(it, lo, hi, color=colour, alpha=0.28, linewidth=0)
    ax.plot(it, med, color=colour, linewidth=1.4, label=label, zorder=3)

ax.set_yscale("log")
ax.set_xlim(0, nit - 1)
ax.set_xlabel("CG iteration")
ax.set_ylabel(r"relative solution error $\|x_k-x^\star\|_2/\|x^\star\|_2$")
ax.set_title("%s, %d iterations, %d seeds" % (MATRIX, nit, nseeds), fontsize=12, pad=10)
ax.grid(axis="both", which="major", linestyle="-", linewidth=0.4, color="#cccccc")
ax.grid(axis="y", which="minor", linestyle=":", linewidth=0.3, color="#e8e8e8")
ax.set_axisbelow(True)
ax.legend(frameon=False, loc="lower left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches="tight")
print("wrote %s (%d seeds, %d iterations)" % (OUT, nseeds, nit))

# The numbers the caption and the surrounding paragraph quote.
print("\nmedian relative solution error")
for mark in (2000, nit - 1):
    vals = ", ".join("%s %.1e" % (label, np.median(np.vstack(series[k])[:, mark]))
                     for k, label, _ in ARMS)
    print("  iteration %-5d %s" % (mark, vals))
Q = np.vstack(series["solerr_p32q"])[:, -1]
N = np.vstack(series["solerr_p32n"])[:, -1]
print("\nat the last iteration: quire smaller on %d of %d seeds, naive on %d"
      % (int((Q < N).sum()), nseeds, int((N < Q).sum())))
