#!/usr/bin/env python3
"""Generate SLURM job lists and a submission manifest for the full experiment set.

Run from the repo root:   python3 hpc/make_joblists.py

Writes hpc/jobs/<program>_<class>.txt and hpc/jobs/manifest.tsv. submit.sh
reads the manifest; nothing else needs editing when a program or matrix changes.

COST MODEL. Per-seed wall time on one core, fitted to measurements on an
i7-12800H and scaled by CLUSTER for a fully packed EPYC Milan node at 2.3 GHz:

  cg_compare_seeded   t = K * (2.245e-6*nnz + 9.279e-5*n)      [5 arms, breaks at
                                                                conv_iter+20]
  generic_ladder      t = 300 * (2.395e-5*n + 9.10e-7*nnz)     [fixed 300 iters]
  qmv / qmv_plain     0.87 / 0.88 of cg_compare_seeded          [measured]
  refinement (x3)     0.55 of cg_compare_seeded                 [measured FLOOR --
                      these do not use the conv_iter+20 break, so on matrices that
                      converge early the true cost is higher]

Treat every figure as +/- 30%. The --time requests carry 1.4x headroom on top,
because Slurm kills a job that exceeds its limit by one second and the seed's work
is lost.
"""
import os
from math import ceil

# Slowdown vs the i7-12800H the cost model was fitted on, measured by probe.sbatch:
#   epyc-genoa (EPYC 9654P)  5.1 s -> 1.70
#   epyc-milan (EPYC 7643)   5.9 s -> 1.97
#   epyc       (EPYC 7643)   6.0 s -> 2.00   <- default; superset of epyc-milan
#   skylake    (Gold 6140)   7.8 s -> 2.60
#   haswell    (E5-2680 v3)  9.7 s -> 3.23
CLUSTER = float(os.environ.get("CLUSTER", 2.0))
HEADROOM = 1.4         # safety factor on --time
SEEDS = 50

# Longest --time the scheduler will accept, in hours. Confirm with `sinfo` and
# override if epyc differs:  MAX_WALL_H=48 python3 hpc/make_joblists.py
# A task whose padded estimate exceeds this is dropped rather than submitted to be
# killed: these programs do not checkpoint, so a TIMEOUT loses the whole seed.
MAX_WALL_H = float(os.environ.get("MAX_WALL_H", 168))

# Billing rate, dollars per core-hour. ARCSIM has quoted three different numbers
# (page prose $0.005, page table $5.10/1000, and $0.008), and the Internal/Federal
# tier differs from Other Academic and External. Confirm before quoting a cost.
RATE = float(os.environ.get("RATE", 0.008))

# name, n, nnz, K (the budget cg_compare_seeded actually reaches)
MATRICES = [
    ("s3dkq4m2", 90449, 4820891, 2000),
    ("s3dkt3m2", 90449, 3753461, 2000),
    ("nasasrb",  54870, 2677324, 2000),
    ("bcsstk37", 25503, 1140977, 2000),
    ("bcsstk36", 23052, 1143140, 2000),
    ("bcsstk38",  8032,  355460, 2000),
    ("nasa4704",  4704,  104756, 2000),
    ("bodyy4",   17546,  121938,  305),
    ("sts4098",   4098,   72356,  735),
    ("bcsstk14",  1806,   63454,  931),
    ("nos2",       957,    4137, 2000),
    ("mhd4800b",  4800,   27520,   61),
    ("bcsstk03",   112,     640,  265),
]
N    = {m[0]: m[1] for m in MATRICES}
NNZ  = {m[0]: m[2] for m in MATRICES}
KBUD = {m[0]: m[3] for m in MATRICES}
ALL  = [m[0] for m in MATRICES]


def t_sweep(m):
    return KBUD[m] * (2.245e-6 * NNZ[m] + 9.279e-5 * N[m]) * CLUSTER


def t_ladder(m):
    return 300 * (2.395e-5 * N[m] + 9.10e-7 * NNZ[m]) * CLUSTER


# MEASURED cost, seconds per seed on one epyc core, from the 2026-10-29 to 11-03
# run of all 5,259 tasks. Medians of 50 seeds, extracted with
# `bash hpc/fit_costs.sh`; the full table including maxima is in
# results/costs.txt.
#
# These supersede the formulas above wherever a matrix appears here, because the
# formulas were fitted on one laptop and are wrong in a way no scalar corrects:
# the old 0.55 x t_sweep multiplier for the refinement programs is right for
# bcsstk38 (0.53 implied) and 2.3x low for bcsstk36 and bcsstk37 (1.21--1.26
# implied). Real cost grew 7.3x between those matrices where t_sweep grows 3.05x.
# Underestimating that way killed 319 tasks at a 4:20 limit.
#
# Medians, not maxima. The maxima are contaminated by node-153, which had a bad
# temperature setting in its management controller and throttled its CPUs about
# 8x; it shows up as max/med ratios of 8--9 in otherwise stable classes. HEADROOM
# below is what covers a slow node, not an inflated central estimate.
MEASURED = {
    "sweep": {
        "s3dkq4m2": 72900, "s3dkt3m2": 62078, "nasasrb": 41568,
        "bcsstk37": 20334, "bcsstk36": 20045, "bcsstk38":  6158,
        "nasa4704":  2005, "bodyy4":     360, "sts4098":    302,
        "bcsstk14":   228, "nos2":       128, "mhd4800b":    24,
        "bcsstk03":     3,
    },
    "ladder": {
        "s3dkq4m2":  1465, "s3dkt3m2":  1465, "nasasrb":   852,
        "bcsstk37":   433, "bcsstk36":   346, "bcsstk38":  137,
        "nasa4704":    82, "bodyy4":     267, "sts4098":    66,
        "bcsstk14":    30, "nos2":        17, "mhd4800b":   41,
        "bcsstk03":     4,
    },
    "qmv": {
        "s3dkq4m2": 73406, "s3dkt3m2": 63173, "nasasrb": 44575,
        "bcsstk37": 18356, "bcsstk36": 18065, "bcsstk38":  6107,
        "nasa4704":  2024, "bodyy4":     368, "sts4098":    298,
        "bcsstk14":   231, "nos2":       131, "mhd4800b":    24,
        "bcsstk03":     3,
    },
    "qmvplain": {
        "s3dkq4m2": 67548, "s3dkt3m2": 61997, "nasasrb": 41504,
        "bcsstk37": 18503, "bcsstk36": 18242, "bcsstk38":  5726,
        "nasa4704":  1866, "bodyy4":     361, "sts4098":    300,
        "bcsstk14":   237, "nos2":       131, "mhd4800b":    24,
        "bcsstk03":     3,
    },
    # The four refinement programs agree within a few percent on a given matrix,
    # so they share one table: the worst median of the four, since the wall limit
    # has to cover the slowest.
    "refine": {
        "s3dkq4m2": 55770, "s3dkt3m2": 50189, "nasasrb": 44499,
        "bcsstk36": 19066, "bcsstk37": 17608, "bcsstk38":  3299,
        "nasa4704":  1700, "bodyy4":     210, "sts4098":    206,
        "bcsstk14":   164, "nos2":       125, "mhd4800b":    12,
        "bcsstk03":     3,
    },
    # maxiter = n. bcsstk36 and bcsstk37 are 57 and 72 hours a seed.
    "maxn": {
        "bcsstk37": 260438, "bcsstk36": 206021, "bcsstk38": 17746,
        "nasa4704":   4679, "bodyy4":       394, "sts4098":   297,
        "bcsstk14":    246, "nos2":          60, "mhd4800b":   24,
        "bcsstk03":      2,
    },
    "budget2k": {
        "s3dkq4m2": 144654, "s3dkt3m2": 133271, "nasasrb": 89728,
    },
}
modelled = []


def t_refine_fallback(m):
    return 0.55 * t_sweep(m)


def measured(prog, fallback):
    """Cost function for `prog`: the measured table, else `fallback`.

    Measurements are per-core seconds on epyc, which is CLUSTER = 2.0, so they
    rescale if CLUSTER changes for a different partition.
    """
    table = MEASURED.get(prog, {})

    def cost(m):
        if m in table:
            return table[m] * (CLUSTER / 2.0)
        if (prog, m) not in modelled:
            modelled.append((prog, m))
        return fallback(m)

    return cost


# program -> (binary, outdir, cost fn, matrices, seeds, extra arg per matrix or None)
#
# `maxn` passes maxiter = n as a fourth argument, which only cg_compare_seeded
# accepts. nasasrb, s3dkt3m2 and s3dkq4m2 are excluded from it: a single seed needs
# 304, 760 and 868 core-hours, past any wall limit and not restartable.
MAXN_MATRICES = [m for m in ALL if m not in ("nasasrb", "s3dkt3m2", "s3dkq4m2")]
BIG3 = ["s3dkq4m2", "s3dkt3m2", "nasasrb"]

PROGRAMS = {
    "sweep": dict(
        binary="bin/cg_compare_seeded", outdir="results/seed_sweep_50",
        cost=measured("sweep", t_sweep),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="5, 6, 8, 9"),
    "ladder": dict(
        binary="bin/generic_ladder_seeded", outdir="results/ladder_50",
        cost=measured("ladder", t_ladder),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="4"),
    # Table 4. A separate program from `ladder` because generic_ladder_seeded
    # seeds the residual from uniform_real(0.5, 1.5) rather than b = A x_true,
    # never records a double reference for r'z, and does not truncate at
    # convergence. See src/ladder_quiregain.cpp.
    "quiregain": dict(
        binary="bin/ladder_quiregain", outdir="results/quiregain_50",
        cost=measured("quiregain", lambda m: 0.4 * t_ladder(m)),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="4"),
    "qmv": dict(
        binary="bin/cg_compare_seeded_q512_qmv", outdir="results/qmv_50/qmv",
        cost=measured("qmv", lambda m: 0.87 * t_sweep(m)),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="10 (quire matvec arm)"),
    "qmvplain": dict(
        binary="bin/qmv_plain", outdir="results/qmv_50/plain",
        cost=measured("qmvplain", lambda m: 0.88 * t_sweep(m)),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="10 (plain matvec arm)"),
    # Program-to-table mapping follows TABLE_MAP.md. Table 11 comes from the v2
    # refinement program; v3_q512 is Table 12's QUIRE arm, paired with noquire as
    # its plain arm. Getting this wrong would leave Table 11 with no source and
    # Table 12 with only one of its two arms.
    "refine11": dict(
        binary="bin/cg_refinement_seeded_v2", outdir="results/refine11_50",
        cost=measured("refine", t_refine_fallback),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="11"),
    "ref12quire": dict(
        binary="bin/cg_refinement_seeded_v3_q512",
        outdir="results/refine12_quire_50",
        cost=measured("refine", t_refine_fallback),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="12 (quire arm)"),
    "ref12plain": dict(
        binary="bin/cg_refinement_noquire", outdir="results/refine12_plain_50",
        cost=measured("refine", t_refine_fallback),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="12 (plain arm)"),
    "refexactb": dict(
        binary="bin/cg_refinement_exactb_q512", outdir="results/refine_exactb_50",
        cost=measured("refine", t_refine_fallback),
        matrices=ALL, seeds=SEEDS, arg=None,
        tables="13"),
    # maxn and budget2k must NOT share seed_sweep_50. They write
    # <matrix>_seed<s>_it<N>.log there, and every aggregation script globs
    # <matrix>_seed*.log -- so their runs get pooled into the main sweep's
    # statistics. The maxn arm iterates to n and the naive arm diverges, which
    # turned Table 5's mean solution-error ratio into 1.1e29 before this was
    # caught.
    "maxn": dict(
        binary="bin/cg_compare_seeded", outdir="results/maxn_50",
        cost=measured("maxn", lambda m: (N[m] / KBUD[m]) * t_sweep(m)),
        matrices=MAXN_MATRICES, seeds=5, arg=lambda m: N[m],
        tables="iteration-budget check"),
    "budget2k": dict(
        binary="bin/cg_compare_seeded", outdir="results/budget2k_50",
        cost=measured("budget2k", lambda m: 2.0 * t_sweep(m)),
        matrices=BIG3, seeds=3,
        arg=lambda m: 4000,
        tables="iteration-budget check"),
}

# Size classes, by per-seed wall time. Separate --time per class so Slurm can
# backfill the cheap tasks into gaps left by the expensive ones.
CLASSES = [
    ("xl", 6 * 3600, float("inf"), 96),
    ("lg", 1800,     6 * 3600,     64),
    ("sm", 0,        1800,         96),
]


def hhmmss(seconds):
    s = int(ceil(seconds / 600.0) * 600)          # round up to 10 minutes
    s = max(s, 1800)                              # never ask for less than 30 min
    return "%02d:%02d:00" % (s // 3600, (s % 3600) // 60)


def mem_for(mx):
    """2 GB a task, which is what ARCSIM asked for after measuring the run.

    Peak MaxRSS across all 27 arrays of the first pass was 207 MB, on qmv_xl;
    187--194 MB for the other xl arrays and under 62 MB for lg. The buckets this
    replaces asked 4, 8 and 16 G by matrix size, up to 80x over, and those
    requests rather than cores were what kept tasks from starting: an epyc node
    has 503 GB, so 16 G a task caps it at 31 tasks however many cores sit idle.

    Scaling below 2 GB would buy nothing -- 128 tasks at 2 GB is half of one
    node's memory -- so the argument is accepted and ignored.
    """
    del mx
    return "2G"


os.makedirs("hpc/jobs", exist_ok=True)
rows, grand, dropped = [], 0.0, []

for prog, cfg in PROGRAMS.items():
    for cname, lo, hi, cap in CLASSES:
        mx = [m for m in cfg["matrices"] if lo <= cfg["cost"](m) < hi]
        too_long = [m for m in mx if cfg["cost"](m) * HEADROOM > MAX_WALL_H * 3600]
        for m in too_long:
            dropped.append((prog, m, cfg["cost"](m) / 3600.0))
        mx = [m for m in mx if m not in too_long]
        if not mx:
            continue
        mx.sort(key=cfg["cost"], reverse=True)    # longest first within the array
        path = "hpc/jobs/%s_%s.txt" % (prog, cname)
        with open(path, "w") as f:
            for m in mx:
                for s in range(1, cfg["seeds"] + 1):
                    extra = "" if cfg["arg"] is None else " %d" % cfg["arg"](m)
                    f.write("%s %d%s\n" % (m, s, extra))
        ntasks = len(mx) * cfg["seeds"]
        longest = max(cfg["cost"](m) for m in mx)
        ch = sum(cfg["cost"](m) for m in mx) * cfg["seeds"] / 3600.0
        grand += ch
        rows.append(dict(cls="%s_%s" % (prog, cname), prog=prog,
                         binary=cfg["binary"], outdir=cfg["outdir"],
                         joblist=path, ntasks=ntasks,
                         time=hhmmss(longest * HEADROOM), mem=mem_for(mx),
                         cap=min(cap, ntasks), coreh=ch, longest=longest / 3600.0))

# Longest-processing-time first: a task cannot be split, so the longest one floors
# the finish time. Submitting it first lets cheap tasks fill the gaps around it.
rows.sort(key=lambda r: -r["longest"])

with open("hpc/jobs/manifest.tsv", "w") as f:
    f.write("#class\tbinary\toutdir\tjoblist\tntasks\ttime\tmem\tcap\n")
    for r in rows:
        f.write("%(cls)s\t%(binary)s\t%(outdir)s\t%(joblist)s\t%(ntasks)d\t"
                "%(time)s\t%(mem)s\t%(cap)d\n" % r)

print("%-18s %6s %10s %8s %7s %9s" %
      ("class", "tasks", "--time", "mem", "cap", "core-h"))
for r in rows:
    print("%-18s %6d %10s %8s %7d %9.1f" %
          (r["cls"], r["ntasks"], r["time"], r["mem"], r["cap"], r["coreh"]))
print("\n%d arrays, %d tasks, %.0f core-hours, about $%.2f at $%g/core-hour"
      % (len(rows), sum(r["ntasks"] for r in rows), grand, grand * RATE, RATE))

if modelled:
    print("\nNot measured, so still using the fitted formula:")
    for prog, m in sorted(modelled):
        print("  %-12s %s" % (prog, m))
    print("  The formulas were fitted on one laptop and were 2.3x low on the")
    print("  mid-size matrices. Measure with: bash hpc/fit_costs.sh, then add")
    print("  the rows to MEASURED above.")

if dropped:
    print("\nDropped: a single seed exceeds the %g h wall limit and these programs"
          "\ndo not checkpoint, so a TIMEOUT would lose the whole run." % MAX_WALL_H)
    for prog, m, h in sorted(dropped, key=lambda d: -d[2]):
        print("  %-10s %-10s %6.1f h per seed" % (prog, m, h))
    print("  Raise with MAX_WALL_H=<hours> if sinfo allows more.")
print("\nTables covered:")
for prog, cfg in PROGRAMS.items():
    print("  %-10s -> %s" % (prog, cfg["tables"]))
