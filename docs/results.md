# Results

Full tables corresponding to the CoNGA26 submission "Exact Accumulation Improves the Inner Product but Not the Attainable Accuracy of Conjugate Gradient." Table numbers match the paper.

## Table 1 — Test suite (13 SuiteSparse matrices)

| Matrix | n | Diagonal ratio | cond(A) |
|---|---|---|---|
| bcsstk03 | 112 | 1.52e6 | 6.79e6 |
| bcsstk14 | 1,806 | 8.94e9 | 1.19e10 |
| bcsstk36 | 23,052 | 1.74e9 | 7.43e11 |
| bcsstk37 | 25,503 | 9.61e8 | 6.87e13 |
| bcsstk38 | 8,032 | 9.26e12 | 5.52e16 |
| bodyy4 | 17,546 | 2.45e2 | 8.06e2 |
| mhd4800b | 4,800 | 3.73e12 | 8.17e13 |
| nasa4704 | 4,704 | 1.52e5 | 4.17e7 |
| nasasrb | 54,870 | 2.55e5 | 5.58e8 |
| nos2 | 957 | 1.23e5 | 5.10e9 |
| s3dkq4m2 | 90,449 | 1.44e7 | 1.90e11 |
| s3dkt3m2 | 90,449 | 2.52e7 | 3.63e11 |
| sts4098 | 4,098 | 6.03e7 | 2.17e8 |

## Table 3 — Quire gain, per-iteration (mean, 50-seed sweep, James's pipeline)

Source: external/james_replication/ladder_seeded_summary.csv

| Matrix | pAp gain (mean) | r'z gain (mean) |
|---|---|---|
| bcsstk03 | 13.35 | 12.25 |
| bcsstk14 | 30.96 | 74.66 |
| bcsstk36 | 190.03 | 422.69 |
| bcsstk37 | 480.46 | 822.76 |
| bcsstk38 | 35.58 | 64.79 |
| bodyy4 | 81.40 | 338.86 |
| mhd4800b | 370.52 | 542.02 |
| nasa4704 | 259.62 | 249.31 |
| nasasrb | 2430.56 | 2840.64 |
| nos2 | 27.60 | 23.79 |
| s3dkq4m2 | 2819.43 | 3846.79 |
| s3dkt3m2 | 1609.77 | 2861.99 |
| sts4098 | 443.38 | 430.01 |

Median across matrices: **260x (pAp)**, **423x (r'z)**.

## Table 4 — Non-transfer: solution-error ratio (quire/naive), James's pipeline

Source: external/james_replication/ladder_seeded_summary.csv

| Matrix | Seeds | Mean | Std | Median |
|---|---|---|---|---|
| bcsstk03 | 50 | 1.1974 | 0.8726 | 0.9761 |
| bcsstk14 | 26 | 1.0268 | 0.1959 | 1.0121 |
| bcsstk36 | 10 | 1.0051 | 0.1027 | 1.0264 |
| bcsstk37 | 10 | 1.0565 | 0.1218 | 1.0682 |
| bcsstk38 | 12 | 0.9915 | 0.1161 | 0.9735 |
| bodyy4 | 25 | 0.9921 | 0.0726 | 0.9772 |
| mhd4800b | 50 | 1.4245 | 1.2045 | 1.0501 |
| nasa4704 | 18 | 1.0445 | 0.3836 | 0.9673 |
| nos2 | 40 | 1.0940 | 0.6017 | 0.9161 |
| sts4098 | 27 | 0.9745 | 0.2670 | 0.9689 |
| s3dkq4m2 | 8 | 1.0022 | 0.0065 | 1.0011 |
| s3dkt3m2 | 8 | 1.0024 | 0.0056 | 1.0038 |
| nasasrb | 10 | 0.9372 | 0.1523 | 0.9562 |

## Table 5 — Pre-convergence pAp gain (this repo's own 2000-iter pipeline)

Source: results/csv/ (this repo's seed-sweep logs, restricted to pre-1e-6-residual iterations)

| Matrix | Seeds | Pre-conv gain (mean) | Pre-conv gain (std) |
|---|---|---|---|
| bcsstk03 | 50 | 6.64 | 8.44 |
| bcsstk14 | 26 | 28.63 | 25.48 |
| bcsstk36 | 5 | 84.45 | 48.66 |
| bcsstk37 | 4 | 263.30 | 273.14 |
| bcsstk38 | 12 | 21.04 | 20.09 |
| bodyy4 | 25 | 163.91 | 303.45 |
| mhd4800b | 50 | 44.97 | 150.85 |
| nasa4704 | 5 | 26.01 | 23.07 |
| nos2 | 29 | 5.50 | 8.77 |
| s3dkq4m2 | 8 | 219.20 | 68.96 |
| s3dkt3m2 | 8 | 1032.50 | 1121.61 |
| sts4098 | 27 | 177.98 | 231.36 |

nasasrb excluded — does not reach 1e-6 relative residual within 2000 iterations for any seed.

## Table 6 — Independent-pipeline TOST equivalence (6-matrix subset)

Source: external/james_replication/stats_equivalence.csv

| Matrix | Pairs | Ratio | 95% CI | TOST p (10%) |
|---|---|---|---|---|
| bcsstk03 | 50 | 0.9898 | [0.844, 1.162] | 0.145 |
| bcsstk14 | 50 | 1.0262 | [0.982, 1.072] | 0.0013 |
| mhd4800b | 50 | 1.0804 | [0.836, 1.396] | 0.444 |
| nasa4704 | 50 | 0.9959 | [0.995, 0.997] | ~0 |
| nos2 | 50 | 1.0013 | [0.999, 1.004] | ~0 |
| sts4098 | 50 | 1.0585 | [0.997, 1.124] | 0.1024 |
| **Pooled (6 matrices)** | 300 | 1.0248 | [0.974, 1.078] | 0.0031 |

## Table 7 — Matvec-quire solution-error ratio (plain/quire matvec, inner products fixed at quire)

**Data status: not yet committed to this repo — see TODO below.**

| Matrix | Seeds | Plain err. (mean) | Quire err. (mean) | Ratio |
|---|---|---|---|---|
| bcsstk03 | 50 | 0.01687 | 0.01128 | 1.50 |
| bcsstk14 | 20 | 0.000263 | 0.000168 | 1.56 |
| bodyy4 | 20 | 2.60e-7 | 1.99e-7 | 1.31 |
| mhd4800b | 50 | 0.000708 | 0.000380 | 1.86 |
| nasa4704 | 20 | 0.13808 | 0.13803 | 1.00 |
| nos2 | 50 | 0.32779 | 0.32702 | 1.00 |
| sts4098 | 20 | 0.000556 | 0.0000817 | 6.81 |

TODO: commit quired_matvec.csv to external/james_replication/ (currently missing).

## Table 8 — Iterative refinement improvement factor (posit32 correction solve)

**Data status: code exists (src/cg_refinement_seeded_v2.cpp), result CSV not yet committed — see TODO below.**

| Matrix | Seeds | Imp. (mean) | Imp. (std) | Seeds worse |
|---|---|---|---|---|
| bcsstk03 | 20 | 3.33 | 2.34 | 2 |
| bcsstk14 | 20 | 36.58 | 8.32 | 0 |
| bcsstk36 | 20 | 2.17 | 0.29 | 0 |
| bcsstk37 | 20 | 4.20 | 0.41 | 0 |
| bcsstk38 | 20 | 3.32 | 1.12 | 0 |
| bodyy4 | 20 | 209.99 | 18.19 | 0 |
| mhd4800b | 20 | 18,559.72 | 12,381.65 | 0 |
| nasa4704 | 20 | 3.99 | 1.30 | 0 |
| nasasrb | 10 | 3.16 | 1.16 | 0 |
| nos2 | 20 | 1.48 | 0.50 | 2 |
| s3dkq4m2 | 12 | 3.03 | 0.53 | 0 |
| s3dkt3m2 | 10 | 3.07 | 0.34 | 0 |
| sts4098 | 20 | 38.66 | 26.16 | 0 |

228 of 232 seed-runs improve; only bcsstk03 and nos2 show occasional regressions.

TODO: run src/cg_refinement_seeded_v2.cpp end to end (or locate existing output logs) and commit the resulting CSV to results/csv/refinement_seeded_summary.csv, then remove this TODO note.
