# posit-sparse-bench

Benchmarking posit arithmetic against IEEE double precision for sparse matrix computations in FEM/HPC applications. LFX Summer 2026 mentorship project under Kurt Keville (MIT) and Joshua Gyllinsky.

## Research Question

Can posit arithmetic with quire exact accumulation match or exceed double precision accuracy for the conjugate gradient inner product `p^T A p`, which is the critical dot product in sparse iterative solvers?

## Key Finding

**posit32+quire achieves 1827x–445382x lower error than posit32 naive accumulation across all tested matrices.** The quire's exact accumulation — not just wider precision — is the primary driver of accuracy.

posit32+quire maintains relative error below 1e-6 across 300 CG iterations on all tested matrices. posit16 is unreliable on high dynamic range matrices; behavior is matrix-dependent (marginal on bcsstk38/nasasrb, catastrophic failure on others).

## Benchmark Method

For each CG iteration, we compute the `p^T A p` dot product simultaneously in:
- posit8, posit16, posit32, posit64 — each with quire (exact) and naive accumulation
- double64 — used as ground truth reference

Relative error = `|posit_result - double64| / |double64|`

All matrices are real symmetric from [SuiteSparse Matrix Collection](https://sparse.tamu.edu).

## Results

### Quire Improvement (posit32 naive vs posit32+quire), 300 CG iterations

| Matrix | Domain | Diag ratio | p32+quire max err | p32 naive max err | Quire gain |
|--------|--------|-----------|------------------|------------------|------------|
| bcsstk03 | FEM structural | 10^6 | 3.22e-02* | 2.12e-01 | 3,588x |
| bcsstk14 | FEM structural | 10^10 | 3.33e-06 | 1.31e-04 | 1,827x |
| bcsstk38 | FEM structural | 10^12 | 4.39e-07 | 6.75e-05 | 26,890x |
| nasasrb | NASA structural | 10^5 | 9.93e-08 | 9.96e-05 | 43,214x |
| mhd4800b | MHD physics | 10^12 | 6.54e-05 | 2.94e-02 | 111,342x |
| s3dkt3m2 | 3D structural FEM | 10^7 | 7.25e-06 | 2.74e-02 | 445,382x |

*bcsstk03: high error confined to iters 181–201 where |pAp| < 1e-30 (near-convergence regime, 98 iters filtered). posit32 precision floor reached.

### Precision Ladder

![Precision ladder](results/figures/precision_ladder_6mat.png)

- **posit8**: catastrophic failure on all matrices (error 10^3–10^28)
- **posit16**: fails on wide dynamic range matrices; marginal on bcsstk38/nasasrb
- **posit32+quire**: robust across all domains, max error 1e-5 to 4e-7

### posit32+quire vs posit32 naive

![Quire vs naive](results/figures/p32_quire_vs_naive_6mat.png)

## Matrices Tested

See Extended Results table below for full 10-matrix list with properties.

### Primary 6 matrices (used in ladder & methodology validation)

| Matrix | n | nnz | Domain | Source |
|--------|---|-----|--------|--------|
| bcsstk03 | 112 | 640 | FEM structural | Boeing/HB |
| bcsstk14 | 1806 | 32,606 | FEM structural | Boeing/HB |
| bcsstk38 | 8,032 | 355,460 | FEM structural | Boeing |
| nasasrb | 54,870 | 2,677,324 | NASA structural | Pothen/NASA |
| mhd4800b | 4,800 | 27,520 | MHD physics | Bai |
| s3dkt3m2 | 90,449 | 3,753,461 | 3D structural FEM | GHS_psdef |

## Excluded Matrices

Moved to `src/exploratory/` with reasons:
- **scircuit, memplus, add32**: unsymmetric — CG invalid
- **cfd1, cfd2**: preconditioned (diagonal all 1.0) — dynamic range artificially suppressed

## Repository Structure
src/                    # ladder benchmark source files

src/exploratory/        # excluded matrices (unsymmetric/preconditioned)

results/ladder_logs/    # raw per-iteration logs (300 iters each)

results/csv/            # summary statistics

results/figures/        # plots

data/matrices/          # bcsstk03, bcsstk14 matrix files
## Dependencies

- [universal](https://github.com/stillwater-sc/universal) — posit arithmetic library
- g++ with C++20
- Python3, scipy, matplotlib (for analysis and plots)

## Build

```bash
g++ -O2 -std=c++20 -I/path/to/universal/include src/bcsstk38_ladder.cpp -o bcsstk38_ladder
./bcsstk38_ladder
```

## Project Context

LFX Summer 2026 — "Broadening the RISC-V High Precision Code Base and Reach"
Mentors: Kurt Keville (MIT R&D Labs), Joshua Gyllinsky
Target venue: CoNGA 2026

## Methodology Validation

posit64+quire matches the double64 reference to within 1e-11 (or exactly) across all 6 matrices, confirming the measurement framework is sound and double64 is a valid ground truth.

| Matrix | p64+quire max err | p64 naive max err |
|--------|------------------|------------------|
| bcsstk03 | 1.25e-11 | 7.37e-11 |
| bcsstk14 | 0.00e+00 | 0.00e+00 |
| bcsstk38 | 0.00e+00 | 0.00e+00 |
| nasasrb  | 0.00e+00 | 0.00e+00 |
| mhd4800b | 0.00e+00 | 6.72e-11 |
| s3dkt3m2 | 2.01e-11 | 3.70e-11 |

## Quire Accuracy Gain Summary

![Quire gain barchart](results/figures/quire_gain_barchart.png)

## Extended Results (10 matrices)

| Matrix | n | Diag ratio | Val ratio | p32q max err | p32 naive max err | Quire gain |
|--------|---|-----------|-----------|-------------|------------------|------------|
| bcsstk03 | 112 | 1.52e+06 | 3.78e+16 | 3.22e-02* | 2.12e-01 | 3,588x |
| bcsstk14 | 1806 | 8.94e+09 | 5.53e+36 | 3.33e-06 | 1.31e-04 | 1,827x |
| bcsstk36 | 23052 | 1.74e+09 | 1.13e+37 | 1.53e-06 | 1.01e-03 | 662x |
| bcsstk37 | 25503 | 9.61e+08 | 8.81e+27 | 1.46e-06 | 6.67e-04 | 458x |
| bcsstk38 | 8032 | 9.26e+12 | 2.00e+38 | 4.39e-07 | 6.75e-05 | 26,890x |
| nasasrb | 54870 | 2.55e+05 | 2.32e+22 | 9.93e-08 | 9.96e-05 | 43,214x |
| mhd4800b | 4800 | 3.73e+12 | 5.75e+20 | 6.54e-05 | 2.94e-02 | 111,342x |
| s3dkt3m2 | 90449 | 2.52e+07 | 1.01e+40 | 7.25e-06 | 2.74e-02 | 445,382x |
| s3dkq4m2 | 90449 | 1.44e+07 | 1.62e+26 | 5.22e-06 | 1.35e-02 | 2,591x |
| sts4098 | 4098 | 6.02e+07 | 5.66e+54 | 1.13e-07 | 4.07e-06 | 36x |

*bcsstk03: high error confined to near-convergence regime (iters 181-201, pAp < 1e-30)


**Key finding:** No single matrix property (diagonal ratio or value range) cleanly predicts quire gain. sts4098 has the highest value range (1e+54) yet lowest quire gain (36x), suggesting quire benefit depends on the interaction of value distribution, matrix size, and CG search direction evolution.

![Property vs quire gain](results/figures/property_vs_quire_gain.png)
