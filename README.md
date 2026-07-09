# posit-sparse-bench

Benchmarking posit arithmetic against IEEE double precision for sparse matrix computations in FEM/HPC applications. LFX Summer 2026 mentorship project under Kurt Keville (MIT) and Joshua Gyllinsky.

## Research Question

Can posit arithmetic with quire exact accumulation match or exceed double precision accuracy for the conjugate gradient inner product `p^T A p`, which is the critical dot product in sparse iterative solvers?

## Key Finding

**posit32+quire achieves 6.6x–4,531x lower error than posit32 naive accumulation across all tested matrices.** The quire's exact accumulation — not just wider precision — is the primary driver of accuracy.

posit32+quire maintains relative error below 4e-2 across 300 CG iterations on all tested matrices (below 1e-2 on 11 of 13; two matrices near precision floor reach ~2-3e-2); below 1e-6 on well-conditioned matrices. posit16 is unreliable on high dynamic range matrices; behavior is matrix-dependent (marginal on bcsstk38/nasasrb, catastrophic failure on others).

## Benchmark Method

For each CG iteration, we compute the `p^T A p` dot product simultaneously in:
- posit8, posit16, posit32, posit64 — each with quire (exact) and naive accumulation
- double64 — used as ground truth reference

Relative error = `|posit_result - double64| / |double64|`

All matrices are real symmetric from [SuiteSparse Matrix Collection](https://sparse.tamu.edu).


## Implementation Details

**Quire configuration:** `quire<N,ES,2>` throughout — capacity parameter 2 gives a 482-bit quire for posit32, derived directly from the Universal library source (qbits = escale*(4*nbits-8) + capacity, where escale = 2^es). Capacity=2 formally guarantees exact accumulation of up to 4 worst-case-magnitude (maxpos^2) terms; in practice, real matrix entries are far below maxpos^2 in magnitude, so overflow does not occur even with hundreds of thousands of accumulations (e.g., bcsstk38's 355,460 nonzeros).

**CG solver:** Jacobi-preconditioned CG — diagonal preconditioner `z[i] = r[i]/diagA[i]`. 300 iterations per run.

**es values:** Per the 2022 Posit Standard (es=2 for all sizes), confirmed directly with John Gustafson — thank you to Prof. Gustafson for clarifying this. Our ladder results below include both configurations: the legacy variable-es convention (posit8 es=0, posit16 es=1, posit32/64 es=2) used in earlier posit literature, and the standard-compliant es=2 uniform configuration, to show how the fixed exponent size affects posit8/16 behavior relative to the older convention.
## Results

### Quire Improvement (posit32 naive vs posit32+quire), 300 CG iterations

| Matrix | Domain | Diag ratio | p32+quire max err | p32 naive max err | Quire gain |
|--------|--------|-----------|------------------|------------------|------------|
| bcsstk03 | FEM structural | 10^6 | 3.22e-02* | 2.12e-01 | 6.6x |
| bcsstk14 | FEM structural | 10^10 | 3.33e-06 | 1.31e-04 | 39x |
| bcsstk38 | FEM structural | 10^12 | 3.47e-08 | 1.95e-06 | 56x |
| nasasrb | NASA structural | 10^5 | 9.28e-09 | 7.37e-07 | 79x |
| mhd4800b | MHD physics | 10^12 | 1.78e-02 | 4.01e-01 | 23x |
| s3dkt3m2 | 3D structural FEM | 10^7 | 1.12e-07 | 1.72e-04 | 1,535x |

*bcsstk03: high error confined to iters 181–201 where |pAp| < 1e-30 (near-convergence regime, 98 iters filtered). posit32 precision floor reached.

### Precision Ladder

![Precision ladder](results/figures/precision_ladder_6mat.png)

- **posit8**: catastrophic failure on all matrices (error 10^3–10^28)
- **posit16**: fails on wide dynamic range matrices; marginal on bcsstk38/nasasrb
- **posit32+quire**: robust across all domains, max error 9.28e-09 to 3.22e-02 (two matrices near the posit32 precision floor account for the upper end; well-conditioned matrices stay below 1e-6)

### posit32+quire vs posit32 naive

![Quire vs naive](results/figures/p32_quire_vs_naive_6mat.png)

## Matrices Tested

See Extended Results table below for full 13-matrix list with properties.

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

results/                # raw per-iteration logs (300 iters each, gitignored — regenerate via run_all.sh)

results/csv/            # summary statistics

results/figures/        # plots

data/matrices/          # 13 tested matrices (see Excluded Matrices for others present but not used)
## Dependencies

- [universal](https://github.com/stillwater-sc/universal) — posit arithmetic library
- g++ with C++20
- Python3, scipy, matplotlib (for analysis and plots)

## Build

```bash
g++ -O2 -std=c++20 -I/path/to/universal/include src/generic_ladder.cpp -o generic_ladder
./generic_ladder data/matrices/bcsstk38.mtx results/bcsstk38_ladder.log
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

## Extended Results (13 matrices)

| Matrix | n | Diag ratio | Val ratio | p32q max err | p32 naive max err | Quire gain |
|--------|---|-----------|-----------|-------------|------------------|------------|
| bcsstk03 | 112 | 1.52e+06 | 3.78e+16 | 3.22e-02* | 2.12e-01 | 6.6x |
| bcsstk14 | 1806 | 8.94e+09 | 5.53e+36 | 3.33e-06 | 1.31e-04 | 39x |
| bcsstk36 | 23052 | 1.74e+09 | 1.13e+37 | 1.56e-08 | 1.73e-06 | 111x |
| bcsstk37 | 25503 | 9.61e+08 | 8.81e+27 | 3.51e-08 | 3.27e-06 | 93x |
| bcsstk38 | 8032 | 9.26e+12 | 2.00e+38 | 3.47e-08 | 1.95e-06 | 56x |
| nasasrb | 54870 | 2.55e+05 | 2.32e+22 | 9.28e-09 | 7.37e-07 | 79x |
| mhd4800b | 4800 | 3.73e+12 | 5.75e+20 | 1.78e-02 | 4.01e-01 | 23x |
| s3dkt3m2 | 90449 | 2.52e+07 | 1.01e+40 | 1.12e-07 | 1.72e-04 | 1,535x |
| s3dkq4m2 | 90449 | 1.44e+07 | 1.62e+26 | 1.11e-07 | 5.05e-04 | 4,531x |
| sts4098 | 4098 | 6.02e+07 | 5.66e+54 | 3.65e-07 | 1.43e-05 | 39x |
| nasa4704 | 4704 | 1.52e+05 | 2.60e+20 | 4.67e-08 | 3.85e-07 | 8x |
| nos2 | 957 | 1.23e+05 | 2.46e+05 | 7.99e-07 | 7.72e-06 | 10x |
| bodyy4 | 17546 | 2.45e+02 | 1.84e+19 | 6.23e-03 | 7.96e-01 | 128x |

*bcsstk03: high error confined to near-convergence regime (iters 181-201, pAp < 1e-30)


**Key finding:** No single matrix property (diagonal ratio or value range) cleanly predicts quire gain. sts4098 has the highest value range (1e+54) yet relatively low quire gain (39x), suggesting quire benefit depends on the interaction of value distribution, matrix size, and CG search direction evolution.

**posit16 failure predictor — negative result:** We tested whether posit16 (es=1 or es=2) failure could be predicted from value ratio, diagonal ratio, matrix size (n), or value-ratio-per-n. None separate failing matrices (bcsstk03, mhd4800b, bodyy4, bcsstk14) from passing matrices (bcsstk36, bcsstk37, bcsstk38, nasasrb, s3dkt3m2, s3dkq4m2, sts4098, nasa4704, nos2) cleanly. Notably, bcsstk38 has the highest value-ratio-per-n (2.49e+34) of any tested matrix yet does not fail, while mhd4800b fails with a value-ratio-per-n three orders of magnitude lower (1.19e+17). This extends the quire-gain finding above: arithmetic reliability under posit16 depends on the interaction of value distribution shape and CG search direction evolution across iterations, not a static summary statistic of the matrix.

### posit16 saturation mechanism (bcsstk38 vs mhd4800b)

To probe this further, we plotted the actual pAp trajectory (ground-truth double64
value vs. posit16+quire vs. posit16 naive) per iteration for both matrices:

![posit16 saturation comparison](results/figures/p16_saturation_compare.png)

bcsstk38's posit16 tracks double64 closely across all iterations — pAp stays within
posit16's representable range. mhd4800b's posit16 (both quire and naive) saturates
at a fixed ceiling (~2.68e+08) from iteration 0 onward, while the true pAp value
climbs to ~6e+14. This is a hard dynamic-range ceiling, not a rounding-error effect:
posit16's max representable magnitude (~2^28 with es=2) is far below mhd4800b's
actual pAp magnitudes, so it saturates immediately rather than degrading gradually.


![Property vs quire gain](results/figures/property_vs_quire_gain.png)

## Full CG Solver Convergence

Beyond measuring a single inner product, we ran complete CG solvers in double64, float32, and posit32+quire simultaneously across 8 matrices, tracking residual norm per iteration.

![CG convergence comparison](results/figures/cg_convergence_compare.png)

Iteration-to-converge (residual < 1e-10), verified per-iteration from raw logs:

| Matrix | double64 | float32 | posit32+quire |
|---|---|---|---|
| bcsstk03 | converges iter 198 | converges iter 382 | enters bounded precision-floor regime (residual 1e-10 to 1e-9, never above 1e-9 after iter 450) from ~iter 450 through 1000 iterations (extended run) — no divergence, no further improvement |
| mhd4800b | converges iter 55 | converges iter 69 | converges iter 69 — exact match with float32 |
| bcsstk14, bcsstk36, bcsstk37, bcsstk38, nasasrb, sts4098 | does not converge (500 iters) | does not converge | does not converge — posit32+quire does not perform worse than double64 or float32 on ill-conditioned matrices |

Key observations:
- **bcsstk14**: posit32+quire convergence curve is visually indistinguishable from double64 — a drop-in replacement result
- **sts4098**: float32 diverges from double64 after iter 200; posit32+quire tracks double64 to the end — posit32+quire beats float32
- **mhd4800b**: all three converge in 70 iterations; posit32+quire reaches 1e-10 vs double64's 1e-13
- **bcsstk03**: double64 converges to 1e-38; posit32+quire enters a bounded precision-floor regime (no divergence), confirmed stable through an extended 1000-iteration run
- **bcsstk38**: ill-conditioned, none converge — but posit32+quire does not make behavior worse

posit32+quire matches or exceeds float32 behavior across all tested matrices in full solver context.

## Reproducing Results (Docker)

Requirements: Docker, git

```bash
git clone https://github.com/Gurleen-kansray/posit-sparse-bench
cd posit-sparse-bench
docker build -t posit-bench .
docker run --rm posit-bench bash run_all.sh
```

Environment: Ubuntu 22.04, g++ 11, Universal v3.80, quire<N,ES,2>, 300 CG iterations per matrix.

## Divergence Analysis

Per-iteration relative error tracking (posit32 quire vs naive, against double64 reference) across 300 CG iterations for all 13 test matrices. Divergence point defined as naive error exceeding quire error by >10x, sustained for 5+ iterations.

| Matrix | Divergence Iter | Max Err (Quire) | Max Err (Naive) | Gain (Max) |
|---|---|---|---|---|
| bcsstk03 | none (floors out) | 3.22e-02 | 2.12e-01 | 6.6x |
| bcsstk14 | 32 | 3.33e-06 | 1.31e-04 | 39.2x |
| bcsstk36 | 6 | 1.56e-08 | 1.73e-06 | 110.7x |
| bcsstk37 | 0 | 3.51e-08 | 3.27e-06 | 93.2x |
| bcsstk38 | 0 | 3.47e-08 | 1.95e-06 | 56.3x |
| bodyy4 | 0 | 6.23e-03 | 7.96e-01 | 127.8x |
| mhd4800b | 6 | 1.77e-02 | 4.01e-01 | 22.6x |
| nasa4704 | 0 | 4.67e-08 | 3.85e-07 | 8.2x |
| nasasrb | 0 | 9.28e-09 | 7.37e-07 | 79.4x |
| nos2 | 2 | 7.99e-07 | 7.72e-06 | 9.7x |
| s3dkq4m2 | 0 | 1.11e-07 | 5.05e-04 | 4531.5x |
| s3dkt3m2 | 0 | 1.12e-07 | 1.72e-04 | 1535.0x |
| sts4098 | 12 | 3.65e-07 | 1.43e-05 | 39.1x |

## Divergence Mechanism (mhd4800b) - Confirmed via Controlled Isolation

We investigated why naive posit32 (no quire) fails to match float32 on mhd4800b's full CG solver (converges at iteration 79 vs float32/posit32+quire at 69), through four successive tests - two hypotheses ruled out, one confirmed causally, not just correlationally.

**Hypothesis 1 (ruled out): individual term-magnitude precision loss.** We measured what fraction of individual p[i]*Ap[i] terms fall outside posit32's precision-favorable magnitude zone (empirically ~3.16e-5 to ~1e5, from an analytical posit32-vs-float32 relative-precision sweep, src/posit_precision_curve.cpp). At the iterations where naive posit32 first diverges (15-25), only 0.04-2.9% of individual terms fall outside this zone - this cannot explain the gap (src/term_probe.cpp).

**Hypothesis 2 (ruled out): catastrophic cancellation.** We measured the ratio of the largest running partial sum to the final accumulated result during pAp accumulation. This ratio stayed ~1.0 for both float32 and posit32-naive at iterations 20-24 - no significant cancellation was occurring (src/cancellation_probe.cpp).

**Confirmed mechanism: magnitude-dependent rounding of the accumulated pAp scalar, compounding through CG's recurrence.** Tracking the relative L2 drift of the p search-direction vector against a double64 ground truth, from iteration 0 onward (src/trajectory_probe.cpp), showed posit32-naive's trajectory deviates from ground truth more than float32's at every iteration from the start - not a sudden fork at iteration 19, but a gradual, compounding divergence that only becomes visible in the residual norm once it's accumulated enough.

To isolate the cause, we ran a controlled hybrid experiment (src/hybrid_probe.cpp): three otherwise-identical double64 CG solvers, differing ONLY in how the scalar pAp is rounded each iteration (unrounded control; rounded via float32; rounded via posit32). With everything else held at double64, posit32's rounding of pAp alone produced 30-45x more downstream trajectory deviation than float32's rounding of the same value, sustained across iterations 0-18.

Finally, we directly compared single-value rounding error of the true pAp scalar (from an unperturbed double64 trajectory) through posit32 vs float32 at its actual observed magnitude each iteration (src/pAp_rounding_probe.cpp), and correlated this against the precision-curve zone from Hypothesis 1's analysis:

- Iterations 0-16: pAp magnitude ranges 1e9-1e13, OUTSIDE posit32's favorable zone (upper bound ~1e5). Posit32's rounding error is 10-90x larger than float32's here.
- Iterations 17-18: pAp crosses into the zone boundary (~1e5-1e6); rounding errors converge (ratio ~1.0).
- Iterations 19+: pAp drops to 10^2-10^4, INSIDE posit32's favorable zone. Posit32 becomes MORE accurate than float32 per-step (ratio 0.001-0.22) - but by this point the earlier compounded error has already been baked into the trajectory.

**Summary:** naive posit32's disadvantage is not caused by individual-term precision loss or cancellation, but by the ACCUMULATED pAp scalar sitting outside posit32's precision-favorable magnitude range during CG's early iterations (when pAp is largest), producing measurably larger single-step rounding error than float32 during exactly this window. This early-iteration error compounds through CG's own recurrence (each iteration's search direction depends on the previous iteration's rounding error), producing the observed ~10-iteration convergence lag - even though posit32 would out-precision float32 later in the same run, once pAp's magnitude shrinks into its favorable zone. Quire eliminates this entirely by accumulating pAp exactly regardless of magnitude, which is why posit32+quire matches float32/double64 from iteration 0.

Full per-iteration data: results/posit_precision_curve.log, results/term_probe_mhd4800b.log, results/term_probe_bcsstk38.log, results/trajectory_mhd4800b.log, results/hybrid_mhd4800b.log, results/pAp_rounding_mhd4800b.log
