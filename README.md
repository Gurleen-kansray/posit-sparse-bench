# posit-sparse-bench

Benchmarking posit arithmetic (with quire exact accumulation) against IEEE double and float32 in sparse conjugate gradient (CG) solvers, on real symmetric matrices from the SuiteSparse collection. LFX Summer 2026 mentorship project, targeting CoNGA-Q@SC26.

## Research Question

Does exact accumulation (the posit quire) at each of CG's three reduction sites — the curvature term p^T A p, the residual inner product r^T z, and the sparse matrix-vector product — improve the solver's final solution accuracy, and if not everywhere, where does it actually pay off?

## Headline Findings

1. **Quire vs naive posit32 (inner products):** posit32+quire achieves large per-iteration accuracy gains at both CG inner products — median 260x lower error on p^T A p and 423x on r^T z across 13 SuiteSparse matrices (50-seed sweep), against a double-precision reference.
2. **Non-transfer to the solution.** Despite those per-iteration gains, the final solution-error ratio (posit32+quire / posit32 naive) clusters near 1.0 across all 13 matrices (median 0.92-1.05). Independently confirmed via a formal TOST equivalence test on a 6-matrix subset (pooled ratio 1.025, p=0.003, 300 paired runs) — see external/james_replication/.
3. **Matrix-vector accumulation behaves differently.** Quiring the sparse matvec (holding inner-product accumulation fixed at quire) gives a real, if smaller and less uniform, solution-accuracy gain — a 1.0x-1.9x band across 6 of 7 matrices tested, with one outlier (sts4098) at ~6.8x.
4. **Iterative refinement recovers the accuracy the inner products lose.** Applying the quire once to a refinement residual (d = b - Ax0, solved once rather than every iteration) — a scheme suggested by Prof. John Gustafson — improves the solution in 228 of 232 seed-runs, with gains ranging ~2x to ~18,500x depending on matrix. This is the paper's central practical result: exact accumulation isn't worth its cost applied in bulk inside the CG loop, but earns its cost applied once at a refinement step.
5. **Formal quire behavior:** quire eliminates accumulation-rounding error entirely (one exact rounding at final readout). It does not eliminate input-quantization error (casting p_i, Ap_i to posit32 before entry), and does not give an nnz-independent error bound. An earlier claim of rel_err <= u was tested and retracted (see docs/quire_error_bound.md).

## Methodology

- CG solver: Jacobi-preconditioned, 300-2000 iterations depending on experiment
- Quire config: quire<N,ES,2> (482-bit for posit32)
- es=2 uniformly, per the 2022 Posit Standard, confirmed directly with Prof. John Gustafson, who identified that early results used a pre-ratified variable-es convention. Correcting to es=2 improved posit16 accuracy by up to 1,043x on some matrices.
- Ground truth: posit64, cross-validated against double64 (agrees to 1e-11 or exactly, across all matrices)
- 13 matrices tested (see docs/results.md for the full properties table); unsymmetric (add32, scircuit, memplus) and artificially-preconditioned (cfd1, cfd2) matrices excluded as CG-invalid, moved to src/exploratory/

**Two separate, independently-run pipelines back this repo's numbers — do not merge their outputs:**

| Pipeline | Location | Settings | Backs |
|---|---|---|---|
| James Quinlan's replication pipeline | external/james_replication/ | MAXITER=300, tol=1e-10 | Tables 3, 4, 6 (pAp/r'z gain, non-transfer, TOST) |
| This repo's own sweep | results/csv/, src/cg_refinement_seeded_v2.cpp | MAXITER=2000, tol=1e-6 | Tables 5, 7, 8 (pre-convergence gain, matvec-quire, refinement) |

## Prof. Quinlan's Three-Part Static/Dynamic Conditioning Extension

- **Part A (static):** posit8/16 fail Cholesky factorization outright under quantization on nearly all matrices; posit32 tracks posit64's condition estimate closely. bcsstk37 is an open anomaly (posit64 itself fails factorization, a structural property, not a precision effect).
- **Part B (dynamic):** p-vector saturation is exactly 0.0 at every iteration, every matrix. Rules out saturation as the divergence mechanism.
- **Part C (correlation):** divergence-onset iteration does not correlate with any static conditioning metric (Spearman rho=0.164, not significant). Divergence depends on CG's dynamic trajectory, not static matrix properties.

## Divergence Mechanism (mhd4800b case study)

Naive posit32 lags float32/quire in full-solver convergence because pAp's magnitude sits outside posit32's precision-favorable zone during early CG iterations, when pAp is largest. This early rounding error compounds through CG's own recurrence. Confirmed causally via controlled hybrid isolation (src/hybrid_probe.cpp), not just correlation. Full derivation in docs/methodology.md.

## Status / Open Work

- cond(x,y) correlation analysis (exploratory, not part of the CoNGA26 submission's main results) is still inconsistent across aggregation methods and should not be read as a settled finding.
- Citations / CITATION.cff: pending final reference list.

## Reproducing Results (Docker)

    git clone https://github.com/Gurleen-kansray/posit-sparse-bench
    cd posit-sparse-bench
    docker build -t posit-bench .
    docker run --rm posit-bench bash run_all.sh

Environment: Ubuntu 22.04, g++ 11, Universal v3.80, quire<N,ES,2>.

## Acknowledgments

Mentors: Kurt Keville (MIT), Joshua Gyllinsky, Prof. John Gustafson (ASU, posit inventor, es=2 standard correction, suggested the iterative-refinement scheme in Finding 4), Theodore Omtzigt (Stillwater Supercomputing), Prof. James Quinlan (University of Maine, three-part conditioning extension, alpha-metric cross-validation, independent TOST replication). Paul Sherman (RISC-V International) for RISC-V Open Lab silicon access and the summation-order framing.

## Full Documentation

- docs/results.md — full 13-matrix result tables (properties, quire gain, solution non-transfer, matvec-quire, refinement)
- docs/methodology.md — CG accumulation-site design, metric definitions, divergence mechanism derivation
- docs/quire_error_bound.md — formal quire error bound analysis and retraction
