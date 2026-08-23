# posit-sparse-bench

Benchmarking posit arithmetic (with quire exact accumulation) against IEEE double and float32 in sparse conjugate gradient (CG) solvers, on real symmetric matrices from the SuiteSparse collection. LFX Summer 2026 mentorship project, targeting CoNGA-Q@SC26.

## Research Question

Can posit arithmetic with quire exact accumulation match or exceed double-precision accuracy for the CG inner product p^T A p — the critical dot product in sparse iterative solvers — and does that accuracy carry through to final solution quality?

## Headline Findings

1. **Quire vs naive posit32:** posit32+quire achieves 6.6x-4,531x lower pAp error than naive posit32 across all 13 tested matrices. Exact accumulation, not just wider precision, drives this.
2. **Formal quire behavior:** quire eliminates accumulation-rounding error entirely (one exact rounding at final readout). It does not eliminate input-quantization error (casting p_i, Ap_i to posit32 before entry), and does not give an nnz-independent error bound. An earlier claim of rel_err <= u was tested and retracted (see docs/quire_error_bound.md).
3. **pAp accuracy gain does not transfer to solution accuracy.** Despite pAp gains up to 4,531x, solution-error ratio (naive/quire) clusters near 1.0 across all 13 matrices (median 0.88-1.02) in a 50-seed sweep, independently confirmed by Prof. James Quinlan via a formal TOST equivalence test (pooled ratio 1.02, p=0.003, 6 matrices). This is the paper's central reframing.
4. **Practical convergence win:** on sts4098, posit32+quire converges in 706 iterations vs float32's 800 (real iteration-count speedup, not just per-step accuracy).

## Methodology

- CG solver: Jacobi-preconditioned, 300-2000 iterations depending on experiment
- Quire config: quire<N,ES,2> (482-bit for posit32)
- es=2 uniformly, per the 2022 Posit Standard, confirmed directly with Prof. John Gustafson, who identified that early results used a pre-ratified variable-es convention. Correcting to es=2 improved posit16 accuracy by up to 1,043x on some matrices.
- Ground truth: posit64, cross-validated against double64 (agrees to 1e-11 or exactly, across all matrices)
- 13 matrices tested (see docs/results.md for full properties table); unsymmetric (add32, scircuit, memplus) and artificially-preconditioned (cfd1, cfd2) matrices excluded as CG-invalid, moved to src/exploratory/

## Prof. Quinlan's Three-Part Static/Dynamic Conditioning Extension

- **Part A (static):** posit8/16 fail Cholesky factorization outright under quantization on nearly all matrices; posit32 tracks posit64's condition estimate closely. bcsstk37 is an open anomaly (posit64 itself fails factorization, a structural property, not a precision effect).
- **Part B (dynamic):** p-vector saturation is exactly 0.0 at every iteration, every matrix. Rules out saturation as the divergence mechanism.
- **Part C (correlation):** divergence-onset iteration does not correlate with any static conditioning metric (Spearman rho=0.164, not significant). Divergence depends on CG's dynamic trajectory, not static matrix properties.

## Divergence Mechanism (mhd4800b case study)

Naive posit32 lags float32/quire in full-solver convergence (79 vs 69 iterations) because pAp's magnitude sits outside posit32's precision-favorable zone (~3.16e-5 to ~1e5) during early CG iterations, when pAp is largest. This early rounding error compounds through CG's own recurrence. Confirmed causally via controlled hybrid isolation (src/hybrid_probe.cpp), not just correlation. Full derivation in docs/methodology.md.

## Status / Open Work

- Iterative refinement (residual correction via quire) and permutation/summation-order effects are active, unfinished experiments, not yet in this README pending completion of the current sweep.
- Citations / CITATION.cff: pending final reference list.

## Reproducing Results (Docker)
git clone https://github.com/Gurleen-kansray/posit-sparse-bench
cd posit-sparse-bench
docker build -t posit-bench .
docker run --rm posit-bench bash run_all.sh

Environment: Ubuntu 22.04, g++ 11, Universal v3.80, quire<N,ES,2>.

## Acknowledgments

Mentors: Kurt Keville (MIT), Joshua Gyllinsky, Prof. John Gustafson (ASU, posit inventor, es=2 standard correction), Theodore Omtzigt (Stillwater Supercomputing), Prof. James Quinlan (University of Maine, three-part conditioning extension, alpha-metric cross-validation, independent TOST replication). Paul Sherman (RISC-V Open Lab) for silicon access and the summation-order framing.

## Full Documentation

- docs/results.md — full 13-matrix result tables
- docs/methodology.md — ladder design, metric definitions, divergence mechanism derivation
- docs/quire_error_bound.md — formal quire error bound analysis and retraction
