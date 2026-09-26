# posit-sparse-bench

Code and data for **"Exact Accumulation Improves the Inner Product but Not the Attainable Accuracy of Conjugate Gradient"** (Kaur, Omtzigt, Quinlan, Sherman — CoNGA26). LFX Summer 2026 mentorship project.

## What this paper shows

The posit standard's quire lets a dot product be summed exactly and rounded once. Conjugate gradient (CG) computes two inner products and a sparse matrix-vector product per iteration — natural candidates for exact accumulation. Across 13 SuiteSparse symmetric positive definite matrices, we find:

1. **The quire delivers its designed gain at the reduction level.** Median 260x lower error on the curvature term p^T A p and 423x on the residual inner product r^T z, against a double-precision reference (50-seed sweep).
2. **That gain does not transfer to the final solution.** The solution-error ratio (posit32+quire vs. posit32 naive) clusters near 1.0 across all 13 matrices (median 0.92-1.05), confirmed by a formal TOST equivalence test on an independently replicated 6-matrix subset (pooled ratio 1.025, p=0.003, 300 paired runs).
3. **The matrix-vector product is the exception.** Quiring the sparse matvec (inner products held fixed at quire) gives a real solution-accuracy gain — a 1.0x-1.9x band across 6 of 7 matrices, with sts4098 an outlier at ~6.8x.
4. **Iterative refinement is where the quire earns its cost.** Applying it once to a refinement residual (d = b - Ax0, solved once rather than every iteration — a scheme suggested by Prof. John Gustafson) improves the solution in 228 of 232 seed-runs, with gains from ~2x to ~18,500x depending on matrix.

Section 6 of the paper connects this to Greenbaum's finite-precision CG theory: rounding error acts as a spectral perturbation on the problem CG actually solves, which offers a candidate mechanism for why per-iteration accuracy doesn't reliably reach the final solution.

## Repository Structure

`docs/results.md` contains the full per-matrix result tables matching the paper's Tables 1-8, with matrix properties, quire gain figures, the solution non-transfer result, the matvec-quire ablation, and the refinement improvement factors.

Two pipelines produced these results: an independent replication pipeline built by co-author James Quinlan (`external/james_replication/`, MAXITER=300, tol=1e-10), and this repo's own solver sweep (`results/csv/`, `src/cg_refinement_seeded_v2.cpp`, MAXITER=2000, tol=1e-6). Both are documented in `docs/methodology.md`, and their outputs should not be merged, as they use different seeds and convergence criteria.

## Methodology

- CG solver: Jacobi-preconditioned, 300-2000 iterations depending on experiment
- Quire config: quire<N,ES,2> (482-bit for posit32), es=2 uniformly per the 2022 Posit Standard
- Ground truth: posit64, cross-validated against double64
- 13 matrices tested; unsymmetric and artificially-preconditioned matrices excluded as CG-invalid

Full metric definitions and the three-accumulation-site design are in `docs/methodology.md`.

## Reproducing results (Docker)

    git clone https://github.com/Gurleen-kansray/posit-sparse-bench
    cd posit-sparse-bench
    docker build -t posit-bench .
    docker run --rm posit-bench bash run_all.sh

Environment: Ubuntu 22.04, g++ 11, Universal v3.80, quire<N,ES,2>.

## Open items

- Table 7's matvec-quire script (James Quinlan's) not yet added to this repo.
- CITATION.cff pending final reference list.

## Acknowledgments

Prof. John Gustafson (ASU, posit inventor) for the es=2 standard correction and for suggesting the iterative-refinement scheme (Finding 4). Prof. James Quinlan (University of Maine) for the independent TOST replication and alpha-metric cross-validation. Theodore Omtzigt (Stillwater Supercomputing). Kurt Keville (MIT) and Joshua Gyllinsky, mentors. Paul Sherman (RISC-V International) for RISC-V Open Lab access.

## Full documentation

- docs/results.md — full 13-matrix result tables matching the paper's Tables 1-8
- docs/methodology.md — accumulation-site design, metric definitions
- docs/quire_error_bound.md — formal quire error bound analysis and retraction
