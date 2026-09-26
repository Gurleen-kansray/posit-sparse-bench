# Methodology

## Three accumulation sites

Each CG iteration has three reduction sites whose accumulation policy (naive posit32 vs quire posit32) is varied independently:

1. **p^T A p** — the curvature term, denominator of alpha_k. Recomputed fresh each iteration; errors here do not directly carry forward.
2. **r^T z** — the preconditioned residual inner product, numerator of alpha_k and feeds beta_k, which scales the previous search direction. Errors here compound forward through every later search direction.
3. **The sparse matrix-vector product A p_k** — each row is itself a dot product.

Naive posit32 and quire posit32 use the identical posit32 format and identical unit round-off (u = 2^-27 ~= 7.45e-9); only the accumulation policy at a given site differs. This isolates the quire's contribution from any general precision difference between number formats (the confound identified in Buoncristiani et al. 2020).

## Divergence mechanism (mhd4800b case study)

Naive posit32 lags float32/quire in full-solver convergence because pAp's magnitude sits outside posit32's precision-favorable zone during early CG iterations, when pAp is largest. This early rounding error compounds through CG's own recurrence rather than being corrected by later iterations. Confirmed via controlled hybrid isolation (src/hybrid_probe.cpp): forcing quire accumulation only during the early high-magnitude iterations reproduces the accuracy gain seen from quiring throughout, and forcing it only later does not.

## Metric definitions

- **Quire gain** (Table 3): ratio of maximum relative error under naive posit32 accumulation to that under quire accumulation, against a double-precision reference, computed per seed then averaged.
- **Solution-error ratio** (Table 4): relative solution error of posit32+quire divided by that of posit32 naive, both measured against the known ground-truth solution vector x*.
- **Pre-convergence gain** (Table 5): same as quire gain, restricted to iterations before the double-precision reference reaches relative residual 1e-6.
- **TOST equivalence** (Table 6): formal two-one-sided-tests procedure testing whether posit32 quire and naive solution errors are statistically equivalent within a stated 10% margin.
- **Refinement improvement factor** (Table 8): relative solution error of the unrefined posit32 solve divided by that of the quire-refined solve.

## Pipelines — do not merge

See the table in the top-level README. James Quinlan's replication pipeline (external/james_replication/, MAXITER=300/tol=1e-10) and this repo's own sweep (results/csv/, MAXITER=2000/tol=1e-6) use different codebases, seeds, and convergence criteria, and back different tables. Treat any cross-pipeline number comparison as informal corroboration, not a merged dataset.

*(Working skeleton — expand with full ladder/probe design and ILU(0)/Jacobi comparison before camera-ready.)*
