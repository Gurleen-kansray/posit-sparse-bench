# Lemma: Quire-Accumulated Inner Product Error Bound

## Setup

Let A be an n×n sparse symmetric matrix with nnz nonzero entries, and let
p, Ap ∈ R^n be posit32-representable vectors. Define the inner product

    s = sum_{i=1}^{nnz} p_i * (Ap)_i

We compare two accumulation strategies:

1. **Naive posit32 accumulation**: each partial sum is rounded to posit32
   after every addition.
2. **Quire accumulation**: all nnz products are accumulated exactly in a
   fixed-point register (the quire) with no intermediate rounding; a single
   rounding step occurs only when the final sum is read out to posit32.

## Quire capacity

For posit32 with es=2, the standard quire configuration used here is
quire<32,2,2> (482 bits, capacity parameter c=2). The quire's fixed-point
range is sized to exactly represent the sum of up to 2^c · (maxpos)^2
worst-case-magnitude terms without overflow, i.e.

    sum |p_i * (Ap)_i| < 2^c * maxpos^2   =>   exact accumulation guaranteed

## Lemma

If no overflow occurs (satisfied whenever the above magnitude bound holds,
which we confirm empirically holds with 23–33 orders of magnitude of margin
for all 13 tested matrices), then

    rel_err(s_quire) = |s_quire - s_exact| / |s_exact| <= u

where u = 2^-(nbits - es - 2) is the posit32 unit roundoff (u ≈ 2^-28 for
es=2), **independent of nnz**.

## Proof sketch

Each product p_i * (Ap)_i is formed exactly in the quire's fixed-point
representation (no rounding at multiply time, since the quire has enough
bits to hold the full product exactly — this is the defining property of
capacity ≥ 1). Accumulation of these exact terms in fixed-point arithmetic
introduces no rounding error at any step, since fixed-point addition is
exact as long as the accumulator doesn't overflow (guaranteed by the
capacity bound above). Therefore the quire holds s_exact bit-for-bit
before readout.

The only rounding operation in the entire computation is the final
quire-to-posit32 conversion, which is a single correctly-rounded cast.
By definition of posit32's unit roundoff u, any single correctly-rounded
cast satisfies relative error ≤ u. Hence:

    rel_err(s_quire) <= u,  for all nnz.

This is in contrast to naive sequential posit32 (or float32) accumulation,
where each of the nnz-1 intermediate additions is independently rounded.
Under the classical floating-point error accumulation model (Wilkinson),
worst-case relative error for sequential summation grows as:

    rel_err(s_naive) <= (nnz - 1) * u + O(u^2)

i.e. naive accumulation error scales linearly with term count, while
quire's does not scale with nnz at all — it is fixed at a single
rounding-unit bound regardless of how many terms are summed.

## Empirical confirmation

Across all 13 SuiteSparse test matrices (nnz range: 640 to 4,824,880 — a
7,500x spread), observed rel_err(pAp_quire) stayed within the single-cast
bound u in every iteration, with no upward trend as nnz increases. This is
the direct empirical signature predicted by the lemma: naive error grows
with matrix size/density (consistent with quire gains scaling up to 4,531x
on the largest, densest matrices — s3dkq4m2, s3dkt3m2), while quire error
does not.

## Open items / caveats for review

- This bound assumes no overflow; the 23–33 order-of-magnitude margin
  observed empirically should be stated as an empirical safety margin, not
  a universal guarantee — pathological matrices with extreme value ranges
  could in principle approach the capacity bound.
- The O(u^2) term in the naive bound is dropped for tightness; should
  confirm this is negligible at posit32 precision for our tested magnitude
  ranges (likely fine, but worth a sentence in the paper acknowledging it).
- This is a *deterministic worst-case* bound (Wilkinson-style), not a
  probabilistic/average-case one. Given Gustafson's background, may be
  worth noting whether a probabilistic bound (à la Higham's stochastic
  rounding analysis) would give a tighter expected-case comparison — could
  preempt a question from him about worst-case vs typical-case framing.
