# posit-sparse-bench

**Precision ladder experiment: posit8/16/32/64+quire vs naive on real SuiteSparse sparse matrices**

LFX Mentorship — RISC-V High Precision | IEEE HPEC 2025 Theme  
Gurleen Kaur | Mentors: Kurt Keville (MIT), Joshua Gyllinsky

---

## Key Finding

**Quire exactness is necessary but not sufficient — posit precision must match matrix dynamic range.**

posit16+quire fails on high-dynamic-range FEM matrices (bcsstk03: rel error 159, bcsstk14: rel error 0.05).  
posit32+quire succeeds on all matrices tested.  
posit64+quire = exact double (zero measured error) on all three matrices.

This provides empirical boundary conditions on Gustafson's HPEC 2025 claim that  
*"16-bit posits can replace 64-bit floats in stable computations via exact dot products."*  
The claim does not hold for FEM structural matrices with dynamic range ~10^16.

---

## Precision Ladder Results (bcsstk14, FEM structural)

| Precision | Max rel error (quire) | Max rel error (naive) | Verdict |
|-----------|----------------------|----------------------|---------|
| posit8    | 6.87×10^5            | 1.36×10^5            | overflow — unusable |
| posit16   | 5.26×10^-2           | 3.68×10^-1           | 7x quire wins — NOT engineering viable |
| posit32   | 3.65×10^-8           | 7.83×10^-7           | 21x quire wins ✓ viable |
| posit64   | 0.000                | 0.000                | exact double ✓ |

## Precision Ladder Results (bcsstk03, FEM structural)

| Precision | Max rel error (quire) | Max rel error (naive) | Verdict |
|-----------|----------------------|----------------------|---------|
| posit8    | 6.71×10^8            | 2.82×10^10           | overflow — unusable |
| posit16   | 1.59×10^2            | 6.39×10^2            | NOT viable |
| posit32   | 5.44×10^-7           | 7.37×10^-6           | 14x quire wins ✓ viable |
| posit64   | 0.000                | 0.000                | exact double ✓ |

## Precision Ladder Results (add32, Circuit/SPICE)

| Precision | Max rel error (quire) | Max rel error (naive) | Verdict |
|-----------|----------------------|----------------------|---------|
| posit8    | 1.00×10^0            | 1.00×10^0            | overflow — unusable |
| posit16   | 3.07×10^-1           | 1.18×10^0            | 4x quire wins — NOT engineering viable |
| posit32   | 9.28×10^-8           | 9.35×10^-6           | 101x quire wins ✓ viable |
| posit64   | 0.000                | 0.000                | exact double ✓ |

---

## Three-way comparison (posit32+quire vs posit32 naive vs double64)

| Matrix | n_used | Max rel — naive | Max rel — quire | Improvement |
|--------|--------|----------------|----------------|-------------|
| add32  | 11/50  | 9.35×10^-6     | 9.28×10^-8     | **101x**    |
| bcsstk03 | 115/300 | 7.37×10^-6  | 5.44×10^-7     | **14x**     |
| bcsstk14 | 64/300  | 7.83×10^-7  | 3.65×10^-8     | **21x**     |

**posit64+quire on bcsstk14: 4.14×10^-15** — matches double64 machine epsilon exactly.

---

## Matrices tested (all verified from .mtx files)

| Matrix | Domain | Size | Dynamic range | Diag ratio (proxy for conditioning) |
|--------|--------|------|---------------|--------------------------------------|
| add32 | Circuit/SPICE (Motorola) | 4,960×4,960 | ~10^4.5 | ~10^2 (condition ~137) |
| bcsstk03 | FEM Structural | 112×112 | ~10^16.6 | ~4.55×10^8 |
| bcsstk14 | FEM Structural | 1,806×1,806 | ~10^15.7 | ~8.94×10^9 |

Download matrices: https://sparse.tamu.edu

---

## Method

- Jacobi-preconditioned Conjugate Gradient
- 50 iterations (add32) / 300 iterations (bcsstk03, bcsstk14)
- At each iteration: pAp computed in double64 (reference), posit+quire, posit naive
- Filter: iterations where |pAp_d| < 1e-6 × max|pAp_d| excluded (near-zero denominator near convergence)
- CG runs entirely in double64 — posit computations logged for comparison only
- Precisions: posit8 (es=0), posit16 (es=1), posit32 (es=2), posit64 (es=2)

---

## Build

Requires [Universal Numbers Library](https://github.com/stillwater-sc/universal) v3.80+

```bash
make INCLUDES=-I/path/to/universal/include
make run_all
```

---

## Repository structure
src/           — C++ harness (three-way comparison + precision ladder)

results/logs/  — raw per-iteration logs

results/csv/   — per-iteration relative error trends

data/matrices/ — .mtx files (Matrix Market, from SuiteSparse)

docs/          — presentation slides
---

## Implications (all supported by log data)

1. **Dynamic range is the limiting factor** — posit16+quire fails on matrices with dynamic range >10^15, regardless of quire exactness
2. **posit32+quire is the minimum viable precision** for FEM structural matrices tested here
3. **posit16+quire is not engineering viable on any matrix tested** — even circuit matrices show 0.31 max relative error at posit16
4. **posit64+quire matches double64 exactly** — zero measured error across all three matrices
5. **Quire beats naive by 14x–101x at posit32** — accumulation error is real and eliminatable by quire

---

## Target publication

HPEC 2026 — *"Precision requirements for posit arithmetic in sparse iterative solvers: a domain-specific empirical study"*  
Combining accuracy results with RISC-V (Milk-V) execution data for §4.3 bitwise reproducibility.

Mentors: Kurt Keville (MIT R&D Labs) · Joshua Gyllinsky

---

## Real-codebase validation: SLFFEA

The precision-ladder results above use standalone benchmark harnesses on extracted
SuiteSparse matrices. As a follow-up, posit32+quire dotX was integrated into the
conjugate gradient solver of SLFFEA v1.5 (beam module) — a real, otherwise-unmodified
open-source FEM package, not a synthetic test case. On a cantilever test model
(6 elements, 864 dof), it converges in 12 CG iterations with posit32+quire vs double
absolute differences in the 1e-9 to 1e-10 range, consistent with the precision-ladder
findings above. Full integration details, patch, and logs in external/slffea-beam/.

---

## Performance overhead: software-emulated posit vs native double

Measured on x86 (WSL2, -O2) using time-bounded microbenchmark; same vector sizes
as precision ladder. Each function timed for 1 second per case.

| matrix (n)       | double (ns) | posit32+quire (ns) | quire/double |
|------------------|-------------|--------------------|--------------|
| bcsstk03 (112)   | 80.9        | 577,264            | 7,132x       |
| bcsstk14 (1806)  | 1,042.6     | 8,244,849          | 7,908x       |
| add32 (4960)     | 2,473.1     | 26,422,731         | 10,684x      |

quire/naive ratio ≈ 1.0 across all cases: the overhead is posit conversion
(software emulation), not quire accumulation specifically. On hardware with a
native posit FPU (the motivation for RISC-V posit work), this overhead disappears.

---

## Convergence boundary: where quire matters for solver correctness

Beyond accuracy differences, quire provides **solver-level correctness** where naive
posit32 fails entirely. Experiment: scale bcsstk03's largest diagonal entries by
factor S to push dynamic range from baseline ~1.5e6 up through 1.5e18, then run
preconditioned CG with all three methods (double, posit32+quire, posit32+naive).

| scale | diag dynamic range | double (iters) | quire (iters) | naive (iters) |
|-------|--------------------|----------------|---------------|---------------|
| 1e+00 | 1.52e+06           | 173            | 263           | 271           |
| 1e+08 | 1.52e+14           | 148            | 265           | 343           |
| 1e+10 | **1.52e+16**       | 156            | 292           | **DIV**       |
| 1e+12 | **1.52e+18**       | 177            | 430           | **DIV**       |

**Finding:** naive posit32 diverges at dynamic range ~1.52e+16. Posit32+quire
converges at the same condition. The quire's exact dot product accumulation
preserves the CG descent direction where naive rounding errors corrupt it
beyond recovery. This is the first empirical demonstration of this boundary
on a real FEM matrix.

Full log and source in results/logs/bcsstk03_convergence_boundary.log and
src/bcsstk03_convergence_boundary.cpp.

---

## Full posit32 solver: does quire robustness survive posit32 matvec?

Previous convergence boundary experiment used double-precision matvec with posit32
dot products only. This experiment runs the full CG solver in posit32 — both matvec
and dot product — to model what a real RISC-V posit deployment would look like.

| scale | diag_range | dbl/dbl | dbl-mv/quire | dbl-mv/naive | p32-mv/quire | p32-mv/naive |
|-------|-----------|---------|-------------|-------------|-------------|-------------|
| 1e+08 | 1.52e+14  | 148     | 265         | 343         | 299         | **DIV**     |
| 1e+09 | 1.52e+15  | 153     | 257         | 280         | 330         | 313         |
| 1e+10 | 1.52e+16  | 156     | 292         | **DIV**     | 367         | **DIV**     |
| 1e+12 | 1.52e+18  | 177     | 430         | **DIV**     | 439         | **DIV**     |

**Findings:**
1. Posit32 matvec pulls the naive divergence boundary down by ~one decade (1.52e+16 → 1.52e+14)
2. Quire survives full posit32 arithmetic all the way to 1.52e+18 — same ceiling as hybrid mode
3. Quire's exact dot product compensates for posit32 matvec rounding, not just naive dot product errors
4. Anomaly at scale=1e+09: p32-mv/naive converges (313 iters) between two diverging scale points — noted for further investigation

For a real RISC-V posit32 deployment, quire is not just more accurate — it is what
keeps the solver converging at all under extreme dynamic range conditions.

### Correction to scale=1e+09 result

The p32-mv/naive "convergence" at scale=1e+09 (313 iterations) is a **false positive**.
Per-iteration residual tracing (src/bcsstk03_anomaly_diag.cpp) shows all three
scales around the anomaly (1e+08, 1e+09, 1e+10) produce identical chaotic residual
behavior after ~iteration 150: erratic oscillation between ~1e-2 and ~1e-6 with no
monotone descent. The scale=1e+09 run happens to cross the 1e-6 threshold once at
iter 312 during this wandering — not genuine convergence.

**Corrected interpretation:** naive posit32 enters a chaotic residual regime at
dynamic range ~1.52e+14. In this regime, convergence checks are unreliable —
a solver could report "converged" while solution quality is actually poor.
Quire avoids this regime entirely by maintaining exact accumulation.

---

## add32: posit32+quire matches double exactly; posit16 fails on value range

add32 is a circuit simulation matrix (n=4960) with diagonal dynamic range of only
5.6 — superficially ideal for posit16. However, full matrix value range spans
1.976e-38 to 4.232e-02, a ratio of ~2.1e+36, far exceeding posit16's representable
range (~9 decades). Result using actual add32_b.mtx right-hand side:

| method        | iters | converged | final rel residual |
|---------------|-------|-----------|--------------------|
| double        | 68    | YES       | 5.245e-09          |
| posit32+quire | 68    | YES       | 5.245e-09          |
| posit16+quire | 303   | NO        | 7.003e+79          |
| posit16+naive | 3     | NO        | 6.650e+32          |

**posit32+quire matches double exactly** — same iteration count, same residual.
posit16 fails regardless of accumulation method because matrix values below
posit16's minimum representable number (~1e-9) are flushed to zero.

**Methodological finding:** diagonal dynamic range is an insufficient metric
for posit precision selection. Full matrix value range must be considered.
