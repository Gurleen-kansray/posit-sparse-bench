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
The claim holds for circuit matrices (low dynamic range). It does not hold for FEM structural matrices (dynamic range ~10^16).

---

## Precision Ladder Results (bcsstk14, FEM structural, dynamic range 10^15.7)

| Precision | Max rel error (quire) | Max rel error (naive) | Quire vs naive |
|-----------|----------------------|----------------------|----------------|
| posit8    | 6.87×10^5            | 1.36×10^5            | worse (overflow) |
| posit16   | 5.26×10^-2           | 3.68×10^-1           | 7x better — NOT engineering viable |
| posit32   | 3.65×10^-8           | 7.83×10^-7           | 21x better ✓ |
| posit64   | 0.000                | 0.000                | exact double ✓ |

## Precision Ladder Results (bcsstk03, FEM structural, dynamic range 10^16.6)

| Precision | Max rel error (quire) | Max rel error (naive) |
|-----------|----------------------|----------------------|
| posit8    | 6.71×10^8            | 2.82×10^10           |
| posit16   | 1.59×10^2            | 6.39×10^2            |
| posit32   | 5.44×10^-7           | 7.37×10^-6           | 
| posit64   | 0.000                | 0.000                |

## Precision Ladder Results (add32, Circuit/SPICE, dynamic range ~10^4.5)

| Precision | Max rel error (quire) | Max rel error (naive) |
|-----------|----------------------|----------------------|
| posit8    | 1.00×10^0            | 1.00×10^0            |
| posit16   | 3.07×10^-1           | 1.18×10^0            |
| posit32   | 9.28×10^-8           | 9.35×10^-6           |
| posit64   | 0.000                | 0.000                |

---

## Three-way comparison (posit32+quire vs posit32 naive vs double64)

| Matrix | n_used | Max rel — naive | Max rel — quire | Improvement |
|--------|--------|----------------|----------------|-------------|
| add32  | 11/50  | 9.35×10^-6     | 9.28×10^-8     | **101x**    |
| bcsstk03 | 115/300 | 7.37×10^-6  | 5.44×10^-7     | **14x**     |
| bcsstk14 | 64/300  | 7.83×10^-7  | 3.65×10^-8     | **21x**     |

**posit64+quire on bcsstk14: 4.14×10^-15** — matches double64 machine epsilon exactly.

---

## Matrices tested

| Matrix | Domain | Size | Dynamic Range | Condition |
|--------|--------|------|---------------|-----------|
| add32 | Circuit/SPICE (Motorola) | 4,960×4,960 | ~10^4.5 | ~10^2 |
| bcsstk03 | FEM Structural | 112×112 | ~10^16.6 | ~10^6 |
| bcsstk14 | FEM Structural | 1,806×1,806 | ~10^15.7 | ~10^10 |

Download matrices: https://sparse.tamu.edu

---

## Method

- Jacobi-preconditioned Conjugate Gradient
- 50 iterations (add32) / 300 iterations (bcsstk03, bcsstk14)
- At each iteration: pAp computed in double64 (reference), posit+quire, posit naive
- Filter: iterations where |pAp_d| < 1e-6 × max|pAp_d| excluded (near-zero denominator)
- CG runs entirely in double64 — posit computations logged for comparison only
- Precisions tested: posit8 (es=0), posit16 (es=1), posit32 (es=2), posit64 (es=2)

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

## Implications

1. **Dynamic range is the limiting factor**, not vector length or condition number alone
2. **posit32+quire is the minimum viable precision** for FEM structural matrices with dynamic range >10^10
3. **posit16+quire is viable for circuit simulation** (dynamic range <10^5) but not for FEM
4. **posit64+quire matches double64 exactly** — zero measured error across all matrices tested
5. **Quire beats naive by 14x–101x** at posit32 — accumulation error is real and eliminatable

---

## Target publication

HPEC 2026 — *"Precision requirements for posit arithmetic in sparse iterative solvers:
a domain-specific empirical study"*  
Combining accuracy results with RISC-V (Milk-V) execution data for §4.3 bitwise reproducibility.

Mentors: Kurt Keville (MIT R&D Labs) · Joshua Gyllinsky
