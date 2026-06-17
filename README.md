# posit-sparse-bench

**Posit32+quire vs posit32 naive vs double64 on real SuiteSparse sparse matrices**

LFX Mentorship — RISC-V High Precision | IEEE HPEC 2025 Theme  
Mentors: Kurt Keville (MIT), Joshua Gyllinsky

---

## What this is

We instrument the pAp dot product inside a Conjugate Gradient solver three ways:
- **double64** — reference
- **posit32+quire** — exact accumulation, single final rounding (§4.2, 2022 Posit Standard)
- **posit32 naive** — round after every multiply-add (standard float behavior)

Every iteration is logged. Results show quire beats naive by **14x–101x** on relative error.

---

## Matrices tested

| Matrix | Domain | Size | Dynamic Range | Source |
|--------|--------|------|---------------|--------|
| add32 | Circuit/SPICE | 4,960×4,960 | ~10^4.5 | SuiteSparse / Hamm |
| bcsstk03 | FEM Structural | 112×112 | ~10^16.6 | SuiteSparse / HB |
| bcsstk14 | FEM Structural | 1,806×1,806 | ~10^15.7 | SuiteSparse / HB |

Download matrices from: https://sparse.tamu.edu

---

## Key results (Jacobi-preconditioned CG, dynamic threshold filter)

| Matrix | n_used | Max rel — naive | Max rel — quire | Improvement |
|--------|--------|----------------|----------------|-------------|
| add32 | 11/50 | 9.35e-6 | 9.28e-8 | **101x** |
| bcsstk03 | 115/300 | 7.37e-6 | 5.44e-7 | **14x** |
| bcsstk14 | 64/300 | 7.83e-7 | 3.65e-8 | **21x** |

**posit64+quire on bcsstk14: 4.14×10^-15** — matches double64 machine epsilon exactly.

---

## Method

- Jacobi preconditioner (diagonal scaling)
- 50 iterations (add32) / 300 iterations (bcsstk03, bcsstk14)
- Filter: iterations where |pAp_d| < 1e-6 × max|pAp_d| excluded
  (denominator too close to zero near convergence — standard CG behavior)
- CG runs entirely in double64; posit computations logged for comparison only

---

## Build

Requires [Universal Numbers Library](https://github.com/stillwater-sc/universal) (v3.80+).

```bash
# Point UNIVERSAL_INC to your universal/include path
make INCLUDES=-I/path/to/universal/include
make run_all
```

---

## Repository structure
src/          — C++ harness for each matrix

results/logs/ — raw per-iteration logs (iter, pAp_d, pAp_p, pAp_pn, diffs)

results/csv/  — per-iteration relative error trends (naive vs quire)

data/matrices/— .mtx files (Matrix Market format, from SuiteSparse)

docs/         — presentation slides, speech notes
---

## Next steps

- RISC-V execution on Milk-V (§4.3 bitwise reproducibility check)
- Quire vs Ozaki comparison on same matrices
- posit64+quire on bcsstk03
- Larger circuit matrix (G2_circuit, 150K×150K)

---

## Citation

If you use this work, please cite:  
Gurleen Kaur, Kurt Keville, Joshua Gyllinsky.  
*"posit32+quire for sparse iterative solvers: FEM and circuit simulation."*  
LFX Mentorship / IEEE HPEC 2025 Theme, June 2026.
