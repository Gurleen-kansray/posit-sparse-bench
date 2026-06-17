# SLFFEA beam module — posit32+quire dotX integration

Validates posit32+quire dotX (see ../../README.md for the synthetic precision-ladder
results) inside a real, otherwise-unmodified FEM solver: SLFFEA v1.5 (slffea.com),
beam module, conjugate gradient solver (bmconj.c).

## What changed

bmconj.c calls a new dotX_posit() alongside the existing double-precision dotX()
inside the CG loop and prints both results plus their absolute difference every
iteration; no change to control flow or convergence criteria. fembeam.c has
numnp_LUD_max lowered from 500 to 5 so the small cantilever test model (6 elements,
864 dof) routes through conjugate gradient rather than direct LU decomposition,
since CG is where dotX is exercised. dotx_posit.cpp (new file) implements
dotX_posit() using the Universal library's quire<32,2,2> for exact fused
dot-product accumulation with a single final rounding to posit32.

See bmconj_fembeam.patch for the exact diff against upstream SLFFEA v1.5.

## Result

Cantilever test case, 12 CG iterations to convergence. Full log in
cantilever_dotx_run.log. Differences land in the 1e-9 to 1e-10 absolute range
on values of order 0.03-0.4 — consistent with the posit32+quire relative-error
magnitudes already established on synthetic SuiteSparse matrices in the main
precision ladder experiment. This is the first validation of that result inside
an actual FEM solver codebase rather than a standalone test harness.

## Provenance

Local-only fork for posit experimentation. Not pushed to SLFFEA's upstream
repository — that origin remains San Le's official project.
