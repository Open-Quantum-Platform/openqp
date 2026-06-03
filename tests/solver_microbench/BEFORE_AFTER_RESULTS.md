# Before/after benchmark: `origin/main` vs `claude/Davidson-ZVector`

Both built libint-free (`-DUSE_LIBINT=OFF`, native Rys integrals — the CI config)
with gfortran 13.3, OpenBLAS64 ILP64. `origin/main` is the pre-everything
baseline (before the Davidson/z-vector stability hardening, the PCG
optimization, the MINRES solver, and the readability refactor — 43 commits).

Test cases (H2O):
- MRSF-TDDFT gradient, BHHLYP/6-31G* (`examples/MRSF-TDDFT/...GRADIENT.inp`)
- TDDFT (RPA) energy, B3LYP5/... (`examples/TDDFT/H2O_B3LYP5-TDDFT_ENERGY.inp`)

## Accuracy — bit-identical

| quantity | origin/main | HEAD |
|---|---|---|
| MRSF Davidson iterations | 3 | 3 |
| MRSF excitation energies (eV) | −7.714810 / 1.232686 / 2.818148 | identical |
| TDDFT Davidson iterations | 4 | 4 |
| TDDFT excitation energies (eV) | 8.114387 / 10.192553 / 10.556486 | identical |
| MRSF z-vector CG gradient (O,z) | −0.18299119 | identical |

The Davidson stability guards (reorthogonalize-on-drift, denominator-floor
preconditioner, non-finite-residual guards) do **not** change iteration count or
energies on well-behaved problems — they only activate on NaN/Inf/near-singular
inputs. Same for the z-vector hardening: default CG output is unchanged.

## Performance

| metric | result |
|---|---|
| Full MRSF gradient wall time (min of 7) | main 0.962 s -> HEAD 0.972 s (~1%, within noise) |
| PCG solver micro-opt, isolated (`run.sh`) | OLD 26.05 s -> NEW 15.70 s = **1.66x**, identical iters/residual |
| z-vector iterations by solver | CG 4 / MINRES 6 / GMRES 7 (MINRES & GMRES converge ~10x tighter than CG's squared-norm test) |

Conclusion: the branch adds robustness (NaN/Inf fail-closed across z-vector and
Davidson), a faster preconditioned-CG inner loop, and a MINRES option for
indefinite operators, with **no measurable cost** to accuracy or wall time on
normal runs.

Reproduce the isolated solver numbers with `./run.sh` in this directory; the
full-calc numbers need a built OpenQP (see the repo build instructions).
