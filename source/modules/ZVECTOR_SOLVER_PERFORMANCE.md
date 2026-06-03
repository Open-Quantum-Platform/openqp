# Z-vector / response solver performance notes

Scope: the iterative solvers behind the CPHF/CPKS "z-vector" gradient equations
and the TDHF/SF/MRSF Davidson eigensolvers — `source/pcg.F90`,
`source/modules/tdhf_z_vector.F90`, `source/modules/tdhf_sf_z_vector.F90`,
`source/modules/tdhf_mrsf_z_vector.F90`.

These notes accompany the numerical-stability hardening on this branch and record
where the *time* actually goes, what has been optimized, and what the remaining
high-value (but compiler-gated) opportunities are.

> Build status: `gfortran` is available, but the **full** OpenQP build is not
> reproducible here — Libint is fetched from `qchemlab.knu.ac.kr`, which the
> network policy blocks (`403 host_not_allowed`). So the integrated z-vector
> files (which pull in the Libint-backed ERI engine, DFT grids, tagarray, ...)
> cannot be compiled or run in this environment; only the self-contained solver
> modules (`pcg.F90`, `minres.F90`) can. Those were compiled and benchmarked via
> `tests/solver_microbench/` (see results below). Anything touching the
> integrated files still needs CI / a Fortran-capable machine to validate.

## Measured results (gfortran 13.3, -O2, tests/solver_microbench)

PCG optimization (carried `rz` + scalar-reduction guards), SPD system,
n = 200000, 300 solves:

| build | iters/solve | final residual | wall time |
|---|---|---|---|
| old (recompute rz + full per-iter scans) | 18 | 2.69556e-11 | 19.9 s |
| new (carried rz + scalar-reduction guards) | 18 | 2.69556e-11 | 12.3 s |

→ ~1.6x faster with **identical iteration count and residual** (equivalence
preserved). This is the cheap-matvec upper bound; real z-vector solves are
Fock-build-bound, so the wall-clock effect there is small — the value is the
proven equivalence + preserved fail-closed stability (NaN/Inf -> PCG_BREAKDOWN,
verified by injection).

MINRES validation (`source/minres.F90`):

| test | result |
|---|---|
| SPD system | matches PCG, max\|Δx\| = 4.7e-10 |
| indefinite (Laplacian − 2I) | converges, ‖b−Ax‖/‖b‖ = 1.9e-12 |
| A = diag(+1,−1), b = (1,1) | **MINRES solves in 2 iters; CG breaks down (pᵀAp=0)** |

## 1. Where the cost is

Each iteration of every z-vector solver performs **one full Fock/ERI build** via
`int2_driver%run(...)` (plus the XC response `utddft_fxc` for DFT). That operator
apply is orders of magnitude more expensive than all of the surrounding vector
algebra (dot products, `axpy`, norms, finiteness scans) combined.

Two consequences:

1. For the **z-vector use case**, wall-time ≈ `n_iter × cost(Fock build)`.
   The dominant levers are therefore **(a) reducing `n_iter`** and **(b) never
   doing a redundant Fock build** — not shaving vector bookkeeping.
2. `pcg_mod` is also a **generic** reusable solver (cheap-matvec callers do
   exist), where per-iteration vector overhead *does* matter. The implemented
   change below helps that case directly and is free for the z-vector case.

## 2. Implemented on this branch (low-risk)

`source/pcg.F90`, `pcg_step` / `pcg_init`:

- **Carry the CG numerator `rz = r·M⁻¹r`** as solver state (`pcg_t%rz`), seeded in
  `pcg_init` and refreshed at the end of each `pcg_step`. The previous code
  recomputed `dot_product(r, y)` at the top of every iteration; that value is
  identical to the prior iteration's `rz_new`, so the recompute was a redundant
  O(n) reduction. **Saves one dot product per iteration**, result unchanged.
- **Removed the redundant per-iteration finiteness scans.** The old `pcg_step`
  ran a blanket `ieee_is_finite` scan over `x, b, p, r, y` on entry and then
  rescanned several of the same vectors again later (e.g. `p` twice, `x` twice,
  `y` twice, and `b` — which never changes after `init`). Fail-closed detection
  now rides on the scalar reductions the algorithm already computes: a NaN/Inf
  anywhere in `p`/`Ap`/`r`/`y` propagates into `pap`, `error = ‖r‖`, or `rz_new`,
  each of which is still checked. The full solution vector is scanned **once**,
  at the moment it is about to be returned as converged (so a finite residual
  can't mask a non-finite `x` entry that slipped through a zero in `Ap`).
  **Net: ~1 dot product + several O(n) passes removed per iteration**, with the
  fail-closed breakdown semantics preserved.

These are equivalence-preserving (same iterates, same converged solution); only
redundant work is removed. Regression test:
`tests/test_zvector_solver_stability.py::
test_pcg_step_carries_rz_instead_of_recomputing_each_iteration`.

## 3. "Can GMRES be made as fast as PCG?"

Not in general — they are different algorithm classes:

| | CG / PCG | GMRES(m) |
|---|---|---|
| recurrence | fixed 3-term | full Arnoldi (orthogonalize vs all prior basis vectors) |
| work / iter | O(n) + 1 matvec | O(j·n) + 1 matvec at inner step j |
| memory | ~5 vectors | m+1 stored basis vectors |
| requires | symmetric positive-definite A | general A |

GMRES carries an inherent O(m²·n) orthogonalization + O(m·n) storage penalty per
restart cycle and loses convergence speed at each restart. It cannot match CG's
per-iteration cost. It is only *needed* when A is not SPD.

**Key point for this code:** the z-vector `(A+B)` operator solved here is
symmetric (that is exactly why the default `z_solver`=CG path works and is the
default; `tdhf_mrsf_z_vector.F90`). GMRES is the non-default robustness fallback.
So the right answer to "make GMRES as fast as CG" is usually "don't use GMRES
here — use a symmetric short-recurrence solver," see §5.

## 4. Concrete GMRES redundancy (compiler-gated follow-up)

In `gmres_solve` (`tdhf_mrsf_z_vector.F90`) the operator is applied (a Fock build)
**once more per restart cycle than necessary**:

- The cycle top recomputes `r = b − A·x` via `apply_operator` (around the
  `! Compute initial residual r = b - A*x` block).
- But `x` has not changed since the previous cycle's
  `recompute_gmres_true_residual` (end of cycle) already computed `A·x` and
  `r = b − A·x` for that same `x`. For the **first** cycle the pre-loop initial
  residual (the `! Compute initial residual` block before the `do iter` loop)
  is likewise recomputed.

Reusing the already-computed `r` / `true_residual` at the cycle top removes **one
Fock build per restart** — the single biggest GMRES cost here.

Why it is *not* implemented on this branch: it restructures the restart control
flow that the stability work deliberately added (the "use the true `‖b−Ax‖` at
each restart, not the preconditioned Arnoldi seed norm" guarantee, pinned by
`test_mrsf_gmres_restart_convergence_uses_true_residual_before_preconditioned_beta`).
Done carefully it is numerically identical (deterministic operator, unchanged
`x`), but it must be validated with a compiler and the corresponding stability
tests updated in lockstep.

## 5. Recommended better solvers / levers (priority order)

1. **Reduce iterations — the highest-value lever.** Each iteration is a Fock
   build, so halving `n_iter` ≈ halves the solve.
   - *Better preconditioner.* The current preconditioner is the diagonal
     orbital-energy-difference inverse (`sanitize_*_preconditioner`). A stronger
     approximate-diagonal / level-shifted preconditioner reduces iteration count.
   - *Warm start.* All callers start from `x = 0`. A diagonal warm start
     `x₀ = M⁻¹b` (one preconditioner apply, no extra Fock build) or reusing the
     previous geometry's z-vector along a trajectory typically removes several
     iterations. Converges to the same solution.
2. **MINRES instead of GMRES for symmetric-indefinite A.** If the only reason a
   case needs GMRES is indefiniteness (not non-symmetry), MINRES gives a CG-style
   short 3-term recurrence (constant memory, ~CG cost per iteration) while
   remaining robust where plain CG breaks down. This is the closest thing to
   "GMRES-robustness at CG speed" and is the recommended fallback.
   *Status:* implemented, integrated, and numerically validated.
   `source/minres.F90` is wired into the MRSF z-vector selection through the same
   `apply_z_operator` / `apply_z_precond` used by CG and GMRES (via the
   `minres_apply_op` / `minres_apply_pc` wrappers). Selection is now
   `z_solver` = 0 → CG (default), **1 → MINRES**, 2 → GMRES (`oqpdata.py`).

   Validated with a libint-free local build (`-DUSE_LIBINT=OFF`, native Rys
   integrals — the same configuration CI uses), H2O BHHLYP/6-31G* MRSF-TDDFT
   gradient, state 3, O z-component:

   | z_solver | gradient (O, z) | iters | residual |
   |---|---|---|---|
   | 0 CG (default) | −0.18299119 | — | unchanged vs baseline |
   | 1 MINRES | −0.18299934 | 6 | 6.3e-7 |
   | 2 GMRES | −0.18299937 | 7 | 1.6e-7 |

   MINRES agrees with GMRES to ~3e-8; the default CG path is byte-identical to
   the pre-change baseline. (The ~8e-6 CG vs MINRES/GMRES offset is pre-existing:
   the CG path stops on the *squared* residual `‖r‖²<zvconv`, i.e. a looser
   `‖r‖<√zvconv`, whereas MINRES/GMRES converge to `‖r‖<zvconv`.)
3. **BiCGStab for genuinely non-symmetric A.** Short recurrence, constant memory,
   ~2 matvecs/iter, no growing orthogonalization — far cheaper per iteration than
   restarted GMRES, when a truly non-symmetric operator must be solved.
4. **Trim GMRES overhead when it must be used:** eliminate the redundant Fock
   build (§4), and raise the restart dimension `m` (memory permitting) to cut the
   restart penalty.

Suggested sequencing: land the safe PCG change (done) → add diagonal warm start +
benchmark iteration counts → evaluate MINRES as the symmetric-indefinite fallback
in place of GMRES → only then revisit GMRES internals.
