# UMRSF Gradient Z-Vector Handoff

Date: 2026-07-01

Private branch target: `codex/umrsf-gradient-zvector-20260701`

Private repository: `git@github.com:karmachoi/openqp-private.git`

Source checkout: `/Users/cheolhochoi/Documents/Code Checking/oqp-uhf-grad-plan`

## Status

This branch preserves the UMRSF gradient/z-vector performance work from the local
`uhf-grad-plan` checkout and merges in the newer upstream analytic M1 overlap
term from `VladimirMakhnev/oqp` commit `6c875b8`.

The Vladimir remote branch was not pushed from this checkout after the user asked
to publish the work to the private repository instead.

## Main Changes

- UMRSF analytic gradient path now uses the coupled alpha/beta response data
  instead of falling back to an expensive numerical response-gradient path.
- UMRSF keeps a response-gradient cache so the closed-form response terms can be
  reused during gradient assembly.
- MRSF-family single-gradient runs can stop Davidson convergence on the requested
  target root for `compute_grad`, making MRSF and UMRSF timing comparisons fairer.
- UMRSF M1 alignment Pulay response now defaults to the analytic
  `umrsf_m1_analytic` path merged from upstream, with `UMRSF_M1FD=1` retained as
  the numerical finite-difference oracle.
- Several hot response contractions were changed to BLAS-backed matrix products,
  and the MRSF integral buffer was guarded against too-small defaults.
- Source-level tests were added for the new target-root and response-gradient
  control logic.

## Validation

- Rebuilt successfully with the existing local OpenQP build tree:
  `cmake --build build --target install --parallel 4`
- Focused regression tests passed:
  `python -m unittest tests.test_davidson_solver_stability tests.test_umrsf_energy_regression tests.test_zvector_solver_stability tests.test_opentrustregion_linalg_config`
- H2O UMRSF gradient sanity check: single-target and full-root runs produced the
  same state-1 energy and gradient at printed precision.
- Thymine SG2, 8-core timing comparison:
  - MRSF target-root gradient, state 3: total about 10.6 s.
  - MRSF target-root gradient, state 1: total about 10.1 s.
  - UMRSF target-root gradient, state 1: total about 23.8 s.

## Notes For Next Developer

- The UMRSF energy stage is now close to MRSF when both are run in equivalent
  single-target-gradient mode. The remaining slowdown is mainly in the UMRSF
  gradient response stage, not in the SCF.
- The analytic M1 path should remain the default. Use `UMRSF_M1FD=1` only as a
  diagnostic oracle because it is intentionally much slower.
- If further speed work is needed, inspect the response/z-vector solve and the
  smooth-basis tracking stage first.
