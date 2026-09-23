# UMRSF Gradient Performance Handoff

Date: 2026-07-02

Private branch target: `codex/umrsf-gradient-zvector-20260701`

Private repository: `git@github.com:karmachoi/openqp-private.git`

Source checkout: `/Users/cheolhochoi/Documents/Code Checking/oqp-uhf-grad-plan`

Commit pushed: `21285da` (on top of `80b75b3`, the previous
`2026-07-01-umrsf-gradient-zvector` handoff state).

## Status

The UMRSF analytic gradient was ~4x slower than MRSF in the **gradient phase**.
This session brought it to **2.77x** (gradient phase) / **1.86x** (total run)
through surgical, byte-identical changes to `source/modules/tdhf_umrsf_gradient.F90`.
Every change was verified to reproduce the previous gradient to machine
precision (thymine 6-31G* + butadiene, each compared against its own dense /
Davidson oracle: `max|new - oracle| = 0` to printed precision).

Only `tdhf_umrsf_gradient.F90` was committed. The pre-existing (not-ours)
working-tree edit to `pyoqp/oqp/utils/input_checker.py` was left uncommitted;
it is what fails 12 input-validation tests (with it reverted to HEAD, all 48
focused tests pass).

## CRITICAL build note (this explained a 2.5x phantom slowdown)

The wall-clock numbers depend entirely on the BLAS the library links:

- **Correct (fast): Apple Accelerate, LP64.** Build the venv with
  `-C cmake.define.LINALG_LIB_INT64=OFF -C cmake.define.BLA_VENDOR=Apple`.
  `otool -L .../oqp/lib/liboqp.dylib` must show `Accelerate.framework`.
- **Wrong (slow): scikit-build default is `LINALG_LIB_INT64=ON` (ILP64)**,
  which on macOS cannot use Accelerate's LP64 interface and silently falls back
  to the **bundled netlib reference BLAS** (statically linked, no external BLAS
  in `otool -L`). That runs ~2.5x slower across every phase (SCF, XC, Davidson,
  gradient) and was the entire cause of an apparent "UMRSF is 25.8s vs codex 11s"
  discrepancy — same source, same machine, only the BLAS differed. Verify with
  `/usr/bin/time` user-seconds (netlib ~189 vs Accelerate ~74 for the MRSF run).

Full working install command:
```
CC=gcc-15 CXX=g++-15 FC=gfortran-15 .venv/bin/python -m pip install \
  --no-build-isolation --force-reinstall --no-deps \
  -C cmake.define.LINALG_LIB_INT64=OFF -C cmake.define.BLA_VENDOR=Apple .
```
The manual `.../build` tree (used via `cmake --build build --target install`)
is already configured Accelerate LP64. Inputs must set `[scf] save_molden=False`
and `[guess] save_mol=false` (an LP64 molden-writer path crashes otherwise;
irrelevant to the gradient).

## Main Changes (all in `tdhf_umrsf_gradient.F90`)

1. **Alignment adjoint solve.** The default `de_m1` / M1 path
   (`umrsf_m1_analytic`) solved `H^T mu = lambda` with a full dense `dgelss`
   **SVD over the whole npair system** — for thymine (npair ~ 6969) this did not
   finish in 2 minutes even on Accelerate. Now: the per-segment block solver
   `umrsf_solve_alignment_adjoint_blocks` (H is exactly block-diagonal), and for
   large blocks a **matrix-free CG -> MINRES -> dense** solve. H is symmetric
   (`max|H-H^T| ~ 1e-14`, verified) and diagonally dominant after Jacobi
   preconditioning, so it converges in ~4 iterations. Env: `UMRSF_ALIGN_CG=0`
   forces MINRES, `UMRSF_ALIGN_DENSE=1` forces the dense `dgesv`/`dgelss` oracle,
   `UMRSF_ALIGN_SVD=1` forces SVD.

2. **de_m1 / dG^f dedup.** `umrsf_genfock_analytic` (the aligned-to-canonical
   dG^f) and `umrsf_m1_analytic` (de_m1) were building the *identical* `hmat`,
   solving the *identical* system, and computing the *identical* `T-bar`.
   `umrsf_genfock_analytic` now returns `T-bar` (optional `tbar_out` arg) and the
   default M1 path uses the new `umrsf_m1_from_tbar` (just the cheap
   `W_m1 = C_a T-bar C_b^T` -> `-Tr(W_m1 S^x)` tail). The oracle branches
   (`UMRSF_GFFD`, `UMRSF_M1FD`) still use the standalone `umrsf_m1_analytic`.

3. **Smooth-basis amplitude-tracking Davidson removed by default.** The gradient
   used to re-diagonalize A in the smooth-aligned basis (`umrsf_track_amplitude_dav`,
   ~6 response matvecs, ~4s) to get a "genuine eigenvector". But per RULES §17 the
   energy step uses the SAME aligner, so the stored `bvec` is already the
   smooth-basis eigenvector (overlap deficit -> 0). Default now uses it directly:
   `xamp = bvec(:,tstate)/||bvec||`, `omega_eig = td_en(tstate)`, no matvecs.
   Re-diagonalization kept as an opt-in oracle: `UMRSF_TRACK=dav|dense|cmp`
   (legacy `UMRSF_TRKDENSE=1`/`UMRSF_TRKCMP=1` map to `dense`/`cmp`).
   MRSF never does this (single ROHF orbital set).

4. **P_eff parallelized.** `umrsf_build_peff` probed `umrsf_orb_matvec` with n^2
   unit matrices (O(n^4), serial). The probes are independent (each writes a
   unique element; `umrsf_orb_matvec` uses only local scratch), so the loops are
   now OMP-parallel with per-thread private scratch. 0.99s -> 0.14s, byte-identical.

5. **Hot `genfock_z` transforms -> BLAS.** The per-z-vector-matvec congruence
   transforms `C z C^T` and `C^T Y C` were `matmul` intrinsics; now explicit
   `dgemm`. (Consistency with MRSF style; no measurable speed change — see below.)

## Validation

- Rebuilt (Accelerate LP64) cleanly.
- Byte-identical gradients: thymine 6-31G*/BHHLYP/SG2 (Max 0.1502110 /
  0.0941788 depending on geometry) and butadiene 6-31G* (Max 0.0633842) each
  reproduce their `UMRSF_ALIGN_DENSE=1` / `UMRSF_TRACK=dav` oracle to
  `max|diff| = 0`.
- Focused tests pass: `tests.test_davidson_solver_stability`,
  `tests.test_zvector_solver_stability` (35 tests OK). The 12 failures in the
  full focused set are from the uncommitted `input_checker.py` change, not ours.
- Timing (thymine 6-31G*/BHHLYP/SG2, nbf=147, Accelerate LP64, 8 cores):

  | Phase    | MRSF  | UMRSF | ratio |
  |----------|-------|-------|-------|
  | SCF      |  3.1s |  2.8s | 0.91x |
  | Energy   |  2.8s |  4.4s | 1.57x |
  | Gradient |  4.1s | 11.3s | 2.77x |
  | Total    | 10.0s | 18.6s | 1.86x |

  (UMRSF gradient phase was 16.3s at the start of this session -> 11.3s now.)

## Notes For Next Developer

- **matmul -> dgemm is NOT a speed lever here.** Three independent tests
  confirmed it: `-fexternal-blas` (routes every matmul to dgemm) gave 0 change;
  the hot `genfock_z` conversion gave 0 change. At nbf~147 the dense products are
  small and the genuinely hot response transforms (`umrsfmntoia`/`umrsfcbc` in
  `tdhf_mrsf_lib`) already use `dgemm`. The remaining ~92 `matmul` calls in the
  UMRSF gradient are one-shot O(n^3) setup/diagnostic code — converting them is a
  consistency cleanup with zero speedup and real risk; do it separately if at all.

- **The remaining 2.77x is mostly the inherent 2x** from UHF's two spin channels:
  the coupled alpha/beta z-vector (~3.9s), the reference UHF gradient (~4.2s, =
  MRSF's own), and the doubled Fock builds. This is not algorithmic waste.

- **The one remaining big lever is an XC-kernel cache for the z-vector**
  (analysis: `.claude/.../workflows/` run this session; verdict CONTAINED-FEASIBLE).
  Each z-vector matvec re-evaluates the f_xc kernel at the FIXED reference density
  on the grid (`utddft_fxc` -> `run_xc`: `compAOs` + `compXC` are identical across
  matvecs; only the trial-density contraction changes). A "compute-once /
  apply-per-matvec" cache would save ~1.0-1.4s (2.77x -> ~2.5x) at ~0.5 GB memory
  (O(N^2) scaling, needs a large-system fallback). It was NOT implemented: the
  reference state to cache is a large coupled set (aoV, aoG1, reference dRho, the
  libxc 1st/2nd-derivative kernel arrays in `xce%XCLib`, weights, pruning maps),
  and correctly caching+restoring it into the shared engine (or reimplementing the
  GGA f_xc contraction in `dft_gridint_fxc.F90` UUpdate, lines ~379-438) is a
  large, high-iteration change on the most numerically delicate DFT code, gated by
  bit-exactness. Judged not worth the effort/risk for a ~10% gradient gain; revisit
  if UMRSF gradients become a production bottleneck or for much larger systems
  (where the XC re-evaluation dominates and the cache memory is affordable).

- **Accuracy knobs left at default:** `UMRSF_TRKTOL` (amplitude, now moot) and
  `UMRSF_ZTOL` (z-vector). Loosening `UMRSF_ZTOL` to 1e-4 cuts the z-vector from
  17 to 12 matvecs and was byte-identical to printed precision on thymine, but it
  is a genuine accuracy tradeoff (changes the relaxation-density convergence), so
  the default 1e-6 was kept.
