# CASSCF orbital-converger framework — implementation notes

Date: 2026-08-01, branch `openqp-fci-option`.
Scope: native CASSCF/SA-CASSCF orbital optimization (`pyoqp/oqp/library/casscf.py`),
new converger framework (`pyoqp/oqp/library/casscf_convergers.py`),
tests (`tests/test_casscf_convergers.py`).

## What OpenTrustRegion turned out to support

- The compiled core ships an OTR bridge, `source/otr_interface.F90`, plus the
  native TRAH machinery in `source/trah_converger.F90` / `source/scf_converger.F90`.
- The OpenTrustRegion library itself exposes a *generic* callback solver
  (`solver(update_orbs, obj_func, n_param, error, settings)` with objective /
  gradient / Hessian-diagonal / Hessian-vector callbacks), so the upstream API is
  not intrinsically SCF-only.
- OpenQP's bridge, however, is **SCF-only in practice**: `update_orbs`, `obj_func`
  and `hess_x_cb` are hard-wired to the SCF converger object (`trah_converger`),
  `calc_fock`, `get_ab_initio_density` and `conv%calc_g_h`/`calc_h_op` — i.e. the
  Fock/density machinery for RHF/UHF/ROHF (`scftype` 1/2/3).  There is no
  Fortran- or Python-facing entry point that accepts an arbitrary cost function
  (an MCSCF orbital-rotation objective) — no `oqp.` symbol exports the OTR solver.
- This build was configured with `ENABLE_OPENTRAH=OFF` (see `source/CMakeLists.txt`:
  the interface file is dropped and nothing links `libopentrustregion`), so even
  the SCF bridge is not linked here.
- Conclusion: a native Python trust-region augmented Hessian was implemented for
  CASSCF.  A clearly-marked optional hook is documented in the
  `casscf_convergers.py` module docstring: if a generic OTR binding (objective +
  gradient + Hessian-vector callbacks) is ever exposed to Python, it can replace
  `_ah_model_step`/`_ah_inner` behind the same `converger=ah` key without
  touching the dispatch or the tests' contracts.

## Design

- `casscf._optimize` gained a **single dispatch point**: it reads
  `mol.config.get('casscf', {}).get('converger', 'twophase')`; anything in
  `{'', 'twophase', 'two-phase', 'default'}` falls through to the *unchanged*
  production code (the framework module is not even imported), everything else
  goes to `casscf_convergers.run_converger(...)`, which returns the same
  `(C, energies, coeffs, history, converged, total_iterations)` contract.
  Verified: a default-path run after the change produces a log byte-identical to
  the pre-change code (timestamps/build-hash/CPU-time lines aside), with the
  same energy and the same macroiteration count (pinned in test (c)).
- The finite-difference Hessian builder was extracted from `_newton_step` into
  `casscf._fd_orbital_hessian` (pure refactor, identical float operations) so
  the AH converger reuses exactly the Hessian the default Newton step builds.
- Convergers:
  - `twophase` — the existing eigenvalue-floored Newton + curvature-gated saddle
    escape, called through the same functions (`_floored_newton_loop`,
    `_escape_saddles`).
  - `ah` — trust-region augmented Hessian: per macroiteration one gradient
    evaluation + one FD Hessian, then the level-shifted AH eigenproblem
    `[[0, a g^T], [a g, H]]` solved in the Hessian eigenbasis with
    microiterations bisecting the scale `a` until the step matches the trust
    radius (`(H - lambda) x = -g`, `lambda <= min(0, w_min)`, so the shifted
    Hessian is PSD and the step is downhill through negative curvature).
    Trust radius adapts on the actual/predicted ratio (shrink < 0.25,
    grow > 0.75 at the boundary); uphill trial steps are rejected and re-solved
    at a smaller radius at zero CI cost; a gradient-converged point with
    curvature below `-ah_saddle_curv_tol` is kicked along the lowest Hessian
    mode (both signs line-searched) and re-converged with the same AH loop,
    accepting only a strict energy gain.  State-averaged CASSCF works unchanged
    because the objective/gradient closure already carries the SA weights.
  - `diis` — Pulay DIIS on (accumulated-rotation, orbital-gradient) pairs over
    the production quasi-Newton (or Powell, following `[casscf] optimizer`)
    step; the extrapolated point is *evaluated* and kept only when its true
    energy beats the plain step, so each iteration is at worst the default
    iteration (one extra CI solve per extrapolated macroiteration).  The
    accumulated rotation is the BCH-truncated sum of accepted steps — exactness
    is not required since every candidate is re-evaluated variationally.
    Ill-conditioned B matrices drop the oldest vectors.
  - `auto` — `ah` with a stagnation watchdog (no objective decrease beyond
    `energy_decrease_tol` for `auto_stagnation` consecutive macroiterations, or
    a collapsed trust region / macroiteration cap).  On stagnation it restarts
    the unchanged two-phase converger from the AH loop's final point.  AH
    accepts only non-increasing steps, so that point is never above the start;
    if AH made no progress the fallback trajectory is exactly today's default.
    `auto` therefore never diverges and at worst matches the default result
    (the reported iteration count then includes the abandoned AH attempts and
    can exceed `max_macro_iterations`).
- Log surface: non-default convergers append `PyOQP converger: ...` trace lines
  (name, DIIS extrapolation tally / AH stagnation / auto fallback, CI-evaluation
  count) to the CASSCF banner block via a `mol._casscf_converger_trace`
  attribute; the default path never sets it, keeping the default log
  byte-identical.

## Honest cost statement

There is no analytic orbital Hessian in this validation-grade CASSCF; the AH
converger reuses the symmetrized finite-difference Hessian from the analytic
gradient: `2*n_par` gradient evaluations = `2*n_par` active-space CI solves
(+ generalized-Fock builds) per macroiteration — the *same* per-iteration cost
as the default Newton converger, plus one CI solve per trial step.  AH
microiterations and rejected-step re-solves act on the fixed quadratic model
and cost no CI solves.  Fine for small/medium active spaces (this path's
target); large systems need an analytic / DM-direct Hessian-vector product
before `ah` scales.

## `[casscf]` keys read by the framework (for the schema/checker)

All read with `dict.get` defaults from `mol.config['casscf']` — no schema entry
required for programmatic use; the orchestrator can add them to
`oqpdata.py` / `input_checker.py` as follows:

| key | type | default | meaning |
| --- | --- | --- | --- |
| `converger` | str | `twophase` | `twophase` \| `ah` \| `diis` \| `auto` (aliases: `two-phase`/`default` -> twophase; `trah`/`augmented-hessian` -> ah) |
| `ah_start_trust_radius` | float | `0.2` | initial trust radius (clipped to the ceiling) |
| `ah_max_trust_radius` | float | `max(max_rotation_norm, 1e-3)` | trust-radius ceiling |
| `ah_min_trust_radius` | float | `1.0e-6` | floor; collapse below it = stagnation |
| `ah_max_micro` | int | `32` | AH scale-bisection microiterations per step |
| `ah_max_rejects` | int | `6` | uphill-step rejections per macroiteration |
| `ah_saddle_curv_tol` | float | `2.5e-2` | deep-negative-curvature escape threshold |
| `ah_saddle_egain_tol` | float | `1.0e-3` | strict gain required to accept an escape |
| `diis_space` | int | `8` | stored rotation/gradient pairs |
| `diis_start` | int | `2` | pairs required before extrapolating |
| `auto_stagnation` | int | `3` | stalled macroiterations before the fallback |

Existing keys keep their meaning and are shared: `max_macro_iterations`,
`gradient_norm_tol`, `energy_decrease_tol`, `step_norm_tol`,
`max_rotation_norm` (also the default AH trust ceiling), `optimizer`
(`diis` accelerates whichever step it selects; `ah`/`auto` replace the step
engine and ignore it), `level_shift` (two-phase only).

Schema/checker integration (2026-08-01, orchestrator): all keys above now have
rows in `OQP_CONFIG_SCHEMA` (`oqpdata.py`) and validation in
`input_checker.py` (`converger` choice set incl. `trah`/`augmented-hessian`
aliases; positive-value checks for the tolerances/counters), so **input files
carry them directly** (`[casscf] converger=ah`).  `ah_max_trust_radius` uses a
`0.0 = auto` sentinel so the schema default does not override the dynamic
`max(max_rotation_norm, 1e-3)` ceiling.  The input-file path is pinned by
`test_converger_key_through_input_schema`; the runtime-injection path remains
covered by the original tests.

## Measured iteration counts (this build, STO-3G, hcore guess, tolerances of `tests/test_casscf.py`)

Hard case (test (b)) — stretched LiH (R = 3.0 A), CAS(2,2), fc=1:

| converger | macroiterations | E(CASSCF) / Eh |
| --- | --- | --- |
| twophase (default) | 7 | -7.7983384275 |
| ah | 7 | -7.7983384275 |
| diis | 6 | -7.7983384275 |
| auto | 7 | -7.7983384275 |

H2O, CAS(4,4), fc=3 (equilibrium; the most differentiating probed case):

| converger | macroiterations | CI evaluations | E(CASSCF) / Eh |
| --- | --- | --- | --- |
| twophase (default) | 20 | (2*n_par+1)/iter + backtracks ~ 500 | -75.0085688882 |
| ah | 16 | 376 | -75.0085688882 |
| diis | 21 | 520 | -75.0085688882 |
| auto | 16 | 376 | -75.0085688882 |

H4 chain (0.740 A spacing), CAS(2,2), fc=1 — two close CASSCF solutions:

| converger | macroiterations | E(CASSCF) / Eh |
| --- | --- | --- |
| twophase (default) | 3 | -2.1153228671 |
| ah | 7 | **-2.1164680613** |
| diis | 3 | -2.1153228671 |
| auto | 7 | **-2.1164680613** |

The AH step follows the negative-curvature direction from the RHF guess into a
deeper CASSCF solution (1.1 mEh below the default's stationary point; both are
variational, below the CASCI energy -2.1147334880).  The default's phase-2
escape does not trigger because the shallow solution's curvature is above the
`2.5e-2` escape threshold.  This is intended AH behaviour, so the tests pin
"never above the default" for H4 rather than energy equality; the two
cross-converger *agreement* systems (<= 1e-8 Eh) are H2O CAS(4,4) and stretched
LiH, plus SA2-CASSCF(2,2)/H4 for the state-averaged check
(twophase -2.1124780053 vs ah -2.1124780069, 3 macroiterations each).

## Test / verification status

- `tests/test_casscf_convergers.py`: 12 passed (~42 s) — agreement (a),
  hard-case counts (b), default-path regression (c) (energy AND iteration count
  pinned to the pre-change baselines recorded on this build: H2O
  -75.0085688882 / 20 iters, H4 -2.1153228671 / 3 iters), H4 never-above
  invariant, SA-CASSCF ah=default, unknown-converger error.
- `tests/test_casscf.py`: 1 passed; `tests/test_casci.py`: 2 passed.
- Default-path log diff vs pre-change baseline: identical except
  timestamps/build-hash/CPU-time lines.

## Unfinished / caveats

- Schema/checker rows for the keys above are intentionally **not** added here
  (concurrent edits own `oqpdata.py` / `input_checker.py`); until they land,
  `[casscf] converger=...` in an *input file* is rejected by the parser.
- The AH Hessian is finite-difference (cost above); a Davidson micro-solver on
  Hessian-vector products only pays off once an analytic Hv product exists —
  the dense (n_par+1) AH eigensolve is exact and CI-free at current sizes.
- `auto`'s fallback grants the two-phase stage a fresh `max_macro_iterations`
  budget, so its *reported* total iteration count can exceed the nominal cap
  (never-diverge is prioritized over the cap's letter).
- Trust-ratio thresholds/factors (0.25/0.75, x0.5/x2, reject x0.25) are fixed
  constants, deliberately not config keys.
- Iteration-count pins in test (c) are change detectors recorded on this build;
  a different BLAS/architecture could legitimately shift them.

## Analytic MCSCF orbital-rotation Hessian (`[casscf] hessian = analytic`)

Date: 2026-08-01, branch `openqp-fci-option`.
Scope: new `pyoqp/oqp/library/casscf_hessian.py`, pluggable-Hessian slot in
`casscf.py` / `casscf_convergers.py`, tests `tests/test_casscf_hessian.py`.
This section supersedes the "no analytic orbital Hessian yet" statements
above: the FD Hessian is now only the (unchanged) default backend.

### What it computes

`casscf._fd_orbital_hessian` differentiates the *variational* orbital
gradient (the CI is re-solved at every displaced point), so its symmetrized
limit is the second derivative of `E(x) = sum_I w_I E_I(C exp(K(x)))` with
`E_I` eigenvalues of the core-folded active Hamiltonian.  The analytic
Hessian therefore has two parts (both implemented; a naive fixed-CI orbital
Hessian alone does NOT reproduce the FD matrix):

1. **Fixed-CI (orbital-orbital) part** from the full-space state-averaged
   1-/2-RDMs `D`, `G` (chemist convention of `casscf._full_rdms`):
   `H_fix = 1/2 (B + B^T)`, `B_kl = Z(h^(l), g^(l))_{p_k q_k} - Z(..)_{q_k p_k}`
   with one-index derivative integrals `h^(l) = [h, K^(l)]`,
   `g^(l) = sum_pos K^(l)-contracted (pq|rs)`, and the directional-derivative
   intermediate
   `Z(t,T)_mn = (Dt)_nm - (tD)_nm + 1/2 sum [G_nqrs T_mqrs + 3 more positions]`
   (for `(t,T) = (h,g)` the pair contraction reproduces the production
   gradient `2(F_qp - F_pq)`; asserted in the tests).  Algebraically this is
   the standard MCSCF orbital Hessian (Siegbahn/Almloef/Heiberg/Roos,
   JCP 74 (1981) 2384; Helgaker/Jorgensen/Olsen MEST Sec. 10.8) written as
   derivative-integral/generalized-Fock ("Y intermediate") contractions so a
   single code path covers ai/ti/ta blocks and the SA case (weighted RDMs).
2. **CI-relaxation part**
   `2 sum_I w_I sum_{J != I} <I|dH/dx_k|J><J|dH/dx_l|I> / (E_I - E_J)`, with
   `dH/dx_k` the effective Hamiltonian of the core-folded derivative
   integrals (folding is linear) and `J` running over the COMPLETE
   determinant-space spectrum from one dense `_symmetric_eigh` of the active
   Hamiltonian, assembled with precomputed spin-summed `E_tu` excitation
   matrices (cached per active space).  Couplings inside the averaged set
   carry `(w_I - w_J)` and cancel exactly for equal-weight SA; a genuine root
   degeneracy with non-zero coupling raises (the objective is non-smooth
   there, and the FD Hessian is equally undefined).  Spin-forbidden partners
   couple to exactly zero and are skipped.

Cost per build: zero `evaluate` calls (no CI solves); one dense active-space
diagonalization + `n_par` derivative-integral transforms (`O(n_par nbf^4)`)
-- roughly one CI-solve-equivalent, reported in the trace as
`analytic hessian builds`.  Memory guard: the `(nact^2, ndet^2)` excitation
stack must stay under 2 GiB, else a clear error advises `hessian = fd`.
Explicit (non-sequential) `cas.active_orbital_indices` /
`core_orbital_indices` selections are rejected at provider construction.

### Wiring (default byte-identical)

- `casscf._optimize` resolves `[casscf] hessian` (dict.get, default `fd`;
  accepted: `fd`/`finite-difference`/`finite_difference`/`numerical`/
  `default` -> FD, `analytic`/`exact` -> analytic, anything else ->
  `ValueError`) into an optional `hess_fn(C, ci_coeffs)` from
  `casscf_hessian.make_hessian_provider` and threads it through
  `_floored_newton_loop` / `_newton_step(hess=...)` / `_escape_saddles` and
  `run_converger(..., hess_fn=...)`.  With the key absent `hess_fn is None`
  and every call site executes the previous FD code path unchanged (pinned by
  the pre-existing default-regression tests, which still pass: energy AND
  iteration counts).
- All convergers consume it: `ah`/`auto` replace the per-macroiteration FD
  build; `twophase`/`diis` feed it to the production Newton step.  The CI
  vectors of the current `evaluate` are reused, so the provider never solves
  a CI itself.
- Trace lines (only when analytic is active, so default logs are unchanged):
  `PyOQP orbital hessian: analytic` and `PyOQP analytic hessian builds: N`.

### Validation (tests/test_casscf_hessian.py, 20 tests)

- Synthetic exactness (development gate, scratch script): random symmetric
  h / 8-fold-symmetric g, FD of the eigenvalue objective in the symmetric
  `exp(sum x_k K_k)` parametrization vs the analytic matrix; the FD error
  scales as h^2 with ratio 9.00 for step ratio 3 (1.6e-3 at h=1e-3, 1.6e-5
  at h=1e-4) for SS root 0/1, SA equal and SA unequal weights -- the
  signature of an exact analytic Hessian.
- 2-RDM source: `make_rdm1_spatial` / `make_rdm2_spatial` (already present
  in `oqp.library.rdm`; no new builder was needed) through
  `casscf._full_rdms`; identity
  `E = sum h*D1 + 0.5 sum (pq|rs)*D2 + Enuc = E_CASCI` asserted to 1e-10 per
  root on all three systems.
- Machinery anchors: `_fold_active` == `fci._active_space` folding to 1e-12;
  dense active Hamiltonian eigenvalues == CI solver energies to 1e-9; Z
  intermediate == production gradient to 1e-10.
- FD agreement gate (`atol 1e-6`, both at the initial RHF orbitals and at
  the converged CASSCF orbitals; residuals are FD truncation, h=1e-4):

  | system | n_par | max|analytic - fd| initial | converged |
  | --- | --- | --- | --- |
  | LiH (3.0 A) STO-3G CAS(2,2) | 11 | 4.3e-7 | 3.7e-7 |
  | H2O STO-3G CAS(4,4) | 12 | 5.2e-7 | 7.6e-7 |
  | H4 STO-3G SA2-CAS(2,2) | 5 | 6.4e-8 | 6.5e-8 |

### Measured 'ah' cost, fd vs analytic hessian (this build, STO-3G)

| case | macroiters fd -> analytic | CI evals fd -> analytic | hessian builds | E (both) / Eh |
| --- | --- | --- | --- | --- |
| LiH (3.0 A) CAS(2,2) | 7 -> 7 | 139 -> 7 | 6 | -7.7983384275 |
| H2O CAS(4,4) | 16 -> 15 | 376 -> 15 | 14 | -75.0085688882 |
| H4 SA2-CAS(2,2) | 3 -> 3 | 23 -> 3 | 2 | -2.1124780069 |

Energies agree to <= 1e-8 Eh (here: to all printed digits); per
macroiteration the CI cost drops from `2*n_par + 1` solves to 1 solve plus
one dense diagonalization.  The one-iteration difference on H2O comes from
the ~1e-7 FD truncation in the fd-Hessian trajectory, not from the model.
The two-phase default with `hessian = analytic` converges to the same LiH
solution (also tested), so the default Newton path can reuse it as-is.

### Schema/checker rows for the orchestrator

`hessian` is read at runtime with a `dict.get` default (tests inject it via
`runner.mol.config`); to carry it in input files add:

| key | type | default | meaning |
| --- | --- | --- | --- |
| `hessian` | str | `fd` | orbital-Hessian backend: `fd` (finite-difference, default; aliases `finite-difference`/`finite_difference`/`numerical`/`default`) \| `analytic` (exact, CI-solve-free; alias `exact`) |

Checker guidance: choice validation over the alias sets above; a warning
(not an error) when `hessian=analytic` is combined with explicit
`cas.active_orbital_indices`, since the provider raises for non-sequential
active spaces at runtime.

### Caveats

- The CI-relaxation term uses the complete dense spectrum: exact and cheap
  for the validation-grade active spaces this path targets, `O(ndet^3)` /
  2 GiB-guarded beyond.  Scaling further needs iterative response solves
  (`(E_I - H) z = P dH/dx |I>`) and a density-direct fixed-CI part -- the
  documented upgrade hook inside `casscf_hessian.py`.
- Near an exact root degeneracy with symmetry-allowed coupling the SA/SS
  objective is non-smooth; the module raises instead of silently returning a
  wrong curvature.
- `_excitation_matrices` is cached (`lru_cache(8)`) keyed by
  `(nact, nalpha, nbeta)`; the stack is orbital-independent so the cache is
  exact across macroiterations and geometries.
