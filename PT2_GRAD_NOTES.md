# PT2 Numerical Gradients (QDPT/CASPT2 family) — Implementation Notes

Branch: `openqp-fci-option`. Scope: give the native PT2 family
(`method = caspt2 / ms-caspt2 / xms-caspt2 / mrmp2 / mcqdpt2 / xmcqdpt2`)
central-difference numerical gradients and wire them into `runtype=grad`,
`runtype=optimize`, and `runtype=meci` (penalty-function).

## What already existed (reuse findings)

* **Energy path**: `SinglePoint.energy()` (pyoqp/oqp/library/single_point.py)
  dispatches PT2 methods to `oqp.library.caspt2_dyall.native_caspt2_energy`,
  which puts the PT2 total energies of the computed roots into `mol.energies`
  (ascending; single-state methods produce one entry). Repeated in-process
  `mol.update_system(x); sp.energy()` cycles are the established pattern all
  optimizers use, and are **bitwise deterministic** for a fixed guess (verified:
  fresh-process vs in-process round trip differ by exactly 0.0 on H2 mrmp2).
* **Gradient consumers**: every gradient-driven runtype funnels through
  `Gradient.gradient()`:
  * `runtype=grad` → `compute_grad` (runfunc.py);
  * `runtype=optimize` → `StateSpecificOpt.one_step` (libscipy.py), reused
    verbatim by the `lib=oqp` (`OQPOpt`, default) and `lib=geometric` backends —
    the optimizers consume `(energies[istate], grads[istate])` only;
  * `runtype=meci` → `MECIOpt.one_step` with `meci_search=penalty` (default),
    `ubp`, or `hybrid` — **all three consume state energies and state gradients
    only; no NACs** (the penalty method of J. Phys. Chem. B 2008, 112, 405).
* **Existing numerical-derivative machinery**:
  * `Hessian.numerical_hess` (single_point.py) = FD **of gradients** via
    subprocess workers (`grad_wrapper`) — needs a per-displacement *gradient*
    runtype, so it is not directly reusable for an energy-only method (it now
    composes with this work: PT2 numerical Hessian = FD of FD, untested, and
    blocked by the input checker for subprocesses — see "Checker rows").
  * `Hessian._compute_vibrational_intensities` — FD of dipoles/polarizabilities
    with displacement 1.0e-3 Bohr via fresh `Runner`s; the same displacement
    convention is adopted here.
  * `NAC.numerical_nac` — FD of wavefunction overlaps; MRSF-specific.
  * No optimizer backend accepts energy-only callables (scipy is called with
    `jac=True`; the oqp engine consumes `(e, g)`), so supplying gradients at the
    `Gradient` seam — rather than building a second optimizer — is the minimal
    wiring that makes *all* existing drivers (optimize/meci/mecp/ts/mep/neb/irc)
    PT2-capable.

## What was added

1. **`pyoqp/oqp/library/pt2_numgrad.py`** (new): `pt2_numerical_gradient(mol,
   grad_list, sp=None)` — central differences `g_i = [E(x+h e_i) − E(x−h e_i)]/2h`
   over all 3N Cartesians around the **current** geometry, re-running the full
   SCF + (CASCI/CASSCF reference) + PT2 pipeline per displaced point through
   `SinglePoint.energy(do_init_scf=False)` on the live molecule (serial loop,
   2·3N energies). All computed PT2 states share the displaced energies, so the
   full `(nstate, natom, 3)` array is returned at the cost of one state
   (a 2-state MECI gradient pair costs the same 6N energies as one state).
   After the loop the central geometry is restored **by recomputation** (one
   extra energy) and compared to the entry energies (drift > 1e-10 Eh logs a
   closure warning — catches SCF-solution drift mid-loop honestly).
2. **`Gradient.gradient()` dispatch** (single_point.py, surgical): PT2-family
   methods route to `pt2_numerical_gradient`; `hf`/`tdhf` and the unknown-method
   `ValueError` are byte-identical (re-verified post-edit: HF analytic gradient
   and the error text unchanged).
3. **`tests/test_pt2_grad.py`** (new): 4 tests, ~3.3 s total (see below).
4. This notes file.

## State-selection convention

Gradient state indices (`[properties] grad`, `[optimize] istate/jstate/kstate`)
index `mol.energies` **directly**, exactly as for hf/tdhf. Which roots are in
`mol.energies` follows the existing `[pt2]` convention: `root` (alias `state`)
for the single-state methods, `target_roots`/`nroot` for multistate, energies
ascending. Examples:

* ground-state optimize, mrmp2: `[optimize] istate=0` (**required** — the
  schema default `istate=1` addresses TDHF excited states and is out of range
  for a single-state PT2 run; `pt2_numgrad` raises a self-explanatory error);
* S1 optimize, mcqdpt2 2 roots: `[pt2] nroot=2` (or `target_roots=0,1`) +
  `[optimize] istate=1`;
* S0/S1 MECI, mcqdpt2: `[pt2] nroot=2`, `[optimize] istate=0 jstate=1`.

## Config keys (read via `dict.get` with defaults in pt2_numgrad.py)

| key | default | meaning |
|---|---|---|
| `[pt2] grad_step` | `1.0e-3` (Bohr) | central-difference half-step |
| `[pt2] grad_guess` | `cold` | `cold`: each displaced energy re-runs the configured `[guess] type` — bit-for-bit identical to a fresh single point at that geometry. `warm`: temporarily sets guess type `json`, which for the RHF references PT2 enforces keeps the in-memory MOs/density as the SCF start (no file I/O). |
| `[pt2] grad_gap_warn` | `1.0e-5` (Eh) | floor for the root-swap warning (below) |

**Step-size choice**: the pipeline is converged to ~1e-9..1e-10 Eh and bitwise
reproducible, so FD noise ≈ eps/h ≈ 1e-6 Eh/Bohr at h=1e-3, balancing the
O(h²) truncation ≈ h²·|E'''| ≈ 1e-6 Eh/Bohr (optimum h ~ eps^(1/3) ~ 1e-3).
Empirically h=5e-4 vs h=1e-3 gradients differ by 5e-7 Eh/Bohr on H4 mcqdpt2 —
consistent with that estimate. 1e-3 Bohr also matches the default `[hess] dx`
and the vibrational-intensity displacement already in the code base.

**Warm vs cold**: warm starts inherit the *semicanonical* orbitals PT2 commits
to `OQP::VEC_MO_A` (the active block straddles occ/virt, so those MOs are not
an RHF solution), but SCF re-converges through the untouched stored density;
measured cold-vs-warm gradient difference on H4 mcqdpt2 is 6.8e-7 Eh/Bohr
(SCF-tolerance re-convergence noise). `cold` is the default because it is
exactly reproducible against independent single points (the FD-vs-FD test
matches to 0.0); on these small systems warm gives no measurable speedup.

## Root-tracking limitation (honest)

Multistate roots are **energy-ordered**; near a crossing, displaced geometries
can reorder the adiabatic states, silently differentiating "state k by energy
order" and contaminating the FD gradient with the neighbouring surface. The PT2
kernels expose no CI vectors after the run, so cross-displacement overlap
re-identification is not possible. `pt2_numgrad` **detects and warns** (log)
whenever any displaced point's smallest adjacent gap falls below
`max(grad_gap_warn, 2 × the largest displacement-induced state shift)` — the
regime where reordering is possible. Swaps with roots *outside* the computed
set (e.g. single-state PT2 crossing an uncomputed root) are undetectable.
Penalty-function MECI keeps a finite gap and is the recommended way to approach
degeneracies.

## MECI feasibility verdict

**Feasible, gradient-only, using the native driver unchanged.** The native
MECI drivers (`MECIOpt` + its `lib=oqp`/`geometric` wrappers) implement
penalty (default), UBP, and hybrid searches, and all consume only
`energies[i], energies[j], grads[i], grads[j]` — no NAC vectors (the UBP
branching-plane estimate is built from gradient differences across iterations).
No penalty re-formulation was needed. `runtype=meci` with mcqdpt2 S0/S1 on H4
executes end-to-end (test (d), 2 iterations, ~0.8 s). Caveat: the endgame of a
tight MECI convergence walks into the near-degenerate regime where the FD
root-swap warning fires and FD gradients of *both* states degrade; the penalty
parameters (`pen_alpha`, `energy_gap`) keep the working gap finite.

## Tests (tests/test_pt2_grad.py; all sto-3g; 4 passed in 3.32 s)

Run:
```
cd /Users/cheolhochoi/Documents/claude/openqp-dev-private
OPENQP_ROOT=$PWD/.venv/lib/python3.11/site-packages/oqp .venv/bin/python -m pytest tests/test_pt2_grad.py -q
```

* **(a)** H2 mrmp2 (CAS(2,2) = full space, so E2=0 exactly and the PT2 energy
  equals CASCI — the full PT2 code path still executes): 3-point fit
  (0.720/0.735/0.750 Å) locates r_min = 0.735105 Å; the production gradient
  there is < 1e-3 Eh/Bohr axially (< 1e-6 transverse), and matches an
  independent FD (fresh `Runner` per displaced point, exact Bohr coordinates
  via `update_system`) **elementwise to 1e-8** (measured: 0.0 — bitwise).
* **(b)** `runtype=optimize` (default `lib=oqp` engine) from 0.760 Å converges
  to r = 0.734876 Å vs fit 0.735105 Å (|Δ| = 2.3e-4 Å < 1e-3 Å).
  **Measured wall time of the optimize run: 0.50 s** (0.97 s for the whole
  test including the scan fit).
* **(c)** H4 mcqdpt2 2-state `runtype=grad` returns finite `(2, 4, 3)`
  gradients; ascending energies; S1 gradient distinct from S0.
* **(d)** H4 mcqdpt2 S0/S1 penalty-MECI smoke: 2 iterations of the native
  driver execute; log shows the optimization steps and the PT2 numerical
  gradient banner.

## Schema/checker rows for the orchestrator (files owned elsewhere)

The runtime wiring is complete, but two declarative gates remain in files this
change does not own. Until they land, programmatic runs must (as the tests do)
construct the `Runner` with `runtype=energy` and flip
`mol.config['input']['runtype']` afterwards, and inject `grad_*` keys into
`mol.config['pt2']` post-construction.

1. **`pyoqp/oqp/molecule/oqpdata.py`** — add to the `'pt2'` schema block
   (the parser rejects unknown input-file keys without these):
   ```python
   'grad_step': {'type': float, 'default': '1.0e-3'},
   'grad_guess': {'type': string, 'default': 'cold'},
   'grad_gap_warn': {'type': float, 'default': '1.0e-5'},
   ```
2. **`pyoqp/oqp/utils/input_checker.py`** — `_check_casci` currently emits
   `ERROR "CASCI is currently implemented only for energy calculations"` for
   any `runtype != energy` for the whole CASCI/CASSCF/PT2 family. Relax for
   the PT2 family to the runtypes now wired (numerical gradients), e.g.:
   ```python
   PT2_GRAD_RUNTYPES = {"energy", "grad", "optimize", "meci", "mecp", "ts", "mep", "neb", "irc"}
   # in _check_casci: if method is in the PT2 family, accept PT2_GRAD_RUNTYPES
   # (keep the energy-only error for casci/casscf/sa-casscf);
   ```
   plus (optional, mirrors the tdhf checks) validate `grad_step > 0`,
   `grad_guess in {cold, warm}`, and that requested state indices
   (`properties.grad`, `optimize.istate/jstate`) are `< max(1, pt2.nroot,
   len(pt2.target_roots))` with a hint that PT2 indices address `mol.energies`.
3. Note for a future checker row: `runtype=hess` (numerical Hessian) for PT2
   composes as FD-of-FD through subprocess workers that themselves run
   `runtype=grad` input files — it additionally needs row 2 to land so the
   worker inputs pass the checker; left untested here.

## Known boundaries

* Serial displacement loop (2·3N+1 energies per gradient call). The
  `numerical_hess`-style multiprocessing fan-out is a natural extension once
  the checker rows land (workers are external `runtype=grad` inputs).
* Displaced runs use `do_init_scf=False` (the `[scf] init_scf` helper, if
  configured, runs only for the central/first energies — same convention as
  the existing optimizer steps).
* PT2 requires closed-shell RHF references (enforced upstream in
  `native_caspt2_energy`); all of the above inherits that boundary.
