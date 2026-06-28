# Gradient derivative-2e screening lever (`OQP_GRAD_CUTOFF`)

Performance push, derivative side. Companion to PR #236 (MRSF response:
digestion + FP32 + `OQP_MRSF_RESP_CUTOFF`) and PR #238 (SCF progressive
screening). Goal: speed up nuclear gradients (HF, DFT ground state,
MRSF-TDDFT excited state) by tolerating looser settings on work that is
recomputed at the **converged** density, validated against gradient accuracy
(Hartree/Bohr), which is stricter than energy.

**Default `OQP_GRAD_CUTOFF = 1.0d-8`** (the historic Schwarz cutoff was an
unusually tight `1.0d-10`). The env var overrides it; **`=1.0d-10` restores the
historic tight screening byte-for-byte**. `1e-8` is the size-robust choice
(max|ΔG| ≤ ~1e-6 a.u. through the 36-atom systems tested) and sits well below
the regression-suite tolerance (≈5e-5, max-element deviation rounded to 4
decimals), so it does not perturb the small-molecule example references.

## STEP 0 — where the gradient time goes (cc-pVDZ, `OMP_NUM_THREADS=8`)

Per-stage wall (and CPU) from the in-tree `measure_time` markers
(1e-grad / XC-grad / 2e-deriv build):

| system (atoms)      | method | 1e (s) | XC (s) | **2e-deriv (s)** | 2e share |
|---------------------|--------|-------:|-------:|-----------------:|---------:|
| benzene (12)        | HF     | 0.006  |   –    | **1.41**         | 99.6%    |
| benzene (12)        | DFT    | 0.006  | 0.061  | **1.43**         | 94.4%    |
| benzene (12)        | MRSF   | 0.006  | 0.121  | **1.55**         | 92.4%    |
| naphthalene (18)    | HF     | 0.016  |   –    | **5.87**         | 99.7%    |
| naphthalene (18)    | DFT    | 0.016  | 0.127  | **6.01**         | 97.7%    |
| naphthalene (18)    | MRSF   | 0.017  | 0.293  | **6.87**         | 95.7%    |
| pentacene (36)      | HF     | 0.103  |   –    | **41.5**         | 99.7%    |

**Conclusion:** the derivative two-electron (gradient-Fock) build dominates the
gradient for all three methods (>92%, growing with size). The XC-gradient
quadrature is 4–7%; 1e/overlap-derivative is negligible. Every method funnels
through the same `grd2_driver_gen` (HF: 1 density; DFT: hybrid-scaled; MRSF:
`mrsf_2e_grad` with 11 densities), so **one lever in `grd2_driver_gen` helps all
three.**

For MRSF the *whole gradient* (1e+XC+2e) is only ~23% of the wall clock; the
preceding SCF + response + **Z-vector** solve dominates — but that is a separate
work item (Z-vector chip), not this task.

## The lever — `OQP_GRAD_CUTOFF`

`source/integrals/grd2.F90:grd2_driver_gen` screens shell quartets with a
Schwarz block cutoff (and density-weighted fine screen `cutoff2 = cutoff/2`),
now **defaulting to `1.0d-8`** (was a hard-coded, unusually tight `1.0d-10`).
The converged-density gradient tolerates the looser cutoff. `OQP_GRAD_CUTOFF`
overrides it; **`=1.0d-10` restores the historic tight screening exactly**.

## Accuracy (max |Δgradient| vs the tight 1e-10 baseline, Hartree/Bohr)

Exact (gradient values, independent of timing). Consistent across methods;
HF is the most sensitive (full exact exchange):

| cutoff  | benzene HF | benzene DFT | benzene MRSF | verdict |
|---------|-----------:|------------:|-------------:|---------|
| 1e-8    | 5.8e-7     | 4.8e-7      | 5.1e-7       | ≤1e-6 ✓ (conservative) |
| 1e-7    | 3.0e-6     | 3.2e-6      | 3.5e-6       | ≤1e-5 ✓ (**recommended**) |
| 3e-7    | 1.4e-5     | 1.9e-5      | 1.7e-5       | ~1e-5 borderline |
| 1e-6    | 1.8e-4     | 4.8e-5      | 4.6e-5       | ✗ exceeds 1e-5 |
| 1e-5    | 8.8e-4     | 3.1e-4      | 2.9e-4       | ✗ |

The error grows steeply past 1e-7 (derivative integrals amplify the dropped
contributions by ~exponent, so the integral-magnitude Schwarz bound
under-protects the gradient). The **default is `1e-8`** (size-robust, ≤~1e-6);
`OQP_GRAD_CUTOFF=1e-7` is an aggressive override (~3.5e-6 on small systems but
growing past 1e-5 for larger HF — see below); `=1e-10` restores the tight
baseline.

## Speedup — scales with system size

CPU = total compute work (`OMP_WAIT_POLICY=passive`, so OpenMP barrier-wait does
not inflate it); runs serial, no concurrent jobs. Wall tracks CPU (e.g. pentacene
HF cut1e-7: wall 39.4→31.4 s = 1.25×). `OMP_NUM_THREADS=8`.

2e-deriv-build speedup (×) vs the tight baseline:

| cutoff | benzene (12 at) | naphthalene (18) | pentacene (36) | max\|ΔG\| trend |
|--------|----------------:|-----------------:|---------------:|----------------|
| 1e-8   | 1.05            | 1.07             | **1.12**       | 5.8e-7 → 1.3e-6 (size-stable, ≤~1e-6 ✓) |
| 1e-7   | 1.10            | 1.16             | **1.25**       | 3e-6 → 2.5e-5 (grows; >1e-5 past ~18 at, HF) |

HF baseline 2e CPU: benzene 10.9 s, naphthalene 45.7 s, pentacene 310 s
(2e-deriv build is ~constant share, so this is the gradient speedup).

The table is HF; **DFT and MRSF track it within ~0.03×** (same `grd2_driver_gen`):
at cut1e-7, naphthalene HF 1.16× / DFT 1.16× / MRSF 1.13×. DFT/MRSF are also
slightly *less* error-sensitive (HF-exchange scale ≤0.5), so cut1e-7 stays within
1e-5 for them up to ~18 atoms where HF already breaches it.

**Settings:**
- **Default `1e-8`** — size-robust (max\|ΔG\| ≤ ~1.3e-6 a.u. across all sizes &
  methods), **~5–12 %** faster (grows with size). Below the regression-suite
  tolerance (≈5e-5), so the small-molecule example references are unaffected.
- `OQP_GRAD_CUTOFF=1e-7` — **~10–25 %** faster, but max\|ΔG\| reaches ~2.5e-5 on
  ~36-atom HF (exceeds the 1e-5 target). Safe within 1e-5 only for ≲18-atom HF;
  DFT/MRSF (HF-exchange scale ≤0.5) tolerate it better (≤~7e-6 at naphthalene).
- `OQP_GRAD_CUTOFF=1e-10` — restore the historic tight screening exactly.

## Ceiling / honest limits

- Screening reduces quartet **count** far more than **work**: at cut1e-5
  benzene HF drops 55% of quartets but only 32% of CPU — the screened quartets
  are disproportionately cheap (diffuse, low angular momentum). The expensive
  (compact, high-am) quartets have large Schwarz integrals and are never
  screenable, capping any screening lever at ~30% even at absurd cutoffs, and
  ~10% within the gradient-accuracy budget on benzene.

## Geometry-optimization safety

Benzene HF/cc-pVDZ optimize from the (slightly off-eq) acene geometry: baseline
and `OQP_GRAD_CUTOFF=1e-7` both converge in **5 geometry steps** (identical
trajectory), to final energy −230.72234960 Ha (Δ 6e-10) and the same step-size
convergence (max step 0.000136 vs 0.000133). The loosened gradient does not
perturb the optimizer path or the converged minimum.

## Finite-difference cross-check (absolute correctness)

Benzene HF/cc-pVDZ, analytic vs central-difference energy gradient (h=0.005 Å):
**max|analytic − FD| = 2.4e-5 Ha/Bohr**, every component agreeing to 3 sig figs
with a uniform offset consistent with FD truncation at this step (energies logged
to 1e-8). Confirms the analytic baseline gradient is correct; the lever's
deviation from baseline (3e-6 at cut1e-7) is below the FD-truncation floor, so the
loosened gradient also tracks the true gradient. (`fd_grad.py`.)

## Derivative-aware Schwarz screen — TESTED, strictly worse (not shipped)

The natural next idea (and the one an adversarial review flagged): screen on the
*derivative* bound `Schwarz_ij * Schwarz_kl * amp`, with `amp = 2·sqrt(max
exponent over the 4 shells)` approximating |∂_A(ab|cd)| ~ 2α·(raised ERI), to
protect compact (high-derivative) quartets and kill the error cliff. Implemented
and measured (benzene HF) — it is **strictly worse than the plain integral
cutoff on the error-vs-work frontier**:

| computed quartets ≈ | plain cutoff max\|ΔG\| | deriv-aware max\|ΔG\| |
|--------------------:|-----------------------:|----------------------:|
| 800k                | 3.0e-6 (cut 1e-7)      | 8.5e-5 (dcut 3e-6) — 28× worse |
| 730k                | 1.4e-5 (cut 3e-7)      | 4.4e-4 (dcut 1e-5) — 30× worse |

Reason: the integral magnitude (with the existing density-weighted fine screen)
already predicts the actual gradient contribution better than the exponent-
amplified bound — the `amp` weighting wastes work protecting compact quartets
that don't need it while over-screening diffuse quartets that *do* contribute to
the valence gradient. No `amp` recalibration helps (the whole frontier is worse),
so the code was reverted. A *rigorous* derivative-Schwarz (evaluate the raised-
angular-momentum pair `(a+1,b|a+1,b)`) might do better but is a large new kernel
for, at best, a few percent — not pursued.

## Dead ends (measured, not assumed)

- **`OQP_GRAD_PRIM_TOL`** (loosen primitive-pair / Rys tolerance `tol_int`
  20→12→10): **no speedup** (≤1% CPU) — the cc-pVDZ contractions are short, so
  the Rys evaluation cost is insensitive to `dtol`. (tol 20→12 also gives
  *zero* gradient change; tol 10 gives ~9e-6.) Not shipped.
- **FP32 / mixed precision for the 2e-deriv contraction:** the build is
  **eval-bound, not digestion-bound** — naphthalene MRSF 2e (11 densities) is
  only ~17% more CPU than HF 2e (1 density), so `get_density`/contraction is a
  small fraction; the cost is the FP64 Rys kernel. FP32-izing that kernel is a
  large, risky rewrite very likely to exceed the 1e-6 gradient budget (response
  lesson: all-FP32 gave ~6% at ~µeV energy). Not pursued.
- **XC-gradient grid/Φ reuse:** XC-grad is only 4–7% of the gradient, and the
  gradient needs 2nd-derivative AO collocation (`aoG2`) that an SCF-Vxc Φ-cache
  (≤1st derivative) would not hold. Low ceiling, high complexity. Not pursued.

## Reproduce

Build (this Mac, gcc-15): see brief. Run via staged `OPENQP_ROOT`.
`OQP_GRAD_STATS=1` prints the coarse/fine skip + computed quartet counts.
Geometry/inputs + sweep harness: session scratchpad `gen_acenes.py`,
`sweep.py`, `fd_grad.py`.
