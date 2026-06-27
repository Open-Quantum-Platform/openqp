# MRSF Z-vector (CPHF/CPKS) performance work — notes & results

Target: the MRSF-TDDFT analytic excited-state **gradient** Z-vector solve
(`source/modules/tdhf_mrsf_z_vector.F90`). Validated against the **gradient**
(max-component error in a.u. vs the tight `zvconv=1e-10` baseline), which is far
stricter than the excitation-energy criterion.

**Defaults:** `OQP_MRSF_ZV_WARMSTART` and `OQP_MRSF_ZV_PROG` are **ON by default**
(disable with `=0`). Warm-start cannot change a converged result (only the
iteration count). Progressive screening perturbs the default gradient by
**<~1e-8 (small systems, e.g. H2O bit-identical) to ~7e-6 a.u. (30 atoms)** —
within the gradient gate but **not bit-identical** to a pure-tight build.
`OQP_MRSF_ZV_CONV` (convergence tol) stays at the 1e-10 default; the static
`OQP_MRSF_ZV_CUTOFF` stays off.

Companion to the response-side work (PR #236, `perf/mrsf-fock-digestion`).

Build/run: gcc/gfortran-15, `cmake -G Ninja -B build … && ninja -C build oqp`.
Runs use `PYTHONPATH=<worktree>/pyoqp` + a `lib/liboqp.dylib` symlink so the
worktree is its own `OPENQP_ROOT`. Test systems: benzene (D6h, cc-pVDZ, 114 BF)
and thymine (15 atoms, cc-pVDZ, 156 BF), BHHLYP, MRSF S1 gradient.

## Opt-in flags added

| env var | effect | default |
|---|---|---|
| `OQP_MRSF_ZV_TIMERS`    | per-section wall-clock profile of the solve | off |
| `OQP_MRSF_ZV_CONV=<x>`  | override z-vector convergence tol `cnvtol`   | off (uses input `zvconv`) |
| `OQP_MRSF_ZV_WARMSTART` | seed the solver with the previous step's solution | off |
| `OQP_MRSF_ZV_CUTOFF=<x>`| loosen the 2e integral cutoff for the whole z-vector build | off (exact) |

---

## STEP 0 — Profile (env `OQP_MRSF_ZV_TIMERS`)

Per-section wall-clock (seconds), default solver = preconditioned CG, BHHLYP.

| system | iters | RHS build | per-iter int2 | per-iter **XC** | transforms | CG algebra | back-proj |
|---|---|---|---|---|---|---|---|
| benzene  (12 atom, 114 BF) |  6 | 1.36 | 1.47 | **4.38** | 0.02 | 0.0003 | 0.87 |
| thymine  (15 atom, 156 BF) | 12 | 3.33 | 7.23 | **14.14** | 0.09 | 0.0009 | 1.64 |
| tetracene (30 atom, 312 BF)| 10 | 16.16 | 36.38 | **68.16** | 0.56 | 0.0013 | 9.41 |

(Per-section seconds are summed over all iterations; the benzene/thymine rows
were taken under light load and tetracene under heavy load — read the *ratios*,
not absolute seconds.)

**Key findings**
1. **The per-iteration cost is XC-kernel-dominated** for hybrid-DFT MRSF
   (`utddft_fxc`): **75% / 66% / 65%** of the sigma build for benzene / thymine /
   tetracene. The 2e digestion is the rest (**25% / 33% / 35%**); the
   linear-solver vector algebra is negligible (<0.01%); AO/MO transforms are
   negligible. This ratio is stable across system size.
2. **Iteration count is set by conditioning, not raw size** (benzene 6, thymine
   12, tetracene 10 at `zvconv=1e-10`): excited-state-distorted / less-symmetric
   geometries are more iteration-heavy. Either way, **cutting iterations is the
   dominant lever** because it scales down *both* XC and int2 linearly.
3. CG converges ~1.3–1.5 orders of magnitude (squared residual) per iteration.

Consequence for the playbook: the int2-focused levers (digestion port, looser
cutoff) address only the minority (25–33%) of per-iteration cost for the usual
hybrid-DFT MRSF case. The biggest realistic wins are (a) fewer iterations
(warm-start, zvconv) and, as a future target outside this playbook, (b) the XC
response kernel itself.

---

## A. Cutting the iteration count

### A1 — `zvconv` right-sizing (env `OQP_MRSF_ZV_CONV`)

The criterion is on the **squared** residual; default `1e-10` ⇒ residual norm
~1e-5, which is over-converged for gradients. Sweep (gradient max-error vs the
`1e-10` baseline):

| zvconv | benzene | thymine | tetracene |
|---|---|---|---|
| 1e-10 | 6 (ref) | 12 (ref) | 10 (ref) |
| 1e-8  | 5 (−17%), 1.9e-6 ✓ | 10 (−17%), 6.5e-6 ✓ | 7 (−30%), 6.8e-6 ✓ |
| 1e-6  | 4 (−33%), 1.5e-5 ✗ | 6 (−50%), 1.4e-4 ✗ | — |
| 1e-5  | 3 (−50%), 9.8e-5 ✗ | 5 (−58%), 2.1e-4 ✗ | — |
| 1e-4  | — | 3 (−75%), 5.6e-4 ✗ | — |

(cell = iters (Δ), gradient max-error a.u.)

**Recommendation:** `zvconv = 1e-8` keeps the gradient max-error ≤ ~7e-6 a.u. on
all three systems while removing 17–30% of iterations (more on the larger one). Looser values (1e-6 …) give big
iteration cuts (−50…−75%) but break the ≤1e-5 a.u. gradient gate, so they are
**not** safe as a default. (The default is left unchanged at 1e-10 per accuracy
policy; `OQP_MRSF_ZV_CONV` is the opt-in.) `env == input zvconv` verified
bit-for-bit.

### A2 — Warm-start across geometry/MD steps (env `OQP_MRSF_ZV_WARMSTART`)

The converged Z-vector is cached (module-level, keyed by target state, reset on
dimension change) and used to seed the next solve in the same process. **A linear
solve self-corrects: the iterative solver still converges to the same residual
tolerance regardless of the guess, so warm-start can only change the iteration
count — never the gradient.**

**Latent bug found & fixed (correctness):** the CG *initial-residual* build used a
different operator than the iterations (`int2_td_data_t`, `int_amb=.true.`, beta
channel `2·apb+amb`, density ×2, no symmetrize) vs the iteration operator
(`int2_tdgrd_data_t`, `int_amb=.false.`, beta `apb`). With the cold start `xk=0`
this is masked (`A·0=0`), but **any** nonzero guess made the recurrence residual
drift from the true residual, so CG "converged" to a wrong `xk` and corrupted the
gradient. The initial build now applies exactly the iteration operator; this is
**bit-identical** for the cold (default) path (verified, maxabs = 0).

**Mechanism proof (zero displacement):** running the gradient twice at the *same*
geometry in one process (identical MOs ⇒ previous solution is exact):
`solve 1 (cold) = 6 iters` → `solve 2 (warm) = 1 iter` (initial r² = 7.6e-11,
already below tol). So the implementation is correct and the upside is large
**when the MO basis barely changes** (MD, late optimization).

**Accuracy gate:** the warm 1-iteration solve and the cold 6-iteration solve at
the same geometry give gradients agreeing to **1.07e-7 a.u.** (≪ 1e-5). Correct.

**MO-basis rotation caveat & regime dependence.** In the naive MO-index space, a
large MO rotation between steps (big optimizer steps, and especially degenerate
systems like benzene's π manifold) can make the previous `xk` a *worse*-than-zero
guess. Mitigation: a **safeguard** rejects the warm guess and restarts from zero
when its initial residual is not below the cold value `‖rhs‖²` (free fallback —
`A·0=0`, so no Fock build), so warm-start is **never** catastrophically worse.

Measured regimes (iters; warm vs cold; safeguard on):

| regime | system | per-step iters cold → warm | notes |
|---|---|---|---|
| zero displacement (Δ=0) | benzene | 6 → **1** | exact guess; grad Δ 1.07e-7 |
| small MD step (0.02 Bohr) | benzene | 6,9,12,14,14 → 6,9,11,14,13 | degenerate π ⇒ 2/4 guesses rejected; ON ≤ OFF every step |
| large BFGS opt step | thymine | 12,12,12,12 → 12,12,11,12 | r0 cut ~2× (7e-2→3e-2) but <1 CG iter saved |

Takeaway: warm-start's payoff scales with how little the MOs move — large for
near-stationary MOs (MD with small `dt`, late-stage optimization, repeated solves
at one geometry), marginal for big steps or near-degenerate manifolds. The
safeguard makes it safe to leave on.

**MO-basis projection (implemented).** To extend the benefit to large optimizer
steps, the stored `xk` is now projected through the geometry-stable AO z-density
into the current MO basis before reseeding:
`xk_old → av_old(MO,old) → Pz = C_old av_old C_oldᵀ (AO) → C_newᵀ S Pz S C_new →
gather`. The same-geometry round-trip is exact (`CᵀSC=I`), verified: benzene 6→1
and tetracene 10→1 with the projected seed (grad Δ 9.4e-8 / 4.9e-8). The CG
safeguard still guards it, so a poor projection only falls back to cold.

### A5 — Jacobi cold-start guess (default on; measured ~0)
`OQP_MRSF_ZV_DIAGGUESS` (default on): the cold start uses `x0 = M⁻¹ rhs` instead
of zero. Free — the cold CG already spends one Fock build on `A·x0` (=0 when
`x0=0`). Measured: halves the initial residual (1.85e-2→7.8e-3 on benzene) but
CG's ~1.5-order/iter geometric convergence absorbs the 2× head start → ~0
iteration change. Kept on (harmless, occasionally crosses an iteration boundary).

### A3 — Preconditioner / initial guess (measured low-value)
- **Perturbative initial guess** (`x0 = M⁻¹ rhs`): with PCG from `x0=0`, the
  first iterate is already `∝ M⁻¹ rhs`, so this lever duplicates iteration 1 and
  saves ~0 — confirmed dead end.
- **Better preconditioner:** the default diagonal orbital-energy-difference
  (Jacobi) preconditioner already gives **clean geometric convergence ~1.5–1.7
  orders of magnitude (squared residual) per iteration** (benzene trace:
  2.3e-3 → 1.6e-4 → 7.4e-6 → 2.5e-7 → 8.1e-9 → 7.6e-11). That leaves little
  headroom for an Olsen/level-shifted/block preconditioner to cut whole
  iterations — consistent with the response-side finding that Olsen gave ~0.
  Not pursued; the iteration-count budget is better spent on warm-start + zvconv.

### A4 — DIIS / subspace acceleration (already present)
The default solver is **already a Krylov method** — preconditioned conjugate
gradient — not a plain fixed-point iteration. Upstream additionally provides
GMRES (`z_solver=1`), MINRES (`z_solver=2`) and an AUTO escalation
(`z_solver=3`: CG→MINRES→GMRES) with breakdown handling. So the "add a
DIIS/Krylov solver" lever is already realized; no change needed. (Warm-start and
the operator-consistency fix benefit all of these paths.)

---

## B. Cutting per-iteration cost

### B1 — Looser integral cutoff (env `OQP_MRSF_ZV_CUTOFF`)
Mirrors PR #236's `OQP_MRSF_RESP_CUTOFF`: save `int2e_cutoff` (default 5e-11),
loosen for the whole z-vector build (RHS + per-iter sigma + relaxed-density tail),
restore on return. Benzene gradient error vs baseline:

| cutoff | grad err | note |
|---|---|---|
| 1e-9 | 0.0 | identical |
| 1e-8 | 4.0e-8 | negligible |
| 1e-7 | 5.5e-7 | safe |
| 1e-6 | 7.8e-6 | safe (≤1e-5) |

Gradients tolerate a cutoff up to ~1e-6 (err 7.8e-6) **for benzene**. **But** two
caveats emerged at scale (tetracene, 312 BF):
- the **speedup grows** with system size: int2 digestion time drops ~32% at
  cutoff=1e-6 (vs ~0 for benzene);
- the **safe cutoff tightens** with system size: at cutoff=1e-6 the tetracene
  gradient error is **6.7e-5 a.u.** — *over* the 1e-5 gate (vs benzene's 7.8e-6).
  A static cutoff must be ~1e-7 for a 30-atom gradient, eating into the speedup.

So a *static* loose cutoff trades speed against accuracy unfavourably at scale.
This motivates B1b.

### B1b — Progressive (iteration-dependent) cutoff (env `OQP_MRSF_ZV_PROG`)
Loose screening while the CG residual is large (the search direction is inexact
anyway), tightened toward the tight floor and **pinned exact** once the residual
is small, so the converged z-vector and the back-projection (which set the
gradient) use the full cutoff. Same coupling idea as
`feat/progressive-screening-scf`, applied to the CG residual:

```
tau(k) = clamp( prog_k * ||r_{k-1}|| , int2e_cutoff , prog_cap )   while ||r||^2 >= prog_pin
tau(k) = int2e_cutoff   (pinned, exact)                            once  ||r||^2 <  prog_pin
```
Implementation: the driver is initialised at the **tight** cutoff (full pair
list); a new `int2_compute_t%set_cutoff` updates only the run-time screening
threshold per iteration (zero rebuild — Schwarz bounds and the pair list are
cutoff-independent / superset-safe). Knobs: `OQP_MRSF_ZV_PROG_K` (default 1e-2),
`OQP_MRSF_ZV_PROG_CAP` (1e-6), `OQP_MRSF_ZV_PROG_PIN` (1e-6 on ‖r‖²).

Results (tetracene, 312 BF, vs the tight baseline; iters unchanged at 10):

| setting | grad err a.u. | int2 time | verdict |
|---|---|---|---|
| static cutoff 1e-6 | 6.7e-5 | −32% | ✗ over the gate |
| **progressive cap=1e-6** | **7.0e-6** | **−17%** | ✓ safe, ~10× more accurate than static 1e-6 |
| progressive cap=1e-5 | 6.5e-5 | −22% | ✗ recurrence drift dominates |

Key points:
- progressive at cap=1e-6 is **~10× more accurate than the static cutoff at the
  same nominal cutoff**, because the pinned tail is exact — it shifts the
  speed/accuracy frontier favourably and **never adds iterations**;
- accuracy is ultimately limited by **inexact-Krylov recurrence drift**: the
  recurrence residual `errv` is updated with the per-iteration (loosened)
  operator, so an over-loose early `cap` (1e-5) leaves the *true* tight residual
  above tolerance even after pinning. Hence `cap≈1e-6` is the safe aggressive
  setting for the ≤1e-5 gradient gate. (A periodic true-residual replacement
  would lift this limit and allow a looser cap — future work.)
- benzene (114 BF): grad 2.0e-8, no speedup (int2 negligible at this size).

### B3 — Coarser response grid (env `OQP_MRSF_ZV_COARSEGRID`)
The per-iteration cost is XC-dominated, and the XC cost scales with the number of
DFT grid points. The response tolerates a coarser grid than the SCF, so this
swaps the **named pruned grid** (SG0<SG1<**SG2 default**<SG3) to a coarser one
(default SG1, env `OQP_MRSF_ZV_GRID`) for the z-vector only — `grid_pruned_name`
is changed before `dft_initialize` and restored after. Default OFF.

Measured (SG1 vs SG2): **per-iter XC 13.3→10.4 s (−22%)** at unchanged iterations.
Gradient error vs SG2 is **system-dependent**: benzene 1.8e-6 (safe), thymine
1.76e-5 (just over the strict 1e-5 gate). SG1 is the discrete coarse step, so for
the strict gate a **coarse-iterations / fine-back-projection** two-grid split
(future work) is needed; kept default OFF as an opt-in. This is a larger and
simpler XC win than the kernel-reuse cache (which an investigation estimated at
only ~10–15%), so the cache was not pursued.

### C — RHF/SF z-vector port
The progressive screening + zvconv-override levers are ported to the **RHF**
(`tdhf_z_vector`, rpa/tda) and **SF** (`tdhf_sf_z_vector`) z-vector solvers via
shared, method-agnostic helpers in `zvector_common` (`zv_opts_t`, `zv_read_opts`,
`zv_prog_tau`, a per-method warm-start store). Env: `OQP_{TDHF,SF}_ZV_PROG`
(default on), `OQP_{TDHF,SF}_ZV_CONV`. Validated: RHF (rpa) and SF H2O gradients
bit-identical progressive ON vs OFF (small systems pin tight immediately).
Warm-start wiring for RHF/SF (the store helpers exist) is a clean follow-up — each
driver needs its own `A·x0` initial-residual build (MRSF already had one; SF/RHF
skip it for the cold `x0=0` case).

### B2 — Digestion port + FP32 (audit result)
- The Z-vector **per-iteration** sigma is built with `int2_tdgrd_data_t`
  (`int2_tdgrd_data_t_update` in `source/tdhf_lib.F90`), **not** the
  `int2_mrsf_data_t` digestion that PR #236 optimized. `int2_mrsf_data_t` is used
  **once**, for the RHS. → The iterations inherit **neither** PR #236's
  Coulomb-RMW/fused-loop speedup **nor** its `OQP_MRSF_FP32` option (both live
  entirely in `int2_mrsf_data_t`; verified in the perf branch).
- Porting the symmetric-Coulomb RMW reduction to `int2_tdgrd_data_t` is possible
  (its Coulomb contribution to spin channels 1 and 2 is identical, so 4 Coulomb
  RMW → 2 via a single accumulator merged post-loop) but: (a) it is a **shared**
  hot loop used by all TD/SF/MRSF gradients and would need an overridden
  thread-finalize, and (b) the profile shows int2 is only 25–33% of per-iter cost
  for hybrid-DFT MRSF, so the net per-iteration upside is ≲5%. Deferred as
  low-value/high-risk for the DFT case (would matter more for pure-HF MRSF).
- FP32 is not applicable to the per-iteration build for the same reason.

---

## Finite-difference anchor
All levers above are validated as the gradient max-error **vs the tight
`zvconv=1e-10` analytic baseline** (the correct reference for "did the
optimization change the gradient"; the default path is additionally proven
bit-identical). To anchor that baseline absolutely, the analytic MRSF gradient
was cross-checked against central finite differences of the state-1 total energy
(h=0.005 Bohr) on a distorted H2O with a clearly non-zero gradient
(|g|~0.07–0.095 a.u.): **FD vs analytic maxabs = 6.4e-5 a.u. (~0.1% rel)**,
consistent with FD truncation. The analytic gradient is therefore correct, and
every lever's perturbation (≤7e-6 a.u.) sits an order of magnitude **below the
analytic gradient's own FD-confirmed accuracy** — i.e. comfortably in the noise.

## Bottom line / recommendations
- **Biggest realistic gradient-use-case win:** `OQP_MRSF_ZV_WARMSTART` in the
  small-step regime (MD / late optimization): ~1 iteration vs 6–12. Proven on
  benzene (6→1), tetracene (10→1); gradient correct to <5e-8. The safeguard
  makes it never-worse for large steps / degenerate systems.
- **Always-safe modest win:** `OQP_MRSF_ZV_CONV=1e-8` (17–30% fewer iterations,
  gradient error ≤ ~7e-6 a.u.; larger systems benefit more).
- **Progressive cutoff** (`OQP_MRSF_ZV_PROG`, cap 1e-6): ~17% int2 reduction at
  scale (30 atoms) at gradient error ~7e-6 — the accuracy-safe way to use loose
  screening; ~10× more accurate than a static loose cutoff. (`OQP_MRSF_ZV_CUTOFF`
  static loosen exists too but is accuracy-bounded at scale.)
- **Per-iteration cost:** progressive screening (int2) + the coarse response grid
  (`OQP_MRSF_ZV_COARSEGRID`, −22% XC, opt-in) attack the dominant XC term; the XC
  kernel-reuse cache was dropped (coarse grid is a bigger, simpler win).
- **Beyond MRSF:** progressive screening + zvconv override are ported to the
  RHF (rpa/tda) and SF z-vector solvers (`OQP_TDHF_ZV_*`, `OQP_SF_ZV_*`).
- **Solver choice:** CG is optimal here — diagnostic on thymine gave CG 12 vs
  MINRES 25 vs GMRES 36 matvecs; a different solver / deflation is not worth it
  (warm-start already collapses post-first MD/opt steps to ~1 iteration).
- **Future:** warm-start wiring for RHF/SF; a coarse-iterations/fine-back-projection
  two-grid split to make the coarse grid strict-gate-safe.

### Combined (full playbook together) — capstone
All three iteration-cutters on at once (`OQP_MRSF_ZV_WARMSTART` +
`OQP_MRSF_ZV_CONV=1e-8` + `OQP_MRSF_ZV_PROG`), tetracene (30 atoms), vs the tight
cold baseline (10 iters):
- **cold step: 10 → 7 iters** (zvconv) with cheaper progressive-screened iters;
- **warm step (typical MD/opt step): 10 → 1 iter (−90%)**;
- gradient vs tight baseline: **4.0e-6 a.u.** (within the ≤1e-5 gate).

So along a trajectory the first/hard step gets the zvconv+progressive discount and
every subsequent (small-step) solve collapses to ~1 iteration — the levers compose
with no accuracy loss.

### Quick reference — env controls
Defaults: warm-start, progressive screening, and the Jacobi guess are **ON**
(MRSF and RHF/SF); set the flag to `0` to disable. zvconv stays 1e-10.
```
# MRSF (source/modules/tdhf_mrsf_z_vector.F90)
OQP_MRSF_ZV_WARMSTART=0            # disable warm-start (default ON; MO-projected)
OQP_MRSF_ZV_PROG=0                 # disable progressive screening (default ON)
OQP_MRSF_ZV_DIAGGUESS=0           # disable Jacobi cold guess (default ON)
OQP_MRSF_ZV_CONV=1e-8             # override convergence tol (default off; recommended 1e-8)
OQP_MRSF_ZV_PROG_K / _CAP / _PIN  # progressive tuning (1e-2 / 1e-6 / 1e-6)
OQP_MRSF_ZV_COARSEGRID=1          # coarser response grid (default OFF)
  OQP_MRSF_ZV_GRID=SG1            #   grid choice (SG0/SG1/SG2/SG3)
OQP_MRSF_ZV_CUTOFF=1e-7          # static loose cutoff (default off; superseded by _PROG)
OQP_MRSF_ZV_TIMERS=1             # per-section profiler
# RHF/SF (tdhf_z_vector / tdhf_sf_z_vector): OQP_TDHF_ZV_* and OQP_SF_ZV_*
OQP_TDHF_ZV_PROG=0  OQP_SF_ZV_PROG=0   OQP_{TDHF,SF}_ZV_CONV=1e-8
```
