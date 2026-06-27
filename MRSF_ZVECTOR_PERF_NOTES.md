# MRSF Z-vector (CPHF/CPKS) performance work — notes & results

Target: the MRSF-TDDFT analytic excited-state **gradient** Z-vector solve
(`source/modules/tdhf_mrsf_z_vector.F90`). Everything here is **default-OFF**,
gated by environment variables, and validated against the **gradient**
(max-component error in a.u. vs the tight `zvconv=1e-10` baseline), which is far
stricter than the excitation-energy criterion.

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
| benzene (114 BF) |  6 | 1.36 | 1.47 | **4.38** | 0.02 | 0.0003 | 0.87 |
| thymine (156 BF) | 12 | 3.33 | 7.23 | **14.14** | 0.09 | 0.0009 | 1.64 |

**Key findings**
1. **The per-iteration cost is XC-kernel-dominated** for hybrid-DFT MRSF
   (`utddft_fxc`): 75% (benzene) / 66% (thymine) of the sigma build. The 2e
   digestion is only 25–33%; the linear-solver vector algebra is negligible
   (<0.01%); AO/MO transforms are negligible.
2. **Iteration count scales with system size / floppiness** (benzene 6 → thymine
   12 at `zvconv=1e-10`); larger and excited-state-distorted geometries are more
   iteration-heavy, so **cutting iterations is the dominant lever** because it
   scales down *both* XC and int2 linearly.
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

| zvconv | benzene iters | benzene grad err | thymine iters | thymine grad err |
|---|---|---|---|---|
| 1e-10 | 6 | (ref) | 12 | (ref) |
| 1e-8  | 5 (−17%) | 1.9e-6 ✓ | 10 (−17%) | 6.5e-6 ✓ |
| 1e-6  | 4 (−33%) | 1.5e-5 ✗ | 6 (−50%) | 1.4e-4 ✗ |
| 1e-5  | 3 (−50%) | 9.8e-5 ✗ | 5 (−58%) | 2.1e-4 ✗ |
| 1e-4  |  — | — | 3 (−75%) | 5.6e-4 ✗ |

**Recommendation:** `zvconv = 1e-8` keeps the gradient max-error ≤ ~6e-6 a.u. on
both systems while removing ~17% of iterations. Looser values (1e-6 …) give big
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
safeguard makes it safe to leave on. (A rigorous MO-basis projection of the
stored `xk` into the new MO basis — via the MO overlap — would extend the benefit
to large steps; left as a future enhancement.)

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

Gradients tolerate a cutoff up to ~1e-6 (err 7.8e-6). **But** the speedup is
small for ≤156 BF (few integrals are screened out at these thresholds); the
benefit grows with system size (consistent with the response-side observation).

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

## Bottom line / recommendations
- **Biggest realistic gradient-use-case win:** `OQP_MRSF_ZV_WARMSTART` in the
  small-step regime (MD / late optimization): up to ~1 iteration vs ~6–12.
- **Always-safe modest win:** `OQP_MRSF_ZV_CONV=1e-8` (~17% fewer iterations,
  gradient error ≤ ~6e-6 a.u.).
- **Cutoff** is safe to ~1e-6 but only pays off on large systems.
- **Highest-value future target (out of scope here):** the per-iteration XC
  response kernel (`utddft_fxc`), which is 66–75% of the per-iteration cost.
