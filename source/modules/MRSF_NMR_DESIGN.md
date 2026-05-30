# MRSF-TDDFT NMR Magnetic Shielding — Theory & Architecture Design

> **DRAFT — design study, not final implementation doctrine.** The formulation
> below (especially the symmetry argument in §0 and the MRSF Lagrangian
> second-derivative in §1–§4) is a working hypothesis to be *confirmed
> numerically* before any implementation. The central symmetry claim — that the
> imaginary antisymmetric first-order magnetic density makes the Coulomb and
> local/semi-local XC density responses vanish for pure functionals — is
> plausible but unproven here; and the **exact-exchange response survives for
> HF/hybrids and must be included**. §9 defines the numerical validation gates
> that must pass before the roadmap (§10) is acted on. Treat all "reuse X
> directly" statements as hypotheses pending those gates.

**Status:** design study only (no code). Prerequisite — the native RHF/pure-DFT
common-gauge-origin (CGO) NMR shielding prototype — is implemented and validated
(`source/modules/nmr_shielding.F90`, branch `feat/dft-nmr`). Note that this
prototype currently uses the **uncoupled** paramagnetic response; for HF/hybrids
that is an approximation (see §9, gate 4) that the coupled exact-exchange response
must correct before it is layered into MRSF.

**Goal:** determine the formally correct expression for the nuclear magnetic
shielding of an MRSF-TDDFT state, choose a formulation, and lay out an
implementation roadmap that maximally reuses the existing MRSF energy / Z-vector
/ gradient machinery.

Throughout, atomic units are used; `α` is the fine-structure constant; `O` is the
gauge origin; `R_N` the nucleus; `r_O = r − O`, `r_N = r − R_N`. Magnetic
perturbations are the uniform external field `B` and the nuclear magnetic dipole
`m_N`.

---

## 0. Notation and the magnetic perturbations

Minimal coupling `p → p + A` with `A = A_B + A_{m}`,
`A_B = ½ B × r_O`, `A_{m} = α² (m_N × r_N)/r_N³`, expands the electronic energy of
a state `|I⟩` to second order in `(B, m_N)`. The shielding tensor is the bilinear
coefficient

```
σ^I_{ts}(N) = d²E_I / (dB_s dm_{N,t})  |_{B=0, m=0}.
```

Three one-electron operators appear (all already implemented and validated in the
RHF/DFT prototype):

| Operator | Expression | Role | Symmetry of real matrix |
|----------|-----------|------|--------------------------|
| Orbital Zeeman `h^B_s` | `½ (L_O)_s = −(i/2)[(r_O)×∇]_s` | `∂h/∂B_s` | antisymmetric (imaginary) |
| PSO `h^{m}_t` | `(L_N)_t/r_N³ = −i[(r_N)×∇]_t/r_N³` | `∂h/∂m_{N,t}` | antisymmetric (imaginary) |
| Diamagnetic `h^{Bm}_{ts}` | `α²[(r_O·r_N)δ_{ts} − (r_N)_t(r_O)_s]/r_N³` | `∂²h/∂B_s∂m_t` | symmetric (real) |

**Central fact that shapes everything below:** `h^B` and `h^{m}` are *imaginary*
(spin-free, spatial) operators, so the first-order density response to `B` is
*imaginary and antisymmetric*. Consequences:

- Coulomb response `J[δP^B] = 0` (J of an antisymmetric density vanishes).
- For non-current-dependent functionals the XC kernel response `f_xc[δP^B] = 0`
  (the real density change is zero to first order).
- Only the **exact-exchange** response `K[δP^B]` survives, and only for hybrids.

So for **pure functionals the magnetic orbital response is uncoupled** in the
2-electron sense (identical to the ground-state pure-DFT CGO result), and the only
genuine coupling that remains in MRSF is through the **amplitude/relaxation**
structure of the correlated state. This is the single most important design
driver.

---

## 1. Formally correct shielding expression for an MRSF-TDDFT state

An MRSF state `|I⟩` is obtained from a high-spin **ROHF/ROKS triplet reference**
`Φ` (two singly-occupied orbitals O1,O2, both α) by a spin-flip response in the
`occ_α → vir_β` particle–hole space, with the **mixed-reference** correction that
averages the response built on the `M_s=+1` and `M_s=−1` reference determinants to
restore the closed-shell and doubly-excited configurations and remove spin
contamination. The state energy is `E_I = E_ref + ω_I`, with `ω_I, X_I` the MRSF
eigenpair (Tamm–Dancoff in the current code, `rpaeig`/`rparedms`).

Writing the shielding as a second derivative of the (made-stationary) MRSF energy
Lagrangian `L_I(B, m)` and using the interchange theorem (need the first-order
response to only **one** perturbation):

```
σ^I_{ts}(N) =  Tr[ γ_I · h^{Bm}_{ts} ]                         (diamagnetic)
            +  d/dB_s ⟨ ∂L_I/∂m_t ⟩ |_0
            =  Tr[ γ_I · h^{Bm}_{ts} ]                         (diamagnetic)
            +  Tr[ γ^{B_s}_I · h^{m}_t ]                        (paramagnetic)
```

where

- `γ_I` is the **relaxed one-particle density matrix of state I** (the same object
  the MRSF analytic gradient already builds: unrelaxed difference density `T_ij`,
  `T_ab` plus the Z-vector orbital-relaxation and the recovered closed/open-shell
  contributions). The diamagnetic term is therefore a *pure expectation value* —
  **no response is required**, only the relaxed density (already available) and the
  diamagnetic operator (already implemented).

- `γ^{B_s}_I` is the **first-order response of the relaxed density of state I to
  the magnetic field `B_s`**. This is the only object that requires new solver
  work. It contains (i) the orbital response `U^{B_s}` (CPKS-like, imaginary), and
  (ii) the response of the MRSF amplitudes `X^{B_s}` and the Z-vector multipliers,
  all driven by `B_s`.

Equivalently, in sum-over-states form (rigorous, complete-manifold limit):

```
σ^{I,para}_{ts}(N) = α² Σ_{J≠I} [ ⟨I|h^B_s|J⟩⟨J|h^{m}_t|I⟩ + ⟨I|h^{m}_t|J⟩⟨J|h^B_s|I⟩ ]
                                   / (E_I − E_J).
```

Because `h^B` and `h^{m}` are **spin-free**, `⟨I|h|J⟩` vanishes unless `I` and `J`
have the same spin multiplicity component (see §6). This restricts the SOS sum and
makes it natural to evaluate within the MRSF manifold (one diagonalization yields
many `J`).

**For the MRSF ground state `S0`** (the lowest MRSF root) this reduces, in the
single-reference limit, to the ordinary CGO ground-state expression — a useful
internal check.

---

## 2. Which formulation? (coupled-perturbed / Z-vector / quadratic response / SOS / hybrid)

| Formulation | What it computes | Pros | Cons |
|-------------|------------------|------|------|
| **SOS (within MRSF manifold)** | para from MRSF transition densities + energy denominators | trivial to assemble from existing MRSF roots; ideal validation oracle; gives excited-state shieldings cheaply | truncation/incompleteness error (MRSF manifold ≠ complete state space); not size-consistent for production |
| **Coupled-perturbed MRSF (CP-MRSF)** | solve response of orbitals + amplitudes to `B`; contract with PSO | exact within the MRSF model; 3 linear solves (one per `B` component); production-scalable | requires deriving/implementing the CP-MRSF response (amplitude + orbital relaxation) |
| **Z-vector / Lagrangian** | analytic 2nd derivative of the MRSF Lagrangian; first-order response to `B`, contract with `m` | reuses the *existing* MRSF gradient Lagrangian + Z-vector solver; avoids per-`m`-perturbation solves (interchange) | same algebra as CP-MRSF; the perturbation is imaginary/antisymmetric (new RHS + density symmetry) |
| **Quadratic response** | single residue of `⟨⟨h^{m}; h^B, V⟩⟩` at `ω→ω_I` | formally complete TDDFT object | full QR machinery is heavy; for TDA-MRSF it *reduces to* the Lagrangian approach |

**Recommendation (hybrid):**

1. **Production path = Z-vector/Lagrangian analytic second derivative**, realized as
   a **CP-MRSF response to `B`** (imaginary RHS) contracted with the PSO operator,
   plus the diamagnetic expectation over the relaxed density. For a TDA-type model
   this *is* the single-residue quadratic response — formally correct.
2. **Validation/first-PoC path = SOS within the MRSF manifold.** Cheap, reuses MRSF
   roots and transition densities, and provides an independent number to check the
   CP-MRSF implementation and to deliver excited-state shieldings early.

This mirrors how the RHF prototype was validated (analytic CPKS-equivalent result
cross-checked against an independent reference).

---

## 3. Additional derivatives required, starting from the MRSF gradient

The MRSF gradient (`tdhf_mrsf_gradient.F90`) assembles `dE_I/dR` from: the relaxed
1-PDM `γ_I` and energy-weighted density `W`, contracted with `∂h/∂R`,
`∂(μν|λσ)/∂R`, `∂S/∂R` (Pulay/W term), and XC grid derivatives. For the shielding
(a second derivative w.r.t. `B` and `m`) the substitutions/additions are:

**Reused as-is:**
- The relaxed density `γ_I` and `W` machinery (the Z-vector solution at zero field).
- The orbital-Hessian / Fock-response action `apply_z_operator` (`sfrogen`/`sfrolhs`
  + `int2_driver%run` + `utddft_fxc`) — **but** invoked on an *antisymmetric*
  density (see §4).

**New (CGO):**
- `∂h/∂B = ½ L_O` — imaginary 1e operator (have it: `angular_momentum_integrals`).
- `∂h/∂m = h^{PSO}` — imaginary 1e operator (have it: `pso_integrals`).
- `∂²h/∂B∂m = h^{dia}` — real 1e operator (have it: `nmr_dia_shielding`).
- A **Fock build from an imaginary/antisymmetric density** for the hybrid-exchange
  part of the `B`-response (pure functionals: skip — `J=0`, `f_xc=0`).

**New (GIAO, deferred):**
- `∂S/∂B` (imaginary "London" overlap derivative),
- `∂h_core/∂B` (GIAO kinetic + nuclear derivatives),
- `∂(μν|λσ)/∂B` (GIAO 2-electron derivatives) — **not present in OpenQP**; the
  largest single piece of future work.

**Not needed (relative to gradients):** the nuclear-coordinate 2e integral
derivatives `∂(μν|λσ)/∂R`. For CGO there are *no* 2e derivative integrals at all;
the 2e response enters only through the Fock build acting on the imaginary
`B`-response density.

---

## 4. Can the MRSF Z-vector equations be reused?

**Yes, with a well-defined extension.** The LHS of the response is unchanged; the
RHS and the density symmetry change.

- **LHS operator (reused):** `apply_z_operator` in `tdhf_mrsf_z_vector.F90` builds
  the action of the MRSF/orbital Hessian `(A±B)` via `sfrogen → orthogonal_transform
  → int2_driver%run → utddft_fxc → sfrolhs`, solved by `pcg.F90` (CG) or GMRES with
  the diagonal preconditioner `sfromcal`. The same solver and preconditioner apply.

- **RHS (new):** replace the gradient RHS (`sfrorhs`, built from a *real symmetric*
  nuclear perturbation) with the **magnetic RHS**: the orbital-space gradient of the
  MRSF Lagrangian w.r.t. `B_s`, i.e. the `½L_O` perturbation projected into the
  occ–vir (and active) rotation space, plus the indirect terms from `T_ij`, `T_ab`.
  Three RHS vectors (s = x,y,z).

- **Density symmetry (the real subtlety):** the gradient Z-vector propagates a
  *real, symmetric* density; the magnetic response density `δP^{B}` is *imaginary
  and antisymmetric*. Therefore inside the reused operator:
  - the **Coulomb** contribution must be **zero** (it would be for an exact
    antisymmetric density, but the existing builder assumes symmetric storage — this
    path must be guarded/branched);
  - the **XC kernel** `utddft_fxc` contribution must be **zero** for pure functionals
    (real density change is second order) — skip it;
  - only the **exchange** `K[δP^B]` survives (hybrids). The Fock builder must accept
    an antisymmetric density and return the antisymmetric exchange.

  Practically: add an "imaginary/antisymmetric" mode to the response build (or a
  dedicated thin wrapper) rather than reusing the symmetric path unmodified. The
  *solver, preconditioner, MO/AO transforms, and orbital-space bookkeeping are all
  reused.*

- **Interchange/2n+1:** solve only for the `B`-response (3 solves); the `m`
  dependence is obtained by contracting with the PSO operator (no `m`-perturbation
  solve). This is the efficient analytic-derivative choice and matches the existing
  Z-vector philosophy.

> Note: the exploration suggested the Z-vector could be reused "without
> modification." That is true for the *solver and orbital-Hessian action*, but the
> RHS construction and the antisymmetric/imaginary density handling are genuinely
> new and are where the MRSF-specific derivation effort sits.

---

## 5. Required one-electron operators (inventory)

| Operator | Needed for | Status in OpenQP |
|----------|-----------|------------------|
| Angular momentum `L_O` about gauge origin | `∂h/∂B`, SOS, CP-RHS | **Done & validated** (`comp_amom_int1_prim`) |
| Diamagnetic `[(r_O·r_N)δ − r_N r_O]/r_N³` | `σ_dia` | **Done & validated** (`comp_nmr_dia_int1_prim`) |
| PSO `(r_N×∇)/r_N³` | `∂h/∂m`, SOS, para contraction | **Done & validated** (`comp_pso_int1_prim`) |
| Magnetic-field perturbation `½L_O` (+ GIAO `∂h_core/∂B`, `∂S/∂B`) | CP-MRSF RHS / GIAO | `½L_O` available; GIAO pieces **absent** |
| GIAO 2e derivative `∂(μν|λσ)/∂B` | GIAO response | **Absent** (largest future item) |

For the CGO PoC, **no new integrals are required** — the three validated operators
plus an antisymmetric-density Fock build suffice.

---

## 6. Spin-flip structure: vanishing and new terms

The magnetic operators (`L_O`, PSO) are **spin-free / spin-conserving**; the MRSF
excitations are **spin-flip** (`occ_α → vir_β`, `Δm_s = ±1`). This produces clean
selection rules:

- **Vanishing — direct amplitude driving:** the one-electron matrix element of a
  spin-free operator between the reference and a pure spin-flip configuration is
  zero (`⟨φ_α|h|φ_β⟩ ∝ ⟨α|β⟩ = 0`). Hence the *direct* `B`-derivative of the
  spin-flip amplitude equations vanishes; the amplitudes do **not** couple to `B`
  through the spin-flip block.

- **Vanishing — SOS cross-spin terms:** in the SOS expression, `⟨I|h^B|J⟩ = 0`
  unless `I` and `J` share the same spin-multiplicity component, so spin-flip
  intermediates of the "wrong" multiplicity drop out — the para sum runs only over
  spin-allowed `J`.

- **Surviving / new — indirect orbital channel:** `B` enters through the
  **spin-conserving orbital response** `U^B` (`α→α`, `β→β`) of *both* the `M_s=+1`
  and `M_s=−1` references; this propagates into the relaxed density `γ^B_I` and the
  amplitude/Z-vector relaxation via the MRSF Lagrangian. This is the physical
  paramagnetic channel.

- **New vs ordinary (spin-conserving) TDDFT:**
  - the relaxed density carries the **7-component MRSF structure** (closed `C`,
    open-shell `O1,O2`, virtual `V` blocks; cf. `mrsf_density`), so the diamagnetic
    expectation and the PSO contraction must use the full active-space density, not a
    closed-shell density;
  - the **mixed-reference average** introduces contributions from both reference
    determinants' orbital responses and their coupling;
  - the diamagnetic term acquires open-shell/active contributions absent in
    closed-shell RHF NMR.

A formal corollary: because the direct spin-flip coupling to `B` vanishes, the
**leading MRSF NMR working equations are structurally close to a (state-relaxed)
spin-conserving CPKS on the ROHF reference, dressed by the MRSF relaxed density** —
which is exactly why the existing relaxed-density + Z-vector machinery is the right
substrate.

---

## 7. Complexity vs. accuracy

| Variant | New code | 2e coupling handled | Accuracy | Use |
|---------|----------|---------------------|----------|-----|
| **SOS-MRSF** | low (reuse roots + transition densities + 1e ops + denominators) | none (state energies already include it) | limited by manifold completeness; great for trends + validation + excited states | validation oracle; first PoC; excited-state shieldings |
| **Uncoupled MRSF** | medium-low (orbital response from denominators + diamagnetic over relaxed density; amplitudes frozen) | none | **exact orbital part for pure functionals/CGO**; error from neglected amplitude relaxation | quick production estimate for pure functionals |
| **Fully coupled MRSF** | high (CP-MRSF: `U^B` + amplitude/Z-vector relaxation + antisymmetric-exchange Fock build) | full (incl. hybrid exchange) | production quality | final method |

Key point from §0: for **pure functionals + CGO**, "uncoupled" and "fully coupled"
differ *only by the amplitude/relaxation response* — the 2e orbital coupling to the
imaginary density vanishes. For **hybrids** (e.g. BHHLYP, the MRSF workhorse) the
imaginary-exchange response is nonzero and the fully coupled path is required for
quantitative results. This argues for supporting the antisymmetric-exchange build
relatively early if BHHLYP results are the publication target.

---

## 8. Minimum proof-of-concept for publishable results

**Target capability:** common-gauge-origin MRSF-TDDFT isotropic shielding for the
MRSF **ground state S0** and a few **low-lying excited states** of small molecules.

**Minimum technical content:**
1. Diamagnetic term = expectation of `h^{dia}` over the MRSF relaxed density `γ_I`
   (reuse the gradient's relaxed density; **no response**).
2. Paramagnetic term via **two independent routes that must agree**:
   - **SOS within the MRSF manifold** (cheap, reuses roots/transition densities);
   - **CP-MRSF response to `B`** (orbital response; amplitude/Z-vector relaxation),
     CGO, antisymmetric density.
3. Functional: start with a **pure GGA** (uncoupled-exact, simplest) to establish
   correctness; then add the **imaginary-exchange build for BHHLYP** to make the
   numbers representative of production MRSF.

**Why it is publishable:**
- **First excited-state NMR shieldings from MRSF-TDDFT** (a capability few methods
  offer and MRSF is well-suited to, given balanced treatment of near-degeneracies).
- **Ground-state shieldings for multireference/diradical/bond-breaking systems**
  where single-reference DFT NMR fails — MRSF's correct static-correlation treatment
  should give systematically better shieldings. A compelling test set: stretched
  H₂/N₂ along dissociation, ortho/meta/para-benzyne, a carbene, or points near a
  conical intersection.

**Validation for the PoC:** S0 reduces to RHF/DFT CGO NMR in the single-reference
limit (direct check against the existing prototype); CP-MRSF vs SOS internal
agreement; absolute S0 numbers vs CCSD(T)/experiment for small closed-shell
molecules; multireference cases vs CASSCF/MRCI NMR or experiment.

---

## 9. Numerical validation gates (prerequisite to implementation)

The §0 symmetry argument is the load-bearing assumption of the whole design. It
**must be verified numerically on the existing ground-state prototype before any
MRSF code is written.** All gates below use only present infrastructure: the
`int2` Coulomb/exchange Fock builders, the RHF/DFT CGO shielding module, and PySCF
as an external oracle. MRSF generalization (gate 5) is *blocked* until gates 0–4
pass.

Define `P^B` ≡ the first-order density response to a uniform field component `B`
(claimed imaginary, antisymmetric in the AO basis), `P` ≡ the SCF ground-state
density, `J(·)`/`K(·)` ≡ the Coulomb / exact-exchange Fock matrices built from a
given density.

### Gate 0 — First-order magnetic density is antisymmetric (foundational)
- **Compute:** `P^B` (the first-order response of the AO density matrix to a field
  component `B`) and the scalar `‖P^B + (P^B)^T‖_F`.
- **Expect:** `P^B` is antisymmetric (the entire §0 symmetry argument — and hence
  gates 1–4 — depends on this; if it fails, the rest is moot).
- **Pass:** `‖P^B + (P^B)^T‖_F < 1e-10`.
- **Tools:** the orbital response built from `½L_O` (antisymmetric 1e operator) and
  the occupied–virtual amplitudes; this is the same object whose Fock image is
  probed in gates 1–2.

### Gate 1 — Coulomb response vanishes (pure-DFT magnetic response)
- **Compute:** `J(P^B)` from an antisymmetric `P^B`, and the scalar `Tr[J(P^B) P]`.
- **Expect:** `J(P^B) = 0` elementwise and `Tr[J(P^B) P] = 0` to machine precision
  (because `(μν|λσ)` is symmetric in `λ↔σ` while `P^B` is antisymmetric).
- **Pass:** `‖J(P^B)‖_F < 1e-10` and `|Tr[J(P^B) P]| < 1e-10`.
- **Tools:** `int2` Coulomb build on an explicitly antisymmetrized test density.

### Gate 2 — Exact-exchange response is nonzero (HF/hybrid)
- **Compute:** `K(P^B)` from the same antisymmetric `P^B`, scaled by the
  exact-exchange fraction `c_x`, for **three exchange fractions explicitly**:
  **HF** (`c_x=1.0`), **BHHLYP** (`c_x=0.5`, the MRSF workhorse), and **PBE0**
  (`c_x=0.25`). The point is to confirm the response behavior is governed by `c_x`,
  not to assume one value generalizes.
- **Expect:** `c_x·K(P^B) ≠ 0` and antisymmetric for each; magnitude should scale
  ~linearly with `c_x` (same `K(P^B)`, different prefactor).
- **Pass:** for each functional `‖c_x·K(P^B)‖_F` is `O(c_x)` and clearly nonzero,
  and `‖K(P^B)+K(P^B)^T‖_F < 1e-10` (antisymmetry preserved).
- **Note:** also confirm the local/semi-local **XC kernel** response `f_xc[P^B]`
  is ~0 for the semi-local part of each functional (real density change is second
  order), i.e. the only surviving 2e channel is exact exchange. Report the residual
  `‖f_xc[P^B]‖_F` for HF (trivially 0), BHHLYP, and PBE0 — if it is *not* ~0 for a
  given functional, that is itself a finding that changes the design.

### Gate 3 — Pure-DFT: is coupled ≈ uncoupled? (hypothesis, not assumption)
- **Hypothesis (to be tested, not assumed):** for a pure functional the coupled
  and uncoupled magnetic shieldings coincide, because `J(P^B)=0` and `f_xc[P^B]=0`
  and there is no exact exchange.
- **Compute:** the isotropic shielding for a pure GGA (e.g. PBE/BLYP) two ways:
  (a) the current **uncoupled** prototype; (b) a **coupled** CPKS solve (full 2e +
  `f_xc` response). Report the per-atom difference `Δ = σ_coupled − σ_uncoupled`
  for H₂O/STO-3G and one larger basis (e.g. 6-31G*).
- **Acceptance:** **let the numbers set the tolerance.** Tabulate `Δ`; if the
  hypothesis holds it should sit at the level of solver/grid/integral-screening
  noise. Establish that empirical noise floor first (e.g. by tightening solver and
  grid), then declare the pass threshold from it rather than hard-coding `1e-4`. A
  `Δ` materially above the noise floor *falsifies* the pure-functional symmetry
  claim and must be explained (residual `f_xc`, current dependence, implementation
  bug) before proceeding.
- **Significance:** this is the end-to-end test of the symmetry claim for pure
  functionals; its outcome is a measurement, not a foregone conclusion.

### Gate 4 — HF/hybrid: uncoupled ≠ coupled shielding
- **Compute:** the isotropic shielding for **HF** (and **BHHLYP**) two ways as in
  gate 3.
- **Expect:** they **differ** by the exact-exchange response; and the **coupled**
  result matches PySCF's coupled common-gauge NMR (the real reference), whereas the
  uncoupled one does not.
- **Pass:** `|σ_coupled − σ_uncoupled|` is clearly nonzero (≫ noise) AND
  `|σ_coupled − σ_PySCF(coupled, common-gauge)| < 5e-2` ppm per atom for
  H₂O/STO-3G.
- **Significance:** quantifies the known limitation of the current HF prototype and
  proves the exact-exchange response is correctly built. This gate is the
  prerequisite "Phase 0" deliverable in §10.

### Gate 5 — Generalize to MRSF (only after gates 0–4 pass)

Split into a ground-state limit (5a) and the genuinely new excited-state
capability (5b). **5a must pass before any excited-state (5b) development begins.**

#### Gate 5a — MRSF ground-state (S0) limit
- **Compute, in order:** (i) MRSF diamagnetic over the relaxed S0 density (§10
  Phase 1); (ii) S0 paramagnetic via both SOS-MRSF and CP-MRSF, using the validated
  antisymmetric Coulomb-zero / exchange-nonzero / `f_xc`-zero(pure) response.
- **Expect / pass:**
  - MRSF S0 reduces to the gate-3/4 ground-state result in the single-reference
    limit (closed-shell, single-determinant-dominated molecules) — per-atom
    agreement to the gate-3/4 noise floor;
  - CP-MRSF(S0) agrees with SOS-MRSF(S0) within manifold completeness;
  - the magnetic response density is antisymmetric to machine precision (reuse the
    PSO `max|A+Aᵀ|`-style diagnostic).
- **Blocking rule:** do not begin Phase 3 (CP-MRSF) until gate 4 demonstrates a
  correct coupled exact-exchange magnetic response on the ground state; do not
  begin gate 5b until 5a passes.

#### Gate 5b — Excited-state MRSF shielding
- **Compute:** isotropic shielding for low-lying MRSF excited states (S1, T1, …)
  via CP-MRSF, cross-checked against SOS-MRSF over the manifold.
- **Expect / pass:** CP-MRSF and SOS-MRSF agree within manifold-completeness
  expectations; the state-relaxed density and its magnetic response remain
  antisymmetric; results are stable w.r.t. number of MRSF roots retained in the SOS
  cross-check.
- **Note:** this is the novel publishable capability (§8); it is only meaningful
  once 5a establishes that the machinery reproduces the ground-state limit.

**Recommended artifacts:** capture gates 0–4 as standalone numerical checks (small
driver + reference numbers, in the spirit of `tests/test_nmr_shielding.py`) so the
symmetry assumptions are regression-protected before MRSF work begins.

---

## 10. Implementation roadmap

Estimated modules, dependencies, and validation. Each phase ends with a concrete,
checkable result. **Phase 0 is gated by §9 (gates 0–4); MRSF phases by gate 5
(5a before 5b).**

### Phase 0 — Foundations (mostly done)
- **Have:** validated CGO `L_O`, `h^{dia}`, `h^{PSO}` integrals + drivers
  (`int1.F90`), RHF/DFT CGO shielding (`nmr_shielding.F90`), MRSF energy + Z-vector
  + gradient.
- **New:** an **antisymmetric/imaginary-density Fock-response helper** (Coulomb=0,
  `f_xc`=0 for pure, exchange-only) — prototype against the *ground-state* hybrid
  NMR CPHF first (extends the current RHF prototype from uncoupled to coupled for
  hybrids). *Module:* extend `int2`/`tdhf_lib` or add `nmr_response_lib`.
  *Validation:* **this phase is exactly §9 gates 0–4** — antisymmetric response
  density (gate 0), Coulomb-zero (gate 1), exchange-nonzero across HF/BHHLYP/PBE0
  (gate 2), pure-DFT coupled≈uncoupled hypothesis test (gate 3), and HF/BHHLYP
  coupled vs PySCF coupled common-gauge (gate 4). Do not proceed to MRSF phases
  until all pass.

### Phase 1 — MRSF diamagnetic term
- **Module:** `source/modules/mrsf_nmr_shielding.F90` (new); Python dispatch
  (`td_prop=nmr` in `runfunc.py`, whitelist in `input_checker.py`).
- **Work:** read MRSF relaxed density `γ_I` (`OQP_td_p` + `OQP_td_mrsf_density` +
  reference density), contract with `nmr_dia_shielding`.
- **Validation:** S0 diamagnetic matches the RHF/DFT prototype in the
  single-reference limit; open-shell/active contributions sane.

### Phase 2 — SOS-MRSF paramagnetic (validation oracle + excited states)
- **Work:** assemble `⟨I|L_O|J⟩`, `⟨I|PSO|J⟩` from MRSF transition densities
  (`OQP_td_bvec_mo`, `mrsfcbc`) in the MO basis; energy-denominator sum with
  spin-selection (§6).
- **Dependencies:** MRSF roots/energies/transition densities (exist); 1e ops (have).
- **Validation:** internal consistency vs Phase 3; excited-state trends; S0 vs
  ground-state prototype.

### Phase 3 — CP-MRSF paramagnetic (production, CGO)
- **Work:** build magnetic RHS (`½L_O` projected to rotation space + `T_ij`,`T_ab`
  indirect terms) — *new* analog of `sfrorhs`; solve with the reused
  CG/GMRES + `sfromcal` preconditioner using the **antisymmetric-density** operator
  (Phase 0 helper inside `apply_z_operator`'s Fock build); form `γ^B_I`; contract
  with PSO; add diamagnetic.
- **Reuse:** `pcg.F90`, `apply_z_operator` solver shell, `sfrogen`/`sfrolhs`,
  MO/AO transforms.
- **Validation:** matches SOS (Phase 2) within manifold-completeness expectations;
  S0 matches coupled ground-state NMR; pure-functional uncoupled limit matches
  Phase-2-with-frozen-amplitudes.

### Phase 4 — Productionization & paper
- Per-atom tensors + isotropic shielding output; examples; regression tests
  (antisymmetry of magnetic response density; S0/excited reference values).
- Validation set: dissociating H₂/N₂, benzynes, carbene; vs CASSCF/MRCI/CCSD(T)/exp.

### Phase 5 — Future extensions (out of PoC scope)
- **GIAO/LAO** gauge treatment: `∂S/∂B`, GIAO core-Hamiltonian and **2e GIAO
  derivative integrals** (large; gauge-origin independence, basis-set convergence).
- **Coupled CPHF/CPKS hybrid exchange** response made fully general.
- **Hybrid functionals** beyond BHHLYP; current-dependent functionals (if ever).
- **Open-shell / UMRSF** references; **MRSF spin–spin coupling**; magnetizabilities.

### Module/dependency summary

| New/changed | Kind | Depends on |
|-------------|------|-----------|
| `mrsf_nmr_shielding.F90` | new module (assembly + output) | int1 ops, MRSF density/Z-vector, response helper |
| antisymmetric-density Fock-response helper | new lib (or `int2`/`tdhf_lib` branch) | `int2` exchange build |
| magnetic RHS builder (analog of `sfrorhs`) | new routine in MRSF z-vector path | MRSF Lagrangian quantities |
| Python dispatch (`td_prop=nmr`) | `runfunc.py`, `input_checker.py` | C symbol in `oqp.h` |
| GIAO integral suite (Phase 5) | new integrals in `mod_1e_primitives`/`int2` | Rys/GIAO recurrences |

### Cross-cutting validation strategy
1. **Limit check:** MRSF S0 → RHF/DFT CGO NMR (single-reference molecules).
2. **Two-route agreement:** CP-MRSF vs SOS-MRSF.
3. **Operator-level:** `L_O`, PSO, dia already validated vs PySCF; magnetic response
   density must be antisymmetric (diagnostic like the PSO `max|A+Aᵀ|` check).
4. **Absolute accuracy:** vs CCSD(T)/experiment (single-ref) and CASSCF/MRCI (multi-ref).
5. **Gauge dependence:** document CGO origin dependence; quantify and motivate the
   eventual GIAO work.

---

## Appendix — key risks / open theory items

> **HIGHEST-RISK DERIVATION — construction of the magnetic RHS for the MRSF
> Lagrangian, and its reduction to the single-reference limit.** This is the single
> item most likely to be wrong, hardest to validate piecewise, and on the critical
> path for the entire production method. Concretely: the `B`-derivative of the MRSF
> stationarity conditions (the analog of `sfrorhs`) must be assembled from the
> orbital response, the amplitude/Z-vector relaxation, the mixed-reference
> (`M_s=±1`) coupling, and the open-shell active-space terms — for an *imaginary,
> antisymmetric* perturbation. Two non-negotiable checks on the derivation:
> (i) it must reduce *analytically and numerically* to the validated single-reference
> ground-state coupled CPHF/CPKS magnetic RHS in the closed-shell limit (gate 5a);
> (ii) the spin-selection vanishings of §6 (no direct spin-flip coupling to `B`)
> must fall out of the algebra rather than being imposed by hand. Recommended
> approach: derive the RHS by differentiating the *existing* MRSF gradient
> Lagrangian symbolically, term by term, swapping the nuclear perturbation for
> `½L_O`; verify each term against a finite-difference magnetic response on a tiny
> system before assembling the full RHS. Until this derivation is checked, treat
> the §1–§4 "reuse" claims as provisional.

- **Antisymmetric-density 2e build:** verify the `int2` exchange builder yields the
  correct antisymmetric exchange for an imaginary density; confirm Coulomb and
  `f_xc` are exactly zero (pure) or correctly handled (hybrid).
- **Mixed-reference bookkeeping:** ensure both `M_s=±1` reference orbital responses
  are included consistently in `γ^B_I`.
- **Spin-selection assumptions (§6):** confirm numerically (e.g., that the direct
  spin-flip `B`-coupling is zero to working precision).
