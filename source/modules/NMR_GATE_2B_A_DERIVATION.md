# Gate 2b-A — Amplitude-Coupled MRSF-TDDFT NMR Paramagnetic Response (Formal Derivation)

> **DRAFT — design/derivation only. No code, no integrals, no tests changed.**
> This document is a *formal handoff* for the next NMR gate. It (i) records the
> exact formal status of the validated ground-state RHF/DFT base, and (ii) derives
> the structure of the *amplitude-coupled* paramagnetic response that Gate 2b adds
> on top of the orbital-only response (Gate 2a). Every "reuse X" / "term Y
> survives" statement is a **hypothesis to be confirmed numerically** by the
> acceptance checks in §10 before any production code is written. The load-bearing
> symmetry results (J(P^B)=0, semi-local f_xc[P^B]=0, exchange-only) are *proven and
> numerically confirmed* on the ground state (`NMR_MAGNETIC_SYMMETRY_NOTE.md`,
> Phase-0 gates 0–4); their MRSF generalization is asserted in
> `MRSF_NMR_DESIGN.md` §7 and is **not** re-derived here.

**Branch:** `feat/dft-nmr` (foundational RHF/DFT CGO base; private).
**Prereq (done & validated):** `source/modules/nmr_shielding.F90` — RHF/closed-shell
DFT, common gauge origin, diamagnetic + uncoupled + Phase-0 coupled CPHF/CPKS
paramagnetic. Phase-0 gates 0–4/6 green (`tests/test_nmr_coupled.py`).
**This gate:** Gate **2b-A** — the formal amplitude-coupled derivation. The coding
follow-up is Gate **2b-B** (not in scope here).

Atomic units throughout. `α` = fine-structure constant. `O` = gauge origin (center
of mass). `R_N` = position of nucleus `N`. `r_O = r − O`, `r_N = r − R_N`. The two
magnetic perturbations are the uniform external field `B` (components `B_s`) and
the nuclear magnetic dipole `m_N` (components `m_{N,t}`).

---

## 0. Where Gate 2b-A sits in the gate ladder

| Gate | Content | State |
|------|---------|-------|
| 0–4 (Phase 0) | ground-state RHF/DFT coupled CPHF/CPKS magnetic response; J=0/f_xc=0/exchange-only proven | **DONE & validated** (this branch) |
| 5a-1 | MRSF S0 single-reference reduction limit | downstream (MRSF branch) |
| 5a-2 | MRSF diamagnetic over the relaxed density | downstream (MRSF branch) |
| **2a** | **orbital-only** MRSF paramagnetic (reference orbital response `U^B`; amplitudes & Z-vector frozen) | experimental (MRSF branch) |
| **2b-A** | **formal amplitude-coupled derivation** (this document) | **design-only** |
| 2b-B | implementation of the amplitude/Z-vector response | future |

Gate 2b-A answers one question precisely: *what are the amplitude- and
Z-vector-response terms that Gate 2a omits, what equations produce them, and what
must they reduce to in the single-reference limit?* It deliberately does **not**
re-derive the diamagnetic term (5a-2) or the orbital channel (2a) except where
needed to define the boundary of the new work.

---

## 1. Current formal status (the substrate 2b-A extends)

This is the validated ground-state base. Everything here is implemented in
`nmr_shielding.F90` and matched to the PySCF common-gauge oracle.

### 1.1 Gauge: Common Gauge Origin (CGO) at the center of mass
A single global gauge origin `O = Σ_A m_A R_A / Σ_A m_A` is used for the vector
potential `A_B = ½ B × r_O`. There are **no** gauge-including (London/GIAO) atomic
orbitals; consequently there are **no** `B`-derivatives of the overlap or of the
two-electron integrals, and **no** 2e derivative integrals of any kind. The 2e
response enters *only* through Fock builds acting on the first-order density
(§1.5). Gauge-origin dependence is the known limitation; GIAO is out of scope for
every gate discussed here.

### 1.2 Decomposition `σ = σ_dia + σ_para`
The isotropic nuclear shielding of nucleus `N` is the bilinear coefficient

```
    σ_{ts}(N) = ∂²E / (∂m_{N,t} ∂B_s) |_{B=0,m=0}  =  σ^dia_{ts}(N) + σ^para_{ts}(N),
    σ_iso(N)  = (1/3) Σ_t σ_{tt}(N) · 1e6   [ppm].
```

(Index convention: `t` = nuclear-moment component, `s` = field component; see §2.
The reported quantity is isotropic, so the `t↔s` labeling does not affect output,
but the full-tensor code path fixes a definite ordering — preserve it.)

### 1.3 Diamagnetic term — full tensor (expectation value, no response)
```
    σ^dia_{ts}(N) = (α²/2) [ δ_{ts} Tr(g^N) − g^N_{s,t} ],
    g^N_{ab}      = Σ_{μν} P_{μν} ⟨ μ | (r_O)_a (r_N)_b / r_N³ | ν ⟩,
```
with `P` the SCF (closed-shell) density. Operator `h^{Bm}_{ts} = α²[(r_O·r_N)δ_{ts}
− (r_O)_s (r_N)_t]/r_N³` is **real symmetric**; the term is a pure expectation
value over `P`. Matches PySCF `dia()` to ~6 figures (O 411.418, H 28.062 ppm).
Code: `nmr_dia_shielding` (`int1.F90`) + assembly in `nmr_shielding.F90`.
**Do not touch the α²/2 prefactor or the `(r_O)·(r_N)/r_N³` field factor.**

### 1.4 Paramagnetic term — uncoupled (sum-over-states)
```
    σ^para,unc_{ts}(N) = 2α² Σ_{i∈occ, a∈vir}  L^O_t(a,i) · PSO^N_s(i,a) / (ε_a − ε_i),
```
`L^O = C^T A_O C` (orbital Zeeman about `O`, full antisymmetric AO matrix),
`PSO^N = C^T A_PSO^N C`. The `2α²` collects the `½` of the orbital-Zeeman operator
`h^B = ½ L_O` and the factor 2 from the two sum-over-states cross terms. Matches
the PySCF *uncoupled* reference for both atoms (O −113.63, H 1.785 ppm). Code:
`sum_ov` in `nmr_shielding.F90`. **PSO stays exactly antisymmetric `A=(M−Mᵀ)/2`;
do not touch the `(r−c)/r³` field factor or the α² prefactor.**

### 1.5 Paramagnetic term — Phase-0 coupled CPHF/CPKS
The coupled response replaces the bare denominator quotient by the solution of the
magnetic CPHF/CPKS equation for each field component `s = x,y,z`:
```
    (ε_a − ε_i) R^{(s)}_{ai}  +  c_x · K[ P^B(R^{(s)}) ]_{ai}  =  b^{(s)}_{ai},
    b^{(s)}_{ai} = (½ L_O)_s |_{ai}  (orbital-Zeeman occ–vir block),
    P^B(R) = C (R − Rᵀ) Cᵀ   (imaginary/antisymmetric first-order density),
```
solved by fixed-point iteration; then
```
    σ^para,cpl_{ts}(N) = − 2α² Σ_{ai} R^{(s)}_{ai} · PSO^N_t(i,a).
```
`c_x = HFscale` (1.0 HF, 0.5 BHHLYP, 0.25 PBE0, 0.0 pure GGA). `K` is the
exact-exchange image of the antisymmetric density, built via the `int2` **A−B**
path (`int2_td_data_t(int_amb=.true.)`, `mntoia`). For `c_x = 0` the loop is
skipped and the result is identical to §1.4. HF matches the PySCF coupled oracle
exactly (O para −230.63). Code: `compute_coupled_para` in `nmr_shielding.F90`.

### 1.6 The current approximation: exchange-only, `J(P^B)=0`, `f_xc[P^B]=0`
Because `P^B` is imaginary/antisymmetric in a real AO basis, the first-order charge
density vanishes pointwise (`ρ^B(r) ≡ 0`), so:
- Coulomb response `J(P^B) = 0` — `(μν|λσ)` is symmetric in `λ↔σ`, contracts to 0;
- semi-local XC-kernel response `f_xc[P^B] = 0` — driven by `δρ, ∇δρ, δτ`, all 0;
- only exact exchange `K(P^B) ≠ 0` survives, scaled by `c_x`.

Proven element-wise and confirmed numerically (`NMR_MAGNETIC_SYMMETRY_NOTE.md`,
gates 0–2; gate-3 PBE coupled ≡ uncoupled to machine precision). The single caveat
is the omitted **current-density** XC dependence — exact within the conventional
(non-current) functional class, which is the standard DFT-NMR approximation.

**This §1.6 symmetry is a statement about the two-electron *orbital* coupling
only.** It is the boundary Gate 2b crosses: the MRSF *amplitude* coupling (§5.2) is
a different channel and is **not** killed by this symmetry (see §8, the key
nuance).

---

## 2. Conventions (fixed for the whole derivation)

**MO index spaces** (MRSF / ROHF-ROKS triplet reference `Φ`, two singly-occupied
α orbitals `O1,O2`):

| Symbol | Space |
|--------|-------|
| `i,j,k,l` | occupied (doubly-occupied closed `C` ∪ singly-occupied open `O1,O2`) |
| `a,b,c,d` | virtual `V` |
| `p,q,r`   | general MO |
| blocks `C / O / V` | closed / open(active) / virtual — the MRSF 7-block density structure (`mrsf_density`) |

**Spin:** `σ,τ ∈ {α,β}`. The MRSF excitations are spin-flip `occ_α → vir_β`
(`Δm_s = ±1`); the magnetic operators are **spin-free / spin-conserving**.

**Cartesian:** `s,t,u,v ∈ {x,y,z}`; `s` = field `B`, `t` = nuclear moment `m_N`.
**AO:** `μ,ν,λ,σ_AO` (the AO `σ` is written `σ_AO` to avoid clashing with the
shielding tensor). Occ–vir response vectors are packed `k = (a−1)·n_occ + i` to
match `iatogen`/`mntoia` (as in `nmr_shielding.F90`).

**One-electron operators** (all validated AO matrices already exist):

| Operator | Definition | Real-matrix symmetry | Source |
|----------|-----------|----------------------|--------|
| Orbital Zeeman `h^B_s` | `½(L_O)_s = −(i/2)[r_O×∇]_s` | antisymmetric | `angular_momentum_integrals` |
| PSO `h^m_t` | `(L_N)_t/r_N³ = −i[r_N×∇]_t/r_N³` | antisymmetric | `pso_integrals` |
| Diamagnetic `h^{Bm}_{ts}` | `α²[(r_O·r_N)δ_{ts} − (r_O)_s(r_N)_t]/r_N³` | symmetric | `nmr_dia_shielding` |

**MRSF state:** `E_I = E_ref + ω_I`, eigenpair `(ω_I, X_I)` of the (Tamm–Dancoff)
spin-flip response matrix `A` on the mixed-reference (`M_s = +1` and `M_s = −1`)
average; `rpaeig`/`rparedms`. The orbital-rotation generator is `κ` (anti-Hermitian
occ–vir + active rotations); `U^ξ ≡ dκ/dξ` is the orbital response to perturbation
`ξ`.

---

## 3. Gate 2b-A objective

**Mathematical objective.** Derive the *amplitude-coupled* contribution to the
first-order response of the **relaxed** one-particle density matrix of an MRSF
state to the magnetic field,

```
    γ^{B_s}_I  =  γ^{B_s, orb}_I  +  γ^{B_s, amp}_I ,
```

where `γ^{B_s, orb}_I` (the reference-orbital response channel) is Gate 2a and is
*not* the deliverable here, and

```
    γ^{B_s, amp}_I  ≡  ∂γ_I/∂(amplitudes) · X^{B_s}  +  ∂γ_I/∂(Z-vector) · z^{B_s}
```

is the **Gate 2b** object: the response of the MRSF amplitudes `X^{B_s}` and the
Z-vector orbital-relaxation multipliers `z^{B_s}` to `B_s`. Gate 2b-A delivers the
equations that define `X^{B_s}`, `z^{B_s}`, and `γ^{B_s,amp}_I`, with their
single-reference reduction.

The paramagnetic shielding is then obtained by the interchange theorem (§5.4) as a
single contraction with the PSO operator:

```
    σ^{I,para}_{ts}(N) = Tr[ γ^{B_s}_I · h^m_t ]
                       = Tr[ γ^{B_s,orb}_I h^m_t ]  +  Tr[ γ^{B_s,amp}_I h^m_t ].
                          \_______ Gate 2a _______/    \_______ Gate 2b _______/
```

Gate 2b-A is correct iff (a) the equations below reproduce the ground-state coupled
CPHF/CPKS of §1.5 in the single-reference closed-shell limit, and (b) the
spin-selection vanishings of §7 fall out of the algebra rather than being imposed.

---

## 4. The MRSF energy Lagrangian and its magnetic second derivative

Make the MRSF state energy stationary in all variational and constrained
parameters by the gradient Lagrangian (the same object the analytic MRSF gradient
already builds; `tdhf_mrsf_gradient.F90`):

```
    L_I[X, ω, κ, z, w]  =  E_ref[κ]                                 (reference energy)
                         +  ⟨X_I | (A[κ] − ω_I) | X_I⟩              (amplitude stationarity)
                         +  Σ_{ai} z_{ai} · F_{ai}[κ]               (orbital Brillouin / Z-vector)
                         +  Σ_{pq} w_{pq} (⟨p|q⟩ − δ_{pq}),         (orthonormality)
```

with `L_I` stationary: `∂L_I/∂X = ∂L_I/∂ω = ∂L_I/∂κ = ∂L_I/∂z = 0`. At
stationarity `γ_I = ∂L_I/∂h` is the relaxed 1-PDM (unrelaxed difference density
`T_ij,T_ab` from `X_I` + the Z-vector orbital relaxation + reference density), and
the *first* derivative w.r.t. either magnetic perturbation is a pure expectation
value:

```
    ∂L_I/∂m_{N,t} = Tr[ γ_I · h^m_t ],     ∂L_I/∂B_s = Tr[ γ_I · h^B_s ].
```

Because the magnetic operators carry no `B`/`m` cross dependence beyond the
diamagnetic 1e term, the mixed second derivative splits as

```
    σ^I_{ts}(N) = ∂²L_I/∂m_{N,t}∂B_s
                = Tr[ γ_I · h^{Bm}_{ts} ]      (diamagnetic; γ_I at zero field — Gate 5a-2)
                + Tr[ γ^{B_s}_I · h^m_t ].      (paramagnetic; needs the B-response of γ_I)
```

By the interchange / `2n+1` principle we choose to solve the response to `B`
(3 components total) rather than to `m` (3 components per nucleus); `h^m_t` carries
no `B`-dependence in CGO, so it enters only through the final contraction. This is
exactly the philosophy of the existing Z-vector gradient code and of §1.5.

`γ^{B_s}_I` therefore requires the first-order response of **all** variational
parameters to `B_s`: the orbitals (`U^{B_s}`), the amplitudes (`X^{B_s}`), and the
Z-vector multipliers (`z^{B_s}`). The orbital piece is Gate 2a; the amplitude +
multiplier pieces are Gate 2b.

---

## 5. Response equations

### 5.1 Orbital response `U^{B_s}` (Gate 2a — substrate, stated for completeness)
The spin-conserving CPKS response of the ROHF/ROKS reference orbitals to the
orbital-Zeeman perturbation, with the §1.6 antisymmetric-density coupling:

```
    (ε_a − ε_i) U^{B_s}_{ai} + c_x · K[ P^B(U^{B_s}) ]_{ai} = (½ L_O)_s |_{ai},
```

solved on **both** the `M_s=+1` and `M_s=−1` reference determinants and averaged
(mixed-reference). For pure functionals `c_x = 0` ⇒ uncoupled orbital response;
for hybrids the exchange image survives (§1.6). This equation is *structurally
identical* to §1.5 — that is the design's central reuse claim.

### 5.2 Amplitude response `X^{B_s}` (the Gate 2b deliverable)
Differentiate the amplitude stationarity condition `(A − ω_I) X_I = 0`,
`⟨X_I|X_I⟩ = 1`, with respect to `B_s`:

```
    (A − ω_I) X^{B_s}_I  =  − ( A^{B_s} − ω_I^{B_s} ) X_I ,
    ω_I^{B_s} = ⟨X_I | A^{B_s} | X_I⟩ ,     ⟨X_I | X^{B_s}_I⟩ = 0   (intermediate normalization).
```

The crucial content is the form of `A^{B_s} ≡ dA/dB_s`. The MRSF response matrix
`A` depends on `B` **only through the orbitals** (CGO: no explicit `B` in the 2e
integrals, no London phase), and the *direct* one-electron coupling of the
spin-free `½L_O` operator into the spin-flip block vanishes (§7). Hence

```
    A^{B_s}  =  Σ_{ai} (∂A/∂κ_{ai}) · U^{B_s}_{ai}        (orbital-channel only),
```

i.e. **`B` enters the amplitude equation purely indirectly, through `U^{B_s}`.**
This is what makes Gate 2b strictly downstream of Gate 2a: solve `U^{B_s}` first,
build `A^{B_s}` as the orbital-relaxation derivative of the MRSF response matrix,
then solve the linear CP-MRSF equation above for `X^{B_s}`. The LHS operator
`(A − ω_I)` is the *same* orbital/amplitude Hessian action already implemented for
the MRSF gradient Z-vector (`apply_z_operator` → `sfrogen`/`sfrolhs` +
`int2_driver%run` + `utddft_fxc`), solved with the same CG/GMRES + `sfromcal`
preconditioner. **Only the RHS and the antisymmetric/imaginary density handling
are new.**

### 5.3 Z-vector response `z^{B_s}`
The relaxed density carries the orbital-relaxation multiplier `z` (solution of the
zero-field Z-vector equation). Its first-order response solves a linear equation
with the same Hessian LHS and an RHS built from `U^{B_s}` and `X^{B_s}`:

```
    H · z^{B_s}  =  − ∂/∂B_s ( RHS_z )  =  − [ (∂RHS_z/∂κ) U^{B_s} + (∂RHS_z/∂X) X^{B_s} ],
```

where `H` is the orbital Hessian (LHS of the gradient Z-vector) and `RHS_z` its
zero-field right-hand side. This is the analytic-second-derivative analog of
`sfrorhs`, now differentiated w.r.t. the *imaginary, antisymmetric* magnetic
perturbation. It is the **single highest-risk piece** (cf. `MRSF_NMR_DESIGN.md`
appendix) and must be assembled term-by-term and finite-difference-checked.

### 5.4 Assembly via the interchange theorem
With `{U^{B_s}, X^{B_s}, z^{B_s}}` in hand, build the response density
`γ^{B_s}_I = γ^{B_s,orb}_I + γ^{B_s,amp}_I` (the latter being the `X^{B_s}`- and
`z^{B_s}`-driven pieces of the relaxed-density map), then

```
    σ^{I,para}_{ts}(N) = Tr[ γ^{B_s}_I · h^m_t ]    contracted once with PSO, for all N,t.
```

Three `B`-solves (s = x,y,z) suffice for all nuclei and all moment components.

---

## 6. Building `γ^{B_s, amp}_I`

The relaxed 1-PDM is `γ_I = γ_ref + T[X_I] + γ_z[z]`, with `T` the unrelaxed
amplitude difference density (`T_ij = −Σ_a X_{ia}X_{ja}`-type, `T_ab = Σ_i
X_{ia}X_{ib}`-type, MRSF 7-block) and `γ_z` the orbital-relaxation density from the
Z-vector. Linearizing in `B_s`:

```
    γ^{B_s,amp}_I  =  (∂T/∂X) · X^{B_s}_I        (amplitude-relaxation density)
                   +  (∂γ_z/∂z) · z^{B_s}        (Z-vector-relaxation density)
                   +  (∂γ_z/∂κ) · U^{B_s} | restricted to the parts not already in γ^{orb}.
```

The bookkeeping boundary between `γ^{orb}` (Gate 2a) and `γ^{amp}` (Gate 2b) must
be fixed once and documented in code; the natural split is: `γ^{orb}` = everything
obtained with `X^{B_s}=0, z^{B_s}=0` (orbitals respond, amplitudes/multipliers
frozen); `γ^{amp}` = the remainder. `γ^{B_s}_I` must be **antisymmetric** to
machine precision (same diagnostic as the PSO `max|A+Aᵀ|` and the Phase-0 gate-0
`max|P^B+P^Bᵀ|`).

---

## 7. Spin-flip selection rules (which terms survive)

The magnetic operators are spin-free; the MRSF excitations are spin-flip. The
following must emerge from the algebra (not be imposed):

- **Direct amplitude driving vanishes.** `⟨φ_α|h^B|φ_β⟩ ∝ ⟨α|β⟩ = 0`, so the
  *explicit* `B`-derivative of the spin-flip amplitude equations is zero — hence
  `A^{B_s}` has the orbital-channel form of §5.2 (no direct 1e term).
- **SOS cross-spin terms vanish.** In the sum-over-states picture, `⟨I|h^B|J⟩ = 0`
  unless `I,J` share the spin-multiplicity component; the para sum runs only over
  spin-allowed intermediates.
- **Surviving channel = spin-conserving orbital response of both references.** `B`
  enters via `U^{B_s}` (`α→α`, `β→β`) on the `M_s=±1` references; this propagates
  into `X^{B_s}`, `z^{B_s}` and thence `γ^{B_s}_I`. This is the physical
  paramagnetic channel, and it is exactly why the existing relaxed-density +
  Z-vector machinery is the right substrate.

**Numerical check (gate 2b-A acceptance):** the direct spin-flip `B`-coupling must
be zero to working precision; the mixed-reference contributions must be included
symmetrically from both `M_s=±1` responses.

---

## 8. Inside vs. outside the current approximation

**Inside (carried over from the validated base, §1.6):**
- CGO at COM; no GIAO; no `∂S/∂B`, no 2e `B`-derivative integrals.
- `J(P^B) = 0`, semi-local `f_xc[P^B] = 0`; **two-electron orbital coupling is
  exact exchange only**, scaled by `c_x`. Holds for the *orbital* response
  `U^{B_s}` (§5.1) in MRSF identically to the ground state.
- Diamagnetic = expectation value over the relaxed density (no response).
- Interchange: solve `B`-response only; contract with PSO.

**Outside (the genuinely new Gate-2b content, beyond Gate 2a):**
- The **amplitude response `X^{B_s}`** (§5.2) and **Z-vector response `z^{B_s}`**
  (§5.3), and their density images `γ^{B_s,amp}_I` (§6).
- The **mixed-reference (`M_s=±1`) coupling** of the orbital/amplitude responses.
- The open-shell/active-space contributions to the relaxed density and its
  response (the 7-block structure), absent in closed-shell RHF NMR.

**KEY NUANCE — do not over-extend the §1.6/gate-3 result.** The "coupled ≡
uncoupled for pure functionals" result (gate 3) is a statement about the
*two-electron orbital* coupling: it vanishes for pure functionals because `J=0`,
`f_xc=0`, `c_x=0`. The **amplitude coupling of §5.2 is a different channel** — it is
the relaxation of the correlated wavefunction, not a 2e-density-symmetry effect —
and it is generally **nonzero even for pure GGA functionals**. Therefore Gate 2b
must be exercised for pure functionals too; one must *not* assume "pure ⇒ amplitude
coupling = 0." The only limit in which `γ^{B_s,amp}_I → 0` is the single-reference
closed-shell limit, where the MRSF amplitudes and Z-vector themselves degenerate
away (§9).

**Still outside everything here:** GIAO/London gauge invariance; current-density
functionals; open-shell UMRSF references; spin–spin coupling; magnetizabilities.

---

## 9. Single-reference reduction (the non-negotiable check)

In the closed-shell, single-reference limit (a molecule well described by one
determinant, no static correlation), the MRSF relaxed density `γ_I → P` (the
closed-shell SCF density) and the amplitude/Z-vector relaxation contributions
vanish: `X^{B_s} → 0`, `z^{B_s} → 0`, hence `γ^{B_s,amp}_I → 0`. The full response
collapses to the orbital channel `γ^{B_s}_I → γ^{B_s,orb}_I = U^{B_s}`-only, which
is *exactly* the ground-state coupled CPHF/CPKS of §1.5:

```
    (ε_a − ε_i) U^{B_s}_{ai} + c_x K[P^B(U^{B_s})]_{ai} = (½L_O)_s|_{ai},
    σ^para,cpl_{ts}(N) = −2α² Σ_{ai} U^{B_s}_{ai} PSO^N_t(i,a).
```

**This is the load-bearing acceptance criterion for Gate 2b-A's derivation:** the
amplitude-coupled equations of §5.2–§6 must reduce — analytically and numerically —
to §1.5, with `γ^{B_s,amp}_I` measured to be zero (to the gate-4 tolerance) on a
closed-shell single-reference molecule. A nonzero residual there falsifies the
derivation.

---

## 10. Implementation handoff (Gate 2b-B)

### 10.1 Files / modules likely touched
| Item | Kind | Notes |
|------|------|-------|
| `source/modules/mrsf_nmr_shielding.F90` | **new** | assembly + output; the MRSF analog of `nmr_shielding.F90` |
| magnetic RHS builder (analog of `sfrorhs`) | **new** routine in the MRSF Z-vector path | `½L_O` projected to rotation space + `T_ij,T_ab` indirect terms; the §5.3 highest-risk piece |
| `tdhf_mrsf_z_vector.F90` | reuse | `apply_z_operator` (LHS Hessian action), `sfrogen`/`sfrolhs`, `sfromcal` preconditioner |
| `pcg.F90` | reuse | CG/GMRES solver shell for §5.2/§5.3 |
| `tdhf_lib` | reuse + extend | `iatogen`, `mntoia`, `int2_td_data_t` (A−B path); add antisymmetric/imaginary-density mode |
| `int2_compute` | reuse | exchange Fock build on the antisymmetric density (Coulomb branched to 0; `f_xc` skipped for pure) |
| `tdhf_mrsf_gradient.F90`, `mrsf_density` | reuse | relaxed density `γ_I`, energy-weighted `W`, 7-block structure |
| `int1.F90` (`angular_momentum_integrals`, `pso_integrals`, `nmr_dia_shielding`) | reuse | all three 1e ops validated — **do not modify** |
| `pyoqp/.../runfunc.py`, `input_checker.py` | extend | dispatch `td_prop=nmr`; whitelist |
| `include/oqp.h` | extend | C symbol for the MRSF NMR entry point |

No new integrals are required for CGO. The antisymmetric-density exchange build is
already prototyped in `nmr_shielding.F90::compute_coupled_para` (the `int2`
`int_amb=.true.` path) and should be promoted to a shared helper.

### 10.2 Tests / oracles needed
- **SOS-MRSF cross-check (internal oracle).** Assemble `⟨I|L_O|J⟩`, `⟨I|PSO|J⟩`
  from MRSF transition densities + energy denominators (spin-selected, §7). CP-MRSF
  (§5) must agree with SOS-MRSF within manifold-completeness; this is the primary
  oracle for the genuinely-MRSF amplitude term (PySCF offers no MRSF reference).
- **Single-reference reduction (§9)** vs the existing RHF/DFT coupled prototype and
  the PySCF coupled CGO oracle (`tests/fixtures/nmr/pyscf_cgo_reference.json`),
  reusing `tests/test_nmr_coupled.py` tolerances (`TOL_HF_ABS = 5e-2`,
  `TOL_DELTA = 0.12`).
- **Finite-difference magnetic response** on a tiny system: validate the §5.3 RHS
  *term by term* against a numerical `∂/∂B` of the MRSF stationarity conditions
  before trusting the assembled RHS (per the appendix's recommended approach).
- **Antisymmetry diagnostics:** `max|γ^{B_s}_I + (γ^{B_s}_I)ᵀ|` and the PSO
  `max|A+Aᵀ|`/`max|diag|` checks, as regression guards.
- New fixtures: an MRSF SOS reference JSON for a small molecule (e.g. the H₂O test
  geometry, or a diradical for the multireference selling point); the SR-limit case
  reuses the existing oracle.

### 10.3 Risks / gotchas
- **Magnetic RHS derivation (§5.3) is the critical-path risk.** Differentiate the
  *existing* gradient Lagrangian symbolically, swapping the nuclear perturbation for
  `½L_O`; FD-check each term. Do not hand-port `sfrorhs`.
- **Amplitude coupling ≠ 0 for pure functionals** (§8 key nuance). The gate-3
  "coupled ≡ uncoupled" shortcut applies only to the 2e orbital channel; do not skip
  Gate 2b for pure GGA.
- **Antisymmetric/imaginary density inside `apply_z_operator`.** The Coulomb path
  must be guarded to zero (the symmetric-storage builder would otherwise pollute it);
  `utddft_fxc` skipped for pure functionals; only exact exchange retained. Add an
  explicit "antisymmetric density" mode rather than reusing the symmetric path.
- **Mixed-reference (`M_s=±1`) bookkeeping** must include both reference orbital
  responses consistently in `U^{B_s}`, `X^{B_s}`, and `γ^{B_s}_I`.
- **Index packing** `k=(a−1)·n_occ+i` must match `iatogen`/`mntoia` (as in the
  ground-state prototype) throughout the RHS/solve/contraction.
- **Coupling sign** `kappa = 1` was validated against PySCF on the ground state;
  re-confirm the sign of the exchange image after promoting the helper.
- **Tensor index convention** (`t = moment`, `s = field`, §2) must be fixed
  consistently across RHS, contraction, and output; isotropic output hides a swap
  but the full tensor does not.
- **Do not alter** PSO antisymmetry `A=(M−Mᵀ)/2`, the α²/α²·2 prefactors, or the
  `(r_O)·(r_N)/r_N³` field factor.

### 10.4 Acceptance criteria for Gate 2b-B (the next coding gate)
1. **SR reduction (blocking):** on a closed-shell single-reference molecule,
   `γ^{B_s,amp}_I → 0` and the total MRSF S0 shielding matches §1.5 / the PySCF
   coupled CGO oracle to the gate-4 tolerance (≤ 5e-2 ppm/atom for HF; Δ-metric ≤
   0.12 ppm for hybrids).
2. **Two-route agreement:** CP-MRSF (§5) vs SOS-MRSF agree within
   manifold-completeness expectations and are stable w.r.t. the number of MRSF roots
   retained in the SOS cross-check.
3. **Antisymmetry:** `max|γ^{B_s}_I + (γ^{B_s}_I)ᵀ| < 1e-9` for every field
   component and functional.
4. **Finite-difference:** the assembled magnetic RHS (§5.3) reproduces a numerical
   `∂/∂B` of the MRSF stationarity conditions term-by-term on a tiny system.
5. **Spin selection (§7):** the direct spin-flip `B`-coupling is numerically zero;
   both `M_s=±1` contributions are present.
6. **Pure-functional amplitude coupling:** demonstrated nonzero (≫ noise) for a
   multireference/diradical case with a pure GGA, confirming the §8 nuance and that
   the amplitude channel is exercised independently of `c_x`.

Only after 1–6 pass does Gate 2b-B proceed to excited-state shieldings (the
downstream gate-5b capability).

---

## 11. Cross-references
- `source/modules/NMR_SHIELDING_STATUS.md` — validated RHF/DFT base (§1 here).
- `source/modules/NMR_MAGNETIC_SYMMETRY_NOTE.md` — proofs of `J(P^B)=0`,
  `f_xc[P^B]=0`, exchange-only (§1.6 here).
- `source/modules/MRSF_NMR_DESIGN.md` — full MRSF design study; §1–§4 (Lagrangian,
  reuse), §6 (spin-flip structure), §9 (gates), appendix (highest-risk RHS). This
  document is the formal expansion of that appendix's amplitude-RHS item.
- `source/modules/nmr_shielding.F90` — the ground-state implementation (§1.3–§1.5).
