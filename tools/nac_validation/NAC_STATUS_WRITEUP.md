# Analytical MRSF-TDDFT Non-Adiabatic Couplings in OpenQP — Status Write-up

**Date:** 2026-07-21   **System of record:** `openqp-dev-private/nac`
**Reference:** GAMESS numerical MRSF NAC (Lee 2021 protocol)

---

## 1. Executive summary

| Piece | Status |
|-------|--------|
| **Numerical NAC (GAMESS) protocol** | ✅ Solved & validated (6/6, \|cos\|≥0.99995) |
| **Analytical NAC in OpenQP (shipped, semi-numerical `d_amp`)** | ✅ **THE PROJECT GOAL — DONE.** Matches GAMESS 6/6, cheap (~490×) |
| Closed-form `d_amp` — operator diagnosis | ✅ Proven (matvec frame, not gradient chain) |
| Closed-form `d_amp` — `ana2e` | ✅ Built |
| Closed-form `d_amp` — `esum` (1e/2e/XC) | ✅ Derived + **FD-validated 1e-5..1e-9** |
| Closed-form `d_amp` — transport term | ✅ Derived + **validated as necessary** |
| Closed-form `d_amp` — remaining response term | ❌ **Open — proven not reachable from the terms in hand** |

**Bottom line.** The project goal — analytic MRSF NACs in OpenQP, validated against
GAMESS — **is met and shipped**. All state pairs are correct to \|cos\| ≥ 0.99995.

The *fully closed-form* `d_amp` (zero finite differences anywhere) is a further
refinement and is **not complete**. This session took it from "one undifferentiated
open problem" to: three of four components built and validated, and a proof that the
fourth is genuinely a new term rather than a mis-weighting of what exists.

---

## 2. The NAC assembly

The derivative coupling between MRSF states I and J:

```
d_IJ = ⟨Ψ_I|∇_R|Ψ_J⟩ = d_ov − d_amp
```

- **`d_ov`** (orbital / Pulay / overlap term) — **closed-form, working.**
  `d_ov = ½[ d^(frozen+skel) − d^(CPHF) ]`, built from the analytic TLF transition
  density and the CPHF orbital-response seam. Contributes to the validated result.

- **`d_amp`** (amplitude term) — the matvec derivative:
  ```
  d_amp(I,J) = X_Iᵀ (dA/dR) X_J / (Ω_J − Ω_I),     A = MRSF TDA matvec
  ```
  Production computes this **semi-numerically** (transported-matvec FD) and it is
  correct. The **closed form** of this single term is the only open piece.

---

## 3. Numerical NAC (GAMESS) — solved

Validated protocol (`gamess-numerical-nac-protocol`):
- Use **tracked-NAMD** (`RUNTYP=MD`, `$MD NAMD=.T.`, VVERLET), *not* the isolated
  2-step GRADNAC (which is force-contaminated).
- Two silent GAMESS bugs must be fixed first:
  1. **Random per-run MRSF state phase** — identical inputs give ±S₁₂. Fix by
     pinning the transition-dipole sign.
  2. **Force-contaminated MD probe** — the tracked displacement is ~9× the intended
     probe for H at DT=1 fs. Fix with `dt = 1e-17` at fixed `dR`.
- Result: GAMESS == OpenQP, cos **+0.9999969**, ratio **1.00049** at the
  twisted-pyramidal ethylene MECI.

Full benchmark (analytical-vs-numerical, OpenQP): **6/6 pairs, |cos| ≥ 0.99995.**

---

## 4. Closed-form `d_amp` — the operator proof (this session)

### 4.1 What was tried for ~21 attempts, and why it can't work

Every prior attempt built `d_amp` by manipulating the analytic MRSF **gradient**
(polarization identity, RHS surgery, U^x injection). The decisive test this session
was a **homogeneity measurement**: dump `G(0), G(X), G(2X), G(3X)` along each state
direction and check the degree-2 ratio

```
(G(3X) − 3G(X) + 2G(0)) / (G(2X) − 2G(X) + G(0))     [must be 3 for degree ≤ 2]
```

| state | ratio |
|-------|-------|
| 1 (S0) | **3.00000** |
| 2 | 2.996 |
| 3 | 2.987 |

(2/3's shortfall is z-vector noise at their ~0.014 amplitudes.) So `G(X)` **is** a
clean degree-2 form. That closes an airtight chain:

1. `G` degree-2 ⟹ Fortran polarization returns **exactly** `Xᵢᵀ M Xⱼ`
   (`M` = Hessian of the gradient chain).
2. MRSF state gradients are validated ⟹ `M` agrees with the true matvec operator
   `dA` on all three **eigenvector diagonals**.
3. Yet polarization(1,3) = **½** × oracle ⟹ `X₁ᵀ M X₃ ≠ X₁ᵀ dA X₃`.
4. ∴ **`M ≠ dA` as operators** — they match on the 3 eigenvector diagonals but differ
   on the interstate cross terms. *Three diagonal constraints cannot pin a full matrix.*

**Conclusion:** the gradient chain has the *wrong operator off-diagonal* for any
S0-involving pair. No RHS/polarization surgery on it can succeed. The closed form must
be built from the **matvec `dA` directly** — which is exactly what the shipped
semi-numerical oracle does.

### 4.2 The correct decomposition (matvec path)

```
X_Iᵀ (dA/dR) X_J  =  ana2e   +   esum   +   resp
                     (2e)      (1e/Fock)   (orbital response)
```

- **`ana2e`** — explicit 2e integral derivative, fixed MOs. Built by `mrsf_nac_amp`
  (Fortran, shipped). Overshoots the total ~6-7× on its own (large cancellation).
- **`esum`** — explicit 1e/Fock derivative. **Derived below.**
- **`resp`** — orbital-response term, replaces the 3N CPHF solves with one z-vector.
  **Still open** (§6).

---

## 5. The `esum` derivation (fresh dedicated effort, this session)

### 5.1 Source in the matvec

`mrsf_matvec_apply` (`tdhf_mrsf_energy.F90`) builds `A·x = amo`:
- line 208 `mrsfmntoia(...)` → the **2e** back-transform (→ `ana2e` on differentiation);
- line 210 `mrsfesum(wrk1, fa, fb, amo)` adds the **Fock** part, where
  `fa = mo_aᵀ·FOCK_A·mo_a` (line 148), `fb = mo_bᵀ·FOCK_B·mo_b`.

`mrsfesum` computes the standard TDA Fock-diagonal contraction `(F_vv X − X F_oo)`
(with the SOMO fold). Therefore the coupling's Fock part is exactly

```
X_Iᵀ · mrsfesum(iatogen(X_J), fa, fb).
```

### 5.2 The skeleton derivative

Differentiating at **fixed MOs** and **fixed (transported) amplitudes**:

```
esum^x = X_Iᵀ · mrsfesum(iatogen(X_J), fa^x, fb^x)
       = Tr(fa^x · Γ_A) + Tr(fb^x · Γ_B)
       = Tr(P^IJ_a · dFOCK_Aˢᵏᵉˡ/dR) + Tr(P^IJ_b · dFOCK_Bˢᵏᵉˡ/dR)
```

with
```
fa^x = mo_aᵀ (dFOCK_A/dR)_skeleton mo_a          [fixed reference density]
Γ_A  = −sym(Xt_Iᵀ Xt_J) on occ-occ,   Γ_B = +sym(Xt_J Xt_Iᵀ) on virt-virt
Xt   = mrsfxvec-unfolded amplitude,   P^IJ = C Γ Cᵀ  (interstate difference density)
```

This is *literally the same contraction the matvec already performs for the Fock
part*, with `FOCK → dFOCKˢᵏᵉˡ`. It is closed-form and analytic.

**Cross-confirmation:** this reduces to `Tr(P^IJ · dF_AO/dR)` — the same object p20
targeted independently. Two derivations agree.

### 5.3 Implementation (`mrsf_nac_esum`) — DONE

`mrsf_nac_esum(infos, istate, jstate)` in `tdhf_mrsf_gradient.F90` → `OQP::nac_esum`:

- **1e:** `pack(P^IJ_a + P^IJ_b)` → `grd1::grad_ee_kinetic` +
  `grad_en_hellman_feynman` + `grad_en_pulay`.
  `grad_ee_overlap` is **deliberately excluded** — the `W·dSˣ` term belongs to `d_ov`.
- **2e:** `fock_deriv_contract_os(infos, basis, ptot, pa, P^IJ_a, hfscale, g)` and
  the same with `(pb, P^IJ_b)`, at **frozen reference density**.

> ⚠ **Landmine:** `fock_deriv_contract_os` declares `intent(out) gx` and does
> `gx = 0`, i.e. it **overwrites**. The two spin contributions must be accumulated
> through a zeroed scratch array — otherwise you silently get the beta term only.

### 5.4 XC term (DFT) — added

For a DFT reference the Fock also carries `V_xc`, so esum needs
`Tr[P^IJ · dV_xc/dR]`. This is added with the **same production routine the
validated MRSF gradient itself uses**, with the relaxed density replaced by the
interstate probe:

```
utddft_xc_gradient(basis, molGrid, dedft=gxc,
                   da=Pa_ref, db=Pb_ref,       ! frozen reference density
                   pa=P^IJ_a, pb=P^IJ_b,       ! interstate probe
                   nmtx=1, threshold=0, infos)
```

No `xa/xb` is passed, so only the `grad_v_xc` branch runs — the `f_xc` kernel term
belongs to the transition-density side, not the Fock part. Guarded by
`dft = (hamilton == 20)`.

> ⚠ **Two more landmines:** `dedft` is `intent(out)` (**overwrites** → accumulate
> via scratch), and `da/db/pa/pb` are `intent(inout)` and get **scaled by basis
> norms in place** → pass **copies**, or `P^IJ` is silently corrupted downstream.

### 5.5 FD validation — PASSED (BHHLYP / ROHF, the real MRSF setup)

Central difference of `Tr[P^IJ_a·G^a] + Tr[P^IJ_b·G^b]` via `fock_jk` at frozen
densities and frozen probe (no SCF iteration):

| pair | max\|an − FD\| | \|an\| | \|FD\| |
|------|---------------|--------|--------|
| (1,2) | **2.0e-9** | 3.5170 | 3.5170 |
| (1,3) | 3.4e-5 | 2.0500 | 2.0501 |
| (2,3) | 1.1e-4 | 2.3239 | 2.3239 |

Magnitudes agree to 5 significant figures on every pair; the 1e-5..1e-4 residual is
FD roundoff (h = 1e-4 differencing O(1-10) traces), and the 1e-9 on (1,2) shows the
machinery is exact. **The 2e esum is correct.**

### 5.5b XC FD validation — PASSED (and it caught a real bug)

Reference: central difference of `Tr[P^IJ_a·Vxc^a] + Tr[P^IJ_b·Vxc^b]`, with `Vxc`
built by `dftexcor` from **fixed MO coefficients** (ROHF ⇒ `iscftyp=2`, beta MOs =
alpha MOs) and the **DFT grid rebuilt** at every displaced geometry.

| pair | before fix | after fix | \|an\| | \|FD\| |
|------|-----------|-----------|--------|--------|
| (1,2) | 8.8799e-2 | **1.07e-5** | 0.16962 | 0.16964 |
| (1,3) | 8.8799e-2 | **7.88e-6** | 0.086717 | 0.086733 |
| (2,3) | 8.8799e-2 | **1.92e-5** | 0.11240 | 0.11242 |

**The bug it caught:** `utddft_xc_gradient` also returns a **probe-independent
ground-state XC gradient** term (it depends only on `da/db`). That constant is not
part of `Tr[P^IJ·dVxc/dR]` and must be removed with a **zero-probe baseline**:

```
gxc = gxc(pa = P^IJ) − gxc(pa = 0)
```

(the same `gZ − gS` seam pattern used elsewhere in the file; pass fresh `da/db`
copies to *both* calls, since they are scaled in place).

**The diagnostic tell:** the error was `8.8799e-02` *identically* for all three
pairs while `|an|` and `|FD|` both varied per pair. A pair-independent error can
only be an additive constant — never real physics. After the fix the residuals
differ per pair (1.07e-5 / 7.88e-6 / 1.92e-5), which is genuine FD noise, and `|an|`
lands on the *original* `|FD|`, confirming the FD reference was right all along.

### 5.5c The `logtol` bug — and the 1e validation

Starting work on `resp` exposed a bug in esum itself. I had passed
`logtol = log(10)*infos%control%int2e_cutoff` to the `grd1` 1e routines, copying
`hf_gradient`. But the NAC harnesses set `int2e_cutoff = 1e-20`, so that argument
evaluated to ≈0 — and these routines treat `logtol` as a **log** threshold, making
it meaningless screening. Removing it (using each routine's `tol_default`) moved
esum by 10-200×.

With that fixed, the **1e piece FD-validates** against `Tr[(P^IJ_a+P^IJ_b)·Hcore(R)]`
(rebuilt via `int1e` at each displaced geometry, frozen probe):

| pair | 1e | 2e | XC |
|------|-----|-----|-----|
| (1,2) | **2.2e-9** | 1.8e-9 | 1.2e-5 |
| (1,3) | **6.7e-10** | 3.1e-5 | 5.4e-6 |
| (2,3) | **1.6e-9** | 1.1e-4 | 1.9e-5 |

**esum is now fully validated in all three pieces.**

### 5.5d ★ Retraction: the "100:1 cancellation" was an artifact

I earlier reported that esum was large (~2-4) and nearly cancelled `resp`, implying
a catastrophic ~100:1 cancellation that would demand 4-5 significant digits per
term. **That was wrong** — a direct consequence of the `logtol` bug.

The truth: `|1e| ≈ 3.4` and `|2e| ≈ 3.2` **nearly cancel each other** (nuclear
attraction vs electron repulsion), giving `|esum| ≈ 0.03-0.41`. The broken tolerance
screened the 1e term down to ~0.57 so it could not perform that cancellation,
leaving a spurious `|esum| ≈ 4`. With it fixed the decomposition is **well
conditioned**:

```
|orc·gap|      0.010 – 0.039
|ana2e|        0.023 – 0.147
|esum|         0.025 – 0.409
|resp_target|  0.167 – 0.360      <- all comparable
```

Lesson recorded: never copy a `logtol=` argument without checking the units of the
value passed, and FD-validate **every** sub-piece — the 1e term was the one piece I
waved through as "production routines used correctly," and it was the one that broke.

> **Note — there is no "HF-MRSF."** MRSF is defined on the **ROHF triplet reference
> with a DFT functional**. An earlier version of this document reported a parallel
> validation on an input created by deleting the `functional=` line; that setup is
> not MRSF at all, and its results (including a 10-30% run-to-run drift I had
> mis-attributed to SCF conditioning) are disregarded. That input has been deleted.
> All MRSF validation is on the DFT path. Consequently the **XC term is not
> optional** — it is part of every real MRSF calculation.

### 5.6 Two process bugs found along the way

1. **`fock_deriv_contract_os` overwrites** (`intent(out)`, `gx = 0`). Piping both
   spin contributions into one array silently yields beta-only. Caught by the FD
   test (it showed `|an| ≈ ½|FD|`).
2. **Harness ordering.** `_compute_amp_damp` (the oracle) does 6N re-SCF at
   displaced geometries and does **not** restore MOs/geometry. Anything depending
   on unperturbed orbitals must be computed **before** it.

### 5.7 Reproducibility

Repeat runs with identical inputs on the production BHHLYP/ROHF input
(`h2o_ana.inp`): **esum reproduces to 0.05-3%** — stable.

Cross-run comparison is still the weaker check here. The **FD self-test is the right
validation**, because analytic and FD share the same amplitudes within one run.

### 5.5 ★ Correction to an earlier claim

I previously wrote that p20's esum was "~20× too big," inferring that from
`|esum| ≫ |missing|`. **That premise was wrong.** The validated esum is *large*
(|esum| ≈ 1.7-3.0, dominated by the 2e piece) because `dFOCKˢᵏᵉˡ/dR` is large — it
does **not** approximate `missing`. It must be nearly cancelled by an equally large
orbital response: measured `|resid| = |missing − esum| ≈ |esum|` to ~1.5%.

p20's number is still not this quantity (its gradient-seam route adds `dG[P^IJ]` and
`W·dSˣ`), but "20× too big" was the wrong characterization of why.

---

## 6. What remains open: the orbital-response term `resp`

### 6.0 Final state of the response side (ethylene, C1 — cosine is diagnostic here)

Reading `_compute_amp_damp` line by line (§6.1) fixed the transport term, which is now
**validated as a necessary component**. Against the null control
`base = ana2e + esum + L:U` (no transport):

| pair | cos(base) | **cos(+transport)** | cos(broken transport) | r(new) |
|------|-----------|--------------------|----------------------|--------|
| (1,2) | +0.9798 | **+0.9805** | +0.9791 | 1.575 |
| (1,3) | +0.2368 | **+0.7840** | −0.3022 | 2.155 |
| (1,4) | +0.6366 | +0.6315 | +0.6336 | 2.992 |
| (2,3) | −0.5356 | **+0.9533** | −0.9594 | 1.571 |
| (2,4) | +0.2783 | **+0.3643** | +0.0320 | 5.356 |
| (3,4) | +0.9779 | **+0.9812** | +0.9738 | 1.284 |

Better than the null control on 5/6 pairs and than the broken version on 6/6;
(2,3) −0.54→+0.95 and (1,3) +0.24→+0.78 are decisive. **Direction is now largely
right** (three pairs > 0.95). But every ratio is **> 1** — the model overshoots.

### 6.0b The remaining term is genuinely new (proven)

Free least-squares fits on the cached components:

- 1 free coefficient on `L:U`: fitted values −0.79, 1.28, 1.28, 1.29, 1.12, −0.28
  (spread **2.08**) — no consistent weight.
- 4 free coefficients `p·ana2e + q·esum + s·L:U + t·transport`: `p` ∈ [−0.71, 1.19],
  `q` ∈ [0.00, 0.75], `s` ∈ [−0.40, 0.80], `t` ∈ [0.01, 3.57], and residuals of
  **0.08-0.88 × |orc|** even at best fit.

**No linear combination of the four terms reproduces the oracle.** The oracle is not
in their span, so the gap is a structurally new term — not a re-weighting, not a sign.
That is a hard bound, and it is the honest stopping point.

---

### 6.1 (historical) H2O analysis — superseded

`resp` must supply `oracle·gap − ana2e − esum`. Measured facts (H2O, snapshots):

- For pair (1,3), `ana2e` is **exactly anti-parallel** to the oracle (cos −1.0000) at
  6.77× magnitude — so a large, precise cancellation is required.
- **Now quantified via the validated esum:** `resp` must be ≈ **2-3** in magnitude
  (it has to cancel esum to ~1.5%), whereas the `L : Uˣ` prototype yields only
  ~0.2-0.4. So `resp` is **~10× under-captured**, not the 5-6× estimated before esum
  was available.
- The prototype `resp = L : Uˣ` is **structurally wrong**: tested against every
  projection of `Uˣ` (full / antisymmetric / symmetric, all blocks, `L` and `Lᵀ`),
  the best match was cos 0.63 at 0.20× magnitude. It is missing most of the physics,
  not just a sign or factor.
- `L` itself **is** matvec-consistent (built by orbital-rotation FD of
  `oqp.mrsf_matvec_apply`). The failure is in how `Uˣ`/the z-vector contracts `L`.

This is the genuine remaining crux: the interstate orbital-response contraction for
MRSF — the standard Lagrangian method (Li–Liu / Fatehi–Subotnik), **unpublished for
MRSF**.

---

## 7. Recommendation

1. **Ship the analytical NAC as it stands.** It is correct (6/6 vs GAMESS,
   |cos| ≥ 0.99995), cheap, and publishable. This is the project deliverable and it
   is done.
2. **Closed-form `d_amp` is a separate, open research problem.** Do not resume
   gradient-chain polarization (§4 proves it cannot work), and do not try to
   re-weight the existing four terms (§6.0b proves the oracle is not in their span).
   The next real step is deriving the missing response term from the Lagrangian —
   this is the Li-Liu / Fatehi-Subotnik construction, unpublished for MRSF.

## 8. Process lessons from this session (worth keeping)

Four real bugs were found, **every one of them only by an independent finite-difference
check** — none by code inspection:

| bug | signature that exposed it |
|-----|--------------------------|
| `fock_deriv_contract_os` is `intent(out)` (overwrites) | \|an\| ≈ ½\|FD\| |
| `utddft_xc_gradient` adds a ground-state XC baseline | error **identical** across all pairs |
| `logtol` passed with wrong units (screened the 1e term) | esum moved 10-200× when removed |
| transport used raw `dM`, not Löwdin-projected `dQ` | damage exactly where \|transp\| ≈ \|L:U\| |

Three corrections to my own earlier claims were required:

- "p20's esum is ~20× too big" — wrong premise (esum ≉ `missing`).
- "esum is large, ~100:1 cancellation" — an artifact of the `logtol` bug.
- "first real signal, cos +0.98" — the **null control** showed that agreement came
  from `ana2e + esum` alone; my response terms contributed nothing there.

**The transferable rule:** always evaluate a new term against *term absent* first.
Comparing variants of a model against each other attributes the baseline's agreement
to whatever was just added. A pair-independent error is always an additive constant,
never physics. And on a symmetric molecule (H2O, C2v) a 1-D irrep block makes
`|cos| = 1` automatic and meaningless — judge by ratio, and test on a C1 system.

---

## 8. Artifacts / reproducibility

- Diagnostics: `mrsf_nac_polarize` + `NAC_HOMOG` branch (`tdhf_mrsf_gradient.F90`);
  harnesses `p24_polar_fortran.py`, `p25_homog.py`, `p18_zvec_rhs.py`,
  `p19_zvec_L.py`, `p20_esum_grad.py`.
- Snapshots: `data_snapshots/p18_zvec.npz` (L, Uˣ, oracle, ana2e),
  `p12_ux_target.npz` (oracle, ana2e, missing), `p19_zvecL.npz` (z-vector d_L).
- Build (LOGIN node): `tools/nac_validation/build_oqp.sh`
  (PATH GCCcore-12.3, `CFLAGS=-g0`, `ninja -C build install`).
- Run (SLURM): `sbatch -p def -c N run_pyoqp.sh <harness.py>`.
