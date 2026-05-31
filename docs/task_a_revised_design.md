# Task A — Native Rys single-center Coulomb second derivative (l <= 3)

Status: **Gate 2 implemented and validated** (see Section 9). The production
basis-basis second derivative uses **angular-momentum (AM) shift identities** and
therefore does **not** differentiate the Rys roots/weights. The Gate-1
`rys_deriv.F90` root/weight X-derivative layer is retained as a validated
auxiliary/diagnostic layer (and as the basis for any future *direct*
charge-center differentiated-quadrature formulation), but it is **not** on the
production `p_AA / p_AB / p_BB` path. Sections 1-2 (the exact `nroots_der2`
formula and the s/p/d/f <= 5 root proof) remain in force; the
root-X-derivative-based assembly originally sketched in Sections 3-4 is
superseded by the AM-shift formulation for the basis-basis blocks.

## 0. Superseded baseline (do not use)

The libint `int1e_ipipnuc` / `int1e_ipnucip` baseline is **invalid** and is
withdrawn. Two settled reasons:
1. `grd1.F90::hess_en` scatters charge-center blocks 5-9 **inside** `do ic = 1,
   nat`, to Hessian atom-blocks indexed by the nucleus `ic`. It therefore needs
   **per-nucleus** single-center second derivatives `p_AA(C)`, `p_AC(C)`,
   `p_BC(C)`; the summed `ipipnuc`/`ipnucip` cannot provide them.
2. OpenQP libint has **only** ERI bindings (`libint2_build_eri/2eri/3eri`); there
   are no 1-body libint bindings. The native single-center Coulomb
   (`comp_coulomb_int1_prim` / `QGaussRys(c, znuc)`) already evaluates `1/|r-c|`
   per nucleus, so the native path is the correct vehicle.

Task A is the native Rys root-derivative for the per-nucleus operator `-Z_C/|r-C|`.

---

## 1. Exact nroots formula for the 1e nuclear-attraction second derivative

**Conventions (from code).** `basis%am` is the angular momentum L with `s=0,
p=1, d=2, f=3` (`basis_tools.F90`); `shi%ang = basis%am(shid) = L`. The
shell-pair sets, for the *energy/gradient* path
(`mod_shell_tools.F90:170,253`):

```
nroots_energy = (Li + Lj + 1)/2 + 1     (Fortran integer division)
```

A Gauss-Rys rule with `n` roots integrates the Rys-variable polynomial exactly to
degree `2n-1`. A nuclear-attraction integral of total angular momentum
`lambda = Li + Lj` needs `n_min = floor(lambda/2) + 1`. OpenQP's `+1` makes
`nroots_energy` cover **one** derivative (lambda -> lambda+1) for free, which is
why the existing gradient path reuses it.

**Second derivative.** Two derivatives raise the total angular momentum by 2
(either `d2/dA2`: bra +2, or `d2/dA dB`: bra +1 and ket +1). So
`lambda_eff = Li + Lj + 2`, and the minimal requirement is

```
nroots_der2 = floor((Li + Lj + 2)/2) + 1                       [Eq. 1]
```

**Consequence (must be enforced in code):** `nroots_der2 > cp%nroots` in general,
so the second-derivative routine **must set the root count from Eq. 1** and not
inherit the stored `cp%nroots`. (The current WIP `comp_coulomb_der2_braC` reuses
`cp%nroots`; for high L that **under-integrates** -- a real defect of the WIP,
independent of the fixed-root error.) The bra angular-momentum buildup is already
handled: `QGaussRys(..., igrd=2)` extends the xyz tables to bra index `iang+2`
(`mod_1e_primitives.F90:1854`, the `igrd1` path).

---

## 2. Proof that all s/p/d/f shell pairs need nroots <= 5

Apply Eq. 1 over `Li, Lj in {0,1,2,3}` (only `Li+Lj` matters):

| Li+Lj | example pair      | lambda_eff = Li+Lj+2 | nroots_der2 = floor(lambda_eff/2)+1 | root routine |
|-------|-------------------|----------------------|-------------------------------------|--------------|
| 0     | (s,s)             | 2                    | 2                                   | rys_rt2      |
| 1     | (s,p)             | 3                    | 2                                   | rys_rt2      |
| 2     | (s,d),(p,p)       | 4                    | 3                                   | rys_rt3      |
| 3     | (s,f),(p,d)       | 5                    | 3                                   | rys_rt3      |
| 4     | (p,f),(d,d)       | 6                    | 4                                   | rys_rt4      |
| 5     | (d,f)             | 7                    | 4                                   | rys_rt4      |
| 6     | (f,f)             | 8                    | **5**                               | **rys_rt5**  |

Maximum is **nroots = 5** at the (f,f) pair. Every s/p/d/f shell pair therefore
uses a **closed-form** root routine (`rys_rt1..rys_rt5`); none reaches
`rys_general` (nroots >= 6, the Golub-Welsch eigenvalue path). `MAX_NROOTS = 7`
and `mxrys = 13` leave ample array headroom. (Cross-check: OpenQP's conservative
`(Li+Lj+3)/2+1` style also caps at 5 for `Li+Lj=6`.)

This is the technical foundation of Option 1: the hard ingredient -- root/weight
X-derivatives -- is needed only for the five closed-form routines, which are
analytically differentiable by hand.

---

## 3. Root/weight X-derivative validation plan (vs FD of rys_rt1..rys_rt5)

The new analytic derivatives `du_t/dX`, `dw_t/dX` (and the second derivatives
`d2u_t/dX2`, `d2w_t/dX2` needed for the full Hessian term) are obtained by
differentiating the **existing** branch expressions of each `rys_rtN` in place
(same Chebyshev/polynomial + `exp(-X)` branches; the branch boundaries
`X <= 3e-7, 1, 3, 5, ...` and the asymptotic large-X forms are inherited
unchanged). They are validated **before** any integral assembly, against central
finite differences of the unmodified `rys_rtN`:

- **Targets:** for each `N = 1..5`, every root `u_t(X)` and weight `w_t(X)`.
- **Reference:** `(rys_rtN(X+h) - rys_rtN(X-h)) / (2h)` for first derivatives;
  the standard 2nd-difference stencil for `d2/dX2`. Non-iterative ->
  truncation-limited.
- **Sampling:** >= 3 X points strictly inside **every** branch interval, plus
  points bracketing each branch boundary (to confirm the derivative is correct on
  both sides; the value is continuous, the derivative branch must match), plus
  large-X asymptotic samples.
- **Step protocol:** central differences at `h in {1e-3, 5e-4, 2.5e-4}` (in X);
  report the error ratio, require ~4.0 (O(h^2)); identify the round-off floor.
- **Pass/fail:** max-abs `< 1e-9` in the truncation-limited regime, O(h^2) ratio
  confirmed. A flat plateau or non-O(h^2) ratio fails (this is precisely the test
  the fixed-root scheme never had; it is the primary guard against repeating it).
- **Boundary stress:** explicitly test X within `1e-10` of each branch switch;
  if a derivative is discontinuous there, the branch expressions disagree and the
  term is rejected.

Only after the root-derivative unit tests pass does work proceed to Section 4.

---

## 4. Revised per-nucleus integral equations

For a shell pair (bra shell on atom A, ket shell on atom B) and a single nucleus
C with charge `Z_C`, the operator is `v_C(r) = -Z_C / |r - C|`. The Rys integral
in 1D (per Cartesian axis, per root t) is `I_t(n_i, n_j)` built by `QGaussRys`
with argument `X = pp%aa * sum((pp%r - C)**2)`, `pp%r = (a_i A + a_j B)/(a_i+a_j)`.
Define the elementary, exact geometric factors (X is quadratic in the centers):

```
dX/dA = 2 pp%aa (a_i/(a_i+a_j)) (pp%r - C)
dX/dB = 2 pp%aa (a_j/(a_i+a_j)) (pp%r - C)
dX/dC = -2 pp%aa (pp%r - C)
d2X/dA dB, d2X/dA2, ... : constant tensors from the same quadratic form
```

The **full** second derivative of the contracted integral w.r.t. two centers
P, Q in {A, B} is (schematically, summed over roots t and the polynomial xyz
factorization):

```
d2 I / dP dQ = sum_t [  (d2_poly/dP dQ) w_t                        (poly-poly, frozen root: the only term the fixed-root scheme kept)
                      + (d_poly/dP)(dw_t/dX)(dX/dQ) + (P<->Q)      (poly-weight cross)
                      + poly ( (d2w_t/dX2)(dX/dP)(dX/dQ)
                             + (dw_t/dX)(d2X/dP dQ) )              (weight curvature + geometric)
                      + (identical structure with root u_t and du_t/dX,
                         since the polynomial recurrence coefficients depend on u_t) ]
```

`d_poly/dP` etc. are the existing angular-momentum shifts (`der_coul_xyz` and the
`igrd=2` tables); `du_t/dX, dw_t/dX, d2u_t/dX2, d2w_t/dX2` are the Section-3
quantities. The four blocks the assembly needs:

- **p_AA(C)** = `d2 I / dA2` (both derivatives on the bra center A).
- **p_AB(C)** = `d2 I / dA dB` (one on bra A, one on ket B).
- **p_BB(C)** = `d2 I / dB2` (from the swapped shell pair, i.e. the same routine
  with bra<->ket).
- **p_AC(C)**, **p_BC(C)** (bra-charge, ket-charge): obtained from translational
  invariance of the single-center operator
  `(d/dA + d/dB + d/dC) v_C = 0`:

```
p_AC(C) = -(p_AA(C) + p_AB(C))            [Eq. 2a]
p_BC(C) = -(p_AB(C)^T + p_BB(C))          [Eq. 2b]
p_CC(C) =  p_AA(C) + p_AB(C) + p_AB(C)^T + p_BB(C)   [Eq. 2c]
```

These are exact for the **single-center** operator (Section 6 of the prior review
proved no extra terms arise; `v_C` depends only on `A-C` and `B-C`). This is why
the per-nucleus operator is mandatory: TI holds per nucleus, not for the summed
`nuc`.

Each `p_AA(C)`, `p_AB(C)` is contracted with the density block over the shell-pair
AOs (matching `bfnrm` normalization and the `-Z_C` charge factor, exactly as
`comp_coulomb_der1` / `comp_coulomb_helfeyder1` do for the gradient) before the
3x3 block is returned to `hess_en`. The contraction convention must equal the one
the validated gradient routines use (verified by test 5 below).

**Note on the WIP factorization.** The current `comp_coulomb_der2_braC` uses a
bra/Hellmann-Feynman split that was shaped around the (wrong) fixed-root scheme.
Per open question, the corrected `p_AA`/`p_AB` must be **re-derived** from
Eq. (the d2 I expansion) with the root-derivative terms, not patched onto the WIP;
`p_AC`/`p_BC` then follow from Eq. 2. Whether the bra/charge split survives is an
implementation decision to make from the derivation, not from the existing code.

---

## 5. Low-symmetry PySCF oracle protocol (with_rinv_at_nucleus)

The integral-level oracle is the **primary** correctness gate (it caught the
fixed-root error last time). It must run on a **low-symmetry, all-distinct-atom**
geometry so that no charge-center cross block is zero or degenerate by symmetry.

- **Molecule:** a C1 system with all atoms inequivalent, e.g. **HOF** (H-O-F bent,
  no symmetry) or a deliberately distorted/rotated geometry; small basis
  (STO-3G or 6-31G) with shells up to the L being tested. H2 / H2O are forbidden
  here -- they zero or degenerate the very blocks under test.
- **Per nucleus C** (loop over all nuclei):
  ```
  with mol.with_rinv_at_nucleus(C_index):
      ipiprinv = mol.intor('int1e_ipiprinv', comp=9).reshape(3,3,nao,nao)  # d2/dA2
      iprinvip = mol.intor('int1e_iprinvip', comp=9).reshape(3,3,nao,nao)  # d2/dA dB
  p_AA_ref(C) = -Z_C * ipiprinv ;  p_AB_ref(C) = -Z_C * iprinvip
  ```
  (The `-Z_C` charge factor is applied explicitly; the bare intors are chargeless.)
- **AO matching:** build an **explicit AO-permutation + normalization map** from
  the PySCF `cart=True` Cartesian ordering to OpenQP's shell/`bfnrm` layout. No
  Frobenius-norm or trace shortcuts -- compare **element by element** after the
  map. (A consistent permutation error otherwise passes silently.)
- **Component check:** confirm the `(a_bra, b_ket)` index assignment by comparing
  `p_AB_ref[a,b]` vs `p_AB_ref[b,a]` for `a != b` to the OpenQP block; this
  catches a component transpose that the (a,b)-symmetric `ipipnuc` would hide.
- **Mixed-block cross-check (breaks TI circularity):** compute `p_AC` two ways --
  (i) Eq. 2a from OpenQP `p_AA,p_AB`, and (ii) directly from PySCF's per-nucleus
  `hcore_generator`-style charge-center block -- and require agreement `< 1e-9`.
- **Pass/fail:** every element max-abs `< 1e-9` (cart, normalized), on the
  low-symmetry molecule, for **each** nucleus separately.

Contracted-level confirmation (secondary, mandatory): the `hess1_selftest`
"nucattr" line -- analytic `hess_en` vs central FD of
`grad_en_pulay + grad_en_hellman_feynman` via the separated-FD harness -- on the
same low-symmetry molecule, h in {8e-3,4e-3,2e-3,1e-3} bohr, O(h^2) ratio ~4.0,
abs `< 1e-6`. This catches charge-factor and density-convention errors. Plus the
**subsystem regression guard** (re-run `fockx`/`cphf_f0x`/`cphf_dpdx`; any change
fails Task A).

---

## 6. High-L policy (g and above: clean abort)

Shells with `L >= 4` (g, h, ...) on either center can require `nroots >= 5`
already at (f,f)=5 but exceed the closed-form regime for `Li+Lj >= 7`
(nroots `>= 5`; `Li+Lj = 7` -> nroots 5 still closed-form; `Li+Lj >= 8`
-> nroots `>= 6` -> `rys_general`). Since the analytic root derivatives are
implemented only for `rys_rt1..rys_rt5`, the native Hessian must **refuse** any
shell pair that would need `nroots >= 6`:

- **Guard location:** at the entry of the native nuclear-attraction Hessian
  routine, compute `nroots_der2` (Eq. 1) per shell pair; if `nroots_der2 > 5`
  (equivalently `Li + Lj >= 8`, i.e. any g-or-higher combination beyond the
  closed-form reach), **abort cleanly** with a specific diagnostic:
  `"Native analytic Hessian: nuclear-attraction 2nd derivative for shells with
   L>=4 (nroots>=6, rys_general) is not implemented; use a basis with l<=f or the
   PySCF bridge / numerical Hessian."`
- **Input-checker mirror:** `pyoqp/oqp/utils/input_checker.py` rejects
  `[hess] type=analytical` (native) when the basis contains `L >= 4` functions,
  with the same message, so the failure is reported before the SCF rather than
  mid-integral.
- **No silent fallback:** consistent with the data contract; the abort is
  explicit. (The PySCF bridge remains available for l >= 4 as a separate,
  user-selected path.)
- This bounds Task A to the analytically-tractable closed-form regime and defers
  the Golub-Welsch root-derivative (or FD-root) work to a clearly separated,
  later task.

---

## 7. Files

### To touch
- `source/integrals/rys.F90` -- add analytic X-derivative routines for
  `rys_rt1..rys_rt5` (new, e.g. `rys_rt1_d`, ... returning `du/dX, dw/dX` and the
  2nd derivatives). **Do not modify** the existing `rys_rtN` value routines.
- `source/integrals/mod_1e_primitives.F90` -- replace the fixed-root
  `der2_coul_xyz` and re-derive `comp_coulomb_der2_braC` to assemble the full
  second derivative using the root-derivative terms (Section 4); set
  `nroots_der2` from Eq. 1 (do not inherit `cp%nroots`); keep the per-nucleus
  `(c, znuc)` interface that `hess_en` calls.
- `source/integrals/grd1.F90` -- `hess_en`: add the `nroots_der2 > 5` /
  `L >= 4` guard (Section 6); update the WIP `@warning` once validated. The
  `do ic` loop and the 9-block scatter are unchanged in structure but must be
  **re-validated** end-to-end against the correct per-nucleus integrals on the
  low-symmetry molecule (the scatter is currently provisional, never run against
  correct inputs).
- `source/modules/hess1_selftest.F90` -- re-enable the "nucattr" PASS/FAIL line;
  add the per-nucleus integral oracle, component-transpose, and mixed-block
  cross-check tests driven from a low-symmetry geometry.
- `pyoqp/oqp/utils/input_checker.py` -- L>=4 rejection for the native analytic
  Hessian (Section 6).
- `docs/task_a_revised_design.md`, `docs/second_derivative_integral_design.md` --
  status updates.

### Forbidden to touch (read-only; frozen validated response subsystem)
- `source/modules/fock_deriv.F90`, `source/modules/cphf.F90`
- `grd1.F90` first-derivative *matrix* builders `der_overlap_matrix`,
  `der_kinetic_matrix`, `der_nucattr_matrix`, and the `comp_*_der1_block` helpers
- `source/modules/fock_deriv_selftest.F90`, `cphf_nuclear_selftest.F90`,
  `cphf_dpdx_selftest.F90`
- the validated overlap/kinetic 2nd-derivative path: `der2_kinovl_xyz`,
  `comp_overlap_der2`, `comp_kinetic_der2`, `hess_ee_overlap`, `hess_ee_kinetic`
- the existing `rys_rtN` value routines and `rys_general` (read for
  differentiation; not modified)
- `source/modules/hf_hessian.F90` -- stays guarded; Task A does not unguard it
  (final assembly, after Tasks A+B, is a separate step).
- all libint files (`libint_f.F90`, `int_libint.F90`) -- Option 1 needs no libint.

---

## 8. Open questions to resolve before coding
1. Confirm the exact `d2/dA2` vs `d2/dA dB` angular-momentum raising in the
   existing xyz factorization, so the polynomial-shift factors multiplying
   `du/dX, dw/dX` are placed correctly (re-derive; do not infer from the WIP).
2. Decide whether `comp_coulomb_der2_braC`'s bra/charge split is retained or
   replaced after the root-derivative re-derivation.
3. Confirm `nroots_der2` (Eq. 1) and the `xyzin` array sizing for the +2 bra
   angular momentum at the largest closed-form case (f,f, nroots=5).
4. Re-validate the hess_en 9-block scatter against correct per-nucleus integrals
   on the low-symmetry molecule before removing its "provisional" status.

---

## 9. Gate 2 — implemented formulation and validation (AM-shift)

**Decision (supersedes the "Gate 2 must use rys_deriv" framing).** The production
basis-basis second derivative is built from angular-momentum raising/lowering
identities, *not* from differentiated Rys roots/weights:

- `d/dA chi_A`, `d2/dA2 chi_A`, `d/dB chi_B`, `d2/dB2 chi_B` are exact linear
  combinations of ordinary single-center Coulomb integrals with shifted angular
  momentum. The Rys roots/weights are simply **recomputed for the shifted
  integral class** at the corrected count
  `nroots_der2 = floor((Li + Lj + 2)/2) + 1` (explicitly set, **not** inherited
  from `cp%nroots`).
- **The production basis-basis second derivative uses AM-shift identities and
  therefore does not differentiate Rys roots/weights.** `rys_deriv.F90` stays as
  a validated auxiliary/diagnostic layer (and a basis for a future direct
  charge-center differentiated-quadrature route); it is not used here.

**Code.** `mod_1e_primitives.F90`:
- `comp_coulomb_der2_blocks(cp, c, znuc, pAA, pAB)` — shared kernel returning the
  uncontracted per-AO blocks `pAA(a,b,i,j) = d2/dA_a dA_b <i|znuc/|r-c||j>` and
  `pAB(a,b,i,j) = d2/dA_a dB_b <i|...|j>` (bra-bra via `der2_coul_xyz`; ket and
  mixed via the explicit index recurrences). Aborts cleanly for `nroots_der2 > 5`
  (any g-or-higher pair).
- `comp_coulomb_der2_braC(...)` — thin density-contraction wrapper over the kernel
  returning `p_XX` and `p_XC = -(p_XX + p_AB)`.

Per-nucleus blocks via single-center translational invariance (after the direct
basis-basis blocks are correct):
```
p_AC(C) = -(p_AA(C) + p_AB(C))
p_BC(C) = -(p_BB(C) + p_AB(C)^T)
p_CC(C) = -(p_AC(C)^T + p_BC(C)^T) = p_AA + p_AB + p_AB^T + p_BB
```

**Validation arbiter = the per-nucleus integral oracle** (not whether
`rys_deriv.F90` is used, and not TI/Hessian-symmetry, which the TI-built blocks
satisfy by construction). `tests/test_hess_nuc_oracle.py`, on a low-symmetry C1
all-distinct-atom molecule (HOF / 6-31G*), checks **element-wise** (no
Frobenius-norm shortcuts):
- `p_AA(C)` vs PySCF `int1e_ipiprinv` at nucleus C (sign/component pinned by FD:
  `int1e_ipiprinv[a,b] = +d2/dA_a dA_b`);
- `p_AB(C)` vs PySCF `int1e_iprinvip` (`[a,b] = +d2/dA_a dB_b`), including
  explicit `a != b` component-order checks;
- `p_BB(C)` is covered by the full `nbf x nbf` sweep (bra on every atom);
- an **explicit AO permutation** (proven on the overlap correlation matrix) plus
  an **explicit per-AO normalization scale** (OpenQP `xx,yy,zz,xy,xz,yz` vs PySCF
  `xx,xy,xz,yy,yz,zz` for d; libcint cart order generally);
- an explicit `-Z_C` charge-factor check.
Pass tolerance `< 1e-9` per element, per nucleus.

**Secondary (contracted) confirmation**, run only after the integral oracle
passes: (a) analytic `hess_en` nuclear-attraction total vs central FD of
`grad_en_pulay + grad_en_hellman_feynman` on the same molecule, required to be at
parity with the validated overlap/kinetic second-derivative FD errors; (b) a
**non-circular** charge-charge check — the TI-built `p_CC` (hess_en block 9) vs an
*independent* FD of the HF charge gradient w.r.t. the charge position (`fd_cc`),
which never uses TI.

**Status:** integral oracle + both secondary confirmations pass; Gate-1
`rys_deriv` regression and the H2O overlap/kinetic selftest remain green. The
`hf_hessian` kernel stays guarded (final assembly is a later, separate step).
