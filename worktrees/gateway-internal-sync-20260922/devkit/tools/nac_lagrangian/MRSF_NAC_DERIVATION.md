# Fully Analytic MRSF-TDDFT Nonadiabatic Couplings — Derivation

**Final production statement (2026-08-01).**  OpenQP evaluates the static
singlet MRSF-TDA interstate derivative coupling

```
d_IJ = <Psi_I | d/dx | Psi_J>,      h_IJ = (E_J - E_I) d_IJ
```

with one **adjoint Z-vector RHS per ordered state pair**, not `3N` forward
CPHF right-hand sides.  Every scientific operation in the production path is
resident Fortran: the exact state-overlap metric, amplitude/Fock bilinear,
ROHF/ROKS orbital source, one-RHS MINRES solve, HF/JK/Pulay and XC adjoint
contractions, ordered-pair accumulation, HST projection, and gap scaling.
`pyoqp/oqp/library/nac_analytic.py` validates the public scope, invokes the
single `mrsf_nac_lagrangian` entry point, and reshapes the final records only.
Nuclear and orbital-generator finite differences remain independent
diagnostics and are not part of production.

Code/legacy compact conventions in Secs. 0--7: `x` = one nuclear
coordinate; `p,q,r,s` MO; `mu,nu` AO; `i,j` occupied-alpha
(doc+socc), `a,b` virtual-beta (socc+virt).  When the Lee et al. paper
is quoted below its notation is used instead: `i,j in C` (doubly
occupied), `x,y in O` (the two singly occupied orbitals), and `a,b in
V` (virtual).  In all sections, `d_IJ = -d_JI`, `h_IJ = +h_JI`, and
`gap = E_J - E_I`.  Symbols called an "ordered kernel" below denote a
one-sided HST derivative before the state-index projection.  Only the returned
observable is constrained by `d_IJ=-d_JI` and `h_IJ=h_JI`.

This document is also a chronological forensic log.  Sections 0--7.67 retain
failed hypotheses, intermediate implementations, and old accuracy snapshots.
They are deliberately preserved as an audit trail, but are **non-normative**.
The published anchor immediately below and Section 8 onward are the final
specification whenever an older statement conflicts with them.

## Published energy-gradient anchor and its exact scope

The primary reference is S. Lee, E. E. Kim, H. Nakata, S. Lee, and
C. H. Choi, *J. Chem. Phys.* **150**, 184111 (2019),
doi:[10.1063/1.5086895](https://doi.org/10.1063/1.5086895),
"Efficient implementations of analytic energy gradient for
mixed-reference spin-flip time-dependent density functional theory
(MRSF-TDDFT)."  The paper and its supplementary derivation were checked
against the authors' local PDF/TeX, not reconstructed from code comments.

The paper fixes the diagonal (`I=J`) Lagrangian normalization that the
implementation must reproduce:

```
Eq. (3.3):
L^k = G^k
    + 2 sum_{ia,s} Z^k_{ia,s} Fbar_{ia,s}
    +   sum_{ix,s} Z^k_{ix,s} Fbar_{ix,s}
    +   sum_{xa,s} Z^k_{xa,s} Fbar_{xa,s}
    - 2 sum_{p<=q,s} W^k_{pq,s}(S_{pq,s}-delta_pq)

Eq. (3.6):       Jbar Zbar = -Rbar
Eq. (3.7):       Zbar_ix = Z_ix,beta
                  Zbar_xa = Z_xa,alpha
                  Zbar_ia = Z_ia,alpha = Z_ia,beta
```

Equations (3.8a--f) give the ROHF orbital Hessian `Jbar`.  The paper
states explicitly that this left-hand side is identical to the existing
SF-TDDFT one.  The MRSF-specific state dependence is on the right:
Eq. (3.9) gives `Rbar` in terms of `H+[T]` and `H[X,X]`, Eq. (3.10)
defines `H+`, and Eqs. (3.11--3.15) define `T`, the expanded `U X`, and
the spin-pairing additions.  After the solve, Eq. (3.16) is `P=T+Z`;
Eq. (3.18) determines the energy-weighted density `W`; and Eq. (3.21)
contracts the relaxed one-particle density, overlap multiplier, and
two-particle density with nuclear integral derivatives:

```
Omega^R = sum h^R P - sum S^R W + sum (mu nu|kappa lambda)^R Gamma.
```

The supplementary derivation fixes the same normalization from the
coefficient-stationarity side.  SI Eqs. (S25a--i) are the paired equations for
the energy-weighted multiplier `Wbar` after inserting `Zbar`.  Subtracting
S25a from S25b, S25c from S25d, and S25e from S25f removes `Wbar` and gives SI
Eqs. (S26a--c).  Their left sides define `Jbar Zbar`, their right sides define
`-Rbar`, and substituting the solution back into S25 determines `Wbar`.  Thus
S25/S26 independently lock the C--O, C--V, and O--V half factors; they are not
license to add a second factor after converting to OpenQP's native tangent and
dual coordinates.

The corresponding OpenQP gradient map is:

```
mrsfcbc / mrsfmntoia   <-> U expansion and U^T projection, Eqs. (2.5),(2.8)
SPC channel scaling    <-> spin pairing, Eqs. (2.13)--(2.15)
mrsfxvec               <-> expanded U X used in T, Eq. (3.11)
mrsfsp                 <-> pairing part of H[X,X], Eqs. (3.12)--(3.15)
sfrorhs                <-> -2 Rbar in the code's internal normalization
sfrolhs                <-> Jbar action, Eqs. (3.6),(3.8), SI Eq. (S26)
sfropcal               <-> xk/2 -> the paper's one-sided Zbar density
mrsfrowcal             <-> W, Eqs. (3.18)--(3.19)
```

The final interstate driver reuses these validated kernels but does not route
its computational adjoint through the legacy `sfropcal` half-density seam.
The exact sign and coordinate conversion are derived in Section 8.5.

There are therefore **two different Hessians** in this code.  The
MRSF Davidson Hessian is the response-state operator assembled through
`mrsfcbc -> int2 -> SPC -> mrsfmntoia`.  The orbital-stationarity
Hessian is the SF/ROHF `Jbar` used by the Z-vector equation.  They must
not be identified with each other.

One typesetting defect in the published paper must not be copied:
Eq. (3.8f) prints `-1/c_H` in the singly-occupied--virtual block.
Equation (3.10), SI Eq. (S26), direct differentiation, and the working
implementation all give `-c_H/2` for that exchange pair.

Most importantly, this is an **energy-gradient paper**.  It proves the
state-diagonal Lagrangian and provides a mandatory `I=J` limit, but it
does not derive an interstate NAC Lagrangian.  In particular it does
not prove the off-diagonal polarization/order of every quadratic
`T[X,X]`, `H[X,X]`, pairing `Gamma`, or `W` term.  The interstate
bilinearization in the remainder of this document is a new derivation,
validated numerically where stated; it must not be described as a
theorem of Lee et al.

---

## Historical archive boundary

Everything from Section 0 through Section 7.67 records the route to the final
formula.  Terms such as "current", "production", "remaining", and "final" in
that archive describe the campaign stage at which they were written.  In
particular, the raw-determinant metric, sign-scanned TLF kernel, same-space
projection, antisymmetric-`dD` fold hypothesis, Python pair algebra, forward
`3N` CPHF interpretation, and irreducible XC-grid-floor claims were all
superseded.  Section 8 is the only normative interstate derivation.

## 0. Representation choice (the lesson of Phase 11)

Work **exclusively in the unfolded determinant space**. Let `V` be the
constant orthogonal fold matrix (`mrsfxvec`): for the singlet sector it maps
the two SOMO ground-pair configurations onto `(lr1 - lr2)/sqrt(2)` (triplet:
`+`), identity elsewhere. Define once, at the reference geometry:

```
Xt_I = V X_I        (theta-independent, frozen tensors)
At   = V A V^T      (symmetric; same spectrum, eigenvectors Xt_I)
```

All derivatives below act on `At`'s **integrals only**. The spin-adaptation
factors (1/sqrt2 on the ground pair, sqrt2 in `get_trans_den`'s SOMO cross
terms) are entries of the constant `V` — they are never differentiated.

*Why this matters:* Phase 11 proved (`PHASE11_fd_findings.md`) that the
production matvec is the composite `fold ∘ A_det ∘ unfold` with the fold
applied on the *output* side; its orbital-rotation derivative picks up the
derivative **of the fold itself**, which is multiplied by the huge ground-
pair amplitude (x_lr ≈ 0.998) and swamps every alpha rotation block. In the
unfolded representation that contamination does not exist by construction.
The production z-vector RHS (`sfrorhs`) already lives on the input-fold
side, which is why it is the structurally consistent partner.

---

## 1. Exact splitting of d_IJ

The MRSF state is `|Psi_J> = sum_m Xt^J_m |Phi_m(C(x))>` over spin-flip
determinants `Phi_m = Phi_i^a` (alpha-occ i → beta-virt a on the ROHF
triplet reference). Differentiate the ket:

```
d_IJ = sum_m Xt^I_m (d Xt^J_m/dx)                       [amplitude part]
     + sum_mn Xt^I_m Xt^J_n <Phi_m | d/dx Phi_n>        [orbital part]
     =: d_amp + d_orb
```

This split is exact and mirrors the branch's `d_IJ = d_ov - d_amp` (their
sign absorbed into the definition below).

---

## 2. The orbital part d_orb — orbital and state antisymmetries

### 2.1 The one-particle reduction

`d/dx phi_p = sum_q U^x_qp phi_q + (moving-AO term)`, so

```
T_pq := <phi_p | d/dx phi_q> = U^x_pq + Sk_pq ,
Sk_pq := sum_{mu,nu} C_mu_p C_nu_q <chi_mu | d/dx chi_nu>     (ket-half overlap derivative)
```

Orthonormality `d/dx <phi_p|phi_q> = 0` gives **T_qp = -T_pq exactly**:
the full orbital-derivative matrix is antisymmetric. (Equivalently
`U^x_pq + U^x_qp = -S^x_pq` with `S^x = Sk + Sk^T`.) This fixes the
orbital-index contraction.  It does **not** require either one-sided ordered
HST kernel to be antisymmetric in the state labels.  The observable evaluated
by both the numerical HST referee and the analytic implementation is

```
d_IJ = 1/2 (Dord_IJ - Dord_JI).
```

That projection is part of the HST definition, not an independent theory
gate and not evidence that each ordered leg is complete.  In the H2O v64
debug artifact the unprojected symmetric parts are sizeable; only comparisons
of the projected observable with an independent numerical reference test the
assembled formula.

Slater–Condon on `<Phi_m| d/dx Phi_n>` (determinants differing by at most
one spin-orbital, and the derivative is a one-particle substitution) yields

```
d_orb = sum_{pq} gamma^IJ_pq T_pq = sum_{p>q} (gamma^IJ_pq - gamma^IJ_qp) T_pq
```

with the **interstate CIS-type transition density** (unfolded amplitudes):

```
occ-occ :  gamma^IJ_ij = - sum_a Xt^I_ia Xt^J_ja        (alpha spin)
vir-vir :  gamma^IJ_ab = + sum_i Xt^I_ia Xt^J_ib        (beta spin)
socc cross terms: generated automatically because the socc orbitals are
   BOTH alpha-occupied and beta-virtual; the sqrt2 factors of
   get_trans_den are the V-entries appearing here, fixed, not fitted.
```

Only the orbital-index antisymmetric part
`gamma_a := (gamma - gamma^T)/2` survives its contraction with `T`.  This is
distinct from the state-index HST projection
`(Dord_IJ-Dord_JI)/2`; neither operation may be substituted for the other.

**Correction to the branch:** same-space (doc-doc / socc-socc / virt-virt)
blocks of `gamma_a` are NOT automatically "pure gauge". The doc-doc and
virt-virt blocks of gamma are symmetric in (I<->J composite) only when the
amplitude products commute; the socc-socc block in particular carries the
ground-pair structure. Keep all blocks of `gamma_a`; let the algebra kill
what vanishes. (The branch zeroes them: `single_point.py:1802` — re-derive,
do not copy.)

### 2.2 Split T into solvable + explicit

```
d_orb = sum gamma_a . U^x            [-> goes to the Lagrangian, Sec. 4]
      + sum gamma_a . Sk             [explicit; ket-half overlap derivative]
```

The second term is the branch's `mrsf_nac_overlap` frozen term — correctly
built there with the real derivative integrals (`der_overlap_matrix_ket`,
`grd1.F90:694`, from `comp_overlap_der1_block`). Keep it, but contract with
the *derived* `gamma_a`, not the sign-scanned TLF kernel.

---

## 3. The amplitude part d_amp — first-order perturbation

`(At - Om_J) Xt_J = 0`. Differentiate, project on `Xt_I` (I ≠ J,
`Xt_I . Xt_J = 0`):

```
Xt_I^T (d Xt_J/dx) = Xt_I^T (dAt/dx) Xt_J / (Om_J - Om_I)
```

Since both states share the reference SCF energy, `Om_J - Om_I = E_J - E_I
= gap`. Now chain-rule the operator:

```
dAt/dx = (dAt/dx)|_skel                       [MO coefficients frozen; integral derivatives]
       + sum_{p>q} (dAt/dtheta_pq) U^x_pq     [orbital response]
```

so, defining the two ingredients

```
h_skel := Xt_I^T (dAt/dx)|_skel Xt_J
L_pq   := Xt_I^T (dAt/dtheta_pq) Xt_J        [interstate orbital gradient of the response operator]
```

```
d_amp = [ h_skel + sum_{p>q} L_pq U^x_pq ] / gap
```

### 3.1 h_skel — already built and FD-validated on the branch

`At = E + G`: `E` = the Fock-diagonal contraction (`mrsfesum`), `G` = the
7-channel 2e kernel (`mrsfmntoia`/`mrsfcbc` + spc scalings, sign `mrst`).

* 2e piece: the bilinear symmetrized density
  `D2_IJ = 1/2 [f_I(a) g_J(b) + f_J(a) g_I(b)]` contracted with real ERI
  derivatives through `grd2_driver` — this is exactly
  `grd2_mrsf_nac_compute_data_t` (`mrsf_nac_amp`), which collapses
  bit-for-bit to the production gradient at I=J (opt-in selftest exists).
* 1e/Fock piece: `Tr[P^IJ dF_skel/dx]` with
  `P^IJ = C Gam^IJ C^T`, `Gam^IJ` = the symmetrized interstate difference
  density (`mrsf_interstate_tden`) — this is `mrsf_nac_esum` (1e + 2e + XC),
  FD-validated to 1e-9..1e-5 on the branch. **Reuse as-is.**

### 3.2 L — analytic, in the unfolded space

`theta_pq` rotates MOs: `dC/dtheta_pq = C K_pq` (antisymmetric generator,
spin-resolved per ROHF block: doc-socc rotates beta only, doc-virt both,
socc-virt alpha only — `sfrogen`, confirmed in Phase 11 findings #4).

* From `E`: `dE/dtheta` inserts `[K, F_MO]` plus the frozen-density Fock
  response; contraction with `Xt_I ⊗ Xt_J` gives the `Gam.F + F.Gam`
  commutator structure — the matrix `M_A = Gam_A F_oo + F_oo Gam_A`,
  `M_B = Gam_B F_vv + F_vv Gam_B` already coded in the `nac_wsx` block
  (`tdhf_mrsf_gradient.F90:1372`ff).
* From `G`: `dG/dtheta` rotates one MO index of each channel kernel;
  bilinear contraction gives
  `L^2e_pq = [G_MO(Xt_I ⊗ Xt_J) + G_MO(Xt_J ⊗ Xt_I)]_{pq}` antisymmetrized
  — the same "input-fold analytic RHS" that Phase 11 verified equals
  `sfrorhs`'s hxa/hxb structure to 1e-15 at I=J (`fcac55a`, `90943f8`).
* Plus the derivative of the frozen AO Fock **through the density**: at
  fixed AO-Fock (the matvec's convention) this term is absent from At but
  reappears in the z-vector LHS; see 4.2 — this bookkeeping choice is why
  the naive `L : U^x` prototype came out at 0.2x magnitude: part of the
  response lives in the *Fock rebuild*, i.e. in `A^orb`, not in L.

Validation gate for L: rotation-FD of the **unfolded** operator
`Xt_I^T At(C e^{theta K}) Xt_J` per generator — free of the output-fold
artifact — must match analytic L to FD floor, per spin block.

---

## 4. The Lagrangian: one z-vector per pair

Collect every U^x dependence from Secs. 2.2 and 3:

```
d_IJ = [ h_skel + sum_{p>q} Ltot_pq U^x_pq ] / gap  +  sum gamma_a . Sk
with  Ltot_pq := L_pq + gap * (gamma_a)_pq
```

### 4.1 Eliminate U^x

Split U^x into the independent rotations (the three ROHF blocks) and the
dependent part fixed by orthonormality `U^x_sym = -1/2 S^x`:

```
sum_{p>q} Ltot_pq U^x_pq
   = sum_{ai in rot-space} Ltot_ai U^x_ai            [independent -> z-vector]
   - 1/2 sum_{pq} Ltot^sym-partner . S^x_pq          [dependent -> -Tr(W^IJ S^x)]
```

The second piece defines the **interstate energy-weighted density**
`W^IJ`; its Fock-commutator part is exactly the exported-but-unused
`OQP::nac_wsx` object. It is LARGE (carries core orbital energies ~ -20 Ha)
and is what cancels the large `esum` — the measured
`|resid| = |missing - esum| ≈ |esum|` to 1.5% on the branch says this
cancellation is real. **It must be added to the assembly** (the branch
computes it and then does not use it).

### 4.2 Interchange (Handy–Schaefer)

`U^x` solves the coupled-perturbed ROHF equations
`A^orb U^x = B^x` (B^x = the CPHF RHS built from skeleton derivatives).
Then

```
sum_ai Ltot_ai U^x_ai = sum_ai z^IJ_ai B^x_ai,     A^orb z^IJ = Ltot|_rot
```

ONE symmetric-solve per pair (GMRES/MINRES on the existing `sfrolhs` LHS —
untouched), then `z^IJ . B^x` is contracted for all 3N coordinates through
the production gradient code, exactly like a normal excited-state gradient.
The plumbing already exists: `OQP::nac_orbgrad_L` RHS hook at
`tdhf_mrsf_z_vector.F90:1800`, relaxed-density contraction via the
`set_mrsf_nac_cphf` seam.

Bookkeeping caveat (3.2): whether the frozen-Fock response sits in `Ltot`
or in `A^orb` must match how B^x is defined on the gradient side. Use the
gradient's own convention (sfrorhs/sfrolhs pairing) — the interchange is
only valid for a consistent (RHS, LHS, B^x) triple. This is the single
most likely place to lose a factor; gate it with the identity
`sum L U^x == z . B^x` checked numerically for ONE displaced geometry.

### 4.3 Final assembly

```
h_IJ  = h_skel(2e bilinear + esum)  +  z^IJ . B^x  -  Tr[ W^IJ S^x ]
d_IJ  = h_IJ / gap  +  Tr[ gamma_a Sk ]
                                     (+ optional ETF correction, Sec. 5)
```

No fitted coefficient appears: the orbital 1/2s and signs are fixed by
orthonormality and the perturbation identity.  Compute both ordered kernels
independently and apply the same HST projection as the numerical observable.
The resulting `d` antisymmetry and `h=gap*d` symmetry are output contracts;
their zero residuals alone are not independent validation evidence.

---

## 5. Conventions that must be locked before comparing to anything

1. **Gap orientation:** `gap = E_J - E_I`; `nacv[i,j] = +h_IJ` symmetric.
   Fix `numerical_nac` (axes swapped at `single_point.py:2273`) and the
   analytic stamp (`:2196`) FIRST, freeze oracle references WITH SIGN.
2. **HST factor:** `(S - S^T)/(2 dt)` + central average — already correct.
3. **ETF:** the FD oracle computes the raw coupling (no electron
   translation factors). The analytic result above is also raw. State it.
   Translational sum rule for the raw coupling (the gate):
   `sum_A d^A_IJ = sum_pq gamma^IJ_a,pq <phi_p| d/dr |phi_q>`
   (rigid translation of all AOs = minus electronic momentum contraction);
   both sides computable — this catches any dropped Pulay piece.
4. **Multiplicity guard:** everything above is per-`mrst`; the Python path
   must branch (or abort) on `mult != 1` until the triplet fold is wired.

---

## 6. Implementation map (what exists / what to write)

| Term | Status | Where |
|---|---|---|
| `h_skel` 2e bilinear | EXISTS, selftested at I=J | `mrsf_nac_amp` + `grd2_mrsf_nac_compute_data_t` |
| `h_skel` 1e/Fock/XC (`esum`) | EXISTS, FD-validated | `mrsf_nac_esum` |
| `gamma_a` (derived, all blocks) | REWRITE (derive; drop sign-scan + cross-only mask) | replaces `_build_nac_gamma_tlf` |
| `Sk` contraction | EXISTS | `mrsf_nac_overlap` / `der_overlap_matrix_ket` |
| `L` analytic | PARTIAL (2e verified at I=J; E-part = nac_wsx commutators; assemble bilinear I≠J) | Phase-11 exports + new assembly |
| `W^IJ . S^x` | COMPUTED, UNUSED — wire in | `nac_wsx` block, `tdhf_mrsf_gradient.F90:1372` |
| z-vector interstate solve | EXISTS (RHS hook + GMRES) | `OQP::nac_orbgrad_L`, `set_mrsf_nac_cphf` |
| `z . B^x` contraction | EXISTS | gradient seam (gZ - gS) |
| FD oracle sign fix | TODO first | `single_point.py:2273, 2196` |
| Gates: L rotation-FD (unfolded), interchange identity, HST projected C1 reference, raw-translation identity, frozen numeric refs with assert | PARTIAL/UPDATED; see 7.61 onward | `tests/`, `tools/nac_lagrangian/` |

### Sequencing (non-negotiable, per the audit's lessons)

1. Fix + freeze the FD oracle (signs included). Commit reference values.
2. Derive & unit-test `gamma_a` against determinant-overlap linearization
   (no sign scan).
3. Analytic `L` vs unfolded rotation-FD, per block.
4. Interchange identity on one geometry.
5. Full assembly vs oracle: signed cos (no alignment), both (I,J) and
   (J,I), H2O pairs + one C1 system (ethylene MECI), sum rule.
6. Only then: benchmark table, GAMESS cross-check.

---

## 7. Gate results (2026-07-31, chc3, H2O/BHHLYP/6-31G*)

### 7.1 Oracle convention fix CONFIRMED (commit 28eaaa7)
Post-fix probe: all three pairs `cos(nacv, (E_j-E_i)*d) = +1.000000, ratio
1.000000`; `|d+d^T| = 0` and `|h-h^T| = 0` exactly. Canonical reference
magnitudes (to freeze): |d(1,2)| = 0.126363, |d(1,3)| = 0.036786,
|d(2,3)| = 0.528614.

CAUTION (infrastructure): `_run_oqp_external` prefers `shutil.which('openqp')`
over the venv that launched the master — on the Mac a stale global
`/opt/homebrew/bin/openqp` silently served all displaced workers (2x-off
|d| from its missing HST factor). Cluster runs were clean (PATH fallback to
`sys.executable -m oqp.pyoqp`). Pin PATH (or patch the resolver) in every
harness.

### 7.2 gamma^IJ: DERIVED closed form, machine-precision confirmed
With the operator phase convention |i->a> = a+_{a beta} a_{i alpha}|ref>
(amplitude-to-sorted-determinant parity (-1)^(noca-1-i); beta parity
det-independent, cancels in bilinears):

```
gamma_alpha[j,i] = ov*delta_ji - sum_a Xt^I[i,a] Xt^J[j,a]   (occ x occ)
gamma_beta [a,b] = ov*delta_ab(core) + sum_i Xt^I[i,a] Xt^J[i,b] (virt x virt)
ov = sum(Xt^I * Xt^J)  (=0 for I != J)
```

Gate A: closed form == exact Slater-Condon to 5e-20/1e-16 on ALL pairs.
Gate D (code-independent): FD of the exact biorthogonal determinant overlap
== sum gamma_pq K_pq to 1e-10..1e-23 on EVERY rotation block. gamma is
exactly the first-order orbital response of the true SF wavefunction
overlap. Structurally: gamma_dv = 0 exactly (no doc-virt block).

### 7.3 The branch TLF kernel is wrong beyond a calibration
vs the exact gamma_a: cross blocks (ds, sv) uniformly ratio ~sqrt(2)
(1.4068..1.4139) and globally sign-flipped (cos = -1.0000 every cross
block); ss/vv/dd blocks MISSING from TLF (zeroed as "gauge") while exact
has them — for the S0-touching pair (1,2) the ss block DOMINATES
(|ss| = 1.003 vs cross ~0.05; exact-overlap ss response = +1.699). The
branch's global 1/2 calibration and sign alignment absorbed sqrt(2)/sign
partially, which is why it could "validate" while wrong.

### 7.4 TLF overlap response != exact overlap response (gate C vs D)
The Fortran TLF (compute_states_overlap) rotation response disagrees with
the exact-overlap response at FIRST order: spurious dv response (0.08-0.54
where exact is 0 exactly), ss response 0 where exact is 1.699, sv wrong
sign/magnitude for some pairs. Consequences:
1. The TLF-based numerical NAC oracle inherits O(1)-per-theta response
   errors UNLESS they cancel in geometric displacements — must quantify
   (check ndtlf default/convergence).
2. We now have an EXACT-overlap numerical oracle for free: the biorthogonal
   det(M) machinery evaluated with the cross-geometry MO overlap. Use IT to
   freeze reference NACs; TLF remains the production fast path.
3. The nonzero exact ss response with frozen amplitudes is the transport/
   gauge coupling (amplitudes must co-rotate with socc); it does not drop —
   the Lagrangian's U^x_ss handles it. Do NOT zero same-space blocks.

### 7.5 Gate C follow-up: tlf=2 is the default -- the mismatch is not truncation
`infos%tddft%tlf` defaults to 2 ("most accurate"), so gate C already ran at
the highest TLF order. Working hypothesis for the C-vs-D discrepancy: the
get_states_overlap PIPELINE (AO-overlap -> s_mo -> possible Jacobi/
alignment preprocessing) computes the overlap response in an ALIGNED gauge
-- consistent with the observed ss response = 0 (alignment absorbs
same-space rotations) and apparent block transfer (dv response where the
raw wavefunction has none). If so this is a GAUGE choice, not a bug; the
analytic assembly must either (a) match that gauge (transport term), or
(b) validate against the raw exact-overlap oracle (gate D machinery) with
NO alignment, and treat the aligned pipeline as production-only. Decide by
reading compute_states_overlap's s_mo preprocessing; do not guess.

### 7.6 KNU-GAMESS cross-check + the broken tlf=0 exact path (2026-07-31)
- GAMESS (namd.src) aligns MO phase/order via `sffase` BEFORE the TLF
  formula; OpenQP's get_structures_ao_overlap does no alignment (raw
  C_old^T S C_new). Gauge difference to keep in mind when comparing
  numerical NACs across the two codes.
- GAMESS `ndtlf=0` = OVEXACT = EXACT hole-pair minor determinants ("no
  tlf"); OpenQP ported it (`ov_exact`, case(0)) but the header comment is
  inverted ("less accurate") and the path was NEVER exercised (default
  tlf=2).
- OpenQP `tlf=0` CRASHES: heap corruption (munmap/SIGSEGV in free()),
  gdb -> mrsf_tlf, valgrind -> TWO invalid 8-byte writes inside ov_exact
  (one landing 8 bytes BEFORE ddet => iipp <= 0, one past a
  compute_states_overlap 720-byte block = s_ia). Static port-diff vs
  GAMESS ityp=1 shows faithful structure; bounds-checked debug build
  (-fcheck=bounds) queued to name the exact line. IMPORTANT: ityp=3
  (s_ia) runs for ALL tlf values -- must verify the production tlf=2 path
  is clean under bounds checking too.
- Davidson state phases are RANDOM per run (3-run experiment: magnitudes
  reproduce to 1e-10, pair signs form a consistent state gauge, product
  of the three pair signs = +1). Regression test resolves the gauge and
  enforces the product rule. Consider pinning phases canonically
  (largest-|X| component > 0) in the production Davidson for NAMD
  stability.
- Infrastructure: `_run_oqp_external` prefers PATH `openqp` over the
  launching venv -- a stale global install silently serves all displaced
  workers (observed on Mac: 2x-off |d| from the old pre-HST-fix global).
  Pin PATH or resolve the console script next to sys.executable.

### 7.7 tlf=0 fixed + THE REPRESENTATION RESOLUTION (2026-07-31)
- After the diagonal-minor fix (2242bff) tlf=0 runs clean (bounds-check
  silent, S diag = 1 exactly).
- Gate C at tlf=0 gives the SAME first-order responses as tlf=2 -- TLF
  truncation differences are higher order in the off-diagonals, so the
  response mismatch vs the raw-determinant gamma is NOT truncation error.
- CONCLUSION: the GAMESS/OpenQP state-overlap FORMULA (7-term contraction
  with sqrt2 SOMO weights) is a SPIN-ADAPTED (CSF-aware) representation of
  the MRSF states, not the raw 90-determinant expansion. Its first-order
  orbital response therefore differs STRUCTURALLY from the raw-det gamma:
  nonzero dv response (partner-determinant contributions outside the
  single-SF space), zero ss response, and the sqrt(2) scale that the
  branch's sign-scanned kernel exhibited (ratio 1.407-1.414). BOTH gammas
  are internally consistent objects:
    raw gamma_a  <-> raw SF-CI representation   (machine-verified, gate D)
    gamma^formula <-> the published MRSF NACME definition (TLF/NAMD)
- The production numerical NAC measures d in the FORMULA representation,
  so the analytic NAC deliverable must use gamma^formula -- DERIVED from
  the 7-term structure with exact minors (tlf=0), not sign-scanned. The
  raw-gamma machinery stays as the independent self-consistency gate.
- NEXT: differentiate compute_states_overlap's contraction analytically:
  S_IJ = sum over (alpham,betam,gammam,deltam) x {s_ij, s_ab, s_ia}
  with s_* = exact hole-pair minors; d(s_*)/dtheta are cofactor
  contractions (Jacobi's formula on the minors). This yields the exact
  formula-consistent gamma^IJ_pq in closed form, including the dv block
  and the sqrt2 weights, with no fitted signs.

### 7.8 gamma^formula derivation campaign (in progress)
- gammaTLF (branch kernel) scored against the formula-FD: ss EXACT (their
  gauge claim holds in the formula representation), but dv response MISSING
  entirely and ds/sv wrong per pair (20%..100% errors, spurious (1,3) sv).
  The 1e-12 branch validation must have probed a restricted set. So the
  formula kernel must be DERIVED; neither raw gamma nor gammaTLF is it.
- Formula structure (compute_states_overlap): 6 explicit terms + column
  NORMALIZATION (s_st(:,i)/norm2) -- part of the formula's definition.
  s_ia enters BILINEARLY in the sqrt2 terms => first-order dv response
  comes from ds_ia(dv) x s_ia0(socc baseline). Confirmed mechanism.
- numpy replica of the full contraction built on exact minors: close but
  NOT exact -- first-order mismatch in ALL rotation blocks (linear in
  theta), meaning a minor-convention detail (transpose/index order/sign)
  differs. Unit-amplitude probe (contraction collapses to
  s_ij*s_ab + s_ia*s_ia) launched to extract the exact conventions from
  the Fortran directly, testing all transpose/role-swap variants.

### 7.9 CONVERGENCE: replica exact + THE ORIENTATION VERDICT (2026-07-31)

**The numpy<->Fortran transpose trap (root cause of a whole findings
family).** Every 2-D tagarray reaches numpy in C order = the TRANSPOSE of
the Fortran matrix. Consequences found and fixed:
1. My gate-C staging rotated MOs in the wrong space (fixed:
   W -> expm(-th K) @ W). The earlier "formula has dv response /
   TLF-vs-exact structural mismatch" conclusions (7.4) were artifacts of
   this. With correct staging the formula's dv response is 0 EXACTLY.
2. My replica-vs-Fortran S comparison read S transposed (SP == SF^T
   element-for-element).
3. PRODUCTION BUG: NACME.nacme() reads td_states_overlap without .T, so
   dc_python held d_ji = -d_ij. Verified ABSOLUTELY (not by inspection)
   with a code-independent exact biorthogonal-overlap oracle at a real
   geometry displacement: dc_py sign-opposed to <I|dJ> on every pair.
   Fixed at the storage boundary (commit 60c5412), keeping 28eaaa7's
   canonical gap axes -- exactly the split the audit gate demanded.

**Formula replica: EXACT.** Literal ov_exact minors (incl. case(3)
overwrite semantics) + 8-term contraction + column normalization
reproduce Fortran S to 1.9e-16 at finite rotations. The formula is now a
transparent Python object that can be differentiated term by term.

**gamma^formula: extracted and gated.** Generator sweep of the exact
replica (Richardson, /2 slot normalization) gives gamma[I,J,p,q] =
dS_IJ/dtheta_pq; the gate sum(gamma o K) == Fortran FD passed on every
block up to the slot-normalization fix (all residuals were exactly 2x
before it, machine-small after). Structure at the reference:
  dv = 0 exactly; ds ~= the branch TLF kernel (1e-4); ss real and equal
  in magnitude to the raw-determinant gamma with opposite sign; sv a
  distinct structure matching neither prior kernel.
Saved as gamma_formula_h2o.npz; wiring target: OQP::nac_gamma_tlf.

**Remaining assembly (unchanged plan):** d_orb = gamma^formula . (U^x + Sk)
via the Lagrangian (z-vector on Ltot = L + gap*gamma_a; -Tr[W S^x] partner);
d_amp = [h_skel + L:U^x]/gap with the validated esum/bilinear pieces.
Validation against the exact-overlap numerical oracle (tlf=0, now fixed)
with signed cos, emergent antisymmetry, C1 system, sum rule.

### 7.10 Campaign close-out (2026-07-31 evening)
- Post boundary-fix package: full suite 9 passed (flipped REF_D + the
  gauge-product rule -- a global d-sign error cannot pass it).
- gamma^formula final gate PASSED (slot normalization /2); replica exact at
  2.1e-16. gamma_formula_h2o.npz saved on chc3.
- EXACT-overlap (tlf=0) numerical references frozen; KEY RESULT: tlf=0 and
  tlf=2 numerical NACs agree to 1e-10..1e-12 (gauge-resolved) at dx=1e-3 --
  the TLF truncation error is NEGLIGIBLE in geometric FD, so the production
  oracle is effectively exact and both flavors are now available and fixed.
- orientation_gate harness: phase-alignment hygiene still crude (S_py axis
  mix-up in the aligned rerun); the definitive orientation evidence is the
  run-1 clean verdict + the suite's gauge-product proof. Tidy later.

NEXT SESSION ENTRY POINT (the remaining Lagrangian assembly):
1. Wire gamma^formula into OQP::nac_gamma_tlf -> d_orb frozen+skeleton via
   mrsf_nac_overlap (replaces the sign-scanned kernel).
2. Ltot = L + gap*gamma_a^formula -> interstate z-vector (hook exists);
   W^IJ.S^x partner (nac_wsx, computed-but-unused).
3. d_amp: h_skel (esum + bilinear 2e, both validated) + L:U^x via the same
   z-vector solve.
4. Assemble d_IJ; gate vs the frozen exact-overlap references with SIGNED
   cos, emergent antisymmetry (compute (I,J),(J,I) independently), C1
   molecule (ethylene MECI), translational sum rule.

### 7.11 Assembly gate v1/v2 (2026-07-31 night): direction PROVEN, bookkeeping open
d_pred = Xt_I.dXt_J + gamma^formula:T vs production d_num, in-process:
- SIGNED cos = +0.984 / |1.0000| / |0.999| on all pairs -- the decomposition
  DIRECTION and the storage-boundary sign fix are confirmed end to end.
- Magnitudes off per pair (0.9x / 1.6-2.7x / 3.2x). Diagnostics:
  * The formula's amplitude metric is trivial along the actual dX directions
    (exact directional derivative == plain dot; metric correction changed
    nothing) -- the excess sits in the gamma:T term.
  * For (2,3), d_num ~= the amplitude term ALONE (0.51 vs 0.59): the
    gamma^formula:T contraction largely double-counts response already
    carried by the ALIGNED displaced amplitudes (the amplitude derivative
    includes transport/co-rotation content). The clean split needs the
    amplitude derivative taken at FROZEN orbital gauge (transported dX)
    with gamma:T carrying all orbital response -- i.e., the transport
    convention must be fixed ONE way in both terms (branch Phase-10 lesson,
    now precisely quantifiable with exact tools).
  * Harness hygiene: run d_num FIRST on the pristine state (assembly runs
    showed (1,3) |d| instability 0.037<->0.065 caused by mol-state mutation
    before numerical_nac; the three pristine freeze runs agreed to 1e-10).
NEXT-SESSION plan (precise):
1. Reorder assembly_gate: d_num first, then sweeps; keep ONE gauge.
2. Transported amplitude derivative: dXt_J^transported = Q^T Xt_J(x) Q-style
   per-block Loewdin transport (branch _compute_amp_damp machinery,
   validated) -> amp term = formula-metric derivative of TRANSPORTED
   amplitudes; gamma:T keeps the FULL T. Gate again -- expect magnitudes
   to close; any residual = the same-space canonical U^x sector.
3. Then term-by-term analytic replacement (z-vector, esum, W.S^x).

### 7.12 *** THE MASTER DECOMPOSITION IS PROVEN (2026-07-31, late) ***
Assembly gate v4 (three-way, antisymmetrized), H2O all pairs:
  LEG1 replica==chain: cos +1.00000000 (1e-6..1e-7)  -- chain rule exact
  LEG2 replica==num:   cos +0.999994..+0.9999998, max|diff| 3.6e-5..5.8e-4
  |d_chain| vs |d_num|: 0.126614/0.126363, 0.036848/0.036786,
                        0.528997/0.528614  (~0.1%, h^2-limited)

  d_num[I,J] = antisym_IJ[ dS/dXt . dXt_J/dx + gamma^formula_IJ : T ],
  T = dM/dx (M = column-normalized C0^T S_AO C(x))

Keys that closed it:
- the production pipeline's (S - S^T) numerator means all comparisons are
  on the ANTISYMMETRIZED combination (the formula's raw dS carries a real
  symmetric part; antisymmetry of d is EMERGENT in this assembly);
- run d_num on the pristine state (harness order);
- the formula amplitude metric is trivial along actual dX directions.

ANALYTIC REPLACEMENT LADDER (each step gated against this scaffold):
 A. dXt_J/dx -> first-order PT on the matvec: (At-Om_J)dXt_J = -(dAt-dOm)Xt_J
    projected: Xt_I.dXt_J = [Xt_I (dAt) Xt_J]/(Om_J-Om_I) -- h_skel already
    validated (esum + bilinear 2e engines); the L:U^x part shares the same
    T machinery below.
 B. T -> Sk (der_overlap_matrix_ket, exists) + U^x (interstate z-vector via
    OQP::nac_orbgrad_L + same-space canonical terms; W^IJ.S^x partner).
 C. gamma^formula -> closed form by Jacobi cofactors of the exact minors
    (mechanical; the replica defines every term).
 D. Wire gamma^formula + Sk into mrsf_nac_overlap (tag interface ready);
    ship d_amp via one z-vector per pair. Benchmarks: H2O + C1 ethylene
    MECI + sum rule; freeze references with the exact-overlap oracle.

### 7.13 Ladder A resolved in principle: why TRANSPORT is mandatory (final)
Diagnostics (chc3): apply_A (mrsf_matvec_apply) is exactly linear (1e-15)
and eigen-consistent (Rayleigh quotients match Omega to 8 digits); the
reference and displaced amplitude bases overlap to |dot| ~ 0.9999. Yet the
raw-frame PT numerator X0_I^T dA X0_J is numerically DESTROYED: the O(h)
eigenvector mismatch delta (||delta||^2 ~ 1.7e-4) couples through the
matvec's HUGE spectral radius (~300 Ha, core excitations) giving
second-order Rayleigh contamination ~0.05 on diagonals and swamping the
~1e-5 off-diagonal signal entirely. THIS is the first-principles,
quantitative justification of the branch's Loewdin-transported matvec:
transport removes the frame mismatch before the spectral radius can
amplify it. Ladder A's analytic path is therefore the TRANSPORTED-frame
PT -- whose skeleton pieces (esum 1e/2e/XC + bilinear 2e) the branch
already FD-validated -- plus the transport-gauge term and L:U^x via the
interstate z-vector. The v4 scaffold remains the exact referee for each
replacement step.

### 7.14 Ladder A2: the response target quantified (2026-07-31, closing)
Scaffold amplitude term vs the ANALYTIC Fortran skeletons (mrsf_nac_amp
bilinear 2e + mrsf_nac_esum 1e/2e/XC), H2O:
- (2,3): amp2e/gap anti-parallel at cos -0.9986, 4.4x -- the giant
  cancellation structure; esum small here.
- S0 pairs (1,2)/(1,3): skeleton/gap is 100-2000x the scaffold amp -- the
  response term must cancel it almost exactly (matches branch Sec.6
  quantitatively, now in the proven convention frame).
- Residual-vs-wsx cos = +-1.000000 is SYMMETRY-FORCED (C2v irrep blocks
  are 1-D; branch landmine #4) -- only ratios are diagnostic: 2.84 / 72.9
  / 1.79. wsx alone does not close the residual; the antisymmetric part
  requires the interstate z-vector, as the Lagrangian derivation says.
Structural conclusion: the amplitude residual = L:U^x/gap with the SAME
U^x entering gamma:T. One z-vector infrastructure closes both. THE single
remaining derivation: the interstate orbital gradient L of the matvec
(assemble from the Phase-11-verified ingredients G_MO/fa/fb in the
unfolded representation), then per-pair z-vector via the existing hook,
then C1 validation (ethylene MECI) where cosines are meaningful.

Frozen per-pair z-vector targets (this run's gauge):
  (1,2) |resid| = 1.157*gap;  (1,3) 0.529*gap;  (2,3) 3.097*gap.

### 7.15 Ladder A3: response direction PROVEN, same-space L is the last piece
Polarized production-RHS L (exact for the bilinear RHS; NAC_DUMP_RHS hook,
no code changes) contracted with U^x = T - Sk over the THREE rotation
blocks, vs the frozen residuals:
  -L.U^x/gap: cos +1.000000 / +1.000000 / +0.998317 (sign = MINUS, matching
  sfrorhs's global -(...)); magnitude coverage 13% / 28% / 71%.
The shortfall is structural, not a factor: the production RHS lives ONLY
in the 86-dim inter-space rotation space, so polarization cannot deliver
L's SAME-SPACE blocks -- and the S0 pairs are ss-dominated (the same
structure gamma^formula showed). THE remaining derivation, now uniquely
identified: assemble L's same-space blocks directly from the verified
matvec exports (G_MO / gchan / fa/fb, unfolded representation), contract
with U^x_same-space (canonical-orbital response), and the amplitude term
closes. Then wire the z-vector for the inter-space part (interchange) and
ship. All targets frozen; every step has an exact referee.

### 7.16 Ladder A4 (session close): the bookkeeping theorem, empirically
Full-block L (matvec rotation-FD, frozen-Fock, antisym 171 + sym 190
generators) contracted with the measured full U^x does NOT close the
amplitude residual (coverage 15-45%; per-pair best signs disagree).
Anatomy frozen: (2,3) sv-dominated (1.404/gap), (1,2) ss (0.256) + sv
(0.151), dd/vv negligible. CONCLUSION (= derivation Sec. 4.2's warning,
now measured): the FOCK-REBUILD response dominates the remainder, and
frozen-L x U^x is an INCONSISTENT pairing -- the response enters
consistently only through the interchange triple (frozen-L RHS,
sfrolhs LHS with the Fock response, B^x contraction), i.e. the
interstate z-vector, plus the same-space canonical F^x term.

NEXT-SESSION SPEC (single, concrete):
1. Inter-space: solve sfrolhs z with RHS = polarized frozen-L
   (A3 machinery, exact) restricted to the rotation space; response
   contribution = z . B^x per coordinate (B^x via the existing gradient
   seam gZ-gS). Gate vs the frozen residuals.
2. Same-space: canonical response U^ss_pq = F^x_pq/(eps_q - eps_p)
   (F^x from the CPHF-consistent Fock derivative; start with the
   skeleton F^x from esum's machinery) x L^ss (measured here; analytic
   assembly from G_MO/gchan/fa/fb after gating).
3. Assemble, close the amplitude term, then full d vs the v4 scaffold,
   then C1 (ethylene MECI) + sum rule.
All referees frozen: v4 scaffold, per-pair residuals, A4 block anatomy.

### 7.17 Ladder A5 (session end): the two decompositions reconciled
Branch transported damp vs our raw scaffold amp (correct same-scale
comparison, no extra /gap): (2,3) 0.520 vs 0.588 at cos +0.9947; S0 pairs:
raw amp ~1e-4 (states follow orbitals rigidly) vs transported 0.028-0.119.
CONCLUSION: two EQUIVALENT decompositions exist --
  ours  : d = raw-amp + gamma^formula : T(full, incl. same-space)   [v4-proven]
  branch: d = transported-amp + cross-only d_ov
differing by where the same-space co-rotation content lives. The branch's
fitted 0.5 factors patched exactly this migration. Our v4 pairing is the
internally-consistent one (proven to 0.1% with zero fitted constants).
Production path stands as spec'd in 7.16; the transported-damp machinery
remains a valid semi-numeric reference for the amplitude term.

### 7.18 Ladder A6: THE INTERCHANGE PATH IS OPERATIONAL (session true end)
The production-grade route -- polarized L^{IJ} pushed through the
OQP::nac_orbgrad_L hook, ONE sfrolhs z-vector solve per pair, z.B^x
extracted by the gZ-gS gradient seam -- reproduces the direct L.U^x
contraction to ~1% (0.1503/0.1491, 0.1508/0.1507, 2.1876/2.1852) and
aligns with the frozen residuals at cos +0.9998..+1.0000 (sign +zB).
This verifies BOTH the Handy-Schaefer interchange through the actual
Fortran machinery AND that the measured U^x satisfies the CPHF equation
consistently. Coverage 13/28/71% as established; the sole remaining term
is the same-space canonical response (magnitudes and blocks frozen in
7.14/7.16). The analytic MRSF NAC is now: proven scaffold + verified
skeleton engines + OPERATIONAL z-vector response + one spec'd remaining
term -- assembly, not research.

### 7.19 Ladder A8: THE AMPLITUDE TERM CLOSES -- physical, scan-free
With the symmetric-generator staging corrected ((C e^{tS})^T = e^{+tS} W),
the SINGLE physical combination

    resp = z.B^x + (La.Ua + Ls.Us)_same-space + (Ls.Us)_inter-space-sym

closes the frozen amplitude residuals at 95.0% / 64.3% / 93.3% with ALL
signs = +1 exactly as derived (diagnostic sign scan finds no better
combination for (1,2)/(2,3)). Absolute misses 0.058/0.189/0.207 -- the
stacked-FD noise floor plus the small inter-space-symmetric Fock
response. Every term is production-implementable:
  t1: OPERATIONAL (A6) -- one sfrolhs solve + gZ-gS seam;
  t2: same-space L from the matvec exports (NO Fock rebuild there --
      same-space rotations leave D invariant, proven);
  t3: the orthonormality/W.S^x family (skeleton machinery exists).
The analytic MRSF NAC is COMPLETE as a specification: proven scaffold
(0.1%), verified skeleton engines, operational z-vector response, and a
closed, parameter-free response formula. Remaining: Fortran residency of
t2/t3, the gamma^formula cofactor closed form, C1 validation, ship.

### 7.20 *** C1 VALIDATION PASSED (ethylene, twisted, no symmetry) ***
The master decomposition on a C1 molecule (twisted ethylene, 18 coords,
S0/S1/S2), where cosines are fully diagnostic:
  LEG1 chain rule: cos +1.00000000 (1e-6) -- exact.
  LEG2 vs production numerical NAC (signed, per component):
    (1,2) cos +0.99998, |d| 0.2175/0.2149 (1.2%)
    (1,3) cos +1.00000000, |d| 0.9550/0.9552 (0.02%), max comp 7.8e-5
    (2,3) cos +0.99961, |d| 3.715/3.680 (0.9%, near-degenerate pair)
Component-level signed agreement (-0.30775/-0.30775, 0.42583/0.42576...).
The amp/orb split differs completely from H2O (amp-dominant here),
confirming the decomposition across regimes. Audit criterion 14 exceeded
(cos > 0.999, magnitudes 0.02-1.2%). The campaign's validation program is
COMPLETE: conventions proven, formula replicated, gamma derived,
decomposition proven on C2v AND C1, amplitude closed scan-free, response
operational through the production z-vector. Remaining work is Fortran
residency and packaging only.

### 7.21 Section C1 DONE + C2 partial (ethylene response closure)
C1 -- gamma^formula CLOSED FORM: cofactor-sensitivity assembly (adjugate
of the reference minors x multilinear contraction partials) matches the
generator sweep at MACHINE PRECISION (1.7e-13/5.2e-13/6.3e-14) after
fixing the s_ia sensitivity index map to the LITERAL overwrite layout
(clean-list rows were 6-long for socc-row minors -- the same trap as
7.9). The O(nbf^2)-sweep is replaced by one linear-algebra pass;
production-scalable.

C2 -- ethylene amplitude closure (restricted sweep): PARTIAL.
(1,2) 78.8% physical (93.6% with s1=-1), (1,3) 22-31%, (2,3) fails
(near-degenerate; zB overshoots 2x with wrong sign). Open items, precise:
 (a) unrestricted same-space sweep (vv/dd may matter in the amp-dominant
     regime, unlike H2O);
 (b) near-degenerate pair handling in the seam (CPHF conditioning +
     /gap amplification);
 (c) zB sign bookkeeping per pair/target ordering in the gradient seam.
Referees: v4-ethylene scaffold (7.20, 0.02-1.2%) stands regardless --
the DECOMPOSITION is C1-proven; only the response-term SPLIT needs the
remaining bookkeeping.

### 7.22 C2 item (a) RESOLVED-NEGATIVE; (b)/(c) sharpened (A10, session end)
The UNRESTRICTED ethylene L sweep changes nothing material vs the
restricted one (76.7/22.4/-301% vs 78.8/22.1/-298%): vv/dd blocks are
NOT the missing ethylene content. Sharpened diagnostics:
 (c) the diagnostic best signs want s1 = -1 for ALL ethylene pairs while
     H2O wanted +1 for all -- a MOLECULE-LEVEL sign in the zB seam chain
     (suspects: the polarized-RHS push/pack orientation vs the seam's
     gradient sign, state/target ordering in set_mrsf_nac_cphf). With
     s1=-1: (1,2) closes 93.4%.
 (b) (2,3) zB magnitude 2.2x the residual regardless of sign --
     near-degenerate conditioning (small gap, /gap amplification) needs
     regularized treatment.
 (1,3) at ~31% keeps a genuinely missing ~0.35 piece -- revisit after
     (b)/(c) are fixed, with the A10 npz (all L/U/zB data saved,
     no re-sweep needed).

### 7.23 C2 FINAL VERDICTS (A11) -- the empirical program closes
(c) Target ordering RULED OUT: zB(tgt=I) == zB(tgt=J) to all digits on
both molecules. The A9/A10 "uniform ethylene s1=-1" was FIT DEGENERACY
among comparable-size terms, not a real sign: standalone-zB best signs
are mixed with weak cosines (0.72/0.90/0.97). The derived all-plus
convention (proven on H2O: 95/93% scan-free) stands.
(b) CONFIRMED SHARP: ethylene gap(2,3) = 10.2 mHa; the /gap
amplification blows stacked-FD noise up to O(1) in the term SPLIT while
the TOTAL d stays at 0.9% (v4). The split is ill-conditioned near
degeneracy by construction.
PRODUCTION PRESCRIPTION (the deep lesson): do not assemble the response
term-by-term. Solve ONE z-vector with the COMBINED RHS
Ltot = L + gap*gamma_a (derivation 4.2) and contract once -- the
split-conditioning problem never arises, and every ingredient of Ltot is
now individually verified. The term-split gates (A2-A11) exist to certify
the ingredients, not to be the production path.

### 7.24 The production-candidate gates: skeleton CERTIFIED, six-term scatter localized
First production-candidate assembly (analytic_nac_gate.py: skel/gap + zL
+ zG + ov) gave right DIRECTIONS but wrong magnitudes (H2O 5-14x, ETH
0.5-1.5x). Term diagnostics (analytic_nac_diag/diag2.py) then localized
everything:

(1) gamma push/read conventions CERTIFIED PERFECT: the engine's own
    OQP::nac_trden_mo echo == pushed gamma^closed blockwise to 0.0e0 on
    both molecules. Python-created 3-D tagarrays keep their dims verbatim
    on the Fortran side (numpy (nbf^2,ns,ns) -> Fortran (nbf^2,ns,ns);
    the reversal happens only for FORTRAN-created records read into
    numpy). New landmine entry.
(2) OQP::nac_wsx = -Tr[W^IJ S^x] exists (exported by mrsf_nac_esum,
    Fock-weighted, the large-esum canceller) and was MISSING from the
    first assembly.
(3) Six-term unit sums [amp,esum,wsx,zL,zG,ov] still do NOT close
    (12-66%), but the LSQ fit residual on C1 ethylene is 0.2-1.9% =>
    the six engines SPAN d_num almost completely; the term DEFINITIONS
    mix content (coefficients scatter per pair).

FROZEN-MO SKELETON GATE (skel_gate.py) -- the decisive split. Freeze C,
X, D at the reference; displace geometry; 1-iter SCF from the reference
json guess (=> FOCK = F(x')[D_ref]); push frozen C; matvec the bilinear
E_IJ(x') = X_I^T A(x') X_J. Sanity: E_ref == diag(omega) to 8 digits.
VERDICT (H2O): FD_skel == amp2e + esum, cos +1.000000 ALL pairs,
maxdiff ~1e-5 = FD truncation. amp2e+esum is EXACTLY the skeleton
derivative, Gamma conventions included; wsx does NOT belong to the
skeleton (adding it breaks the match). ALL remaining error lives in the
response/CSF bookkeeping.

### 7.25 Exact-term assembly v3 (no fitted constants, no unknown objects)
The interstate Hellmann-Feynman identity fixes every term:

  d_IJ = [skel_IJ + M^IJ : U^x]/gap_IJ + gamma^IJ : (Sk_MO + U^x)

with M^IJ the FULL unrestricted orbital derivative of E_IJ, available
from the frozen A8/A10 generator sweeps: M = La + Ls (La antisym part,
Ls symmetric part incl diagonal; conventions verified against the sweep
staging line by line). PT identity: <X_I|dX_J> = X_I^T dA X_J/gap
EXACTLY (project the first-order PT sum onto X_I: only K=I survives).
Same-space antisym responses need no special handling here because the
exact FD U^x carries them; production replaces U^x contractions with
(i) ONE z-vector, combined antisym RHS [M + gap*gamma]_a via the
orbgrad hook (the hook antisymmetrizes internally), (ii) the symmetric
part -1/2 Tr[(M_sym)S^x_MO] via dSfull, (iii) same-space antisym via
the canonical chain U_ss = (eps S^x - Phi)/(deps) with Phi = F^x_skel +
G[dD] interchanged into one extra z-solve, (iv) gamma-side bookkeeping
via the certified ov engine + gamma_a z-RHS. New export NAC_DUMP_DS
(OQP::dbg_dsket/dbg_dsfull, bfnrm applied) supplies Sk/S^x matrices.
Gate = assembly_v3.py: C1 sym(U^x)==-1/2 S^x_MO; C2 amplitude part vs
frozen d_amp; C3 full d vs d_num, SIGNED.

### 7.26 v3d/v3e: the assembly WIRING is machine-certified; gauges identified
v3d/v3e (v3d_h2o.py, v3e_h2o.py) triangulated the v3b/v3c failures:

(1) T3: my ampdir == A8's amp_directional to 1e-11 (metric replication).
(2) CROSS-SESSION GAUGE LANDMINE: the frozen A8/A10 npz objects (La, Ls,
    Ux, d_amp) live in THAT session's orbital-sign and state-phase gauge.
    After a rebuild, pairing them singly with current-session objects
    (gamma:Ux, w_rot along Ux, C1 tests) is INVALID -- max|Ux_npz -
    Ux_now| = 2.99. Only gauge-invariant PAIRINGS (M:Ux both from one
    session) survive across sessions.
(3) T4a (v3e): the v4 identity REPRODUCED in-session at machine level:
    antisym[ampdir(dX_FD) + gamma:T_FD] vs d_num =
    (2,3): cos +1.00000000, maxdiff 1.1e-7 (EXACT);
    (1,2),(1,3): |pred| == |d_num| to all printed digits, cos -1 = the
    known inter-run Davidson phase (phi_1 = -1), pair-sign product +1.
    THE ASSEMBLY WIRING IS PROVEN. The missing step in v3c/v3d was the
    v4 ANTISYMMETRIZATION d = (dS_IJ - dS_JI)/2 (the formula dS carries
    a symmetric part).
(4) T4c == T4b to 1e-7: gamma:(Sk_analytic + Ux) == gamma:T_FD exactly.
    The dSket export's same-center omission (T1' forensics: on-site AO
    blocks exported as 0, FD nonzero) is BENIGN for the gamma
    contraction: C gamma C^T is antisymmetric while the on-site dSket
    content is symmetric -> contracts to zero. Convention documented.
(5) Open: dX_PT != dX_FD (T2'), traced to non-physical gauge content in
    the FD orbital response Ux entering through same-space blocks of
    the response bilinears. v3f probes it by an h-consistency map
    (physical U^x is h-independent; arbitrary degenerate-sector
    rotations are not) and re-judges the continuous-gauge assembly with
    a cleaned U.

### 7.27 v3f/v3g: THE ANALYTIC ASSEMBLY CLOSES ON H2O (machine level)
v3f killed the degenerate-gauge hypothesis: Ux is h-CONSISTENT
(|Ux(1e-3)-Ux(5e-4)| at FD-truncation level, cleaning mask empty), and
J3 proved gamma:(Sk_an + Ux) == gamma:T_FD at machine level inside the
total. The only wrong object was my w = dA X_J composition.

v3g found and fixed it:
(W3) SUM RULE: X_J^T w must equal d(omega_J)/dx. The exact referee
     w_ref (FULL displaced-SCF matvec on sg-transported reference
     amplitudes, transported back) satisfies it to 5e-7; my
     w_skel + w_rot violates it by up to 0.39. The missing channel is
     the Fock DENSITY-RESPONSE G[dD]: mrsf_matvec_apply is frozen-Fock
     (rotating C never rebuilds F_AO[D(C)]), so neither w_skel (frozen
     D) nor w_rot (rotated MO transform of the stored F_AO) contains
     G[dD]. Its signature: the w_ref - w_mine difference peaks on the
     core-Fock eigenvectors (w ~ 19-20 Ha).
(W4) PT FORMULA CERTIFIED: dX = sum_{k != J} V_k (V_k^T w_ref)/
     (om_J - w_k) reproduces the actual FD amplitude derivative to
     cos +0.999999..+0.99999, maxdiff 1.7e-4..1.5e-3 (FD truncation).
(W5) FULL ANALYTIC-STRUCTURE ASSEMBLY vs d_num (H2O):
       (1,3) cos +1.00000000 maxdiff 4.2e-9   (MACHINE EXACT)
       (1,2) |pred| == |d_num| all digits     (phase gauge phi1=-1)
       (2,3) 0.12%                            (phi2*phi3=-1; product +1)

THE CERTIFIED FORMULA (all ingredients defined, each with a referee):

  d_IJ = antisym_IJ[ ampdir_J(dX_J) + gamma^IJ : (Sk + U^x) ]
  dX_J = (om_J - A)^{-1}|_{perp J} (dA X_J)
  dA X_J = [AO-derivative skeleton, frozen C,D]
         + [MO-rotation response along U^x, frozen Fock]
         + [Fock density response G[dD(U^x)]]        <- the 7.27 term
  U^x = CPHF orbital response; Sk = C^T dSket_AO C (export, benign
        on-site omission); gamma = closed form; ampdir = the formula-
        metric directional derivative (linear in dX).

PRODUCTION ADJOINT FORM (per pair, no per-coordinate solves):
  ampdir_J(dX)[I] = G_met[I,J] . dX  (G_met analytic from the kernel)
  => amp term^c = ytil . (dA X_J)^c,  ytil = (om_J - A)^{-1} G_met|perp
     ONE amplitude-space solve per pair (MINRES on the matvec), then
     ytil.(dA X_J)^c decomposes into gradient-type contractions:
     skeleton -> esum/amp-engine pattern on the (ytil, X_J) pair;
     U^x-linear parts (rot + G[dD] + gamma:U^x) -> ONE z-vector with
     the combined RHS (the 7.23 prescription, now with every term
     derived and referee'd). This is the final production math.

### 7.28 ETHYLENE (C1) CLOSES -- the analytic assembly is PROVEN on both systems
v3g on ethylene (same script, ETH inputs), SIGNED, no phase flips this run:

  pair (1,2): cos +0.99999979  maxdiff 9.6e-5   (|d| 0.2148/0.2149)
  pair (1,3): cos +1.00000000  maxdiff 1.2e-5   (|d| 0.9552)
  pair (2,3): cos +0.99999997  maxdiff 4.6e-4   (|d| 3.6799/3.6805)

The (2,3) pair -- gap 10.2 mHa, the near-degeneracy that poisoned every
term-split approach (7.22-7.23) -- closes at 0.012% relative with the
exact resolvent denominators: the stacked-FD amplification problem is
GONE in the analytic structure. W4: PT == FD-dX to 3e-3 even across the
10 mHa gap. W3: w_ref sum rule 5e-6..1.5e-5 (exact); the G[dD] channel
on ethylene is a 1-2% correction (vs up to 60% on H2O) -- system-
dependent, never optional.

This exceeds the original audit criterion 14 (C1 cos > 0.999) by ~5
orders of magnitude. The formula of 7.27 is final. Remaining work is
IMPLEMENTATION ONLY (production adjoint form, end of 7.27): G_met
extraction, one MINRES amplitude solve per pair, the (ytil, X_J)
skeleton contraction engines, and ONE combined-RHS z-vector for all
U^x-linear content; then rewire analytical_nac() and refreeze the suite.

### 7.29 THE PER-PAIR ADJOINT FORM IS CERTIFIED ON BOTH MOLECULES
v4a (adjoint gate) + v3h (in-process merge; v4a's cross-process pairing
tripped the 7.26 phase landmine on H2O -- gate construction, not math):

  ADJOINT IDENTITY  ytil.w == G_met.PT(w):  1.7e-16 (H2O), 2.0e-15 (ETH)
  ADJOINT-FORM d vs d_num (signed, in-process phases):
      H2O: (1,2) 1.2e-6   (1,3) 4.2e-9   (2,3) 3.4e-4
      ETH: (1,2) 9.6e-5   (1,3) 1.2e-5   (2,3) 4.6e-4
  == the per-coordinate numbers exactly. ONE amplitude-space adjoint
  solve per pair (ytil = (om_J - A)^{-1}|perp G_met) replaces all
  per-coordinate amplitude work. G_met extracted by unit sweep passes a
  random-direction linearity check to FD precision.

STATUS: the derivation + validation program of this campaign is
COMPLETE. d = antisym[ ytil.(dA X_J)^c + gamma:(Sk + U^x)^c ] with every
object certified. Remaining work is ENGINEERING ONLY:
  (i)   MINRES the ytil solve on the matvec (no full A);
  (ii)  analytic (dA X_J): skeleton vector engine on the (ytil, X_J)
        pair (amp/esum generalization), U^x-linear content (rotation +
        G[dD] + gamma:U^x) via ONE combined-RHS z-vector per pair;
  (iii) closed-form G_met from the kernel (replace the unit sweep);
  (iv)  rewire analytical_nac(), refreeze references, extend the suite.

### 7.30 Production engineering: chain-polarization RULED OUT; the v7 plan
Toward replacing the FD ingredients (w_U, U^x channels) analytically:

(a) v5 (v5_prod_h2o.py): MINRES ytil == eigen to 2.6e-9 (production
    solve path OK). SLOT INJECTION works: pushing ytil into a bvec slot
    makes the existing bilinear engines (amp2e/esum/wsx) compute the
    (ytil, X_J) pair objects with NO new Fortran -- certified against
    ytil.w_skel at 1e-5 (gate G-A). The v5 piecewise assembly
    [engines + combined-RHS seam + gamma:Sk] reaches (1,2) 5% but
    misses the 2e-W-sym and same-space channels elsewhere.
(b) v5b/v5c/v6/v6b: polarizing the WHOLE production gradient chain
    (4-term bilinear extraction; +1/2 g(0) constant is REQUIRED -- the
    3-term form leaves -g0/2) REPRODUCES the exact bilinear for the
    (2,3) pair (2.4%) but FAILS for S0 pairs: the channel split shows
    a spurious net from imperfect cancellation between the z/P channel
    (+0.59) and the W channel (-0.44). The chain's internal elimination
    bookkeeping is STATE-SPECIFIC (normalization/eigen assumptions),
    not bilinear-safe. W build (mrsfrowcal) itself has no omega terms.
    VERDICT: do not polarize the chain; assemble the interstate terms
    explicitly (the v5 route) and complete the two missing channels.
(c) The missing analytic tool is G[P] for an arbitrary symmetric
    density. NO new Fortran needed for the gate: OQP_SCF_GUESS_D
    (the #154 native density guess) + 1-iteration SCF gives
    F[D + eps P] and hence G[P] by FD in eps. Production later adds a
    clean bind(C) wrapper around fock_jk.
(d) v7 plan: (i) sym channel -1/2 Mt_sym:S^x per coordinate via TWO
    frozen-Fock matvec directional derivatives along V_c = -1/2 S^x_MO
    PLUS the G-channel Tr[dD(V_c) G[P_F(ytil,X)]] (one G build/pair);
    (ii) same-space antisym canonical chain: Q_ss = [Mt+gamma]_a/deps;
    direct (eps S^x - F^x_skel) contraction (F^x_skel from the skel
    workers' FOCK export) + one extra seam solve with RHS = G[Q_ss];
    (iii) judge on H2O, then ethylene.

NOTE (session hygiene): a parallel session checked out another branch
in the primary clone; nac-lagrangian work continues in the dedicated
worktree repo_nac/. All refs intact; chc3 mirror unaffected.

### 7.31 v7a: full-sym channel helps (1,3)/(2,3), regresses (1,2) -- OPEN
v7a replaced wsx by T3 = [staged directional matvec along V_c = -1/2
S^x_MO] + [G-channel via the in-process G[P] build: DM records <- D +
eps*Pij (NAC_DUMP_PIJ export), oqp.hf_energy at control.maxit=1, G =
(F_eps - F0)/eps -- includes the XC kernel response automatically; all
records saved/restored]. Results (H2O, sign-resolved):
    (1,2) 8.6e-2 (v5-wsx: 6.2e-3  <- REGRESSION)
    (1,3) 3.7e-3 (v5-wsx: 6.0e-2  <- large improvement)
    (2,3) 2.1e-1 (v5-wsx: 3.1e-1  <- improvement)
T3 magnitudes: staged 0.40/0.20/0.85, G-ch 0.19/0.003/0.57 per pair.
OPEN QUESTION: the precise decomposition boundary between T3, the seam
(T2), and the same-space channel (T4, still unimplemented). The pair-
dependent sign of the v7a-vs-v5 shift says the current T3 either
double-counts part of the seam's internal W(z)S^x or mis-handles a
block weighting. NEXT: implement the U-channel terms EXACTLY as listed
in Sec 4 of this document (term-by-term, no shortcuts), reusing the
now-certified tools: MINRES ytil, slot-injection engines, staged
directional matvec, in-process G[P] build, seam solves, F^x_skel from
the FOCK-exporting skel workers, and judge each term against the
in-process FD channel referees (v3g machinery).

### 7.32 UNIT-L SEAM CALIBRATION + the calibrated bookkeeping (v7e/v7f)
Gauge invariance (v7c probe): [Mt_a + 2 gamma_a]_same-space == 0
numerically (socc pair: +1.418146 vs -1.418189, sum 4.3e-5) => T4 == 0
IDENTICALLY. The same-space channel cancels between the amplitude
response and gamma:U. (v7b's T4 explosion was a maxit-leak execution
bug: the G-build's control.maxit=1 leaked into the displaced SCFs.)

UNIT-L SEAM CALIBRATION (v7e_seamcal.py): inject L = e_pq into the
orbgrad hook; the seam returns
    seam(e_pq) = - U^x_FULL(p,q)      (r-fits -0.9994..-1.0000)
i.e. MINUS the FULL orbital response entry (sym content included), and
the hook antisymmetrizes L internally (e_qp -> -seam to 1e-9).

The algebra then fixes the assembly uniquely. With X := Mt_full + gamma
and U(lo,hi) = -U(hi,lo) - S^x eliminated:

  sum_pq X_pq U_pq = -seam(X)                                [T2 = -seam]
      - sum_cross-pairs X_(lo,hi) S^x                        \
      - 1/2 sum_same-space X_sym S^x                          } = X : V
  V[r,s] = -S^x[r,s] * (1/2 same-space; 1 if space(r)<space(s); else 0)

THE V-MASK IS EXACTLY THE OV-ENGINE'S trden_ss WEIGHTING (1/2 same-
space, 1 on isp<jsq) -- the Fortran comment said it all along; now it
is derived AND calibrated. gamma's own V-contraction (the cross
single-triangle -sum gamma_(lo,hi) S^x) is nonzero and included.
v7f: d = T1 - seam(Lt+gamma) + [staged(V) + Tr[dD(V) G[P~]] + gamma:V]
       + gamma:Sk,   with dD(V) = C(V D0 + D0 V^T)C^T for general V.

### 7.33 v7g U-channel forensics: gamma-side bookkeeping VERIFIED; the
remaining open item is the Mt-side discrete identity
v7g (v7g_h2o.py) compared, per pair per coordinate, the EXACT U-channel
referee  exact_U = ytil.w_ref - T1 + gamma:Ux_FD  against the
production candidate  prod_U = -seam(Lt+gamma) + T3''.

VERIFIED: the gamma-side bookkeeping closes where the same-space
channel is small: pair (1,3): [-seam(gamma) + gamma:V] = -0.19385 vs
exact gamma:U = -0.19027 (2%). The calibrated seam semantics
(seam(e_pq) = -U_full) and the ov-weight V-mask are CORRECT for the
gamma side.

OPEN: (a) pairs with a large same-space channel ((1,2): gamma_a(socc)
= -0.709, gamma:U_ss ~ -0.31) rely on the exact cancellation
[Mt + 2gamma]_a,ss = 0 (proven at generator level, 4.3e-5) carrying
over to the OMITTED-channel assembly -- but the assembled totals still
miss: (1,2) prod_U = +0.1767 vs exact +0.3555 (a clean 1/2 signature),
(1,3) prod_U = -0.0760 vs exact +0.0745 (clean -1), (2,3) wild. The
Mt-side pieces (staged directional matvec + G[P]-channel) carry a
remaining bookkeeping error with half/sign signatures.
(b) NEXT (fresh session): derive on paper the exact DISCRETE identity
sum_pq X_pq U_pq = -seam(X) + X:V + [ss-channel] with the CALIBRATED
seam semantics, then verify EACH equality against the v7g-style
referees (rerun v7g with np.savez of all pieces for offline algebra).
All tools certified and in place: MINRES ytil, slot-injection engines,
staged directional matvec, in-process G[P] build (hf_energy maxit=1),
double seams, FD referees, gauge-invariance probes.

### 7.34 v7h + offline algebra: the discrete identity is MAPPED
v7h harvested the FULL Mt matrices (single-element staging, 361 dirs x
6 pairs) + closed-form G-channel Mt_G[p,q] = 2 sum_s occ^s_q [C^T G_s
C]_pq + w_ref/Ux/double seams, all saved (H2O_energy_tlf0_v7h.npz).
Offline algebra on the npz established, at machine level:

(1) Ux_FD structure: OFF-DIAGONALS satisfy orthonormality U+U^T=-S^x
    exactly; DIAGONAL is zeroed by the sg-alignment gauge (the whole
    v3f "C1 failure" mystery resolved). Decomposition identity checks
    to 1e-8.
(2) gamma-side seam identity MACHINE-VERIFIED:
    -seam(gamma) == sum_pairs pack(gamma).U(hi,lo)  (2e-7..7e-4).
(3) ss-antisym channel of X = Mt+gamma auto-cancels (|ss_asym| ~ 1e-5;
    the gauge invariance Mt_a,ss = -gamma_a,ss inside the matrices).
(4) THE SEAM MISMATCH DIAGNOSED: the polarized-RHS L != pack(Mt)|rot --
    sfrorhs ADDS the "2*F*T orbital-relaxation terms that the Davidson
    matvec does NOT have" (its own comment/env NAC_ZERO_2FT). The
    production T2 must inject Lmat = Mt + gamma DIRECTLY into the hook
    (or use NAC_ZERO_2FT + reconcile), not the polarized RHS.
(5) Primary gate Mt-completeness: (Mt+Mt_G):Ux vs [ytil.w_ref - T1]:
    cos +1.000000 all pairs; magnitudes -2.9%/-3.3%/+12.7%. The full-d
    offline assembly T1 + X:Ux + gamma:Sk closes to 1.9e-2/8.7e-3/
    2.1e-1 (sign-resolved). Diag channel contributes zero.
OPEN (last layer): the 3%/13% Mt-channel deficit -- candidates: FD
referee precision (Ux ~0.3-1.5% rel; near-degenerate (8,9) worse) vs a
small missing Mt channel. Next: tighten the referee (smaller h,
Richardson) or cross-check Mt:U against ytil.(w_ref - w_skel) with
w_skel from the FOCK-exporting workers, then inject Lmat = Mt+gamma
into the hook (v7i, 6 solves) for the production seam form.

### 7.35 v7i: the DIRECT-INJECTION seam is certified; precision layer open
(1) PRODUCTION T2 CONFIRMED: injecting X = Mt + Mt_G + gamma directly
    into the orbgrad hook gives -seam(X) == sum_pairs pack(X).U(hi,lo)
    at 1.5e-5..1.2e-2 abs (vs |pack| up to 2.8) on ALL pairs. Together
    with 7.32/7.34 this certifies the complete production mechanism for
    the cross-space channel: ONE z-solve per pair on the X matrix.
(2) The h=5e-4 sweep of v7i is internally anomalous (Mt:U collapses to
    ~40-60% of the h=1e-3 value; Richardson made things worse, while
    the w_ref side stayed h-consistent). Sweep-2 state contamination or
    a degenerate-sector sg issue -- diagnose offline from
    H2O_energy_tlf0_v7i.npz (Ux1/Ux2/w1/w2 saved).
(3) Best current closures (h=1e-3 referees): offline assembly
    T1 + X:Ux + gamma:Sk = 1.9e-2 / 8.7e-3 / 2.1e-1 (sign-resolved).
    The remaining discrepancy is the Mt-channel deficit (-3%/-3%/+13%
    against [ytil.w_ref - T1]) whose origin (referee precision vs a
    small missing channel) is the LAST open question.

STATUS SUMMARY of the production program:
  CERTIFIED: T1 slot-injection engines; gamma:Sk; the seam transfer
  function (unit-L, -U_full) and the direct-injection T2; the V-mask
  weights; ss-antisym auto-cancellation; G[P] in-process build; MINRES
  ytil; the full-Mt harvest method; the decomposition identity (1e-8).
  OPEN: the few-percent Mt-channel deficit; the near-degenerate pair's
  FD-referee quality; then ethylene + rewiring.

### 7.36 THE DISPLACED-FRAME LANDMINE (session close-out)
Direct comparison of two same-build processes' displaced sweeps:
|Ux(v7h) - Ux(v7i)| up to 1.94 (entries in near-degenerate virtual /
space-boundary clusters); |w_ref| norms identical but entries differ
by 0.085 (frame rotations). THE DISPLACED-GEOMETRY ORBITAL GAUGE IS
RUN-NONDETERMINISTIC: within near-degenerate clusters the SCF mixing
at x +- h varies run to run; sign-alignment (sg) cannot fix continuous
mixing. Consequences, all consistent with every observation of this
campaign:
  - FD sweep data (Ux, w_ref, dX) is valid ONLY in-process. Cross-
    process pairings are invalid even at fixed phases (this, not the
    rebuild, also explains part of 7.26).
  - The certified machine-level closures (v3g/v3h, and the offline
    v7h total) work because the TOTAL is frame-gauge invariant: the
    amplitude channel and gamma:U carry the frame content coherently.
  - The seam identity survived cross-frame (1e-3) because rot-block U
    entries are frame-stable.
  - The v7h in-process Mt-completeness deficit (-3%/-3%/+13%) is the
    one standing quantitative open item; part of it may be frame/FD
    precision of the SPLIT (not missing physics), since the production
    seam-based form needs no FD-U at all.
NEXT SESSION RECIPE: (1) run ONE in-process gate combining: engines T1,
direct-injection seam T2, S^x-contraction terms (elim + ss-sym from the
X matrices), gamma:Sk, and judge vs d_num IN THE SAME PROCESS with the
Mt matrices harvested in that process (the v7h+v7i union script);
(2) if the few-% gap persists, compute the T3-channel exactly via
in-process FD (staged V-directions vs actual U-sym channel) to decide
referee-vs-physics; (3) ethylene; (4) rewire analytical_nac().

### 7.37 v7j UNIFIED GATE: the production form is CERTIFIED SELF-CONSISTENT;
the residual is a REAL small Mt channel
One process, one frame, no cross pairings (v7j_h2o.py):
  J4 seam identity: -seam(X) == pack(X).U at 2.7e-5..7.8e-4 (machine).
  J2 PRODUCTION FORM == J3 FD-U TOTAL to 3 digits:
     production: 1.896e-2 / 8.68e-3 / 2.06e-1  (sign-resolved maxdiff)
     FD-U total: 1.936e-2 / 8.71e-3 / 2.05e-1
  => the no-FD production assembly [T1 - seam(X) + elim(X) + ss_sym(Mt)
     + gamma:Sk] faithfully implements the current best decomposition.
  J1 Mt-completeness deficit REPRODUCES EXACTLY across processes:
     -2.9% / -3.3% / +12.8%  =>  a REAL missing channel in
     Mt = staged-frozen + G[P_Fock], not referee noise.
CANDIDATES for the missing Mt channel (next forensic layer):
  (a) the matvec may read MO energies (E_MO records) for its diagonal
      instead of C^T F C -- staged C-rotations would then miss the
      eps-diagonal response channel; check mrsf_matvec_apply's fa/fb/
      diagonal construction;
  (b) G-build fidelity for the response Fock (SCF hfscale vs response
      hfscale channels in the matvec's Fock usage);
  (c) w_ref referee frame content under near-degeneracy ((2,3) is an
      OVERSHOOT while normal pairs UNDERSHOOT -- two sources).
Everything else in the production chain is now certified. Final d gap
vs theory (machine): 1.9e-2 / 8.7e-3 / 2.1e-1 -- entirely attributable
to the J1 channel.

### 7.38 Candidate eliminations (v7k/v7l/v7m) -- the J1 channel search
narrows to the referee's displaced-Fock fidelity
- (a) eps/get_jacobi channel: RULED OUT -- eps-consistent staging
  (push E_MO = diag(C'^T F C') per staged matvec) leaves every J1..J4
  number identical to all printed digits.
- DM/SM direct dependence of the sigma: RULED OUT -- record
  perturbation probes give |dAx|/eps ~ 6e-12 (machine zero).
- (b') G-build contamination: RULED OUT at machine level -- the
  in-process hf_energy(maxit=1) G[P] passes eps-linearity 1.3e-10,
  self-adjointness 2.9e-11, G[0] = 7e-11. The G[P] tool is certified.
- Diagonal channel: M~ diagonals are EXACTLY zero in the harvest
  (frozen AND G) -- empirically no diagonal-scaling response of the
  bilinear; note MTG diag zero warrants one code-glance (GM diag).
REMAINING HYPOTHESES for the reproducible -3%/-3%/+13% J1 deficit:
  (i) the w_skel REFEREE side: the skel workers' 1-iteration SCF at
      DISPLACED geometries runs the full SinglePoint pipeline --
      possible vshift/DIIS handling of the stored FOCK at x' (the
      reference-geometry hf_energy G-build is clean, but the worker
      path was never given the same fidelity probe);
  (ii) finite-h second-order content in w_ref - w_skel;
  (iii) a genuinely quadratic-in-U cross term (h*U) not representable
      as M~:U at first order.
NEXT: probe (i) by rebuilding one displaced F(x')[D_ref] via the CLEAN
hf_energy route inside a worker and diffing the FOCK records; then (ii)
by an h-refined skel worker pair. All tools exist.

### 7.39 h-SCALING VERDICT: the J1 deficit is a REAL first-order channel
v7o (one process, dual-h sweeps): the deficit is h-INDEPENDENT
  (1,2): 0.03355 -> 0.03341 (h -> h/2);  (2,3): 0.3968 -> 0.3952.
NOT central-FD truncation. Combined with 7.38's eliminations and the
machine-level identity d_num == ytil.w_ref + gamma:(Sk+Ux) (v3g W5 +
the adjoint identity), the FULL COHERENCE CHAIN is:
  v3g total(machine) - J3 total == J1 deficit exactly
  => the deficit is the AMP-CHANNEL LINEARIZATION GAP:
     ytil.(w_ref - w_skel)  is NOT fully representable as  Mt : Ux
     with Mt = [C-rotation staging + G[dD]-channel], even though every
     record dependence of the matvec (VEC_MO, FOCK, E_MO, DM, SM) is
     accounted for. The missing content is h-independent, reproducible,
     and lives disproportionately on the near-degenerate pair.
Working hypothesis for its origin: frame/cluster mixing content in the
displaced-frame objects that cancels in the gauge-invariant TOTAL but
not in the (Mt, Ux) split; alternatively a genuine cross term outside
the (theta at x=0) linearization. Either way, PRODUCTION IMPLICATION:
do not build the amp channel from Mt:U. The two clean production
routes for the amp channel are
  (A) the derivative-sigma VECTOR engine: assemble w = (dA/dx) X
      directly in Fortran (derivative-ERI/Fock sigma build; the
      existing amp/esum engines are its TRACES) and contract with ytil
      -- exact by construction, no U-decomposition at all; or
  (B) keep the certified FD w_ref for validation and derive (A)
      against it.
The gamma channel is fully production-certified (seam + V-mask + Sk).
Data for the Delta-structure hunt: H2O_energy_tlf0_v7o.npz (Ux1/Ux2/
w1/w2 + all matrices, one frame).

### 7.40 Delta-structure hunt: H2O fits are C2v ARTIFACTS; ethylene run launched
Offline Delta characterization on H2O:
- Delta is TRANSLATIONALLY CLEAN (sums ~2e-5) => genuine electronic
  channel, well-defined.
- The seductive exact fit Delta(1,3) = c*gamma:Sk (residual 3e-10,
  c = -0.0624586 ~ -1/16) is a C2v ARTIFACT: that pair's irrep pattern
  is ONE-dimensional, so ANY channel in it is "exactly proportional" to
  any other. The old landmine holds: judge structure only on C1.
- The amplitude frame-transport hypothesis (Delta = G_met.(Sk-action on
  X)) FAILS: the transport vectors vanish identically for (1,3) while
  Delta does not; (2,3) fits at only ~21-24%.
=> the Delta-structure identification REQUIRES the ethylene dataset.
Launched: v7o on ETH (full Mt harvest 1444 dirs x 6 pairs + dual-h
sweeps + all matrices saved) -- the C1 fit basis for the channel hunt.

### 7.41 ETHYLENE Delta fits (C1): no clean single channel; PATCHING CLOSED
ETH v7o harvest (h-independent deficits confirmed again: (1,2) 0.017,
(1,3) 0.372, (2,3) 0.541) + phase-anchored ctx (w-probe alignment) +
offline C1 fits against {gamma:Sk, transposes, d_num, occ/virt Sk- and
S^x-transport channels}:
  (1,2): best single cos -0.65; LSQ 14% with degenerate coefficients.
  (1,3): several candidates all |cos| ~ 0.99 -- the pair's response is
         single-mode dominated, so even C1 leaves partial collinearity.
  (2,3): best single sxA (occ-side S^x transport) cos -0.957 (29%).
sxA recurs suggestively across pairs but with NO universal coefficient.
VERDICT: the Delta channel is not a simple one-term patch of Mt:U.
The candidate-patching program is CLOSED. The production completion is
route (A) of 7.39: the derivative-sigma amp channel built as a WHOLE
(interstate gradient-type contraction of the (ytil, X_J) pair with full
relaxation), gated per coordinate against the certified w_ref referee.
All datasets for that work are frozen: H2O_energy_tlf0_{v7h,v7i,v7o}.npz
and ETH_energy_{v7o,ctx,dnum}.npz (local + chc3).

### 7.42 v8 w-VECTOR FINGERPRINT: the Delta channel is the SOCC/FOLD
sector's orbital response (route-A target LOCALIZED)
v8_wprobe (H2O, c0=5, J=3): [w_ref - w_skel] vs [staged(Ux) + G-vector]
per raw amplitude slot: |lhs|=0.651, |rhs|=0.663, |diff|=0.089 (~14%).
The residual is STRUCTURED, not noise:
  - worst slot = LR1 (the +-RS spin-folded SOCC-pair slot): staged
    OVERSHOOTS 10x (lhs +0.006 vs rhs +0.067);
  - systematic socc2-row (i=6) misses across many virtuals
    (a=7,8,9,12,13,15,17);
  - doc rows only on the LUMO column (i=3,4 @ a=7).
=> the missing channel is the C-response of the MRSF SPIN-PAIRING /
FOLD machinery itself (the U-matrix pairing + spc scalings + the +-RS
fold conventions; JCP 158,194105 structures) -- a structured response
that naive single-element C-staging mis-weights on the socc/lr slots.
ROUTE-A first derivation item is therefore: the fold/pairing-sector
orbital derivative (kernel socc rows + lr slot), to be implemented in
mrsf_nac_wpair and gated per slot against this probe (v8_wprobe.py,
extendable to any coordinate/pair). Bonus: the G-channel VECTOR
(FOCK-record perturbation + matvec) measured here validates the
Mt_G closed form at the vector level.

### 7.43 v8/v8b/v8c/v9: the elimination program's terminus
- v8b == v8 to every digit (global phase only): the sg_apply transport
  and sign-fixed-C referee are EQUIVALENT (the matvec's exact lr2-axis
  annihilation makes the fold transport safe). Referee methodology
  exonerated.
- socc orbitals DO sign-flip asymmetrically at +-h displacements (sg
  probe) -- harmless given the above.
- v8c: the pure C-channel [A(C')-A(C)] at converged F matches the
  staged response to |diff| = 0.0197 (3%); the v8 "socc fingerprint"
  was G-vector misplacement in the probe wiring.
- v9: the MEASURED vector G-channel (per-coordinate G[dD] build +
  FOCK-perturbed matvec) reproduces the closed-form Mt_G contraction
  TO ALL DIGITS -- the closed form is CORRECT (another certification).
RESULT: the invariant residual (1.94e-2 / 8.7e-3 / 2.06e-1 on the
final d; ~3%/13% on the U-channel) survives every representation:
closed-G == measured-G; staged-C == true-C to 3%; h-independent;
process-independent; structured. The ONLY remaining unverified link in
the entire chain is the numerical STAGING as a proxy for the analytic
d/d(theta) of the sigma build. ROUTE-A PROPER therefore reduces to ONE
question: derive d(sigma)/d(theta) term-by-term FROM THE SIGMA SOURCE
(mrsfmntoia/mrsfcbc channels, spc scalings, fold application) and
compare against the staged directional derivative -- wherever they
differ is the final term. All probes, referees, and datasets for that
comparison are in place.

### 7.44 v10: THE PARADOX, precisely stated (session terminus)
v10 (single coordinate, in-process):
- The C-channel residual [A(C',INTS_d) - A(C0,INTS_d)]/FD minus the
  staged directional response is h-INDEPENDENT: 0.01965 (h) vs 0.01968
  (h/2). The invariant residual LIVES IN THE C-CHANNEL.
- (The in-process 1-iter "worker emulation" is broken -- trueG blew up
  to 32.8; discard. The F-channel accounting stands via 7.43's
  closed-G == measured-G.)
THE PARADOX: every link is individually certified --
  (i)   Ux == C0^{-1} dC_raw (orthonormality, exact);
  (ii)  the staging measures the code's own directional C-derivative
        (t = 1e-5, linearity certified);
  (iii) cross terms (d2A/dtheta dx, d2A/dtheta2, U(h)-curvature) cancel
        at EVEN order in the central FD;
yet the true-path C-response differs from the staged response by an
h-independent 3% (13% on the near-degenerate pair). One of the hidden
premises is false. Candidate premises to attack in the route-A paper
derivation: the smoothness/single-valuedness of the displaced SCF
solution C(x) in the near-degenerate virtual sector (is C(+h) vs
C(-h) on the SAME smooth branch?); the exactness of the sg-branch
choice; whether the matvec output depends on pieces of C beyond the
occupied+active blocks that the M-frame Ux does not constrain (virtual
-virtual rotations of C(x') are NOT fixed by any condition used here
-- THE VIRTUAL-GAUGE: Ux_vv is whatever the displaced SCF returns,
but the staged and true paths agree on it BY CONSTRUCTION since Ux is
measured from the same runs... unless sign-fixing branch-cuts).
NEXT SESSION: start exactly here, with v10_probe.py as the instrument
(fix its worker branch by using SEPARATE PROCESSES for the 1-iter
skel), and the derivation target of 7.43.

### 7.45 v11/v12: branch and curvature ELIMINATED; two named measurements left
- v11 B1: NO branch kink -- u+ vs u- agree to 6e-3 (vv worst), path
  smooth. B2: the ONE-SIDED FD differs from the staged response along
  its own direction by the same 0.0195 (== the central value).
- v12 t-scan: the staged response is LINEAR in t to 2e-4 over
  t = 1e-5..1e-3. No curvature.
Combined with 7.44: smooth path + linear response + (u+ defines C+
reconstruction) should force agreement -- so one of exactly TWO things
remains: (a) the reconstruction C0(1 + h u+) != C+_signfixed beyond
O(h^2) (measure |C0(1+h u+) - C+| directly), or (b) the MIXED second
derivative d2A/dx dtheta: the lhs C-channel is measured AT DISPLACED
INTEGRALS while the staging runs at reference integrals -- an O(h)
term that would explain 0.0195 IFF the pure-C-channel h-scaling is
NOT constant. The v10 h/2 "constancy" is untrustworthy (same-run
maxit-hack contamination, trueG blowup); the ONLY clean h-test so far
is v7o's J1 (which mixes channels). NEXT-SESSION MEASUREMENT #1:
fresh-process pure-C-channel at h and h/2 (no maxit hacks) + the
reconstruction check. If (b) confirms: the referee/staging mismatch is
a benign O(h) cross term and the TRUE production closure should be
judged only through h->0-extrapolated or integral-consistent referees
-- potentially dissolving the entire "invariant residual" as a
referee artifact at finite h, with d_num (production, its own FD) the
final arbiter.

### 7.46 v13: both named measurements done; the residual is an h^0
STATE-DEPENDENCE of the response operator
Clean process (no maxit hacks), one-sided, h and h/2:
  M1 reconstruction |C0(1+h u) - C+_sf| = 2.65e-4 -> 1.33e-4 (EXACT
     halving: pure O(h^2)); u is faithful. Premise (i) closed.
  M2 |lhs - staged| = 0.01952 -> 0.01957 (h-INDEPENDENT, clean).
     The mixed d2A/dxdtheta hypothesis (predicting halving) is DEAD.
Since every input difference between the lhs pair and the staged pair
(C: h*u certified; integrals: O(h); FOCK/DM/eps records: O(h)) scales
with h, an h-independent response difference forces the conclusion:
the matvec's response operator differs at O(h^0) between the
DISPLACED-state evaluation (lhs pair, run inside the displaced-SCF
process state) and the REFERENCE-state evaluation (staged pair) --
i.e., some piece of process state that the record set does not
capture (int2 driver internals / screening tables / module-level
caches?) shifts the C-response.
NEXT MEASUREMENT (#1 of next session): the STATE-SWAP bisection --
evaluate the staged pair [A(C0), A(C0(1+h u))] INSIDE the displaced
state (before restore) and compare with lhs (expected == to ~3e-4);
then bisect which state element (int2 re-init at reference vs
displaced, grid, screening) carries the h^0 shift. The answer defines
the final production term or exonerates the assembly entirely (with
d_num as the only valid referee).

### 7.47 v14: state-dependence ELIMINATED; the terminus sharpens to the
entrywise scaling of the reconstruction error
State-swap bisection: IN the displaced state |lhs - staged| = 0.019559;
in the reference state 0.019453; the state shift itself is only 6.2e-4.
NOT state-dependence. Therefore the entire residual is A's response to
e(h) = C+_signfixed - C0(1 + h u), whose MAX-norm scales as h^2 (M1)
while the response scales as h^1 (M2 constancy) => SOME ENTRIES of
e(h) must scale as h^1 (masked by larger h^2 entries in the max norm),
with O(1) gain. NEXT-SESSION MEASUREMENT #1 (final): print e(h) and
e(h/2) ENTRYWISE; the entries scaling as h^1 name the missing direction
-- their (p,q) orbital labels ARE the final term of the derivation.
(Candidate origin of an h^1 reconstruction error: the inverse-metric
correction (C0^T S_cross)^{-1} vs C0 in u's definition, second-order
fold/sg interplay, or SCF-convergence tails in specific blocks.)

### 7.48 v15: THE FINAL TERM IS MEASURED (and my 7.46 scaling misread
corrected)
CORRECTION: in 7.46 I misread "2.65e-4 -> 1.33e-4" as O(h^2); halving
is h^1. The entrywise map settles it:
  - 161/192 entries of e(h) = C+_sf - C0(1+h u) scale EXACTLY as h^1
    (ratio 0.500); only 31 as h^2.
  - RESPONSE ATTRIBUTION: perturbing C0 by the h^1-part of e alone
    reproduces THE ENTIRE residual: |dA|/h = 0.01959 vs target 0.0196;
    the h^2-part gives 0.00028. THE MISSING DIRECTION IS MEASURED:
    V := e(h)/h (an O(1) matrix, virtual-row dominated: rows 16,17,6
    + socc 5).
  - Consequence: u = (M - Sk)/h MISSES the O(1) direction V, despite
    the seemingly exact identity M - Sk = C0^T S_cross dC. The identity
    chain (with completeness C0 C0^T = S0^{-1}) predicts e = O(h^2) --
    so ONE of its inputs is not what it seems. NEXT-SESSION #1 (pure
    numerics, 5 min): verify M - Sk == C0^T S_cross_record dC entry by
    entry from the run's own overlap_ao record; wherever it fails names
    the convention (basis_overlap orientation / bfnrm / old-basis copy)
    behind V.
PRODUCTION MEANING: the U-channel contraction must use (Ux + V) -- and
once V's analytic form is identified from the convention fix, the
Mt-channel closes: Mt:(Ux+V) == exact was verified IMPLICITLY by the
response attribution. The end of the campaign is one convention-
identification away.

### 7.49 v16: THE CULPRIT -- the Sk record's DIAGONAL
Reconstruction from the run's own records:
  M record == C0^T S_cross Cd to 2.7e-7  (M construction exact).
  Sk record != C0^T S_cross C0: deviation 1.73e-4, DIAGONAL-
  CONCENTRATED, and the worst rows are EXACTLY V's dominant rows
  (17,17), (16,16), (6,6), (5,5)...
MECHANISM (complete): the spurious O(h)-ish diagonal in Sk =>
u = (M - Sk)/h carries a spurious O(1) DIAGONAL => e = h*(diag-err
rows of C0) -- exactly the measured V pattern (rows 16/17/6/5 spread
over AO columns). The ENTIRE invariant residual (2-13% channel, the
1.9e-2/8.7e-3/0.21 production gap) traces to THIS single diagonal
convention gap in the frozen-C cross-overlap.
REMAINING (one check + one decision): (a) is the Sk-record deviation
O(h) (a genuine convention error in get_structures' second-call path
-- then FIX it and u becomes exact) or O(h^2)*large-coefficient
(<phidot|phidot> ~ 170 for core functions -- then the TRUE Sk-diag
must use the product C0^T S_cross C0, not the record)? Measure
D(h) vs D(h/2). (b) Either way the PRODUCTION fix is identical:
build u from the RAW PRODUCT C0^T S_cross_record C0 (bypassing the
routine's Sk output), i.e., compute Sk in numpy from overlap_ao --
one-line change in every sweep/probe -- then re-run v7j/J1: expected
theory-level closure. THE CAMPAIGN'S LAST STEP IS A ONE-LINE FIX.

### 7.50 *** THE ONE-LINE FIX CLOSES THE U-CHANNEL *** (v17)
With Sk built as the RAW PRODUCT C0^T S_cross_record C0 (bypassing the
routine's second-call output):
  M1 reconstruction: 9.1e-7 (h) -> 2.3e-7 (h/2)  [pure O(h^2): exact]
  M2 |lhs - staged|: 0.00040 -> 0.00020          [collapsed 50x from
     0.0195; clean O(h) = the benign mixed truncation term]
THE ENTIRE "invariant residual" chased through 7.31-7.49 was the Sk
record's diagonal. The U-channel is CLOSED to FD-truncation level.
v18 (v7j + the fix in the sweeps) is running as the full re-judge:
expected theory-level J1/J3 closure on the complete assembly. The fix
must be propagated to every sweep-based referee (A8-lineage, v3-v15)
and noted as a LANDMINE: "OQP::overlap_mo_non_orthogonal" from a
second same-geometry call is NOT C0^T S_cross C0 on the diagonal --
always rebuild Sk from OQP::overlap_ao_non_orthogonal directly.

### 7.51 v18: J1 unchanged -- and rightly so; the deficit is CORNERED
INTO THE F-CHANNEL
v18 (full gate + the 7.49 fix): J1/J3/J2 identical to v7j to all
digits. This is CONSISTENT, not a failure: Mt's diagonal is exactly
zero, so the diagonal Sk correction cannot enter Mt:Ux. The v13-v17
discovery was a VECTOR-norm sub-mystery (a real landmine + real fix,
but ytil-projection-blind).
THE DECISIVE CONSEQUENCE of v17 stands: with the fixed u, the pure
C-channel VECTOR identity [A(C')-A(C0)]|sameF == staged(u) closed to
4e-4 (O(h) truncation). Therefore the remaining J1 deficit
(1.94e-2 / 8.7e-3 / 2.06e-1) can only live in the F-CHANNEL:
   [A(C0, F_d[D']) - A(C0, F_d[D_ref])] X    (true)
   vs  MTG:Ux == gvec(dD_model(Ux))          (model)
i.e., THE DENSITY-RESPONSE MODEL dD = C(U D0 + D0 U^T)C^T vs the true
Delta-D. Everything else in the chain is now closed at the vector
level. MEASUREMENT (next session #1): the worker-style 1-iter skel in
SEPARATE PROCESSES (the in-process emulation is broken, 7.44) gives
w_skel(F[D_ref]); [w_ref - w_skel - stagedC] = the true F-channel;
diff against gvec(dD(Ux)) names the dD-model error (candidates: the
ROHF occupation model, the XC response in gbuild vs the sigma's Fock
usage, spin-resolution of dD).

### 7.52 v19: the F-channel is 89% closed; the LAST OBJECT is the
fold-sector Fock response beyond the linear FOCK-record channel
Separate-process 1-iter workers (this process's phases) + fixed-Sk Ux:
  |w_ref|=0.066, |w_skel|=0.669, |stagedC|=0.651
  trueF = w_ref - w_skel - stagedC: |trueF|=0.0871
  gvec(dD_model): |gvec|=0.0914;  |trueF - gvec| = 0.0102 (11%)
  The residual concentrates ONCE MORE on LR1 (+0.0072) and the socc2
  rows -- the fold sector.
ACCOUNTING NOW COMPLETE TO ~1%: U-channel = stagedC (closed, 4e-4)
+ F-channel (dD-model + measured-G, 89%) + [0.010 fold-sector
residual]. The last named object: the sigma's SOCC/fold-sector Fock
usage responds to the density beyond the linear FOCK_A/FOCK_B-record
channel (candidates: the open-shell spc coupling's fa/fb combinations,
XC response pieces reaching the socc rows differently than gbuild's
record perturbation represents). This is a ~1%-of-channel effect
(~2e-2 on d for the worst pair), fold-sector-local, and is THE
remaining derivation item for theory-level closure. Everything else
in the U-channel is closed at FD-truncation level.

### 7.53 The antisym-dD derivation: proposed, tested, FALSIFIED (with
partial signal) -- the fold-term needs the sigma-source reading
PROPOSED (paper): if the worker skeleton were Lowdin-orthonormalized,
its density baseline D_g = D_ref - (h/2) C{S^x,D0}C^T would exactly
cancel the sym-U part of the dD model, leaving dD = C[U_a, D0]C^T.
TESTED offline (v7h npz): J1 with MTG:U_antisym --
  (1,3) improves 8.7e-3 -> 7.0e-3; assembly (2,3) 0.205 -> 0.120;
  BUT (1,2) worsens 1.9e-2 -> 1.3e-1 and J1(2,3) -> 0.54. FALSIFIED
as the global answer. Consistently: the G-A gate (engines == worker
w_skel at 1e-5..3e-4) already bounds any orthonormalization shift far
below the residual -- the premise was wrong; the sym-dD channel is
REAL under the actual referee convention.
STATE AT CLOSE: the fold-sector residual (7.52: 0.010 vector-level,
LR1+socc rows; -> 1.9e-2/8.7e-3/0.21 on d) remains THE item. Its
derivation requires reading the sigma source's SOCC/fold machinery
(mrsfmntoia / mrsfcbc spc channels, the U-matrix pairing of JCP 158,
194105, the fold application) -- a focused fresh-context task. All
referees (v19 slot-resolved, v7h/v7o matrices) are frozen and synced.

### 7.54 v20: the fold term is the missing XC-kernel Fock response;
derived and closed analytically

The sigma-source reading removes the ambiguity left in 7.52. In
`mrsf_matvec_apply`, `mrsfcbc -> int2_driver -> mrsfmntoia` constructs
the trial-vector two-electron part. Its SPC-scaled channels and the
`trans` pairing contain no ground-state `DM_A/B` or `FOCK_A/B`
dependency. The ONLY ground-state-Fock dependency of the sigma build is

```
call mrsfesum(infos, wrk1, fa, fb, amo, ivec)
```

Therefore the 7.52 residual cannot be a new `mrsfcbc`, `mrsfmntoia`,
SPC, or pairing channel. It must be a missing component of the response
Fock passed through the already-existing `mrsfesum` fold.

Let the reference AO densities be

```
P^s = C O_s C^T,                    s in {alpha,beta},
C^x = C U^x.
```

The first-order density used by the fold is consequently

```
P^{s,x} = C (U^x O_s + O_s U^{x,T}) C^T.
```

For a hybrid KS reference the complete response Fock is

```
F^{s,x}_resp = G^s_JK[P^{alpha,x},P^{beta,x}]
             + sum_t f_xc^{s,t} P^{t,x}.
```

The v19 `gbuild` probe perturbed only `OQP::DM_A/B` and called
`hf_energy`. Source inspection shows why that is incomplete:
`calc_jk_xc` builds J/K from `DM_A/B`, but `calc_dft_xc -> dftexcor`
reconstructs the XC density from `VEC_MO_A/B`. A DM-only perturbation
therefore differentiates J/K while leaving the XC potential unchanged.
The exact missing object is

```
Delta F^s_xc = sum_t f_xc^{s,t} P^{t,x},
Delta sigma_xc = R_X(Delta F^alpha_xc, Delta F^beta_xc),
```

where `R_X` is not a new operator: it is exactly the linear Fock map
implemented by `mrsfesum`. In its generic block it forms
`Xtilde F_beta - F_alpha Xtilde`; for the singlet LR1/LR2 sector it
adds the four `xlr/sqrt(2)` row/column corrections and overwrites LR1
with the corresponding signed open-shell contractions and diagonal
combination. Substituting `Delta F_xc` for `fa/fb` therefore predicts
both the magnitude and the LR1+socc localization seen in v19.

v20 (`v20_fold_audit.py`, H2O coordinate 5, private workers) separates
the channels and gates the analytic result:

```
|trueF|                                      = 8.704022e-2
actual displaced Fock -> mrsfesum:
  |trueF - g_actual|                        = 2.2264e-7
v19 DM-only JK response:
  |trueF - g_JK|                            = 1.016826e-2
XC correction alone:
  |g_xc|                                    = 1.016837e-2
  |(trueF-g_JK) - g_xc|                     = 2.6537e-7
analytic get_response_packed vs MO+DM FD:
  |G_alpha^analytic-G_alpha^FD|             = 2.0554e-10
  |G_beta^analytic-G_beta^FD|               = 2.0471e-10
analytic JK+XC -> mrsfesum:
  |trueF - g_analytic|                      = 2.6539e-7
```

This closes the entire F channel to the finite-difference/noise floor.
The largest old residual slot remains LR1 (`+0.00723675`) before the XC
correction, exactly as the fold formula predicts. A diagnostic C entry
point, `mrsf_nac_response`, now exposes the existing
`scf_addons:get_response_packed` JK+XC kernel through packed
`OQP::nac_dm1_{a,b}` -> `OQP::nac_v1_{a,b}` records. This is a gated
candidate path, not yet propagated into production `nac_analytic.py`.

BUILD NOTE: rebuilding exposed a pre-existing symbol collision between
the external Fortran routine `mrsf_nac_wpair` and
`bind(C,name="mrsf_nac_wpair")`. Renaming only the internal scaffold to
`mrsf_nac_wpair_impl` restores a clean build without changing behavior.

### 7.55 Full-response propagation and the converged-referee reversal

`mrsf_nac_response` was wired into `nac_analytic.py`: for every ordered
pair, the interstate alpha/beta densities are passed to
`get_response_packed`, transformed to MO space, and folded into
`MT_response`. The resulting reference formula is

```
d_IJ = antisym[T1 - seam(X) + X:V + gamma:Sk],
X = MT_frozen + MT_response + gamma.
```

Against the old H2O numerical referee (`dx=1e-3`) the apparent v21
closure was 5.370e-4, 5.354e-5, and 7.804e-4 for pairs (1,2), (1,3),
and (2,3). This verdict was premature: new `dx=5e-4` and `2.5e-4`
freezes agree pairwise to 3.08e-7, 1.32e-9, and 9.23e-6, whereas the old
`dx=1e-3` referee differs from `2.5e-4` by 3.223e-4, 5.229e-5, and
2.771e-3.

Against the converged `dx=2.5e-4` referee, v2 gives

```
pair       production seam     same-process direct U
(1,2)      5.8134e-4            3.1921e-5
(1,3)      1.9583e-6            1.9922e-5
(2,3)      3.4193e-3            2.6446e-4
```

Thus the full amplitude/F-response expression is sound; the remaining
error is introduced by the no-FD replacement of direct U with the
z-vector interchange seam. The decisive J4 ordered-pair mismatch is
6.551e-3 for (3,2), which becomes the 3.419e-3 antisymmetrized (2,3)
error. Forcing NAC-only MINRES and a 1e-10 residual changes none of the
pair errors, falsifying the earlier loose-CG hypothesis. The MINRES
change remains a useful precision guard, not the structural fix.

### 7.56 Resident Fortran wpair reference engine

The production-scale Python orbital-generator loop was moved into the
resident Fortran `mrsf_nac_wpair_impl` routine. The C/tagarray interface
accepts `nac_ytil` and `nac_xstate` and exports `nac_mt_frozen`; Python
now performs only one thin call per ordered pair. The rebuilt CFFI path
and H2O v25 gate reproduce the previous Python-harvest v2 coupling to
the following state-pair maxima:

```
(1,2)  1.2741e-9
(1,3)  1.1012e-10
(2,3)  3.7301e-10
```

This is a resident Fortran reference implementation, not the final
analytic kernel: it still evaluates central orbital generators inside
Fortran. The next implementation step is to replace that internal
O(nbf^2) harvest with the closed-form bilinear `mrsfcbc/mrsfmntoia`
adjoint while preserving the tested external interface. Before the
ethylene or Acrolein verdict, however, the larger seam/interchange
error isolated in 7.55 must be corrected and regated against the
converged small-displacement H2O reference.

### 7.57 Lee-gradient audit: the old seam mixed two normalizations

The paper/SI audit summarized at the start of this document changes the
interpretation of 7.55--7.56.  Let the independent ROHF rotation be

```
kappa = (kappa_SD, kappa_DV, kappa_SV),
SD = docc-socc, DV = docc-virt, SV = socc-virt.
```

The physical reference-stationarity residual used by native OpenQP is

```
r = (F_beta_SD, F_alpha_DV + F_beta_DV, F_alpha_SV).
```

Lee et al.'s `Fbar`, coefficient stationarity, and the one-half factors
in Eq. (3.8)/SI Eq. (S26) imply

```
H_native := dr/dkappa = 2 Jbar,
l_native := dG/dkappa = 2 Rbar.
```

The Lee Lagrange multiplier therefore obeys

```
H_native Z_Lee = -l_native,
Z_Lee = Zbar.
```

OpenQP's nuclear response is written `H_native U^R=B^R` with
`B^R=-partial_R r`.  Production therefore stores the sign-reversed
computational adjoint `zeta=-Z_Lee`, solves

```
H_native zeta = +l_native,
zeta^T B^R = l_native^T U^R,
```

and adds `zeta^T B^R` to the ordered derivative.  The code record is named
`nac_rohf_z`, but its sign is `zeta=-Zbar` in Lee's multiplier convention.
There is no extra *factor* of two.  In the legacy gradient solver the same
accounting is distributed differently: `sfrorhs` builds `-2 Rbar`,
`sfrolhs` applies `Jbar`, its solved vector is `xk=2 Zbar`, and `sfropcal`
inserts `xk/2`.  Feeding the computational `zeta=-Zbar` through that legacy
half-density seam, or feeding `xk=2Zbar` into a native adjoint contraction,
mixes sign and normalization conventions.
The v3 production path therefore bypasses `sfropcal/mrsfrowcal` for the
interstate adjoint and contracts the native solution directly.

The native tangent embedding and dual projection are also not mutual
array inverses.  In block form,

```
rohf_unpack_trial(kappa):
  SD -> beta only
  DV -> alpha = kappa and beta = kappa
  SV -> alpha only

rohf_pack_trial(g_alpha,g_beta):
  SD <- g_beta
  DV <- g_alpha + g_beta
  SV <- g_alpha
```

Thus the induced coordinate metric is

```
E^T E = diag(1,2,1) = pack(unpack(1))
```

over `(SD,DV,SV)`.  `unpack` is the tangent embedding `E`; `pack` is
the cotangent projection `E^T`.  In particular the DV factor two is
already present on both sides of the native stationarity equation and
must not be "corrected" afterward.  The old source comment calling the
two routines inverses is mathematically misleading.

The audit also fixed a non-square TagArray trap in the forward gate.
A Fortran record `A(ltot,ncart)` exposed through OQPData must be
recovered in Python as

```
raw.reshape(ncart, ltot).T
```

not as a bare `raw.T`.  Any earlier native-U result read with the latter
layout is void.

### 7.58 Native Z-vector production path and the DFT lifecycle bug

The production numerical route is now resident in Fortran, with Python
only orchestrating one call sequence per ordered state pair:

```
Xmat = MT_frozen + MT_response + gamma
Ldual = E^T (Xmat - transpose(Xmat))            [pack_rohf_dual]
H_native zeta = Ldual                          [one ROHF/ROKS solve]
zB_HF = zeta^T B^R_HF/JK/Pulay                 [analytic Fortran adjoint]
zB_XC = zeta^T B^R_XC                          [analytic Fortran adjoint]
d_IJ = T1 + zB_HF + zB_XC + Xmat:V + gamma:Sk
```

The entry points are `mrsf_nac_rohf_solve`,
`mrsf_nac_rohf_hf_adjoint`, and `mrsf_nac_xc_adjoint` in
`source/modules/mrsf_nac_interchange.F90`.  The solver uses the same
`cphf_apbx_rohf` operator and coordinates as `hf_hessian_rohf`, requires
a residual-converged solve, and performs one adjoint solve rather than
`3N` forward CPHF solves.  The HF routine transposes the nuclear
one-electron, two-electron, JK-response, and Pulay contractions.  The XC
routine evaluates the moving-grid mixed derivative and adds the
occupied-space reorthonormalization response through one
`f_xc[P_z]` build.

For the full symmetric physical rotation density

```
delta P_z^s = C_v^s z^s C_o^{s,T} + transpose,
```

the exact XC normalization is

```
Tr[V_xc delta P_z] = 2 z^T r_xc,
z^T b_xc = -1/2 d_R Tr[V_xc delta P_z].
```

Therefore the final `-0.5` in `mrsf_nac_xc_adjoint` is required.  An
early gate appeared to make the analytic XC result twice the forward
one, but the factor was not algebraic.  `dft_initialize` appends
functionals to process-global XC state; `cphf_solve_rohf` initialized
the grid and returned without `dftclean`, after which the XC adjoint
initialized it again.  BHHLYP was consequently present twice.  The
solver now brackets its DFT state with `dft_initialize/dftclean` (as do
the reusable RHF/UHF CPHF drivers), and the `-1/2` factor is retained.
Do not revive the rejected `-1/4` workaround.

The `mrsf_nac_wpair_impl` status also changed after 7.56: its central
orbital-generator harvest has been replaced by a closed-form
`mrsfcbc/mrsfmntoia` adjoint.  It uses a two-vector batched ERI build,
adds both spin-pair sides through `mrsfsp`, evaluates the frozen-Fock
term in MO space, and preserves the tested C/tagarray interface.  Thus
the expensive production kernels are Fortran-resident; Python retains
only state-pair orchestration and gate logic.

### 7.59 The apparent 7.55 seam failure was primarily loose SCF/TDHF

The former `3.419e-3` H2O verdict refined the nuclear displacement but
left both SCF and TDHF at `1e-6`.  It therefore tested a stable finite
difference of incompletely stationary orbitals, not the Lagrangian
limit.  This is especially damaging for the near-degenerate `(3,2)`
pair, whose large adjoint amplifies a small stationarity residual.

The same-process v33/v34 comparison isolates this effect.  At
`scf.conv=tdhf.conv=1e-6`, native CPHF `U` differed from the independently
reconverged finite-difference `U` by

```
SD 8.974e-3,  DV 1.158e-3,  SV 6.699e-3,
L.U max difference = 6.579e-3.
```

At `scf.conv=tdhf.conv=1e-10`, with the same geometry and algebra, these
became

```
SD 8.262e-6,  DV 3.164e-6,  SV 8.462e-6,
L.U max difference = 1.181e-5.
```

This three-orders-of-magnitude collapse identifies loose electronic
convergence as the main source of the old large "seam" error.  It also
explains why merely tightening the Z solver did not repair 7.55: the Z
equation was already solved, while its reference orbitals and numerical
referee were not stationary enough.  `nac_analytic.py` now refuses
`scf.conv > 1e-8` and recommends `1e-10` near a crossing.

Independent component gates on ordered pair `(3,2)` give, at tight
convergence,

```
native direct L.U vs native adjoint z.B       1.151e-5
analytic HF/JK/Pulay vs forward RHS           3.116e-8
analytic XC vs finite-difference XC RHS       1.70490e-4
full analytic Z contraction vs forward RHS    1.70474e-4
```

The remaining `1.7e-4` is the XC mixed-derivative/forward-FD gate level;
it is not evidence for a missing factor or a second MRSF sigma channel.

### 7.60 H2O v36 production verdict (tight reference)

The first end-to-end gate of the pairwise one-Z-vector production path
uses H2O/BHHLYP/6-31G*, singlet MRSF, `tlf=0`,
`scf.conv=tdhf.conv=1e-10`.  The analytic artifact is

```
~/nac_audit/probe/H2O_energy_tlf0_tight_v36_z.npz
```

and the numerical reference is

```
~/nac_audit/probe/H2O_energy_tlf0_tight_dx25e5_dnum.npz
dx = 2.5e-4 Angstrom.
```

Gauge-resolved maximum component errors are

```
pair (1,2)   2.10885792e-5
pair (1,3)   1.50539330e-5
pair (2,3)   2.10869153e-4
overall      2.10869153e-4   PASS (criterion 3.0e-4)
```

The tight numerical references at `dx=5.0e-4` and `2.5e-4` agree to
`2.49864104e-7` overall, so displacement truncation is well below the
production error.  One artifact trap must be recorded: the
`dcv_reference` array embedded in `H2O_energy_tlf0_tight_v36_z.npz`
and the accompanying v36 log still point to the older loose-convergence
freeze and display a misleading `3.27076e-3`.  The certified comparison
is the artifact's `dcv` against the separate tight `dx25e5` file named
above.

This is a real H2O production pass, not a general completion claim.
Tight-convergence ethylene/C1, Acrolein, translational/sum-rule gates,
and upstream duplicate-fix review remain.  More fundamentally, the Lee
energy-gradient paper does not prove the interstate bilinearization;
the H2O result validates this implementation and its diagonal
normalization on one system.  Broader systems and the explicit `I=J`
collapse remain the evidence required before calling the MRSF NAC
Lagrangian generally established.

### 7.61 Former v3 ordered-HST formulation (superseded by Section 8)

Sections 7.55--7.60 describe the route by which the implementation was
localized, but their H2O-only accuracy verdict is superseded below.  The
production observable is the central/Hammes--Schiffer--Tully (HST) derivative
of a state-overlap matrix.  Write the one-sided analytic derivative before the
state projection as `Dord_IJ`.  The definition shared by the numerical referee
and production code is

```
d_IJ = 1/2 (Dord_IJ - Dord_JI),
h_IJ = (E_J-E_I) d_IJ.
```

The zero residuals of `d+d^T` and `h-h^T` are therefore output-contract
checks, not independent proof of the interstate Lagrangian.  In particular,
the unprojected H2O v64 kernels have non-negligible symmetric parts.

For a fixed ket state `J`, differentiation of the MRSF state-overlap formula
gives the ordered kernel

```
Dord_IJ = Gmet_IJ^T X_J^R + gamma_IJ : (Sk^R + U^R).
```

Here `Gmet_IJ` is the derivative of the exact `tlf=0` overlap expression with
respect to the ket amplitudes; `gamma_IJ` is its orbital-overlap derivative;
`Sk` is the ket-half AO overlap derivative; and `U^R` is the MO response.  On
the physical folded configuration space, define the amplitude adjoint by

```
(Omega_J - A) y_IJ = Q_J Gmet_IJ,
Q_J = 1 - X_J X_J^T.
```

Symmetry of the MRSF/TDA response operator then eliminates the coordinate-wise
amplitude response:

```
Gmet_IJ^T X_J^R = y_IJ^T A^R X_J.
```

The derivative of this bilinear response contraction is split in the same way
as an analytic gradient: a frozen-orbital skeleton plus an orbital-gradient
source.  In the implementation,

```
T1   = [mrsf_nac_amp + mrsf_nac_esum](y_IJ, X_J),
M_IJ = MT_frozen(y_IJ,X_J) + MT_response(P_IJ) + gamma_IJ,
```

and the ordered kernel is

```
Dord_IJ = T1 + z_IJ^T B^R_HF/JK/Pulay + z_IJ^T B^R_XC
              + M_IJ : V^R + gamma_IJ : Sk^R.
```

`V^R` is the overlap-fixed symmetric/reorthonormalization part of the MO
derivative.  The independent ROHF rotations are eliminated by the *adjoint*
stationarity equation

```
H_native^T zeta_IJ = E^T (M_IJ-M_IJ^T),
```

where `E` is `rohf_unpack_trial` and `E^T` is the native dual projection used
by `pack_rohf_dual`.  The ROHF orbital Hessian is symmetric in these physical
tangent/dual coordinates, so the same matrix-action implementation used by a
CPHF driver may solve the adjoint equation.  This does not turn the production
algorithm into a `3N` forward CPHF calculation: production supplies exactly
one state-pair RHS and evaluates `z^T B^R` for all coordinates analytically.
The public production entry is consequently named
`mrsf_nac_rohf_zvector`; `cphf_solve_rohf` is only the historical name of the
reused orbital-Hessian linear solver.

The Lee-gradient normalization fixes the remaining apparent factor ambiguity.
For native OpenQP stationarity coordinates,

```
H_native = 2 Jbar,
l_native = E^T (M_II-M_II^T) = +2 Rbar,
(2 Jbar) zeta = +2 Rbar  =>  zeta = -Zbar.
```

The legacy diagonal chain instead solves `Jbar xk=-2Rbar`, obtains
`xk=2Zbar`, and inserts `xk/2` in `sfropcal`.  The sign-reversed native
computational adjoint `zeta=-Zbar` must never be sent through that legacy
half-density seam.  Likewise,
`pack(unpack(kappa))=diag(1,2,1)` over `(SD,DV,SV)` is the correct induced
metric, not a missing factor that should be repaired after the solve.

### 7.62 The actual residual source: XC fuzzy-cell moving-grid response

After tight SCF/TDHF and the native Z-vector seam had been established, the
remaining ordered-pair error localized completely to the XC part of
`mrsf_nac_esum`: the one-electron and two-electron skeletons agreed with
finite differences to `5e-8` and `2e-8`, while the XC skeleton differed by
`1.7974e-3`.  A zero-probe subtraction could remove an accidentally included
ground-state term, but it could not differentiate the atom-centred quadrature
itself.

For a linear interstate probe `P`, OpenQP evaluates a finite quadrature

```
E_xc[P;R] = sum_g w_g(R) q_P(r_g(R);R).
```

The consistent derivative contains three pieces:

```
dE_xc = sum_g w_g dq_P|basis
      + sum_g q_P dw_g|partition
      + sum_g w_g grad_r(q_P) . d r_owner.
```

The last two terms are a pair: every atom-centred slice point moves with its
owner, so retaining only the normalized partition derivative or only the
integrand owner motion produces a large spurious result.  For fuzzy cells
`p_O=c_O/sum_K c_K`, the implemented logarithmic form is

```
d log p_O = d log c_O - sum_K p_K d log c_K.
```

Only the point owner and the two atoms of a Becke pair affect a given `mu_ij`.
Precomputed `R_ij`/unit vectors and thread-local logarithmic derivatives reduce
the new work from `O(Ngrid*Natom^3)` to `O(Ngrid*Natom^2)` and avoid inner-loop
allocation.  Production invokes the linear-probe branch with
`include_ground_state=.false.` and `include_weight_derivative=.true.`.  The
routine now accumulates into `dedft` (`intent(inout)`), so the ordinary
state-diagonal MRSF gradient is not overwritten.

The decisive ethylene pair `(2,3)` XC check changed from `1.7974e-3` to
`6.5165e-8`.  This supersedes the `1.7e-4` "XC floor" in 7.59: that value was
not an unavoidable grid error but a missing analytic moving-grid term.

### 7.63 Historical v65/v60/v62 numerical reference gates

All gates below use independently reconverged displaced workers, one process
per displacement, `scf.conv=tdhf.conv=1e-10`, and a central step of
`2.5e-4 Angstrom`.  One state gauge is solved globally; pair signs are not
chosen independently.

```
system / artifact                    pair       max component error
H2O v65                              (1,2)      4.58748968e-8
                                     (1,3)      2.15320715e-7
                                     (2,3)      1.12841328e-6
ethylene v60                         (1,2)      2.70646099e-6
                                     (1,3)      1.25137596e-7
                                     (2,3)      7.06993091e-6
Acrolein v62, tlf=2                  (1,2)      2.13691099e-6
                                     (1,3)      2.30623982e-6
                                     (2,3)      2.22002911e-6
```

The corresponding maximum energy mismatches are `0`, `2.416e-13`, and
`2.274e-13 Hartree`.  Every pair passes the deliberately loose publication
gate `3e-4`; the observed production errors are two orders of magnitude
smaller.  The certified artifacts are

```
~/nac_audit/probe/H2O_production_optimized_v65.{npz,out}
~/nac_audit/probe/ETH_production_movinggrid_v60.{npz,out}
~/nac_audit/probe/Acrolein_production_movinggrid_v62.{npz,out}
```

The frozen references are the matching `*_dx25e5_dnum.npz` files and their
worker directories.  `nac_reference_gate.py` is the reproducible comparison;
it also checks the returned-array contracts, but those construction identities
must not be counted as separate theory evidence.

### 7.64 Lifecycle and diagonal-gradient regressions

The ordered pair engine must be independent of what ran earlier in the same
process.  `mrsf_nac_amp` previously reused a resident `OQP::td_p` created by a
gradient call; it now owns a zero two-particle-response slot for the NAC
skeleton.  Together with the paired `dft_initialize/dftclean` lifecycle, this
gives exact history invariance:

```
H2O v64: energy -> NAC  versus  energy -> gradient(root 3) -> NAC
maximum gauge-resolved difference over every pair = 0.0.
```

The ordinary Lee MRSF gradient is unchanged by the new accumulation and grid
metadata plumbing.  Against the frozen H2O reference, the final binary gives

```
ground/total energy max difference = 7.1054e-14
TD energy max difference           = 1.2698e-14
gradient max difference            = 1.6048e-13.
```

`OQP_NAC_SELFTEST` additionally proves that the bilinear 2e derivative engine
collapses bit-for-bit to the production quadratic 2e engine at `I=J`
(`max difference=0`, production norm `9.6123e-1`).  Its scope is deliberately
narrow: it does not by itself prove the full esum, native RHS/Z normalization,
HF/XC adjoints, or `Vmask` diagonal continuation.

The meaningful Lee limit is not the returned `d_II`, which is zero by the
real-state/HST convention.  Let `Delta_IJ=Omega_J-Omega_I`.  The overlap
derivative gives `Gmet_IJ -> X_I`, hence

```
Delta_IJ y_IJ -> X_I,
Delta_IJ Dord_IJ = Q_R[X_I,X_J] + Delta_IJ A_R[gamma_IJ].
```

The diagonal continuation is `Q_R[X_I,X_I]=d Omega_I/dR`, with
`gamma_II=0`; the total Lee gradient is `dE_0/dR+dOmega_I/dR`.  Frozen
off-diagonal artifacts satisfy `max|Delta*y-X_I|=1.10e-6` for H2O and
`6.16e-6` for the ethylene `(2,3)` pair, consistent with the current
`EPSA=1e-6` state-space unit sweep.  This algebra establishes the required
limit.  A duplicated-slot H2O v71 source gate additionally sets `y=X_I`,
`gamma=0`, and evaluates the pair source with two distinct but identical
amplitude slots.  It verifies the sign-correct Lee relation

```
ell_pair + rhs_legacy = 0,
ell_pair = +2 Rbar,       rhs_legacy = -2 Rbar,
```

with state maxima `3.22e-8`, `1.60e-7`, and `3.93e-7`
(`H2O_diagonal_rhs_all_v71.npz`).  This closes the source normalization and
the fact that the computational native adjoint is `zeta=-Zbar`.  A full
duplicated-slot *value* comparison against `dOmega_I/dR` remains a useful
future diagnostic; the present 2e self-test and v71 source gate must not be
described as that full Eq. (3.21) comparison.

### 7.65 Raw translation rule, not an ETF zero-sum rule

The current result is a raw electronic derivative coupling and contains no
electron translation factor (ETF).  Therefore `sum_A d_IJ^A=0` is the wrong
gate.  The correct implementation identity is

```
sum_A d_IJ^A = antisym_IJ sum_A (gamma_IJ : Sk_A).
```

On the H2O v64 debug decomposition, the atom sums of antisymmetrized `T1`,
`z_HF`, `z_XC`, and `Vmask` are each below `3.1e-14`; the full raw sum equals
the `gamma:Sk` sum to `4.1e-14`.  The raw atom sums themselves are nonzero:

```
H2O       1.1084e-1
ethylene  4.9231e-2
Acrolein  1.7962e-1.
```

`translation_gate.py` encodes this structural rule and deliberately never
subtracts an atomic mean or imposes zero.  Such a subtraction is not an ETF.
An ETF-corrected observable would require a separately derived mode and
separate references.

### 7.66 Partition-function audit and performance recheck

Adding the analytic weight derivative exposed dormant derivative/mapping
errors outside the default SSF path.  The Fortran type IDs are
`SSF=0, ERF=1, BECKE4=2`; Python had ERF and Becke reversed.  The ERF
derivative lacked one factor `1/(1-x^2)`, smoothstep orders 2--5 omitted the
chain factor `-SCALEF`, and the smoothstep-5 dispatch limit was `0.74` instead
of its definition's `0.73`.  These are corrected before the weight derivative
is exposed to production.

An independent H2O ERF grid gate gives maximum pair errors
`4.5806e-8`, `2.1533e-7`, and `1.1294e-6`, with zero energy mismatch
(`H2O_erf_dx25e5_dnum_v67.npz` versus `H2O_erf_production_v68.npz`).  The
optimized default H2O result reproduces v61.  The optimized Acrolein v66 rerun
reproduces v62 to `5.12e-11` maximum component and independently retains the
`2.30624e-6` numerical-reference error; its wall time was 1328.64 s on chc3.

### 7.67 Historical production boundary before resident finalization

At this historical checkpoint the following statement was supported:
**static singlet MRSF-TDDFT analytic raw NAC was implemented and independently
gated for H2O, C1 ethylene, and Acrolein.**  It was not correct to call the
entire OpenQP SOC-NAMD stack production-ready:

* analytic NAC currently rejects multiplicities other than singlet;
* `md` runtype is not implemented;
* trajectory-continuous state permutation/near-degenerate subspace rotation
  is not yet applied consistently to NAC, SOC, gradients, and coefficients;
* that checkpoint's driver still retained Python state-pair algebra, the
  `G_met`/`gamma` construction, MINRES, and small contractions.  Section 8.6
  supersedes this item: the final production pair path is resident Fortran.

For any future OpenQP numerical work, implement the numerical kernel in
Fortran first and keep Python to API orchestration and validation.  A future
NAMD shakedown must replay consecutive Acrolein frames, establish one causal
state/subspace gauge, compare `tau_IJ=sum_A v_A.d_IJ^A` with the time-overlap
coupling at two timesteps, remove COM velocity for the raw convention (or add
a separately validated ETF mode), and only then exercise hopping dynamics.

---

## 8. Final normative Lagrangian and production contract

This section supersedes every formula and implementation-status statement in
Sections 0--7.67.  It separates three operations that were repeatedly mixed
during the campaign:

1. differentiation of the one-sided, column-normalized MRSF state overlap;
2. symmetric polarization of the Lee diagonal excitation-gradient source;
3. adjoint elimination of the ROHF/ROKS orbital response.

The first makes an ordered HST leg, the second constructs its interstate MRSF
source, and the third is the one-RHS Z-vector.  Only after both ordered legs
have been evaluated is the HST state-index projection applied.

### 8.1 Ordered HST derivative and state response

Let `S_IJ(R0,R)` be the exact MRSF state-overlap formula with the bra at the
reference geometry and the ket at `R`.  Define its one-sided analytic
derivative

```
Dord_IJ^R = [d S_IJ(R0,R) / dR]_(R=R0).
```

The observable shared by the numerical referee and the analytic code is

```
d_IJ^R = 1/2 (Dord_IJ^R - Dord_JI^R),
h_IJ^R = (Omega_J-Omega_I) d_IJ^R.
```

Thus `d+d^T=0` and `h-h^T=0` are returned-array contracts.  They are not
independent evidence that either ordered leg is complete.

At the reference geometry the normalized MRSF/TDA eigenvectors obey

```
A X_K = Omega_K X_K,                 X_I^T X_J = delta_IJ.
```

For `I != J`, differentiating the eigenproblem and projecting on `X_I` gives

```
X_I^T X_J^R = X_I^T A^R X_J / Delta_IJ,
Delta_IJ = Omega_J-Omega_I.
```

The general state-space adjoint notation is

```
(Omega_J-A)^T y_IJ = Q_J Gmet_IJ,       Q_J=1-X_J X_J^T.
```

The exact overlap metric has `Gmet_IJ=X_I` at the identity.  Since `A` is
symmetric, `I != J`, and `Q_J X_I=X_I`, the spectral solution reduces exactly
to the production state-response vector

```
y_IJ = X_I / Delta_IJ
```

with the redundant folded slot set to zero.  No coordinate-wise Davidson
response solve is needed.  The small-gap division is nevertheless physical:
the driver rejects a zero or numerically unresolved `Delta_IJ`; it does not
pretend that an isolated-state NAC is well defined at an exact degeneracy.

The one-sided chain rule is then

```
Dord_IJ^R = y_IJ^T A^R X_J + gamma_IJ : (Sk^R + U^R),
```

where `Sk^R=<chi|d_R chi>` is the ket-half AO-overlap derivative in the MO
basis, `U^R` is the MO coefficient response, and `gamma_IJ` is the orbital
derivative of the exact state-overlap formula.  `gamma_IJ` is ordered in the
state labels and must not be replaced by `-gamma_JI`.

### 8.2 Symmetric interstate polarization

Lee Eq. (3.21) supplies a diagonal excitation-gradient functional.  If its
homogeneous amplitude-quadratic part is denoted by

```
Q_R[X] = X^T A^R X,
```

the unique real symmetric bilinear continuation used here is

```
B_R[y,X] = 1/2 { Q_R[y+X] - Q_R[y] - Q_R[X] }.
```

Equivalently, if a diagnostic evaluates the *full* gradient
`G_R[X]=G_R[0]+Q_R[X]`, it must use

```
B_R[y,X] = 1/2 { G_R[y+X] - G_R[y] - G_R[X] + G_R[0] }.
```

Production does not obtain this bilinear by subtracting four complete
gradients.  Its Fortran kernels construct the two-vector terms directly.  In
particular, every quadratic MRSF channel is symmetrized between the two slots;
the closed `wpair` source contains

```
1/2 [ H(y,KX) + H(X,Ky) ],
```

and the two-electron and spin-pair density channels use the corresponding
`1/2(yX+Xy)` products.  A left-slot-only or right-slot-only continuation is
not the implemented Lagrangian.

This interstate polarization is a new derivation.  Lee et al. prove the
diagonal functional and its stationarity, not this off-diagonal continuation.
The mandatory diagonal identity is `B_R[X,X]=Q_R[X]`.  The literal real-state
self-overlap derivative has `gamma_II=0` and the returned HST `d_II=0`; the
Lee value test therefore uses two distinct, duplicated amplitude slots and
then takes their diagonal continuation.  Those are different statements and
must not be conflated.

### 8.3 Exact streamed state-overlap metric

The production metric is the analytic first derivative at the identity of the
literal `ndtlf=0` MRSF overlap formula.  Before normalization, its seven
determinant-block contractions are products of the `s_ij`, `s_ab`, and `s_ia`
minor families, including the two SOMO spin-pair blocks.  For raw overlaps
`R_KJ`, the formula normalizes each ket column:

```
n_J  = [sum_K R_KJ^2]^(1/2),
S_IJ = R_IJ / n_J,

dS_IJ/dR_KJ = delta_IK/n_J - R_IJ R_KJ/n_J^3.
```

Production first reverses that normalization, accumulates the analytic
product-rule sensitivities of all seven contractions, and contracts them with
the exact determinant cofactors at the identity.  Direct cofactors are
required: some relevant minors are singular, so an inverse-based Jacobi
formula is not valid.  For the independent generator

```
K_pq=+1, K_qp=-1  (p>q),
```

half of the directional derivative is placed in each orbital slot,

```
gamma_IJ(p,q) =  1/2 dS_IJ/dtheta_pq,
gamma_IJ(q,p) = -gamma_IJ(p,q).
```

Orbital-slot antisymmetry is exact; state-index antisymmetry is deliberately
not imposed.  Same-space blocks are retained.  A fixed ket column `J` shares
one normalization denominator, so `mrsf_nac_metric_column` streams only
`O(nstate*nbf^2)` storage rather than materializing
`O(nstate^2*nbf^2)` data.

The independent gates are:

```
resident streamed column vs closed cofactor oracle     5.219e-13
exact-generator residual by block:
  doc-doc  2.385e-18      doc-socc 6.476e-13      doc-virt 0
  socc-socc 1.110e-15     socc-virt 8.327e-13     virt-virt 5.551e-13
orbital-slot antisymmetry residual                     0
same-space signal (must not be projected away)         7.0910e-1
ordered-state non-antisymmetry signal                  7.0111e-1
```

The old raw-determinant `gamma`, sign-scanned TLF kernel, and same-space-zero
claims therefore have no production role.

### 8.4 Pair orbital source and the complete ordered formula

Symmetric polarization of `y_IJ^T A^R X_J` produces an explicit nuclear
skeleton and an orbital derivative.  In the current Fortran decomposition,

```
T1_IJ^R = mrsf_nac_amp(y_IJ,X_J)^R
        + mrsf_nac_esum(y_IJ,X_J)^R,

M_IJ = MT_frozen(y_IJ,X_J)
     + MT_response(PairDensity_IJ)
     + gamma_IJ.
```

`MT_frozen` contains the closed `mrsfcbc`/ERI/`mrsfsp` bilinear and the
frozen-Fock derivative.  `MT_response` is the full JK plus XC response to the
pair density and includes the semilocal `f_xc` response through
`get_response_packed`.  It is not reconstructed from Python AO records.

Let `V^R` denote the dependent symmetric/reorthonormalization part of the MO
response fixed by `U^R+U^{R,T}=-S^R`.  Let `E` be the native ROHF tangent
embedding and `E^T` its dual projection.  The independent orbital source is

```
ell_IJ = E^T (M_IJ-M_IJ^T).
```

After adjoint elimination (Section 8.5), the complete ordered derivative is

```
Dord_IJ^R = T1_IJ^R
          + zeta_IJ^T B_HF/JK/Pulay^R
          + zeta_IJ^T B_XC^R
          + M_IJ : V^R
          + gamma_IJ : Sk^R.
```

There are no fitted factors or post-solve block repairs in this expression.
The resident accumulator stores this as `OQP::nac_dp_ordered`; the finalizer
forms

```
OQP::nac_dcv  = 1/2 (Dord-Dord^T_state),
OQP::nac_nacv = (Omega_J-Omega_I) OQP::nac_dcv.
```

### 8.5 Why this is a Z-vector, not `3N` CPHF

For each nuclear coordinate, the forward ROHF/ROKS response convention is

```
H_native U^R = B^R,                 B^R = -partial_R r,
```

where the independent stationarity residual is

```
r = (F_beta_SD, F_alpha_DV+F_beta_DV, F_alpha_SV).
```

A forward CPHF implementation would solve this equation for every one of the
`3N` columns and then form `ell_IJ^T U^R`.  The adjoint interchange instead
solves once per ordered state pair,

```
H_native^T zeta_IJ = ell_IJ,
```

and evaluates all coordinates as

```
ell_IJ^T U^R = zeta_IJ^T B^R.
```

The physical ROHF tangent/dual Hessian is symmetric, so production reuses the
same matrix action and the symmetric-indefinite MINRES route.  The internal
routine name `cphf_solve_rohf` is historical; `nrhs=1` and the surrounding
algorithm are an adjoint Z-vector calculation.  The only `3N` forward solve
is an explicit diagnostic gate.

Lee Eqs. (3.6)--(3.10) and SI S25/S26 fix the sign.  In OpenQP native
coordinates,

```
H_native = 2 Jbar,
ell_II   = +2 Rbar,
Jbar Zbar_Lee = -Rbar,

H_native zeta = ell  =>  zeta = -Zbar_Lee.
```

This computational `zeta` is what contracts with `B^R`.  The legacy diagonal
chain distributes the same normalization differently:

```
sfrorhs = -2 Rbar,
Jbar xk = -2 Rbar  =>  xk = 2 Zbar_Lee,
sfropcal inserts xk/2 = Zbar_Lee.
```

Consequently the diagonal cross-seam closure is

```
zeta + xk/2 = 0,
```

not `zeta-xk/2=0`.  Passing `zeta` through `sfropcal`, or adding another
factor to the DV block, is wrong.  The native maps already satisfy

```
pack(unpack(kappa)) = diag(1,2,1) kappa       over (SD,DV,SV),
```

which is the induced tangent/dual metric rather than a missing normalization.

### 8.6 Resident Fortran implementation map

The final scientific call graph is:

```
pyoqp analytic_nac
  -> mrsf_nac_lagrangian                         [one C call]
     -> mrsf_nac_metric_column                   [exact streamed gamma]
     -> mrsf_nac_wpair_impl                      [closed frozen source]
     -> mrsf_nac_amp + mrsf_nac_esum             [explicit skeleton]
     -> mrsf_nac_response                        [JK+fxc pair response]
     -> mrsf_nac_rohf_pair_overlap               [ell, V, gamma:Sk]
     -> mrsf_nac_rohf_zvector                     [one pair RHS]
     -> mrsf_nac_rohf_hf_adjoint + _xc_adjoint   [all coordinates]
     -> mrsf_nac_pair_accumulate + _finalize     [HST and gap]
```

The owning files are:

```
source/modules/mrsf_nac_driver.F90
source/modules/mrsf_nac_metric_data.F90
source/modules/mrsf_nac_interchange.F90
source/modules/tdhf_mrsf_gradient.F90
source/modules/tdhf_mrsf_energy.F90
pyoqp/oqp/library/nac_analytic.py                [scope/API/reshape only]
```

Production never materializes the all-pair metric, never sweeps orbital
generators, never performs pair algebra in Python, and never solves forward
CPHF for every nuclear coordinate.  Python scripts under
`tools/nac_lagrangian/` are diagnostic referees only.

### 8.7 XC moving-grid derivative for a linear relaxed-density probe

This was the last substantive diagonal-Lagrangian error.  Let the reference
spin densities be `D=(D_a,D_b)` and let the relaxed linear probe be
`P=(P_a,P_b)`.  On a GGA/meta-GGA grid point define

```
rho_P   = (rho_Pa, rho_Pb),

sigma_P = ( 2 grad(rho_a).grad(rho_Pa),
            2 grad(rho_b).grad(rho_Pb),
              grad(rho_Pa).grad(rho_b)
            + grad(rho_a).grad(rho_Pb) ),

q_P = e_rho.rho_P + e_sigma.sigma_P + e_tau.tau_P.
```

There is no extra `1/2` in `sigma_P`.  For LDA the sigma/tau terms vanish;
for GGA tau vanishes.  In the standard Lee MRSF gradient the probe is the
relaxed one-particle density

```
P = T + Z                                      [Lee Eq. (3.16)].
```

Here `Z` is Lee's relaxed-density contribution `Zbar_Lee`, inserted by the
legacy diagonal path as `xk/2`.  It is **not** the sign-reversed computational
adjoint `zeta`; in the matched diagonal convention
`Zbar_Lee=xk/2=-zeta`.

OpenQP calls `utddft_xc_gradient` for this MRSF term without `xa/xb`, hence
`doFxc=.false.`.  It is the linear `grad_v_xc` branch.  The separate
transition-density `grad_f_xc` branch is not part of this call.  The required
`f_xc` response of the ROHF orbital Hessian and interstate pair source is
already evaluated in `utddft_fxc/get_response_packed`; inventing an additional
`X f_xc X` moving-grid term in the MRSF gradient would double count a different
object.  The current moving-grid extension intentionally aborts if `xa/xb`
are present because that general third-derivative case has not been derived.

For an atom-centred slice owned by atom `O`, write its finite quadrature weight
as `w_g=w_g^base p_O`, with normalized fuzzy cell

```
p_O = c_O / sum_K c_K.
```

The exact derivative of the *discrete quadrature used by the energy* is

```
D_A E_xc[P] = (fixed AO/basis contribution)
            + sum_g w_g [ delta_AO grad_r(q_P)
                         + q_P D_A log(p_O) ],

D_A log(p_O) = D_A log(c_O) - sum_K p_K D_A log(c_K).
```

The owner-motion term and normalized-partition term are inseparable.  Keeping
only one gives a large spurious contribution.  Only the point owner and the
two atoms in each Becke pair enter a local `D_A mu_ij`; logarithmic cell
derivatives therefore give `O(Ngrid*Natom^2)` work and no inner-loop
allocation.

The production use is deliberately split:

* interstate `mrsf_nac_esum` and `mrsf_nac_xc_adjoint` request the complete
  linear-probe derivative with `include_ground_state=.false.` and moving-grid
  response enabled;
* the ordinary diagonal gradient first performs its established ground plus
  fixed-grid relaxed-density sweep, then performs a second resident
  `weight_derivative_only` sweep using a copy of `P=T+Z`, with the ground state
  excluded.  That second sweep adds only owner motion plus partition response
  and cannot double count the AO/basis term.

Allowing the correction while `include_ground_state=.true.` is forbidden.
The naive mixed call translates a ground-state integrand without its matching
ground partition derivative and was the v89 failure.

### 8.8 v84--v90 diagonal-value evidence

The following H2O/BHHLYP diagonal continuation compares the duplicated-slot
pair Lagrangian with the Lee analytic excitation gradient.  Values are maximum
Cartesian-component errors in hartree/bohr for states 1--3.  This diagnostic
sets `gamma=0`, as required by the literal diagonal metric limit; consequently
it does not test the exact metric of Section 8.3.  The phrase "exact metric" in
the v84 row identifies the binary/campaign stage only.

```
gate    change                                      state 1       state 2       state 3
v84     exact metric + native MINRES                 2.8941616e-5  8.4319324e-6  1.2237472e-5
v85     force exact zeta=-xk/2                      ~2.8938569e-5 ~8.4575500e-6 ~1.2217980e-5
v86     tighten solve/electronic thresholds         ~2.8941210e-5 ~8.4380600e-6 ~1.2264900e-5
v87     disable only pair-esum grid derivative       7.0430381e-6  1.4347500e-6  2.3365700e-6
v88     pure HF, removing XC entirely                3.9063379e-6  6.3663000e-8  3.5826000e-7
v89     naive ground+probe moving-grid call           approximately 5.33e-1 for every state (REJECTED)
v90     correction-only moving grid for P=T+Z        6.6297107e-7  3.6219095e-8  6.6150102e-8
```

The unchanged v84--v86 errors exclude the Z solver and ordinary convergence
as the limiting cause; the exact metric is absent from this gate by
construction.  v87 and the pure-HF v88 gate localize the defect to XC
quadrature response.  v89 proves that a ground and probe mixture is not the
derivative being sought.  v90 implements the correction-only `P=T+Z`
derivative and passes the `5e-6` diagonal-value gate.
Its `max|zeta+xk/2|` is `4.20265769e-7`; the pure-HF v88 closure is
`4.61e-11`.

The v85/v86 entries are the rounded campaign log values; the v84, v87, v88,
and v90 entries are the retained gate values used for the final decisions.
The final pointer-hardened binary rerun is retained as
`H2O_diagonal_value_v94_final.npz`; it reproduces the worst Lee-value error
`6.62971067e-7` and native/legacy Z closure `4.20265769e-7`, both below the
final `5e-6` gates.

The diagnosis is therefore precise: the ROHF/ROKS Hessian was not the main
problem, and no missing antisymmetric MRSF fold term was required.  The error
sat at the interface between Lee's relaxed density `P=T+Z` and the moving
atom-centred XC quadrature.

### 8.9 Final independent NAC reference gates

Independent references reconverge every displaced geometry in a separate
process, use one global state gauge, and never choose signs pair by pair.  The
current resident formula gives the following maximum component errors.  The
`tlf` label in this table belongs to the finite-displacement HST *reference*
calculation.  The analytic production metric itself always uses the exact
`ndtlf=0` identity derivative described in Section 8.3; it has no separate
`tlf=2` production branch.

```
system / final gate              pair (1,2)      pair (1,3)      pair (2,3)
H2O v94, tlf=0                   7.68935260e-8   2.09421460e-7   2.30391559e-6
C1 ethylene v94, tlf=0           3.21820655e-6   1.97306408e-7   7.55836774e-6
Acrolein v94, tlf=2              2.13672772e-6   2.30452788e-6   2.21028518e-6
```

H2O and C1 ethylene pass their final `1e-5` component gates; Acrolein passes
its final `2e-5` component gate.  The final artifacts are
`H2O_final_fortran_v94.npz`, `ETH_final_fortran_v94.npz`, and
`Acrolein_final_fortran_v94.npz`.  The earlier v80 resident accumulator and
finalizer reproduced the v76 H2O assembly exactly, proving that moving the
last pair algebra from Python to Fortran did not change the observable.  The
frozen NPZ files and their separate-process worker directories under
`~/nac_audit/probe/` remain the numerical referees; historical v65/v60/v62
values in Section 7.63 and the v76--v78 closeout labels are superseded by the
final-binary v94 runs.

For the raw electronic coupling, translation is not a zero-sum constraint.
Without an electron-translation factor (ETF), the exact code identity is

```
sum_A d_IJ^A = antisym_IJ sum_A gamma_IJ : Sk_A.
```

Subtracting an atomic mean would hide an error and is not an ETF.

### 8.10 Diagonal Lee checks and what they prove

The final diagonal evidence is layered:

* the bilinear two-electron engine collapses bit-for-bit to the quadratic
  Lee engine at duplicated equal slots;
* the v71 source test gives `ell_pair=+2Rbar` and
  `rhs_legacy=-2Rbar`, with state maxima `3.22e-8`, `1.60e-7`, and
  `3.93e-7`;
* the native/legacy multiplier test gives `zeta=-xk/2=-Zbar_Lee` within its
  stated solver/source tolerance;
* v90 identifies and closes the complete duplicated-slot excitation-gradient
  value, including the previously missing XC moving-grid response of `P=T+Z`;
  the final binary repeat `H2O_diagonal_value_v94_final.npz` passes at
  `6.62971067e-7` worst value error and `4.20265769e-7` Z closure.

These checks establish the computational sign, source normalization, the
diagonal normalization of the chosen symmetric continuation, and its full
excitation-gradient diagonal value for the tested H2O/BHHLYP case.  They do
not by themselves establish the off-diagonal continuation or turn the Lee
energy-gradient paper into an off-diagonal NAC proof; that support comes from
the derivation in Sections 8.1--8.5 and the independent H2O, C1 ethylene, and
Acrolein interstate gates in Section 8.9.

### 8.11 Supported scope and remaining limitations

The supportable statement is:

> OpenQP contains a resident-Fortran, static, singlet, two-SOMO
> ROHF/ROKS MRSF-TDA analytic raw electronic NAC, using one adjoint Z-vector
> RHS per ordered state pair, independently gated on H2O, C1 ethylene, and
> Acrolein, with its Lee diagonal continuation gated on H2O/BHHLYP.

The boundaries are equally important:

* UMRSF and multiplicities other than singlet are rejected by the production
  driver, even though some lower-level metric/bilinear kernels have broader
  diagnostic support.
* This is a raw electronic derivative coupling.  No ETF, SOC-NAC combination,
  `md` runtime, hopping dynamics, or trajectory-continuous state/subspace
  gauge is supplied by this derivation.
* The isolated-pair `1/Delta_IJ` form is not a near-degenerate subspace
  treatment.  Exactly or numerically unresolved gaps are rejected; physically
  clustered states require a separately derived subspace connection.
* The exact metric is the first derivative at the identity of the literal
  column-normalized `ndtlf=0` overlap formula.  It contains no production
  nuclear finite difference, but it is not a finite-displacement all-order
  overlap propagator.  The Acrolein v94 gate establishes end-to-end
  compatibility with a `tlf=2` numerical HST reference for that tested system;
  it is neither a separate analytic `tlf=2` implementation nor a general proof
  of every finite-displacement tracking mode.
* The implemented XC moving-grid extension is restricted to the linear-probe
  `doFxc=.false.` branch.  A future `xa/xb` third-functional-derivative moving
  grid needs its own derivation and gates.
* The final moving-grid value evidence is the named BHHLYP/H2O GGA gate (with
  the historical H2O ERF partition check in Section 7.66).  The LDA/GGA/meta-
  GGA linear-probe algebra is implemented generically, but this document does
  not claim an independent diagonal-value reference for every functional,
  quadrature, or basis.
* Lee et al. establish the diagonal energy gradient.  The interstate
  symmetric polarization and ordered exact-overlap metric are new work and
  must be cited and validated as such.
* Production numerical kernels belong in Fortran.  Python remains a thin API,
  reshaping, orchestration, and independent validation layer only.

### 8.12 Superseded claims retained only for forensic value

The following archived conclusions must not re-enter the implementation:

```
raw determinant gamma is the production metric                 false
the exact TLF metric has no same-space content                  false
gamma_IJ may be replaced by -gamma_JI before HST projection    false
an antisymmetric-dD fold term closes the F channel              falsified by v19/7.53
production requires orbital-generator finite differences       false
production requires 3N forward CPHF                             false
the remaining XC error is an irreducible fixed-grid floor       false
the ROHF/ROKS Z-vector Hessian is the dominant residual source  false
Python performs the production state-pair algebra               false
```

The durable lessons of the v3--v19 campaign are instead the exact ordered HST
decomposition, strict process-isolated numerical workers, one global state
gauge, the Fortran/Python TagArray transpose distinction, the `Sk` record
diagonal bug, the necessity of full JK+XC response, and the need to
differentiate both the owner motion and normalized fuzzy-cell weights of the
finite XC quadrature.
