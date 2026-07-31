# Fully Analytic MRSF-TDDFT Nonadiabatic Couplings — Derivation

**Goal.** A closed-form interstate derivative coupling

```
d_IJ = <Psi_I | d/dx | Psi_J>,      h_IJ = (E_J - E_I) d_IJ
```

for MRSF-TDA states, with ONE z-vector solve per state pair (no finite
differences anywhere), reusing OpenQP's validated gradient machinery. The
structure is the interstate generalization of the Handy–Schaefer / Fatehi–
Subotnik Lagrangian, specialized to the MRSF three-space (doc/socc/virt)
ROHF reference and the spin-pair fold.

Conventions: `x` = one nuclear coordinate; `p,q,r,s` MO; `mu,nu` AO;
`i,j` occupied-alpha (doc+socc), `a,b` virtual-beta (socc+virt);
antisymmetric `d_IJ = -d_JI`, symmetric `h_IJ = +h_JI`, `gap = E_J - E_I`.

---

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

## 2. The orbital part d_orb — antisymmetry is emergent

### 2.1 The one-particle reduction

`d/dx phi_p = sum_q U^x_qp phi_q + (moving-AO term)`, so

```
T_pq := <phi_p | d/dx phi_q> = U^x_pq + Sk_pq ,
Sk_pq := sum_{mu,nu} C_mu_p C_nu_q <chi_mu | d/dx chi_nu>     (ket-half overlap derivative)
```

Orthonormality `d/dx <phi_p|phi_q> = 0` gives **T_qp = -T_pq exactly**:
the full orbital-derivative matrix is antisymmetric. (Equivalently
`U^x_pq + U^x_qp = -S^x_pq` with `S^x = Sk + Sk^T`.) This is where the
antisymmetry of d_IJ comes from — it must NOT be imposed by hand; if the
assembled d is not antisymmetric to solver tolerance, a term is missing.

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

Only the antisymmetric part `gamma_a := (gamma - gamma^T)/2` survives —
this is the first-principles derivation of the sign set that the branch
obtained by scanning `(sg_ij, sg_ab, sg_ia) = (1,-1,-1)` against an oracle,
and of the antisymmetrization `kernel_pair(I,J) - kernel_pair(J,I)`.

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

No fitted coefficient anywhere: the 1/2s and signs are fixed by
orthonormality and the perturbation identity. `d` antisymmetry and `h`
symmetry are consequences (gamma_a antisym, T antisym, gap flips sign under
I<->J) — compute (I,J) and (J,I) independently and CHECK, never stamp.

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
| Gates: L rotation-FD (unfolded), interchange identity, antisymmetry-emergent, C1 system, sum rule, frozen numeric refs with assert | TODO | new `tests/` files |

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
