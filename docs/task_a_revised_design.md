# Task A (revised) — 1e Nuclear-Attraction Second Derivative: Native Rys vs new Libint 1-body bindings

Status: **design only, no code.** Supersedes the Task A section of
`second_derivative_integral_design.md`, whose `ipipnuc`/`ipnucip` baseline is
**invalid** (proven below).

## 0. Why the previous baseline is invalid (settled)

`grd1.F90::hess_en` loops over nuclei `do ic = 1, nat` and, **inside that loop**,
scatters charge-center blocks 5-9 to Hessian atom-block coordinates that contain
`ic` as a row/column index:

```
block 5/6 -> (A, ic)/(ic, A) = p_AC(C=ic)
block 7/8 -> (B, ic)/(ic, B) = p_BC(C=ic)
block 9   -> (ic, ic)        = p_AA+bAB+bAB^T+p_BB, all at C=ic
```

Each is a physically distinct Hessian element `d2E/dA_a dC_b` per nucleus C.
Therefore `hess_en` requires the **per-nucleus single-center** second derivatives
of `-Z_C/|r-C|`: `p_AA(C)`, `p_AC(C)`, `p_BC(C)` (and `p_BB(C)` from the swapped
pair). The summed operators `int1e_ipipnuc`/`int1e_ipnucip` (= sum_C) cannot
populate blocks 5-9: the sum collapses the index `ic` the scatter uses as an atom
coordinate.

The required quantities are the OpenQP equivalents of PySCF `int1e_ipiprinv` /
`int1e_iprinvip` evaluated nucleus-by-nucleus (`with_rinv_at_nucleus(C)`).

Two findings constrain how to get them:
- **OpenQP libint is ERI-only.** `libint_f.F90` exposes only
  `libint2_build_eri{,1,2}`/`2eri`/`3eri` (four-center 2e). There are **no**
  1-body libint bindings (no `elecpot`, no nucleus-centered Coulomb). Grep for
  `rinv|nuclear|elecpot|onebody|overlap|kinetic` in `int_libint.F90` and
  `libint_f.F90` returns nothing.
- **OpenQP already has a native single-center Coulomb at an arbitrary center.**
  `int1.F90` calls `comp_coulomb_int1_prim(cntp, ig, c(:), -znuc, blk)` per
  nucleus; `comp_coulomb_*` / `QGaussRys(ryscomp, cp, id, c, znuc, ...)` build
  `1/|r-c|` via Rys for any `c`. The WIP `comp_coulomb_der2_braC` already takes
  `coord(:,ic), -zq(ic)` per nucleus. Its first derivative is correct; only its
  *second*-derivative scheme (fixed-root) is wrong.

So Task A is a choice between (1) fixing the native Rys single-center second
derivative, or (2) building 1-body libint bindings from scratch.

---

## 1. The mathematical core both options must get right

`X = pp%aa * sum((pp%r - c)**2)` with `pp%r = (a_i r_i + a_j r_j)/(a_i+a_j)` the
bra-ket product center. Rys roots `u_t(X)`, weights `w_t(X)` depend on `X`, which
depends on the basis centers and on `c`. The exact second derivative is

```
d2 I / dA dB = sum_t [ (d2 poly_t/dA dB) w_t
                     + (d poly_t/dA)(dw_t/dX)(dX/dB) + (A<->B)
                     + poly_t ( (d2w_t/dX2)(dX/dA)(dX/dB) + (dw_t/dX)(d2X/dA dB) )
                     + (root-shift cross terms u_t and du_t/dX analogous) ]
```

The fixed-root recursion-twice keeps only the first term (`d2 poly` at frozen
roots) and drops every `du_t/dX`, `dw_t/dX` term -- which is nonzero because X
depends on the differentiated centers. That is the confirmed root cause of the
on-axis failure. `dX/dA`, `d2X/dA dB` are elementary (X is quadratic in centers);
the missing pieces are `du_t/dX`, `dw_t/dX`, `d2u_t/dX2`, `d2w_t/dX2`.

---

## 2. Option 1 -- Native Rys root-derivative (extend comp_coulomb_* / QGaussRys)

### 2a. Differentiability of the roots (decisive cost driver)

`rys.F90` computes roots/weights two ways (verified):
- **nroots 1-5: closed-form** piecewise polynomial + `exp(-x)` branches
  (`rys_rt1..rys_rt5`). These are **analytically differentiable in X by hand**;
  `du/dX`, `dw/dX`, second derivatives are obtained by differentiating the same
  Chebyshev/poly branches. Bounded, self-contained work.
- **nroots >= 6: `rys_general`** -- discretized Stieltjes + Golub-Welsch
  **eigenvalue** quadrature (plus an asymptotic Hermite branch). No closed-form
  X-derivative; differentiating it means differentiating an eigenvalue solve
  (Hellmann-Feynman-style dlambda/dX) or finite-differencing the root finder.

Practical reach: `BAS_MXANG = 6`, `mxrys = 13`, `MAX_NROOTS = 7`. nroots for the
1e nuclear-attraction *second* derivative ~ `(l_i + l_j + 2)/2 + 1`. For shells up
to f (l <= 3) this is nroots <= 5 -> entirely closed-form. Only g-and-higher on
both centers (l >= 4) reach `rys_general`. So Option 1 is closed-form for the
common case and only needs the hard eigenvalue-derivative path for high L.

### 2b. Scope of work
- Add `du/dX`, `dw/dX` (and 2nd derivatives) for `rys_rt1..rys_rt5` by
  differentiating the existing branch polynomials. (New code, but mechanical and
  unit-testable against finite differences of the existing root routines.)
- Decide high-L (nroots >= 6) policy: either implement the Golub-Welsch root
  derivative, or restrict the native Hessian to l <= 3 and reject l >= 4 in the
  input checker (documented limitation), or finite-difference the roots there.
- Rewrite `comp_coulomb_der2_braC` (and helpers) to assemble the full second
  derivative using `du/dX,dw/dX` + `dX/dA,d2X/dAdB` + the existing poly shifts.
  Keep the per-nucleus `(c, znuc)` interface `hess_en` already calls -- the
  driver and the 9-block scatter are unchanged.

### 2c. Pros / cons
- + No new external dependency; buildable in any environment (incl. this sandbox).
- + Reuses QGaussRys, the per-nucleus interface, and the validated hess_en scatter.
- + Closed-form for l <= 3 (covers most production bases).
- - The root-derivative algebra is the exact thing that broke once; demands the
    strongest validation (Section 4).
- - High-L (l >= 4) needs the eigenvalue-derivative path or a documented cap.

---

## 3. Option 2 -- New Libint 1-body derivative bindings (from scratch)

OpenQP has **zero** 1-body libint bindings. To use `int1e_ipiprinv`-equivalent
integrals would require building the entire 1-body libint layer:
- Fortran `bind(C)` interfaces for the libint 1-body engines (the
  `libint2_build_<onebody>` family) -- none exist; `libint_f.F90` has only the
  eri/2eri/3eri families.
- Operator-center control: libint's nuclear/`rinv` operator is configured through
  the `libint2_params`/engine `set` API in C++; OpenQP's wrapper has no path to
  set an operator center per nucleus. This is the `with_rinv_at_nucleus`
  equivalent and must be plumbed through new C-ABI calls.
- `init`/`cleanup`/buffer handling for the 1-body engine (separate from the eri
  engine already wired).
- Regenerate `oqp-libint` with `INCLUDE_ONEBODY >= 2` AND confirm the chosen
  tarball actually compiles the 1-body deriv-2 targets.
- A new compute/scatter path feeding per-nucleus `ipiprinv`/`iprinvip` into
  hess_en.

### 3c. Pros / cons
- + Integrals themselves are reference-quality (libint handles root derivatives
    internally); no Rys-root algebra to get wrong.
- - Largest surface of *new* code: an entire 1-body libint binding layer OpenQP
    has never had, including C++-side operator-center control.
- - Not buildable in the current sandbox (libint host unreachable); blocks all
    validation here.
- - Component-ordering / normalization / AO-permutation matching vs OpenQP bfnrm
    is a fresh, untested integration surface (the §4 permutation risks apply in
    full).
- - Per-nucleus operator-center control via the libint C ABI is non-trivial and
    may require upstream `oqp-libint` changes, not just Fortran.

---

## 4. Validation tests (apply to whichever option; designed to catch the prior failure classes)

Mandatory, in order; a term is not accepted until all pass.

1. **Root-derivative unit test (Option 1 only).** `du/dX`, `dw/dX` for
   `rys_rt1..5` vs central finite difference of the existing `rys_rtN` routines in
   X. Pass: O(h^2), abs < 1e-9 across >= 3 X values spanning all branch
   boundaries (x <= 3e-7, <=1, <=3, <=5, >5). This is the test the fixed-root
   scheme never had.
2. **Per-nucleus integral oracle (the fixed-root guard).** For each nucleus C,
   compare `p_AA(C)`, `p_AB(C)` element-wise to PySCF
   `with_rinv_at_nucleus(C); int1e_ipiprinv / int1e_iprinvip` (times -Z_C),
   matching `cart=True` ordering and OpenQP bfnrm via an **explicit AO-permutation
   map** (no Frobenius/trace shortcuts). Pass: max-abs < 1e-9. Run on a
   **low-symmetry, all-distinct-atom** geometry (e.g. distorted HOF or
   off-axis-rotated water) so no block is zero by symmetry.
3. **Component-transpose check.** Explicitly compare `H_AB[a,b]` vs `H_AB[b,a]`
   for `a != b` against the `iprinvip` oracle, to catch a (a_bra,b_ket)
   permutation that the (a,b)-symmetric `ipipnuc` would hide.
4. **Mixed-block independent cross-check (breaks the TI circularity).** Compute
   `H_AC(C)` two ways: (i) via translational invariance from `(p_AA,p_AB)`, and
   (ii) directly from the per-nucleus oracle (PySCF hcore_generator block for that
   nucleus). Require agreement < 1e-9. Without this, TI and symmetry are satisfied
   by construction and cannot falsify the assembly.
5. **Contracted FD vs production gradients (mandatory, not optional).** hess1e
   "nucattr" line: analytic `hess_en` vs central FD of
   `grad_en_pulay + grad_en_hellman_feynman`, using the separated-FD harness
   (basis vs charge positions independently). Pass: O(h^2) ratio ~4.0 over
   h in {8e-3,4e-3,2e-3,1e-3} bohr, into the SCF/FD floor; abs < 1e-6. Catches
   missing/double `-Z_C` and density-convention errors.
6. **Charge-factor check.** Verify the `-Z_C` factor by comparing a one-nucleus
   contribution scaled vs PySCF's charged `hcore_generator` block (not the bare
   intor).
7. **Symmetry / TI as CONSISTENCY checks only** (demoted): pre-symmetrization
   asymmetry < 1e-12 (asserted, not just recorded); per-atom 3x3 column sum over
   all atoms < 1e-10. These cannot validate a TI-built assembly on their own
   (they pass by construction); they are guards against gross bugs.
8. **Subsystem regression guard (boundary enforcement).** Re-run
   `fockx_selftest`, `cphf_f0x_selftest`, `cphf_dpdx_selftest` before/after; any
   change in their results fails Task A (proves the frozen response subsystem was
   not perturbed).
9. **Low-symmetry full-block coverage.** Require the validation molecule to have
   all `(3N,3N)` atom blocks nonzero so blocks 5-9 are independently constrained;
   H2/H2O under-determine the charge-center cross blocks.

---

## 5. Risk comparison

| Risk | Option 1 (native Rys) | Option 2 (libint 1-body) |
|---|---|---|
| Root-derivative algebra wrong (repeat of fixed-root) | HIGH -- but directly testable (test 1) here | none (libint internal) |
| New external/C-ABI surface | none | HIGH -- whole 1-body layer + operator-center control |
| Buildable/validatable in current sandbox | YES | NO (libint host unreachable) |
| Component/AO-permutation/normalization mismatch | low (same engine as production 1e) | HIGH (fresh integration) |
| High-L (l>=4) coverage | needs eigenvalue-deriv or documented cap | native to libint |
| Touches frozen subsystem | no | no |
| Effort | medium (closed-form l<=3) + high (l>=4) | high (binding layer) + upstream libint |

---

## 6. Recommendation

**Option 1 (native Rys root-derivative), scoped to l <= 3 first.** Rationale:
- It is the only option **buildable and validatable in the current environment**;
  Option 2 cannot even be compiled here.
- For l <= 3 the roots are **closed-form and hand-differentiable**, so the hard
  part (test 1) is bounded and unit-testable -- precisely the guard the fixed-root
  attempt lacked.
- It reuses the per-nucleus interface and the (provisional) hess_en scatter
  unchanged; no new dependency, no C-ABI work.
- High-L (l >= 4) is handled as a *separate, later* decision (eigenvalue-root
  derivative vs documented cap), not a blocker for the first validated term.

Option 2 remains the fallback if the native root derivative proves intractable at
the required accuracy, but it should not be the first attempt: it maximizes new,
unvalidatable-here code and a fresh integration surface.

---

## 7. Exact files

### To touch (Option 1)
- `source/integrals/rys.F90` -- add root/weight X-derivatives for `rys_rt1..5`
  (new routines, e.g. `rys_rtN_deriv`); do not alter the existing root routines.
- `source/integrals/mod_1e_primitives.F90` -- rewrite `comp_coulomb_der2_braC`
  (and remove/replace the fixed-root `der2_coul_xyz` helper) to assemble the full
  second derivative with the root-derivative terms; keep the per-nucleus
  `(c, znuc)` signature `hess_en` calls.
- `source/integrals/grd1.F90` -- `hess_en` body only if the block-builder
  signature changes; the `do ic` loop and 9-block scatter stay. Update the WIP
  `@warning` once validated; re-validate the scatter end-to-end (it is currently
  provisional, never run against correct integrals).
- `source/modules/hess1_selftest.F90` -- re-enable the "nucattr" PASS/FAIL line;
  add tests 2-4 (per-nucleus oracle, transpose, mixed cross-check) driven from a
  low-symmetry geometry.
- `docs/second_derivative_integral_design.md`, this note -- status updates.

### Not to touch (read-only; frozen subsystem)
- `source/modules/fock_deriv.F90`, `source/modules/cphf.F90`
- `grd1.F90` first-derivative *matrix* builders `der_overlap_matrix`,
  `der_kinetic_matrix`, `der_nucattr_matrix`, and `comp_*_der1_block`
- `source/modules/fock_deriv_selftest.F90`, `cphf_nuclear_selftest.F90`,
  `cphf_dpdx_selftest.F90`
- the validated overlap/kinetic 2nd-deriv path (`der2_kinovl_xyz`,
  `comp_overlap_der2`, `comp_kinetic_der2`, `hess_ee_overlap`, `hess_ee_kinetic`)
- `source/modules/hf_hessian.F90` -- stays guarded; Task A does not unguard it.

### Libint files (Option 2 only; not recommended first)
- `source/integrals/libint_f.F90`, `int_libint.F90` -- new 1-body bind(C)
  interfaces, init/cleanup, operator-center control, buffer handling.
- `external/CMakeLists.txt` and the `oqp-libint` tarball -- `INCLUDE_ONEBODY>=2`
  plus possible upstream operator-center support.

---

## 8. Open questions to resolve before coding (Option 1)
1. Exact nroots needed for the *second* derivative at each (l_i, l_j); confirm the
   l <= 3 -> nroots <= 5 closed-form boundary and the xyzin array sizing for the
   extra angular momentum.
2. High-L policy: implement Golub-Welsch root derivative, finite-difference the
   roots for nroots >= 6, or cap at l = 3 with an input-checker rejection.
3. Whether `comp_coulomb_der2_braC`'s current bra/ket/charge factorization is the
   right decomposition for the full (root-derivative-complete) second derivative,
   or whether it must be restructured (re-derive from the per-nucleus operator,
   do not patch).
4. Confirm the hess_en 9-block scatter against a correct per-nucleus integral on a
   low-symmetry molecule **before** declaring it validated.
