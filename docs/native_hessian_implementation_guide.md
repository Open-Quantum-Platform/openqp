# Native HF/DFT Analytic Hessian — Implementation Guide

A runbook for completing the **native** (PySCF-free at runtime) analytic HF/DFT
Hessian in OpenQP, to be executed on a **libint-capable build machine**. It
consolidates the design decisions and validated building blocks already on this
branch. For the investigation narrative and rationale, see
`analytic_hessian_design.md`.

---

## 1. What is already done (validated, on this branch)

- **Working analytic HF/DFT Hessians via the PySCF bridge** (RHF/UHF + DFT),
  `pyoqp/oqp/library/external.py::analytic_hessian_from_pyscf`. This is the
  supported path today and the **element-wise validation oracle** for the native
  work. Regression test: `tests/test_pyscf_hessian_bridge.py`.
- **Native, validated Hessian terms:** `grd1::hess_nn` (nuclear repulsion),
  `grd1::hess_ee_overlap`, `grd1::hess_ee_kinetic` (1e overlap/kinetic).
- **Complete one-electron CPHF RHS derivative matrices:**
  `grd1::der_overlap_matrix` (`dS/dR`), `der_kinetic_matrix` (`dT/dR`),
  `der_nucattr_matrix` (`dV/dR`), built from the non-contracting primitives
  `comp_overlap_der1_block`, `comp_kinetic_der1_block`,
  `comp_coulomb_der1_block`, `comp_coulomb_helfeyder1_block`.
- **Validation harnesses:** `source/modules/hess1_selftest.F90` (a `bind(C)`
  harness that finite-differences the production gradient routines), the Fortran
  recursion test `tests/fortran/test_der2_recursion.F90`, and the NumPy
  reference `tests/test_nuclear_repulsion_hessian.py`.
- **Driver bookkeeping** (`hess_en`'s nine-block scatter via translational
  invariance) is validated in structure; only the second-derivative *integrals*
  feeding it were wrong (fixed-root `der2_coul` is invalid for Rys integrals --
  see `analytic_hessian_design.md`).

Anything not listed above (1e nuclear-attraction 2nd derivative, all 2e 2nd
derivatives, the CPHF solve, DFT XC 2nd derivatives) is **not** implemented; the
native `hf_hessian` kernel stays guarded.

---

## 2. Build prerequisites (one-time, maintainer task)

1. **Regenerate `oqp-libint` with second derivatives.** The custom tarball
   referenced in `external/CMakeLists.txt`
   (`oqp-libint-v2.7.1.1-am4`) must be generated with derivative order 2
   (`INCLUDE_ERI=2`, `INCLUDE_ONEBODY=2` in libint's config). Publish it and
   point `external/CMakeLists.txt` at the new URL.
2. **Build with `USE_LIBINT=ON`** (the CMake default). The Rys and libint 2e
   engines coexist and are selected at runtime via `libint2_active`
   (`int2.F90`), so enabling libint does not remove the Rys path. Energy/gradient
   behavior is unchanged.
3. Keep `USE_LIBINT=OFF` as the lighter build: in it the native analytic Hessian
   is unavailable and `[hess] type=analytical` must fall back explicitly to the
   PySCF bridge / numerical Hessian (no silent fallback).

The libint wrapper already has the deriv-2 ERI code path (`int_libint.F90`,
`case (2): libint2_build_eri2`, gated by `#if INCLUDE_ERI >= 2`).

---

## 3. Implementation tasks (dependency order)

### Tasks A & B — second-derivative integrals (1e nuclear-attraction, 2e ERI)

> **Correction (libint is the pragmatic path, not a mathematical necessity).**
> Native second-derivative nuclear-attraction and ERI integrals are *feasible*
> without libint, but require extending the Rys engine with **root/weight
> derivatives**. Code-grounded reason: the Rys argument is center-dependent --
> 1e: `X = pp%aa*sum((pp%r-c)**2)` (`mod_1e_primitives.F90:383`); 2e:
> `X = rho*sum((P-Q)**2)` with `P,Q` the Gaussian-product centers
> (`grd2_rys.F90:309`). The roots `u(X)`/weights `w(X)` therefore depend on the
> differentiated centers. First derivatives work with *frozen* roots because the
> polynomial angular-momentum shift alone suffices; the **second** derivative
> picks up a `(dX/dA)*d(roots)/dX` term that a fixed-root double-shift drops
> (this is exactly why the `der2_coul` fixed-root attempt failed -- it is a
> mathematical fact about the scheme, not an implementation bug). A correct
> native route implements `du/dX`, `dw/dX` and chains through (as GAMESS/HONDO
> do), or uses a parallel Obara-Saika/McMurchie-Davidson auxiliary-integral
> scheme. **Libint deriv-2 is chosen here only because it is lower-effort**, and
> because the deriv-2 ERI path (`libint2_build_eri2`) is already wired
> (`int_libint.F90`, gated by `#if INCLUDE_ERI >= 2`).

### Task A — 1e nuclear-attraction second derivative  (libint 1-body deriv-2)
- Compute libint 1-body deriv-2 integrals: `ipipnuc` (d2/dA2), `ipnucip`
  (d2/dA dB), and per-nucleus `ipiprinv`/`iprinvip` (charge center).
- Feed them into the existing `grd1::hess_en` block bookkeeping (already
  validated) in place of the removed `comp_coulomb_der2_braC` path. The PySCF
  `hcore_generator` assembly (basis-center `nuc` blocks + per-nucleus `rinv`
  blocks, combined by translational invariance) is the blueprint.
- **Validate:** `hess1_selftest` "nucattr" line should reach ~1e-10; element-wise
  vs PySCF `mol.intor('int1e_ipipnuc')`/`int1e_ipnucip` and
  `with_rinv_at_nucleus`+`int1e_ipiprinv`.

### Task B — 2e second derivative  (libint deriv-2 ERI)
- Use `libint2_build_eri2` (already wired) for `d2(uv|ls)`; contract
  `1/2 sum P_uv P_ls d2(uv|ls)` with Coulomb and exchange (apply `HFscale` for
  hybrids). Write a `hess_2e` driver paralleling `grd2_driver`/`grd2_rys`.
- **Validate:** finite difference of the existing 2e gradient (`grd2_driver`);
  element-wise cross-check vs PySCF `int2e_ipip1`/`int2e_ipvip1`.

### Task C — CPHF orbital response

> **Key enabler (native, no libint): the A-matrix builder already exists.**
> `scf_addons.F90::get_response_packed(basis, infos, molGrid, mo_a, dm1_tri,
> v1_tri, mo_b)` builds the AO linear-response vector `v1 = J/K(dm1) + fxc(dm1)`
> from a first-order density, for RHF/UHF/ROHF, using the native Rys 2e engine
> (`fock_jk`) plus the DFT XC kernel (`tddft_fxc`/`utddft_fxc`). This *is* the
> CPHF/CPKS A-matrix action. So the CPHF solver needs no libint: drive `pcg`
> with an `update` callback that (i) `iatogen`: MO trial vector U -> AO density,
> (ii) `get_response_packed`: AO density -> AO response Fock, (iii) `mntoia`:
> back to MO occ-vir, (iv) add the orbital-energy diagonal `(e_a - e_i) U_ia`.
> Preconditioner = `1/(e_a - e_i)`. This mirrors `tdhf_z_vector::compute_apbx`
> but calls the public `get_response_packed` instead of the private TD path.
>
> **Validate the solver independently of geometry derivatives** via the static
> **dipole polarizability**: build the dipole-integral RHS (AO dipole matrices
> exist in `electric_moments`), solve `A U = -mu_ai`, contract to get
> `alpha = -2 sum_ia mu_ai U_ia`, and compare to PySCF
> `mf.Polarizability().polarizability()` (or finite difference of the SCF dipole
> under a field). This exercises the exact A-matrix + PCG the Hessian uses, with
> no libint and no geometry-derivative integrals -- so the solver can be written
> and fully validated in any environment. Only the *nuclear* RHS (the 2e
> `dJ/dx[P]` piece, Task B's libint deriv-1) is environment-gated.

- **CPHF solver `cphf_solve(infos, nrhs, bvec, uvec)`** — DONE
  (`source/modules/cphf.F90`), validated vs PySCF dipole polarizability.
- **RHS `B^x` (per nuclear coordinate x):** the 1e parts are done
  (`der_overlap_matrix`, `der_kinetic_matrix`, `der_nucattr_matrix`). The
  remaining piece is `F^x` (the 2e derivative-Fock, below).

#### Task C1 — `F^x` 2e derivative-Fock builder (native, no libint)

> **Correction (driver rewrite, not a wrapper).** The code-grounded reality:
> `grd2_rys.F90` *fuses* the 4-index density contraction into its innermost
> assembly (`compute_der_ijkl`, ~line 773), emitting only the scalar gradient
> `gdat%fd(3,4)`; the bare `d(uv|ls)/dx` block is never materialized. So `F^x`
> reuses ~70% of `grd2_rys` UNCHANGED (`compute_rys_rw`, `compute_coefficients`,
> `compute_xyz_p0q0`, `compute_xyz_ijkl`, `compute_der_xyz_ijkl` -- the derivative
> shift `fi = g[m+1]*2*alpha - g[m-1]*(i-1)`), but the final assembly must be a
> NEW routine that *scatters* into an `(nbf,nbf,3,natom)` Fock derivative
> (`F^x_uv += d(uv|ls)/dx * D_ls` Coulomb, minus exchange * `HFscale`) instead of
> contracting all four indices to a scalar. Plus a 2-index density supplier in
> place of the 4-index `get_density`. ~150-250 new LOC, no new integrals, no
> libint. This is what makes the full native CPHF chain runnable for the nuclear
> case.

- Transform to the MO occupied-virtual RHS:
  `B^x_{ia} = (C^T (dh/dx + dG^x[P]) C)_{ia} - e_i (C^T dS/dx C)_{ia}`,
  where `dh/dx = dT/dx + dV/dx` (done) and `dG^x[P]` is `F^x`.
- **Solve** `A U^x = B^x` with `cphf_solve` (done).
- Assemble the response contribution to the Hessian from `U^x`.

#### Validation protocol for `F^x` and the full CPHF chain

The chain `F^x -> B^x -> U^x -> dC/dx -> dP/dx` is validated bottom-up, each
stage against a reference of known quality. The two-stage design avoids the
"reference-limited comparison" trap: the integral stages are checked against
*exact/non-iterative* references (truncation-limited, ~1e-9), and only the SCF
response stage uses converged-SCF references (whose quality is then *proven* by
h^2 scaling rather than asserted).

1. **`F^x` AO derivative-Fock (exact, non-iterative reference).**
   - *Trace identity:* for any fixed symmetric matrix `M`,
     `sum_uv M_uv F^x_uv` must equal the existing, already-validated 2e gradient
     contraction `grd2_driver` for the density product `M (x) P`. Equivalently,
     contracting `F^x[P]` with `P` reproduces the 2e part of the energy gradient
     `dE_2e/dx` (a single scalar per coordinate) to ~1e-9.
   - *Element-wise FD:* `F^x_uv = d/dx ( G_uv[P] )` at FROZEN density `P`, where
     `G[P] = fock_jk(P)`. Central difference:
     `F^x_uv ~ (G_uv[P; R+he_x] - G_uv[P; R-he_x]) / (2h)`, `h ~ 1e-4` bohr,
     `P` held fixed at the R-geometry converged density. This is the cleanest
     check: no SCF iteration in the reference, so accuracy is truncation-limited
     (expect ~1e-8..1e-9, scaling as h^2). Implemented in the `hess1_selftest`
     style (the existing harness already FDs `fock_jk`-type builds).

2. **`B^x` MO right-hand side.** Assemble
   `B^x_{ia} = (C^T(dh/dx + F^x[P])C)_{ia} - e_i (C^T dS/dx C)_{ia}`; check that
   each constituent (`dh/dx`, `dS/dx` already ~1e-14; `F^x` from step 1) is
   validated before combination. No separate FD needed beyond steps 1 and the
   existing 1e checks.

3. **`U^x` / `dC/dx` (converged-SCF reference, quality proven by scaling).**
   - The CPHF solution `U^x` gives the orbital response; the occupied MO
     coefficient derivative is `dC_{mu i}/dx = sum_a U^x_{ai} C_{mu a}`
     (occ-vir rotation; the occ-occ block is fixed by `dS/dx`).
   - *Reference:* converge SCF at displaced geometries `R +- h e_x`, extract the
     converged MO coefficients `C(R +- h)`, and form the central-difference
     `dC/dx`. **Phase/degeneracy fix:** align each displaced `C` to the reference
     by `sign`/Procrustes on the occupied block (MOs are defined up to sign and
     occupied-subspace rotation), then compare the gauge-invariant
     `dP/dx = sum_i (dC_i C_i^T + C_i dC_i^T)` density derivative rather than raw
     `dC` to sidestep gauge ambiguity.
   - *Primary comparison:* `dP/dx` from CPHF (`dP/dx = C U^x C^T + h.c.` on the
     occ-vir block) vs central-difference of the converged AO density
     `P(R +- h)`. This is gauge-free and the natural physical quantity.
   - *Expected accuracy & the anti-self-deception step:* do NOT report a single
     number. Run the FD at `h = 4e-3, 2e-3, 1e-3, 5e-4` and confirm the
     CPHF-vs-FD difference **decreases ~4x per halving (O(h^2))** until it hits
     the SCF-convergence floor. Tighten SCF (`scf.conv` to 1e-10/1e-11) so the
     floor is below the truncation curve over at least two h values. A genuine
     match shows the O(h^2) trend; a *reference-limited* comparison shows a flat
     plateau -- reporting the trend is the proof of correctness.

4. **Native-only.** All of steps 1-3 use only OpenQP (native Rys 2e, the SCF
   driver for displaced geometries, `cphf_solve`). PySCF is an OPTIONAL
   cross-check (`mf.Hessian()` includes the same CPHF response), never required.

Final (after second-derivative integrals land): full HF Hessian (direct + CPHF)
vs central FD of the analytic HF gradient on H2/H2O/formaldehyde, and optional
PySCF `mf.Hessian()` cross-check.

### Task D — DFT exchange-correlation second derivatives
- Add the XC-kernel grid second-derivative terms for RKS/UKS. Gate functional
  classes in `pyoqp/oqp/utils/input_checker.py` (reject meta-GGA / unsupported
  grid derivatives / unavailable ECP 2nd derivatives rather than proceeding).

### Task E — Integrate and unguard
- Assemble all terms in the native `hf_hessian` kernel; route
  `[hess] type=analytical` to native when available (add a `[hess] backend`
  selector), keeping the PySCF bridge as an automatic cross-check.
- Run the acceptance set (H2 / H2O / formaldehyde): report max-abs, RMS,
  symmetry error, displacement size, per `analytic_hessian_design.md`.
- Remove the guard in `source/modules/hf_hessian.F90`.

---

## 4. Validation playbook (use throughout)

- **PySCF element-wise oracle** (test-time only, not a runtime dependency):
  build the same molecule in PySCF and compare the analytic integrals/blocks
  directly -- `int1e_ipipnuc`, `int1e_ipnucip`, `int1e_ipiprinv`, `int2e_ipip1`,
  and the full `mf.Hessian().kernel()`. The bridge already wires PySCF.
- **Finite-difference harness:** follow `hess1_selftest` -- FD the production
  gradient routines and compare to the analytic second derivative, contracting
  with a fixed matrix so basis normalization is handled consistently. The
  separated-FD trick (FD basis vs charge positions independently, since
  `grad_en_*` take `coord` separately from `basis%atoms%xyz`) localizes failures
  by sub-block.
- **Acceptance:** central finite differences of analytic gradients on
  H2/H2O/formaldehyde; small basis first, then BHHLYP/6-31G* style cases.

---

## 5. Key file map

| Concern | File / routine |
|---|---|
| Hessian dispatch (Python) | `pyoqp/oqp/library/single_point.py::Hessian.analytical_hess` |
| PySCF bridge + oracle | `pyoqp/oqp/library/external.py::analytic_hessian_from_pyscf` |
| Native Hessian terms | `source/integrals/grd1.F90` (`hess_nn`, `hess_ee_overlap/kinetic`, `hess_en`, `der_*_matrix`) |
| 1e primitives / blocks | `source/integrals/mod_1e_primitives.F90` (`comp_*_der1`, `comp_*_der1_block`, `der_*_xyz`) |
| libint deriv ERI wrapper | `source/integrals/int_libint.F90` (`libint2_build_eri2`, `INCLUDE_ERI`) |
| 2e driver / engines | `source/integrals/int2.F90` (`int2_compute_t`, `libint2_active`), `int2e_rys`, `int2e_libint` |
| CPHF solve template | `source/modules/tdhf_z_vector.F90` (`compute_apbx`, `precond`, `tdhf_cg_data`), `source/pcg.F90` |
| Response Fock build | `source/tdhf_lib.F90` (`int2_td_data_t`, `iatogen`, `mntoia`) |
| MO data (tagarray) | `OQP::VEC_MO_A`, `OQP::E_MO_A`, `OQP::DM_A`, `OQP::FOCK_A`, `OQP::SM` |
| Self-test harness | `source/modules/hess1_selftest.F90` (`hess1_selftest` bind(C)) |
| Guarded native kernel | `source/modules/hf_hessian.F90` |
