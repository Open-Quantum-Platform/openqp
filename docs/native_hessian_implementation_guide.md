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

- **RHS `B^x` (per nuclear coordinate x):** the 1e parts are done
  (`der_overlap_matrix`, `der_kinetic_matrix`, `der_nucattr_matrix`). Add the 2e
  response `dJ/dx[P]`, `dK/dx[P]` (libint deriv-1 Fock build, or the existing
  `int2_td_data_t` path driven with derivative integrals). Transform to the MO
  occupied-virtual block:
  `B^x_{ia} = (C^T (dh/dx + dG/dx[P]) C)_{ia} - e_i (C^T dS/dx C)_{ia}`.
  The MO transform path is validated (`C^T S C = I` check in `hess1_selftest`).
- **Solve `A U^x = B^x`:** reuse the Z-vector CG block from
  `tdhf_z_vector.F90` *verbatim* with `B^x` as the right-hand side -- the static
  CPHF A-matrix is the same orbital Hessian `(A+B)` that `compute_apbx` already
  applies (response Fock via `int2_td_data_t` + `int2_compute_t%run` + DFT
  `fxc`), with diagonal preconditioner `precond` (`1/(e_a-e_i)`). Expose/refactor
  `compute_apbx`, `precond`, `tdhf_cg_data` into a shared module so CPHF can call
  them. The A-matrix uses *undifferentiated* 2e integrals (Rys, no libint deriv).
- Assemble the response contribution to the Hessian from `U^x`.
- **Validate:** `U^x` vs finite differences of the converged MO coefficients;
  the full HF Hessian (direct + CPHF) vs FD of the HF gradient and vs the PySCF
  bridge (which already includes CPHF).

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
