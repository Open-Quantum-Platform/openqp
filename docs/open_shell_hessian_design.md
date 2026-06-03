# Open-shell (UHF / ROHF) analytic HF/DFT Hessian — design & status

This note tracks extending the validated closed-shell (RHF/RKS) analytic
Hessian (`source/modules/hf_hessian.F90`) to open-shell UHF and ROHF
references. The closed-shell kernel reads only the alpha density/MOs and treats
`nocc` as doubly occupied, so it must never run on an open-shell SCF; a
`scftype >= 2` guard in `hf_hessian` aborts rather than returning a wrong
matrix, and the Python input checker reports open-shell HF Hessians as
`unsupported_scf_type` until the response assembly lands.

The Hessian decomposes as

    H = E_nn'' + H_skeleton(1e + 2e, fixed density) + H_response(CPHF) .

`E_nn''` (`grd1::hess_nn`) is density-independent. The two open-shell-specific
subsystems are the skeleton and the CPHF orbital response.

## Stage 1 — open-shell skeleton (DONE, validated)

The fixed-density skeleton is spin-agnostic once it is given the correct
**total** density and the correct **energy-weighted (Lagrangian)** density:

  * total AO density `P = P_alpha + P_beta`  (`OQP_DM_A + OQP_DM_B`);
  * `W` from `grd1::eijden`, whose `scftype >= 2` branch already builds
    `W = -(Fa.Pa + Fb.Pb)` for UHF/ROHF;
  * the two-electron contraction via `grd2_uhf_compute_data_t` (Coulomb from
    the total density, exchange from each spin density) — the same object the
    open-shell production gradient `hf_2e_grad` uses.

Validated by `hess_skel_open_selftest` (bind C, open-shell analog of
`hess_skel_selftest`): analytic skeleton vs central FD of the open-shell
frozen-density gradient. OH radical / 6-31g*: UHF and ROHF both give
`max|analytic - FD| ~ 6e-9`. See `tests/test_hess_skel_open_selftest.py`.

## Stage 2 — open-shell (UHF) CPHF solver (DONE, HF validated)

`cphf_mod::cphf_solve_uhf` solves `M U = B` with the unknown/RHS vectors laid
out as the alpha occ-vir block (`la = nocca*nvira`) followed by the beta
occ-vir block (`lb = noccb*nvirb`), each in `iatogen`/`mntoia` occ-major order.
The UHF orbital-Hessian action is

    (M U)^s_ia = (e^s_a - e^s_i) U^s_ia + [ C^s^T dF^s C^s ]_ia ,
    dF^s = J[dP^a + dP^b] - c_x K[dP^s]   (+ spin-resolved f_xc for UKS),
    dP^s_mn = sum_ia ( C^s_mi U^s_ia C^s_na + C^s_ma U^s_ia C^s_ni ).

The 2e response reuses `scf_addons::fock_jk` (scftype>=2 -> `int2_urohf`:
Coulomb from the spin-summed trial density, exchange same-spin) and
`mod_dft_gridint_fxc::utddft_fxc` for the UKS kernel. The PCG solver is driven
per RHS with a diagonal `1/(e_a-e_i)` preconditioner.

Validated by `cphf_uhf_polarizability_selftest`: for a closed-shell molecule
run as UHF (mult 1, `Pa = Pb`) the UHF static dipole polarizability must equal
the validated closed-shell `cphf_static_polarizability`. H2O/6-31g: UHF
isotropic 4.02713729 vs RHF 4.02713607 (full tensor agrees to ~1e-6),
confirming the spin coupling and the `alpha_pq = -2 sum_s sum_ia mu U`
normalization. See `tests/test_cphf_uhf_polarizability.py`. The UKS `f_xc`
coupling compiles via `utddft_fxc` but is **not yet** finite-difference
validated.

## Stage 3 — UHF CPHF response assembly (DONE, HF validated)

`hf_hessian_mod::hf_hessian_uhf` (dispatched from `hf_hessian` for
`scftype == 2`, HF only) mirrors the validated closed-shell response block per
spin and sums over `s in {alpha, beta}` with **single** (not doubled)
occupation factors. Per spin `s`:

  * `S^x,s`, `h^x,s` MO transforms with `C^s`; occ-occ blocks `s1oo^s`.
  * CPHF RHS `B^s_ia = -(h^x_ia + G^{s,x}[P]_ia) + e^s_i S^x,s_ia - G^s[d0]_ia`,
    where `G^{s,x}[P]_ia` is the 2e response-Fock skeleton (probe
    `C^s_a C^s_i^T`), `d0^s = -sum_ij S^x,s_ij C^s_i C^s_j^T`
    (reorthonormalization, factor 1, not the closed-shell 2), and
    `G^s[d0] = J[d0^a + d0^b] - c_x K[d0^s]` from `fock_jk` (2-column).
  * Solve all `3N` RHS with `cphf_solve_uhf` -> `U^s`.
  * Relaxed `dC^s`, `dP^s` (occupation 1), `moe1^s`; assemble
    ```
    H_resp_xy = sum_s [ Tr[dP^s,y h^x] + Tr[dP^s,y G^{s,x}[P]] ]
              - 2 sum_s sum_i e^s_i (dC^s,y_i . S^x . C^s_i)
              -   sum_s sum_kl s1oo^s,x_kl moe1^s,y_kl
              -   sum_s Tr[Mi^s,x G^{s,y}[P]]
    ```
    symmetrized.  The closed-shell global factors {4,4,2} on
    {Tr[dP h^x]+A2, eps-overlap, occ-occ} become the per-spin {1,2,1} sums here
    (the closed-shell 2x came from double occupancy).

**Key building block (DONE):** the *two-density* derivative-Fock contraction
`fock_deriv_mod::fock_deriv_contract_os` (type `grd2_fockprobe_os_data_t`). It
takes a separate Coulomb density (total `Pa+Pb`) and exchange density (spin
`P^s`) and uses the full (not `1/2`) open-shell exchange factor, returning the
genuine `g_x = sum_uv M_uv ( J^x[pcoul] - c_x K^x[pexch] )`. Validated by
`fockx_os_selftest` against a frozen-density finite difference of
`Tr[M . fock_jk_spin]` (OH / 6-31g*, `max|an-FD| ~ 6.5e-9`, ratio 1.0); see
`tests/test_fockx_os_selftest.py`.

Validation of the full UHF Hessian: the OpenQP numerical Hessian (central FD of
the validated UHF analytic gradient).  Agreement at the FD-truncation level with
an exactly symmetric analytic matrix:
  * OH / STO-3G (dx=0.005): `max|analytic - numerical| ~ 6.7e-5`,
    `H_zz ~ 0.7545` (bond axis);
  * bent H2O+ / STO-3G (doublet): `~ 6.0e-5`;
  * bent H2O+ / cc-pVDZ (d functions, bfnrm != 1): `~ 1.3e-4`.
See `tests/test_uhf_hessian.py`.  Open-shell DFT (UKS) is still gated off: the
UKS `f_xc` response compiles but is not finite-difference validated.

## Stage 4 — ROHF CPHF response (PARTIAL: solver done, Hessian socc-block WIP)

ROHF uses a single MO set with a doubly-occupied / singly-occupied / virtual
partition, so the rotation space and coupling differ from UHF; it is **not** a
relabelled UHF path. The skeleton (Stage 1) already covers ROHF; only the
response differs.

### Stage 4a — ROHF CPHF solver (DONE, validated)

`cphf_mod::cphf_solve_rohf` solves `H theta = B` over the docc/socc/virt rotation
space (`rohf_pack_trial`/`rohf_unpack_trial`, the layout of
`scf_converger::pack_rohf_trial`: socc-docc, virt-docc, virt-socc blocks). The
orbital-Hessian action replicates the validated TRAH ROHF operator
(`scf_converger::calc_h_op`): per spin `Fvv x - x Foo` with `Foo`/`Fvv` the
occ-occ / vir-vir blocks of the converged spin Fock matrices (`OQP_FOCK_A/B`) in
the MO basis (full blocks, so non-canonical orbitals are handled), plus the
response Fock from `scf_addons::get_response_packed`. Diagonal
orbital-energy-gap preconditioner.

Validated by `cphf_rohf_polarizability_selftest`: a closed-shell molecule run as
ROHF (offset = n_socc = 0) reduces the rotation space to virt-docc and the ROHF
orbital Hessian to (twice) the RHF one; the ROHF isotropic polarizability is
4.02713548 vs the validated RHF 4.02713607 (~6e-7). See
`tests/test_cphf_rohf_polarizability.py`.

### Stage 4b — ROHF analytic Hessian (WIP, gated)

`hf_hessian_mod::hf_hessian_rohf` (dispatched for `scftype == 3`, HF only)
assembles the Hessian as `E_nn'' + skeleton + response`. The ROHF energy has the
same functional form as UHF in (Pa, Pb), so the response is evaluated
**semi-numerically**, reusing the validated analytic open-shell gradient: with
the CPHF-relaxed alpha/beta orbital derivatives `dCa`/`dCb` (built UHF-style from
the unpacked rotation `xa`/`xb`; the socc-docc rotation lives in `xb` since socc
is beta-virtual, relaxing Pb and leaving Pa invariant), the full electronic
Hessian column is the central finite difference, over BOTH geometry and the
relaxed orbital path, of the electronic gradient
(`grad_ee_overlap(W') + grad_ee_kinetic(P') + grad_en(P') + grad_2e(Pa',Pb')`,
with `W' = -(Pa' Fa' Pa' + Pb' Fb' Pb')`, `Fa'/Fb'` rebuilt as `Hcore' +
fock_jk` at the displaced geometry). Nuclear repulsion is added analytically
(`hess_nn`). The non-canonical CPHF RHS replaces the orbital energies by the full
Fock occ-occ blocks (`Foo^s`), reducing to the validated UHF RHS in the
canonical limit.

**Status / validation.** The closed-shell (offset=0) limit is EXACT: H2O run as
ROHF reproduces the RHF numerical Hessian to ~3e-5 (so the virt-docc RHS,
response, `dCa`/`dCb` build and the geometry+orbital FD are all correct). For a
genuine open shell the singly-occupied rotation blocks of the **CPHF
right-hand side** (socc-docc and virt-socc, via the non-canonical Foo-coupling /
reorthonormalization terms) still carry a residual error: OH/STO-3G `max|an -
num| ~ 3e-3` (localized to the degenerate-pi in-plane component) and bent
H2O+/STO-3G `~ 8e-3` (~1%, in-plane block; the bond-axis/z block is exact to
~1e-5). The bug is isolated to the socc-block RHS terms (the solver and the
offset=0 path are exact). Until those are corrected and finite-difference
validated, the Python input checker gates ROHF analytic Hessians to
`type=numerical`; the kernel remains reachable (dispatched) for development.

Remaining work: derive/validate the exact non-canonical CPHF RHS for the
socc-docc and virt-socc blocks (the Foo-coupling and reorthonormalization-density
terms), then re-enable the input-checker capability and add a `test_rohf_hessian`
finite-difference regression mirroring `test_uhf_hessian`.
