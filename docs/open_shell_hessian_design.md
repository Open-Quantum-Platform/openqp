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

## Stage 4 — ROHF CPHF response (DONE, validated)

ROHF uses a single MO set with a doubly-occupied / singly-occupied / virtual
partition, so the rotation space and coupling differ from UHF; it is **not** a
relabelled UHF path. The skeleton (Stage 1) already covers ROHF; only the
response differs.

### Stage 4a — ROHF CPHF solver (DONE, validated)

`cphf_mod::cphf_solve_rohf` solves `H theta = B` over the docc/socc/virt rotation
space (`rohf_pack_trial`/`rohf_unpack_trial`, the layout of
`scf_converger::pack_rohf_trial`: socc-docc, virt-docc, virt-socc blocks). The
Fock-transform part of the orbital-Hessian action is the EXACT commutator
`[F^s_MO, K]_vo` per spin, where `K` is the antisymmetric MO rotation built from
the packed trial vector (vir-occ_alpha from `xa`, socc-docc occ-occ from `xb`)
and `F^s_MO` the converged spin Fock in the MO basis (full matrix, `OQP_FOCK_A/B`),
plus the response Fock from `scf_addons::get_response_packed`. The canonical form
`Fvv x - x Foo` (used by the TRAH operator `scf_converger::calc_h_op`) is NOT
sufficient: the raw spin-Fock vir-occ blocks are nonzero (~0.03 Ha) for the
non-canonical ROHF orbitals -- even though the *packed* orbital gradient is zero
at convergence -- and their commutator coupling to the socc rotations is the term
the canonical form drops. Validated in PySCF against a finite difference of the
packed orbital gradient: canonical form `max err ~7e-2`, commutator form `~3e-6`.
Diagonal orbital-energy-gap preconditioner.

Validated by `cphf_rohf_polarizability_selftest`: a closed-shell molecule run as
ROHF (offset = n_socc = 0) reduces the rotation space to virt-docc and the ROHF
orbital Hessian to (twice) the RHF one; the ROHF isotropic polarizability is
4.02713548 vs the validated RHF 4.02713607 (~6e-7). See
`tests/test_cphf_rohf_polarizability.py`.

### Stage 4b — ROHF analytic Hessian (DONE, validated)

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
(`hess_nn`).

The non-canonical CPHF right-hand side is the Pulay form
`B^s_ai = -(h^x + G2e + G[d0])_ai + sum_{j in occ} ( S^x_aj F^s_ji + F^s_aj S^x_ji )`,
i.e. the *occupied-projected anticommutator* of the overlap derivative and the
spin Fock. The first sum is the usual `eps_i S^x_ai` in the canonical limit; the
second (`F^s_aj` is a vir-occ Fock element, zero canonically) is the
non-canonical correction the socc rotations need. Both the operator commutator
and this RHS term are required: with the canonical operator + canonical RHS the
ROHF Hessian is off by ~3e-3; the commutator operator alone leaves ~6e-4; both
together reach the finite-difference floor.

**Status / validation.** Validated against the OpenQP numerical Hessian (central
FD of the analytic ROHF gradient), exactly symmetric, at the FD-truncation level:
  * OH / STO-3G (degenerate pi socc): `max|analytic - numerical| ~ 2.5e-5`;
  * bent H2O+ / STO-3G (doublet): `~ 1.3e-4`;
  * bent H2O+ / cc-pVDZ (d functions): `~ 1.4e-4`.
A per-perturbation density-derivative isolation (non-degenerate H2O+/STO-3G,
calibrated against re-SCF `dPa`/`dPb`) shows the CPHF-relaxed `dPa`/`dPb` match to
~2e-5 for every Cartesian perturbation. See `tests/test_rohf_hessian.py`. The
Python input checker enables ROHF (HF) analytic Hessians. Open-shell DFT (ROKS)
remains gated (the ROKS `f_xc` response is not finite-difference validated).

Note: PySCF (2.13) has NO analytic ROHF Hessian (`mf.Hessian()` raises
`NotImplementedError`) -- mature codes compute ROHF frequencies numerically -- so
this native analytic ROHF HF Hessian goes beyond the common reference
implementations. The key correctness tools were (a) a PySCF finite-difference
validation of the orbital-Hessian *operator* (which exposed the missing
commutator coupling) and (b) the per-perturbation re-SCF density-derivative
isolation (which localized and then confirmed the fixes).

## Stage 5 — open-shell DFT (UKS done/validated; ROKS pending)

The XC contribution is added on top of the open-shell HF Hessian, mirroring the
closed-shell RKS structure for two spins (open-shell `dftexcor`/`derexc_blk` with
`urohf=.true.`; the HF-exchange fraction stays in `hfscale`; the CPHF operator
already carries the open-shell `f_xc` kernel via `utddft_fxc`/`get_response_packed`):

  * RHS XC: one central FD of the spin XC Fock along the combined
    geometry + reorthonormalization path, giving the XC skeleton `dVxc/dR` and
    `f_xc[d0]`, subtracted from the CPHF right-hand side.
  * `dHse`: XC skeleton + density response, central FD of the analytic open-shell
    XC gradient (`derexc_blk`) along `R +/- h`, `P^s +/- h dP^s`.
  * `dHt3`: the XC part of the energy-weighted term,
    `- sum_s sum_kl s1oo^s,x_kl (vxc^s,y + fxc[dP^s,y])_kl`, with coefficient
    **-1 per spin** (the UHF occupation factor; the closed-shell RKS value is -2
    for double occupancy -- getting this factor right was the key fix).

### UKS (DONE, validated)

`hf_hessian_uhf` uses the analytic UHF response assembly, whose energy-weighting
already carries Vxc through the KS orbital energies (`eps^s`), so the three XC
pieces above compose cleanly. Validated against the OpenQP numerical Hessian,
exactly symmetric, at the FD floor:
  * offset=0 (closed-shell run as UKS) vs RKS analytic: `~7e-6`;
  * OH / 6-31G bhhlyp (hybrid): `~2.1e-5`; bent H2O+ / 6-31G bhhlyp: `~1.9e-5`;
  * OH / 6-31G PBE (pure GGA): `~2.8e-5`.
See `tests/test_uks_hessian.py`. The input checker enables UHF/UKS analytic
Hessians.

### ROKS (DONE, validated)

`hf_hessian_rohf` evaluates the HF response semi-numerically (`resp_grad`). The
XC is folded into the SAME finite difference, with NO separate `dHt3` term (a
separate analytic XC energy-weighting would double-count the Vxc already carried
by W'):
  * the spin Fock used for the energy-weighted density `W'` is the full KS Fock
    (Hcore + J - c_x K + Vxc), Vxc from the open-shell `dftexcor` at the displaced
    geometry/orbitals;
  * the explicit open-shell XC gradient (`derexc_blk`, `urohf=.true.`) at the
    displaced density is added to the gradient.
The CPHF right-hand side carries the XC skeleton `dVxc/dR` + `f_xc[d0]` (FD of the
spin XC Fock along the geometry+reorthonormalization path), and the CPHF operator
already includes the open-shell `f_xc` kernel (`get_response_packed`). A grid
warm-up flushes stale state left by the CPHF solver before the moving-grid FD.

Validated against the OpenQP numerical Hessian, exactly symmetric, at the FD
floor: OH/6-31G bhhlyp `~3.6e-5`; bent H2O+/6-31G bhhlyp `~1.2e-4`; OH/6-31G PBE
`~4.4e-4`. See `tests/test_roks_hessian.py`. The input checker enables ROHF/ROKS
analytic Hessians.

NB: the earlier "gross error" attempts combined Vxc-in-W' with a separate analytic
`dHt3` (double-counting the XC energy-weighting). The correct semi-numerical ROKS
decomposition is the full-KS `W'` (Vxc included) plus the explicit `derexc`
gradient, and no `dHt3` -- the resp_grad finite difference already contains the
W-derivative.

## Stage 6 — robustness sweep & feature gates

Broad analytic-vs-numerical validation across functionals, bases, ECPs and
references (OH, H2O+, HBr+; UHF/UKS/ROHF/ROKS), `max|analytic - numerical|`:

| class | examples | result |
|-------|----------|--------|
| LDA / GGA / hybrid / meta-GGA | slater, PBE, b3lyp5, bhhlyp, **m06-2x** | 3e-5 - 4e-4 ✅ |
| f-functions (cc-pVTZ) | UHF/ROHF/UKS, OH | 2e-5 - 6e-5 ✅ |
| polyatomic | H2O+ ROKS/PBE | 7e-5 ✅ |
| **range-separated (CAM/LC)** | cam-b3lyp | ~1e-1 ❌ |
| **ECP** | HBr+ / LANL2DZ | 1e2 - 1e18 ❌ |

The two failures are **pre-existing limitations of the OpenQP analytic Hessian**,
not the open-shell extension: closed-shell RHF/RKS fail identically (RKS
cam-b3lyp ~1.2e-1, RHF/HBr/LANL2DZ ~3e2). The native derivative-integral
machinery (`der_nucattr`/`hess_en`, the Rys 2e second-derivative kernel,
`fock_deriv_contract`) does not include ECP second derivatives, and the 2e
derivative integrals are not range-separation (CAM) screened.

These are now **gated to the numerical Hessian** rather than silently returning a
wrong matrix:
  * runtime guard in `hf_hessian` (all references): aborts with a clear message if
    `any(ecp_zn_num /= 0)` or (DFT) `infos%dft%cam_flag`;
  * Python input checker: reports range-separated functionals (`cam*`, `dtcam*`,
    `stg*`, `lc-*`, `wb97`) and known ECP basis sets (`lanl*`, `sbkjc`, `crenb`,
    `stuttgart`) as `unsupported_feature` -> use `[hess] type=numerical`.

Supported analytic-Hessian scope (validated): RHF/RKS, UHF/UKS, ROHF/ROKS, with
LDA/GGA/hybrid/meta-GGA functionals and standard (incl. f-function) basis sets.

## Stage 7 — ECP (effective core potential) second derivatives

ECP is no longer a limitation: the analytic Hessian now contracts the ECP
second derivatives directly from **libecpint** (which already implements them up
to derivative order 2 — only the C wrapper was missing).

Wiring:
  * `source/wrapper/libecpint_wrapper.cpp` — new `compute_second_derivs` C entry
    point (mirrors `compute_first_derivs`), concatenating the
    `3N(3N+1)/2` packed `ncart x ncart` derivative matrices;
  * `ecpint.F90` — Fortran interface for `compute_second_derivs`;
  * `ecp.F90`
    - `add_ecphess(basis, coord, denab, hess)` — the ECP **skeleton**
      `sum_uv P_uv d^2 V_ECP_uv/dR_I dR_J`, decoded from libecpint's atom-pair
      packing (`H_START(I,J,N)`; diagonal blocks store 6 components
      `{xx,xy,xz,yy,yz,zz}`, off-diagonal blocks 9 row-major `{xx..zz}`) and
      scattered symmetrically into the `(3N,3N)` Hessian;
    - `ecp_deriv_ints(basis, coord, dVecp)` — the ECP **first**-derivative
      integrals `dV_ECP/dR` (uncontracted), added into the core-Hamiltonian
      derivative so the ECP drives the CPHF RHS and the orbital relaxation.

Assembly per reference (`hf_hessian.F90`):
  * **RHF/UHF** (analytic CPHF): `add_ecphess` supplies the skeleton; `dVecp` is
    added into `dVa` (the nuclear-attraction derivative tensor) so the ECP enters
    the RHS, the energy-weighted term and `F^x` — the exact parallel of how
    point-charge nuclear attraction is split into `hess_en` (skeleton) + `dVa`
    (response);
  * **ROHF/ROKS** (semi-numerical `resp_grad`): the ECP gradient `add_ecpder` is
    folded into the displaced gradient, so the central FD over geometry+orbital
    path yields both the ECP skeleton and response. The ECP centre `ecp_coord`
    (sized `3*num_ecps`, indexed per ECP centre, **not** per atom) is moved in
    lockstep with its atom via an `iecp_atom` map; the ECP is **not** added to the
    Fock used for `W'` (that would double-count the ECP Pulay term already in
    `add_ecpder`).

Latent bug fixed along the way: the CPHF response built `dVa` from
`der_nucattr_matrix(..., basis%atoms%zn, ...)` — the **full** nuclear charge —
while the SCF, `hess_en` and `hess_nn` use the **ECP-screened** charge
`zn - ecp_zn_num`. Harmless when `ecp_zn_num = 0`, this made the response see a
+Z core the valence never feels and blew the ECP Hessian up by ~1e2-1e25. All
three kernels now pass the screened charge.

Validation (`max|analytic - numerical|`, LANL2DZ, exactly symmetric):

| system | reference | result |
|--------|-----------|--------|
| HBr | RHF, RKS/bhhlyp | 3.8e-5 / 1.5e-4 ✅ |
| HBr+ | UHF, UKS/bhhlyp | 1.2e-5 / 3.0e-5 ✅ |
| HBr+ | ROHF, ROKS/bhhlyp | 1.9e-5 / 1.2e-5 ✅ |
| NaCl (two ECP centres) | RHF, ROHF | 1.1e-6 / 8.7e-7 ✅ |
| BrH2+ (3 atoms) | RHF | 3.1e-5 ✅ |

Test: `tests/test_ecp_hessian.py`. The ECP gate is removed from the runtime guard
and the Python input checker; **range-separated (CAM/LC) remains gated** to the
numerical Hessian (Stage 8 candidate: split the 2e derivative-integral assembly
into long-range Coulomb + short-range erfc-exchange passes).

## Stage 8 — range-separated (CAM/LC) functionals

CAM is no longer a limitation either. CAM-type functionals split `1/r12` into a
long-range part (full Coulomb `J` + exchange `K[1/r]` scaled by `cam_alpha`) and
a short-range part (erfc-attenuated exchange `K[erf(mu r)/r]` scaled by
`cam_beta`). OpenQP's 2e derivative-integral engine is already erfc-attenuation
capable (`grd2_rys_compute`/`grd2_rys_hess_compute` take `mu2`), so every Fock
build in the analytic Hessian only had to run the established two-pass split:

| Hessian term | routine | CAM handling |
|--------------|---------|--------------|
| 2e second-derivative skeleton | `grd2_hess_driver` | **new** two-pass wrapper (this stage), mirroring `grd2_driver` |
| CPHF response-Fock derivative `G[P]^x` / `F^x` | `fock_deriv_contract(_os)` -> `grd2_driver` | already auto-detects `cam_flag` |
| reorthonormalization `G[d0]`, relaxed-density `G[dP]` | `fock_jk` -> `int2_driver` | already CAM-aware |
| CPHF A-matrix (orbital Hessian) | `cphf_apbx`/`_uhf`/`_rohf` | already CAM-aware |
| XC skeleton/kernel (`dHse`, RHS `dVxc`, `f_xc`) | `dftexcor` / `derexc_blk` | functional's own XC; unaffected by the range split |

The only code change was splitting `grd2_hess_driver` into a thin CAM two-pass
wrapper plus `grd2_hess_driver_gen` (the former body), exactly parallel to
`grd2_driver`/`grd2_driver_gen`. Everything else was already in place — the
robustness sweep's ~1e-1 CAM error was entirely the un-split 2e skeleton.

Validation (`max|analytic - numerical|`, 6-31g OH/OH+, exactly symmetric):

| functional | RHF/RKS | UHF/UKS | ROHF/ROKS |
|------------|---------|---------|-----------|
| cam-b3lyp | 1.5e-4 | 4.2e-5 | 1.9e-5 |
| wb97x | 2.7e-5 | 2.7e-5 | — |
| lc-blyp | 2.7e-5 | 4.1e-5 | — |

CAM + ECP compose: CAM-B3LYP / LANL2DZ / HBr (RHF) = 2.2e-4. Test:
`tests/test_cam_hessian.py`. The CAM gate is removed from the runtime guard and
the input checker.

**Final analytic-Hessian scope (all finite-difference validated):** RHF/RKS,
UHF/UKS, ROHF/ROKS, with LDA/GGA/hybrid/meta-GGA **and range-separated (CAM/LC)**
functionals, standard (incl. f-function) **and ECP** basis sets. No remaining
feature gates for the ground-state HF/DFT analytic Hessian.
