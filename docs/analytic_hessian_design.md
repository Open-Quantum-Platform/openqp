# OpenQP Analytic Hessian Design Note

This note defines the staged analytic-Hessian contract for HF/DFT, conventional TDDFT, SF-TDDFT, and MRSF-TDDFT in OpenQP. It is a design and acceptance document, not a support claim for the unimplemented scientific kernels.

## Units and data contract

OpenQP Hessian drivers must return a real Cartesian nuclear Hessian with shape `(3N, 3N)` in the same flattened atom-major coordinate order consumed by the existing `normal_mode()` and thermochemistry path: `x1, y1, z1, x2, y2, z2, ...`. The matrix elements are derivatives of the electronic energy with respect to OpenQP nuclear coordinates in the units already expected by the numerical Hessian/frequency driver.

The final Python driver may symmetrize a Hessian before frequency analysis, but it must record or expose the pre-symmetrization maximum asymmetry. Symmetrization is a presentation/output safeguard, not a way to hide missing analytic terms. Analytic Hessian dispatch must use no silent numerical fallback: if a method/state is not implemented, `[hess] type=analytical` must fail explicitly with a method-specific diagnostic.

## HF analytic Hessian

The ground-state HF Hessian is the first native target. The implementation should reuse the existing SCF and HF-gradient data layout and accumulate:

- one-electron second derivative terms for overlap, kinetic, nuclear attraction, and available ECP contributions;
- two-electron Coulomb/exchange second derivative terms;
- orbital-response terms from CPHF/CPKS-like response equations;
- nuclear-repulsion second derivatives; and
- final packing into the common `(3N, 3N)` Hessian storage.

Initial validation must compare H2/H2O/formaldehyde HF Hessians against central finite differences of analytic gradients. Report max absolute error, RMS error, symmetry error, and displacement size.

### Native implementation progress

The native HF Hessian is being built term-by-term. Each term is validated in isolation before the full `[hess] type=analytical` HF path is unguarded, so an incomplete Hessian is never returned.

- **Done — nuclear-repulsion second derivatives.** `grd1::hess_nn` accumulates the `E_nn` second derivatives directly into the `(3N, 3N)` atom-major Hessian, reusing the effective-charge convention of `grd1::grad_nn`. Each pair `(k,l)` contributes the `3x3` block `Zk*Zl * (3 p_a p_b / r^5 - delta_ab / r^3)` to the `(k,k)`/`(l,l)` diagonal blocks and its negative to the `(k,l)`/`(l,k)` off-diagonal blocks. Validated by finite difference of `grd1::grad_nn` (Fortran link test `tests/fortran/test_hess_nn.F90`) and by a dependency-light NumPy reference (`tests/test_nuclear_repulsion_hessian.py`), which also checks Hessian symmetry and translational invariance.
- **In progress — one-electron `der2` integral terms** (`sum P * d2h - sum W * d2S`):
  - **Done — overlap and kinetic second derivatives.** `mod_1e_primitives::der2_kinovl_xyz` provides the 1D bra-center second-derivative recursion (`d2[j,i] = 4 ai^2 [j,i+2] - 2 ai (2i+1) [j,i] + i(i-1) [j,i-2]`), and `comp_overlap_der2` / `comp_kinetic_der2` assemble the 3x3 Cartesian blocks. The drivers `grd1::hess_ee_overlap` / `hess_ee_kinetic` scatter them into the `(3N,3N)` Hessian, reusing translational invariance of the two-center integrals and the same factor-2 convention as `grad_ee_overlap` / `grad_ee_kinetic`. Validated by (a) `tests/fortran/test_der2_recursion.F90` (the closed-form recursion equals composing the first-derivative recursion twice, to 1e-15), and (b) `hess1_selftest` (a `bind(C)` harness that finite-differences the production gradient routines on a real H2O/6-31G* basis: overlap and kinetic max|analytic-FD| ~ 1e-10, symmetric to ~1e-14).
  - **In progress — nuclear-attraction second derivatives.** The 1/r operator uses Rys quadrature and depends on three centers: bra A, ket B, and the nuclear charge C. The design avoids any second-derivative Rys-root machinery by reducing everything to the three basis-center blocks:
    - `d2/dA2`, `d2/dB2`, `d2/dA dB` are pure polynomial recursions on the Rys-expanded 1D integrals (raise/lower the bra and/or ket index per root), needing no derivatives of the Rys roots/weights.
    - the charge-center and mixed blocks follow from translational invariance `d/dC = -(d/dA + d/dB)`: e.g. `d2/dC2 = d2/dA2 + 2 d2/dA dB + d2/dB2`, `d2/dA dC = -(d2/dA2 + d2/dA dB)`.

    **Done so far:** `der2_coul_xyz` provides the 1D bra-center second-derivative recursion for the Rys integrals (same closed form as `der2_kinovl_xyz`), validated in `tests/fortran/test_der2_recursion.F90` (matches composing `der_coul_xyz` twice to ~1e-15). The Rys build must request `igrd=2` (so `xyzin` reaches bra index `iang+2`) and a correspondingly enlarged `xyzin` array.

    **Investigated and rejected approach (naive ket raising).** A first attempt computed the ket-center derivatives by directly raising the ket (first) index of `xyzin` with the ket exponent (`der_coul_xyz_j`), mirroring the validated bra recursion, then assembled the AA/AB/BB blocks and scattered the C/mixed blocks via translational invariance. A `bind(C)` FD harness on H2/STO-3G localized a bug: the **transverse (perpendicular-to-axis) Hessian elements matched FD exactly, but the on-axis (bond/charge-axis) elements were wrong** (e.g. zz analytic `-0.196` vs FD `-0.488`). Root cause: `QGaussRys` builds the integral asymmetrically (VRR in the ket index, then HRR to the bra), so the post-HRR ket index is *not* a clean ket angular momentum that the `2*aj` raising operator differentiates correctly. The bra raising is fine (it is exactly the production Pulay path), but naive ket raising is not a valid physical `d/dB`. The algebraic recursion self-check (`der2_j == der_j(der_j)`) passes but does not catch this, because it only verifies internal consistency, not physical correctness.

    **Corrected approach (implemented, WIP).** `d2/dB2` is obtained from the bra second derivative of the *swapped* shell pair (bra<->ket), reusing the validated bra path. The mixed `d2/dA dB` uses translational invariance `d2/dA dB = -(d2/dA2 + d2/dA dC)`, where `d2/dA dC` is the bra derivative of the Hellmann-Feynman charge term. Implemented as: `DQGaussRys` extended with an `igrd` parameter (so the charge-derivative build reaches bra index `iang+2`); `comp_coulomb_der2_braC` returns `p_AA = d2/dA2` (standard build) and `p_AC = d2/dA dC` (`DQGaussRys` + `der_helfey_xyz` + `der_coul_xyz`); and `hess_en` calls it twice per ordered pair (standard and swapped) and scatters all nine atom blocks.

    **Open issue (root cause found).** A separated finite-difference diagnostic (`hess1_selftest`, writing `/tmp/hess1_sep.out`) FDs the basis (Pulay) and charge (Hellmann-Feynman) gradient pieces against the basis and charge positions *independently* (possible because `grad_en_*` take the charge `coord` separately from `basis%atoms%xyz`). On H2/STO-3G this isolates the failure: the basis-basis **off-diagonal** (mixed bra-ket) block matches FD, but the **diagonal** `d2/dA2` blocks are too small (e.g. analytic 1.09 vs FD 1.49). So `der2_coul_xyz` -- the fixed-root angular-momentum recursion applied twice -- does **not** give the true second derivative of the nuclear-attraction (Rys, 1/r) integral, even though `der_coul_xyz` gives the exact first derivative. Bumping the Rys root count moves the answer *away* from the FD value, confirming this is not under-integration but an invalid scheme: the Rys roots/weights depend on the basis center (through `X = aa|P-c|^2`), and that root dependence contributes to the second derivative but is absent from a fixed-root double application. Overlap/kinetic worked precisely because they have no Rys roots.

    **Consequence.** The nuclear-attraction (and, by extension, the two-electron) second derivatives require a genuine second-derivative Rys scheme (root/weight derivatives, i.e. extending the Rys engine), not the fixed-root recursion-twice shortcut. The WIP `comp_coulomb_der2_braC` / `hess_en` and the `DQGaussRys` `igrd` extension are committed but remain unvalidated and unwired; `hf_hessian` stays guarded. This is the same Rys-engine work already flagged as required for the 2e term.

    **Reference blueprint (PySCF / libcint).** PySCF's `hessian/rhf.py` confirms this: it builds the core-Hamiltonian and nuclear-attraction Hessian from dedicated second-derivative integrals, never from first-derivative recursion applied twice. The relevant libcint intors and their OpenQP targets:
    - `int1e_ipipnuc` (comp=9) = `d2/dA2` nuclear attraction (both derivatives on the bra center) -- the term `der2_coul_xyz` fails to reproduce.
    - `int1e_ipnucip` (comp=9) = `d2/dA dB` (one derivative on bra, one on ket).
    - `int1e_ipipkin` / `int1e_ipkinip`, `int1e_ipipovlp` / `int1e_ipovlpip` = the kinetic and overlap second derivatives (OpenQP already has working equivalents for these via `der2_kinovl_xyz`, since they carry no Rys roots).
    - `int1e_ipiprinv` / `int1e_iprinvip` = second derivatives of `1/|r-C|` about a nucleus (the charge-center / Hellmann-Feynman terms), obtained per nucleus with `with_rinv_at_nucleus`.

    Assembly: `get_hcore` returns the basis-center blocks (`h1aa`=ipipnuc+ipipkin, `h1ab`=ipnucip+ipkinip); `hcore_generator(iatm,jatm)` adds the per-nucleus `rinv` blocks scaled by the nuclear charge, combining the basis-center and charge-center pieces with the translational-invariance structure (`hcore = -rinv2aa - rinv2ab` plus bra/ket add-backs). OpenQP needs the equivalents of `ipipnuc` / `ipnucip` / `ipiprinv` / `iprinvip` -- i.e. real second-derivative Rys integrals -- after which the `hess_en` driver's block bookkeeping can be reused.

    **Chosen integral backend: libint (deriv order 2).** Rather than extend OpenQP's own Rys engine for second derivatives or add libcint, the native Hessian will use the libint integrals OpenQP already links. The libint wrapper (`int_libint.F90`) already exposes a `deriv_order` interface and a `libint2_build_eri2` (second-derivative ERI) code path gated by `#if INCLUDE_ERI >= 2`, and the Rys and libint 2e engines already coexist and are selected at runtime via `libint2_active` (see `int2.F90`), so turning libint on does not remove the Rys path. Plan: (a) regenerate the custom `oqp-libint` tarball with deriv order 2 and build with `INCLUDE_ERI=2`; (b) use libint's second-derivative 1e and 2e integrals (`ipipnuc`/`ipnucip`/`ipiprinv` and `eri2`) to feed the already-validated `hess_en`/Hessian block bookkeeping; (c) keep `USE_LIBINT=OFF` as a lighter default build in which the native analytic Hessian is unavailable and dispatch falls back, with an explicit error, to the PySCF bridge or numerical Hessian (no silent fallback). PySCF/libcint remains the *test-time* validation oracle (element-wise vs `int1e_ipipnuc`, `int2e_ipip1`), not a runtime dependency.
- **Deferred — two-electron `d2(uv|ls)` terms and the CPHF orbital-response term.** The two-electron second-derivative integrals require extending the Rys engine (`grd2_rys`) beyond its current first-derivative-only (`nder = 1`) form; the CPHF solve can reuse the `pcg` solver following the `tdhf_z_vector` pattern.

Until the electronic terms above are implemented and validated, the native `hf_hessian` kernel remains an explicit guarded scaffold (it aborts rather than returning a partial matrix), and analytical HF/DFT Hessian requests are served by the clearly-labeled external PySCF bridge with no silent numerical fallback.

## External PySCF bridge (supported delivery path)

The external PySCF bridge (`pyoqp/oqp/library/external.py::analytic_hessian_from_pyscf`) is the supported way to obtain analytic HF and DFT Hessians today. It builds a guarded PySCF mean-field object from the OpenQP molecule, runs the PySCF analytic Hessian, and returns it in OpenQP `(3N, 3N)` atom-major convention; `[hess] type=analytical` dispatch routes ground-state HF/DFT here with no silent numerical fallback. Coverage: `scf.type` `rhf` and `uhf` (RKS/UKS for DFT). `rohf` raises explicitly because PySCF provides no ROHF analytic Hessian; TDDFT/SF/MRSF remain `NotImplementedError`.

Three latent bugs that made this path non-functional were fixed and the result is validated element-wise against direct PySCF (`tests/test_pyscf_hessian_bridge.py`, RHF and UHF, `max|dH| < 1e-8`):
- the bridge was gated off by `mol.usempi`, which is hard-coded `True` even in serial; the guard now checks the MPI manager's real `use_mpi` flag.
- the geometry was double-converted -- `mol.get_system()` already returns Bohr, and the `ANGSTROM_TO_BOHR` constant is actually the Bohr->Angstrom factor (0.529177); the coordinates are now passed straight through to PySCF's `unit='Bohr'`.
- the molden frequency writer crashed on a `(-1, 1)`-reshaped `freqs` (`moldenwriter.write_frequency`); it now flattens before formatting.

## DFT analytic Hessian

DFT builds on the HF Hessian structure but adds exchange-correlation grid derivatives and Kohn-Sham response terms. The first accepted scope should explicitly name supported functional classes. If meta-GGA, unsupported grid derivatives, or unavailable ECP second derivatives are not ready, input checking must reject them rather than proceeding with incomplete terms.

DFT validation mirrors HF validation with small molecules and stable functionals, starting from a tiny basis and later checking canonical user-facing cases such as BHHLYP/6-31g* only after the small finite-difference matrix is clean.

## TDDFT analytic Hessian

Conventional TDDFT/TDA/RPA Hessians are state-specific excited-state Hessians. They must route through `[input] method=tdhf`, `[tdhf] type=tda` or `rpa`, and `[hess] state=N`, where state 0 remains the reference/ground state and excited states use positive indices.

The implementation should mirror the TDDFT gradient/Z-vector flow: response amplitudes, relaxed density/intermediate construction, response-vector derivative terms, XC response-kernel terms for DFT, and state-specific Hessian accumulation. TDA should be validated before full RPA so X/Y coupling terms cannot be accidentally skipped.

## MRSF-TDDFT analytic Hessian

MRSF-TDDFT Hessians must be treated as a later, separately validated tier rather than as a consequence of the HF/DFT or conventional TDDFT scaffolds. The response-root convention is OpenQP/MRSF-specific: root 0 is the high-spin ROHF reference, root 1 maps to physical S0, root 2 maps to physical S1, and so on. Hessian tests, diagnostics, and output must preserve that mapping and must report `<S^2>` for the tracked target root.

The MRSF implementation must start from a validated MRSF gradient/Z-vector baseline. Before enabling a native MRSF Hessian kernel, finite-difference checks must show stable analytic gradients for multiple molecules and roots, with no-fix controls for known SPC/Z-vector density terms. The Hessian accumulation should reuse the established spin-adapted density channels and MRSF Z-vector intermediates exactly; it must not reinterpret ROHF-MRSF as a general UHF path or substitute `mo_a`/`mo_b` changes without separate UMRSF validation.

Initial MRSF Hessian support remains explicitly disabled until those gradient and root-tracking checks are clean. A source-level scaffold may expose ABI and dispatch guardrails, but runtime `[hess] type=analytical` for `tdhf.type=mrsf` must fail explicitly rather than returning zeros or falling back to numerical displacements.

## Unsupported first-release cases

The first release should explicitly reject or defer unsupported features, including unavailable ECP second derivatives, unvalidated functional classes, state crossings/root flips, noncollinear or spin-orbit Hessians, DFTB+ Hessians, UMRSF gradients/Hessians, and any TD/SF case that lacks finite-difference validation evidence.

## Acceptance evidence

Each staged PR must include tests that fail before implementation, pass after the minimal change, and remain dependency-light where possible. Scientific support claims require live finite-difference validation against analytic gradients; source-level scaffolding and static tests only prove routing and guardrails, not Hessian correctness.
