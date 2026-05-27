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
