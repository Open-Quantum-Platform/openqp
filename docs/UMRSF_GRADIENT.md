# UMRSF-TDDFT Analytic Energy Gradient

Analytic nuclear energy gradients for UHF-referenced MRSF-TDDFT
(UMRSF-TDDFT) in OpenQP.

## Overview
Excited-state analytic energy gradients for a spin-unrestricted (UHF)
MRSF-TDDFT reference, via a Lagrangian / Z-vector formulation. The HF
limit is validated. Hybrid and pure-GGA paths include the XC response
terms, but remain under finite-difference validation.

## Validation
HF-limit gradients agree with finite-difference references to
<=1e-5 Ha/Bohr for the validated low-lying singlet states.

For near-degenerate states, gradients follow state character
(transition-density overlap + spin) rather than energy ordering.

Current DFT/XC status: the coupled alpha/beta Z-vector solve is active
and agrees with the dense oracle on H2O/BHHLYP, but the total analytic
gradient still shows a small XC-related finite-difference residual on
that test case. Treat DFT UMRSF gradients as experimental until that
residual is closed.

## Implementation
Response/Z-vector module: source/modules/tdhf_umrsf_z_vector.F90

Gradient assembly module: source/modules/tdhf_umrsf_gradient.F90

The standard TD-gradient pipeline first calls the UMRSF Z-vector entry
point, which solves the spin-coupled alpha/beta UMRSF response and
caches the response-gradient contribution. The alpha and beta MO
rotations are separate unknown blocks in one coupled solve; the
mean-field/XC response couples those spin blocks. The gradient entry
point then evaluates the reference UHF-triplet gradient and adds the
cached response term.

A UMRSF-TDDFT gradient is requested through the standard OpenQP input
deck. The reference is the UHF triplet (`[scf] type=uhf`,
`multiplicity=3`), the response is `[tdhf] type=umrsf`, and the gradient
of a chosen excited state is selected with `runtype=grad` plus
`[properties] grad=`:

    [input]
    runtype=grad
    method=tdhf
    functional=bhhlyp        # omit for the HF limit; any hybrid or pure-GGA name
    basis=6-31g*
    system=
       6   0.000000   0.000000   0.000000
       # ... Z  x  y  z (Angstrom), one line per atom ...

    [scf]
    type=uhf
    multiplicity=3

    [tdhf]
    type=umrsf
    nstate=10

    [properties]
    grad=1                   # 1-based index of the excited state to differentiate

`nstate` must be large enough to include the requested state. Multiple
states may be listed under `grad=` (comma-separated). Gradients are
returned in Hartree/Bohr.
