# UMRSF-TDDFT Analytic Energy Gradient

Analytic nuclear energy gradients for UHF-referenced MRSF-TDDFT
(UMRSF-TDDFT) in OpenQP.

## Overview
Excited-state analytic energy gradients for a spin-unrestricted (UHF)
MRSF-TDDFT reference, via a Lagrangian / Z-vector formulation. Supports
HF, hybrid, and pure-GGA functionals.

## Validation
Gradients agree with reference to <=1e-5 Ha/Bohr across
{HF, hybrid, pure-GGA} x {common Pople/Dunning basis sets}, for the
low-lying singlet states (S1/S2/S3), up to production size (nbf ~150).

For near-degenerate states, gradients follow state character
(transition-density overlap + spin) rather than energy ordering.

## Implementation
Main module: source/modules/tdhf_umrsf_gradient.F90

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
