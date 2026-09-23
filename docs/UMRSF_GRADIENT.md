# UMRSF-TDDFT Analytic Energy Gradient

Analytic nuclear energy gradients for UHF-referenced MRSF-TDDFT
(UMRSF-TDDFT) in OpenQP.

## Overview
Excited-state analytic energy gradients for a spin-unrestricted (UHF)
MRSF-TDDFT reference, via a Lagrangian / Z-vector formulation. The
implemented equations cover the HF limit and full-range LDA/GGA
functionals, including conventional global hybrids. DFT gradients include
the reference UKS XC-kernel response and the fixed-weight explicit XC
skeleton derivative.

## Validation
HF-limit gradients agree with character-followed finite-difference
references to <=1e-5 Ha/Bohr for the validated low-lying singlet states.
The SOMO-corrected response density also closes the former S2 structural
error.

On 2026-08-02, C1-distorted H2CO/6-31G* with a UHF-triplet reference and
singlet UMRSF target root 1 (`nstate=3`) was compared with two-point central
finite differences at `dx=1e-3` Bohr. The calculation used a 96-radial by
302-angular per-atom grid with quadrature and AO pruning disabled. All 12
Cartesian components were nonzero; the maximum component errors were
1.39e-5 Ha/Bohr for BHHLYP and 2.52e-5 Ha/Bohr for BLYP. These two cases do
not establish support for every LibXC functional in the implemented classes.
The reproducible commands are:

    python tools/validate_gradients.py umrsf-bhhlyp --nstate 1 --dx 1e-3 --tol 1e-4
    python tools/validate_gradients.py umrsf-blyp   --nstate 1 --dx 1e-3 --tol 1e-4

Diagnostic internal oracles on the BHHLYP case give a maximum Cartesian
component discrepancy of 2.01e-9 Ha/Bohr for the frozen-density 11-channel
two-electron derivative and a largest alpha/beta elementwise discrepancy of
3.92e-11 for the analytic one-sided two-electron generalized Fock. The full
coupled alpha/beta Z-vector stationarity residual was 3.38e-10 against the
requested 1e-9 tolerance.

An analytic gradient is evaluated for the selected response eigenvector at
one geometry. Across displaced geometries, energy-index labels can exchange
character. Finite-difference validation of near-degenerate or higher roots
therefore requires overlap/character following; the supplied reproducible
DFT regressions intentionally check only separated root 1.

Range-separated CAM/LRC, meta-GGA, double-hybrid, and explicitly overridden
or functional-specific `spc_coco`/`spc_ovov`/`spc_coov` parameterizations are
outside this analytic-gradient formulation and are rejected for
gradient-driven UMRSF runtypes. They remain available for energy calculations.

OpenQP's XC-gradient consumers currently omit derivatives of the
atom-partition moving-grid quadrature weights. The electronic XC response is
included, but the omitted weight derivative leaves a grid-dependent
analytic-versus-finite-difference floor. Use a dense unpruned grid with AO
pruning disabled for validation. UMRSF Hessians are not implemented; both
analytical and numerical `runtype=hess` requests are currently unsupported.

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

The dedicated UMRSF solver first eliminates the diagonal occupied-occupied
and virtual-virtual blocks, then applies a matrix-free PCG trial to the
symmetric occupied-virtual block with automatic MINRES fallback. `maxit_zv`
is the total PCG-plus-MINRES iteration budget. `zvconv` is checked against an
explicit Euclidean relative residual for both the reduced solve and the full
coupled alpha/beta equation; an unconverged response aborts rather than being
cached. `[tdhf] z_solver` is not consulted by this dedicated path and should
be left at its default value.

With the implementation's one-sided derivative
`G_pq = d omega / d eta_pq` for `C_q += eta_pq C_p`, the corresponding
orthogonal-angle gradient is `g = -antisym(G)`. The Z-vector operator emitted
by the generalized Fock of the multiplier scalar is already the adjoint
canonicality action `M = -H^T`; therefore the implemented `M z = -R` solve is
the Lagrangian equation `H^T z = -g`, not a forward-Jacobian solve.

A UMRSF-TDDFT gradient is requested through the standard OpenQP input
deck. The reference is the UHF triplet (`[scf] type=uhf`,
`multiplicity=3`), the response is `[tdhf] type=umrsf`, and the gradient
of a chosen excited state is selected with `runtype=grad` plus
`[properties] grad=`:

    [input]
    runtype=grad
    method=tdhf
    functional=bhhlyp        # omit for HF; use semilocal or global-hybrid LDA/GGA
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
