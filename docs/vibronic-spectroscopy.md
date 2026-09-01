# Harmonic vibronic spectroscopy Python API

OpenQP provides a post-processing model for Franck--Condon (FC) and
first-order Herzberg--Teller (HT) vibronic spectra from two harmonic surfaces.
It also evaluates double-harmonic excited-state infrared intensities and
nonresonant Raman activities when the corresponding target-state property
derivatives are supplied explicitly.

The implementation does not represent MRSF states as determinants.  Electronic
structure enters only through state-specific energies, Hessians, transition
dipoles, and property derivatives supplied by the caller.

## Harmonic FC model

The two sets of mass-weighted normal coordinates obey

\[
Q_e = J Q_g + K .
\]

`HarmonicOverlapEngine` evaluates the exact multidimensional generating
function for the two harmonic oscillator bases.  It includes frequency change,
Duschinsky mixing, and displacement.  It does not use an independent-mode
product approximation.

```python
from oqp.library import (
    HarmonicVibronicModel,
    harmonic_vibronic_spectrum,
)

model = HarmonicVibronicModel.create(
    ground_frequencies_cm1=[1595.0, 3657.0, 3756.0],
    excited_frequencies_cm1=[1510.0, 3510.0, 3650.0],
    J=[[0.9982005399, -0.0599640065, 0.0],
       [0.0599640065,  0.9982005399, 0.0],
       [0.0,           0.0,          1.0]],
    K=[0.12, -0.04, 0.02],
    coordinate_unit="sqrt(amu)*bohr",
    coordinate_phase_convention="OpenQP Duschinsky aligned normal modes",
    orthogonality_tolerance=1.0e-8,
)

spectrum = harmonic_vibronic_spectrum(
    model,
    electronic_origin_cm1=32000.0,
    origin_kind="zero_zero",
    max_final_quanta=7,
    temperature_kelvin=300.0,
    max_initial_quanta=4,
    broadening_kind="gaussian",
    fwhm_cm1=120.0,
    normalization="sum",
)
```

Positive frequencies are required.  `J` must be square and orthogonal within
the stated tolerance.  The implementation reports and rejects a noncanonical
transformation; it does not replace `J` by a nearest orthogonal matrix.

At nonzero temperature, `max_initial_quanta` is mandatory.  Initial-state
populations are normalized by the exact harmonic partition function, rather
than by the truncated state list.  `retained_thermal_population` therefore
measures the population retained by the selected basis.  By default the
calculation fails if this value is below 0.999.  The final-state truncation is
reported separately as `franck_condon_completeness`.  By default the
calculation also fails when that completeness is below 0.999, before any
unit-sum normalization could conceal omitted vibronic strength.  An explicitly
lower `minimum_franck_condon_completeness` permits a deliberately partial band.

`origin_kind="zero_zero"` interprets `electronic_origin_cm1` as the 0--0 band
origin.  `origin_kind="adiabatic_minima"` interprets it as the energy difference
between the two potential-energy minima and adds the two zero-point energies.
Gaussian and Lorentzian profiles are normalized to unit area.  With
`normalization="sum"`, the final numerical spectrum also integrates to one on
the returned grid.  These are relative vibronic strengths, not absorption
cross sections.

## First-order HT term

The first-order transition moment is

\[
\boldsymbol{M}_{v_e v_g} =
\boldsymbol{\mu}_0\langle v_e|v_g\rangle +
\sum_k \left(\frac{\partial\boldsymbol{\mu}}{\partial Q_{g,k}}\right)_0
\langle v_e|Q_{g,k}|v_g\rangle .
\]

Both the equilibrium transition dipole and its derivative carry the electronic
state label and electronic-state phase convention.  The derivative separately
records the normal-coordinate phase convention, which must match the harmonic
model.  It must use the ground-state normal coordinates and the same coordinate
unit as `K`.

```python
from oqp.library import (
    ElectronicTransitionMoment,
    ExcitedStatePropertyDerivative,
)

transition = ElectronicTransitionMoment.create(
    [0.03, 0.00, 0.01],
    electronic_state="MRSF S1",
    electronic_phase_convention="overlap-aligned MRSF S1",
)

transition_gradient = ExcitedStatePropertyDerivative.create(
    transition_dipole_derivative,       # shape (3, nmode)
    property_kind="transition_dipole",
    coordinate_basis="ground_normal",
    coordinate_unit="sqrt(amu)*bohr",
    property_unit="e*bohr",
    electronic_state="MRSF S1",
    electronic_phase_convention="overlap-aligned MRSF S1",
    coordinate_phase_convention="OpenQP Duschinsky aligned normal modes",
    provenance="central difference of tracked MRSF transition dipoles",
    electronic_state_role="target_excited_state",
)
```

Passing an HT derivative without an equilibrium transition dipole is rejected.
For an electronically forbidden transition, supply an explicit zero dipole so
that the intended Condon limit is unambiguous.  A mismatch in state label,
coordinate unit, electronic phase, or normal-coordinate phase is also rejected.

## State-tracked numerical property derivatives

`finite_difference_excited_state_property` forms a central difference from
property values already calculated at the plus and minus displacements.  It
requires state-overlap diagnostics for every displacement.  Transition-dipole
differences additionally require explicit unit-modulus phase factors for both
geometries.

Each supplied phase factor multiplies the corresponding raw transition moment
before differencing and must align that complete transition moment with the
central-geometry convention.  It must therefore include every electronic
phase relevant to the bra and ket; an excited-state overlap phase alone is not
sufficient when the reference-state phase has not already been fixed.

This adapter checks the supplied isolated-state overlap magnitudes, but does
not itself establish root identity.  An MRSF calculation must additionally
retain agreement in energy, spatial symmetry, spin sector, and the complete
spin-adapted response vector.  A degenerate or near-degenerate manifold
requires projector tracking and is outside this isolated-state adapter.

The adapter performs no electronic-structure calculation.  In particular, it
cannot and does not replace a missing MRSF target-state property by the ROHF or
SCF reference property.

## Analytic MRSF property-derivative hook

`mrsf_analytic_property_derivative_from_provider` is the fail-closed entry for
an OpenQP analytic first-order property implementation.  The provider must
implement `get_mrsf_analytic_property_derivative(property_kind=...,
electronic_state=...)` and return the derivative tensor together with the
electronic method (`MRSF-TDDFT` or `MRSF-TDHF`), derivative order, target-state
role, units, coordinate basis, electronic phase, normal-coordinate phase, and
analytic provenance.

The hook never queries a generic dipole, polarizability, ROHF, or SCF property
accessor.  Consequently, an OpenQP build without the dedicated MRSF analytic
property method raises `TypeError`; it cannot silently label a reference-state
property as an MRSF result.  The present code defines and verifies this public
contract, but it does not claim that the underlying analytic MRSF dipole or
polarizability response equations are already available.

## Excited-state IR and nonresonant Raman

`excited_state_ir_intensities` accepts a target-state permanent-dipole
derivative in excited-state normal coordinates, with shape `(3, nmode)`, and
returns double-harmonic intensities in km mol^-1.  The required derivative unit
is `e*bohr` per `sqrt(amu)*bohr`, numerically equivalent to
`e/sqrt(amu)`; its intensity conversion factor is 974.88011 km mol^-1.

`excited_state_raman_activities` accepts a target-state static-polarizability
derivative with shape `(3, 3, nmode)`.  It returns

\[
S_k = 45\bar{\alpha}_k'^2 + 7\gamma_k'^2
\]

in `bohr^4/amu`, together with the polarized and unpolarized depolarization
ratios.  The nonresonant polarizability derivative must be real and symmetric.

These two functions report fundamental-mode intensities in the double-harmonic
or Placzek limit.  A mode is IR inactive when its target-state dipole derivative
vector is zero, and nonresonant Raman inactive when both polarizability
invariants vanish.  Overtones, combination bands, anharmonic intensity
borrowing, and resonance enhancement are outside these two functions.

Both functions require the explicit declaration
`electronic_state_role="target_excited_state"`.  The literal roles `reference_state`, `ROHF`,
`SCF`, and unspecified state properties are rejected.

These functions evaluate spectra from supplied property derivatives.  They do
not constitute analytic MRSF excited-state dipole or polarizability derivatives.

## Resonance-Raman proof of principle

`resonance_raman_fc_ht` evaluates the resonant term of a harmonic
Kramers--Heisenberg--Dirac sum over intermediate excited-state vibrational
levels.  It requires the electronic origin, laser wavenumber, damping, phase-
consistent transition dipole, and optional first-order HT derivative.  The
result is a complex scattering tensor for each requested final vibrational
state in `(e*bohr)^2/cm^-1`.  `damping_cm1` is the positive imaginary-energy
parameter used directly in the resonant denominator, not a separately
converted full width at half maximum.

The intermediate-state truncation is compared with the exact closure sum for
the supplied constant-plus-linear transition dipole.  By default the
calculation fails if the retained transition strength is below 0.999.  The
returned `intermediate_transition_completeness` records this fraction.  This
closure condition does not replace an explicit convergence comparison of the
complex tensor with respect to `max_intermediate_quanta`, particularly near a
vibronic resonance.

This deliberately limited model omits the antiresonant term, higher-order
non-Condon terms, coordinate dependence beyond first order, and absolute
cross-section prefactors.  It must not be described as a complete resonance-
Raman calculation.

## JSON command-line calculation

A small FC calculation can be evaluated without an OpenQP input keyword:

```bash
python -m oqp.library.vibronic input.json output.json
```

The record in `examples/VIBRONIC/HARMONIC_FC.json` gives the complete input
schema.  The output contains stick factors, broadened intensities, thermal
population retention, and FC completeness.

## Scope boundary

This API provides:

- multidimensional harmonic FC factors and relative absorption spectra;
- first-order HT transition moments from an explicit transition-dipole
  derivative;
- target-excited-state double-harmonic IR intensities from an explicit dipole
  derivative;
- target-excited-state nonresonant Raman activities from an explicit static-
  polarizability derivative; and
- a resonant-term harmonic FC/HT Raman tensor for method assessment.

It does not provide anharmonic vibrational states, an internally calculated
analytic MRSF transition-dipole derivative, an analytic target-state MRSF
polarizability derivative, absolute absorption or Raman cross sections, or a
complete resonance-Raman treatment.

The derivative and IR, Raman, and resonance-Raman result objects provide
`to_dict()` records that retain units, model labels, and the stated scope.  The
records are directly serializable with `json.dumps`.
