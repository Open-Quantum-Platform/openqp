# Lagged Baeck-An time-derivative coupling for FSSH

## Observable and approximation

This option approximates the scalar time-derivative coupling used in the
electronic propagation of same-spin FSSH trajectories.  It does not construct
a nuclear-coordinate nonadiabatic-coupling vector and does not alter the MRSF
response space, state energies, or active-state gradient.

For the adiabatic energy gap

\[
\Delta E_{ij}(t)=E_i(t)-E_j(t),
\]

the magnitude at the centre of three consecutive energy points is

\[
\left|\tau^{\mathrm{BA}}_{ij}(t_n)\right|
=\frac{1}{2}\sqrt{
\frac{\mathrm d^2\Delta E_{ij}(t_n)/\mathrm dt^2}
     {\Delta E_{ij}(t_n)}}.
\]

The nonuniform three-point curvature already used by the resident Fortran
diagnostic is retained.  A pair is set to zero when the radicand is nonpositive,
the centre gap is zero, or the centre gap exceeds `ba_gap_max`.

## Time alignment and sign

The centred curvature at `t_n` becomes available only after the energies at
`t_(n+1)` have been evaluated.  Causal dynamics therefore uses it during the
electronic propagation at `t_(n+1)`.  The method is explicitly one nuclear step
lagged; it is not described as an instantaneous Baeck-An coupling.

Energies alone do not determine the wavefunction gauge.  The magnitude above
receives the pairwise sign of the phase-tracked, centred overlap TDC.  A pair
whose overlap sign is exactly indeterminate is set to zero.  The first interval
and an interval following a discontinuous history use the overlap NPI TDC as a
warm-up value.  Dense trajectory records distinguish the sources:

- `tdc_source=1`: overlap NPI warm-up;
- `tdc_source=3`: lagged Baeck-An magnitude with overlap-transported sign.

The production comparison uses `tdc=baeck_an` and `rescale=isotropic` so that it
tests the approximate electronic coupling without adding an analytic NAC-vector
calculation at a hop.  Energy-conserving isotropic momentum adjustment remains
part of the accepted-hop treatment.

## Scope and limitations

- same-spin FSSH only;
- no full `3N` coupling vector or direction-specific momentum adjustment;
- no claim that the energy-only expression supplies a Berry phase or a signed
  electronic-state gauge;
- one-step time lag and the `ba_gap_max` pair selection are part of the named
  approximation and must be preserved in comparisons and restart signatures.
