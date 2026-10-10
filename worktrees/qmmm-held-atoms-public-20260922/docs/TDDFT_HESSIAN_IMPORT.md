# Import of the closed-shell TDHF/TDDFT nuclear Hessian

## Provenance and scope

The present analytic TDHF/TDDFT nuclear-Hessian methodology and its
implementation effort were initiated by Hiroya Nakata. This OpenQP development
continues that work and records his originating contribution explicitly.

This development is based on the supplied GAMESS implementation and its
verification record. The principal TDHF/TDDFT Hessian routines are in
`tddgrd.src`, with connections and second-derivative terms in `hess.src`,
`cphf.src`, `hss1a.src`, `hss2a.src`, and related source files. The reference
documents and sources describe another program and were used only as
scientific and historical material for the independent OpenQP implementation.

The verified GAMESS limit defines the initial OpenQP scope:

- closed-shell restricted reference;
- singlet TDHF or TDDFT;
- full linear response, not the Tamm--Dancoff approximation;
- Hartree--Fock, LDA, GGA, and hybrid-GGA exchange-correlation models;
- one MPI rank.

Open-shell references, triplet targets, Tamm--Dancoff response, meta-GGA
functionals, and distributed evaluation remain unsupported until their
response equations and numerical criteria are established independently.

## Mathematical decomposition

For one excited state with energy

\[
E_I(\mathbf R) = E_0(\mathbf R) + \omega_I(\mathbf R),
\]

the Cartesian Hessian is

\[
\frac{d^2 E_I}{dR_K dR_L} =
\frac{d^2 E_0}{dR_K dR_L} +
\frac{d^2 \omega_I}{dR_K dR_L}.
\]

Stationarity of the excited-state Lagrangian separates the excitation-energy
term into symmetric second-derivative contractions and directional response
rows,

\[
\frac{d^2 \omega_I}{dR_K dR_L} =
H^{1e}_{KL} + H^{2e}_{KL} + H^{xc}_{KL}
+ \frac{G_{KL}+G_{LK}}{2}.
\]

`assemble_tdhf_cartesian_hessian` accepts the already combined total
fixed-density contribution, \(H^{fixed}=d^2E_0+H^{1e}+H^{2e}\), and adds
\(H^{xc}\) and the symmetrized response rows. This interface prevents the
ground-state term from being added twice. The directional matrix `G` must
remain unsymmetrized until final assembly: its
antisymmetric part is a convergence and completeness diagnostic for the
coupled response equations.

## Correspondence between the GAMESS derivation and OpenQP

| Quantity | GAMESS description | OpenQP starting point |
| --- | --- | --- |
| \(d^2E_0\) | ground-state FCM | `hf_hessian_mod::hf_hessian` |
| relaxed difference density and energy-weighted density | `PEFF`, `WEFF` | `OQP::td_p`, `OQP::WAO`, and `tdhf_z_vector_mod` |
| transition amplitudes | paired \(U=X+Y\), \(V=X-Y\) | `OQP::td_xpy`, `OQP::td_xmy` |
| one-electron second derivatives | overlap, kinetic, nuclear attraction | `grd1` Hessian routines |
| two-electron second derivatives | relaxed two-particle contraction | `grd2_hess_driver` with a TDDFT density provider |
| ground-state orbital response | CPKS | `cphf_mod::cphf_solve` |
| amplitude response | CP-TDHF/TDDFT | new paired-response solve using the existing TD sigma operation |
| derivative Z-vector | derivative of \((A+B)Z=R\) | new response solve using `tdhf_z_vector_mod` conventions |
| XC second derivatives and response rows | moving-grid \(f_{xc}\), \(k_{xc}\), and fourth functional derivatives | extend `mod_dft_gridint_tdxc_grad` after the TDHF terms agree |

## Order of implementation and verification

1. Preserve the current analytic ground-state Hessian as \(d^2E_0\).
2. Form the closed-shell TDHF one- and two-electron second-derivative
   contractions with fixed relaxed densities. The initial implementation is
   `tdhf_hessian_fixed_density_mod::build_tdhf_fixed_density_hessian`; it returns
   the total TDHF fixed-density contribution and explicitly rejects DFT until
   the XC quadrature term is present.
3. Solve the ground-state orbital response, paired amplitude response, and
   derivative Z-vector for every Cartesian perturbation. The projected paired
   solve and analytic coordinate derivatives of both response operators are
   implemented in the TDHF Hessian response modules.
4. Assemble and symmetrize the directional response rows only at the final
   sum. Record their maximum antisymmetric element before symmetrization.
5. Compare every Cartesian element of the TDHF Hessian with a central
   difference of the analytic OpenQP excited-state gradient at several steps.
6. Add LDA, then GGA and hybrid-GGA quadrature terms. Verify translation and
   rotation sum rules and compare against the supplied GAMESS/PySCF reference
   data.

Agreement of frequencies alone is insufficient: all Cartesian elements,
directional-row asymmetry, and reference-free translation and rotation sums
must satisfy stated numerical thresholds.
