# SC-NEVPT2 Nuclear Gradient

OpenQP computes the nuclear gradient of the single-state, strongly contracted
NEVPT2 (SC-NEVPT2) total energy analytically. The implemented energy is the
Dyall-Hamiltonian, strongly contracted correction to a state-specific CASSCF
reference. Set `[input] method=caspt2`, `[pt2] h0=dyall`,
`contraction=strong`, and use a gradient-driven run type.

Unlike CASSCF, SC-NEVPT2 is not variational with respect to the reference
orbitals or CI coefficients. Its derivative therefore includes coupled orbital
and CI response, the response of the semicanonical inactive and virtual
orbitals, and the derivative-integral contraction of the resulting relaxed
one- and two-particle densities. The public energy and gradient calls share one
CASSCF reference, one second-order evaluation, and one coupled response solve
at each geometry.

## Selecting the route

`[pt2] gradient` controls the nuclear-gradient route:

| Value | Behavior |
| --- | --- |
| `auto` | Use the analytic derivative when the input and current reference satisfy its scope; otherwise report the reason and use central differences. This is the default. |
| `analytic` | Require the analytic derivative. An unsupported input, non-stationary reference, ill-conditioned semicanonical response, or unsolved coupled response is an error. |
| `numerical` | Force Cartesian central differences of the converged SC-NEVPT2 energy. |

The analytic route applies only when all of the following are true:

- `method=caspt2`, `variant=caspt2`, `h0=dyall`, and
  `contraction=strong` select single-state SC-NEVPT2;
- `reference=casscf` selects a state-specific, stationary CASSCF reference;
- the reference is a closed-shell singlet from RHF orbitals, with no density
  functional and a contiguous inactive/active/virtual partition;
- the native determinant-space integral backend and compiled
  `nevpt2_gradient` entry point are available;
- no IPEA, real, imaginary, or QDPT energy-denominator shift is requested; and
- the configured `[cas]` or `[pt2] max_memory` can hold the active-space,
  orbital-space, and native Cartesian density intermediates.

`auto` does not hide an implementation failure. It falls back only for a
declared scope or reference-applicability condition; an unexpected native or
Python failure remains an error.

## Input example

The concise `.oqp` request — the SC-NEVPT2 route is the `caspt2` model with the
Dyall zeroth-order Hamiltonian and strong contraction selected in `pt2(...)`:

```text
caspt2/sto-3g opt(S0) cas(active_electrons=2,active_orbitals=2,frozen_core=1)
pt2(h0=dyall,contraction=strong,gradient=analytic)
geom="""
Li  0.0  0.0  0.0
H   0.0  0.0  1.6
"""
```

The equivalent legacy `.inp` spelling is:

```ini
[input]
system=
  Li  0.0  0.0  0.0
  H   0.0  0.0  1.6
basis=sto-3g
method=caspt2
runtype=optimize

[cas]
active_electrons=2
active_orbitals=2
frozen_core=1

[ci]
target_spin=singlet

[pt2]
reference=casscf
h0=dyall
contraction=strong
gradient=analytic

[optimize]
istate=0
```

The runnable input is
[`examples/WF_methods/LiH_SC-NEVPT2_optimize.inp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/LiH_SC-NEVPT2_optimize.inp).
Single-state SC-NEVPT2 publishes one energy and one gradient row, so
`[optimize] istate=0` is the public selector. Use `[pt2] target_roots` to select
the physical CASSCF root corrected by the single-state calculation.

The equivalent Python helper exposes the same route:

```python
from oqp.openqp import OpenQP

job = OpenQP("lih_sc_nevpt2_opt", silent=1)
job.molecule("Li 0 0 0; H 0 0 1.6", charge=0, multiplicity=1)
job.nevpt2(active_electrons=2, active_orbitals=2, frozen_core=1,
           contraction="strong", gradient="analytic",
           basis="sto-3g", runtype="optimize",
           casscf={"gradient_norm_tol": "1.0e-10"})
mol = job.run()
```

## Stationarity and response diagnostics

The SC-NEVPT2 Lagrangian assumes a stationary CASSCF reference. Its gradient
error is first order in the residual CASSCF orbital-rotation gradient, so
`[casscf] gradient_norm_tol` controls the useful precision of the nuclear
gradient. Each analytic evaluation reports:

- the CASSCF orbital-gradient norm;
- the smallest semicanonical inactive or virtual orbital-energy gap;
- the number of exactly degenerate pairs treated as free orbital gauges;
- the smallest Dyall denominator;
- the relative residual of the coupled orbital/CI response solve; and
- the maximum asymmetry of the total orbital Lagrangian.

An exactly degenerate inactive or virtual pair is accepted only when the
corresponding Lagrangian numerator also vanishes at working precision. The
rotation then is a free gauge and its multiplier is zero. A merely near-
degenerate pair remains ill-conditioned and is refused by the analytic route.
This distinction permits symmetry-required degeneracies, such as the two
virtual pi orbitals of linear LiH, without assigning an unstable finite
multiplier.

The native derivative-integral contraction returns the complete relaxed-
density gradient. It is not projected through the integral-symmetry petite
list, because a state-specific or excited CASSCF root need not furnish a
totally symmetric contracted density.

## Verification and optimization reference

The analytic derivative is tested against five-point central differences of
the SC-NEVPT2 total energy. For LiH/STO-3G with CAS(2,2), a frozen Li 1s core,
and exactly degenerate virtual pi orbitals, the largest analytic-minus-finite-
difference component is `3.4e-11` Hartree/Bohr.

The LiH optimization converges to
`E(SC-NEVPT2) = -7.8819849582` Hartree and
`r(Li-H) = 2.9270761417` Bohr. Repeating the calculation reproduces the bond
length to `3e-13` Bohr; starting from 1.75 Angstrom reaches the same minimum to
`2.5e-7` Bohr.

The former linear H4/STO-3G CAS(2,2) optimization is not a suitable coordinate
reference. On that surface H4 dissociates into two H2 fragments separated by
about 7.4 Bohr, and their relative translation is nearly free. Independent
runs can have the same energy to 12 decimal places while their Cartesian
coordinates differ by about 0.067 Bohr. That behavior is a property of the
nearly flat dissociation surface, not evidence of a nuclear-gradient defect.

## Cost and memory

One analytic evaluation replaces the `6 * natom` displaced SC-NEVPT2 energies
needed by a two-point Cartesian central difference. Its dense memory estimate
includes the orbital-space fourth-rank tensors, active-space RDM and adjoint
tensors, the native `nbf^4` density copy, and, for a harmonic basis, the
simultaneous `nbf_cart^2 * nbf^2` and `nbf_cart^4` Cartesian expansions. If
this estimate exceeds the configured memory ceiling, `gradient=auto` selects
the numerical route and `gradient=analytic` reports the required memory.
