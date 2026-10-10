# CASSCF Nuclear Gradient

OpenQP computes analytic nuclear gradients for both state-specific CASSCF and
state-averaged CASSCF (SA-CASSCF). State-specific `method=casscf`
differentiates the optimized `[casscf] root`. Dedicated `method=sa-casscf`
differentiates either the weighted objective or an individual averaged root,
as selected by `[casscf] gradient_state`.

The analytic path replaces the `6 * natom` displaced CASSCF calculations a
central difference needs with one CI solve at the converged orbitals plus one
pass of derivative integrals, so its cost stops growing with the number of
nuclei. It is also free of a failure mode the numerical route has: each
displaced CASSCF run is an independent chance to converge onto a different
solution branch, and nothing in a numerical gradient says that it did.

## Energy expression

With `x` a nuclear coordinate and the superscript `x` a derivative of the AO
integrals at fixed MO coefficients,

```text
dE/dx = sum D^AO h^x + 1/2 sum Gamma^AO (..|..)^x - sum X^AO S^x + dV_NN/dx
```

where `D^AO = C D C^T` is the spin-summed one-particle density, `X^AO = C F C^T`
is the energy-weighted density built from the generalized Fock matrix, and
`Gamma^AO` is the two-particle density.

There is **no orbital-response or CI-response term** for a state-specific
CASSCF energy or for the weighted SA-CASSCF objective. That is a consequence
of full variationality. For state-specific CASSCF it holds because a converged
solution is stationary with respect to every wavefunction parameter:

- the CI vector is an eigenvector of the active Hamiltonian at the current
  orbitals — true of an excited root as well, since a stationary point need not
  be a minimum;
- the non-redundant orbital-rotation gradient is what the optimizer drove to
  zero;
- rotations *within* the inactive, active and virtual blocks leave the energy
  invariant, the active block because the CI spans the full CAS.

This is a different object from the orbital-rotation gradient the optimizer
reports as `|g_orb|`. That gradient vanishing is the **precondition** for the
expression above, not an approximation of it. An individual state in an
SA-CASSCF calculation is not stationary with respect to the common optimized
orbitals; its effective densities therefore include a coupled orbital and CI
Z-vector response.

## Convergence matters more than usual

The gradient error at a non-stationary point is **first order** in `|g_orb|`, so
a loosely converged orbital optimization gives a plausible-looking but wrong
gradient. Converge tighter than you would for an energy:

```ini
[casscf]
gradient_norm_tol=1.0e-7
```

Every gradient run prints its own stationarity diagnostics. The driver refuses
to return a gradient when `|g_orb|` exceeds
`max(1e-4, 100 x [casscf] gradient_norm_tol)`, or when the generalized-Fock
asymmetry is nonfinite or exceeds `1e-6` Hartree:

```text
 Differentiated root                        0
 Orbital-rotation gradient norm        2.02955556E-08
 Generalized Fock asymmetry            6.17901000E-09
 Active 2-RDM correction vectors           10
```

`Generalized Fock asymmetry` is `max |F_pq - F_qp|`; the generalized Fock matrix
is symmetric at a stationary point, so this is an independent acceptance
condition for CI and active-active stationarity. It remains informative for a
full active space, where no non-redundant orbital rotations exist and
`|g_orb|` is identically zero. `Active 2-RDM correction vectors` is the rank of
the all-active two-particle correction (see
[Implementation notes](#implementation-notes)).

## Running it

### `.oqp`

```text
casscf/sto-3g grad guess(type=hcore) scf(conv=1e-10,maxit=80,forced_attempt=3,save_molden=false)
cas(active_electrons=4,active_orbitals=4,frozen_core=3,max_memory=512)
casscf(max_macro_iterations=50,gradient_norm_tol=1e-07) ci(solver=dense)
geom="../geometries/H2O-95fa99c614ed.xyz"
```

Runnable inputs:
[`examples/WF_methods/H2O_CASSCF_CAS44_grad.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/H2O_CASSCF_CAS44_grad.oqp)
and
[`examples/WF_methods/H4_CASSCF_CAS22_ROOT1_grad.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/H4_CASSCF_CAS22_ROOT1_grad.oqp),
each with a same-stem legacy `.inp` companion.

### Python

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_casscf_grad", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.casscf(active_electrons=4, active_orbitals=4, frozen_core=3,
           basis="sto-3g", runtype="grad", gradient_norm_tol=1.0e-7)

mol = job.run()
print("gradient:", mol.get_results()["gradient"])
```

`root=1` selects an excited root; the helper widens `[ci] nroot` itself so the
requested root is solved for. The returned gradient still occupies public
array slot `0`; `[properties] grad` is an array selector, not the physical CI
root number.

### Excited roots

`[casscf] root` selects the objective the orbitals are optimized for, and the
gradient follows it. On a symmetry-free H4 chain at STO-3G, the two roots give
genuinely different energies **and** genuinely different gradients, each of
which reproduces its own finite difference:

| `[casscf] root` | Energy (Ha) | `max abs(analytic - 5-point O(h^4) FD)` |
| ---: | ---: | ---: |
| `0` | `-2.1174290487` | `6.9e-10` |
| `1` | `-1.7712058509` | `4.4e-09` |

## SA-CASSCF analytic gradients

An SA-CASSCF calculation optimizes
`L = sum_I w_I E_I`. Two different derivatives can therefore be requested:

| `[casscf] gradient_state` | Differentiated quantity | Response |
| --- | --- | --- |
| `averaged` (default) | `dL/dx`, the weighted objective | None; the averaged density matrices are variational. |
| A root in `[state_average] target_roots` | `dE_J/dx`, one averaged root | Coupled orbital and CI Z-vector. |

For an individual root, the CI vector remains an eigenvector and is stationary
with respect to CI variation, but its energy is not stationary with respect to
the common state-averaged orbitals. OpenQP solves

```text
H_SA z = -g_J
```

where `H_SA` is the analytic, CI-relaxed Hessian of the weighted objective and
`g_J` is the orbital gradient of state `J`. This is the Schur complement of the
coupled orbital and CI response equations. The resulting orbital and CI
multipliers are incorporated into relaxed one- and two-particle density
matrices before the derivative-integral contraction.

The selector and state-array entry must agree. For example:

```ini
[input]
method=sa-casscf
runtype=grad

[casscf]
gradient_state=1
gradient_norm_tol=1.0e-7
zvector_tol=1.0e-8
zvector_degeneracy_tol=1.0e-8

[properties]
grad=1

[state_average]
enabled=true
nstate=2
equal_weights=true
target_roots=0,1
```

Use `gradient_state=averaged` and `grad=0` for the weighted-objective
derivative. Runnable examples are
[`LiH_SA-CASSCF_ANALYTIC_grad.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/LiH_SA-CASSCF_ANALYTIC_grad.oqp)
and
[`LiH_SA-CASSCF_ROOT1_grad.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/LiH_SA-CASSCF_ROOT1_grad.oqp).
The older
[`LiH_SA-CASSCF_grad.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/LiH_SA-CASSCF_grad.oqp)
uses `method=casscf` with state averaging enabled and remains a numerical
compatibility example.

## Geometry optimization

Because `runtype=grad` is analytic, gradient-driven optimizers run on CASSCF
without a numerical-gradient inner loop. Use `runtype=optimize` with the same
`[cas]`, `[ci]` and `[casscf]` sections. Keep `gradient_norm_tol` tight for
every point, not only the first. State-specific CASSCF publishes one gradient
row, so the optimizer selector must be zero even for an excited physical root:

```ini
[casscf]
root=1

[optimize]
istate=0
```

For `grad`, `optimize`, `ts`, `mep`, and `irc`, public energy slot `0` contains
the same selected physical root as public gradient slot `0`. The complete
physical-root energy list remains available in `OQP::CASSCF_ENERGIES` for
diagnostics and analysis.

A runnable ground-state example is
[`examples/WF_methods/LiH_CASSCF_optimize.inp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/LiH_CASSCF_optimize.inp).

An individual SA-CASSCF root also supports `optimize`, `ts`, `mep`, and `irc`
without a numerical-gradient inner loop. The target roots must be the
contiguous sequence `0,1,...`, and `[casscf] gradient_state` must equal
`[optimize] istate`, so the optimizer receives the energy and gradient of the
same physical root. The weighted-objective derivative is currently restricted
to `runtype=grad`: the optimizer interface reports state energies and has no
energy entry for `sum_I w_I E_I`.

## Scope and what is refused

Implemented on a closed-shell RHF reference with the contiguous
inactive/active/virtual partition used by the native CASSCF driver:

- the state-specific gradient of `[casscf] root` for `method=casscf`;
- the weighted-objective gradient for dedicated `method=sa-casscf`;
- the coupled orbital and CI Z-vector gradient of an averaged root for
  dedicated `method=sa-casscf`.

The analytic path refuses the following cases rather than silently applying an
invalid formula:

- **UHF or ROHF references.** The derivation and derivative-integral
  contraction currently require a closed-shell RHF reference.
- **An individual SA root at a degeneracy.** The adiabatic state energy is not
  differentiable at a crossing. A gap smaller than
  `zvector_degeneracy_tol` is refused.
- **An undefined or inaccurate Z-vector.** The individual-root path refuses a
  right-hand side with significant numerical-null-space leakage, an excessive
  residual, or an asymmetric relaxed generalized Fock matrix.
- **A nonstationary SA solution.** Both SA derivatives require the weighted
  orbital gradient to satisfy the acceptance threshold.
- **Weighted-objective geometry optimization.** The objective is not one of
  the state-indexed energies reported to the optimizer; use `runtype=grad`.
- **Unequal weights or a noncontiguous averaged root set in molecular input.**
  Root tracking between orbital macroiterations currently uses energy order,
  so a root interchange could change the physical objective. The response
  equations themselves support unequal weights, but the molecular calculation
  is refused until overlap tracking is implemented.
- **`[properties] grad` other than `0`.** State-specific CASSCF publishes one
  array row. Select the physical state with `[casscf] root`, while keeping the
  gradient-array selector at zero.
- **A non-Hartree-Fock Hamiltonian.** Leave `[input] functional` empty.
- **A non-contiguous active space.**
- **`method=casci`.** Its orbitals are not optimized, so it needs the
  orbital-response term this derivation drops.

Full active spaces with no non-redundant orbital rotations are supported: their
orbital-gradient norm is exactly zero by construction, while the independent
generalized-Fock-asymmetry condition still verifies CI and active-active
stationarity. An excited-root Davidson calculation with an explicit subspace
is also supported even when its energy optimization used the Python fallback;
gradient eligibility is assessed independently of native energy-driver
eligibility.

For `optimize`, `ts`, `mep`, and `irc`, use `[cas] orbital_source=rhf`.
JSON-sourced coefficients are fixed in the AO basis at one geometry and are no
longer orthonormal in the AO metric after the nuclei move.

## Implementation notes

`Gamma^AO` splits into a separable part plus an all-active correction. The
separable part, symmetrized over `mu <-> nu`, is exactly `4x` the Hartree-Fock
two-particle density expression evaluated at the CASSCF AO density, so it reuses
that kernel term for term.

The correction is not separable, and a `nbf^4` AO tensor is neither necessary
nor affordable. Symmetrized over the eight-element permutation group of
`(tu|vw)` it is a symmetric `nact^2 x nact^2` matrix, so its eigendecomposition

```text
dGamma_{mu nu la si} = sum_k lambda_k A^k_{mu nu} A^k_{la si},   A^k = C_act V^k C_act^T
```

turns it into at most `nact(nact+1)/2` separable products of ordinary symmetric
AO matrices — the rank printed as `Active 2-RDM correction vectors`. Storage is
`O(nact(nact+1)/2 x nbf_cart^2)` instead of `nbf_cart^4`.

The same factorized AO contraction is used for the weighted SA objective. The
individual-state response additionally builds the analytic SA orbital Hessian,
the Z-vector, and relaxed orbital and CI density contributions. It requires the
MO electron-repulsion tensor and several `nbf^4` contractions for the effective
generalized Fock matrix, so its memory requirement is `O(nbf^4)`, matching the
validation-grade scope of `[casscf] hessian=analytic`. Large active spaces are
refused by the existing determinant and excitation-stack memory limits rather
than approximated.

The symmetry petite list is deliberately **not** used. It is valid only for a
totally symmetric contracted density, which the state-specific two-particle
density of an arbitrary root does not guarantee — a spatially degenerate root in
an abelian subgroup is the counterexample. The cost of declining is speed; the
cost of opting in wrongly would be a wrong gradient.

### Cost

Measured for state-specific CASSCF on an M3 Ultra with the default thread
count, the analytic gradient
adds 3-10% to the wall time of the corresponding `runtype=energy` run and
1-21 MB to its peak resident set; the CASSCF orbital optimization dominates
both. Within the gradient itself, roughly 70-76% is the shared derivative-ERI
kernel and 20-23% is OpenMP load imbalance in the derivative-integral driver,
so the two-particle density assembly specific to CASSCF is a single-digit
percentage even at CAS(10,10).

The all-active correction is linear in its rank: at a fixed cc-pVDZ basis the
two-electron gradient step grows about 4.8 ms per correction vector, which is
1.3% of that step at CAS(4,4) and 6.8% at CAS(10,10).

The weighted SA objective has the same derivative-integral cost as the
state-specific gradient. An individual SA root also requires one analytic SA
orbital Hessian, one eigendecomposition in the non-redundant orbital-rotation
space, and two transition-density builds per averaged state.

## Verification

The implementation is checked against a **five-point O(h^4)** finite-difference
stencil at `h = 1e-3` bohr. The usual two-point O(h^2) difference is not
accurate enough here: along a bond-stretching coordinate the third energy
derivative is large enough that its truncation error alone reaches ~1e-5, which
is indistinguishable from a real defect.

| System | Differentiated quantity | `max abs(analytic - FD)` | `max abs(sum_A dE/dR_A)` |
| --- | --- | ---: | ---: |
| H2O/STO-3G CAS(4,4) | State-specific root 0 | `4.1e-09` | `2.3e-14` |
| H4/STO-3G CAS(2,2), symmetry-free | State-specific root 0 | `6.9e-10` | `1.4e-15` |
| H4/STO-3G CAS(2,2), symmetry-free | State-specific root 1 | `4.4e-09` | `5.0e-16` |
| H2O/cc-pVDZ CAS(6,6) | State-specific root 0 | `6.2e-08` | `1.2e-13` |
| H2O/cc-pVDZ CAS(2,2), spherical d shell | SA weighted objective | `<1.0e-06` | `<1.0e-08` |
| H2O/cc-pVDZ CAS(2,2), spherical d shell | SA root 1 Z-vector | `<1.0e-06` | `<1.0e-08` |

The state-specific H2O/cc-pVDZ CAS(6,6) row is converged to
`gradient_norm_tol=1e-6` rather than `1e-7`; since the gradient error is first
order in `|g_orb|`, `6e-8` is the corresponding accuracy.

The SA tests also verify the exact identity
`dL/dx = sum_I w_I dE_I/dx`, recovery of the independent state-specific
gradient in a one-state average, energy and gradient pairing in an actual
root-1 optimizer step, and fail-closed UHF and ROHF input validation. An
abstract integral model covers unequal weights, near-degenerate roots, and the
expected five-point convergence before molecular input is exercised.

Translational invariance, `max abs(sum_A dE/dR_A)`, is the sharpest cheap check
available on the Pulay and two-particle terms — the Hellmann-Feynman part alone
would not satisfy it.

!!! warning "A numerical CASSCF reference can silently mix solution branches"

    When cross-checking against a finite difference yourself, sample the energy
    along the coordinate and confirm it is smooth before treating a disagreement
    as a gradient defect. Three separate geometries produced finite-difference
    "failures" during development that were not gradient errors at all: a
    displaced point converging 0.024 Ha higher, a symmetric geometry sitting on
    a branch 8.3e-6 Ha above the one the displaced points fall onto, and an
    out-of-plane displacement that reordered nearly degenerate reference
    orbitals across the `[cas] sort_orbitals=energy` active-space boundary and
    so selected a different active space.

## See also

- [Wavefunction methods](../keywords/wavefunction.md) — `[cas]`, `[ci]`,
  `[casscf]` and `[state_average]` keywords.
- [Optimization](optimization.md) — gradient-driven geometry optimization.
- [Hessian and Frequencies](hessian.md) — second derivatives.
