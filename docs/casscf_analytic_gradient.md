# Analytic state-specific CASSCF nuclear gradient

`runtype=grad` (and every gradient-driven runtype: `optimize`, `ts`, `mep`,
`irc`) with `method=casscf` computes the nuclear gradient analytically.
One CI solve at the converged orbitals plus one pass of derivative integrals
replaces the `6 * natom` full CASSCF runs a central difference needs, and the
cost no longer grows with the number of nuclei.

## Two different objects called "the CASSCF gradient"

OpenQP contains both, and they must not be confused.

| Quantity | Where | What it is |
| --- | --- | --- |
| `casscf_gfock_grad` | `source/modules/casscf_kernel.F90` | derivative with respect to the **orbital rotation** parameters `kappa`; the orbital optimizer drives it to zero |
| `casscf_gradient` | `source/modules/casscf_gradient.F90` | derivative with respect to the **nuclear coordinates**; what this page documents |

The first vanishing is a *precondition* for the second, not an approximation of
it.

## Conventions

Chemist notation, `(pq|rs)` grouping the charge distributions `pq` and `rs`.
Spin-summed density matrices,

```
D_pq       = sum_sigma <a+_{p sigma} a_{q sigma}>
Gamma_pqrs = sum_{sigma,tau} <a+_{p sigma} a+_{r tau} a_{s tau} a_{q sigma}>
```

so that

```
E = sum_pq D_pq h_pq  +  1/2 sum_pqrs Gamma_pqrs (pq|rs)  +  V_NN.
```

With inactive orbitals `i,j`, active `t,u,v,w` and virtual `a,b`, CASSCF has
`D` equal to 2 on the inactive diagonal and the CI 1-RDM `gamma_tu` on the
active block, and `Gamma` equal to the separable form

```
Gamma^sep_pqrs = D_pq D_rs - 1/2 D_ps D_rq
```

everywhere except the all-active block, which is the CI 2-RDM. These are
exactly the conventions of `casscf.py::_full_rdms` and of the generalized-Fock
engine `casscf_kernel.F90`, whose

```
F_pq = sum_r D_pr h_qr + sum_rst Gamma_prst (qr|st)
```

is reused rather than re-derived.

## The gradient

For a nuclear coordinate `x`, with a superscript `x` denoting a derivative of
the AO integrals at fixed MO coefficients,

```
dE/dx = sum_{mu nu} D^AO_{mu nu} h^x_{mu nu}
      + 1/2 sum_{mu nu la si} Gamma^AO_{mu nu la si} (mu nu|la si)^x
      - sum_{mu nu} X^AO_{mu nu} S^x_{mu nu}
      + dV_NN/dx
```

with `D^AO = C D C^T`, `X^AO = C F C^T` and `Gamma^AO` the four-index
transform of `Gamma`. There is **no orbital-response and no CI-response
(Z-vector) term**.

### Why no response term

Three stationarity conditions, all of which hold at a converged
*state-specific* CASSCF solution:

1. **CI stationarity.** The CI vector is an eigenvector of the active
   Hamiltonian at the current orbitals, so the energy is stationary with
   respect to every CI parameter. This holds for an excited root as well — a
   stationary point need not be a minimum, and the derivation only needs
   stationarity.
2. **Non-redundant orbital stationarity.** The inactive-active,
   inactive-virtual and active-virtual rotation gradient is what the optimizer
   converged to zero.
3. **Redundant orbital stationarity.** Rotations *within* the inactive, active
   and virtual blocks leave the energy invariant — for the active block because
   the CI space spans the full active CAS — so the energy is stationary there
   too, automatically.

Together these make `F` symmetric, and `X^AO = C F C^T` is then the ordinary
energy-weighted density. As a check on the convention: for a closed-shell RHF
reference `F_ij = 2 eps_i delta_ij`, and the Pulay term reduces to the familiar
`-2 sum_i eps_i S^x_ii` that `grd1::eijden` builds for HF.

Because orbital and CI stationarity are numerical conditions, the driver
reports `|g_orb|` and `max |F_pq - F_qp|` with every gradient. PyOQP refuses a
gradient when `|g_orb|` exceeds `max(1e-4, 100 x [casscf]
gradient_norm_tol)`, or when the generalized-Fock asymmetry is nonfinite or
exceeds `1e-6` Hartree. The second condition also applies to a full active
space, where the absence of non-redundant orbital rotations makes `|g_orb|`
identically zero but does not by itself establish CI stationarity.

## The two-particle density in the AO basis

The derivative-ERI driver `grd2` asks for a two-particle density block per
shell quartet and assumes it carries the full eight-fold permutational symmetry
of `(mu nu|la si)`. Splitting `Gamma^AO` into its separable part and an
all-active correction,

```
Gamma^AO      = Gamma^sep(D^AO) + dGamma
dGamma_{tuvw} = Gamma_tuvw - Gamma^sep_tuvw      (active MO indices only)
```

the separable part, symmetrized over `mu <-> nu`, is

```
Gamma^sep = D_{mu nu} D_{la si}
          - 1/4 (D_{mu la} D_{nu si} + D_{nu la} D_{mu si}),
```

which is `4x` the Hartree-Fock `get_density` expression at `hfscale = 1`. So
the CASSCF two-particle density is the Hartree-Fock one evaluated at the CASSCF
AO density, plus a correction confined to the active space.

That correction is not separable, and an `nbf^4` AO tensor is neither necessary
nor affordable. After symmetrization `dGamma` is a symmetric matrix over the
composite index `(tu)`, of dimension `nact^2`, so its eigendecomposition

```
M_{(tu),(vw)} = sum_k lambda_k V^k_{tu} V^k_{vw}
```

turns it into a short sum of separable products of ordinary AO matrices,

```
dGamma_{mu nu la si} = sum_k lambda_k A^k_{mu nu} A^k_{la si},
A^k = C_act V^k C_act^T,
```

with at most `nact(nact+1)/2` nonzero terms — the symmetrization annihilates
the antisymmetric half of the composite space. Each `A^k` is a symmetric AO
matrix, so it goes through the same `bfnrm` folding and spherical-to-Cartesian
expansion (`build_cart_density`) the density does, and each term is manifestly
eight-fold symmetric. Storage is `O(nact(nact+1)/2 x nbf_cart^2)` instead of
`nbf_cart^4` — two arrays of that size, `A^k` and the pre-scaled `lambda_k A^k`
the per-quartet contraction loop reads. The rank actually used is reported in
the log as `Active 2-RDM correction vectors`.

Against the `nbf^4` MO-integral buffer the CASSCF driver already holds, that
footprint is not observable: over the benchmark set below the whole gradient
adds 1–21 MB to the peak resident set of the corresponding energy run.

## Scope and what is refused

Implemented: the **state-specific** gradient of the root selected by `[casscf]
root`, for `method=casscf` with a closed-shell RHF reference and the contiguous
inactive/active/virtual partition the native CASSCF driver runs on.

Refused, with a message rather than a silently wrong number:

- **`method=sa-casscf` / `[state_average] enabled=true`.** A state-averaged run
  optimizes `sum_I w_I E_I`, and neither of its two derivatives is this
  formula: the *averaged objective* needs the weight-averaged RDMs, and an
  *individual averaged state* is not variational at all — only the weighted sum
  is stationary against orbital rotation, so it requires a coupled orbital+CI
  Z-vector response. Both live in `docs/sa_casscf_gradients.md` and
  `pyoqp/oqp/library/casscf_sa_gradient.py`, selected by `[casscf]
  gradient_state`; this entry point refuses rather than reusing its own
  expression on averaged quantities.
- **A non-Hartree-Fock Hamiltonian.** The energy expression above has no
  exchange-correlation term; CASSCF already requires `[input] functional` unset.
- **A non-stationary starting point**, according to either the `|g_orb|` or
  generalized-Fock-asymmetry condition above.
- **A `[properties] grad` selector other than `0`.** State-specific CASSCF
  publishes one array row; `[casscf] root` selects which physical root occupies
  that public slot. Gradient-driven calculations place both the selected-root
  energy and its gradient in slot 0; the complete physical-root energy list
  remains available as `OQP::CASSCF_ENERGIES`.

`method=casci` is not covered: CASCI orbitals are not optimized, so its
gradient needs the orbital-response term this derivation drops.

## Symmetry

The symmetry petite list is deliberately *not* used for the CASSCF two-electron
gradient. It is valid only when the contracted density is totally symmetric,
and the state-specific two-particle density of an arbitrary `[casscf] root`
carries no such guarantee — a spatially degenerate root in an abelian subgroup
is the obvious counterexample. The cost of declining is speed; the cost of
opting in wrongly is a wrong gradient.

## Examples

- `examples/WF_methods/H2O_CASSCF_CAS44_grad.inp` — H2O/STO-3G CAS(4,4),
  ground root.
- `examples/WF_methods/H4_CASSCF_CAS22_ROOT1_grad.inp` — linear H4/STO-3G
  CAS(2,2) with `[casscf] root=1`, pinning root selection.
- `examples/WF_methods/LiH_CASSCF_optimize.inp` — LiH/STO-3G CAS(2,2)
  geometry optimization using the analytic gradient with
  `[optimize] istate=0`.

## Python API

```python
from oqp.openqp import OpenQP

job = OpenQP()
job.system(...)
job.casscf(active_electrons=4, active_orbitals=4, frozen_core=3,
           runtype='grad')
job.run()
```

`job.casscf(root=1, runtype='grad')` differentiates the selected excited root.
