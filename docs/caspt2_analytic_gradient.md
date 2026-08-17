# Analytic CASPT2 / XMS-CASPT2 nuclear gradient

`runtype=grad` (and every gradient-driven runtype: `optimize`, `ts`, `mep`,
`irc`) computes the CASPT2-family nuclear gradient analytically for the variants
listed under [Scope](#scope).  One PT2 evaluation plus one pass of derivative
integrals replaces the `6 * natom` displaced PT2 energies a central difference
needs, and the cost no longer grows with the number of nuclei.

`[pt2] gradient` selects the route:

| value | behaviour |
| --- | --- |
| `auto` (default) | analytic where the variant has one, central differences otherwise |
| `analytic` | refuse rather than fall back |
| `numerical` | always central differences |

## Where the work happens

| Piece | File | Role |
| --- | --- | --- |
| Theory | `pyoqp/oqp/library/caspt2_gradient.py` | rebuilds the CASPT2 state, solves the amplitude, CI, XMS-rotation and orbital-response equations, reduces the Lagrangian to three AO objects |
| Integrals | `source/modules/caspt2_gradient.F90` | contracts those three objects with the derivative one- and two-electron integrals |

The split follows the CASPT2 energy path itself, where the driver is Python and
liboqp computes.  Derivative integrals exist only in Fortran, so that is the
part that has to be Fortran; everything else is linear algebra on objects the
energy path already forms.

## Why the determinant formulation makes this tractable

OpenQP's CASPT2 is **uncontracted**: the first-order interacting space `Q` is the
set of external determinants, an orthonormal basis whose labels do not depend on
the geometry.  There is no internally-contracted overlap metric, no linear
dependence removal and no perturber normalization -- the three things that make
an internally-contracted CASPT2 gradient hard.  What is left is

```
E   = E_ref + E2 ,   E_ref = <Psi_0|H|Psi_0> + V_NN
E2  = V . T ,        V = (P_Q H |Psi_0>) ,   T = -G(D)^-1 V
D   = P_Q H0 P_Q - E0 ,                      E0 = <Psi_0|H0|Psi_0>
```

with `G` the (possibly regularized) denominator function.

## Zeroth-order Hamiltonian, in invariant form

The energy path builds `H0 = diag(eps)` with `eps_p = f_pp` and `f` the
closed+active Fock `h + J(D_sa) - 1/2 K(D_sa)`.  That representation is tied to
the semicanonical basis.  The gradient uses instead the **block-diagonal Fock
operator**

```
H0^inv = sum_{p,q in the same orbital block} f_pq E_pq
```

which is numerically IDENTICAL at the semicanonical orbitals -- the within-block
off-diagonals of `f` vanish there, and the module checks it and refuses
otherwise -- but is invariant under rotations inside the inactive, active and
virtual blocks.  Those rotations then need no constraint and no multiplier.

The IPEA shift biases the active *diagonal* of `H0` in a particular basis and so
breaks that invariance; it is refused rather than differentiated approximately.

## Amplitudes, denominators and shifts

`E2 = V.T` is not stationary in `T` once a level shift is applied, so the
amplitudes carry an explicit multiplier.  Solving the multiplier equation gives
`lambda = T` exactly for a single state -- the Lagrangian collapses to the
*shifted* Hylleraas functional -- and

```
lambda_I = -u_I G(D_I)^-1 (sum_J u_J V_J)
```

in the multistate case, where `u` is the effective-Hamiltonian mixing vector.

The derivative of `G(D)` is the derivative of a **matrix function**, given by the
Daleckii-Krein divided-difference formula

```
lambda . dG(D) . T = Tr[Omega dD] ,
Omega~_ij = G^[1](w_i, w_j) lambda~_i T~_j      (symmetrized)
```

with `G^[1](a,b) = (G(a)-G(b))/(a-b)`.  The option validator in
`caspt2_dyall._caspt2_options` admits only three regularization families (real
shift; real + imaginary; GAMESS ISA `edshft`), and all three have

```
G(w) = d + e/d ,  d = w + level_shift ,  e in {0, imaginary_shift^2, edshft}
```

so `G^[1](a,b) = 1 - e/(d_a d_b)` and `Omega` is exactly rank <= 2 (rank <= 4 in
the multistate case):

```
Omega = sym(lambda T^T) - e sym(lambda' T'^T) ,   X' = (D + s)^-1 X
```

No `n_ext x n_ext` object is ever formed, and the shift derivative is **exact**
rather than neglected.  This is what the `imaginary_shift` and `edshft` rows of
the test matrix measure: dropping the term would be a first-order error, not a
rounding one.

## Reference relaxation

The reference orbitals and CI vectors are not variational for the PT2 energy, so
both carry multipliers.

**CI.**  One projected linear solve per reference root,

```
(H_model - E_r) y_r = -(1 - c_r c_r^T) g_r ,     y_r . c_r = 0
```

Converting `y_r` into an extra transition density means its coupling back into
the orbital equation needs no separate bookkeeping: the orbital gradient is
computed from the augmented densities.

**Orbitals.**  The constraint is whatever actually determines the reference
orbitals as a function of the geometry:

| `[pt2] reference` | constraint |
| --- | --- |
| `casscf` | the (state-averaged) CASSCF orbital-rotation gradient vanishes on the non-redundant pairs |
| `casci` with `[cas] orbital_source=rhf` | the orbitals diagonalize the RHF Fock |

Both Jacobians are closed-form.  With `kappa_A` an antisymmetric generator and
the rotation-derivative rule `dM = [M, kappa]` applied to every tensor index,

```
dF^RHF/dt_A = [F^RHF, kappa_A] - G([D^RHF, kappa_A])
dF^gen/dt_A = [F^gen,  kappa_A] - F^gen(dD_A, dGamma_A)
```

the first of which reduces to the textbook CPHF `A` matrix.  **No finite
difference appears anywhere in the implementation.**

### Which rotations are parameters

The PT2 energy depends on the reference orbitals only through the **spans** of
the three orbital blocks -- including the frozen core, whose orbitals are the
lowest eigenvectors of the closed+active Fock restricted to the inactive span,
which a within-inactive rotation leaves unchanged.  Rotations inside a block are
therefore redundant *for the energy*.

They are not automatically redundant for the *constraint*, and the distinction
matters:

* `reference=casscf` -- a within-block rotation leaves the CASSCF orbital
  gradient at zero on the inter-block pairs, so it is a true gauge freedom and
  is excluded.  The parameter set is the non-redundant pairs.
* `reference=casci` -- the RHF occupied space is the first `nocc` reference
  orbitals, and the CASPT2 **active block straddles that boundary**.  An
  active-active rotation therefore moves an orbital across it, changes `D^RHF`,
  and breaks the constraint on the inter-block pairs as well.  Those rotations
  are parameters, constrained by RHF canonicality inside the active block.
  Within-inactive and within-virtual rotations lie wholly inside the RHF
  occupied and virtual spaces, leave `D^RHF` untouched, and are excluded -- their
  Jacobian columns are identically zero, which would make the system singular
  for a symmetry-degenerate pair such as the 2p shell of Li.

The whole orbital response is assembled and solved densely.  That is
proportionate: the module is validation-grade in the same sense the CASPT2
energy path is -- the external determinant space is enumerated in full -- so the
orbital-pair and CI dimensions are small by construction.

## XMS

For `xms-caspt2` / `xmcqdpt2` the references are pre-rotated by the eigenvectors
`R` of the state-averaged Fock in the model space.  `R` depends on the geometry;
its response is eliminated with first-order eigenvector perturbation theory,

```
sum_ki (dL/dR_ki) dR_ki = Tr[Z_R dM] ,
Z_R = sym(R Phi^T R^T) ,   Phi_ji = (R^T dL/dR)_ji / (lam_i - lam_j)
```

so the whole rotation response is one extra weight matrix on the Fock operator
plus one extra term in the reference-vector gradient.  A near-degenerate
model-space Fock makes `dR/dx` undefined and is refused.

The effective-Hamiltonian eigenvectors need **no** response -- Hellmann-Feynman
for an eigenvalue -- but the gap to the neighbouring root is checked before the
derivative is taken.

## The gradient itself

```
dE/dx = sum_{mu nu}      D^AO_{mu nu}    h^x_{mu nu}
      + 1/2 sum_{mu nu la si} Gamma^AO   (mu nu|la si)^x
      - sum_{mu nu}      X^AO_{mu nu}    S^x_{mu nu}
      + dV_NN/dx
```

`X` is the symmetric generalized Fock of the relaxed densities.  The
two-particle density is handed to the derivative-ERI driver factorized,

```
Gamma^AO_{mu nu la si} = sum_k lam_k A^k_{mu nu} A^k_{la si}
```

from the eigendecomposition of the eight-fold-symmetrized `Gamma` over the
composite index `(mu nu)`.  That is an exact rewriting -- a change of storage,
not a fit -- and every `A^k` is an ordinary symmetric AO matrix, so it takes the
same `bfnrm` folding and spherical-to-Cartesian expansion the density takes.
Unlike the CASSCF gradient, no Hartree-Fock-like separable part is split off:
the CASPT2 two-particle density is not separable in any orbital block, because
the first-order amplitudes correlate inactive, active and virtual space at once.

## Scope

| supported | not supported |
| --- | --- |
| `caspt2`, `mrmp2` (single state) | `ms-caspt2` -- the **multi-set** construction (per-state orbitals, per-state full-Fock-matrix H0, inter-state Loewdin-minor rotations); use `xms-caspt2` |
| `mcqdpt2` (single-set multistate) | `h0=dyall` (NEVPT2) |
| `xms-caspt2`, `xmcqdpt2` | `contraction=strong` (SC-NEVPT2) |
| `level_shift`, `imaginary_shift`, `edshft` -- exact | `ipea_shift != 0` |
| `reference=casci` (RHF orbitals), `reference=casscf` (state-specific or state-averaged) | `[cas] orbital_source` reading orbitals from a file |
| the PT2 frozen core (`[pt2] frozen`) | |

Beyond the variant scope, the gradient is refused rather than approximated when
the reference cannot support it: an unconverged CASSCF, orbitals that are not
semicanonical, a singular orbital-response system, a degenerate
effective-Hamiltonian root, or a reconstruction that does not reproduce the
reported energy.  Each writes the offending number into the log.

## Validation

Measured on an off-symmetry H4 (and LiH for the frozen core), against a
five-point stencil of independently computed total energies at three step sizes.
`max |analytic - 5-point|` in Hartree/Bohr, at `h = 1e-3` Bohr:

| case | residual |
| --- | --- |
| `caspt2`, CASCI reference | 6e-11 |
| + `level_shift=0.15` | 7e-11 |
| + `imaginary_shift=0.20` | 7e-11 |
| + `edshft=0.05` | 6e-11 |
| `caspt2`, CASSCF reference | 5e-8 |
| `mcqdpt2`, 2 roots | 6e-11 / 2e-10 |
| `xms-caspt2`, 2 roots | 5e-11 / 2e-10 |
| `caspt2`, 6-31G | 1e-9 |
| `caspt2`, LiH with the default frozen core | 1e-11 |

The CASSCF-reference row is looser because the finite-difference side inherits
the CASSCF convergence of every displaced point; its residual falls from 3.5e-6
at `h = 4e-3` to 5e-8 at `h = 1e-3`, which is the behaviour of the reference,
not of the quantity being tested.

Translational invariance `|sum_A g_A|` and rotational invariance
`|sum_A r_A x g_A|` hold to 1e-15 and 1e-13.  The numbers above were reproduced
independently on macOS/arm64 with Accelerate ILP64 and on Linux/x86-64 with
OpenBLAS ILP64, agreeing to the quoted digit.

`tests/test_caspt2_analytic_grad.py` carries the same checks in a form the CI
can run, together with the scope refusals and an internal consistency test that
the Lagrangian is stationary along every constrained orbital rotation.
