# SA-CASSCF nuclear gradients: weighted objective and individual states

Two different derivatives live under the name "SA-CASSCF gradient", and they
are not approximations of one another.

| Quantity | Stationary? | Response needed | This page |
| --- | --- | --- | --- |
| `dL/dx`, `L = sum_I w_I E_I` | yes, by construction | none | Path A |
| `dE_J/dx`, one averaged root | **no** | coupled orbital + CI Z-vector | Path B |

`docs/casscf_analytic_gradient.md` covers the third, older object: the
*state-specific* gradient of a run that optimized one root's energy. That
formula must not be reused for either quantity here, and the code refuses
rather than doing so.

Throughout: chemist notation `(pq|rs)`, spin-summed densities

```
D_pq       = sum_sigma <a+_{p sigma} a_{q sigma}>
Gamma_pqrs = sum_{sigma tau} <a+_{p sigma} a+_{r tau} a_{s tau} a_{q sigma}>
E          = sum_pq D_pq h_pq + 1/2 sum_pqrs Gamma_pqrs (pq|rs) + V_NN
```

with inactive `i,j`, active `t,u,v,w`, virtual `a,b`, and `nc`/`na` the
inactive/active counts of the contiguous partition the native driver runs on.
These are the conventions of `casscf.py::_full_rdms`, `casscf_kernel.F90` and
`casscf_hessian.py`, reused rather than re-derived.

## 1. Parameters, and what the SA solution actually satisfies

Orbitals are rotated `C(kappa) = C exp(K)`, `K = sum_k kappa_k (e_{p_k q_k} -
e_{q_k p_k})` over the non-redundant pairs of `casscf._nonredundant_pairs`
(active-inactive, virtual-inactive, virtual-active). CI vectors are the
eigenvectors of the core-folded active Hamiltonian `Heff(kappa)` in the
determinant basis; write its complete spectrum as `(eps_M, |M>)` and let the
averaged roots be the subset `{|I>}` with weights `w_I`, `sum_I w_I = 1`.

CI variations are parameterized as an orthogonal rotation of that eigenbasis,
`|K(u)> = sum_M |M> (exp u)_{MK}`, `u` antisymmetric.

The converged SA-CASSCF solution satisfies exactly two families of conditions:

```
(O)   g^SA_k(kappa, u; x) = d/dkappa_k [ sum_I w_I E_I ] = 0     (npar of them)
(C)   r_{MK}(kappa, u; x) = <M| Heff |K> = 0,  M != K, K averaged
```

**(C) is not implied by (O).** With equal weights the SA objective only fixes
the *span* of the averaged CI vectors; which basis of that span is used is
fixed by (C), not by stationarity. Writing the Lagrangian with (C) as an
explicit constraint (rather than reusing `dL_SA/du = 0`) is what makes the
equal-weight case come out right; both routes give the same orbital equation,
and Section 4 shows where they differ and why the implemented choice is the
numerically safer one.

## 2. Path A -- gradient of the weighted objective

`L = sum_I w_I E_I` is stationary with respect to every wavefunction parameter:
(O) is its orbital stationarity, and (C) plus the eigenvalue conditions give
its CI stationarity. Redundant rotations (inactive-inactive,
virtual-virtual, active-active) leave it invariant. Hence no response term,
and the state-specific expression applies verbatim to the *averaged* densities:

```
dL/dx = sum_{mu nu} D^SA,AO_{mu nu} h^x_{mu nu}
      + 1/2 sum Gamma^SA,AO (mu nu|la si)^x
      - sum_{mu nu} X^SA,AO_{mu nu} S^x_{mu nu}
      + dV_NN/dx
```

with `D^SA = sum_I w_I D^I`, `Gamma^SA = sum_I w_I Gamma^I`, `X^SA = C F^SA C^T`
and `F^SA` the weight-averaged generalized Fock the energy path already builds
(`casscf_gfock_grad_w`).

One structural point makes this cheap. `Gamma^sep(D) = D_pq D_rs - 1/2 D_ps
D_rq` is quadratic in `D`, so `sum_I w_I Gamma^sep(D^I) != Gamma^sep(D^SA)` in
general -- but the only block in which two root-dependent factors meet is the
all-active one, and that block is overwritten by the averaged active 2-RDM
anyway. Every surviving block is affine in `gamma^I`. So

```
Gamma^SA = Gamma^sep(D^SA) + dGamma^SA,
dGamma^SA_{tuvw} = sum_I w_I Gamma^I_{tuvw} - Gamma^sep(gamma^SA)_{tuvw}
```

exactly, and Path A is the existing state-specific machinery fed
weight-averaged RDMs -- same `nact(nact+1)/2` factorization of the active
correction, same `grd2` contraction, no `nbf^4` AO tensor.

Path A is the derivative a future *state-averaged geometry optimization*
should follow. It is not the gradient of any physical state. The current
optimizer API exposes only state-indexed energies, so it refuses Path A rather
than pairing this objective gradient with an individual-state energy.

## 3. Path B -- gradient of one averaged root

### 3.1 Lagrangian

```
L_J = E_J(theta; x)
    + sum_k        z_k    g^SA_k(theta; x)
    + sum_{M<K}    y_{MK} r_{MK}(theta; x)
```

`theta = (kappa, u)`. Requiring `dL_J/dtheta = 0` removes every implicit
derivative, after which `dE_J/dx = dL_J/dx` at fixed `theta`.

### 3.2 The blocks

At the solution `Heff` is diagonal in `{|M>}`, which makes every second
derivative below elementary:

```
dE_J/dkappa_k = g^J_k                                (state-J orbital gradient)
dE_J/du       = 0                                    (c_J is an eigenvector)
dr_{MK}/dkappa_k = A_k^{MK} = <M| dHeff/dkappa_k |K>
dr_{MK}/du_{PQ}  = delta_{(PQ),(MK)} (eps_M - eps_K)
dg^SA_k/du_{PQ}  = 2 (w_Q - w_P) A_k^{PQ}
dg^SA_k/dkappa_l = H^oo_{kl}                         (fixed-CI orbital Hessian)
```

`dE_J/du = 0` is worth stating plainly: an individual averaged state *is*
stationary against CI variation (it is an eigenvector), and *is not*
stationary against orbital rotation. Only the second failure has to be
repaired -- but repairing it drags the CI in, because `g^SA` depends on the CI
vectors.

### 3.3 The CI multipliers eliminate in closed form

The `u_{PQ}` row of `dL_J/dtheta = 0` reads

```
2 (w_Q - w_P) (A^{PQ} . z) + y_{PQ} (eps_P - eps_Q) = 0
  =>  y_{PQ} = 2 (w_Q - w_P) (A^{PQ} . z) / (eps_Q - eps_P)
```

and `y_{PQ} = 0` whenever `w_P = w_Q` (both external, or an equal-weight pair
of averaged roots): the constraint is then decoupled from the orbital
multipliers. A vanishing denominator `eps_P = eps_Q` with `w_P != w_Q` is a
genuine singularity and is **refused**, not regularized -- it is the same
condition `casscf_hessian.py` already refuses when building the SA Hessian.

### 3.4 The orbital equation is the CI-relaxed SA Hessian

Substituting `y` into the `kappa` row,

```
[ H^oo + sum_{M<K} 2 (w_K - w_M) A^{MK} (A^{MK})^T / (eps_K - eps_M) ] z = -g^J
```

The bracket is exactly the matrix `casscf_hessian.analytic_orbital_hessian`
(and its Fortran twin `casscf_anhess.F90`) already assembles: `H^oo` is its
"part 1" fixed-CI Hessian, and the sum is its "part 2" CI relaxation, whose
per-ordered-visit coefficient `0.5 (w_I - w_m)` for averaged pairs and `w_I`
for external partners reproduces `2 (w_K - w_M)/(eps_K - eps_M)` term for term.
So the Z-vector is solved against **the exact second derivative of the
state-averaged objective**, and no separate Hessian implementation exists:

```
    H~ z = -g^J,     H~ = the SA orbital Hessian with CI relaxation folded in
```

`H~` is the Schur complement of the full `(orbital, CI)` Hessian. Solving in
the `npar`-dimensional orbital space instead of `npar + nstate*ndet` is not an
approximation.

### 3.5 Conditioning and the null space

`H~` is symmetric and generally indefinite (a converged SA solution need not
be a minimum of the SA objective in every direction). It is solved through a
symmetric eigendecomposition,

```
z = sum_{|lambda_i| > tau} v_i (v_i . rhs) / lambda_i,
tau = max([casscf] zvector_tol * max_i |lambda_i|, 1e-12)
```

with two fail-closed checks:

* if a **discarded** direction carries a non-negligible RHS component
  `|v_i . g^J|`, the individual-state gradient is not defined: the SA objective
  is flat along `v_i` while `E_J` is not, so `E_J` depends on which member of a
  continuum of equivalent SA solutions the optimizer happened to stop at.
  Refused.
* the residual `|H~ z + g^J|` is computed and reported, and refused above a
  threshold.

Directions that are exact redundancies of the SA objective are harmless: along
them `dL_SA/dtheta` vanishes identically, so `d2 L_SA/dtheta dx` vanishes too
and an arbitrary null component of `z` cannot change the gradient. The check
above is what separates those from accidental near-null directions.

### 3.6 CI relaxation vectors

Reducing the multiplier term back to determinant space (the derivation is in
Section 4) gives one vector per averaged root:

```
zeta_I = sum_{M != I, w_M != w_I} |M> (A^{MI} . z) / (E_I - eps_M)
T_CI   = 2 sum_I w_I <zeta_I| Heff |c_I>
```

`A^{MI}_k = <M| dHeff/dkappa_k |c_I>` is precisely the `amp[k, M]` array
`casscf_hess_amp` already produces for the Hessian, so `zeta_I` costs one
`(npar, ndet)` matrix-vector product per averaged root on top of the Hessian
build. `<c_I|zeta_I> = 0` by construction (and is re-projected numerically).

### 3.7 Effective densities

`dE_J/dx = dL_J/dx` at fixed `theta` is again a Hellmann-Feynman-plus-Pulay
expression, now over *relaxed* densities. Three contributions:

**(i) the state itself.** `D^J`, `Gamma^J` -- the ordinary full-space RDMs of
root `J`.

**(ii) orbital relaxation.** `sum_k z_k g^SA_k` is the directional derivative
of the fixed-density energy functional along `K(z)`, and such a derivative is
the same functional evaluated on *one-index transformed* densities:

```
D^z      = [K(z), D^SA]
Gamma^z  = one-index transform of Gamma^SA by K(z)
         = Gamma^cross(D^z, D^SA) + oit_z(dGamma^SA)
Gamma^cross(A,B)_pqrs = A_pq B_rs + B_pq A_rs - 1/2 (A_ps B_rq + B_ps A_rq)
```

using that the one-index transform is a derivation, so it acts on
`Gamma^sep(D) = D (x) D` by the product rule. For the low-rank active
correction `dGamma^SA = sum_k lambda_k M^k (x) M^k` (with `M^k` the active
eigenvector of the symmetrized correction embedded in the full space) the same
product rule gives `sum_k lambda_k (N^k (x) M^k + M^k (x) N^k)`, `N^k = [K(z),
M^k]`.

**(iii) CI relaxation.** `dT_CI/dx` needs the transition densities between
`zeta_I` and `c_I`:

```
gamma~_tu     = sum_I 2 w_I * sym <zeta_I| E_tu |c_I>
Gamma~_tuvw   = sum_I 2 w_I * sym <zeta_I| e_tuvw |c_I>
```

Because `<zeta_I|c_I> = 0` the core-energy term drops, so these map into the
full orbital space *without* the inactive occupancy:

```
D^ci      = gamma~ on the active block, zero elsewhere
Gamma^ci  = Gamma^cross(D^ci, D^inact) + Gamma~ on the all-active block,
D^inact   = 2 on the inactive diagonal, zero elsewhere
```

which is exactly the inactive-Fock folding `f_tu = h_tu + sum_i [2 (tu|ii) -
(ti|iu)]` written as a two-particle density. Symmetrized transition RDMs are
obtained by polarization of the ordinary RDM builder,
`sym T(a,b) = 1/2 [ T(a+b) - T(a) - T(b) ]`, so no transition-RDM kernel is
needed and the validated `rdm1_spatial` / `rdm2_spatial` engines are reused.

**Totals.**

```
D^eff     = D^J + D^z + D^ci
Gamma^eff = Gamma^J + Gamma^z + Gamma^ci
F^eff_pq  = sum_r D^eff_pr h_qr + sum_rst Gamma^eff_prst (qr|st)
```

### 3.8 The Pulay term, and why `F^eff` must come out symmetric

Orbital orthonormality `C^T S C = 1` fixes the symmetric part of the MO-basis
coefficient response, `U + U^T = -C^T S^x C`; the antisymmetric part is a free
orbital rotation and is annihilated by the stationarity of `L_J`. Since
`dE~/dU_mp = 2 F^eff_pm`, the surviving piece is

```
- sum_{mu nu} X^eff_{mu nu} S^x_{mu nu},    X^eff = C * sym(F^eff) * C^T
```

exactly as in the state-specific case, and the antisymmetric part of `F^eff`
must vanish. It does, in all three redundant blocks, but for three different
reasons:

* **inactive-inactive and virtual-virtual** -- the SA energy functional is
  identically invariant there, so moving along such a rotation maps one SA
  solution onto another; the whole non-redundant gradient stays zero, hence
  `H^oo_{(pq),(ij)} = 0` and `dE_J/dkappa_ij = 0` separately.
* **active-active** -- `E_J` at *fixed CI coefficients* is not invariant, and
  neither is `z . g^SA`. What saves it is that an active-active orbital
  rotation is equivalent to a CI transformation, so its effect on `L_J` is
  `-dL_J/du` along the induced CI direction, and that is zero by the `u` row of
  the stationarity conditions. This is a real dependence of the Pulay term on
  having solved the CI multipliers: drop them and the Pulay term is wrong,
  not merely incomplete.

`max |F^eff_pq - F^eff_qp|` is therefore a live correctness probe, not
decoration, and is reported with every Path B gradient.

## 4. Why the CI constraint is written as (C) and not as `dL_SA/du = 0`

Using `dL_SA/du = 0` as the CI constraint gives multipliers `z^u_{MK}` related
to the `y` above by `y_{MK} = 2 (w_K - w_M) z^u_{MK}`, and an identical orbital
equation. The two formulations differ only in what they do with an
equal-weight pair of averaged roots `I, K` (`w_I = w_K`, `eps_I != eps_K`):

* `y_{IK} = 0` -- the pair contributes nothing;
* `z^u` keeps two terms, one in `zeta_I` and one in `zeta_K`, which cancel
  exactly in `dT_CI/dx` because `<c_K| Heff^x |c_I>` is symmetric.

They agree, but the second form evaluates a cancelling pair of terms with
denominators `+-(E_I - E_K)`, which is a catastrophic cancellation as soon as
two equal-weight averaged roots approach each other -- precisely the regime
state averaging exists for. The implementation therefore skips them, which is
also literally the `coef == 0` skip `casscf_hess_relax` performs when building
`H~`, so one rule covers both.

## 5. What is refused

Path B is refused, with the reason named, when:

* the requested state is not one of the averaged roots;
* two roots with different weights are degenerate to within
  `zvector_degeneracy_tol` and have non-vanishing orbital coupling (the SA
  objective is not smooth there and `H~` does not exist) -- inherited from the
  Hessian builder;
* the requested root is degenerate with any other CI root to within the same
  tolerance: at an exact crossing the individual adiabatic energy is not
  differentiable, whatever the weights. The gap is reported so a
  near-degenerate but resolved case can be told apart from a real crossing;
* the Z-vector RHS has a component in the numerical null space of `H~`
  (Section 3.5), or the solved residual is above tolerance;
* `|F^eff - F^eff,T|` is above tolerance (Section 3.8);
* `|g^SA|` at the reported orbitals is above the acceptance limit -- the whole
  construction assumes (O) holds;
* the Hamiltonian is not Hartree-Fock, the active space is not contiguous, or
  the excitation stack exceeds the dense-spectrum memory guard.

## 6. How it is verified

Verification is layered, cheapest and sharpest first, in
`tests/test_casscf_sa_gradient.py`.

**An abstract model, before any molecule.** A one-parameter family of integrals
`(h(x), g(x))` in a fixed orthonormal basis defines a complete SA-CASSCF
problem with no basis set, no derivative integrals and no Pulay term. The
orbitals are optimized to stationarity at each `x` and the root energies are
differentiated by a five-point O(h^4) difference of the *re-converged*
solution, against

```
dE_J/dx = sum D^eff (dh/dx)^C + 1/2 sum Gamma^eff (dg/dx)^C.
```

A failure there is a defect in the Lagrangian/Z-vector algebra and in nothing
else, which is why it runs first. Observed agreement, 1e-10 to 3e-10 (the
finite-difference roundoff floor) over: equal and unequal weights, two and
three averaged states, with and without inactive orbitals, non-contiguous
averaged roots, and a larger active space. Step-size refinement at
h = 8e-3 .. 1e-3 shows the expected O(h^4) approach to the analytic value until
roundoff takes over at ~1e-10 -- so the residual is truncation error, not a
constant offset.

**Two identities that carry no truncation error at all.**

* `dL/dx = sum_I w_I dE_I/dx` holds to 1e-15..1e-16. It is exact because the
  weights are constants, so it fails sharply if either path is wrong, and it
  ties Path A and Path B together without a finite difference anywhere. On a
  molecule it additionally covers the AO transform, the factorized
  two-particle density and the derivative-integral contraction.
* `max |F^eff - F^eff,T|` comes out at ~1e-14.

**Negative controls.** Removing the response must break the answer, or the
tests cannot distinguish "applied" from "computed and discarded":

| what is removed | effect |
| --- | --- |
| the whole Z-vector | gradient error 8e-3 .. 9e-2, i.e. 7-9 orders above the agreement |
| the CI multipliers only | `max |F^eff - F^eff,T|` jumps from ~1e-14 to 4.9e-2 and the run is refused |

The second row is the quantitative form of Section 3.8: the active-active block
of `F^eff` is annihilated by CI stationarity and by nothing else, so the
symmetry check detects that specific omission on its own, with no reference
value needed.

**Near-degeneracy.** With unequal weights the CI response denominators *are*
the root gap. Searching the model for an avoided crossing and differentiating
on the way in: at gaps of 1.2e-1, 1.4e-2 and 1.7e-3 Eh the agreement is
3e-10, 1.4e-9 and 2.9e-9 (with the finite-difference step shrunk to match),
while the sum rule stays at 1e-15. At the crossing itself the run is refused.

**A one-state "average" must reproduce the state-specific gradient.** With one
weight of 1.0 both paths must return exactly what `casscf_gradient.F90`
computes through its own independent AO density and contraction. Path B is the
stronger of the two: it runs the whole Z-vector machinery, which must then find
`g^J = g^SA = 0`, solve `z = 0` and reduce to the same number.

**Molecules**, against a five-point O(h^4) difference of the SA-CASSCF energy,
plus translational invariance -- the sharpest cheap check on the Pulay and
two-particle terms, which the Hellmann-Feynman part alone would satisfy
anyway. Each finite difference also asserts that the sampled energies are
smooth: a state average has many stationary points and a displaced run can land
on a different one, which is not a gradient defect and must not be reported as
one.

On LiH/STO-3G CAS(2,2) averaged over the two lowest singlets, the run reports

```
SA orbital gradient |g_SA|      1.968e-09     converged
state orbital gradient |g_J|    2.063e-01     root 1 is NOT stationary
Z-vector norm                   1.250e-01     Z-vector residual  2.2e-16
Hessian null directions         0             null-space leakage 0.0
CI relaxation |zeta|            5.216e-03
effective Fock asymmetry        4.980e-11
2-particle factorization check  4.441e-16
```

`|g_J| = 0.21` against `|g_SA| = 2e-9` is the whole justification for this page
in one pair of numbers: the averaged objective is stationary and the individual
root is not, by eight orders of magnitude. The two examples differentiate the
same geometry to `-9.78e-4` (objective) and `+1.54e-2` (root 1).

Unequal weights are implemented and are covered by the abstract protocol, but a
molecular input cannot currently request them: the SA-CASSCF preflight refuses
unequal weights because the averaged roots are followed by energy order and
overlap tracking is not implemented, so a root flip between macroiterations
would move the weights onto different physical states. That is a pre-existing
constraint on the energy path, not a limitation of the response.

Water/cc-pVDZ supplies a spherical d shell and therefore exercises the
spherical-basis derivative transformation absent from the STO-3G cases. Both
Path A and the root-1 Path B gradient
agree with five-point energy differences within `1e-6` Eh/bohr and satisfy
translational invariance within `1e-8` Eh/bohr. An optimizer integration test
also checks that the root-1 energy and gradient are selected from matching
state-array entries. UHF and ROHF requests are rejected by input validation
because this derivation and its AO contraction are restricted to a closed-shell
RHF reference.

## 7. Cost

Path A is the cost of the state-specific gradient.

Path B adds one analytic SA orbital Hessian (`O(ndet^3)` for the dense active
spectrum, `O(nact^2 ndet^2)` for the excitation stack), one `npar^3`
eigendecomposition, `2 nstate` extra RDM builds, and a handful of `nbf^4`
contractions for `F^eff`. It is `O(nbf^4)` in memory because the MO ERI tensor
is needed for the Hessian -- the same validation-grade scope as
`[casscf] hessian = analytic`, and the reason large active spaces are declined
rather than approximated.

## 8. Input selection

```
[casscf]
  gradient_state = averaged        # Path A (default), the weighted objective
  gradient_state = 1               # Path B, the averaged root with CI index 1
  zvector_tol = 1.0e-8             # null-space / residual conditioning
  zvector_degeneracy_tol = 1.0e-8  # gap below which a crossing is refused
```

Anything else is an error, including a `gradient_state` naming a root that is
not in `[state_average] target_roots`. The selector is deliberately NOT
`[properties] grad`: that list names states, and the weighted objective is not
a state. A `[properties] grad` entry that disagrees with `gradient_state` is
rejected rather than silently redirecting the calculation.

From the Python API:

```python
job.sa_casscf(active_electrons=2, active_orbitals=2, frozen_core=1,
              nstate=2, gradient_state=1, runtype="grad")
```

Examples: `examples/WF_methods/LiH_SA-CASSCF_ANALYTIC_grad.inp` (Path A) and
`examples/WF_methods/LiH_SA-CASSCF_ROOT1_grad.inp` (Path B). Implementation:
`pyoqp/oqp/library/casscf_sa_gradient.py` and the derivative-integral
contraction `source/modules/casscf_ao_gradient.F90`.
