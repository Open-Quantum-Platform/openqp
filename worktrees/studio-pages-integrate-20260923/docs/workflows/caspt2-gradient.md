# CASPT2 Nuclear Gradient

OpenQP computes the CASPT2-family nuclear gradient analytically for the variants
listed under [Scope](#scope), and by central differences for the rest. Set
`runtype=grad` (or any gradient-driven run type: `optimize`, `ts`, `mep`, `irc`)
and select the route with `[pt2] gradient`.

| `[pt2] gradient` | Unsupported variant | Failed precondition at this geometry |
| --- | --- | --- |
| `auto` (default) | central differences, logged | central differences, logged |
| `analytic` | refuse | refuse |
| `numerical` | central differences | central differences |

The second column is explained under
[When the analytic route does not apply](#when-the-analytic-route-does-not-apply).

The analytic path is one PT2 evaluation plus one pass of derivative integrals,
so its cost stops growing with the number of nuclei; the central difference
needs `2 * 3 * natom` displaced PT2 energies.

## Why an uncontracted CASPT2 makes this reachable

OpenQP's CASPT2 is **uncontracted**: the first-order interacting space is the
external determinant space, an orthonormal basis whose labels do not depend on
the geometry. There is no internally-contracted overlap metric, no
linear-dependence removal, and no perturber normalization — the three things
that make an internally-contracted CASPT2 gradient hard. What remains is

```text
E   = E_ref + E2 ,   E_ref = <Psi_0|H|Psi_0> + V_NN
E2  = V . T ,        V = (P_Q H |Psi_0>) ,   T = -G(D)^-1 V
```

with `G` the (possibly regularized) denominator function.

## What is differentiated

```text
dE/dx = sum D^AO h^x + 1/2 sum Gamma^AO (..|..)^x - sum X^AO S^x + dV_NN/dx
```

with `D^AO`, `Gamma^AO` and the energy-weighted density `X^AO` built from the
**relaxed** densities. Unlike a converged state-specific CASSCF, a CASPT2 energy
is stationary in none of its parameters by itself, so every one of them carries
a multiplier:

- **the first-order amplitudes.** `E2 = V.T` is not stationary in `T` once a
  level shift is applied. Solving the multiplier equation gives `lambda = T`
  exactly for a single state — the Lagrangian collapses to the *shifted*
  Hylleraas functional — and a genuine solve in the multistate case.
- **the reference CI vector**, through one projected linear solve per root.
- **the XMS state rotation**, whose eigenvectors depend on the geometry.
- **the reference orbitals**, through a Z-vector against whatever actually fixes
  them: CASSCF stationarity for `reference=casscf`, RHF canonicality for
  `reference=casci`.

Both orbital Jacobians are closed form (the RHF one reduces to the textbook CPHF
`A` matrix). No finite difference appears anywhere in the implementation.

The effective-Hamiltonian eigenvectors of a multistate run need **no** response —
Hellmann-Feynman applies to an eigenvalue — but OpenQP checks the gap to the
neighbouring root before differentiating, because a degenerate root has no
well-defined mixing vector.

## Denominator shifts are exact

A shifted amplitude is not stationary in the Hylleraas functional, so the shift
reaches the derivative through the derivative of the denominator *function*, not
just through the amplitudes. OpenQP evaluates it exactly, using the
divided-difference (Daleckii–Krein) derivative of the matrix function `G(D)`. For
every regularization the input validator admits,

```text
G(w) = d + e/d ,  d = w + level_shift ,  e in {0, imaginary_shift^2, edshft}
```

so the derivative weight has rank at most two and no `n_ext x n_ext` object is
ever formed. `level_shift`, `imaginary_shift` and `edshft` therefore all carry an
analytic gradient of the same quality as the unshifted case.

`ipea_shift` is different: it biases the active *diagonal* of the zeroth-order
Hamiltonian in a particular orbital basis, which is not invariant under rotations
inside the active block. A nonzero IPEA shift is refused rather than
differentiated approximately.

## Scope

| Supported | Not supported |
| --- | --- |
| `caspt2`, `mrmp2` (single state) | `ms-caspt2` — the **multi-set** construction (per-state orbitals, per-state full-Fock-matrix `H0`, inter-state Löwdin-minor rotations); use `xms-caspt2` |
| `mcqdpt2` (single-set multistate) | `h0=dyall` (NEVPT2) |
| `xms-caspt2`, `xmcqdpt2` | `contraction=strong` (SC-NEVPT2) |
| `level_shift`, `imaginary_shift`, `edshft` | `ipea_shift` other than `0.0` |
| `reference=casci` (RHF orbitals) and `reference=casscf` (state-specific or state-averaged) | `[cas] orbital_source` reading orbitals from a file — imported orbitals are not a differentiable function of the geometry |
| the PT2 frozen core, `[pt2] frozen` | |

Under `gradient=auto` an unsupported combination falls back to central
differences and records why in the log. Under `gradient=analytic` it raises.

## When the analytic route does not apply

The derivation rests on conditions that hold almost everywhere and can fail at
one particular geometry. Each message carries the offending number:

- the CASCI reference orbitals no longer diagonalize the RHF Fock, or the CASSCF
  reference is not stationary (`g_orb` above tolerance);
- the orbitals are not semicanonical, so the zeroth-order Hamiltonian the
  gradient differentiates is not the one the energy used;
- the orbital-response system is singular, meaning two reference orbitals are
  degenerate under the reference condition and the multipliers are not
  determined;
- the requested effective-Hamiltonian root is degenerate with a neighbour, or
  the XMS model-space Fock has two degenerate eigenvalues;
- the gradient module's own reconstruction of the energy disagrees with the
  reported one, which would mean the gradient does not belong to the printed
  energy.

These are preconditions of the **route**, not of the energy, so `auto` treats
them exactly like an unsupported variant: it falls back to central differences
and logs the condition that failed. One caveat is fundamental rather than a
matter of speed. When the failed precondition is a *degenerate* (or
near-degenerate) effective-Hamiltonian root, the central difference is taken of
the ascending *sorted* energies, and a displaced geometry can cross the seam
and switch branches. The sorted energy is not differentiable at the crossing,
so the result is step- and coordinate-dependent and is **not** the gradient of
one physical state — the log says a sorted difference was used, and it should
be treated as a search direction, not as a state gradient. A penalty-function
MECI search walks into exactly this case by construction, which it tolerates
because its objective is built from the average and the gap of the sorted
pair; on a two-iteration H4 `mcqdpt2` search, 98 steps take the analytic
gradient and 5 fall back at root gaps between 7e-7 and 9e-7 Eh. Anything that
needs the derivative of one diabatic-following state at a seam needs
state-following machinery that this fallback does not provide.

`gradient=analytic` refuses instead, naming the condition.

Three conditions are **not** routed and stop the run whatever `[pt2] gradient`
says, because they are about the caller or the build rather than the geometry:
no PT2 energy on the molecule, a `liboqp` without the `caspt2_gradient` entry
point, and a nonzero status out of the gradient kernel.

## Selecting the state

`[properties] grad` (or `[optimize] istate`) addresses the published PT2 states.
Single-state CASPT2 publishes one, so the only valid selector is `0`. A
multistate or XMS run publishes the diagonalized effective-Hamiltonian roots in
ascending order.

## Example

The complete runnable request, in concise `.oqp` form:

```text
caspt2/sto-3g grad(S0) pt2(gradient=analytic)
cas(active_electrons=2,active_orbitals=2,frozen_core=1)
geom="h4.xyz"
```

The equivalent legacy `.inp` sections are (add an `[input] system=` geometry
block to run it):

```ini
[input]
runtype=grad
basis=sto-3g
method=caspt2

[properties]
grad=0

[cas]
active_electrons=2
active_orbitals=2
frozen_core=1
orbital_source=rhf

[ci]
nroot=1
solver=dense
integral_backend=native

[pt2]
reference=casci
gradient=analytic
```

The shipped examples are `examples/WF_methods/H4_CASPT2_grad.inp` (analytic),
`H4_CASPT2_numgrad.inp` (the central-difference companion), and
`H4_XMS-CASPT2_grad.inp` (multistate/XMS).

## Accuracy

Against a five-point finite difference of independently computed total energies,
at `h = 1e-3` Bohr on an off-symmetry H4 (LiH for the frozen-core row):

| Case | `max abs(analytic - 5-point)` |
| --- | --- |
| `caspt2`, CASCI reference | `6e-11` |
| + `level_shift=0.15` | `6e-11` |
| + `imaginary_shift=0.20` | `6e-11` |
| + `edshft=0.05` | `7e-11` |
| `caspt2`, state-specific CASSCF reference | `4e-9` |
| `xms-caspt2`, SA-CASSCF reference, 2 roots | `3e-10` / `2e-10` |
| `mcqdpt2`, 2 roots | `6e-11` / `2e-10` |
| `xms-caspt2`, 2 roots | `6e-11` / `1e-10` |
| `mrmp2` / `mcqdpt2` / `xmcqdpt2` on their default direct engine | `7e-11` / `2e-10` |
| `caspt2`, 6-31G | `1e-9` |
| `caspt2`, LiH with the default frozen core | `1e-11` |
| `caspt2`, LiH/cc-pVDZ (first case with `d` shells) | `2e-9` |
| `caspt2`, BeH2/6-31G, frozen core splitting the inactive block | `2e-10` |
| `caspt2`, LiH/cc-pVTZ (`f` shells, 44 basis functions), `h = 2e-3` | `1e-7` |
| `xms-caspt2`, H2O/STO-3G, three centres, `h = 2e-3` | `6e-10` / `6e-9` |

The CASSCF rows were measured with the reference converged to
`[casscf] gradient_norm_tol = 1e-9`. That matters: at the default `1e-6` they
read `5e-8` and `5e-7`, and tightening the reference *alone* moves them to
`4e-9` and `3e-10` without touching the gradient. The residual there is the
finite-difference side inheriting the CASSCF convergence of every displaced
point, not the analytic derivative. The cc-pVTZ row is looser for the same kind
of reason: at 44 basis functions and total energies near `-8` Ha, the
reference's convergence divided by `12h` puts the floor near `1e-8`.

Translational invariance `max abs(sum_A dE/dR_A)` and rotational invariance hold
to `1e-15` and `1e-13`. The numbers were reproduced to the quoted digit on
macOS/arm64 with Accelerate ILP64, Linux/x86-64 with OpenBLAS ILP64, and KNU
chc4 through Slurm with GCC 12.3.0 and MKL ILP64.

!!! warning "Comparing against a finite difference yourself"

    The same caution as for CASSCF applies, and more strongly: each displaced
    PT2 run repeats the whole SCF/CASSCF/semicanonicalization pipeline, so a
    disagreement can come from a displaced point landing on a different solution
    branch or reordering nearly degenerate orbitals across the
    `[cas] sort_orbitals=energy` active-space boundary. Sample the energy along
    the coordinate and confirm it is smooth before treating a disagreement as a
    gradient defect.

## See also

- [Wavefunction methods](../keywords/wavefunction.md) — `[cas]`, `[ci]`,
  `[pt2]`, `[casscf]` and `[state_average]` keywords.
- [CASSCF Nuclear Gradient](casscf-gradient.md) — the state-specific CASSCF
  derivative this one builds its reference on.
- [Optimization](optimization.md) — gradient-driven geometry optimization.
