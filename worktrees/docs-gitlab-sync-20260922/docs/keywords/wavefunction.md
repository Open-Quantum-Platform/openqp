# Wavefunction methods

OpenQP provides native determinant-space FCI, CASCI, CASSCF, state-averaged
CASSCF, and second-order multireference perturbation methods. These workflows
require an RHF reference: leave `[input] functional` empty and use
`[scf] type=rhf` with `multiplicity=1`.

All of them run `runtype=energy`. State-specific `method=casscf` and dedicated
`method=sa-casscf` use analytic nuclear gradients, and so do single-state
CASPT2/MRMP2, MCQDPT2, XMS-CASPT2/XMCQDPT2, and the single-state strongly
contracted NEVPT2 variant (SC-NEVPT2). SA-CASSCF provides the gradient of
either the weighted objective or an individual averaged root; the latter
includes the coupled orbital and CI response through a Z-vector. The
compatibility spelling `method=casscf` with `[state_average] enabled=true`,
multi-set MS-CASPT2, and partially contracted NEVPT2 use Cartesian central
differences of converged total energies. The supported analytic and numerical
derivatives also connect to the gradient-driven run types described below. See
[CASSCF Nuclear Gradient](../workflows/casscf-gradient.md),
[CASPT2 Nuclear Gradient](../workflows/caspt2-gradient.md), and
[SC-NEVPT2 Nuclear Gradient](../workflows/sc-nevpt2-gradient.md). CASCI and FCI
remain energy-only.

The method selects which sections are read:

| `[input] method` | Required controls | Optional controls |
| --- | --- | --- |
| `fci` | `[fci]` | — |
| `casci` | `[cas]`, `[ci]` | `[state_average]` |
| `casscf`, `sa-casscf` | `[cas]`, `[ci]`, `[casscf]` | `[state_average]` |
| `caspt2`, `ms-caspt2`, `xms-caspt2` | `[cas]`, `[ci]`, `[pt2]` | `[casscf]`, `[state_average]` |
| `mrmp2`, `mcqdpt2`, `xmcqdpt2` | `[cas]`, `[ci]`, `[pt2]` | `[casscf]`, `[state_average]` |

## Minimal inputs

FCI uses its self-contained section:

```ini
[input]
system=
  H  0.0  0.0  0.0
  H  0.0  0.0  0.74
method=fci
runtype=energy
basis=sto-3g
functional=

[scf]
type=rhf
multiplicity=1

[fci]
nroot=1
solver=auto
```

CASCI and CASSCF separate the active-space, CI, and orbital-optimization
controls:

```ini
[input]
system=
  H  0.0  0.0  0.0
  H  0.0  0.0  0.74
method=casscf
runtype=energy
basis=sto-3g
functional=

[cas]
active_electrons=2
active_orbitals=2
frozen_core=0

[ci]
nroot=2
solver=auto
target_spin=singlet

[casscf]
root=0
converger=trah
hessian=fd
```

For equal-weight two-state SA-CASSCF, add:

```ini
[state_average]
enabled=true
nstate=2
equal_weights=true
target_roots=0,1
```

For CASPT2 on the converged reference, add `[pt2]` and set
`[input] method=caspt2`:

```ini
[pt2]
reference=casscf
h0=fock
contraction=uncontracted
ipea_shift=0.0
```

Runnable examples are in
[`examples/WF_methods`](https://github.com/Open-Quantum-Platform/openqp/tree/main/examples/WF_methods).

## `[fci]`

`[fci]` is used only by `method=fci`. It combines the active-space and solver
controls that CAS methods keep in separate sections.

| Keyword | Default | Meaning |
| --- | --- | --- |
| `nroot` | `1` | Number of CI roots returned. |
| `active_electrons`, `active_orbitals` | `0`, `0` | Optional restricted active space; zeros mean full CI. |
| `frozen_core` | `0` | Number of doubly occupied orbitals folded into the core. |
| `solver` | `auto` | `auto`, `dense`, or `davidson`. |
| `max_det` | `5000` | Determinant-count safety limit. |
| `max_memory` | `2048` | Combined working-memory ceiling in MiB. |
| `eig_tol` | `1.0e-10` | Eigenpair residual tolerance. |
| `davidson_maxiter`, `davidson_subspace` | `100`, `0` | Davidson iteration and subspace limits; zero subspace selects the automatic size. |
| `target_spin` | `any` | `any`, a multiplicity name such as `singlet`, or an integer multiplicity. |
| `irrep` | `any` | Return the lowest root of this spatial irrep instead of the lowest root overall. See [Selecting a root by spatial symmetry](#selecting-a-root-by-spatial-symmetry). |
| `irrep_min_purity` | `0.5` | Fraction of a root's weight that must sit in its dominant irrep before the root is accepted as belonging to it. |
| `print_ci_vectors`, `save_ci_vectors`, `save_rdm` | `false` | Optional CI diagnostics and artifacts. |

The only integral backend is `integral_backend=native`. Use
`integral_cutoff` to screen small determinant-Hamiltonian contributions.

## `[cas]`

`[cas]` defines the orbital partition shared by CASCI, CASSCF, and PT2.

| Keyword | Default | Meaning |
| --- | --- | --- |
| `active_electrons` | `0` | Total active electrons; required by the compact helpers. |
| `active_orbitals` | `0` | Number of active spatial orbitals. |
| `frozen_core` | `0` | Number of frozen doubly occupied orbitals. |
| `active_orbital_indices`, `core_orbital_indices` | empty | Explicit 1-based orbital selections. |
| `orbital_source` | `rhf` | `rhf`, `json`, or another supported orbital source. |
| `orbital_file` | empty | File used by a non-RHF orbital source. |
| `sort_orbitals` | `energy` | Ordering used before the active partition is formed. |
| `localize` | `none` | Optional orbital localization. |
| `max_det`, `max_memory` | `5000`, `2048` | Determinant and combined memory safety limits. |

Explicit scattered orbital selections are supported by CASCI. CASSCF and PT2
currently require a contiguous core/active partition.

## `[ci]`

`[ci]` controls the active-space eigensolver for CAS methods.

| Keyword | Default | Meaning |
| --- | --- | --- |
| `nroot` | `1` | Number of reference roots to solve. |
| `solver` | `auto` | `auto`, `dense`, or `davidson`. |
| `eig_tol` | `1.0e-10` | Eigenpair residual tolerance. |
| `davidson_maxiter`, `davidson_subspace` | `100`, `0` | Davidson limits. |
| `target_spin` | `any` | Requested spin multiplicity. |
| `irrep` | `any` | Return the lowest root of this spatial irrep instead of the lowest root overall. See [Selecting a root by spatial symmetry](#selecting-a-root-by-spatial-symmetry). |
| `irrep_min_purity` | `0.5` | Fraction of a root's weight that must sit in its dominant irrep before the root is accepted as belonging to it. |
| `root_tracking` | `energy` | Root ordering/tracking policy. |
| `print_ci_vectors`, `save_ci_vectors`, `save_rdm` | `false` | Optional diagnostics and artifacts. |

`nroot` must cover every CASSCF root, state-average target, and PT2 target.

### Selecting a root by spatial symmetry

`irrep` asks a correlated method for a state by its spatial symmetry rather
than by its position in the energy ordering. It is the only way to do so: the
solvers otherwise return the lowest roots, and an excited state of a given
symmetry can sit anywhere in that list.

```
[ci]
nroot=1
irrep=b1

[symmetry]
enabled=true
```

`[symmetry] enabled=true` is required, because the request is resolved against
the per-MO irrep labels that symmetry detection produces. The name must be an
irrep of the detected abelian subgroup; the run reports which subgroup was
detected and lists the available names when one does not match.

Determinants are classified by the direct product of the irreps of their
occupied active orbitals, which gives each root a weight per irrep. A root is
then accepted on two conditions, in this order: the irrep holding the most of
that weight -- the root's dominant irrep -- must be the one requested, and it
must hold at least `irrep_min_purity` of the total. Lowering the threshold
therefore never admits a root whose dominant irrep is some other one; it only
relaxes how dominant the requested irrep has to be. Determinants that could not
be classified stay in the total, so they lower the purity rather than quietly
vanishing from it.

The default of `0.5` is the weakest threshold that still identifies a dominant
irrep; raise it to demand a cleaner symmetry eigenstate. Values outside
`0 < value <= 1` are rejected rather than clamped -- below zero the filter would
accept anything, above one it would accept nothing, and `NaN` would disable the
test silently.

The request is refused, rather than partly applied, in three cases:

* **CASSCF and SA-CASSCF.** The orbital optimizer does not follow a
  symmetry-selected root across macroiterations, so the request would be
  honoured by the final CI and ignored by everything that produced the
  orbitals. Use `method=casci`.
* **An unclassified active orbital.** Symmetry detection leaves an orbital it
  cannot assign -- a near-degenerate pair rotated by the SCF, for instance --
  unlabelled. Such an orbital would be coded as totally symmetric, so every
  determinant containing it would be *mislabelled* rather than declined.
* **A solve the native engine declines.** The Python fallback driver has no
  irrep classification, and the two paths disagreeing about which root a run
  returns is worse than not offering the feature on the fallback.

Runnable examples:
[`H2O_CASCI_IRREP_B1`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/H2O_CASCI_IRREP_B1.inp)
for `[ci] irrep`, and
[`H2_FCI_IRREP_B1U`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/WF_methods/H2_FCI_IRREP_B1U.inp)
for `[fci] irrep` combined with `target_spin`. Both predict their reference
from the full spectrum of the same deck, so a filter that silently did nothing
would be caught.

## `[casscf]`

| Keyword | Default | Meaning |
| --- | --- | --- |
| `max_macro_iterations` | `20` | Orbital-optimization iteration limit. |
| `root` | `0` | Zero-based state-specific objective root. |
| `converger` | `trah` | `trah`, `ah`, `twophase`, `diis`, or `auto`. |
| `hessian` | `fd` | Finite-difference (`fd`) or analytic orbital Hessian. |
| `gradient_norm_tol` | `1.0e-6` | Orbital-gradient convergence threshold. |
| `energy_decrease_tol` | `1.0e-10` | Energy-change threshold. |
| `step_norm_tol` | `1.0e-8` | Orbital-step threshold. |
| `max_rotation_norm` | `0.2` | Maximum orbital-rotation norm. |
| `canonicalize` | `true` | Canonicalize the converged orbitals. |
| `gradient_state` | `averaged` | For dedicated `method=sa-casscf`, differentiate the weighted objective (`averaged`) or the named zero-based CI root. |
| `zvector_tol` | `1.0e-8` | Relative cutoff for the individual-state Z-vector null space, leakage, and residual tests. |
| `zvector_degeneracy_tol` | `1.0e-8` | Root-gap threshold in Hartree below which an individual-state derivative is refused. |
| `grad_step` | `1.0e-3` | Central-difference half-step in Bohr for a nuclear gradient. |
| `grad_guess` | `cold` | `cold` restores the same central-point SCF data before every displacement; `warm` follows the preceding solution serially. |
| `grad_gap_warn` | `1.0e-5` | Energy-gap threshold in Hartree used in the displaced-root ordering warning. |
| `grad_ranks_per_group` | `0` | MPI ranks assigned to one displaced energy; zero selects the automatic distribution. |

The `ah_*`, `diis_*`, and diagnostic keywords provide detailed control of the
selected converger. `max_macro_iterations=0` runs the fixed-orbital CASCI
scaffold and performs no orbital optimization.

`gradient_norm_tol` and `hessian` refer to orbital optimization. In
particular, `hessian=analytic` selects an analytic orbital Hessian; it does not
select the nuclear derivative. State-specific CASSCF has a separate analytic
nuclear-gradient kernel.

`root` also selects which state `runtype=grad` differentiates. Tighten
`gradient_norm_tol` to `1.0e-7` for gradient runs: the analytic gradient is
valid at a stationary point, and its error is first order in the converged
orbital-gradient norm. See
[CASSCF Nuclear Gradient](../workflows/casscf-gradient.md).

The analytic state-specific result has one reported row: use `[properties]
grad=0` for a direct gradient and `[optimize] istate=0` for `optimize`, `ts`,
`mep`, or `irc`. The physical CI state is selected only by `[casscf] root`.

For dedicated `method=sa-casscf`, `gradient_state=averaged` selects the
stationary weighted objective, while an integer selects the corresponding
physical CI root and activates its coupled orbital and CI Z-vector response. For
a direct individual-root gradient, `[properties] grad` must equal that physical
root. Gradient-driven optimization also requires `[optimize] istate` to match
the root and `target_roots` to be the contiguous sequence beginning at zero.
The weighted objective is currently available only for `runtype=grad`, because
it is not one of the reported state energies used by the optimizers.

The `grad_*` controls in the table apply only to the numerical compatibility
path, `method=casscf` with state averaging enabled. They do not replace or tune
either analytic derivative.

## `[state_average]`

| Keyword | Default | Meaning |
| --- | --- | --- |
| `enabled` | `false` | Enable a state-averaged CAS objective. |
| `nstate` | `0` | Number of averaged states. |
| `target_roots` | empty | Zero-based CI roots included in the average. |
| `equal_weights` | `true` | Assign equal normalized weights. |
| `weights` | `1.0` | Explicit weights when `equal_weights=false`. |
| `root_tracking` | `overlap` | Accepted for compatibility; see below — roots are in fact followed by energy order. |

The number of weights must equal the number of target roots. Weights are
normalized internally and must contain a positive total.

Dedicated `method=sa-casscf` uses analytic nuclear gradients. The weighted
objective is stationary with respect to the optimized orbital and CI
parameters, so its derivative uses weight-averaged density matrices without a
response term. An individual root is not stationary with respect to the common
state-averaged orbitals; its analytic derivative therefore includes the
coupled orbital and CI response through a Z-vector.

The current molecular input requires equal weights over a contiguous root
block. Despite the `root_tracking=overlap` spelling, no overlap matching is
performed between orbital macroiterations: each CI solve is consumed in energy
order. Equal weights over a contiguous block are invariant when nearby roots
interchange, which is why unequal weights or a noncontiguous `target_roots`
subset — where an interchange would move weight onto different physical
states — are refused at validation.

`method=casscf` with `enabled=true` remains an explicit compatibility route
that uses central differences. It includes orbital relaxation at displaced
geometries and retains the `grad_*` controls above.

## `[pt2]`

| Keyword | Default | Meaning |
| --- | --- | --- |
| `reference` | `casscf` | `casscf` or fixed-orbital `casci` reference. |
| `variant` | `auto` | Single-state, multistate, or extended-multistate variant inferred from the method by default. |
| `h0` | `fock` | `fock` for CASPT2/QDPT or `dyall` for NEVPT2. |
| `contraction` | `uncontracted` | `uncontracted` or `strong` (SC-NEVPT2). |
| `target_roots`, `nroot` | empty, `0` | PT2 target model space; zero inherits `[ci] nroot`. |
| `ipea_shift`, `level_shift`, `imaginary_shift`, `edshft` | `0.0` | Denominator regularization controls. |
| `frozen` | `auto` | PT2 frozen core: freeze the standard deep-core orbitals automatically, use `0` to correlate all inactive orbitals, or give an explicit count. |
| `engine` | `auto` | Automatic, direct/Fortran, or dense QDPT engine. |
| `nproc` | `0` | Worker count; zero selects the automatic value. |
| `max_terms` | `30000000` | Direct-engine streamed-term limit. |
| `max_memory` | `2048` | PT2 memory ceiling in MiB, combined with the tighter `[cas]` ceiling. |
| `semi_canonical` | `true` | Semicanonicalize the reference orbitals. |
| `gradient` | `auto` | Nuclear-gradient route: `auto` (analytic where it applies — including the SC-NEVPT2 analytic scope — and central differences otherwise, both for an unsupported variant and for a geometry where a precondition of the derivation fails), `analytic` (refuse rather than fall back), or `numerical`. See [CASPT2 Nuclear Gradient](../workflows/caspt2-gradient.md) and [SC-NEVPT2 Nuclear Gradient](../workflows/sc-nevpt2-gradient.md). |
| `grad_step` | `1.0e-3` | Central-difference half-step in Bohr for a nuclear gradient. |
| `grad_guess` | `cold` | Displaced-geometry SCF starting-data policy (`cold` or `warm`). |
| `grad_gap_warn` | `1.0e-5` | Energy-gap threshold in Hartree used in the root-ordering warning. |
| `grad_ranks_per_group` | `0` | MPI ranks assigned to one displaced energy; zero selects the automatic distribution. |

Select NEVPT2 through `[input] method=caspt2` with `h0=dyall`; add
`contraction=strong` for strongly contracted NEVPT2. The Python
`job.nevpt2(...)` convenience helper lowers to this same sectioned form rather
than introducing a separate `method=nevpt2` input value. The
QDPT family (`mrmp2`, `mcqdpt2`, and `xmcqdpt2`) uses the diagonal-Fock direct
engine and accepts `edshft`; `edshft` is mutually exclusive with real or
imaginary level shifts.

## Nuclear gradients and geometry calculations

| Methods | Derivative | Supported run types |
| --- | --- | --- |
| State-specific `method=casscf` | Analytic derivative of `[casscf] root` | `grad`, `optimize`, `ts`, `mep`, `irc` |
| Dedicated `method=sa-casscf`, individual `gradient_state` | Analytic coupled orbital and CI Z-vector derivative | `grad`, `optimize`, `ts`, `mep`, `irc` |
| Dedicated `method=sa-casscf`, `gradient_state=averaged` | Analytic derivative of the weighted objective | `grad` |
| `method=casscf` with state averaging enabled | Central difference | `grad`, `optimize`, `ts`, `mep`, `irc` |
| CASPT2, MRMP2, MCQDPT2, XMS-CASPT2, XMCQDPT2 | Analytic (`[pt2] gradient`, `auto` falls back to central differences) | `grad`, `optimize`, `ts`, `mep`, `irc` |
| SC-NEVPT2 (strong contraction, Dyall `h0`) | Analytic when in scope (`[pt2] gradient`) | `grad`, `optimize`, `ts`, `mep`, `irc` |
| Multi-set MS-CASPT2, partially contracted NEVPT2 | Central difference | `grad`, `optimize`, `ts`, `mep`, `irc` |

FCI and CASCI remain energy-only. `meci`, `mecp`, and `neb` are not connected
to these wavefunction-gradient calculations.

For state-specific CASSCF, the physical CI state is selected by `[casscf] root`,
while `[properties] grad` and `[optimize] istate` remain zero because the result
has one reported row. For analytic SA-CASSCF, `gradient_state` is either
`averaged` or a physical root in `[state_average] target_roots`. A direct
individual-root gradient and its optimizer selector must name that same root.
The compatibility numerical route instead indexes the state list published by
`[state_average] target_roots`.

The SC-NEVPT2 analytic route requires `[pt2] h0=dyall`,
`contraction=strong`, `reference=casscf`, and a single state-specific root on a
closed-shell RHF reference. With `gradient=auto`, an out-of-scope input or a
geometry-specific applicability failure is reported and routed to central
differences. With `gradient=analytic`, the same condition is an error rather
than a silent method change. See
[SC-NEVPT2 Nuclear Gradient](../workflows/sc-nevpt2-gradient.md) for the full
contract and diagnostics.

The central-difference calculations currently retain energy ordering without
matching CI vectors between geometries. OpenQP warns when two reported roots become
closer than the displacement-induced energy change, but it cannot detect a
crossing with an uncomputed root.

## Python helpers

The fluent API writes the same sections and validates the required active
space:

```python
from oqp.openqp import OpenQP

job = OpenQP("h4_xms")
job.molecule(geometry="h4.xyz", charge=0)
job.caspt2(
    active_electrons=2,
    active_orbitals=2,
    variant="xms-caspt2",
    nroot=2,
    basis="sto-3g",
)
mol = job.run()
print(mol.get_results()["energies"])
```

Available helpers include `fci`, `casci`, `casscf`, `sa_casscf`, `caspt2`,
`nevpt2`, and `qdpt2`. A helper call resets method-owned state from an earlier
call on the same builder unless that state is explicitly supplied again.

For a gradient-driven state-specific CASSCF calculation,
`job.casscf(root=1, runtype="optimize")` selects physical root 1 and reports it
in result row zero. For analytic SA-CASSCF,
`job.sa_casscf(nstate=2, gradient_state=1, runtype="grad")` differentiates
physical CI root 1 and sets `[properties] grad=1`; omit `gradient_state` to
differentiate the weighted objective. The older `state=1` argument remains a
positional alias for the second entry of `target_roots`. The CASSCF helpers
also accept `grad_step`, `grad_guess`, `grad_gap_warn`, and
`grad_ranks_per_group` for the numerical compatibility route.
