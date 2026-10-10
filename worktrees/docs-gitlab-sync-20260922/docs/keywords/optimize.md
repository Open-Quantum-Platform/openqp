# `[optimize]`

The `[optimize]` section selects the optimization backend, target states, and
convergence thresholds for geometry and reaction-path workflows.

!!! note "Two input surfaces"

    Concise `.oqp` geometry drivers always use the native OpenQP engine and do
    not accept `lib`, `optimizer`, `step_size`, `step_tol`, or `mep_maxit`.
    Those five keys remain documented here for traditional sectioned `.inp`
    files and the compatible Python workflow API.

In Python scripts, use `job.workflow.optimize(...)` for the optimization runtype, common
optimization options, and backend options. For example,
`job.workflow.optimize(istate=0)` uses the native backend and its automatic DLC
coordinates. Supplying `coordsys` routes an explicit coordinate-system override
to the native optimizer while keeping state and convergence options in
`[optimize]`. The lower-level `job.optimize(...)` section helper remains
available for existing scripts.

## Keywords

### `lib`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `oqp` |
| Values | `oqp`, `geometric`, `scipy` |
| Used by | optimizer backend selection |

Traditional `.inp` and Python API backend selector. The native `oqp` backend
supports `optimize`, `ts`, `meci`, `mecp`, `tci`, `neb`, `irc`, and `mep`.
geomeTRIC supports `optimize`, `meci`, `mecp`, `ts`, `irc`, and `neb`; it is
retained as an optional compatibility path for advanced constraints beyond
native frozen distances. SciPy is a legacy backend for `optimize`, `meci`,
`mecp`, and `mep`.

DL-FIND is not a current user-facing optimizer backend.

### `optimizer`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `bfgs` |
| Values | `bfgs`, `cg`, `l-bfgs-b`, `newton-cg` |
| Used by | SciPy backend |

Selects the SciPy optimizer when `lib=scipy`.

### `maxit`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `30` |
| Used by | optimization iterations |

Maximum number of geometry steps. The shipped multistate MECI regression uses
`maxit=1` only as a short continuous-integration smoke; production searches
normally require a substantially larger limit.

### `states`

| Field | Value |
| --- | --- |
| Type | integer array |
| Default | empty |
| Used by | BaekA multistate MECI |

Ordered internal response roots for `[optimize] meci_search=baeka`. Supply at
least two distinct, ascending, consecutive values. In concise `.oqp`, users
write physical labels directly, for example
`meci(S0,S1,S2,S3,algorithm=baeka)`; OpenQP creates this internal list.

### `istate`, `jstate`, `kstate`

| Field | Value |
| --- | --- |
| Type | integer |
| Defaults | `1`, `2`, `3` |
| Used by | state-specific and crossing searches |

Target state indices. HF/DFT ground-state optimization uses `istate=0`.
TDHF/MRSF state-specific workflows use positive state indices. MECI uses
`istate` and `jstate` for the established two-state algorithms. The legacy TCI
compatibility route uses all three and requires `istate < jstate < kstate`;
new multistate inputs use [`states`](#states) with `meci_search=baeka`.

### `imult`, `jmult`

| Field | Value |
| --- | --- |
| Type | integer |
| Defaults | `1`, `3` |
| Used by | MECP multiplicities |

Multiplicity labels for MECP-style crossing searches. MECP requires different
values.

### `energy_shift`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0e-6` |
| Used by | optimization convergence |

Energy-change convergence threshold. BaekA evaluates this as a same-current-
`sigma` objective change: the previous geometry's mean-energy and penalty
terms are recombined using the current penalty weight before `delta_F` is
formed.

### `energy_gap`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0e-5` |
| Used by | crossing-point convergence |

Target energy-gap threshold for crossing searches. For BaekA, OpenQP tests the
outer span `E_last - E_first`, not the largest individual adjacent gap. The
canonical `algorithm=baeka` default is `1.0e-4` Hartree. Traditional sectioned
inputs retain the longstanding global `1.0e-5` default unless they explicitly
set `energy_gap=1.0e-4`. See the
[BaekA multistate MECI workflow](../workflows/baeka-multistate-meci.md).

### `meci_search`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` |
| Values | `auto`, `auglag`, `penalty`, `ubp`, `hybrid`, `baeka` |
| Used by | MECI search algorithm |

Selects the MECI search strategy. `auglag`, `penalty`, `ubp`, and `hybrid` use
the traditional `istate`/`jstate` pair. `baeka` accepts the two-or-more-root
[`states`](#states) list and selects the additive BaekA independent-gap
objective. Each maps to `algorithm=<name>` inside `meci(...)`, so BaekA is
`algorithm=baeka` and the augmented Lagrangian is `algorithm=auglag`.

`auto` resolves to `auglag` for two states and to `baeka` for three or more.
With the native optimizer, BaekA is also held back as a rescue on the recovery
budget if `auglag` does not meet the convergence criteria.

`auglag` is an augmented Lagrangian on top of the usual branching-plane update,
so it needs no derivative coupling. A least-squares multiplier removes the
mean-gradient component along the gap direction, leaving a gap term and a
projected mean gradient that are orthogonal. Because orthogonal terms cannot
cancel, a vanishing total gradient forces the energy gap to zero rather than
balancing it against the descent. `penalty` has no such multiplier: its
stationary point sits at a finite gap, so on its own it stalls above
[`energy_gap`](#energy_gap). The strength of the gap term is set by
[`gap_sigma`](#gap_sigma).

### `mecp_search`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` |
| Values | `auto`, `sqp`, `auglag`, `penalty`, `quad` |
| Used by | MECP (spin-crossing) search algorithm |

Selects the MECP search strategy. The two states differ in spin multiplicity,
so the derivative coupling vanishes and the branching space is the single
gradient-difference direction.

`auto` selects `sqp` with the native optimizer and `auglag` on the backends
that supply their own, since SQP replaces the optimizer rather than providing
an objective to it.

`sqp` solves the KKT equations of the constrained problem

```
minimise 0.5 (E_i + E_j)   subject to   E_j - E_i = 0
```

for the step and the multiplier together, using a damped BFGS approximation to
the Hessian of the Lagrangian started from the native model Hessian. The
multiplier is therefore a result rather than a formula, and there is no penalty
parameter: [`gap_sigma`](#gap_sigma) is not read. It works in delocalized
internal coordinates, so [`coordsys`](oqp.md) applies to the native optimizer
and not to this search, and it requires `lib=oqp`.

`auglag` is the same augmented Lagrangian described above; with `gap_sigma=1`
its converged form is the plain Bearpark gradient projection. It is what `auto`
selects on the scipy and geomeTRIC backends. `penalty` is the Levine-Martinez
smooth penalty.

`quad` is the legacy fixed-weight quadratic penalty
`0.5 (E_i + E_j) + gap_weight (E_j - E_i)^2`. Its stationary condition balances
the mean gradient against the gap term, and because the two are not orthogonal
they cancel at a residual gap of order `1/gap_weight`. It therefore cannot
satisfy a tight [`energy_gap`](#energy_gap) and is retained only to reproduce
runs made before the converging objectives were added.

### `gap_sigma`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `10.0` |
| Used by | `auglag` crossing searches |

Scales the gap term of the `auglag` objective relative to the projected mean
gradient. `1.0` reproduces the plain Bearpark gradient projection; the larger
default reaches the seam faster and keeps the quasi-Newton history better
conditioned. Must be positive: zero removes the gap term entirely, and a
negative value drives the search away from the seam.

### `pen_sigma`, `pen_alpha`

| Field | Value |
| --- | --- |
| Type | float |
| Global schema defaults | `1.0`, `0.0` |
| BaekA method defaults | `1.0`, `0.02` Hartree |
| Used by | penalty-function crossing searches |

For `meci_search=baeka`, `pen_sigma` is the initial current dimensionless
penalty weight and `pen_alpha` is a positive energy-valued smoothing parameter.
Canonical `algorithm=baeka` injects `alpha=0.02` Hartree and rejects explicit
`alpha=0`. A traditional BaekA input inherits the global `pen_alpha=0.0`
sentinel when the key is omitted; BaekA interprets that sentinel as the fixed
`0.02` Hartree method default. It does not derive `alpha` from a gradient RMS,
whose units are incompatible.

### `pen_incre`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0` |
| Used by | legacy penalty-function multiplier update |

Historical multiplicative update control used by the established `tci`
workflow and older penalty implementation:

```text
sigma_next = sigma_current * pen_incre
```

BaekA does not use this key. Its additive controls are `pen_delta` and
`pen_jump`.

### `pen_delta`, `pen_jump`

| Field | Value |
| --- | --- |
| Types | float; float array |
| Defaults | `0.025`; `10,10,25,25,100,100,1000,1000,3000` |
| Used by | BaekA adaptive sigma updates |

`pen_delta` is the small dimensionless additive increment after an ordinary
nonstationary BaekA macro-step:

```text
sigma_next = sigma_current + pen_delta
```

`pen_jump` is an ordered schedule of larger additive increments. OpenQP
consumes one entry only when the same-sigma projected stationary tests pass but
the outer energy span remains above `energy_gap`. Schedule entries are neither
multipliers nor thresholds. In concise `.oqp`, use the public names
`delta_beta` and `beta_schedule`. See the
[BaekA workflow](../workflows/baeka-multistate-meci.md).

### `gap_weight`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0` |
| Used by | crossing-search objective |

Weight applied to a crossing-search penalty objective. The public BaekA route
uses `sigma`, `alpha`, `delta_beta`, `beta_schedule`, and `gap`; for BaekA this
lower-level key is fixed at `1.0`, and `sigma` is the published penalty
multiplier. Other legacy crossing algorithms retain their established use of
`gap_weight`.

### `rmsd_grad`, `max_grad`

| Field | Value |
| --- | --- |
| Type | float |
| Defaults | `1.0e-4`, `3.0e-4` |
| Used by | gradient convergence |

RMS and maximum gradient convergence thresholds for the established
optimizers. BaekA reuses the legacy-named `rmsd_grad` value as the threshold
for the Euclidean norms `||g_parallel||_2 / sigma` and
`||g_perpendicular||_2`; those BaekA quantities are not RMS values. `max_grad`
remains a logged generic diagnostic and is not a BaekA termination test.
Concise BaekA input rejects an explicitly supplied `max_grad` to avoid
implying otherwise.

### `rmsd_step`, `max_step`

| Field | Value |
| --- | --- |
| Type | float |
| Defaults | `1.0e-3`, `2.0e-3` |
| Used by | step convergence |

RMS and maximum geometry-step convergence thresholds. They remain available
as generic optimization diagnostics but are not BaekA termination tests;
concise BaekA input rejects explicit `rmsd_step` and `max_step` controls.

### `step_size`, `step_tol`

| Field | Value |
| --- | --- |
| Type | float |
| Defaults | `0.1`, `1.0e-2` |
| Used by | SciPy/lightweight optimization paths |

Step-size controls used by legacy optimization paths.

These are SciPy-backend compatibility keys for traditional `.inp` files and
the Python API. They are not accepted in concise `.oqp` geometry drivers.

### `mep_maxit`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `10` |
| Used by | MEP workflow |

Maximum number of MEP iterations.

This is a legacy SciPy MEP key and is not accepted in concise `.oqp`. Use
`mep(points=N,step=...)`; `points` lowers to the native `maxit` limit and
`step` lowers to [`[oqp] mep_step`](oqp.md#mep_step).

### `init_scf`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | geometry-step SCF handling |

Requests initial SCF cycles during geometry steps.
