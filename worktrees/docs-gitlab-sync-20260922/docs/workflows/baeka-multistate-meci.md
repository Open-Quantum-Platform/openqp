# BaekA Multistate MECI

`BaekA` is OpenQP's adaptive penalty-function algorithm for minimum-energy
conical intersections involving two or more electronic states. Select it as an
algorithm of the ordinary `meci(...)` driver; it is not a separate three-state
calculation type.

Baek, Lee, Filatov, and Choi introduced and validated the independent-gap
adaptive construction for a minimum-energy three-state conical intersection.
OpenQP exposes the same mathematical construction for `N >= 2` states. The
two-state case is the one-gap specialization, while support for more than three
states is an OpenQP generalization. The cited paper validates the three-state
case; neither non-three-state case is directly validated by those benchmarks.
See the [method reference](../references.md#baeka-multistate-meci).

## Canonical `.oqp` Commands

The same public command covers two, three, and more states:

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1,algorithm=baeka)
geom="guess.xyz"
```

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1,S2)
geom="guess.xyz"
```

```text
mrsf(nstate=6)/bhhlyp/6-31g* meci(S0,S1,S2,S3)
geom="guess.xyz"
```

A two-state `meci(S0,S1)` call defaults to `auglag`, the augmented Lagrangian,
so write `algorithm=baeka` when BaekA is intended for two states. With three or
more states OpenQP infers `baeka`, because the other MECI algorithms are
two-state only; the option is then redundant and the canonical rendering omits
it (`meci(S0,S1,S2,algorithm=baeka)` and `meci(S0,S1,S2)` are the same request).

All selected states must be consecutive response roots in the same physical
spin manifold. OpenQP normalizes their order. Thus
`meci(S0,S1,S2,algorithm=baeka)` and
`meci(T0,T1,T2,T3,algorithm=baeka)` have valid shapes, whereas skipping an
intervening root or mixing singlets and triplets is rejected. A crossing
between different spin manifolds is an MECP problem.

A three-state production-style route is:

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1,S2,maxit=100,trust=0.15)
geom="guess.xyz"
```

Three or more states select BaekA by themselves, and the public BaekA controls
default to `sigma=1.0`, `alpha=0.02`, `delta_beta=0.025`,
`beta_schedule="10,10,25,25,100,100,1000,1000,3000"`, `gap=1e-4`,
`energy_shift=1e-6`, `rmsd_grad=1e-4`, `coordsys=auto`, and `trust_max=0.5`.
Write one only to change it, for example `meci(S0,S1,S2,gap=2e-4)`.

The concise route selects the native OpenQP optimizer automatically. Do not add
`lib=oqp`. MRSF also supplies its internal high-spin working reference
automatically; labels in `meci(...)` are the physical target surfaces.

## Objective and Independent Constraints

For `N` energy-ordered states, define the `N - 1` adjacent gaps
`d_a = E_(a+1) - E_a`. BaekA minimizes the mean state energy plus one penalty
for each independent gap:

```text
F = (1/N) sum_a E_a + sigma sum_a P(d_a)
P(d) = d^2 / (d + alpha)
```

Using adjacent gaps removes redundant pairwise conditions: making all `N - 1`
adjacent gaps zero makes every selected state degenerate. The objective
gradient uses the corresponding adjacent gradient differences. `alpha` and the
energy gaps have Hartree units; `sigma` and its additive increments are
dimensionless.

BaekA tests degeneracy with the outer span

```text
E_last - E_first <= gap
```

not with the largest individual adjacent gap. The default `gap=1e-4` Hartree is
the tight outer-span criterion used for this method.

## Adaptive Sigma State Machine

`sigma` is the current dimensionless penalty weight. The update depends on the
current optimization state:

1. After an ordinary, nonstationary macro-step, apply the small additive update
   `sigma <- sigma + delta_beta`.
2. If the projected stationary tests pass but the outer span is still above
   `gap`, consume the next entry `beta` from `beta_schedule` and apply the
   larger additive update `sigma <- sigma + beta`. The strengthened objective
   is re-formed immediately at the same geometry.
3. If the projected stationary tests and outer-span test all pass, the BaekA
   search is converged.

The `beta_schedule` entries are additive increments. They are neither
multipliers nor gap thresholds, and they are consumed only by case 2. Exhausting
the schedule while another large jump is required is an error rather than a
silent continuation with an invented value.

!!! important "BaekA does not use `pen_incre`"

    `pen_incre` is the historical multiplicative control for the older
    two-state penalty implementation. BaekA instead uses the dedicated
    additive controls `delta_beta` and `beta_schedule`, which lower to
    `pen_delta` and `pen_jump`. The earlier simplified native TCI prototype
    multiplied `sigma` by `pen_incre`; that behavior is not the exact BaekA
    algorithm.

## Projected Stationary Tests

The `N - 1` adjacent gap-gradient vectors define the constraint subspace.
OpenQP obtains its independent numerical basis with an SVD, then decomposes the
BaekA objective gradient into:

- a component parallel to the gap-constraint subspace, divided by the current
  `sigma`; and
- a perpendicular component along the intersection seam.

A local stationary point requires `||g_parallel||_2 / sigma` and
`||g_perpendicular||_2` to be at or below the legacy-named `rmsd_grad`
threshold. These are Euclidean norms, not RMS values. It also requires the
objective change to be at or below `energy_shift`. Because `sigma` evolves,
OpenQP recombines the previous
geometry's sigma-independent mean-energy and penalty terms using the **current
same sigma** before forming this `delta_F`; a mere change of penalty weight
therefore cannot masquerade as geometry progress or destroy a valid
same-objective comparison.

The complete BaekA convergence test is the same-sigma `delta_F`, the two
projected gradient tests, and the outer energy span. Generic `max_grad`,
`rmsd_step`, and `max_step` values may still appear in optimization diagnostics,
but they are not the BaekA termination criteria. Concise `.oqp` therefore
rejects those three controls when they are written explicitly in a BaekA call.

## Public BaekA Controls

| `.oqp` control | Legacy key | Default | Meaning |
| --- | --- | --- | --- |
| `sigma` | `pen_sigma` | `1.0` | Current initial dimensionless penalty weight. |
| `alpha` | `pen_alpha` | `0.02` Hartree | Positive energy smoothing parameter in `P(d)`. |
| `delta_beta` | `pen_delta` | `0.025` | Small dimensionless additive increment after an ordinary nonstationary step. |
| `beta_schedule` | `pen_jump` | `10,10,25,25,100,100,1000,1000,3000` | Ordered larger additive increments used only when locally stationary with an open outer span. |
| `gap` | `energy_gap` | `1e-4` Hartree | Maximum permitted outer span `E_last - E_first`. |

`alpha` must represent a positive energy. In concise `.oqp`, omitting it uses
`0.02` Hartree and explicitly writing `alpha=0` is rejected. The traditional
schema retains a global `pen_alpha=0.0` sentinel because the older two-state
penalty algorithm uses that value differently; under `meci_search=baeka`, an
omitted/sentinel `pen_alpha=0.0` means the fixed BaekA default `0.02` Hartree.
It does **not** request a gradient-derived automatic value. A gradient RMS has
different units and is not a valid replacement for energy-valued `alpha`.

Common native controls remain available: `maxit`, `energy_shift`,
`rmsd_grad`, `coordsys`, `trust`, and `trust_max`.

Use the public `gap` spelling for BaekA. The internal compatibility spelling
`energy_gap` represents the same setting, so supplying both is rejected as a
duplicate rather than choosing one silently.

`gap_weight` is fixed at `1.0` for BaekA. Adjust `sigma` instead; introducing a
second independent weight would change both the published objective and its
projected convergence interpretation.

## State Selection

- Supply at least two distinct, consecutive states of one multiplicity.
- Set `nstate` high enough to contain the highest requested state. Computing a
  few extra roots can help when inspecting nearby-state character.
- Start from a geometry reasonably close to the desired intersection seam.
- The optimizer follows selected response-root indices. It does not track
  electronic character or relabel roots during the geometry search.
- `auglag` is the two-state default, and `penalty`, `ubp`, and `hybrid` remain
  available. Select `algorithm=baeka` when the additive independent-gap
  formulation is intended, including for a two-state calculation.

## Python and Legacy `.inp` Configuration

The concise list of physical labels lowers to an ordered internal state list.
The Python workflow accepts the same public names:

```python
job.workflow.meci(
    algorithm="baeka",
    states=[1, 2, 3],
    sigma=1.0,
    delta_beta=0.025,
    beta_schedule=[10, 25],
    gap=1e-4,
)
```

The equivalent legacy `.inp` control surface is:

```ini
[input]
runtype=meci
method=tdhf
functional=bhhlyp
basis=6-31g*

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=6

[optimize]
lib=oqp
meci_search=baeka
states=1,2,3,4
maxit=100
pen_sigma=1.0
pen_alpha=0.02
pen_delta=0.025
pen_jump=10,10,25,25,100,100,1000,1000,3000
energy_gap=1e-4

[oqp]
coordsys=auto
trust=0.15
trust_max=0.5
```

`[optimize] states` is authoritative for BaekA and must contain two or more
distinct, ascending, consecutive internal roots. The older `istate`/`jstate`
pair remains the traditional two-state spelling for the other MECI algorithms.

## `tci` Compatibility

Existing `tci(S0,S1,S2,...)` and sectioned `runtype=tci` inputs remain
supported with their pre-existing three-state adaptive-penalty algorithm and
multiplicative `pen_incre` update. They are not aliases for BaekA. New BaekA
inputs should use the explicit general form:

```text
meci(S0,S1,S2,algorithm=baeka)
```

The shipped `C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.inp` remains a regression for that
legacy behavior. The separate `C2H4_BHHLYP-MRSFTDDFT_BAEKA_OQP.inp` regression
exercises the new MECI algorithm. Their small iteration limits keep continuous
integration short; they are smoke tests, not converged scientific templates.
