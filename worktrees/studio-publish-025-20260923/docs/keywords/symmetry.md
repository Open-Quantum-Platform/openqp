# `[symmetry]`

The `[symmetry]` section controls point-group metadata, the guarded abelian
integral reduction, and optional experimental symmetry tiers. Labeling can be
useful even when reductions are disabled.

!!! note "Detection is on by default"

    `enabled` defaults to `true`. Detection and labeling are active unless you
    turn them off. The two *reduction* keywords —
    [`use_integral_symmetry`](#use_integral_symmetry) and
    [`use_response_symmetry`](#use_response_symmetry) — remain `False` by
    default, so the geometry is not reoriented, the full integral list is still
    used, and reference-SCF or ordinary HF/DFT ground-state energies are
    unchanged. An MRSF physical `S0` remains a response root and can change
    when symmetry information repairs response-block coverage.

    Detection is **not** purely cosmetic, though. The detected irreps are what
    let the excited-state solver reach every symmetry block, so excited-state
    results can differ from a run with detection off — and where they do, the
    difference is a correction. See
    [Molecular symmetry](../workflows/symmetry.md).

## Keywords

### `enabled`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `true` |
| Values | `false`, `true`, `auto` |
| Used by | symmetry detection |

Controls symmetry detection. When on, OpenQP detects the point group and its
Abelian subgroup, labels orbitals, states and normal modes, and makes the
per-pair irrep table available to the excited-state solver.

What it does **not** do: it does not reorient the molecule and does not reduce
the integral list. Those follow `use_integral_symmetry` only.

Set `enabled=false` to restore the pre-detection behavior — no labels, no irrep
table, and the excited-state trial vectors chosen purely by orbital-energy gap
(see [Molecular symmetry](../workflows/symmetry.md) for why that can matter).

### `point_group`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` |
| Used by | requested point group |

Requested point group, or `auto` for automatic detection.

### `subgroup`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` |
| Used by | Abelian subgroup selection |

Requested Abelian subgroup, or `auto` for automatic choice.

### `label_mo`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | MO labeling |

Prints the molecular-orbital symmetry table. This is a **display** control
only: setting it to `False` suppresses the printed table but does not disable
the orbital irreps themselves, which the excited-state solver needs regardless.

### `label_states`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | response-state labeling |

Labels excited states by symmetry where possible. When `[tdhf] multiplicity` is
set, the labels are reported as spin-resolved terms (`1A1`, `3B2`, ...) rather
than bare irreps.

### `label_modes`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | vibrational-mode labeling |

Labels vibrational modes by symmetry where possible. The labels appear as a
`Symmetry` column in the printed frequency table:

```
   Mode  Symmetry     Frequency(cm-1)      IR(km/mol)        Raman(activity)
      1        a1           2043.11         0.492662           104.268372
      2        a1           4488.55         1.622292           558.762685
      3        b2           4790.80         0.866146           220.085777
```

The column is present whenever mode labels could be assigned, including for a
Hessian restored with `[hess] read=true`.

### `use_integral_symmetry`

| Field | Value |
| --- | --- |
| Type | string/boolean-like |
| Default | `True` |
| Used by | integral symmetry reduction |

Enables the guarded abelian petite-list/skeleton-Fock reduction for supported
`energy` and `grad` calculations. The default route verifies the AO operator
against the overlap matrix, verifies density invariance before each reduced
build, and falls back to the complete integral path when a safety condition is
not met. ROHF calculations using pFON also stay on the complete C1 integral
path because fractional-occupation convergence can be sensitive to the changed
summation order. Set this to `False` to disable integral reduction. The value
`full` requests the experimental non-abelian tier and should be validated
against a C1 reference run.

### `move_to_standard_frame`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | integral symmetry reduction |

Controls whether a calculation with `use_integral_symmetry=true` is rotated and
translated to a standard symmetry frame before the integral reduction is
staged. The default `False` keeps the molecule in its input frame; OpenQP then
uses the corresponding input-frame symmetry operators for the integral and
gradient reductions. The no-move route is available for `energy` and `grad`
calculations and uses the abelian integral-symmetry tier. The non-abelian
`use_integral_symmetry=full` tier requires `move_to_standard_frame=true`; OpenQP
rejects `full` together with `move_to_standard_frame=false` instead of silently
using a smaller group. As with every integral-symmetry calculation, validate new
systems against a run with `use_integral_symmetry=false`.

```ini
[symmetry]
enabled=true
use_integral_symmetry=true
move_to_standard_frame=false
```

**This is the keyword that reorients the molecule** into the symmetry standard
orientation. Turning symmetry detection on does not.

Accepts `true` for the Abelian subgroup (machine-exact) or `full` for the
complete point group (a larger reduction, accurate to roughly 1e-7).

!!! warning "`full` is declined for DFT"

    The `full` tier reduces the two-electron integrals over the complete point
    group, while the exchange-correlation grid necessarily stays on the Abelian
    operations — Lebedev angular grids are not invariant under C3/C6 rotations.
    The two halves would then reduce over different groups, which is a measured
    error of about 3e-04 Hartree.

    When a functional is set, OpenQP therefore falls back to the exact Abelian
    tier and says so in the symmetry log block. You still get the reduction,
    just the tier that is exact. Remove the `[input] functional` if you need
    the full group — OpenQP selects DFT with `method=hf` *plus* a functional,
    so `method=hf` alone is not what distinguishes the two.

The fallback is never silent: if the reduction is requested but does not engage
for any reason, the log states that the run used the full (C1) integral list.

Accepted values are `false`, `true` and `full`. `true` selects the Abelian
subgroup; `full` requests the complete point group.

#### What the reduction covers

Enabling this keyword reduces two separate parts of the calculation:

* **The two-electron integral list.** Symmetry-equivalent shell quartets are
  computed once and reweighted, and the resulting skeleton Fock matrix is
  projected back.
* **The exchange–correlation quadrature.** Grid slices belonging to
  symmetry-equivalent atoms are computed for one representative atom per orbit
  and scaled by the size of that orbit; the images are skipped.

The XC reduction engages only when the molecule actually has an orbit to
exploit — that is, when at least one atom is a symmetry image of another. For a
C1 molecule, or a geometry in which every atom sits on a special position, every
atom is its own orbit, so nothing is skipped and the reduction reports itself
inactive rather than doing identity work.

The XC reduction is **not compatible with the cross-iteration Φ cache**. When it
engages, `[scf] xc_phi_cache` is disabled automatically for that run: a cached
grid block built under the reduction would otherwise be replayed unreduced and
give a silently wrong `E_xc`.

#### Geometry and run types

Enabling the reduction **reorients the molecule** into the symmetry standard
orientation, and all outputs — geometry, gradients, orbitals — are reported in
that frame.

The reduction is applied only to `runtype = energy` and `runtype = grad`. Other
run types (optimisation, numerical Hessian, NEB, MEP, properties) are excluded,
because they either displace the geometry themselves or consume an externally
supplied one, and a per-step frame change would be assembled inconsistently.

#### `full` with a functional

`use_integral_symmetry = full` is declined for any run that sets `[input]
functional`, and the exact Abelian tier is used instead. The two halves of the
Fock matrix would otherwise reduce over different groups: Lebedev angular grids
are invariant under the axis-aligned octahedral operations but not under C3 or
C6, so the XC half can only ever be Abelian-symmetric while the JK half would be
forced symmetric under the full group. The mismatch is a measured error rather
than a rounding difference. The log states when this fallback happens; it is not
silent.

Hartree–Fock runs are unaffected and may use `full`.

### `use_response_symmetry`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | response solver symmetry reduction |

Enables irrep-blocked response updates — the experimental projection that
confines Davidson residuals to a single irrep. The checker marks this as
experimental and recommends validating excitation energies against an unblocked
run.

This is **separate** from the symmetry-block coverage of the initial trial
vectors, which is not experimental and follows detection alone. The symmetry log
block distinguishes the two:

| `response blocking` | Meaning |
| --- | --- |
| `pair_table_staged` | Irrep table built for trial-vector coverage; the residual projection is **off** |
| `active` | The experimental residual projection is running |

### `tolerance`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0e-5` |
| Used by | symmetry detection |

Geometry tolerance for symmetry detection. Must be positive.

### `strict`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | requested/detected group matching |

Requires stricter agreement between requested and detected symmetry labels.
