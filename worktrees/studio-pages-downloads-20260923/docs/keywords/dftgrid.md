# `[dftgrid]`

The `[dftgrid]` section controls DFT quadrature and optional overrides for
hybrid or range-separated functional parameters. Most users should leave these
defaults unchanged unless reproducing a benchmark or testing a functional.

## Keywords

### `hfscale`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `-1.0` |
| Used by | exact-exchange scale override |

Overrides the Hartree-Fock exchange fraction. Negative values mean that OpenQP
uses the value associated with the selected functional.

### `cam_flag`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | range-separated functional setup |

Enables explicit CAM/range-separated parameter handling.

### `cam_alpha`, `cam_beta`, `cam_mu`

| Field | Value |
| --- | --- |
| Type | float |
| Defaults | `-1.0`, `-1.0`, `-1.0` |
| Used by | CAM/range-separated functional setup |

Override the CAM alpha, beta, and mu parameters. Negative values mean that
OpenQP uses the selected functional default.

### `rad_type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `ta` |
| Used by | radial quadrature |

Selects the radial grid family.

### `rad_npts`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `96` |
| Used by | radial quadrature |

Sets the number of radial points.

### `ang_npts`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `302` |
| Used by | angular quadrature |

Sets the number of angular grid points.

### `partfun`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `ssf` |
| Used by | molecular grid partitioning |

Selects the partition function used to assign atomic grid weights.

### `pruned`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `SG2` |
| Used by | pruned DFT grids |

Selects a pruned grid preset: `SG0`, `SG1`, `SG2`, or `SG3`. Write
`pruned=none` (in `.inp` input an empty value, `pruned=`, means the same) to
use the unpruned `rad_npts` × `ang_npts` grid at every radius; some regression
examples do this so their reference numbers do not depend on a pruning
scheme.

### `grid_ao_pruned`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | AO grid pruning |

Enables AO screening on grid points.

### `grid_ao_threshold`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0e-15` |
| Used by | AO grid pruning |

AO values below this threshold can be treated as negligible during grid pruning.

### `grid_ao_sparsity_ratio`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.9` |
| Used by | AO grid pruning |

Controls the sparsity threshold used by AO grid pruning.
