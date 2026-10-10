# `[pcm]`

The `[pcm]` section controls implicit-solvent settings. The current production
path is an energy-only reference-SCF PCM/ddX workflow for RHF/ROHF references.

## Background

PCM replaces explicit solvent molecules with a dielectric continuum surrounding
the solute cavity. In OpenQP, the production path uses ddX for
domain-decomposition continuum models and currently couples that reaction field
to the reference SCF energy. This page documents the implemented input
contract; state-specific excited-state PCM, PCM gradients, and PCM
optimizations are future workflow extensions. See
[References](../references.md#pcm-and-ddx) for PCM and domain-decomposition
ddPCM literature.

## Minimal Example

`.oqp`:

```text
hf/6-31g* pcm(epsilon=78.3553)
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("pcm_keywords")
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.hf(basis="6-31g*")
job.workflow.pcm(epsilon=78.3553)
```

Legacy `.inp`:

```ini
[pcm]
enabled=true
backend=ddx
mode=reference_scf
model=ddpcm
epsilon=78.3553
```

## Keywords

### `enabled`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | PCM runtime activation |

Enables the PCM reaction-field contribution. When true, the input checker
requires the current production scope: `backend=ddx`, `mode=reference_scf`,
`runtype=energy`, and RHF/ROHF.

### `backend`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `ddx` |
| Values | `ddx`, `pcmsolver` |
| Used by | PCM backend selection |

Selects the PCM backend. Only `ddx` has an implemented runtime energy path in
the current production scope.

### `mode`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `reference_scf` |
| Values | `reference_scf`, `reference_scf_plus_post_state`, `post_state_correction` |
| Used by | PCM coupling mode |

Selects the coupling mode. Only `reference_scf` is currently implemented.

### `model`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `ddpcm` |
| Values | `ddcosmo`, `ddpcm`, `ddlpb`, `iefpcm`, `cpcm` |
| Used by | continuum model |

Selects the continuum model. With `backend=ddx`, use `ddcosmo`, `ddpcm`, or
`ddlpb`.

### `solvent`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `water` |
| Used by | solvent label |

Solvent label for readability and future solvent-model support.

### `epsilon`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `78.3553` |
| Used by | dielectric constant |

Static dielectric constant. Must be numeric and greater than `1`.

### `radii`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `uff` |
| Used by | cavity radii model |

Selects the cavity radii model.
