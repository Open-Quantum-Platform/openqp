# `[guess]`

The `[guess]` section controls initial orbitals, JSON restart files, and manual
MO swaps. Most users can start with the default Huckel guess.

## Minimal Example

```ini
[guess]
type=huckel
```

## Keywords

### `type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `huckel` |
| Values | `huckel`, `modhuckel`, `hcore`, `json`, `auto`, `sap`, `minao` |
| Used by | initial density/orbital construction |

Selects the initial guess mode.

| Value | Meaning |
| --- | --- |
| `huckel` | Native extended-Huckel style guess. |
| `modhuckel` | Modified Huckel guess. |
| `hcore` | Core-Hamiltonian guess. |
| `json` | Read orbitals and molecular data from a JSON restart file. |
| `auto` | Read `file` if it exists; otherwise fall back to Huckel. |
| `sap` | Superposition-of-atomic-potentials guess. |
| `minao` | Projected minimal-basis density guess. |

Use `json` for exact restarts and `auto` for workflows that may or may not have
a previous restart file.

### `file`

| Field | Value |
| --- | --- |
| Type | string |
| Default | empty |
| Used by | `type=json`, `type=auto` |

Path to a JSON restart file. With `type=json`, the file must exist. With
`type=auto`, OpenQP uses it when present and otherwise computes a fresh guess.

### `file2`

| Field | Value |
| --- | --- |
| Type | string |
| Default | empty |
| Used by | NACME and previous-step workflows |

Path to a second JSON restart file. NACME can use `file2` for previous-step
information. If it is omitted, NACME can instead use `[input] system2`.

### `save_mol`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | restart, orbital export, and inspection workflows |

Requests saving molecular data for later restart or inspection. Supported
bases also add portable AO basis, molecular-orbital, vibrational, and
state-specific Dyson data for visualization. See
[Orbital and Vibrational Output](../workflows/orbital-vibrational-output.md).

### `continue_geom`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | JSON restart geometry |

When true, reuses geometry from the JSON restart file. The input checker warns
if `continue_geom=true` is set without `type=json`.

### `swapmo`

| Field | Value |
| --- | --- |
| Type | comma-separated integer list |
| Default | empty |
| Used by | manual orbital reordering before SCF |

Swaps molecular orbitals in pairs. The list must contain an even number of
indices:

```ini
[guess]
swapmo=2,7,8,21
```

This swaps `2` with `7` and `8` with `21`. Use it only when you know the target
orbital ordering.
