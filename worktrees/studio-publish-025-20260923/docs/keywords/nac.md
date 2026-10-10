# `[nac]`

The `[nac]` section controls NAC and NACME state-pair workflows. NAC and NACME
currently require `[input] method=tdhf` and `[tdhf] type=mrsf`.

## Keywords

### `type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `numerical` |
| Used by | NAC dispatch |

Selects the NAC implementation.

### `states`

| Field | Value |
| --- | --- |
| Type | state-pair list |
| Default | `1 2` |
| Used by | NAC and NACME state pairs |

Lists pairs of excited states. Each state index must be at least `1`.

Example:

```ini
[nac]
states=1 2,2 3
```

### `dt`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1` |
| Used by | time-derivative coupling style inputs |

Time step for workflows that use time-separated geometries.

### `dx`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.0001` |
| Used by | finite-difference NAC |

Finite-difference displacement.

### `bp`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | branching-plane style workflows |

Enables branching-plane behavior where wired.

### `nproc`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | NAC workers |

Worker count. Must be at least `1`.

### `restart`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | NAC restart |

Continues a NAC workflow where supported.

### `clean`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | temporary NAC files |

Removes temporary NAC files where supported.

### `align`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `reorder` |
| Values | `reorder`, `no` |
| Used by | NACME alignment |

Controls alignment of states/geometries in the Python layer. The input checker
warns for values other than `reorder` and `no`.
