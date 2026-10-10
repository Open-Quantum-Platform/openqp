# `[neb]`

The `[neb]` section defines the product endpoint and image count for nudged
elastic band calculations. Use it with `[input] runtype=neb`. All native NEB
algorithm controls live in [`[oqp]`](oqp.md); concise `.oqp` input routes them
there automatically.

## Keywords

### `product`

| Field | Value |
| --- | --- |
| Type | string |
| Default | empty |
| Used by | NEB product endpoint |

Path to the product endpoint XYZ file. The reactant geometry comes from
`[input] system`; the product geometry comes from `product`.

Example:

```ini
[neb]
product=HCN_RHF-DFT_NEB_OQP_product.xyz
nimage=7
```

### `nimage`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `5` |
| Used by | NEB image count |

Number of NEB images. The input checker requires at least `3`.

For concise input, write the endpoint and native controls together:

```text
dft/pbe0/def2-svp neb(S0,product="product.xyz",images=7,output="path.xyz")
geom="reactant.xyz"
```

The concise `images` spelling lowers to `nimage`, and `output` lowers to
[`[oqp] neb_output`](oqp.md#neb_output).

The remaining six keys are retained only for traditional `.inp` geomeTRIC NEB
compatibility. Concise `.oqp` rejects them and uses native `spring`, `fmax`,
`frms`, boolean `climb`, `align`, and `opt_ends` instead. The native replacements
have different semantics or units, so old values should not be copied
mechanically.

### `k`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0` |
| Used by | legacy geomeTRIC NEB spring strength |

### `maxg`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.1` |
| Used by | legacy geomeTRIC NEB maximum-gradient convergence |

### `avgg`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.05` |
| Used by | legacy geomeTRIC NEB average-gradient convergence |

### `climb`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.5` |
| Used by | legacy geomeTRIC NEB climbing threshold |

This float is distinct from the boolean [`[oqp] climb`](oqp.md#climb) used by
native NEB.

### `align`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | legacy geomeTRIC NEB image alignment |

### `optep`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | legacy geomeTRIC NEB endpoint optimization |
