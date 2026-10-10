# `[mp2]`

The `[mp2]` section controls spin-component scaling for standalone
`[input] method=mp2` calculations. Ordinary MP2 does not require this section;
the defaults are conventional unscaled MP2.

MP2 requires an HF reference, so `[input] functional` must be empty. RHF
references support analytic nuclear gradients and gradient-driven geometry
calculations; UHF and ROHF references remain energy-only. In concise input,
select the reference with
`mp2(reference=rhf|rohf|uhf)/BASIS`; the value lowers to `[scf] type` and is
not a keyword in this section.

## Minimal Examples

Conventional MP2 in `.oqp`:

```text
mp2/6-31g
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mp2")
job.molecule(geometry="water")
job.theory.mp2(basis="6-31g")
mol = job.run()
```

Analytic RHF-MP2 gradient in `.oqp`:

```text
mp2/6-31g grad(S0)
geom="h2o.xyz"
```

The spin-scaling preset is applied to both the energy and its analytic
gradient. For example, replace the first line with
`mp2(variant=scs-mp2)/6-31g` for SCS-MP2.

Legacy `.inp`:

```ini
[input]
method=mp2
runtype=energy
functional=
```

For spin-component-scaled MP2, add an exact `.oqp` section call:

```text
mp2/6-31g mp2(variant=scs-mp2)
geom="h2o.xyz"
```

The legacy section is:

```ini
[mp2]
variant=scs-mp2
```

Custom scaling:

```ini
[mp2]
variant=custom
same_spin_scale=0.50
opposite_spin_scale=1.10
```

## Keywords

### `variant`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `mp2` |
| Values | `mp2`, `conventional`, `scs-mp2`, `scs`, `sos-mp2`, `sos`, `os-mp2`, `os`, `opposite-spin`, `ss-mp2`, `ss`, `same-spin`, `sss-mp2`, `sss`, `scs-mi-mp2`, `scs-mi`, `custom` |
| Used by | standalone MP2 energy and RHF analytic gradient |

Selects the MP2 spin-scaling preset. The scales multiply same-spin and
opposite-spin correlation components:

```text
E(MP2) = c_ss * (E_aa + E_bb) + c_os * E_ab
```

| Variant | Same-spin scale `c_ss` | Opposite-spin scale `c_os` |
| --- | --- | --- |
| `mp2`, `conventional` | `1.0` | `1.0` |
| `scs-mp2`, `scs` | `1.0 / 3.0` | `1.2` |
| `sos-mp2`, `sos` | `0.0` | `1.3` |
| `os-mp2`, `os`, `opposite-spin` | `0.0` | `1.0` |
| `ss-mp2`, `ss`, `same-spin`, `sss-mp2`, `sss` | `1.0` | `0.0` |
| `scs-mi-mp2`, `scs-mi` | `1.29` | `0.40` |
| `custom` | `same_spin_scale` | `opposite_spin_scale` |

Named presets set both scale factors internally. If you provide
`same_spin_scale` or `opposite_spin_scale` with a named preset, those explicit
numbers are not used; choose `variant=custom` instead.

### `same_spin_scale`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0` |
| Used by | `variant=custom` |

Scale factor for the same-spin MP2 contribution, `E_aa + E_bb`, when
`variant=custom`. The input checker requires finite numeric values for custom
scales.

### `opposite_spin_scale`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0` |
| Used by | `variant=custom` |

Scale factor for the opposite-spin MP2 contribution, `E_ab`, when
`variant=custom`. The input checker requires finite numeric values for custom
scales.

## Notes

- RHF MP2 accepts `runtype=energy`, `grad`, `optimize`, `ts`, `mep`, and `irc`.
- UHF and ROHF MP2 accept `runtype=energy` only; derivative requests stop
  before the electronic-structure calculation.
- `method=mp2` rejects non-empty `[input] functional` values.
- The reference is selected through `[scf] type`; RHF, UHF, and ROHF references
  are supported for energies.
- SOS and OS-only variants set the same-spin scale to zero and skip same-spin
  work.

See the [MP2 workflow](../workflows/mp2.md) for complete input and Python
examples.
