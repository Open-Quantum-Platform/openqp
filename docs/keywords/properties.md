# `[properties]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `scf_prop` | `el_mom,mulliken` | SCF properties. Add `nmr` for NMR shielding. |
| `nmr_gauge` | `cgo` | `cgo` or `giao`. |
| `td_prop` | `False` | Response-state properties. |
| `grad` | `0` | Gradient state list. |
| `nac` | empty | NAC state-pair request where used. |
| `export` | `False` | Export selected computed data. |
| `title` | empty | Output title. |
| `back_door` | `False` | Developer/debug option. |

## Gradients

`grad=0` requests the reference-state gradient for HF/DFT workflows. Response
workflows use response-state numbering and require `[tdhf] nstate` to cover the
requested state.

## NMR

Use:

```ini
[properties]
scf_prop=nmr
nmr_gauge=cgo
```

See [NMR, IR, and Raman](../workflows/nmr-ir-raman.md).
