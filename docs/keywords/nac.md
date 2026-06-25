# `[nac]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `type` | `numerical` | NAC implementation selector. |
| `dt` | `1` | Time step for time-derivative coupling style inputs. |
| `dx` | `0.0001` | Finite-difference displacement. |
| `bp` | `False` | Branching-plane style control where used. |
| `nproc` | `1` | Worker count. |
| `restart` | `False` | Continue a NAC workflow. |
| `clean` | `False` | Clean temporary NAC files. |
| `states` | `1 2` | State pairs such as `1 2,2 3`; indices are excited states. |
| `align` | `reorder` | `reorder` or `no` for NACME alignment in the Python layer. |

NAC and NACME workflows require `[input] method=tdhf` and `[tdhf] type=mrsf`.
For NACME, also provide previous-step information through `[guess] file2` or
`[input] system2`.
