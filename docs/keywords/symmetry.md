# `[symmetry]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `enabled` | `false` | `false`, `true`, or `auto`. |
| `point_group` | `auto` | Requested point group or automatic detection. |
| `subgroup` | `auto` | Requested Abelian subgroup or automatic choice. |
| `label_mo` | `True` | Label molecular orbitals. |
| `label_states` | `True` | Label response states. |
| `label_modes` | `True` | Label vibrational modes. |
| `use_integral_symmetry` | `False` | Enable integral symmetry reduction when validated for the workflow. |
| `use_response_symmetry` | `False` | Enable response symmetry reduction when validated for the workflow. |
| `tolerance` | `1.0e-5` | Geometry tolerance for symmetry detection. |
| `strict` | `False` | Strict requested/detected-group matching. |

Symmetry metadata can be useful even when reductions are disabled. Keep
reduction flags off unless the specific workflow has been validated with them.
