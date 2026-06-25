# `[geometric]`

The `[geometric]` section controls the external geomeTRIC backend selected by:

```ini
[optimize]
lib=geometric
```

| Keyword | Default | Notes |
| --- | --- | --- |
| `coordsys` | `tric` | Coordinate system passed to geomeTRIC. |
| `trust` | `0.1` | Initial trust radius. |
| `tmax` | `0.3` | Maximum trust radius. |
| `convergence_set` | `GAU` | geomeTRIC convergence preset. |
| `prefix` | `geometric` | File prefix for geomeTRIC artifacts. |
| `hessian` | `never` | geomeTRIC Hessian option. |
| `irc_direction` | `forward` | IRC direction. |
| `constraints_file` | empty | Constraint file path. |
| `enforce` | `0.0` | Constraint enforcement strength. |
| `conmethod` | `0` | Constraint method selector. |

geomeTRIC is currently connected to state-specific optimization, MECI, MECP, TS,
IRC, and NEB workflows.
