# `[oqp]`

The `[oqp]` section controls the native OpenQP optimizer selected by:

```ini
[optimize]
lib=oqp
```

| Keyword | Default | Notes |
| --- | --- | --- |
| `coordsys` | `tric` | Coordinate system; `tric` is the usual default. |
| `trust` | `0.2` | Initial trust radius. |
| `trust_max` | `0.5` | Maximum trust radius. |
| `follow` | `0` | Mode-following selector for transition-state style steps. |
| `spring` | `0.05` | NEB spring constant. |
| `climb` | `True` | Enable climbing image in NEB. |
| `fmax` | `2.0e-3` | NEB force threshold. |
| `climb_fmax` | `0.05` | Climbing-image threshold. |
| `neb_dt` | `0.5` | NEB integration step. |
| `maxmove` | `0.2` | Maximum NEB image move. |
| `opt_ends` | `True` | Optimize NEB endpoints. |
| `end_fmax` | `1.0e-3` | Endpoint force threshold. |
| `irc_step` | `0.1` | IRC step size. |
| `irc_direction` | `forward` | IRC direction. |
| `mep_step` | `0.1` | MEP step size. |

See [Optimization](../workflows/optimization.md).
