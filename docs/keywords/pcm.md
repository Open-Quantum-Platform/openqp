# `[pcm]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `enabled` | `False` | Enable the PCM reaction-field contribution. |
| `backend` | `ddx` | Runtime production path is `ddx`. |
| `mode` | `reference_scf` | Production mode for the first PCM path. |
| `model` | `ddpcm` | Continuum model. |
| `solvent` | `water` | Solvent label. |
| `epsilon` | `78.3553` | Static dielectric constant. |
| `radii` | `uff` | Cavity radii model. |

## Scope

The documented production path is:

```ini
[pcm]
enabled=true
backend=ddx
mode=reference_scf
model=ddpcm
```

Use it with `[input] runtype=energy` and RHF/ROHF references. Other backend and
mode names may be recognized for validation/planning, but they are not the
current production runtime path.
