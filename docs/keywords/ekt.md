# `[ekt]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `ip` | `True` | Run the ionization-potential channel. |
| `ea` | `False` | Run the electron-affinity channel. |

Use `[input] runtype=ekt`, `[input] method=tdhf`, `[tdhf] type=mrsf`, and at
least one of `ip=true` or `ea=true`.

See [MRSF-EKT](../workflows/ekt.md).
