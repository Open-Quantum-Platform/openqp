# `[ekt]`

The `[ekt]` section selects MRSF-EKT ionization-potential and electron-affinity
channels. Use it with `[input] runtype=ekt`, `[input] method=tdhf`, and
`[tdhf] type=mrsf`.

## Keywords

### `ip`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | MRSF-EKT ionization-potential channel |

Runs the IP channel.

### `ea`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | MRSF-EKT electron-affinity channel |

Runs the EA channel.

At least one of `ip` or `ea` must be true:

```ini
[input]
runtype=ekt
method=tdhf

[tdhf]
type=mrsf
nstate=10

[ekt]
ip=true
ea=false
```
