# `[json]`

The `[json]` section carries advanced restart metadata used by internal JSON
interfaces. Most users should configure the calculation through the normal
route and [`[guess]`](guess.md) section instead.

### `scf_type`, `basis`, `library`, `do_init`

| Keyword | Type | Default | Meaning |
| --- | --- | --- | --- |
| `scf_type` | string | empty | Saved/imported SCF reference metadata used by a JSON guess. |
| `basis` | string | empty | Saved/imported basis metadata used by a JSON guess. |
| `library` | string | empty | Saved/imported basis-library mapping used by a JSON guess. |
| `do_init` | string | `no` | Internal `no`/`yes` marker set when a JSON basis/library mismatch requires updating `scf.init_*`; not a normal user switch. |

These keys are compatibility metadata, not an alternative to the `.oqp`
method, basis, and reference syntax.
