# `[neb]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `product` | empty | Product endpoint XYZ file. Required for NEB. |
| `nimage` | `5` | Number of NEB images. |

NEB is wired through the native `oqp` optimizer and the geomeTRIC backend. The
reactant comes from `[input] system`; the product endpoint comes from
`[neb] product`.
