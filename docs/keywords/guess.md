# `[guess]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `type` | `huckel` | `huckel`, `modhuckel`, `hcore`, `json`, `auto`, `sap`, or `minao`. |
| `file` | empty | JSON restart file for `type=json` or `type=auto`. |
| `file2` | empty | Previous-step restart file for NACME-style workflows. |
| `save_mol` | `False` | Save molecular data for restart/debug workflows. |
| `continue_geom` | `False` | Reuse geometry from a JSON restart. Only applies with `type=json`. |
| `swapmo` | empty | Optional MO swap list. |

## Restart Behavior

Use `type=json` with `file=<restart.json>` when a calculation must restart from
stored orbitals or molecular data. Use `type=auto` when OpenQP should read the
restart file if it exists and otherwise fall back to a Huckel guess.

For NACME, provide either `[guess] file2` or `[input] system2` so the overlap can
use previous-step information.
