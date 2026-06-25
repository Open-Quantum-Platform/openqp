# `[optimize]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `lib` | `oqp` | `oqp`, `geometric`, or `scipy`. |
| `optimizer` | `bfgs` | SciPy optimizer selector. |
| `maxit` | `30` | Maximum geometry steps. |
| `istate` | `1` | First target state. Use `0` for HF/DFT reference-state minima. |
| `jstate` | `2` | Second target state for crossing searches. |
| `kstate` | `3` | Third target state for TCI. |
| `imult` | `1` | First multiplicity for MECP-style workflows. |
| `jmult` | `3` | Second multiplicity for MECP-style workflows. |
| `energy_shift` | `1.0e-6` | Energy-change threshold. |
| `energy_gap` | `1.0e-5` | Crossing-gap threshold. |
| `meci_search` | `penalty` | `penalty`, `ubp`, or `hybrid`. |
| `gap_weight` | `1.0` | Gap contribution weight. |
| `init_scf` | `False` | Run initial SCF cycles during each geometry step. |

## Optimizer Backends

- `oqp`: native optimizer. Supports `optimize`, `ts`, `meci`, `mecp`, `tci`,
  `neb`, `irc`, and `mep`.
- `geometric`: external geomeTRIC backend. Supports `optimize`, `meci`, `mecp`,
  `ts`, `irc`, and `neb`.
- `scipy`: legacy lightweight backend for `optimize`, `meci`, `mecp`, and `mep`.

DL-FIND is not a current user-facing optimizer backend in this repository.
