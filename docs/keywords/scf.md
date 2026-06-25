# `[scf]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `type` | `rhf` | `rhf`, `rohf`, or `uhf`. |
| `multiplicity` | `1` | Must be consistent with the reference type. |
| `maxit` | `30` | Maximum SCF iterations. |
| `maxdiis` | `7` | DIIS subspace size. |
| `diis_type` | `cdiis` | `none`, `cdiis`, `ediis`, `adiis`, or `vdiis`. |
| `conv` | `1.0e-6` | SCF convergence threshold. |
| `incremental` | `True` | Incremental Fock build. |
| `pfon` | `False` | Pseudo-fractional occupation number method. |
| `mom` | `False` | Maximum overlap method. |
| `rstctmo` | `False` | Restrict MO ordering. |
| `converger_type` | `diis` | `diis`, `soscf`, `trah`, `auto`, or `ml`. |
| `init_converger` | `None` | Optional initial converger. |
| `alternative_scf` | `trah` | Fallback converger for repeated attempts. |
| `escalation` | empty | Optional comma-separated fallback ladder such as `soscf,trah`. |
| `trh_impl` | `auto` | Native TRAH by default; explicit `otr` falls back when OpenTRAH is absent. |

## TRAH

`converger_type=trah` selects trust-region augmented Hessian SCF convergence.
The package build skips external OpenTRAH and uses the native implementation.
Source builds may enable the external OpenTRAH library, but native TRAH remains
the documented default behavior.

## Manager Modes

`converger_type=auto` and `converger_type=ml` are Python-side manager modes.
They resolve to concrete native convergers before native dispatch.
