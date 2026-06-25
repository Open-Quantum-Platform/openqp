# `[tdhf]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `type` | `rpa` | `rpa`, `tda`, `sf`, `mrsf`, or `umrsf`. |
| `nstate` | `1` | Number of response states. |
| `target` | `1` | Target response state for selected workflows. |
| `multiplicity` | `1` | Response-state multiplicity. |
| `maxit` | `50` | Response solver maximum iterations. |
| `maxit_zv` | `50` | Z-vector maximum iterations. |
| `conv` | `1.0e-6` | Response convergence threshold. |
| `zvconv` | `1.0e-6` | Z-vector convergence threshold. |
| `nvdav` | `50` | Davidson subspace dimension. |
| `z_solver` | `0` | Z-vector solver selector. |
| `gmres_dim` | `50` | GMRES dimension when used. |

## Types

- `rpa` and `tda`: ordinary TDHF/TDDFT.
- `sf`: spin-flip TDDFT.
- `mrsf`: mixed-reference spin-flip TDDFT.
- `umrsf`: UHF-reference MRSF energy path.

MRSF-EKT should use `[input] runtype=ekt`, `[tdhf] type=mrsf`, and the `[ekt]`
section rather than legacy `mrsf_ekt_ip` or `mrsf_ekt_ea` TDHF types.
