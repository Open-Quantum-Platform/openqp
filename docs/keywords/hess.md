# `[hess]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `type` | `numerical` | `numerical` or `analytical`. |
| `state` | `0` | Hessian target state. HF/DFT uses `0`; TDHF workflows use excited-state indices. |
| `dx` | `0.01` | Finite-difference step for numerical Hessians. |
| `nproc` | `1` | Worker count for numerical Hessians. |
| `read` | `False` | Read existing Hessian data where supported. |
| `restart` | `False` | Continue a numerical Hessian workflow. |
| `temperature` | `298.15` | Thermochemistry temperature list in Kelvin. |
| `clean` | `False` | Clean temporary Hessian files. |

## Analytical Hessians

`type=analytical` is supported for validated HF/DFT ground-state cases. The
input checker will request `type=numerical` when the method, state, or basis
falls outside the current analytical path.
