# `[dftgrid]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `hfscale` | `-1.0` | Exact-exchange scale override. Negative means use the functional default. |
| `cam_flag` | `False` | Enable CAM/range-separated parameters explicitly. |
| `cam_alpha` | `-1.0` | CAM alpha override. Negative means use the functional default. |
| `cam_beta` | `-1.0` | CAM beta override. Negative means use the functional default. |
| `cam_mu` | `-1.0` | CAM mu override. Negative means use the functional default. |
| `rad_type` | `ta` | Radial quadrature family. |
| `rad_npts` | `96` | Radial grid size. |
| `ang_npts` | `302` | Angular grid size. |
| `partfun` | `ssf` | Molecular partition function. |
| `pruned` | `SG2` | Pruned grid selector. |
| `grid_ao_pruned` | `True` | Enable AO grid pruning. |
| `grid_ao_threshold` | `1.0e-15` | AO pruning threshold. |
| `grid_ao_sparsity_ratio` | `0.9` | AO sparsity threshold for pruning. |

Leave the hybrid and CAM overrides negative unless you are deliberately testing
a non-standard functional parameterization.
