# `[tdhf]`

The `[tdhf]` section controls TDHF, TDDFT, spin-flip TDDFT, MRSF-TDDFT, and
UMRSF-TDDFT response calculations. Use `[input] method=tdhf` to activate these
workflows.

## Background

MRSF-TDDFT is OpenQP's main multistate response method. It starts from an
open-shell high-spin reference, builds spin-flip response spaces, and mixes the
reference density information so target states are less affected by ordinary
spin-flip spin contamination. In practice this makes the same response
machinery useful for multiconfigurational ground-state surfaces, excited-state
surfaces, conical-intersection work, gradients, NACME, SOC, and MRSF-EKT. See
[References](../references.md#mrsf-tddft) for the original theory papers and
recent overview articles.

## Minimal MRSF-TDDFT Example

`.oqp`:

```text
mrsf(nstate=5)/bhhlyp/6-31g*
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("mrsf_keywords")
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=5)
```

Legacy `.inp`:

```ini
[input]
method=tdhf

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=5
```

## Keywords

### `type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `rpa` |
| Values | `rpa`, `tda`, `sf`, `mrsf`, `umrsf`, `qmrsf_dk`, `mrsf_ekt_ip`, `mrsf_ekt_ea` |
| Used by | response model selection |

Selects the response model. Use `mrsf` for production MRSF-TDDFT workflows.
`sf` selects ordinary spin-flip TDDFT, and `umrsf` selects the unrestricted
MRSF energy path. `qmrsf_dk` selects the quintet-reference dressed-kernel
method described in [QMRSF-DK](../workflows/qmrsf-dk.md). The legacy
`mrsf_ekt_ip` and `mrsf_ekt_ea` values are
energy-only; the current EKT workflow should use `[input] runtype=ekt`,
`[tdhf] type=mrsf`, and the `[ekt]` section.

MRSF and SF workflows require an ROHF reference in the current code path.
UMRSF-TDDFT requires a UHF reference. QMRSF-DK requires a quintet
(`[scf] multiplicity=5`) ROHF/ROKS reference and `[input] runtype=energy`.

### `nstate`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | number of response roots |

Number of excited states to compute. It must be at least as large as the highest
state requested by gradients, optimizations, NACME, SOC, Hessians, or EKT.

### `nstate_s`, `nstate_t`

| Field | Value |
| --- | --- |
| Type | integer |
| Defaults | `0`, `0` |
| Used by | unequal singlet/triplet SOC spaces |

Optional numbers of singlet and triplet roots for SOC. Zero means to use
`nstate` for that common count. In a `.oqp` file, do not set these
internal selectors directly; write `soc(ns=3,nt=5)`. Both counts must be
provided together, and that form must not be combined with route `nstate`.

### `target`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | target-state workflows |

Target response state for workflows that read a single TDHF/MRSF state.

### `multiplicity`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | response-state spin selection |

Requested response-state multiplicity. For MRSF-TDDFT, this is the target spin
multiplicity after spin flip, not necessarily the same as the high-spin ROHF
reference multiplicity.
For SOC, do not set this as a single target multiplicity; the SOC workflow
computes singlet and triplet response roots internally.

### `maxit`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `50` |
| Used by | Davidson/response solver |

Maximum number of response-solver iterations.

### `conv`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0e-6` |
| Used by | response convergence |

Convergence threshold for the response solver.

### `nvdav`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `50` |
| Used by | Davidson subspace |

Maximum Davidson subspace dimension. The input checker warns when `nvdav` is
smaller than `nstate`.

### `maxit_zv`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `50` |
| Used by | Z-vector solver |

Maximum number of Z-vector iterations for gradient/property workflows.

### `zvconv`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0e-6` |
| Used by | Z-vector convergence |

Convergence threshold for Z-vector equations.

### `z_solver`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` |
| Values | `0`, `1`, `2`, `3` |
| Used by | Z-vector linear solver |

Selects the Z-vector solver. In the source comments, `0` is CG, `1` is legacy
GMRES, `2` is MINRES, and `3` is AUTO.

### `gmres_dim`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `50` |
| Used by | GMRES Z-vector solver |

Subspace dimension for GMRES when that solver is selected.

### `tlf`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `2` |
| Used by | MRSF state-overlap minor determinants |

Selects how the MRSF state overlap between consecutive geometries,
`<Psi_I(t-dt)|Psi_J(t)>`, is evaluated. The reference determinant is shared, so
the overlap factorizes into a contraction of the response amplitudes with three
classes of minor determinants of the MO overlap matrix: `s_ij` one-hole
occupied minors, `s_ab` particle minors, and `s_ia` mixed minors. `tlf` selects
the treatment of `s_ij` and `s_ab`; `s_ia` is always exact.

| `tlf` | Minors | Notes |
| --- | --- | --- |
| `0` (`notlf`, `exact`) | Exact Gaussian-elimination minors, no truncation | Invariant to orbital rotations between steps. This is *not* the zeroth-order TLF(0) of the paper, which is not implemented. |
| `1` | First-order truncated Leibniz formula, TLF(1) | JCTC **15**, 882 (2019) |
| `2` | Second-order truncated Leibniz formula, TLF(2) | Most accurate TLF approximation; KNU-GAMESS `ndtlf=2` |

The truncated Leibniz formula assumes the MOs of consecutive steps are nearly
orthonormal, i.e. that the MO overlap matrix is close to diagonal. When
near-degenerate doubly occupied orbitals rotate into each other within one
nuclear step -- a 45-degree mixing of two occupied orbitals has been observed in
hot uracil trajectories -- the diagonal MO overlaps fall to about 0.7 and TLF(2)
returns a collapsed state overlap (diagonal elements around 0.3-0.4) even though
the SCF solution and the MRSF surfaces are continuous. Norm-preserving
interpolation then turns that collapse into a large spurious time-derivative
coupling.

The exact minors are invariant to such rotations, and for molecules the size of
uracil (30 occupied alpha orbitals, 6-31G*) they cost the same wall time as
TLF(2). For large systems where the `nvir^2` particle minors dominate, the
recommended route is Jacobi's complementary-minor identity -- all one- and
two-hole minors from one LU factorization of the occupied block -- rather than
truncation.

The state-overlap section of the log states which evaluation was used
(`state-overlap minors: exact minor determinants (tlf=0, default; ...)` or
`TLF(n) truncated-Leibniz minors; ...`). NAMD additionally warns when every
column norm of the retained state overlap falls below 0.5.

### `hfscale`, `cam_alpha`, `cam_beta`, `cam_mu`

| Field | Value |
| --- | --- |
| Type | float |
| Defaults | `-1.0`, `-1.0`, `-1.0`, `-1.0` |
| Used by | response functional parameter overrides |

Override exact-exchange and CAM/range-separated parameters for response
calculations. Negative values mean use the selected functional defaults.

For `type=qmrsf_dk` these set the exchange carried by the dressed kernel,
independently of the reference. With a global hybrid the kernel uses
`hfscale`; with a range-separated reference it uses
`cam_alpha*K + cam_beta*K(erf(cam_mu*r)/r)`, and any of the three that is left
negative is inherited from `[dftgrid]`.

### `spc_coco`, `spc_ovov`, `spc_coov`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `-1.0` |
| Used by | spin-purification correction parameters |

Advanced MRSF/SF response parameters. Leave negative unless following a specific
validated protocol.

### `conf_threshold`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `5.0e-2` |
| Used by | configuration analysis/output |

Threshold for reporting or using response configurations.

### `ixcore`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `-1` |
| Used by | core-level response workflows |

Core-orbital selector for core-level response calculations such as X-ray
absorption workflows.

### `resp_cutoff`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `auto` (follows `perf`; `1e-8` baseline) |
| Used by | MRSF-TDDFT — response 2e-integral cutoff |

2e-integral cutoff for the MRSF response build. `1e-8` is exact to ≪ µeV;
looser values (e.g. `1e-6`) trade a few µeV for speed. Never tighter than the
SCF integral cutoff. See [Performance](../performance.md).

### `fp32`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` (follows `perf`; no preset enables it) |
| Values | `on`, `off`, `auto` |
| Used by | MRSF-TDDFT — single-precision response digestion |

Single-precision MRSF response Fock digestion. Non-reproducible and can flip
near-degenerate states; net-slower than FP64 on CPU. Opt-in only.

### `zv_warmstart`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` (follows `perf`; on at `perf` ≥ 1) |
| Values | `on`, `off`, `auto` |
| Used by | MRSF-TDDFT gradients — z-vector (CPHF) warm-start |

Reuse the previous geometry step's z-vector as the CG/GMRES initial guess
(exact; the solve still converges to the same tolerance). Most effective across
many nearby geometries (optimization, MD).
