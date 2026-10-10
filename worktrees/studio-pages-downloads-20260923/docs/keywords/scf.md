# `[scf]`

The `[scf]` section controls the reference wavefunction and SCF convergence
strategy. It is required for HF, DFT, MP2, TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT,
SOC, NACME, and EKT workflows.

## Minimal Examples

Closed-shell singlet:

```ini
[scf]
type=rhf
multiplicity=1
```

Triplet ROHF reference for MRSF-TDDFT:

```ini
[scf]
type=rohf
multiplicity=3
```

## Keywords

### `type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `rhf` |
| Values | `rhf`, `rohf`, `uhf` |
| Used by | reference selection |

Selects the SCF reference. Use `rhf` for closed-shell singlets, `rohf` for the
high-spin reference commonly used by SF/MRSF workflows, and `uhf` for
unrestricted references including UMRSF-TDDFT. Standalone MP2 supports RHF,
UHF, and ROHF references.

The input checker rejects `type=rhf` with `multiplicity>1`.

### `multiplicity`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | spin state and electron occupation |

Sets the spin multiplicity, `2S+1`. SF/MRSF examples usually use a triplet ROHF
reference with `multiplicity=3`.
In the high-level Python API, `job.theory.mrsf(...)` supplies this
usual MRSF reference multiplicity implicitly; set `multiplicity` in
`job.molecule(...)` for ordinary HF/DFT references or when overriding the
default deliberately.

### `maxit`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `30` |
| Used by | SCF iteration loop |

Maximum number of SCF iterations before failure or fallback handling.

### `conv`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0e-6` |
| Used by | SCF convergence test |

Energy/density convergence threshold for the SCF cycle. Tighten this for
high-precision properties or difficult state comparisons.

### `maxdiis`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `7` |
| Used by | DIIS subspace |

Maximum number of vectors kept in the DIIS extrapolation space. Set
`maxdiis <= 13` whenever an active DIIS stage uses `ediis`, `adiis`, `vdiis`,
or any nonzero `vshift`; those paths use the deterministic simplex solver,
whose exact supported history limit is 13 states. The rule applies to the
primary DIIS/`auto`/`ml` converger, to an enabled `init_scf` stage when its
effective `init_converger` is DIIS, and to any DIIS stage selected through
`escalation` or `alternative_scf`. A calculation using only `soscf` or `trah`,
or unshifted `cdiis`, may use a larger history.

### `diis_type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `cdiis` |
| Values | `none`, `cdiis`, `ediis`, `adiis`, `vdiis` |
| Used by | DIIS convergence accelerator |

Selects the DIIS variant. `cdiis` is the default. `vdiis` can be useful for
difficult cases with level shifting. Active `ediis`, `adiis`, and `vdiis`
stages require `maxdiis <= 13`.

### `diis_reset_mod`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `10` |
| Used by | DIIS reset logic |

Controls periodic DIIS reset behavior.

### `diis_reset_conv`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.005` |
| Used by | DIIS reset logic |

Convergence threshold used by the DIIS reset heuristic.

### `cdiis_switch`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.3` |
| Used by | CDIIS switching heuristic |

Threshold for switching into CDIIS-style extrapolation.

### `vdiis_vshift_switch`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.003` |
| Used by | VDIIS level-shift behavior |

Threshold for VDIIS/level-shift switching.

### `vshift`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.0` |
| Used by | virtual-orbital level shift |

Applies a level shift during SCF. Larger values can stabilize convergence but
may slow final convergence. For an active DIIS stage, any nonzero value selects
the C-DIIS/E-DIIS/C-DIIS cascade regardless of `diis_type`, so `maxdiis` must
not exceed 13.

### `incremental`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | incremental Fock build |

Enables incremental Fock builds where supported.

### `mom`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | maximum-overlap method |

Enables MOM occupation tracking for challenging electronic states.

### `mom_switch`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.003` |
| Used by | MOM activation heuristic |

Threshold controlling when MOM behavior is activated.

### `pfon`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | pseudo-fractional occupation number method |

Enables pFON smoothing for difficult SCF cases.

### `pfon_start_temp`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `2000.0` |
| Used by | pFON schedule |

Initial pFON electronic temperature.

### `pfon_cooling_rate`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `50.0` |
| Used by | pFON schedule |

Cooling rate used as pFON occupations are annealed.

### `pfon_nsmear`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `5.0` |
| Used by | pFON occupation smoothing |

Controls the smearing width/count used by the pFON method.

### `converger_type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `diis` |
| Values | `diis`, `soscf`, `trah`, `auto`, `ml` |
| Used by | SCF convergence manager |

Selects the main SCF convergence strategy. `auto` and `ml` are Python-side
manager modes that choose a concrete native converger before dispatch.

### `init_converger`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `None` |
| Values | empty, `none`, `diis`, `soscf`, `trah`, `auto`, `ml` |
| Used by | initial SCF phase |

Overrides the converger used during an optional initial SCF stage.

### `alternative_scf`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `trah` |
| Values | `diis`, `soscf`, `trah`, `auto`, `ml` |
| Used by | fallback SCF attempts |

Fallback converger after repeated SCF difficulty.

### `escalation`

| Field | Value |
| --- | --- |
| Type | comma-separated string |
| Default | empty |
| Values | `diis`, `soscf`, `trah` stages |
| Used by | custom fallback ladder |

Overrides the default escalation ladder. Example:

```ini
[scf]
escalation=soscf,trah
```

### `forced_attempt`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `2` |
| Used by | repeated SCF attempts |

Controls how many attempts the SCF manager can force before declaring failure
or escalating.

### `init_scf`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `no` |
| Values | `no`, `rhf`, `uhf`, `rohf`, `rks`, `uks`, `roks` |
| Used by | initial low-cost SCF stage |

Requests an initial SCF calculation before the main SCF. KS-style initial SCF
values (`rks`, `uks`, `roks`) require `[input] functional`.

### `init_basis`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `none` |
| Used by | initial SCF basis selection |

Basis used for the initial SCF stage. `library` can be used with
`init_library`.

### `init_library`

| Field | Value |
| --- | --- |
| Type | multiline string |
| Default | empty |
| Used by | tagged basis mapping for `init_basis=library` |

Maps atom tags to basis names for the initial SCF basis.

### `init_it`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `15` |
| Used by | initial SCF stage |

Maximum number of iterations in the initial SCF stage.

### `init_conv`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.001` |
| Used by | initial SCF stage |

Convergence threshold for the initial SCF stage.

### `save_molden`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `True` |
| Used by | Molden output |

Writes Molden orbital files when supported. Hessian calculations combine SCF
orbitals and normal modes in the frequency Molden file. MRSF-EKT calculations
also write a state-labeled Dyson Molden file after the EKT kernel completes.
See [Orbital and Vibrational Output](../workflows/orbital-vibrational-output.md).

### `rstctmo`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | MO ordering |

Restricts molecular-orbital ordering in workflows that need stable state or
orbital labels.

### `scal_rel`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` |
| Values | `0`, `1`, `2` |
| Used by | scalar-relativistic one-electron Hamiltonian |

Selects the scalar relativistic Douglas-Kroll-Hess correction applied to the
spin-free one-electron core Hamiltonian:

| Value | Meaning |
| --- | --- |
| `0` | Disabled. |
| `1` | First-order Douglas-Kroll-Hess correction. |
| `2` | First- and second-order Douglas-Kroll-Hess correction. |

This is not the SOC switch. `scal_rel` changes the scalar, spin-free
Hamiltonian; SOC is selected separately with `[input] runtype=soc` and
`[input] soc_2e`. See [References](../references.md#scalar-relativistic-correction)
for the DKH background.

The input-file default is `scal_rel=0`. In the high-level Python API,
`job.workflow.soc(...)` sets `scal_rel=2` by default for the SOC workflow
because the documented SOC path uses a DKH2 spin-free Hamiltonian before the
SOC matrix is formed. Pass `scal_rel=0`, `1`, or `2` to override that helper
default.

### `stability`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | SCF stability analysis |

Requests stability analysis where wired.

### `soscf_lvl_shift`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0` |
| Used by | SOSCF convergence |

Level shift used by the SOSCF converger.

### `verbose`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | SCF logging |

Controls SCF verbosity.

## TRAH Keywords

TRAH is the trust-region augmented Hessian SCF converger. The package build uses
the native TRAH implementation by default when TRAH is requested.

### `trh_impl`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` |
| Used by | TRAH implementation selection |

Selects the TRAH implementation. `auto` uses the available implementation and
falls back to native behavior when external OpenTRAH is absent.

### `trh_stab`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | TRAH stability mode |

Enables TRAH stability behavior.

### `trh_ls`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | TRAH line search |

Enables line search inside TRAH.

### `trh_sub_solver`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `davidson` |
| Used by | TRAH micro-iteration solver |

Selects the subsystem solver used by TRAH.

### `trh_nrtv`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | TRAH random trial vectors |

Number of random trial vectors used in the TRAH subsystem.

### `trh_r0`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.4` |
| Used by | TRAH trust radius |

Initial TRAH trust radius.

### `trh_jd_start`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `30` |
| Used by | TRAH Jacobi-Davidson switch |

Iteration threshold for the Jacobi-Davidson path in TRAH.

### `trh_nmic`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `50` |
| Used by | TRAH micro-iterations |

Maximum number of TRAH micro-iterations.

### `trh_gred`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.001` |
| Used by | TRAH global reduction |

Global reduction factor for the TRAH microproblem.

### `trh_lred`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.0001` |
| Used by | TRAH local reduction |

Local reduction factor for the TRAH microproblem.

### `xc_c2f`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` (follows `perf`) |
| Values | `on`, `off`, `auto` |
| Used by | DFT SCF — coarse-to-fine XC grid during descent |

Use a coarse XC grid during the SCF descent and the full grid in the convergence
tail. Exact at convergence. See [Performance](../performance.md).

### `xc_phi_cache`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` (follows `perf`; no preset enables it) |
| Values | `on`, `off`, `auto` |
| Used by | DFT SCF — collocation-Φ cache |

Cache the collocation Φ across SCF iterations (exact). Opt-in; on the CPU
benchmark it was net-negative, so no preset enables it.

### `xc_incdft`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `auto` (follows `perf`; no preset enables it) |
| Values | `on`, `off`, `auto` |
| Used by | DFT SCF — incremental DFT (experimental) |

Incremental XC from density differences. Experimental; typically slows SCF
convergence. Opt-in only.

### `grad_cutoff`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `auto` (follows `perf`; `1.0d-10` baseline) |
| Used by | gradients — 2e-derivative Schwarz cutoff |

Schwarz block cutoff for the 2e-derivative gradient build. Looser values trade a
small gradient error for speed (`1.0d-8` ≈ 5×10⁻⁷ a.u.).

## Progressive Screening Keywords

Progressive screening is an opt-in SCF acceleration path. During the early
descent it may use looser integral and XC-grid thresholds; once the DIIS error
falls below `pscreen_tight`, it returns to the normal converged thresholds. No
performance preset enables it automatically.

### `pscreen`, `pscreen_k`, `pscreen_cap`, `pscreen_tight`

| Keyword | Type | Default | Meaning |
| --- | --- | --- | --- |
| `pscreen` | boolean | `False` | Enable progressive integral/grid screening. |
| `pscreen_k` | float | `1.0e-2` | Safety factor in the iteration-dependent cutoff, approximately `pscreen_k * DIIS error`. |
| `pscreen_cap` | float | `1.0e-8` | Loosest permitted integral cutoff. Values looser than `1e-8` can risk changing the SCF state. |
| `pscreen_tight` | float | `1.0e-4` | DIIS-error threshold for returning to tight cutoffs; the effective value is at least `10 * scf.conv`. |

### `pscreen_xc_dcut`, `pscreen_xc_aocut`, `pscreen_grid_rad`, `pscreen_grid_ang`

| Keyword | Type | Default | Meaning |
| --- | --- | --- | --- |
| `pscreen_xc_dcut` | float | `0.0` | Loose XC density cutoff during descent; zero disables this ramp. |
| `pscreen_xc_aocut` | float | `0.0` | Loose XC AO-pruning threshold during descent; zero disables it. |
| `pscreen_grid_rad` | integer | `0` | Coarse radial-grid size during descent; both grid sizes must be positive to enable this progressive-screening override. |
| `pscreen_grid_ang` | integer | `0` | Coarse angular-grid size during descent; zero disables this override, independent of the separate `xc_c2f` schedule. |

Keep `pscreen_tight` safely above the final SCF convergence threshold so the
calculation reaches the normal integral and grid accuracy before convergence.
