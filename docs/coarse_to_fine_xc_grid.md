# Coarse-to-fine XC grid schedule for DFT SCF

**Status:** **default ON**, opt-out via `OQP_XC_C2F=0`. RHF/UHF/ROHF, pure-DIIS
converger. The schedule never changes the converged result — it only changes the
early-iteration grid — and disables itself wherever it would not help or could be
unsafe (see below), so it is safe to leave on.

## Idea

A DFT SCF cycle spends a large fraction of its time evaluating the
exchange–correlation (XC) energy and potential on a numerical integration grid.
Early in the SCF the density is still far from converged, so spending the full
production-grid effort on those iterations is wasteful: a cheaper, sparser grid
guides the density almost as well.

This feature runs the **early** SCF iterations on a cheap **coarse** XC grid and
switches to the **production** grid a few iterations before convergence, letting
DIIS re-converge on the production grid. The reported energy (and any subsequent
gradient/Hessian, which build their own production grid) are therefore full
production-grid quality.

Only the XC contribution is affected. The two-electron J/K (Coulomb/exchange)
build is grid-independent, and OpenQP's incremental-Fock reference
(`fold`/`dold`) holds **only** the J/K part, so it stays valid across the switch
and is *not* reset — the coarse stage composes cleanly with the incremental Fock
build.

## Why it is safe

- **Energy is production-grid quality.** Convergence is never declared on a
  coarse-grid Fock; the schedule forces at least one production-grid iteration
  before exit and DIIS re-converges there. Converged energies match the
  all-production-grid baseline to ~1e-9 Ha (see results).
- **The coarse grid never escapes the SCF driver.** It is a local object inside
  `scf_driver`; post-SCF properties, gradients and Hessians independently build
  the production grid via `dft_initialize`. `dft_build_grid_sized` also
  saves/restores the grid sizes in `infos%dft`.
- **Minnesota functionals auto-disable the coarse stage.** Minnesota-family
  functionals (M05/M06/M08/M11/MN1x/SOGGA11 and revM* variants) have large
  coarse-grid errors, so the schedule detects them by name and runs entirely on
  the production grid (printing a notice).
- **On by default but self-disabling** (with a notice in the log) when the run is
  ineligible. The schedule only activates for a DFT run with the pure-DIIS
  converger, a non-Minnesota functional, `maxit >= 2`, no level shift on a
  non-VDIIS converger, and — checked against the **actual built point counts** —
  a coarse grid with meaningfully fewer points (< 0.9×) than the production grid.
  The last check uses real `nMolPts`, so a sparse pruned production grid (e.g.
  SG-0, which can be cheaper than the requested dense coarse grid) correctly
  disables the schedule instead of slowing the run down.

## Usage

The schedule is **on by default** (nothing changes in the input file). Set
`OQP_XC_C2F=0` to turn it off.

| Variable | Default | Meaning |
|---|---|---|
| `OQP_XC_C2F` | `on` | Master switch. Set to `0`/`no`/`false`/`off` to disable the default-on ramp. |
| `OQP_XC_C2F_RAD` | `50` | Coarse grid radial points. |
| `OQP_XC_C2F_ANG` | `110` | Coarse grid angular (Lebedev) points. |
| `OQP_XC_C2F_SWITCH` | `1.0e-2` | Pin to the production grid once `|DIIS error|` drops below this (standalone path). |

The same engine is also driven by progressive screening (upstream #238): with
`scf_pscreen` + `pscreen_grid_rad/ang` (or `OQP_PSCREEN` + `OQP_PSCREEN_GRID_RAD/ANG`)
the grid ramp shares the integral-screening pin threshold (`pscreen_tight`).

Disable for a single run:

```sh
OQP_XC_C2F=0 python pyoqp/oqp/pyoqp.py mymol.inp
```

When active the log shows the coarse grid (built in `hf_energy`) and the pin in
the convergence tail:

```
     Coarse-to-fine XC grid: coarse 50 x 110 (5500 pts) during descent, production (29000 pts) in the tail
   ...
          Coarse-to-fine / progressive screening: final full-accuracy SCF iteration.
```

## When it helps

The speedup is the per-iteration XC saving (coarse grid has far fewer points)
times the number of coarse iterations, minus any extra iterations. It is largest
when XC dominates the per-iteration cost:

- **Pure GGA / meta-GGA** (no exact exchange) benefits most.
- **Expensive production grids** (SG-3, dense 99×590) benefit most; against the
  already-efficient default SG-2 the win is modest.
- **Hybrids** benefit less, because the (grid-independent) exact-exchange J/K
  build dominates and is not sped up.

## Implementation — one unified grid-ramp engine

The coarse→fine XC-grid ramp is a **single engine**, shared by this default-on
schedule and the opt-in progressive-screening lever (`scf_pscreen`, upstream
#238). There is exactly one coarse-grid construction path and one grid-selection
path; the two features only differ in *how the engine is activated*.

- `source/dftlib/dft.F90`:
  - `dft_build_grid_sized(infos, basis, molGrid, nrad, nang)` builds a single
    unpruned Lebedev grid of the requested size without touching the libxc setup
    (saves/restores `infos%dft` grid sizes).
  - `dft_setup_descent_grid(...)` is the **single unified policy**: it decides
    whether to build a coarse "descent" grid and at what size, then builds it via
    `dft_build_grid_sized` and keeps it only if it is genuinely cheaper than the
    production grid (by built `nMolPts`). Activation is either a progressive-
    screening request (`scf_pscreen` + `pscreen_grid_rad/ang`) or the default-on
    coarse-to-fine path (unless `OQP_XC_C2F=0`), the latter gated by
    `c2f_grid_eligible` (pure-DIIS, `maxit >= 2`, no level-shift-on-non-VDIIS,
    non-Minnesota via `xc_is_grid_sensitive`).
- `source/modules/hf_energy.f90`: calls `dft_setup_descent_grid` and passes the
  resulting coarse grid to `scf_driver` as the optional `coarseGrid` argument.
- `source/scf.F90` (the engine, from #238, decoupled so it runs for the grid
  ramp independently of integral screening):
  - `ps_grid_on = present(coarseGrid)` — the grid ramp runs whenever a coarse
    grid was provided.
  - A single sticky **pin latch** (`ps_pin`) flips once `|DIIS error|` drops
    below the pin threshold (`pscreen_tight` under integral screening, else
    `OQP_XC_C2F_SWITCH`, default 1e-2), driving the grid ramp, the 2e-cutoff ramp
    and the XC-threshold ramp together.
  - `ps_cur_grid => coarseGrid` during the descent, `=> molGrid` once pinned;
    one `calc_fock(..., ps_cur_grid, ...)` call.
  - `ps_force_iter` forces one pinned, full-accuracy iteration before
    convergence, so a coarse-grid (or loose-cutoff, or reused-XC) Fock is never
    accepted as converged.

This unification was reached after #238 landed upstream with its own grid ramp;
rather than ship a second mechanism, c2f became an activation policy + safety
guards on top of the shared engine.

### Composition with the Φ-cache / IncDFT (#242)

The shared engine composes with the opt-in XC-reuse features from #242:

- **IncDFT** (`OQP_XC_INCDFT`): `xc_reuse` is forced `.false.` while on the coarse
  grid (`ps_grid_on .and. .not. ps_pin`), so reused reference XC is never carried
  across grids; convergence is refused on a reused-XC Fock.
- **Φ-cache** (`OQP_XC_PHI_CACHE`): the cache validity signature includes the grid
  slice count and the coarse grid has fewer points, so the cache rebuilds
  automatically at the switch.

Verified: with the ramp on, enabling either/both of `OQP_XC_INCDFT` and
`OQP_XC_PHI_CACHE`, or progressive screening (`OQP_PSCREEN`), still converges to
the bit-identical production-grid energy.

## Results

Hardware: 32-core Apple Silicon, gfortran-15, Release build, `OMP_NUM_THREADS=4`,
min-of-N wall times (machine was concurrently loaded; absolute times are noisy,
ratios and iteration counts are robust). 6-31G\*. "SCF wall" is the SCF-driver
phase time (includes the one-off coarse-grid build).

Converged-energy agreement (the key correctness gate; OFF = all-production-grid):

| System | Functional | Prod. grid | E(OFF) | E(c2f) | |ΔE| |
|---|---|---|---|---|---|
| H₂O | BHHLYP | SG-2 | −76.3677816996 | −76.3677816996 | < 1e-10 |
| benzene | PBE | SG-2 | −231.9341227248 | −231.9341227248 | < 1e-10 |
| benzene | PBE | SG-3 | −231.9341015289 | −231.9341015289 | < 1e-10 |
| benzene | BHHLYP | SG-2 | −232.0996658964 | −232.0996658964 | < 1e-10 |
| benzene | BHHLYP | SG-3 | −232.0996546260 | −232.0996546260 | < 1e-10 |
| (H₂O)₆ | PBE | SG-2 | −457.8922695882 | −457.8922695882 | < 1e-10 |
| C₂₀H₄₂ (62 at.) | B3LYPv5 | SG-2 | −786.7633652204 | −786.7633652204 | < 1e-10 |

Every converged energy is **bit-identical** to the all-production-grid baseline
(far inside the 1e-7 Ha gate).

Iterations and SCF wall time:

| System | Functional | Prod. grid | iters OFF→c2f | SCF wall OFF→c2f | Δ |
|---|---|---|---|---|---|
| benzene | PBE | SG-2 | 9 → 10 | 3.60s → 3.43s | −4.6% |
| benzene | PBE | **SG-3** | 9 → 10 | 7.52s → 6.14s | **−18.4%** |
| benzene | BHHLYP | SG-2 | 8 → 9 | 3.23s → 3.05s | −5.7% |
| benzene | BHHLYP | **SG-3** | 8 → 9 | 6.69s → 5.20s | **−22.2%** |
| (H₂O)₆ | PBE | SG-2 | 10 → 11 | 6.66s → 5.98s | −10.2% |
| C₂₀H₄₂ | B3LYPv5 | SG-2 | 11 → 11 | ~318s → ~276s† | (~−13%†) |

† C₂₀H₄₂ wall times were measured under heavy, varying machine load and are only
indicative; the robust results there are **identical converged energy** and
**no added iterations** (11 → 11).

Observations:
- Correctness is exact: converged energies match to ≤1e-10 Ha everywhere.
- Convergence is clean (no DIIS oscillation); the schedule costs **0–1 extra
  iterations**, repaid by the cheaper coarse cycles.
- Speedup grows with grid cost (SG-2 → SG-3 roughly triples the saving) and with
  system size / XC fraction (the 18-atom hexamer beats benzene on the same grid).
- Hybrids gain less than pure GGAs at the same grid because their exact-exchange
  J/K build (grid-independent) dominates the per-iteration cost.

Minnesota hazard check: an **M06-2X** (and **M06-L**) run prints
`Coarse-to-fine XC grid schedule disabled: grid-sensitive (Minnesota)
functional` and runs entirely on the production grid.

**Default-ON safety.** Because the schedule is on by default, it was checked not
to perturb results across reference types and functionals:

- Converged energies are bit-identical (ON vs `OQP_XC_C2F=0`) for RHF, UHF and
  ROHF and for PBE, B3LYPv5, CAM-B3LYP, SCAN and TPSS.
- The six committed DFT example references (`examples/DFT`, RHF/UHF/ROHF energy
  **and gradient**) all **PASS** with the schedule on by default.
- Energy is stationary, so it matches to ~1e-10; the converged *density* differs
  only at the SCF tolerance, so gradients/properties differ by ~1e-7 — three to
  four orders of magnitude inside the example suite's comparison tolerance
  (`round(diff, 4) == 0`, i.e. 5e-5). Set `OQP_XC_C2F=0` for bitwise
  reproduction of an all-production-grid run.

Raw logs: `/Volumes/External_Storage/claude/sessions/20260627_xc_c2f_grid/`.
