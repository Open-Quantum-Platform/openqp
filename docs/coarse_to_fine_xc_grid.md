# Coarse-to-fine XC grid schedule for DFT SCF

**Status:** opt-in, default OFF. Environment-variable gated. RHF/UHF/ROHF, pure-DIIS
converger only.

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
- **Off by default**, and disabled (with a notice in the log) when the run is
  ineligible. The schedule only activates for a DFT run with the pure-DIIS
  converger, a non-Minnesota functional, `maxit >= 2`, no level shift on a
  non-VDIIS converger, and — checked against the **actual built point counts** —
  a coarse grid with meaningfully fewer points (< 0.9×) than the production grid.
  The last check uses real `nMolPts`, so a sparse pruned production grid (e.g.
  SG-0, which can be cheaper than the requested dense coarse grid) correctly
  disables the schedule instead of slowing the run down.

## Usage

Enable via environment variable (nothing changes in the input file):

| Variable | Default | Meaning |
|---|---|---|
| `OQP_XC_C2F` | (off) | Master switch. Set to `1`/`yes`/`true`/`on` to enable. |
| `OQP_XC_C2F_RAD` | `50` | Coarse grid radial points. |
| `OQP_XC_C2F_ANG` | `110` | Coarse grid angular (Lebedev) points. |
| `OQP_XC_C2F_SWITCH` | `1.0e-2` | Switch to the production grid once `|DIIS error|` drops below this. |
| `OQP_XC_C2F_RESET` | (off) | If set, clear the DIIS subspace at the switch (more defensive; costs ~1 iteration). Default keeps the subspace. |

Example:

```sh
OQP_XC_C2F=1 python pyoqp/oqp/pyoqp.py mymol.inp
```

The schedule prints a banner when active, and a line at the grid switch:

```
     Coarse-to-fine XC grid schedule: ON
       coarse grid  = 50 x 110 (radial x angular)
       switch to production grid when |DIIS error| <  1.00E-02
   ...
          Switching XC grid: coarse -> production (DIIS error  3.43E-03).
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

## Implementation

- `source/dftlib/dft.F90`: `dft_build_grid_sized(infos, basis, molGrid, nrad, nang)`
  builds a single unpruned Lebedev grid of the requested size without touching
  the libxc functional setup (saves/restores `infos%dft` grid sizes).
- `source/scf.F90`:
  - `c2f_get_config` reads the environment, checks eligibility, returns the
    coarse sizes / switch threshold / reset flag.
  - `xc_is_grid_sensitive` detects Minnesota-family functionals.
  - The SCF loop builds the coarse grid once, compares `coarse_grid%nMolPts`
    against `molGrid%nMolPts` (disabling the schedule if the coarse grid is not
    cheaper), uses the coarse grid via the `on_coarse` branch of the `calc_fock`
    call, and switches to the production grid (`molgrid`) once
    `|DIIS error| < OQP_XC_C2F_SWITCH` (or, as a fail-safe, a few iterations
    before `maxit`). `coarse_this_iter` guards the convergence/stall checks so a
    coarse-grid Fock can never be accepted as converged.

> **Note (cosmetic):** the "Delta E" column on the first production-grid row
> after the switch spans the coarse→production grid change, so that single value
> is not a true SCF energy step (the energy column itself is correct). The
> "Switching XC grid" banner is printed immediately above it as a cue.

This change was reviewed with a multi-agent adversarial pass; the gating guards
above (maxit, level-shift, and the point-count cheapness check) close the edge
cases it surfaced.

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

Minnesota hazard check: with `OQP_XC_C2F=1`, an **M06-2X** run prints
`Coarse-to-fine XC grid schedule requested but disabled: grid-sensitive
(Minnesota) functional` and runs entirely on the production grid.

Raw logs: `/Volumes/External_Storage/claude/sessions/20260627_xc_c2f_grid/`.
