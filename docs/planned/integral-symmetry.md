# Petite-list integral reduction: it works, and it was silently off

## Correction

An earlier revision of this note recorded the reduction as **measured and worth
nothing** (benzene/cc-pVTZ, 40.1 s off vs 40.5 s on). That conclusion was wrong.
The benchmark compared C1 against C1: at cc-pVTZ the reduction never engaged.

## What it is actually worth

Freshly built binary, 8 OpenMP threads, D2h, RHF, wall clock of the SCF module
(the log's `Total ... Wall time`, not process wall clock -- process time is
dominated by Python startup and the molden write, which is how the effect was
missed the first time).

| system | nbf | `enabled=false` | petite, \|G\| = 8 | speedup |
|---|---|---|---|---|
| benzene / cc-pVTZ | 264 | 39.414 s | 10.937 s | **3.60x** |
| anthracene / 6-31G* | 230 | 7.694 s | 2.295 s | **3.37x** |
| naphthalene / 6-31G* | 166 | 2.921 s | 0.849 s | **3.40x** |

Energies agree to 1e-10 or better with identical iteration counts. The
non-abelian `full` tier reaches **8.06x** on benzene/cc-pVTZ (24 D6h operations).

The ceiling for a **planar** molecule is 4x, not 8x: sigma_h fixes every shell,
so it is in the stabilizer of every quartet and `q4 <= |G|/2`. The measured
3.4-3.6x is 85-90% of that. A non-planar D2h system is the case that could
approach 8x; nothing here measures one.

## Why the original benchmark saw nothing

`stage_integral_symmetry_maps` bailed with `skipped_basis_mismatch`. It built
the shell list as `(atom, l)` pairs with no purity flag, so `_normalize_shells`
defaulted every shell to Cartesian and `maps['n_ao']` came out 300 against
`nbf = 264`.

**The blast radius was much wider than cc-pVTZ.** Shell purity comes from the
Basis Set Exchange per-shell tag, and the 6-311G family is spherical too --
measured: naphthalene/6-311++G** gives nbf = 276, not the Cartesian 286. The
reduction worked **only** on bases whose l >= 2 shells BSE tags Cartesian, i.e.
6-31G*/6-31G** and below. Every cc-pVXZ, aug-cc-pVXZ, def2 and 6-311G run
silently fell back to C1.

## The fix

Shell purity is now exported from the library rather than guessed:
`oqp_get_basis_spherical` (`source/c_interop.F90`) returns the **effective**
per-shell flag, `HARMONIC_ACTIVE .and. harmonic(i) == 1 .and. am(i) >= 2`,
matching `c2s_ncomp` in `source/integrals/cart2sph.F90`.

Deriving it outside that one place had already gone wrong twice, in opposite
directions, and the AO total cannot arbitrate between them because Cartesian and
spherical sizes agree for l <= 1:

| shells passed to `build_reduction_maps` | n_ao vs nbf | max abs(T S T^T - S) | |
|---|---|---|---|
| all Cartesian (the old staging code) | 300 vs 264 | -- | bails, silently C1 |
| all pure (the naive fix) | 264 vs 264 | **1.2** | wrong maps, would corrupt the Fock |
| effective, `l >= 2` | 264 vs 264 | **3.1e-16** | exact |

Verified against the real overlap matrix from live benzene runs at cc-pVDZ and
cc-pVTZ. The existing `T S T^T = S` guard remains as the run-time check.

Measured per-shell flags: 6-31G* reports 0 spherical shells of 48 (so the
staging input is byte-identical to before on Cartesian bases); cc-pVTZ reports
24 of 96 -- 18 d and 6 f, with s and p correctly Cartesian.

## The fallback is no longer silent

`_dump_symmetry_log` was called only on the success path, so a run that asked
for the reduction and fell back to C1 printed nothing -- and the fallback is
energy-identical, so nothing else could reveal it. It is now emitted on every
exit, and a non-active status prints:

```
   integral reduction   : skipped_runtype_hess
   *** the petite-list reduction is NOT active: this run used the full (C1) integral list ***
```

`skipped_basis_mismatch` additionally reports the two AO counts.

## Still open

- **`tests/smoke_petite_*.py` are not collected.** `testpaths = ["tests"]` but
  pytest's default `python_files` never matches `smoke_*.py`, and no workflow or
  doc invokes them. `smoke_petite_benchmark.py` already asserts
  `status == 'active'` on two cc-pVTZ cases -- it would have caught this bug the
  day it was written. Waking it is the single highest-value follow-up, and an
  energy-only assertion cannot substitute: the fail-safe path gives the
  identical energy.
- **`benzene_full_dft` fails that benchmark**: dE = 3.14e-04 against a
  documented `full`-tier tolerance of 5e-7. Pre-existing and unrelated to this
  change (it is a 6-31G* case, where the staging input is unchanged), but it is
  a real defect in the non-abelian tier with DFT and now finally visible.
- **Do not enable by default yet.** Reorientation rewrites the geometry into the
  standard frame, and `coord` is a required, exactly-compared key in the example
  regression suite.
- **`runtype=grad` returns the gradient in the standard frame** with no
  back-rotation (`input_to_standard` is written and consumed nowhere). See
  `input-to-standard.md`.
- **Finite-difference child jobs reorient independently.** Reproduced here: a
  `runtype=hess` parent is correctly excluded (`skipped_runtype_hess`) while its
  children stage into *different* frames -- c2v, cs, and one
  `skipped_orientation_not_converged` -- within a single run. See
  `hessian-reorient.md`. This fix widens engagement to spherical bases, so it
  makes that defect bite on far more runs; it must land first.
- **`grd2` per-caller opt-in must land first** as well, for the same reason.
  See `grd2-optin.md`.
- `nint` and `thr_nshq` are in the `reduction(+:...)` clause at
  `source/integrals/int2.F90` but never initialised before the parallel region
  (only `nschwz` is), so their post-region values are indeterminate.
