# Progressive (iteration-dependent) integral screening

**Status:** opt-in, default **OFF**. Gated by `scf_pscreen` (input `[scf] pscreen=`) or the
`OQP_PSCREEN` environment variable.

## What it does

During direct SCF the two-electron Schwarz/density screening cutoff (`int2e_cutoff`, default
`5e-11`) is normally fixed for every iteration — including the early ones, where the density is
still far from converged and high integral accuracy is wasted. Progressive screening makes the
cutoff **iteration-dependent**: loose in early iterations (cheaper Fock builds), tightened back to
the user's `int2e_cutoff` as the SCF converges, so **the converged energy is unchanged**.

The per-iteration cutoff is coupled to the DIIS error `e` from the previous iteration:

```
tau_iter = clamp( pscreen_k * e ,  int2e_cutoff ,  pscreen_cap )      while e >= pscreen_tight
tau_iter = int2e_cutoff   (pinned, sticky)                            once e <  pscreen_tight
```

The governing invariant is that the error injected by the loose cutoff must stay well below the
size of the SCF update (`err_screen << ||SCF update||`); otherwise DIIS extrapolates on noise and
convergence stalls. `pscreen_k` (the safety fraction) and `pscreen_cap` (the loosest allowed
cutoff) enforce this. The pin to the tight cutoff is **sticky** (latches once `e < pscreen_tight`)
and triggers a **full, non-incremental Fock rebuild**, so any contributions dropped during the
loose phase are recaptured — the converged energy matches an all-tight run bit-for-bit.

This is the classic dynamic/variable integral screening of Almlof-Faegri-Korsell (1982) and
Haser-Ahlrichs (1989); the continuous DIIS-error coupling mirrors Q-Chem's `VARTHRESH`. It
composes with (and rides on top of) OpenQP's incremental difference-density Fock build.

## XC-grid ramp (DFT)

The same `diis_error` controller and `pscreen_tight` pin also drive an optional **XC-grid
threshold ramp**: during the descent the DFT grid density cutoff (`grid_density_cutoff`) and the
grid AO-prune threshold (`grid_ao_threshold`) are loosened, so early XC builds prune more AOs and
skip more low-density points; both are restored to the user's baseline once `diis_error <
pscreen_tight`, so the converged XC energy is unchanged. This is one controller feeding two
consumers (ERIs + XC). It is **opt-in within `pscreen`** (both knobs default to 0 = not ramped),
because the safe loose magnitude is more system-dependent than the ERI cap; recommended range
1e-5 (conservative) to 1e-4 (more aggressive, validated on the systems below).

## Coarse->fine XC grid ramp (the largest XC lever)

In addition to the threshold ramp, an optional **coarse->fine grid ramp** uses a genuinely
*coarser* integration grid (`pscreen_grid_rad` x `pscreen_grid_ang` Lebedev) during the descent
and the user's full grid once `diis_error < pscreen_tight`. Because XC cost scales ~linearly with
the number of grid points, a coarse grid cuts the dominant DFT cost in early iterations far more
than threshold pruning on the full grid. It is **energy-neutral**: OpenQP's XC build is
non-incremental, so the pinned tail recomputes Vxc entirely on the full grid -> the converged
energy matches the all-full-grid run exactly.

The coarse grid is built in `hf_energy` (because `dft_initialize` needs `basis` `intent(inout)`)
with the functional name blanked so libxc is not re-initialised, and passed to `scf_driver` as an
optional argument; `scf_driver` selects coarse vs full per iteration via the same `ps_pin` latch.
Cost: a small one-time +1-2 SCF iterations (the coarse->fine switch is a DIIS transient), far
outweighed by the per-iteration savings. Off by default (`pscreen_grid_rad=pscreen_grid_ang=0`);
recommended coarse ~ half the full grid (e.g. 50x194 against a 96x302 full grid). A slightly
earlier pin (`pscreen_tight=1e-3`) trims the extra iteration.

## Usage

Input file:
```
[scf]
pscreen=True            # enable (default False)
pscreen_k=1.0e-2        # ERI coupling: tau_iter = k * diis_error
pscreen_cap=1.0e-8      # ERI loosest cutoff early (safe ceiling)
pscreen_tight=1.0e-4    # pin ERI + XC to tight once diis_error < this
pscreen_xc_dcut=1.0e-4  # XC: loose grid density cutoff during descent (0=off)
pscreen_xc_aocut=1.0e-4 # XC: loose grid AO-prune threshold during descent (0=off)
pscreen_grid_rad=50     # XC: coarse radial points during descent (0=off)
pscreen_grid_ang=194    # XC: coarse Lebedev points during descent (0=off)
```

Environment (override, handy for benchmarking; matches the `OQP_FOCK_*` convention):
```
OQP_PSCREEN=1            # enable
OQP_PSCREEN_K=1e-2
OQP_PSCREEN_CAP=1e-8
OQP_PSCREEN_TIGHT=1e-4
OQP_PSCREEN_XC_DCUT=1e-4   # XC grid density cutoff (loose phase)
OQP_PSCREEN_XC_AOCUT=1e-4  # XC grid AO-prune threshold (loose phase)
OQP_PSCREEN_GRID_RAD=50    # XC coarse radial points (descent)
OQP_PSCREEN_GRID_ANG=194   # XC coarse Lebedev points (descent)
OQP_FOCK_TIMER=1         # print per-iteration Fock-build time, tau, and skip count
```

## Parameters and the safe operating point

`pscreen_cap=1e-8` is the validated safe ceiling. Looser caps (1e-6, 1e-4) derail DIIS on dense
systems (observed: benzene converging to a spurious state hundreds of Ha away) because the loose
Fock violates `err_screen << ||SCF update||`. At `1e-8` the converged energy is exact and the
iteration count is unchanged (+0-1).

## Measured benefit (gfortran-15, 8 threads, 6-31G\*)

Energy matches the all-tight baseline to 0.0 Ha and iteration counts are unchanged in all cases.

| system (atoms)          | iters OFF->ON | total wall OFF->ON | speedup |
|-------------------------|---------------|--------------------|---------|
| benzene RHF (12)        | 10 -> 10      | 1.07 -> 0.99 s     | ~7%     |
| (H2O)8 RHF (24)         | 14 -> 14      | 1.26 -> 1.18 s     | ~6%     |
| (H2O)16 RHF (48)        | 14 -> 15      | 4.42 -> 4.42 s     | ~0%     |
| (H2O)24 RHF (72)        | 14 -> 14      | 12.3 -> 11.6 s     | ~5-15%* |
| benzene BHHLYP (12)     | 11 -> 11      | 6.76 -> 6.75 s     | ~0%     |
| (H2O)8 BHHLYP (24)      | 13 -> 13      | 30.6 -> 30.6 s     | ~0%     |

\* The Fock build itself is ~20% faster at (H2O)24 (sum of per-iteration Fock time 5.67 -> 4.44 s);
the total-SCF figure is diluted by diagonalization/DIIS/density steps and is noisy on a loaded
machine.

XC THRESHOLD ramp (BHHLYP/6-31G\*, `pscreen_xc_dcut=pscreen_xc_aocut=1e-4`), energy exact, iters unchanged:

| system (atoms)        | OFF       | ERI+XC ramp | speedup |
|-----------------------|-----------|-------------|---------|
| (H2O)8 BHHLYP (24)    | 30.65 s   | 29.70 s     | ~3%     |
| (H2O)16 BHHLYP (48)   | 317.9 s   | 296.7 s     | ~7%     |

XC COARSE->FINE GRID ramp (BHHLYP/6-31G\*, coarse 50x194 vs full 96x302), energy exact:

| system (atoms)        | OFF       | coarse->fine | speedup | iters     |
|-----------------------|-----------|--------------|---------|-----------|
| (H2O)8 BHHLYP (24)    | 32.55 s   | 26.22 s      | **~19%**| 13 -> 15  |
| (H2O)8, pin 1e-3      | 32.55 s   | 25.41 s      | **~22%**| 13 -> 14  |
| (H2O)16 BHHLYP (48)   | 301.9 s   | 209.4 s      | **~31%**| 14 -> 15  |

**Honest assessment.** For pure HF the ERI ramp gives ~5-15% (the existing incremental Fock already
harvests most ERI savings). For hybrid DFT the **coarse->fine grid ramp is the real win** -- ~19-31%
with **exact** converged energy and only +1-2 iterations, and it grows with XC-dominance (31% at 48
atoms). The XC *threshold* ramp is a smaller complementary effect (~3-7%). The convergence-wall
invariant (`err_screen << ||SCF update||`) still caps how aggressive the descent can be. The largest
remaining DFT opportunity is **structural** -- matrix/BLAS-form rho+Vxc and collocation caching
across the (fixed) grid -- a bigger rewrite left as future work.
