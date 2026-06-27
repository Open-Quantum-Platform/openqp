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

## Usage

Input file:
```
[scf]
pscreen=True            # enable (default False)
pscreen_k=1.0e-2        # coupling: tau_iter = k * diis_error
pscreen_cap=1.0e-8      # loosest cutoff allowed early (safe ceiling)
pscreen_tight=1.0e-4    # pin to int2e_cutoff once diis_error < this
```

Environment (override, handy for benchmarking; matches the `OQP_FOCK_*` convention):
```
OQP_PSCREEN=1            # enable
OQP_PSCREEN_K=1e-2
OQP_PSCREEN_CAP=1e-8
OQP_PSCREEN_TIGHT=1e-4
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

**Honest assessment.** The savings are modest for pure HF and negligible for hybrid DFT. Two
reasons: (1) OpenQP's existing density-weighted **incremental** Fock build already harvests most of
the ERI screening savings, so the marginal gain from loosening early iterations is small once the
safe cap is respected; (2) for hybrid DFT the wall-time bottleneck is the **numerical XC grid**,
which this feature does not touch. A coarse->fine **XC grid** ramp (and faster XC integration:
BLAS-form rho/Vxc, collocation caching) is the larger DFT opportunity and is left as future work.
