# DFT XC-grid "compute-once / reuse" optimizations

**Default-OFF**, env-gated optimizations of the SCF DFT cost. OpenQP already
reuses the J/K (HF) part of the Fock matrix incrementally across SCF iterations
(`scf.F90` `dold`/`fold`); these target the XC build, which previously rebuilt
everything from the density every iteration.

| Env var | Default | Effect |
|---|---|---|
| `OQP_XC_PHI_CACHE=1` | off | **Opt 1** — cache the collocation matrix Φ across SCF iterations (exact) |
| `OQP_XC_INCDFT=1` | off | **Opt 2** — incremental XC reuse (experimental; *measured non-improvement*, see below) |
| `OQP_XC_TIMING=1` | off | per-iteration timers: `[SCFTIME]` (J/K vs XC wall split) and `[XCTIME]` (geom-Φ vs density-driven) |

Measurements: Apple-silicon macOS, gfortran-15, Release, 4 OpenMP threads,
6-31G(d), default Becke grid. The box was heavily loaded, so **per-iteration
ratios and SCF iteration counts are the robust metrics**; total-wall numbers are
indicative.

---

## Opt 1 — Collocation-Φ cache (exact)

`Φ[μ,i] = φ_μ(r_i)` (basis values + grid derivatives) depends only on
geometry + basis + grid + integration threshold, **not** on the density, yet the
native grid loop (`run_xc` → `compAOs`/`pruneAOs`) recomputed it on every Fock
build. Opt 1 stores the post-pruning Φ block (plus weights and significant-AO
metadata) per grid slice on the first build of a geometry and replays it on later
SCF iterations.

* New module `source/dftlib/dft_gridint_phi_cache.F90`; build/replay branch in
  `run_xc`; `resetPrunedPointers` gains a `gather=.false.` mode (replay re-does
  only the density-dependent wavefunction compression, skips the geometry-only AO
  gather). Integrates with the symmetry-reduced weight path (the symw-scaled
  weights are geometry-only and are cached as-is).
* Opt-in only from the density-driven SCF Fock build (`dmatd_density_blk`, the
  path `calc_fock` takes). It is deliberately **not** set in `dmatd_blk`, which
  one-shot callers reach (e.g. the finite-difference Hessian via `dftexcor`, at a
  fresh displaced geometry each time) where the per-slice cache would be built but
  never replayed — pure memory/copy overhead. Gradients/response use their own
  consumers and never cache. Validity keyed by a geometry hash + derivative order
  + DFT threshold ⇒ rebuilds transparently on any change.
* Replayed Φ is **bit-for-bit identical**, so the converged energy is unchanged.

**Validation** (converged energy, `OQP_XC_PHI_CACHE` 1 vs 0):

| system / functional | ΔE (Ha) | iters | build → replay |
|---|---|---|---|
| benzene / B3LYPV5 | 0.0e+00 | 9 → 9 | 1 build, 8 replays |
| benzene / M06-2X (meta-GGA) | 0.0e+00 | 9 → 9 | 1 build, 8 replays |
| benzene / PBE | 0.0e+00 | 9 → 9 | 1 build, 8 replays |
| C20H42 / PBE | 0.0e+00 | 12 → 12 | 1 build, 11 replays |

**Speedup.** The cache removes ~85–90 % of the collocation (`compAOs`) cost; the
total benefit depends on how large a share of the XC build that collocation is —
set by grid density and by whether the build is MO- or density-driven (the SCF
path is density-driven, where the per-point ρ/libxc/assembly work is relatively
larger, so collocation is a smaller slice):

| case | XC build / iter (off → replay) | XC frac of Fock | total wall (off → on) |
|---|---|---|---|
| benzene PBE | 0.046 → 0.034 s (−26 %) | ~39 % | 1.73 → 1.62 s |
| C20H42 PBE | 1.07 → 0.93 s (−13 %) | ~34 % | 37.7 → 36.5 s |

The geom-Φ phase per iteration (aggregate thread-seconds, `[XCTIME]`) drops on
replay, confirming `compAOs` is eliminated, e.g. C20H42 PBE: 0.81 → 0.15 s.

**Cost / when it helps.**
* Memory: the cache is `Σ_slices numAOs_pruned × numPts × numAOVecs × 8 B`
  (C20H42 PBE on the default grid ≈ 4.8 GB), logged as `cacheMB` under
  `OQP_XC_TIMING`. Off by default precisely because of this memory/recompute
  trade.
* The win grows with coarser-relative-to-basis grids, meta-GGA, MO-driven builds,
  and whenever the 2e J/K cost is reduced (it composes with 2e-integral screening
  — as J/K shrinks, XC and this optimization become a larger share of SCF time).

---

## Opt 2 — Incremental XC (IncDFT): experimental, **measured non-improvement**

Implemented as whole-matrix XC reuse: store `V_xc[D_ref]`/`E_xc[D_ref]` from the
last full build and reuse it on iterations inside a late-SCF DIIS-error window,
with a periodic forced full rebuild and a return to full builds near convergence
(tied to the existing `scf.F90` incremental refresh, `diis_error < 1e-4`). The SCF
convergence test additionally refuses to accept a converged iteration whose XC was
reused, so the converged Fock/energy always come from a full XC build (exact final).
Gated `OQP_XC_INCDFT=1` (`source/dftlib/dft_incdft.F90`; plumbed `scf.F90` →
`calc_fock` → `calc_jk_xc`).

**Result: it does not help — it hurts convergence and is a net loss.** The
converged energy stays within the 1e-7 gate, but the iteration count blows up:

| system | baseline iters | IncDFT iters | reuses | ΔE (Ha) |
|---|---|---|---|---|
| benzene / PBE | 9 | 18 | 7 | 0.0e+00 |
| C20H42 / PBE | 12 | 31 | 16 | −7e-10 |

**Each reused XC build costs ≈ one extra SCF iteration.** Because XC is nonlinear
and the reused matrix is *global*, freezing it produces a Fock inconsistent with
the still-evolving density, which DIIS must work off; and one extra SCF iteration
(which includes a J/K build — always ≥ the XC build here) costs more than the XC
build it skipped. So whole-matrix XC reuse can only lose.

**The correct approach is per-batch screening** (Q-Chem-style incremental DFT):
keep `V_xc[D]` *accurate* by recomputing only the grid batches where `ΔP` is
significant and reusing the rest, via reference subtraction
`V = V_ref + Σ_active (V_s[D] − V_s[D_ref])`. Because the assembled Fock stays
correct to the screening threshold, convergence is **not** penalized, and the
saving is the screened-out batch fraction (grows late in SCF as `ΔP` sparsifies).
This needs the engine to evaluate a reference density per active batch (cheap
given the Opt-1 Φ cache) and per-slice `ΔP` screening — a larger change left as
future work. The current `OQP_XC_INCDFT` path is retained, off by default and
clearly experimental, as the gating/reference-store scaffolding for that work and
to make the negative result reproducible (`[SCFTIME] ... xc_reused=T/F`).

## Reproduce

```sh
# build (note LP64 BLAS flag, required on macOS/Accelerate):
cmake -G Ninja -B build -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15 \
  -DCMAKE_Fortran_COMPILER=gfortran-15 -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON -DENABLE_OPENMP=ON -DENABLE_PYTHON=OFF -DUSE_LIBINT=OFF \
  -DLINALG_LIB=auto -DLINALG_LIB_INT64=OFF
ninja -C build oqp

# run with timing, toggling the optimizations:
OQP_XC_PHI_CACHE=1 OQP_XC_TIMING=1  pyoqp input.inp     # Opt 1
OQP_XC_INCDFT=1    OQP_XC_TIMING=1  pyoqp input.inp     # Opt 2 (experimental)
```
