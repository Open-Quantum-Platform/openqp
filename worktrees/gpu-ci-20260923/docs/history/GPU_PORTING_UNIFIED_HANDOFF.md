# `feat/gpu-porting-unified` — Handoff

**Branch:** `feat/gpu-porting-unified`
**Worktree:** `/Users/cheolhochoi/Documents/claude/openqp-private-opeqp-GPU`
**Canonical source repo:** `/Users/cheolhochoi/Documents/claude/gradient-fix/openqp`
**Base:** `origin/main` @ `1890fef` (Fix post-v1.1 test failures, #159)
**Remote home:** `private` only (`https://github.com/karmachoi/openqp-private`). **Never push to public `origin`.**
**Status at handoff:** clean worktree from `origin/main`; this handoff doc is the first commit. No GPU code merged yet.

---

## Goal

One long-term GPU control plane for OpenQP. **First executable target = one MRSF-TDDFT / TDHF response path with deterministic CPU/GPU parity.** HF/DFT (J/K, quadrature, ERI) is *designed-for* but **not built** in this milestone.

## Branch inventory (vs `origin/main` @ 1890fef)

| Tag | Branch | HEAD | Role |
|---|---|---|---|
| XCR | `feat/gpu-xc-response-rebase` | 7f9eda0 | Clean additive XC-response scaffold + `ENABLE_CUDA` build option |
| WSM | `feat/gpu-workspace-manager` | 6413215 | Unified `GpuWorkspaceManager` (residency/cache). **Subsumes PMB + XCC.** |
| MRT | `feat/gpu-metc-regression-test` | fc92929 | Real METC CUDA kernel + CPU/GPU parity test. **Strict superset of `private/feat/gpu-metc`.** |
| DVT | `perf/tdhf-davidson-timers` | 6064b5b | `OQP_TIMER` timing/observability layer |
| PMB | `perf/gpu-metc-persistent-buffers` | fb8309e | Superseded by WSM — **do not merge directly** |
| XCC | `perf/tdhf-xc-response-cache` | 7a704d1 | Superseded by WSM — **do not merge directly** |
| METC | `private/feat/gpu-metc` | 5e4b23e | Ancestor of MRT; **do not merge** (catch-up merge `5e4b23e` drags DFTB/JSON) |

## Source of truth per subsystem

| Subsystem | Owner | Reconciliation |
|---|---|---|
| GPU config (`GpuConfig`) | MRT `gpu.py` ⊕ XCR's `xc_response` target | one config; `target` multi-valued (see decision 4) |
| Workspace / residency / cache | **WSM** `gpu_workspace.py` | dock onto unified `GpuConfig` |
| CUDA build block | XCR | fold MRT's block in; one centralized block (decision 7) |
| METC kernel | MRT | rewrite host wrapper onto workspace pointers |
| XC-response scaffold | XCR | wire `.cu` to device pointers from manager |
| Timing / observability | DVT | generalise Davidson manifest → namespaced registry (`metc.*`, `xc_response.*`) |

## Approved consolidation sequence

```
origin/main
→ feat/gpu-xc-response-rebase            (XCR — scaffold + ENABLE_CUDA)
→ feat/gpu-workspace-manager             (WSM — Stage-1 schema foundation; carries PMB+XCC)
→ METC via curated cherry-picks from feat/gpu-metc-regression-test (NOT private/feat/gpu-metc, NOT merge 5e4b23e)
→ perf/tdhf-davidson-timers              (DVT — only if still clean/orthogonal)
→ kernel residency refit                 (new work: d3→RESIDENT_INPUT, f3→RESIDENT_ACCUM)
```

Do **not** merge `perf/gpu-metc-persistent-buffers` or `perf/tdhf-xc-response-cache` independently (subsumed by WSM).
Do **not** merge `private/feat/gpu-metc` directly (use MRT; avoid its DFTB/JSON catch-up merge).

### Curated METC commits from MRT (drop catch-up merge `5e4b23e`)
- `5524e67`(+`61bfe9e`,`4a0a7e0`) config/build scaffold → resolve into the ONE CUDA block + unified `GpuConfig`
- `f70c5ed`(+`89a1cd6`,`d945e6f`) kernel + MRSF wiring → **rewrite** onto residency manager
- `46b9fd0`(+`6fe6869`,`43594ce`,`fc92929`) benchmark + bench timing
- `0ff7785` CPU-vs-GPU parity test (MRT-only) — the parity gate

## Decision resolutions (locked)

1. **Remote home:** `private` only. `git push -u private feat/gpu-porting-unified`. Never public `origin`.
2. **Route:** WSM supersedes PMB + XCC for direct merge.
3. **Sequence:** as above (XCR → WSM → METC-via-MRT → DVT → refit).
4. **`GpuConfig.target` multi-valued** with backward-compatible singleton coercion:
   `"metc" → ("metc",)`, `"xc_response" → ("xc_response",)`, `("metc","xc_response") → both`.
   Workspace manager still namespaces buffers per concrete target string.
5. **`f3` accumulator:** `PER_THREAD` is the first implementation policy (mirrors per-OpenMP-thread `this%f3`).
   Keep `SINGLE_ATOMIC` as a supported schema mode but **not** the default. First refit preserves the CPU threading model.
6. **History style:** curated cherry-picks for METC/MRT. Never merge whole private branches with unrelated catch-up commits.
7. **CMake:** one centralized CUDA integration path. No competing CUDA blocks. Top-level CMake baseline stays controlled until the XCR CUDA block is deliberately integrated.

## Risks / gotchas

- METC kernel currently does `cudaMalloc/H2D(full d3,f3)/kernel/D2H/cudaFree` **every Davidson flush** — the perf thesis is making `d3` RESIDENT_INPUT (upload once/pass) and `f3` RESIDENT_ACCUM (zero once, download once).
- `atomicAdd` ⇒ FP order non-deterministic ⇒ parity must stay **tolerance-based** (`rtol=1e-9, atol=1e-11`), never bit-exact.
- `this%f3` is per-OpenMP-thread inside the integral-buffer flush callback ⇒ device residency must be `PER_THREAD` (`nthreads×T`) for the first refit.
- Conflicts to hand-resolve: two `ENABLE_CUDA` blocks, two `GpuConfig`, two `_check_gpu` (XCR vs MRT). DVT touches different Fortran files (`util.F90`, `tdhf_mrsf_energy.F90`) so low conflict risk vs METC's `tdhf_mrsf_lib.F90`.

## Operational notes

- Branch was created from `origin/main` and initially auto-tracked `origin/main`. The `git push -u private` retargets upstream to `private` (closing the "bare push hits public origin" footgun).
- `private` remote is HTTPS (`https://github.com/karmachoi/openqp-private`); non-interactive auth may require a credential helper / PAT. SSH form `git@github.com:karmachoi/openqp-private.git` avoids this if preferred.
- Full audit + rationale: `GPU_PORTING_CONSOLIDATION_PLAN.md` and `GPU_HFDFT_ARCHITECTURE.md` in `gradient-fix/`.
