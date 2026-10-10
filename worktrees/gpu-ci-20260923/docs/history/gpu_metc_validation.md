# GPU METC validation & benchmark runbook (C1c)

This document describes how to validate and benchmark the GPU MRSF/UMRSF METC
residency path **on a real Fortran/CUDA machine**. It also records the
correctness contract that the C1c-prep hardening enforces.

> **Status:** C1b (resident workspace refit) is complete and accepted. C1c —
> the actual `ENABLE_CUDA=ON` build, CPU/GPU parity validation, and performance
> benchmark — **has not been run yet**. No GPU performance is claimed anywhere
> in this branch. The scripts below are prepared so C1c can be executed
> immediately once GPU hardware is available.

## Correctness contract: all-or-nothing GPU sessions

Integrals are streamed to the contraction flush by flush and discarded after
each flush. Once the resident **device** `f3` becomes the live accumulator,
earlier flushes cannot be replayed on the host. Therefore a GPU failure that
happens *after* the session is active cannot be recovered by falling back to the
CPU — doing so would mix a partial device `f3` with a host `f3` and silently
return an incorrect Fock matrix.

The backend (`source/gpu_backend.F90`) models this with an explicit three-state
machine:

| State              | Meaning                                                        | On failure |
|--------------------|---------------------------------------------------------------|------------|
| `GPU_METC_INACTIVE`| No device accumulation yet; host `f3` is pristine.            | Safe **soft** fallback: the host owns the whole accumulation. Made sticky so a session can never start late and mix results. |
| `GPU_METC_ACTIVE`  | Session live; resident device `f3` holds (partial) sums.      | **Fatal abort** — never fall back (would mix partial device + host `f3`). |
| `GPU_METC_FAILED`  | Sticky, poisoned.                                             | Caller must abort; `finalize` never downloads a poisoned `f3`. |

Key invariants (locked by `tests/test_gpu_metc_no_mixed_fallback.py`):

- `gpu_backend_metc_try_contract` aborts (via `gpu_backend_metc_fatal`, which
  tears down the session and `error stop`s) on an **active-session** flush
  failure — it never returns quietly to let the caller mix in a host result.
- A **begin-time** failure is a legitimate soft CPU fallback by default, but is
  made sticky (`metc_soft_disabled`) so the GPU is not retried mid-accumulation.
- `gpu_backend_metc_finalize` downloads the device `f3` **only** when the state
  is `ACTIVE`.
- `gpu_backend_metc_reset` is called at the start of every fresh accumulation
  (`int2_mrsf_data_t_parallel_start`, `cur_pass == 1`) so a new build may
  attempt the GPU again.
- Both the MRSF (`.false.`) and UMRSF (`.true.`) update paths route through the
  same guarded `try_contract`, so neither can mix.

### Strict mode

`OQP_GPU_METC_STRICT=1` makes even a *begin-time* failure fatal (no CPU
fallback) — for users who explicitly require the GPU. An active-session failure
is always fatal regardless of this flag. The flag is read identically by the
Fortran backend and by `pyoqp/oqp/utils/gpu.py` (`gpu_metc_strict`).

For validation-only smoke tests, `OQP_GPU_METC_INDUCE_BEGIN_FAIL=1` forces
`oqp_gpu_metc_session_begin` to fail before any device accumulation. With
`OQP_GPU_METC_STRICT=1`, this must produce a nonzero exit and the fatal
all-or-nothing diagnostic; without strict mode it remains a soft begin-time CPU
fallback. Do not use this hook for benchmarks or production calculations.

## Environment variables

| Variable                 | Effect                                                          |
|--------------------------|----------------------------------------------------------------|
| `OQP_GPU_METC=1`         | Enable the GPU METC residency path.                            |
| `OQP_GPU_METC_STRICT=1`  | Make a begin-time GPU failure fatal instead of CPU fallback.   |
| `OQP_GPU_METC_TIMING=1`  | Emit a separated-region timing line to stderr at session end.  |
| `OQP_GPU_METC_INDUCE_BEGIN_FAIL=1` | Validation-only hook: force METC session begin failure for strict-mode abort smoke tests. |

## 1. Build with CUDA

```bash
cmake -S . -B build-cuda -DENABLE_CUDA=ON
cmake --build build-cuda -j
```

There is a single `ENABLE_CUDA` option in the build; do not add a second one.

## 2. CPU/GPU parity validation

```bash
python tools/gpu/validate_metc_cuda.py \
    --run-cmd openqp \
    --example examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.json \
    --tol 1e-7
```

The harness runs the case with `OQP_GPU_METC=0` and `OQP_GPU_METC=1` and
compares scalar results within a relative tolerance.

- Exit `0`: parity validated (correctness only — **no** speedup implied).
- Exit `1`: parity mismatch or a run/build error.
- Exit `2`: SKIPPED — no CUDA toolchain (this is **not** a pass).

The scalar extractor (`extract_scalars`) is intentionally simple; adapt
`SCALAR_KEYWORDS` to the exact OpenQP result fields you want to compare on the
target machine. Parity must be green **before** any benchmarking.

## 3. Performance benchmark

```bash
python tools/gpu/benchmark_metc.py \
    --run-cmd openqp \
    --example examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.json \
    --repeat 5 --warmup 1
```

The harness measures end-to-end wall time for the CPU and GPU runs and parses
the `OQP_GPU_METC_TIMING` line (upload `d3` / zero `f3` / download `f3`) so the
device-transfer regions are reported separately from kernel/host time. The
per-flush contraction is deliberately **not** wall-clock timed (it is the hot
path), so the kernel region is derived at the run level.

**Honesty contract:** the benchmark only prints a speedup when it has actually
measured a GPU run on a CUDA machine. With no GPU it reports CPU-only numbers
and explicitly claims no speedup. It never extrapolates or fabricates numbers.

## Definition of done for C1c

1. `ENABLE_CUDA=ON` build succeeds.
2. `validate_metc_cuda.py` exits `0` on a deterministic MRSF (and UMRSF) case.
3. `benchmark_metc.py` produces measured CPU vs GPU timings with the separated
   device-transfer regions, on real hardware.

Only after steps 1–3 may any GPU performance number be reported.
