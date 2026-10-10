# METC GPU benchmark harness

This directory contains the standalone correctness/performance harness for the experimental OpenQP MRSF/UMRSF CUDA METC contraction path.

## Build/run on a CUDA node

```bash
nvcc -O3 -allow-unsupported-compiler -std=c++14 \
  benchmarks/metc_gpu/metc_cuda_benchmark.cu \
  -o metc_cuda_benchmark
./metc_cuda_benchmark > metc_cuda_benchmark.csv
python3 benchmarks/metc_gpu/analyze_metc_benchmark.py \
  metc_cuda_benchmark.csv --outdir metc_analysis
```

The benchmark emits deterministic MRSF and UMRSF cases.  The CUDA result is compared with a CPU reference implementation of the same METC update equations.

### Two timing columns — read the right one

- `gpu_total_ms` / `speedup_total`: the integrated `oqp_gpu_metc_contract` path, i.e. per-call `cudaMalloc` + full H2D transfer of `d3`/`f3` + kernel + D2H + `cudaFree`. This is **dominated by data movement** and must **not** be quoted as kernel/GPU performance.
- `kernel_ms` / `speedup_kernel`: the **compute ceiling** — only the kernel launch, with device buffers pre-allocated and pre-populated, 5-launch warmup, averaged over 50 `cudaEvent`-timed launches. Quote this when judging the kernel itself.

The gap between the two columns is the headroom the planned device-residency redesign would recover (upload `d3` once, keep `f3` resident and accumulated on device across buffer flushes, stream only `(ids, ints)` batches, per-thread CUDA streams, download `f3` once).

Note that the METC contraction has very low arithmetic intensity (~2 flops per f3 cell touched, against an 8-byte `d3` read plus an atomic RMW), so even `kernel_ms` is bandwidth/atomic-throughput bound rather than FLOP bound; the "compute ceiling" here is effectively an L2-atomic-throughput ceiling.

## Acceptance criteria for manuscript-quality correctness

- `ierr == 0` for all cases.
- `max_abs_diff < 1e-9` for all isolated contractions.
- Final OpenQP MRSF/UMRSF excitation energies should agree with CPU reference to within `1e-6` eV before performance claims are made.
