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

The benchmark emits deterministic MRSF and UMRSF cases.  The CUDA result is compared with a CPU reference implementation of the same METC update equations.  Reported GPU timings include the current prototype's full host/device tensor copies and per-call allocations; these are intentionally conservative and should not be used as final optimized speedup claims.

## Acceptance criteria for manuscript-quality correctness

- `ierr == 0` for all cases.
- `max_abs_diff < 1e-9` for all isolated contractions.
- Final OpenQP MRSF/UMRSF excitation energies should agree with CPU reference to within `1e-6` eV before performance claims are made.
