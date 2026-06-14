# BLAS/LAPACK Integer ABI Benchmarks

This directory records the BLAS/LAPACK comparison that motivated the macOS
LP64 exception in OpenQP's build system.

## Policy recommendation

- Use ILP64 BLAS/LAPACK by default.
- Allow LP64 only on macOS for native Accelerate, where it is the practical
  Apple Silicon path and the fastest tested option.
- Treat Homebrew OpenBLAS LP64 as unsupported for normal use: it is slower on
  these benchmarks and failed the local TDDFT cases.
- If OpenBLAS is built from source, build and use the ILP64 interface.
- Use ILP64 MKL on Intel systems where MKL is available.
- Keep Linux/HPC builds ILP64-only.

## Benchmark setup

The benchmark molecule is thymine, using PubChem CID 1135 3D coordinates. Each
BLAS variant was run three times with:

- BHHLYP/6-31G* DFT, TDDFT, and MRSF-TDDFT energy inputs.
- BHHLYP/cc-pVDZ DFT, TDDFT, and MRSF-TDDFT energy inputs.
- `OMP_NUM_THREADS=8`.
- Single-thread BLAS where the BLAS library exposes a direct thread control
  environment variable.

The exact input decks are in `inputs/`.

## Results summary

Mean OpenQP wall times are in `thymine_bhhlyp_blas_means.tsv`. Failed cases are
listed in `thymine_bhhlyp_blas_failures.tsv`.

On Apple Silicon, Accelerate LP64 is fastest for every successful thymine case.
Source OpenBLAS ILP64 is stable but slower. Homebrew OpenBLAS LP64 is slower and
fails the local TDDFT cases with LAPACK/Davidson breakdowns.

On the Intel Mac test system, MKL ILP64 is fastest overall. Source OpenBLAS
ILP64 is consistently faster than Accelerate and much faster than Homebrew
OpenBLAS LP64. Homebrew OpenBLAS LP64 is not competitive on the larger thymine
workloads.

On the Linux/HPC test system, MKL ILP64 was tested through the site module
environment and provides the Linux/HPC reference numbers in the summary table.
