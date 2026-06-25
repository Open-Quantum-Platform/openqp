# Installation

## Recommended Install

Use the Python package when possible:

```bash
pip install openqp
```

For a local source checkout:

```bash
git clone https://github.com/Open-Quantum-Platform/openqp.git
cd openqp
pip install .
```

The top-level package build installs the Python package, native library, header
files, and data files together. Normal command-line use does not require
`OPENQP_ROOT` after installation.

## Requirements

- Python 3.9 or newer
- GCC, G++, and Gfortran
- CMake 3.25 or newer
- BLAS/LAPACK
- `cffi`, NumPy, SciPy, and geomeTRIC
- Ninja, recommended for source builds
- OpenMPI or another MPI implementation, only when building with MPI

## Source Build

The default source install is:

```bash
pip install .
```

For development builds where you want to inspect the native build directory:

```bash
cmake -B build -G Ninja \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_INSTALL_PREFIX=. \
  -DENABLE_OPENMP=ON
ninja -C build install
cd pyoqp
pip install .
```

On macOS, prefer Homebrew GCC and the native Accelerate BLAS/LAPACK stack:

```bash
cmake -B build -G Ninja \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran-15 \
  -DCMAKE_INSTALL_PREFIX=. \
  -DENABLE_OPENMP=ON \
  -DLINALG_LIB=auto \
  -DLINALG_LIB_INT64=OFF
ninja -C build install
cd pyoqp
pip install .
```

Adjust the compiler suffix to match the Homebrew GCC version installed on the
machine.

## Common CMake Options

| Option | Default | Meaning |
| --- | --- | --- |
| `-DENABLE_MPI=ON` | `OFF` | Enable MPI support. Use an MPI Fortran compiler wrapper such as `mpif90`. |
| `-DENABLE_OPENMP=ON` | `OFF` in CMake, `ON` for Python package builds | Enable OpenMP parallel sections. |
| `-DUSE_LIBINT=ON` | `ON` in CMake, `OFF` for Python package builds | Use Libint for ERIs instead of the native Rys path. |
| `-DLINALG_LIB=<vendor>` | `auto` | Select BLAS/LAPACK provider. |
| `-DLINALG_LIB_INT64=ON` | `ON` | Use ILP64 BLAS/LAPACK. |
| `-DENABLE_OPENTRAH=OFF` | `ON` in CMake, `OFF` for Python package builds | Skip the external OpenTrustRegion library and use native TRAH. |
| `-DOQP_REUSE_EXTERNALS=OFF` | `ON` | Disable reusable bundled-external build caches. |

ILP64 BLAS/LAPACK is the normal build mode. LP64
(`-DLINALG_LIB_INT64=OFF`) is supported only on macOS, mainly for a consistent
native Accelerate build.

## Runtime Files

Installed packages resolve runtime files package-locally first. Source-tree
development layouts are also detected when the native library has been installed
into the checkout. Keep `OPENQP_ROOT` only as a compatibility fallback for
custom layouts where Python and the OpenQP runtime tree are separated.

## OpenMP Threads

OpenQP accepts the OpenMP thread count from the command line:

```bash
openqp h2o.inp --omp 16
```

or from the input file:

```ini
[input]
omp_threads=16
```

Precedence is `--omp`, then `[input] omp_threads`, then `OMP_NUM_THREADS`, then
the built-in default.

## Test

```bash
openqp --run_tests all
```

For a smaller first check:

```bash
openqp examples/HF/H2O_RHF-HF_ENERGY.inp
```
