# Build Options

This page lists the OpenQP CMake options exposed by the source tree. It is meant
for source builds, package maintainers, and developers who need to tune native
dependencies.

Most users should install OpenQP with:

```bash
pip install openqp
```

For source checkouts, the top-level Python package build uses scikit-build and
passes a few CMake overrides from `pyproject.toml`. A direct CMake build uses
the defaults in `CMakeLists.txt`.

## Passing Options

For direct CMake builds:

```bash
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=. \
  -DENABLE_OPENMP=ON \
  -DLINALG_LIB=auto
ninja -C build install
```

For `pip install .`, pass CMake options with `CMAKE_ARGS`:

```bash
CMAKE_ARGS="-DENABLE_MPI=ON -DLINALG_LIB=OpenBLAS" pip install .
```

The package build already sets several defaults for ordinary command-line use;
see [Python Package Build Defaults](#python-package-build-defaults).

## Core Options

| Option | Direct CMake default | Python package default | Values | Use |
| --- | --- | --- | --- | --- |
| `CMAKE_BUILD_TYPE` | `Release` | `Release` | `Release`, `Debug`, `RelWithDebInfo` | Select optimization and debug flags for single-configuration generators such as Ninja. |
| `CMAKE_INSTALL_PREFIX` | CMake default | `.` | Path | Install root for `lib`, `include`, and `share` runtime files. Package builds install into the wheel layout. |
| `CMAKE_C_COMPILER` | Auto | Auto | Compiler path | C compiler. Set with matching C++, Fortran compilers for reproducible builds. |
| `CMAKE_CXX_COMPILER` | Auto | Auto | Compiler path | C++ compiler, required for C++ wrappers and Libint verification. |
| `CMAKE_Fortran_COMPILER` | Auto | Auto | Compiler path | Fortran compiler. Use an MPI wrapper when building with MPI. |
| `BUILD_SHARED_LIBS` | `ON` | `ON` | `ON`, `OFF` | Build shared libraries when `ON`; `OFF` forces static-library lookup and static executable link flags for GNU. |
| `BUILD_TESTING` | `ON` | `ON` | `ON`, `OFF` | Enable CTest targets and optional smoke tests. |
| `ENABLE_PYTHON` | `ON` | `ON` | `ON`, `OFF` | Build the Python CFFI extension and install Python runtime files. |
| `ENABLE_OPENMP` | `OFF` | `ON` | `ON`, `OFF` | Enable OpenMP parallel kernels. |
| `ENABLE_MPI` | `OFF` | `OFF` | `ON`, `OFF` | Enable MPI support. Requires MPI and usually MPI compiler wrappers. |

## Integral, Optimizer, and Solvent Backends

| Option | Direct CMake default | Python package default | Values | Use |
| --- | --- | --- | --- | --- |
| `USE_LIBINT` | `ON` | `OFF` | `ON`, `OFF` | Build and use the external Libint two-electron integral path. `OFF` uses OpenQP's native Rys path. |
| `ENABLE_OPENTRAH` | `ON` | `OFF` | `ON`, `OFF` | Build and link the external OpenTrustRegion solver. `OFF` drops `otr_interface.F90` and uses native TRAH fallback behavior. |
| `ENABLE_DDX` | `OFF` | `OFF` | `ON`, `OFF` | Enable the optional ddX continuum-solvation backend. |
| `DDX_ROOT` | Unset | Unset | Path or environment variable | When `ENABLE_DDX=ON`, use a prebuilt ddX install. If unset, OpenQP builds ddX v0.8.0 with matching BLAS integer width. |

## BLAS and LAPACK

| Option | Default | Values | Use |
| --- | --- | --- | --- |
| `LINALG_LIB` | `auto` | `auto`, `netlib`, `FLAME`, `Intel10_64ilp`, `Intel10_64ilp_seq`, `Intel10_64_dyn`, `OpenBLAS` | Select the BLAS/LAPACK provider. |
| `CMAKE_FIND_ROOT_PATH_MODE_PACKAGE` | CMake default | `BOTH`, `ONLY`, `NEVER` | Useful for wheel/cross builds where CMake must search both host and package roots. |

`LINALG_LIB=auto` first tries a system BLAS/LAPACK. If nothing suitable is
found, OpenQP falls back to the bundled NetLib build. `LINALG_LIB=none` is
listed internally but rejected because OpenQP and OpenTrustRegion require
BLAS/LAPACK.

Every build is ILP64 -- one 8-byte BLAS/LAPACK integer model on every platform,
with no option to change it. `LINALG_LIB_INT64` was removed, and passing
`-DLINALG_LIB_INT64=OFF` fails the configure on purpose so that a stale cache
or an old build script cannot silently produce a mixed-width build.

Use ILP64 OpenBLAS or MKL on Linux, MKL ILP64 on Windows, and Accelerate's
`$NEWLAPACK$ILP64` interface on macOS 13.3 or newer.

## External Dependency Cache

OpenQP builds several bundled external dependencies when needed, including
Libint, NLopt, Libxc, tagarray, libecpint, NetLib LAPACK/BLAS, OpenTrustRegion,
fprettify, and ddX.

| Option | Default | Values | Use |
| --- | --- | --- | --- |
| `OQP_REUSE_EXTERNALS` | `ON` | `ON`, `OFF` | Reuse bundled external builds across fresh OpenQP build directories. |
| `OQP_EXTERNALS_ROOT` | Platform cache path | Path | Override the reusable external dependency cache root. |
| `OQP_EXTERNAL_ROOT` | Unset | Path | Backward-compatible alias for `OQP_EXTERNALS_ROOT`. |
| `OQP_LIBXC_C_FLAGS_RELEASE` | Empty | Compiler flags | Advanced: extra release C flags passed only to the external Libxc build. |

Default cache roots are:

| Platform | Default cache root |
| --- | --- |
| macOS | `~/Library/Caches/openqp/externals` |
| Linux with `XDG_CACHE_HOME` | `$XDG_CACHE_HOME/openqp/externals` |
| Linux with `HOME` | `~/.cache/openqp/externals` |
| Fallback | `<build>/external-cache` |

Disable reuse when debugging dependency rebuilds:

```bash
cmake -B build -G Ninja -DOQP_REUSE_EXTERNALS=OFF
```

Use a shared cache location on a build machine:

```bash
cmake -B build -G Ninja -DOQP_EXTERNALS_ROOT=/opt/openqp-cache/externals
```

## Developer and Diagnostic Options

| Option | Default | Values | Use |
| --- | --- | --- | --- |
| `ENABLE_ASAN` | `OFF` | `ON`, `OFF` | Enable AddressSanitizer flags for GNU builds. |
| `ENABLE_LSAN` | `OFF` | `ON`, `OFF` | Enable LeakSanitizer flags for GNU builds. |
| `ENABLE_UBSAN` | `OFF` | `ON`, `OFF` | Enable UndefinedBehaviorSanitizer flags for GNU builds. |
| `ENABLE_TSAN` | `OFF` | `ON`, `OFF` | Enable ThreadSanitizer flags for GNU builds. |
| `ENABLE_MSAN` | `OFF` | `ON`, `OFF` | Declared but not currently used. |
| `ENABLE_Formatter` | `ON` | `ON`, `OFF` | Add the Fortran formatting target using fprettify. |

Example debug build with runtime checks:

```bash
cmake -B build-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_ASAN=ON \
  -DENABLE_UBSAN=ON
ninja -C build-debug
```

## Python Package Build Defaults

The top-level `pip install .` build uses these CMake overrides from
`pyproject.toml`:

| CMake define | Package value | Reason |
| --- | --- | --- |
| `USE_LIBINT` | `OFF` | Prefer the native Rys path for package builds. |
| `ENABLE_OPENMP` | `ON` | Enable parallel production kernels by default. |
| `ENABLE_OPENTRAH` | `OFF` | Avoid the external OpenTrustRegion dependency in ordinary package builds. |
| `CMAKE_INSTALL_PREFIX` | `.` | Keep native runtime files package-local. |

Platform wheel builds add a few more environment-level settings:

| Platform | Wheel build settings |
| --- | --- |
| Linux | `LINALG_LIB=auto`, `CMAKE_FIND_ROOT_PATH_MODE_PACKAGE=BOTH` |
| macOS | Homebrew GCC compilers, `LINALG_LIB=auto` |

## Common Recipes

Linux development build:

```bash
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=. \
  -DUSE_LIBINT=OFF \
  -DENABLE_OPENMP=ON \
ninja -C build install
```

macOS Homebrew GCC build:

```bash
cmake -B build -G Ninja \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran-15 \
  -DCMAKE_INSTALL_PREFIX=. \
  -DENABLE_OPENMP=ON \
  -DLINALG_LIB=auto \
ninja -C build install
```

MPI build:

```bash
cmake -B build-mpi -G Ninja \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DENABLE_MPI=ON \
  -DENABLE_OPENMP=ON
ninja -C build-mpi install
```

ddX-enabled build:

```bash
cmake -B build-ddx -G Ninja \
  -DENABLE_DDX=ON \
  -DDDX_ROOT=/path/to/ddx/install
ninja -C build-ddx install
```

If `DDX_ROOT` is omitted, OpenQP builds the matching ddX dependency itself.
