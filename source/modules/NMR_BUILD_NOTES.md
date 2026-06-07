# NMR build notes

These notes record the build/install recipes verified for the gated NMR
CGO/GIAO development branch. They are build-environment notes only; they do not
claim native GIAO shielding support.

> **Caveat (historical):** the LP64 recipes below (`-DLINALG_LIB_INT64=OFF`,
> `BLA_SIZEOF_INTEGER=4`, `INTEGER_SIZE=4`) predate the upstream change that
> made ILP64 BLAS/LAPACK mandatory (current CMake rejects LP64 with a
> FATAL_ERROR). They are kept for the historical record of the ABI debugging
> only — build current sources with ILP64 (`LINALG_LIB_INT64=ON`, the default).

## Verified macOS Apple Silicon pip install

Use Homebrew GCC/GFortran and native Apple Accelerate BLAS/LAPACK. Clear stale
BLAS/LAPACK environment variables before invoking pip.

```bash
env -u BLAS_LIBRARIES -u LAPACK_LIBRARIES \
  CC=/opt/homebrew/bin/gcc-15 \
  CXX=/opt/homebrew/bin/g++-15 \
  FC=/opt/homebrew/bin/gfortran-15 \
  CMAKE_ARGS='-DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran-15 -DLINALG_LIB=auto -DLINALG_LIB_INT64=OFF -DBLA_VENDOR=Apple' \
  /opt/homebrew/bin/python3.11 -m pip install . -v
```

Verified evidence from the successful macOS install:

- BLAS: `/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework`
- LAPACK: `/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework`
- OpenTrustRegion and OpenQP were configured with `BLA_SIZEOF_INTEGER=4` / `INTEGER_SIZE=4`.
- The wheel installed as `PyOpenQP-0.1.0` for `/opt/homebrew/bin/python3.11`.

## Verified Lima/Linux pip install

Use the OpenQP Linux venv and OpenBLAS with the LP64 integer path.

```bash
env CC=gcc CXX=g++ FC=gfortran \
  CMAKE_ARGS="-DLINALG_LIB=OpenBLAS -DLINALG_LIB_INT64=OFF" \
  /home/cheolhochoi.guest/venvs/openqp311/bin/python -m pip install . -v
```

Verified evidence from the successful Lima install:

- BLAS: `/usr/lib/aarch64-linux-gnu/libopenblas.so`
- LAPACK: `/usr/lib/aarch64-linux-gnu/libopenblas.so`
- OpenTrustRegion and OpenQP were configured with `BLA_SIZEOF_INTEGER=4` / `INTEGER_SIZE=4`.
- The wheel installed as `PyOpenQP-0.1.0` for `/home/cheolhochoi.guest/venvs/openqp311/bin/python`.

## Pitfalls and boundaries

- On Apple Silicon, the Homebrew prefix is `/opt/homebrew`, not `/usr/local`.
- Stale `BLAS_LIBRARIES` or `LAPACK_LIBRARIES` values pointing at
  `/usr/local/opt/...` can poison macOS builds and produce misleading link
  failures such as unresolved BLAS symbols from `blas_wrap.F90.o`.
- The validated pip/install path here uses `LINALG_LIB_INT64=OFF`, i.e.
  `BLA_SIZEOF_INTEGER=4` / `INTEGER_SIZE=4`.
- Do not solve future ABI/build issues by hiding them behind global compiler,
  BLAS/LAPACK, or integer-size changes. Diagnose build-system and ABI problems
  separately from NMR/GIAO physics.
- For live tests in a repo shared between macOS and Linux, require the
  platform-specific library suffix (`liboqp.dylib` on Darwin, `liboqp.so` on
  Linux) before treating a repo-root build as usable.
