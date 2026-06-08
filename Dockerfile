# Use a base Linux image with Python support.
# ubuntu:24.04 (noble) ships gcc-14, matching the CI gcc-build matrix, so the
# image and CI compile with the same toolchain (avoids skew like the old
# focal/gcc-9 rejecting flags that gcc-14 accepts).
FROM ubuntu:24.04

# Set timezone to Asia/Seoul to avoid timezone prompt
ENV DEBIAN_FRONTEND=noninteractive
# Ubuntu 24.04 marks the system Python as externally managed (PEP 668); allow
# the image's pip installs into the system environment.
ENV PIP_BREAK_SYSTEM_PACKAGES=1
RUN ln -fs /usr/share/zoneinfo/Asia/Seoul /etc/localtime && \
    echo "Asia/Seoul" > /etc/timezone

# Build toolchain and runtime deps.  NOTE: no libblas-dev/liblapack-dev here --
# OpenQP requires an ILP64 (8-byte integer) BLAS/LAPACK, which the distro
# reference packages are not; with LINALG_LIB=auto CMake silently falls back to
# compiling the bundled (unoptimized) NetLib reference library.  We build
# OpenBLAS as ILP64 below instead, so the image links an optimized BLAS/LAPACK.
# (`make`/`perl` are needed by the OpenBLAS build.)
RUN apt-get update && apt-get install -y --fix-missing \
    gcc-14 g++-14 gfortran-14 cmake ninja-build make perl \
    openmpi-bin libopenmpi-dev \
    python3-pip wget git \
    && apt-get clean

# Install Python dependencies
RUN pip3 install cffi

# CMake comes from apt (24.04 ships 3.28 >= the required 3.25). The previous
# pinned download of a linux-x86_64 CMake tarball is dropped: it could not run on
# the arm64 leg of the multi-platform image (wrong-arch ELF under QEMU).

# Build OpenBLAS as ILP64 (INTERFACE64=1) with UNSUFFIXED symbols, so it exports
# `dgemm_` etc. with the 8-byte integer width OpenQP calls.  The distro's
# libopenblas64 instead exports `dgemm_64_`-suffixed symbols, which OpenQP does
# not call -- hence building from source here.  DYNAMIC_ARCH=1 keeps the image
# portable (OpenBLAS selects the kernel for the host CPU at run time, not the
# build machine's ISA), which also matters for the multi-platform image.
# Build the library targets only (libs netlib shared) -- `netlib` bundles LAPACK
# into libopenblas; we skip the `tests` target that `make all` runs, since a few
# edge-case extension self-tests can fail on some arches and would abort the build.
WORKDIR /tmp
ARG OPENBLAS_VERSION=0.3.28
RUN wget -q https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OPENBLAS_VERSION}/OpenBLAS-${OPENBLAS_VERSION}.tar.gz && \
    tar -xzf OpenBLAS-${OPENBLAS_VERSION}.tar.gz && \
    cd OpenBLAS-${OPENBLAS_VERSION} && \
    make -j"$(nproc)" INTERFACE64=1 USE_OPENMP=1 DYNAMIC_ARCH=1 \
         CC=gcc-14 FC=gfortran-14 libs netlib shared && \
    make PREFIX=/opt/openblas INTERFACE64=1 install && \
    ln -sf libopenblas.so /opt/openblas/lib/libopenblas64.so && \
    ln -sf libopenblas.a  /opt/openblas/lib/libopenblas64.a && \
    cd /tmp && rm -rf "OpenBLAS-${OPENBLAS_VERSION}" "OpenBLAS-${OPENBLAS_VERSION}.tar.gz"
# CMake's FindBLAS searches for an `openblas64`-named library when
# BLA_SIZEOF_INTEGER=8 (ILP64).  Our OpenBLAS uses UNSUFFIXED symbols (matching
# OpenQP's unsuffixed `dgemm` calls), so expose it under the openblas64 name via
# symlinks above -- without them, find_package(BLAS) misses it and OpenQP
# silently falls back to the bundled NetLib reference library.
ENV LD_LIBRARY_PATH="/opt/openblas/lib:${LD_LIBRARY_PATH}"

# Copy and compile the checked-out OpenQP source.  GitHub Actions has already
# checked out the branch/PR being tested, so do not clone main again here.
# LINALG_LIB=OpenBLAS + CMAKE_PREFIX_PATH point the ILP64 BLAS/LAPACK search at
# the OpenBLAS built above (LAPACK is bundled in libopenblas).
COPY . /opt/openqp
WORKDIR /opt/openqp
RUN cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 -DCMAKE_Fortran_COMPILER=gfortran-14 -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=ON -DLINALG_LIB=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/openblas
# NINJA_JOBS caps compile parallelism for memory-constrained builders (the
# Fortran modules are RAM-heavy); empty = use all cores (the CI default).
ARG NINJA_JOBS=
RUN ninja ${NINJA_JOBS:+-j${NINJA_JOBS}} -C build install
RUN cd pyoqp && pip3 install .

# Install dftd4
WORKDIR /opt
RUN wget https://github.com/dftd4/dftd4/releases/download/v3.3.0/dftd4-3.3.0-linux-x86_64.tar.xz
RUN tar -xf dftd4-3.3.0-linux-x86_64.tar.xz && rm dftd4-3.3.0-linux-x86_64.tar.xz
RUN mv dftd4-3.3.0 /opt/dftd4
ENV PATH="/opt/dftd4/bin:$PATH"
ENV DFTD4_BIN=/opt/dftd4/bin/dftd4

# Install dftd4 via pip
RUN pip3 install dftd4

# Set environment variables, initializing LD_LIBRARY_PATH if not defined
ENV OPENQP_ROOT=/opt/openqp
ENV OMP_NUM_THREADS=4

# Run a lightweight install smoke test.  The full example suite is covered by
# the regular CI workflow; Docker image builds should only check that the
# installed launcher can execute a bundled input, not gate on reference-sensitive
# example comparisons.
RUN openqp /opt/openqp/examples/other/h2o_rhf_6-31g_hf.inp

# Set entrypoint if required
ENTRYPOINT ["bash"]
