# OpenQP build environment: every prerequisite of the OpenQP image build that
# does not depend on OpenQP sources (compiler toolchain, ILP64 OpenBLAS,
# Python dependencies). Published as openqp/openqp-buildenv:<tag> by
# .github/workflows/docker-build.yml, which falls back to building it locally
# when the registry image is missing. Bump the tag there and in the root
# Dockerfile FROM line together whenever this file changes.

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

# Python dependencies: cffi for the pyoqp build, dftd4 for the dispersion
# correction (pyoqp consumes dftd4 through its Python API only).
RUN pip3 install cffi dftd4

# Build OpenBLAS as ILP64 (INTERFACE64=1) with UNSUFFIXED symbols, so it exports
# `dgemm_` etc. with the 8-byte integer width OpenQP calls.  The distro's
# libopenblas64 instead exports `dgemm_64_`-suffixed symbols, which OpenQP does
# not call -- hence building from source here.  DYNAMIC_ARCH=1 keeps the image
# portable (OpenBLAS selects the kernel for the host CPU at run time, not the
# build machine's ISA).
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
