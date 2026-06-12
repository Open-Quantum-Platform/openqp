# OpenQP image: install the checked-out sources with the standard
# `pip install .` route (scikit-build-core builds the native library and the
# Python package together) on top of the prebuilt build environment (compiler
# toolchain, ILP64 OpenBLAS, Python dependencies -- see
# docker/buildenv.Dockerfile), so image builds do not reinstall prerequisites.
# Keep the tag below in sync with .github/workflows/docker-build.yml and
# docker/buildenv.Dockerfile.
FROM openqp/openqp-buildenv:1

# Copy and install the checked-out OpenQP source.  GitHub Actions has already
# checked out the branch/PR being tested, so do not clone main again here.
# USE_LIBINT=OFF, ENABLE_OPENMP=ON, and LINALG_LIB_INT64=ON are the pyproject
# defaults; CC/CXX/FC select the gcc-14 toolchain and CMAKE_ARGS points the
# ILP64 BLAS/LAPACK search at the OpenBLAS in the build environment (LAPACK is
# bundled in libopenblas).  OPENQP_ROOT is intentionally not set: the
# pip-installed package locates itself, and pointing it at the source tree
# would be wrong.
COPY . /opt/openqp
WORKDIR /opt/openqp
ENV CC=gcc-14 CXX=g++-14 FC=gfortran-14
ENV CMAKE_ARGS="-DLINALG_LIB=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/openblas"
# NINJA_JOBS caps compile parallelism for memory-constrained builders (the
# Fortran modules are RAM-heavy); empty = use all cores (the CI default).
ARG NINJA_JOBS=
RUN ${NINJA_JOBS:+env CMAKE_BUILD_PARALLEL_LEVEL=${NINJA_JOBS}} pip3 install .

ENV OMP_NUM_THREADS=4

# Run a lightweight install smoke test.  The full example suite is covered by
# the regular CI workflow; Docker image builds should only check that the
# installed launcher can execute a bundled input, not gate on reference-sensitive
# example comparisons.
RUN openqp /opt/openqp/examples/other/h2o_rhf_6-31g_hf.inp

# Set entrypoint if required
ENTRYPOINT ["bash"]
