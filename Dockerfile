# OpenQP image: compile the checked-out sources on top of the prebuilt build
# environment (compiler toolchain, ILP64 OpenBLAS, Python dependencies -- see
# docker/buildenv.Dockerfile), so image builds do not reinstall prerequisites.
# Keep the tag below in sync with .github/workflows/docker-build.yml and
# docker/buildenv.Dockerfile.
FROM openqp/openqp-buildenv:1

# Copy and compile the checked-out OpenQP source.  GitHub Actions has already
# checked out the branch/PR being tested, so do not clone main again here.
# LINALG_LIB=OpenBLAS + CMAKE_PREFIX_PATH point the ILP64 BLAS/LAPACK search at
# the OpenBLAS in the build environment (LAPACK is bundled in libopenblas).
COPY . /opt/openqp
WORKDIR /opt/openqp
RUN cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 -DCMAKE_Fortran_COMPILER=gfortran-14 -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=ON -DLINALG_LIB=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/openblas
# NINJA_JOBS caps compile parallelism for memory-constrained builders (the
# Fortran modules are RAM-heavy); empty = use all cores (the CI default).
ARG NINJA_JOBS=
RUN ninja ${NINJA_JOBS:+-j${NINJA_JOBS}} -C build install
RUN cd pyoqp && pip3 install .

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
