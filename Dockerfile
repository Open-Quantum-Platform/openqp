# syntax=docker/dockerfile:1.7
# OpenQP image: install the checked-out sources with the standard
# `pip install .` route (scikit-build-core builds the native library and the
# Python package together) on top of the prebuilt build environment (compiler
# toolchain, ILP64 OpenBLAS, Python dependencies -- see
# docker/buildenv.Dockerfile), so image builds do not reinstall prerequisites.
# Keep the tag below in sync with .github/workflows/docker-build.yml and
# docker/buildenv.Dockerfile.
FROM openqp/openqp-buildenv:1

ARG OPENQP_VERSION=dev
ARG OPENQP_REVISION=unknown
LABEL org.opencontainers.image.title="OpenQP" \
      org.opencontainers.image.description="Open Quantum Platform research software" \
      org.opencontainers.image.source="https://github.com/Open-Quantum-Platform/openqp" \
      org.opencontainers.image.documentation="https://open-quantum-platform.github.io/openqp-docs/" \
      org.opencontainers.image.version="${OPENQP_VERSION}" \
      org.opencontainers.image.revision="${OPENQP_REVISION}" \
      org.opencontainers.image.licenses="LicenseRef-OpenQP-Research-1.0"

# Install Python build/runtime dependencies in a source-independent layer.
# Changes to pyproject.toml invalidate this layer automatically; source-only
# PRs then reuse it instead of re-downloading scipy, basis_set_exchange, etc.
COPY pyproject.toml /tmp/openqp-pyproject.toml
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    python3 - <<'PY'
import subprocess
import sys
import tomllib

with open("/tmp/openqp-pyproject.toml", "rb") as handle:
    project = tomllib.load(handle)

requirements = []
requirements.extend(project.get("build-system", {}).get("requires", []))
requirements.extend(project.get("project", {}).get("dependencies", []))

deduped = []
for requirement in requirements:
    if requirement not in deduped:
        deduped.append(requirement)

subprocess.check_call([sys.executable, "-m", "pip", "install", *deduped])
PY

# Copy and install the checked-out OpenQP source.  GitHub Actions has already
# checked out the branch/PR being tested, so do not clone main again here.
# USE_LIBINT=OFF and ENABLE_OPENMP=ON are the pyproject
# defaults; CC/CXX/FC select the gcc-14 toolchain and CMAKE_ARGS points the
# ILP64 BLAS/LAPACK search at the OpenBLAS in the build environment (LAPACK is
# bundled in libopenblas).  OPENQP_ROOT is intentionally not set: the
# pip-installed package locates itself, and pointing it at the source tree
# would be wrong.
COPY . /opt/openqp
COPY LICENSE LICENSING.md SUSTAINABILITY.md THIRD_PARTY_NOTICES.md /usr/share/licenses/openqp/
COPY licenses/third_party/ /usr/share/licenses/openqp/third_party/
WORKDIR /opt/openqp
ENV CC=gcc-14 CXX=g++-14 FC=gfortran-14
ENV CMAKE_ARGS="-DLINALG_LIB=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/openblas -DOQP_EXTERNALS_ROOT=/root/.cache/openqp/externals"
# NINJA_JOBS caps compile parallelism for memory-constrained builders (the
# Fortran modules are RAM-heavy); empty = use all cores (the CI default).
ARG NINJA_JOBS=
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    --mount=type=cache,target=/root/.cache/openqp/externals,sharing=locked \
    ${NINJA_JOBS:+env CMAKE_BUILD_PARALLEL_LEVEL=${NINJA_JOBS}} \
    pip3 install --no-build-isolation --no-deps .

ENV OMP_NUM_THREADS=4

# Run a lightweight install smoke test.  The full example suite is covered by
# the regular CI workflow; Docker image builds should only check that the
# installed launcher can execute a bundled input, not gate on reference-sensitive
# example comparisons.
RUN openqp /opt/openqp/examples/other/h2o_rhf_6-31g_hf.inp

# The external-build cache mount above is intentionally ephemeral.  Verify that
# the installed application itself retained the replaceable DFT-D4 libraries
# and their complete corresponding source before this image can be published.
RUN python3 - <<'PY'
from pathlib import Path
import sys

from oqp.runtime import resolve_oqp_root

root, _ = resolve_oqp_root()
root = Path(root)
if sys.platform == "darwin":
    libraries = [
        "libdftd4.3.dylib", "libmulticharge.0.dylib", "libmctc-lib.0.dylib"
    ]
else:
    libraries = ["libdftd4.so.3", "libmulticharge.so.0", "libmctc-lib.so.0"]
required = [root / "lib" / name for name in libraries]
source = root / "share" / "corresponding-source" / "dftd4-stack"
required.extend([
    source / "README.md",
    source / "openqp-external-build.cmake",
    source / "mctc-lib-0.4.2" / "LICENSE",
    source / "multicharge-0.3.0" / "LICENSE",
    source / "dftd4-3.7.0" / "COPYING",
    source / "dftd4-3.7.0" / "COPYING.LESSER",
])
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise SystemExit(f"DFT-D4 runtime/corresponding-source files missing: {missing}")
print(f"DFT-D4 shared libraries and corresponding source retained under {root}")
PY

# Set entrypoint if required
ENTRYPOINT ["bash"]
