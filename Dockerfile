# syntax=docker/dockerfile:1.7
#
# Release-candidate container for OpenQP.  Both inputs are immutable registry
# digests recorded in docker/base-images.lock.json.  The build environment is
# used only to compile and verify a wheel; no compiler, source checkout, or
# external-build cache is inherited by the final runtime image.
#
# The current openqp-buildenv:1 image is Ubuntu 24.04/Python 3.12 and is
# linux/amd64-only.  The runtime therefore deliberately uses the same Python
# minor version and a newer glibc baseline.  Moving the container to Python
# 3.11 requires publishing the proposed Python-3.11 builder first; see
# docs/openqp-dev-buildenv-hardening-plan.md.
FROM openqp/openqp-buildenv:1@sha256:83a2ba2108bc2bb1123e7dfac31320833817242681ee41f36e5ee950dfb22a3e AS builder

ARG OPENQP_VERSION=1.3.0
ARG OPENQP_REVISION

ENV CC=gcc-14 \
    CXX=g++-14 \
    FC=gfortran-14 \
    CMAKE_ARGS="-DUSE_LIBINT=OFF -DENABLE_OPENMP=ON -DENABLE_OPENTRAH=OFF -DENABLE_MPI=OFF -DENABLE_DDX=OFF -DENABLE_OPENQP_DFTB=OFF -DOQP_ALLOW_REFERENCE_BLAS=OFF -DLINALG_LIB=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/openblas -DOQP_EXTERNALS_ROOT=/root/.cache/openqp/externals -DOQP_SOURCE_REVISION=${OPENQP_REVISION}"
ARG NINJA_JOBS=

WORKDIR /opt/openqp

# Download exact source-independent build and runtime wheelhouses. Build tools
# are installed only from the recorded local wheels, so the candidate inventory
# identifies the exact distributions and hashes that created the OpenQP wheel.
# The final stage has no networked package installation.
COPY pyproject.toml /tmp/openqp-pyproject.toml
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    python3 - <<'PY'
import subprocess
import sys
import tomllib
from pathlib import Path

project = tomllib.loads(Path("/tmp/openqp-pyproject.toml").read_text())
build_requirements = list(project.get("build-system", {}).get("requires", []))
runtime_requirements = list(project.get("project", {}).get("dependencies", []))

Path("/opt/openqp-build-wheelhouse").mkdir(parents=True, exist_ok=True)
subprocess.check_call(
    [
        sys.executable,
        "-m",
        "pip",
        "download",
        "--only-binary=:all:",
        "--dest=/opt/openqp-build-wheelhouse",
        *dict.fromkeys(build_requirements),
    ]
)
subprocess.check_call(
    [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--force-reinstall",
        "--no-index",
        "--find-links=/opt/openqp-build-wheelhouse",
        *dict.fromkeys(build_requirements),
    ]
)
Path("/opt/openqp-wheelhouse").mkdir(parents=True, exist_ok=True)
subprocess.check_call(
    [
        sys.executable,
        "-m",
        "pip",
        "download",
        "--only-binary=:all:",
        "--dest=/opt/openqp-wheelhouse",
        *dict.fromkeys(runtime_requirements),
    ]
)
PY

# Build exactly one wheel from the checked-out source.  BuildKit cache mounts
# are ephemeral and never copied into the runtime image.
COPY . /opt/openqp
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    --mount=type=cache,target=/root/.cache/openqp/externals,sharing=locked \
    ${NINJA_JOBS:+env CMAKE_BUILD_PARALLEL_LEVEL=${NINJA_JOBS}} \
    python3 -m pip wheel --no-build-isolation --no-deps \
      --wheel-dir=/opt/openqp-wheelhouse .

RUN python3 .github/scripts/record_wheelhouse.py /opt/openqp-wheelhouse \
      --build-wheelhouse=/opt/openqp-build-wheelhouse \
      --expected-openqp-version="${OPENQP_VERSION}" \
      --output=/opt/openqp-wheelhouse/wheelhouse-manifest.json

# Install from the offline wheelhouse in the builder and run the same strong
# numerical/ABI/D4 canonical-replaceability gate used for release wheels.
RUN python3 -m pip install --no-index --find-links=/opt/openqp-wheelhouse \
      "OpenQP==${OPENQP_VERSION}" && \
    OQP_WHEEL_SMOKE_EXTERNAL_RUNTIME_PATH=/opt/openblas/lib \
      python3 .github/scripts/wheel_smoke_test.py /opt/openqp

# Copy only the non-glibc shared-library closure required by the installed
# OpenQP wheel.  The collector fails on unresolved dependencies or colliding
# SONAMEs and records package versions, hashes, and license files.
RUN python3 .github/scripts/collect_docker_runtime.py \
      --output=/opt/openqp-runtime \
      --base-lock=docker/base-images.lock.json


FROM python:3.12.13-slim-trixie@sha256:229a2c5bfa27522db7815ea81f9bed70af17ccb9de9fc7ad142b1877b5830d36 AS runtime

ARG OPENQP_VERSION=1.3.0
ARG OPENQP_REVISION
LABEL org.opencontainers.image.title="OpenQP" \
      org.opencontainers.image.description="Open Quantum Platform research software" \
      org.opencontainers.image.source="https://github.com/Open-Quantum-Platform/openqp" \
      org.opencontainers.image.documentation="https://open-quantum-platform.github.io/openqp-docs/" \
      org.opencontainers.image.version="${OPENQP_VERSION}" \
      org.opencontainers.image.revision="${OPENQP_REVISION}" \
      org.opencontainers.image.licenses="LicenseRef-OpenQP-Research-1.0"

ENV OMP_NUM_THREADS=4 \
    LD_LIBRARY_PATH=/opt/openblas/lib:/opt/openqp-runtime/lib \
    HOME=/work \
    XDG_CACHE_HOME=/work/.cache

# The builder-produced wheelhouse is mounted read-only for this instruction;
# it is not copied into an image layer.  No compiler or package index is used.
RUN --mount=type=bind,from=builder,source=/opt/openqp-wheelhouse,target=/tmp/wheelhouse,ro \
    python -m pip install --no-cache-dir --no-index \
      --find-links=/tmp/wheelhouse "OpenQP==${OPENQP_VERSION}"

# MRSF MINRES callbacks are module procedures, so packaged OpenQP libraries do
# not require stack trampolines. Verify the installed wheel without rewriting
# its ELF headers or wheel RECORD.
RUN --mount=type=bind,from=builder,source=/opt/openqp/.github/scripts/normalize_elf_stack.py,target=/tmp/normalize_elf_stack.py,ro \
    python /tmp/normalize_elf_stack.py --openqp-package

COPY --from=builder /opt/openqp-runtime/openblas/ /opt/openblas/lib/
COPY --from=builder /opt/openqp-runtime/lib/ /opt/openqp-runtime/lib/
COPY --from=builder /opt/openqp-runtime/licenses/ /usr/share/licenses/openqp/runtime-packages/
COPY --from=builder /opt/openqp-runtime/runtime-library-manifest.json /usr/share/licenses/openqp/
COPY --from=builder /opt/openqp-wheelhouse/wheelhouse-manifest.json /usr/share/licenses/openqp/
COPY LICENSE LICENSING.md SUSTAINABILITY.md THIRD_PARTY_NOTICES.md /usr/share/licenses/openqp/
COPY docker/THIRD_PARTY_RUNTIME.md /usr/share/licenses/openqp/
COPY licenses/third_party/ /usr/share/licenses/openqp/third_party/

# Verify the actual minimal final stage, including an end-to-end calculation,
# exact DFT-D4 load paths and removal tests, corresponding source, shared
# dependency closure, legal files, and absence of build tools/static archives.
RUN --mount=type=bind,from=builder,source=/opt/openqp/.github/scripts/container_runtime_smoke.py,target=/tmp/container_runtime_smoke.py,ro \
    OPENQP_EXPECTED_REVISION="${OPENQP_REVISION}" \
      python /tmp/container_runtime_smoke.py

RUN mkdir -p /work/.cache && chown -R 65532:65532 /work
WORKDIR /work
USER 65532:65532

# Preserve the historical interactive-container interface.  Users may invoke
# OpenQP inside the shell or override the entrypoint for batch execution.
ENTRYPOINT ["bash"]
