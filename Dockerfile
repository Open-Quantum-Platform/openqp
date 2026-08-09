# syntax=docker/dockerfile:1.7
# OpenQP image: install the checked-out sources with the standard
# `pip install .` route (scikit-build-core builds the native library and the
# Python package together) on top of the prebuilt build environment (compiler
# toolchain, ILP64 OpenBLAS, Python dependencies -- see
# docker/buildenv.Dockerfile), so image builds do not reinstall prerequisites.
# Keep the tag below in sync with .github/workflows/docker-build.yml and
# docker/buildenv.Dockerfile.
FROM openqp/openqp-buildenv:1

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
import hashlib
import json
from pathlib import Path
import sys

from oqp.runtime import resolve_oqp_root

root, _ = resolve_oqp_root()
root = Path(root)
if sys.platform == "darwin":
    libraries = {
        "dftd4": "libdftd4.3.dylib",
        "multicharge": "libmulticharge.0.dylib",
        "mctc-lib": "libmctc-lib.0.dylib",
    }
else:
    libraries = {
        "dftd4": "libdftd4.so.3",
        "multicharge": "libmulticharge.so.0",
        "mctc-lib": "libmctc-lib.so.0",
    }
required = [root / "lib" / name for name in libraries.values()]
source = root / "share" / "corresponding-source" / "dftd4-stack"
required.extend([
    source / "README.md",
    source / "BUILD-INFO.json",
    source / "apply-patch.cmake",
    source / "generate-build-info.cmake",
    source / "openqp-external-build.cmake",
    source / "patches" / "mctc-lib-0.4.2-disable-tests.patch",
    source / "patches" / "dftd4-3.7.0-disable-tests-and-mstore.patch",
    source / "mctc-lib-0.4.2" / "LICENSE",
    source / "multicharge-0.3.0" / "LICENSE",
    source / "dftd4-3.7.0" / "COPYING",
    source / "dftd4-3.7.0" / "COPYING.LESSER",
])
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise SystemExit(f"DFT-D4 runtime/corresponding-source files missing: {missing}")
build_info_path = source / "BUILD-INFO.json"
build_info_text = build_info_path.read_text(encoding="utf-8")
build_info = json.loads(build_info_text)
if (
    build_info.get("schema") != "org.open-quantum-platform.dftd4-build-info"
    or build_info.get("schema_version") != 1
    or build_info.get("canonical_runtime_names") != libraries
):
    raise SystemExit("DFT-D4 BUILD-INFO.json schema/runtime names are invalid")
components = {item["name"]: item for item in build_info.get("components", [])}
expected_components = {
    "mctc-lib": (
        "0.4.2",
        "https://github.com/grimme-lab/mctc-lib/archive/refs/tags/v0.4.2.tar.gz",
        "c7aa45c0a3e6f96e3316e15fc6cdbe48b15234940d3773927a57bb7bfe9744ac",
        "Apache-2.0",
    ),
    "multicharge": (
        "0.3.0",
        "https://github.com/grimme-lab/multicharge/archive/refs/tags/v0.3.0.tar.gz",
        "2fcc1f80871f404f005e9db458ffaec95bb28a19516a0245278cd3175b63a6b2",
        "Apache-2.0",
    ),
    "dftd4": (
        "3.7.0",
        "https://github.com/dftd4/dftd4/archive/refs/tags/v3.7.0.tar.gz",
        "f00b244759eff2c4f54b80a40673440ce951b6ddfa5eee1f46124297e056f69c",
        "LGPL-3.0-or-later",
    ),
}
for name, expected in expected_components.items():
    component = components.get(name, {})
    actual = tuple(
        component.get(field)
        for field in ("version", "archive_url", "sha256", "license")
    )
    if actual != expected:
        raise SystemExit(f"DFT-D4 BUILD-INFO component is invalid: {name}")
if len(build_info.get("components", [])) != len(expected_components):
    raise SystemExit("DFT-D4 BUILD-INFO has duplicate/extra components")
build = build_info.get("build", {})
blas = build.get("blas", {})
if (
    not build.get("cmake_version")
    or not build.get("generator")
    or not build.get("system", {}).get("name")
    or not build.get("system", {}).get("processor")
    or build.get("build_type") != "Release"
    or build.get("build_shared_libs") is not True
    or not isinstance(build.get("openmp"), bool)
    or not blas.get("requested_provider")
    or not blas.get("resolved_provider")
    or blas.get("integer_bytes") != 8
    or not blas.get("resolved_blas_libraries")
    or not blas.get("resolved_lapack_libraries")
):
    raise SystemExit("DFT-D4 BUILD-INFO.json resolved build data are incomplete")
for library_kind in ("resolved_blas_libraries", "resolved_lapack_libraries"):
    if any("/" in name or "\\" in name for name in blas[library_kind]):
        raise SystemExit("Absolute BLAS/LAPACK path leaked into DFT-D4 BUILD-INFO")
for language in ("c", "cxx", "fortran"):
    compiler = build.get("compilers", {}).get(language, {})
    if (
        not compiler.get("id")
        or not compiler.get("version")
        or compiler.get("executable") != Path(compiler.get("executable", "")).name
    ):
        raise SystemExit(f"DFT-D4 compiler record is invalid: {language}")
forwarded_flags = build.get("forwarded_flags", {})
if set(forwarded_flags) != {
    "c", "c_release", "fortran", "fortran_release"
} or not all(isinstance(value, str) for value in forwarded_flags.values()):
    raise SystemExit("DFT-D4 forwarded flag record is incomplete")
patches = {item["file"]: item for item in build_info.get("patches", [])}
expected_patches = {
    "mctc-lib-0.4.2-disable-tests.patch": "mctc-lib",
    "dftd4-3.7.0-disable-tests-and-mstore.patch": "dftd4",
}
if len(build_info.get("patches", [])) != len(expected_patches):
    raise SystemExit("DFT-D4 BUILD-INFO has duplicate/extra patches")
for patch_name, component in expected_patches.items():
    patch_path = source / "patches" / patch_name
    digest = hashlib.sha256(patch_path.read_bytes()).hexdigest()
    if patches.get(patch_name) != {
        "component": component,
        "file": patch_name,
        "sha256": digest,
        "strip": 1,
    }:
        raise SystemExit(f"DFT-D4 patch record is invalid: {patch_name}")
for forbidden in (
    "/private/tmp/", "/tmp/", "/private/var/folders/", "/root/.cache/",
    "/.cache/openqp/", "/host-cache/", "/Library/Caches/openqp/", "/opt/openqp",
):
    if forbidden in build_info_text:
        raise SystemExit(f"Transient path leaked into DFT-D4 BUILD-INFO: {forbidden}")
print(f"DFT-D4 shared libraries and corresponding source retained under {root}")
PY

# Set entrypoint if required
ENTRYPOINT ["bash"]
