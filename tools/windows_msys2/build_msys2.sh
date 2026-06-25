#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
artifact_root="${OQP_MSYS2_ARTIFACT_ROOT:-$repo_root/dist/msys2-ucrt64}"
build_dir="${OQP_MSYS2_BUILD_DIR:-$repo_root/build/msys2-ucrt64}"
prefix="${OQP_MSYS2_PREFIX:-$artifact_root/prefix}"
wheelhouse="${OQP_MSYS2_WHEELHOUSE:-$artifact_root/wheelhouse}"
jobs="${OQP_MSYS2_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)}"

if [[ "${MSYSTEM:-}" != "UCRT64" ]]; then
  echo "Run this script from an MSYS2 UCRT64 shell." >&2
  echo "Current MSYSTEM=${MSYSTEM:-unset}" >&2
  exit 2
fi

if [[ "${OQP_MSYS2_INSTALL_DEPS:-0}" == "1" ]]; then
  pacman --needed --noconfirm -S \
    git \
    mingw-w64-ucrt-x86_64-gcc \
    mingw-w64-ucrt-x86_64-gcc-fortran \
    mingw-w64-ucrt-x86_64-python \
    mingw-w64-ucrt-x86_64-python-pip \
    mingw-w64-ucrt-x86_64-python-numpy \
    mingw-w64-ucrt-x86_64-python-scipy \
    mingw-w64-ucrt-x86_64-cmake \
    mingw-w64-ucrt-x86_64-ninja \
    mingw-w64-ucrt-x86_64-pkgconf \
    zip
fi

export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"
export FC="${FC:-gfortran}"
export CMAKE_BUILD_PARALLEL_LEVEL="$jobs"
export CMAKE_ARGS="${CMAKE_ARGS:-} -DUSE_LIBINT=OFF -DENABLE_OPENMP=ON -DENABLE_OPENTRAH=OFF -DLINALG_LIB=netlib -DLINALG_LIB_INT64=ON -DOQP_EXTERNALS_ROOT=$repo_root/.cache/openqp/externals"

python -m pip install --user --break-system-packages --upgrade cffi scikit-build-core build setuptools wheel "jsonschema<4.18" geometric basis_set_exchange

rm -rf "$build_dir" "$prefix" "$wheelhouse"
mkdir -p "$prefix" "$wheelhouse"

python -m build --wheel --no-isolation --skip-dependency-check --outdir "$wheelhouse" "$repo_root"

wheel_path="$(find "$wheelhouse" -maxdepth 1 -type f -name '*.whl' | head -n 1)"
if [[ -z "$wheel_path" ]]; then
  echo "No wheel was created in $wheelhouse" >&2
  exit 1
fi

python -m pip install --no-deps --prefix "$prefix" "$wheel_path"

site_packages="$(
  OQP_PREFIX="$prefix" python - <<'PY'
import os
import sysconfig

prefix = os.environ["OQP_PREFIX"]
print(sysconfig.get_path("purelib", vars={"base": prefix, "platbase": prefix}))
PY
)"

export PYTHONPATH="$site_packages${PYTHONPATH:+:$PYTHONPATH}"
export PATH="$prefix/bin:$PATH"

cat > "$prefix/activate-openqp.sh" <<EOF
# Source this from an MSYS2 UCRT64 shell before running OpenQP from this prefix.
export PYTHONPATH="$site_packages\${PYTHONPATH:+:\$PYTHONPATH}"
export PATH="$prefix/bin:\$PATH"
EOF

python - <<'PY'
import oqp
from pathlib import Path

root = Path(oqp.oqp_root)
required = [
    root / "include" / "oqp.h",
    root / "lib" / "liboqp.dll",
    root / "share" / "basis_sets",
]
missing = [str(path) for path in required if not path.exists()]
if missing:
    raise SystemExit("Missing OpenQP runtime files: " + ", ".join(missing))
print("OpenQP runtime root:", root)
PY

example="$repo_root/examples/other/h2o_rhf_6-31g_hf.inp"
openqp "$example" --omp 2

archive="$artifact_root/openqp-msys2-ucrt64-prefix.zip"
rm -f "$archive"
(cd "$artifact_root" && zip -qr "$archive" "$(basename "$prefix")" wheelhouse)
echo "Created $archive"
echo "To use this install later, run: source $prefix/activate-openqp.sh"
