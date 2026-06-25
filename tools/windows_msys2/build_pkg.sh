#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
pkgbuild_dir="$repo_root/tools/windows_msys2/mingw-w64-openqp"
artifact_root="${OQP_MSYS2_ARTIFACT_ROOT:-$repo_root/dist/msys2-ucrt64}"
package_out="$artifact_root/packages"

if [[ "${MSYSTEM:-}" != "UCRT64" ]]; then
  echo "Run this script from an MSYS2 UCRT64 shell." >&2
  echo "Current MSYSTEM=${MSYSTEM:-unset}" >&2
  exit 2
fi

if [[ "${OQP_MSYS2_INSTALL_DEPS:-0}" == "1" ]]; then
  pacman --needed --noconfirm -S \
    base-devel \
    git \
    mingw-w64-ucrt-x86_64-gcc \
    mingw-w64-ucrt-x86_64-gcc-fortran \
    mingw-w64-ucrt-x86_64-python \
    mingw-w64-ucrt-x86_64-python-pip \
    mingw-w64-ucrt-x86_64-python-cffi \
    mingw-w64-ucrt-x86_64-python-numpy \
    mingw-w64-ucrt-x86_64-python-scipy \
    mingw-w64-ucrt-x86_64-python-build \
    mingw-w64-ucrt-x86_64-python-installer \
    mingw-w64-ucrt-x86_64-python-scikit-build-core \
    mingw-w64-ucrt-x86_64-python-setuptools \
    mingw-w64-ucrt-x86_64-python-wheel \
    mingw-w64-ucrt-x86_64-cmake \
    mingw-w64-ucrt-x86_64-ninja \
    mingw-w64-ucrt-x86_64-pkgconf
fi

# These two runtime dependencies are not currently available as MSYS2 packages.
# Install them only for local validation until separate PKGBUILDs are added.
if [[ "${OQP_MSYS2_INSTALL_PYPI_DEPS:-1}" == "1" ]]; then
  python -m pip install --user --break-system-packages --upgrade "jsonschema<4.18" basis_set_exchange geometric
fi

rm -rf "$package_out"
mkdir -p "$package_out"

if [[ -z "${OQP_MSYS2_EXTERNALS_ROOT:-}" ]]; then
  if [[ -n "${RUNNER_TEMP:-}" ]]; then
    export OQP_MSYS2_EXTERNALS_ROOT="$(cygpath -u "$RUNNER_TEMP")/oqp-ext"
  else
    export OQP_MSYS2_EXTERNALS_ROOT="${TMPDIR:-/tmp}/oqp-ext"
  fi
fi
rm -rf "$OQP_MSYS2_EXTERNALS_ROOT"
mkdir -p "$OQP_MSYS2_EXTERNALS_ROOT"

cd "$pkgbuild_dir"
# GitHub's Windows checkout may apply CRLF line endings depending on
# repository settings. makepkg sources these files directly and rejects CRLF.
sed -i 's/\r$//' PKGBUILD ./*.install
rm -f ./*.pkg.tar.*
makepkg-mingw --noconfirm --syncdeps --cleanbuild --clean

pkg_file="$(find "$pkgbuild_dir" -maxdepth 1 -type f -name 'mingw-w64-ucrt-x86_64-openqp-*.pkg.tar.*' | head -n 1)"
if [[ -z "$pkg_file" ]]; then
  echo "No MSYS2 package was created in $pkgbuild_dir" >&2
  exit 1
fi

cp "$pkg_file" "$package_out/"
pacman -U --noconfirm "$pkg_file"

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

openqp "$repo_root/examples/other/h2o_rhf_6-31g_hf.inp" --omp 2
echo "Created package artifact under $package_out"
