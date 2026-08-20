#!/usr/bin/env bash
# Build the self-contained Linux OpenQP bundle.  Runs INSIDE a manylinux_2_28
# container (see .github/workflows/linux-bundle.yml): linking against the
# runner's glibc would produce an archive that refuses to start on the RHEL and
# Rocky 8 machines a lot of clusters still run.
#
# Produces /src/dist/openqp-<version>.tar.gz containing a directory that runs
# with nothing installed:  ./openqp job.inp
set -euo pipefail

PY=/opt/python/cp311-cp311/bin/python
CACHE=/host-cache/openqp
BLAS="$CACHE/openblas64"
export XDG_CACHE_HOME=/host-cache

echo "::group::toolchain"
"$PY" --version
gcc --version | head -1
gfortran --version | head -1
echo "::endgroup::"

# --- ILP64 OpenBLAS -------------------------------------------------------
# The same BLAS the Linux wheels bundle: 8-byte integers, plain symbol names,
# runtime CPU dispatch so one archive serves every machine of this
# architecture.  MKL would add hundreds of megabytes for no portability gain.
if [ ! -e "$BLAS/lib/libopenblas.so" ]; then
  echo "::group::build OpenBLAS"
  git clone --depth 1 --branch v0.3.30 https://github.com/OpenMathLib/OpenBLAS.git /tmp/OpenBLAS
  make -C /tmp/OpenBLAS -j"$(nproc)" INTERFACE64=1 USE_OPENMP=1 DYNAMIC_ARCH=1 NO_STATIC=1 libs netlib shared
  make -C /tmp/OpenBLAS PREFIX="$BLAS" INTERFACE64=1 USE_OPENMP=1 DYNAMIC_ARCH=1 NO_STATIC=1 install
  echo "::endgroup::"
else
  echo "reusing cached OpenBLAS in $BLAS"
fi
export PKG_CONFIG_PATH="$BLAS/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="$BLAS/lib:${LD_LIBRARY_PATH:-}"

# --- OpenQP ---------------------------------------------------------------
echo "::group::build OpenQP"
cd /src
"$PY" -m pip install --upgrade pip
"$PY" -m pip install cmake ninja numpy cffi scikit-build-core pyinstaller
CMAKE_ARGS="-DCMAKE_FIND_ROOT_PATH_MODE_PACKAGE=BOTH -DLINALG_LIB=OpenBLAS -DENABLE_MPI=OFF -DUSE_LIBINT=OFF -DENABLE_OPENTRAH=OFF -DBUILD_TESTING=OFF" \
  "$PY" -m pip install -v . --no-build-isolation
"$PY" -m pip install scipy "numpy<2.2" basis_set_exchange geometric
echo "::endgroup::"

# `import oqp` prints "Failed to import mpi4py" on a non-MPI build, so take
# only the last line -- the path.
PKG=$("$PY" -c "import os, oqp; print(os.path.dirname(oqp.__file__))" | tail -n1)
VER=$("$PY" -c "import tomllib,pathlib; print(tomllib.loads(pathlib.Path('pyproject.toml').read_text())['project']['version'])")
echo "package: $PKG   version: $VER"

# --- freeze ---------------------------------------------------------------
echo "::group::freeze"
mkdir -p /tmp/oqpwork
cat > /tmp/oqpwork/openqp_entry.py <<'PY_ENTRY'
import sys
from oqp.pyoqp import main
if __name__ == "__main__":
    sys.exit(main())
PY_ENTRY

# liboqp and the DFT-D4 libraries are opened through cffi at run time and are
# invisible to PyInstaller's static analysis; add them by name.  OpenBLAS and
# the GCC Fortran runtime travel too -- the user has neither.
add=()
for so in "$PKG"/lib/*.so*; do add+=(--add-binary "$so:oqp/lib"); done
for so in "$BLAS"/lib/libopenblas*.so*; do
  [ -L "$so" ] || add+=(--add-binary "$so:oqp/lib")
done
for name in libgfortran.so.5 libquadmath.so.0 libgomp.so.1 libstdc++.so.6; do
  path=$(gcc -print-file-name="$name")
  if [ -f "$path" ]; then add+=(--add-binary "$(readlink -f "$path"):oqp/lib"); else echo "not found: $name"; fi
done
printf 'bundling: %s\n' "${add[@]}" | head -40

"$PY" -m PyInstaller --noconfirm --clean --onedir --name openqp \
  --distpath /tmp/oqpdist --workpath /tmp/oqpwork --specpath /tmp/oqpwork \
  --collect-data oqp --collect-submodules oqp \
  --collect-data basis_set_exchange --collect-submodules basis_set_exchange \
  --collect-data geometric --collect-submodules geometric \
  --collect-submodules scipy --collect-submodules numpy \
  --recursive-copy-metadata basis_set_exchange \
  --copy-metadata OpenQP --copy-metadata geometric \
  "${add[@]}" /tmp/oqpwork/openqp_entry.py
echo "::endgroup::"

ROOT=/tmp/oqpdist/openqp
cp -R "$PKG/share" "$ROOT/share"
cp -R /src/examples "$ROOT/examples"

cat > "$ROOT/README.txt" <<EOF
OpenQP for Linux -- standalone build
====================================

Nothing to install.  Unpack anywhere and run, for example:

    ./openqp examples/HF/H2O_RHF-HF_ENERGY.inp

Results are written next to the input file (.log, .json).

Threads: set OMP_NUM_THREADS to the number of physical cores, e.g.

    OMP_NUM_THREADS=8 ./openqp myjob.inp

This bundle carries its own Python, an ILP64 OpenBLAS with runtime CPU
dispatch, the GCC Fortran runtime and the DFT-D4 libraries.  It neither reads
nor needs any system installation of them.

Built inside manylinux_2_28, so it runs on any distribution with glibc 2.28 or
newer (RHEL/Rocky 8 and later, Ubuntu 18.10 and later).

OpenQP version : $VER
Architecture   : ${ARCH:-unknown}
Commit         : ${GITHUB_SHA:-unknown}
Built          : $(date -u '+%Y-%m-%d %H:%M:%S') UTC
Workflow run   : ${GITHUB_SERVER_URL:-}/${GITHUB_REPOSITORY:-}/actions/runs/${GITHUB_RUN_ID:-}
EOF
tail -6 "$ROOT/README.txt"

# --- smoke test: nothing on PATH, nothing preloaded -----------------------
echo "::group::smoke test"
mkdir -p /tmp/oqpsmoke
cp /src/examples/HF/H2O_RHF-HF_ENERGY.inp /tmp/oqpsmoke/hf.inp
cd /tmp/oqpsmoke
env -i PATH=/usr/bin:/bin HOME=/root OMP_NUM_THREADS=2 "$ROOT/openqp" hf.inp --nompi
grep "TOTAL energy" hf.log
grep -q -- "-76\.01074651" hf.log || { echo "energy does not match the reference"; exit 1; }
echo "::endgroup::"

mkdir -p /src/dist
tar -czf "/src/dist/openqp-${VER}.tar.gz" -C /tmp/oqpdist openqp
ls -lh "/src/dist/openqp-${VER}.tar.gz"
