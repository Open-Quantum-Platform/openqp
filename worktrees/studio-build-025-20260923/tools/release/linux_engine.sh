#!/usr/bin/env bash
# AlmaLinux 8 / glibc 2.28, Python 3.11 and GCC 13, with MKL ILP64.
# Run in a task-owned container on an admitted Linux build lane.
set -euo pipefail
dnf install -y -q gcc-toolset-13-gcc gcc-toolset-13-gcc-c++ gcc-toolset-13-gcc-gfortran \
  python3.11 python3.11-devel python3.11-pip git make patch which zip
source /opt/rh/gcc-toolset-13/enable
export CMAKE_BUILD_PARALLEL_LEVEL=8
export CC=/opt/rh/gcc-toolset-13/root/usr/bin/gcc
export CXX=/opt/rh/gcc-toolset-13/root/usr/bin/g++
export FC=/opt/rh/gcc-toolset-13/root/usr/bin/gfortran
"$CC" --version
"$FC" --version
/usr/bin/python3.11 -m venv /out/engine-venv
PY=/out/engine-venv/bin/python
"$PY" -m pip install --upgrade pip
"$PY" -m pip install cmake ninja scikit-build-core pyinstaller==6.19.0 'numpy<2.2' \
  cffi mkl-devel==2026.1.0
# PyPI MKL runtime names carry a SONAME suffix; expose linker names only
# inside this isolated venv, without touching a system installation.
"$PY" - <<'PYLINK'
from pathlib import Path
for lib in Path('/out/engine-venv/lib').glob('libmkl*.so.*'):
    alias = lib.with_name(lib.name.split('.so.')[0] + '.so')
    if not alias.exists():
        alias.symlink_to(lib.name)
PYLINK
export MKLROOT=/out/engine-venv
export CMAKE_PREFIX_PATH="$MKLROOT"
export LD_LIBRARY_PATH="$MKLROOT/lib:${LD_LIBRARY_PATH:-}"
export OQP_EXTERNALS_ROOT=/cache
# The mounted checkout is owned by the host user, not container root.
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory GIT_CONFIG_VALUE_0=/src
export CMAKE_ARGS="-DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC -DLINALG_LIB=Intel10_64ilp_seq -DENABLE_MPI=OFF -DENABLE_DDX=OFF -DENABLE_OPENMP=ON -DUSE_LIBINT=OFF -DENABLE_OPENTRAH=OFF -DOQP_REUSE_EXTERNALS=ON -DOQP_EXTERNALS_ROOT=/cache"
"$PY" -m pip install -v /src --no-build-isolation
"$PY" -m pip install 'numpy<2.2' basis_set_exchange
"$PY" /studio/tools/release/freeze_engine.py --source /src --output /out/frozen/openqp
"$PY" /studio/tools/release/qualify_engine.py /out/frozen/openqp /out/qualification
