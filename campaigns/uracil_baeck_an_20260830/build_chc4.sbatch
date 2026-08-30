#!/bin/bash
#SBATCH --job-name=oqp-ba-build
#SBATCH --partition=chc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00

set -euo pipefail

campaign_root=${OQP_BA_CAMPAIGN_ROOT:?set OQP_BA_CAMPAIGN_ROOT}
source_commit=${OQP_SOURCE_COMMIT:?set OQP_SOURCE_COMMIT}
bundle_path=${OQP_SOURCE_BUNDLE:?set OQP_SOURCE_BUNDLE}
source_root="$campaign_root/source/openqp-$source_commit"
attempt_root="$campaign_root/build/attempt-001"
external_root=/bighome/cheolho.choi/.cache/openqp/externals
gcc_root=/opt/soft/install/GCCcore/12.3.0
mkl_root=/opt/soft/install/imkl/2023.1.0/mkl/2023.1.0
python_root=/opt/soft/install/Python/3.11.3-GCCcore-12.3.0

test ! -e "$attempt_root"
mkdir -p "$campaign_root/build/logs" "$attempt_root" "$external_root"

export MODULEPATH=/etc/scl/modulefiles:/etc/modulefiles:/usr/share/modulefiles:/opt/nvidia/hpc_sdk/modulefiles:/opt/soft/modulefiles/all
set +u
source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load GCC/12.3.0 Python/3.11.3-GCCcore-12.3.0 imkl/2023.1.0
set -u

export PATH="$gcc_root/bin:$python_root/bin:$PATH"
export LD_LIBRARY_PATH="$gcc_root/lib64:$python_root/lib:$mkl_root/lib/intel64:${LD_LIBRARY_PATH:-}"
export MKLROOT="$mkl_root"
export MKL_INTERFACE_LAYER=ILP64
export MKL_THREADING_LAYER=GNU
export CMAKE_BUILD_PARALLEL_LEVEL="$SLURM_CPUS_PER_TASK"
export CFLAGS="${CFLAGS:-} -gdwarf-2"
export CXXFLAGS="${CXXFLAGS:-} -gdwarf-2"
export TMPDIR="$attempt_root/tmp"
mkdir -p "$TMPDIR" "$attempt_root/wheelhouse"

{
  date -Ins
  hostname
  uname -a
  lscpu
  printf 'SLURM_JOB_ID=%s\nSLURM_CPUS_PER_TASK=%s\n' "$SLURM_JOB_ID" "$SLURM_CPUS_PER_TASK"
  "$gcc_root/bin/gcc" --version | head -n 1
  "$gcc_root/bin/gfortran" --version | head -n 1
  "$python_root/bin/python3.11" -V
  sha256sum "$bundle_path"
  git init --bare "$attempt_root/bundle-verify.git"
  git -C "$attempt_root/bundle-verify.git" bundle verify "$bundle_path"
  ls -l "$mkl_root/lib/intel64/libmkl_rt.so"
} > "$attempt_root/machine-toolchain.txt" 2>&1

git clone "$bundle_path" "$source_root"
git -C "$source_root" checkout --detach "$source_commit"
test "$(git -C "$source_root" rev-parse HEAD)" = "$source_commit"
test -z "$(git -C "$source_root" status --porcelain)"

"$python_root/bin/python3.11" -m venv "$attempt_root/venv"
"$attempt_root/venv/bin/python" -m pip install --upgrade pip
cd "$source_root"
export CMAKE_ARGS="-DCMAKE_C_COMPILER=$gcc_root/bin/gcc -DCMAKE_CXX_COMPILER=$gcc_root/bin/g++ -DCMAKE_Fortran_COMPILER=$gcc_root/bin/gfortran -DUSE_LIBINT=OFF -DENABLE_OPENMP=ON -DLINALG_LIB=Intel10_64ilp_seq -DENABLE_OPENTRAH=OFF -DOQP_REUSE_EXTERNALS=ON -DOQP_EXTERNALS_ROOT=$external_root -DMKLROOT=$mkl_root -DCMAKE_C_FLAGS=-gdwarf-2 -DCMAKE_CXX_FLAGS=-gdwarf-2 -DCMAKE_Fortran_FLAGS=-gdwarf-2"
"$attempt_root/venv/bin/python" -m pip wheel . --no-deps \
  --wheel-dir "$attempt_root/wheelhouse" > "$attempt_root/build.log" 2>&1

wheel=$(find "$attempt_root/wheelhouse" -maxdepth 1 -type f -name 'openqp-*.whl' -print -quit)
test -n "$wheel"
"$attempt_root/venv/bin/python" -m pip install \
  'numpy<2.2' scipy cffi basis_set_exchange geometric pytest
"$attempt_root/venv/bin/python" -m pip install --no-deps --force-reinstall "$wheel"

if grep -Eqi 'Building external.*libint|libint2-build' "$attempt_root/build.log"; then
  echo 'FAIL: Libint build detected' >&2
  exit 1
fi
if grep -Eq 'OQP_BLAS_INT[= :]4|blasint4' "$attempt_root/build.log"; then
  echo 'FAIL: LP64 BLAS integer width detected' >&2
  exit 1
fi
package_root=$("$attempt_root/venv/bin/python" -c \
  'import importlib.util; print(next(iter(importlib.util.find_spec("oqp").submodule_search_locations)))')
native_lib=$(find "$package_root" -type f -name 'liboqp*.so' -print -quit)
test -n "$native_lib"
ldd "$native_lib" > "$attempt_root/ldd.txt"
if grep -Eq 'libmkl_(intel|gf)_lp64' "$attempt_root/ldd.txt"; then
  echo 'FAIL: LP64 MKL interface linked' >&2
  exit 1
fi

export OPENQP_ROOT="$package_root"
"$attempt_root/venv/bin/python" -m pytest -q \
  tests/test_namd_baeck_an_provider.py \
  tests/test_namd_baeck_an.py \
  tests/test_namd_analytic_directional.py \
  tests/test_oqp_input.py > "$attempt_root/regression.log" 2>&1
if grep -Eq '[0-9]+ failed|[0-9]+ skipped' "$attempt_root/regression.log"; then
  cat "$attempt_root/regression.log" >&2
  exit 1
fi

{
  printf 'status=verified\nsource_commit=%s\nwheel=%s\n' "$source_commit" "$wheel"
  sha256sum "$wheel" "$bundle_path" "$attempt_root/build.log" \
    "$attempt_root/regression.log" "$native_lib"
  grep -E 'Intel10_64ilp|MKL|blasint8|Reusing cached external|Building external' \
    "$attempt_root/build.log" | tail -n 160
  date -Ins
} > "$attempt_root/VERIFIED"

echo "BUILD_VERIFIED $attempt_root"
