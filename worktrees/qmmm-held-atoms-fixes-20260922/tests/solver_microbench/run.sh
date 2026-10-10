#!/usr/bin/env bash
# Standalone micro-benchmark + validation for the z-vector iterative solvers.
#
# It compiles the REAL source/pcg.F90 and source/minres.F90 against a synthetic
# symmetric system (no BLAS/LAPACK/Libint, so it runs anywhere gfortran exists)
# and reports:
#   * PERF   : old vs new pcg_step wall time on a cheap-matvec SPD system
#   * STAB   : fail-closed behaviour on NaN/Inf injection
#   * MINRES : correctness on SPD (matches PCG) + robustness on indefinite A
#              (incl. the textbook A=diag(+1,-1) case where CG breaks down)
#
# Usage:  cd tests/solver_microbench && ./run.sh
set -euo pipefail
cd "$(dirname "$0")"

FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -ffree-line-length-none}
ROOT=$(git rev-parse --show-toplevel)

# Baseline pcg.F90 = the commit just before the "Carry PCG rz ..." optimization.
OLD_REF=${OLD_REF:-f601f83}

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

cp "$ROOT/source/pcg.F90"    "$work/pcg_new.F90"
cp "$ROOT/source/minres.F90" "$work/minres.F90"
git -C "$ROOT" show "$OLD_REF:source/pcg.F90" > "$work/pcg_old.F90"
cp precision.F90 messages.F90 bench_common.F90 perf_driver.F90 minres_check.F90 "$work/"

build() {  # $1 = variant tag, $2 = pcg file, $3 = main, rest = extra srcs
  local tag=$1 pcg=$2 main=$3; shift 3
  local d="$work/build_$tag"; mkdir -p "$d"; ( cd "$d"
    $FC $FFLAGS -c ../precision.F90 ../messages.F90
    $FC $FFLAGS -c "../$pcg"
    for s in "$@"; do $FC $FFLAGS -c "../$s"; done
    $FC $FFLAGS -c "../$main"
    $FC $FFLAGS -o run ./*.o )
}

echo "== building =="
build perfold pcg_old.F90 perf_driver.F90 bench_common.F90
build perfnew pcg_new.F90 perf_driver.F90 bench_common.F90
build minres  pcg_new.F90 minres_check.F90 minres.F90 bench_common.F90

echo; echo "================= PERFORMANCE: OLD pcg.F90 ================="
"$work/build_perfold/run"
echo; echo "================= PERFORMANCE: NEW pcg.F90 ================="
"$work/build_perfnew/run"
echo; echo "================= MINRES VALIDATION ======================="
"$work/build_minres/run"
