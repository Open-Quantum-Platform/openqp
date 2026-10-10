#!/usr/bin/env bash
set -euo pipefail
mkdir -p reports
{
  printf 'Source: %s\n' "${CI_COMMIT_SHA:-local}"
  nvcc --version
  cmake --version
  c++ --version
} > reports/cuda-toolchain.txt

# This standalone library currently passes C int to host BLAS/LAPACK.
# Its ABI is distinct from the OpenQP engine's mandatory ILP64 ABI.
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_INSTALL_PREFIX="$PWD/build/install" \
  -DOPENQP_GPU_METC=ON -DOPENQP_GPU_METC_ONLY=OFF \
  -DOPENQP_GPU_DFT=ON -DOPENQP_GPU_GRAD=ON -DOPENQP_GPU_DFBUILD=ON \
  -DBLA_VENDOR=OpenBLAS -DBLA_SIZEOF_INTEGER=4 \
  2>&1 | tee reports/cuda-configure.log
cmake --build build --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-2}" \
  2>&1 | tee reports/cuda-build.log
cmake --install build 2>&1 | tee reports/cuda-install.log

for name in libopenqp_gpu.so libopenqp_gpu_metc.so libopenqp_gpu_df.so; do
  test -s "build/$name"
  ldd "build/$name" > "reports/$name.linkage.txt"
  if grep -q 'not found' "reports/$name.linkage.txt"; then
    cat "reports/$name.linkage.txt"
    exit 1
  fi
done
test -x build/openqp_gpu_build_df
(cd build && sha256sum libopenqp_gpu.so libopenqp_gpu_metc.so libopenqp_gpu_df.so) \
  > reports/cuda-libraries.sha256

# Also exercise the documented BLAS-free standalone configuration.
cmake -S . -B build/metc-only -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DOPENQP_GPU_METC_ONLY=ON 2>&1 | tee reports/metc-only-configure.log
cmake --build build/metc-only --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-2}" \
  2>&1 | tee reports/metc-only-build.log
