// Fused spherical-harmonic rotation (Kronecker-apply) kernel -- CUDA C++.
//
// Standalone, nvcc-buildable port of the validated CuPy RawKernel
// (spherical_eri/cuda_kernel.py).  One thread block per shell quartet; the
// reference block and the four R^l matrices are staged in shared memory and the
// four contractions are done in a single launch with no global temporaries:
//
//   out[i,j,k,l] = sum_{M,N,O,P} R1[i,M] R2[j,N] R3[k,O] R4[l,P] block[M,N,O,P]
//
// Templated on the element type and on the (compile-time) per-shell dimensions
// Nk = 2*lk+1, so the inner sums unroll.  Use apply_quartet_fused<T>(...) to
// dispatch on (l1,l2,l3,l4) at runtime.

#pragma once
#include <cuda_runtime.h>
#include <cstdio>
#include <cstddef>

namespace sph_eri {

template <typename T, int N1, int N2, int N3, int N4>
__global__ void apply_fused(const T* __restrict__ blocks,
                            const T* __restrict__ R1g,
                            const T* __restrict__ R2g,
                            const T* __restrict__ R3g,
                            const T* __restrict__ R4g,
                            T* __restrict__ out, int B) {
  constexpr int NOUT = N1 * N2 * N3 * N4;
  const int q = blockIdx.x;
  if (q >= B) return;

  extern __shared__ char smem_raw[];
  T* a  = reinterpret_cast<T*>(smem_raw);
  T* b  = a  + NOUT;
  T* r1 = b  + NOUT;
  T* r2 = r1 + N1 * N1;
  T* r3 = r2 + N2 * N2;
  T* r4 = r3 + N3 * N3;

  for (int e = threadIdx.x; e < NOUT;   e += blockDim.x) a[e]  = blocks[(size_t)q * NOUT + e];
  for (int e = threadIdx.x; e < N1 * N1; e += blockDim.x) r1[e] = R1g[(size_t)q * N1 * N1 + e];
  for (int e = threadIdx.x; e < N2 * N2; e += blockDim.x) r2[e] = R2g[(size_t)q * N2 * N2 + e];
  for (int e = threadIdx.x; e < N3 * N3; e += blockDim.x) r3[e] = R3g[(size_t)q * N3 * N3 + e];
  for (int e = threadIdx.x; e < N4 * N4; e += blockDim.x) r4[e] = R4g[(size_t)q * N4 * N4 + e];
  __syncthreads();

  // step 1: b[i,N,O,P] = sum_M r1[i,M] a[M,N,O,P]
  for (int e = threadIdx.x; e < NOUT; e += blockDim.x) {
    int P = e % N4, t = e / N4, O = t % N3; t /= N3; int Nn = t % N2, i = t / N2;
    T s = 0;
    for (int M = 0; M < N1; ++M) s += r1[i * N1 + M] * a[(((size_t)M * N2 + Nn) * N3 + O) * N4 + P];
    b[e] = s;
  }
  __syncthreads();
  // step 2: a[i,j,O,P] = sum_N r2[j,N] b[i,N,O,P]
  for (int e = threadIdx.x; e < NOUT; e += blockDim.x) {
    int P = e % N4, t = e / N4, O = t % N3; t /= N3; int j = t % N2, i = t / N2;
    T s = 0;
    for (int Nn = 0; Nn < N2; ++Nn) s += r2[j * N2 + Nn] * b[(((size_t)i * N2 + Nn) * N3 + O) * N4 + P];
    a[e] = s;
  }
  __syncthreads();
  // step 3: b[i,j,k,P] = sum_O r3[k,O] a[i,j,O,P]
  for (int e = threadIdx.x; e < NOUT; e += blockDim.x) {
    int P = e % N4, t = e / N4, k = t % N3; t /= N3; int j = t % N2, i = t / N2;
    T s = 0;
    for (int O = 0; O < N3; ++O) s += r3[k * N3 + O] * a[(((size_t)i * N2 + j) * N3 + O) * N4 + P];
    b[e] = s;
  }
  __syncthreads();
  // step 4: out[i,j,k,l] = sum_P r4[l,P] b[i,j,k,P]
  for (int e = threadIdx.x; e < NOUT; e += blockDim.x) {
    int l = e % N4, t = e / N4, k = t % N3; t /= N3; int j = t % N2, i = t / N2;
    T s = 0;
    for (int Pp = 0; Pp < N4; ++Pp) s += r4[l * N4 + Pp] * b[(((size_t)i * N2 + j) * N3 + k) * N4 + Pp];
    out[(size_t)q * NOUT + e] = s;
  }
}

template <typename T, int N1, int N2, int N3, int N4>
inline cudaError_t launch(const T* blocks, const T* R1, const T* R2, const T* R3,
                          const T* R4, T* out, int B, int threads = 128) {
  constexpr int NOUT = N1 * N2 * N3 * N4;
  size_t shbytes = (size_t)(2 * NOUT + N1 * N1 + N2 * N2 + N3 * N3 + N4 * N4) * sizeof(T);
  auto k = apply_fused<T, N1, N2, N3, N4>;
  if (shbytes > 48u * 1024u) {
    cudaFuncSetAttribute(k, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shbytes);
  }
  k<<<B, threads, shbytes>>>(blocks, R1, R2, R3, R4, out, B);
  return cudaGetLastError();
}

// Thread-per-quartet variant: one thread evaluates a whole quartet with the
// direct sum, no shared memory.  For tiny quartets (L=0, NOUT=1) this avoids the
// one-block-per-quartet overhead and is ~40-80x faster (see bench/test_lowl); for
// L>=1 the block-per-quartet `apply_fused` is far better, so this is only routed
// for (00|00) below.
template <typename T, int N1, int N2, int N3, int N4>
__global__ void apply_tpq(const T* __restrict__ blocks, const T* __restrict__ R1g,
                          const T* __restrict__ R2g, const T* __restrict__ R3g,
                          const T* __restrict__ R4g, T* __restrict__ out, int B) {
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= B) return;
  constexpr int NOUT = N1 * N2 * N3 * N4;
  const T* bk = blocks + (size_t)q * NOUT;
  const T* r1 = R1g + (size_t)q * N1 * N1;
  const T* r2 = R2g + (size_t)q * N2 * N2;
  const T* r3 = R3g + (size_t)q * N3 * N3;
  const T* r4 = R4g + (size_t)q * N4 * N4;
  T* o = out + (size_t)q * NOUT;
  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j)
      for (int k = 0; k < N3; ++k)
        for (int l = 0; l < N4; ++l) {
          T s = 0;
          for (int M = 0; M < N1; ++M)
            for (int Nn = 0; Nn < N2; ++Nn)
              for (int O = 0; O < N3; ++O)
                for (int P = 0; P < N4; ++P)
                  s += r1[i * N1 + M] * r2[j * N2 + Nn] * r3[k * N3 + O] *
                       r4[l * N4 + P] * bk[(((size_t)M * N2 + Nn) * N3 + O) * N4 + P];
          o[(((size_t)i * N2 + j) * N3 + k) * N4 + l] = s;
        }
}

template <typename T, int N1, int N2, int N3, int N4>
inline cudaError_t launch_tpq(const T* blocks, const T* R1, const T* R2,
                              const T* R3, const T* R4, T* out, int B,
                              int threads = 256) {
  int blocks_n = (B + threads - 1) / threads;
  apply_tpq<T, N1, N2, N3, N4><<<blocks_n, threads>>>(blocks, R1, R2, R3, R4, out, B);
  return cudaGetLastError();
}

// ---- single-buffer in-place variant (general, runtime dims) ----------------
// apply_fused stages TWO NOUT buffers in shared memory, which caps it near
// (44|44) on a 164 KB A100.  Here each thread owns whole axis-lines, so the
// per-axis transform is race-free *in place* and needs only ONE buffer: that
// doubles the angular momentum that fits in shared memory (reaching (55|55)),
// and the global-resident form (block kept in `out`, only the small R^l in
// shared) removes the shared-memory ceiling entirely (e.g. (66|66) and beyond).
// Slower than the curated compile-time-sized fused kernels, so it is used as the
// dispatch fallback for the high-L / uncurated classes only.
template <typename T>
__device__ inline void transform_axis(T* buf, const T* R, int n, int outer, int inner) {
  int lines = outer * inner;
  for (int lid = threadIdx.x; lid < lines; lid += blockDim.x) {
    int o = lid / inner, b = lid % inner;
    size_t base = (size_t)o * n * inner + b;
    T v[24];  // n = 2l+1; supports l <= 11
    for (int k = 0; k < n; ++k) v[k] = buf[base + (size_t)k * inner];
    for (int i = 0; i < n; ++i) {
      T s = 0;
      for (int m = 0; m < n; ++m) s += R[i * n + m] * v[m];
      buf[base + (size_t)i * inner] = s;
    }
  }
}

template <typename T>
__global__ void apply_inplace_smem(const T* __restrict__ blocks, const T* __restrict__ R1g,
                                   const T* __restrict__ R2g, const T* __restrict__ R3g,
                                   const T* __restrict__ R4g, T* __restrict__ out, int B,
                                   int N1, int N2, int N3, int N4) {
  int q = blockIdx.x; if (q >= B) return;
  int NOUT = N1 * N2 * N3 * N4;
  extern __shared__ char sm[];
  T* buf = reinterpret_cast<T*>(sm);
  T* r1 = buf + NOUT; T* r2 = r1 + N1*N1; T* r3 = r2 + N2*N2; T* r4 = r3 + N3*N3;
  for (int e = threadIdx.x; e < NOUT;  e += blockDim.x) buf[e] = blocks[(size_t)q*NOUT + e];
  for (int e = threadIdx.x; e < N1*N1; e += blockDim.x) r1[e] = R1g[(size_t)q*N1*N1 + e];
  for (int e = threadIdx.x; e < N2*N2; e += blockDim.x) r2[e] = R2g[(size_t)q*N2*N2 + e];
  for (int e = threadIdx.x; e < N3*N3; e += blockDim.x) r3[e] = R3g[(size_t)q*N3*N3 + e];
  for (int e = threadIdx.x; e < N4*N4; e += blockDim.x) r4[e] = R4g[(size_t)q*N4*N4 + e];
  __syncthreads();
  transform_axis<T>(buf, r1, N1, 1,        N2*N3*N4); __syncthreads();
  transform_axis<T>(buf, r2, N2, N1,       N3*N4);    __syncthreads();
  transform_axis<T>(buf, r3, N3, N1*N2,    N4);       __syncthreads();
  transform_axis<T>(buf, r4, N4, N1*N2*N3, 1);        __syncthreads();
  for (int e = threadIdx.x; e < NOUT; e += blockDim.x) out[(size_t)q*NOUT + e] = buf[e];
}

template <typename T>
__global__ void apply_inplace_gmem(const T* __restrict__ blocks, const T* __restrict__ R1g,
                                   const T* __restrict__ R2g, const T* __restrict__ R3g,
                                   const T* __restrict__ R4g, T* __restrict__ out, int B,
                                   int N1, int N2, int N3, int N4) {
  int q = blockIdx.x; if (q >= B) return;
  int NOUT = N1 * N2 * N3 * N4;
  T* buf = out + (size_t)q * NOUT;  // block resident in global; `out` doubles as scratch
  extern __shared__ char sm[];
  T* r1 = reinterpret_cast<T*>(sm);
  T* r2 = r1 + N1*N1; T* r3 = r2 + N2*N2; T* r4 = r3 + N3*N3;
  for (int e = threadIdx.x; e < NOUT;  e += blockDim.x) buf[e] = blocks[(size_t)q*NOUT + e];
  for (int e = threadIdx.x; e < N1*N1; e += blockDim.x) r1[e] = R1g[(size_t)q*N1*N1 + e];
  for (int e = threadIdx.x; e < N2*N2; e += blockDim.x) r2[e] = R2g[(size_t)q*N2*N2 + e];
  for (int e = threadIdx.x; e < N3*N3; e += blockDim.x) r3[e] = R3g[(size_t)q*N3*N3 + e];
  for (int e = threadIdx.x; e < N4*N4; e += blockDim.x) r4[e] = R4g[(size_t)q*N4*N4 + e];
  __syncthreads();
  transform_axis<T>(buf, r1, N1, 1,        N2*N3*N4); __syncthreads();
  transform_axis<T>(buf, r2, N2, N1,       N3*N4);    __syncthreads();
  transform_axis<T>(buf, r3, N3, N1*N2,    N4);       __syncthreads();
  transform_axis<T>(buf, r4, N4, N1*N2*N3, 1);        __syncthreads();
}

template <typename T>
inline cudaError_t launch_inplace(const T* blocks, const T* R1, const T* R2, const T* R3,
                                  const T* R4, T* out, int N1, int N2, int N3, int N4,
                                  int B, int threads = 128) {
  size_t shR = (size_t)(N1*N1 + N2*N2 + N3*N3 + N4*N4) * sizeof(T);
  size_t shFull = shR + (size_t)(N1*N2*N3*N4) * sizeof(T);
  if (shFull <= 160u * 1024u) {           // whole block fits in one shared buffer
    if (shFull > 48u * 1024u)
      cudaFuncSetAttribute(apply_inplace_smem<T>,
                           cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shFull);
    apply_inplace_smem<T><<<B, threads, shFull>>>(blocks, R1, R2, R3, R4, out, B, N1, N2, N3, N4);
  } else {                                 // block too big: keep it global-resident
    apply_inplace_gmem<T><<<B, threads, shR>>>(blocks, R1, R2, R3, R4, out, B, N1, N2, N3, N4);
  }
  return cudaGetLastError();
}

// Runtime dispatch on (l1,l2,l3,l4) for a curated set of quartet types.
template <typename T>
inline bool apply_quartet_fused(int l1, int l2, int l3, int l4,
                                const T* blocks, const T* R1, const T* R2,
                                const T* R3, const T* R4, T* out, int B,
                                int threads = 0) {
  // Adaptive block size by output size NOUT (one block per quartet): small NOUT
  // wants few warps so more blocks/SM run concurrently; large NOUT wants more
  // warps/block to recover occupancy under its big working set.  Measured optima
  // on A100 (FP64): NOUT<=200 ->64, <=1000 ->128, <=4000 ->256, else ->768
  // (~2-2.4x at (44|44)-(66|66); low-L unregressed, (11|11) slightly improved).
  // Pass threads>0 to override.
  if (threads <= 0) {
    long nout = (long)(2*l1+1) * (2*l2+1) * (2*l3+1) * (2*l4+1);
    threads = nout <= 200 ? 64 : nout <= 1000 ? 128 : nout <= 4000 ? 256 : 768;
  }
  // (00|00): NOUT=1, so map one thread per quartet rather than one block.
  if (l1 == 0 && l2 == 0 && l3 == 0 && l4 == 0) {
    launch_tpq<T, 1, 1, 1, 1>(blocks, R1, R2, R3, R4, out, B);
    return true;
  }
#define SPH_CASE(a, b, c, d)                                                    \
  if (l1 == (a) && l2 == (b) && l3 == (c) && l4 == (d)) {                       \
    launch<T, 2 * (a) + 1, 2 * (b) + 1, 2 * (c) + 1, 2 * (d) + 1>(              \
        blocks, R1, R2, R3, R4, out, B, threads);                              \
    return true;                                                               \
  }
  SPH_CASE(1, 1, 1, 1)
  SPH_CASE(2, 2, 1, 1)
  SPH_CASE(2, 2, 2, 2)
  SPH_CASE(3, 3, 2, 2)
  SPH_CASE(4, 4, 3, 3)
  SPH_CASE(5, 5, 4, 4)
  // diagonal classes of the LibintX eri4 benchmark that still fit the two-buffer
  // kernel (up to (44|44) ~ 105 KB on a 164 KB A100):
  SPH_CASE(2, 0, 2, 0)
  SPH_CASE(2, 2, 2, 0)
  SPH_CASE(3, 3, 3, 3)
  SPH_CASE(4, 4, 4, 4)
#undef SPH_CASE
  // Everything else (incl. (55|55), (66|66), and any uncurated class) falls back
  // to the general single-buffer in-place kernel, which has no compile-time table
  // and no shared-memory ceiling (global-resident above ~160 KB).
  return launch_inplace<T>(blocks, R1, R2, R3, R4, out,
                           2 * l1 + 1, 2 * l2 + 1, 2 * l3 + 1, 2 * l4 + 1, B,
                           threads) == cudaSuccess;
}

}  // namespace sph_eri
