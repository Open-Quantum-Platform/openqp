#include <climits>
#include <cuda_runtime.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

namespace {

static long g_metc_kernel_launches = 0;
static double g_metc_kernel_ms = 0.0;
static double g_metc_sync_s = 0.0;

static const int OQP_METC_PROFILE_KERNELS = 5;
static const int OQP_METC_PROFILE_MAX_SAMPLES = 8192;
static const int OQP_METC_BLOCK_ACCUM_SLOTS = 1024;

struct KernelProfile {
  const char* name;
  long count;
  double total_ms;
  double min_ms;
  double max_ms;
  long total_blocks;
  long total_threads_per_block;
  long max_blocks;
  long max_threads_per_block;
  long long total_elements;
  long long min_elements;
  long long max_elements;
  long long total_ncur;
  long long min_ncur;
  long long max_ncur;
  long long total_est_flops;
  long long total_est_bytes;
  long long total_est_atomics;
  double samples[OQP_METC_PROFILE_MAX_SAMPLES];
};

static KernelProfile g_kernel_profiles[OQP_METC_PROFILE_KERNELS] = {
    {"mrsf_metc_kernel", 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0.0}},
    {"umrsf_metc_kernel", 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0.0}},
    {"umrsf_metc_combined_coulomb_kernel", 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0.0}},
    {"umrsf_metc_combined_coulomb_warp_kernel", 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0.0}},
    {"umrsf_metc_two_phase_accum_kernel", 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0.0}},
};

static double host_now_sec() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return static_cast<double>(ts.tv_sec) + static_cast<double>(ts.tv_nsec) * 1.0e-9;
}

static bool profile_enabled() {
  const char* v = getenv("OQP_GPU_PROFILE");
  if (!v) v = getenv("OQP_GPU_METC_PROFILE");
  return v && (v[0] == '1' || v[0] == 't' || v[0] == 'T' ||
               v[0] == 'y' || v[0] == 'Y' || v[0] == 'o' || v[0] == 'O');
}

static bool env_flag_enabled(const char* name) {
  const char* v = getenv(name);
  return v && (v[0] == '1' || v[0] == 't' || v[0] == 'T' ||
               v[0] == 'y' || v[0] == 'Y' || v[0] == 'o' || v[0] == 'O');
}

static int umrsf_variant_id() {
  const char* v = getenv("OQP_GPU_METC_VARIANT");
  if (!v || !*v) v = getenv("OQP_GPU_ACCUM_VARIANT");
  if (v && (strcmp(v, "combined_coulomb_warp") == 0 ||
            strcmp(v, "warp_agg") == 0 ||
            strcmp(v, "accum_v2") == 0)) {
    return 2;
  }
  if (v && (strcmp(v, "two_phase_accum") == 0 ||
            strcmp(v, "block_private") == 0 ||
            strcmp(v, "accum_v3") == 0)) {
    return 3;
  }
  if ((v && (strcmp(v, "combined_coulomb") == 0 || strcmp(v, "accum_v1") == 0)) ||
      env_flag_enabled("OQP_GPU_METC_COMBINE_COULOMB")) {
    return 1;
  }
  return 0;
}

static int cmp_double(const void* a, const void* b) {
  const double da = *static_cast<const double*>(a);
  const double db = *static_cast<const double*>(b);
  return (da > db) - (da < db);
}

static long long active_terms_per_nf_m(int is_umrsf, int cur_pass) {
  if (is_umrsf) {
    if (cur_pass == 1) return 152;  // m 0-7: 16 terms; m 8-9: 8 terms; m 10: 8 terms
    if (cur_pass == 2) return 8;    // m 10 exchange-only continuation
    return 0;
  }
  if (cur_pass == 1) return 88;     // m 0-3: 16 terms; m 4-6: 8 terms
  if (cur_pass == 2) return 8;      // m 6 exchange-only continuation
  return 0;
}

static void estimate_work(int profile_id, int ncur, int nf, int nmatrix, int cur_pass,
                          long long* est_flops, long long* est_bytes,
                          long long* est_atomics) {
  const int is_umrsf = (profile_id == 1 || profile_id == 2 || profile_id == 3 || profile_id == 4);
  const long long total_elements = static_cast<long long>(ncur) * nf * nmatrix;
  const long long d3_reads_per_nf = active_terms_per_nf_m(is_umrsf, cur_pass);
  long long atomic_updates_per_nf = d3_reads_per_nf;
  if ((profile_id == 2 || profile_id == 3 || profile_id == 4) && cur_pass == 1) {
    // UMRSF combined_coulomb variant: m<8 Coulomb terms keep the same 8 d3
    // reads but combine duplicate f3 targets, so those 8 atomics become 4.
    // Exchange terms and m=8/9/10 remain unchanged: 8*12 + 2*8 + 1*8 = 120.
    atomic_updates_per_nf = 120;
  }
  const long long d3_reads = static_cast<long long>(ncur) * nf * d3_reads_per_nf;
  const long long atomics = static_cast<long long>(ncur) * nf * atomic_updates_per_nf;
  // Approximate only: each active term does one multiply and one atomic add in FP64.
  // Atomic update is counted as read+write of f3 plus one d3 read. ids/int metadata
  // are counted once per launched logical element.
  *est_flops = 2LL * d3_reads + 2LL * total_elements;
  *est_bytes = d3_reads * 8LL + atomics * 16LL + total_elements * 12LL;
  *est_atomics = atomics;
}

static void record_kernel_profile(int id, float ms, int blocks, int threads,
                                  int ncur, int nf, int nmatrix, int cur_pass) {
  if (id < 0 || id >= OQP_METC_PROFILE_KERNELS) return;
  KernelProfile& p = g_kernel_profiles[id];
  const double dms = static_cast<double>(ms);
  const long long elements = static_cast<long long>(ncur) * nf * nmatrix;
  long long est_flops = 0;
  long long est_bytes = 0;
  long long est_atomics = 0;
  estimate_work(id, ncur, nf, nmatrix, cur_pass, &est_flops, &est_bytes, &est_atomics);
  if (p.count == 0) {
    p.min_ms = p.max_ms = dms;
    p.min_elements = p.max_elements = elements;
    p.min_ncur = p.max_ncur = ncur;
  } else {
    if (dms < p.min_ms) p.min_ms = dms;
    if (dms > p.max_ms) p.max_ms = dms;
    if (elements < p.min_elements) p.min_elements = elements;
    if (elements > p.max_elements) p.max_elements = elements;
    if (ncur < p.min_ncur) p.min_ncur = ncur;
    if (ncur > p.max_ncur) p.max_ncur = ncur;
  }
  if (blocks > p.max_blocks) p.max_blocks = blocks;
  if (threads > p.max_threads_per_block) p.max_threads_per_block = threads;
  p.count += 1;
  p.total_ms += dms;
  p.total_blocks += blocks;
  p.total_threads_per_block += threads;
  p.total_elements += elements;
  p.total_ncur += ncur;
  p.total_est_flops += est_flops;
  p.total_est_bytes += est_bytes;
  p.total_est_atomics += est_atomics;
  if (p.count <= OQP_METC_PROFILE_MAX_SAMPLES) {
    p.samples[p.count - 1] = dms;
  }
}

__device__ __forceinline__ int idx4(int f, int m, int row, int col, int nf, int nmatrix, int nbf) {
  return f + nf * (m + nmatrix * (row + nbf * col));
}

__device__ __forceinline__ double atomic_add_double(double* address, double value) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
  unsigned long long int* address_as_ull = reinterpret_cast<unsigned long long int*>(address);
  unsigned long long int old = *address_as_ull;
  unsigned long long int assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(value + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
#else
  return atomicAdd(address, value);
#endif
}

__device__ __forceinline__ void add4(double* f3, int f, int m, int row, int col,
                                     int nf, int nmatrix, int nbf, double value) {
  atomic_add_double(&f3[idx4(f, m, row, col, nf, nmatrix, nbf)], value);
}

__device__ __forceinline__ void atomic_add_double_warp_agg(double* address, double value) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700
  const unsigned active = __activemask();
  const unsigned long long key = reinterpret_cast<unsigned long long>(address);
  const unsigned same_key = __match_any_sync(active, key);
  const int lane = threadIdx.x & 31;
  const int leader = __ffs(same_key) - 1;
  if (lane == leader) {
    double sum = 0.0;
    for (int src = 0; src < 32; ++src) {
      if (same_key & (1u << src)) {
        sum += __shfl_sync(active, value, src);
      }
    }
    atomic_add_double(address, sum);
  }
#else
  atomic_add_double(address, value);
#endif
}

__device__ __forceinline__ void add4_warp(double* f3, int f, int m, int row, int col,
                                          int nf, int nmatrix, int nbf, double value) {
  atomic_add_double_warp_agg(&f3[idx4(f, m, row, col, nf, nmatrix, nbf)], value);
}

__device__ __forceinline__ unsigned int hash_key_int(int key) {
  unsigned int x = static_cast<unsigned int>(key);
  x ^= x >> 16;
  x *= 0x7feb352dU;
  x ^= x >> 15;
  x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

__device__ __forceinline__ void add4_block_private(int* keys, double* vals,
                                                   double* f3, int f, int m, int row, int col,
                                                   int nf, int nmatrix, int nbf, double value) {
  const int key = idx4(f, m, row, col, nf, nmatrix, nbf);
  unsigned int slot = hash_key_int(key) & (OQP_METC_BLOCK_ACCUM_SLOTS - 1);
  for (int probe = 0; probe < 16; ++probe) {
    int old = atomicCAS(&keys[slot], -1, key);
    if (old == -1 || old == key) {
      atomic_add_double(&vals[slot], value);
      return;
    }
    slot = (slot + 1) & (OQP_METC_BLOCK_ACCUM_SLOTS - 1);
  }
  atomic_add_double(&f3[key], value);
}

__device__ __forceinline__ double get4(const double* d3, int f, int m, int row, int col,
                                       int nf, int nmatrix, int nbf) {
  return d3[idx4(f, m, row, col, nf, nmatrix, nbf)];
}

__global__ void mrsf_metc_kernel(const int* ids, const double* ints, int ncur,
                                 double* f3, const double* d3, int nf, int nmatrix,
                                 int nbf, int cur_pass, double scale_exchange,
                                 double scale_coulomb) {
  int linear = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncur * nf * nmatrix;
  if (linear >= total) return;

  int m = linear % nmatrix;
  int tmp = linear / nmatrix;
  int f = tmp % nf;
  int n = tmp / nf;

  int i = ids[4 * n + 0] - 1;
  int j = ids[4 * n + 1] - 1;
  int k = ids[4 * n + 2] - 1;
  int l = ids[4 * n + 3] - 1;
  double val = ints[n];
  double xval = val * scale_exchange;
  double cval = val * scale_coulomb;

  if (cur_pass == 1) {
    if (m < 4) {
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
    }
    if (m < 7) {
      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
  } else if (cur_pass == 2 && m == 6) {
    add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
    add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
    add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
    add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
  }
}

__global__ void umrsf_metc_kernel(const int* ids, const double* ints, int ncur,
                                  double* f3, const double* d3, int nf, int nmatrix,
                                  int nbf, int cur_pass, double scale_exchange,
                                  double scale_coulomb) {
  int linear = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncur * nf * nmatrix;
  if (linear >= total) return;

  int m = linear % nmatrix;
  int tmp = linear / nmatrix;
  int f = tmp % nf;
  int n = tmp / nf;

  int i = ids[4 * n + 0] - 1;
  int j = ids[4 * n + 1] - 1;
  int k = ids[4 * n + 2] - 1;
  int l = ids[4 * n + 3] - 1;
  double val = ints[n];
  double xval = val * scale_exchange;
  double cval = val * scale_coulomb;

  if (cur_pass == 1) {
    if (m < 8) {
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));

      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
    if (m == 8 || m == 9) {
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
  }

  if ((cur_pass == 1 || cur_pass == 2) && m == 10) {
    add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
    add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
    add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
    add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
  }
}


__global__ void umrsf_metc_combined_coulomb_kernel(const int* ids, const double* ints, int ncur,
                                  double* f3, const double* d3, int nf, int nmatrix,
                                  int nbf, int cur_pass, double scale_exchange,
                                  double scale_coulomb) {
  int linear = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncur * nf * nmatrix;
  if (linear >= total) return;

  int m = linear % nmatrix;
  int tmp = linear / nmatrix;
  int f = tmp % nf;
  int n = tmp / nf;

  int i = ids[4 * n + 0] - 1;
  int j = ids[4 * n + 1] - 1;
  int k = ids[4 * n + 2] - 1;
  int l = ids[4 * n + 3] - 1;
  double val = ints[n];
  double xval = val * scale_exchange;
  double cval = val * scale_coulomb;

  if (cur_pass == 1) {
    if (m < 8) {
      // Algorithm-preserving accumulation variant: combine duplicate Coulomb
      // writes within one (n,f,m) thread before one atomic add per f3 target.
      // This keeps all inputs and target indices identical to the reference
      // implementation, but cuts the Coulomb atomic updates for m<8 from 8 to 4.
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * (get4(d3, f, m, k, l, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, l, k, nf, nmatrix, nbf)));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * (get4(d3, f, m, i, j, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, j, i, nf, nmatrix, nbf)));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * (get4(d3, f, m, k, l, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, l, k, nf, nmatrix, nbf)));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * (get4(d3, f, m, i, j, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, j, i, nf, nmatrix, nbf)));

      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
    if (m == 8 || m == 9) {
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
  }

  if ((cur_pass == 1 || cur_pass == 2) && m == 10) {
    add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
    add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
    add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
    add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
  }
}

__global__ void umrsf_metc_combined_coulomb_warp_kernel(const int* ids, const double* ints, int ncur,
                                  double* f3, const double* d3, int nf, int nmatrix,
                                  int nbf, int cur_pass, double scale_exchange,
                                  double scale_coulomb) {
  int linear = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncur * nf * nmatrix;
  if (linear >= total) return;

  int m = linear % nmatrix;
  int tmp = linear / nmatrix;
  int f = tmp % nf;
  int n = tmp / nf;

  int i = ids[4 * n + 0] - 1;
  int j = ids[4 * n + 1] - 1;
  int k = ids[4 * n + 2] - 1;
  int l = ids[4 * n + 3] - 1;
  double val = ints[n];
  double xval = val * scale_exchange;
  double cval = val * scale_coulomb;

  if (cur_pass == 1) {
    if (m < 8) {
      // Algorithm-preserving warp aggregation variant: combine duplicate Coulomb
      // writes within one (n,f,m) thread before one atomic add per f3 target.
      // This keeps all inputs and target indices identical to the combined-coulomb
      // variant, then lets lanes in a warp with the same f3 address aggregate before atomicAdd.
      add4_warp(f3, f, m, i, j, nf, nmatrix, nbf, cval * (get4(d3, f, m, k, l, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, l, k, nf, nmatrix, nbf)));
      add4_warp(f3, f, m, k, l, nf, nmatrix, nbf, cval * (get4(d3, f, m, i, j, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, j, i, nf, nmatrix, nbf)));
      add4_warp(f3, f, m, j, i, nf, nmatrix, nbf, cval * (get4(d3, f, m, k, l, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, l, k, nf, nmatrix, nbf)));
      add4_warp(f3, f, m, l, k, nf, nmatrix, nbf, cval * (get4(d3, f, m, i, j, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, j, i, nf, nmatrix, nbf)));

      add4_warp(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4_warp(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4_warp(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4_warp(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4_warp(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4_warp(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4_warp(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4_warp(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
    if (m == 8 || m == 9) {
      add4_warp(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4_warp(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4_warp(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4_warp(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4_warp(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4_warp(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4_warp(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4_warp(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
  }

  if ((cur_pass == 1 || cur_pass == 2) && m == 10) {
    add4_warp(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
    add4_warp(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
    add4_warp(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
    add4_warp(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
    add4_warp(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
    add4_warp(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
    add4_warp(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
    add4_warp(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
  }
}

__global__ void umrsf_metc_two_phase_accum_kernel(const int* ids, const double* ints, int ncur,
                                  double* f3, const double* d3, int nf, int nmatrix,
                                  int nbf, int cur_pass, double scale_exchange,
                                  double scale_coulomb) {
  int linear = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncur * nf * nmatrix;
  extern __shared__ unsigned char smem[];
  int* skeys = reinterpret_cast<int*>(smem);
  double* svals = reinterpret_cast<double*>(&skeys[OQP_METC_BLOCK_ACCUM_SLOTS]);
  for (int slot = threadIdx.x; slot < OQP_METC_BLOCK_ACCUM_SLOTS; slot += blockDim.x) {
    skeys[slot] = -1;
    svals[slot] = 0.0;
  }
  __syncthreads();

  if (linear < total) {

  int m = linear % nmatrix;
  int tmp = linear / nmatrix;
  int f = tmp % nf;
  int n = tmp / nf;

  int i = ids[4 * n + 0] - 1;
  int j = ids[4 * n + 1] - 1;
  int k = ids[4 * n + 2] - 1;
  int l = ids[4 * n + 3] - 1;
  double val = ints[n];
  double xval = val * scale_exchange;
  double cval = val * scale_coulomb;

  if (cur_pass == 1) {
    if (m < 8) {
      // Algorithm-preserving accumulation variant: combine duplicate Coulomb
      // writes within one (n,f,m) thread before one atomic add per f3 target.
      // This keeps all inputs and target indices identical to the reference
      // implementation, but cuts the Coulomb atomic updates for m<8 from 8 to 4.
      add4_block_private(skeys, svals, f3, f, m, i, j, nf, nmatrix, nbf, cval * (get4(d3, f, m, k, l, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, l, k, nf, nmatrix, nbf)));
      add4_block_private(skeys, svals, f3, f, m, k, l, nf, nmatrix, nbf, cval * (get4(d3, f, m, i, j, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, j, i, nf, nmatrix, nbf)));
      add4_block_private(skeys, svals, f3, f, m, j, i, nf, nmatrix, nbf, cval * (get4(d3, f, m, k, l, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, l, k, nf, nmatrix, nbf)));
      add4_block_private(skeys, svals, f3, f, m, l, k, nf, nmatrix, nbf, cval * (get4(d3, f, m, i, j, nf, nmatrix, nbf) +
                                                      get4(d3, f, m, j, i, nf, nmatrix, nbf)));

      add4_block_private(skeys, svals, f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
    if (m == 8 || m == 9) {
      add4_block_private(skeys, svals, f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4_block_private(skeys, svals, f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
  }

  if ((cur_pass == 1 || cur_pass == 2) && m == 10) {
    add4_block_private(skeys, svals, f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
    add4_block_private(skeys, svals, f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
    add4_block_private(skeys, svals, f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
    add4_block_private(skeys, svals, f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
    add4_block_private(skeys, svals, f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
    add4_block_private(skeys, svals, f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
    add4_block_private(skeys, svals, f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
    add4_block_private(skeys, svals, f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
  }

  }

  __syncthreads();
  for (int slot = threadIdx.x; slot < OQP_METC_BLOCK_ACCUM_SLOTS; slot += blockDim.x) {
    const int key = skeys[slot];
    const double value = svals[slot];
    if (key >= 0 && value != 0.0) {
      atomic_add_double(&f3[key], value);
    }
  }
}

}  // namespace

extern "C" void oqp_gpu_metc_cuda_profile_reset(void) {
  g_metc_kernel_launches = 0;
  g_metc_kernel_ms = 0.0;
  g_metc_sync_s = 0.0;
  for (int i = 0; i < OQP_METC_PROFILE_KERNELS; ++i) {
    const char* name = g_kernel_profiles[i].name;
    memset(&g_kernel_profiles[i], 0, sizeof(KernelProfile));
    g_kernel_profiles[i].name = name;
  }
}

extern "C" void oqp_gpu_metc_cuda_profile_snapshot(long *launches, double *kernel_ms,
                                                    double *sync_s) {
  if (launches) *launches = g_metc_kernel_launches;
  if (kernel_ms) *kernel_ms = g_metc_kernel_ms;
  if (sync_s) *sync_s = g_metc_sync_s;
}

extern "C" void oqp_gpu_metc_cuda_profile_dump(void) {
  for (int i = 0; i < OQP_METC_PROFILE_KERNELS; ++i) {
    KernelProfile& p = g_kernel_profiles[i];
    if (p.count <= 0) continue;
    const long nsample = (p.count < OQP_METC_PROFILE_MAX_SAMPLES) ? p.count : OQP_METC_PROFILE_MAX_SAMPLES;
    double sorted[OQP_METC_PROFILE_MAX_SAMPLES];
    for (long j = 0; j < nsample; ++j) sorted[j] = p.samples[j];
    qsort(sorted, static_cast<size_t>(nsample), sizeof(double), cmp_double);
    const long median_idx = nsample / 2;
    long p95_idx = static_cast<long>(0.95 * static_cast<double>(nsample - 1));
    if (p95_idx < 0) p95_idx = 0;
    if (p95_idx >= nsample) p95_idx = nsample - 1;
    const double median_ms = sorted[median_idx];
    const double p95_ms = sorted[p95_idx];
    const double mean_ms = p.total_ms / static_cast<double>(p.count);
    const double mean_blocks = static_cast<double>(p.total_blocks) / static_cast<double>(p.count);
    const double mean_elements = static_cast<double>(p.total_elements) / static_cast<double>(p.count);
    const double intensity = (p.total_est_bytes > 0)
        ? static_cast<double>(p.total_est_flops) / static_cast<double>(p.total_est_bytes)
        : 0.0;
    fprintf(stderr,
            "OQP_GPU_KERNEL_PROFILE schema=1 kernel=%s call_count=%ld total_s=%.9f"
            " mean_s=%.9f median_s=%.9f p95_s=%.9f min_s=%.9f max_s=%.9f"
            " block_dim=%ld grid_blocks_mean=%.3f grid_blocks_max=%ld ncur_min=%lld ncur_max=%lld"
            " elements_total=%lld elements_mean=%.3f elements_min=%lld elements_max=%lld"
            " est_flops=%lld est_bytes=%lld est_atomics=%lld est_arith_intensity=%.9f sample_count=%ld"
            " timing_note=kernel_event_and_sync_overlap\n",
            p.name, p.count, p.total_ms * 1.0e-3, mean_ms * 1.0e-3,
            median_ms * 1.0e-3, p95_ms * 1.0e-3, p.min_ms * 1.0e-3,
            p.max_ms * 1.0e-3, p.max_threads_per_block, mean_blocks,
            p.max_blocks, p.min_ncur, p.max_ncur, p.total_elements,
            mean_elements, p.min_elements, p.max_elements, p.total_est_flops,
            p.total_est_bytes, p.total_est_atomics, intensity, nsample);
  }
}

// Borrowed-pointer kernel launcher (METC-C1b resident hot path).
extern "C" int oqp_gpu_metc_launch(const int* d_ids, const double* d_ints, int ncur,
                                   double* d_f3, const double* d_d3, int nf,
                                   int nmatrix, int nbf, int cur_pass,
                                   double scale_exchange, double scale_coulomb,
                                   int is_umrsf) {
  if (ncur < 0 || nf <= 0 || nbf <= 0 ||
      nmatrix != (is_umrsf ? 11 : 7) || (cur_pass != 1 && cur_pass != 2))
    return (int)cudaErrorInvalidValue;
  if (nf > INT_MAX / nmatrix || nbf > INT_MAX / (nf * nmatrix) ||
      nbf > INT_MAX / (nf * nmatrix * nbf) ||
      ncur > (INT_MAX - 255) / (nf * nmatrix))
    return (int)cudaErrorInvalidValue;
  if (ncur == 0) return 0;

  int threads = 256;
  int total = ncur * nf * nmatrix;
  int blocks = (total + threads - 1) / threads;
  const bool prof = profile_enabled();
  cudaEvent_t start, stop;
  if (prof) {
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
  }
  const int variant = is_umrsf ? umrsf_variant_id() : 0;
  if (variant == 3) {
    const size_t shared_bytes = static_cast<size_t>(OQP_METC_BLOCK_ACCUM_SLOTS) *
        (sizeof(int) + sizeof(double));
    umrsf_metc_two_phase_accum_kernel<<<blocks, threads, shared_bytes>>>(
        d_ids, d_ints, ncur, d_f3, d_d3, nf, nmatrix, nbf, cur_pass,
        scale_exchange, scale_coulomb);
  } else if (variant == 2) {
    umrsf_metc_combined_coulomb_warp_kernel<<<blocks, threads>>>(
        d_ids, d_ints, ncur, d_f3, d_d3, nf, nmatrix, nbf, cur_pass,
        scale_exchange, scale_coulomb);
  } else if (variant == 1) {
    umrsf_metc_combined_coulomb_kernel<<<blocks, threads>>>(
        d_ids, d_ints, ncur, d_f3, d_d3, nf, nmatrix, nbf, cur_pass,
        scale_exchange, scale_coulomb);
  } else if (is_umrsf) {
    umrsf_metc_kernel<<<blocks, threads>>>(d_ids, d_ints, ncur, d_f3, d_d3, nf,
                                           nmatrix, nbf, cur_pass, scale_exchange,
                                           scale_coulomb);
  } else {
    mrsf_metc_kernel<<<blocks, threads>>>(d_ids, d_ints, ncur, d_f3, d_d3, nf,
                                          nmatrix, nbf, cur_pass, scale_exchange,
                                          scale_coulomb);
  }

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    if (prof) { cudaEventDestroy(start); cudaEventDestroy(stop); }
    return static_cast<int>(err);
  }
  double t_sync = host_now_sec();
  if (prof) cudaEventRecord(stop, 0);
  err = cudaDeviceSynchronize();
  if (prof) {
    g_metc_sync_s += host_now_sec() - t_sync;
    cudaEventSynchronize(stop);
    float ms = 0.0f;
    cudaEventElapsedTime(&ms, start, stop);
    g_metc_kernel_ms += static_cast<double>(ms);
    g_metc_kernel_launches += 1;
    const int profile_id = (variant == 3) ? 4 : ((variant == 2) ? 3 : ((variant == 1) ? 2 : (is_umrsf ? 1 : 0)));
    record_kernel_profile(profile_id, ms, blocks, threads, ncur, nf, nmatrix, cur_pass);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  }
  return static_cast<int>(err);
}

// Legacy owning wrapper (METC-B).  Preserved only as an optional fallback / comparison path.
extern "C" int oqp_gpu_metc_contract_owning(const int* ids, const double* ints, int ncur,
                                      double* f3, const double* d3, int nf,
                                      int nmatrix, int nbf, int cur_pass,
                                      double scale_exchange, double scale_coulomb,
                                      bool is_umrsf) {
  if (ncur < 0 || nf <= 0 || nbf <= 0 ||
      nmatrix != (is_umrsf ? 11 : 7) || (cur_pass != 1 && cur_pass != 2))
    return (int)cudaErrorInvalidValue;
  if (nf > INT_MAX / nmatrix || nbf > INT_MAX / (nf * nmatrix) ||
      nbf > INT_MAX / (nf * nmatrix * nbf) ||
      ncur > (INT_MAX - 255) / (nf * nmatrix))
    return (int)cudaErrorInvalidValue;
  if (ncur == 0) return 0;

  if (!ids || !ints || !f3 || !d3) return (int)cudaErrorInvalidValue;
  for (size_t i = 0; i < static_cast<size_t>(ncur) * 4; ++i)
    if (ids[i] < 1 || ids[i] > nbf) return (int)cudaErrorInvalidValue;
  size_t ids_bytes = static_cast<size_t>(4) * ncur * sizeof(int);
  size_t ints_bytes = static_cast<size_t>(ncur) * sizeof(double);
  size_t tensor_count = static_cast<size_t>(nf) * nmatrix * nbf * nbf;
  size_t tensor_bytes = tensor_count * sizeof(double);

  int* d_ids = nullptr;
  double* d_ints = nullptr;
  double* d_f3 = nullptr;
  double* d_d3 = nullptr;

  cudaError_t err = cudaMalloc(&d_ids, ids_bytes);
  if (err != cudaSuccess) return static_cast<int>(err);
  err = cudaMalloc(&d_ints, ints_bytes);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMalloc(&d_f3, tensor_bytes);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMalloc(&d_d3, tensor_bytes);
  if (err != cudaSuccess) goto cleanup;

  err = cudaMemcpy(d_ids, ids, ids_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMemcpy(d_ints, ints, ints_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMemcpy(d_f3, f3, tensor_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMemcpy(d_d3, d3, tensor_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;

  // Preserve the caller's initial accumulator, as in the CPU update.

  err = static_cast<cudaError_t>(oqp_gpu_metc_launch(d_ids, d_ints, ncur, d_f3, d_d3, nf,
                                                     nmatrix, nbf, cur_pass, scale_exchange,
                                                     scale_coulomb, is_umrsf ? 1 : 0));
  if (err != cudaSuccess) goto cleanup;

  err = cudaMemcpy(f3, d_f3, tensor_bytes, cudaMemcpyDeviceToHost);

cleanup:
  if (d_d3) cudaFree(d_d3);
  if (d_f3) cudaFree(d_f3);
  if (d_ints) cudaFree(d_ints);
  if (d_ids) cudaFree(d_ids);
  return static_cast<int>(err);
}

// Compatibility entry point for callers of the original METC kernel.
extern "C" int oqp_gpu_metc_contract(const int* ids, const double* ints, int ncur,
    double* f3, const double* d3, int nf, int nmatrix, int nbf, int cur_pass,
    double scale_exchange, double scale_coulomb, bool is_umrsf) {
  return oqp_gpu_metc_contract_owning(ids, ints, ncur, f3, d3, nf, nmatrix,
      nbf, cur_pass, scale_exchange, scale_coulomb, is_umrsf);
}
extern "C" int oqp_gpu_metc_device_count(void) {
  int count = 0;
  return cudaGetDeviceCount(&count) == cudaSuccess ? count : 0;
}
