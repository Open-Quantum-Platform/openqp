/* Minimal real workspace runtime for the unified GPU workspace bridge (METC-C1a).
 *
 * This backs the C ABI declared in source/gpu_workspace_bridge.F90 with a real
 * contiguous-arena allocator plus pointer lookup by byte offset.  Under
 * OQP_CUDA_ENABLE the arena is a CUDA device allocation; otherwise it is a host
 * allocation so non-CUDA builds and source-level tests stay meaningful.
 *
 * METC-C1a boundary: this file owns ONLY allocation, pointer arithmetic, and
 * release.  It contains no CUDA kernels and no METC contraction logic.  The
 * METC contraction wrapper is NOT rewired to these pointers yet (that is C1b).
 */

/* Expose clock_gettime / CLOCK_MONOTONIC on strict-libc builds. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <limits.h>

#ifdef OQP_CUDA_ENABLE
#include <cuda_runtime.h>
#endif

/* ==========================================================================
 * Per-region timing instrumentation (C1c benchmark harness preparation).
 *
 * Host wall-clock accumulators for the once-per-pass / once-per-accumulation
 * device transfer regions, so the benchmark can report SEPARATED timing
 * regions (upload d3 / zero f3 / download f3) instead of one opaque total.
 * The per-flush contraction hot path is deliberately NOT timed here -- timing
 * every flush would pollute the hot loop and the residency invariants.
 *
 * These are host wall times around the device calls, NOT a substitute for
 * on-device CUDA-event timing.  In a non-CUDA build they simply measure the
 * host memcpy fallbacks.  No speedup is implied or reported from this file.
 * ========================================================================== */
#define OQP_GPU_METC_T_UPLOAD_D3    0
#define OQP_GPU_METC_T_ZERO_F3      1
#define OQP_GPU_METC_T_DOWNLOAD_F3  2
#define OQP_GPU_METC_T_UPLOAD_BATCH 3
#define OQP_GPU_METC_T_ALLOC        4
#define OQP_GPU_METC_T_FREE         5
#define OQP_GPU_METC_NTIMERS        6

static double g_metc_secs[OQP_GPU_METC_NTIMERS];
static long   g_metc_calls[OQP_GPU_METC_NTIMERS];
static long long g_metc_h2d_bytes;
static long long g_metc_d2h_bytes;
static long long g_metc_device_bytes_live;
static long long g_metc_device_bytes_peak;
static long g_metc_h2d_calls;
static long g_metc_d2h_calls;
static long g_metc_malloc_calls;
static long g_metc_free_calls;

static int oqp_gpu_metc_env_flag(const char *name) {
  const char *v = getenv(name);
  return v && (strcmp(v, "1") == 0 || strcmp(v, "true") == 0 ||
               strcmp(v, "TRUE") == 0 || strcmp(v, "on") == 0 ||
               strcmp(v, "ON") == 0 || strcmp(v, "yes") == 0 ||
               strcmp(v, "YES") == 0);
}

static int oqp_gpu_metc_profile_requested(void) {
  return oqp_gpu_metc_env_flag("OQP_GPU_PROFILE") ||
         oqp_gpu_metc_env_flag("OQP_GPU_METC_PROFILE");
}

static double oqp_ws_now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void oqp_ws_timing_add(int idx, double t0) {
  g_metc_secs[idx] += oqp_ws_now_sec() - t0;
  g_metc_calls[idx] += 1;
}

/* Copy up to n region times (seconds) and call counts into the caller's
 * arrays.  Returns the number of entries written.  Either pointer may be
 * NULL.  Index order: upload_d3, zero_f3, download_f3. */
int oqp_gpu_metc_get_timings(double *secs_out, long *count_out, int n) {
  int m = (n < OQP_GPU_METC_NTIMERS) ? n : OQP_GPU_METC_NTIMERS;
  int i;
  if (m < 0) m = 0;
  for (i = 0; i < m; ++i) {
    if (secs_out)  secs_out[i]  = g_metc_secs[i];
    if (count_out) count_out[i] = g_metc_calls[i];
  }
  return m;
}

/* Reset all region timers/counters (call before a benchmark region). */
void oqp_gpu_metc_reset_timings(void) {
  int i;
  for (i = 0; i < OQP_GPU_METC_NTIMERS; ++i) {
    g_metc_secs[i] = 0.0;
    g_metc_calls[i] = 0;
  }
  g_metc_h2d_bytes = 0;
  g_metc_d2h_bytes = 0;
  g_metc_h2d_calls = 0;
  g_metc_d2h_calls = 0;
  g_metc_malloc_calls = 0;
  g_metc_free_calls = 0;
  g_metc_device_bytes_live = 0;
  g_metc_device_bytes_peak = 0;
}

#ifdef OQP_CUDA_ENABLE
extern void oqp_gpu_metc_cuda_profile_reset(void);
extern void oqp_gpu_metc_cuda_profile_snapshot(long *launches, double *kernel_ms,
                                               double *sync_s);
extern void oqp_gpu_metc_cuda_profile_dump(void);
#endif

static void oqp_gpu_metc_profile_reset_all(void) {
  oqp_gpu_metc_reset_timings();
#ifdef OQP_CUDA_ENABLE
  oqp_gpu_metc_cuda_profile_reset();
#endif
}

static void oqp_gpu_metc_note_alloc(int64_t bytes, double t0) {
  oqp_ws_timing_add(OQP_GPU_METC_T_ALLOC, t0);
  g_metc_malloc_calls += 1;
  g_metc_device_bytes_live += (long long)bytes;
  if (g_metc_device_bytes_live > g_metc_device_bytes_peak) {
    g_metc_device_bytes_peak = g_metc_device_bytes_live;
  }
}

static void oqp_gpu_metc_note_free(int64_t bytes, double t0) {
  oqp_ws_timing_add(OQP_GPU_METC_T_FREE, t0);
  g_metc_free_calls += 1;
  g_metc_device_bytes_live -= (long long)bytes;
  if (g_metc_device_bytes_live < 0) g_metc_device_bytes_live = 0;
}

/* Status codes (0 == success), shared in spirit with the Fortran bridge ierr. */
#define OQP_GPU_WS_OK         0
#define OQP_GPU_WS_ERR_HANDLE 1
#define OQP_GPU_WS_ERR_RANGE  2
#define OQP_GPU_WS_ERR_ARGS   3

#define OQP_GPU_WS_MAX_SLOTS 64

typedef struct {
  void   *base;
  int64_t total_bytes;
  int     in_use;
  int     is_device;
} oqp_gpu_ws_slot_t;

/* Fixed-size handle table.  Handles are 1-based slot indices; 0 == invalid. */
static oqp_gpu_ws_slot_t g_slots[OQP_GPU_WS_MAX_SLOTS];

static oqp_gpu_ws_slot_t *oqp_gpu_ws_slot_of(int handle) {
  if (handle <= 0 || handle > OQP_GPU_WS_MAX_SLOTS) {
    return NULL;
  }
  oqp_gpu_ws_slot_t *slot = &g_slots[handle - 1];
  if (!slot->in_use) {
    return NULL;
  }
  return slot;
}

/* acquire/allocate: reserve a contiguous arena of total_bytes for `target`
 * (0 = metc, 1 = xc_response).  Returns a positive handle on success, 0 on
 * failure (bad args, table full, or out of memory). */
int oqp_gpu_ws_acquire(int target, int64_t total_bytes) {
  if (total_bytes <= 0) {
    return 0;
  }
  if (target != 0 && target != 1) {
    return 0;
  }

  int slot_index = -1;
  for (int i = 0; i < OQP_GPU_WS_MAX_SLOTS; ++i) {
    if (!g_slots[i].in_use) {
      slot_index = i;
      break;
    }
  }
  if (slot_index < 0) {
    return 0;
  }

  double t_alloc = oqp_ws_now_sec();
  void *base = NULL;
  int is_device = 0;
#ifdef OQP_CUDA_ENABLE
  if (cudaMalloc(&base, (size_t)total_bytes) != cudaSuccess) {
    return 0;
  }
  is_device = 1;
  oqp_gpu_metc_note_alloc(total_bytes, t_alloc);
#else
  base = malloc((size_t)total_bytes);
  if (base == NULL) {
    return 0;
  }
  oqp_gpu_metc_note_alloc(total_bytes, t_alloc);
#endif

  g_slots[slot_index].base = base;
  g_slots[slot_index].total_bytes = total_bytes;
  g_slots[slot_index].in_use = 1;
  g_slots[slot_index].is_device = is_device;
  return slot_index + 1; /* 1-based handle */
}

/* validate: confirm [offset, offset+nbytes) lies within the handle's arena. */
int oqp_gpu_ws_validate(int handle, int64_t offset, int64_t nbytes) {
  oqp_gpu_ws_slot_t *slot = oqp_gpu_ws_slot_of(handle);
  if (slot == NULL) {
    return OQP_GPU_WS_ERR_HANDLE;
  }
  if (nbytes <= 0 || offset < 0) {
    return OQP_GPU_WS_ERR_ARGS;
  }
  if (offset > slot->total_bytes || nbytes > slot->total_bytes - offset) {
    return OQP_GPU_WS_ERR_RANGE;
  }
  return OQP_GPU_WS_OK;
}

/* ptr/borrow: return base + offset, or NULL on invalid handle / out-of-range. */
void *oqp_gpu_ws_ptr(int handle, int64_t offset) {
  oqp_gpu_ws_slot_t *slot = oqp_gpu_ws_slot_of(handle);
  if (slot == NULL) {
    return NULL;
  }
  if (offset < 0 || offset >= slot->total_bytes) {
    return NULL;
  }
  return (void *)((char *)slot->base + offset);
}

/* total_bytes accessor: arena size for the handle, or -1 if invalid. */
int64_t oqp_gpu_ws_total_bytes(int handle) {
  oqp_gpu_ws_slot_t *slot = oqp_gpu_ws_slot_of(handle);
  if (slot == NULL) {
    return -1;
  }
  return slot->total_bytes;
}

/* release: free the arena.  Double release is safe: a released (or never
 * allocated) handle is not in_use, so this returns OQP_GPU_WS_ERR_HANDLE. */
int oqp_gpu_ws_release(int handle) {
  oqp_gpu_ws_slot_t *slot = oqp_gpu_ws_slot_of(handle);
  if (slot == NULL) {
    return OQP_GPU_WS_ERR_HANDLE;
  }
  double t_free = oqp_ws_now_sec();
#ifdef OQP_CUDA_ENABLE
  if (slot->is_device) {
    int rc = (int)cudaFree(slot->base);
    if (rc != 0) return rc;
  } else {
    free(slot->base);
  }
#else
  free(slot->base);
#endif
  oqp_gpu_metc_note_free(slot->total_bytes, t_free);
  slot->base = NULL;
  slot->total_bytes = 0;
  slot->in_use = 0;
  slot->is_device = 0;
  return OQP_GPU_WS_OK;
}

/* ==========================================================================
 * METC residency session (METC-C1b).
 *
 * A METC session owns one contiguous workspace arena (via oqp_gpu_ws_acquire)
 * laid out as four regions:
 *
 *   [ d3 | f3[0..nthreads) | ids[0..nthreads) | ints[0..nthreads) ]
 *
 *   - d3   : shared, read-only density/input tensor (per-pass refresh).
 *   - f3   : PER_THREAD accumulator (zeroed once at accumulation start,
 *            downloaded once at the end).
 *   - ids  : PER_THREAD index scratch, sized for max_ncur (concurrency-safe).
 *   - ints : PER_THREAD integral scratch, sized for max_ncur.
 *
 * Per-thread regions give every concurrent OpenMP `update` call its own scratch
 * and f3 slice, so no resident buffer is shared-written across threads.  The
 * hot path performs NO cudaMalloc/cudaFree; it borrows these pointers and only
 * copies the current flush's ids/ints into the calling thread's scratch.
 * ========================================================================== */

#define OQP_GPU_METC_MAX_SESSIONS 64
#define OQP_GPU_METC_ERR_SESSION 5
#define OQP_GPU_METC_ERR_THREAD  6
#define OQP_GPU_METC_ERR_CAP     7
#define OQP_GPU_METC_ERR_F3_CHECK 8

#ifdef OQP_CUDA_ENABLE
/* Borrowed-pointer kernel launcher, implemented in gpu_metc_cuda.cu.  It does
 * NOT allocate, copy, or free device memory -- it only launches the kernel. */
extern int oqp_gpu_metc_launch(const int *d_ids, const double *d_ints, int ncur,
                               double *d_f3, const double *d_d3, int nf,
                               int nmatrix, int nbf, int cur_pass,
                               double scale_exchange, double scale_coulomb,
                               int is_umrsf);
#endif

typedef struct {
  int     in_use;
  int     ws_handle;     /* underlying workspace arena handle */
  int     nthreads;
  int     max_ncur;      /* resident scratch capacity (integrals per flush) */
  int     nf, nmatrix, nbf;
  int64_t bytes_tensor;  /* bytes of one d3/f3 tensor slice */
  int64_t ids_cap_bytes; /* bytes of one thread's ids scratch */
  int64_t ints_cap_bytes;/* bytes of one thread's ints scratch */
  int64_t off_d3, off_f3, off_ids, off_ints;
  int64_t total_bytes;
  int     f3_check_enabled;
  double *f3_check_ref;
  double *f3_check_d3;
  double  session_start_s;
  double  contract_wall_s;
  long    contract_calls;
} oqp_gpu_metc_session_t;

static oqp_gpu_metc_session_t g_metc[OQP_GPU_METC_MAX_SESSIONS];


static int oqp_gpu_metc_f3_check_requested(void) {
  return oqp_gpu_metc_env_flag("OQP_GPU_F3_CHECK") ||
         oqp_gpu_metc_env_flag("OQP_GPU_METC_F3_CHECK");
}

static double oqp_gpu_metc_env_double(const char *name, double fallback) {
  const char *v = getenv(name);
  char *end = NULL;
  if (!v || !*v) return fallback;
  double x = strtod(v, &end);
  return (end && end != v) ? x : fallback;
}

static inline int idx4_host(int f, int m, int row, int col,
                            int nf, int nmatrix, int nbf) {
  return f + nf * (m + nmatrix * (row + nbf * col));
}

static inline void add4_host(double *f3, int f, int m, int row, int col,
                             int nf, int nmatrix, int nbf, double value) {
  f3[idx4_host(f, m, row, col, nf, nmatrix, nbf)] += value;
}

static inline double get4_host(const double *d3, int f, int m, int row, int col,
                               int nf, int nmatrix, int nbf) {
  return d3[idx4_host(f, m, row, col, nf, nmatrix, nbf)];
}

static void oqp_gpu_metc_f3_check_accumulate(oqp_gpu_metc_session_t *s, int thread,
                                             const int *ids, const double *ints,
                                             int ncur, int cur_pass,
                                             double scale_exchange,
                                             double scale_coulomb, int is_umrsf) {
  if (!s || !s->f3_check_enabled || !s->f3_check_ref || !s->f3_check_d3) return;
  double *f3 = s->f3_check_ref + (int64_t)thread * (s->bytes_tensor / (int64_t)8);
  const double *d3 = s->f3_check_d3;
  int nf = s->nf, nmatrix = s->nmatrix, nbf = s->nbf;
  for (int n = 0; n < ncur; ++n) {
    int i = ids[4 * n + 0] - 1;
    int j = ids[4 * n + 1] - 1;
    int k = ids[4 * n + 2] - 1;
    int l = ids[4 * n + 3] - 1;
    double val = ints[n];
    double xval = val * scale_exchange;
    double cval = val * scale_coulomb;
    if (i < 0 || j < 0 || k < 0 || l < 0 || i >= nbf || j >= nbf || k >= nbf || l >= nbf) continue;
    for (int f = 0; f < nf; ++f) {
      for (int m = 0; m < nmatrix; ++m) {
        if (!is_umrsf) {
          if (cur_pass == 1) {
            if (m < 4) {
              add4_host(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4_host(d3, f, m, k, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4_host(d3, f, m, i, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4_host(d3, f, m, l, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4_host(d3, f, m, i, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4_host(d3, f, m, k, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4_host(d3, f, m, j, i, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4_host(d3, f, m, l, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4_host(d3, f, m, j, i, nf, nmatrix, nbf));
            }
            if (m < 7) {
              add4_host(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, i, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, i, nf, nmatrix, nbf));
            }
          } else if (cur_pass == 2 && m == 6) {
            add4_host(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, l, nf, nmatrix, nbf));
            add4_host(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, j, nf, nmatrix, nbf));
            add4_host(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, k, nf, nmatrix, nbf));
            add4_host(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, j, nf, nmatrix, nbf));
            add4_host(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, l, nf, nmatrix, nbf));
            add4_host(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, i, nf, nmatrix, nbf));
            add4_host(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, k, nf, nmatrix, nbf));
            add4_host(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, i, nf, nmatrix, nbf));
          }
        } else {
          if (cur_pass == 1) {
            if (m < 8) {
              add4_host(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4_host(d3, f, m, k, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4_host(d3, f, m, i, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4_host(d3, f, m, l, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4_host(d3, f, m, i, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4_host(d3, f, m, k, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4_host(d3, f, m, j, i, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4_host(d3, f, m, l, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4_host(d3, f, m, j, i, nf, nmatrix, nbf));
              add4_host(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, i, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, i, nf, nmatrix, nbf));
            }
            if (m == 8 || m == 9) {
              add4_host(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, i, nf, nmatrix, nbf));
              add4_host(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, i, nf, nmatrix, nbf));
            }
            if (m == 10) {
              add4_host(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, j, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, l, nf, nmatrix, nbf));
              add4_host(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, i, nf, nmatrix, nbf));
              add4_host(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, k, nf, nmatrix, nbf));
              add4_host(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, i, nf, nmatrix, nbf));
            }
          } else if (cur_pass == 2 && m == 10) {
            add4_host(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, l, nf, nmatrix, nbf));
            add4_host(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, j, nf, nmatrix, nbf));
            add4_host(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, j, k, nf, nmatrix, nbf));
            add4_host(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, j, nf, nmatrix, nbf));
            add4_host(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, l, nf, nmatrix, nbf));
            add4_host(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, l, i, nf, nmatrix, nbf));
            add4_host(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, i, k, nf, nmatrix, nbf));
            add4_host(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4_host(d3, f, m, k, i, nf, nmatrix, nbf));
          }
        }
      }
    }
  }
}

static int oqp_gpu_metc_f3_check_compare(oqp_gpu_metc_session_t *s, const double *gpu_f3) {
  if (!s || !s->f3_check_enabled) return OQP_GPU_WS_OK;
  if (!s->f3_check_ref || !gpu_f3) {
    fprintf(stderr, "OQP_GPU_METC_F3_CHECK passed=false notes=missing_cpu_or_gpu_f3\n");
    return OQP_GPU_METC_ERR_F3_CHECK;
  }
  int64_t n = (int64_t)s->nthreads * (s->bytes_tensor / (int64_t)8);
  double norm_cpu2 = 0.0, norm_gpu2 = 0.0, norm_delta2 = 0.0;
  double max_abs = -1.0;
  int64_t max_idx = -1;
  for (int64_t q = 0; q < n; ++q) {
    double c = s->f3_check_ref[q];
    double g = gpu_f3[q];
    if (!isfinite(c) || !isfinite(g)) {
      fprintf(stderr, "OQP_GPU_METC_F3_CHECK passed=false notes=nan_or_inf index=%" PRId64 " cpu=%.17g gpu=%.17g\n", q, c, g);
      return OQP_GPU_METC_ERR_F3_CHECK;
    }
    double d = g - c;
    double ad = fabs(d);
    norm_cpu2 += c * c;
    norm_gpu2 += g * g;
    norm_delta2 += d * d;
    if (ad > max_abs) { max_abs = ad; max_idx = q; }
  }
  double norm_cpu = sqrt(norm_cpu2);
  double norm_gpu = sqrt(norm_gpu2);
  double norm_delta = sqrt(norm_delta2);
  double rms = n > 0 ? sqrt(norm_delta2 / (double)n) : 0.0;
  double atol = oqp_gpu_metc_env_double("OQP_GPU_METC_F3_CHECK_ATOL", 1.0e-7);
  double rtol = oqp_gpu_metc_env_double("OQP_GPU_METC_F3_CHECK_RTOL", 1.0e-7);
  int passed = (max_abs <= atol + rtol * fmax(1.0, norm_cpu));
  int64_t per_thread = s->bytes_tensor / (int64_t)8;
  int thread = (per_thread > 0 && max_idx >= 0) ? (int)(max_idx / per_thread) : -1;
  int64_t local = (per_thread > 0 && max_idx >= 0) ? (max_idx % per_thread) : -1;
  int f = -1, m = -1, row = -1, col = -1;
  if (local >= 0) {
    f = (int)(local % s->nf);
    int64_t tmp = local / s->nf;
    m = (int)(tmp % s->nmatrix);
    tmp = tmp / s->nmatrix;
    row = (int)(tmp % s->nbf);
    col = (int)(tmp / s->nbf);
  }
  fprintf(stderr,
          "OQP_GPU_METC_F3_CHECK nf=%d nmatrix=%d nbf=%d nthreads=%d n_elements=%" PRId64
          " norm_cpu=%.17g norm_gpu=%.17g norm_delta=%.17g max_abs_error=%.17g"
          " rms_abs_error=%.17g max_error_index=%" PRId64 ":%d:%d:%d:%d:%d"
          " atol=%.3g rtol=%.3g passed=%s notes=%s\n",
          s->nf, s->nmatrix, s->nbf, s->nthreads, n, norm_cpu, norm_gpu,
          norm_delta, max_abs, rms, max_idx, thread, f, m, row, col,
          atol, rtol, passed ? "true" : "false", passed ? "ok" : "exceeds_tolerance");
  return passed ? OQP_GPU_WS_OK : OQP_GPU_METC_ERR_F3_CHECK;
}

static oqp_gpu_metc_session_t *oqp_gpu_metc_of(int session) {
  if (session <= 0 || session > OQP_GPU_METC_MAX_SESSIONS) {
    return NULL;
  }
  oqp_gpu_metc_session_t *s = &g_metc[session - 1];
  if (!s->in_use) {
    return NULL;
  }
  return s;
}

/* Compute the contiguous arena layout for a METC session.  Pure arithmetic so
 * it can be cross-checked against the Python manifest without any allocation. */
void oqp_gpu_metc_layout(int nf, int nmatrix, int nbf, int nthreads,
                         int max_ncur, int64_t *off_d3, int64_t *off_f3,
                         int64_t *off_ids, int64_t *off_ints,
                         int64_t *total_bytes) {
  /* The kernel uses 32-bit tensor offsets; reject shapes before multiplication. */
  if (nf <= 0 || nmatrix <= 0 || nbf <= 0 || nthreads <= 0 || max_ncur <= 0 ||
      nf > INT_MAX / nmatrix || nbf > INT_MAX / (nf * nmatrix) ||
      nbf > INT_MAX / (nf * nmatrix * nbf) ||
      (int64_t)nthreads + 1 > INT64_MAX / ((int64_t)nf * nmatrix * nbf * nbf * 8) ||
      (int64_t)nthreads > INT64_MAX / ((int64_t)max_ncur * 24)) {
    if (off_d3) *off_d3 = -1;
    if (off_f3) *off_f3 = -1;
    if (off_ids) *off_ids = -1;
    if (off_ints) *off_ints = -1;
    if (total_bytes) *total_bytes = -1;
    return;
  }
  int64_t bytes_tensor = (int64_t)nf * nmatrix * nbf * nbf * (int64_t)8;
  int64_t ids_cap = (int64_t)4 * max_ncur * (int64_t)4;   /* 4 int32 per integral */
  int64_t ints_cap = (int64_t)max_ncur * (int64_t)8;      /* 1 double per integral */
  if (bytes_tensor * ((int64_t)nthreads + 1) >
      INT64_MAX - (int64_t)nthreads * (ids_cap + ints_cap)) {
    if (total_bytes) *total_bytes = -1;
    return;
  }
  int64_t d3 = 0;
  int64_t f3 = d3 + bytes_tensor;
  int64_t ids = f3 + (int64_t)nthreads * bytes_tensor;
  int64_t ints = ids + (int64_t)nthreads * ids_cap;
  int64_t total = ints + (int64_t)nthreads * ints_cap;
  if (off_d3) *off_d3 = d3;
  if (off_f3) *off_f3 = f3;
  if (off_ids) *off_ids = ids;
  if (off_ints) *off_ints = ints;
  if (total_bytes) *total_bytes = total;
}

/* begin: acquire the arena and record the layout.  Returns a positive session
 * id on success, 0 on failure.  Does NOT upload d3 or zero f3. */
int oqp_gpu_metc_session_begin(int nf, int nmatrix, int nbf, int nthreads,
                               int max_ncur) {
  const char *induce_begin_fail = getenv("OQP_GPU_METC_INDUCE_BEGIN_FAIL");
  if (induce_begin_fail && strcmp(induce_begin_fail, "1") == 0) {
    return 0;
  }
  if (nf <= 0 || nmatrix <= 0 || nbf <= 0 || nthreads <= 0 || max_ncur <= 0) {
    return 0;
  }
  oqp_gpu_metc_profile_reset_all();
  int slot = -1;
  for (int i = 0; i < OQP_GPU_METC_MAX_SESSIONS; ++i) {
    if (!g_metc[i].in_use) {
      slot = i;
      break;
    }
  }
  if (slot < 0) {
    return 0;
  }

  oqp_gpu_metc_session_t *s = &g_metc[slot];
  oqp_gpu_metc_layout(nf, nmatrix, nbf, nthreads, max_ncur, &s->off_d3,
                      &s->off_f3, &s->off_ids, &s->off_ints, &s->total_bytes);

  if (s->total_bytes <= 0) return 0;
  int ws = oqp_gpu_ws_acquire(0 /* metc target */, s->total_bytes);
  if (ws == 0) {
    return 0;
  }

  s->ws_handle = ws;
  s->nthreads = nthreads;
  s->max_ncur = max_ncur;
  s->nf = nf;
  s->nmatrix = nmatrix;
  s->nbf = nbf;
  s->bytes_tensor = (int64_t)nf * nmatrix * nbf * nbf * (int64_t)8;
  s->ids_cap_bytes = (int64_t)4 * max_ncur * (int64_t)4;
  s->ints_cap_bytes = (int64_t)max_ncur * (int64_t)8;
  s->f3_check_enabled = oqp_gpu_metc_f3_check_requested();
  s->f3_check_ref = NULL;
  s->f3_check_d3 = NULL;
  s->session_start_s = oqp_ws_now_sec();
  s->contract_wall_s = 0.0;
  s->contract_calls = 0;
  if (s->f3_check_enabled) {
    size_t ref_n = (size_t)((int64_t)nthreads * s->bytes_tensor);
    s->f3_check_ref = (double *)calloc(1, ref_n);
    s->f3_check_d3 = (double *)malloc((size_t)s->bytes_tensor);
    if (!s->f3_check_ref || !s->f3_check_d3) {
      free(s->f3_check_ref);
      free(s->f3_check_d3);
      s->f3_check_ref = NULL;
      s->f3_check_d3 = NULL;
      oqp_gpu_ws_release(ws);
      return 0;
    }
  }
  s->in_use = 1;
  return slot + 1;
}

int64_t oqp_gpu_metc_session_total_bytes(int session) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  return s ? s->total_bytes : (int64_t)-1;
}

int oqp_gpu_metc_session_capacity(int session) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  return s ? s->max_ncur : -1;
}

/* Shared d3 input pointer (same for every thread). */
void *oqp_gpu_metc_d3_ptr(int session) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return NULL;
  return oqp_gpu_ws_ptr(s->ws_handle, s->off_d3);
}

/* Per-thread f3 accumulator slice (thread in [0, nthreads)). */
void *oqp_gpu_metc_f3_ptr(int session, int thread) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s || thread < 0 || thread >= s->nthreads) return NULL;
  return oqp_gpu_ws_ptr(s->ws_handle, s->off_f3 + (int64_t)thread * s->bytes_tensor);
}

/* Per-thread ids scratch slice. */
void *oqp_gpu_metc_ids_ptr(int session, int thread) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s || thread < 0 || thread >= s->nthreads) return NULL;
  return oqp_gpu_ws_ptr(s->ws_handle, s->off_ids + (int64_t)thread * s->ids_cap_bytes);
}

/* Per-thread ints scratch slice. */
void *oqp_gpu_metc_ints_ptr(int session, int thread) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s || thread < 0 || thread >= s->nthreads) return NULL;
  return oqp_gpu_ws_ptr(s->ws_handle, s->off_ints + (int64_t)thread * s->ints_cap_bytes);
}

/* Reject a flush whose integral count exceeds the resident scratch capacity. */
int oqp_gpu_metc_session_check_ncur(int session, int ncur) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return OQP_GPU_METC_ERR_SESSION;
  if (ncur < 0 || ncur > s->max_ncur) return OQP_GPU_METC_ERR_CAP;
  return OQP_GPU_WS_OK;
}

static int oqp_gpu_ws_memset_like(void *dst, int value, size_t n, int is_device) {
#ifdef OQP_CUDA_ENABLE
  if (is_device) return (int)cudaMemset(dst, value, n);
#else
  (void)is_device;
#endif
  memset(dst, value, n);
  return 0;
}

static int oqp_gpu_ws_copy_h2d(void *dst, const void *src, size_t n, int is_device) {
  g_metc_h2d_calls += 1;
  g_metc_h2d_bytes += (long long)n;
#ifdef OQP_CUDA_ENABLE
  if (is_device) return (int)cudaMemcpy(dst, src, n, cudaMemcpyHostToDevice);
#else
  (void)is_device;
#endif
  memcpy(dst, src, n);
  return 0;
}

static int oqp_gpu_ws_copy_d2h(void *dst, const void *src, size_t n, int is_device) {
  g_metc_d2h_calls += 1;
  g_metc_d2h_bytes += (long long)n;
#ifdef OQP_CUDA_ENABLE
  if (is_device) return (int)cudaMemcpy(dst, src, n, cudaMemcpyDeviceToHost);
#else
  (void)is_device;
#endif
  memcpy(dst, src, n);
  return 0;
}

static int oqp_gpu_metc_is_device(int session) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return 0;
  oqp_gpu_ws_slot_t *slot = oqp_gpu_ws_slot_of(s->ws_handle);
  return slot ? slot->is_device : 0;
}

/* Zero the whole per-thread f3 region once at the start of accumulation. */
int oqp_gpu_metc_session_zero_f3(int session) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return OQP_GPU_METC_ERR_SESSION;
  void *f3 = oqp_gpu_ws_ptr(s->ws_handle, s->off_f3);
  if (!f3) return OQP_GPU_WS_ERR_HANDLE;
  size_t n = (size_t)((int64_t)s->nthreads * s->bytes_tensor);
  double t0 = oqp_ws_now_sec();
  int rc = oqp_gpu_ws_memset_like(f3, 0, n, oqp_gpu_metc_is_device(session));
  if (rc != 0) return rc;
  if (s->f3_check_ref) memset(s->f3_check_ref, 0, n);
  oqp_ws_timing_add(OQP_GPU_METC_T_ZERO_F3, t0);
  return OQP_GPU_WS_OK;
}

/* Upload/refresh the shared d3 tensor (once per pass). */
int oqp_gpu_metc_session_upload_d3(int session, const double *d3_host) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return OQP_GPU_METC_ERR_SESSION;
  if (!d3_host) return OQP_GPU_WS_ERR_ARGS;
  void *d3 = oqp_gpu_metc_d3_ptr(session);
  if (!d3) return OQP_GPU_WS_ERR_HANDLE;
  double t0 = oqp_ws_now_sec();
  int rc = oqp_gpu_ws_copy_h2d(d3, d3_host, (size_t)s->bytes_tensor, oqp_gpu_metc_is_device(session));
  if (rc != 0) return rc;
  if (s->f3_check_enabled && s->f3_check_d3) {
    memcpy(s->f3_check_d3, d3_host, (size_t)s->bytes_tensor);
  }
  oqp_ws_timing_add(OQP_GPU_METC_T_UPLOAD_D3, t0);
  return OQP_GPU_WS_OK;
}

/* Copy one flush's ids/ints into the calling thread's resident scratch. */
int oqp_gpu_metc_session_upload_batch(int session, int thread, const int *ids,
                                      const double *ints, int ncur) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return OQP_GPU_METC_ERR_SESSION;
  if (thread < 0 || thread >= s->nthreads) return OQP_GPU_METC_ERR_THREAD;
  int cap = oqp_gpu_metc_session_check_ncur(session, ncur);
  if (cap != OQP_GPU_WS_OK) return cap;
  if (ncur == 0) return OQP_GPU_WS_OK;
  if (!ids || !ints) return OQP_GPU_WS_ERR_ARGS;
  int is_device = oqp_gpu_metc_is_device(session);
  double t0 = oqp_ws_now_sec();
  int rc = oqp_gpu_ws_copy_h2d(oqp_gpu_metc_ids_ptr(session, thread), ids,
                      (size_t)((int64_t)4 * ncur * (int64_t)4), is_device);
  if (rc != 0) return rc;
  rc = oqp_gpu_ws_copy_h2d(oqp_gpu_metc_ints_ptr(session, thread), ints,
                      (size_t)((int64_t)ncur * (int64_t)8), is_device);
  if (rc != 0) return rc;
  oqp_ws_timing_add(OQP_GPU_METC_T_UPLOAD_BATCH, t0);
  return OQP_GPU_WS_OK;
}

/* Hot-path contract for one flush on one thread: copy ids/ints into the
 * thread's resident scratch, then launch the kernel against borrowed d3/f3
 * pointers.  No cudaMalloc/cudaFree.  Under the host fallback the kernel is not
 * executed (OQP_CUDA_ENABLE undefined); the scratch copy still happens so the
 * data-movement bookkeeping is testable. */
int oqp_gpu_metc_session_contract(int session, int thread, const int *ids,
                                  const double *ints, int ncur, int cur_pass,
                                  double scale_exchange, double scale_coulomb,
                                  int is_umrsf) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return OQP_GPU_METC_ERR_SESSION;
  if (s->nmatrix != (is_umrsf ? 11 : 7) || (cur_pass != 1 && cur_pass != 2))
    return OQP_GPU_WS_ERR_ARGS;
  if (ncur > 0 && (!ids || !ints)) return OQP_GPU_WS_ERR_ARGS;
  int capacity = oqp_gpu_metc_session_check_ncur(session, ncur);
  if (capacity != 0) return capacity;
  for (int64_t i = 0; i < (int64_t)ncur * 4; ++i)
    if (ids[i] < 1 || ids[i] > s->nbf) return OQP_GPU_WS_ERR_ARGS;
  double t_contract = oqp_ws_now_sec();
  int up = oqp_gpu_metc_session_upload_batch(session, thread, ids, ints, ncur);
  if (up != OQP_GPU_WS_OK) return up;
  if (ncur == 0) return OQP_GPU_WS_OK;
#ifdef OQP_CUDA_ENABLE
  int rc = oqp_gpu_metc_launch(
      (const int *)oqp_gpu_metc_ids_ptr(session, thread),
      (const double *)oqp_gpu_metc_ints_ptr(session, thread), ncur,
      (double *)oqp_gpu_metc_f3_ptr(session, thread),
      (const double *)oqp_gpu_metc_d3_ptr(session), s->nf, s->nmatrix, s->nbf,
      cur_pass, scale_exchange, scale_coulomb, is_umrsf);
  s->contract_wall_s += oqp_ws_now_sec() - t_contract;
  s->contract_calls += 1;
  if (rc == OQP_GPU_WS_OK) {
    oqp_gpu_metc_f3_check_accumulate(s, thread, ids, ints, ncur, cur_pass,
                                     scale_exchange, scale_coulomb, is_umrsf);
  }
  return rc;
#else
  (void)cur_pass;
  (void)scale_exchange;
  (void)scale_coulomb;
  (void)is_umrsf;
  return 9; /* no CUDA backend: never report a no-op as a contraction */
#endif
}

/* Download the full per-thread f3 region to host once at accumulation end.
 * The host f3 array keeps the existing (nfocks,nmatrix,nbf,nbf,nthreads)
 * layout, so the existing Fortran thread reduction still applies. */
int oqp_gpu_metc_session_download_f3(int session, double *f3_host) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return OQP_GPU_METC_ERR_SESSION;
  if (!f3_host) return OQP_GPU_WS_ERR_ARGS;
  void *f3 = oqp_gpu_ws_ptr(s->ws_handle, s->off_f3);
  if (!f3) return OQP_GPU_WS_ERR_HANDLE;
  size_t n = (size_t)((int64_t)s->nthreads * s->bytes_tensor);
  double t0 = oqp_ws_now_sec();
  int rc = oqp_gpu_ws_copy_d2h(f3_host, f3, n, oqp_gpu_metc_is_device(session));
  if (rc != 0) return rc;
  oqp_ws_timing_add(OQP_GPU_METC_T_DOWNLOAD_F3, t0);
  return oqp_gpu_metc_f3_check_compare(s, f3_host);
}

/* Release the arena and free the session slot.  Double end is safe. */
int oqp_gpu_metc_session_end(int session) {
  oqp_gpu_metc_session_t *s = oqp_gpu_metc_of(session);
  if (!s) return OQP_GPU_METC_ERR_SESSION;
  int rc = oqp_gpu_ws_release(s->ws_handle);
  s->ws_handle = 0;
  /* Optional separated-region timing dump for the benchmark harness. */
  if (getenv("OQP_GPU_METC_TIMING")) {
    fprintf(stderr,
            "OQP_GPU_METC_TIMING upload_d3_s=%.6f zero_f3_s=%.6f download_f3_s=%.6f"
            " upload_d3_n=%ld zero_f3_n=%ld download_f3_n=%ld\n",
            g_metc_secs[OQP_GPU_METC_T_UPLOAD_D3],
            g_metc_secs[OQP_GPU_METC_T_ZERO_F3],
            g_metc_secs[OQP_GPU_METC_T_DOWNLOAD_F3],
            g_metc_calls[OQP_GPU_METC_T_UPLOAD_D3],
            g_metc_calls[OQP_GPU_METC_T_ZERO_F3],
            g_metc_calls[OQP_GPU_METC_T_DOWNLOAD_F3]);
  }
  if (oqp_gpu_metc_profile_requested()) {
    long kernel_launches = 0;
    double kernel_ms = 0.0;
    double sync_s = 0.0;
#ifdef OQP_CUDA_ENABLE
    oqp_gpu_metc_cuda_profile_snapshot(&kernel_launches, &kernel_ms, &sync_s);
#endif
    double session_wall_s = oqp_ws_now_sec() - s->session_start_s;
    double contract_host_overhead_s = s->contract_wall_s - sync_s;
    if (contract_host_overhead_s < 0.0) contract_host_overhead_s = 0.0;
    fprintf(stderr,
            "OQP_GPU_PROFILE schema=1 target=metc nf=%d nmatrix=%d nbf=%d nthreads=%d max_ncur=%d"
            " total_bytes=%" PRId64 " device_bytes_peak=%lld h2d_bytes=%lld d2h_bytes=%lld"
            " h2d_calls=%ld d2h_calls=%ld cuda_malloc_calls=%ld cuda_free_calls=%ld"
            " alloc_s=%.6f free_s=%.6f upload_d3_s=%.6f upload_batch_s=%.6f zero_f3_s=%.6f download_f3_s=%.6f"
            " session_wall_s=%.6f contract_wall_s=%.6f contract_calls=%ld contract_host_overhead_s=%.6f"
            " upload_d3_n=%ld upload_batch_n=%ld zero_f3_n=%ld download_f3_n=%ld"
            " kernel_launches=%ld kernel_s=%.6f sync_s=%.6f f3_check=%d"
            " timing_semantics=kernel_event_time_overlaps_sync_wait\n",
            s->nf, s->nmatrix, s->nbf, s->nthreads, s->max_ncur,
            s->total_bytes, g_metc_device_bytes_peak, g_metc_h2d_bytes,
            g_metc_d2h_bytes, g_metc_h2d_calls, g_metc_d2h_calls,
            g_metc_malloc_calls, g_metc_free_calls,
            g_metc_secs[OQP_GPU_METC_T_ALLOC], g_metc_secs[OQP_GPU_METC_T_FREE],
            g_metc_secs[OQP_GPU_METC_T_UPLOAD_D3],
            g_metc_secs[OQP_GPU_METC_T_UPLOAD_BATCH],
            g_metc_secs[OQP_GPU_METC_T_ZERO_F3],
            g_metc_secs[OQP_GPU_METC_T_DOWNLOAD_F3],
            session_wall_s, s->contract_wall_s, s->contract_calls,
            contract_host_overhead_s,
            g_metc_calls[OQP_GPU_METC_T_UPLOAD_D3],
            g_metc_calls[OQP_GPU_METC_T_UPLOAD_BATCH],
            g_metc_calls[OQP_GPU_METC_T_ZERO_F3],
            g_metc_calls[OQP_GPU_METC_T_DOWNLOAD_F3],
            kernel_launches, kernel_ms * 1.0e-3, sync_s, s->f3_check_enabled);
#ifdef OQP_CUDA_ENABLE
    oqp_gpu_metc_cuda_profile_dump();
#endif
  }
  free(s->f3_check_ref);
  free(s->f3_check_d3);
  s->f3_check_ref = NULL;
  s->f3_check_d3 = NULL;
  s->f3_check_enabled = 0;
  s->in_use = 0;
  s->ws_handle = 0;
  s->total_bytes = 0;
  return rc;
}
