/*
 * Runtime control of BLAS-library internal threading.
 *
 * Grid-based DFT (and other slice-parallel code) calls many small BLAS
 * operations from inside OpenMP parallel regions.  If the linked BLAS
 * library runs its own thread pool (e.g. pthread builds of OpenBLAS),
 * every such call is dispatched to that pool, which serializes on the
 * pool lock and oversubscribes the machine.  The fix is to switch the
 * BLAS library to single-threaded mode around such regions.
 *
 * The BLAS library is only known at link/run time, so the setter is
 * looked up dynamically with dlsym().  Supported: OpenBLAS, MKL, BLIS.
 * If no known symbol is found (e.g. reference BLAS, Apple Accelerate),
 * the calls are no-ops and we report "unknown" (-1).
 */

#define _GNU_SOURCE   /* RTLD_DEFAULT in <dlfcn.h> on glibc < 2.34 */

#include <stdint.h>

#if defined(_WIN32)

/* No dlsym on Windows -- compile to no-ops. */
int64_t oqp_blas_thread_count(void) { return -1; }
void oqp_blas_thread_set(int64_t n) { (void)n; }

#else

#include <dlfcn.h>
#include <pthread.h>

typedef void (*set_fn_t)(int);
typedef int (*get_fn_t)(void);

static set_fn_t set_fn;
static get_fn_t get_fn;
static pthread_once_t resolved = PTHREAD_ONCE_INIT;

static void
resolve_syms(void)
{
    /* OpenBLAS */
    set_fn = (set_fn_t) dlsym(RTLD_DEFAULT, "openblas_set_num_threads");
    get_fn = (get_fn_t) dlsym(RTLD_DEFAULT, "openblas_get_num_threads");
    if (set_fn && get_fn) return;

    /* MKL */
    set_fn = (set_fn_t) dlsym(RTLD_DEFAULT, "MKL_Set_Num_Threads");
    get_fn = (get_fn_t) dlsym(RTLD_DEFAULT, "MKL_Get_Max_Threads");
    if (set_fn && get_fn) return;

    /* BLIS */
    set_fn = (set_fn_t) dlsym(RTLD_DEFAULT, "bli_thread_set_num_threads");
    get_fn = (get_fn_t) dlsym(RTLD_DEFAULT, "bli_thread_get_num_threads");
    if (set_fn && get_fn) return;

    set_fn = 0;
    get_fn = 0;
}

/*
 * @brief Get the current BLAS thread count
 * @return thread count, or -1 if the BLAS library is not recognized
 */
int64_t
oqp_blas_thread_count(void)
{
    pthread_once(&resolved, resolve_syms);
    if (!get_fn) return -1;
    return get_fn();
}

/*
 * @brief Set the BLAS thread count; no-op for n < 1 or unknown BLAS
 */
void
oqp_blas_thread_set(int64_t n)
{
    pthread_once(&resolved, resolve_syms);
    if (!set_fn || n < 1) return;
    set_fn((int) n);
}

#endif

/*
 * @brief Whether liboqp was compiled with OpenMP support (1) or not (0).
 *
 * Lets the Python frontend warn that an `omp_threads` / OMP_NUM_THREADS request
 * has no effect on a serial build.  Compiled with the project's C flags, which
 * include -fopenmp when ENABLE_OPENMP is on, so _OPENMP reflects the build.
 */
int
oqp_have_openmp(void)
{
#ifdef _OPENMP
    return 1;
#else
    return 0;
#endif
}

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * @brief Set the OpenMP thread count at runtime (no-op for n < 1 or serial build).
 *
 * Lets the frontend honour an `omp_threads` request that arrives too late for the
 * OMP_NUM_THREADS environment variable -- in particular the programmatic
 * Runner(input_dict=...) path, where the value is only known after liboqp has
 * loaded.  omp_set_num_threads takes effect for subsequent parallel regions.
 */
void
oqp_omp_set_num_threads(int n)
{
#ifdef _OPENMP
    if (n >= 1) omp_set_num_threads(n);
#else
    (void) n;
#endif
}
