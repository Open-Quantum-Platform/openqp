/* Direct-integral METC backend; independent of the density-fitting engine.
 * Arrays use f + nf*(m + nmatrix*(row + nbf*col)); ids are 1-based int32.
 * All calls must be serialized by the caller. Session thread indices select
 * private accumulator slices; they do not make the global profiler thread-safe.
 * A session is bound to the current CUDA device from begin through end.
 * Negative/positive nonzero status means failure; never treat it as CPU fallback.
 */
#ifndef OPENQP_GPU_METC_H
#define OPENQP_GPU_METC_H
#include <stdint.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
int oqp_gpu_metc_device_count(void);
/* f3 is in/out: adds the contraction to its initial value. */
int oqp_gpu_metc_contract(const int *ids, const double *ints, int ncur,
    double *f3, const double *d3, int nf, int nmatrix, int nbf, int cur_pass,
    double scale_exchange, double scale_coulomb, bool is_umrsf);
int oqp_gpu_metc_get_timings(double *secs_out, long *count_out, int n);
void oqp_gpu_metc_reset_timings(void);
int oqp_gpu_ws_acquire(int target, int64_t total_bytes);
int oqp_gpu_ws_validate(int handle, int64_t offset, int64_t nbytes);
void *oqp_gpu_ws_ptr(int handle, int64_t offset);
int64_t oqp_gpu_ws_total_bytes(int handle);
int oqp_gpu_ws_release(int handle);
void oqp_gpu_metc_layout(int nf, int nmatrix, int nbf, int nthreads,
                         int max_ncur, int64_t *off_d3, int64_t *off_f3,
                         int64_t *off_ids, int64_t *off_ints,
                         int64_t *total_bytes);
int oqp_gpu_metc_session_begin(int nf, int nmatrix, int nbf, int nthreads,
                               int max_ncur);
int64_t oqp_gpu_metc_session_total_bytes(int session);
int oqp_gpu_metc_session_capacity(int session);
void *oqp_gpu_metc_d3_ptr(int session);
void *oqp_gpu_metc_f3_ptr(int session, int thread);
void *oqp_gpu_metc_ids_ptr(int session, int thread);
void *oqp_gpu_metc_ints_ptr(int session, int thread);
int oqp_gpu_metc_session_check_ncur(int session, int ncur);
int oqp_gpu_metc_session_zero_f3(int session);
int oqp_gpu_metc_session_upload_d3(int session, const double *d3_host);
int oqp_gpu_metc_session_upload_batch(int session, int thread, const int *ids,
                                      const double *ints, int ncur);
int oqp_gpu_metc_session_contract(int session, int thread, const int *ids,
                                  const double *ints, int ncur, int cur_pass,
                                  double scale_exchange, double scale_coulomb,
                                  int is_umrsf);
int oqp_gpu_metc_session_download_f3(int session, double *f3_host);
int oqp_gpu_metc_session_end(int session);
#ifdef __cplusplus
}
#endif
#endif
