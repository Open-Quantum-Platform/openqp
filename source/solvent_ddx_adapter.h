#ifndef OQP_SOLVENT_DDX_ADAPTER_H
#define OQP_SOLVENT_DDX_ADAPTER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct oqp_ddx_smoke_result_s {
  double energy;
  double x_norm;
  double s_norm;
  double xi_norm;
  double first_cavity_value;
  int nbasis;
  int ncav;
} oqp_ddx_smoke_result_t;

/*
 * Run a minimal ddX host-code lifecycle using a small point-charge model.
 *
 * This is intentionally not OpenQP SCF coupling yet. It verifies that the
 * OpenQP-owned ddX adapter can allocate a model, build a source term, solve,
 * retrieve forward/adjoint/cavity outputs, validate basic numerical sanity,
 * and clean up. Returns 0 on success and nonzero on failure.
 */
int oqp_ddx_run_point_charge_smoke(oqp_ddx_smoke_result_t* result,
                                   char* message,
                                   int message_len);

#ifdef __cplusplus
}
#endif

#endif /* OQP_SOLVENT_DDX_ADAPTER_H */
