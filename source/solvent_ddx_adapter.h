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
  double q_cav_norm;
  double q_cav_fd_derivative;
  double q_cav_fd_direct_abs_error;
  double q_cav_fd_abs_error;
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

/*
 * Run the same point-charge model through the explicit host-code PCM path:
 * build host-supplied psi and phi_cav, call ddx_pcm_setup/solve/solve_adjoint,
 * and retrieve ddX's cavity-projected ddPCM state%q through ddx_get_xi.
 * This is the next seam needed before OpenQP can provide psi/phi_cav from
 * AO densities in the SCF loop.
 */
int oqp_ddx_run_explicit_pcm_smoke(oqp_ddx_smoke_result_t* result,
                                   char* message,
                                   int message_len);

/*
 * Explicit PCM smoke that also copies the candidate reaction-field handoff
 * arrays for OpenQP's external_charge_potential seam: cavity_xyz is a
 * 3*ncav Bohr coordinate array and q_cav is the cavity projection of ddPCM
 * state%q returned by ddx_get_xi. Caller owns the output buffers.
 */
int oqp_ddx_run_explicit_pcm_reaction_field_smoke(
    oqp_ddx_smoke_result_t* result,
    double* cavity_xyz,
    double* q_cav,
    int max_cav,
    int* ncav_written,
    char* message,
    int message_len);

#ifdef __cplusplus
}
#endif

#endif /* OQP_SOLVENT_DDX_ADAPTER_H */
