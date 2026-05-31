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

/*
 * Production PCM reaction-field seam for the SCF loop (real ddX only).
 *
 * Two-phase, stateless handshake so the OpenQP Fortran SCF can supply phi_cav
 * built from the AO density:
 *
 *   1. oqp_ddx_pcm_cavity(): build the ddPCM cavity from solute geometry and
 *      return the cavity-point Cartesian coordinates (Bohr, 3*ncav) and count.
 *   2. The caller evaluates phi_cav (electronic + nuclear) at those points.
 *   3. oqp_ddx_pcm_solve(): rebuild the identical model, run the ddPCM
 *      forward/adjoint solve for the supplied phi_cav, and return the
 *      cavity-projected adjoint charge q_cav (ddx_get_xi) plus the ddX
 *      solvation energy in esolv.
 *
 * The cavity is a deterministic function of (geometry, radii, discretization),
 * so the two calls agree on ncav and cavity-point ordering. Both return 0 on
 * success, 2 when OpenQP was built without OQP_ENABLE_DDX, and 1 on a
 * ddX/runtime error (with a human-readable message). Atomic cavity radii are
 * derived from the nuclear charges via a small built-in van der Waals table.
 *
 * NOTE: the sign/scale convention of q_cav relative to OpenQP's
 * external_charge_potential seam is PROVISIONAL and not yet validated against a
 * trusted PCM reference. See docs/solvent_ddx_scf_integration_seam.md.
 */
int oqp_ddx_pcm_cavity(int natom, const double* xyz_bohr, const double* charges,
                       double epsilon, int max_cav, int* ncav_out,
                       double* cav_xyz_out, char* message, int message_len);

int oqp_ddx_pcm_solve(int natom, const double* xyz_bohr, const double* charges,
                      double epsilon, int ncav, const double* phi_cav,
                      double* q_cav_out, double* esolv_out, char* message,
                      int message_len);

#ifdef __cplusplus
}
#endif

#endif /* OQP_SOLVENT_DDX_ADAPTER_H */
