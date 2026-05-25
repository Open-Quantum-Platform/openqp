#include <stdio.h>
#include <stdlib.h>

#include "solvent_ddx_adapter.h"

int main(void) {
  oqp_ddx_smoke_result_t result;
  char message[512];

  const int status = oqp_ddx_run_point_charge_smoke(
      &result, message, (int)sizeof(message));
  if (status != 0) {
    fprintf(stderr, "ddX adapter smoke test failed: %s\n", message);
    return status;
  }

  printf("%s\n", message);
  printf("energy=%.16f nbasis=%d ncav=%d x_norm=%.8e s_norm=%.8e xi_norm=%.8e q_cav_norm=%.8e\n",
         result.energy, result.nbasis, result.ncav, result.x_norm,
         result.s_norm, result.xi_norm, result.q_cav_norm);

  const int explicit_status = oqp_ddx_run_explicit_pcm_smoke(
      &result, message, (int)sizeof(message));
  if (explicit_status != 0) {
    fprintf(stderr, "explicit ddX PCM smoke test failed: %s\n", message);
    return explicit_status;
  }

  printf("%s\n", message);
  printf("explicit_energy=%.16f nbasis=%d ncav=%d x_norm=%.8e s_norm=%.8e xi_norm=%.8e q_cav_norm=%.8e\n",
         result.energy, result.nbasis, result.ncav, result.x_norm,
         result.s_norm, result.xi_norm, result.q_cav_norm);

  const int max_cav = 2000;
  double* cavity_xyz = (double*)calloc((size_t)3 * max_cav, sizeof(double));
  double* q_cav = (double*)calloc((size_t)max_cav, sizeof(double));
  int ncav_written = 0;
  if (cavity_xyz == NULL || q_cav == NULL) {
    fprintf(stderr, "failed to allocate reaction-field smoke buffers\n");
    free(cavity_xyz);
    free(q_cav);
    return 1;
  }

  const int reaction_status = oqp_ddx_run_explicit_pcm_reaction_field_smoke(
      &result, cavity_xyz, q_cav, max_cav, &ncav_written, message,
      (int)sizeof(message));
  if (reaction_status != 0) {
    fprintf(stderr, "explicit ddX PCM reaction-field smoke test failed: %s\n", message);
    free(cavity_xyz);
    free(q_cav);
    return reaction_status;
  }

  printf("%s\n", message);
  printf("reaction_field_ncav=%d q_cav_norm=%.8e first_q_cav=%.8e first_cavity_xyz=(%.8e,%.8e,%.8e)\n",
         ncav_written, result.q_cav_norm, q_cav[0], cavity_xyz[0],
         cavity_xyz[1], cavity_xyz[2]);
  free(cavity_xyz);
  free(q_cav);
  return 0;
}
