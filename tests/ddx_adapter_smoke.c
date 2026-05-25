#include <stdio.h>

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
  printf("energy=%.16f nbasis=%d ncav=%d x_norm=%.8e s_norm=%.8e xi_norm=%.8e\n",
         result.energy, result.nbasis, result.ncav, result.x_norm,
         result.s_norm, result.xi_norm);
  return 0;
}
