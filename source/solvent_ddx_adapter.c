#include "solvent_ddx_adapter.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef OQP_ENABLE_DDX
#include "ddx.h"
#endif

static void set_message(char* message, int message_len, const char* text) {
  if (message == NULL || message_len <= 0) {
    return;
  }
  snprintf(message, (size_t)message_len, "%s", text == NULL ? "" : text);
}

#ifdef OQP_ENABLE_DDX
static int fail_with_ddx_error(void* error, char* message, int message_len) {
  char ddx_message[2000];
  ddx_message[0] = '\0';
  ddx_get_error_message(error, ddx_message, (int)sizeof(ddx_message));
  set_message(message, message_len, ddx_message);
  return 1;
}

static int check_ddx_error(void* error, char* message, int message_len) {
  if (ddx_get_error_flag(error) != 0) {
    return fail_with_ddx_error(error, message, message_len);
  }
  return 0;
}
#endif

int oqp_ddx_run_point_charge_smoke(oqp_ddx_smoke_result_t* result,
                                   char* message,
                                   int message_len) {
  if (result != NULL) {
    memset(result, 0, sizeof(*result));
  }

#ifndef OQP_ENABLE_DDX
  set_message(message, message_len, "OpenQP was built without OQP_ENABLE_DDX");
  return 2;
#else
  const int nsph = 12;
  const double bohr_per_angstrom = 1.0 / 0.5291772109;
  const double charges[12] = {
      -0.04192, -0.04192, -0.04198, -0.04192, -0.04192, -0.04198,
       0.04193,  0.04193,  0.04197,  0.04193,  0.04193,  0.04197};
  const double radii_angstrom[12] = {
      4.00253, 4.00253, 4.00253, 4.00253, 4.00253, 4.00253,
      2.99956, 2.99956, 2.99956, 2.99956, 2.99956, 2.99956};
  const double x_angstrom[12] = {
      0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
      0.00103, 0.00103, 0.00000, -0.00103, -0.00103, 0.00000};
  const double y_angstrom[12] = {
      2.29035, 2.29035, 0.00000, -2.29035, -2.29035, 0.00000,
      4.05914, 4.05914, 0.00000, -4.05914, -4.05914, 0.00000};
  const double z_angstrom[12] = {
      1.32281, -1.32281, -2.64562, -1.32281, 1.32281, 2.64562,
      2.34326, -2.34326, -4.68652, -2.34326, 2.34326, 4.68652};

  int status = 1;
  void* error = NULL;
  void* model = NULL;
  void* electrostatics = NULL;
  void* state = NULL;
  double* solute_multipoles = NULL;
  double* psi = NULL;
  double* forces = NULL;
  double* x = NULL;
  double* s = NULL;
  double* xi = NULL;
  double* cavity = NULL;

  double radii[12];
  double centres[36];
  for (int i = 0; i < nsph; ++i) {
    radii[i] = radii_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 0] = x_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 1] = y_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 2] = z_angstrom[i] * bohr_per_angstrom;
  }

  error = ddx_allocate_error();
  if (error == NULL) {
    set_message(message, message_len, "Failed to allocate ddX error object");
    goto cleanup;
  }

  const int model_pcm = 2;
  const int enable_forces = 1;
  const double epsilon = 78.3553;
  const double kappa = 0.0;
  const double eta = 0.1;
  const double shift = 0.0;
  const int lmax = 8;
  const int n_lebedev = 302;
  const int incore = 0;
  const int maxiter = 100;
  const int jacobi_n_diis = 20;
  const int enable_fmm = 1;
  const int fmm_multipole_lmax = 7;
  const int fmm_local_lmax = 6;
  const int n_proc = 1;
  const int length_logfile = 0;
  char logfile[1] = {'\0'};

  model = ddx_allocate_model(
      model_pcm, enable_forces, epsilon, kappa, eta, shift, lmax, n_lebedev,
      incore, maxiter, jacobi_n_diis, enable_fmm, fmm_multipole_lmax,
      fmm_local_lmax, n_proc, nsph, centres, radii, length_logfile, logfile,
      error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  const int nbasis = ddx_get_n_basis(model);
  const int ncav = ddx_get_n_cav(model);
  if (nbasis <= 0 || ncav <= 0) {
    set_message(message, message_len, "Invalid ddX dimensions");
    goto cleanup;
  }

  const int nmultipoles = 1;
  solute_multipoles = (double*)calloc((size_t)nmultipoles * nsph, sizeof(double));
  psi = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  forces = (double*)calloc((size_t)3 * nsph, sizeof(double));
  x = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  s = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  xi = (double*)calloc((size_t)ncav, sizeof(double));
  cavity = (double*)calloc((size_t)3 * ncav, sizeof(double));
  if (!solute_multipoles || !psi || !forces || !x || !s || !xi || !cavity) {
    set_message(message, message_len, "Allocation failure in ddX adapter");
    goto cleanup;
  }

  for (int i = 0; i < nsph; ++i) {
    solute_multipoles[i] = charges[i] / sqrt(4.0 * M_PI);
  }

  electrostatics = ddx_allocate_electrostatics(model, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  ddx_multipole_electrostatics(model, nsph, nmultipoles, solute_multipoles,
                               electrostatics, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  ddx_multipole_psi(model, nbasis, nsph, nmultipoles, solute_multipoles, psi,
                    error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  state = ddx_allocate_state(model, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  const double tol = 1.0e-9;
  const int read_guess = 0;
  const double energy = ddx_ddrun(model, state, electrostatics, nbasis, nsph,
                                  psi, tol, forces, read_guess, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  ddx_get_x(state, nbasis, nsph, x);
  ddx_get_s(state, nbasis, nsph, s);
  ddx_get_xi(state, model, ncav, xi);
  ddx_get_cavity(model, ncav, cavity);

  double x_norm = 0.0;
  double s_norm = 0.0;
  double xi_norm = 0.0;
  for (int i = 0; i < nbasis * nsph; ++i) {
    x_norm += x[i] * x[i];
    s_norm += s[i] * s[i];
  }
  for (int i = 0; i < ncav; ++i) {
    xi_norm += xi[i] * xi[i];
  }
  x_norm = sqrt(x_norm);
  s_norm = sqrt(s_norm);
  xi_norm = sqrt(xi_norm);

  const double reference_energy = -0.00017974013712832552;
  if (fabs(energy - reference_energy) > 1.0e-6) {
    set_message(message, message_len, "Unexpected ddX energy");
    goto cleanup;
  }
  if (x_norm <= 0.0 || s_norm <= 0.0 || xi_norm <= 0.0 || cavity[0] <= 0.0) {
    set_message(message, message_len, "Unexpected ddX solution norms");
    goto cleanup;
  }

  if (result != NULL) {
    result->energy = energy;
    result->x_norm = x_norm;
    result->s_norm = s_norm;
    result->xi_norm = xi_norm;
    result->q_cav_norm = xi_norm;
    result->first_cavity_value = cavity[0];
    result->nbasis = nbasis;
    result->ncav = ncav;
  }
  set_message(message, message_len, "ddX adapter smoke test passed");
  status = 0;

cleanup:
  if (state != NULL && error != NULL) ddx_deallocate_state(state, error);
  if (electrostatics != NULL && error != NULL) ddx_deallocate_electrostatics(electrostatics, error);
  if (model != NULL && error != NULL) ddx_deallocate_model(model, error);
  free(solute_multipoles);
  free(psi);
  free(forces);
  free(x);
  free(s);
  free(xi);
  free(cavity);
  return status;
#endif
}

int oqp_ddx_run_explicit_pcm_smoke(oqp_ddx_smoke_result_t* result,
                                   char* message,
                                   int message_len) {
  if (result != NULL) {
    memset(result, 0, sizeof(*result));
  }

#ifndef OQP_ENABLE_DDX
  set_message(message, message_len, "OpenQP was built without OQP_ENABLE_DDX");
  return 2;
#else
  const int nsph = 12;
  const double bohr_per_angstrom = 1.0 / 0.5291772109;
  const double charges[12] = {
      -0.04192, -0.04192, -0.04198, -0.04192, -0.04192, -0.04198,
       0.04193,  0.04193,  0.04197,  0.04193,  0.04193,  0.04197};
  const double radii_angstrom[12] = {
      4.00253, 4.00253, 4.00253, 4.00253, 4.00253, 4.00253,
      2.99956, 2.99956, 2.99956, 2.99956, 2.99956, 2.99956};
  const double x_angstrom[12] = {
      0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
      0.00103, 0.00103, 0.00000, -0.00103, -0.00103, 0.00000};
  const double y_angstrom[12] = {
      2.29035, 2.29035, 0.00000, -2.29035, -2.29035, 0.00000,
      4.05914, 4.05914, 0.00000, -4.05914, -4.05914, 0.00000};
  const double z_angstrom[12] = {
      1.32281, -1.32281, -2.64562, -1.32281, 1.32281, 2.64562,
      2.34326, -2.34326, -4.68652, -2.34326, 2.34326, 4.68652};

  int status = 1;
  void* error = NULL;
  void* model = NULL;
  void* state = NULL;
  double* solute_multipoles = NULL;
  double* psi = NULL;
  double* phi_cav = NULL;
  double* x = NULL;
  double* s = NULL;
  double* xi = NULL;
  double* cavity = NULL;

  double radii[12];
  double centres[36];
  for (int i = 0; i < nsph; ++i) {
    radii[i] = radii_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 0] = x_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 1] = y_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 2] = z_angstrom[i] * bohr_per_angstrom;
  }

  error = ddx_allocate_error();
  if (error == NULL) {
    set_message(message, message_len, "Failed to allocate ddX error object");
    goto cleanup;
  }

  const int model_pcm = 2;
  const int enable_forces = 1;
  const double epsilon = 78.3553;
  const double kappa = 0.0;
  const double eta = 0.1;
  const double shift = 0.0;
  const int lmax = 8;
  const int n_lebedev = 302;
  const int incore = 0;
  const int maxiter = 100;
  const int jacobi_n_diis = 20;
  const int enable_fmm = 1;
  const int fmm_multipole_lmax = 7;
  const int fmm_local_lmax = 6;
  const int n_proc = 1;
  const int length_logfile = 0;
  char logfile[1] = {'\0'};

  model = ddx_allocate_model(
      model_pcm, enable_forces, epsilon, kappa, eta, shift, lmax, n_lebedev,
      incore, maxiter, jacobi_n_diis, enable_fmm, fmm_multipole_lmax,
      fmm_local_lmax, n_proc, nsph, centres, radii, length_logfile, logfile,
      error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  const int nbasis = ddx_get_n_basis(model);
  const int ncav = ddx_get_n_cav(model);
  if (nbasis <= 0 || ncav <= 0) {
    set_message(message, message_len, "Invalid ddX dimensions");
    goto cleanup;
  }

  const int nmultipoles = 1;
  solute_multipoles = (double*)calloc((size_t)nmultipoles * nsph, sizeof(double));
  psi = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  phi_cav = (double*)calloc((size_t)ncav, sizeof(double));
  x = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  s = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  xi = (double*)calloc((size_t)ncav, sizeof(double));
  cavity = (double*)calloc((size_t)3 * ncav, sizeof(double));
  if (!solute_multipoles || !psi || !phi_cav || !x || !s || !xi || !cavity) {
    set_message(message, message_len, "Allocation failure in explicit ddX PCM smoke");
    goto cleanup;
  }

  for (int i = 0; i < nsph; ++i) {
    solute_multipoles[i] = charges[i] / sqrt(4.0 * M_PI);
  }

  ddx_multipole_psi(model, nbasis, nsph, nmultipoles, solute_multipoles, psi,
                    error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  ddx_get_cavity(model, ncav, cavity);
  for (int icav = 0; icav < ncav; ++icav) {
    const double rx = cavity[3 * icav + 0];
    const double ry = cavity[3 * icav + 1];
    const double rz = cavity[3 * icav + 2];
    double potential = 0.0;
    for (int isph = 0; isph < nsph; ++isph) {
      const double dx = rx - centres[3 * isph + 0];
      const double dy = ry - centres[3 * isph + 1];
      const double dz = rz - centres[3 * isph + 2];
      const double distance = sqrt(dx * dx + dy * dy + dz * dz);
      if (distance > 1.0e-14) {
        potential += charges[isph] / distance;
      }
    }
    phi_cav[icav] = potential;
  }

  state = ddx_allocate_state(model, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  const double tol = 1.0e-9;
  ddx_pcm_setup(model, state, ncav, nbasis, nsph, psi, phi_cav, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_guess(model, state, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_solve(model, state, tol, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_guess_adjoint(model, state, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_solve_adjoint(model, state, tol, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  const double energy = ddx_pcm_energy(model, state, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  ddx_get_x(state, nbasis, nsph, x);
  ddx_get_s(state, nbasis, nsph, s);
  ddx_get_xi(state, model, ncav, xi);

  double x_norm = 0.0;
  double s_norm = 0.0;
  double xi_norm = 0.0;
  for (int i = 0; i < nbasis * nsph; ++i) {
    x_norm += x[i] * x[i];
    s_norm += s[i] * s[i];
  }
  for (int i = 0; i < ncav; ++i) {
    xi_norm += xi[i] * xi[i];
  }
  x_norm = sqrt(x_norm);
  s_norm = sqrt(s_norm);
  xi_norm = sqrt(xi_norm);

  if (x_norm <= 0.0 || s_norm <= 0.0 || xi_norm <= 0.0) {
    set_message(message, message_len, "Unexpected explicit ddX PCM solution norms");
    goto cleanup;
  }

  if (result != NULL) {
    result->energy = energy;
    result->x_norm = x_norm;
    result->s_norm = s_norm;
    result->xi_norm = xi_norm;
    result->q_cav_norm = xi_norm;
    result->first_cavity_value = cavity[0];
    result->nbasis = nbasis;
    result->ncav = ncav;
  }
  set_message(message, message_len, "explicit ddX PCM smoke test passed");
  status = 0;

cleanup:
  if (state != NULL && error != NULL) ddx_deallocate_state(state, error);
  if (model != NULL && error != NULL) ddx_deallocate_model(model, error);
  free(solute_multipoles);
  free(psi);
  free(phi_cav);
  free(x);
  free(s);
  free(xi);
  free(cavity);
  return status;
#endif
}

