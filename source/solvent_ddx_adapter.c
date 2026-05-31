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


int oqp_ddx_run_explicit_pcm_reaction_field_smoke(
    oqp_ddx_smoke_result_t* result,
    double* cavity_xyz_out,
    double* q_cav_out,
    int max_cav,
    int* ncav_written,
    char* message,
    int message_len) {
  if (result != NULL) {
    memset(result, 0, sizeof(*result));
  }
  if (ncav_written != NULL) {
    *ncav_written = 0;
  }

#ifndef OQP_ENABLE_DDX
  set_message(message, message_len, "OpenQP was built without OQP_ENABLE_DDX");
  return 2;
#else
  if (cavity_xyz_out == NULL || q_cav_out == NULL || max_cav <= 0) {
    set_message(message, message_len, "Invalid reaction-field output buffers");
    return 1;
  }

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
  void* state_plus = NULL;
  void* state_minus = NULL;
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
  if (max_cav < ncav) {
    set_message(message, message_len, "Reaction-field output buffers are too small");
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
    set_message(message, message_len, "Allocation failure in reaction-field ddX PCM smoke");
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
  double q_cav_norm = 0.0;
  for (int i = 0; i < nbasis * nsph; ++i) {
    x_norm += x[i] * x[i];
    s_norm += s[i] * s[i];
  }
  for (int icav = 0; icav < ncav; ++icav) {
    q_cav_norm += xi[icav] * xi[icav];
    q_cav_out[icav] = xi[icav];
    cavity_xyz_out[3 * icav + 0] = cavity[3 * icav + 0];
    cavity_xyz_out[3 * icav + 1] = cavity[3 * icav + 1];
    cavity_xyz_out[3 * icav + 2] = cavity[3 * icav + 2];
  }
  x_norm = sqrt(x_norm);
  s_norm = sqrt(s_norm);
  q_cav_norm = sqrt(q_cav_norm);

  if (x_norm <= 0.0 || s_norm <= 0.0 || q_cav_norm <= 0.0) {
    set_message(message, message_len, "Unexpected reaction-field ddX PCM solution norms");
    goto cleanup;
  }

  const double finite_difference_delta = 1.0e-6;
  state_plus = ddx_allocate_state(model, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  phi_cav[0] += finite_difference_delta;
  ddx_pcm_setup(model, state_plus, ncav, nbasis, nsph, psi, phi_cav, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_guess(model, state_plus, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_solve(model, state_plus, tol, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  const double energy_plus = ddx_pcm_energy(model, state_plus, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  state_minus = ddx_allocate_state(model, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  phi_cav[0] -= 2.0 * finite_difference_delta;
  ddx_pcm_setup(model, state_minus, ncav, nbasis, nsph, psi, phi_cav, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_guess(model, state_minus, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  ddx_pcm_solve(model, state_minus, tol, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  const double energy_minus = ddx_pcm_energy(model, state_minus, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  phi_cav[0] += finite_difference_delta;

  const double q_cav_fd_derivative =
      (energy_plus - energy_minus) / (2.0 * finite_difference_delta);
  const double q_cav_fd_direct_abs_error = fabs(q_cav_fd_derivative - q_cav_out[0]);
  const double q_cav_fd_abs_error = fabs(q_cav_fd_derivative + 0.5 * q_cav_out[0]);

  if (result != NULL) {
    result->energy = energy;
    result->x_norm = x_norm;
    result->s_norm = s_norm;
    result->xi_norm = q_cav_norm;
    result->q_cav_norm = q_cav_norm;
    result->q_cav_fd_derivative = q_cav_fd_derivative;
    result->q_cav_fd_direct_abs_error = q_cav_fd_direct_abs_error;
    result->q_cav_fd_abs_error = q_cav_fd_abs_error;
    result->first_cavity_value = cavity[0];
    result->nbasis = nbasis;
    result->ncav = ncav;
  }
  if (ncav_written != NULL) {
    *ncav_written = ncav;
  }
  set_message(message, message_len, "explicit ddX PCM reaction-field smoke test passed");
  status = 0;

cleanup:
  if (state_minus != NULL && error != NULL) ddx_deallocate_state(state_minus, error);
  if (state_plus != NULL && error != NULL) ddx_deallocate_state(state_plus, error);
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

#ifdef OQP_ENABLE_DDX
/*
 * Bondi-style van der Waals radii (Angstrom) for the common main-group
 * elements, indexed by nuclear charge Z. Used only to construct the ddPCM
 * cavity. These are provisional defaults, not an OpenQP-curated radius set, and
 * carry the same "unvalidated convention" caveat as q_cav.
 */
static double vdw_radius_angstrom(int z) {
  static const double table[] = {
      1.50,                               /* 0: dummy */
      1.20, 1.40,                         /* H  He */
      1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, /* Li..Ne */
      2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88  /* Na..Ar */
  };
  const int n = (int)(sizeof(table) / sizeof(table[0]));
  if (z < 0) z = 0;
  if (z >= n) return 1.80; /* provisional fallback for heavier elements */
  return table[z];
}

/*
 * Allocate a ddPCM model for the given solute geometry. Centres are supplied in
 * Bohr; per-atom cavity radii are derived from the integer nuclear charges.
 * Returns the model (caller deallocates) or NULL on error with message set.
 * The model is a deterministic function of (geometry, radii, discretization),
 * so oqp_ddx_pcm_cavity and oqp_ddx_pcm_solve agree on ncav and point order.
 */
static void* build_pcm_model(int natom, const double* xyz_bohr,
                             const double* charges, double epsilon, void* error,
                             char* message, int message_len) {
  const double bohr_per_angstrom = 1.0 / 0.5291772109;
  const int model_pcm = 2;
  const int enable_forces = 1;
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

  double* centres = (double*)calloc((size_t)3 * natom, sizeof(double));
  double* radii = (double*)calloc((size_t)natom, sizeof(double));
  if (!centres || !radii) {
    set_message(message, message_len, "Allocation failure building ddX cavity");
    free(centres);
    free(radii);
    return NULL;
  }
  for (int i = 0; i < natom; ++i) {
    centres[3 * i + 0] = xyz_bohr[3 * i + 0];
    centres[3 * i + 1] = xyz_bohr[3 * i + 1];
    centres[3 * i + 2] = xyz_bohr[3 * i + 2];
    const int z = (int)(charges[i] + 0.5);
    radii[i] = vdw_radius_angstrom(z) * bohr_per_angstrom;
  }

  void* model = ddx_allocate_model(
      model_pcm, enable_forces, epsilon, kappa, eta, shift, lmax, n_lebedev,
      incore, maxiter, jacobi_n_diis, enable_fmm, fmm_multipole_lmax,
      fmm_local_lmax, n_proc, natom, centres, radii, length_logfile, logfile,
      error);
  free(centres);
  free(radii);
  if (ddx_get_error_flag(error) != 0) {
    fail_with_ddx_error(error, message, message_len);
    if (model != NULL) ddx_deallocate_model(model, error);
    return NULL;
  }
  return model;
}
#endif

int oqp_ddx_pcm_cavity(int natom, const double* xyz_bohr, const double* charges,
                       double epsilon, int max_cav, int* ncav_out,
                       double* cav_xyz_out, char* message, int message_len) {
  if (ncav_out != NULL) {
    *ncav_out = 0;
  }
#ifndef OQP_ENABLE_DDX
  (void)natom;
  (void)xyz_bohr;
  (void)charges;
  (void)epsilon;
  (void)max_cav;
  (void)cav_xyz_out;
  set_message(message, message_len, "OpenQP was built without OQP_ENABLE_DDX");
  return 2;
#else
  if (natom <= 0 || xyz_bohr == NULL || charges == NULL || ncav_out == NULL ||
      cav_xyz_out == NULL || max_cav <= 0) {
    set_message(message, message_len, "Invalid arguments to oqp_ddx_pcm_cavity");
    return 1;
  }

  int status = 1;
  void* error = NULL;
  void* model = NULL;
  double* cavity = NULL;

  error = ddx_allocate_error();
  if (error == NULL) {
    set_message(message, message_len, "Failed to allocate ddX error object");
    return 1;
  }

  model = build_pcm_model(natom, xyz_bohr, charges, epsilon, error, message,
                          message_len);
  if (model == NULL) goto cleanup;

  const int ncav = ddx_get_n_cav(model);
  if (ncav <= 0) {
    set_message(message, message_len, "Invalid ddX cavity size");
    goto cleanup;
  }
  if (ncav > max_cav) {
    set_message(message, message_len,
                "ddX cavity larger than caller-provided buffer");
    goto cleanup;
  }

  cavity = (double*)calloc((size_t)3 * ncav, sizeof(double));
  if (cavity == NULL) {
    set_message(message, message_len, "Allocation failure in oqp_ddx_pcm_cavity");
    goto cleanup;
  }
  ddx_get_cavity(model, ncav, cavity);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  for (int i = 0; i < 3 * ncav; ++i) {
    cav_xyz_out[i] = cavity[i];
  }
  *ncav_out = ncav;
  set_message(message, message_len, "ddX PCM cavity built");
  status = 0;

cleanup:
  if (model != NULL && error != NULL) ddx_deallocate_model(model, error);
  free(cavity);
  return status;
#endif
}

int oqp_ddx_pcm_solve(int natom, const double* xyz_bohr, const double* charges,
                      double epsilon, int ncav, const double* phi_cav,
                      double* q_cav_out, double* esolv_out, char* message,
                      int message_len) {
  if (esolv_out != NULL) {
    *esolv_out = 0.0;
  }
#ifndef OQP_ENABLE_DDX
  (void)natom;
  (void)xyz_bohr;
  (void)charges;
  (void)epsilon;
  (void)ncav;
  (void)phi_cav;
  (void)q_cav_out;
  set_message(message, message_len, "OpenQP was built without OQP_ENABLE_DDX");
  return 2;
#else
  if (natom <= 0 || xyz_bohr == NULL || charges == NULL || ncav <= 0 ||
      phi_cav == NULL || q_cav_out == NULL) {
    set_message(message, message_len, "Invalid arguments to oqp_ddx_pcm_solve");
    return 1;
  }

  int status = 1;
  void* error = NULL;
  void* model = NULL;
  void* state = NULL;
  double* solute_multipoles = NULL;
  double* psi = NULL;
  double* phi_cav_copy = NULL;
  double* xi = NULL;

  error = ddx_allocate_error();
  if (error == NULL) {
    set_message(message, message_len, "Failed to allocate ddX error object");
    return 1;
  }

  model = build_pcm_model(natom, xyz_bohr, charges, epsilon, error, message,
                          message_len);
  if (model == NULL) goto cleanup;

  const int nbasis = ddx_get_n_basis(model);
  const int model_ncav = ddx_get_n_cav(model);
  if (nbasis <= 0 || model_ncav <= 0) {
    set_message(message, message_len, "Invalid ddX dimensions in pcm_solve");
    goto cleanup;
  }
  if (model_ncav != ncav) {
    set_message(message, message_len,
                "ddX cavity size disagrees with caller phi_cav length");
    goto cleanup;
  }

  const int nsph = natom;
  const int nmultipoles = 1;
  solute_multipoles =
      (double*)calloc((size_t)nmultipoles * nsph, sizeof(double));
  psi = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  phi_cav_copy = (double*)calloc((size_t)ncav, sizeof(double));
  xi = (double*)calloc((size_t)ncav, sizeof(double));
  if (!solute_multipoles || !psi || !phi_cav_copy || !xi) {
    set_message(message, message_len, "Allocation failure in oqp_ddx_pcm_solve");
    goto cleanup;
  }

  /*
   * psi is the solute source in spherical harmonics. For this first energy
   * gate it is built from nuclear monopoles only (charges as Z), matching the
   * branch's explicit-PCM smoke. Folding the electronic density into psi is a
   * deferred refinement; see docs/solvent_ddx_scf_integration_seam.md.
   */
  for (int i = 0; i < nsph; ++i) {
    solute_multipoles[i] = charges[i] / sqrt(4.0 * M_PI);
  }
  ddx_multipole_psi(model, nbasis, nsph, nmultipoles, solute_multipoles, psi,
                    error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  for (int i = 0; i < ncav; ++i) {
    phi_cav_copy[i] = phi_cav[i];
  }

  state = ddx_allocate_state(model, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

  const double tol = 1.0e-9;
  ddx_pcm_setup(model, state, ncav, nbasis, nsph, psi, phi_cav_copy, error);
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

  /*
   * q_cav is the cavity-projected ddPCM adjoint charge (state%q via
   * ddx_get_xi). Its sign/scale relative to OpenQP's external_charge_potential
   * seam is PROVISIONAL: the branch only finite-difference-validated
   * dE/dphi_cav ~= -0.5*q_cav for the energy, not the Fock-operator charge.
   */
  ddx_get_xi(state, model, ncav, xi);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  for (int i = 0; i < ncav; ++i) {
    q_cav_out[i] = xi[i];
  }
  if (esolv_out != NULL) {
    *esolv_out = energy;
  }
  set_message(message, message_len, "ddX PCM solve complete");
  status = 0;

cleanup:
  if (state != NULL && error != NULL) ddx_deallocate_state(state, error);
  if (model != NULL && error != NULL) ddx_deallocate_model(model, error);
  free(solute_multipoles);
  free(psi);
  free(phi_cav_copy);
  free(xi);
  return status;
#endif
}

int oqp_ddx_pcm_solve_multipole_source(
    int natom, const double* xyz_bohr, const double* cavity_charges,
    int nmultipoles, const double* source_multipoles, double epsilon, int ncav,
    double* phi_source_out, double* q_cav_out, double* esolv_out, char* message,
    int message_len) {
  if (esolv_out != NULL) {
    *esolv_out = 0.0;
  }
#ifndef OQP_ENABLE_DDX
  (void)natom;
  (void)xyz_bohr;
  (void)cavity_charges;
  (void)nmultipoles;
  (void)source_multipoles;
  (void)epsilon;
  (void)ncav;
  (void)phi_source_out;
  (void)q_cav_out;
  set_message(message, message_len, "OpenQP was built without OQP_ENABLE_DDX");
  return 2;
#else
  if (natom <= 0 || xyz_bohr == NULL || cavity_charges == NULL ||
      nmultipoles <= 0 || source_multipoles == NULL || ncav <= 0 ||
      phi_source_out == NULL || q_cav_out == NULL) {
    set_message(message, message_len,
                "Invalid arguments to oqp_ddx_pcm_solve_multipole_source");
    return 1;
  }
  int mmax = 0;
  while ((mmax + 1) * (mmax + 1) < nmultipoles) ++mmax;
  if ((mmax + 1) * (mmax + 1) != nmultipoles) {
    set_message(message, message_len,
                "nmultipoles must be a perfect square for ddX multipoles");
    return 1;
  }

  int status = 1;
  void* error = NULL;
  void* model = NULL;
  void* electrostatics = NULL;
  void* state = NULL;
  double* solute_multipoles = NULL;
  double* psi = NULL;
  double* forces = NULL;
  double* xi = NULL;

  error = ddx_allocate_error();
  if (error == NULL) {
    set_message(message, message_len, "Failed to allocate ddX error object");
    return 1;
  }

  model = build_pcm_model(natom, xyz_bohr, cavity_charges, epsilon, error,
                          message, message_len);
  if (model == NULL) goto cleanup;

  const int nbasis = ddx_get_n_basis(model);
  const int model_ncav = ddx_get_n_cav(model);
  if (nbasis <= 0 || model_ncav <= 0) {
    set_message(message, message_len,
                "Invalid ddX dimensions in multipole-source pcm solve");
    goto cleanup;
  }
  if (model_ncav != ncav) {
    set_message(message, message_len,
                "ddX cavity size disagrees with caller q_cav length");
    goto cleanup;
  }

  const int nsph = natom;
  solute_multipoles =
      (double*)calloc((size_t)nmultipoles * nsph, sizeof(double));
  psi = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  forces = (double*)calloc((size_t)3 * nsph, sizeof(double));
  xi = (double*)calloc((size_t)ncav, sizeof(double));
  if (!solute_multipoles || !psi || !forces || !xi) {
    set_message(message, message_len,
                "Allocation failure in oqp_ddx_pcm_solve_multipole_source");
    goto cleanup;
  }

  memcpy(solute_multipoles, source_multipoles,
         (size_t)nmultipoles * (size_t)nsph * sizeof(double));

  ddx_multipole_electrostatics_0(model, nsph, ncav, nmultipoles,
                                 solute_multipoles, phi_source_out, error);
  if (check_ddx_error(error, message, message_len)) goto cleanup;

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

  ddx_get_xi(state, model, ncav, xi);
  if (check_ddx_error(error, message, message_len)) goto cleanup;
  for (int i = 0; i < ncav; ++i) {
    q_cav_out[i] = xi[i];
  }
  if (esolv_out != NULL) {
    *esolv_out = energy;
  }
  set_message(message, message_len, "ddX PCM multipole-source solve complete");
  status = 0;

cleanup:
  if (state != NULL && error != NULL) ddx_deallocate_state(state, error);
  if (electrostatics != NULL && error != NULL) ddx_deallocate_electrostatics(electrostatics, error);
  if (model != NULL && error != NULL) ddx_deallocate_model(model, error);
  free(solute_multipoles);
  free(psi);
  free(forces);
  free(xi);
  return status;
#endif
}

