#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ddx.h"

static void print_error(void* error) {
  char message[2000];
  ddx_get_error_message(error, message, 2000);
  fprintf(stderr, "ddX error: %s\n", message);
}

static int check_error(void* error) {
  if (ddx_get_error_flag(error) != 0) {
    print_error(error);
    return 1;
  }
  return 0;
}

int main(void) {
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

  double radii[12];
  double centres[36];
  for (int i = 0; i < nsph; ++i) {
    radii[i] = radii_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 0] = x_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 1] = y_angstrom[i] * bohr_per_angstrom;
    centres[3 * i + 2] = z_angstrom[i] * bohr_per_angstrom;
  }

  void* error = ddx_allocate_error();
  if (error == NULL) {
    fprintf(stderr, "Failed to allocate ddX error object\n");
    return 1;
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

  void* model = ddx_allocate_model(
      model_pcm, enable_forces, epsilon, kappa, eta, shift, lmax, n_lebedev,
      incore, maxiter, jacobi_n_diis, enable_fmm, fmm_multipole_lmax,
      fmm_local_lmax, n_proc, nsph, centres, radii, length_logfile, logfile,
      error);
  if (check_error(error)) return 1;

  const int nbasis = ddx_get_n_basis(model);
  const int ncav = ddx_get_n_cav(model);
  if (nbasis <= 0 || ncav <= 0) {
    fprintf(stderr, "Invalid ddX dimensions nbasis=%d ncav=%d\n", nbasis, ncav);
    return 1;
  }

  const int nmultipoles = 1;
  double* solute_multipoles = (double*)calloc((size_t)nmultipoles * nsph, sizeof(double));
  double* psi = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  double* forces = (double*)calloc((size_t)3 * nsph, sizeof(double));
  double* x = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  double* s = (double*)calloc((size_t)nbasis * nsph, sizeof(double));
  double* xi = (double*)calloc((size_t)ncav, sizeof(double));
  double* cavity = (double*)calloc((size_t)ncav, sizeof(double));
  if (!solute_multipoles || !psi || !forces || !x || !s || !xi || !cavity) {
    fprintf(stderr, "Allocation failure in ddX adapter smoke test\n");
    return 1;
  }

  for (int i = 0; i < nsph; ++i) {
    solute_multipoles[i] = charges[i] / sqrt(4.0 * M_PI);
  }

  void* electrostatics = ddx_allocate_electrostatics(model, error);
  if (check_error(error)) return 1;

  ddx_multipole_electrostatics(model, nsph, nmultipoles, solute_multipoles,
                               electrostatics, error);
  if (check_error(error)) return 1;

  ddx_multipole_psi(model, nbasis, nsph, nmultipoles, solute_multipoles, psi,
                    error);
  if (check_error(error)) return 1;

  void* state = ddx_allocate_state(model, error);
  if (check_error(error)) return 1;

  const double tol = 1.0e-9;
  const int read_guess = 0;
  const double energy = ddx_ddrun(model, state, electrostatics, nbasis, nsph,
                                  psi, tol, forces, read_guess, error);
  if (check_error(error)) return 1;

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
    fprintf(stderr, "Unexpected ddX energy %.16f (reference %.16f)\n", energy,
            reference_energy);
    return 1;
  }
  if (x_norm <= 0.0 || s_norm <= 0.0 || xi_norm <= 0.0 || cavity[0] <= 0.0) {
    fprintf(stderr, "Unexpected ddX solution norms x=%g s=%g xi=%g cavity0=%g\n",
            x_norm, s_norm, xi_norm, cavity[0]);
    return 1;
  }

  printf("ddX adapter smoke test passed\n");
  printf("energy=%.16f nbasis=%d ncav=%d x_norm=%.8e s_norm=%.8e xi_norm=%.8e\n",
         energy, nbasis, ncav, x_norm, s_norm, xi_norm);

  ddx_deallocate_state(state, error);
  ddx_deallocate_electrostatics(electrostatics, error);
  ddx_deallocate_model(model, error);
  free(solute_multipoles);
  free(psi);
  free(forces);
  free(x);
  free(s);
  free(xi);
  free(cavity);
  return 0;
}
