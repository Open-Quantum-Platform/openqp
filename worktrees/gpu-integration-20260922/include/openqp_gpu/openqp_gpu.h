/* openqp-gpu : GPU-accelerated density-fitting SCF for OpenQP.
 *
 * A self-contained CUDA library (cuBLAS + cuSOLVER only; no rotation-method or
 * other integral-engine dependency).  The density-fitting tensor B, the
 * one-electron matrices H and S, and the initial guess density are supplied by
 * the caller; everything from there — J/K builds, Fock assembly, diagonalization,
 * DIIS — runs device-resident on one GPU.
 *
 * All matrices are packed lower-triangular in OpenQP order  t = i*(i+1)/2 + j
 * (i>=j), unless noted.  Integer arguments are passed by pointer for Fortran
 * ISO_C_BINDING compatibility.  info returns 0 on success, non-zero on error.
 *
 * Environment (read by the library):
 *   OPENQP_GPU_B     path to the density-fitting tensor B (required)
 *   CUDA_VISIBLE_DEVICES  selects the GPU
 * Behaviour knobs default to the tuned production values; see README.
 */
#ifndef OPENQP_GPU_H
#define OPENQP_GPU_H

#ifdef __cplusplus
extern "C" {
#endif

/* ---- whole-SCF drivers -------------------------------------------------- */

/* Restricted (closed-shell) density-fitting HF SCF.
 *   h,s      : packed H (core Hamiltonian) and S (overlap), length nbf*(nbf+1)/2
 *   nbf,nocc : basis size, number of doubly-occupied orbitals
 *   enuc     : nuclear repulsion energy
 *   sc,se    : Coulomb and exchange scale (1.0,1.0 = HF; hybrid-DFT sets these)
 *   conv_e,conv_d : energy and density/gradient convergence thresholds
 *   maxit    : iteration cap
 *   e_out    : converged total energy (scalar)
 *   d_out,c_out,eps_out : packed density, MO coefficients, orbital energies
 *   ncyc,info: cycles taken; status
 * The initial guess is read from OQP_SCF_GUESS_D (packed density file) if set,
 * else an on-device generalized-Wolfsberg-Helmholtz guess, else core.
 */
void routec_scf_solve(const double* h, const double* s, const int* nbf,
                      const int* nocc, const double* enuc,
                      const double* sc, const double* se,
                      const double* conv_e, const double* conv_d,
                      const int* maxit,
                      double* e_out, double* d_out,
                      double* c_out, double* eps_out,
                      int* ncyc, int* info);

/* Unrestricted (open-shell) density-fitting HF SCF.  da0,db0 are optional
 * packed alpha/beta guess densities (pass NULL for the built-in guess). */
void routec_scf_solve_uhf(const double* h, const double* s, const int* nbf,
                          const int* nalpha, const int* nbeta,
                          const double* enuc,
                          const double* sc, const double* se,
                          const double* conv_e, const double* conv_d,
                          const int* maxit,
                          const double* da0, const double* db0,
                          double* e_out, double* da_out, double* db_out,
                          double* ca_out, double* cb_out,
                          double* epsa_out, double* epsb_out,
                          int* ncyc, int* info);

/* ---- per-iteration J/K seam (for driving from an external SCF loop) ------ */

/* Build the closed/open-shell Fock matrix contribution J/K from density d.
 * nfocks = 1 (RHF/RKS) or 2 (UHF/UKS, alpha|beta stacked). */
void routec_fock_jk(const double* d, double* f, const int* nbf,
                    const int* nfocks, const double* scale_exch,
                    const double* scale_coul, int* info);

/* Persistent-tensor J/K apply (init once, apply many, free once): for building
 * J and K of a set of trial vectors without re-uploading B each call. */
int  routec_jkm_init(const int* nbf);
void routec_jkm_free(void);
void routec_jkm_apply(const double* x, const int* nbf, const int* nvec,
                      double* jm, double* km, int* info);
/* Low-rank variant (rank-compressed trial densities; e.g. response). */
void routec_jkm_apply_lr(const double* u, const double* v, const int* nbf,
                         const int* nout, const int* ranks,
                         double* jm, double* km, int* info);

/* ---- exchange-correlation (DFT) : provided by the XC module ------------- */
/* Adds the XC potential matrix fxc for density d on the supplied grid.
 * (Present in src/xc.cu; wired into the RKS/UKS driver — see README status.) */
void routec_vxc(const double* d, double* fxc, const int* nbf,
                const int* functional, double* exc, int* info);

/* ---- nuclear gradient (2-electron part) : provided by the gradient module */
/* de += the density-fitted two-electron nuclear gradient for density d at
 * geometry xyz (natm atoms).  hfscale/coulscale set the exact-exchange and
 * Coulomb weights (HF = 1,1; hybrid DFT scales exchange).  The one-electron
 * and (for DFT) exchange-correlation gradient terms are added by the caller.
 * (src/grad.cu; end-to-end validation pending — see README status.) */
void routec_grad2(const double* d, const double* xyz, double* de,
                  const int* nbf, const int* natm, const double* hfscale,
                  const double* coulscale, int* info);

#ifdef __cplusplus
}
#endif

#endif /* OPENQP_GPU_H */
