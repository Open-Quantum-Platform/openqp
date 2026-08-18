#include <stdint.h>
#include <stdbool.h>

/* Internal deterministic simplex-QP test ABI. H/forbidden are column-major. */
void oqp_simplex_qp_solve(int64_t n, const double *h, const double *g,
                          double *x, double *value, int *status);
void oqp_simplex_qp_solve_avoid(int64_t n, const double *h, const double *g,
                                int64_t nforbidden, const double *forbidden,
                                int64_t forbid_vertices_before, double *x,
                                double *value, int *status);

typedef double xyz_t[3];

typedef struct oqp_handle_t {
    void *inf;
    xyz_t * xyz;
    double *qn;
    double *mass;
    xyz_t * grad;
    struct molecule *mol_prop;
    struct energy_results *mol_energy;
    struct dft_parameters *dft;
    struct tddft_parameters *tddft;
    struct control_parameters *control;
    struct mpi_communicator *mpiinfo;
    struct electron_shell *elshell;
    struct trah_control *trah;
} oqp_handle_t;

struct Cstring{
  int64_t length;
  char *string;
};

struct molecule {
    int64_t   natom;
    int64_t   charge;
    int64_t   nelec;
    int64_t   nelec_A;
    int64_t   nelec_B;
    int64_t   mult;
    int64_t   nvelec;
    int64_t   nocc;
};

struct energy_results {
    double   energy;
    double   enuc;
    double   psinrm;
    double   ehf1;
    double   vee;
    double   nenergy;
    double   etot;
    double   vne;
    double   vnn;
    double   vtot;
    double   tkin;
    double   virial;
    double   excited_energy;
    bool     SCF_converged;
    bool     Davidson_converged;
    bool     Z_Vector_converged;
};

struct dft_parameters {
    char XC_functional_name[20];
    double hfscale;
    double cam_alpha;
    double cam_beta;
    double cam_mu;
    double MP2SS_Scale;
    double MP2OS_Scale;
    bool cam_flag;
    bool dh_flag;
    bool grid_pruned;
    bool grid_ao_pruned;
    double grid_ao_threshold;
    double grid_ao_sparsity_ratio;
    char grid_pruned_name[16];
    int64_t grid_num_ang_grids;
    int64_t grid_rad_size;
    int64_t grid_ang_size;
    double  grid_density_cutoff;
    int64_t dft_partfun;
    int64_t rad_grid_type;
    int64_t dft_bfc_algo;
    bool dft_wt_der;
};

struct tddft_parameters {
    int64_t nstate;
    int64_t target_state;
    int64_t maxvec;
    int64_t mult;
    double cnvtol;
    double zvconv;
    bool debug_mode;
    bool tda;
    int64_t tlf;
    double hfscale;
    double cam_alpha;
    double cam_beta;
    double cam_mu;
    double spc_coco;
    double spc_ovov;
    double spc_coov;
    int32_t* ixcore;
    int64_t ixcore_len;
    int64_t z_solver;
    int64_t gmres_dim;
    bool umrsf;
};

struct control_parameters {
    int64_t   hamilton;
    int64_t   scftype;
    char      runtype[20];
    int64_t   guess;
    int64_t   active_basis;
    int64_t   maxit;
    int64_t   maxit_dav;
    int64_t   maxit_zv;
    int64_t   maxdiis;
    int64_t   diis_reset_mod;
    double    diis_reset_conv;
    int64_t   diis_type;
    double    cdiis_switch;
    double    vdiis_vshift_switch;
    double    vshift;
    bool      mom;
    bool      pfon;
    double    mom_switch;
    double    pfon_start_temp;
    double    pfon_cooling_rate;
    double    pfon_nsmear;
    double    conv;
    int64_t   scf_incremental;
    double    int2e_cutoff;
    int64_t   scf_pscreen;
    double    pscreen_k;
    double    pscreen_cap;
    double    pscreen_tight;
    double    pscreen_xc_dcut;
    double    pscreen_xc_aocut;
    int64_t   pscreen_grid_rad;
    int64_t   pscreen_grid_ang;
    int64_t   esp;
    int64_t   resp_target;
    double    esp_constr;
    bool      basis_set_issue;
    double    conf_print_threshold;
    bool      rstctmo;
    int64_t   scal_rel;
    int64_t   soc_2e;
    int64_t   converger_type;
    double    soscf_lvl_shift;
    int64_t   verbose;
    bool      trh_stab;
    bool      trh_ls;
    int64_t   trh_sub_solver;
    int64_t   trh_nrtv;
    double    trh_r0;
    int64_t   trh_jd_start;
    int64_t   trh_nmic;
    double    trh_gred;
    double    trh_lred;
    int64_t   trh_impl;
    bool      sd_scf;
    bool      pcm_enabled;
    double    pcm_epsilon;
    /* Performance knobs -- keep in sync with control_parameters in source/types.F90 */
    int64_t   xc_c2f;
    int64_t   xc_phi_cache;
    int64_t   xc_incdft;
    double    grad_cutoff;
    double    mrsf_resp_cutoff;
    int64_t   mrsf_fp32;
    int64_t   mrsf_zv_warmstart;
    bool      qmmm_flag;
    /* Coupled cluster -- keep in sync with control_parameters in source/types.F90 */
    int64_t   cc_maxit;
    double    cc_conv;
    int64_t   cc_ndiis;
    int64_t   cc_nfzc;
    int64_t   cc_triples;
    int64_t   cc_cholesky;
    double    cc_cholesky_tol;
    int64_t   cc_cholesky_direct;
};

struct mpi_communicator {
        int32_t comm;
        bool debug_mode;
        bool usempi;
};

struct electron_shell {
        int id;
	int element_id;
	int32_t ang_mom;
	int32_t harmonic;
	int32_t ecp_nam;
	int* num_expo;
	double* expo;
	double* coef;
        int* ecp_am;
        int* ecp_rex;
	double* ecp_coord;
	int* ecp_zn;
};

oqp_handle_t *oqp_init();
int oqp_clean(oqp_handle_t * c_handle);
int oqp_have_openmp(void);

/* 1 when the memory budget is a flat OQP_MEMORY_LIMIT_GB that still has to
   cover memory already allocated; 0 when it is the live probe, which reports
   what remains and has therefore already subtracted it. */
int oqp_memory_budget_includes_resident(void);
void oqp_omp_set_num_threads(int n);
int64_t oqp_get(struct oqp_handle_t *c_handle, char *code,
        int32_t *type_id, int32_t *ndims, int64_t *dims, void **v);
int64_t oqp_alloc(struct oqp_handle_t *c_handle, char *code,
        int32_t *type_id, int32_t *ndims, int64_t *dims, void **v);
int64_t oqp_del(struct oqp_handle_t *c_handle, char *code);
int64_t oqp_get_nbf(struct oqp_handle_t *c_handle);
int64_t oqp_get_basis(struct oqp_handle_t *c_handle,
        int64_t *nsh, int64_t *nprim, int64_t *nbf,
        int64_t **bt, int64_t **at, int64_t **cdeg, double **ex, double **cc);

/* Per-shell "stored as spherical harmonics" flag, 1 = spherical.
   This is the EFFECTIVE flag: a shell is spherical only when the run is
   harmonic, the shell is tagged pure, and l >= 2 (s and p stay Cartesian
   even in a spherical basis). Needed to know the AO layout; do not try to
   reconstruct it from the angular momenta alone. */
int64_t oqp_get_basis_spherical(struct oqp_handle_t *c_handle,
        int64_t *nsh, int64_t **spherical);

/* `mass` is optional, pass NULL if not needed */
int oqp_set_atoms(struct oqp_handle_t * c_handle, int64_t natoms, double * x, double * y, double * z, double * q, double * mass);
void oqp_set_harmonic_active(bool flag);
void oqp_banner(struct oqp_handle_t *inf);

void apply_basis(struct oqp_handle_t *inf);

void append_shell(struct oqp_handle_t *inf);
void append_ecp(struct oqp_handle_t *inf);

void int1e(struct oqp_handle_t *inf);

void int2e(struct oqp_handle_t *inf);

void guess_hcore(struct oqp_handle_t *inf);
void guess_huckel(struct oqp_handle_t *inf);
void guess_modhuckel(struct oqp_handle_t *inf);
void guess_json(struct oqp_handle_t *inf);
void guess_sap(struct oqp_handle_t *inf);
void guess_minao(struct oqp_handle_t *inf);
void proj_dm_newbas(struct oqp_handle_t *inf);

void hf_energy(struct oqp_handle_t *inf);
void hf_gradient(struct oqp_handle_t *inf);
void fci_ao_integrals(struct oqp_handle_t *inf);
/* Matrix-free QDPT2 streaming kernel (diagonal-Fock H0): accumulates the
 * per-state first-order couplings V_I(Phi) and diagonal E0(Phi) over the
 * external singles/doubles of the reference determinants into a hash table.
 * Handle-free: pure array in/out; cap must be a power of two; returns the
 * number of unique external determinants, or -1 on hash overflow. */
int64_t qdpt2_stream_kernel(int32_t norb, int32_t ncore, int32_t nact,
    int64_t nsup, int32_t nstate, int32_t nthreads, int64_t cap,
    const int64_t *sup_a, const int64_t *sup_b, const double *cvec,
    const double *h1e, const double *eri, const double *eps,
    int64_t *out_ka, int64_t *out_kb, double *out_e0, double *out_v);
/* AO -> MO integral transformation (mo_transform.F90). The four-index
 * transform is four DGEMMs with transposed output, so the index roll between
 * quarter transformations is free and no nbf^4 temporary is transposed.
 * mo_transform_eri returns 0, or -1 when its one internal nbf^4 buffer could
 * not be allocated (the caller falls back to the Python einsum). */
void mo_transform_h1e(int32_t nbf, const double *hcore, const double *coeff,
    double *h1e);
int64_t mo_transform_eri(int32_t nbf, const double *eri_ao, const double *coeff,
    double *eri_mo);
/* Active-space setup (fci_setup.F90): the spin-orbital expansion of the
 * integrals and the frozen-core fold. Both are exact copies/sums, so they
 * reproduce the Python references bit for bit. `active` and `core` are 0-based
 * MO indices; `eref` receives the single inactive-energy scalar. */
void fci_spin_orbital_integrals(int32_t norb, const double *h1e,
    const double *eri, double *hspin, double *gspin, int32_t nthreads);
void fci_fold_core(int32_t norb, int32_t nact, int32_t ncore,
    const int32_t *active, const int32_t *core, const double *h1e,
    const double *eri, double *h_act, double *eref);
/* Working-set ceiling for the blocked scratch of the string-driven sigma and
 * RDM (fci_sigma_strings.F90), in bytes.  Set from the `[cas]`/`[fci]`/`[pt2]`
 * max_memory budget; without it the kernels use a built-in 512 MiB default.
 * The value is clamped to a 16 MiB floor, and a non-positive value restores
 * the default.  Both kernels block their alpha-string range to stay under it,
 * so this bounds scratch rather than deciding whether they can run at all. */
void fci_set_work_bytes_cap(int64_t bytes);
int64_t fci_get_work_bytes_cap(void);

/* <S^2> per root from the determinant CI vectors (fci_setup.F90). `dets` is in
 * CI order, so it is sorted internally with a permutation back to CI position
 * rather than searched directly. Returns 0, or -1 on allocation failure. */
int64_t fci_spin_square(int32_t norb, int64_t ndet, int32_t nvec,
    const int64_t *dets, const double *civec, int32_t ms2, double *s2,
    int32_t nthreads);
/* Determinant-CI Fortran engine (fci_hamiltonian.F90): ILP64-safe dense
 * symmetric eigensolver, dense Hamiltonian build, diagonal, and matrix-free
 * block application Y = 0.5(H+H^T) X for the Davidson solver. */
int64_t oqp_dsyevd(int32_t n, double *a, double *w);
void fci_dense_hamiltonian(int32_t nspin, int64_t ndet, const int64_t *dets,
    const double *hspin, const double *gspin, double cutoff, double *hmat,
    int32_t nthreads);
void fci_hamiltonian_diag(int32_t nspin, int64_t ndet, const int64_t *dets,
    const double *hspin, const double *gspin, double *diag);
int64_t fci_hamiltonian_matvec(int32_t nspin, int64_t ndet, const int64_t *dets,
    const double *hspin, const double *gspin, double cutoff, int32_t nvec,
    const double *x, double *y, int32_t nthreads);
/* ------------------------------------------------------------------------
 * One entry point for a complete CI solve (fci_driver.F90).
 *
 * Runs active-space setup, the frozen-core fold, the spin-orbital expansion,
 * determinant generation, the dense/Davidson dispatch, the Davidson iteration,
 * the spin diagnostics and root selection without returning to the caller --
 * the kernels above are what it is built from, and remain callable on their
 * own for the tests that pin them.
 *
 * Handle-free by design (like qdpt2_stream_kernel): the dominant caller is the
 * CASSCF microiteration, whose MO integrals are a line-search trial point that
 * does not live in the handle, and the eventual casscf_energy(inf) is itself
 * Fortran holding h1e/eri in local arrays.
 *
 * `h1e`/`eri` are the FULL MO-basis integrals over iopt[FCI_I_NORB] orbitals,
 * C-order; `active`/`core` are 0-based MO indices (pass a one-element dummy for
 * `core` when ncore == 0).  `energies` and `s2` hold nroot entries, `civecs` is
 * C-order [ndet, nroot] so root k is civecs[:, k].  `s2` is written only when
 * iopt[FCI_I_WANT_S2] is set.
 *
 * Returns the number of roots written (== nroot), or a negative status.  Every
 * negative value means "could not do it here"; the caller falls back to the
 * Python solve_fci, which remains the numerical pin and owns all user-facing
 * validation and error messages.
 *
 * OPTION SCHEMA -- authoritative copy is the parameter block in
 * source/modules/fci_driver.F90; pyoqp/oqp/library/fci.py mirrors it and
 * tests/test_fci_solve.py checks that the mirror is exact.  Only what the
 * compute path reads is packed here: the ~111 [cas]/[ci]/[casscf] keys are
 * Python-side sections with no representation in source/types.F90, and they do
 * not belong on the shared control_parameters struct. */
enum {
  FCI_I_NORB      = 0,  /* MOs spanned by h1e/eri                         */
  FCI_I_NACT      = 1,  /* active orbitals                                */
  FCI_I_NCORE     = 2,  /* frozen doubly-occupied orbitals                */
  FCI_I_NALPHA    = 3,  /* active alpha electrons                         */
  FCI_I_NBETA     = 4,  /* active beta electrons                          */
  FCI_I_NROOT     = 5,  /* roots requested and returned                   */
  FCI_I_SOLVER    = 6,  /* 0 auto, 1 dense, 2 davidson                    */
  FCI_I_MAXITER   = 7,  /* Davidson iteration cap                         */
  FCI_I_SUBSPACE  = 8,  /* Davidson subspace cap, 0 = auto                */
  FCI_I_MULT      = 9,  /* target spin multiplicity, 0 = any              */
  FCI_I_MAXMEMORY = 10, /* working-set budget, MiB                        */
  FCI_I_NTHREADS  = 11, /* OpenMP threads for the kernels                 */
  FCI_I_WANT_S2   = 12, /* 1 = also return <S^2> per returned root        */
  FCI_I_GUESS     = 13, /* 1 = civecs holds nroot Davidson start vectors  */
  FCI_NIOPT       = 14
};
enum {
  FCI_D_ECORE     = 0,  /* scalar added to every returned root            */
  FCI_D_EIG_TOL   = 1,  /* eigenpair residual tolerance                   */
  FCI_D_CUTOFF    = 2,  /* integral screening cutoff                      */
  FCI_NDOPT       = 3
};
int64_t fci_solve(const int32_t *iopt, const double *dopt,
    const int32_t *active, const int32_t *core,
    const double *h1e, const double *eri,
    double *energies, double *civecs, double *s2);
/* Spin-orbital determinant RDMs (rdm_kernel.F90). rdm2_spinorb returns 0 on
 * success and -1 when `cap` was too small for the reachable intermediates,
 * in which case the caller falls back to the Python enumeration. */
void rdm1_spinorb(int32_t nspin, int64_t ndet, const int64_t *dets,
    const double *civec, double *d1);
int64_t rdm2_spinorb(int32_t nspin, int64_t ndet, const int64_t *dets,
    const double *civec, int64_t cap, double *d2, int32_t nthreads);
/* Spin-summed spatial RDMs. Same build as the spin-orbital pair above, but the
 * [2n,...] spin-orbital tensor is never materialised: the spin sum is applied
 * while unpacking, so rdm2_spatial writes n^4 instead of 16 n^4. */
void rdm1_spatial(int32_t norb, int64_t ndet, const int64_t *dets,
    const double *civec, double *d1);
int64_t rdm2_spatial(int32_t norb, int64_t ndet, const int64_t *dets,
    const double *civec, int64_t cap, double *d2, int32_t nthreads);
/* String-factorized spatial d1+d2 of one CI vector on the canonical CAS
 * product determinant list (fci_sigma_strings.F90); nonzero status = the list
 * is not the canonical product, fall back to the walking kernels above. */
int64_t rdm12_strings_c(int32_t norb, int64_t ndet, const int64_t *dets,
    const double *civec, double *d1, double *d2, int32_t nthreads);
/* Spin-free dm1..dm4 (PySCF make_dm1234 convention) from the determinant CI
 * vector; `upto` selects how many are produced and the caller must have
 * allocated every array it requests. Returns 0, or -1 on allocation failure. */
int64_t nevpt2_make_rdms(int32_t norb, int64_t ndet, const int64_t *dets,
    const double *civec, int32_t upto, double *dm1, double *dm2,
    double *dm3, double *dm4);
/* CASSCF generalized Fock and orbital gradient (casscf_kernel.F90). `gamma`
 * and `gamma2` carry one active 1-/2-RDM per state-average root; the weighted
 * sum is accumulated internally. `grad` is filled over the non-redundant
 * rotation pairs in casscf.py `_nonredundant_pairs` order. */
void casscf_gfock_grad(int32_t nbf, int32_t ncore, int32_t nact, int32_t nroot,
    const double *weights, const double *gamma, const double *gamma2,
    const double *h1e, const double *eri, double *fock, double *grad);
/* CI-relaxation amplitudes of the analytic CASSCF orbital Hessian
 * (casscf_hess_kernel.F90): the per-pair derivative operator applied to the
 * reference CI vector and projected on the active eigenbasis, batched over
 * the rotation pairs. */
void casscf_hess_amp(int32_t nact, int64_t ndet, int32_t npar,
    const double *stack, const double *fder, const double *gder,
    const double *wmat, const double *vecs, double *amp);
/* W_tua = (E_tu|c>)_a for one CI vector (casscf_hess_kernel.F90). */
void casscf_hess_wmat(int32_t nact, int64_t ndet, const double *stack,
    const double *civec, double *wmat);
/* CI-relaxation accumulation for one averaged root (casscf_hess_kernel.F90).
 * Returns 0, or j+1 when eigenstate j is degenerate with the reference and
 * still genuinely coupled -- the caller raises, since no orbital Hessian
 * exists there. */
int64_t casscf_hess_relax(int32_t npar, int64_t ndet, int32_t navg, int32_t iavg,
    const double *ovl, const double *weights, const double *eps, double e_i,
    double degen_tol, double noise_tol, const double *amp, double *hess);
/* Fixed-CI half of the analytic CASSCF orbital Hessian (casscf_hess_bmat.F90):
 * the B columns and the folded active derivative integrals, assembled from the
 * sparse one-index derivative slabs so no nbf^4 temporary is formed. */
void casscf_hess_bmat(int32_t nbf, int32_t ncore, int32_t nact, int32_t npar,
    const int32_t *pairs, const double *dmat, const double *rdm2,
    const double *h1e, const double *eri, double *bmat, double *fder,
    double *gder);
/* Spin-summed excitation matrices (E_tu)_{row,col} over the determinant list
 * (casscf_exc_stack.F90). Returns 0, or -1 when 2*nact exceeds the 62-bit
 * determinant encoding, in which case the caller falls back to Python. */
int64_t casscf_excitation_stack(int32_t nact, int64_t ndet, const int64_t *dets,
    double *stack);
/* Determinant-space bookkeeping and mean-field Fock of the PT2 path
 * (pt2_kernel.F90).  Determinant keys are the fci.py integers, so these need
 * 2*norb <= 62.  pt2_external_indices returns the number of external
 * determinants written; pt2_occupation_blocks returns the block count and
 * fills `starts` with nblock+1 boundaries into the `order` permutation. */
void pt2_effective_fock(int32_t nbf, const double *h1e, const double *eri,
    const double *dmat, double *fock);
void pt2_h0_dyall_active(int32_t nbf, int32_t ncore, int32_t nact,
    const double *h1e, const double *eri, double *hact);
int64_t pt2_external_indices(int32_t norb, int32_t ncore, int32_t nact,
    int64_t ndet, const int64_t *dets, int64_t *ext);
void pt2_diagonal_h0(int32_t norb, int64_t ndet, const int64_t *dets,
    const double *eps, double *diag);
int64_t pt2_occupation_blocks(int32_t norb, int32_t ncore, int32_t nact,
    int64_t next, const int64_t *dets, const int64_t *ext, int64_t *order,
    int64_t *starts);
/* Loewdin minors of the MS/XMS orbital rotation, T[t,s] = det(R[t_occ,s_occ]).
 * `occ` is the [nstr,nelec] table of occupied orbitals per single-spin string,
 * in _string_basis order. */
void pt2_minor_transform(int32_t norb, int32_t nelec, int64_t nstr,
    const double *r, const int64_t *occ, double *tmat);
/* Dominant Koopmans intermediates of the strongly contracted NEVPT2
 * (nevpt2_koopmans.F90).  `h2e` is the physicist-ordered active ERI tensor and
 * dm3/dm4 the spin-free RDMs in the PySCF make_dm1234 convention. */
void nevpt2_f3ca_f3ac(int32_t nact, const double *h2e, const double *dm4,
    double *f3ca, double *f3ac);
void nevpt2_a16(int32_t nact, const double *h1e, const double *h2e,
    const double *dm3, const double *f3ca, const double *f3ac, double *a16);
void nevpt2_a22(int32_t nact, const double *h1e, const double *h2e,
    const double *dm2, const double *dm3, const double *f3ca,
    const double *f3ac, double *a22);
/* Sij / Srs intermediates.  `hdm1` is an argument of the Python _hdm3 and _a9
 * but never appears in their bodies, so it is not passed.  nevpt2_a7 returns
 * the reduced 2-RDM alongside a7, as the Python does. */
void nevpt2_hdm3(int32_t nact, const double *dm1, const double *dm2,
    const double *dm3, const double *hdm2, double *hdm3);
void nevpt2_a9(int32_t nact, const double *h1e, const double *h2e,
    const double *hdm2, const double *hdm3, double *a9);
void nevpt2_a7(int32_t nact, const double *h1e, const double *h2e,
    const double *dm1, const double *dm2, const double *dm3,
    double *rm2, double *a7);
/* Sir intermediates.  a12's RDM arguments are indexed [q,p,...] while its
 * output is [p,q,a,b]; a13's carry only q as a spectator. */
void nevpt2_a12(int32_t nact, const double *h1e, const double *h2e,
    const double *dm2, const double *dm3, double *a12);
void nevpt2_a13(int32_t nact, const double *h1e, const double *h2e,
    const double *dm1, const double *dm2, const double *dm3, double *a13);
/* Closed+active mean-field Fock h + J - K/2 used to canonicalize the CASSCF
 * orbitals (casscf_kernel.F90); shares its J/K builder with the generalized
 * Fock above. */
void casscf_effective_fock(int32_t nbf, const double *dmat, const double *h1e,
    const double *eri, double *fout);
/* CASSCF orbital rotation C <- C exp(K) (casscf_orbrot.F90). `pairs` is the
 * non-redundant rotation-pair list in casscf.py `_nonredundant_pairs` order;
 * the antisymmetric K is built from it and exponentiated by scaling-and-
 * squaring with the degree-13 Pade approximant. Both return 0, or the LAPACK
 * info of a failed Pade solve. */
int32_t casscf_orbital_rotate(int32_t nbf, int32_t npar, const int32_t *pairs,
    const double *vec, const double *cin, double *cout);
int32_t casscf_expm(int32_t n, const double *kin, double *eout);
/* CASSCF orbital-optimization convergers (casscf_ah.F90). casscf_ah_model_step
 * returns 0 on success, 1 when the augmented-Hessian eigenvector has no
 * reference component (the caller steps along the lowest mode instead), and 2
 * when LAPACK failed. casscf_diis_coeffs writes `nused` coefficients for the
 * LAST `nused` stored gradients; nused == 0 means no stable extrapolation. */
int32_t casscf_ah_model_step(int32_t npar, const double *grad, const double *w,
    const double *umat, double trust, int32_t max_micro, double v0tol,
    double *step, double *shift, double *pred, int32_t *nmic);
void casscf_lowest_mode_step(int32_t npar, const double *grad, const double *w,
    const double *umat, double trust, double *step, double *pred);
void casscf_diis_coeffs(int32_t nvec, int32_t npar, const double *gmat,
    double condmax, double *coef, int32_t *nused);
/* ------------------------------------------------------------------------
 * One entry point for a complete CASSCF / SA-CASSCF run (casscf_driver.F90).
 *
 * The wavefunction stack's counterpart to hf_energy(inf) / mp2_energy(inf):
 * the whole default two-phase orbital optimizer runs inside liboqp -- the
 * macroiteration loop with its backtracking line search, the 2*n_par
 * finite-difference orbital Hessian, the level-shifted Newton step, the
 * curvature-gated saddle escape, canonicalization, the final CI and its spin
 * diagnostics -- built out of fci_solve and the CASSCF kernels above.
 * `OQP::Hcore`, `OQP::AO_ERI` and `OQP::VEC_MO_A` are read from the handle and
 * `OQP::VEC_MO_A` is overwritten with the optimized orbitals.
 *
 * `weights`/`roots` carry the resolved state-average plan (NSTATE entries;
 * `roots` are 0-based CI root indices).  `energies` and `s2` receive NROOT
 * entries.  `history` is the macroiteration table, C-order [MAXHIST, 5] as
 * (it, E, dE, |g_orb|, |step|), so Python formats the log unchanged; `stats`
 * receives 10 counters -- (rows written, macroiterations, converged,
 * evaluations, analytic Hessian builds, DIIS extrapolations accepted, DIIS
 * extrapolations attempted, `ah` stagnation flag, `auto` outcome, `auto`
 * macroiterations before falling back) -- from which Python formats the
 * converger trace.
 *
 * The `twophase` (default), `ah`, `diis` and `auto` convergers all run here,
 * with either orbital-Hessian builder (`CAS_I_HESSIAN`).
 *
 * Returns 0, or a negative status meaning "could not do it here" -- the caller
 * then runs the Python optimizer, which remains the numerical pin and owns
 * every user-facing message.  Two negative statuses are refusals rather than
 * fallbacks and the Python re-raises their messages: -8, a root degeneracy
 * with non-zero orbital coupling (the state-averaged objective is not smooth,
 * so no orbital Hessian exists), and -9, an excitation stack past the
 * dense-spectrum memory guard.
 *
 * OPTION SCHEMA -- authoritative copy is the parameter block in
 * source/modules/casscf_driver.F90; pyoqp/oqp/library/casscf.py mirrors it and
 * tests/test_casscf_energy.py checks that the mirror is exact.  The CI half of
 * the run reuses the FCI_I_* / FCI_D_* layout above verbatim. */
enum {
  CAS_I_NCORE     = 0,  /* inactive doubly-occupied orbitals               */
  CAS_I_NACT      = 1,  /* active orbitals                                 */
  CAS_I_NALPHA    = 2,  /* active alpha electrons                          */
  CAS_I_NBETA     = 3,  /* active beta electrons                           */
  CAS_I_NSTATE    = 4,  /* averaged roots in weights/roots                 */
  CAS_I_NROOT     = 5,  /* CI roots solved and returned                    */
  CAS_I_SOLVER    = 6,  /* 0 auto, 1 dense, 2 davidson                     */
  CAS_I_MAXITER   = 7,  /* Davidson iteration cap                          */
  CAS_I_SUBSPACE  = 8,  /* Davidson subspace cap, 0 = auto                 */
  CAS_I_MULT      = 9,  /* target spin multiplicity, 0 = any               */
  CAS_I_MAXMEMORY = 10, /* CI working-set budget, MiB                      */
  CAS_I_NTHREADS  = 11, /* OpenMP threads for the kernels                  */
  CAS_I_MAXMACRO  = 12, /* macroiteration cap                              */
  CAS_I_OPTIMIZER = 13, /* 0 newton, 1 powell                              */
  CAS_I_CANONICAL = 14, /* 1 = canonicalize the converged orbitals         */
  CAS_I_MAXESCAPE = 15, /* saddle escapes, 0 disables phase 2              */
  CAS_I_MAXHIST   = 16, /* rows the caller's history buffer holds          */
  CAS_I_CONVERGER = 17, /* 0 twophase, 1 ah, 2 diis, 3 auto                */
  CAS_I_HESSIAN   = 18, /* 0 finite difference, 1 analytic                 */
  CAS_I_AH_MICRO  = 19, /* AH scale-bisection microiterations              */
  CAS_I_AH_REJECT = 20, /* uphill-step rejections per macroiteration       */
  CAS_I_DIIS_SPACE= 21, /* stored rotation/gradient pairs                  */
  CAS_I_DIIS_START= 22, /* pairs required before extrapolating             */
  CAS_I_AUTO_STAG = 23, /* stalled macroiterations before the fallback     */
  CAS_NIOPT       = 24
};
enum {
  CAS_D_ENUC      = 0,  /* nuclear repulsion, added to every root          */
  CAS_D_EIG_TOL   = 1,  /* CI eigenpair residual tolerance                 */
  CAS_D_CUTOFF    = 2,  /* CI integral screening cutoff                    */
  CAS_D_GRAD_TOL  = 3,  /* |g_orb| convergence threshold                   */
  CAS_D_ENER_TOL  = 4,  /* macroiteration energy-decrease threshold        */
  CAS_D_STEP_TOL  = 5,  /* macroiteration step-norm threshold              */
  CAS_D_MAXROT    = 6,  /* trust cap on |step|                             */
  CAS_D_SHIFT     = 7,  /* Hessian eigenvalue floor (level shift)          */
  CAS_D_FD_STEP   = 8,  /* finite-difference Hessian displacement          */
  CAS_D_SADDLE_C  = 9,  /* deep-negative-curvature threshold, phase 2      */
  CAS_D_SADDLE_E  = 10, /* strict energy gain to accept an escape          */
  CAS_D_AH_START  = 11, /* AH initial trust radius                         */
  CAS_D_AH_MAXTR  = 12, /* AH trust-radius ceiling                         */
  CAS_D_AH_MINTR  = 13, /* AH trust-radius floor / stagnation              */
  CAS_D_AH_SADC   = 14, /* AH deep-negative-curvature threshold            */
  CAS_D_AH_SADE   = 15, /* AH strict energy gain to accept an escape       */
  CAS_NDOPT       = 16
};
int64_t casscf_energy(struct oqp_handle_t *inf, const int32_t *iopt,
    const double *dopt, const double *weights, const int32_t *roots,
    double *energies, double *s2, double *history, int32_t *stats);
/* ------------------------------------------------------------------------
 * Analytic state-specific CASSCF NUCLEAR gradient (casscf_gradient.F90).
 *
 * Not to be confused with casscf_gfock_grad above, which is the derivative
 * with respect to the ORBITAL ROTATION parameters.  That gradient vanishing is
 * a precondition for this one: the expression implemented here is the fully
 * variational
 *
 *   dE/dx = sum D h^x + 1/2 sum Gamma (..|..)^x - sum X S^x + dV_NN/dx
 *
 * with no orbital-response or CI-response (Z-vector) term, which is valid only
 * at a converged state-specific solution.
 *
 * Takes the SAME iopt/dopt schema as casscf_energy, so a gradient is requested
 * with the arrays that configured the energy it differentiates; the
 * optimizer-only entries are ignored.  `weights` must be the single entry 1.0
 * and `roots` the single 0-based root the orbitals were optimized for
 * ([casscf] root): a state-averaged run is REFUSED (status -20) rather than
 * given the state-specific formula.  Individual SA state gradients need a
 * Lagrangian/Z-vector response that is not implemented.
 *
 * Reads OQP::Hcore, OQP::AO_ERI and OQP::VEC_MO_A; writes infos%atoms%grad.
 *
 * `info` receives CAS_G_NINFO diagnostics, of which the first two let the
 * caller reject a gradient taken at a non-stationary point.
 *
 * Returns 0, or a negative status: -1..-9 are the casscf_energy codes,
 * -20 state-averaged run, -21 non-Hartree-Fock Hamiltonian, -22 allocation,
 * -23 eigensolver failure, -24 the CI / generalized-Fock evaluation failed. */
enum {
  CAS_G_GNORM  = 0,  /* |g_orb|, the orbital-rotation gradient norm      */
  CAS_G_ENERGY = 1,  /* energy of the differentiated root                */
  CAS_G_FASYM  = 2,  /* max |F_pq - F_qp|, stationarity probe            */
  CAS_G_NVEC   = 3,  /* rank of the factorized active 2-RDM correction   */
  CAS_G_NINFO  = 4
};
int64_t casscf_gradient(struct oqp_handle_t *inf, const int32_t *iopt,
    const double *dopt, const double *weights, const int32_t *roots,
    double *info);
/* ------------------------------------------------------------------------
 * Nuclear gradient from a FACTORIZED AO two-particle density
 * (casscf_ao_gradient.F90).
 *
 *   dE/dx = sum D h^x + 1/2 sum Gamma (..|..)^x - sum X S^x + dV_NN/dx
 *
 * with X the energy-weighted density of the term as written (the sign the
 * overlap-derivative kernel wants is applied inside), and
 *
 *   Gamma_{mu nu la si} =
 *       sum_t sepc[t] [ P^t_{mu nu} P^t_{la si}
 *                       - 1/4 (P^t_{mu la} P^t_{nu si}
 *                            + P^t_{nu la} P^t_{mu si}) ]
 *     + sum_k lam[k] A^k_{mu nu} A^k_{la si}.
 *
 * The first family is the Hartree-Fock separable form at an arbitrary
 * symmetric matrix and carries every separable and BILINEAR piece by
 * polarization; the second is the Coulomb-type product an all-active
 * eigendecomposition produces.  Both are eight-fold symmetric and expand to
 * the Cartesian-effective basis matrix by matrix, which a dense nbf^4 AO
 * tensor does not.
 *
 * `dmat`, `xmat`, `sepm` and `avm` are C-order [nbf,nbf] (or [n,nbf,nbf])
 * buffers; every matrix must be symmetric, and is symmetrized on entry.
 * `nsep >= 1` is required.  Writes infos%atoms%grad.
 *
 * Returns 0, or -30 (sizes disagree with the handle) / -31 (allocation). */
enum {
  CAS_AG_NSEP  = 0,  /* separable two-particle factors contracted        */
  CAS_AG_NVEC  = 1,  /* low-rank two-particle factors contracted         */
  CAS_AG_NINFO = 2
};
int64_t casscf_ao_gradient(struct oqp_handle_t *inf, int32_t nbf,
    const double *dmat, const double *xmat,
    int32_t nsep, const double *sepc, const double *sepm,
    int32_t nvec, const double *lam, const double *avm,
    double *info);
void hf_hessian(struct oqp_handle_t *inf);
void hess1_selftest(struct oqp_handle_t *inf);
void grd2_hess_selftest(struct oqp_handle_t *inf);
void cholesky_eri_selftest(struct oqp_handle_t *inf);
void hess_skel_selftest(struct oqp_handle_t *inf);
void hess_skel_open_selftest(struct oqp_handle_t *inf);
void electric_dipole_au(struct oqp_handle_t *inf, double *dipole);
void cphf_static_polarizability(struct oqp_handle_t *inf, double *alpha);
void vibrational_intensities_native(struct oqp_handle_t *inf, int64_t nmode, int64_t ncoord,
        double *modes, double *dipole_derivatives, double *polarizability_derivatives,
        double *infrared_intensities, double *mode_dipole_derivatives,
        double *raman_activities, double *mode_polarizability_derivatives);
void cphf_polarizability_selftest(struct oqp_handle_t *inf);
void cphf_uhf_polarizability_selftest(struct oqp_handle_t *inf);
void cphf_rohf_polarizability_selftest(struct oqp_handle_t *inf);
void fockx_selftest(struct oqp_handle_t *inf);
void fockx_os_selftest(struct oqp_handle_t *inf);

void tdhf_energy(struct oqp_handle_t *inf);
void tdhf_z_vector(struct oqp_handle_t *inf);
void tdhf_gradient(struct oqp_handle_t *inf);
void tdhf_hessian(struct oqp_handle_t *inf);

void tdhf_sf_energy(struct oqp_handle_t *inf);
void tdhf_sf_z_vector(struct oqp_handle_t *inf);
void tdhf_sf_gradient(struct oqp_handle_t *inf);
void tdhf_sf_hessian(struct oqp_handle_t *inf);

void tdhf_mrsf_energy(struct oqp_handle_t *inf);
void tdhf_umrsf_energy(struct oqp_handle_t *inf);
void tdhf_mrsf_ekt_ip(struct oqp_handle_t *inf);
void tdhf_mrsf_ekt_ea(struct oqp_handle_t *inf);
void tdhf_mrsf_z_vector(struct oqp_handle_t *inf);
void tdhf_mrsf_gradient(struct oqp_handle_t *inf);

void mp2_energy(struct oqp_handle_t *inf);
void mp2_gradient(struct oqp_handle_t *inf);
void ccsd_t_energy(struct oqp_handle_t *inf);

void electric_moments(struct oqp_handle_t *inf);
void electric_moments_excited(struct oqp_handle_t *inf);
void get_structures_ao_overlap(struct oqp_handle_t *inf);
void get_states_overlap(struct oqp_handle_t *inf);
void mrsf_namd_hop(struct oqp_handle_t *inf);
double oqp_namd_counter_random(int64_t seed, int64_t stream, int64_t step);
void oqp_namd_counter_normal_fill(int64_t seed, int64_t stream, int64_t count,
        double *values);
int oqp_namd_baeck_an_tdc(int64_t nstate, double dt_left, double dt_right,
        double gap_max, const double *energies_old,
        const double *energies_center, const double *energies_current,
        double *tdc_row_major);
int oqp_namd_nacme_gate(int64_t nstate, const double *candidate,
        const double *reference, const int *reference_mask, int compare_mode,
        double invariant_tol, double abs_tol, double rel_tol, double *metrics,
        int64_t *counts);
/* Native ODP CV types: 1 distance, 2 d(a,b)-d(c,d), 3 angle a-b-c. */
int oqp_odp_umbrella_evaluate(int64_t natom, int64_t ncv,
        const int *cv_types, const int64_t *cv_atoms,
        const double *cv_scales, const double *reference_r,
        const double *reference_p, double center, double k_parallel,
        double k_perpendicular, const double *coordinates, double *energy,
        double *force, double *xi, double *cv_raw, double *cv_scaled,
        double *cv_perpendicular, double *perpendicular_norm,
        double *energy_parallel, double *energy_perpendicular);
int oqp_namd_droplet_boundary(int64_t natom, int64_t ngroup,
        const double *coordinates, const double *masses,
        const int64_t *group_index, const double *center, double radius,
        double buffer, double force_constant, double max_penetration_limit,
        double *energy, double *forces, double *max_penetration,
        int64_t *active_count);
int oqp_namd_com_restraint(int64_t natom, const double *coordinates,
        const double *masses, const int64_t *selected, const double *center,
        double force_constant, double *energy, double *forces,
        double *displacement_norm);
int oqp_namd_langevin_thermostat(int64_t natom, double dt,
        double temperature, double friction, int64_t seed, int64_t stream,
        int64_t step, const double *masses, double *velocities, double *heat);
int oqp_maximum_overlap_assignment(int n, const double *overlap_row_major,
        int *assignment, double *signs, double *matched, double *margins);
int oqp_diagonal_phase_tracking(int n, const double *overlap_row_major,
        double *signs, double *matched, double *margins);
void resp_charges(struct oqp_handle_t *inf);
void mulliken(struct oqp_handle_t *inf);
void mulliken_excited(struct oqp_handle_t *inf);
void lowdin(struct oqp_handle_t *inf);

/* Native DFT-D4 dispersion (source/dftd4_interface.F90). Self-contained:
   takes atomic numbers + coordinates (Bohr), not the oqp handle.  The legacy
   entry point preserves its historical neutral-system behavior. */
void oqp_dftd4_disp(int nat, const int *z, const double *xyz,
                    const char *func, int lfunc, int do_grad,
                    double *energy, double *grad, int *ier);

/* Charge-aware ABI. param_mode=0 selects damping by functional name;
   param_mode=1 reads damping as [s6, s8, s9, a1, a2, alp]. */
void oqp_dftd4_disp_v2(int nat, const int *z, const double *xyz,
                       double total_charge, const char *func, int lfunc,
                       int param_mode, const double *damping, int do_grad,
                       double *energy, double *grad, int *ier);

void soc_mrsf(struct oqp_handle_t *inf);
void dk_scalar(struct oqp_handle_t *inf);
void nmr_shielding(struct oqp_handle_t *inf);
void nmr_giao_shielding_debug(struct oqp_handle_t *inf);
void nmr_giao_shielding(struct oqp_handle_t *inf);

void espf_op_corr(struct oqp_handle_t *inf);
void form_esp_charges(struct oqp_handle_t *inf);
void grad_esp_qmmm(struct oqp_handle_t *inf);
void grad_esp_qmmm_excited(struct oqp_handle_t *inf);
