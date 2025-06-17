#include <stdint.h>
#include <stdbool.h>

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
    double    diis_method_threshold;
    int64_t   diis_type;
    double    vdiis_cdiis_switch;
    double    vdiis_vshift_switch;
    double    vshift_cdiis_switch;
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
    int64_t   esp;
    int64_t   resp_target;
    double    esp_constr;
    bool      basis_set_issue;
    double    conf_print_threshold;
    bool      rstctmo;
    // SOSCF parameters
    int64_t   soscf_type;
    double    soscf_lvl_shift;
    int64_t   soscf_reset_mod;
    int64_t   verbose;
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
int64_t oqp_get(struct oqp_handle_t *c_handle, char *code,
        int32_t *type_id, int32_t *ndims, int64_t *dims, void **v);
int64_t oqp_alloc(struct oqp_handle_t *c_handle, char *code,
        int32_t *type_id, int32_t *ndims, int64_t *dims, void **v);
int64_t oqp_del(struct oqp_handle_t *c_handle, char *code);
int64_t oqp_get_nbf(struct oqp_handle_t *c_handle);
int64_t oqp_get_basis(struct oqp_handle_t *c_handle,
        int64_t *nsh, int64_t *nprim, int64_t *nbf,
        int64_t **bt, int64_t **at, int64_t **cdeg, double **ex, double **cc);

/* `mass` is optional, pass NULL if not needed */
int oqp_set_atoms(struct oqp_handle_t * c_handle, int64_t natoms, double * x, double * y, double * z, double * q, double * mass);
void oqp_banner(struct oqp_handle_t *inf);

void apply_basis(struct oqp_handle_t *inf);

void append_shell(struct oqp_handle_t *inf);
void append_ecp(struct oqp_handle_t *inf);

void int1e(struct oqp_handle_t *inf);

void guess_hcore(struct oqp_handle_t *inf);
void guess_huckel(struct oqp_handle_t *inf);
void guess_json(struct oqp_handle_t *inf);
void proj_dm_newbas(struct oqp_handle_t *inf);

void hf_energy(struct oqp_handle_t *inf);
void hf_gradient(struct oqp_handle_t *inf);

void tdhf_energy(struct oqp_handle_t *inf);
void tdhf_z_vector(struct oqp_handle_t *inf);
void tdhf_gradient(struct oqp_handle_t *inf);

void tdhf_sf_energy(struct oqp_handle_t *inf);
void tdhf_sf_z_vector(struct oqp_handle_t *inf);
void tdhf_sf_gradient(struct oqp_handle_t *inf);

void tdhf_mrsf_energy(struct oqp_handle_t *inf);
void tdhf_mrsf_z_vector(struct oqp_handle_t *inf);
void tdhf_mrsf_gradient(struct oqp_handle_t *inf);

void electric_moments(struct oqp_handle_t *inf);
void get_structures_ao_overlap(struct oqp_handle_t *inf);
void get_states_overlap(struct oqp_handle_t *inf);
void resp_charges(struct oqp_handle_t *inf);
void mulliken(struct oqp_handle_t *inf);
void lowdin(struct oqp_handle_t *inf);

