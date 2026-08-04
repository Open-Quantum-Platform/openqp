!> @brief One bind(C) entry point for a complete CASSCF / SA-CASSCF run.
!>
!> Every other method in OpenQP exposes a single C entry point per method --
!> `hf_energy(inf)`, `mp2_energy(inf)` -- with everything inside it Fortran and
!> Python only parsing options and reading results.  The wavefunction stack
!> instead exposed ~20 fine-grained kernels that a Python loop orchestrated:
!> for one CASSCF macroiteration the boundary was crossed once per AO->MO
!> transform, once per CI solve, twice per state-average root for the RDMs,
!> once per generalized-Fock build and once per orbital rotation -- and the
!> default finite-difference orbital Hessian runs 2*n_par of those per
!> macroiteration, so a 20-macroiteration CAS(6,6)/cc-pVDZ run crossed the
!> boundary ~2.4e5 times.  Considerable effort went into making that crossing
!> cheap (cached backend lookups, no-op array conversions removed); this
!> removes the seam instead.  `casscf_energy` crosses it ONCE per run.
!>
!> What runs inside
!> ----------------
!> Every `[casscf] converger`, over either orbital-Hessian builder:
!>
!>   * `twophase` (default): the macroiteration loop
!>     `casscf.py::_floored_newton_loop` -- level-shifted (eigenvalue-floored)
!>     Newton with the 12-step halving backtrack -- followed by the
!>     curvature-gated saddle-escape phase `_escape_saddles`, including its
!>     six-amplitude kick line search and its strict-energy-gain acceptance;
!>   * `ah`: the trust-region augmented-Hessian loop of
!>     `casscf_convergers.py::_ah_inner`, its rejection/radius bookkeeping and
!>     its own `_curvature_escape` (which re-converges with the AH loop, and
!>     stamps rerun history rows differently from `_escape_saddles` -- both
!>     conventions are reproduced);
!>   * `diis`: `_diis_optimize`, i.e. Pulay extrapolation of the accumulated
!>     rotation over the production step, with every candidate re-evaluated
!>     variationally so acceleration can only match or improve;
!>   * `auto`: `_auto_optimize`, `ah` under a stagnation watchdog with an
!>     automatic fallback to `twophase` from the best point reached;
!>   * the orbital Hessian itself, either the symmetrized finite-difference one
!>     (2*n_par gradient evaluations, which dominate the run) or the exact
!>     analytic one (`casscf_anhess.F90`, zero extra CI solves);
!>   * canonicalization (`_canonicalize` + `_effective_fock`);
!>   * the final CI on the optimized orbitals and its spin diagnostics;
!>   * the commit of the optimized orbitals into `OQP::VEC_MO_A`.
!>
!> Each of these is built out of the kernels that already existed and are still
!> individually callable and individually pinned: `fci_solve` (fci_driver.F90),
!> `mo_transform_h1e`/`mo_transform_eri`, `rdm1_spatial`/`rdm2_spatial`,
!> `casscf_gfock_grad`, `casscf_effective_fock`, `casscf_orbital_rotate`,
!> `oqp_dsyevd_f`.  Calling exactly those, in exactly the order the Python
!> driver called them, is what makes the native path reproduce the Python one
!> to the last bit everywhere the two use the same LAPACK.
!>
!> This is precisely the in-process Fortran caller `fci_solve` was shaped
!> handle-free for: the microiteration's MO integrals are a line-search trial
!> point held in local arrays here, never in `infos%dat`, so routing them
!> through the handle would have cost an nbf^4 copy per evaluation.
!>
!> What stays in Python, and why
!> -----------------------------
!>   * option parsing and validation -- integer bookkeeping executed once per
!>     run, and the source of every user-facing error message;
!>   * the state-average plan, resolved to explicit `weights`/`roots` arrays;
!>   * the log, including the macroiteration table and the converger trace: the
!>     driver returns the history rows and the counters and lets `_write_log`
!>     format them, so the printed output is byte-identical to the Python
!>     path's;
!>   * the converger option parsing, including which spellings are accepted --
!>     an unrecognized `converger` or `hessian` makes Python decline the driver
!>     so that `casscf_convergers.run_converger` / `_hessian_provider` raise the
!>     established message.
!>
!> The NumPy implementations in `casscf_convergers.py` and `casscf_hessian.py`
!> stay in the tree as the numerical pins and as the fallback for every
!> declined run -- the same `_gfock_backend()` pattern the rest of the stack
!> uses.
!>
!> Any negative status means "could not do it here" and the caller re-runs the
!> Python optimizer, which is still the numerical pin (see
!> tests/test_casscf_energy.py, which compares the two arms in the same
!> process so both use identical integrals and identical CI solves).  Two of
!> those statuses are REFUSALS, not fallbacks -- `CAS_ERR_DEGEN` (a root
!> degeneracy with live orbital coupling: the state-averaged objective is not
!> smooth and no orbital Hessian is defined) and `CAS_ERR_STACK` (an excitation
!> stack past the dense-spectrum memory guard).  Returning a finite Hessian in
!> the first case would be wrong, so the Python path re-runs and raises its own
!> message rather than the driver inventing one.
!>
!> Option schema
!> -------------
!> `[cas]`/`[casscf]`/`[ci]`/`[state_average]` have no representation in
!> `source/types.F90` and carry 111 keys between them; putting a method's
!> private options on the shared `control_parameters` struct would put them in
!> every method's namespace.  `fci_driver.F90` already answered this: pack only
!> what the compute path reads into small fixed-layout `iopt`/`dopt` arrays
!> whose authoritative schema is the parameter block below.  `include/oqp.h`
!> documents it, `pyoqp/oqp/library/casscf.py` mirrors it, and
!> `tests/test_casscf_energy.py::test_option_schema_matches_fortran` parses
!> THIS FILE to assert the mirror rather than trust it.
!>
!> The CI half of the run is configured with the *same* `iopt`/`dopt` layout
!> `fci_solve` already defines; the constants are imported from
!> `fci_driver_mod` rather than re-declared, so that schema cannot drift here.
!>
!> Array conventions follow the rest of the stack: everything crosses in the
!> caller's C order.  Note the one place that matters -- `OQP::VEC_MO_A` is
!> stored AO-fastest, i.e. the Fortran view `mo(ao, mo)`, while every kernel
!> here wants the *C-order* coefficient buffer `coeff[ao][mo]` (MO fastest).
!> The two are transposes and the matrix is not symmetric, so the driver
!> transposes explicitly on the way in and on the way out rather than
!> reinterpreting the buffer.
module casscf_driver_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use fci_driver_mod, only: fci_solve, build_determinants, &
                            FCI_I_NORB, FCI_I_NACT, FCI_I_NCORE, &
                            FCI_I_NALPHA, FCI_I_NBETA, FCI_I_NROOT, &
                            FCI_I_SOLVER, FCI_I_MAXITER, FCI_I_SUBSPACE, &
                            FCI_I_MULT, FCI_I_MAXMEMORY, FCI_I_NTHREADS, &
                            FCI_I_WANT_S2, FCI_I_GUESS, FCI_NIOPT, &
                            FCI_D_ECORE, FCI_D_EIG_TOL, FCI_D_CUTOFF, &
                            FCI_NDOPT, FCI_MAX_NSPIN
  use fci_hamiltonian_mod, only: oqp_dsyevd_f
  use mo_transform_mod, only: mo_transform_h1e, mo_transform_eri
  use casscf_kernel_mod, only: casscf_gfock_grad, casscf_effective_fock, &
                               casscf_gfock_grad_w
  use casscf_fastints_mod, only: cas_partial_eri, cas_ao_density, &
                                 cas_wmo_from_ao, cas_dmo_to_ao, &
                                 cas_act_exchange_tensor, cas_act_fock
  use fci_sigma_strings_mod, only: rdm12_strings
  use casscf_orbrot_mod, only: casscf_orbital_rotate
  use casscf_ah_mod, only: casscf_ah_model_step, casscf_lowest_mode_step, &
                           casscf_diis_coeffs
  use casscf_anhess_mod, only: cas_anhess_ctx_t, cas_anhess_init, &
                               cas_anhess_build, CAS_HESS_ERR_DEGEN, &
                               CAS_HESS_ERR_STACK
  use rdm_kernel_mod, only: rdm1_spatial, rdm2_spatial
  use trah_core_mod, only: trah_provider_t, trah_params_t, trah_result_t, trah_run
  ! DGEMM/DGEMV through the ILP64 wrapper layer (AGENTS.md rule 1); see the
  ! same note in fci_driver.F90 for why `only:` cannot be used here.
  use oqp_linalg
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: casscf_energy

  ! ------------------------------------------------------------------ schema
  ! AUTHORITATIVE OPTION SCHEMA -- mirrored in include/oqp.h and in
  ! pyoqp/oqp/library/casscf.py (_CAS_IOPT / _CAS_DOPT), and checked by
  ! tests/test_casscf_energy.py.  Indices are 0-based, as seen from C.
  !
  ! nbf, the AO ERIs, Hcore, the MO coefficients and the nuclear repulsion all
  ! come from the handle, so none of them appears here.
  !
  ! iopt (int32):
  integer, parameter :: CAS_I_NCORE     = 0   ! inactive doubly-occupied orbitals
  integer, parameter :: CAS_I_NACT      = 1   ! active orbitals
  integer, parameter :: CAS_I_NALPHA    = 2   ! active alpha electrons
  integer, parameter :: CAS_I_NBETA     = 3   ! active beta electrons
  integer, parameter :: CAS_I_NSTATE    = 4   ! averaged roots in weights/roots
  integer, parameter :: CAS_I_NROOT     = 5   ! CI roots solved and returned
  integer, parameter :: CAS_I_SOLVER    = 6   ! 0 auto, 1 dense, 2 davidson
  integer, parameter :: CAS_I_MAXITER   = 7   ! Davidson iteration cap
  integer, parameter :: CAS_I_SUBSPACE  = 8   ! Davidson subspace cap, 0 = auto
  integer, parameter :: CAS_I_MULT      = 9   ! target multiplicity, 0 = any
  integer, parameter :: CAS_I_MAXMEMORY = 10  ! CI working-set budget, MiB
  integer, parameter :: CAS_I_NTHREADS  = 11  ! OpenMP threads for the kernels
  integer, parameter :: CAS_I_MAXMACRO  = 12  ! macroiteration cap
  integer, parameter :: CAS_I_OPTIMIZER = 13  ! 0 newton, 1 powell
  integer, parameter :: CAS_I_CANONICAL = 14  ! 1 = canonicalize the orbitals
  integer, parameter :: CAS_I_MAXESCAPE = 15  ! saddle escapes, 0 disables phase 2
  integer, parameter :: CAS_I_MAXHIST   = 16  ! rows the history buffer holds
  integer, parameter :: CAS_I_CONVERGER = 17  ! 0 twophase, 1 ah, 2 diis, 3 auto
  integer, parameter :: CAS_I_HESSIAN   = 18  ! 0 finite difference, 1 analytic
  integer, parameter :: CAS_I_AH_MICRO  = 19  ! AH scale-bisection microiterations
  integer, parameter :: CAS_I_AH_REJECT = 20  ! uphill-step rejections per macro
  integer, parameter :: CAS_I_DIIS_SPACE= 21  ! stored rotation/gradient pairs
  integer, parameter :: CAS_I_DIIS_START= 22  ! pairs required before extrapolating
  integer, parameter :: CAS_I_AUTO_STAG = 23  ! stalled macros before the fallback
  integer, parameter :: CAS_NIOPT       = 24
  !
  ! dopt (double):
  integer, parameter :: CAS_D_ENUC      = 0   ! nuclear repulsion, added to every root
  integer, parameter :: CAS_D_EIG_TOL   = 1   ! CI eigenpair residual tolerance
  integer, parameter :: CAS_D_CUTOFF    = 2   ! CI integral screening cutoff
  integer, parameter :: CAS_D_GRAD_TOL  = 3   ! |g_orb| convergence threshold
  integer, parameter :: CAS_D_ENER_TOL  = 4   ! macroiteration energy-decrease threshold
  integer, parameter :: CAS_D_STEP_TOL  = 5   ! macroiteration step-norm threshold
  integer, parameter :: CAS_D_MAXROT    = 6   ! trust cap on |step|
  integer, parameter :: CAS_D_SHIFT     = 7   ! Hessian eigenvalue floor (level shift)
  integer, parameter :: CAS_D_FD_STEP   = 8   ! finite-difference Hessian displacement
  integer, parameter :: CAS_D_SADDLE_C  = 9   ! deep-negative-curvature threshold
  integer, parameter :: CAS_D_SADDLE_E  = 10  ! strict energy gain to accept an escape
  integer, parameter :: CAS_D_AH_START  = 11  ! AH initial trust radius
  integer, parameter :: CAS_D_AH_MAXTR  = 12  ! AH trust-radius ceiling
  integer, parameter :: CAS_D_AH_MINTR  = 13  ! AH trust-radius floor / stagnation
  integer, parameter :: CAS_D_AH_SADC   = 14  ! AH deep-negative-curvature threshold
  integer, parameter :: CAS_D_AH_SADE   = 15  ! AH strict gain to accept an escape
  integer, parameter :: CAS_NDOPT       = 16

  ! Status codes.  Success is 0; every negative value means "the caller should
  ! run the Python optimizer", which owns the user-facing messages.
  integer(i8), parameter :: CAS_ERR_INPUT       = -1_i8  ! unsupported/invalid sizes
  integer(i8), parameter :: CAS_ERR_ALLOC       = -2_i8  ! out of memory
  integer(i8), parameter :: CAS_ERR_TAGS        = -3_i8  ! handle is missing a record
  integer(i8), parameter :: CAS_ERR_CI          = -4_i8  ! fci_solve declined
  integer(i8), parameter :: CAS_ERR_EIGEN       = -5_i8  ! LAPACK failure
  integer(i8), parameter :: CAS_ERR_TRANSFORM   = -6_i8  ! AO->MO buffer allocation
  integer(i8), parameter :: CAS_ERR_ROTATE      = -7_i8  ! expm/Pade failure
  integer(i8), parameter :: CAS_ERR_DEGEN       = -8_i8  ! non-smooth SA objective
  integer(i8), parameter :: CAS_ERR_STACK       = -9_i8  ! active space too large

  ! Backtracking budget and acceptance tolerance of `_floored_newton_loop`,
  ! and the escape kick amplitudes of `_escape_saddles`.  Fixed constants in
  ! the Python too (not config keys), so they are fixed here.
  integer, parameter :: CAS_MAX_BACKTRACK = 12
  real(dp), parameter :: CAS_ACCEPT_TOL = 1.0e-12_dp
  real(dp), parameter :: CAS_ESCAPE_AMPS(6) = &
      [0.3_dp, 0.2_dp, 0.1_dp, -0.1_dp, -0.2_dp, -0.3_dp]

  ! Fixed algorithmic constants of `casscf_convergers.py` -- standard
  ! trust-region ratio thresholds and factors, the AH reference-component
  ! cutoff and the DIIS conditioning ceiling.  Deliberately not config keys
  ! there, so deliberately not options here.
  real(dp), parameter :: CAS_RHO_SHRINK   = 0.25_dp
  real(dp), parameter :: CAS_RHO_GROW     = 0.75_dp
  real(dp), parameter :: CAS_GROW_FACTOR  = 2.0_dp
  real(dp), parameter :: CAS_REJECT_FACTOR = 0.25_dp
  real(dp), parameter :: CAS_PRED_FLOOR   = 1.0e-13_dp
  real(dp), parameter :: CAS_V0_TOL       = 1.0e-10_dp
  real(dp), parameter :: CAS_DIIS_CONDMAX = 1.0e14_dp
  integer, parameter :: CAS_DIIS_BACKTRACK = 12

  !> `stats` layout returned to the caller (all int32):
  !>   0 history rows written        5 DIIS extrapolations accepted
  !>   1 macroiterations             6 DIIS extrapolations attempted
  !>   2 converged flag              7 `ah` stagnation flag
  !>   3 energy/gradient evaluations 8 `auto` outcome (0 converged, 1 stalled,
  !>   4 analytic Hessian builds       2 macroiteration cap)
  !>                                 9 `auto` macroiterations before falling back
  integer, parameter :: CAS_NSTATS = 10

  !> Resolved trust-region parameters of the `ah` converger.  Python owns the
  !> option parsing (including the "<= 0 means follow max_rotation_norm"
  !> default for the ceiling), so these arrive already resolved.
  type :: cas_ah_par_t
    real(dp) :: start_trust = 0.2_dp
    real(dp) :: max_trust = 0.5_dp
    real(dp) :: min_trust = 1.0e-6_dp
    integer :: max_micro = 32
    integer :: max_rejects = 6
    real(dp) :: saddle_curv_tol = 2.5e-2_dp
    real(dp) :: saddle_egain_tol = 1.0e-3_dp
  end type cas_ah_par_t

  !> Everything one energy/gradient evaluation needs, so `cas_evaluate` can be
  !> the plain procedure Fortran gives us in place of the Python closure.
  type :: cas_ctx_t
    integer :: n = 0            !< orbitals (== nbf)
    integer :: nc = 0           !< inactive
    integer :: na = 0           !< active
    integer :: nstate = 0       !< state-average roots
    integer :: nroot = 0        !< CI roots solved per evaluation
    integer :: npar = 0         !< non-redundant rotation pairs
    integer :: nthreads = 1
    integer(i8) :: ndet = 0
    integer :: ncall = 0        !< energy/gradient evaluations, for the trace
    integer(c_int32_t) :: iopt_ci(0:FCI_NIOPT-1) = 0_c_int32_t
    real(dp) :: dopt_ci(0:FCI_NDOPT-1) = 0.0_dp
    integer(c_int32_t), allocatable :: active(:), core(:), pairs(:)
    real(dp), allocatable :: hcore(:)          !< C-order [n,n], symmetric
    real(dp), pointer :: eri_ao(:) => null()   !< handle record, nbf^4
    real(dp), allocatable :: weights(:)
    integer, allocatable :: roots(:)           !< 0-based CI root indices
    integer(i8), allocatable :: dets(:)        !< CI ordering, not key ordering
    ! reusable work buffers -- one allocation per run, not per evaluation
    real(dp), allocatable :: h1e(:), eri(:), civ(:), enci(:), s2w(:)
    real(dp), allocatable :: g1(:), g2(:), fock(:), cvec(:)
    ! ---- fast-evaluation buffers (AO J/K + partial transforms; see
    !      casscf_fastints.F90).  The full nbf^4 MO tensor `eri` above is only
    !      written at the once-per-run call sites (final CI, canonicalization)
    !      and by the analytic-Hessian build, never per gradient evaluation.
    real(dp), allocatable :: h1e_act(:)        !< folded active h, C-order [na,na]
    real(dp), allocatable :: eri_act(:)        !< active ERIs, C-order [na^4]
    real(dp), allocatable :: eact(:)           !< (j t|u v), Fortran (n, na^3)
    real(dp), allocatable :: tp1(:), tp2(:), tp3(:)  !< transform scratch
    real(dp), allocatable :: ytr(:)            !< exchange-slice shuffle scratch
    real(dp), allocatable :: yten(:)           !< (a v|c w) slice, [na^2 n^2]
    real(dp), allocatable :: dcore(:,:), dwork(:,:)  !< AO densities
    real(dp), allocatable :: js(:,:), ks(:,:)        !< AO J/K scratch
    real(dp), allocatable :: fao(:,:), tmw(:,:)      !< AO Fock / GEMM scratch
    real(dp), allocatable :: wcore(:,:)              !< core mean field in MO
    real(dp), allocatable :: wmo(:)            !< per-state W, [nstate*n^2]
    real(dp), allocatable :: tmpna(:)          !< (na, n) scratch
    integer(c_int32_t) :: iopt_act(0:FCI_NIOPT-1) = 0_c_int32_t
    real(dp) :: dopt_act(0:FCI_NDOPT-1) = 0.0_dp
    integer(c_int32_t), allocatable :: active_act(:), core_act(:)
    logical :: warm = .false.  !< ctx%civ holds the previous point's CI vectors
    ! ---- orbital-Hessian backend (`[casscf] hessian`)
    integer :: hessmode = 0     !< 0 finite difference, 1 analytic
    integer :: nhess = 0        !< analytic Hessian builds, for the trace
    type(cas_anhess_ctx_t) :: ah
    !> CI vectors AT THE CURRENT POINT.  The analytic Hessian is a function of
    !> the orbitals and of the CI solution at them, and `ctx%civ` is clobbered
    !> by every trial evaluation (2*n_par of them on the FD path, one per
    !> backtrack and one per DIIS extrapolation elsewhere), so the accepted
    !> point's vectors are snapshotted here the moment they are known.
    real(dp), allocatable :: civ_ref(:)
    !> Integrals for the Hessian build.  `hess_fn(C, coeffs)` re-transforms at
    !> C in the Python, and must here too: after a rejected DIIS extrapolation
    !> `ctx%h1e`/`ctx%eri` belong to the *rejected* point, not to C.
    real(dp), allocatable :: h1e_h(:), eri_h(:)
  end type cas_ctx_t

  !> `[casscf] converger = trah`: CASSCF's half of the shared trust-region
  !> augmented-Hessian optimizer in `source/trah_core.F90`.
  !>
  !> The core needs four things from a method; this supplies them out of
  !> `cas_evaluate` and `cas_rotate`, so the whole converger is the trust-region
  !> machinery SCF already uses driven by CASSCF energies and gradients.
  !>
  !> The Hessian-vector product (`mode = 0`, the default) never assembles the
  !> orbital Hessian: it is the central difference of the ANALYTIC orbital
  !> gradient along the requested direction,
  !>
  !>     J.v = [ g(C exp(K(t v))) - g(C exp(-K(t v))) ] / 2t,   |t v| = fd_step
  !>
  !> i.e. exactly the object `fd_orbital_hessian` builds column by column, taken
  !> along one direction instead of along all `npar` of them.  Two CI solves per
  !> product rather than `2*npar` per Hessian, no `npar^2` matrix is ever
  !> stored, and `casscf_hess_bmat` is never called at all.  Because the CI is
  !> re-solved at each displaced point, the CI relaxation is included exactly;
  !> nothing is held fixed.
  !>
  !> What this does and does not buy (measured, macOS/ARM/Accelerate, GCC 15):
  !> against the DEFAULT finite-difference Hessian -- what `twophase` and `ah`
  !> run unless the user opts in -- it is an order of magnitude fewer CI
  !> evaluations and correspondingly faster.  Against the ASSEMBLED ANALYTIC
  !> Hessian it is still SLOWER at every size reachable here (H2O/cc-pVTZ
  !> CAS(6,6): 13.0 s matrix-free vs 6.6 s assembled), because each product is
  !> two full energy/gradient evaluations and each of those contains an
  !> O(nbf^5) AO->MO transform.  CG now asks for 4-5 products per
  !> macroiteration rather than the 14 it once did (see `cas_trah_precond` and
  !> `trah_optimize`), but the remaining budget is close to irreducible: at
  !> cc-pVTZ the run spends 72 products in CG and another 32 in the one
  !> stability check, and 104 products is 208 evaluations against the assembled
  !> arm's 29.  The matrix-free cost per macroiteration is independent of
  !> `npar` while the assembled one is `npar * O(nbf^4)`, so the crossover moves
  !> with the size of the rotation space, not with the basis alone.  Both arms
  !> are kept.
  !>
  !> One correction is needed and it is not optional.  `g` is the gradient in
  !> the DISPLACED frame -- the derivative with respect to a *further* rotation
  !> -- so `J = dg/dkappa` is the second derivative of E only up to the
  !> Baker-Campbell-Hausdorff term of `exp(K(kappa)) exp(K(delta))`.  That term
  !> is exactly antisymmetric in the two directions and of order |g|, which is
  !> why `fd_orbital_hessian` and the analytic builder both symmetrize
  !> explicitly.  Matrix-free there is no transpose to average with, so the
  !> antisymmetric part is subtracted in closed form:
  !>
  !>     H.v = J.v - 1/4 ( R[q_l,p_l] - R[p_l,q_l] ),   R = K(v) G - G K(v)
  !>
  !> with G the FULL antisymmetric gradient matrix `G[p,q] = 2(F[q,p]-F[p,q])`
  !> from the generalized Fock -- full, not just the non-redundant pairs, since
  !> the commutator leaks into the active-active block, whose gradient is
  !> non-zero for SA-CASSCF.  One nbf^3 commutator per product; pinned against
  !> the assembled Hessian in tests/test_casscf_trah.py.
  !>
  !> `mode = 1` (`[casscf] hessian = analytic`) instead assembles the exact
  !> Hessian once per macroiteration and multiplies against it.  That is the
  !> "keep the Hessian" option, kept so the two can be compared in one binary.
  type, extends(trah_provider_t) :: cas_trah_provider_t
    type(cas_ctx_t), pointer :: ctx => null()
    real(dp), allocatable :: cbuf(:,:)    !< current orbitals, C-order [n,n]
    real(dp), allocatable :: ctry(:,:)    !< trial orbitals
    real(dp), allocatable :: gcur(:)      !< gradient at the current point
    real(dp), allocatable :: gp(:), gm(:) !< displaced gradients / scratch
    real(dp), allocatable :: sv(:)        !< scaled displacement
    real(dp), allocatable :: gam(:,:)     !< full antisymmetric gradient matrix
    real(dp), allocatable :: kap(:,:), w1(:,:), w2(:,:)
    real(dp), allocatable :: hess(:,:)    !< assembled Hessian (mode 1 only)
    real(dp), allocatable :: dscr(:), fscr(:)  !< 1-RDM / effective Fock scratch
    real(dp) :: fdstep = 1.0e-4_dp
    real(dp) :: e_init = 0.0_dp           !< objective at the starting point
    real(dp) :: g_init = 0.0_dp           !< |g| at the starting point
    integer  :: mode = 0                  !< 0 matrix-free, 1 assembled
    integer  :: nmatvec = 0
    logical  :: seeded = .false.
    integer(i8) :: status = 0_i8
  contains
    procedure :: grad_hdiag   => cas_trah_grad_hdiag
    procedure :: hess_vec     => cas_trah_hess_vec
    procedure :: trial_energy => cas_trah_trial_energy
    procedure :: apply_step   => cas_trah_apply_step
  end type cas_trah_provider_t

contains

  !> C-bound entry point: one call per CASSCF run.
  !>
  !> @param[in]     c_handle  the OQP handle (`OQP::Hcore`, `OQP::AO_ERI` and
  !>                          `OQP::VEC_MO_A` are read; `OQP::VEC_MO_A` is
  !>                          overwritten with the optimized orbitals)
  !> @param[in]     iopt      integer options, CAS_NIOPT entries (schema above)
  !> @param[in]     dopt      real options, CAS_NDOPT entries (schema above)
  !> @param[in]     weights   state-average weights, [NSTATE]
  !> @param[in]     roots     0-based CI root indices averaged over, [NSTATE]
  !> @param[out]    energies  the NROOT final CI energies
  !> @param[out]    s2        <S^2> per final root, [NROOT]
  !> @param[out]    history   macroiteration table, C-order [MAXHIST, 5] as
  !>                          (it, E, dE, |g_orb|, |step|)
  !> @param[out]    stats     [rows written, macroiterations, converged,
  !>                           energy/gradient evaluations]
  !> @return        0, or a negative status meaning "run the Python optimizer"
  function casscf_energy_C(c_handle, iopt, dopt, weights, roots, energies, &
                           s2, history, stats) result(status) &
      bind(C, name="casscf_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    integer(c_int32_t), intent(in) :: iopt(0:*)
    real(dp), intent(in) :: dopt(0:*)
    real(dp), intent(in) :: weights(0:*)
    integer(c_int32_t), intent(in) :: roots(0:*)
    real(dp), intent(inout) :: energies(0:*), s2(0:*), history(0:*)
    integer(c_int32_t), intent(inout) :: stats(0:*)
    integer(i8) :: status
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    status = casscf_energy(inf, iopt, dopt, weights, roots, energies, s2, &
                           history, stats)
  end function casscf_energy_C

  !> The run itself.  Kept separate from the bind(C) shell so an in-process
  !> Fortran caller (a future geometry optimizer, say) can reach it too.
  function casscf_energy(infos, iopt, dopt, weights, roots, energies, s2, &
                         history, stats) result(status)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_Hcore, OQP_AO_ERI, &
                                   OQP_VEC_MO_A
    type(information), target, intent(inout) :: infos
    integer(c_int32_t), intent(in) :: iopt(0:*)
    real(dp), intent(in) :: dopt(0:*)
    real(dp), intent(in) :: weights(0:*)
    integer(c_int32_t), intent(in) :: roots(0:*)
    real(dp), intent(inout) :: energies(0:*), s2(0:*), history(0:*)
    integer(c_int32_t), intent(inout) :: stats(0:*)
    integer(i8) :: status

    type(cas_ctx_t), target :: ctx
    real(dp), contiguous, pointer :: hcore_p(:) => null()
    real(dp), contiguous, pointer :: eri_p(:) => null()
    real(dp), contiguous, pointer :: mo_p(:,:) => null()
    real(dp), allocatable :: cbuf(:,:), curv_w(:), curv_u(:,:)
    integer :: n, nc, na, nstate, nroot, maxmacro, optimizer, canonical
    integer :: maxescape, maxhist, ierr, k, i, j, nhist, niter, nvirt
    integer :: idx, p, q, ndet_ok, converger, hessmode
    logical :: converged, have_curv, stagnated
    real(dp) :: gradtol, enertol, steptol, maxrot, shift, fdstep
    real(dp) :: saddle_c, saddle_e, objective
    type(cas_ah_par_t) :: ahpar
    integer :: n_used, n_extrap, auto_code, auto_it

    status = 0_i8
    do k = 0, CAS_NSTATS - 1
      stats(k) = 0_c_int32_t
    end do
    n_used = 0
    n_extrap = 0
    auto_code = 0
    auto_it = 0
    stagnated = .false.

    n  = infos%basis%nbf
    nc = int(iopt(CAS_I_NCORE))
    na = int(iopt(CAS_I_NACT))
    nstate    = int(iopt(CAS_I_NSTATE))
    nroot     = int(iopt(CAS_I_NROOT))
    maxmacro  = int(iopt(CAS_I_MAXMACRO))
    optimizer = int(iopt(CAS_I_OPTIMIZER))
    canonical = int(iopt(CAS_I_CANONICAL))
    maxescape = int(iopt(CAS_I_MAXESCAPE))
    maxhist   = int(iopt(CAS_I_MAXHIST))
    converger = int(iopt(CAS_I_CONVERGER))
    hessmode  = int(iopt(CAS_I_HESSIAN))

    ahpar%max_micro       = max(1, int(iopt(CAS_I_AH_MICRO)))
    ahpar%max_rejects     = max(0, int(iopt(CAS_I_AH_REJECT)))
    ahpar%start_trust     = dopt(CAS_D_AH_START)
    ahpar%max_trust       = dopt(CAS_D_AH_MAXTR)
    ahpar%min_trust       = dopt(CAS_D_AH_MINTR)
    ahpar%saddle_curv_tol = dopt(CAS_D_AH_SADC)
    ahpar%saddle_egain_tol = dopt(CAS_D_AH_SADE)

    gradtol  = dopt(CAS_D_GRAD_TOL)
    enertol  = dopt(CAS_D_ENER_TOL)
    steptol  = dopt(CAS_D_STEP_TOL)
    maxrot   = dopt(CAS_D_MAXROT)
    shift    = dopt(CAS_D_SHIFT)
    fdstep   = dopt(CAS_D_FD_STEP)
    saddle_c = dopt(CAS_D_SADDLE_C)
    saddle_e = dopt(CAS_D_SADDLE_E)

    if (n <= 0 .or. na <= 0 .or. nc < 0 .or. nc + na > n) then
      status = CAS_ERR_INPUT
      return
    end if
    if (2 * na > FCI_MAX_NSPIN) then
      status = CAS_ERR_INPUT
      return
    end if
    if (nstate < 1 .or. nroot < 1 .or. maxhist < 2 .or. fdstep <= 0.0_dp) then
      status = CAS_ERR_INPUT
      return
    end if
    if (converger < 0 .or. converger > 4 .or. hessmode < 0 .or. hessmode > 1) then
      status = CAS_ERR_INPUT
      return
    end if
    if (converger == 1 .or. converger == 3 .or. converger == 4) then
      if (ahpar%start_trust <= 0.0_dp .or. ahpar%max_trust <= 0.0_dp .or. &
          ahpar%min_trust <= 0.0_dp) then
        status = CAS_ERR_INPUT
        return
      end if
    end if
    do k = 0, nstate - 1
      if (roots(k) < 0 .or. int(roots(k)) >= nroot) then
        status = CAS_ERR_INPUT
        return
      end if
    end do

    ! ---- handle records
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore_p)
    call tagarray_get_data(infos%dat, OQP_AO_ERI, eri_p)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_p)
    if (.not. associated(hcore_p) .or. .not. associated(eri_p) .or. &
        .not. associated(mo_p)) then
      status = CAS_ERR_TAGS
      return
    end if
    if (size(hcore_p) < n * (n + 1) / 2 .or. &
        size(eri_p) < int(n, i8)**4 .or. size(mo_p) < n * n) then
      status = CAS_ERR_TAGS
      return
    end if

    ! ---- context
    ctx%n = n
    ctx%nc = nc
    ctx%na = na
    ctx%nstate = nstate
    ctx%nroot = nroot
    ctx%nthreads = max(1, int(iopt(CAS_I_NTHREADS)))
    nvirt = n - nc - na
    ctx%npar = na * nc + nvirt * nc + nvirt * na
    if (ctx%npar <= 0) then
      status = CAS_ERR_INPUT
      return
    end if
    ctx%eri_ao => eri_p

    ! `[cas]` may in principle select a non-sequential active space, but the
    ! CASSCF optimizer itself is written around contiguous
    ! inactive/active/virtual blocks (`_nonredundant_pairs`, `_full_rdm1`), so
    ! the driver only ever sees the sequential plan; Python declines otherwise.
    allocate(ctx%active(0:na-1), ctx%core(0:max(nc, 1)-1), &
             ctx%pairs(0:2*ctx%npar-1), ctx%weights(0:nstate-1), &
             ctx%roots(0:nstate-1), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if
    do k = 0, na - 1
      ctx%active(k) = int(nc + k, c_int32_t)
    end do
    ctx%core = 0_c_int32_t
    do k = 0, nc - 1
      ctx%core(k) = int(k, c_int32_t)
    end do
    do k = 0, nstate - 1
      ctx%weights(k) = weights(k)
      ctx%roots(k) = int(roots(k))
    end do

    ! `_nonredundant_pairs`: (active, inactive), (virtual, inactive),
    ! (virtual, active), each with the first index outermost.
    idx = 0
    do p = nc, nc + na - 1
      do q = 0, nc - 1
        ctx%pairs(2*idx) = int(p, c_int32_t)
        ctx%pairs(2*idx + 1) = int(q, c_int32_t)
        idx = idx + 1
      end do
    end do
    do p = nc + na, n - 1
      do q = 0, nc - 1
        ctx%pairs(2*idx) = int(p, c_int32_t)
        ctx%pairs(2*idx + 1) = int(q, c_int32_t)
        idx = idx + 1
      end do
    end do
    do p = nc + na, n - 1
      do q = nc, nc + na - 1
        ctx%pairs(2*idx) = int(p, c_int32_t)
        ctx%pairs(2*idx + 1) = int(q, c_int32_t)
        idx = idx + 1
      end do
    end do

    ! ---- CI option block, in fci_solve's own schema
    ctx%iopt_ci = 0_c_int32_t
    ctx%iopt_ci(FCI_I_NORB)      = int(n, c_int32_t)
    ctx%iopt_ci(FCI_I_NACT)      = int(na, c_int32_t)
    ctx%iopt_ci(FCI_I_NCORE)     = int(nc, c_int32_t)
    ctx%iopt_ci(FCI_I_NALPHA)    = iopt(CAS_I_NALPHA)
    ctx%iopt_ci(FCI_I_NBETA)     = iopt(CAS_I_NBETA)
    ctx%iopt_ci(FCI_I_NROOT)     = int(nroot, c_int32_t)
    ctx%iopt_ci(FCI_I_SOLVER)    = iopt(CAS_I_SOLVER)
    ctx%iopt_ci(FCI_I_MAXITER)   = iopt(CAS_I_MAXITER)
    ctx%iopt_ci(FCI_I_SUBSPACE)  = iopt(CAS_I_SUBSPACE)
    ctx%iopt_ci(FCI_I_MULT)      = iopt(CAS_I_MULT)
    ctx%iopt_ci(FCI_I_MAXMEMORY) = iopt(CAS_I_MAXMEMORY)
    ctx%iopt_ci(FCI_I_NTHREADS)  = int(ctx%nthreads, c_int32_t)
    ctx%iopt_ci(FCI_I_WANT_S2)   = 0_c_int32_t
    ctx%dopt_ci(FCI_D_ECORE)   = dopt(CAS_D_ENUC)
    ctx%dopt_ci(FCI_D_EIG_TOL) = dopt(CAS_D_EIG_TOL)
    ctx%dopt_ci(FCI_D_CUTOFF)  = dopt(CAS_D_CUTOFF)

    ctx%ndet = binom(na, int(iopt(CAS_I_NALPHA))) * &
               binom(na, int(iopt(CAS_I_NBETA)))
    if (ctx%ndet < 1_i8 .or. int(nroot, i8) > ctx%ndet) then
      status = CAS_ERR_INPUT
      return
    end if

    call cas_alloc(ctx, ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if
    call build_determinants(na, int(iopt(CAS_I_NALPHA)), &
                            int(iopt(CAS_I_NBETA)), ctx%ndet, ctx%dets)

    ! ---- Hcore: the handle stores the lower triangle row-wise, as
    ! `_unpack_lower_triangle` reads it.  The unpacked matrix is symmetric, so
    ! its C-order buffer and its column-major view coincide.
    idx = 0
    do i = 0, n - 1
      do j = 0, i
        ctx%hcore(int(i, i8)*int(n, i8) + int(j, i8)) = hcore_p(idx + 1)
        ctx%hcore(int(j, i8)*int(n, i8) + int(i, i8)) = hcore_p(idx + 1)
        idx = idx + 1
      end do
    end do

    ! ---- MO coefficients.  `mo_p(u, p)` is AO-fastest; the kernels want the
    ! C-order buffer coeff[u][p], i.e. the transpose.  Not symmetric: this is a
    ! real transpose, never a reinterpretation.
    allocate(cbuf(0:n-1, 0:n-1), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if
    do i = 0, n - 1
      do j = 0, n - 1
        cbuf(j, i) = mo_p(i + 1, j + 1)
      end do
    end do

    ! ---- orbital-Hessian backend (`[casscf] hessian = fd | analytic`)
    ctx%hessmode = hessmode
    if (hessmode == 1) then
      allocate(ctx%civ_ref(0:ctx%ndet*int(nroot, i8) - 1_i8), &
               ctx%h1e_h(0:int(n, i8)**2 - 1_i8), &
               ctx%eri_h(0:int(n, i8)**4 - 1_i8), stat=ierr)
      if (ierr /= 0) then
        status = CAS_ERR_ALLOC
        return
      end if
      status = cas_anhess_init(ctx%ah, n, nc, na, int(iopt(CAS_I_NALPHA)), &
                               int(iopt(CAS_I_NBETA)), ctx%npar, nstate, &
                               ctx%nthreads, ctx%ndet, ctx%dets)
      if (status < 0_i8) then
        ! The 2 GiB excitation-stack guard and the encoding limit are refusals,
        ! not fallbacks: hand them back so Python raises its own message.
        if (status == CAS_HESS_ERR_STACK) status = CAS_ERR_STACK
        return
      end if
    end if

    ! ============================================== orbital-optimization phase
    nhist = 0
    select case (converger)
    case (1)          ! ah: trust-region augmented Hessian + curvature escape
      call ah_converge(ctx, cbuf, maxmacro, gradtol, enertol, steptol, &
                       fdstep, ahpar, 0, maxescape, history, maxhist, nhist, &
                       converged, niter, objective, stagnated, status)
      if (status < 0_i8) return

    case (2)          ! diis: Pulay acceleration over the production step
      call diis_optimize(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                         steptol, maxrot, shift, fdstep, saddle_c, saddle_e, &
                         maxescape, int(iopt(CAS_I_DIIS_SPACE)), &
                         int(iopt(CAS_I_DIIS_START)), history, maxhist, nhist, &
                         converged, niter, objective, n_used, n_extrap, status)
      if (status < 0_i8) return

    case (3)          ! auto: ah with a two-phase fallback on stagnation
      call auto_optimize(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                         steptol, maxrot, shift, fdstep, saddle_c, saddle_e, &
                         maxescape, ahpar, max(1, int(iopt(CAS_I_AUTO_STAG))), &
                         history, maxhist, nhist, converged, niter, objective, &
                         auto_code, auto_it, stagnated, status)
      if (status < 0_i8) return

    case (4)          ! trah: matrix-free trust-region augmented Hessian
      call trah_optimize(ctx, cbuf, maxmacro, gradtol, fdstep, ahpar, hessmode, &
                         history, maxhist, nhist, converged, niter, objective, &
                         status)
      if (status < 0_i8) return

    case default      ! twophase: damped Newton, then the saddle escape
      call newton_loop(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                       steptol, maxrot, shift, fdstep, history, maxhist, &
                       nhist, .true., converged, niter, curv_w, curv_u, &
                       have_curv, objective, status)
      if (status < 0_i8) return

      if (converged .and. optimizer == 0 .and. maxescape > 0 .and. have_curv) then
        call escape_saddles(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                            steptol, maxrot, shift, fdstep, saddle_c, saddle_e, &
                            maxescape, history, maxhist, nhist, niter, &
                            curv_w, curv_u, objective, status)
        if (status < 0_i8) return
        converged = .true.
      end if
    end select

    ! ==================================================== canonicalization
    if (canonical /= 0) then
      call canonicalize(ctx, cbuf, status)
      if (status < 0_i8) return
    end if

    ! ---- commit the optimized orbitals (the Fortran core reads this buffer)
    do i = 0, n - 1
      do j = 0, n - 1
        mo_p(i + 1, j + 1) = cbuf(j, i)
      end do
    end do

    ! ================================ final CI on the optimized orbitals
    call mo_transform_h1e(int(n, c_int32_t), ctx%hcore, cbuf, ctx%h1e)
    if (mo_transform_eri(int(n, c_int32_t), ctx%eri_ao, cbuf, ctx%eri) /= 0_i8) then
      status = CAS_ERR_TRANSFORM
      return
    end if
    ctx%iopt_ci(FCI_I_NROOT) = int(nroot, c_int32_t)
    ctx%iopt_ci(FCI_I_WANT_S2) = 1_c_int32_t
    ndet_ok = int(fci_solve(ctx%iopt_ci, ctx%dopt_ci, ctx%active, ctx%core, &
                            ctx%h1e, ctx%eri, ctx%enci, ctx%civ, ctx%s2w))
    ctx%iopt_ci(FCI_I_WANT_S2) = 0_c_int32_t
    if (ndet_ok < 0) then
      status = CAS_ERR_CI
      return
    end if
    do k = 0, nroot - 1
      energies(k) = ctx%enci(k)
      s2(k) = ctx%s2w(k)
    end do

    stats(0) = int(nhist, c_int32_t)
    stats(1) = int(niter, c_int32_t)
    stats(2) = merge(1_c_int32_t, 0_c_int32_t, converged)
    stats(3) = int(ctx%ncall, c_int32_t)
    stats(4) = int(ctx%nhess, c_int32_t)
    stats(5) = int(n_used, c_int32_t)
    stats(6) = int(n_extrap, c_int32_t)
    stats(7) = merge(1_c_int32_t, 0_c_int32_t, stagnated)
    stats(8) = int(auto_code, c_int32_t)
    stats(9) = int(auto_it, c_int32_t)
    ! `ctx` is a local, so its allocatable components -- including the
    ! excitation stack inside `ctx%ah` -- are released on return, on the error
    ! paths as well as this one.
  end function casscf_energy

  ! =================================================================== work
  subroutine cas_alloc(ctx, ierr)
    type(cas_ctx_t), intent(inout) :: ctx
    integer, intent(out) :: ierr
    integer(i8) :: n2, n4, na2, na4
    integer(i8) :: n_i8, na_i8
    integer :: k

    n2 = int(ctx%n, i8)**2
    n4 = int(ctx%n, i8)**4
    na2 = int(ctx%na, i8)**2
    na4 = int(ctx%na, i8)**4
    n_i8 = int(ctx%n, i8)
    na_i8 = int(ctx%na, i8)
    allocate(ctx%hcore(0:n2-1), ctx%h1e(0:n2-1), ctx%eri(0:n4-1), &
             ctx%civ(0:ctx%ndet*int(ctx%nroot, i8)-1), &
             ctx%enci(0:ctx%nroot-1), ctx%s2w(0:ctx%nroot-1), &
             ctx%g1(0:int(ctx%nstate, i8)*na2-1), &
             ctx%g2(0:int(ctx%nstate, i8)*na4-1), &
             ctx%fock(0:n2-1), ctx%cvec(0:ctx%ndet-1), &
             ctx%dets(ctx%ndet), stat=ierr)
    if (ierr /= 0) return
    allocate(ctx%h1e_act(0:na2-1), ctx%eri_act(0:na4-1), &
             ctx%eact(0:n_i8*na_i8**3-1), &
             ctx%tp1(0:n_i8**3*na_i8-1), ctx%tp2(0:n_i8**2*na_i8**2-1), &
             ctx%tp3(0:n_i8*na_i8**3-1), &
             ctx%ytr(0:n_i8**3*na_i8-1), ctx%yten(0:n_i8**2*na_i8**2-1), &
             ctx%dcore(0:ctx%n-1, 0:ctx%n-1), ctx%dwork(0:ctx%n-1, 0:ctx%n-1), &
             ctx%js(0:ctx%n-1, 0:ctx%n-1), ctx%ks(0:ctx%n-1, 0:ctx%n-1), &
             ctx%fao(0:ctx%n-1, 0:ctx%n-1), ctx%tmw(0:ctx%n-1, 0:ctx%n-1), &
             ctx%wcore(0:ctx%n-1, 0:ctx%n-1), &
             ctx%wmo(0:int(ctx%nstate, i8)*n2-1), &
             ctx%tmpna(0:na_i8*n_i8-1), &
             ctx%active_act(0:max(ctx%na, 1)-1), ctx%core_act(0:0), stat=ierr)
    if (ierr /= 0) return
    do k = 0, ctx%na - 1
      ctx%active_act(k) = int(k, c_int32_t)
    end do
    ctx%core_act(0) = 0_c_int32_t
  end subroutine cas_alloc

  !> One energy/gradient evaluation at the orbitals `cbuf`.
  !>
  !> Reproduces `casscf.py::_optimize`'s `evaluate` closure -- one CI solve, the
  !> per-root spatial RDMs, then the weighted generalized Fock and orbital
  !> gradient -- but never forms the full nbf^4 MO ERI tensor.  The CI runs on
  !> the frozen-core-folded ACTIVE problem (h from the AO-basis inactive Fock,
  !> ERIs from the shared partial transformation) and the Fock/gradient build
  !> gets its per-root mean field W from AO-basis J/K contractions plus a
  !> one-index back-transform.  Same numbers to rounding, at
  !> O(nbf^4 nact) + O(nbf^4) per evaluation instead of O(nbf^5) -- and the
  !> TRAH converger calls this twice per Hessian-vector product, which is what
  !> made the full transform the dominant cost of every CASSCF run.
  !>
  !> After the first solve the CI restarts from the previous point's vectors
  !> (`ctx%warm`): inside one orbital optimization consecutive CI problems
  !> differ by O(step), so a warm Davidson converges in a handful of
  !> applications where a cold start pays the full price.  The warm start only
  !> replaces the `auto` solver choice -- an explicit `[ci] solver=dense` is
  !> honored -- and only when no spin filter is active (the filtered window
  !> logic wants its own root growth).
  subroutine cas_evaluate(ctx, cbuf, objective, grad, ierr)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(in) :: cbuf(0:,0:)
    real(dp), intent(out) :: objective
    real(dp), contiguous, intent(out) :: grad(0:)
    integer(i8), intent(out) :: ierr

    integer :: k, r, na, n, nc, t, u
    integer(i8) :: j, cap, rc, off1, off2, a, b
    real(dp) :: efold

    ierr = 0_i8
    ctx%ncall = ctx%ncall + 1
    na = ctx%na
    n = ctx%n
    nc = ctx%nc

    ! ---- inactive Fock in AO, folded core energy, active h
    call cas_ao_density(n, nc, 0, cbuf, ctx%g1, ctx%dcore, ctx%tmpna)
    call cas_wmo_from_ao(n, ctx%eri_ao, ctx%hcore, cbuf, ctx%dcore, &
                         ctx%js, ctx%ks, ctx%fao, ctx%tmw, ctx%wcore)
    efold = 0.0_dp
    do b = 0_i8, int(n, i8) - 1_i8
      do a = 0_i8, int(n, i8) - 1_i8
        efold = efold + ctx%dcore(a, b) &
                      * (ctx%hcore(a * int(n, i8) + b) + ctx%fao(a, b))
      end do
    end do
    efold = 0.5_dp * efold
    do t = 0, na - 1
      do u = 0, na - 1
        ctx%h1e_act(int(t, i8) * int(na, i8) + int(u, i8)) = &
            ctx%wcore(nc + t, nc + u)
      end do
    end do

    ! ---- active ERIs and the (j t|u v) slice, one shared partial transform
    call cas_partial_eri(n, na, nc, ctx%eri_ao, cbuf, ctx%tp1, ctx%tp2, &
                         ctx%tp3, ctx%eri_act, ctx%eact)

    ! ---- CI on the folded active problem
    ctx%iopt_act = ctx%iopt_ci
    ctx%iopt_act(FCI_I_NORB) = int(na, c_int32_t)
    ctx%iopt_act(FCI_I_NACT) = int(na, c_int32_t)
    ctx%iopt_act(FCI_I_NCORE) = 0_c_int32_t
    ctx%iopt_act(FCI_I_NROOT) = int(ctx%nroot, c_int32_t)
    if (ctx%warm .and. ctx%iopt_ci(FCI_I_SOLVER) == 0_c_int32_t .and. &
        ctx%iopt_ci(FCI_I_MULT) == 0_c_int32_t) then
      ctx%iopt_act(FCI_I_SOLVER) = 2_c_int32_t
      ctx%iopt_act(FCI_I_GUESS) = 1_c_int32_t
    end if
    ctx%dopt_act = ctx%dopt_ci
    ctx%dopt_act(FCI_D_ECORE) = ctx%dopt_ci(FCI_D_ECORE) + efold
    if (fci_solve(ctx%iopt_act, ctx%dopt_act, ctx%active_act, ctx%core_act, &
                  ctx%h1e_act, ctx%eri_act, ctx%enci, ctx%civ, &
                  ctx%s2w) < 0_i8) then
      ierr = CAS_ERR_CI
      return
    end if
    ctx%warm = .true.

    cap = ctx%ndet * int(2 * na, i8)**2 + 1_i8
    do k = 0, ctx%nstate - 1
      r = ctx%roots(k)
      ! `civ` is C-order [ndet, nroot]; root r is the strided column r.
      do j = 0_i8, ctx%ndet - 1_i8
        ctx%cvec(j) = ctx%civ(j * int(ctx%nroot, i8) + int(r, i8))
      end do
      off1 = int(k, i8) * int(na, i8)**2
      off2 = int(k, i8) * int(na, i8)**4
      ! string-factorized RDMs (one DGEMM); the walking kernels remain the
      ! numerical pin and the fallback for non-product determinant lists
      if (rdm12_strings(na, ctx%ndet, ctx%dets, ctx%cvec, ctx%g1(off1), &
                        ctx%g2(off2), ctx%nthreads) /= 0) then
        call rdm1_spatial(int(na, c_int32_t), ctx%ndet, ctx%dets, ctx%cvec, &
                          ctx%g1(off1))
        rc = rdm2_spatial(int(na, c_int32_t), ctx%ndet, ctx%dets, ctx%cvec, &
                          cap, ctx%g2(off2), 0_c_int32_t)
        if (rc /= 0_i8) then
          ierr = CAS_ERR_ALLOC
          return
        end if
      end if
    end do

    ! ---- per-state mean field W = C^T [F_core + J(D_act) - K(D_act)/2] C.
    ! The active-density J/K come from the partial-transform slices (t2 and
    ! the exchange-shaped yten), so nothing past the core build and the shared
    ! pass-1 transform ever streams the nbf^4 tensor again -- per-state cost
    ! is O(nbf^2 nact^2) regardless of how many states are averaged.
    call cas_act_exchange_tensor(n, na, nc, cbuf, ctx%tp1, ctx%ytr, ctx%yten)
    do k = 0, ctx%nstate - 1
      off1 = int(k, i8) * int(na, i8)**2
      call cas_act_fock(n, na, ctx%tp2, ctx%yten, ctx%g1(off1), ctx%fao, &
                        ctx%js, ctx%ks, ctx%dwork)
      ! W_k = C^T F_k C
      call dgemm('N', 'N', n, n, n, 1.0_dp, cbuf(0, 0), n, ctx%dwork(0, 0), &
                 n, 0.0_dp, ctx%tmw(0, 0), n)
      call dgemm('N', 'T', n, n, n, 1.0_dp, ctx%tmw(0, 0), n, cbuf(0, 0), n, &
                 0.0_dp, ctx%wmo(int(k, i8) * int(n, i8)**2), n)
    end do
    call casscf_gfock_grad_w(n, nc, na, ctx%nstate, ctx%weights, ctx%g1, &
                             ctx%g2, ctx%wmo, ctx%eact, ctx%fock, grad)

    objective = 0.0_dp
    do k = 0, ctx%nstate - 1
      objective = objective + ctx%weights(k) * ctx%enci(ctx%roots(k))
    end do
  end subroutine cas_evaluate

  !> C <- C exp(K(vec)) over the non-redundant pairs.
  subroutine cas_rotate(ctx, cin, vec, cout, ierr)
    type(cas_ctx_t), intent(in) :: ctx
    real(dp), contiguous, intent(in) :: cin(0:,0:), vec(0:)
    real(dp), contiguous, intent(inout) :: cout(0:,0:)
    integer(i8), intent(out) :: ierr
    integer(c_int32_t) :: info

    info = casscf_orbital_rotate(int(ctx%n, c_int32_t), &
                                 int(ctx%npar, c_int32_t), ctx%pairs, vec, &
                                 cin, cout)
    ierr = merge(0_i8, CAS_ERR_ROTATE, info == 0_c_int32_t)
  end subroutine cas_rotate

  ! ==================================================================== trah
  !> `[casscf] converger = trah`: run the shared TRAH core over CASSCF.
  !>
  !> The history table is filled the way `newton_loop` fills it -- a seed row at
  !> iteration 0 and one row per accepted macroiteration -- so `_write_log`
  !> formats it unchanged.
  subroutine trah_optimize(ctx, cbuf, maxmacro, gradtol, fdstep, par_ah, &
                           hessmode, history, maxhist, nhist, converged, niter, &
                           objective, status)
    type(cas_ctx_t), target, intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, maxhist, hessmode
    real(dp), intent(in) :: gradtol, fdstep
    type(cas_ah_par_t), intent(in) :: par_ah
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist
    logical, intent(out) :: converged
    integer, intent(out) :: niter
    real(dp), intent(out) :: objective
    integer(i8), intent(inout) :: status

    type(cas_trah_provider_t) :: prov
    type(trah_params_t) :: par
    type(trah_result_t) :: tres
    integer, allocatable :: hit(:)
    real(dp), allocatable :: he(:), hde(:), hg(:), hs(:)
    integer :: n, npar, ierr, nh, r

    n = ctx%n
    npar = ctx%npar
    converged = .false.
    niter = 0
    objective = 0.0_dp

    prov%nparam = npar
    prov%ctx => ctx
    prov%mode = hessmode
    prov%fdstep = fdstep
    allocate(prov%cbuf(0:n-1, 0:n-1), prov%ctry(0:n-1, 0:n-1), &
             prov%gcur(0:npar-1), prov%gp(0:npar-1), prov%gm(0:npar-1), &
             prov%sv(0:npar-1), prov%gam(0:n-1, 0:n-1), prov%kap(0:n-1, 0:n-1), &
             prov%w1(0:n-1, 0:n-1), prov%w2(0:n-1, 0:n-1), &
             prov%dscr(0:int(n, i8)**2 - 1_i8), &
             prov%fscr(0:int(n, i8)**2 - 1_i8), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if
    if (hessmode == 1) then
      allocate(prov%hess(npar, npar), stat=ierr)
      if (ierr /= 0) then
        status = CAS_ERR_ALLOC
        return
      end if
    end if
    prov%cbuf = cbuf

    allocate(hit(maxhist), he(maxhist), hde(maxhist), hg(maxhist), hs(maxhist), &
             stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    ! The micro-solver is the Steihaug-Toint truncated CG rather than the
    ! augmented-Hessian Davidson the SCF path defaults to.  Both live in the
    ! core; the choice follows the cost of one Hessian-vector product.  For SCF
    ! that product is one Fock-like contraction, so the Davidson can afford ~40
    ! of them per macroiteration; here it is two CI solves, and CG's relative
    ! residual test stops after a handful.  Negative curvature is still handled
    ! head-on -- CG walks to the trust boundary along it -- and the core's
    ! stability check (a Davidson for the lowest Hessian eigenpair) still runs
    ! at apparent convergence, which is what catches the spurious CASSCF saddles
    ! the two-phase converger needs a whole escape phase for.
    par%nmac = maxmacro
    par%nmic = par_ah%max_micro
    par%conv_tol = gradtol
    par%r0 = par_ah%start_trust
    ! Trust-radius ceiling: the core default max(4, 8*r0), NOT
    ! `max_rotation_norm`.  `max_rotation_norm` is a FIXED cap designed for the
    ! level-shifted Newton step of the two-phase and `ah` convergers, which have
    ! no other way to keep a step sane; a trust region does that job adaptively
    ! from the actual/predicted energy ratio, and clamping it to 0.2 costs
    ! roughly twice the macroiterations for nothing -- measured 12 -> 30 on
    ! H2O/cc-pVDZ CAS(6,6) and 14 -> 28 on stretched H2O/6-31G CAS(6,6), at the
    ! same energy.  Set `ah_start_trust_radius` to control it.
    par%dmax = 0.0_dp
    par%sub_solver = 2
    par%nrtv = 1
    par%deterministic = .false.
    par%rms_gnorm = .false.        ! [casscf] gradient_norm_tol is a 2-norm
    par%verbose = .false.          ! Python owns the CASSCF log
    par%want_history = .true.
    ! Two settings the SCF path deliberately does NOT take (its Hessian-vector
    ! product is one Fock-like contraction; this one is two CI solves, so the
    ! trade is the other way round here).  Both are measured in the comment on
    ! `cas_trah_precond`.
    par%cg_res_floor = gradtol     ! never solve tighter than |g| has to be
    par%stab_seed_diag = .true.    ! start the stability check at the soft modes

    call trah_run(prov, par, tres, hit, he, hde, hg, hs, nh)
    if (tres%ierr < 0) then
      status = prov%status
      if (status >= 0_i8) status = CAS_ERR_INPUT
      return
    end if

    cbuf = prov%cbuf
    objective = tres%energy
    converged = tres%converged
    niter = tres%iter

    call push_history(history, maxhist, nhist, 0, prov%e_init, 0.0_dp, &
                      prov%g_init, 0.0_dp)
    do r = 1, nh
      call push_history(history, maxhist, nhist, hit(r), he(r), hde(r), hg(r), hs(r))
    end do
  end subroutine trah_optimize

  !> Gradient, preconditioner diagonal and objective at the provider's orbitals.
  subroutine cas_trah_grad_hdiag(this, g, hdiag, e, ierr)
    class(cas_trah_provider_t), intent(inout) :: this
    real(dp), intent(out) :: g(:), hdiag(:), e
    integer,  intent(out) :: ierr
    integer(i8) :: rc
    integer :: n, p, q

    ierr = 0
    n = this%ctx%n
    call cas_evaluate(this%ctx, this%cbuf, e, g, rc)
    if (rc < 0_i8) then
      this%status = rc
      ierr = -1
      return
    end if
    call cas_snap_ci(this%ctx)
    this%gcur = g
    if (.not. this%seeded) then
      this%e_init = e
      this%g_init = sqrt(sum(g*g))
      this%seeded = .true.
    end if

    ! Full antisymmetric gradient matrix from the generalized Fock, in the same
    ! convention `casscf_gfock_grad` uses for the non-redundant components:
    ! G[p,q] = 2 (F[q,p] - F[p,q]).  The active-active block is NOT redundant
    ! for a state-averaged objective, and the BCH correction in `hess_vec`
    ! contracts against it, so the whole matrix is formed.
    do p = 0, n - 1
      do q = 0, n - 1
        this%gam(p, q) = 2.0_dp * (this%ctx%fock(int(q, i8)*int(n, i8) + int(p, i8)) &
                                 - this%ctx%fock(int(p, i8)*int(n, i8) + int(q, i8)))
      end do
    end do

    call cas_trah_precond(this, hdiag)

    if (this%mode == 1) then
      call build_hessian(this%ctx, this%cbuf, this%fdstep, this%hess, this%status)
      if (this%status < 0_i8) ierr = -1
    end if
  end subroutine cas_trah_grad_hdiag

  !> Preconditioner diagonal: the orbital-energy-difference approximation to
  !> the orbital Hessian diagonal,
  !>
  !>     hdiag[p,q] ~ 2 (n_q - n_p) (eps_p - eps_q),
  !>
  !> with `eps` the diagonal of the closed+active mean-field Fock
  !> (`casscf_effective_fock`, the same matrix canonicalization uses) and `n`
  !> the state-averaged occupations.  For a virtual/inactive rotation this is
  !> the familiar 4(eps_a - eps_i).
  !>
  !> This is a PRECONDITIONER, not a Hessian: it only sets the CG's metric, so
  !> an approximation costs iterations and never accuracy.  Building the exact
  !> diagonal would cost an assembled Hessian, which is the whole thing this
  !> converger exists to avoid.
  !>
  !> Regularization, and why it is not optional
  !> ------------------------------------------
  !> The approximation is worst exactly where it is used hardest.  For an
  !> active/virtual rotation the occupation factor is `n_t`, so a weakly
  !> occupied active orbital (n_t ~ 3e-3 in H2O/cc-pVDZ CAS(6,6)) produces a
  !> diagonal of ~5e-3 while the true Hessian diagonal there is ~1e-3 and
  !> sometimes NEGATIVE -- measured against the assembled analytic Hessian, the
  !> ratio reaches 13.  Those are 96 of the 140 rotations at cc-pVDZ and 300 of
  !> 412 at cc-pVTZ, and dividing the residual by a near-zero denominator is
  !> what makes CG grind.  A Levenberg-style shift `sigma = 1e-4 max(hdiag)`
  !> damps them; it is relative, so it carries across systems, and it changes
  !> only the metric, never the fixed point.
  !>
  !> Measured together with the CG residual floor of `trah_optimize`, on an
  !> idle Xeon 8368 (Linux/OpenBLAS/GCC 11), best of 3 alternating runs, CG
  !> products per macroiteration and CI evaluations before -> after:
  !>
  !>   H2O/cc-pVDZ CAS(6,6)        10.6 -> 4.4   382 -> 176   3.93 -> 1.70 s
  !>   H2O/cc-pVTZ CAS(6,6)        14.3 -> 5.1   488 -> 237  26.44 -> 13.01 s
  !>   H2O 1.5x /6-31G CAS(6,6)    11.6 -> 6.6   468 -> 244   7.30 -> 3.29 s
  !>
  !> at converged energies identical to the last printed digit in all three.
  !>
  !> The size of the shift matters and is not a free parameter: at
  !> `5e-4 max(hdiag)` and above the run converges to a DIFFERENT, higher
  !> CASSCF solution (-76.0833 instead of -76.1134 on cc-pVDZ), and at 2e-4 it
  !> still finds the right one but takes more macroiterations than it saves.
  !> What did NOT work, and was measured rather than assumed: the exact
  !> diagonal frozen from one analytic Hessian build at the starting geometry
  !> is WORSE than this approximation (135 vs 120 CG products) -- the diagonal
  !> has to track the orbitals, a stale exact one is not an improvement.
  subroutine cas_trah_precond(this, hdiag)
    class(cas_trah_provider_t), intent(inout) :: this
    real(dp), intent(out) :: hdiag(:)
    integer :: n, nc, na, i, j, k, l, p, q
    real(dp) :: val, shift

    n = this%ctx%n
    nc = this%ctx%nc
    na = this%ctx%na

    this%dscr = 0.0_dp
    do i = 0, nc - 1
      this%dscr(int(i, i8)*int(n, i8) + int(i, i8)) = 2.0_dp
    end do
    do k = 0, this%ctx%nstate - 1
      do i = 0, na - 1
        do j = 0, na - 1
          this%dscr(int(nc + i, i8)*int(n, i8) + int(nc + j, i8)) = &
              this%dscr(int(nc + i, i8)*int(n, i8) + int(nc + j, i8)) &
              + this%ctx%weights(k) &
                * this%ctx%g1(int(k, i8)*int(na, i8)**2 + int(i, i8)*int(na, i8) &
                              + int(j, i8))
        end do
      end do
    end do

    ! The mean field of `dscr` used to come from `casscf_effective_fock` over
    ! `ctx%h1e`/`ctx%eri` -- which cas_evaluate no longer fills (and which,
    ! after a rejected extrapolation, belonged to the wrong point anyway).
    ! Build it from the AO integrals at the CURRENT orbitals instead: same
    ! matrix, no nbf^4 MO tensor, no staleness.
    call cas_dmo_to_ao(n, this%cbuf, this%dscr, this%ctx%tmw, this%ctx%dwork)
    call cas_wmo_from_ao(n, this%ctx%eri_ao, this%ctx%hcore, this%cbuf, &
                         this%ctx%dwork, this%ctx%js, this%ctx%ks, &
                         this%ctx%fao, this%ctx%tmw, this%fscr)

    do l = 0, this%ctx%npar - 1
      p = int(this%ctx%pairs(2*l))
      q = int(this%ctx%pairs(2*l + 1))
      val = 2.0_dp &
          * (this%dscr(int(q, i8)*int(n, i8) + int(q, i8)) &
           - this%dscr(int(p, i8)*int(n, i8) + int(p, i8))) &
          * (this%fscr(int(p, i8)*int(n, i8) + int(p, i8)) &
           - this%fscr(int(q, i8)*int(n, i8) + int(q, i8)))
      hdiag(l + 1) = max(val, 1.0e-3_dp)
    end do

    shift = 1.0e-4_dp * maxval(hdiag(1:this%ctx%npar))
    hdiag(1:this%ctx%npar) = hdiag(1:this%ctx%npar) + shift
  end subroutine cas_trah_precond

  !> hv = H.v.  See the type's documentation for the derivation of the BCH
  !> correction that makes the matrix-free product symmetric.
  subroutine cas_trah_hess_vec(this, v, hv, ierr)
    class(cas_trah_provider_t), intent(inout) :: this
    real(dp), intent(in)  :: v(:)
    real(dp), intent(out) :: hv(:)
    integer,  intent(out) :: ierr
    integer(i8) :: rc
    integer :: n, npar, l, p, q
    real(dp) :: nv, t, obj

    ierr = 0
    n = this%ctx%n
    npar = this%ctx%npar

    if (this%mode == 1) then
      ! the assembled Hessian is symmetric, so its C-order buffer and its
      ! column-major view are the same matrix
      call dgemv('N', npar, npar, 1.0_dp, this%hess, npar, v, 1, 0.0_dp, hv, 1)
      return
    end if

    nv = sqrt(sum(v*v))
    if (nv <= 0.0_dp) then
      hv = 0.0_dp
      return
    end if
    t = this%fdstep / nv
    this%sv = t * v

    call cas_rotate(this%ctx, this%cbuf, this%sv, this%ctry, rc)
    if (rc >= 0_i8) call cas_evaluate(this%ctx, this%ctry, obj, this%gp, rc)
    if (rc < 0_i8) then
      this%status = rc
      ierr = -1
      return
    end if
    this%sv = -this%sv
    call cas_rotate(this%ctx, this%cbuf, this%sv, this%ctry, rc)
    if (rc >= 0_i8) call cas_evaluate(this%ctx, this%ctry, obj, this%gm, rc)
    if (rc < 0_i8) then
      this%status = rc
      ierr = -1
      return
    end if
    this%nmatvec = this%nmatvec + 1

    do l = 0, npar - 1
      hv(l + 1) = (this%gp(l) - this%gm(l)) / (2.0_dp * t)
    end do

    ! R = K(v) G - G K(v); subtract the antisymmetric BCH part.
    this%kap = 0.0_dp
    do l = 0, npar - 1
      p = int(this%ctx%pairs(2*l))
      q = int(this%ctx%pairs(2*l + 1))
      this%kap(p, q) = this%kap(p, q) + v(l + 1)
      this%kap(q, p) = this%kap(q, p) - v(l + 1)
    end do
    call dgemm('N', 'N', n, n, n, 1.0_dp, this%kap, n, this%gam, n, 0.0_dp, this%w1, n)
    call dgemm('N', 'N', n, n, n, 1.0_dp, this%gam, n, this%kap, n, 0.0_dp, this%w2, n)
    do l = 0, npar - 1
      p = int(this%ctx%pairs(2*l))
      q = int(this%ctx%pairs(2*l + 1))
      hv(l + 1) = hv(l + 1) - 0.25_dp * ((this%w1(q, p) - this%w2(q, p)) &
                                       - (this%w1(p, q) - this%w2(p, q)))
    end do
  end subroutine cas_trah_hess_vec

  !> Objective at (current orbitals rotated by p); the current point is unchanged.
  subroutine cas_trah_trial_energy(this, p, e, ierr)
    class(cas_trah_provider_t), intent(inout) :: this
    real(dp), intent(in)  :: p(:)
    real(dp), intent(out) :: e
    integer,  intent(out) :: ierr
    integer(i8) :: rc
    ierr = 0
    this%sv = p
    call cas_rotate(this%ctx, this%cbuf, this%sv, this%ctry, rc)
    if (rc >= 0_i8) call cas_evaluate(this%ctx, this%ctry, e, this%gp, rc)
    if (rc < 0_i8) then
      this%status = rc
      ierr = -1
    end if
  end subroutine cas_trah_trial_energy

  !> Commit the rotation.
  subroutine cas_trah_apply_step(this, p, ierr)
    class(cas_trah_provider_t), intent(inout) :: this
    real(dp), intent(in) :: p(:)
    integer,  intent(out) :: ierr
    integer(i8) :: rc
    ierr = 0
    this%sv = p
    call cas_rotate(this%ctx, this%cbuf, this%sv, this%ctry, rc)
    if (rc < 0_i8) then
      this%status = rc
      ierr = -1
      return
    end if
    this%cbuf = this%ctry
  end subroutine cas_trah_apply_step

  ! ============================================================ Newton loop
  !> `casscf.py::_floored_newton_loop`, including the finite-difference
  !> Hessian, the level-shifted Newton step and the halving backtrack.
  !>
  !> `curv_w`/`curv_u` come back as the RAW (unfloored) eigendecomposition of
  !> the Hessian built by the FINAL step, so the saddle-escape phase costs zero
  !> extra CI solves -- the same contract the Python has.
  subroutine newton_loop(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                         steptol, maxrot, shift, fdstep, history, maxhist, &
                         nhist, seed_history, converged, niter, curv_w, &
                         curv_u, have_curv, objective, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, optimizer, maxhist
    real(dp), intent(in) :: gradtol, enertol, steptol, maxrot, shift, fdstep
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist
    logical, intent(in) :: seed_history
    logical, intent(out) :: converged, have_curv
    integer, intent(out) :: niter
    real(dp), allocatable, intent(inout) :: curv_w(:), curv_u(:,:)
    real(dp), intent(out) :: objective
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: grad(:), grad_new(:), step(:), ctry(:,:)
    integer :: npar, it, bt, ierr
    integer(i8) :: rc
    real(dp) :: obj_old, obj_new, gnorm, step_norm
    logical :: first_row

    npar = ctx%npar
    converged = .false.
    have_curv = .false.
    niter = 0
    allocate(grad(0:npar-1), grad_new(0:npar-1), step(0:npar-1), &
             ctry(0:ctx%n-1, 0:ctx%n-1), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    call cas_evaluate(ctx, cbuf, objective, grad, rc)
    if (rc < 0_i8) then
      status = rc
      return
    end if
    call cas_snap_ci(ctx)
    first_row = seed_history
    if (first_row) call push_history(history, maxhist, nhist, 0, objective, &
                                     0.0_dp, vnorm(grad), 0.0_dp)

    do it = 1, maxmacro
      niter = it
      gnorm = vnorm(grad)
      if (gnorm < gradtol) then
        converged = .true.
        exit
      end if

      if (optimizer == 0) then
        call newton_step(ctx, cbuf, grad, fdstep, shift, maxrot, step, &
                         curv_w, curv_u, status)
        if (status < 0_i8) return
        have_curv = .true.
      else
        ! `_powell_step`: scaled steepest descent, no curvature.
        step = -grad
        step_norm = vnorm(step)
        if (step_norm > maxrot) step = step * (maxrot / step_norm)
      end if

      obj_old = objective
      do bt = 1, CAS_MAX_BACKTRACK
        call cas_rotate(ctx, cbuf, step, ctry, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        call cas_evaluate(ctx, ctry, obj_new, grad_new, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        if (obj_new <= obj_old + CAS_ACCEPT_TOL) exit
        step = 0.5_dp * step
      end do

      ! ctx%civ still belongs to the last trial evaluation, which IS the point
      ! being accepted; snapshot before anything else can clobber it.
      call cas_snap_ci(ctx)
      cbuf = ctry
      objective = obj_new
      grad = grad_new
      step_norm = vnorm(step)
      call push_history(history, maxhist, nhist, it, objective, &
                        objective - obj_old, vnorm(grad), step_norm)
      if (abs(objective - obj_old) < enertol .and. step_norm < steptol) then
        converged = vnorm(grad) < gradtol
        exit
      end if
    end do
  end subroutine newton_loop

  !> `casscf.py::_newton_step`, over whichever Hessian `[casscf] hessian`
  !> selects (`build_hessian` dispatches).
  subroutine newton_step(ctx, cbuf, grad, fdstep, shift, maxrot, step, &
                         curv_w, curv_u, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(in) :: cbuf(0:,0:), grad(0:)
    real(dp), intent(in) :: fdstep, shift, maxrot
    real(dp), contiguous, intent(out) :: step(0:)
    real(dp), allocatable, intent(inout) :: curv_w(:), curv_u(:,:)
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: hess(:,:), ge(:), wf(:)
    integer :: npar, k, i, ierr
    real(dp) :: nrm

    npar = ctx%npar
    if (allocated(curv_w)) deallocate(curv_w)
    if (allocated(curv_u)) deallocate(curv_u)
    allocate(hess(npar, npar), curv_w(npar), curv_u(npar, npar), &
             ge(npar), wf(npar), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    call build_hessian(ctx, cbuf, fdstep, hess, status)
    if (status < 0_i8) return

    ! The Hessian is symmetric by construction, so its C-order buffer and its
    ! column-major view coincide and DSYEVD may read it in place.  The columns
    ! of `curv_u` are the eigenvectors, ascending in `curv_w`.
    curv_u = hess
    call oqp_dsyevd_f(npar, curv_u, curv_w, ierr)
    if (ierr /= 0) then
      status = CAS_ERR_EIGEN
      return
    end if

    do i = 1, npar
      wf(i) = merge(curv_w(i), shift, curv_w(i) > shift)
    end do
    ! ge = U^T g, then step = -U (ge / wf).  U is NOT symmetric, so the two
    ! products need opposite DGEMV trans flags; swapping them is the silent
    ! wrong-step failure documented in the porting notes.
    call dgemv('T', npar, npar, 1.0_dp, curv_u, npar, grad, 1, 0.0_dp, ge, 1)
    do i = 1, npar
      ge(i) = ge(i) / wf(i)
    end do
    call dgemv('N', npar, npar, -1.0_dp, curv_u, npar, ge, 1, 0.0_dp, step, 1)

    nrm = vnorm(step)
    if (nrm > maxrot) then
      do k = 0, npar - 1
        step(k) = step(k) * (maxrot / nrm)
      end do
    end if
    deallocate(hess, ge, wf)
  end subroutine newton_step

  !> `casscf.py::_fd_orbital_hessian`: 2*n_par gradient evaluations, then the
  !> explicit symmetrization.  This is where a CASSCF run spends its time.
  subroutine fd_orbital_hessian(ctx, cbuf, hh, hess, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(in) :: cbuf(0:,0:)
    real(dp), intent(in) :: hh
    real(dp), contiguous, intent(out) :: hess(:,:)
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: e(:), gp(:), gm(:), ctry(:,:)
    integer :: npar, k, i, j, ierr
    integer(i8) :: rc
    real(dp) :: obj, tmp

    npar = ctx%npar
    allocate(e(0:npar-1), gp(0:npar-1), gm(0:npar-1), &
             ctry(0:ctx%n-1, 0:ctx%n-1), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    e = 0.0_dp
    do k = 0, npar - 1
      e(k) = hh
      call cas_rotate(ctx, cbuf, e, ctry, rc)
      if (rc < 0_i8) then
        status = rc
        return
      end if
      call cas_evaluate(ctx, ctry, obj, gp, rc)
      if (rc < 0_i8) then
        status = rc
        return
      end if
      e(k) = -hh
      call cas_rotate(ctx, cbuf, e, ctry, rc)
      if (rc < 0_i8) then
        status = rc
        return
      end if
      call cas_evaluate(ctx, ctry, obj, gm, rc)
      if (rc < 0_i8) then
        status = rc
        return
      end if
      e(k) = 0.0_dp
      do i = 0, npar - 1
        hess(i + 1, k + 1) = (gp(i) - gm(i)) / (2.0_dp * hh)
      end do
    end do

    ! 0.5 * (H + H^T), element for element as the NumPy expression forms it
    do j = 1, npar
      do i = j, npar
        tmp = 0.5_dp * (hess(i, j) + hess(j, i))
        hess(i, j) = tmp
        hess(j, i) = tmp
      end do
    end do
    deallocate(e, gp, gm, ctry)
  end subroutine fd_orbital_hessian

  ! ================================================== orbital-Hessian backend
  !> Snapshot the CI vectors of the point just evaluated.
  !>
  !> The analytic Hessian is a function of the orbitals AND of the CI solution
  !> at them (`hess_fn(C, coeffs)` in the Python).  `ctx%civ` is a scratch
  !> buffer that every trial evaluation overwrites, so the accepted point's
  !> vectors have to be copied out the moment they are current.  On the FD path
  !> nothing reads them and the copy is skipped entirely.
  subroutine cas_snap_ci(ctx)
    type(cas_ctx_t), intent(inout) :: ctx
    if (ctx%hessmode == 0) return
    ctx%civ_ref = ctx%civ
  end subroutine cas_snap_ci

  !> `[casscf] hessian`: the finite-difference builder or the analytic one.
  !>
  !> The FD branch is the original code path, untouched.  The analytic branch
  !> re-transforms the AO integrals at `cbuf` first, exactly as the Python
  !> `hess_fn` does -- and not as an accident of tidiness: after a rejected DIIS
  !> extrapolation `ctx%h1e`/`ctx%eri` belong to the rejected trial point, not
  !> to the orbitals the Hessian is wanted at.
  subroutine build_hessian(ctx, cbuf, fdstep, hess, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(in) :: cbuf(0:,0:)
    real(dp), intent(in) :: fdstep
    real(dp), contiguous, intent(out) :: hess(:,:)
    integer(i8), intent(inout) :: status

    integer :: npar, i, j
    integer(i8) :: rc
    real(dp) :: tmp

    if (ctx%hessmode == 0) then
      call fd_orbital_hessian(ctx, cbuf, fdstep, hess, status)
      return
    end if

    npar = ctx%npar
    ctx%nhess = ctx%nhess + 1
    call mo_transform_h1e(int(ctx%n, c_int32_t), ctx%hcore, cbuf, ctx%h1e_h)
    if (mo_transform_eri(int(ctx%n, c_int32_t), ctx%eri_ao, cbuf, ctx%eri_h) &
        /= 0_i8) then
      status = CAS_ERR_TRANSFORM
      return
    end if
    rc = cas_anhess_build(ctx%ah, ctx%pairs, ctx%weights, ctx%roots, ctx%nroot, &
                          ctx%dets, ctx%h1e_h, ctx%eri_h, ctx%civ_ref, hess)
    if (rc < 0_i8) then
      ! A degeneracy refusal must NOT become a number.  It comes back as its
      ! own status so the Python path re-raises the message that explains the
      ! objective is not smooth there.
      if (rc == CAS_HESS_ERR_DEGEN) then
        status = CAS_ERR_DEGEN
      else
        status = CAS_ERR_ALLOC
      end if
      return
    end if
    ! The CI-relaxation half closes with a GEMM, so the two triangles agree
    ! only to rounding.  The Python feeds this matrix to `_symmetric_eigh`,
    ! which symmetrizes before diagonalizing; do the same so both arms pose the
    ! same eigenproblem.  (The FD Hessian is already exactly symmetric, which
    ! is why that branch is left alone.)
    do j = 1, npar
      do i = j + 1, npar
        tmp = 0.5_dp * (hess(i, j) + hess(j, i))
        hess(i, j) = tmp
        hess(j, i) = tmp
      end do
    end do
  end subroutine build_hessian

  ! =============================================== ah: trust-region augmented
  !> `casscf_convergers.py::_ah_inner`: one trust-region AH macroiteration loop.
  !>
  !> Each macroiteration builds the orbital Hessian, solves the level-shifted
  !> augmented-Hessian model in its eigenbasis (`casscf_ah_model_step`, which
  !> is CI-free), evaluates the trial point, and either accepts it or shrinks
  !> the radius and re-solves the same model.  `curv_w`/`curv_u` come back as
  !> the raw eigendecomposition built by the FINAL step, so the escape phase
  !> costs no extra CI solves -- the same contract `newton_loop` has.
  subroutine ah_inner(ctx, cbuf, maxmacro, gradtol, enertol, steptol, fdstep, &
                      par, stagnation_break, history, maxhist, nhist, &
                      converged, niter, curv_w, curv_u, have_curv, objective, &
                      stagnated, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, maxhist, stagnation_break
    real(dp), intent(in) :: gradtol, enertol, steptol, fdstep
    type(cas_ah_par_t), intent(in) :: par
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist
    logical, intent(out) :: converged, have_curv, stagnated
    integer, intent(out) :: niter
    real(dp), allocatable, intent(inout) :: curv_w(:), curv_u(:,:)
    real(dp), intent(out) :: objective
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: grad(:), grad_new(:), step(:), ctry(:,:)
    real(dp), allocatable :: hess(:,:), uc(:,:)
    real(dp) :: trust, gnorm, obj_old, obj_new, de, pred, step_norm, rho
    real(dp) :: best_obj, prd(1), shf(1)
    integer :: npar, it, rej, ierr, st, stall
    integer(c_int32_t) :: nmc(1)
    integer(i8) :: rc
    logical :: accepted

    npar = ctx%npar
    trust = par%start_trust
    converged = .false.
    stagnated = .false.
    have_curv = .false.
    niter = 0
    objective = 0.0_dp
    allocate(grad(0:npar-1), grad_new(0:npar-1), step(0:npar-1), &
             ctry(0:ctx%n-1, 0:ctx%n-1), hess(npar, npar), &
             uc(0:npar-1, 0:npar-1), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if
    if (allocated(curv_w)) deallocate(curv_w)
    if (allocated(curv_u)) deallocate(curv_u)
    allocate(curv_w(npar), curv_u(npar, npar), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    call cas_evaluate(ctx, cbuf, objective, grad, rc)
    if (rc < 0_i8) then
      status = rc
      return
    end if
    call cas_snap_ci(ctx)
    nhist = 0
    call push_history(history, maxhist, nhist, 0, objective, 0.0_dp, &
                      vnorm(grad), 0.0_dp)
    best_obj = objective
    stall = 0

    do it = 1, maxmacro
      niter = it
      gnorm = vnorm(grad)
      if (gnorm < gradtol) then
        converged = .true.
        exit
      end if

      call build_hessian(ctx, cbuf, fdstep, hess, status)
      if (status < 0_i8) return
      curv_u = hess
      call oqp_dsyevd_f(npar, curv_u, curv_w, ierr)
      if (ierr /= 0) then
        status = CAS_ERR_EIGEN
        return
      end if
      have_curv = .true.
      ! The AH engine reads the eigenvectors as the C-order buffer U[i][k];
      ! LAPACK left them in Fortran COLUMNS.  U is not symmetric, so this
      ! transpose is mandatory and omitting it is silent (every shape square).
      uc = transpose(curv_u)

      obj_old = objective
      accepted = .false.
      de = 0.0_dp
      pred = 0.0_dp
      obj_new = objective
      grad_new = grad
      do rej = 0, par%max_rejects
        st = int(casscf_ah_model_step(int(npar, c_int32_t), grad, curv_w, uc, &
                                      trust, int(par%max_micro, c_int32_t), &
                                      CAS_V0_TOL, step, shf, prd, nmc))
        if (st == 1) then
          ! no reference component: the pure negative-curvature case
          call casscf_lowest_mode_step(int(npar, c_int32_t), grad, curv_w, uc, &
                                       trust, step, prd)
        else if (st /= 0) then
          status = CAS_ERR_EIGEN
          return
        end if
        pred = prd(1)
        call cas_rotate(ctx, cbuf, step, ctry, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        call cas_evaluate(ctx, ctry, obj_new, grad_new, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        de = obj_new - obj_old
        if (de <= CAS_ACCEPT_TOL) then
          accepted = .true.
          exit
        end if
        trust = max(CAS_REJECT_FACTOR * trust, par%min_trust)
        if (trust <= par%min_trust) exit
      end do
      if (.not. accepted) then
        stagnated = .true.
        exit
      end if

      ! adaptive radius from the predicted-vs-actual decrease ratio
      step_norm = vnorm(step)
      if (pred < -CAS_PRED_FLOOR) then
        rho = de / pred
        if (rho < CAS_RHO_SHRINK) then
          trust = max(0.5_dp * trust, par%min_trust)
        else if (rho > CAS_RHO_GROW .and. step_norm >= 0.8_dp * trust) then
          trust = min(CAS_GROW_FACTOR * trust, par%max_trust)
        end if
      end if

      call cas_snap_ci(ctx)
      cbuf = ctry
      objective = obj_new
      grad = grad_new
      call push_history(history, maxhist, nhist, it, objective, &
                        objective - obj_old, vnorm(grad), step_norm)
      if (abs(objective - obj_old) < enertol .and. step_norm < steptol) then
        converged = vnorm(grad) < gradtol
        exit
      end if
      if (stagnation_break > 0) then
        if (objective < best_obj - enertol) then
          best_obj = objective
          stall = 0
        else
          stall = stall + 1
          if (stall >= stagnation_break) then
            stagnated = .true.
            exit
          end if
        end if
      end if
    end do
    deallocate(grad, grad_new, step, ctry, hess, uc)
  end subroutine ah_inner

  !> `casscf_convergers.py::_curvature_escape`.
  !>
  !> Same idea as `escape_saddles` -- kick along the lowest mode, line-search
  !> both signs, re-converge, accept only a strict energy gain -- but it
  !> re-converges with the AH loop so the converger stays AH end to end, and it
  !> stamps each rerun row with `total_it + that row's own index` where
  !> `escape_saddles` stamps them all with the running total.  The two differ
  !> in the Python and both are reproduced.
  subroutine ah_curvature_escape(ctx, cbuf, maxmacro, gradtol, enertol, &
                                 steptol, fdstep, par, maxescape, history, &
                                 maxhist, nhist, total_it, curv_w, curv_u, &
                                 have_curv, objective, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, maxescape, maxhist
    real(dp), intent(in) :: gradtol, enertol, steptol, fdstep
    type(cas_ah_par_t), intent(in) :: par
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist, total_it
    real(dp), allocatable, intent(inout) :: curv_w(:), curv_u(:,:)
    logical, intent(in) :: have_curv
    real(dp), intent(inout) :: objective
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: kick(:), grad(:), cbest(:,:), ctry(:,:), c2(:,:)
    real(dp), allocatable :: w2(:), u2(:,:), hrow(:,:)
    integer :: npar, escapes, amp, kmin, ierr, it2, nrow2, r
    integer(i8) :: rc
    real(dp) :: best_obj, on, obj2
    logical :: conv2, have2, stag2

    if (.not. have_curv) return          ! `curv is None`: nothing to inspect

    npar = ctx%npar
    allocate(kick(0:npar-1), grad(0:npar-1), cbest(0:ctx%n-1, 0:ctx%n-1), &
             ctry(0:ctx%n-1, 0:ctx%n-1), c2(0:ctx%n-1, 0:ctx%n-1), &
             hrow(5, maxhist), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    escapes = 0
    do while (escapes < maxescape)
      kmin = minloc(curv_w, 1)
      if (curv_w(kmin) >= -par%saddle_curv_tol) exit

      best_obj = 0.0_dp
      do amp = 1, size(CAS_ESCAPE_AMPS)
        kick = CAS_ESCAPE_AMPS(amp) * curv_u(1:npar, kmin)
        call cas_rotate(ctx, cbuf, kick, ctry, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        call cas_evaluate(ctx, ctry, on, grad, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        if (amp == 1 .or. on < best_obj) then
          best_obj = on
          cbest = ctry
        end if
      end do

      c2 = cbest
      nrow2 = 0
      call ah_inner(ctx, c2, maxmacro, gradtol, enertol, steptol, fdstep, par, &
                    0, hrow, maxhist, nrow2, conv2, it2, w2, u2, have2, obj2, &
                    stag2, status)
      if (status < 0_i8) return
      do r = 2, nrow2
        call push_history(history, maxhist, nhist, total_it + int(hrow(1, r)), &
                          hrow(2, r), hrow(3, r), hrow(4, r), hrow(5, r))
      end do
      total_it = total_it + it2
      if (.not. conv2 .or. obj2 >= objective - par%saddle_egain_tol .or. &
          .not. have2) exit

      cbuf = c2
      objective = obj2
      if (allocated(curv_w)) deallocate(curv_w)
      if (allocated(curv_u)) deallocate(curv_u)
      call move_alloc(w2, curv_w)
      call move_alloc(u2, curv_u)
      escapes = escapes + 1
    end do
    deallocate(kick, grad, cbest, ctry, c2, hrow)
  end subroutine ah_curvature_escape

  !> `casscf_convergers.py::_ah_converge`: the AH loop plus its escape phase.
  subroutine ah_converge(ctx, cbuf, maxmacro, gradtol, enertol, steptol, &
                         fdstep, par, stagnation_break, maxescape, history, &
                         maxhist, nhist, converged, total_it, objective, &
                         stagnated, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, maxhist, stagnation_break, maxescape
    real(dp), intent(in) :: gradtol, enertol, steptol, fdstep
    type(cas_ah_par_t), intent(in) :: par
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist
    logical, intent(out) :: converged, stagnated
    integer, intent(out) :: total_it
    real(dp), intent(out) :: objective
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: curv_w(:), curv_u(:,:)
    integer :: niter
    logical :: have_curv

    call ah_inner(ctx, cbuf, maxmacro, gradtol, enertol, steptol, fdstep, par, &
                  stagnation_break, history, maxhist, nhist, converged, niter, &
                  curv_w, curv_u, have_curv, objective, stagnated, status)
    if (status < 0_i8) return
    total_it = niter
    if (converged) then
      call ah_curvature_escape(ctx, cbuf, maxmacro, gradtol, enertol, steptol, &
                               fdstep, par, maxescape, history, maxhist, nhist, &
                               total_it, curv_w, curv_u, have_curv, objective, &
                               status)
      if (status < 0_i8) return
      converged = .true.
    end if
  end subroutine ah_converge

  ! ================================================================== diis
  !> `casscf_convergers.py::_diis_optimize`.
  !>
  !> Each macroiteration takes the production step (Newton or Powell) with the
  !> production backtracking, records the (accumulated-rotation, gradient)
  !> pair, and once `diis_start` pairs exist EVALUATES the Pulay-extrapolated
  !> point, keeping whichever candidate has the lower true energy.  Because
  !> every candidate is re-evaluated variationally the acceleration can only
  !> match or improve the plain iteration, which is also why the BCH-truncated
  !> accumulated rotation does not need to be exact.
  subroutine diis_optimize(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                           steptol, maxrot, shift, fdstep, saddle_c, saddle_e, &
                           maxescape, nspace_in, nstart_in, history, maxhist, &
                           nhist, converged, total_it, objective, n_used, &
                           n_extrap, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, optimizer, maxhist, maxescape
    integer, intent(in) :: nspace_in, nstart_in
    real(dp), intent(in) :: gradtol, enertol, steptol, maxrot, shift, fdstep
    real(dp), intent(in) :: saddle_c, saddle_e
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist
    logical, intent(out) :: converged
    integer, intent(out) :: total_it, n_used, n_extrap
    real(dp), intent(out) :: objective
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: grad(:), grad_new(:), grad_x(:), step(:), astep(:)
    real(dp), allocatable :: tvec(:), tnew(:), txt(:), coef(:)
    real(dp), allocatable :: c0(:,:), ctry(:,:), cacc(:,:), cx(:,:)
    real(dp), allocatable :: gstore(:,:), tstore(:,:)
    real(dp), allocatable :: curv_w(:), curv_u(:,:)
    integer(c_int32_t) :: nused(1)
    integer :: npar, nspace, nstart, nstore, it, bt, k, ierr, ncf
    integer(i8) :: rc
    real(dp) :: obj_old, obj_new, obj_x, nrm, step_norm
    logical :: have_curv

    npar = ctx%npar
    nspace = max(2, nspace_in)
    nstart = max(2, nstart_in)
    converged = .false.
    have_curv = .false.
    n_used = 0
    n_extrap = 0
    total_it = 0
    objective = 0.0_dp

    allocate(grad(0:npar-1), grad_new(0:npar-1), grad_x(0:npar-1), &
             step(0:npar-1), astep(0:npar-1), tvec(0:npar-1), &
             tnew(0:npar-1), txt(0:npar-1), coef(0:nspace-1), &
             c0(0:ctx%n-1, 0:ctx%n-1), ctry(0:ctx%n-1, 0:ctx%n-1), &
             cacc(0:ctx%n-1, 0:ctx%n-1), cx(0:ctx%n-1, 0:ctx%n-1), &
             gstore(0:npar-1, 0:nspace-1), tstore(0:npar-1, 0:nspace-1), &
             stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    c0 = cbuf
    tvec = 0.0_dp
    nstore = 0
    call cas_evaluate(ctx, cbuf, objective, grad, rc)
    if (rc < 0_i8) then
      status = rc
      return
    end if
    call cas_snap_ci(ctx)
    nhist = 0
    call push_history(history, maxhist, nhist, 0, objective, 0.0_dp, &
                      vnorm(grad), 0.0_dp)

    do it = 1, maxmacro
      total_it = it
      if (vnorm(grad) < gradtol) then
        converged = .true.
        exit
      end if

      if (optimizer == 0) then
        call newton_step(ctx, cbuf, grad, fdstep, shift, maxrot, step, &
                         curv_w, curv_u, status)
        if (status < 0_i8) return
        have_curv = .true.
      else
        step = -grad
        nrm = vnorm(step)
        if (nrm > maxrot) step = step * (maxrot / nrm)
      end if

      obj_old = objective
      astep = step
      do bt = 1, CAS_DIIS_BACKTRACK
        call cas_rotate(ctx, cbuf, astep, ctry, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        call cas_evaluate(ctx, ctry, obj_new, grad_new, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        if (obj_new <= obj_old + CAS_ACCEPT_TOL) exit
        astep = 0.5_dp * astep
      end do
      call cas_snap_ci(ctx)          ! candidate so far: the plain step
      cacc = ctry
      tnew = tvec + astep

      ! FIFO of (accumulated rotation, gradient) pairs, oldest column first --
      ! which is also the layout `casscf_diis_coeffs` wants (its C-order
      ! [nvec, npar] buffer IS this Fortran (npar, nvec) array).
      if (nstore == nspace) then
        do k = 0, nspace - 2
          gstore(:, k) = gstore(:, k + 1)
          tstore(:, k) = tstore(:, k + 1)
        end do
        nstore = nspace - 1
      end if
      gstore(:, nstore) = grad_new
      tstore(:, nstore) = tnew
      nstore = nstore + 1

      if (nstore >= nstart) then
        call casscf_diis_coeffs(int(nstore, c_int32_t), int(npar, c_int32_t), &
                                gstore, CAS_DIIS_CONDMAX, coef, nused)
        ncf = int(nused(1))
        if (ncf > 0) then
          txt = 0.0_dp
          do k = 0, ncf - 1
            txt = txt + coef(k) * tstore(:, nstore - ncf + k)
          end do
          call cas_rotate(ctx, c0, txt, cx, rc)
          if (rc < 0_i8) then
            status = rc
            return
          end if
          call cas_evaluate(ctx, cx, obj_x, grad_x, rc)
          if (rc < 0_i8) then
            status = rc
            return
          end if
          n_extrap = n_extrap + 1
          if (obj_x < obj_new - CAS_ACCEPT_TOL) then
            n_used = n_used + 1
            call cas_snap_ci(ctx)
            cacc = cx
            obj_new = obj_x
            grad_new = grad_x
            tnew = txt
            gstore(:, nstore - 1) = grad_x
            tstore(:, nstore - 1) = txt
          end if
        end if
      end if

      step_norm = vnorm(tnew - tvec)
      cbuf = cacc
      objective = obj_new
      grad = grad_new
      tvec = tnew
      call push_history(history, maxhist, nhist, it, objective, &
                        objective - obj_old, vnorm(grad), step_norm)
      if (abs(objective - obj_old) < enertol .and. step_norm < steptol) then
        converged = vnorm(grad) < gradtol
        exit
      end if
    end do

    if (converged .and. optimizer == 0 .and. maxescape > 0 .and. have_curv) then
      call escape_saddles(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                          steptol, maxrot, shift, fdstep, saddle_c, saddle_e, &
                          maxescape, history, maxhist, nhist, total_it, &
                          curv_w, curv_u, objective, status)
      if (status < 0_i8) return
      converged = .true.
    end if
  end subroutine diis_optimize

  ! ================================================================== auto
  !> `casscf_convergers.py::_auto_optimize`: `ah` with a stagnation watchdog
  !> and an automatic fallback to the two-phase converger.
  !>
  !> The AH loop accepts only non-increasing steps, so its final point is its
  !> best point and the fallback never restarts above the initial one; at worst
  !> `auto` reproduces the default converger's trajectory from there.
  subroutine auto_optimize(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                           steptol, maxrot, shift, fdstep, saddle_c, saddle_e, &
                           maxescape, par, stagnation, history, maxhist, nhist, &
                           converged, total_it, objective, auto_code, auto_it, &
                           stagnated, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, optimizer, maxhist, maxescape, stagnation
    real(dp), intent(in) :: gradtol, enertol, steptol, maxrot, shift, fdstep
    real(dp), intent(in) :: saddle_c, saddle_e
    type(cas_ah_par_t), intent(in) :: par
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist
    logical, intent(out) :: converged, stagnated
    integer, intent(out) :: total_it, auto_code, auto_it
    real(dp), intent(out) :: objective
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: curv_w(:), curv_u(:,:), hrow(:,:)
    integer :: nrow2, it2, offset, r, ierr
    logical :: conv1, conv2, have2

    auto_code = 0
    auto_it = 0
    call ah_converge(ctx, cbuf, maxmacro, gradtol, enertol, steptol, fdstep, &
                     par, stagnation, maxescape, history, maxhist, nhist, &
                     conv1, total_it, objective, stagnated, status)
    if (status < 0_i8) return
    if (conv1) then
      converged = .true.
      return
    end if

    auto_code = merge(1, 2, stagnated)     ! 1 stalled, 2 hit the macro cap
    auto_it = total_it
    allocate(hrow(5, maxhist), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if
    nrow2 = 0
    call newton_loop(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, steptol, &
                     maxrot, shift, fdstep, hrow, maxhist, nrow2, .true., &
                     conv2, it2, curv_w, curv_u, have2, objective, status)
    if (status < 0_i8) return
    offset = total_it
    do r = 2, nrow2
      call push_history(history, maxhist, nhist, offset + int(hrow(1, r)), &
                        hrow(2, r), hrow(3, r), hrow(4, r), hrow(5, r))
    end do
    total_it = offset + it2
    converged = conv2
    deallocate(hrow)

    if (conv2 .and. optimizer == 0 .and. maxescape > 0 .and. have2) then
      call escape_saddles(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                          steptol, maxrot, shift, fdstep, saddle_c, saddle_e, &
                          maxescape, history, maxhist, nhist, total_it, &
                          curv_w, curv_u, objective, status)
      if (status < 0_i8) return
      converged = .true.
    end if
  end subroutine auto_optimize

  ! ========================================================= saddle escape
  !> `casscf.py::_escape_saddles`.  Only a DEEP negative curvature triggers a
  !> step, and the escaped point is kept only when it strictly lowers the
  !> energy, so the phase is a guaranteed no-op wherever phase 1 already found
  !> a minimum.
  subroutine escape_saddles(ctx, cbuf, maxmacro, optimizer, gradtol, enertol, &
                            steptol, maxrot, shift, fdstep, curv_tol, &
                            egain_tol, maxescape, history, maxhist, nhist, &
                            total_it, curv_w, curv_u, objective, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer, intent(in) :: maxmacro, optimizer, maxescape, maxhist
    real(dp), intent(in) :: gradtol, enertol, steptol, maxrot, shift, fdstep
    real(dp), intent(in) :: curv_tol, egain_tol
    real(dp), intent(inout) :: history(0:*)
    integer, intent(inout) :: nhist, total_it
    real(dp), allocatable, intent(inout) :: curv_w(:), curv_u(:,:)
    real(dp), intent(inout) :: objective
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: kick(:), grad(:), cbest(:,:), ctry(:,:)
    real(dp), allocatable :: c2(:,:), w2(:), u2(:,:), hrow(:,:)
    integer :: npar, escapes, amp, kmin, ierr, it2, nrow2, r
    integer(i8) :: rc
    real(dp) :: best_obj, on, obj2
    logical :: conv2, have2

    npar = ctx%npar
    allocate(kick(0:npar-1), grad(0:npar-1), cbest(0:ctx%n-1, 0:ctx%n-1), &
             ctry(0:ctx%n-1, 0:ctx%n-1), c2(0:ctx%n-1, 0:ctx%n-1), &
             hrow(5, maxhist), stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if

    escapes = 0
    do while (escapes < maxescape)
      kmin = minloc(curv_w, 1)
      if (curv_w(kmin) >= -curv_tol) exit

      best_obj = 0.0_dp
      do amp = 1, size(CAS_ESCAPE_AMPS)
        kick = CAS_ESCAPE_AMPS(amp) * curv_u(1:npar, kmin)
        call cas_rotate(ctx, cbuf, kick, ctry, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        call cas_evaluate(ctx, ctry, on, grad, rc)
        if (rc < 0_i8) then
          status = rc
          return
        end if
        if (amp == 1 .or. on < best_obj) then
          best_obj = on
          cbest = ctry
        end if
      end do

      ! re-converge from the kicked point into a PRIVATE history buffer: the
      ! rows are only committed once, and (as in the Python) every row of a
      ! rerun is stamped with the running total, not with its own index
      c2 = cbest
      nrow2 = 0
      call newton_loop(ctx, c2, maxmacro, optimizer, gradtol, enertol, &
                       steptol, maxrot, shift, fdstep, hrow, maxhist, nrow2, &
                       .true., conv2, it2, w2, u2, have2, obj2, status)
      if (status < 0_i8) return
      total_it = total_it + it2
      do r = 2, nrow2
        call push_history(history, maxhist, nhist, total_it, hrow(2, r), &
                          hrow(3, r), hrow(4, r), hrow(5, r))
      end do
      if (.not. conv2 .or. obj2 >= objective - egain_tol .or. .not. have2) exit

      cbuf = c2
      objective = obj2
      if (allocated(curv_w)) deallocate(curv_w)
      if (allocated(curv_u)) deallocate(curv_u)
      call move_alloc(w2, curv_w)
      call move_alloc(u2, curv_u)
      escapes = escapes + 1
    end do
  end subroutine escape_saddles

  ! ======================================================= canonicalization
  !> `casscf.py::_canonicalize`: diagonalize the inactive / active / virtual
  !> blocks of the closed+active mean-field Fock so the inactive and virtual
  !> orbitals carry orbital energies.
  !>
  !> Note the density: the Python uses the ROOT-0 CI density with weight 1.0
  !> even for SA-CASSCF (`_solve_active(..., [1.0], [0])`), so this does too.
  subroutine canonicalize(ctx, cbuf, status)
    type(cas_ctx_t), intent(inout) :: ctx
    real(dp), contiguous, intent(inout) :: cbuf(0:,0:)
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: dmat(:), feff(:), sub(:,:), evals(:), blk(:,:)
    integer :: n, nc, na, i, j, ierr, lo, hi, m, b
    integer(i8) :: rc, jd
    real(dp) :: sym

    n = ctx%n
    nc = ctx%nc
    na = ctx%na

    call mo_transform_h1e(int(n, c_int32_t), ctx%hcore, cbuf, ctx%h1e)
    if (mo_transform_eri(int(n, c_int32_t), ctx%eri_ao, cbuf, ctx%eri) /= 0_i8) then
      status = CAS_ERR_TRANSFORM
      return
    end if
    ctx%iopt_ci(FCI_I_NROOT) = 1_c_int32_t
    rc = fci_solve(ctx%iopt_ci, ctx%dopt_ci, ctx%active, ctx%core, ctx%h1e, &
                   ctx%eri, ctx%enci, ctx%civ, ctx%s2w)
    ctx%iopt_ci(FCI_I_NROOT) = int(ctx%nroot, c_int32_t)
    if (rc < 0_i8) then
      status = CAS_ERR_CI
      return
    end if

    allocate(dmat(0:int(n, i8)**2 - 1_i8), feff(0:int(n, i8)**2 - 1_i8), &
             stat=ierr)
    if (ierr /= 0) then
      status = CAS_ERR_ALLOC
      return
    end if
    ! nroot == 1 here, so the CI vector is already contiguous
    do jd = 0_i8, ctx%ndet - 1_i8
      ctx%cvec(jd) = ctx%civ(jd)
    end do
    call rdm1_spatial(int(na, c_int32_t), ctx%ndet, ctx%dets, ctx%cvec, ctx%g1)

    ! `_full_rdm1`: 2 on the inactive diagonal, the active 1-RDM on the active
    ! block, zero everywhere else.
    dmat = 0.0_dp
    do i = 0, nc - 1
      dmat(int(i, i8)*int(n, i8) + int(i, i8)) = 2.0_dp
    end do
    do i = 0, na - 1
      do j = 0, na - 1
        dmat(int(nc + i, i8)*int(n, i8) + int(nc + j, i8)) = &
            ctx%g1(int(i, i8)*int(na, i8) + int(j, i8))
      end do
    end do

    call casscf_effective_fock(int(n, c_int32_t), dmat, ctx%h1e, ctx%eri, feff)

    do b = 1, 3
      select case (b)
      case (1)
        lo = 0
        hi = nc - 1
      case (2)
        lo = nc
        hi = nc + na - 1
      case default
        lo = nc + na
        hi = n - 1
      end select
      m = hi - lo + 1
      if (m < 2) cycle
      allocate(sub(m, m), evals(m), blk(m, n), stat=ierr)
      if (ierr /= 0) then
        status = CAS_ERR_ALLOC
        return
      end if
      ! Feff is symmetric here, but the Python still forms 0.5*(sub + sub^T)
      ! before diagonalizing, so the same expression is written out.
      do i = 1, m
        do j = 1, m
          sym = 0.5_dp * (feff(int(lo + i - 1, i8)*int(n, i8) + int(lo + j - 1, i8)) &
                        + feff(int(lo + j - 1, i8)*int(n, i8) + int(lo + i - 1, i8)))
          sub(i, j) = sym
        end do
      end do
      call oqp_dsyevd_f(m, sub, evals, ierr)
      if (ierr /= 0) then
        status = CAS_ERR_EIGEN
        return
      end if
      ! Cnew[:, idx] = C[:, idx] @ V, i.e. in the C-order buffer (MO index
      ! first) blk = V^T * cbuf(lo:hi, :).
      call dgemm('T', 'N', m, n, m, 1.0_dp, sub, m, cbuf(lo, 0), n, &
                 0.0_dp, blk, m)
      do i = 1, m
        do j = 0, n - 1
          cbuf(lo + i - 1, j) = blk(i, j + 1)
        end do
      end do
      deallocate(sub, evals, blk)
    end do
    deallocate(dmat, feff)
  end subroutine canonicalize

  ! ================================================================ helpers
  !> Append one macroiteration row to the caller's C-order [maxhist, 5] table.
  subroutine push_history(history, maxhist, nhist, it, e, de, gnorm, snorm)
    real(dp), intent(inout) :: history(0:*)
    integer, intent(in) :: maxhist, it
    integer, intent(inout) :: nhist
    real(dp), intent(in) :: e, de, gnorm, snorm
    integer(i8) :: base

    if (nhist >= maxhist) return
    base = int(nhist, i8) * 5_i8
    history(base)     = real(it, dp)
    history(base + 1) = e
    history(base + 2) = de
    history(base + 3) = gnorm
    history(base + 4) = snorm
    nhist = nhist + 1
  end subroutine push_history

  pure function vnorm(v) result(nrm)
    real(dp), intent(in) :: v(0:)
    real(dp) :: nrm
    nrm = sqrt(sum(v * v))
  end function vnorm

  pure function binom(n, k) result(c)
    integer, intent(in) :: n, k
    integer(i8) :: c
    integer :: i
    c = 0_i8
    if (k < 0 .or. k > n) return
    c = 1_i8
    do i = 1, min(k, n - k)
      c = c * int(n - min(k, n - k) + i, i8) / int(i, i8)
    end do
  end function binom

end module casscf_driver_mod
