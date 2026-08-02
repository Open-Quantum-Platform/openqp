!> @file afqmc_oqp.F90
!>
!> @brief OpenQP entry point for phaseless AFQMC.
!>
!> This is the whole compute path: `[input] method = afqmc` dispatches to the
!> `afqmc` C symbol, which reads the converged SCF reference straight out of the
!> tagarray, builds the Cholesky-factorized molecular Hamiltonian in memory, and
!> runs the propagation.  Nothing is written to or read from disk in between,
!> and no part of the loop returns to Python.
!>
!> It replaces the earlier arrangement, where a Python bridge
!> (`openqp_to_afqmc.py`) ran OpenQP, factorized the ERIs with NumPy, wrote a
!> HAMILTONIAN text file, and spawned a separate `openqp-afqmc-native`
!> executable.
!>
!> The kernels under source/afqmc/ are vendored from
!> Open-Quantum-Platform/AFQMC (see afqmc_host.F90 for the host shims).

module afqmc_oqp_mod

  use precision, only: dp

  implicit none

  private
  public :: afqmc

  character(len=*), parameter :: module_name = "afqmc_oqp_mod"

  !> ---------------------------------------------------------------------
  !> [afqmc] option vector -- the single place this schema is defined on the
  !> Fortran side.  The Python side mirrors it in
  !> pyoqp/oqp/library/afqmc_driver.py (_pack_options); the two must agree.
  !> Everything travels as real64 so one tagarray record carries the lot.
  !> ---------------------------------------------------------------------
  integer, parameter, public :: AFQMC_NOPT       = 18
  integer, parameter :: OPT_TRIAL_KIND           = 1   !< see the trial-kind codes below
  integer, parameter :: OPT_CHOL_TOL             = 2
  integer, parameter :: OPT_NWALKERS             = 3
  integer, parameter :: OPT_NSTEPS               = 4
  integer, parameter :: OPT_TIMESTEP             = 5
  integer, parameter :: OPT_SEED                 = 6
  integer, parameter :: OPT_STABILIZE_EVERY      = 7
  integer, parameter :: OPT_POPCONTROL_EVERY     = 8
  integer, parameter :: OPT_ESTIMATE_EVERY       = 9
  integer, parameter :: OPT_ACCUMULATE_AFTER     = 10
  integer, parameter :: OPT_FORCE_BIAS_BOUND     = 11
  integer, parameter :: OPT_OO_ORBITALS          = 12  !< 1 -> apply the OOORBDAT rotation
  integer, parameter :: OPT_NLOW                 = 13  !< lower states to deflate (LOWSTATE)
  integer, parameter :: OPT_LOW_MAX              = 14
  integer, parameter :: OPT_LOW_CAP              = 15
  integer, parameter :: OPT_LOW_START            = 16
  integer, parameter :: OPT_TRIAL_STATE          = 17  !< 1-based MRSF root for TRIAL_MRSF_STATE
  integer, parameter :: OPT_TRIAL_THRESH         = 18  !< drop trial determinants below this |c|

  !> Trial-wavefunction kinds.  0-2 keep the upstream AFQMC numbering (the
  !> vendored kernels branch on `MRSF_CIS`); 3 is the OpenQP addition that takes
  !> the trial straight from an MRSF-TDDFT root in memory, with no file.
  integer, parameter, public :: TRIAL_MEAN_FIELD = 0
  integer, parameter, public :: TRIAL_SF_FILE    = 1
  integer, parameter, public :: TRIAL_MRSF_FILE  = 2
  integer, parameter, public :: TRIAL_MRSF_STATE = 3

contains

!> C-bound entry point: `[input] method = afqmc` dispatches here.
  subroutine afqmc_C(c_handle) bind(C, name="afqmc")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call afqmc(inf)
  end subroutine afqmc_C

!> @brief Phaseless AFQMC on the converged SCF reference.
  subroutine afqmc(infos)

    use types, only: information
    use io_constants, only: iw
    use basis_tools, only: basis_set
    use printing, only: print_module_info
    use messages, only: show_message, WITH_ABORT
    use oqp_tagarray_driver
    use int2e_mod, only: int2e
    use AFQMC_Module, only: AFQMC_Control_Type
    use IO_Files, only: afqmc_iw => IW
    use afqmc_trial_mod, only: build_mrsf_state_trial, read_trial_header

    implicit none

    interface
      subroutine AFQMC_Run(NVar, NSpin, NOcc, NChol, ENu, Hc, Chol, Cntrl, &
                           MRSF_CIS, MRSF_File, OO_Orb, OO_File, NLow, Low_File, &
                           Low_Max, Low_Cap, Low_Start, Results, Ext_Coeff, Ext_Occ)
        use AFQMC_Module, only: AFQMC_Control_Type
        implicit none
        integer, intent(in) :: NVar, NSpin, NChol, NOcc(1:NSpin)
        double precision, intent(in) :: ENu
        double precision, intent(in) :: Hc(1:NVar*(NVar+1)/2)
        double precision, intent(in) :: Chol(1:NVar,1:NVar,1:NChol)
        type(AFQMC_Control_Type), intent(in) :: Cntrl
        integer, intent(in) :: MRSF_CIS, OO_Orb, NLow, Low_Start
        character(len=8), intent(in) :: MRSF_File, OO_File, Low_File
        double precision, intent(in) :: Low_Max, Low_Cap
        double precision, optional, intent(out) :: Results(1:4)
        double complex, optional, intent(in) :: Ext_Coeff(:)
        integer, optional, intent(in) :: Ext_Occ(:,:,:)
      end subroutine
    end interface

    character(len=*), parameter :: subroutine_name = "afqmc"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    real(kind=dp), contiguous, pointer :: hcore(:), mo_a(:,:), eri_flat(:), opts(:)
    character(len=*), parameter :: tags_required(2) = (/ character(len=80) :: &
      OQP_Hcore, OQP_VEC_MO_A /)

    type(AFQMC_Control_Type) :: cntrl
    real(kind=dp), allocatable :: h1_mo(:), chol_ao(:,:), chol_mo(:,:,:)
    real(kind=dp) :: results(4), residual, enuc
    integer :: nbf, nmo, nchol, nocc(2), nspin, ig
    integer :: trial_kind, oo_orb, nlow, low_start, trial_state, mrsf_flag
    integer :: ndet_hdr
    real(kind=dp) :: low_max, low_cap, trial_thresh
    character(len=8) :: trial_file, oo_file, low_file
    complex(kind=dp), allocatable :: trial_coeff(:)
    integer, allocatable :: trial_occ(:,:,:)

    open (unit=iw, file=infos%log_filename, position="append")
    afqmc_iw = iw

    call print_module_info('AFQMC', &
      'Phaseless auxiliary-field quantum Monte Carlo')

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nmo = nbf

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    call read_options(infos, opts)
    call unpack_controls(opts, cntrl, trial_kind, oo_orb, nlow, low_max, low_cap, &
                         low_start, trial_state, trial_thresh)

!   --- trial wavefunction ----------------------------------------------------
!   A spin-flip trial lives in a different M_S sector than the reference, so the
!   electron counts come from the trial, not from the SCF.
    nspin = 2
    trial_file = '        '
    select case (trial_kind)

    case (TRIAL_MEAN_FIELD)
      nocc(1) = int(infos%mol_prop%nelec_A)
      nocc(2) = int(infos%mol_prop%nelec_B)

    case (TRIAL_SF_FILE, TRIAL_MRSF_FILE)
      trial_file = trial_file_name(trial_kind)
      call read_trial_header(trim(trial_file), ndet_hdr, nocc(1), nocc(2))
      write(iw,'(/1x,"AFQMC trial: ",i0," determinants from ",a)') &
        ndet_hdr, trim(trial_file)

    case (TRIAL_MRSF_STATE)
      call build_mrsf_state_trial(infos, trial_state, trial_thresh, nocc, &
                                  trial_coeff, trial_occ)

    end select

    if (max(nocc(1), nocc(2)) > nmo) &
      call show_message("afqmc: electron count exceeds the orbital count", WITH_ABORT)

!   The vendored kernels branch on MRSF_CIS /= 0 to select the multideterminant
!   overlap/force-bias/local-energy path; an in-memory MRSF trial takes the same
!   route as a file-based one.
    mrsf_flag = trial_kind
    if (trial_kind == TRIAL_MRSF_STATE) mrsf_flag = TRIAL_MRSF_FILE

    enuc = infos%mol_energy%vnn
    if (enuc == 0.0_dp) enuc = infos%mol_energy%nenergy

    write(iw,'(/1x,"Orbitals (nmo)          : ",i0)') nmo
    write(iw,'(1x,"Electrons (alpha/beta)  : ",i0,"/",i0)') nocc(1), nocc(2)
    write(iw,'(1x,"Core (nuclear) energy   : ",f20.12)') enuc
    write(iw,'(1x,"Walkers / steps         : ",i0," / ",i0)') cntrl%NWlk, cntrl%NStep
    write(iw,'(1x,"Time step / seed        : ",f12.6," / ",i0)') cntrl%dT, cntrl%ISeed

!   --- one-electron: h_pq = C^T h_munu C, kept packed like OQP::Hcore ---
    allocate(h1_mo(nmo*(nmo+1)/2))
    call AFQMC_TFT(nmo, nbf, "triangle", hcore, nbf, mo_a, nbf, h1_mo, nmo)

!   --- two-electron: Cholesky-factorize the AO ERIs, then rotate to MOs ------
!   The factorization is done in the AO basis and each vector transformed
!   afterwards.  That is exact, not an approximation: (pq|rs) = sum_g L^g_pq
!   L^g_rs is preserved term by term under L^g -> C^T L^g C, and the AO route
!   avoids ever forming the nmo^4 MO tensor.
    close(iw)
    call int2e(infos)
    open (unit=iw, file=infos%log_filename, position="append")
    call tagarray_get_data(infos%dat, OQP_ERI_AO, eri_flat)

    call pivoted_cholesky(eri_flat, nbf, cntrl%Chol_Tol, chol_ao, nchol, residual)

    write(iw,'(/1x,"Cholesky vectors        : ",i0," (threshold ",es9.2,")")') &
      nchol, cntrl%Chol_Tol
    write(iw,'(1x,"Residual diagonal bound : ",es12.4)') residual
    write(iw,'(1x,"max|h1(MO)| / max|L(AO)|: ",es12.4," / ",es12.4)') &
      maxval(abs(h1_mo)), maxval(abs(chol_ao))
    if (.not. all(h1_mo == h1_mo) .or. .not. all(chol_ao == chol_ao)) &
      call show_message("afqmc: the molecular Hamiltonian is not finite", WITH_ABORT)

    allocate(chol_mo(nmo, nmo, nchol))
    do ig = 1, nchol
      call AFQMC_TFT(nmo, nbf, "square", chol_ao(1,ig), nbf, mo_a, nbf, &
                     chol_mo(1,1,ig), nmo)
    end do
    deallocate(chol_ao)
    write(iw,'(1x,"max|L(MO)|              : ",es12.4)') maxval(abs(chol_mo))
    if (.not. all(chol_mo == chol_mo)) &
      call show_message("afqmc: the MO Cholesky vectors are not finite", WITH_ABORT)

!   --- propagate -------------------------------------------------------------
    oo_file = 'OOORBDAT'
    low_file = 'LOWSTATE'

!   AFQMC_Run writes its per-step table to the same unit, so the log stays open
!   across the propagation.
    if (allocated(trial_coeff)) then
      call AFQMC_Run(nmo, nspin, nocc, nchol, enuc, h1_mo, chol_mo, cntrl, &
                     mrsf_flag, trial_file, oo_orb, oo_file, nlow, low_file, &
                     low_max, low_cap, low_start, results, trial_coeff, trial_occ)
    else
      call AFQMC_Run(nmo, nspin, nocc, nchol, enuc, h1_mo, chol_mo, cntrl, &
                     mrsf_flag, trial_file, oo_orb, oo_file, nlow, low_file, &
                     low_max, low_cap, low_start, results)
    end if

    write(iw,'(/,2x,60("="))')
    write(iw,'(2x,a)') 'AFQMC (phaseless, ph-AFQMC)'
    write(iw,'(2x,60("="))')
    write(iw,'(2x,a,f20.10)') 'E(mixed estimator, final)   = ', results(1)
    write(iw,'(2x,a,f20.10)') 'E(hybrid estimator, final)  = ', results(2)
    write(iw,'(2x,a,f20.10)') 'E(mixed, block average)     = ', results(3)
    write(iw,'(2x,a,f20.10)') 'E(hybrid, block average)    = ', results(4)
    write(iw,'(2x,60("="),/)')

!   The block-averaged mixed estimator is the AFQMC energy; fall back to the
!   final mixed value when no post-equilibration block was accumulated.
    if (cntrl%NStep > cntrl%NStep_Accum) then
      infos%mol_energy%energy = results(3)
    else
      infos%mol_energy%energy = results(1)
    end if
    infos%mol_energy%etot = infos%mol_energy%energy

    deallocate(h1_mo, chol_mo)
    close(iw)

  end subroutine afqmc

!> @brief Fetch the [afqmc] option vector written by the Python driver.
  subroutine read_options(infos, opts)
    use types, only: information
    use oqp_tagarray_driver
    use messages, only: show_message, WITH_ABORT
    type(information), target, intent(inout) :: infos
    real(kind=dp), contiguous, pointer, intent(out) :: opts(:)
    character(len=*), parameter :: tags_opt(1) = (/ character(len=80) :: &
      OQP_afqmc_options /)
    call data_has_tags(infos%dat, tags_opt, module_name, "read_options", WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_afqmc_options, opts)
    if (size(opts) /= AFQMC_NOPT) &
      call show_message("afqmc: OQP::afqmc_options has the wrong length", WITH_ABORT)
  end subroutine read_options

!> @brief Map the option vector onto AFQMC_Control_Type and the trial controls.
  subroutine unpack_controls(opts, cntrl, trial_kind, oo_orb, nlow, &
                             low_max, low_cap, low_start, trial_state, trial_thresh)
    use AFQMC_Module, only: AFQMC_Control_Type
    use messages, only: show_message, WITH_ABORT
    real(kind=dp), pointer, intent(in) :: opts(:)
    type(AFQMC_Control_Type), intent(out) :: cntrl
    integer, intent(out) :: trial_kind, oo_orb, nlow, low_start, trial_state
    real(kind=dp), intent(out) :: low_max, low_cap, trial_thresh

    interface
      subroutine AFQMC_Init(C)
        use AFQMC_Module, only: AFQMC_Control_Type
        implicit none
        type(AFQMC_Control_Type), intent(inout) :: C
      end subroutine
    end interface

    call AFQMC_Init(cntrl)
    trial_kind = 0
    oo_orb = 0
    nlow = 0
    low_max = 0.0_dp
    low_cap = 0.0_dp
    low_start = 0
    trial_state = 1
    trial_thresh = 1.0e-8_dp

    if (associated(opts)) then
      trial_kind          = nint(opts(OPT_TRIAL_KIND))
      cntrl%Chol_Tol      = opts(OPT_CHOL_TOL)
      cntrl%NWlk          = nint(opts(OPT_NWALKERS))
      cntrl%NStep         = nint(opts(OPT_NSTEPS))
      cntrl%dT            = opts(OPT_TIMESTEP)
      cntrl%ISeed         = nint(opts(OPT_SEED))
      cntrl%NStep_Stblz   = nint(opts(OPT_STABILIZE_EVERY))
      cntrl%NStep_PopCntrl = nint(opts(OPT_POPCONTROL_EVERY))
      cntrl%NStep_Estimate = nint(opts(OPT_ESTIMATE_EVERY))
      cntrl%NStep_Accum   = nint(opts(OPT_ACCUMULATE_AFTER))
      cntrl%x_bound       = opts(OPT_FORCE_BIAS_BOUND)
      oo_orb              = nint(opts(OPT_OO_ORBITALS))
      nlow                = nint(opts(OPT_NLOW))
      low_max             = opts(OPT_LOW_MAX)
      low_cap             = opts(OPT_LOW_CAP)
      low_start           = nint(opts(OPT_LOW_START))
      trial_state         = nint(opts(OPT_TRIAL_STATE))
      trial_thresh        = opts(OPT_TRIAL_THRESH)
    end if

!   E_Bound tracks dT, so it must be recomputed after dT is overridden.
    cntrl%E_Bound = sqrt(2.0_dp/cntrl%dT)

    if (cntrl%NWlk < 1 .or. cntrl%NStep < 1 .or. cntrl%dT <= 0.0_dp) &
      call show_message("afqmc: nwalkers, nsteps and timestep must be positive", WITH_ABORT)
    if (min(cntrl%NStep_Stblz, cntrl%NStep_PopCntrl, cntrl%NStep_Estimate) < 1) &
      call show_message("afqmc: the *_every intervals must be positive", WITH_ABORT)
    if (trial_kind < TRIAL_MEAN_FIELD .or. trial_kind > TRIAL_MRSF_STATE) &
      call show_message("afqmc: trial must be mean_field, sf, mrsf or mrsf_state", &
                        WITH_ABORT)
    if (trial_state < 1) &
      call show_message("afqmc: state must be a 1-based MRSF root index", WITH_ABORT)
    if (nlow > 0 .and. trial_kind == TRIAL_MEAN_FIELD) &
      call show_message("afqmc: nlow needs an SF/MRSF-CIS determinant trial", WITH_ABORT)
  end subroutine unpack_controls

!> @brief Determinant-expansion file each trial kind reads from the run directory.
  function trial_file_name(trial_kind) result(name)
    integer, intent(in) :: trial_kind
    character(len=8) :: name
    select case (trial_kind)
    case (1)
      name = 'SFDAT   '
    case (2)
      name = 'MRSFDAT '
    case default
      name = '        '
    end select
  end function trial_file_name

!> @brief Pivoted incomplete Cholesky of the AO ERI supermatrix.
!>
!> Factorizes (mu nu|la si) = sum_g L^g_{mu nu} L^g_{la si}, keeping vectors
!> until the largest residual diagonal falls below @p tol.
!>
!> Only the residual *diagonal* and the already-accepted vectors are carried,
!> so the cost is O(nbf^2 * nchol^2) with the inner update done by DGEMV -- the
!> supermatrix is never copied and never updated in full.  (The Python bridge
!> this replaces diagonalized the whole nbf^2 x nbf^2 supermatrix instead,
!> which is O(nbf^6).)
!>
!> Because the supermatrix is positive semidefinite, the largest residual
!> diagonal at exit bounds every residual element, so @p residual is a rigorous
!> bound on the reconstruction error -- no O(nbf^4 * nchol) check needed.
  subroutine pivoted_cholesky(eri, nbf, tol, vectors, nchol, residual)
    use messages, only: show_message, WITH_ABORT
    real(kind=dp), contiguous, pointer, intent(in) :: eri(:)
    integer, intent(in) :: nbf
    real(kind=dp), intent(in) :: tol
    real(kind=dp), allocatable, intent(out) :: vectors(:,:)  !< (nbf*nbf, nchol)
    integer, intent(out) :: nchol
    real(kind=dp), intent(out) :: residual

    real(kind=dp), allocatable :: diag(:), col(:), pivot_row(:), grown(:,:)
    integer :: npair, capacity, p, pmax, ig
    real(kind=dp) :: vmax, scale

    npair = nbf*nbf
    allocate(diag(npair), col(npair))

!   The supermatrix diagonal is (mu nu|mu nu), i.e. eri at the repeated pair.
    do p = 1, npair
      diag(p) = eri(p + (p-1)*npair)
    end do

    capacity = min(npair, 10*nbf + 50)
    allocate(vectors(npair, capacity), pivot_row(capacity))
    nchol = 0

    do
      pmax = maxloc(diag, dim=1)
      vmax = diag(pmax)
      if (vmax <= tol) exit
      if (vmax < 0.0_dp) &
        call show_message("afqmc: AO ERI supermatrix is not positive semidefinite", &
                          WITH_ABORT)
      if (nchol == npair) exit

      if (nchol == capacity) then
!       Grow geometrically; the final count is not known in advance.
        allocate(grown(npair, min(npair, 2*capacity)))
        grown(:,1:capacity) = vectors
        call move_alloc(grown, vectors)
        capacity = size(vectors, 2)
        deallocate(pivot_row)
        allocate(pivot_row(capacity))
      end if

!     Residual column at the pivot: M(:,pmax) - sum_{g<=nchol} L^g_pmax L^g(:)
      col(1:npair) = eri((pmax-1)*npair + 1 : pmax*npair)
      if (nchol > 0) then
        pivot_row(1:nchol) = vectors(pmax, 1:nchol)
        call dgemv('N', npair, nchol, -1.0_dp, vectors, npair, pivot_row, 1, &
                   1.0_dp, col, 1)
      end if

      nchol = nchol + 1
      scale = 1.0_dp/sqrt(vmax)
      vectors(:, nchol) = col*scale
      diag = diag - vectors(:, nchol)**2
!     Rounding can drive an exhausted direction slightly negative; clamp so the
!     pivot search cannot latch onto noise.
      where (diag < 0.0_dp) diag = 0.0_dp
    end do

    residual = maxval(diag)
    if (nchol == 0) &
      call show_message("afqmc: chol_tol discarded every Cholesky vector", WITH_ABORT)

!   Hand back exactly the accepted vectors.
    allocate(grown(npair, nchol))
    grown = vectors(:,1:nchol)
    call move_alloc(grown, vectors)

    deallocate(diag, col, pivot_row)
  end subroutine pivoted_cholesky

end module afqmc_oqp_mod
