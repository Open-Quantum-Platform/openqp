module qmmm_mod

  use oqp_linalg
  use resp_mod, only: add_atom_grid
  implicit none

  character(len=*), parameter :: module_name = "qmmm_mod"

  private
  public get_mm_energy
  public form_esp_charges
  public print_mm_energy
  public oqp_esp_qmmm
  public grad_esp_qmmm
  public espf_op_corr
  public add_potqm_contributions
  public form_esp_charges_excited

  contains

  subroutine espf_op_corr_C(c_handle) bind(C, name="espf_op_corr")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
!    call form_esp_charges(inf)
    call espf_op_corr(inf)
  end subroutine espf_op_corr_C

  subroutine form_esp_charges_C(c_handle) bind(C, name="form_esp_charges")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call form_esp_charges(inf)
  end subroutine form_esp_charges_C

  subroutine grad_esp_qmmm_C(c_handle) bind(C, name="grad_esp_qmmm")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call grad_esp_qmmm(inf)
  end subroutine grad_esp_qmmm_C

  subroutine grad_esp_qmmm_excited_C(c_handle) bind(C, name="grad_esp_qmmm_excited")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call grad_esp_qmmm_excited(inf)
  end subroutine grad_esp_qmmm_excited_C

  subroutine form_esp_charges_excited(infos)
    use precision, only: dp
    use io_constants, only: iw
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use oqp_tagarray_driver
    use int1, only: omp_qmmm
    use constants, only: tol_int
    use mathlib, only: traceprod_sym_packed
    implicit none

    character(len=*), parameter :: module_name = "form_esp_charges_excited"
    character(len=*), parameter :: subroutine_name = "form_esp_charges_excited"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    integer :: nat, nbf, nbf2, ok, istate
    integer :: npt, nptcur, i
    logical :: urohf, use_relaxed

    real(dp), allocatable :: tmp(:), chg_op(:)
    real(dp), allocatable, target :: xyz(:,:), ttt(:,:)
    real(dp), pointer :: wn(:,:)
    real(dp) :: tol, chg, corr

    real(dp), contiguous, pointer :: dmat_a(:), dmat_b(:), partial_charges(:)
    real(dp), contiguous, pointer :: td_p(:,:), td_abxc(:)

    character(len=*), parameter :: tags_required(5) = (/ character(len=80) :: &
         OQP_DM_A, OQP_DM_B, OQP_TD_P, OQP_TD_ABXC, OQP_SM /)

    character(len=*), parameter :: tags_qmmm(1) = (/ character(len=80) :: &
         OQP_partial_charges /)

    open (unit=IW, file=infos%log_filename, position="append")

    basis => infos%basis
    basis%atoms => infos%atoms

    nat  = ubound(infos%atoms%zn, 1)
    nbf  = basis%nbf
    nbf2 = nbf * (nbf + 1) / 2
    istate = infos%tddft%target_state

    urohf = infos%control%scftype == 2 .or. infos%control%scftype == 3
    use_relaxed = .true.

    ! ESPF_ROHF=1: use the ROHF reference density for ESPF charge fitting instead
    ! of the S1 relaxed density.  This matches GAMESS's ESPF implementation, which
    ! always fits charges from the ROHF reference (not the response density).  With
    ! the hard-pruned GAMESS grid (ESPF_GAMESS=1), the ROHF density has smaller
    ! ESPF fitting residuals at the grid boundary, so the force discontinuity when
    ! points blink in/out is much smaller -- reproducing GAMESS-level conservation.
    block
      character(len=8) :: env_r
      integer :: st_r
      call get_environment_variable('ESPF_ROHF', env_r, status=st_r)
      if (st_r == 0) then
        if (trim(env_r) == '1' .or. trim(env_r) == 'on') use_relaxed = .false.
      end if
    end block

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)

    allocate(tmp(nbf2), stat=ok)
    if (ok /= 0) call show_message('Cannot allocate tmp', WITH_ABORT)
    tmp = 0.0_dp

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    tmp = dmat_a
    if (urohf) then
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      tmp = tmp + dmat_b
    end if

    if (use_relaxed) then
      call tagarray_get_data(infos%dat, OQP_TD_P, td_p)
      tmp = tmp + td_p(:,1) + td_p(:,2)
      write(iw,'(4x,a,i4)') 'Using RELAXED excited-state density for ESPF charges, state ', istate
    else
      write(iw,'(4x,a,i4)') 'Using ROHF reference density for ESPF charges (ESPF_ROHF), state ', istate
    end if

    npt = nat * (132 + 152 + 192 + 350)
    allocate(xyz(npt,3), ttt(nat,npt), chg_op(nbf2), stat=ok)
    if (ok /= 0) call show_message('Cannot allocate ESPF arrays', WITH_ABORT)

    call form_espf_grid(nat, npt, 4, [1.4_dp,1.6_dp,1.8_dp,2.0_dp], &
                        [132,152,192,350], [3,3,3,0], &
                        infos%atoms%zn, infos%atoms%xyz, xyz, ttt, nptcur)

    wn  => ttt(:,:nptcur)
    tol = log(10.0d0) * tol_int

    call infos%dat%reserve_data(OQP_partial_charges, TA_TYPE_REAL64, nat, &
         comment=OQP_partial_charges_comment)
    partial_charges = 0.0_dp
    call data_has_tags(infos%dat, tags_qmmm, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_partial_charges, partial_charges)

    corr=0.0_dp
    do i = 1, nat
      call omp_qmmm(basis, i, transpose(xyz(1:nptcur,:)), wn, chg_op, nat, logtol=tol)
      chg = traceprod_sym_packed(tmp, chg_op, nbf)
      partial_charges(i) = infos%atoms%zn(i) + chg
      corr = corr + partial_charges(i)
    end do

    do i = 1, nat
      partial_charges(i) = partial_charges(i) + (infos%mol_prop%charge-corr)/nat
    end do

    call print_charges(infos, partial_charges, iw)

    deallocate(tmp, xyz, ttt, chg_op)
  end subroutine form_esp_charges_excited


  subroutine espf_op_corr(infos)!,dmat,nbf)
    use precision, only: dp
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring
    use int1, only: omp_qmmm
    use constants, only: tol_int
    use mathlib, only: traceprod_sym_packed
    implicit none

    character(len=*), parameter :: subroutine_name = "espf_op_corr"

    integer :: nbf!, intent(in) :: nbf

    type(information), target, intent(inout) :: infos

    integer :: nat,nelec,nbf2,ok
    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: chg_op(:)
    real(kind=dp), target, allocatable :: ttt(:,:), xyz(:,:)
    real(kind=dp), pointer :: coord(:,:)
    real(kind=dp), pointer :: wn(:,:)

    integer :: npt, nptcur
    integer :: i
    real(kind=dp) :: chg

    integer, parameter :: nlayers = 4
    real(kind=dp), parameter :: &
      layers(nlayers) = [1.4, 1.6, 1.8, 2.0]
    integer, parameter :: npt_layer(nlayers) = [132,152,192,350]
    integer, parameter :: typ_layer(nlayers) = [3,3,3,0]
    real(kind=dp), allocatable :: sum_op(:), corr(:)

    real(kind=dp) :: tol

    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    character(len=*), parameter :: tags_alpha(1) = &
      (/ character(len=80) :: OQP_DM_A /)
    character(len=*), parameter :: tags_beta(1) = &
      (/ character(len=80) :: OQP_DM_B /)
    !==============================================================================
    ! Tag Arrays for Accessing Data
    !==============================================================================
    real(kind=dp), contiguous, pointer :: smat(:), chg_ops_corr(:,:)
    character(len=*), parameter :: tags_general(1) = &
      (/ character(len=80) :: OQP_SM /)
    character(len=*), parameter :: tags(1) = &
      (/ character(len=80) :: OQP_ESPF_CORR /)
!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
!   Allocate memory
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(infos%atoms%zn,1)
    nelec = infos%mol_prop%nelec
    npt = nat * sum(npt_layer)

    allocate(xyz(npt,3), &
             ttt(nat,npt), &
             chg_op(nbf2), &
             sum_op(nbf2), &
             corr(nbf2), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)

    call infos%dat%reserve_data(OQP_ESPF_CORR, TA_TYPE_REAL64, nbf2*nat, (/ nbf2, nat /), comment=OQP_ESPF_CORR_comment)
    call data_has_tags(infos%dat, tags, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_ESPF_CORR, chg_ops_corr)

    call form_espf_grid(nat,npt,nlayers,layers,npt_layer,typ_layer,infos%atoms%zn,infos%atoms%xyz,xyz,ttt,nptcur)

    wn => ttt(:,:nptcur)
    tol = log(10.0d0)*tol_int
    sum_op = 0
    do i = 1, nat
       chg_ops_corr(:,i) = 0.0_dp
       call omp_qmmm( &
            basis,                        &
            i,                            &
            transpose(xyz(1:nptcur,:)),   &
            wn,                           &
            chg_ops_corr(:,i),            &
            nat,                          &
            logtol = tol )

       sum_op(:) = sum_op(:) + chg_ops_corr(:,i)
    end do

    corr(:) = (sum_op(:) + smat(:)) / real(nat, dp)

    do i = 1, nat
       chg_ops_corr(:,i) = chg_ops_corr(:,i) - corr(:)
    end do
    deallocate(xyz, ttt, chg_op, sum_op, corr)

  end subroutine espf_op_corr

  subroutine add_potqm_contributions(infos, dens, dh)
    use precision,  only : dp
    use basis_tools, only : basis_set
    use types,      only : information
    use messages,   only : show_message, with_abort
    use mathlib,    only : traceprod_sym_packed
    use oqp_tagarray_driver
    implicit none

    character(len=*), parameter :: subroutine_name = "add_potqm_contributions"
    character(len=*), parameter :: module_name     = "qmmm_espf"  ! adjust if needed

    type(information), target, intent(inout) :: infos
    real(dp), contiguous, intent(in)        :: dens(:)
    real(dp), contiguous, intent(inout)     :: dh(:)

    type(basis_set), pointer :: basis
    integer :: nbf, nbf2, nat, i, ok

    real(dp), contiguous, pointer :: smat(:), potqm(:,:), chg_ops_corr(:,:)

    real(dp), allocatable :: q(:)
    real(dp), allocatable :: v(:)
!    real(dp), allocatable :: dh(:)

    character(len=*), parameter :: tags(3) = &
         (/ character(len=80) :: OQP_SM, OQP_ESPF_CORR, OQP_POTQM /)

    if(.not.infos%control%qmmm_flag) return

    basis      => infos%basis
    basis%atoms => infos%atoms

    nbf  = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat  = ubound(infos%atoms%zn, 1)

    call data_has_tags(infos%dat, tags, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM,        smat)
    call tagarray_get_data(infos%dat, OQP_ESPF_CORR, chg_ops_corr)
    call tagarray_get_data(infos%dat, OQP_POTQM,     potqm)

    allocate(q(nat), v(nat), stat=ok)
    if (ok /= 0) call show_message("Cannot allocate memory in add_potqm_contributions", WITH_ABORT)

    q(:) = 0.0_dp
    do i = 1, nat
       q(i) = traceprod_sym_packed(dens, chg_ops_corr(:,i), nbf)
    end do

    v = matmul(potqm, q)

    dh(:) = 0.0_dp
    do i = 1, nat
       dh(:) = dh(:) - v(i) * chg_ops_corr(:,i)
    end do
    deallocate(q, v)

  end subroutine add_potqm_contributions


!> @brief Compute the MM and atomic QM/MM contributions to energy
!
!> @detail Classical contributions to energy and the QM atomic charge times MM potential
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Oct, 2024_ Initial release
  subroutine get_mm_energy(infos,emm)
    use precision, only: dp
    use oqp_tagarray_driver
    use types, only: information
    use messages, only: WITH_ABORT
    implicit none

    character(len=*), parameter :: subroutine_name = "get_mm_energy"

    real(kind=dp), intent(inout) :: emm
    real(kind=dp), contiguous, pointer :: mm_potential(:), mm_energy(:)
    character(len=*), parameter :: tags_qmmm(2) = (/ character(len=80) :: &
       OQP_mm_potential, OQP_mm_energy /)

    type(information), target, intent(inout) :: infos

    emm = 0.0d0
    if(.not.infos%control%qmmm_flag) return

    call data_has_tags(infos%dat, tags_qmmm, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_mm_energy, mm_energy)
    call tagarray_get_data(infos%dat, OQP_mm_potential, mm_potential)
    emm=sum(mm_energy)+dot_product(mm_potential,infos%atoms%zn)

    return

  end subroutine get_mm_energy

!--------------------------------------------------------------------------------

!> @brief Print MM energy decomposition in file
!
!> @detail Classical contributions to energy and the QM atomic charge times MM potential
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Oct, 2024_ Initial release
  subroutine print_mm_energy(infos)
    use precision, only: dp
    use oqp_tagarray_driver
    use types, only: information
    use messages, only: WITH_ABORT
    use io_constants, only: iw
    implicit none

    character(len=*), parameter :: subroutine_name = "get_mm_energy"

    real(kind=dp), contiguous, pointer :: mm_potential(:), mm_energy(:),partial_charges(:)
    character(len=*), parameter :: tags_qmmm(3) = (/ character(len=80) :: &
       OQP_mm_potential, OQP_mm_energy, OQP_partial_charges /)

    type(information), target, intent(inout) :: infos

    if(.not.infos%control%qmmm_flag) return

    call data_has_tags(infos%dat, tags_qmmm, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_mm_energy, mm_energy)
    call tagarray_get_data(infos%dat, OQP_mm_potential, mm_potential)
    call tagarray_get_data(infos%dat, OQP_partial_charges, partial_charges)

    write(iw,*)
! QM/MM interaction energy
    write(iw,"(' QM/MM: Electrostatic energy      = ',F20.10)")&
    dot_product(partial_charges,mm_potential)
! Classical forcefields
    if(abs(mm_energy(1)).gt.1.0d-10)&
       write(iw,"(' MM: Bonded energy                = ',F20.10)") mm_energy(1)
    if(abs(mm_energy(2)).gt.1.0d-10)&
       write(iw,"(' MM: Angle bending energy         = ',F20.10)") mm_energy(2)
    if(abs(mm_energy(3)).gt.1.0d-10)&
       write(iw,"(' MM: Periodic Torsion energy      = ',F20.10)") mm_energy(3)
    if(abs(mm_energy(4)).gt.1.0d-10)&
       write(iw,"(' MM: RBTorsion energy             = ',F20.10)") mm_energy(4)
    if(abs(mm_energy(5)).gt.1.0d-10)&
       write(iw,"(' MM: Nonbonded energy             = ',F20.10)") mm_energy(5)

! Custom Classical forcefields
    if(abs(mm_energy(6)).gt.1.0d-10)&
       write(iw,"(' MM: CustomBonded energy          = ',F20.10)") mm_energy(6)
    if(abs(mm_energy(7)).gt.1.0d-10)&
       write(iw,"(' MM: CustomAngle energy           = ',F20.10)") mm_energy(7)
    if(abs(mm_energy(8)).gt.1.0d-10)&
       write(iw,"(' MM: CustomTorsion energy         = ',F20.10)") mm_energy(8)
    if(abs(mm_energy(9)).gt.1.0d-10)&
       write(iw,"(' MM: CustomNonbonded energy       = ',F20.10)") mm_energy(9)
    if(abs(mm_energy(10)).gt.1.0d-10)&
       write(iw,"(' MM: CustomExternal energy        = ',F20.10)") mm_energy(10)
! Implicit solvation
    if(abs(mm_energy(11)).gt.1.0d-10)&
       write(iw,"(' MM: GBSAOBC energy               = ',F20.10)") mm_energy(11)
    if(abs(mm_energy(12)).gt.1.0d-10)&
       write(iw,"(' MM: CustomGB energy              = ',F20.10)") mm_energy(12)

! Amoeba forcefields
    if(abs(mm_energy(13)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaMultipole energy       = ',F20.10)") mm_energy(13)
    if(abs(mm_energy(14)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaVdw energy             = ',F20.10)") mm_energy(14)
    if(abs(mm_energy(15)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaBond energy            = ',F20.10)") mm_energy(15)
    if(abs(mm_energy(16)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaAngle energy           = ',F20.10)") mm_energy(16)
    if(abs(mm_energy(17)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaTorsion energy         = ',F20.10)") mm_energy(17)
    if(abs(mm_energy(18)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaOutOfPlaneBend energy  = ',F20.10)") mm_energy(18)
    if(abs(mm_energy(19)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaPiTorsion energy       = ',F20.10)") mm_energy(19)
    if(abs(mm_energy(20)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaStretchBend energy     = ',F20.10)") mm_energy(20)
    if(abs(mm_energy(21)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaUreyBradley energy     = ',F20.10)") mm_energy(21)
    if(abs(mm_energy(22)).gt.1.0d-10)&
       write(iw,"(' MM: AmoebaVdw14 energy           = ',F20.10)") mm_energy(22)
! Extras
    if(abs(mm_energy(23)).gt.1.0d-10)&
       write(iw,"(' MM: CMMotion energy              = ',F20.10)") mm_energy(23)
    if(abs(mm_energy(24)).gt.1.0d-10)&
       write(iw,"(' MM: CustomCompoundBond energy    = ',F20.10)") mm_energy(24)
    if(abs(mm_energy(25)).gt.1.0d-10)&
       write(iw,"(' MM: MonteCarloBarostad           = ',F20.10)") mm_energy(25)
    write(iw,*)

    call print_charges(infos,partial_charges,iw)

    return

  end subroutine print_mm_energy

!--------------------------------------------------------------------------------

!> @brief Compute the MM and atomic QM/MM contributions to energy
!
!> @detail Classical contributions to energy and the QM atomic charge times MM potential
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Oct, 2024_ Initial release

  subroutine form_esp_charges(infos)!,dmat,nbf)
    use precision, only: dp
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring
    use int1, only: omp_qmmm
    use constants, only: tol_int
    use mathlib, only: traceprod_sym_packed
    implicit none

    character(len=*), parameter :: subroutine_name = "form_esp_charges"

    integer :: nbf!, intent(in) :: nbf
    real(kind=dp), allocatable  :: dmat(:)!(nbf*nbf)!, intent(in) :: dmat(nbf*nbf)

    type(information), target, intent(inout) :: infos

    real(kind=dp), contiguous, pointer :: partial_charges(:)
    character(len=*), parameter :: tags_qmmm(1) = (/ character(len=80) :: &
       OQP_partial_charges /)

    integer :: nat,nelec,nbf2,ok
    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: chg_op(:)
    real(kind=dp), target, allocatable :: ttt(:,:), xyz(:,:)
    real(kind=dp), pointer :: coord(:,:)
    real(kind=dp), pointer :: wn(:,:)

    integer :: npt, nptcur
    integer :: i
    real(kind=dp) :: chg,corr

    integer, parameter :: nlayers = 4
    real(kind=dp), parameter :: &
      layers(nlayers) = [1.4, 1.6, 1.8, 2.0]
    integer, parameter :: npt_layer(nlayers) = [132,152,192,350]
    integer, parameter :: typ_layer(nlayers) = [3,3,3,0]

    real(kind=dp) :: tol

    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    character(len=*), parameter :: tags_alpha(1) = &
      (/ character(len=80) :: OQP_DM_A /)
    character(len=*), parameter :: tags_beta(1) = &
      (/ character(len=80) :: OQP_DM_B /)


!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
!   Allocate memory
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(infos%atoms%zn,1)
    nelec = infos%mol_prop%nelec
    npt = nat * sum(npt_layer)


    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)

    ! Get beta-spin tag arrays if needed
    if (infos%control%scftype > 1) then
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    end if

    allocate(xyz(npt,3), &
             ttt(nat,npt), &
             chg_op(nbf2), &
             dmat(nbf2), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    dmat = dmat_a
    if (infos%control%scftype > 1) then
      dmat = dmat + dmat_b
    end if

    call form_espf_grid(nat,npt,nlayers,layers,npt_layer,typ_layer,infos%atoms%zn,infos%atoms%xyz,xyz,ttt,nptcur)

! Compute integrals and form partial charges
    wn => ttt(:,:nptcur)
    tol = log(10.0d0)*tol_int

    call infos%dat%reserve_data(OQP_partial_charges, TA_TYPE_REAL64, infos%mol_prop%natom, comment=OQP_partial_charges_comment)
    call data_has_tags(infos%dat, tags_qmmm, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_partial_charges, partial_charges)

! Form partial charges and compute the correction for total charge
    corr=0.0_dp!nelec/nat

    do i=1,nat
       call omp_qmmm(basis, i, transpose(xyz(1:nptcur,:)), wn, chg_op, nat, logtol=tol)
       chg = traceprod_sym_packed(dmat,chg_op,nbf)
       partial_charges(i) = infos%atoms%zn(i) + chg
       corr = corr + partial_charges(i)
    end do

! Correct for the total charge
    do i=1,nat
       partial_charges(i) = partial_charges(i) + (infos%mol_prop%charge - corr)/nat
    end do
    deallocate(xyz, ttt, chg_op, dmat)

  end subroutine form_esp_charges

!--------------------------------------------------------------------------------

!> @brief Compute the QM/MM one-electron hamiltonian using ESPF method
!> @param[in]      infos        OQP handle
!> @param[in,out]  hqmmm        QM/MM one-electron hamiltonian
!> @param[in]      mm_potential classical MM potential on QM centers
!> @param[in]      smat         overlap matrix
!> @param[in]      logtol       tolerance threshold for integrals
!
!> @detail This subroutine computes the one-electron hamiltonian for introducing the QM/MM interaction
!>         in the core hamiltonian using the electrostatic potential fitted (ESPF) method.
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Jul, 2024_ Initial release
  subroutine oqp_esp_qmmm(infos, Hqmmm, mm_potential, smat, logtol)
    use precision, only: dp
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring
    use int1, only: omp_qmmm
    use lebedev, only: lebedev_get_grid
    use elements, only: ELEMENTS_VDW_RADII
    implicit none

    character(len=*), parameter :: subroutine_name = "oqp_esp_qmmm"

    type(information), target, intent(inout) :: infos
    real(kind=dp), contiguous, intent(inout) :: Hqmmm(:)
    real(kind=dp), contiguous, intent(in) :: smat(:),mm_potential(:)
    real(kind=dp), intent(in) :: logtol

    integer :: nbf, nbf2, ok
    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: wt(:), chg_op(:)
    real(kind=dp), allocatable :: xyz(:,:)
    real(kind=dp), target, allocatable :: ttt(:,:)
    real(kind=dp), pointer :: coord(:,:)
    real(kind=dp), pointer :: wn(:,:)

    integer :: nat, npt, nptcur
    integer :: i
    real(kind=dp) :: mm_pot_av

    integer, parameter :: nlayers = 4
    real(kind=dp), parameter :: &
      layers(nlayers) = [1.4, 1.6, 1.8, 2.0]
    integer, parameter :: npt_layer(nlayers) = [132,152,192,350]
    integer, parameter :: typ_layer(nlayers) = [3,3,3,0]

    logical :: restr

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

!   Allocate memory
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(infos%atoms%zn, 1)
    npt = nat * sum(npt_layer)

    allocate(xyz(npt,3), &
             wt(npt), &
             ttt(nat,npt), &
             chg_op(nbf2), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

! Compute the numerical Lebedev grid (xyz) and integral weights (TTT)
    call form_espf_grid(nat,npt,nlayers,layers,npt_layer,typ_layer,infos%atoms%zn,infos%atoms%xyz,xyz,ttt,nptcur)

! Compute integrals and form QM/MM hamiltonian
    wn => ttt(:,:nptcur)
    mm_pot_av = sum(mm_potential)/nat
    hqmmm = -smat*mm_pot_av
    do i=1,nat
       call omp_qmmm(basis, i, transpose(xyz(1:nptcur,:)), wn, chg_op, nat, logtol)
       hqmmm = hqmmm + (mm_potential(i)-mm_pot_av)*chg_op
    end do

  end subroutine oqp_esp_qmmm

!--------------------------------------------------------------------------------

!> @brief Compute the gradient contributions from the QM/MM interaction energy using ESPF method
!> @param[in]      infos        OQP handle
!> @param[in]      dens         density matrix
!> @param[in,out]  grad         energy gradient
!> @param[in]      logtol       tolerance threshold for integrals
!
!> @detail This subroutine computes the analytic derivatives of the QM/MM interaction energy
!>         from the electrostatic potential fitted (ESPF) method. Here the polarizable term
!>         (Q^x*mm_potential) and the classical term (Q*mm_potential^x) computed by OpenMM
!>         are added to the QM gradient. The MM atoms gradient is added in the python layer.
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Jul, 2024_ Initial release
  subroutine grad_esp_qmmm(infos)!, dens, grad, logtol)
    use oqp_tagarray_driver
    use precision, only: dp
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring
    use grd1, only: grad_elpot, grad_ee_overlap
    use lebedev, only: lebedev_get_grid
    use elements, only: ELEMENTS_VDW_RADII
    use constants, only: tol_int
    implicit none

    character(len=*), parameter :: subroutine_name = "grad_esp_qmmm"

    type(information), target, intent(inout) :: infos
!    real(kind=dp), intent(inout) :: grad(:,:)
!    real(kind=dp), intent(inout) :: dens(:)
    real(kind=dp) :: logtol

    integer :: nbf, nbf2, ok
    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: wt(:), dens(:)
    real(kind=dp), allocatable :: xyz(:,:)
    real(kind=dp), target, allocatable :: ttt(:,:)

    integer :: nat, npt, nptcur
    integer :: i, j, k

    integer, parameter :: nlayers = 4
    real(kind=dp), parameter :: &
      layers(nlayers) = [1.4, 1.6, 1.8, 2.0]
    integer, parameter :: npt_layer(nlayers) = [132,152,192,350]
    integer, parameter :: typ_layer(nlayers) = [3,3,3,0]

!tagarray
    real(kind=dp), contiguous, pointer :: partial_charges(:), mm_potential(:),&
                                          espf_grad_ta(:,:), &
                                          dmat_a(:), dmat_b(:)
    character(len=*), parameter :: tags_qmmm(4) = (/ character(len=80) :: &
      OQP_partial_charges, OQP_POTMM, OQP_ESPF_GRAD, OQP_DM_A/)

    character(len=*), parameter :: tags_beta(1) = (/ character(len=80) :: &
      OQP_DM_B/)

    real(kind=dp) :: mm_pot_av
    logical :: skip_legacy
    character(len=8) :: env_l
    integer :: st_l

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

!   Replace legacy espf_grad weight term with the complete pseudoinverse
!   derivative (espf_grad_weight); ESPF_LEGACY=1 restores the old term.
    skip_legacy = .true.
    call get_environment_variable('ESPF_LEGACY', env_l, status=st_l)
    if (st_l == 0) then
      if (trim(env_l) == '1' .or. trim(env_l) == 'on') skip_legacy = .false.
    end if

!   Allocate memory
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(infos%atoms%zn, 1)
    npt = nat * sum(npt_layer)

    logtol = log(10.0d0)*tol_int
    allocate(xyz(npt,3), &
             wt(npt), &
             ttt(nat,npt), &
             dens(nbf2), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    call form_espf_grid(nat,npt,nlayers,layers,npt_layer,typ_layer,infos%atoms%zn,infos%atoms%xyz,xyz,ttt,nptcur)

! ESP gradient contribution
!   Tagarray
    call infos%dat%reserve_data(OQP_ESPF_GRAD, TA_TYPE_REAL64, nat*3, (/ 3, nat /), comment=OQP_ESPF_GRAD_comment)
    call data_has_tags(infos%dat, tags_qmmm, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_POTMM, mm_potential)
    call tagarray_get_data(infos%dat, OQP_partial_charges, partial_charges)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_ESPF_GRAD, espf_grad_ta)
   ! Get beta-spin tag arrays if needed
    if (infos%control%scftype > 1) then
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    end if

    espf_grad_ta = 0.0_dp
    dens = dmat_a
    if (infos%control%scftype>=2) dens =  dens + dmat_b
! Compute integrals and form ESP operators
!   Compute the corrected mm potential
    mm_pot_av = sum(mm_potential)/nat
    do i=1,nat
       mm_potential(i)=mm_pot_av-mm_potential(i)
    end do

!   Add integral gradient term, mm_potential*[(T^+T)^-1*T^+]*V^x
    do i=1,nptcur
       wt(i)=-dot_product(ttt(:,i),mm_potential)
       call grad_elpot(basis, xyz(i,:), wt(i), dens, espf_grad_ta)
    end do

!   Add overlap derivative correction for the total charge conservation
    if(abs(mm_pot_av).gt.1.0e-6) then
       dens=-mm_pot_av*dens
       call grad_ee_overlap(basis, dens, espf_grad_ta, logtol)
       dens=-dens/mm_pot_av
    end if

!   Add weights gradient term, -mm_potential*[(T^+T)^-1*T^+]*T^xQ + q*mm_potential^x
    if (.not. skip_legacy) then
    call espf_grad(&
         x=xyz(:nptcur,1),&
         y=xyz(:nptcur,2),&
         z=xyz(:nptcur,3),&
         at=infos%atoms%xyz,&
         wt=wt,&
         zn=infos%atoms%zn,&
         pchg=partial_charges,&
         grad=espf_grad_ta)
    end if
!    espf_grad = grad

!   Complete pseudoinverse-weight derivative dZ/dR (a_i = phi_i-<phi> = -mm_potential)
    call espf_grad_weight(basis, nat, nptcur, infos%atoms%xyz, xyz, &
                          ttt, -mm_potential(1:nat), dens, espf_grad_ta, logtol, &
                          ELEMENTS_VDW_RADII(int(infos%atoms%zn)))

!   Restablish mm potential to the original one
    do i=1,nat
! it seems not to be:  mm_potential(i)= - mm_pot_av-mm_potential(i)
       mm_potential(i)= mm_pot_av-mm_potential(i)
    end do

  end subroutine grad_esp_qmmm

  subroutine grad_esp_qmmm_excited(infos)!, dens, grad, logtol)
    use oqp_tagarray_driver
    use precision, only: dp
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring
    use grd1, only: grad_elpot, grad_ee_overlap
    use lebedev, only: lebedev_get_grid
    use elements, only: ELEMENTS_VDW_RADII
    use constants, only: tol_int
    implicit none

    character(len=*), parameter :: subroutine_name = "grad_esp_qmmm_excited"

    type(information), target, intent(inout) :: infos
!    real(kind=dp), intent(inout) :: grad(:,:)
!    real(kind=dp), intent(inout) :: dens(:)
    real(kind=dp) :: logtol

    integer :: nbf, nbf2, ok
    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: wt(:), dens(:)
    real(kind=dp), allocatable :: xyz(:,:)
    real(kind=dp), target, allocatable :: ttt(:,:)

    integer :: nat, npt, nptcur
    integer :: i, j, k

    integer, parameter :: nlayers = 4
    real(kind=dp), parameter :: &
      layers(nlayers) = [1.4, 1.6, 1.8, 2.0]
    integer, parameter :: npt_layer(nlayers) = [132,152,192,350]
    integer, parameter :: typ_layer(nlayers) = [3,3,3,0]

!tagarray
    real(kind=dp), contiguous, pointer :: partial_charges(:), mm_potential(:),&
                                          espf_grad_ta(:,:), &
                                          dmat_a(:), dmat_b(:)
    real(dp), contiguous, pointer :: td_p(:,:), td_abxc(:)
    character(len=*), parameter :: tags_qmmm(6) = (/ character(len=80) :: &
      OQP_partial_charges, OQP_POTMM, OQP_ESPF_GRAD, OQP_DM_A, &
      OQP_TD_P, OQP_TD_ABXC/)
    character(len=*), parameter :: tags_beta(1) = (/ character(len=80) :: &
      OQP_DM_B/)

    real(kind=dp) :: mm_pot_av
    logical :: use_relaxed
    logical :: skip_legacy
    character(len=8) :: env_l
    integer :: st_l

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms
    use_relaxed = .true.
    ! ESPF_ROHF=1: use ROHF reference density in the ESPF gradient.  The ROHF
    ! density has smaller ESP fitting residuals at grid-boundary points; when the
    ! GAMESS hard-pruned grid (ESPF_GAMESS=1) removes those points the force
    ! discontinuity is therefore much smaller, reproducing GAMESS energy conservation.
    block
      character(len=8) :: env_r2
      integer :: st_r2
      call get_environment_variable('ESPF_ROHF', env_r2, status=st_r2)
      if (st_r2 == 0) then
        if (trim(env_r2) == '1' .or. trim(env_r2) == 'on') use_relaxed = .false.
      end if
    end block
    ! The legacy espf_grad weight term is replaced by the complete pseudoinverse
    ! derivative (espf_grad_weight, ported from GAMESS DVESPF/INIDZ): FD-verified
    ! to cut the QM-atom force-energy inconsistency from RMS 1.6e-3 to 5.2e-4
    ! Ha/bohr at large displacement.  Set ESPF_LEGACY=1 to restore the old term.
    skip_legacy = .true.
    call get_environment_variable('ESPF_LEGACY', env_l, status=st_l)
    if (st_l == 0) then
      if (trim(env_l) == '1' .or. trim(env_l) == 'on') skip_legacy = .false.
    end if
!   Allocate memory
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(infos%atoms%zn, 1)
    npt = nat * sum(npt_layer)

    logtol = log(10.0d0)*tol_int
    allocate(xyz(npt,3), &
             wt(npt), &
             ttt(nat,npt), &
             dens(nbf2), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    call form_espf_grid(nat,npt,nlayers,layers,npt_layer,typ_layer,infos%atoms%zn,infos%atoms%xyz,xyz,ttt,nptcur)

! ESP gradient contribution
!   Tagarray
    call infos%dat%reserve_data(OQP_ESPF_GRAD, TA_TYPE_REAL64, nat*3, (/ 3, nat /), comment=OQP_ESPF_GRAD_comment)
    call data_has_tags(infos%dat, tags_qmmm, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_POTMM, mm_potential)
    call tagarray_get_data(infos%dat, OQP_partial_charges, partial_charges)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_ESPF_GRAD, espf_grad_ta)

   ! Get beta-spin tag arrays if needed
    if (infos%control%scftype > 1) then
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    end if

    espf_grad_ta = 0.0_dp
    dens = dmat_a
    if (infos%control%scftype>=2) dens =  dens + dmat_b
    if (use_relaxed) then
      ! Default: S1 relaxed density = ROHF + td_p (orbital-relaxation correction)
      call tagarray_get_data(infos%dat, OQP_TD_P, td_p)
      dens = dens + td_p(:,1) + td_p(:,2)
    else
      ! ESPF_ROHF=1: stop at the ROHF reference density -- do NOT add td_abxc.
      ! This mirrors GAMESS, which always uses the ROHF density for ESPF fitting
      ! (not the response/relaxed density).
      continue
    end if

! Compute integrals and form ESP operators
!   Compute the corrected mm potential
    mm_pot_av = sum(mm_potential)/nat
    do i=1,nat
       mm_potential(i)=mm_pot_av-mm_potential(i)
    end do

!   Add integral gradient term, mm_potential*[(T^+T)^-1*T^+]*V^x
    do i=1,nptcur
       wt(i)=-dot_product(ttt(:,i),mm_potential)
       call grad_elpot(basis, xyz(i,:), wt(i), dens, espf_grad_ta)
    end do

!   Add overlap derivative correction for the total charge conservation
    if(abs(mm_pot_av).gt.1.0e-6) then
       dens=-mm_pot_av*dens
       call grad_ee_overlap(basis, dens, espf_grad_ta, logtol)
       dens=-dens/mm_pot_av
    end if

!   Add weights gradient term, -mm_potential*[(T^+T)^-1*T^+]*T^xQ + q*mm_potential^x
    if (.not. skip_legacy) then
    call espf_grad(&
         x=xyz(:nptcur,1),&
         y=xyz(:nptcur,2),&
         z=xyz(:nptcur,3),&
         at=infos%atoms%xyz,&
         wt=wt,&
         zn=infos%atoms%zn,&
         pchg=partial_charges,&
         grad=espf_grad_ta)
    end if
!    espf_grad = grad

!   Add the ESPF weight-matrix (pseudoinverse) derivative term dZ/dR.  Here
!   mm_potential currently holds (<phi> - phi_i), so the mean-removed potential
!   a_i = phi_i - <phi> = -mm_potential.
    call espf_grad_weight(basis, nat, nptcur, infos%atoms%xyz, xyz, &
                          ttt, -mm_potential(1:nat), dens, espf_grad_ta, logtol, &
                          ELEMENTS_VDW_RADII(int(infos%atoms%zn)))

!   Restablish mm potential to the original one
    do i=1,nat
! it seems not to be:  mm_potential(i)= - mm_pot_av-mm_potential(i)
       mm_potential(i)= mm_pot_av-mm_potential(i)
    end do

  end subroutine grad_esp_qmmm_excited

!--------------------------------------------------------------------------------

!> @brief Add to the gradient the integral weight derivatives and classical MM contributions
!> @param[in]      infos        OQP handle
!> @param[in]      dens         density matrix
!> @param[in,out]  grad         energy gradient
!> @param[in]      logtol       tolerance threshold for integrals
!
!> @detail This subroutine computes the analytic derivatives of the QM/MM interaction energy
!>         corresponding to the integral weight derivatives and the classical MM contributions
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Jul, 2024_ Initial release
  subroutine espf_grad(x, y, z, at, wt, zn, pchg, grad)
    use precision, only: dp
    implicit none

    real(kind=dp), intent(in) :: x(:), y(:), z(:), at(:,:), zn(:)
    real(kind=dp), intent(inout) :: grad(:,:)
    real(kind=dp), contiguous, intent(in) :: pchg(:)
    real(kind=dp), intent(inout) :: wt(:)
!    real(kind=dp), target, allocatable :: gradient_mm(:,:)

    integer :: i, j, nat, npts, dim1,dim2

    nat = ubound(at, 2)
    npts = ubound(x, 1)


    ! Compute the matrix of the inverse distances
    do i = 1, nat
      do j = 1, npts
        grad(1,i) = grad(1,i) &
            + wt(j)*(pchg(i)-zn(i))*(at(1,i)-x(j))/norm2(at(:,i)-[x(j),y(j),z(j)])**3
        grad(2,i) = grad(2,i) &
            + wt(j)*(pchg(i)-zn(i))*(at(2,i)-y(j))/norm2(at(:,i)-[x(j),y(j),z(j)])**3
        grad(3,i) = grad(3,i) &
            + wt(j)*(pchg(i)-zn(i))*(at(3,i)-z(j))/norm2(at(:,i)-[x(j),y(j),z(j)])**3
      end do
    end do

  end subroutine espf_grad

!--------------------------------------------------------------------------------

!> @brief ESPF weight-matrix (pseudoinverse) derivative contribution to the
!>        QM-atom gradient of the QM/MM electrostatic coupling energy.
!>
!> @detail The ESPF charges are q = Z u (electronic part), Z = (T^T T)^-1 T^T,
!>         T_{i,k} = 1/|R_i - g_k| (atom i, grid point k), u_k = Tr(P V_k^AO) the
!>         electronic ESP at grid point k.  The coupling energy contains
!>         E = sum_i a_i sum_k Z_{i,k} u_k  with a_i = phi_i - <phi> (mean-removed
!>         MM potential).  grad_elpot already differentiates u_k at fixed Z; this
!>         routine adds the term from dZ/dR (the GAMESS DVESPF/INIDZ "DTZ.V"
!>         contribution), with the atom-centred grid treated as fixed (matching
!>         GAMESS).  Closed form (M = T^T T, Minv = Z Z^T):
!>           b = Minv a,  g = Z u,  G = T^T g,  B = T^T b
!>           dE/dR_{J,c} = sum_k dT(J,k,c) [ b_J (u_k - G_k) - g_J B_k ]
!>           dT(J,k,c) = (g_k - R_J)_c / |g_k - R_J|^3
!>
!> @author  ported/derived from GAMESS espf.src (INIDZ/DVESPF), 2026-06
!>
!> @detail  Smooth-switching mode (ESPF_SMOOTH=1): the grid carries smooth
!>          weights s_k (espf_smooth_s) and Z = M^-1 T diag(s), M = T diag(s) T^T.
!>          E = a^T M^-1 T diag(s) u, so dE/dR has, besides the dT term (now with
!>          an s_k factor), an extra weight-derivative term:
!>            dE/dR += sum_k ds_k/dR * B_k (u_k - G_k),
!>          with ds_k/dR_{J} = [prod_{b/=J} S_b] S'_J (1/delta) (R_J-g_k)/|g_k-R_J|.
  subroutine espf_grad_weight(basis, nat, nptcur, at, gxyz, zmat, amm, dens, grad, logtol, rvdw)
    use precision, only: dp
    use basis_tools, only: basis_set
    use int1, only: electrostatic_potential_unweighted
    use messages, only: show_message, with_abort
    implicit none
    type(basis_set), intent(inout)       :: basis
    integer, intent(in)                  :: nat, nptcur
    real(kind=dp), intent(in)            :: at(:,:)        ! (3,nat) atom coords
    real(kind=dp), intent(in)            :: gxyz(:,:)      ! (npt,3) grid coords
    real(kind=dp), intent(in)            :: zmat(:,:)      ! (nat,npt) Z (weighted)
    real(kind=dp), intent(in)            :: amm(:)         ! (nat) mean-removed MM potential a_i
    real(kind=dp), intent(inout)         :: dens(:)        ! packed density for ESP
    real(kind=dp), intent(inout)         :: grad(:,:)      ! (3,nat) accumulated dE/dR
    real(kind=dp), intent(in)            :: logtol
    real(kind=dp), intent(in)            :: rvdw(:)        ! (nat) base VDW radii

    real(kind=dp), allocatable :: gx(:), gy(:), gz(:), u(:)
    real(kind=dp), allocatable :: tmat(:,:), minv(:,:), bvec(:), gvec(:), &
                                  gbig(:), bbig(:), sk(:)
    integer, allocatable :: ipiv(:)
    real(kind=dp) :: dx, dy, dz, r2, r3, scal, fac, sw_delta, sw_scale, d, xx, &
                     sp, sig, pe, coef, wk
    integer :: i, k, j, info
    character(len=32) :: envv
    integer :: status
    logical :: smooth

    ! runtime toggle / scale (default ON, scale 1.0); allows FD calibration
    ! without rebuilds:  ESPF_WDERIV=0 disables, ESPF_WSCALE=<f> rescales.
    call get_environment_variable('ESPF_WDERIV', envv, status=status)
    if (status == 0) then
      if (trim(envv) == '0' .or. trim(envv) == 'off') return
    end if
    scal = 1.0_dp
    call get_environment_variable('ESPF_WSCALE', envv, status=status)
    if (status == 0) then
      if (len_trim(envv) > 0) read(envv,*) scal
    end if
    smooth = .true.
    call get_environment_variable('ESPF_SMOOTH', envv, status=status)
    if (status == 0) then
      if (trim(envv) == '0' .or. trim(envv) == 'off') smooth = .false.
    end if
    sw_delta = 0.7_dp
    call get_environment_variable('ESPF_SWDELTA', envv, status=status)
    if (status == 0) then
      if (len_trim(envv) > 0) read(envv,*) sw_delta
    end if
    sw_scale = 1.8_dp
    call get_environment_variable('ESPF_SWSCALE', envv, status=status)
    if (status == 0) then
      if (len_trim(envv) > 0) read(envv,*) sw_scale
    end if

    allocate(gx(nptcur), gy(nptcur), gz(nptcur), u(nptcur))
    gx = gxyz(1:nptcur,1); gy = gxyz(1:nptcur,2); gz = gxyz(1:nptcur,3)

    ! u_k = Tr(P V_k^AO): electronic ESP at the grid points (same kernel/sign as
    ! the omp_qmmm operator used to build the charges).
    call electrostatic_potential_unweighted(basis, gx, gy, gz, dens, u, logtol)

    ! T_{i,k} = 1/|R_i - g_k|
    allocate(tmat(nat, nptcur), sk(nptcur))
    do k = 1, nptcur
      do i = 1, nat
        dx = at(1,i) - gx(k); dy = at(2,i) - gy(k); dz = at(3,i) - gz(k)
        tmat(i,k) = 1.0_dp / sqrt(dx*dx + dy*dy + dz*dz)
      end do
    end do

    if (smooth) then
      call espf_smooth_s(gxyz(:nptcur,:), nptcur, at, nat, sw_scale*rvdw, sw_delta, sk)
    else
      sk = 1.0_dp
    end if

    allocate(minv(nat,nat), bvec(nat), gvec(nat), gbig(nptcur), bbig(nptcur))
    ! b = M^-1 a ;  g = Z u
    if (smooth) then
      ! M = T diag(s) T^T ; solve M b = a
      do j = 1, nat
        do i = 1, nat
          minv(i,j) = sum(sk(1:nptcur)*tmat(i,1:nptcur)*tmat(j,1:nptcur))
        end do
      end do
      bvec = amm
      allocate(ipiv(nat))
      call dgesv(nat, 1, minv, nat, ipiv, bvec, nat, info)
      if (info /= 0) call show_message('espf_grad_weight: dgesv failed', WITH_ABORT)
      deallocate(ipiv)
    else
      minv = matmul(zmat(1:nat,1:nptcur), transpose(zmat(1:nat,1:nptcur)))  ! ZZ^T=M^-1
      bvec = matmul(minv, amm)
    end if
    gvec = matmul(zmat(1:nat,1:nptcur), u)
    ! G_k = (T^T g)_k ,  B_k = (T^T b)_k
    gbig = matmul(transpose(tmat), gvec)
    bbig = matmul(transpose(tmat), bvec)

    ! dT term: dE/dR_{i,c} += scal * sum_k s_k dT(i,k,c) [ b_i(u_k-G_k) - g_i B_k ]
    do i = 1, nat
      do k = 1, nptcur
        dx = gx(k) - at(1,i); dy = gy(k) - at(2,i); dz = gz(k) - at(3,i)
        r2 = dx*dx + dy*dy + dz*dz
        r3 = r2 * sqrt(r2)
        fac = scal * sk(k) * ( bvec(i)*(u(k) - gbig(k)) - gvec(i)*bbig(k) ) / r3
        grad(1,i) = grad(1,i) + fac*dx
        grad(2,i) = grad(2,i) + fac*dy
        grad(3,i) = grad(3,i) + fac*dz
      end do
    end do

    ! ds term (smooth only): dE/dR_{J,c} += scal * sum_k ds_k/dR_{J,c} B_k(u_k-G_k)
    if (smooth) then
      do k = 1, nptcur
        wk = bbig(k) * (u(k) - gbig(k))            ! B_k (u_k - G_k)
        if (wk == 0.0_dp) cycle
        do j = 1, nat
          dx = gx(k) - at(1,j); dy = gy(k) - at(2,j); dz = gz(k) - at(3,j)
          d = sqrt(dx*dx + dy*dy + dz*dz)
          xx = (d - sw_scale*rvdw(j)) / sw_delta
          sp = espf_dsstep(xx)
          if (sp == 0.0_dp) cycle
          sig = espf_sstep(xx)                     ! in (0,1) where sp/=0
          pe = sk(k) / sig                         ! prod_{b/=j} S_b
          ! ds_k/dR_{J,c} = pe * sp/delta * (R_J - g_k)_c/d  = pe*sp/delta*(-dx)/d
          coef = scal * wk * pe * sp / (sw_delta * d)
          grad(1,j) = grad(1,j) - coef*dx
          grad(2,j) = grad(2,j) - coef*dy
          grad(3,j) = grad(3,j) - coef*dz
        end do
      end do
    end if

    deallocate(gx, gy, gz, u, tmat, sk, minv, bvec, gvec, gbig, bbig)
  end subroutine espf_grad_weight

!--------------------------------------------------------------------------------

!> @brief C1 smoothstep S(x): 0 for x<=0, 1 for x>=1, x^2(3-2x) in between.
  elemental function espf_sstep(x) result(s)
    use precision, only: dp
    real(kind=dp), intent(in) :: x
    real(kind=dp) :: s
    if (x <= 0.0_dp) then
      s = 0.0_dp
    else if (x >= 1.0_dp) then
      s = 1.0_dp
    else
      s = x*x*(3.0_dp - 2.0_dp*x)
    end if
  end function espf_sstep

!> @brief S'(x) for the C1 smoothstep.
  elemental function espf_dsstep(x) result(d)
    use precision, only: dp
    real(kind=dp), intent(in) :: x
    real(kind=dp) :: d
    if (x <= 0.0_dp .or. x >= 1.0_dp) then
      d = 0.0_dp
    else
      d = 6.0_dp*x*(1.0_dp - x)
    end if
  end function espf_dsstep

!> @brief Smooth per-point grid weights s_k = prod_b S((|g_k-R_b|-Rvdw_b)/delta).
!>        Zero inside any atom's VDW sphere, 1 outside, smooth in the shell of
!>        width delta.  Makes the ESPF grid contribution a smooth function of
!>        nuclear geometry (replaces the hard VDW pruning in add_atom_grid).
  subroutine espf_smooth_s(xyz, np, at, nat, rvdw, delta, s)
    use precision, only: dp
    real(kind=dp), intent(in)  :: xyz(:,:)     ! (npt,3)
    integer, intent(in)        :: np, nat
    real(kind=dp), intent(in)  :: at(:,:)      ! (3,nat)
    real(kind=dp), intent(in)  :: rvdw(:)      ! (nat)
    real(kind=dp), intent(in)  :: delta
    real(kind=dp), intent(out) :: s(:)         ! (np)
    integer :: k, b
    real(kind=dp) :: d
    do k = 1, np
      s(k) = 1.0_dp
      do b = 1, nat
        d = norm2(xyz(k,:) - at(:,b))
        s(k) = s(k) * espf_sstep((d - rvdw(b))/delta)
        if (s(k) == 0.0_dp) exit
      end do
    end do
  end subroutine espf_smooth_s

!> @brief Weighted ESPF pseudoinverse Z = M^-1 T diag(s), M = T diag(s) T^T,
!>        T_{i,k}=1/|R_i-g_k|.  Reduces to the unweighted fit when s=1.
  subroutine espf_weights_w(x, y, z, ttt, at, s)
    use precision, only: dp
    use messages, only: show_message, with_abort
    implicit none
    real(kind=dp), intent(in)    :: x(:), y(:), z(:), at(:,:), s(:)
    real(kind=dp), intent(inout) :: ttt(:,:)        ! (nat, npts) <- Z on output
    integer :: i, j, k, nat, npts, info
    real(kind=dp), allocatable :: tmat(:,:), mmat(:,:)
    integer, allocatable :: ipiv(:)
    nat = ubound(at, 2)
    npts = ubound(x, 1)
    allocate(tmat(nat,npts), mmat(nat,nat), ipiv(nat))
    do k = 1, npts
      do i = 1, nat
        tmat(i,k) = 1.0_dp / norm2(at(:,i) - [x(k), y(k), z(k)])
      end do
    end do
    ! M = T diag(s) T^T
    do j = 1, nat
      do i = 1, nat
        mmat(i,j) = sum(s(1:npts)*tmat(i,1:npts)*tmat(j,1:npts))
      end do
    end do
    ! RHS = T diag(s) ; solve M Z = RHS  -> Z = M^-1 T diag(s)
    do k = 1, npts
      ttt(1:nat,k) = tmat(1:nat,k)*s(k)
    end do
    call dgesv(nat, npts, mmat, nat, ipiv, ttt, size(ttt,1), info)
    if (info /= 0) call show_message('espf_weights_w: dgesv failed', WITH_ABORT)
    deallocate(tmat, mmat, ipiv)
  end subroutine espf_weights_w

!--------------------------------------------------------------------------------

!> @brief Computes the integral weights for computing the ESPF integrals
!> @param[in]      x,y,z        Grid Cartesian coordinates
!> @param[in,out]  ttt          Integral weights
!> @param[in]      at           Atom Cartesian coordinates
!
!> @detail Form the electrostatic kernel T(nat,npts) and forms the pseudoinverse
!>         (T^+T)^-1*T^+ (in which + = dagger)
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Jul, 2024_ Initial release
  subroutine espf_weights(x, y, z, ttt, at)
    use io_constants, only: iw
    use precision, only: dp
    implicit none

    real(kind=dp), intent(in) :: x(:), y(:), z(:), at(:,:)
    real(kind=dp), intent(inout) :: ttt(:,:)

    integer :: i, j, nat, npts, lwork, info
    real(kind=dp), allocatable :: b(:,:), r(:,:)
    real(kind=dp), allocatable :: work(:)
    real(kind=dp) :: workk(1)

    nat = ubound(at, 2)
    npts = ubound(x, 1)

    allocate(b(npts,npts), r(nat,npts))

    ! Compute the matrix of the inverse distances
    b=0.0d0
    do j = 1, npts
      b(j,j)=1.0d0
      do i = 1, nat
        r(i,j) = 1/norm2(at(:,i) - [x(j), y(j), z(j)])
      end do
    end do

    ! Workspace query for optimal lwork size
    lwork = -1
    call dgels('T', nat, npts, npts, r, nat, b, npts, workk, lwork, info)
    lwork = int(workk(1))
    allocate(work(lwork))

    ! Solve the least-squares problem T * X ≈ B
    call dgels('T', nat, npts, npts, r, nat, b, npts, work, lwork, info)

    ttt(:nat,:npts)=b(:nat,:npts)

  end subroutine espf_weights

!--------------------------------------------------------------------------------

!> @brief Adds atomic-centered spherical grid to form the molecular grid for ESP calculations
!> @param[in]      nat          Number of QM centers
!> @param[in]      npt          Maximum number of gridpoints
!> @param[in]      nlayers      Number of layers for the grid
!> @param[in]      layers       Radius for the different layers
!> @param[in]      npt_layer    Number of points per layer
!> @param[in]      typ_layer
!> @param[in]      zn           atomic charges
!> @param[in]      atoms_xyz    atomic coordinates
!> @param[in,out]  xyz          grid coordinates
!> @param[in,out]  ttt          grid weights
!> @param[in,out]  nptcur       final number of grid points
!
!> @detail Form the numerical grid for ESPF computations and construct the integral weights
!>         This routine is adapted from oqp/modules/resp.F90
!
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Oct, 2024_ Initial release
  subroutine form_espf_grid(nat,npt,nlayers,layers,npt_layer,typ_layer,zn,atoms_xyz,xyz,ttt,nptcur)
    use precision, only: dp
    use lebedev, only: lebedev_get_grid
    use elements, only: ELEMENTS_VDW_RADII
    use physical_constants, only: ANGSTROM_TO_BOHR
    use messages, only: show_message, with_abort
    implicit none

    integer, intent(in) :: nat,npt,nlayers
    real(kind=dp), intent(in) :: zn(nat),atoms_xyz(3,nat),layers(nlayers)
    integer, intent(in) :: npt_layer(nlayers),typ_layer(nlayers)
    real(kind=dp), intent(inout) :: xyz(npt,3),ttt(nat,npt)
    integer, intent(inout) :: nptcur

    real(kind=dp), allocatable :: leb(:,:), lebw(:)
    real(kind=dp), allocatable :: vdwrad(:), excl_vdw(:), wt(:)
    integer, allocatable :: neigh(:)

    integer :: nadd, nleb, ok
    integer :: i, j, k, layer, iz
    logical :: keepall, smooth, gamess_mode
    character(len=16) :: env_k
    integer :: st_k
    real(kind=dp) :: sw_delta, sw_scale, rmax_i, vdwenv_i, prec1_bohr
    real(kind=dp), allocatable :: gamess_vdw(:)

    !> GAMESS ESPF VDW radii (Å, from LEBGRD in espf.src, Emsley 1991 + Bondi 1964).
    !> Zero entries fall back to OpenQP's ELEMENTS_VDW_RADII at runtime.
    real(kind=dp), parameter :: gamess_rvdw_ang(104) = [ &
      1.20_dp, 1.22_dp,                                                & ! H, He
      0.00_dp, 0.00_dp,                                                & ! Li, Be
      2.08_dp, 1.85_dp, 1.54_dp, 1.50_dp, 1.35_dp, 1.60_dp,          & ! B-Ne
      2.31_dp, 0.00_dp, 2.05_dp, 2.00_dp, 1.90_dp, 1.85_dp,          & ! Na-S
      1.81_dp, 1.91_dp,                                                & ! Cl, Ar
      2.31_dp,                                                         & ! K
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Ca-Fe
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,          & ! Co-Ge
      2.00_dp, 2.00_dp, 1.95_dp, 1.98_dp,                             & ! As, Se, Br, Kr
      2.44_dp,                                                         & ! Rb
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Sr-Ru
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,          & ! Rh-Sn
      2.20_dp, 2.20_dp, 2.15_dp, 0.00_dp,                             & ! Sb, Te, I, Xe
      2.62_dp,                                                         & ! Cs
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Ba-Nd
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Pm-Yb
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Lu-Re
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Os-Pb
      2.40_dp, 0.00_dp, 0.00_dp, 0.00_dp,                             & ! Bi, Po, At, Rn
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Fr-Am
      0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, & ! Cm-Rf
      0.00_dp, 0.00_dp, 0.00_dp                                        & ! Db, Sg, Bh (pad to 104)
    ]

    ! npt may be larger than npt_layer total when GAMESS mode adds 146/shell;
    ! use max(npt, nat*nlayers*146) in the caller, or just pad npt there.
    allocate(leb(max(maxval(npt_layer),146),3), &
             lebw(max(maxval(npt_layer),146)), &
             vdwrad(nat), &
             excl_vdw(nat), &
             gamess_vdw(nat), &
             neigh(nat), &
             wt(npt), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

!   Set up the grid
    nptcur = 0  ! Current number of point in a grid

!   ESPF_KEEPALL=1 keeps every Lebedev point (fixed grid membership) instead of
!   removing points that fall inside neighbour VDW spheres.  The hard removal is
!   a step function of geometry: as atoms move, points blink in/out, making the
!   ESPF energy non-smooth and leaking energy in dynamics.  Keeping all points
!   makes the atom-centred grid translate smoothly with the nuclei.
    keepall = .false.
    call get_environment_variable('ESPF_KEEPALL', env_k, status=st_k)
    if (st_k == 0) then
      if (trim(env_k) == '1' .or. trim(env_k) == 'on') keepall = .true.
    end if
    ! Smooth-switching weighted ESPF grid (DEFAULT): keep all points + smooth
    ! weights for energy-conserving dynamics.  ESPF_SMOOTH=0 restores the hard
    ! VDW-pruned grid (matches the original/GAMESS grid).
    smooth = .true.
    call get_environment_variable('ESPF_SMOOTH', env_k, status=st_k)
    if (st_k == 0) then
      if (trim(env_k) == '0' .or. trim(env_k) == 'off') smooth = .false.
    end if
    sw_delta = 0.7_dp
    call get_environment_variable('ESPF_SWDELTA', env_k, status=st_k)
    if (st_k == 0) then
      if (len_trim(env_k) > 0) read(env_k,*) sw_delta
    end if
    sw_scale = 1.8_dp
    call get_environment_variable('ESPF_SWSCALE', env_k, status=st_k)
    if (st_k == 0) then
      if (len_trim(env_k) > 0) read(env_k,*) sw_scale
    end if
    if (smooth) keepall = .true.

!   GAMESS-identical mode: ESPF_GAMESS=1 reproduces GAMESS LEBGRD exactly
!   (RVDW table from Emsley 1991/Bondi 1964, NANG=146 Lebedev per shell,
!    linear shell placement RINC=IS*(4r-r-0.1)/NRAD, fixed r+0.1 Å exclusion).
    gamess_mode = .false.
    call get_environment_variable('ESPF_GAMESS', env_k, status=st_k)
    if (st_k == 0) then
      if (trim(env_k) == '1' .or. trim(env_k) == 'on') gamess_mode = .true.
    end if
    prec1_bohr = 0.1_dp * ANGSTROM_TO_BOHR  ! GAMESS PREC1 = 0.1 Å
    ! Build GAMESS VDW radius array in bohr; fall back to OpenQP for unknown elements
    do i = 1, nat
      iz = int(zn(i))
      if (iz >= 1 .and. iz <= size(gamess_rvdw_ang) .and. gamess_rvdw_ang(iz) > 0.0_dp) then
        gamess_vdw(i) = gamess_rvdw_ang(iz) * ANGSTROM_TO_BOHR
      else
        gamess_vdw(i) = ELEMENTS_VDW_RADII(iz)
      end if
    end do
    ! GAMESS mode forces hard pruning (no smooth, no keepall)
    if (gamess_mode) then
      smooth  = .false.
      keepall = .false.
    end if

!   Fixed exclusion radii: scaled by the MINIMUM layer factor (layers(1) = 1.4)
!   so the exclusion sphere is independent of the current layer being built.
!
!   Why not base VDW (layer 1.0)?
!     Too permissive: OpenQP shells sit at 1.4–2.0 × r_vdw, so base-VDW
!     exclusion admits ~1400 extra near-atom intermediate-zone points that
!     make the T matrix ill-conditioned and energy noisier.
!
!   Why not layer-scaled exclusion (original)?
!     The outer shells (1.8–2.0×) then sit right at their own exclusion
!     boundary.  A C–H stretch of 0.1 Å flips the outer-shell point toward H
!     from "kept" to "excluded" (margin ≈ 0.10 Å), creating discrete grid
!     membership changes → non-smooth energy surface → drift.
!
!   Why min-layer (1.4 × r_vdw)?
!     • Same inner-shell points are excluded as before (inner shell exclusion
!       is identical to the original 1.4-layer-scaled rule).
!     • Outer shells (1.8–2.0×) now have a margin of ~0.5–0.8 Å before any
!       point flips out; normal MD vibrational amplitudes (< 0.3 Å) are safe.
!     • Intermediate-zone points (1.4–2.0 × r_vdw from neighbours) are admitted
!       on outer shells but at physically reasonable positions (not near singularity).
!   This is the correct fixed-scale complement to OpenQP's 1.4–2.0 shell scheme.
    if (gamess_mode) then
      excl_vdw = gamess_vdw + prec1_bohr        ! GAMESS: r_vdw + 0.1 Å (fixed, Emsley radii)
    else
      excl_vdw = ELEMENTS_VDW_RADII(int(zn)) * layers(1)  ! = 1.4 × r_vdw (innermost layer)
    end if

!   Loop over layers and add spherical grid points on each atom to the molecular grid
    do layer = 1, nlayers

!     Get the atomic radii on which to place new grid layer
      if (gamess_mode) then
        ! GAMESS shell radius: IS*(RMAX-VDWEnv)/NRAD where RMAX=4*r_vdw, VDWEnv=r_vdw+0.1 Å
        ! => layer*(nlayers*gamess_vdw - gamess_vdw - prec1)/nlayers
        do i = 1, nat
          vdwrad(i) = real(layer,dp) * (real(nlayers,dp)*gamess_vdw(i) - gamess_vdw(i) - prec1_bohr) &
                    / real(nlayers,dp)
        end do
      else
        vdwrad = ELEMENTS_VDW_RADII(int(zn))*layers(layer)
      end if

!     Get grid
      if (gamess_mode) then
        nleb = 146
        leb = 0; lebw = 0
        call lebedev_get_grid(nleb, leb, lebw, 0)   ! type 0 = standard Lebedev, NANG=146
      else
        nleb = npt_layer(layer)
        leb = 0
        lebw = 0
        call lebedev_get_grid(nleb, leb, lebw, typ_layer(layer))
      end if

!     Add new grid layer for each atom, remove inner points
      do i = 1, nat
        ! GAMESS: atom IA is included in its own exclusion loop, so IS=1 shell
        ! (radius < r_vdw+0.1) is always fully self-excluded.  add_atom_grid skips
        ! self (cur_atom), so mimic GAMESS by skipping this layer when the shell
        ! radius is within the atom's own exclusion sphere.
        if (gamess_mode .and. vdwrad(i) <= excl_vdw(i)) then
          nadd = 0
          nptcur = nptcur + nadd
          cycle
        end if
        if (keepall) then
          do j = 1, nleb
            xyz(nptcur+j,1) = atoms_xyz(1,i) + vdwrad(i)*leb(j,1)
            xyz(nptcur+j,2) = atoms_xyz(2,i) + vdwrad(i)*leb(j,2)
            xyz(nptcur+j,3) = atoms_xyz(3,i) + vdwrad(i)*leb(j,3)
            wt(nptcur+j)    = lebw(j)
          end do
          nadd = nleb
        else
          call add_atom_grid( &
            x=xyz(nptcur+1:,1), &
            y=xyz(nptcur+1:,2), &
            z=xyz(nptcur+1:,3), &
            wts=wt(nptcur+1:), &
            nadd=nadd, &
            atpts=leb(:nleb,:), &
            atwts=lebw(:nleb), &
            atoms_xyz=atoms_xyz, &
            atoms_rad=vdwrad, &
            cur_atom=i, &
            neighbours=neigh, &
            excl_rad=excl_vdw)
        end if
        nptcur = nptcur + nadd
      end do
    end do

!   Compute the weights
    if (smooth) then
      block
        real(kind=dp) :: rvdw(nat), s(nptcur)
        integer :: kk, nkeep
        ! Core radius = sw_scale*Rvdw: larger scale zeroes (and drops) more inner
        ! points, keeping the grid ~ as sparse as the hard prune for efficiency.
        rvdw = sw_scale * ELEMENTS_VDW_RADII(int(zn))
        call espf_smooth_s(xyz(:nptcur,:), nptcur, atoms_xyz, nat, rvdw, sw_delta, s)
        ! Compact out points with s_k == 0 (deep inside a neighbour VDW core):
        ! they contribute exactly 0 to both energy and gradient (S=0 and S'=0),
        ! so dropping them is discontinuity-free and recovers the keep-all cost.
        nkeep = 0
        do kk = 1, nptcur
          if (s(kk) > 0.0_dp) then
            nkeep = nkeep + 1
            xyz(nkeep,:) = xyz(kk,:)
            s(nkeep) = s(kk)
          end if
        end do
        nptcur = nkeep
        call espf_weights_w(xyz(:nptcur,1), xyz(:nptcur,2), xyz(:nptcur,3), ttt, atoms_xyz, s(:nptcur))
      end block
    else
      call espf_weights( &
           x=xyz(:nptcur,1),  &
           y=xyz(:nptcur,2),  &
           z=xyz(:nptcur,3),  &
           ttt=ttt,           &
           at = atoms_xyz)
    end if

    ! one-shot grid-size report (ESPF_NPRINT=1)
    if (st_k >= 0) then
      call get_environment_variable('ESPF_NPRINT', env_k, status=st_k)
      if (st_k == 0 .and. trim(env_k) == '1') then
        block
          integer :: ud
          open(newunit=ud, file='/tmp/espf_npts.txt', position='append', action='write')
          if (gamess_mode) then
            write(ud,'(a,i7)') 'GAMESS nptcur=', nptcur
          else if (smooth) then
            write(ud,'(a,i7,a,f5.2,a,f5.2)') 'smooth nptcur=', nptcur, &
                  '  scale=', sw_scale, '  delta=', sw_delta
          else
            write(ud,'(a,i7)') 'pruned nptcur=', nptcur
          end if
          close(ud)
        end block
      end if
    end if

  end subroutine form_espf_grid

!--------------------------------------------------------------------------------

!> @brief Print partial ESPF charges
!>
!> @param[in]      infos        OQP handle
!> @param[in]      chg          atomic partial charges
!
!> @author   Miquel Huix-Rotllant
!
!> @detail This is a modified copy of print_charges of oqp/modules/resp.F90
!
!     REVISION HISTORY:
!> @date _Oct, 2024_ Initial release
  subroutine print_charges(infos,partial_charges,iw)
    use precision, only: dp
    use oqp_tagarray_driver
    use elements, only: ELEMENTS_SHORT_NAME
    use types, only: information
    use messages, only: show_message, with_abort
    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), contiguous :: partial_charges(:)
    integer, intent(in) :: iw

    integer :: i, elem, nat

    nat = ubound(infos%atoms%zn, 1)

    write(iw,'(2/)')
    write(iw,'(4x,a)') '=============================='
    write(iw,'(4x,a)') 'ESPF QM/MM charges calculation'
    write(iw,'(4x,a)') '=============================='

    write(iw,'(/,30("^"))')
    write(iw,'(/a8,a8,a14)') '#', 'Name', 'Charge'
    write(iw,'(30("-"))')

    do i = 1, nat
      elem = nint(infos%atoms%zn(i))
      write(iw,'(i8,a8,f14.6)') i, ELEMENTS_SHORT_NAME(elem), partial_charges(i)
    end do

    write(iw,'(30("="))')

  end subroutine print_charges

end module qmmm_mod
