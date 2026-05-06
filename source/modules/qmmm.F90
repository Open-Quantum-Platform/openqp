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
      call tagarray_get_data(infos%dat, OQP_TD_ABXC, td_abxc)
      tmp = tmp + td_abxc
      write(iw,'(4x,a,i4)') 'Using UNRELAXED excited-state density for ESPF charges, state ', istate
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

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

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
    call espf_grad(&
         x=xyz(:nptcur,1),&
         y=xyz(:nptcur,2),&
         z=xyz(:nptcur,3),&
         at=infos%atoms%xyz,&
         wt=wt,&
         zn=infos%atoms%zn,&
         pchg=partial_charges,&
         grad=espf_grad_ta)
!    espf_grad = grad

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

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms
    use_relaxed = .true.
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
      call tagarray_get_data(infos%dat, OQP_TD_P, td_p)
      dens = dens + td_p(:,1) + td_p(:,2)
    else
      call tagarray_get_data(infos%dat, OQP_TD_ABXC, td_abxc)
      dens = dens + td_abxc
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
    call espf_grad(&
         x=xyz(:nptcur,1),&
         y=xyz(:nptcur,2),&
         z=xyz(:nptcur,3),&
         at=infos%atoms%xyz,&
         wt=wt,&
         zn=infos%atoms%zn,&
         pchg=partial_charges,&
         grad=espf_grad_ta)
!    espf_grad = grad

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
    use messages, only: show_message, with_abort
    implicit none

    integer, intent(in) :: nat,npt,nlayers
    real(kind=dp), intent(in) :: zn(nat),atoms_xyz(3,nat),layers(nlayers)
    integer, intent(in) :: npt_layer(nlayers),typ_layer(nlayers)
    real(kind=dp), intent(inout) :: xyz(npt,3),ttt(nat,npt)
    integer, intent(inout) :: nptcur

    real(kind=dp), allocatable :: leb(:,:), lebw(:)
    real(kind=dp), allocatable :: vdwrad(:), wt(:)
    integer, allocatable :: neigh(:)

    integer :: nadd, nleb, ok
    integer :: i, j, k, layer

    allocate(leb(maxval(npt_layer),3), &
             lebw(maxval(npt_layer)), &
             vdwrad(nat), &
             neigh(nat), &
             wt(npt), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

!   Set up the grid
    nptcur = 0  ! Current number of point in a grid

!   Loop over layers and add spherical grid points on each atom to the molecular grid
    do layer = 1, nlayers

!     Get the atomic radii on which to place new grid layer
      vdwrad = ELEMENTS_VDW_RADII(int(zn))*layers(layer)

!     Get grid
      nleb = npt_layer(layer)
      leb = 0
      lebw = 0
      call lebedev_get_grid(nleb, leb, lebw, typ_layer(layer))

!     Add new grid layer for each atom, remove inner points
      do i = 1, nat
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
          neighbours=neigh)
        nptcur = nptcur + nadd
      end do
    end do

!   Compute the weights
    call espf_weights( &
         x=xyz(:nptcur,1),  &
         y=xyz(:nptcur,2),  &
         z=xyz(:nptcur,3),  &
         ttt=ttt,           &
         at = atoms_xyz)

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
