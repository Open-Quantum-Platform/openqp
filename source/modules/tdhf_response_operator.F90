module tdhf_response_operator_mod

  use precision, only: dp

  implicit none

  private
  public :: build_tdhf_response_matrices
  public :: apply_tdhf_ao_operators

contains

!###############################################################################

  subroutine build_tdhf_response_matrices(infos, amb_mo, apb_mo)
    ! Construct dense occupied-virtual representations of A-B and A+B by
    ! applying the same ERI and XC operations used by the TD energy solver to
    ! the complete set of unit vectors.  The Hessian response equations use
    ! these matrices for their projected solves and residual checks.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_E_MO_A, OQP_VEC_MO_A
    use tdhf_lib, only: iatogen, mntoia, esum
    use mathlib, only: orthogonal_transform
    use messages, only: show_message, WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: amb_mo(:,:), apb_mo(:,:)

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: mo(:,:), eps(:)
    real(kind=dp), allocatable, target :: dx(:,:,:)
    real(kind=dp), allocatable :: amb_ao(:,:,:), apb_ao(:,:,:)
    real(kind=dp), allocatable :: unitvec(:,:), mo_density(:,:), scr(:,:)
    integer :: ivec, lexc, nbf, nocc, nvir

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nocc = infos%mol_prop%nocc
    nvir = nbf - nocc
    lexc = nocc*nvir

    if (infos%control%scftype /= 1 .or. infos%mol_prop%mult /= 1 .or. &
        infos%tddft%mult /= 1 .or. infos%tddft%tda) then
      call show_message('TD response matrices require a closed-shell restricted '// &
                        'singlet full-response calculation.', WITH_ABORT)
    end if
    if (any(shape(amb_mo) /= [lexc,lexc]) .or. &
        any(shape(apb_mo) /= [lexc,lexc])) then
      call show_message('TD response-matrix output has the wrong shape.', WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, OQP_E_MO_A, eps)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    allocate(dx(nbf,nbf,lexc), unitvec(lexc,lexc), &
             mo_density(nbf,nbf), scr(nbf,nbf), source=0.0_dp)
    do ivec = 1, lexc
      unitvec(ivec,ivec) = 1.0_dp
      call iatogen(unitvec(:,ivec), mo_density, nocc, nocc)
      call orthogonal_transform('t', nbf, mo, mo_density, dx(:,:,ivec), scr)
    end do

    allocate(amb_ao(nbf,nbf,lexc), apb_ao(nbf,nbf,lexc))
    call apply_tdhf_ao_operators(infos, dx, amb_ao, apb_ao)

    do ivec = 1, lexc
      call mntoia(apb_ao(:,:,ivec), apb_mo(:,ivec), mo, mo, nocc, nocc)
      call esum(eps, apb_mo, unitvec, nocc, ivec)
      call mntoia(amb_ao(:,:,ivec), amb_mo(:,ivec), mo, mo, nocc, nocc)
      call esum(eps, amb_mo, unitvec, nocc, ivec)
    end do

    deallocate(dx, amb_ao, apb_ao, unitvec, mo_density, scr)
  end subroutine build_tdhf_response_matrices

!###############################################################################

  subroutine apply_tdhf_ao_operators(infos, density, amb_ao, apb_ao)
    ! Apply the two-electron and XC parts of A-B and A+B to general AO
    ! transition densities.  Orbital-energy differences are not added here.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use mathlib, only: symmetrize_matrix
    use mod_dft, only: dft_initialize, dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_fxc, only: tddft_fxc
    use messages, only: show_message, WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: density(:,:,:)
    real(kind=dp), intent(out) :: amb_ao(:,:,:), apb_ao(:,:,:)

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(int2_td_data_t), target :: int2_data
    type(dft_grid_t) :: molgrid
    real(kind=dp), contiguous, pointer :: mo(:,:)
    real(kind=dp), allocatable, target :: dx(:,:,:)
    real(kind=dp), pointer :: amb_work(:,:,:), apb_work(:,:,:)
    real(kind=dp) :: scale_exch
    integer :: ivec, nbf, nvec
    logical :: dft

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nvec = size(density,3)
    dft = infos%control%hamilton == 20
    if (size(density,1) /= nbf .or. size(density,2) /= nbf .or. &
        any(shape(amb_ao) /= [nbf,nbf,nvec]) .or. &
        any(shape(apb_ao) /= [nbf,nbf,nvec])) then
      call show_message('AO TD response-operator arrays have incompatible shapes.', WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    allocate(dx, source=density)
    if (dft) call dft_initialize(infos, basis, molgrid)
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%dft%hfscale
    int2_data = int2_td_data_t(d2=dx, int_apb=.true., int_amb=.true., &
      tamm_dancoff=.false., tamm_dancoff_coulomb=.true., &
      scale_exchange=scale_exch)
    call int2_driver%run(int2_data, &
      cam=dft.and.infos%dft%cam_flag, alpha=infos%dft%cam_alpha, &
      beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
    apb_work => int2_data%apb(:,:,:,1)
    amb_work => int2_data%amb(:,:,:,1)
    if (dft) then
      do ivec = 1, nvec
        call symmetrize_matrix(dx(:,:,ivec), nbf)
      end do
      call tddft_fxc(basis=basis, molgrid=molgrid, isvecs=.true., wf=mo, &
        fx=apb_work, dx=dx, nmtx=nvec, threshold=0.0_dp, infos=infos)
    end if
    amb_ao = amb_work
    apb_ao = apb_work
    call int2_driver%clean()
    if (dft) call dftclean(infos)
    deallocate(dx)
  end subroutine apply_tdhf_ao_operators

end module tdhf_response_operator_mod
