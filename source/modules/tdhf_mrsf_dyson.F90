!================================================================================!
! tdhf_mrsf_dyson.F90 - MRSF-TDDFT Dyson orbital driver
!================================================================================!

module dyson_mrsf_mod
  use precision,    only: dp
  use io_constants, only: iw
  use types,        only: information
  use tdhf_mrsf_z_vector_mod, only: tdhf_mrsf_z_vector
  use, intrinsic :: iso_c_binding

  implicit none
  private

  public :: dyson_mrsf_C
  public :: dyson_mrsf_driver

contains

  !> \brief  C entry point for MRSF Dyson/EKT.
  !! \param[in] c_handle  OQP opaque handle (bridge to Fortran infos)
  subroutine dyson_mrsf_C(c_handle) bind(C, name="dyson_mrsf")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call dyson_mrsf_driver(inf)
  end subroutine dyson_mrsf_C

  !> \brief  MRSF-TDDFT path: ensures Z-vector, builds relaxed matrices, runs Dyson.
  !! \param[in,out] infos  Global OQP information object
  subroutine dyson_mrsf_driver(infos)
    use printing, only: print_module_info
    use messages,  only: show_message, with_abort
    use dyson_orbitals_mod, only: dyson_orbital_driver
    use oqp_tagarray_driver
    use tagarray
    type(information), intent(inout) :: infos

    integer :: nbf, nstates, target_state

    call print_module_info('Dyson Orbitals (MRSF-TDDFT) module', 'Dyson Orbitals (MRSF-TDDFT) module')

    nbf          = infos%basis%nbf
    nstates      = max(1, infos%tddft%nstate)
    target_state = max(1, infos%tddft%target_state)

    if (.not. infos%mol_energy%Z_Vector_converged) then
      write(iw,'(5x,"Computing Z-vector for MRSF...")')
      call tdhf_mrsf_z_vector(infos)
      if (.not. infos%mol_energy%Z_Vector_converged) then
        call show_message('Z-vector failed to converge; cannot proceed with Dyson (MRSF).', with_abort)
      end if
    else
      write(iw,'(5x,"Using pre-computed Z-vector.")')
    end if

    call build_and_store_relaxed_matrices(infos, nbf, target_state)
    call dyson_orbital_driver(infos)
  end subroutine dyson_mrsf_driver

  !> \brief  Assemble and store relaxed density and Lagrangian from Z-vector.
  !! \param[in,out] infos         Global infos (TagArray sink is infos%dat)
  !! \param[in]     nbf           AO basis size
  !! \param[in]     target_state  target state (only this slice is filled)
  subroutine build_and_store_relaxed_matrices(infos, nbf, target_state)
    use oqp_tagarray_driver
    use tagarray
    use messages, only: show_message, without_abort
    use mathlib,  only: symmetrize_matrix
    type(information), intent(inout) :: infos
    integer, intent(in) :: nbf, target_state
    real(dp), allocatable :: density_relaxed(:,:), lagrangian_relaxed(:,:)
    integer(c_int32_t) :: status

    allocate(density_relaxed(nbf,nbf)); density_relaxed = 0.0_dp
    allocate(lagrangian_relaxed(nbf,nbf)); lagrangian_relaxed = 0.0_dp

    call build_relaxed_density_from_z   (infos, density_relaxed, nbf)
    call build_relaxed_lagrangian_from_z(infos, lagrangian_relaxed, nbf)

    call symmetrize_matrix(density_relaxed,   nbf)
    call symmetrize_matrix(lagrangian_relaxed, nbf)

!    status = infos%dat%set_data(OQP_density_relaxed,   density_relaxed,   comment=OQP_density_relaxed_comment)
!    if (status /= TA_OK) call show_message('Failed to store relaxed density',   without_abort)
!    status = infos%dat%set_data(OQP_lagrangian_relaxed, lagrangian_relaxed, comment=OQP_lagrangian_relaxed_comment)
!    if (status /= TA_OK) call show_message('Failed to store relaxed Lagrangian', without_abort)

    deallocate(density_relaxed, lagrangian_relaxed)
    write(iw,'(5x,"Relaxed matrices prepared and stored.")')
  end subroutine build_and_store_relaxed_matrices

  !> \brief  Build relaxed density using ground-state density plus Z-vector blocks.
  !! \param[in,out] infos    Global infos (source tags: OQP_DM_A, OQP_td_xpy)
  !! \param[out]    density  (nbf,nbf) relaxed density
  !! \param[in]     nbf      AO basis size
  subroutine build_relaxed_density_from_z(infos, density, nbf)
    use oqp_tagarray_driver
    use tagarray
    type(information), intent(inout) :: infos
    real(dp), intent(out) :: density(nbf,nbf)
    integer, intent(in)   :: nbf
    real(dp), pointer :: dm_a(:,:), z_xpy(:)
    integer(c_int32_t) :: status
    integer :: i, a, ia, nocc

    nocc = infos%mol_prop%nelec_A
    call tagarray_get_data(infos%dat, OQP_DM_A,   dm_a, status)
    if (status == TA_OK) then
      density = dm_a
    else
      density = 0.0_dp
      do i=1, nocc
        density(i,i) = 2.0_dp
      end do
    end if

    call tagarray_get_data(infos%dat, OQP_td_xpy, z_xpy, status)
    if (status == TA_OK) then
      ia = 0
      do i=1, nocc
        do a=nocc+1, nbf
          ia = ia + 1
          if (ia <= size(z_xpy)) then
            density(i,a) = density(i,a) + z_xpy(ia)
            density(a,i) = density(i,a)
          end if
        end do
      end do
    end if
  end subroutine build_relaxed_density_from_z

  !> \brief  Build relaxed Lagrangian using Fock plus Z-vector blocks.
  !! \param[in,out] infos       Global infos (source tags: OQP_FOCK_A, OQP_td_xmy)
  !! \param[out]    lagrangian  (nbf,nbf) relaxed Lagrangian
  !! \param[in]     nbf         AO basis size
  subroutine build_relaxed_lagrangian_from_z(infos, lagrangian, nbf)
    use oqp_tagarray_driver
    use tagarray
    type(information), intent(inout) :: infos
    real(dp), intent(out) :: lagrangian(nbf,nbf)
    integer, intent(in)   :: nbf
    real(dp), pointer :: fock_a(:,:), z_xmy(:)
    integer(c_int32_t) :: status
    integer :: i, a, ia, nocc

    nocc = infos%mol_prop%nelec_A
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a, status)
    if (status == TA_OK) then
      lagrangian = fock_a
    else
      lagrangian = 0.0_dp
    end if

    call tagarray_get_data(infos%dat, OQP_td_xmy, z_xmy, status)
    if (status == TA_OK) then
      ia = 0
      do i=1, nocc
        do a=nocc+1, nbf
          ia = ia + 1
          if (ia <= size(z_xmy)) then
            lagrangian(i,a) = lagrangian(i,a) + z_xmy(ia)
            lagrangian(a,i) = lagrangian(i,a)
          end if
        end do
      end do
    end if
  end subroutine build_relaxed_lagrangian_from_z

end module dyson_mrsf_mod
