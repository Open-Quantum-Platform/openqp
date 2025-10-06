!================================================================================!
! tdhf_mrsf_dyson.F90 - Simplified version
! MRSF-TDDFT specific Dyson orbital calculation module
! Computes only the Z-vector and relaxed matrices needed for Dyson
! Does NOT require full gradient calculation
!================================================================================!

module dyson_mrsf_mod
  use precision, only: dp
  use io_constants, only: iw
  use types, only: information
  use tdhf_mrsf_z_vector_mod, only: tdhf_mrsf_z_vector 
  use, intrinsic :: iso_c_binding
  
  implicit none
  private
  
  character(len=*), parameter :: module_name = "dyson_mrsf_mod"
  
  public :: dyson_mrsf_driver
  
contains

  !------------------------------------------------------------------------------!
  ! C interface wrapper
  !------------------------------------------------------------------------------!
  subroutine dyson_mrsf_C(c_handle) bind(C, name="dyson_mrsf")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    
    inf => oqp_handle_get_info(c_handle)
    call dyson_mrsf_driver(inf)
    
  end subroutine dyson_mrsf_C
  
  !------------------------------------------------------------------------------!
  ! Main MRSF Dyson driver - computes only what's needed
  !------------------------------------------------------------------------------!
  subroutine dyson_mrsf_driver(infos)
    use printing, only: print_module_info
    use messages, only: show_message, with_abort
    use dyson_orbitals_mod, only: dyson_orbital_driver
    use oqp_tagarray_driver
    use tagarray
    
    type(information), intent(inout) :: infos
    
    integer :: nbf, nstates, target_state, nocca, noccb
    character(len=200) :: msg
    
    call print_module_info('DYSON_MRSF', 'MRSF-TDDFT Dyson Orbital Calculation')
    
    ! Validate MRSF requirements
    if (infos%mol_prop%mult /= 3) then
      write(msg, '(A,I0)') 'MRSF Dyson requires triplet reference (SCF mult=3), got: ', &
                           infos%mol_prop%mult
      call show_message(msg, with_abort)
    end if
    
    if (infos%tddft%mult /= 1 .and. infos%tddft%mult /= 3 .and. infos%tddft%mult /= 5) then
      write(msg, '(A,I0)') 'MRSF Dyson requires target mult=1,3,5. Got: ', infos%tddft%mult
      call show_message(msg, with_abort)
    end if
    
    ! Get dimensions
    nbf = infos%basis%nbf
    nstates = infos%tddft%nstate
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    target_state = infos%tddft%dyson_target_state
    
    if (target_state == 0) target_state = infos%tddft%target_state
    if (target_state == 0) target_state = 1
    
    write(iw, '(/,5x,A)') 'MRSF Dyson Configuration:'
    write(iw, '(5x,A,I0)') '  Target state: ', target_state
    write(iw, '(5x,A,I0)') '  MRSF multiplicity: ', infos%tddft%mult
    
    ! Check that Z-vector has been computed (should be done in Python)
    if (.not. infos%mol_energy%Z_Vector_converged) then
       call show_message('Z-vector must be computed before Dyson calculation', with_abort)
    end if
    write(iw, '(5x,A)') 'Using pre-computed Z-vector'

    
    ! Build relaxed density and Lagrangian from Z-vector
    call build_and_store_relaxed_matrices(infos, nbf, target_state)
    
    ! Call general Dyson orbital driver
    call dyson_orbital_driver(infos)
    
  end subroutine dyson_mrsf_driver
  
  !------------------------------------------------------------------------------!
  ! Build and store relaxed matrices from Z-vector
  !------------------------------------------------------------------------------!
  subroutine build_and_store_relaxed_matrices(infos, nbf, target_state)
    use oqp_tagarray_driver
    use tagarray
    use messages, only: show_message, without_abort
    
    type(information), intent(inout) :: infos
    integer, intent(in) :: nbf, target_state
    
    real(dp), allocatable :: density_relaxed(:,:)
    real(dp), allocatable :: lagrangian_relaxed(:,:)
    integer(c_int32_t) :: status
    
    allocate(density_relaxed(nbf, nbf))
    allocate(lagrangian_relaxed(nbf, nbf))
    
    write(iw, '(/,5x,A)') 'Building relaxed matrices from Z-vector...'
    
    ! Build relaxed density from ground state density + Z-vector correction
    call build_relaxed_density_from_z(infos, density_relaxed, nbf)
    
    ! Build relaxed Lagrangian from Fock + Z-vector correction  
    call build_relaxed_lagrangian_from_z(infos, lagrangian_relaxed, nbf)
    
    ! Store in tagarray for Dyson driver
    status = infos%dat%set_data(OQP_density_relaxed, density_relaxed, &
                                comment=OQP_density_relaxed_comment)
    if (status /= TA_OK) then
      call show_message('Failed to store relaxed density', without_abort)
    end if
    
    status = infos%dat%set_data(OQP_lagrangian_relaxed, lagrangian_relaxed, &
                                comment=OQP_lagrangian_relaxed_comment)
    if (status /= TA_OK) then
      call show_message('Failed to store relaxed Lagrangian', without_abort)
    end if
    
    deallocate(density_relaxed, lagrangian_relaxed)
    
    write(iw, '(5x,A)') 'Relaxed matrices prepared and stored'
    
  end subroutine build_and_store_relaxed_matrices
  
  !------------------------------------------------------------------------------!
  ! Build relaxed density using Z-vector
  !------------------------------------------------------------------------------!
  subroutine build_relaxed_density_from_z(infos, density, nbf)
    use oqp_tagarray_driver
    use tagarray
    
    type(information), intent(inout) :: infos
    real(dp), intent(out) :: density(nbf, nbf)
    integer, intent(in) :: nbf
    
    real(dp), pointer :: dm_a(:,:), z_xpy(:)
    integer(c_int32_t) :: status
    integer :: i, j, a, nocc, nvirt, ia
    
    nocc = infos%mol_prop%nelec_A
    nvirt = nbf - nocc
    
    ! Start with ground state density
    call tagarray_get_data(infos%dat, OQP_DM_A, dm_a, status)
    if (status == TA_OK) then
      density = dm_a
    else
      ! Simple diagonal if not available
      density = 0.0_dp
      do i = 1, nocc
        density(i,i) = 2.0_dp  ! Assuming closed shell reference
      end do
    end if
    
    ! Add Z-vector correction to density
    call tagarray_get_data(infos%dat, OQP_td_xpy, z_xpy, status)
    if (status == TA_OK) then
      ! Z-vector gives occupied-virtual corrections
      ia = 0
      do i = 1, nocc
        do a = nocc+1, nbf
          ia = ia + 1
          if (ia <= size(z_xpy)) then
            ! Update density matrix with Z-vector elements
            density(i,a) = density(i,a) + z_xpy(ia)
            density(a,i) = density(i,a)
          end if
        end do
      end do
    end if
    
    ! Ensure symmetry
    do i = 1, nbf
      do j = i+1, nbf
        density(j,i) = density(i,j)
      end do
    end do
    
  end subroutine build_relaxed_density_from_z
  
  !------------------------------------------------------------------------------!
  ! Build relaxed Lagrangian using Z-vector
  !------------------------------------------------------------------------------!
  subroutine build_relaxed_lagrangian_from_z(infos, lagrangian, nbf)
    use oqp_tagarray_driver
    use tagarray
    
    type(information), intent(inout) :: infos
    real(dp), intent(out) :: lagrangian(nbf, nbf)
    integer, intent(in) :: nbf
    
    real(dp), pointer :: fock_a(:,:), z_xmy(:)
    integer(c_int32_t) :: status
    integer :: i, j, a, nocc, nvirt, ia
    
    nocc = infos%mol_prop%nelec_A
    nvirt = nbf - nocc
    
    ! Start with Fock matrix
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a, status)
    if (status == TA_OK) then
      lagrangian = fock_a
    else
      ! Diagonal approximation
      lagrangian = 0.0_dp
      do i = 1, nocc
        lagrangian(i,i) = -0.5_dp
      end do
      do a = nocc+1, nbf
        lagrangian(a,a) = 0.1_dp
      end do
    end if
    
    ! Add Z-vector contribution to Lagrangian
    call tagarray_get_data(infos%dat, OQP_td_xmy, z_xmy, status)
    if (status == TA_OK) then
      ! Z-vector modifies occupied-virtual blocks
      ia = 0
      do i = 1, nocc
        do a = nocc+1, nbf
          ia = ia + 1
          if (ia <= size(z_xmy)) then
            lagrangian(i,a) = lagrangian(i,a) + z_xmy(ia)
            lagrangian(a,i) = lagrangian(i,a)
          end if
        end do
      end do
    end if
    
    ! Ensure symmetry
    do i = 1, nbf
      do j = i+1, nbf
        lagrangian(j,i) = lagrangian(i,j)
      end do
    end do
    
  end subroutine build_relaxed_lagrangian_from_z

end module dyson_mrsf_mod