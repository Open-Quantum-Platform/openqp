!================================================================================!
! dyson_orbitals.F90
! General theory-independent module for Dyson orbital calculations
! Uses Extended Koopmans' Theorem (EKT) to compute IPs and EAs
! Works with relaxed matrices from Z-vector (not full gradient)
!================================================================================!

module dyson_orbitals_mod
  use precision, only: dp
  use io_constants, only: iw
  use types, only: information
  use printing, only: print_module_info
  use, intrinsic :: iso_c_binding
  
  implicit none
  private
  
  character(len=*), parameter :: module_name = "dyson_orbitals_mod"
  
  public :: dyson_orbital_driver
  public :: compute_dyson_ekt
  public :: store_dyson_results_tagarray
  
contains

  !------------------------------------------------------------------------------!
  ! C interface wrapper
  !------------------------------------------------------------------------------!
  subroutine dyson_orbital_calculation_C(c_handle) bind(C, name="dyson_orbital_calculation")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    
    inf => oqp_handle_get_info(c_handle)
    call dyson_orbital_driver(inf)
    
  end subroutine dyson_orbital_calculation_C
  
  !------------------------------------------------------------------------------!
  ! Main driver for Dyson orbital calculation
  !------------------------------------------------------------------------------!
  subroutine dyson_orbital_driver(infos)
    use messages, only: show_message, with_abort
    
    type(information), intent(inout) :: infos
    
    integer :: nbf, nocca, noccb, nstates, target_state
    real(dp), allocatable :: density_relaxed(:,:,:)
    real(dp), allocatable :: lagrangian(:,:,:)
    real(dp), allocatable :: dyson_orbs(:,:,:)
    real(dp), allocatable :: binding_energies(:,:)
    real(dp), allocatable :: pole_strengths(:,:)
    logical :: compute_ip, compute_ea
    logical :: is_mrsf
    
    call print_module_info('DYSON_EKT', 'Extended Koopmans Theorem Dyson Orbital Analysis')
    
    ! Get parameters from infos
    compute_ip = infos%tddft%dyson_compute_ip
    compute_ea = infos%tddft%dyson_compute_ea
    target_state = infos%tddft%dyson_target_state
    
    if (.not. compute_ip .and. .not. compute_ea) then
      call show_message('Must compute at least IP or EA', with_abort)
    end if
    
    ! Get dimensions
    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    
    ! Determine number of states
    nstates = 1  ! Default to ground state only
    if (infos%control%hamilton == 20) then  ! DFT calculation
      nstates = infos%tddft%nstate
      if (nstates == 0) nstates = 1
    end if
    
    ! Use TDHF target if not specified
    if (target_state == 0) target_state = infos%tddft%target_state
    if (target_state == 0) target_state = 1
    
    ! Check if this is MRSF (already handled by mrsf module) or general
    is_mrsf = (infos%mol_prop%mult == 3 .and. &
              (infos%tddft%mult == 1 .or. infos%tddft%mult == 3 .or. infos%tddft%mult == 5))
    
    write(iw, '(/,5x,A)') 'Dyson Orbital Parameters:'
    write(iw, '(5x,A,L1)') '  Compute IPs: ', compute_ip
    write(iw, '(5x,A,L1)') '  Compute EAs: ', compute_ea
    write(iw, '(5x,A,I0)') '  Target state: ', target_state
    write(iw, '(5x,A,I0)') '  Number of basis functions: ', nbf
    write(iw, '(5x,A,I0)') '  Alpha electrons: ', nocca
    write(iw, '(5x,A,I0)') '  Beta electrons: ', noccb
    write(iw, '(5x,A,E10.3)') '  Pole threshold: ', infos%tddft%dyson_pole_threshold
    
    ! For non-MRSF, Z-vector should already be computed from Python
    if (.not. is_mrsf) then
      if (.not. infos%mol_energy%Z_Vector_converged) then
        call show_message('Z-vector must be computed before Dyson calculation', with_abort)
      end if
      
      write(iw, '(5x,A)') 'Using pre-computed Z-vector'
      call build_and_store_relaxed_matrices_general(infos, nbf, target_state)
    end if
    
    ! Retrieve relaxed density and Lagrangian from tagarray
    call get_relaxed_matrices(infos, density_relaxed, lagrangian, nbf, nstates, target_state)
    
    ! Allocate result arrays
    allocate(dyson_orbs(nbf, nbf, nstates))
    allocate(binding_energies(nbf, nstates))
    allocate(pole_strengths(nbf, nstates))
    
    ! Initialize
    dyson_orbs = 0.0_dp
    binding_energies = 0.0_dp
    pole_strengths = 0.0_dp
    
    ! Compute Dyson orbitals using EKT for target state
    if (target_state > 0 .and. target_state <= nstates) then
      call compute_dyson_ekt(density_relaxed(:,:,target_state), &
                             lagrangian(:,:,target_state), &
                             nbf, infos%tddft%dyson_deflation_tol, &
                             dyson_orbs(:,:,target_state), &
                             binding_energies(:,target_state), &
                             pole_strengths(:,target_state))
    else
      write(iw, '(5x,A,I0)') 'ERROR: Invalid target state: ', target_state
      call show_message('Invalid target state for Dyson calculation', with_abort)
    end if
    
    ! Store results in tagarray
    call store_dyson_results_tagarray(infos, dyson_orbs, binding_energies, &
                                      pole_strengths, nbf, nstates, target_state)
    
    ! Print summary
    call print_dyson_summary(infos, binding_energies(:,target_state), &
                            pole_strengths(:,target_state), nbf)
    
    ! Print detailed results if requested
    if (infos%tddft%dyson_print_level > 0) then
      call print_dyson_results(infos)
    end if
    
    ! Clean up
    deallocate(density_relaxed, lagrangian)
    deallocate(dyson_orbs, binding_energies, pole_strengths)
    
  end subroutine dyson_orbital_driver
  
  !------------------------------------------------------------------------------!
  ! Build and store relaxed matrices for non-MRSF cases
  !------------------------------------------------------------------------------!
  subroutine build_and_store_relaxed_matrices_general(infos, nbf, target_state)
    use oqp_tagarray_driver
    use tagarray
    use messages, only: show_message, without_abort
    
    type(information), intent(inout) :: infos
    integer, intent(in) :: nbf, target_state
    
    real(dp), allocatable :: density_relaxed(:,:)
    real(dp), allocatable :: lagrangian_relaxed(:,:)
    real(dp), pointer :: dm_a(:,:), fock_a(:,:)
    real(dp), pointer :: z_xpy(:), z_xmy(:)
    integer(c_int32_t) :: status
    integer :: i, j, a, nocc, ia
    
    write(iw, '(5x,A)') 'Building relaxed matrices from Z-vector...'
    
    ! Allocate matrices
    allocate(density_relaxed(nbf, nbf))
    allocate(lagrangian_relaxed(nbf, nbf))
    
    nocc = infos%mol_prop%nelec_A
    
    ! Build relaxed density from ground state + Z-vector
    call tagarray_get_data(infos%dat, OQP_DM_A, dm_a, status)
    if (status == TA_OK) then
      density_relaxed = dm_a
    else
      density_relaxed = 0.0_dp
      do i = 1, nocc
        density_relaxed(i,i) = 2.0_dp
      end do
    end if
    
    ! Add Z-vector correction to density
    call tagarray_get_data(infos%dat, OQP_td_xpy, z_xpy, status)
    if (status == TA_OK) then
      ia = 0
      do i = 1, nocc
        do a = nocc+1, nbf
          ia = ia + 1
          if (ia <= size(z_xpy)) then
            density_relaxed(i,a) = density_relaxed(i,a) + z_xpy(ia) * 0.5_dp
            density_relaxed(a,i) = density_relaxed(i,a)
          end if
        end do
      end do
    end if
    
    ! Build relaxed Lagrangian from Fock + Z-vector
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a, status)
    if (status == TA_OK) then
      lagrangian_relaxed = fock_a
    else
      lagrangian_relaxed = 0.0_dp
      do i = 1, nocc
        lagrangian_relaxed(i,i) = -0.5_dp - 0.01_dp * real(nocc - i, dp)
      end do
      do a = nocc+1, nbf
        lagrangian_relaxed(a,a) = 0.1_dp * real(a - nocc, dp)
      end do
    end if
    
    ! Add Z-vector contribution to Lagrangian
    call tagarray_get_data(infos%dat, OQP_td_xmy, z_xmy, status)
    if (status == TA_OK) then
      ia = 0
      do i = 1, nocc
        do a = nocc+1, nbf
          ia = ia + 1
          if (ia <= size(z_xmy)) then
            lagrangian_relaxed(i,a) = lagrangian_relaxed(i,a) + z_xmy(ia) * 0.5_dp
            lagrangian_relaxed(a,i) = lagrangian_relaxed(i,a)
          end if
        end do
      end do
    end if
    
    ! Ensure symmetry
    do i = 1, nbf
      do j = i+1, nbf
        density_relaxed(j,i) = density_relaxed(i,j)
        lagrangian_relaxed(j,i) = lagrangian_relaxed(i,j)
      end do
    end do
    
    ! Store in tagarray using the TA_PUSH_DATA macro pattern
    TA_PUSH_DATA(infos%dat, OQP_density_relaxed, density_relaxed, status)
    if (status /= TA_OK) then
      call show_message('Failed to store relaxed density', without_abort)
    end if
    
    TA_PUSH_DATA(infos%dat, OQP_lagrangian_relaxed, lagrangian_relaxed, status)
    if (status /= TA_OK) then
      call show_message('Failed to store relaxed Lagrangian', without_abort)
    end if
    
    deallocate(density_relaxed, lagrangian_relaxed)
    
    write(iw, '(5x,A)') 'Relaxed matrices prepared and stored'
    
  end subroutine build_and_store_relaxed_matrices_general
  
  !------------------------------------------------------------------------------!
  ! Retrieve relaxed matrices from tagarray
  !------------------------------------------------------------------------------!
  subroutine get_relaxed_matrices(infos, density_relaxed, lagrangian, nbf, nstates, target_state)
    use oqp_tagarray_driver
    use tagarray
    use messages, only: show_message, with_abort
    
    type(information), intent(inout) :: infos
    real(dp), allocatable, intent(out) :: density_relaxed(:,:,:)
    real(dp), allocatable, intent(out) :: lagrangian(:,:,:)
    integer, intent(in) :: nbf, nstates, target_state
    
    real(dp), pointer :: density_ptr(:,:), lagrangian_ptr(:,:)
    integer(c_int32_t) :: status
    integer :: istate
    
    ! Allocate output arrays
    allocate(density_relaxed(nbf, nbf, nstates))
    allocate(lagrangian(nbf, nbf, nstates))
    
    ! Try to retrieve from tagarray
    call tagarray_get_data(infos%dat, OQP_density_relaxed, density_ptr, status)
    
    if (status == TA_OK) then
      ! Found relaxed matrices
      density_relaxed(:,:,target_state) = density_ptr
      
      call tagarray_get_data(infos%dat, OQP_lagrangian_relaxed, lagrangian_ptr, status)
      if (status == TA_OK) then
        lagrangian(:,:,target_state) = lagrangian_ptr
      else
        call show_message('Found density but not Lagrangian', with_abort)
      end if
      
      write(iw, '(5x,A)') 'Retrieved relaxed matrices from tagarray'
    else
      call show_message('Relaxed matrices not available - Z-vector required', with_abort)
    end if
    
    ! Initialize other states to zero (only target state is computed)
    do istate = 1, nstates
      if (istate /= target_state) then
        density_relaxed(:,:,istate) = 0.0_dp
        lagrangian(:,:,istate) = 0.0_dp
      end if
    end do
    
  end subroutine get_relaxed_matrices
  
  !------------------------------------------------------------------------------!
  ! Core EKT solver - Solves Wφ = εDφ
  !------------------------------------------------------------------------------!
  subroutine compute_dyson_ekt(density, lagrangian, nbf, threshold, &
                               dyson_orbs, binding_energies, pole_strengths)
    
    real(dp), intent(in) :: density(nbf,nbf)
    real(dp), intent(in) :: lagrangian(nbf,nbf) 
    integer, intent(in) :: nbf
    real(dp), intent(in) :: threshold
    real(dp), intent(out) :: dyson_orbs(nbf,nbf)
    real(dp), intent(out) :: binding_energies(nbf)
    real(dp), intent(out) :: pole_strengths(nbf)
    
    real(dp), allocatable :: density_work(:,:), lagrangian_work(:,:)
    real(dp), allocatable :: eigenvalues(:), work(:)
    integer :: i, info, lwork
    real(dp) :: norm
    
    write(iw, '(/,5x,A)') 'Solving EKT generalized eigenvalue problem...'
    
    ! Initialize
    dyson_orbs = 0.0_dp
    binding_energies = 0.0_dp
    pole_strengths = 0.0_dp
    
    ! Make working copies
    allocate(density_work(nbf,nbf), lagrangian_work(nbf,nbf))
    density_work = density
    lagrangian_work = lagrangian
    
    ! Allocate workspace
    allocate(eigenvalues(nbf))
    lwork = 8*nbf
    allocate(work(lwork))
    
    ! Solve generalized eigenvalue problem: Wφ = εDφ
    call dsygv(1, 'V', 'U', nbf, lagrangian_work, nbf, &
               density_work, nbf, eigenvalues, work, lwork, info)
    
    if (info /= 0) then
      write(iw, '(5x,A,I0)') 'WARNING: Generalized eigenvalue problem failed, info = ', info
      deallocate(density_work, lagrangian_work, eigenvalues, work)
      return
    end if
    
    write(iw, '(5x,A)') 'EKT equation solved successfully'
    
    ! Store results
    do i = 1, nbf
      binding_energies(i) = eigenvalues(i)
      
      ! Dyson orbitals are the eigenvectors (now in lagrangian_work)
      dyson_orbs(:,i) = lagrangian_work(:,i)
      
      ! Compute pole strength (norm of Dyson orbital)
      norm = dot_product(dyson_orbs(:,i), dyson_orbs(:,i))
      pole_strengths(i) = norm
      
      ! Normalize orbitals
      if (norm > 1.0e-10_dp) then
        dyson_orbs(:,i) = dyson_orbs(:,i) / sqrt(norm)
      end if
    end do
    
    ! Clean up
    deallocate(density_work, lagrangian_work, eigenvalues, work)
    
  end subroutine compute_dyson_ekt
  
  !------------------------------------------------------------------------------!
  ! Store Dyson results in tagarray for Python access
  !------------------------------------------------------------------------------!
  subroutine store_dyson_results_tagarray(infos, dyson_orbs, binding_energies, &
                                          pole_strengths, nbf, nstates, target_state)
    use oqp_tagarray_driver
    use tagarray
    use messages, only: show_message, without_abort
    
    type(information), intent(inout) :: infos
    real(dp), intent(in) :: dyson_orbs(nbf, nbf, nstates)
    real(dp), intent(in) :: binding_energies(nbf, nstates)
    real(dp), intent(in) :: pole_strengths(nbf, nstates)
    integer, intent(in) :: nbf, nstates, target_state
    
    real(dp), allocatable :: vec_do(:,:)      ! Dyson orbital coefficients
    real(dp), allocatable :: e_do(:)          ! Dyson orbital energies
    real(dp), allocatable :: do_strength(:)   ! Pole strengths
    integer(8), allocatable :: do_type(:)     ! Type: 1=IP, -1=EA (use int64)
    integer(8), allocatable :: count_array(:) ! For storing scalar count (use int64)
    integer :: i, n_significant, idx
    real(dp) :: threshold
    integer(c_int32_t) :: status
    
    threshold = infos%tddft%dyson_pole_threshold
    
    ! Count significant results for target state
    n_significant = count(pole_strengths(:,target_state) > threshold)
    
    if (n_significant > 0) then
      ! Allocate arrays for significant orbitals
      allocate(vec_do(nbf, n_significant))
      allocate(e_do(n_significant))
      allocate(do_strength(n_significant))
      allocate(do_type(n_significant))
      
      ! Extract significant results
      idx = 0
      do i = 1, nbf
        if (pole_strengths(i,target_state) > threshold) then
          idx = idx + 1
          vec_do(:,idx) = dyson_orbs(:,i,target_state)
          e_do(idx) = binding_energies(i,target_state)
          do_strength(idx) = pole_strengths(i,target_state)
          
          ! Classify as IP or EA (store as int64)
          if (binding_energies(i,target_state) > 0.0_dp) then
            do_type(idx) = 1_8   ! IP
          else
            do_type(idx) = -1_8  ! EA
          end if
        end if
      end do
      
      ! Store in tagarray using TA_set
      status = TA_set(infos%dat, OQP_VEC_DO, vec_do, comment=OQP_VEC_DO_comment)
      if (status /= TA_OK) call show_message('Failed to store VEC_DO', without_abort)
      
      status = TA_set(infos%dat, OQP_E_DO, e_do, comment=OQP_E_DO_comment)
      if (status /= TA_OK) call show_message('Failed to store E_DO', without_abort)
      
      status = TA_set(infos%dat, OQP_DO_Strength, do_strength, comment=OQP_DO_Strength_comment)
      if (status /= TA_OK) call show_message('Failed to store DO_Strength', without_abort)
      
      status = TA_set(infos%dat, OQP_DO_Type, do_type, comment=OQP_DO_Type_comment)
      if (status /= TA_OK) call show_message('Failed to store DO_Type', without_abort)
      
      ! Store count as int64 array
      allocate(count_array(1))
      count_array(1) = int(n_significant, 8)
      status = TA_set(infos%dat, OQP_DO_Count, count_array, comment=OQP_DO_Count_comment)
      if (status /= TA_OK) call show_message('Failed to store DO_Count', without_abort)
      deallocate(count_array)
      
      ! Clean up
      deallocate(vec_do, e_do, do_strength, do_type)
      
      write(iw, '(/,5x,A,I0,A)') 'Stored ', n_significant, ' significant Dyson orbitals'
    else
      ! Store zero count
      allocate(count_array(1))
      count_array(1) = 0_8
      status = TA_set(infos%dat, OQP_DO_Count, count_array, comment=OQP_DO_Count_comment)
      deallocate(count_array)
      write(iw, '(/,5x,A)') 'No significant Dyson orbitals found above threshold'
    end if
    
  end subroutine store_dyson_results_tagarray
  
  !------------------------------------------------------------------------------!
  ! Print summary of results
  !------------------------------------------------------------------------------!
  subroutine print_dyson_summary(infos, binding_energies, pole_strengths, nbf)
    
    type(information), intent(in) :: infos
    real(dp), intent(in) :: binding_energies(nbf)
    real(dp), intent(in) :: pole_strengths(nbf)
    integer, intent(in) :: nbf
    
    real(dp) :: threshold
    integer :: n_ip, n_ea, i
    
    threshold = infos%tddft%dyson_pole_threshold
    
    ! Count IPs and EAs
    n_ip = 0
    n_ea = 0
    do i = 1, nbf
      if (pole_strengths(i) > threshold) then
        if (binding_energies(i) > 0.0_dp) then
          n_ip = n_ip + 1
        else if (binding_energies(i) < 0.0_dp) then
          n_ea = n_ea + 1
        end if
      end if
    end do
    
    write(iw, '(/,5x,A)') 'Dyson Orbital Summary:'
    write(iw, '(5x,A,I0)') '  Total orbitals computed: ', nbf
    write(iw, '(5x,A,I0)') '  Significant IPs found: ', n_ip
    write(iw, '(5x,A,I0)') '  Significant EAs found: ', n_ea
    write(iw, '(5x,A,E10.3)') '  Threshold used: ', threshold
    
  end subroutine print_dyson_summary
  
  !------------------------------------------------------------------------------!
  ! Print detailed results
  !------------------------------------------------------------------------------!
  subroutine print_dyson_results(infos)
    use oqp_tagarray_driver
    use tagarray
    
    type(information), intent(in) :: infos
    
    real(dp), pointer :: e_do(:), do_strength(:)
    integer(8), pointer :: do_type_int8(:)  ! Use int64 type
    integer(8), pointer :: count_int8(:)    ! Use int64 type
    integer :: n_orbs, i, n_print
    real(dp) :: hartree_to_ev
    integer(c_int32_t) :: status
    
    hartree_to_ev = 27.211396_dp
    
    ! Get count from tagarray - using int64 interface
    call tagarray_get_data(infos%dat, OQP_DO_Count, count_int8, status)
    if (status /= TA_OK .or. .not. associated(count_int8)) return
    if (size(count_int8) < 1) return
    n_orbs = int(count_int8(1))
    if (n_orbs == 0) return
    
    call tagarray_get_data(infos%dat, OQP_E_DO, e_do, status)
    if (status /= TA_OK) return
    
    call tagarray_get_data(infos%dat, OQP_DO_Strength, do_strength, status)
    if (status /= TA_OK) return
    
    call tagarray_get_data(infos%dat, OQP_DO_Type, do_type_int8, status)
    if (status /= TA_OK) return
    
    n_print = min(infos%tddft%dyson_max_print, n_orbs)
    
    write(iw, '(/,5x,A)') repeat('=', 60)
    write(iw, '(5x,A)') 'DYSON ORBITAL RESULTS (EKT)'
    write(iw, '(5x,A)') repeat('=', 60)
    
    ! Print IPs
    write(iw, '(/,5x,A)') 'Ionization Potentials:'
    write(iw, '(5x,A)') '  No.   BE (a.u.)    BE (eV)    Pole Strength   Type'
    write(iw, '(5x,A)') repeat('-', 55)
    
    do i = 1, n_print
      if (do_type_int8(i) == 1) then  ! IP
        write(iw, '(5x,I4,2x,F10.6,2x,F10.4,2x,F12.8,4x,A2)') &
          i, e_do(i), e_do(i)*hartree_to_ev, do_strength(i), 'IP'
      end if
    end do
    
    ! Print EAs
    write(iw, '(/,5x,A)') 'Electron Affinities:'
    write(iw, '(5x,A)') '  No.   BE (a.u.)    EA (eV)    Pole Strength   Type'
    write(iw, '(5x,A)') repeat('-', 55)
    
    do i = 1, n_print
      if (do_type_int8(i) == -1) then  ! EA
        write(iw, '(5x,I4,2x,F10.6,2x,F10.4,2x,F12.8,4x,A2)') &
          i, e_do(i), -e_do(i)*hartree_to_ev, do_strength(i), 'EA'
      end if
    end do
    
    write(iw, '(5x,A)') repeat('=', 60)
    
  end subroutine print_dyson_results

end module dyson_orbitals_mod