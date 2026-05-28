module tdhf_mrsf_ekt_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_ekt_mod"

contains

  subroutine tdhf_mrsf_ekt_ip_C(c_handle) bind(C, name="tdhf_mrsf_ekt_ip")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_ekt(inf, .false.)
  end subroutine tdhf_mrsf_ekt_ip_C

  subroutine tdhf_mrsf_ekt_ea_C(c_handle) bind(C, name="tdhf_mrsf_ekt_ea")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_ekt(inf, .true.)
  end subroutine tdhf_mrsf_ekt_ea_C

  subroutine tdhf_mrsf_ekt(infos, electron_affinity)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time
    use precision, only: dp
    use mathlib, only: unpack_matrix, pack_matrix, orthogonal_transform_sym
    use eigen, only: diag_symm_full
    use printing, only: print_module_info
    use tdhf_mrsf_energy_mod, only: tdhf_mrsf_energy

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_ekt"
    type(information), target, intent(inout) :: infos
    logical, intent(in) :: electron_affinity

    type(basis_set), pointer :: basis
    integer :: nbf, nbf2, nroot, ok, i, j, k, ierr
    integer :: nactive
    integer, allocatable :: active(:)
    real(kind=dp), parameter :: metric_tol = 1.0e-8_dp
    real(kind=dp), allocatable :: dens_pack(:), lag_pack(:)
    real(kind=dp), allocatable :: density_mo(:,:), fock_mo(:,:), lagrangian_mo(:,:)
    real(kind=dp), allocatable :: ekt_metric(:,:), ekt_operator(:,:)
    real(kind=dp), allocatable :: op_active(:,:), metric_active(:,:), eig(:), vec_active(:,:)
    real(kind=dp), allocatable :: orbitals(:,:), strengths(:)
    real(kind=dp), contiguous, pointer :: fock_a(:), dmat_a(:), mo_a(:,:)
    real(kind=dp), contiguous, pointer :: density_store(:,:), lagrangian_store(:,:)
    real(kind=dp), contiguous, pointer :: fock_store(:,:), orbital_store(:,:), strength_store(:)
    character(len=*), parameter :: tags_required(3) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_DM_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_alloc(5) = (/ character(len=80) :: &
      OQP_mrsf_ekt_density_mo, OQP_mrsf_ekt_lagrangian_mo, OQP_mrsf_ekt_fock_mo, &
      OQP_mrsf_ekt_orbitals_mo, OQP_mrsf_ekt_strengths /)

    ! Run the parent MRSF calculation first so the selected neutral state and
    ! response vectors are available.  The EKT step then builds the IP/EA
    ! generalized eigenproblem in the MRSF molecular-orbital basis.
    call tdhf_mrsf_energy(infos)

    open(unit=iw, file=infos%log_filename, position="append")
    if (electron_affinity) then
      call print_module_info('MRSF_EKT_EA','Computing MRSF-EKT electron affinities')
    else
      call print_module_info('MRSF_EKT_IP','Computing MRSF-EKT ionization potentials')
    end if

    basis => infos%basis
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nroot = max(1, min(infos%tddft%nstate, nbf))

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    allocate(dens_pack(nbf2), lag_pack(nbf2), density_mo(nbf,nbf), fock_mo(nbf,nbf), &
             lagrangian_mo(nbf,nbf), ekt_metric(nbf,nbf), ekt_operator(nbf,nbf), &
             orbitals(nbf,nroot), strengths(nroot), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    allocate(active(nbf), source=0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    ! MRSF-EKT uses the target-state one-particle density P and the
    ! energy-weighted density / Lagrangian W in the neutral-state MO basis.
    ! This first OpenQP port follows the GAMESS EKT algebra with the currently
    ! available spin-adapted ROHF alpha density and Fock/Lagrangian proxy.
    call orthogonal_transform_sym(nbf, nbf, dmat_a, mo_a, nbf, dens_pack)
    call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, lag_pack)
    call unpack_matrix(dens_pack, density_mo, nbf, 'U')
    call unpack_matrix(lag_pack, fock_mo, nbf, 'U')

    lagrangian_mo = 0.0_dp
    do j = 1, nbf
      do i = 1, nbf
        do k = 1, nbf
          lagrangian_mo(i,j) = lagrangian_mo(i,j) + fock_mo(i,k)*density_mo(k,j)
        end do
      end do
    end do
    lagrangian_mo = 0.5_dp*(lagrangian_mo + transpose(lagrangian_mo))

    if (electron_affinity) then
      ! EKT-EA: (F - W) * x = (I - P) * x * lambda
      ekt_metric = -density_mo
      do i = 1, nbf
        ekt_metric(i,i) = ekt_metric(i,i) + 1.0_dp
      end do
      ekt_operator = fock_mo - lagrangian_mo
    else
      ! EKT: W * x = P * x * lambda
      ekt_metric = density_mo
      ekt_operator = lagrangian_mo
    end if

    nactive = 0
    do i = 1, nbf
      if (ekt_metric(i,i) > metric_tol) then
        nactive = nactive + 1
        active(nactive) = i
      end if
    end do
    if (nactive == 0) call show_message('MRSF-EKT metric has no active subspace', WITH_ABORT)

    allocate(op_active(nactive,nactive), metric_active(nactive,nactive), eig(nactive), &
             vec_active(nactive,nactive), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    do j = 1, nactive
      do i = 1, nactive
        op_active(i,j) = ekt_operator(active(i), active(j))
        metric_active(i,j) = ekt_metric(active(i), active(j))
      end do
    end do

    call solve_symmetric_generalized(op_active, metric_active, eig, vec_active, nactive, ierr)
    if (ierr /= 0) call show_message('(A,I0)', 'MRSF-EKT generalized eigensolver failed; info=', ierr, WITH_ABORT)

    orbitals = 0.0_dp
    do j = 1, min(nroot, nactive)
      do i = 1, nactive
        orbitals(active(i),j) = vec_active(i,j)
      end do
      strengths(j) = sum(orbitals(:,j)**2)
    end do

    call infos%dat%remove_records(tags_alloc)
    call infos%dat%reserve_data(OQP_mrsf_ekt_density_mo, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), &
        comment='MRSF-EKT metric density P in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_lagrangian_mo, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), &
        comment='MRSF-EKT Lagrangian W in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_fock_mo, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), &
        comment='MRSF-EKT Fock matrix in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_orbitals_mo, TA_TYPE_REAL64, nbf*nroot, (/ nbf, nroot /), &
        comment='MRSF-EKT Dyson-like orbital coefficients in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_strengths, TA_TYPE_REAL64, nroot, &
        comment='MRSF-EKT pole strengths')

    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_density_mo, density_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_lagrangian_mo, lagrangian_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_fock_mo, fock_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_orbitals_mo, orbital_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_strengths, strength_store)
    density_store = density_mo
    lagrangian_store = lagrangian_mo
    fock_store = fock_mo
    orbital_store = orbitals
    strength_store = strengths

    if (electron_affinity) then
      write(iw,'(/,2x,"MRSF-EKT electron affinities (hartree)")')
    else
      write(iw,'(/,2x,"MRSF-EKT ionization potentials (hartree)")')
    end if
    write(iw,'(2x,"state",5x,"energy",11x,"strength")')
    do i = 1, min(nroot, nactive)
      write(iw,'(2x,I5,2x,F16.8,2x,F16.8)') i, eig(i), strengths(i)
    end do
    call flush(iw)

    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine tdhf_mrsf_ekt

  subroutine solve_symmetric_generalized(operator, metric, eigenvalues, eigenvectors, n, ierr)
    use precision, only: dp
    use eigen, only: diag_symm_full
    implicit none
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: operator(n,n)
    real(kind=dp), intent(in) :: metric(n,n)
    real(kind=dp), intent(out) :: eigenvalues(n), eigenvectors(n,n)
    integer, intent(out) :: ierr
    integer :: i, j
    real(kind=dp), allocatable :: reduced(:,:)

    allocate(reduced(n,n))
    reduced = operator
    do j = 1, n
      do i = 1, n
        if (metric(i,i) > 1.0e-12_dp .and. metric(j,j) > 1.0e-12_dp) then
          reduced(i,j) = reduced(i,j) / sqrt(metric(i,i)*metric(j,j))
        end if
      end do
    end do
    reduced = 0.5_dp*(reduced + transpose(reduced))
    call diag_symm_full(1, n, reduced, n, eigenvalues, ierr)
    eigenvectors = reduced
  end subroutine solve_symmetric_generalized

end module tdhf_mrsf_ekt_mod
