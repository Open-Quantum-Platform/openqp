module population_analysis

  implicit none

  character(len=*), parameter :: module_name = "population_analysis"

  integer, parameter :: POP_MULLIKEN = 0
  integer, parameter :: POP_LOWDIN = 1

  private
  public mulliken
  public lowdin
  public run_population_analysis
  public POP_MULLIKEN
  public POP_LOWDIN

!--------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------

  subroutine mulliken_C(c_handle) bind(C, name="mulliken")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info, oqp_handle_refresh_ptr
    use strings, only: Cstring
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call mulliken(inf)
  end subroutine mulliken_C

!--------------------------------------------------------------------------------

  subroutine lowdin_C(c_handle) bind(C, name="lowdin")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info, oqp_handle_refresh_ptr
    use strings, only: Cstring
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call lowdin(inf)
  end subroutine lowdin_C

!--------------------------------------------------------------------------------

  subroutine mulliken(infos)
    use precision, only: dp
    use io_constants, only: iw
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: orbital_pop(:)
    real(kind=dp), allocatable :: chg(:)
    integer :: nat, ok

    open (unit=IW, file=infos%log_filename, position="append")

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms
    nat = ubound(infos%atoms%zn,1)

    allocate(orbital_pop(basis%nbf), &
             chg(nat), &
             source=0.0d0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    write(iw,'(2/)')
    write(iw,'(4x,a)') '============================'
    write(iw,'(4x,a)') 'Mulliken population analysis'
    write(iw,'(4x,a)') '============================'
    call flush(iw)

    call run_population_analysis(infos, basis, orbital_pop, chg, POP_MULLIKEN)

    write(iw,'(/,2X,A)') 'Gross AO population (Mulliken)'
    call print_ao_pop(infos, orbital_pop)

    write(iw,'(/,2X,A)') 'Atomic partial charges (Mulliken)'
    call print_charges(infos, chg)

    close(iw)

  end subroutine mulliken

!--------------------------------------------------------------------------------

  subroutine lowdin(infos)
    use precision, only: dp
    use io_constants, only: iw
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: orbital_pop(:)
    real(kind=dp), allocatable :: chg(:)
    integer :: nat, ok

    open (unit=IW, file=infos%log_filename, position="append")

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms
    nat = ubound(infos%atoms%zn,1)

    allocate(orbital_pop(basis%nbf), &
             chg(nat), &
             source=0.0d0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    write(iw,'(2/)')
    write(iw,'(4x,a)') '=========================='
    write(iw,'(4x,a)') 'Lowdin population analysis'
    write(iw,'(4x,a)') '=========================='
    call flush(iw)

    call run_population_analysis(infos, basis, orbital_pop, chg, POP_LOWDIN)

    write(iw,'(/,2X,A)') 'Gross AO population (Lowdin)'
    call print_ao_pop(infos, orbital_pop)

    write(iw,'(/,2X,A)') 'Atomic partial charges (Lowdin)'
    call print_charges(infos, chg)

    close(iw)

  end subroutine lowdin

!--------------------------------------------------------------------------------

!> @brief Run population analysis
!>
!> @param[in]      infos        QOP handle, fortran
!> @param[in]      basis        basis set
!> @param[out]     orbital_pop  AO population
!> @param[out]     chg          atomic partial charges
!> @param[in]      sel          Mulliken(0) or Lowdin(1) analysis selector
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine run_population_analysis(infos, basis, orbital_pop, chg, sel)
    use precision, only: dp
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use mathlib, only: unpack_matrix
    implicit none

    character(len=*), parameter :: subroutine_name = "run_population_analysis"

    type(information), intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), allocatable, intent(inout) :: orbital_pop(:), chg(:)
    integer, intent(in) :: sel

    integer :: nbf, nbf2, ok
    logical :: urohf
    real(kind=dp), allocatable :: dens(:,:), tmp(:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:), smat(:)
    character(len=*), parameter :: tags_alpha(2) = (/ character(len=80) :: &
      OQP_SM, OQP_DM_A /)
    character(len=*), parameter :: tags_beta(1) = (/ character(len=80) :: &
      OQP_DM_B /)

    ! Allocate memory
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    urohf = infos%control%scftype == 2 .or. infos%control%scftype == 3

    allocate(tmp(nbf2), dens(nbf,nbf), source=0.0d0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_SM, smat)

    tmp = dmat_a
    if (urohf) then
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      tmp = tmp + dmat_b
    end if

    ! Make density matrix
    call unpack_matrix(tmp, dens)

    ! Run the analysis
    select case(sel)
    case (POP_MULLIKEN)
      call get_orb_pop_mulliken(smat, dens, orbital_pop)
    case (POP_LOWDIN)
      call get_orb_pop_lowdin(smat, dens, orbital_pop)
    case default
      call show_message('Unknown population analysis method', WITH_ABORT)
    end select

    ! Get electronic populations on aotms
    call get_atomic_pop(basis, orbital_pop, chg)

    ! Compute partial charges
    chg = infos%atoms%zn - chg

  end subroutine

!--------------------------------------------------------------------------------

!> @brief Compute Mulliken's atomic orbital population
!>

!> @brief Compute atomic population from orbital population
!>
!> @param[in]      s            overlap matrix, packed
!> @param[in]      d            density matrix, square
!> @param[out]     pop          AO population
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine get_atomic_pop(basis, ao_pop, at_pop)
    use precision, only: dp
    use basis_tools, only: basis_set
    use constants, only: NUM_CART_BF
    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: ao_pop(:)
    real(kind=dp), intent(out) :: at_pop(:)

    integer :: i, iatom, i0, i1

    do i = 1, basis%nshell
      iatom = basis%origin(i)
      i0 = basis%ao_offset(i)
      i1 = basis%ao_offset(i)+NUM_CART_BF(basis%am(i))-1
      at_pop(iatom) = at_pop(iatom) + sum(ao_pop(i0:i1))
    end do
  end subroutine

!--------------------------------------------------------------------------------

!> @brief Compute Mulliken's atomic orbital population
!>
!> @param[in]      s            overlap matrix, packed
!> @param[in]      d            density matrix, square
!> @param[out]     pop          AO population
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine get_orb_pop_mulliken(s, d, pop)
    use precision, only: dp
    use mathlib, only: unpack_matrix
    implicit none

    real(kind=dp), intent(in) :: s(:), d(:,:)
    real(kind=dp), intent(out) :: pop(:)

    real(kind=dp), allocatable :: stmp(:,:)
    integer :: nbf

    nbf = ubound(d,1)
    allocate(stmp(nbf,nbf))
    call unpack_matrix(s, stmp)

    pop = sum(d*stmp,dim=2)
  end subroutine

!--------------------------------------------------------------------------------

!> @brief Compute orbital population using Lowdin's symmetrically
!>        orthogonalized basis set
!>
!> @param[in]      s            overlap matrix, packed
!> @param[in]      d            density matrix, square
!> @param[out]     pop          AO population
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine get_orb_pop_lowdin(s, d, pop)
    use precision, only: dp
    use eigen, only: diag_symm_packed
    use mathlib, only: unpack_matrix
    use oqp_linalg
    implicit none

    real(kind=dp), intent(in) :: s(:), d(:,:)
    real(kind=dp), intent(out) :: pop(:)

    real(kind=dp), allocatable :: stmp(:), tmp(:,:), &
            eval(:), evec(:,:), ovlsqrt(:,:)
    integer :: ierr
    integer :: nbf, i

    nbf = ubound(d,1)
    allocate(eval(nbf), evec(nbf,nbf), tmp(nbf,nbf), ovlsqrt(nbf,nbf))

    ! Compute S^{1/2}
    stmp = s
    call diag_symm_packed(1, nbf, nbf, nbf, stmp, eval, evec, ierr)

    do i = 1, nbf
      tmp(:,i) = evec(:,i)*sqrt(eval(i))
    end do

    call dgemm('n','t',nbf,nbf,nbf, &
            1.0d0, evec, nbf, &
                   tmp, nbf, &
            0.0d0, ovlsqrt, nbf)

    ! Compute Tr(S^{1/2}*dens*S^{1/2})
    call dgemm('n','t',nbf,nbf,nbf, &
            1.0d0, ovlsqrt, nbf, &
                   d, nbf, &
            0.0d0, tmp, nbf)

    pop = sum(tmp*ovlsqrt,dim=2)

  end subroutine

!--------------------------------------------------------------------------------

!> @brief Print partial charges
!>
!> @param[in]      infos        OQP handle
!> @param[in]      chg          atomic partial charges
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine print_charges(infos, chg)
    use precision, only: dp
    use elements, only: ELEMENTS_SHORT_NAME
    use types, only: information
    type(information), intent(in) :: infos
    real(kind=dp), intent(in) :: chg(:)

    integer :: i, elem, nat
    nat = ubound(infos%atoms%zn, 1)

    write(*,'(/,30("^"))')
    write(*,'(/a8,a8,a14)') '#', 'Name', 'Charge'
    write(*,'(30("-"))')

    do i = 1, nat
      elem = nint(infos%atoms%zn(i))
      write(*,'(i8,a8,f14.6)') i, ELEMENTS_SHORT_NAME(elem), chg(i)
    end do

    write(*,'(30("="))')

  end subroutine print_charges

!--------------------------------------------------------------------------------

!> @brief Print gross AO population
!>
!> @param[in]      infos        OQP handle
!> @param[in]      pop          atomic orbital population
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine print_ao_pop(infos, pop)
    use precision, only: dp
    use types, only: information
    type(information), intent(in) :: infos
    real(kind=dp), intent(in) :: pop(:)

    integer :: i, nbf
    nbf = infos%basis%nbf

    write(*,'(/,34("^"))')
    write(*,'(/a8,a11,a15)') '#', 'A  N  L', 'Population'
    write(*,'(34("-"))')

    do i = 1, nbf
      write(*,'(i8,a12,f14.6)') i, infos%basis%bf_label(i), pop(i)
    end do

    write(*,'(34("="))')

  end subroutine print_ao_pop

!--------------------------------------------------------------------------------

end module population_analysis
