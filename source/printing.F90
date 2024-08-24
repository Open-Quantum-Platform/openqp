module printing

  use precision, only: dp

  implicit none

  character(len=*), parameter :: module_name = "printing"

  private
  public :: print_module_info
  public :: print_sym_labeled
  public :: print_mo_range
  public :: print_eigvec_vals_labeled
  public :: print_square
  public :: print_sympack
  public :: print_sym
  public :: print_ev_sol

contains

!> @brief Print MODULE information
!> @detail Printout the information of each MODULES of OQP
  subroutine print_module_info(module_title,module_info)
    use io_constants, only: iw

    implicit none

    character(len=*), intent(in) :: module_title, module_info

    write(iw,'(/20x,40("+")/&
             &23X,"MODULE: ",A/&
             &23X,A/&
             &20X,40("+"))') module_title, module_info

  end subroutine print_module_info

!> @brief Print symmetric packed matrix `d` of dimension `n`
!> @detail The rows will be labeled with basis function tags
  subroutine print_sym_labeled(d, n, basis)
    use precision, only: dp
    use io_constants, only: iw
    use basis_tools, only: basis_set

    implicit none

    real(kind=dp), intent(in) :: d(*)
    integer, intent(in) :: n
    type(basis_set), intent(in) :: basis

    integer, parameter :: maxcolumns = 5
    integer :: i0, ila, i, j0, j1

    do i0 = 1, n, maxcolumns
      write(iw, '(/,15x,*(4x,i4,3x))') (i, i=i0, min(n, i0+maxcolumns-1))
      write(iw,'(G0)')
      ila = 0
      do i = i0, n
        j0 = i0 + i*(i-1)/2
        j1 = j0 + min(ila, maxcolumns-1)
        write (iw,'(i5,2x,a8,*(f11.6))') i, basis%bflab(i), d(j0:j1)
        ila = ila + 1
      end do
    end do
  end subroutine print_sym_labeled

!> @brief Printing out MOs
  subroutine print_mo_range(basis, infos, mostart, moend)
    use io_constants, only: iw
    use types, only: information
    use basis_tools, only: basis_set
    implicit none

    type(basis_set), intent(in) :: basis
    type(information), intent(inout) :: infos
    integer, intent(in) :: mostart, moend
    integer :: mo0, mo1

    write (iw,fmt="(/&
            &10x, 31('=')/&
            &10x, 'Molecular Orbitals and Energies'/&
            &10x, 31('='))")

    mo0 = max(mostart, 1)
    mo1 = min(moend, basis%nbf)
    call print_eigvec_vals_labeled(basis, infos, mo0, mo1)

  end subroutine print_mo_range
!>
!>    @brief    print eigenvector/values, with MO symmetry labels
  subroutine print_eigvec_vals_labeled(basis, infos, mostart, moend)

    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    implicit none

    character(len=*), parameter :: subroutine_name = "print_eigvec_vals_labeled"

    type(basis_set), intent(in) :: basis
    type(information), intent(inout) :: infos

    integer, intent(in) :: mostart, moend

    integer :: i, imax, imin, j, nmax
    character(len=*), parameter :: fmt1 = '(/15x,10(8x,i4,5x))'
    character(len=*), parameter :: fmt2 = '(15x,10f17.10)'
    character(len=*), parameter :: fmt4 = '(i5,2x,a8,10f17.10)'

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      mo_energy_a(:), mo_energy_b(:), mo_a(:,:), mo_b(:,:)
    character(len=*), parameter :: tags_alpha(2) = (/ character(len=80) :: &
      OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(2) = (/ character(len=80) :: &
      OQP_E_MO_B, OQP_VEC_MO_B /)

!Print out eigendata, with mo symmetry labels
!The rows are labeled with the basis function names.

    nmax = 5

    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    write(iw,'(/,A)') '   -------------- Alpha Orbitals -------------'
    do imin = mostart, moend, nmax
      imax = min(imin+nmax-1, moend)

      write(iw,fmt1) (i, i=imin, imax)
      write(iw,fmt2) (mo_energy_a(i), i=imin, imax)

      do j = 1, basis%nbf
        write(iw,fmt4) j, basis%bflab(j), mo_a(j,imin:imax)
      end do

    end do

    if ((infos%mol_prop%nelec_b /= 0) .and. (infos%control%scftype == 2)) then
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
      write(iw,'(/,A)') '   -------------- Beta Orbitals -------------'
      do imin = mostart, moend, nmax
        imax = min(imin+nmax-1, moend)

        write(iw,fmt1) (i, i=imin, imax)
        write(iw,fmt2) (mo_energy_b(i), i=imin, imax)

        do j = 1, basis%nbf
          write(iw,fmt4) j, basis%bflab(j), mo_b(j,imin:imax)
        end do

      end do
    end if

  end subroutine print_eigvec_vals_labeled

!> @brief Print out a square matrix
!> @param[in] v           rectanbular matrix
!> @param[in] m           number of columns in `V`
!> @param[in] n           number of rows in `V`
!> @param[in] ndim        leading dimension of `V`
!> @param[in] tag         optional, will be printed at the beginning of each line
!> @param[in] maxcolumns  optional, number of columns to wrap printing
  subroutine print_square(v, m, n, ndim, tag, maxcolumns)
    use io_constants, only: iw

    implicit none

    real(kind=dp), intent(in) :: V(NDIM, M)
    integer, intent(in) :: m, n, ndim
    integer, optional, intent(in) :: maxcolumns
    character(*), optional, intent(in) :: tag

    character(:), allocatable :: ttag
    integer :: imin, imax, i, j, mxlen

    mxlen = m
    if (present(maxcolumns)) mxlen = maxcolumns
    if (present(tag)) then
        ttag = trim(tag)
    else
        ttag = ""
    end if

    do imin = 1, m, mxlen
      imax = min(imin+mxlen-1, m)
      write (iw, '(a)') ttag
      write (iw, '(a,6x, *(4x, i4, 4x))') ttag, (i, i=imin, imax)
      write (iw, '(a)') ttag
      do j = 1, n
        write (iw, '(a,i5, 1x, *(f12.7))') ttag, j, (v(j, i), i=imin, imax)
      end do
    end do
  end subroutine

!> @brief Print out a symmetric matrix in packed format
!> @param[in] d           symmetric matrix in packed format
!> @param[in] n           matric dimension
  subroutine print_sympack(d, n)
    use io_constants, only: iw

    implicit none

    real(kind=dp), intent(in) :: d(*)
    integer, intent(in) :: n

    integer, parameter :: mxlen = 5
    integer :: i, j, i0, il, j0, jl

    do i0 = 1, n, mxlen
      write (iw,'(/,6x,*(4x,I4,4x))') (i, i=i0, min(n,i0+mxlen-1))
      write (iw,*)
      do i = i0, n
        il = i-i0
        j0 = i0+(i*i-i)/2
        jl = j0+min(il, mxlen-1)
        write (iw, '(i5,1x,*(f12.7))') i, (d(j), j=j0, jl)
      end do
    end do
  end subroutine

!> @brief Print symmetric matrix `D` in square format
!> @param[in]   d   matrix to print, only lower triangle is referenced
!> @param[in]   n   rank of matrix `D`
!> @param[in]   ld  leading dimension of matrix `D`
  subroutine print_sym(d, n, ld)
    use io_constants, only: iw

    implicit none

    real(kind=dp), intent(in) :: d(ld,*)
    integer, intent(in) :: n, ld

    integer, parameter :: mxlen = 5

    integer :: j, i, i0

    do i0 = 1, n, mxlen
      write (iw,fmt='(/,6x,*(4x,i4,4x))') (i, i=i0, min(n,i0+mxlen-1))
      write (iw,*)
      do i = i0, n
        write (iw,fmt='(i5,x,*(f12.7))') i, (d(j,i), j=i0, min(n,i0+min(i-i0,mxlen-1)))
      end do
    end do
  end subroutine

!> @brief Print the solution of the eigenvalue problem
!> @param[in]   v    matrix of eigenvectors, v(ldv,m)
!> @param[in]   e    array of eigenvectors, e(m)
!> @param[in]   m    dimension of column space
!> @param[in]   n    dimension of row space
!> @param[in]   ldv  leading dimension of `V`
  subroutine print_ev_sol(v, e, m, n, ldv)
    use io_constants, only: iw

    implicit none

    real(kind=dp), intent(in) :: v(ldv,m), e(m)
    integer, intent(in) :: m, n, ldv
    integer, parameter :: mxlen = 5

    integer :: i, j, i0, i1

    do  i0 = 1, m, mxlen
      i1 = min(m, i0+mxlen-1)
      write(iw,'(/,15x,*(4x,i4,3x))') (i, i=i0, i1)
      write(iw,'(/,15x,*(f11.6))') (e(i), i=i0, i1)
      write(iw,*)
      do j = 1, n
        write(iw,'(i5,10x,*(f11.6))') j, (v(j,i), i=i0, i1)
      end do
    end do

  end subroutine

end module printing
