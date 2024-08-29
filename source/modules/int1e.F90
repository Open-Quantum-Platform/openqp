module int1e_mod

  implicit none

  character(len=*), parameter :: module_name = "int1e_mod"

  private
  public int1e

contains

  subroutine int1e_C(c_handle) bind(C, name="int1e")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call int1e(inf)
  end subroutine int1e_C

!> @brief Calculate the basic H, S, and T 1e-integrals
  subroutine int1e(infos)

    use types, only: information
    use oqp_tagarray_driver
    use precision, only: dp
    use io_constants, only: iw
    use constants, only: tol_int
    use int1, only: omp_hst
    use basis_tools, only: basis_set
    use printing, only: print_sym_labeled
    use messages, only: show_message, WITH_ABORT
    use strings, only: Cstring, fstring
    use physical_constants, only: BOHR_TO_ANGSTROM
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "int1e"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    real(kind=dp) :: tol
    integer :: i, nbf, nat, nbf2

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      hcore(:), tmat(:), smat(:)
    character(len=*), parameter :: tags_general(3) = (/ character(len=80) :: &
      OQP_SM, OQP_TM, OQP_Hcore /)

    logical dbg
    dbg = .false.

!   Files open:
!   LOG: Read and Write: Main output file
    open(unit=iw,  file=infos%log_filename, position="append")

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms
!
    call print_module_info('int1e','Computing H, S and T Matrices')

!   Print out the Cartesian Coordinates.
    write(iw, fmt="(&
              &/21X,32('=')&
              &/21X,a&
              &/21X,32('=')&
              &/8X,'ATOM     ZNUC',11X,'X',14X,'Y',14X,'Z'&
              &/6X,62('-'))") "Cartesian Coordinate in Angstrom"

    do i = 1, size(basis%atoms%zn(:))
       write(iw,'(7x,i4,5x,f4.1,3(x,f15.9))') &
               i, basis%atoms%zn(i), basis%atoms%xyz(1:3,i)*BOHR_TO_ANGSTROM
    end do

!   Allocate H, S and T matrices
    nbf2 = basis%nbf*(basis%nbf+1)/2

    call infos%dat%remove_records(tags_general)

    call infos%dat%reserve_data(OQP_SM, TA_TYPE_REAL64, nbf2, comment=OQP_SM_comment)
    call infos%dat%reserve_data(OQP_TM, TA_TYPE_REAL64, nbf2, comment=OQP_TM_comment)
    call infos%dat%reserve_data(OQP_Hcore, TA_TYPE_REAL64, nbf2, comment=OQP_Hcore_comment)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_TM, tmat)
    call tagarray_get_data(infos%dat, OQP_Hcore, Hcore)

!   Create arrays of atomic coordinates and charges for one-electron code
    nbf = basis%nbf
    nat = ubound(infos%atoms%zn,1)

!   Compute conventional H, S, and T integrals
    tol = log(10.0d0)*tol_int
    call omp_hst(basis, infos%atoms%xyz, infos%atoms%zn, hcore, smat, tmat, logtol=tol)

    if (dbg) then
        write(iw,'(/"BARE NUCLEUS HAMILTONIAN INTEGRALS (H=T+V)")')
        call print_sym_labeled(hcore,nbf, basis)

        write(iw,'(/"OVERLAP MATRIX")')
        call print_sym_labeled(Smat,nbf, basis)

        write(iw,'(/"KINETIC ENERGY INTEGRALS")')
        call print_sym_labeled(tmat,nbf, basis)
    end if

    write(iw,"(/1X,'...... End Of One Electron Integrals ......'/)")

    close(iw)

  end subroutine int1e

end module int1e_mod
