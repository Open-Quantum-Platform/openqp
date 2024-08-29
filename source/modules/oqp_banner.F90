!> @brief   The initialization of Open Quantum Platform (OpenQP = OQP in source code level)
!> @details This module initialize entire OQP in Fortran side.
!>          It does:
!>          1) Setting up the log file
!>          2) Printing out author information
!>          3) Printing out the basic information regarding OS, date, HW Specs.
!>
!> @param infos(in,out)     Molecule information
module oqp_banner_mod

  character(len=*), parameter :: module_name = "oqp_banner_mod"

contains

  subroutine oqp_banner_C(c_handle) bind(C, name="oqp_banner")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call oqp_banner(inf)
  end subroutine oqp_banner_C

  subroutine oqp_banner(infos)
    use messages, only: show_message, with_abort
    use types, only: information
!$  use omp_lib, only: omp_get_max_threads
    use oqp_tagarray_driver
    use iso_c_binding, only: c_char
    implicit none
    type(information), intent(inout) :: infos
    integer :: iw, CPU_core, i
    character(len=28) :: cdate, hostname

  ! Section of Tagarray for the log filename
  ! We are getting lot file name from Python via tagarray

    character(len=1,kind=c_char), contiguous, pointer :: log_filename(:)
    character(len=*), parameter :: subroutine_name = "oqp_banner"
    character(len=*), parameter :: tags_general(1) = (/ character(len=80) :: &
          OQP_log_filename /)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_log_filename, log_filename)
    allocate(character(ubound(log_filename,1)) :: infos%log_filename)
    do i = 1, ubound(log_filename,1)
       infos%log_filename(i:i) = log_filename(i)
    end do

    open (newunit=iw, file=infos%log_filename, position="append")

    write(iw, '(/,10x, "***********************************************************")')
    write(iw, '(10x,   "*                                                         *")')
    write(iw, '(10x,   "*             OpenQP: Open Quantum Platform               *")')
    write(iw, '(10x,   "*                                                         *")')
    write(iw, '(10x,   "*                Version: 1.0 Aug, 2024                   *")')
    write(iw, '(10x,   "*                                                         *")')
    write(iw, '(10x,   "***********************************************************")')
    write(iw, '(10x,   "*     The most efficient implementation of MRSF-TDDFT.    *")')
    write(iw, '(10x,   "***********************************************************")')
    write(iw, '(10x,   "*                                                         *")')
    write(iw, '(10x,   "*   OpenQP was initiated by Prof. Cheol Ho Choi in 2012.  *")')
    write(iw, '(10x,   "*                                                         *")')
    write(iw, '(10x,   "*   It has since been developed by:                       *")')
    write(iw, '(10x,   "*   Dr. Vladimir Mironov                                  *")')
    write(iw, '(10x,   "*   Dr. Konstantin Komarov                                *")')
    write(iw, '(10x,   "*   Mr. Igor Gerasimov                                    *")')
    write(iw, '(10x,   "*   Dr. Hiroya Nakata                                     *")')
    write(iw, '(10x,   "*                                                         *")')
    write(iw, '(10x,   "*   In 2024, Prof. Jingbai Li at Hoffmann Institute of    *")')
    write(iw, '(10x,   "*   Advanced Materials began developing PyOQP.            *")')
    write(iw, '(10x,   "*                                                         *")')
    write(iw, '(10x,   "***********************************************************")')

    call fdate(cdate)
    call hostnm(hostname)
    CPU_core = 1
  !$  CPU_core = omp_get_max_threads()

    write(iw,'(/20x,A,"Job starts at ",A/,22x," on the host of ",A/,22x," with ",I4," CPU cores.")') ' ', cdate, hostname, CPU_core

    close (iw)

  end subroutine oqp_banner

end module oqp_banner_mod
