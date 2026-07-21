!> Native DFT-D4 dispersion interface for OpenQP.
!>
!> Thin bind(C) shim over the dftd4 Fortran library (statically linked into
!> liboqp). Replaces the former `dftd4` Python (cffi) package, which is capped
!> at Python <= 3.12. This routine has no Python dependency and is called
!> directly through the existing cffi boundary (see include/oqp.h).
! Kind of the integer that the dftd4 stack (mctc-lib) was compiled with for its
! `num` argument: it tracks OpenQP's BLAS integer size (BLA_SIZEOF_INTEGER) so
! the dftd4 deps are built LP64/ILP64 to match the host. OpenQP compiles this
! file with -fdefault-integer-8, so we must hand mctc the right-kind integer
! explicitly rather than relying on the default. Set in source/CMakeLists.txt.
#ifndef OQP_D4_INT_KIND
#define OQP_D4_INT_KIND 4
#endif

module dftd4_interface
  use iso_c_binding
  use iso_fortran_env, only: wp => real64
  use mctc_io, only: structure_type, new
  use dftd4, only: d4_model, new_d4_model, damping_param, &
                   get_rational_damping, get_dispersion, realspace_cutoff
  implicit none
  private
  public :: oqp_dftd4_disp

contains

  !> Compute the DFT-D4 dispersion energy and (optionally) nuclear gradient.
  !>
  !>  nat      number of atoms
  !>  z        atomic numbers, z(nat)
  !>  xyz      Cartesian coordinates in Bohr, xyz(3, nat)
  !>  func     functional name (e.g. "pbe0"), not null-terminated
  !>  lfunc    number of characters in func
  !>  do_grad  /= 0 to also evaluate the gradient
  !>  energy   dispersion energy in Hartree
  !>  grad     dispersion gradient in Hartree/Bohr, grad(3, nat) (zeroed if no grad)
  !>  ier      0 on success, 1 if the functional has no D4 damping parameters
  subroutine oqp_dftd4_disp(nat, z, xyz, func, lfunc, do_grad, energy, grad, ier) &
       bind(C, name="oqp_dftd4_disp")
    integer(c_int), value :: nat, lfunc, do_grad
    integer(c_int), intent(in) :: z(nat)
    real(c_double), intent(in) :: xyz(3, nat)
    character(kind=c_char), intent(in) :: func(lfunc)
    real(c_double), intent(out) :: energy
    real(c_double), intent(out) :: grad(3, nat)
    integer(c_int), intent(out) :: ier

    type(structure_type) :: mol
    type(d4_model) :: disp
    class(damping_param), allocatable :: param
    character(len=:), allocatable :: fname
    real(wp), allocatable :: g(:, :)
    real(wp) :: e, sig(3, 3)
    integer :: i

    ier = 0
    energy = 0.0_c_double
    grad = 0.0_c_double

    fname = ''
    do i = 1, lfunc
      fname = fname // func(i)
    end do

    call new(mol, num=int(z, OQP_D4_INT_KIND), xyz=real(xyz, wp))
    call new_d4_model(disp, mol)
    call get_rational_damping(trim(fname), param, s9=1.0_wp)
    if (.not. allocated(param)) then
      ier = 1
      return
    end if

    if (do_grad /= 0) then
      allocate(g(3, nat))
      ! NB: dftd4's get_dispersion writes to `sigma` unconditionally whenever a
      ! gradient is requested, despite it being declared optional (see the dftd4
      ! C API in src/dftd4/api.f90). It MUST be supplied here or it SIGBUSes.
      call get_dispersion(mol, disp, param, realspace_cutoff(), e, gradient=g, sigma=sig)
      grad = real(g, c_double)
    else
      call get_dispersion(mol, disp, param, realspace_cutoff(), e)
    end if
    energy = real(e, c_double)
  end subroutine oqp_dftd4_disp

end module dftd4_interface
