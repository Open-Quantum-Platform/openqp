!> Native DFT-D4 dispersion interface for OpenQP.
!>
!> Thin bind(C) shim over the dynamically linked dftd4 Fortran library.
!> Replaces the former `dftd4` Python (cffi) package, which is capped at Python
!> <= 3.12. This routine has no Python dependency and is called directly
!> through the existing cffi boundary (see include/oqp.h).
! Kind of the ordinary integer that mctc-lib uses for its `num` argument.  The
! DFT-D4 stack deliberately keeps four-byte default integers; WITH_ILP64 changes
! only its BLAS interface. OpenQP compiles this file with -fdefault-integer-8,
! so the shim converts `num` explicitly. Set in source/CMakeLists.txt.
#ifndef OQP_D4_INT_KIND
#define OQP_D4_INT_KIND 4
#endif

module dftd4_interface
  use iso_c_binding
  use iso_fortran_env, only: wp => real64
  use mctc_io, only: structure_type, new
  use dftd4, only: d4_model, new_d4_model, damping_param, &
                   get_rational_damping, get_dispersion, realspace_cutoff
  use dftd4_damping_rational, only: rational_damping_param
  implicit none
  private
  public :: oqp_dftd4_disp, oqp_dftd4_disp_v2

contains

  !> Legacy ABI: compute neutral-system DFT-D4 dispersion using parameters
  !> selected by functional name.  Keep this symbol and its numerical behavior
  !> for existing callers; charge-aware callers should use oqp_dftd4_disp_v2.
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

    real(c_double), parameter :: legacy_damping(6) = 0.0_c_double

    call dftd4_disp_impl(nat, z, xyz, 0.0_c_double, func, lfunc, &
                         0_c_int, legacy_damping, do_grad, energy, grad, ier)
  end subroutine oqp_dftd4_disp


  !> Charge-aware DFT-D4 ABI with optional explicit rational damping.
  !>
  !>  total_charge  molecular total charge in units of |e|
  !>  param_mode    0: load parameters by functional name
  !>                1: use damping = [s6, s8, s9, a1, a2, alp]
  !>  ier           0 on success, 1 for unknown functional, 2 for bad mode
  subroutine oqp_dftd4_disp_v2(nat, z, xyz, total_charge, func, lfunc, &
       param_mode, damping, do_grad, energy, grad, ier) &
       bind(C, name="oqp_dftd4_disp_v2")
    integer(c_int), value :: nat, lfunc, param_mode, do_grad
    integer(c_int), intent(in) :: z(nat)
    real(c_double), intent(in) :: xyz(3, nat)
    real(c_double), value :: total_charge
    character(kind=c_char), intent(in) :: func(lfunc)
    real(c_double), intent(in) :: damping(6)
    real(c_double), intent(out) :: energy
    real(c_double), intent(out) :: grad(3, nat)
    integer(c_int), intent(out) :: ier

    call dftd4_disp_impl(nat, z, xyz, total_charge, func, lfunc, &
                         param_mode, damping, do_grad, energy, grad, ier)
  end subroutine oqp_dftd4_disp_v2


  subroutine dftd4_disp_impl(nat, z, xyz, total_charge, func, lfunc, &
       param_mode, damping, do_grad, energy, grad, ier)
    integer(c_int), intent(in) :: nat, lfunc, param_mode, do_grad
    integer(c_int), intent(in) :: z(nat)
    real(c_double), intent(in) :: xyz(3, nat), total_charge
    character(kind=c_char), intent(in) :: func(lfunc)
    real(c_double), intent(in) :: damping(6)
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

    call new(mol, num=int(z, OQP_D4_INT_KIND), xyz=real(xyz, wp), &
             charge=real(total_charge, wp))
    call new_d4_model(disp, mol)
    select case (param_mode)
    case (0)
      call get_rational_damping(trim(fname), param, s9=1.0_wp)
      if (.not. allocated(param)) then
        ier = 1
        return
      end if
    case (1)
      allocate(rational_damping_param :: param)
      select type (param)
      type is (rational_damping_param)
        param = rational_damping_param( &
          s6=real(damping(1), wp), s8=real(damping(2), wp), &
          s9=real(damping(3), wp), a1=real(damping(4), wp), &
          a2=real(damping(5), wp), alp=real(damping(6), wp))
      end select
    case default
      ier = 2
      return
    end select

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
  end subroutine dftd4_disp_impl

end module dftd4_interface
