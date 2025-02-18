module electric_moments_mod

  implicit none

  character(len=*), parameter :: module_name = "electric_moments_mod"

  private
  public electric_moments

contains

  subroutine electric_moments_C(c_handle) bind(C, name="electric_moments")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call electric_moments(inf)
  end subroutine electric_moments_c

  subroutine electric_moments(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring
    use mathlib, only: traceprod_sym_packed, triangular_to_full
    use int1, only: multipole_integrals
    use physical_constants, only: AU_TO_DEBYE, AU_TO_BUCK, AU_TO_OCT
    use xyz_order

    implicit none

    character(len=*), parameter :: subroutine_name = "electric_moments"

    type(information), target, intent(inout) :: infos

    integer :: nbf, nbf2, ok
    logical :: urohf
    type(basis_set), pointer :: basis
    real(kind=8), allocatable :: mints(:,:)
    real(kind=8) :: com(3), dip(3), qxyz(6), quad(3,3), dr(3), z
    real(kind=8) :: oxyz(10), oct(10)
    integer :: nat, i

    real(kind=8), contiguous, pointer :: dmat_a(:), dmat_b(:)
    integer(4) :: status

    urohf = infos%control%scftype == 2 .or. infos%control%scftype == 3

!   3. LOG: Write: Main output file
    open (unit=IW, file=infos%log_filename, position="append")

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

!   Allocate memory
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    write(iw,'(2/)')
    write(iw,'(4x,a)') '========================'
    write(iw,'(4x,a)') 'Electric moment analysis'
    write(iw,'(4x,a)') '========================'
    call flush(iw)

    nat = ubound(basis%atoms%zn,1)
    allocate(mints(nbf2,19), source=0.0d0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
    call check_status(status, module_name, subroutine_name, OQP_DM_A)

    if (urohf) then
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b, status)
      call check_status(status, module_name, subroutine_name, OQP_DM_B)
    endif

    com = 0
    do i = 1, nat
        com = com + basis%atoms%xyz(:,i)*basis%atoms%mass(i)
    end do
    com = com / sum(basis%atoms%mass)

    call multipole_integrals(basis, mints, com, 3)

    dip = 0
    qxyz = 0
    oxyz = 0

    do i = 1, nat
      dr = infos%atoms%xyz(:,i)-com
      z = infos%atoms%zn(i) - infos%basis%ecp_zn_num(i)

      dip = dip + z * dr

      qxyz(XX_) = qxyz(XX_) + z * dr(X__)*dr(X__)
      qxyz(YY_) = qxyz(YY_) + z * dr(Y__)*dr(Y__)
      qxyz(ZZ_) = qxyz(ZZ_) + z * dr(Z__)*dr(Z__)
      qxyz(XY_) = qxyz(XY_) + z * dr(X__)*dr(Y__)
      qxyz(XZ_) = qxyz(XZ_) + z * dr(X__)*dr(Z__)
      qxyz(YZ_) = qxyz(YZ_) + z * dr(Y__)*dr(Z__)

      oxyz(XXX) = oxyz(XXX) + z * dr(X__)*dr(X__)*dr(X__)
      oxyz(YYY) = oxyz(YYY) + z * dr(Y__)*dr(Y__)*dr(Y__)
      oxyz(ZZZ) = oxyz(ZZZ) + z * dr(Z__)*dr(Z__)*dr(Z__)
      oxyz(XXY) = oxyz(XXY) + z * dr(X__)*dr(X__)*dr(Y__)
      oxyz(XXZ) = oxyz(XXZ) + z * dr(X__)*dr(X__)*dr(Z__)
      oxyz(YYX) = oxyz(YYX) + z * dr(Y__)*dr(Y__)*dr(X__)
      oxyz(YYZ) = oxyz(YYZ) + z * dr(Y__)*dr(Y__)*dr(Z__)
      oxyz(ZZX) = oxyz(ZZX) + z * dr(Z__)*dr(Z__)*dr(X__)
      oxyz(ZZY) = oxyz(ZZY) + z * dr(Z__)*dr(Z__)*dr(Y__)
      oxyz(XYZ) = oxyz(XYZ) + z * dr(X__)*dr(Y__)*dr(Z__)
    end do

    do i = 1, 3
        dip(i) = dip(i) - traceprod_sym_packed(mints(:,i), dmat_a, nbf)
    end do

    do i = 1, 6
        qxyz(i) = qxyz(i) - traceprod_sym_packed(mints(:,3+i), dmat_a, nbf)
    end do

    do i = 1, 10
        oxyz(i) = oxyz(i) - traceprod_sym_packed(mints(:,9+i), dmat_a, nbf)
    end do

    if (urohf) then
      do i = 1, 3
        dip(i) = dip(i) - traceprod_sym_packed(mints(:,i), dmat_b, nbf)
      end do
      do i = 1, 6
        qxyz(i) = qxyz(i) - traceprod_sym_packed(mints(:,3+i), dmat_b, nbf)
      end do
    do i = 1, 10
        oxyz(i) = oxyz(i) - traceprod_sym_packed(mints(:,9+i), dmat_b, nbf)
    end do
    end if

    dip = dip*AU_TO_DEBYE

    ! Assemble quadrupole tensor:
    ! Q_ab = 1/2 * ( 3 * r_a * r_b - r^2 * \delta(a,b) )
    quad(X__,X__) = 0.5*(2*qxyz(XX_) - (qxyz(YY_)+qxyz(ZZ_)))
    quad(X__,Y__) = 0.5*(3*qxyz(XY_)                        )
    quad(X__,Z__) = 0.5*(3*qxyz(XZ_)                        )
    quad(Y__,Y__) = 0.5*(2*qxyz(YY_) - (qxyz(XX_)+qxyz(ZZ_)))
    quad(Y__,Z__) = 0.5*(3*qxyz(YZ_)                        )
    quad(Z__,Z__) = 0.5*(2*qxyz(ZZ_) - (qxyz(XX_)+qxyz(YY_)))
    call triangular_to_full(quad,3,'u')
    quad = quad*AU_TO_BUCK

    ! Assemble octopole tensor:
    ! O_abc = 1/2*(5*r_a*r_b*r_c - r^2*(r_a*\delta(b,c)+r_b*\delta(a,c)+r_c*\delta(a,b))
    oct(XXX) = 0.5*(2*oxyz(XXX) - 3*(oxyz(YYX)+oxyz(ZZX)))
    oct(YYY) = 0.5*(2*oxyz(YYY) - 3*(oxyz(XXY)+oxyz(ZZY)))
    oct(ZZZ) = 0.5*(2*oxyz(ZZZ) - 3*(oxyz(YYZ)+oxyz(XXZ)))
    oct(XXY) = 0.5*(4*oxyz(XXY) -   (oxyz(YYY)+oxyz(ZZY)))
    oct(XXZ) = 0.5*(4*oxyz(XXZ) -   (oxyz(YYZ)+oxyz(ZZZ)))
    oct(YYX) = 0.5*(4*oxyz(YYX) -   (oxyz(XXX)+oxyz(ZZX)))
    oct(YYZ) = 0.5*(4*oxyz(YYZ) -   (oxyz(XXZ)+oxyz(ZZZ)))
    oct(ZZX) = 0.5*(4*oxyz(ZZX) -   (oxyz(XXX)+oxyz(YYX)))
    oct(ZZY) = 0.5*(4*oxyz(ZZY) -   (oxyz(XXY)+oxyz(YYY)))
    oct(XYZ) = 0.5*(5*oxyz(XYZ)                          )

    oct = oct*AU_TO_OCT

    write(iw,'(/1x,a)') 'At point (Bohr):'
    write(iw,'(4x,40x,3(a8,7x))') 'X', 'Y', 'Z'
    write(iw,'(4x,40x,3f15.8,sp,f12.5)') com

    write(iw,'(/4x,a)') 'electric charge (a.u.):'
    write(iw,'(4x,34x,a6,sp,f12.5)') 'C_o', real(infos%mol_prop%charge)

    write(iw,'(/4x,a)') 'electric dipole (Debye):'
    write(iw,'(4x,40x,3(a8,7x),a10)') 'X', 'Y', 'Z', 'Norm'
    write(iw,'(4x,34x,a6,4f15.8)') 'D_o', dip, norm2(dip)

    write(iw,'(/4x,a)') 'electric quadrupole (Buckingham):'
    write(iw,'(4x,40x,3(a8,7x))') 'X', 'Y', 'Z'
    write(iw,'(4x,34x,a6,3F15.8)') 'Q_X', quad(:,1)
    write(iw,'(4x,34x,a6,3f15.8)') 'Q_Y', quad(:,2)
    write(iw,'(4x,34x,a6,3f15.8)') 'Q_Z', quad(:,3)

    write(iw,'(/4x,a)') 'electric octoupole (Buckingham*Angstrom):'
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_xxx', oct(XXX)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_xxy', oct(XXY)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_xxz', oct(XXZ)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_xyy', oct(YYX)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_xyz', oct(XYZ)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_xzz', oct(ZZX)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_yyy', oct(YYY)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_yyz', oct(YYZ)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_yyz', oct(ZZY)
    write(iw,'(4x,34x,a6,*(F15.8))') 'O_zzz', oct(ZZZ)

    close(iw)

  end subroutine electric_moments

end module electric_moments_mod
