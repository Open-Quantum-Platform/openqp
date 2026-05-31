!> @brief OpenQP <-> ddX PCM reaction-field bridge for the SCF energy path.
!>
!> @details This module is the Fortran half of the energy-only PCM seam. It
!> declares the iso_c_binding interfaces to the two production C adapter
!> entry points (source/solvent_ddx_adapter.c) and orchestrates the closed
!> reaction-field loop inside a single SCF Fock build:
!>
!>     D  ->  phi_cav  ->  ddX q_cav  ->  V_pcm  ->  Fock / E_pcm
!>
!> The C adapter returns status 2 when OpenQP was built without OQP_ENABLE_DDX,
!> so a PCM-enabled run on a non-ddX build aborts here with a clear message
!> rather than silently producing a vacuum result.
!>
!> PROVISIONAL CONVENTIONS (not yet validated against a trusted PCM reference):
!>   * phi_cav sign:  phi_total = sum_k Z_k/|r-R_k| - phi_elec
!>   * q_cav sign/scale: ddX cavity-projected adjoint charge (ddx_get_xi) used
!>     directly as the external-charge vector for external_charge_potential.
!>   * E_pcm: the ddX solvation energy is reported as the PCM energy term.
!> See docs/solvent_ddx_scf_integration_seam.md.
module solvent_pcm

  use precision, only: dp
  use iso_c_binding, only: c_int, c_double, c_char, c_null_char
  use types, only: information
  use basis_tools, only: basis_set
  use messages, only: show_message, with_abort
  use int1, only: electrostatic_potential_unweighted, external_charge_potential

  implicit none
  private

  public :: add_pcm_reaction_field

  ! Maximum Lebedev cavity points per atom (must match n_lebedev used when the
  ! ddX model is built in build_pcm_model(); used only to size the receive
  ! buffer for the cavity coordinates).
  integer(c_int), parameter :: MAX_CAV_PER_ATOM = 302

  interface
    integer(c_int) function oqp_ddx_pcm_cavity(natom, xyz_bohr, charges, &
        eps, max_cav, ncav_out, cav_xyz_out, message, message_len) &
        bind(C, name="oqp_ddx_pcm_cavity")
      import :: c_int, c_double, c_char
      integer(c_int), value :: natom
      real(c_double), intent(in) :: xyz_bohr(*)
      real(c_double), intent(in) :: charges(*)
      real(c_double), value :: eps
      integer(c_int), value :: max_cav
      integer(c_int), intent(out) :: ncav_out
      real(c_double), intent(out) :: cav_xyz_out(*)
      character(kind=c_char), intent(out) :: message(*)
      integer(c_int), value :: message_len
    end function oqp_ddx_pcm_cavity

    integer(c_int) function oqp_ddx_pcm_solve(natom, xyz_bohr, charges, &
        eps, ncav, phi_cav, q_cav_out, esolv_out, message, message_len) &
        bind(C, name="oqp_ddx_pcm_solve")
      import :: c_int, c_double, c_char
      integer(c_int), value :: natom
      real(c_double), intent(in) :: xyz_bohr(*)
      real(c_double), intent(in) :: charges(*)
      real(c_double), value :: eps
      integer(c_int), value :: ncav
      real(c_double), intent(in) :: phi_cav(*)
      real(c_double), intent(out) :: q_cav_out(*)
      real(c_double), intent(out) :: esolv_out
      character(kind=c_char), intent(out) :: message(*)
      integer(c_int), value :: message_len
    end function oqp_ddx_pcm_solve
  end interface

contains

  !> @brief Add the ddX PCM reaction-field operator to the Fock matrices and
  !>        return the (provisional) PCM energy contribution.
  !> @param[in]    basis   AO basis (read-only)
  !> @param[in]    infos   run information; uses control%pcm_epsilon, atoms, natom
  !> @param[in]    d       packed AO density blocks (nbf_tri, nfocks)
  !> @param[in]    nfocks  number of spin blocks
  !> @param[inout] f       packed AO Fock blocks (nbf_tri, nfocks); V_pcm added
  !> @param[out]   e_pcm   PCM energy contribution (provisional)
  subroutine add_pcm_reaction_field(basis, infos, d, nfocks, f, e_pcm)
    type(basis_set),   intent(in)    :: basis
    type(information), intent(in)    :: infos
    real(dp),          intent(in)    :: d(:,:)
    integer,           intent(in)    :: nfocks
    real(dp),          intent(inout) :: f(:,:)
    real(dp),          intent(out)   :: e_pcm

    integer(c_int) :: natom, ncav, max_cav
    integer :: nbf_tri, ii, iat, icav, rc
    real(dp) :: eps, esolv, phin, dx, dy, dz, r
    real(dp), allocatable :: xyz(:,:), charges(:)
    real(dp), allocatable :: cav_xyz(:), cx(:), cy(:), cz(:)
    real(dp), allocatable :: phi_elec(:), phi_cav(:), q_cav(:)
    real(dp), allocatable :: dtot(:), vpcm(:)
    character(kind=c_char) :: cmsg(256)
    character(len=256) :: fmsg

    e_pcm = 0.0_dp

    natom = int(infos%mol_prop%natom, c_int)
    eps = infos%control%pcm_epsilon
    nbf_tri = size(d, 1)

    allocate(xyz(3, natom), charges(natom))
    xyz(:,:) = infos%atoms%xyz(:, 1:natom)
    charges(:) = infos%atoms%zn(1:natom)

    ! ---- Phase 1: build ddX cavity, retrieve cavity-point coordinates ------
    max_cav = MAX_CAV_PER_ATOM * natom
    allocate(cav_xyz(3 * max_cav))
    ncav = 0
    rc = oqp_ddx_pcm_cavity(natom, xyz, charges, eps, max_cav, ncav, &
                            cav_xyz, cmsg, int(size(cmsg), c_int))
    if (rc /= 0) then
      call c_message_to_fortran(cmsg, fmsg)
      call show_message('PCM (ddX) cavity build failed: '//trim(fmsg), with_abort)
    end if

    allocate(cx(ncav), cy(ncav), cz(ncav))
    do icav = 1, ncav
      cx(icav) = cav_xyz(3*(icav-1) + 1)
      cy(icav) = cav_xyz(3*(icav-1) + 2)
      cz(icav) = cav_xyz(3*(icav-1) + 3)
    end do

    ! ---- Phase 2: total solute electrostatic potential at cavity points ----
    ! Electronic part from the AO density on a temporary copy of the total
    ! density, so the SCF density blocks are not disturbed.
    allocate(dtot(nbf_tri), source=0.0_dp)
    do ii = 1, nfocks
      dtot(:) = dtot(:) + d(:, ii)
    end do

    allocate(phi_elec(ncav), phi_cav(ncav))
    call electrostatic_potential_unweighted(basis, cx, cy, cz, dtot, phi_elec)

    ! phi_total = sum_k Z_k/|r-R_k| - phi_elec   (PROVISIONAL sign convention)
    do icav = 1, ncav
      phin = 0.0_dp
      do iat = 1, natom
        dx = cx(icav) - xyz(1, iat)
        dy = cy(icav) - xyz(2, iat)
        dz = cz(icav) - xyz(3, iat)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        if (r > 1.0e-12_dp) phin = phin + charges(iat) / r
      end do
      phi_cav(icav) = phin - phi_elec(icav)
    end do

    ! ---- Phase 3: ddX solve -> cavity-projected adjoint charge q_cav -------
    allocate(q_cav(ncav))
    rc = oqp_ddx_pcm_solve(natom, xyz, charges, eps, ncav, phi_cav, q_cav, &
                           esolv, cmsg, int(size(cmsg), c_int))
    if (rc /= 0) then
      call c_message_to_fortran(cmsg, fmsg)
      call show_message('PCM (ddX) solve failed: '//trim(fmsg), with_abort)
    end if

    ! ---- Phase 4: V_pcm AO matrix, add to Fock blocks, report E_pcm --------
    allocate(vpcm(nbf_tri))
    call external_charge_potential(basis, vpcm, cx, cy, cz, q_cav)
    do ii = 1, nfocks
      f(:, ii) = f(:, ii) + vpcm(:)
    end do

    ! Provisional PCM energy term reported to the SCF total energy.
    e_pcm = esolv

  end subroutine add_pcm_reaction_field

  !> @brief Copy a NUL-terminated C character buffer into a Fortran string.
  subroutine c_message_to_fortran(cmsg, fmsg)
    character(kind=c_char), intent(in)  :: cmsg(:)
    character(len=*),       intent(out) :: fmsg
    integer :: i
    fmsg = ''
    do i = 1, min(size(cmsg), len(fmsg))
      if (cmsg(i) == c_null_char) exit
      fmsg(i:i) = cmsg(i)
    end do
  end subroutine c_message_to_fortran

end module solvent_pcm
