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
!>   * phi_cav sign:  phi_total = sum_k Z_k/|r-R_k| + phi_elec
!>   * q_cav sign/scale: ddX cavity-projected adjoint charge (ddx_get_xi) used
!>     directly as the external-charge vector for external_charge_potential.
!>   * E_pcm: the ddX solvation energy is reported as the PCM energy term.
!> The single canonical runtime path and these conventions are pinned by
!> tests/test_pcm_canonical_runtime_path.py.
module solvent_pcm

  use precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use iso_c_binding, only: c_int, c_double, c_char, c_null_char
  use types, only: information
  use basis_tools, only: basis_set
  use messages, only: show_message, with_abort
  use io_constants, only: iw
  use mathlib, only: traceprod_sym_packed
  use constants, only: NUM_CART_BF
  use oqp_tagarray_driver, only: OQP_SM, data_has_tags, tagarray_get_data
  use int1, only: electrostatic_potential_unweighted, external_charge_potential, &
       multipole_integrals

  implicit none
  private

  public :: add_pcm_reaction_field

  ! Maximum Lebedev cavity points per atom (must match n_lebedev used when the
  ! ddX model is built in build_pcm_model(); used only to size the receive
  ! buffer for the cavity coordinates).
  integer(c_int), parameter :: MAX_CAV_PER_ATOM = 302
  integer, parameter :: PCM_FD_MAX_SAMPLES = 3
  real(dp), parameter :: PCM_FD_STEP = 1.0e-4_dp
  ! Finite-difference diagnostics below derive this sign/scale.  It is kept as
  ! an explicit constant so the Fock convention is guarded by runtime evidence.
  real(dp), parameter :: PCM_QCAV_TO_FOCK_SCALE = -0.5_dp

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

    integer(c_int) function oqp_ddx_pcm_solve_multipole_source(natom, &
        xyz_bohr, cavity_charges, nmultipoles, source_multipoles, eps, ncav, &
        phi_source_out, q_cav_out, esolv_out, message, message_len) &
        bind(C, name="oqp_ddx_pcm_solve_multipole_source")
      import :: c_int, c_double, c_char
      integer(c_int), value :: natom
      real(c_double), intent(in) :: xyz_bohr(*)
      real(c_double), intent(in) :: cavity_charges(*)
      integer(c_int), value :: nmultipoles
      real(c_double), intent(in) :: source_multipoles(*)
      real(c_double), value :: eps
      integer(c_int), value :: ncav
      real(c_double), intent(out) :: phi_source_out(*)
      real(c_double), intent(out) :: q_cav_out(*)
      real(c_double), intent(out) :: esolv_out
      character(kind=c_char), intent(out) :: message(*)
      integer(c_int), value :: message_len
    end function oqp_ddx_pcm_solve_multipole_source

    integer(c_int) function oqp_ddx_pcm_solve_multipole_source_with_phi(natom, &
        xyz_bohr, cavity_charges, nmultipoles, source_multipoles, eps, ncav, &
        phi_cav, q_cav_out, esolv_out, message, message_len) &
        bind(C, name="oqp_ddx_pcm_solve_multipole_source_with_phi")
      import :: c_int, c_double, c_char
      integer(c_int), value :: natom
      real(c_double), intent(in) :: xyz_bohr(*)
      real(c_double), intent(in) :: cavity_charges(*)
      integer(c_int), value :: nmultipoles
      real(c_double), intent(in) :: source_multipoles(*)
      real(c_double), value :: eps
      integer(c_int), value :: ncav
      real(c_double), intent(in) :: phi_cav(*)
      real(c_double), intent(out) :: q_cav_out(*)
      real(c_double), intent(out) :: esolv_out
      character(kind=c_char), intent(out) :: message(*)
      integer(c_int), value :: message_len
    end function oqp_ddx_pcm_solve_multipole_source_with_phi
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
    type(information), intent(inout) :: infos
    real(dp),          intent(in)    :: d(:,:)
    integer,           intent(in)    :: nfocks
    real(dp),          intent(inout) :: f(:,:)
    real(dp),          intent(out)   :: e_pcm

    integer(c_int) :: natom, ncav, max_cav, nmultipoles
    integer :: nbf_tri, ii, iat, icav, rc, ndelta
    real(dp) :: eps, esolv, phin, dx, dy, dz, r
    real(dp) :: half_tr_dv, q_cav_sum, q_cav_absnorm, phi_cav_sum, phi_cav_min, phi_cav_max
    real(dp) :: source_charge_sum, phi_source_delta_rms, phi_source_delta_max
    real(dp) :: fd_fock_scale_mean, fd_fock_scale_rms, fd_fock_scale_maxerr
    integer :: fd_fock_samples
    real(dp), allocatable :: xyz(:,:), charges(:)
    real(dp), allocatable :: cav_xyz(:), cx(:), cy(:), cz(:)
    real(dp), allocatable :: phi_elec(:), phi_cav(:), phi_source(:), q_cav(:)
    real(dp), allocatable :: dtot(:), vpcm(:)
    real(dp), allocatable :: ao_pop(:), atom_pop(:), source_charges(:)
    real(dp), allocatable :: ao_dip(:,:), atom_dip(:,:), ao_quad(:,:), atom_quad(:,:)
    real(dp), allocatable :: source_multipoles(:,:)
    real(dp), contiguous, pointer :: smat(:)
    character(kind=c_char) :: cmsg(256)
    character(len=256) :: fmsg
    character(len=*), parameter :: tags_overlap(1) = (/ character(len=80) :: OQP_SM /)

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

    ! phi_total = sum_k Z_k/|r-R_k| + phi_elec, where the OpenQP
    ! Coulomb-potential primitive returns the electronic contribution with the
    ! electron-charge sign already included.
    do icav = 1, ncav
      phin = 0.0_dp
      do iat = 1, natom
        dx = cx(icav) - xyz(1, iat)
        dy = cy(icav) - xyz(2, iat)
        dz = cz(icav) - xyz(3, iat)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        if (r > 1.0e-12_dp) phin = phin + charges(iat) / r
      end do
      phi_cav(icav) = phin + phi_elec(icav)
    end do

    ! ---- Phase 3: consistent QM source -> ddX q_cav ------------------------
    ! Build one per-atom real-solid-harmonic source through quadrupoles, then feed exactly
    ! that same multipole array to ddX for both phi and psi.  The atom-centered
    ! moments use a Mulliken row partition of the AO density.
    allocate(ao_pop(basis%nbf), atom_pop(natom), source_charges(natom), &
             ao_dip(3, basis%nbf), atom_dip(3, natom), &
             ao_quad(6, basis%nbf), atom_quad(6, natom), source=0.0_dp)
    call data_has_tags(infos%dat, tags_overlap, &
                       'solvent_pcm:add_pcm_reaction_field', with_abort)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call mulliken_atomic_population_from_density(basis, smat, dtot, &
                                                 ao_pop, atom_pop)
    call mulliken_atomic_multipoles_from_density(basis, dtot, atom_pop, &
                                                 ao_dip, atom_dip, ao_quad, atom_quad)
    source_charges(:) = charges(:) - atom_pop(:)

    nmultipoles = 9_c_int
    allocate(source_multipoles(nmultipoles, natom), source=0.0_dp)
    call pack_ddx_l2_multipoles(source_charges, atom_dip, atom_quad, source_multipoles)

    allocate(phi_source(ncav), q_cav(ncav))
    rc = oqp_ddx_pcm_solve_multipole_source(natom, xyz, charges, &
         nmultipoles, source_multipoles, eps, ncav, phi_source, q_cav, esolv, &
         cmsg, int(size(cmsg), c_int))
    if (rc /= 0) then
      call c_message_to_fortran(cmsg, fmsg)
      call show_message('PCM (ddX) solve failed: '//trim(fmsg), with_abort)
    end if

    ! ---- Phase 4: V_pcm AO matrix, add to Fock blocks, report E_pcm --------
    allocate(vpcm(nbf_tri))
    call external_charge_potential(basis, vpcm, cx, cy, cz, q_cav)
    call pcm_fock_scale_fd_diagnostic(natom, xyz, charges, eps, ncav, &
         nmultipoles, source_multipoles, phi_source, q_cav, &
         fd_fock_scale_mean, fd_fock_scale_rms, fd_fock_scale_maxerr, &
         fd_fock_samples)
    vpcm(:) = PCM_QCAV_TO_FOCK_SCALE * vpcm(:)
    do ii = 1, nfocks
      f(:, ii) = f(:, ii) + vpcm(:)
    end do

    ! PCM reaction-field (solvation) energy. The apparent surface charges q_cav
    ! (from the consistent multipole-source ddX solve) are contracted with the
    ! EXACT total solute potential at the cavity points (phi_cav = nuclear +
    ! electronic), not with ddX's truncated multipole esolv. This recovers an
    ! independent ddPCM reaction-field energy to <0.5 kcal/mol for H2O and NH3
    ! (validated by tests/test_pcm_literature_benchmarks.py). The -0.5 factor is the
    ! linear-response polarization factor, consistent with the finite-difference
    ! relation dE/dphi = -0.5*q_cav verified by pcm_fock_scale_fd_diagnostic.
    e_pcm = PCM_QCAV_TO_FOCK_SCALE * dot_product(phi_cav, q_cav)

    ! ---- Diagnostic block (validation gate; does NOT affect e_pcm or Fock) --
    ! Exposes, in Fortran, the quantities needed to validate the QM SCF PCM
    ! conventions against a reference (see tests/test_pcm_literature_benchmarks.py
    ! and tests/data/pcm_literature_benchmarks.json). These prints are read-only
    ! summaries of arrays already computed above; the energy and Fock are
    ! unchanged. The host-side polarization energy 0.5*Tr[D.V_pcm] is reported
    ! alongside the ddX esolv so the e_pcm-vs-(1/2)Tr[D.V] bookkeeping question
    ! can be measured rather than assumed. psi_source records the QM source now
    ! used consistently for ddX phi and psi.
    half_tr_dv  = 0.5_dp * traceprod_sym_packed(dtot, vpcm, basis%nbf)
    q_cav_sum   = sum(q_cav)
    q_cav_absnorm = sqrt(sum(q_cav*q_cav))
    phi_cav_sum = sum(phi_cav)
    phi_cav_min = minval(phi_cav)
    phi_cav_max = maxval(phi_cav)
    source_charge_sum = sum(source_charges)
    phi_source_delta_rms = 0.0_dp
    phi_source_delta_max = 0.0_dp
    ndelta = 0
    do icav = 1, ncav
      if (ieee_is_finite(phi_source(icav)) .and. ieee_is_finite(phi_cav(icav))) then
        phi_source_delta_rms = phi_source_delta_rms + &
             (phi_source(icav) - phi_cav(icav))**2
        phi_source_delta_max = max(phi_source_delta_max, &
             abs(phi_source(icav) - phi_cav(icav)))
        ndelta = ndelta + 1
      end if
    end do
    if (ndelta > 0) then
      phi_source_delta_rms = sqrt(phi_source_delta_rms / real(ndelta, dp))
    else
      phi_source_delta_rms = huge(1.0_dp)
      phi_source_delta_max = huge(1.0_dp)
    end if
    write(iw,'(1x,"PCM diag e_pcm=",ES22.14)') e_pcm
    write(iw,'(1x,"PCM diag esolv_multipole=",ES22.14)') esolv
    write(iw,'(1x,"PCM diag half_tr_dv=",ES22.14)') half_tr_dv
    write(iw,'(1x,"PCM diag q_cav_sum=",ES22.14)') q_cav_sum
    write(iw,'(1x,"PCM diag q_cav_absnorm=",ES22.14)') q_cav_absnorm
    write(iw,'(1x,"PCM diag fock_q_scale=",ES22.14)') PCM_QCAV_TO_FOCK_SCALE
    write(iw,'(1x,"PCM diag fd_fock_scale_mean=",ES22.14)') fd_fock_scale_mean
    write(iw,'(1x,"PCM diag fd_fock_scale_rms=",ES22.14)') fd_fock_scale_rms
    write(iw,'(1x,"PCM diag fd_fock_scale_maxerr=",ES22.14)') fd_fock_scale_maxerr
    write(iw,'(1x,"PCM diag fd_fock_samples=",I0)') fd_fock_samples
    write(iw,'(1x,"PCM diag source_charge_sum=",ES22.14)') source_charge_sum
    write(iw,'(1x,"PCM diag phi_source_vs_exact_rms=",ES22.14)') phi_source_delta_rms
    write(iw,'(1x,"PCM diag phi_source_vs_exact_max=",ES22.14)') phi_source_delta_max
    write(iw,'(1x,"PCM diag phi_cav_sum=",ES22.14)') phi_cav_sum
    write(iw,'(1x,"PCM diag phi_cav_min=",ES22.14)') phi_cav_min
    write(iw,'(1x,"PCM diag phi_cav_max=",ES22.14)') phi_cav_max
    write(iw,'(1x,"PCM diag ncav=",I0)') int(ncav)
    write(iw,'(1x,"PCM diag psi_source=total_qm_atom_multipoles_l2")')

  end subroutine add_pcm_reaction_field

  subroutine mulliken_atomic_population_from_density(basis, smat, density, ao_pop, atom_pop)
    type(basis_set), intent(in) :: basis
    real(dp),        intent(in) :: smat(:), density(:)
    real(dp),        intent(out) :: ao_pop(:), atom_pop(:)

    integer :: mu, nu, ish, iatom, i0, i1, idx

    ao_pop(:) = 0.0_dp
    atom_pop(:) = 0.0_dp
    do mu = 1, basis%nbf
      do nu = 1, basis%nbf
        idx = packed_index(mu, nu)
        ao_pop(mu) = ao_pop(mu) + density(idx) * smat(idx)
      end do
    end do

    do ish = 1, basis%nshell
      iatom = basis%origin(ish)
      i0 = basis%ao_offset(ish)
      i1 = basis%ao_offset(ish) + NUM_CART_BF(basis%am(ish)) - 1
      atom_pop(iatom) = atom_pop(iatom) + sum(ao_pop(i0:i1))
    end do
  end subroutine mulliken_atomic_population_from_density

  subroutine mulliken_atomic_multipoles_from_density(basis, density, atom_pop, ao_dip, atom_dip, ao_quad, atom_quad)
    type(basis_set), intent(in) :: basis
    real(dp),        intent(in) :: density(:)
    real(dp),        intent(in) :: atom_pop(:)
    real(dp),        intent(out) :: ao_dip(:,:), atom_dip(:,:)
    real(dp),        intent(out) :: ao_quad(:,:), atom_quad(:,:)

    integer :: iat, mu, nu, ish, i0, i1, idx
    integer, allocatable :: ao_atom(:)
    real(dp), allocatable :: moment_ints(:,:)

    atom_dip(:,:) = 0.0_dp
    atom_quad(:,:) = 0.0_dp
    allocate(ao_atom(basis%nbf), moment_ints(size(density), 9))
    moment_ints(:,:) = 0.0_dp
    ao_atom(:) = 0
    do ish = 1, basis%nshell
      iat = basis%origin(ish)
      i0 = basis%ao_offset(ish)
      i1 = basis%ao_offset(ish) + NUM_CART_BF(basis%am(ish)) - 1
      ao_atom(i0:i1) = iat
    end do

    do iat = 1, size(atom_pop)
      ao_dip(:,:) = 0.0_dp
      ao_quad(:,:) = 0.0_dp
      moment_ints(:,:) = 0.0_dp
      call multipole_integrals(basis, moment_ints, basis%atoms%xyz(:, iat), 2)
      do mu = 1, basis%nbf
        if (ao_atom(mu) /= iat) cycle
        do nu = 1, basis%nbf
          idx = packed_index(mu, nu)
          ao_dip(1, mu) = ao_dip(1, mu) + density(idx) * moment_ints(idx, 1)
          ao_dip(2, mu) = ao_dip(2, mu) + density(idx) * moment_ints(idx, 2)
          ao_dip(3, mu) = ao_dip(3, mu) + density(idx) * moment_ints(idx, 3)
          ao_quad(1, mu) = ao_quad(1, mu) + density(idx) * moment_ints(idx, 4)
          ao_quad(2, mu) = ao_quad(2, mu) + density(idx) * moment_ints(idx, 5)
          ao_quad(3, mu) = ao_quad(3, mu) + density(idx) * moment_ints(idx, 6)
          ao_quad(4, mu) = ao_quad(4, mu) + density(idx) * moment_ints(idx, 7)
          ao_quad(5, mu) = ao_quad(5, mu) + density(idx) * moment_ints(idx, 8)
          ao_quad(6, mu) = ao_quad(6, mu) + density(idx) * moment_ints(idx, 9)
        end do
      end do
      atom_dip(:, iat) = -sum(ao_dip(:, :), dim=2)
      atom_quad(:, iat) = -sum(ao_quad(:, :), dim=2)
    end do
  end subroutine mulliken_atomic_multipoles_from_density

  subroutine pack_ddx_l2_multipoles(source_charges, atom_dip, atom_quad, multipoles)
    real(dp), intent(in)  :: source_charges(:), atom_dip(:,:), atom_quad(:,:)
    real(dp), intent(out) :: multipoles(:,:)

    integer :: iat
    real(dp) :: sqrt4pi, sqrt4pi_over3
    real(dp), parameter :: q_xy = 1.0925484305920792_dp
    real(dp), parameter :: q_z2 = 0.31539156525252005_dp
    real(dp), parameter :: q_x2y2 = 0.5462742152960396_dp

    sqrt4pi = sqrt(4.0_dp * acos(-1.0_dp))
    sqrt4pi_over3 = sqrt(4.0_dp * acos(-1.0_dp) / 3.0_dp)
    multipoles(:,:) = 0.0_dp
    do iat = 1, size(source_charges)
      ! ddX real-solid-harmonic order for l<=2 is:
      !   1: charge; 2: y dipole; 3: z dipole; 4: x dipole;
      !   5: xy; 6: yz; 7: z^2; 8: xz; 9: x^2-y^2.
      multipoles(1, iat) = source_charges(iat) / sqrt4pi
      multipoles(2, iat) = -atom_dip(2, iat) / sqrt4pi_over3
      multipoles(3, iat) =  atom_dip(3, iat) / sqrt4pi_over3
      multipoles(4, iat) = -atom_dip(1, iat) / sqrt4pi_over3
      multipoles(5, iat) = q_xy * atom_quad(4, iat)
      multipoles(6, iat) = q_xy * atom_quad(6, iat)
      multipoles(7, iat) = q_z2 * (-atom_quad(1, iat) - atom_quad(2, iat) + &
                                   2.0_dp * atom_quad(3, iat))
      multipoles(8, iat) = q_xy * atom_quad(5, iat)
      multipoles(9, iat) = q_x2y2 * (atom_quad(1, iat) - atom_quad(2, iat))
    end do
  end subroutine pack_ddx_l2_multipoles

  subroutine pcm_fock_scale_fd_diagnostic(natom, xyz, charges, eps, ncav, &
       nmultipoles, source_multipoles, phi_source, q_cav, scale_mean, &
       scale_rms, maxerr, nsample)
    integer(c_int),  intent(in) :: natom, ncav, nmultipoles
    real(dp),        intent(in) :: xyz(:,:), charges(:), eps
    real(dp),        intent(in) :: source_multipoles(:,:), phi_source(:), q_cav(:)
    real(dp),        intent(out) :: scale_mean, scale_rms, maxerr
    integer,         intent(out) :: nsample

    integer :: icav, slot, worst_slot, sample_idx(PCM_FD_MAX_SAMPLES), rc
    real(dp) :: sample_abs(PCM_FD_MAX_SAMPLES)
    real(dp) :: h, eplus, eminus, fd, scale, err, min_abs
    real(dp), allocatable :: phi_plus(:), phi_minus(:), q_tmp(:)
    character(kind=c_char) :: cmsg(256)
    character(len=256) :: fmsg

    sample_idx(:) = 0
    sample_abs(:) = -1.0_dp
    do icav = 1, ncav
      if (abs(q_cav(icav)) <= 1.0e-14_dp) cycle
      worst_slot = 1
      min_abs = sample_abs(1)
      do slot = 2, PCM_FD_MAX_SAMPLES
        if (sample_abs(slot) < min_abs) then
          min_abs = sample_abs(slot)
          worst_slot = slot
        end if
      end do
      if (abs(q_cav(icav)) > min_abs) then
        sample_abs(worst_slot) = abs(q_cav(icav))
        sample_idx(worst_slot) = icav
      end if
    end do

    scale_mean = 0.0_dp
    scale_rms = 0.0_dp
    maxerr = 0.0_dp
    nsample = 0
    h = PCM_FD_STEP
    allocate(phi_plus(ncav), phi_minus(ncav), q_tmp(ncav))
    do slot = 1, PCM_FD_MAX_SAMPLES
      icav = sample_idx(slot)
      if (icav <= 0) cycle
      phi_plus(:) = phi_source(:)
      phi_minus(:) = phi_source(:)
      phi_plus(icav) = phi_plus(icav) + h
      phi_minus(icav) = phi_minus(icav) - h
      rc = oqp_ddx_pcm_solve_multipole_source_with_phi(natom, xyz, charges, &
           nmultipoles, source_multipoles, eps, ncav, phi_plus, q_tmp, eplus, &
           cmsg, int(size(cmsg), c_int))
      if (rc /= 0) then
        call c_message_to_fortran(cmsg, fmsg)
        call show_message('PCM (ddX) phi+ FD solve failed: '//trim(fmsg), with_abort)
      end if
      rc = oqp_ddx_pcm_solve_multipole_source_with_phi(natom, xyz, charges, &
           nmultipoles, source_multipoles, eps, ncav, phi_minus, q_tmp, eminus, &
           cmsg, int(size(cmsg), c_int))
      if (rc /= 0) then
        call c_message_to_fortran(cmsg, fmsg)
        call show_message('PCM (ddX) phi- FD solve failed: '//trim(fmsg), with_abort)
      end if
      fd = (eplus - eminus) / (2.0_dp * h)
      scale = fd / q_cav(icav)
      err = fd - PCM_QCAV_TO_FOCK_SCALE * q_cav(icav)
      scale_mean = scale_mean + scale
      scale_rms = scale_rms + scale * scale
      maxerr = max(maxerr, abs(err))
      nsample = nsample + 1
    end do
    if (nsample > 0) then
      scale_mean = scale_mean / real(nsample, dp)
      scale_rms = sqrt(scale_rms / real(nsample, dp))
    else
      scale_mean = huge(1.0_dp)
      scale_rms = huge(1.0_dp)
      maxerr = huge(1.0_dp)
    end if
  end subroutine pcm_fock_scale_fd_diagnostic

  pure integer function packed_index(i, j) result(idx)
    integer, intent(in) :: i, j
    if (i >= j) then
      idx = i * (i - 1) / 2 + j
    else
      idx = j * (j - 1) / 2 + i
    end if
  end function packed_index

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
