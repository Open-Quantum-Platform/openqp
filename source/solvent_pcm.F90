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
!> SOURCE CONSISTENCY (ddX forward Phi vs adjoint Psi):
!>   * Phi (forward solve RHS): the EXACT total solute potential phi_cav at the
!>     cavity points, built from the full AO density (electrostatic_potential_
!>     unweighted) plus the analytic nuclear term.
!>   * Psi (adjoint solve source): a FULL-DENSITY source. Atom-centered real
!>     solid-harmonic multipoles M_lm are accumulated for l = 0..PCM_PSI_LMAX
!>     (=8) from the full AO density by numerical quadrature over a dedicated
!>     source-projection molecular grid that reproduces the reference
!>     ddCOSMO/ddPCM density partition: per-atom (PARENT-ATOM) point
!>     assignment with Becke-original (3-iteration) fuzzy-cell weights and
!>     Treutler-Ahlrichs sqrt(R_i/R_j) atomic-size shifting over the Becke
!>     Bragg-Slater table (H = 0.35 A), WITH the literature outside-sphere
!>     leak continuation q*rsph^(2l+1)/r^(l+1) for r>rsph
!>     (build_full_density_multipoles / pcm_grid_update). The moments are in
!>     the exact ddX harmonic convention (solvent_pcm_harmonics), then mapped
!>     to psi by the ddX rule psi(lm,isph)=4*pi/((2l+1) rsph^l) M_lm(isph)
!>     using the production cavity radii (oqp_ddx_pcm_radii). The production
!>     q_cav is the ddX adjoint charge from oqp_ddx_pcm_solve_psi(psi,
!>     phi_cav): both the forward RHS and the adjoint source are full-density
!>     quantities. This is recorded by "PCM diag pcm_source_mode=full_density_
!>     multipoles_lmax8_exact_phi" and "PCM diag psi_source=full_density_grid_
!>     multipoles_lmax8_becke3_treutler_parent_atom_leak".
!>   * NOTE: the per-sphere moments are partition-defined integrals
!>     (Becke-original/Treutler cells), the same convention the reference
!>     ddPCM implementations project on their per-atom Becke grids; the two
!>     codes agree in the fine-grid limit, with only quadrature-mesh
!>     differences remaining (it is NOT claimed to be bit-identical to
!>     PySCF's grid-projected psi).
!>   * The legacy l<=2 atom-centered Mulliken multipole solve (Phi AND Psi from
!>     the l<=2 source) is still run as a DIAGNOSTIC only, to report the
!>     source-vs-exact phi residual and the q_cav shift between the old l<=2 psi
!>     and the new full-density psi (PCM diag q_cav_*_vs_*_rms).
!>
!> VALIDATED SCALAR CONVENTIONS (analytic Born-ion/ddX oracle gate):
!>   * phi_cav sign:  phi_total = sum_k Z_k/|r-R_k| + phi_elec
!>   * q_cav sign/scale: ddX cavity-projected adjoint charge (ddx_get_xi) used
!>     directly as the external-charge vector for external_charge_potential.
!>   * E_pcm: -0.5 * dot_product(phi_cav, q_cav), with NO additional dielectric
!>     factor: ddX folds the full dielectric response into its ddPCM R_eps
!>     operators, so -0.5*<phi_cav, q_cav> = ddx_pcm_energy = the PHYSICAL
!>     solvation free energy. Proven by the Born-ion oracle (point charge q
!>     centered in a single sphere of radius R): -0.5*<phi,q_cav> reproduces
!>     -(1/2)(1-1/eps)*q^2/R to machine precision at eps = 78.3553 and eps = 2.
!>     An extra f(eps) = (eps-1)/eps here (as in PySCF's solvent.ddpcm) would
!>     double-count the dielectric scaling, by -1.3% at eps=78 and -50% at eps=2.
!> The single canonical runtime path and these conventions are pinned by
!> tests/test_pcm_canonical_runtime_path.py.
module solvent_pcm

  use precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use iso_c_binding, only: c_int, c_double, c_char, c_null_char, &
       c_int64_t, c_bool
  use types, only: information
  use basis_tools, only: basis_set
  use messages, only: show_message, with_abort
  use io_constants, only: iw
  use mathlib, only: traceprod_sym_packed
  use oqp_tagarray_driver, only: OQP_SM, data_has_tags, tagarray_get_data
  use int1, only: electrostatic_potential_unweighted, external_charge_potential, &
       multipole_integrals
  use mod_dft_molgrid, only: dft_grid_t
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t, xc_options_t, run_grid_aos
  use mod_dft_partfunc, only: PTYPE_BECKE3
  use dft, only: dft_initialize
  use solvent_pcm_harmonics, only: pcm_ylmscale, accumulate_point_multipole, &
       accumulate_point_multipole_leak

  implicit none
  private

  public :: add_pcm_reaction_field

  ! Maximum angular momentum of the full-density adjoint source Psi. The ddX
  ! model is built with lmax = 8 (solvent_ddx_adapter.c::build_pcm_model), so
  ! nbasis = (PCM_PSI_LMAX+1)^2 = 81 must match ddx_get_n_basis().
  integer, parameter :: PCM_PSI_LMAX = 8

  ! Maximum Lebedev cavity points per atom (must match n_lebedev used when the
  ! ddX model is built in build_pcm_model(); used only to size the receive
  ! buffer for the cavity coordinates).
  integer(c_int), parameter :: MAX_CAV_PER_ATOM = 302
  integer, parameter :: PCM_FD_MAX_SAMPLES = 3
  real(dp), parameter :: PCM_FD_STEP = 1.0e-4_dp
  ! Finite-difference diagnostics below derive this sign/scale.  It is kept as
  ! an explicit constant so the Fock convention is guarded by runtime evidence.
  real(dp), parameter :: PCM_QCAV_TO_FOCK_SCALE = -0.5_dp

  ! Grid consumer for the PCM full-density-Psi production path. It integrates
  ! the electronic density on a dedicated Becke-original/Treutler-shifted
  ! molecular grid, assigns each weighted point to its PARENT atom (the atom
  ! whose atomic grid generated the slice, xce%currAtom), and accumulates the
  ! corresponding negative electronic charge into ddX-convention real-solid-
  ! harmonic multipoles through PCM_PSI_LMAX. Nuclear monopoles are added by
  ! the driver after the grid loop. This reproduces the reference
  ! ddCOSMO/ddPCM source projection (per-atom Becke-partitioned moments) up
  ! to quadrature-mesh differences.
  type, extends(xc_consumer_t) :: pcm_psi_grid_consumer_t
    integer :: lmax = 0
    integer :: nbasis = 0
    integer :: natom = 0
    real(dp), pointer :: xyz(:,:) => null()
    real(dp), allocatable :: vscales(:)
    real(dp), allocatable :: radii(:)
    real(dp), allocatable :: multipoles(:,:,:)
  contains
    procedure :: parallel_start => pcm_grid_parallel_start
    procedure :: parallel_stop  => pcm_grid_parallel_stop
    procedure :: update         => pcm_grid_update
    procedure :: postUpdate     => pcm_grid_post_update
    procedure :: clean          => pcm_grid_clean
  end type pcm_psi_grid_consumer_t

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

    integer(c_int) function oqp_ddx_pcm_radii(natom, charges, radii_bohr_out, &
        message, message_len) bind(C, name="oqp_ddx_pcm_radii")
      import :: c_int, c_double, c_char
      integer(c_int), value :: natom
      real(c_double), intent(in) :: charges(*)
      real(c_double), intent(out) :: radii_bohr_out(*)
      character(kind=c_char), intent(out) :: message(*)
      integer(c_int), value :: message_len
    end function oqp_ddx_pcm_radii

    integer(c_int) function oqp_ddx_pcm_solve_psi(natom, xyz_bohr, charges, &
        eps, ncav, nbasis, psi, phi_cav, q_cav_out, esolv_out, message, &
        message_len) bind(C, name="oqp_ddx_pcm_solve_psi")
      import :: c_int, c_double, c_char
      integer(c_int), value :: natom
      real(c_double), intent(in) :: xyz_bohr(*)
      real(c_double), intent(in) :: charges(*)
      real(c_double), value :: eps
      integer(c_int), value :: ncav
      integer(c_int), value :: nbasis
      real(c_double), intent(in) :: psi(*)
      real(c_double), intent(in) :: phi_cav(*)
      real(c_double), intent(out) :: q_cav_out(*)
      real(c_double), intent(out) :: esolv_out
      character(kind=c_char), intent(out) :: message(*)
      integer(c_int), value :: message_len
    end function oqp_ddx_pcm_solve_psi
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

    integer(c_int) :: natom, ncav, max_cav, nmultipoles, nbasis
    integer :: nbf_tri, ii, iat, icav, rc, ndelta
    real(dp) :: eps, f_epsilon, esolv, esolv_source, esolv_l2, phin, dx, dy, dz, r
    real(dp) :: half_tr_dv, q_cav_sum, q_cav_absnorm, phi_cav_sum, phi_cav_min, phi_cav_max
    real(dp) :: source_charge_sum, phi_source_delta_rms, phi_source_delta_max
    real(dp) :: q_cav_shift_rms, q_cav_full_vs_l2_rms, psi_full_norm, mult_full_norm
    real(dp) :: fd_fock_scale_mean, fd_fock_scale_rms, fd_fock_scale_maxerr
    integer :: fd_fock_samples
    real(dp), allocatable :: xyz(:,:), charges(:)
    real(dp), allocatable :: cav_xyz(:), cx(:), cy(:), cz(:)
    real(dp), allocatable :: phi_elec(:), phi_cav(:), phi_source(:), q_cav(:), q_cav_source(:)
    real(dp), allocatable :: q_cav_l2(:)
    real(dp), allocatable :: dtot(:), vpcm(:)
    real(dp), allocatable :: ao_pop(:), atom_pop(:), source_charges(:)
    real(dp), allocatable :: ao_dip(:,:), atom_dip(:,:), ao_quad(:,:), atom_quad(:,:)
    real(dp), allocatable :: source_multipoles(:,:)
    real(dp), allocatable :: multipoles_full(:,:), psi_full(:,:), radii(:)
    real(dp), contiguous, pointer :: smat(:)
    character(kind=c_char) :: cmsg(256)
    character(len=256) :: fmsg
    character(len=*), parameter :: tags_overlap(1) = (/ character(len=80) :: OQP_SM /)

    e_pcm = 0.0_dp

    natom = int(infos%mol_prop%natom, c_int)
    eps = infos%control%pcm_epsilon
    ! DIAGNOSTIC-ONLY dielectric factor f(eps) = (eps-1)/eps. It is reported in
    ! the diagnostics below but is NOT applied to e_pcm or the Fock operator:
    ! ddX folds the COMPLETE dielectric response into its ddPCM R_eps operators,
    ! so its pcm_energy = 0.5*<xs,psi> (= -0.5*<phi_cav,q_cav> by the adjoint
    ! identity) is already the physical solvation free energy. This is proven by
    ! the Born-ion oracle: -0.5*<phi,q_cav> = -(1/2)(1-1/eps)q^2/R to machine
    ! precision. PySCF's solvent.ddpcm applies an extra 0.5*f_eps*<psi,Xvec>
    ! scaling on top of its R_eps solve, which is why it FAILS the same Born
    ! oracle; it must not be imitated here.
    f_epsilon = (eps - 1.0_dp) / eps
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

    ! ---- Phase 3: QM source -> ddX q_cav -----------------------------------
    ! The PRODUCTION adjoint source Psi is FULL-DENSITY (3c): atom-centered real
    ! solid-harmonic multipoles for l = 0..PCM_PSI_LMAX accumulated from the AO
    ! density by parent-atom Becke-partitioned grid quadrature, in the exact ddX harmonic
    ! convention, mapped to psi by the ddX rule. The l<=2 Mulliken-multipole
    ! source (3a/3b) is retained ONLY as a diagnostic baseline.
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

    allocate(phi_source(ncav), q_cav(ncav), q_cav_source(ncav), q_cav_l2(ncav))

    ! (3a) DIAGNOSTIC -- legacy all-multipole solve: both Phi and Psi from the
    ! l<=2 atom-centered source. Exposes phi_source (source-vs-exact phi
    ! residual) and q_cav_source. Does not feed the Fock matrix or e_pcm.
    rc = oqp_ddx_pcm_solve_multipole_source(natom, xyz, charges, &
         nmultipoles, source_multipoles, eps, ncav, phi_source, q_cav_source, &
         esolv_source, cmsg, int(size(cmsg), c_int))
    if (rc /= 0) then
      call c_message_to_fortran(cmsg, fmsg)
      call show_message('PCM (ddX) diagnostic source solve failed: '//trim(fmsg), with_abort)
    end if

    ! (3b) DIAGNOSTIC -- previous production path: exact phi_cav + l<=2 multipole
    ! Psi. Kept so the q_cav shift from upgrading to full-density Psi can be
    ! measured (q_cav_full_vs_l2_rms below).
    rc = oqp_ddx_pcm_solve_multipole_source_with_phi(natom, xyz, charges, &
         nmultipoles, source_multipoles, eps, ncav, phi_cav, q_cav_l2, esolv_l2, &
         cmsg, int(size(cmsg), c_int))
    if (rc /= 0) then
      call c_message_to_fortran(cmsg, fmsg)
      call show_message('PCM (ddX) l2-psi diagnostic solve failed: '//trim(fmsg), with_abort)
    end if

    ! (3c) PRODUCTION solve -- full-density Psi (l = 0..PCM_PSI_LMAX) from the AO
    ! density + the EXACT total cavity potential phi_cav (Phase 2). Both the
    ! forward RHS and the adjoint source are full-density. q_cav is the
    ! cavity-projected adjoint charge (ddx_get_xi) from this solve and is what
    ! drives the Fock matrix and e_pcm.
    nbasis = int((PCM_PSI_LMAX+1)**2, c_int)
    allocate(multipoles_full(nbasis, natom), psi_full(nbasis, natom), radii(natom))
    ! Query the production ddX cavity radii FIRST: the full-density Psi must use
    ! each sphere's rsph for the outside-sphere "leak" continuation (the QM
    ! density tail beyond the small vdW sphere), exactly as PySCF's
    ! cache_fake_multipoles does with (r_vdw/r)^(2l+1).
    rc = oqp_ddx_pcm_radii(natom, charges, radii, cmsg, int(size(cmsg), c_int))
    if (rc /= 0) then
      call c_message_to_fortran(cmsg, fmsg)
      call show_message('PCM (ddX) radii query failed: '//trim(fmsg), with_abort)
    end if
    call build_full_density_multipoles(basis, infos, dtot, charges, radii, xyz, &
                                       int(natom), PCM_PSI_LMAX, multipoles_full, &
                                       mult_full_norm)
    call multipoles_to_psi(multipoles_full, radii, PCM_PSI_LMAX, psi_full)
    psi_full_norm = sqrt(sum(psi_full*psi_full))

    rc = oqp_ddx_pcm_solve_psi(natom, xyz, charges, eps, ncav, nbasis, &
         psi_full, phi_cav, q_cav, esolv, cmsg, int(size(cmsg), c_int))
    if (rc /= 0) then
      call c_message_to_fortran(cmsg, fmsg)
      call show_message('PCM (ddX) full-density psi solve failed: '//trim(fmsg), with_abort)
    end if

    ! ---- Phase 4: V_pcm AO matrix, add to Fock blocks, report E_pcm --------
    allocate(vpcm(nbf_tri))
    call external_charge_potential(basis, vpcm, cx, cy, cz, q_cav)
    ! FD-validate dE/dphi = -0.5*q_cav about the EXACT phi_cav baseline using the
    ! SAME full-density psi as the production solve, so the perturbed re-solves
    ! match the production forward/adjoint source.
    call pcm_fock_scale_fd_diagnostic(natom, xyz, charges, eps, ncav, nbasis, &
         psi_full, phi_cav, q_cav, &
         fd_fock_scale_mean, fd_fock_scale_rms, fd_fock_scale_maxerr, &
         fd_fock_samples)
    ! Fock reaction-field operator. The variational PCM free energy
    !   E_pcm = -0.5 * <phi_cav(D), q_cav(D)>
    ! is quadratic in D (both phi_cav and q_cav are linear in D for the
    ! symmetric ddPCM response), so its derivative dE/dD = -ext(q_cav):
    ! the explicit factor of 1/2 cancels against the two equal D-dependent terms
    ! (the phi-side and psi-side contributions, equal by the symmetry of the
    ! continuum reaction-field kernel). The full coupling is therefore
    ! 2*PCM_QCAV_TO_FOCK_SCALE = -1, NOT the bare -0.5 explicit-phi factor that
    ! pcm_fock_scale_fd_diagnostic verifies for dE/dphi. No dielectric factor is
    ! applied: q_cav already carries the full eps response (see f_epsilon note).
    vpcm(:) = 2.0_dp * PCM_QCAV_TO_FOCK_SCALE * vpcm(:)
    do ii = 1, nfocks
      f(:, ii) = f(:, ii) + vpcm(:)
    end do

    ! PCM reaction-field (solvation) energy. The apparent surface charges q_cav
    ! (from the full-density-Psi exact-phi ddX solve) are contracted with the
    ! EXACT total solute potential at the cavity points (phi_cav = nuclear +
    ! electronic). The -0.5 factor is the linear-response polarization factor.
    ! No additional dielectric factor: -0.5*<phi_cav,q_cav> equals ddX's
    ! pcm_energy and the physical solvation free energy (Born-ion oracle).
    !
    ! FOCK DERIVATIVE SCOPE: the production SCF uses the full linear-dielectric
    ! coupling V_pcm = -external_charge_potential(q_cav), recorded as
    ! fock_mode=ddpcm_physical_full_variational_coupling. The finite-difference
    ! probe above only verifies the explicit dE/dphi relation (-0.5*q_cav);
    ! future analytic gradients/response work should add a dedicated dPsi/dD
    ! check for the grid-projected full-density source.
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
    ! RMS shift in the surface charge from the full-density exact-phi production
    ! solve relative to the legacy all-multipole (l<=2 Phi and Psi) solve.
    q_cav_shift_rms = sqrt(sum((q_cav - q_cav_source)**2) / real(ncav, dp))
    ! RMS shift from upgrading the adjoint source from l<=2 multipole Psi to the
    ! full-density Psi, both with the EXACT phi_cav: the quantitative measure of
    ! what the full-density Psi buys over the previous production path.
    q_cav_full_vs_l2_rms = sqrt(sum((q_cav - q_cav_l2)**2) / real(ncav, dp))
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
    write(iw,'(1x,"PCM diag esolv_full_density_psi=",ES22.14)') esolv
    write(iw,'(1x,"PCM diag esolv_l2_psi_exact_phi=",ES22.14)') esolv_l2
    write(iw,'(1x,"PCM diag esolv_source_multipole=",ES22.14)') esolv_source
    write(iw,'(1x,"PCM diag half_tr_dv=",ES22.14)') half_tr_dv
    write(iw,'(1x,"PCM diag q_cav_sum=",ES22.14)') q_cav_sum
    write(iw,'(1x,"PCM diag q_cav_absnorm=",ES22.14)') q_cav_absnorm
    write(iw,'(1x,"PCM diag fock_q_scale=",ES22.14)') PCM_QCAV_TO_FOCK_SCALE
    write(iw,'(1x,"PCM diag f_epsilon=",ES22.14)') f_epsilon
    write(iw,'(1x,"PCM diag fock_q_coupling=",ES22.14)') &
         2.0_dp * PCM_QCAV_TO_FOCK_SCALE
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
    write(iw,'(1x,"PCM diag q_cav_source_vs_exact_rms=",ES22.14)') q_cav_shift_rms
    write(iw,'(1x,"PCM diag q_cav_full_vs_l2_rms=",ES22.14)') q_cav_full_vs_l2_rms
    write(iw,'(1x,"PCM diag multipoles_full_norm=",ES22.14)') mult_full_norm
    write(iw,'(1x,"PCM diag psi_full_norm=",ES22.14)') psi_full_norm
    block
      integer :: ldiag, mdiag, lmdiag
      real(dp) :: lnorm
      do ldiag = 0, PCM_PSI_LMAX
        lnorm = 0.0_dp
        do mdiag = -ldiag, ldiag
          lmdiag = ldiag*ldiag + ldiag + 1 + mdiag
          lnorm = lnorm + sum(multipoles_full(lmdiag,:)**2)
        end do
        write(iw,'(1x,"PCM diag mult_l",I0,"_norm=",ES22.14)') ldiag, sqrt(lnorm)
      end do
      write(iw,'(1x,"PCM diag atom_q0=",10ES16.8)') &
           multipoles_full(1,:) * sqrt(4.0_dp*acos(-1.0_dp))
    end block
    write(iw,'(1x,"PCM diag pcm_source_mode=full_density_multipoles_lmax8_exact_phi")')
    write(iw,'(1x,"PCM diag psi_source=full_density_grid_multipoles_lmax8_becke3_treutler_parent_atom_leak")')
    write(iw,'(1x,"PCM diag fock_mode=ddpcm_physical_full_variational_coupling")')

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
      i1 = basis%ao_offset(ish) + basis%naos(ish) - 1
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
      i1 = basis%ao_offset(ish) + basis%naos(ish) - 1
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


  subroutine build_full_density_multipoles(basis, infos, density_packed, charges, &
       radii, xyz, natom, lmax, multipoles, mult_norm)
    type(basis_set), intent(in) :: basis
    type(information), intent(inout) :: infos
    real(dp), intent(in) :: density_packed(:)
    real(dp), intent(in), target :: charges(:), xyz(:,:)
    real(dp), intent(in) :: radii(:)
    integer, intent(in) :: natom, lmax
    real(dp), intent(out) :: multipoles(:,:)
    real(dp), intent(out) :: mult_norm

    type(basis_set) :: basis_work
    type(dft_grid_t), target :: molgrid
    type(xc_options_t) :: xc_opts
    type(pcm_psi_grid_consumer_t) :: dat
    real(dp), allocatable, target :: density_full(:,:)
    integer :: nbf, i, j, idx, maxl, nang
    real(dp) :: sqrt4pi
    character(kind=c_char) :: saved_xcname(20)
    integer(c_int64_t) :: saved_partfun, saved_bfc_algo, saved_nrad, &
         saved_nang, saved_radtype
    logical(c_bool) :: saved_pruned

    nbf = basis%nbf
    if (size(multipoles,1) /= (lmax+1)**2 .or. size(multipoles,2) /= natom) then
      call show_message('PCM full-density Psi multipole buffer has wrong shape', with_abort)
    end if

    allocate(density_full(nbf, nbf), source=0.0_dp)
    do i = 1, nbf
      do j = 1, nbf
        idx = packed_index(i, j)
        density_full(i,j) = density_packed(idx) * basis%bfnrm(i) * basis%bfnrm(j)
      end do
    end do

    basis_work = basis
    ! HF/reference-SCF inputs can leave infos%dft%xc_functional_name unset/garbage.
    ! dft_set_options still inspects the C-string even when need_functional=.false.,
    ! so follow the existing SAP-grid pattern: temporarily blank the name while
    ! constructing the quadrature grid, then restore it.
    saved_xcname = infos%dft%xc_functional_name
    infos%dft%xc_functional_name = c_null_char
    ! The PCM source-projection grid is PINNED to the reference ddCOSMO/ddPCM
    ! density-partition convention, independent of any user XC-grid settings:
    !   * Becke's ORIGINAL fuzzy-cell partition (3 softening iterations,
    !     JCP 88, 2547 (1988)) -- PTYPE_BECKE3,
    !   * Treutler-Ahlrichs atomic-size surface shifting chi = sqrt(R_i/R_j)
    !     (JCP 102, 346 (1995)) over the Becke Bragg-Slater table (H = 0.35 A)
    !     -- dft_bfc_algo = 2,
    !   * unpruned 240x302 atomic grids on the standard (MHL) radial map.
    ! Together with the parent-atom point assignment in pcm_grid_update this
    ! makes the per-sphere source moments converge to the SAME partitioned
    ! integrals as the reference ddPCM implementations (e.g. PySCF's
    ! ddcosmo.make_psi_vmat on its Becke-partitioned atomic grids); the grid
    ! mesh itself only controls quadrature accuracy, not the partition limit.
    saved_partfun = infos%dft%dft_partfun
    saved_bfc_algo = infos%dft%dft_bfc_algo
    saved_nrad = infos%dft%grid_rad_size
    saved_nang = infos%dft%grid_ang_size
    saved_radtype = infos%dft%rad_grid_type
    saved_pruned = infos%dft%grid_pruned
    infos%dft%dft_partfun = int(PTYPE_BECKE3, c_int64_t)
    infos%dft%dft_bfc_algo = 2_c_int64_t
    infos%dft%grid_rad_size = 240_c_int64_t
    infos%dft%grid_ang_size = 302_c_int64_t
    infos%dft%rad_grid_type = 0_c_int64_t
    infos%dft%grid_pruned = .false._c_bool
    call dft_initialize(infos, basis_work, molgrid, verbose=.false., need_functional=.false.)
    infos%dft%dft_partfun = saved_partfun
    infos%dft%dft_bfc_algo = saved_bfc_algo
    infos%dft%grid_rad_size = saved_nrad
    infos%dft%grid_ang_size = saved_nang
    infos%dft%rad_grid_type = saved_radtype
    infos%dft%grid_pruned = saved_pruned
    infos%dft%xc_functional_name = saved_xcname

    dat%lmax = lmax
    dat%nbasis = (lmax+1)**2
    dat%natom = natom
    dat%xyz => xyz
    allocate(dat%vscales(dat%nbasis))
    call pcm_ylmscale(lmax, dat%vscales)
    allocate(dat%radii(natom))
    dat%radii(:) = radii(1:natom)
    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    maxl = maxval(basis_work%am)
    nang = maxl + 2
    xc_opts%isGGA = .false.
    xc_opts%needTau = .false.
    xc_opts%hasBeta = .false.
    xc_opts%isWFVecs = .false.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molgrid%maxSlicePts
    xc_opts%limPts = molgrid%maxNRadTimesNAng
    xc_opts%numAtoms = natom
    xc_opts%maxAngMom = nang
    xc_opts%nDer = 0
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => density_full
    xc_opts%molGrid => molgrid
    ! No density-based weight screening for the source projection: the Becke
    ! partition has small but nonzero tail weights on every atom's grid, and
    ! the per-sphere moments should integrate them like the reference does.
    xc_opts%dft_threshold = 0.0_dp
    xc_opts%ao_threshold = infos%dft%grid_ao_threshold
    ! Keep AO pruning enabled using the normal DFT threshold; the consumer uses
    ! xce%compMOs/compRho so both pruned and unpruned paths are handled by the engine.
    xc_opts%ao_sparsity_ratio = infos%dft%grid_ao_sparsity_ratio
    if (infos%dft%grid_pruned) xc_opts%ao_sparsity_ratio = 0.0_dp

    call run_grid_aos(xc_opts, dat, basis_work)

    multipoles(:,:) = dat%multipoles(:,:,1)
    ! Add nuclear point charges exactly at their own atom centers. For c=0 only
    ! the l=0 real harmonic contributes: q * Y_00 = q/sqrt(4*pi).
    sqrt4pi = sqrt(4.0_dp * acos(-1.0_dp))
    do i = 1, natom
      multipoles(1,i) = multipoles(1,i) + charges(i) / sqrt4pi
    end do
    mult_norm = sqrt(sum(multipoles*multipoles))

    call dat%clean()
  end subroutine build_full_density_multipoles

  subroutine multipoles_to_psi(multipoles, radii, lmax, psi)
    real(dp), intent(in) :: multipoles(:,:), radii(:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: psi(:,:)

    integer :: iat, l, m, lm
    real(dp) :: pi4, denom

    pi4 = 4.0_dp * acos(-1.0_dp)
    psi(:,:) = 0.0_dp
    do iat = 1, size(multipoles,2)
      do l = 0, lmax
        denom = real(2*l + 1, dp) * radii(iat)**l
        do m = -l, l
          lm = l*l + l + 1 + m
          psi(lm,iat) = pi4 * multipoles(lm,iat) / denom
        end do
      end do
    end do
  end subroutine multipoles_to_psi

  subroutine pcm_grid_parallel_start(self, xce, nthreads)
    class(pcm_psi_grid_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    if (allocated(self%multipoles)) deallocate(self%multipoles)
    allocate(self%multipoles(self%nbasis, self%natom, nthreads), source=0.0_dp)
  end subroutine pcm_grid_parallel_start

  subroutine pcm_grid_parallel_stop(self)
    class(pcm_psi_grid_consumer_t), intent(inout) :: self
    if (allocated(self%multipoles)) then
      if (size(self%multipoles,3) > 1) then
        self%multipoles(:,:,1) = sum(self%multipoles, dim=3)
      end if
      call self%pe%allreduce(self%multipoles(:,:,1), size(self%multipoles(:,:,1)))
    end if
  end subroutine pcm_grid_parallel_stop

  subroutine pcm_grid_update(self, xce, mythread)
    class(pcm_psi_grid_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(dp), allocatable :: rho(:,:)
    real(dp) :: qpt, c(3)
    integer :: ipt, iown

    ! PARENT-ATOM partition, exactly the reference ddCOSMO/ddPCM source
    ! projection (PySCF ddcosmo.make_psi_vmat): every point of the current
    ! slice belongs to the atom whose atomic grid generated it
    ! (xce%currAtom); the fuzzy-cell share of the molecular density at that
    ! point is already carried by the Becke-original/Treutler-shifted
    ! partition weight inside xce%wts (see build_full_density_multipoles).
    ! The multipole about the owning sphere uses the outside-sphere leak
    ! continuation (q*rsph^(2l+1)/r^(l+1) for r>rsph), exactly as PySCF's
    ! cache_fake_multipoles caps points beyond the vdW sphere, which keeps
    ! the QM density tail from blowing up the bare interior moment q*r^l.
    iown = xce%currAtom
    if (iown < 1 .or. iown > self%natom) then
      call show_message('PCM psi grid consumer: slice parent atom not set', &
                        with_abort)
    end if

    call xce%compMOs()
    allocate(rho(2, xce%numPts), source=0.0_dp)
    call xce%compRho(rho)

    do ipt = 1, xce%numPts
      qpt = -sum(rho(:,ipt)) * xce%wts(ipt)
      if (qpt == 0.0_dp) cycle
      c(:) = xce%xyzw(ipt,1:3) - self%xyz(:,iown)
      call accumulate_point_multipole_leak(c, qpt, self%lmax, self%vscales, &
                                           self%radii(iown), &
                                           self%multipoles(:,iown,mythread))
    end do
  end subroutine pcm_grid_update

  subroutine pcm_grid_post_update(self, xce, mythread)
    class(pcm_psi_grid_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    ! No per-slice postprocessing needed; update accumulates directly.
  end subroutine pcm_grid_post_update

  subroutine pcm_grid_clean(self)
    class(pcm_psi_grid_consumer_t), intent(inout) :: self
    if (allocated(self%vscales)) deallocate(self%vscales)
    if (allocated(self%radii)) deallocate(self%radii)
    if (allocated(self%multipoles)) deallocate(self%multipoles)
    nullify(self%xyz)
  end subroutine pcm_grid_clean

  subroutine pcm_fock_scale_fd_diagnostic(natom, xyz, charges, eps, ncav, &
       nbasis, psi, phi_cav, q_cav, scale_mean, scale_rms, maxerr, nsample)
    integer(c_int),  intent(in) :: natom, ncav, nbasis
    real(dp),        intent(in) :: xyz(:,:), charges(:), eps
    real(dp),        intent(in) :: psi(:,:), phi_cav(:), q_cav(:)
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
      phi_plus(:) = phi_cav(:)
      phi_minus(:) = phi_cav(:)
      phi_plus(icav) = phi_plus(icav) + h
      phi_minus(icav) = phi_minus(icav) - h
      rc = oqp_ddx_pcm_solve_psi(natom, xyz, charges, eps, ncav, nbasis, &
           psi, phi_plus, q_tmp, eplus, cmsg, int(size(cmsg), c_int))
      if (rc /= 0) then
        call c_message_to_fortran(cmsg, fmsg)
        call show_message('PCM (ddX) full-psi phi+ FD solve failed: '//trim(fmsg), with_abort)
      end if
      rc = oqp_ddx_pcm_solve_psi(natom, xyz, charges, eps, ncav, nbasis, &
           psi, phi_minus, q_tmp, eminus, cmsg, int(size(cmsg), c_int))
      if (rc /= 0) then
        call c_message_to_fortran(cmsg, fmsg)
        call show_message('PCM (ddX) full-psi phi- FD solve failed: '//trim(fmsg), with_abort)
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
