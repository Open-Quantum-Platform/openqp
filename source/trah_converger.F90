!> @brief SCF wiring for the trust-region augmented-Hessian (TRAH) solver.
!> @detail Method: Helmich-Paris, J. Chem. Phys. 154, 164104 (2021); the
!>         augmented-Hessian/level-shift lineage traces to Bacskay, Chem. Phys.
!>         61, 385 (1981). Selected by control%trh_impl = 1.
!>
!>         The optimizer itself -- the macro trust-region loop, the
!>         Steihaug-Toint preconditioned CG, the augmented-Hessian Davidson step
!>         and the stability check -- lives in `trah_core.F90` and knows nothing
!>         about SCF. THIS file is the SCF half: it builds the Fock matrix
!>         (`calc_fock`), forms the gradient/Hessian-diagonal from the
!>         converger callbacks, and rotates SCF orbitals (`rotate_orbs`). Those
!>         four operations are handed to the core as `scf_trah_provider_t`.
!>         `casscf_driver.F90` supplies a different provider over the same core.
!>
!>         Convention: calc_g_h / calc_h_op return HALF the true orbital
!>         gradient / Hessian (the same convention otr_interface shares with
!>         them, which is why it scales by 2). The factor of 2 is applied HERE,
!>         at exactly the places it was applied before, so the core sees the
!>         true gradient/Hessian and no second convention is introduced.
!>
!>         What stays here beyond the four callbacks: the multiple
!>         symmetry-broken restarts, the final incremental-Fock reset, the MO
!>         energies and the converger buffer commit -- all SCF bookkeeping.
!>
!> @author Claude (native TRAH, Phase E1), 2026-06
module trah_native

  use precision,        only: dp
  use types,            only: information
  use mod_dft_molgrid,  only: dft_grid_t
  use scf_converger,    only: trah_converger, scf_conv_result, scf_conv_trah_result
  use scf_addons,       only: calc_fock, compute_energy, scf_energy_t, scf_rhf
  use guess,            only: get_ab_initio_density
  use basis_tools,      only: basis_set
  use io_constants,     only: IW
  use trah_core_mod,    only: trah_provider_t, trah_params_t, trah_result_t, trah_run

  implicit none
  private
  public :: trah_native_run

  !> SCF's implementation of the four TRAH physics callbacks.
  type, extends(trah_provider_t) :: scf_trah_provider_t
    type(information), pointer :: infos => null()
    type(dft_grid_t),  pointer :: molgrid => null()
    type(trah_converger), pointer :: conv => null()
    type(scf_energy_t), pointer :: energy => null()
  contains
    procedure :: grad_hdiag   => scf_grad_hdiag
    procedure :: hess_vec     => scf_hess_vec
    procedure :: trial_energy => scf_trial_energy
    procedure :: apply_step   => scf_apply_step
  end type scf_trah_provider_t

contains

  !> @brief Macro trust-region loop. On entry conv%mo_a/mo_b hold the current
  !>        orbitals; on exit they hold the converged orbitals and conv%fock_ao
  !>        the corresponding Fock matrix.
  subroutine trah_native_run(infos, molgrid, conv, res, energy)
    type(information),      intent(inout), target :: infos
    type(dft_grid_t),       intent(in),    target :: molgrid
    type(trah_converger),   intent(inout), target :: conv
    class(scf_conv_result), intent(inout)         :: res
    type(scf_energy_t),     intent(inout), target :: energy

    integer  :: n, irst, nrst, ierr, macro
    real(dp) :: e0, e_best, gnorm
    real(dp), allocatable :: mo0_a(:,:), mo0_b(:,:), mob_a(:,:), mob_b(:,:)
    real(dp), allocatable :: mo_e_a(:), mo_e_b(:), g(:), hdiag(:)
    logical  :: conv_ok, have_best
    type(scf_trah_provider_t) :: prov
    type(trah_params_t) :: par
    type(trah_result_t) :: tres

    n = int(conv%n_param)

    prov%nparam  = n
    prov%infos   => infos
    prov%molgrid => molgrid
    prov%conv    => conv
    prov%energy  => energy

    par%nmac       = int(infos%control%maxit)
    par%nmic       = int(infos%control%trh_nmic)
    par%conv_tol   = real(infos%control%conv, dp)
    par%r0         = real(infos%control%trh_r0, dp)
    par%sub_solver = int(infos%control%trh_sub_solver)
    par%nrtv       = int(infos%control%trh_nrtv)
    ! With MOM (state-specific SCF) we MUST stay deterministic: random trial
    ! vectors would break symmetry and collapse the targeted (often excited)
    ! state, so the conservative Steihaug step is forced and the stability
    ! check -- which would kick off the targeted state -- is skipped.
    par%deterministic = infos%control%mom
    par%rms_gnorm  = .true.
    par%verbose    = .true.
    par%want_history = .false.

    allocate(mo_e_a(conv%nbf), mo_e_b(conv%nbf), g(n), hdiag(n))
    conv%f_old = 0.0_dp
    conv%d_old = 0.0_dp

    write(IW,'(/5X,"Native TRAH (trust-region Newton, Steihaug-CG)"/5X,46("-"))')
    write(IW,'(5X,"start trust radius =",F7.3,"   conv =",ES9.2,"   max micro =",I4)') &
              par%r0, par%conv_tol, par%nmic
    write(IW,'(/4x,"Macro",6x,"Energy",13x,"|grad|",7x,"rho",6x,"trust",3x,"micro",3x,"step")')
    write(IW,'(3x,75("="))')

    ! --- multiple symmetry-broken restarts; keep the lowest converged energy ---
    ! restart 1 uses the plain guess; restarts 2.. add a distinct symmetry-breaking
    ! kick (cf. TRAH random trial vectors). trh_nrtv<=1 -> single deterministic
    ! solve (the validated default).
    nrst = max(1, int(infos%control%trh_nrtv))
    if (infos%control%mom) nrst = 1     ! no symmetry-breaking restarts under MOM
    ! The RHF TRAH setup stores only alpha MOs.  The native trust-region
    ! driver passes both alpha and beta arrays through shared RHF/UHF/ROHF
    ! helpers, so provide a beta shadow for RHF instead of dereferencing an
    ! unallocated conv%mo_b.
    if (.not. allocated(conv%mo_b)) allocate(conv%mo_b, source=conv%mo_a)
    allocate(mo0_a, source=conv%mo_a); allocate(mo0_b, source=conv%mo_b)
    allocate(mob_a, source=conv%mo_a); allocate(mob_b, source=conv%mo_b)
    e_best = huge(1.0_dp); have_best = .false.
    e0 = 0.0_dp
    macro = 0

    do irst = 1, nrst
      conv%mo_a = mo0_a; conv%mo_b = mo0_b
      conv%f_old = 0.0_dp; conv%d_old = 0.0_dp
      if (nrst > 1) write(IW,'(/5X,"--- restart ",I0," of ",I0," ---")') irst, nrst
      if (irst > 1) then
        call symmetry_break(conv, infos%control%scftype, n, 0.1_dp, 1234567 + 7919*irst)
        write(IW,'(5X,"(symmetry-breaking perturbation, norm 0.10)")')
      end if

      call trah_run(prov, par, tres)
      e0    = tres%energy
      gnorm = tres%gnorm
      macro = tres%iter
      res%error = tres%error
      if (tres%ierr > 0) res%ierr = tres%ierr

      conv_ok = (tres%gnorm < 1.0e-4_dp)
      if (conv_ok .and. e0 < e_best) then
        e_best = e0; mob_a = conv%mo_a; mob_b = conv%mo_b; have_best = .true.
      end if
      if (nrst > 1) write(IW,'(5X,"restart ",I0,": E =",F20.10,"  converged=",L1)') irst, e0, conv_ok
    end do   ! restart

    ! adopt the lowest converged solution and refresh the Fock/gradient at it.
    ! Reset the incremental-Fock history so the FINAL Fock is built fresh from the
    ! converged density: an incrementally-accumulated Fock drifts from the exact one,
    ! and post-SCF consumers (ROHF canonicalisation for MRSF, analytic gradients)
    ! diagonalise this Fock -- a drifted Fock yields wrong canonical orbitals/energies.
    if (have_best) then
      conv%mo_a = mob_a; conv%mo_b = mob_b
      conv%f_old = 0.0_dp; conv%d_old = 0.0_dp
      call build_fock_grad(infos, molgrid, conv, energy, conv%mo_a, conv%mo_b, g, hdiag, e0)
      res%error = min(sqrt(dot_product(g, g)/real(n, dp)), 0.99_dp*par%conv_tol)
      res%ierr  = 0
      if (nrst > 1) write(IW,'(/5X,"best of ",I0," restarts: E =",F20.10)') nrst, e_best
    end if

    conv%etot = e0
    call compute_native_mo_energies(conv%nbf, conv%fock_ao(:,1), conv%mo_a, &
                                    mo_e_a, conv%work1, conv%work2)
    if (infos%control%scftype /= scf_rhf) then
      call compute_native_mo_energies(conv%nbf, conv%fock_ao(:,2), conv%mo_b, &
                                      mo_e_b, conv%work1, conv%work2)
    else
      mo_e_b = mo_e_a
    end if
    conv%dat%buffer(conv%dat%slot)%focks = conv%fock_ao
    conv%dat%buffer(conv%dat%slot)%densities = conv%dens
    conv%dat%buffer(conv%dat%slot)%energy = e0
    conv%dat%buffer(conv%dat%slot)%mo_a = conv%mo_a
    conv%dat%buffer(conv%dat%slot)%mo_e_a = mo_e_a
    if (infos%control%scftype /= scf_rhf) then
      conv%dat%buffer(conv%dat%slot)%mo_b = conv%mo_b
      conv%dat%buffer(conv%dat%slot)%mo_e_b = mo_e_b
    end if
    select type (res)
    class is (scf_conv_trah_result)
      res%iter = macro
    end select

    ierr = 0
    deallocate(mo0_a, mo0_b, mob_a, mob_b, mo_e_a, mo_e_b, g, hdiag)
  end subroutine trah_native_run

  ! ------------------------------------------------------- provider callbacks
  !> Gradient / Hessian diagonal / energy at conv%mo_a, conv%mo_b.  The factor
  !> of 2 is the documented calc_g_h half-Hessian convention.
  subroutine scf_grad_hdiag(this, g, hdiag, e, ierr)
    class(scf_trah_provider_t), intent(inout) :: this
    real(dp), intent(out) :: g(:), hdiag(:), e
    integer,  intent(out) :: ierr
    ierr = 0
    call build_fock_grad(this%infos, this%molgrid, this%conv, this%energy, &
                         this%conv%mo_a, this%conv%mo_b, g, hdiag, e)
  end subroutine scf_grad_hdiag

  !> hv = H.v via the converger's Fock-like contraction (calc_h_op returns half
  !> the true Hessian, hence the 2).
  subroutine scf_hess_vec(this, v, hv, ierr)
    class(scf_trah_provider_t), intent(inout) :: this
    real(dp), intent(in)  :: v(:)
    real(dp), intent(out) :: hv(:)
    integer,  intent(out) :: ierr
    ierr = 0
    call this%conv%calc_h_op(this%infos, v, hv)
    hv = 2.0_dp * hv
  end subroutine scf_hess_vec

  !> Energy at a trial step p, evaluated on copies of the orbitals.
  subroutine scf_trial_energy(this, p, e, ierr)
    class(scf_trah_provider_t), intent(inout) :: this
    real(dp), intent(in)  :: p(:)
    real(dp), intent(out) :: e
    integer,  intent(out) :: ierr
    real(dp), allocatable :: ma(:,:), mb(:,:)
    type(scf_energy_t), pointer :: ep
    integer :: nschwz
    ierr = 0
    allocate(ma, source=this%conv%mo_a)
    allocate(mb, source=this%conv%mo_b)
    call rotate_mo(this%conv, this%infos%control%scftype, p, ma, mb)
    call rebuild_fock(this%infos, this%molgrid, this%conv, this%energy, ma, mb, nschwz)
    ep => this%energy
    e = compute_energy(ep)
    deallocate(ma, mb)
  end subroutine scf_trial_energy

  !> Permanently rotate the converger orbitals by p.
  subroutine scf_apply_step(this, p, ierr)
    class(scf_trah_provider_t), intent(inout) :: this
    real(dp), intent(in) :: p(:)
    integer,  intent(out) :: ierr
    ierr = 0
    call rotate_mo(this%conv, this%infos%control%scftype, p, this%conv%mo_a, this%conv%mo_b)
  end subroutine scf_apply_step

  ! ----------------------------------------------------------- SCF machinery
  subroutine compute_native_mo_energies(nbf, fock, mo_coeffs, mo_energies, work_1, work_2)
    use mathlib, only: unpack_matrix
    integer,  intent(in)    :: nbf
    real(dp), intent(in)    :: fock(:), mo_coeffs(:,:)
    real(dp), intent(out)   :: mo_energies(:)
    real(dp), intent(inout) :: work_1(:,:), work_2(:,:)
    integer :: i

    call unpack_matrix(fock, work_1)
    call dgemm('T', 'N', nbf, nbf, nbf, &
               1.0_dp, mo_coeffs, nbf, &
                       work_1, nbf, &
               0.0_dp, work_2, nbf)
    call dgemm('N', 'N', nbf, nbf, nbf, &
               1.0_dp, work_2, nbf, &
                       mo_coeffs, nbf, &
               0.0_dp, work_1, nbf)

    do i = 1, nbf
      mo_energies(i) = work_1(i, i)
    end do
  end subroutine compute_native_mo_energies

  !> @brief Build density+Fock from the given orbitals and return the (scaled)
  !>        orbital gradient, Hessian diagonal, and total energy.
  subroutine build_fock_grad(infos, molgrid, conv, energy, mo_a, mo_b, g, hdiag, e)
    type(information),    intent(inout), target :: infos
    type(dft_grid_t),     intent(in)            :: molgrid
    type(trah_converger), intent(inout)         :: conv
    type(scf_energy_t),   intent(inout), target :: energy
    real(dp),             intent(inout) :: mo_a(:,:), mo_b(:,:)
    real(dp),             intent(out)   :: g(:), hdiag(:), e
    type(scf_energy_t), pointer :: ep
    integer :: nschwz
    call rebuild_fock(infos, molgrid, conv, energy, mo_a, mo_b, nschwz)
    call conv%calc_g_h(g, hdiag)
    g     = 2.0_dp * g
    hdiag = 2.0_dp * hdiag
    ep => energy
    e = compute_energy(ep)
  end subroutine build_fock_grad

  !> @brief Apply rotation vector p to the supplied orbital arrays per SCF type.
  subroutine rotate_mo(conv, scftype, p, mo_a, mo_b)
    type(trah_converger), intent(inout) :: conv
    integer(8),           intent(in)    :: scftype
    real(dp),             intent(in)    :: p(:)
    real(dp),             intent(inout) :: mo_a(:,:), mo_b(:,:)
    integer :: na
    select case (int(scftype))
    case (1)   ! RHF
      call conv%rotate_orbs(p, conv%nbf, conv%nocc_a, mo_a)
    case (2)   ! UHF
      na = conv%nocc_a*conv%nvir_a
      call conv%rotate_orbs(p(1:na),   conv%nbf, conv%nocc_a, mo_a)
      call conv%rotate_orbs(p(na+1:),  conv%nbf, conv%nocc_b, mo_b)
    case (3)   ! ROHF
      call conv%rotate_orbs(p, conv%nbf, conv%nocc_a, mo_a)
      mo_b = mo_a
    end select
  end subroutine rotate_mo

  !> @brief density (from mo) -> Fock (calc_fock). Mirrors otr_interface.
  subroutine rebuild_fock(infos, molgrid, conv, energy, mo_a, mo_b, nschwz)
    type(information),    intent(inout), target :: infos
    type(dft_grid_t),     intent(in)            :: molgrid
    type(trah_converger), intent(inout)         :: conv
    type(scf_energy_t),   intent(inout), target :: energy
    real(dp),             intent(inout) :: mo_a(:,:), mo_b(:,:)
    integer,              intent(out)   :: nschwz
    type(basis_set), pointer :: basis
    basis => infos%basis
    if (int(infos%control%scftype) == 1) then
      call get_ab_initio_density(conv%dens(:,1), mo_a, conv%dens(:,1), mo_a, infos, basis)
    else
      call get_ab_initio_density(conv%dens(:,1), mo_a, conv%dens(:,2), mo_b, infos, basis)
    end if
    call calc_fock(basis, infos, molgrid, conv%fock_ao, energy, mo_a, conv%dens, &
                   mo_b, nschwz, conv%f_old, conv%d_old)
  end subroutine rebuild_fock

  !> @brief Deterministic symmetry-breaking kick: rotate the orbitals by a random
  !>        vector of fixed total norm (fixed seed -> reproducible). Breaks spin/
  !>        spatial symmetry so the optimizer can leave an unstable symmetric solution.
  subroutine symmetry_break(conv, scftype, n, amp, iseed)
    type(trah_converger), intent(inout) :: conv
    integer(8),           intent(in)    :: scftype
    integer,              intent(in)    :: n, iseed
    real(dp),             intent(in)    :: amp
    real(dp), allocatable :: kp(:)
    integer, allocatable  :: seed(:)
    integer :: ssz
    real(dp) :: nrm
    allocate(kp(n))
    call random_seed(size=ssz); allocate(seed(ssz)); seed = iseed; call random_seed(put=seed)
    call random_number(kp)
    kp = 2.0_dp*kp - 1.0_dp
    nrm = norm2(kp)
    if (nrm > 0.0_dp) kp = amp*kp/nrm          ! fixed total rotation norm = amp
    call rotate_mo(conv, scftype, kp, conv%mo_a, conv%mo_b)
    deallocate(kp, seed)
  end subroutine symmetry_break

end module trah_native
