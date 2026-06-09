!> @brief Trust-region augmented-Hessian (TRAH) SCF solver.
!> @detail Method: Helmich-Paris, J. Chem. Phys. 154, 164104 (2021); the
!>         augmented-Hessian/level-shift lineage traces to Bacskay, Chem. Phys.
!>         61, 385 (1981). Selected by control%trh_impl = 0. Reuses the physics callbacks
!>         already provided by the trah_converger (calc_g_h, calc_h_op, rotate_orbs)
!>         plus calc_fock/get_ab_initio_density/compute_energy -- only the
!>         optimization shell is implemented here.
!>
!>         Algorithm: trust-region Newton. Each macro step solves the
!>         subproblem  min_p  g.p + 1/2 p.H p  s.t. |p| <= Delta  by the
!>         Steihaug-Toint preconditioned conjugate-gradient method (preconditioner
!>         M = diag(h_diag)); the step is accepted/rejected and Delta updated from
!>         the ratio of actual to predicted energy reduction. H.p is formed matrix-
!>         free via calc_h_op (a Fock-like contraction); H is never built.
!>
!>         Convention: calc_g_h / calc_h_op return half the true orbital gradient /
!>         Hessian (same as otr_interface, which scales by 2), so we scale by 2.
!>
!> @note  E1 implementation (CG-Steihaug + basic trust control). Hardening of the
!>        micro-solver for pathological/negative-gap cases (Jacobi-Davidson,
!>        restarts, random trial vectors) is Phase-E2.
!> @author Claude (native TRAH, Phase E1), 2026-06
module trah_native

  use precision,        only: dp
  use types,            only: information
  use mod_dft_molgrid,  only: dft_grid_t
  use scf_converger,    only: trah_converger, scf_conv_result, scf_conv_trah_result
  use scf_addons,       only: calc_fock, compute_energy, scf_energy_t
  use guess,            only: get_ab_initio_density
  use basis_tools,      only: basis_set
  use io_constants,     only: IW

  implicit none
  private
  public :: trah_native_run

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

    integer  :: n, macro, nmac, nmic, micro_used, irst, nrst
    real(dp) :: delta, dmax, conv_tol, gnorm, e0, etrial, rho, pred, snorm, e_best
    real(dp), allocatable :: g(:), hdiag(:), p(:), mo0_a(:,:), mo0_b(:,:), mob_a(:,:), mob_b(:,:)
    real(dp), allocatable :: vmin(:)
    real(dp) :: lam
    logical  :: accepted, conv_ok, have_best
    real(dp), parameter :: stab_eig_tol = 1.0e-4_dp   ! Hessian eigenvalue below -this = unstable
    ! near-convergence guards: once the model can no longer predict a meaningful
    ! energy reduction (pred below FP noise) or the trust radius collapses while
    ! the gradient is already small, the energy is converged even if |g| has not
    ! reached the (tight) gradient tolerance. A trust collapse with a large |g| is
    ! instead a genuine stall and is reported as non-convergence.
    real(dp), parameter :: pred_floor = 1.0e-11_dp
    real(dp), parameter :: delta_min  = 1.0e-4_dp
    real(dp), parameter :: gtol_fp    = 1.0e-4_dp
    real(dp), parameter :: stab_step  = 1.0e-3_dp   ! step above this at small |g| = saddle escape

    n        = int(conv%n_param)
    nmac     = int(infos%control%maxit)
    nmic     = int(infos%control%trh_nmic)
    conv_tol = real(infos%control%conv, dp)
    delta    = real(infos%control%trh_r0, dp)
    dmax     = max(4.0_dp, 8.0_dp*delta)

    allocate(g(n), hdiag(n), p(n), vmin(n))
    conv%f_old = 0.0_dp
    conv%d_old = 0.0_dp

    write(IW,'(/5X,"Native TRAH (trust-region Newton, Steihaug-CG)"/5X,46("-"))')
    write(IW,'(5X,"start trust radius =",F7.3,"   conv =",ES9.2,"   max micro =",I4)') &
              delta, conv_tol, nmic
    write(IW,'(/4x,"Macro",6x,"Energy",13x,"|grad|",7x,"rho",6x,"trust",3x,"micro",3x,"step")')
    write(IW,'(3x,75("="))')

    ! --- multiple symmetry-broken restarts; keep the lowest converged energy ---
    ! restart 1 uses the plain guess; restarts 2.. add a distinct symmetry-breaking
    ! kick (cf. TRAH random trial vectors). trh_nrtv<=1 -> single deterministic
    ! solve (the validated default).
    nrst = max(1, int(infos%control%trh_nrtv))
    if (infos%control%mom) nrst = 1     ! no symmetry-breaking restarts under MOM
    allocate(mo0_a, source=conv%mo_a); allocate(mo0_b, source=conv%mo_b)
    allocate(mob_a, source=conv%mo_a); allocate(mob_b, source=conv%mo_b)
    e_best = huge(1.0_dp); have_best = .false.

    do irst = 1, nrst
      conv%mo_a = mo0_a; conv%mo_b = mo0_b
      conv%f_old = 0.0_dp; conv%d_old = 0.0_dp
      delta = real(infos%control%trh_r0, dp)
      if (nrst > 1) write(IW,'(/5X,"--- restart ",I0," of ",I0," ---")') irst, nrst
      if (irst > 1) then
        call symmetry_break(conv, infos%control%scftype, n, 0.1_dp, 1234567 + 7919*irst)
        write(IW,'(5X,"(symmetry-breaking perturbation, norm 0.10)")')
      end if

      ! initial point at the (possibly perturbed) orbitals
      call build_fock_grad(infos, molgrid, conv, energy, conv%mo_a, conv%mo_b, g, hdiag, e0)

    do macro = 1, nmac
      gnorm = sqrt(dot_product(g, g) / real(n, dp))

      ! trust-region subproblem  ->  step p, predicted reduction pred.
      ! Default (davidson/jacobi_davidson) = augmented-Hessian + random trial vectors
      ! (robust; finds broken-symmetry minima, like OpenTrustRegion). trh_sub_solver=tcg
      ! selects the cheaper Steihaug-Toint CG for routine/easy cases.
      ! With MOM (state-specific SCF) we MUST stay deterministic: random trial vectors
      ! would break symmetry and collapse the targeted (often excited) state, so force
      ! the conservative Steihaug step that preserves the MOM-selected occupied space.
      ! The micro-solver runs BEFORE the convergence test so the aug-Hessian eigensolve
      ! can detect a negative-curvature (saddle) direction even when |g|~0 -- this is the
      ! stability check that lets us escape an unstable solution a prior DIIS pass landed on.
      if (infos%control%trh_sub_solver == 2 .or. infos%control%mom) then
        call steihaug_cg(infos, conv, g, hdiag, delta, n, nmic, p, pred, micro_used)
      else
        call aughess_step(infos, conv, g, hdiag, delta, n, nmic, p, pred, micro_used)
      end if
      snorm = sqrt(dot_product(p, p))

      ! converged only at a genuine minimum: small gradient AND no escape step.
      if (gnorm < conv_tol .and. snorm < stab_step) then
        ! Stability check (skip under MOM, which must stay on the targeted state):
        ! find the lowest eigenvalue of the orbital Hessian H. If it is negative the
        ! point is a saddle/unstable solution -- kick along that mode and keep going
        ! to descend into the lower (symmetry-broken) minimum.
        if (.not. infos%control%mom) then
          call lowest_hessian_eig(infos, conv, hdiag, n, nmic, lam, vmin)
          if (lam < -stab_eig_tol) then
            write(IW,'(4x,i4,2x,f20.10,2x,es12.4,3x,"unstable (Hess eig ",es10.2,") - escaping")') &
                  macro, e0, gnorm, lam
            call apply_step(conv, infos%control%scftype, 0.1_dp*vmin)
            call build_fock_grad(infos, molgrid, conv, energy, conv%mo_a, conv%mo_b, g, hdiag, e0)
            delta = real(infos%control%trh_r0, dp)
            cycle
          end if
        end if
        write(IW,'(4x,i4,2x,f20.10,2x,es12.4,3x,"CONVERGED")') macro-1, e0, gnorm
        res%error = gnorm
        exit
      end if

      ! model can no longer predict a meaningful ENERGY reduction (pred ~ FP noise).
      ! The energy is converged, but the orbital GRADIENT may still be loose -- near a
      ! minimum pred ~ 0.5*g.p ~ gnorm^2, so it underflows the floor while gnorm is
      ! still ~1e-7. Post-SCF consumers (ROHF canonicalisation for MRSF, analytic
      ! gradients) need accurate orbitals, so take the (descent) Newton/CG step once to
      ! tighten the gradient quadratically before exiting (OTR likewise steps on to
      ! ~1e-8). The energy ratio is FP-noise-dominated here, so accept unconditionally.
      if (pred <= pred_floor .and. gnorm < gtol_fp .and. snorm < stab_step) then
        do while (gnorm > conv_tol .and. snorm > 0.0_dp .and. pred <= pred_floor)
          call apply_step(conv, infos%control%scftype, p)
          call build_fock_grad(infos, molgrid, conv, energy, conv%mo_a, conv%mo_b, g, hdiag, e0)
          gnorm = sqrt(dot_product(g, g) / real(n, dp))
          if (infos%control%trh_sub_solver == 2 .or. infos%control%mom) then
            call steihaug_cg(infos, conv, g, hdiag, delta, n, nmic, p, pred, micro_used)
          else
            call aughess_step(infos, conv, g, hdiag, delta, n, nmic, p, pred, micro_used)
          end if
          snorm = sqrt(dot_product(p, p))
        end do
        ! If a step revived the model (pred > floor) without reaching conv_tol, resume
        ! the normal trust-region loop; otherwise the orbitals are as tight as FP allows.
        if (pred > pred_floor .and. gnorm >= conv_tol) cycle
        write(IW,'(4x,i4,2x,f20.10,2x,es12.4,3x,"CONVERGED (FP precision)")') macro, e0, gnorm
        ! report error below conv_tol so the SCF driver recognises convergence and
        ! does NOT re-diagonalise the raw Fock (which would corrupt ROHF orbitals)
        res%error = min(gnorm, 0.99_dp*conv_tol)
        exit
      end if

      ! trial energy at trial orbitals (copies; conv%mo_* untouched)
      etrial = trial_energy(infos, molgrid, conv, energy, p)
      if (pred > 0.0_dp) then
        rho = (e0 - etrial) / pred
      else
        rho = -1.0_dp
      end if

      accepted = (rho > 0.1_dp)

      ! Backtracking line search: if the full trust step is rejected, shrink it
      ! along its own direction (energy-only evals, no extra Hessian transforms)
      ! until the energy decreases. The first CG direction is -M^{-1}g (descent),
      ! so a sufficiently short step always lowers the energy. This rescues steps
      ! that overshoot into an ascent region (e.g. small-denominator ROHF
      ! docc-socc rotations) instead of collapsing the trust radius.
      if (.not. accepted .and. pred > pred_floor) then
        block
          real(dp) :: fac, et2
          integer  :: ls
          fac = 0.5_dp
          do ls = 1, 5
            et2 = trial_energy(infos, molgrid, conv, energy, fac*p)
            if (et2 < e0 - 1.0e-12_dp) then
              p = fac*p; snorm = fac*snorm; etrial = et2
              rho = (e0 - etrial) / (pred*fac)   ! reduction vs.\ scaled prediction
              accepted = .true.
              exit
            end if
            fac = 0.5_dp*fac
          end do
        end block
      end if
      write(IW,'(4x,i4,2x,f20.10,2x,es12.4,2x,f7.3,2x,f7.3,3x,i4,3x,a)') &
            macro, merge(etrial, e0, accepted), gnorm, rho, delta, micro_used, &
            merge('acc', 'rej', accepted)
      call flush(IW)

      if (accepted) then
        ! commit the rotation, then rebuild gradient/Hessian-diag at the new point
        call apply_step(conv, infos%control%scftype, p)
        call build_fock_grad(infos, molgrid, conv, energy, conv%mo_a, conv%mo_b, g, hdiag, e0)
      end if

      ! trust-radius update
      if (rho < 0.25_dp) then
        delta = 0.25_dp * delta
      else if (rho > 0.75_dp .and. snorm > 0.8_dp*delta) then
        delta = min(2.0_dp*delta, dmax)
      end if

      ! trust region collapsed: converged (small |g|) or a genuine stall (large |g|)
      if (delta < delta_min) then
        if (gnorm < gtol_fp) then
          write(IW,'(4x,i4,2x,f20.10,2x,es12.4,3x,"CONVERGED (trust radius minimal)")') &
                macro, e0, gnorm
          res%error = min(gnorm, 0.99_dp*conv_tol)
        else
          write(IW,'(5X,"Native TRAH: trust region collapsed without convergence, |g|=",ES10.3)') gnorm
          res%error = gnorm
          res%ierr  = 4
        end if
        exit
      end if

      if (macro == nmac) then
        write(IW,'(5X,"Native TRAH: reached max macro iterations.")')
        res%error = gnorm
        res%ierr  = 4
      end if
      end do   ! macro

      conv_ok = (gnorm < gtol_fp)
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
      res%error = min(sqrt(dot_product(g, g)/real(n, dp)), 0.99_dp*conv_tol)
      res%ierr  = 0
      if (nrst > 1) write(IW,'(/5X,"best of ",I0," restarts: E =",F20.10)') nrst, e_best
    end if

    conv%etot = e0
    select type (res)
    class is (scf_conv_trah_result)
      res%iter = macro
    end select

    deallocate(g, hdiag, p, vmin, mo0_a, mo0_b, mob_a, mob_b)
  end subroutine trah_native_run

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

  !> @brief Energy at a trial step p, evaluated on copies of the orbitals.
  function trial_energy(infos, molgrid, conv, energy, p) result(e)
    type(information),    intent(inout), target :: infos
    type(dft_grid_t),     intent(in)            :: molgrid
    type(trah_converger), intent(inout)         :: conv
    type(scf_energy_t),   intent(inout), target :: energy
    real(dp),             intent(in)    :: p(:)
    real(dp) :: e
    real(dp), allocatable :: ma(:,:), mb(:,:)
    type(scf_energy_t), pointer :: ep
    integer :: nschwz
    allocate(ma, source=conv%mo_a)
    allocate(mb, source=conv%mo_b)
    call rotate_mo(conv, infos%control%scftype, p, ma, mb)
    call rebuild_fock(infos, molgrid, conv, energy, ma, mb, nschwz)
    ep => energy
    e = compute_energy(ep)
    deallocate(ma, mb)
  end function trial_energy

  !> @brief Permanently rotate the converger orbitals by p.
  subroutine apply_step(conv, scftype, p)
    type(trah_converger), intent(inout) :: conv
    integer(8),           intent(in)    :: scftype
    real(dp),             intent(in)    :: p(:)
    call rotate_mo(conv, scftype, p, conv%mo_a, conv%mo_b)
  end subroutine apply_step

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

  !> @brief Steihaug-Toint preconditioned CG for the trust-region subproblem.
  !>        Minimises m(p)=g.p+1/2 p.H p with |p|<=delta. Preconditioner M=diag(hdiag).
  !>        H.x via calc_h_op (scaled by 2). Returns step p and pred = -m(p) >= 0.
  subroutine steihaug_cg(infos, conv, g, hdiag, delta, n, nmic, p, pred, used)
    type(information),    intent(inout), target :: infos
    type(trah_converger), intent(inout) :: conv
    real(dp),             intent(in)    :: g(:), hdiag(:), delta
    integer,              intent(in)    :: n, nmic
    real(dp),             intent(out)   :: p(:), pred
    integer,              intent(out)   :: used
    real(dp), allocatable :: r(:), y(:), d(:), hd(:), hx(:), hp(:)
    real(dp) :: ry, ry_new, curv, alpha, beta, tau, pd, dd, pp, gp, php, rnorm0, rnorm
    integer  :: k

    allocate(r(n), y(n), d(n), hd(n), hx(n), hp(n))
    p = 0.0_dp
    hp = 0.0_dp                    ! hp accumulates H.p alongside p (hd=2*calc_h_op is H.d),
    r = g                          ! residual of (H p + g); at p=0, r=g
    call precond(hdiag, r, y)      ! y = M^-1 r
    d = -y
    ry = dot_product(r, y)
    rnorm0 = sqrt(dot_product(r, r))
    used = 0

    do k = 1, nmic
      used = k
      call conv%calc_h_op(infos, d, hx)
      hd = 2.0_dp * hx
      curv = dot_product(d, hd)
      if (curv <= 0.0_dp) then       ! negative curvature -> go to trust boundary
        call to_boundary(p, d, delta, tau)
        p = p + tau*d
        hp = hp + tau*hd
        exit
      end if
      alpha = ry / curv
      ! check trust-region boundary
      pp = dot_product(p, p); pd = dot_product(p, d); dd = dot_product(d, d)
      if (pp + 2.0_dp*alpha*pd + alpha*alpha*dd >= delta*delta) then
        call to_boundary(p, d, delta, tau)
        p = p + tau*d
        hp = hp + tau*hd
        exit
      end if
      p = p + alpha*d
      hp = hp + alpha*hd
      r = r + alpha*hd
      rnorm = sqrt(dot_product(r, r))
      if (rnorm <= min(0.1_dp, sqrt(rnorm0))*rnorm0 .or. rnorm < 1.0e-10_dp) exit
      call precond(hdiag, r, y)
      ry_new = dot_product(r, y)
      beta = ry_new / ry
      d = -y + beta*d
      ry = ry_new
    end do

    ! predicted reduction = -(g.p + 1/2 p.H p); hp already holds H.p, so the
    ! curvature term p.H.p = p.hp -- no extra Hessian-vector (Fock) build needed.
    gp  = dot_product(g, p)
    php = dot_product(p, hp)
    pred = -(gp + 0.5_dp*php)
    deallocate(r, y, d, hd, hx, hp)
  end subroutine steihaug_cg

  !> @brief y = M^-1 r with M = diag(hdiag), small/negative diagonal floored.
  !> @note  Using |hdiag| was tried to handle indefinite ROHF Hessians but
  !>        destabilised well-behaved UHF cases (e.g. [Fe(H2O)6]2+); the robust
  !>        handling of an indefinite Hessian belongs in the micro-solver
  !>        (augmented-Hessian / level shift), not the preconditioner.
  subroutine precond(hdiag, r, y)
    real(dp), intent(in)  :: hdiag(:), r(:)
    real(dp), intent(out) :: y(:)
    integer :: i
    real(dp) :: di
    do i = 1, size(r)
      di = hdiag(i)
      if (di < 1.0e-6_dp) di = 1.0e-6_dp
      y(i) = r(i) / di
    end do
  end subroutine precond

  !> @brief Positive root tau of |p + tau d| = delta (move to trust boundary).
  subroutine to_boundary(p, d, delta, tau)
    real(dp), intent(in)  :: p(:), d(:), delta
    real(dp), intent(out) :: tau
    real(dp) :: a, b, c, disc
    a = dot_product(d, d)
    b = 2.0_dp*dot_product(p, d)
    c = dot_product(p, p) - delta*delta
    disc = max(b*b - 4.0_dp*a*c, 0.0_dp)
    tau = (-b + sqrt(disc)) / (2.0_dp*a)
  end subroutine to_boundary

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

  !> @brief A.[c;y] for the bordered (augmented) Hessian A=[[0,g^T],[g,H]] (alpha=1).
  !>        v(1)=head, v(2:n+1)=orbital rotation. H.y via calc_h_op (scaled by 2).
  subroutine aug_matvec(infos, conv, g, v, n, av)
    type(information),    intent(inout), target :: infos
    type(trah_converger), intent(inout)         :: conv
    real(dp), intent(in)  :: g(:), v(:)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: av(:)
    real(dp), allocatable :: hx(:)
    allocate(hx(n))
    call conv%calc_h_op(infos, v(2:n+1), hx)
    av(1)      = dot_product(g, v(2:n+1))
    av(2:n+1)  = v(1)*g + 2.0_dp*hx
    deallocate(hx)
  end subroutine aug_matvec

  !> @brief Augmented-Hessian (RFO) trust-region step via Davidson for the lowest
  !>        eigenpair of [[0,g^T],[g,H]]. The lowest eigenvector [c;y] yields the
  !>        level-shifted Newton step kappa = y/c, which follows negative curvature
  !>        into the lower (symmetry-broken) basin -- unlike Steihaug-CG, which
  !>        truncates there. Matrix-free (calc_h_op), preconditioned by [0;hdiag],
  !>        step truncated to the trust radius. This is the defining TRAH micro-solver.
  subroutine aughess_step(infos, conv, g, hdiag, delta, n, nmic, p, pred, used)
    type(information),    intent(inout), target :: infos
    type(trah_converger), intent(inout)         :: conv
    real(dp), intent(in)  :: g(:), hdiag(:), delta
    integer,  intent(in)  :: n, nmic
    real(dp), intent(out) :: p(:), pred
    integer,  intent(out) :: used
    integer :: nn, mmax, m, i, k, info, lwork, n_rtv, j, mw, ssz
    real(dp) :: theta, rnorm, c0, snorm, php, di, nv
    real(dp), allocatable :: V(:,:), W(:,:), Tm(:,:), u(:), au(:), r(:), tc(:)
    real(dp), allocatable :: eig(:), work(:), hx(:)
    integer,  allocatable :: seed(:)

    nn   = n + 1
    mmax = min(max(int(nmic), 6), 40)
    allocate(V(nn,mmax), W(nn,mmax), u(nn), au(nn), r(nn), tc(nn), hx(n))
    allocate(eig(mmax), work(8*mmax + 2*mmax*mmax))
    lwork = size(work)

    ! initial subspace, vector 1: gradient-seeded [1 ; -M^{-1} g]
    V(1,1) = 1.0_dp
    do i = 1, n
      di = hdiag(i); if (di < 1.0e-6_dp) di = 1.0e-6_dp
      V(1+i,1) = -g(i)/di
    end do
    V(:,1) = V(:,1)/norm2(V(:,1))
    m = 1

    ! random trial vectors (head 0, random tail): give the eigensolver a generic
    ! component along symmetry-breaking negative-curvature modes that the gradient
    ! cannot see at a symmetric point -- this is how the solver locates broken-
    ! symmetry minima in a single run. Deterministic fixed seed -> reproducible.
    n_rtv = min(max(int(infos%control%trh_nrtv), 1), mmax-1)
    call random_seed(size=ssz); allocate(seed(ssz)); seed = 20260607; call random_seed(put=seed)
    do j = 1, n_rtv
      tc(1) = 0.0_dp
      call random_number(tc(2:nn)); tc(2:nn) = 2.0_dp*tc(2:nn) - 1.0_dp
      do i = 1, m
        tc = tc - dot_product(V(:,i), tc)*V(:,i)
      end do
      nv = norm2(tc)
      if (nv > 1.0e-8_dp) then
        m = m + 1; V(:,m) = tc/nv
      end if
    end do
    deallocate(seed)

    used = 0; mw = 0
    do k = 1, mmax
      ! augmented matrix-vector products for any new basis columns
      do while (mw < m)
        mw = mw + 1
        call aug_matvec(infos, conv, g, V(:,mw), n, W(:,mw)); used = used + 1
      end do
      ! Rayleigh matrix Tm = V^T W  (m x m, symmetric); BLAS over the large nn dim
      allocate(Tm(m,m))
      call dgemm('T','N', m, m, nn, 1.0_dp, V, nn, W, nn, 0.0_dp, Tm, m)
      call dsyev('V','U', m, Tm, m, eig(:m), work, lwork, info)
      theta = eig(1)                       ! lowest eigenvalue
      call dgemv('N', nn, m, 1.0_dp, V, nn, Tm(:,1), 1, 0.0_dp, u, 1)   ! Ritz vector
      call dgemv('N', nn, m, 1.0_dp, W, nn, Tm(:,1), 1, 0.0_dp, au, 1)  ! A u
      deallocate(Tm)
      r = au - theta*u                     ! residual
      rnorm = norm2(r)
      if (rnorm < 1.0e-6_dp .or. m == mmax) exit
      ! preconditioned correction t = (D - theta)^{-1} r,  D = [0 ; hdiag]
      tc(1) = r(1)/sign(max(abs(0.0_dp-theta),1.0e-6_dp), 0.0_dp-theta)
      do i = 1, n
        di = hdiag(i) - theta
        if (abs(di) < 1.0e-6_dp) di = sign(1.0e-6_dp, di)
        tc(1+i) = r(1+i)/di
      end do
      ! orthonormalize against current basis (modified Gram-Schmidt)
      do i = 1, m
        tc = tc - dot_product(V(:,i), tc)*V(:,i)
      end do
      nv = norm2(tc)
      if (nv < 1.0e-8_dp) exit
      m = m + 1
      V(:,m) = tc/nv
    end do

    ! step kappa = y/c from the lowest Ritz vector (head normalized positive)
    c0 = u(1)
    if (c0 < 0.0_dp) then; u = -u; c0 = -c0; end if
    p = u(2:nn)/max(c0, 1.0e-8_dp)
    snorm = norm2(p)
    if (snorm > delta) p = p*(delta/snorm)

    call conv%calc_h_op(infos, p, hx)
    php  = 2.0_dp*dot_product(p, hx)
    pred = -(dot_product(g, p) + 0.5_dp*php)
    deallocate(V, W, u, au, r, tc, eig, work, hx)
  end subroutine aughess_step

  !> @brief Lowest eigenpair of the orbital Hessian H (n x n, NOT bordered), by
  !>        Davidson with several random starts (so it reliably finds a negative
  !>        symmetry-breaking mode at a converged point, where the gradient is ~0
  !>        and the bordered form degenerates). Matrix-free: H.x = 2*calc_h_op.
  subroutine lowest_hessian_eig(infos, conv, hdiag, n, nmic, lam, vmin)
    type(information),    intent(inout), target :: infos
    type(trah_converger), intent(inout)         :: conv
    real(dp), intent(in)  :: hdiag(:)
    integer,  intent(in)  :: n, nmic
    real(dp), intent(out) :: lam, vmin(:)
    integer :: mmax, m, i, k, info, lwork, mw, ssz
    real(dp) :: rnorm, di, nv, theta
    real(dp), allocatable :: V(:,:), W(:,:), Tm(:,:), u(:), r(:), tc(:), eig(:), work(:), hx(:)
    integer,  allocatable :: seed(:)

    mmax = min(max(int(nmic), 8), 40)
    allocate(V(n,mmax), W(n,mmax), u(n), r(n), tc(n), hx(n))
    allocate(eig(mmax), work(8*mmax + 2*mmax*mmax)); lwork = size(work)

    call random_seed(size=ssz); allocate(seed(ssz)); seed = 987654321; call random_seed(put=seed)
    m = 0
    do i = 1, min(4, mmax)                 ! a few random trial vectors
      call random_number(tc); tc = 2.0_dp*tc - 1.0_dp
      do k = 1, m; tc = tc - dot_product(V(:,k), tc)*V(:,k); end do
      nv = norm2(tc); if (nv > 1.0e-8_dp) then; m = m + 1; V(:,m) = tc/nv; end if
    end do
    deallocate(seed)

    theta = 0.0_dp; mw = 0
    do k = 1, mmax
      do while (mw < m)
        mw = mw + 1
        call conv%calc_h_op(infos, V(:,mw), hx); W(:,mw) = 2.0_dp*hx
      end do
      allocate(Tm(m,m))
      call dgemm('T','N', m, m, n, 1.0_dp, V, n, W, n, 0.0_dp, Tm, m)
      call dsyev('V','U', m, Tm, m, eig(:m), work, lwork, info)
      theta = eig(1)
      call dgemv('N', n, m, 1.0_dp, V, n, Tm(:,1), 1, 0.0_dp, u, 1)
      call dgemv('N', n, m, 1.0_dp, W, n, Tm(:,1), 1, 0.0_dp, r, 1)
      r = r - theta*u
      deallocate(Tm)
      rnorm = norm2(r)
      if (rnorm < 1.0e-5_dp .or. m == mmax) exit
      do i = 1, n
        di = hdiag(i) - theta
        if (abs(di) < 1.0e-6_dp) di = sign(1.0e-6_dp, di)
        tc(i) = r(i)/di
      end do
      do i = 1, m; tc = tc - dot_product(V(:,i), tc)*V(:,i); end do
      nv = norm2(tc); if (nv < 1.0e-8_dp) exit
      m = m + 1; V(:,m) = tc/nv
    end do

    lam = theta
    nv = norm2(u); if (nv > 0.0_dp) u = u/nv
    vmin = u
    deallocate(V, W, u, r, tc, hx, eig, work)
  end subroutine lowest_hessian_eig

end module trah_native
