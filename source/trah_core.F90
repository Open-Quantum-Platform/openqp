!> @brief Method-independent core of the trust-region augmented-Hessian (TRAH)
!> optimizer: the macro trust-region loop and its three micro-solvers.
!>
!> `trah_converger.F90` implemented TRAH wired straight into SCF -- its macro
!> loop called `build_fock_grad` to build the Fock itself and `rotate_mo` to
!> rotate SCF orbitals, so nothing else could reuse it.  Everything the
!> algorithm actually needs from the physics is four operations, and they are
!> the four deferred bindings of `trah_provider_t`:
!>
!>   * `grad_hdiag(g, hdiag, e)` -- gradient, Hessian diagonal and energy at the
!>     current point;
!>   * `hess_vec(v, hv)`         -- hv = H.v, matrix-free;
!>   * `trial_energy(p, e)`      -- energy at the current point stepped by p,
!>                                  without committing the step;
!>   * `apply_step(p)`           -- commit the step.
!>
!> SCF supplies these from `build_fock_grad` / `calc_h_op` / `rotate_orbs`
!> (`trah_converger.F90`); CASSCF supplies its own from `cas_evaluate` /
!> `cas_rotate` (`modules/casscf_driver.F90`).  Nothing here knows about Fock
!> matrices, molecular grids, or CI vectors.
!>
!> Sign and scale convention
!> -------------------------
!> The SCF callbacks `calc_g_h` / `calc_h_op` return HALF the true orbital
!> gradient / Hessian, which is the convention `otr_interface` shares with them
!> (it scales by 2).  That convention is NOT re-invented here and it is NOT
!> pushed onto other methods: a provider returns the TRUE gradient, Hessian
!> diagonal and Hessian-vector product, and the SCF provider applies the
!> documented factor of 2 at exactly the places `trah_converger.F90` applied it
!> before.  So the core sees one convention, and the SCF path sees the
!> unchanged one.
!>
!> Algorithm (unchanged from the SCF implementation this was extracted from):
!> trust-region Newton.  Each macro step solves
!>   min_p  g.p + 1/2 p.H p   s.t. |p| <= Delta
!> either by the Steihaug-Toint preconditioned CG (`sub_solver = 2`) or by the
!> augmented-Hessian/RFO Davidson step (the default, and the defining TRAH
!> micro-solver, which follows negative curvature instead of truncating at it).
!> The step is accepted or rejected from the ratio of actual to predicted energy
!> reduction, with a logarithmic (halving) line search rescuing steps that
!> overshoot, and the trust radius is updated from that ratio.  At a point that
!> looks converged the lowest Hessian eigenpair is computed and a negative
!> eigenvalue is escaped along its mode -- the stability check.
!>
!> Method: Helmich-Paris, J. Chem. Phys. 154, 164104 (2021); the
!> augmented-Hessian/level-shift lineage traces to Bacskay, Chem. Phys. 61, 385
!> (1981).
module trah_core_mod

  use precision,     only: dp
  use io_constants,  only: IW
  ! DGEMM/DGEMV/DSYEV through the ILP64 wrapper layer (AGENTS.md rule 1).
  use oqp_linalg

  implicit none
  private

  public :: trah_provider_t, trah_params_t, trah_result_t
  public :: trah_run, trah_micro_step, trah_lowest_hessian_eig

  !> The physics behind one TRAH optimization.  A method implements these four
  !> and gets the whole trust-region machinery.
  type, abstract :: trah_provider_t
    integer :: nparam = 0     !< length of the rotation vector
  contains
    procedure(trah_gh_i), deferred :: grad_hdiag
    procedure(trah_hv_i), deferred :: hess_vec
    procedure(trah_te_i), deferred :: trial_energy
    procedure(trah_as_i), deferred :: apply_step
  end type trah_provider_t

  abstract interface
    !> Gradient, Hessian diagonal and energy at the provider's current point.
    subroutine trah_gh_i(this, g, hdiag, e, ierr)
      import :: trah_provider_t, dp
      class(trah_provider_t), intent(inout) :: this
      real(dp), intent(out) :: g(:), hdiag(:), e
      integer,  intent(out) :: ierr
    end subroutine trah_gh_i

    !> hv = H.v with the TRUE orbital Hessian.  Never assembles H unless the
    !> implementation chooses to.
    subroutine trah_hv_i(this, v, hv, ierr)
      import :: trah_provider_t, dp
      class(trah_provider_t), intent(inout) :: this
      real(dp), intent(in)  :: v(:)
      real(dp), intent(out) :: hv(:)
      integer,  intent(out) :: ierr
    end subroutine trah_hv_i

    !> Energy at (current point stepped by p); the current point is unchanged.
    subroutine trah_te_i(this, p, e, ierr)
      import :: trah_provider_t, dp
      class(trah_provider_t), intent(inout) :: this
      real(dp), intent(in)  :: p(:)
      real(dp), intent(out) :: e
      integer,  intent(out) :: ierr
    end subroutine trah_te_i

    !> Commit the step p.
    subroutine trah_as_i(this, p, ierr)
      import :: trah_provider_t, dp
      class(trah_provider_t), intent(inout) :: this
      real(dp), intent(in) :: p(:)
      integer,  intent(out) :: ierr
    end subroutine trah_as_i
  end interface

  !> Control parameters.  The SCF driver fills these from
  !> `infos%control%trh_*`; the CASSCF driver from its own `[casscf]` block.
  type :: trah_params_t
    integer  :: nmac = 30           !< macroiteration cap
    integer  :: nmic = 32           !< micro-solver iteration cap
    integer  :: sub_solver = 1      !< 1 augmented-Hessian Davidson, 2 Steihaug-CG
    integer  :: nrtv = 1            !< random trial vectors in the AH subspace
    real(dp) :: r0 = 0.4_dp         !< initial trust radius
    !> Trust-radius ceiling.  <= 0 keeps the SCF default max(4, 8*r0); the
    !> CASSCF driver passes the resolved `[casscf] ah_max_trust_radius`, which
    !> itself defaults to `max_rotation_norm`, so the documented cap on |step|
    !> keeps its meaning under this converger too.
    real(dp) :: dmax = 0.0_dp
    real(dp) :: conv_tol = 1.0e-6_dp!< |g| convergence threshold
    !> Deterministic mode (SCF/MOM): forces Steihaug-CG and disables the
    !> stability check, so a targeted state is never symmetry-broken.
    logical  :: deterministic = .false.
    !> |g| measured as an RMS norm (SCF) or as a plain 2-norm (CASSCF).  The
    !> two drivers' `conv` keys mean different things, so the norm is a
    !> parameter rather than a hard-coded choice.
    logical  :: rms_gnorm = .true.
    logical  :: verbose = .true.    !< write the macroiteration table to IW
    !> Optional history callback context: the CASSCF driver records one row per
    !> accepted macroiteration, the SCF driver does not.
    logical  :: want_history = .false.
  end type trah_params_t

  !> What the macro loop reports back.
  type :: trah_result_t
    real(dp) :: energy = 0.0_dp
    real(dp) :: gnorm = 0.0_dp
    real(dp) :: error = 0.0_dp      !< value the SCF driver stores in res%error
    real(dp) :: step_norm = 0.0_dp
    real(dp) :: de = 0.0_dp
    integer  :: iter = 0
    integer  :: ierr = 0            !< 4 = not converged, <0 = provider failure
    logical  :: converged = .false.
  end type trah_result_t

  ! Near-convergence guards: once the model can no longer predict a meaningful
  ! energy reduction (pred below FP noise) or the trust radius collapses while
  ! the gradient is already small, the energy is converged even if |g| has not
  ! reached the (tight) gradient tolerance.  A trust collapse with a large |g|
  ! is instead a genuine stall and is reported as non-convergence.
  real(dp), parameter :: stab_eig_tol = 1.0e-4_dp  !< Hessian eig below -this = unstable
  real(dp), parameter :: pred_floor   = 1.0e-11_dp
  real(dp), parameter :: delta_min    = 1.0e-4_dp
  real(dp), parameter :: gtol_fp      = 1.0e-4_dp
  real(dp), parameter :: stab_step    = 1.0e-3_dp  !< step above this at small |g| = saddle escape

contains

  !> One TRAH macro trust-region loop from the provider's current point.
  !>
  !> On entry the provider is positioned at the starting point; on exit at the
  !> converged (or stalled) one.  `hist_it/hist_e/hist_de/hist_g/hist_s` are
  !> filled for `nhist` accepted macroiterations when `par%want_history`.
  subroutine trah_run(prov, par, res, hist_it, hist_e, hist_de, hist_g, hist_s, nhist)
    class(trah_provider_t), intent(inout) :: prov
    type(trah_params_t),    intent(in)    :: par
    type(trah_result_t),    intent(inout) :: res
    integer,  intent(out), optional :: hist_it(:)
    real(dp), intent(out), optional :: hist_e(:), hist_de(:), hist_g(:), hist_s(:)
    integer,  intent(out), optional :: nhist

    integer  :: n, macro, micro_used, ierr, nh
    real(dp) :: delta, dmax, gnorm, e0, etrial, rho, pred, snorm, lam, obj_old
    real(dp), allocatable :: g(:), hdiag(:), p(:), vmin(:)
    logical  :: accepted

    n     = prov%nparam
    delta = par%r0
    dmax  = merge(par%dmax, max(4.0_dp, 8.0_dp*delta), par%dmax > 0.0_dp)
    nh    = 0
    snorm = 0.0_dp
    res%ierr = 0
    res%converged = .false.

    allocate(g(n), hdiag(n), p(n), vmin(n))

    call prov%grad_hdiag(g, hdiag, e0, ierr)
    if (ierr /= 0) then
      res%ierr = ierr
      return
    end if
    gnorm = gnorm_of(g, n, par%rms_gnorm)
    res%energy = e0
    res%gnorm  = gnorm

    do macro = 1, par%nmac
      res%iter = macro
      gnorm = gnorm_of(g, n, par%rms_gnorm)

      ! trust-region subproblem -> step p, predicted reduction pred.  The
      ! micro-solver runs BEFORE the convergence test so the aug-Hessian
      ! eigensolve can detect a negative-curvature (saddle) direction even when
      ! |g| ~ 0 -- this is the stability check that lets us escape an unstable
      ! solution an earlier converger landed on.
      call trah_micro_step(prov, par, g, hdiag, delta, n, p, pred, micro_used, ierr)
      if (ierr /= 0) then
        res%ierr = ierr
        return
      end if
      snorm = sqrt(dot_product(p, p))

      ! converged only at a genuine minimum: small gradient AND no escape step.
      if (gnorm < par%conv_tol .and. snorm < stab_step) then
        if (.not. par%deterministic) then
          call trah_lowest_hessian_eig(prov, hdiag, n, par%nmic, lam, vmin, ierr)
          if (ierr /= 0) then
            res%ierr = ierr
            return
          end if
          if (lam < -stab_eig_tol) then
            if (par%verbose) write(IW, &
              '(4x,i4,2x,f20.10,2x,es12.4,3x,"unstable (Hess eig ",es10.2,") - escaping")') &
              macro, e0, gnorm, lam
            call prov%apply_step(0.1_dp*vmin, ierr)
            if (ierr == 0) call prov%grad_hdiag(g, hdiag, e0, ierr)
            if (ierr /= 0) then
              res%ierr = ierr
              return
            end if
            delta = par%r0
            cycle
          end if
        end if
        if (par%verbose) write(IW,'(4x,i4,2x,f20.10,2x,es12.4,3x,"CONVERGED")') &
              macro-1, e0, gnorm
        res%error = gnorm
        res%converged = .true.
        exit
      end if

      ! model can no longer predict a meaningful ENERGY reduction (pred ~ FP
      ! noise).  The energy is converged, but the orbital GRADIENT may still be
      ! loose -- near a minimum pred ~ 0.5*g.p ~ gnorm^2, so it underflows the
      ! floor while gnorm is still ~1e-7.  Post-SCF consumers need accurate
      ! orbitals, so take the (descent) step once to tighten the gradient
      ! quadratically before exiting.  The energy ratio is FP-noise-dominated
      ! here, so accept unconditionally.
      if (pred <= pred_floor .and. gnorm < gtol_fp .and. snorm < stab_step) then
        do while (gnorm > par%conv_tol .and. snorm > 0.0_dp .and. pred <= pred_floor)
          call prov%apply_step(p, ierr)
          if (ierr == 0) call prov%grad_hdiag(g, hdiag, e0, ierr)
          if (ierr /= 0) then
            res%ierr = ierr
            return
          end if
          gnorm = gnorm_of(g, n, par%rms_gnorm)
          call trah_micro_step(prov, par, g, hdiag, delta, n, p, pred, micro_used, ierr)
          if (ierr /= 0) then
            res%ierr = ierr
            return
          end if
          snorm = sqrt(dot_product(p, p))
        end do
        if (pred > pred_floor .and. gnorm >= par%conv_tol) cycle
        if (par%verbose) write(IW, &
              '(4x,i4,2x,f20.10,2x,es12.4,3x,"CONVERGED (FP precision)")') macro, e0, gnorm
        ! report error below conv_tol so the SCF driver recognises convergence
        ! and does NOT re-diagonalise the raw Fock (which would corrupt ROHF
        ! orbitals)
        res%error = min(gnorm, 0.99_dp*par%conv_tol)
        res%converged = .true.
        exit
      end if

      ! trial energy at the trial point (the provider's point is untouched)
      call prov%trial_energy(p, etrial, ierr)
      if (ierr /= 0) then
        res%ierr = ierr
        return
      end if
      if (pred > 0.0_dp) then
        rho = (e0 - etrial) / pred
      else
        rho = -1.0_dp
      end if

      accepted = (rho > 0.1_dp)

      ! Backtracking line search: if the full trust step is rejected, shrink it
      ! along its own direction (energy-only evaluations, no extra Hessian
      ! transforms) until the energy decreases.  The first CG direction is
      ! -M^{-1}g (descent), so a sufficiently short step always lowers the
      ! energy.  This rescues steps that overshoot into an ascent region instead
      ! of collapsing the trust radius.
      if (.not. accepted .and. pred > pred_floor) then
        block
          real(dp) :: fac, et2
          integer  :: ls
          fac = 0.5_dp
          do ls = 1, 5
            call prov%trial_energy(fac*p, et2, ierr)
            if (ierr /= 0) then
              res%ierr = ierr
              return
            end if
            if (et2 < e0 - 1.0e-12_dp) then
              p = fac*p; snorm = fac*snorm; etrial = et2
              rho = (e0 - etrial) / (pred*fac)
              accepted = .true.
              exit
            end if
            fac = 0.5_dp*fac
          end do
        end block
      end if
      if (par%verbose) then
        write(IW,'(4x,i4,2x,f20.10,2x,es12.4,2x,f7.3,2x,f7.3,3x,i4,3x,a)') &
              macro, merge(etrial, e0, accepted), gnorm, rho, delta, micro_used, &
              merge('acc', 'rej', accepted)
        call flush(IW)
      end if

      if (accepted) then
        obj_old = e0
        call prov%apply_step(p, ierr)
        if (ierr == 0) call prov%grad_hdiag(g, hdiag, e0, ierr)
        if (ierr /= 0) then
          res%ierr = ierr
          return
        end if
        if (par%want_history .and. present(nhist)) then
          if (nh < size(hist_e)) then
            nh = nh + 1
            hist_it(nh) = macro
            hist_e(nh)  = e0
            hist_de(nh) = e0 - obj_old
            hist_g(nh)  = gnorm_of(g, n, par%rms_gnorm)
            hist_s(nh)  = snorm
          end if
        end if
      end if

      ! trust-radius update
      if (rho < 0.25_dp) then
        delta = 0.25_dp * delta
      else if (rho > 0.75_dp .and. snorm > 0.8_dp*delta) then
        delta = min(2.0_dp*delta, dmax)
      end if

      ! trust region collapsed: converged (small |g|) or a genuine stall.
      ! `gnorm` is deliberately the value from the top of this macroiteration.
      if (delta < delta_min) then
        if (gnorm < gtol_fp) then
          if (par%verbose) write(IW, &
            '(4x,i4,2x,f20.10,2x,es12.4,3x,"CONVERGED (trust radius minimal)")') macro, e0, gnorm
          res%error = min(gnorm, 0.99_dp*par%conv_tol)
          res%converged = .true.
        else
          if (par%verbose) write(IW, &
            '(5X,"Native TRAH: trust region collapsed without convergence, |g|=",ES10.3)') gnorm
          res%error = gnorm
          res%ierr  = 4
        end if
        exit
      end if

      if (macro == par%nmac) then
        if (par%verbose) write(IW,'(5X,"Native TRAH: reached max macro iterations.")')
        res%error = gnorm
        res%ierr  = 4
      end if
    end do

    res%energy = e0
    res%gnorm  = gnorm_of(g, n, par%rms_gnorm)
    res%step_norm = snorm
    if (present(nhist)) nhist = nh
    deallocate(g, hdiag, p, vmin)
  end subroutine trah_run

  !> |g|, as an RMS norm (SCF `conv`) or a plain 2-norm (CASSCF `grad_tol`).
  pure function gnorm_of(g, n, rms) result(val)
    real(dp), intent(in) :: g(:)
    integer,  intent(in) :: n
    logical,  intent(in) :: rms
    real(dp) :: val
    val = dot_product(g, g)
    if (rms) then
      val = sqrt(val / real(n, dp))
    else
      val = sqrt(val)
    end if
  end function gnorm_of

  !> Trust-region subproblem: whichever micro-solver the parameters select.
  !>
  !> The augmented-Hessian Davidson step is the default and is what makes this
  !> TRAH rather than plain trust-region Newton -- it follows negative curvature
  !> into the lower basin instead of truncating there.  With a state-targeted
  !> run (`deterministic`) the random trial vectors would break symmetry and
  !> collapse the targeted state, so the conservative Steihaug step is forced.
  subroutine trah_micro_step(prov, par, g, hdiag, delta, n, p, pred, used, ierr)
    class(trah_provider_t), intent(inout) :: prov
    type(trah_params_t),    intent(in)    :: par
    real(dp), intent(in)  :: g(:), hdiag(:), delta
    integer,  intent(in)  :: n
    real(dp), intent(out) :: p(:), pred
    integer,  intent(out) :: used, ierr
    if (par%sub_solver == 2 .or. par%deterministic) then
      call steihaug_cg(prov, g, hdiag, delta, n, par%nmic, p, pred, used, ierr)
    else
      call aughess_step(prov, par, g, hdiag, delta, n, par%nmic, p, pred, used, ierr)
    end if
  end subroutine trah_micro_step

  !> @brief Steihaug-Toint preconditioned CG for the trust-region subproblem.
  !>        Minimises m(p)=g.p+1/2 p.H p with |p|<=delta.  Preconditioner
  !>        M=diag(hdiag).  H.x via the provider.  Returns p and pred = -m(p).
  subroutine steihaug_cg(prov, g, hdiag, delta, n, nmic, p, pred, used, ierr)
    class(trah_provider_t), intent(inout) :: prov
    real(dp), intent(in)  :: g(:), hdiag(:), delta
    integer,  intent(in)  :: n, nmic
    real(dp), intent(out) :: p(:), pred
    integer,  intent(out) :: used, ierr
    real(dp), allocatable :: r(:), y(:), d(:), hd(:), hp(:)
    real(dp) :: ry, ry_new, curv, alpha, beta, tau, pd, dd, pp, gp, php, rnorm0, rnorm
    integer  :: k

    ierr = 0
    allocate(r(n), y(n), d(n), hd(n), hp(n))
    p  = 0.0_dp
    hp = 0.0_dp                    ! hp accumulates H.p alongside p
    r  = g                         ! residual of (H p + g); at p=0, r=g
    call precond(hdiag, r, y)      ! y = M^-1 r
    d = -y
    ry = dot_product(r, y)
    rnorm0 = sqrt(dot_product(r, r))
    used = 0

    do k = 1, nmic
      used = k
      call prov%hess_vec(d, hd, ierr)
      if (ierr /= 0) return
      curv = dot_product(d, hd)
      if (curv <= 0.0_dp) then       ! negative curvature -> go to trust boundary
        call to_boundary(p, d, delta, tau)
        p = p + tau*d
        hp = hp + tau*hd
        exit
      end if
      alpha = ry / curv
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
    ! curvature term costs no extra Hessian-vector product.
    gp  = dot_product(g, p)
    php = dot_product(p, hp)
    pred = -(gp + 0.5_dp*php)
    deallocate(r, y, d, hd, hp)
  end subroutine steihaug_cg

  !> @brief y = M^-1 r with M = diag(hdiag), small/negative diagonal floored.
  !> @note  Using |hdiag| was tried to handle indefinite ROHF Hessians but
  !>        destabilised well-behaved UHF cases; the robust handling of an
  !>        indefinite Hessian belongs in the micro-solver (augmented-Hessian /
  !>        level shift), not the preconditioner.
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

  !> @brief A.[c;y] for the bordered (augmented) Hessian A=[[0,g^T],[g,H]].
  subroutine aug_matvec(prov, g, v, n, av, ierr)
    class(trah_provider_t), intent(inout) :: prov
    real(dp), intent(in)  :: g(:), v(:)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: av(:)
    integer,  intent(out) :: ierr
    real(dp), allocatable :: hx(:)
    allocate(hx(n))
    call prov%hess_vec(v(2:n+1), hx, ierr)
    av(1)     = dot_product(g, v(2:n+1))
    av(2:n+1) = v(1)*g + hx
    deallocate(hx)
  end subroutine aug_matvec

  !> @brief Augmented-Hessian (RFO) trust-region step via Davidson for the
  !>        lowest eigenpair of [[0,g^T],[g,H]].  The lowest eigenvector [c;y]
  !>        yields the level-shifted Newton step kappa = y/c, which follows
  !>        negative curvature into the lower (symmetry-broken) basin -- unlike
  !>        Steihaug-CG, which truncates there.  Matrix-free, preconditioned by
  !>        [0;hdiag], step truncated to the trust radius.  This is the defining
  !>        TRAH micro-solver.
  subroutine aughess_step(prov, par, g, hdiag, delta, n, nmic, p, pred, used, ierr)
    class(trah_provider_t), intent(inout) :: prov
    type(trah_params_t),    intent(in)    :: par
    real(dp), intent(in)  :: g(:), hdiag(:), delta
    integer,  intent(in)  :: n, nmic
    real(dp), intent(out) :: p(:), pred
    integer,  intent(out) :: used, ierr
    integer :: nn, mmax, m, i, k, info, lwork, n_rtv, j, mw, ssz
    real(dp) :: theta, rnorm, c0, snorm, php, di, nv
    real(dp), allocatable :: V(:,:), W(:,:), Tm(:,:), u(:), au(:), r(:), tc(:)
    real(dp), allocatable :: eig(:), work(:), hx(:)
    integer,  allocatable :: seed(:)

    ierr = 0
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

    ! random trial vectors (head 0, random tail): give the eigensolver a
    ! generic component along symmetry-breaking negative-curvature modes that
    ! the gradient cannot see at a symmetric point.  Deterministic fixed seed.
    n_rtv = min(max(par%nrtv, 1), mmax-1)
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
      do while (mw < m)
        mw = mw + 1
        call aug_matvec(prov, g, V(:,mw), n, W(:,mw), ierr)
        if (ierr /= 0) return
        used = used + 1
      end do
      ! Rayleigh matrix Tm = V^T W  (m x m, symmetric); BLAS over the large nn dim
      allocate(Tm(m,m))
      call dgemm('T','N', m, m, nn, 1.0_dp, V, nn, W, nn, 0.0_dp, Tm, m)
      call dsyev('V','U', m, Tm, m, eig(:m), work, lwork, info)
      theta = eig(1)
      call dgemv('N', nn, m, 1.0_dp, V, nn, Tm(:,1), 1, 0.0_dp, u, 1)   ! Ritz vector
      call dgemv('N', nn, m, 1.0_dp, W, nn, Tm(:,1), 1, 0.0_dp, au, 1)  ! A u
      deallocate(Tm)
      r = au - theta*u
      rnorm = norm2(r)
      if (rnorm < 1.0e-6_dp .or. m == mmax) exit
      ! preconditioned correction t = (D - theta)^{-1} r,  D = [0 ; hdiag]
      tc(1) = r(1)/sign(max(abs(0.0_dp-theta),1.0e-6_dp), 0.0_dp-theta)
      do i = 1, n
        di = hdiag(i) - theta
        if (abs(di) < 1.0e-6_dp) di = sign(1.0e-6_dp, di)
        tc(1+i) = r(1+i)/di
      end do
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

    call prov%hess_vec(p, hx, ierr)
    if (ierr /= 0) return
    php  = dot_product(p, hx)
    pred = -(dot_product(g, p) + 0.5_dp*php)
    deallocate(V, W, u, au, r, tc, eig, work, hx)
  end subroutine aughess_step

  !> @brief Lowest eigenpair of the orbital Hessian H (n x n, NOT bordered), by
  !>        Davidson with several random starts (so it reliably finds a negative
  !>        symmetry-breaking mode at a converged point, where the gradient is
  !>        ~0 and the bordered form degenerates).  Matrix-free.
  subroutine trah_lowest_hessian_eig(prov, hdiag, n, nmic, lam, vmin, ierr)
    class(trah_provider_t), intent(inout) :: prov
    real(dp), intent(in)  :: hdiag(:)
    integer,  intent(in)  :: n, nmic
    real(dp), intent(out) :: lam, vmin(:)
    integer,  intent(out) :: ierr
    integer :: mmax, m, i, k, info, lwork, mw, ssz
    real(dp) :: rnorm, di, nv, theta
    real(dp), allocatable :: V(:,:), W(:,:), Tm(:,:), u(:), r(:), tc(:), eig(:), work(:)
    integer,  allocatable :: seed(:)

    ierr = 0
    mmax = min(max(int(nmic), 8), 40)
    allocate(V(n,mmax), W(n,mmax), u(n), r(n), tc(n))
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
        call prov%hess_vec(V(:,mw), W(:,mw), ierr)
        if (ierr /= 0) return
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
    deallocate(V, W, u, r, tc, eig, work)
  end subroutine trah_lowest_hessian_eig

end module trah_core_mod
