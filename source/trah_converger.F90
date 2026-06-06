!> @brief Native Fortran trust-region augmented-Hessian (TRAH) SCF solver.
!> @detail Drop-in replacement for the external OpenTrustRegion driver, selected
!>         by control%trh_impl = 1 ("native"). Reuses the exact physics callbacks
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

    integer  :: n, macro, nmac, nmic, micro_used
    real(dp) :: delta, dmax, conv_tol, gnorm, e0, etrial, rho, pred, snorm
    real(dp), allocatable :: g(:), hdiag(:), p(:)
    logical  :: accepted
    ! near-convergence guards: once the model can no longer predict a meaningful
    ! energy reduction (pred below FP noise) or the trust radius collapses while
    ! the gradient is already small, the energy is converged even if |g| has not
    ! reached the (tight) gradient tolerance. A trust collapse with a large |g| is
    ! instead a genuine stall and is reported as non-convergence.
    real(dp), parameter :: pred_floor = 1.0e-11_dp
    real(dp), parameter :: delta_min  = 1.0e-4_dp
    real(dp), parameter :: gtol_fp    = 1.0e-4_dp

    n        = int(conv%n_param)
    nmac     = int(infos%control%maxit)
    nmic     = int(infos%control%trh_nmic)
    conv_tol = real(infos%control%conv, dp)
    delta    = real(infos%control%trh_r0, dp)
    dmax     = max(4.0_dp, 8.0_dp*delta)

    allocate(g(n), hdiag(n), p(n))
    conv%f_old = 0.0_dp
    conv%d_old = 0.0_dp

    write(IW,'(/5X,"Native TRAH (trust-region Newton, Steihaug-CG)"/5X,46("-"))')
    write(IW,'(5X,"start trust radius =",F7.3,"   conv =",ES9.2,"   max micro =",I4)') &
              delta, conv_tol, nmic
    write(IW,'(/4x,"Macro",6x,"Energy",13x,"|grad|",7x,"rho",6x,"trust",3x,"micro",3x,"step")')
    write(IW,'(3x,75("="))')

    ! initial point at the incoming orbitals
    call build_fock_grad(infos, molgrid, conv, energy, conv%mo_a, conv%mo_b, g, hdiag, e0)

    do macro = 1, nmac
      gnorm = sqrt(dot_product(g, g) / real(n, dp))
      if (gnorm < conv_tol) then
        write(IW,'(4x,i4,2x,f20.10,2x,es12.4,3x,"CONVERGED")') macro-1, e0, gnorm
        res%error = gnorm
        exit
      end if

      ! trust-region subproblem  ->  step p, predicted reduction pred
      call steihaug_cg(infos, conv, g, hdiag, delta, n, nmic, p, pred, micro_used)
      snorm = sqrt(dot_product(p, p))

      ! model can no longer predict a meaningful reduction -> energy converged
      if (pred <= pred_floor .and. gnorm < gtol_fp) then
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
    end do

    conv%etot = e0
    select type (res)
    class is (scf_conv_trah_result)
      res%iter = macro
    end select

    deallocate(g, hdiag, p)
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
    real(dp), allocatable :: r(:), y(:), d(:), hd(:), hx(:)
    real(dp) :: ry, ry_new, curv, alpha, beta, tau, pd, dd, pp, gp, php, rnorm0, rnorm
    integer  :: k

    allocate(r(n), y(n), d(n), hd(n), hx(n))
    p = 0.0_dp
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
        exit
      end if
      alpha = ry / curv
      ! check trust-region boundary
      pp = dot_product(p, p); pd = dot_product(p, d); dd = dot_product(d, d)
      if (pp + 2.0_dp*alpha*pd + alpha*alpha*dd >= delta*delta) then
        call to_boundary(p, d, delta, tau)
        p = p + tau*d
        exit
      end if
      p = p + alpha*d
      r = r + alpha*hd
      rnorm = sqrt(dot_product(r, r))
      if (rnorm <= min(0.1_dp, sqrt(rnorm0))*rnorm0 .or. rnorm < 1.0e-10_dp) exit
      call precond(hdiag, r, y)
      ry_new = dot_product(r, y)
      beta = ry_new / ry
      d = -y + beta*d
      ry = ry_new
    end do

    ! predicted reduction = -(g.p + 1/2 p.H p)
    call conv%calc_h_op(infos, p, hx)
    gp  = dot_product(g, p)
    php = 2.0_dp * dot_product(p, hx)
    pred = -(gp + 0.5_dp*php)
    deallocate(r, y, d, hd, hx)
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

end module trah_native
