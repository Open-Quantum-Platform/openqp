!> @brief Fortran engine for the CASSCF orbital-optimization convergers of
!> pyoqp casscf_convergers.py.
!>
!> Three bind(C) entry points replace the NumPy in the converger module (the
!> Python implementations remain as the numerical pin and as the fallback when
!> the symbols are unavailable):
!>
!>  * casscf_ah_model_step     -- trust-region augmented-Hessian model step
!>  * casscf_lowest_mode_step  -- saddle-escape step along the lowest mode
!>  * casscf_diis_coeffs       -- Pulay extrapolation coefficients
!>
!> Augmented-Hessian model step
!> ----------------------------
!> Working in the Hessian eigenbasis H = U diag(w) U^T, for a scale `alpha` the
!> lowest eigenpair (lambda, v) of the bordered matrix
!>
!>     [[0, alpha g_e^T], [alpha g_e, diag(w)]],      g_e = U^T g
!>
!> gives x = v[1:] / (alpha v[0]), which solves (H - lambda) x = -g with
!> lambda <= min(0, w_min): a positive-definite shifted Hessian, hence a
!> descent direction even through negative curvature.  Larger alpha means a
!> more negative shift and a shorter step, so the microiterations first grow
!> alpha geometrically until |x| is inside the trust radius, then bisect in
!> log-alpha until |x| lands in [0.8, 1.0] * trust.  The control flow, the
!> microiteration budget and the fall-through cases are transcribed from the
!> Python one-for-one, including the hard rescale when the budget runs out
!> while the step is still outside the region.
!>
!> The bordered eigenproblem is solved DENSELY, with LAPACK DSYEVD, exactly as
!> the Python does through `fci._symmetric_eigh` -> `oqp_dsyevd`.  This is
!> deliberate and must stay that way.  The bordered matrix is an *arrowhead*
!> matrix, so its lowest eigenpair is formally the root of a scalar secular
!> equation and the O(n^3) eigensolve collapses to O(n) per Newton step; that
!> shortcut was prototyped and rejected before this engine was written.  Over
!> 720 randomized cases it got the eigenvalue right to 2.7e-15 relative but the
!> step only to ~3e-4, because x_i divides by (lambda - w_i) ~ 1e-11 and an
!> eigenvalue accurate to 1e-15 absolute is not a step accurate to much.  In
!> the regime the optimizer actually converges into -- small gradient, one deep
!> negative mode -- the root lies nearer a Hessian eigenvalue than double
!> precision resolves and the iteration lands on the pole: 115 of 180 such
!> cases returned a finite but WRONG step where the dense path correctly
!> reports no reference component and the caller falls back to the lowest-mode
!> step.  Doing it properly needs LAPACK dlaed4's pole-offset treatment (solve
!> for d = lambda - w_k so the small denominator is never formed by
!> cancellation).  Until someone writes that, the dense solve stays: it is the
!> correct, better-tested code.
!>
!> Arrays cross in the caller's C order (last index fastest).  The eigenvector
!> matrix is NOT symmetric, so its C-order buffer is not its column-major view;
!> it is gathered index-explicitly into ut(i,k) = U[i,k] and the two
!> projections then differ in transpose mode -- ge = U^T g is the transposed
!> DGEMV, step = U x the untransposed one.  Getting that pair the wrong way
!> round is silent: both shapes are square, so it produces a plausible,
!> completely wrong step rather than a crash.
module casscf_ah_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: casscf_ah_model_step, casscf_lowest_mode_step, casscf_diis_coeffs

contains

  !> Lowest eigenpair of the bordered matrix [[0, alpha ge^T], [alpha ge, W]].
  !>
  !> Fills `x(0:n-1)` with v[1:]/(alpha v[0]) and returns `lam`.  `status` is
  !> 0 on success, 1 when the reference component v[0] is below `v0tol` (the
  !> pure negative-curvature case the caller handles by stepping along the
  !> lowest Hessian mode), and 2 when LAPACK failed.
  subroutine ah_solve(n, ge, w, alpha, v0tol, lam, x, status)
    integer, intent(in) :: n
    real(dp), intent(in) :: ge(0:n-1), w(0:n-1), alpha, v0tol
    real(dp), intent(out) :: lam, x(0:n-1)
    integer, intent(out) :: status

    real(dp), allocatable :: a(:,:), evals(:), work(:)
    integer, allocatable :: iwork(:)
    real(dp) :: qwork(1), v0
    integer :: m, i, info, lwork, liwork, iqwork(1)

    m = n + 1
    status = 0
    lam = 0.0_dp
    x = 0.0_dp

    allocate(a(0:m-1, 0:m-1), evals(0:m-1))
    a = 0.0_dp
    do i = 0, n - 1
      a(0, i + 1) = alpha * ge(i)
      a(i + 1, 0) = alpha * ge(i)
      a(i + 1, i + 1) = w(i)
    end do

    ! Workspace query, then the solve.  DSYEVD returns eigenvalues ascending
    ! and overwrites `a` with the eigenvectors as columns.
    call dsyevd('V', 'U', m, a, m, evals, qwork, -1, iqwork, -1, info)
    if (info /= 0) then
      deallocate(a, evals)
      status = 2
      return
    end if
    lwork = int(qwork(1))
    liwork = int(iqwork(1))
    allocate(work(lwork), iwork(liwork))
    call dsyevd('V', 'U', m, a, m, evals, work, lwork, iwork, liwork, info)
    deallocate(work, iwork)
    if (info /= 0) then
      deallocate(a, evals)
      status = 2
      return
    end if

    lam = evals(0)
    v0 = a(0, 0)
    if (abs(v0) < v0tol) then
      status = 1
      deallocate(a, evals)
      return
    end if
    do i = 0, n - 1
      x(i) = a(i + 1, 0) / (alpha * v0)
    end do
    deallocate(a, evals)
  end subroutine ah_solve

  !> Quadratic model value ge.x + 1/2 sum_i w_i x_i^2.
  pure function ah_model(n, ge, w, x) result(val)
    integer, intent(in) :: n
    real(dp), intent(in) :: ge(0:n-1), w(0:n-1), x(0:n-1)
    real(dp) :: val
    integer :: i
    val = 0.0_dp
    do i = 0, n - 1
      val = val + ge(i) * x(i) + 0.5_dp * w(i) * x(i) * x(i)
    end do
  end function ah_model

  pure function vnorm(n, x) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(0:n-1)
    real(dp) :: r
    integer :: i
    r = 0.0_dp
    do i = 0, n - 1
      r = r + x(i) * x(i)
    end do
    r = sqrt(r)
  end function vnorm

  !> Trust-region augmented-Hessian model step.
  !>
  !> @param[in]  npar      number of orbital-rotation parameters
  !> @param[in]  grad      orbital gradient, [npar]
  !> @param[in]  w         Hessian eigenvalues (ascending), [npar]
  !> @param[in]  umat      Hessian eigenvectors, C-order [npar,npar] (columns
  !>                       are the modes, i.e. umat[i][k] = U[i,k])
  !> @param[in]  trust     trust radius
  !> @param[in]  max_micro microiteration budget
  !> @param[in]  v0tol     reference-component cutoff
  !> @param[out] step      full-space step U x, [npar]
  !> @param[out] shift     the AH level shift lambda
  !> @param[out] pred      predicted model decrease
  !> @param[out] nmic      microiterations consumed
  !> @return     0 on success; 1 when the AH eigenvector has no reference
  !>             component (caller steps along the lowest mode instead);
  !>             2 when LAPACK failed (caller falls back to the Python path)
  function casscf_ah_model_step(npar, grad, w, umat, trust, max_micro, v0tol, &
                                step, shift, pred, nmic) result(status) &
      bind(C, name="casscf_ah_model_step")
    integer(c_int32_t), value :: npar, max_micro
    real(dp), value :: trust, v0tol
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: grad(0:*), w(0:*), umat(0:*)
    real(dp), intent(inout) :: step(0:*), shift(0:*), pred(0:*)
    integer(c_int32_t), intent(inout) :: nmic(0:*)
    integer(c_int32_t) :: status

    real(dp), allocatable :: ge(:), x(:), x_in(:), ut(:,:)
    real(dp) :: lam, s_in, a_lo, a_hi, a_mid, xn
    integer :: n, i, j, st, mic, budget
    logical :: have_in

    n = int(npar)
    budget = max(1, int(max_micro))
    status = 0_c_int32_t
    shift(0) = 0.0_dp
    pred(0) = 0.0_dp
    nmic(0) = 1_c_int32_t
    if (n <= 0) return

    allocate(ge(0:n-1), x(0:n-1), x_in(0:n-1), ut(0:n-1, 0:n-1))

    ! Gather the eigenvectors index-explicitly so ut(i,k) = U[i,k] in the
    ! ordinary mathematical sense, whatever the buffer order.
    do i = 0, n - 1
      do j = 0, n - 1
        ut(i, j) = umat(int(i, i8)*int(n, i8) + int(j, i8))
      end do
    end do
    ! ge = U^T g:  ge(k) = sum_i U[i,k] g(i) = sum_i ut(i,k) g(i), i.e. the
    ! TRANSPOSED matvec; the step below is the untransposed one.
    call dgemv('T', n, n, 1.0_dp, ut, n, grad, 1, 0.0_dp, ge, 1)

    mic = 1
    call ah_solve(n, ge, w, 1.0_dp, v0tol, lam, x, st)
    shift(0) = lam
    nmic(0) = int(mic, c_int32_t)
    if (st /= 0) then
      status = int(st, c_int32_t)
      deallocate(ge, x, x_in, ut)
      return
    end if
    if (vnorm(n, x) <= trust) then
      pred(0) = ah_model(n, ge, w, x)
      call dgemv('N', n, n, 1.0_dp, ut, n, x, 1, 0.0_dp, step, 1)
      deallocate(ge, x, x_in, ut)
      return
    end if

    ! ---- step longer than the radius: grow alpha until inside, then bisect
    a_lo = 1.0_dp
    a_hi = 1.0_dp
    have_in = .false.
    s_in = lam
    do while (mic < budget)
      a_hi = a_hi * 4.0_dp
      mic = mic + 1
      call ah_solve(n, ge, w, a_hi, v0tol, lam, x, st)
      shift(0) = lam
      nmic(0) = int(mic, c_int32_t)
      if (st /= 0) then
        status = int(st, c_int32_t)
        deallocate(ge, x, x_in, ut)
        return
      end if
      if (vnorm(n, x) <= trust) then
        x_in = x
        s_in = lam
        have_in = .true.
        exit
      end if
      a_lo = a_hi
    end do

    if (.not. have_in) then
      ! budget exhausted while still outside: hard rescale
      xn = vnorm(n, x)
      if (xn > 0.0_dp) x = x * (trust / xn)
      shift(0) = lam
      pred(0) = ah_model(n, ge, w, x)
      call dgemv('N', n, n, 1.0_dp, ut, n, x, 1, 0.0_dp, step, 1)
      nmic(0) = int(mic, c_int32_t)
      deallocate(ge, x, x_in, ut)
      return
    end if

    do while (mic < budget .and. vnorm(n, x_in) < 0.8_dp * trust)
      a_mid = sqrt(a_lo * a_hi)
      mic = mic + 1
      call ah_solve(n, ge, w, a_mid, v0tol, lam, x, st)
      if (st == 2) then
        status = 2_c_int32_t
        deallocate(ge, x, x_in, ut)
        return
      end if
      if (st == 1) exit
      if (vnorm(n, x) <= trust) then
        x_in = x
        s_in = lam
        a_hi = a_mid
      else
        a_lo = a_mid
      end if
    end do

    shift(0) = s_in
    pred(0) = ah_model(n, ge, w, x_in)
    call dgemv('N', n, n, 1.0_dp, ut, n, x_in, 1, 0.0_dp, step, 1)
    nmic(0) = int(mic, c_int32_t)
    deallocate(ge, x, x_in, ut)
  end function casscf_ah_model_step

  !> Saddle-escape trial step along the lowest Hessian mode (downhill sign).
  !>
  !> @param[in]  npar  number of orbital-rotation parameters
  !> @param[in]  grad  orbital gradient, [npar]
  !> @param[in]  w     Hessian eigenvalues (ascending), [npar]
  !> @param[in]  umat  Hessian eigenvectors, C-order [npar,npar]
  !> @param[in]  trust trust radius
  !> @param[out] step  trial step, [npar]
  !> @param[out] pred  predicted model change
  subroutine casscf_lowest_mode_step(npar, grad, w, umat, trust, step, pred) &
      bind(C, name="casscf_lowest_mode_step")
    integer(c_int32_t), value :: npar
    real(dp), value :: trust
    real(dp), intent(in) :: grad(0:*), w(0:*), umat(0:*)
    real(dp), intent(inout) :: step(0:*), pred(0:*)

    integer :: n, i
    real(dp) :: overlap, sgn, gdotstep

    n = int(npar)
    pred(0) = 0.0_dp
    if (n <= 0) return

    ! u = U[:,0]; the C-order buffer gives U[i,0] at umat(i*n + 0).
    overlap = 0.0_dp
    do i = 0, n - 1
      overlap = overlap + grad(i) * umat(int(i, i8)*int(n, i8))
    end do
    sgn = 1.0_dp
    if (overlap > 0.0_dp) sgn = -1.0_dp

    do i = 0, n - 1
      step(i) = (sgn * trust) * umat(int(i, i8)*int(n, i8))
    end do
    ! step = sgn*trust*u, so g.step is sgn*trust*(g.u) -- reuse the overlap
    ! rather than summing a second time.
    gdotstep = sgn * trust * overlap
    pred(0) = gdotstep + 0.5_dp * w(0) * trust * trust
  end subroutine casscf_lowest_mode_step

  !> Pulay (DIIS) coefficients minimizing |sum_i c_i g_i| subject to sum c_i = 1.
  !>
  !> Drops the oldest stored vector while the bordered B matrix is
  !> ill-conditioned, exactly as the Python does, and returns the coefficients
  !> for the LAST `nused` entries.
  !>
  !> @param[in]  nvec    number of stored gradient vectors
  !> @param[in]  npar    length of each gradient vector
  !> @param[in]  gmat    gradients, C-order [nvec,npar] (row i is g_i)
  !> @param[in]  condmax conditioning ceiling (the Python uses 1e14)
  !> @param[out] coef    coefficients for the last `nused` vectors, [nvec]
  !> @param[out] nused   number of coefficients written; 0 means "no stable
  !>                     extrapolation exists" (the Python returns None)
  subroutine casscf_diis_coeffs(nvec, npar, gmat, condmax, coef, nused) &
      bind(C, name="casscf_diis_coeffs")
    integer(c_int32_t), value :: nvec, npar
    real(dp), value :: condmax
    real(dp), intent(in) :: gmat(0:*)
    real(dp), intent(inout) :: coef(0:*)
    integer(c_int32_t), intent(inout) :: nused(0:*)

    real(dp), allocatable :: gram(:,:), b(:,:), bs(:,:), rhs(:), sv(:), work(:)
    integer, allocatable :: ipiv(:)
    real(dp) :: acc, qwork(1), cond
    integer :: nv, np, start, n, m, i, j, k, info, lwork
    logical :: ok

    nv = int(nvec)
    np = int(npar)
    nused(0) = 0_c_int32_t
    if (nv < 2 .or. np <= 0) return

    ! Full Gram matrix once; every retry is a trailing sub-block of it.
    !
    ! gram(i,j) = sum_k gmat[i,k] gmat[j,k] is G^T G, and the caller's C-order
    ! [nvec, npar] buffer viewed column-major is exactly the (npar x nvec)
    ! matrix whose columns are the error vectors -- so this is one DSYRK rather
    ! than a scalar triple loop.  BLAS fills one triangle; the mirror is
    ! explicit because the retry loop below reads both.
    allocate(gram(0:nv-1, 0:nv-1))
    call dsyrk('U', 'T', nv, np, 1.0_dp, gmat, np, 0.0_dp, gram, nv)
    do i = 0, nv - 1
      do j = i + 1, nv - 1
        gram(j, i) = gram(i, j)
      end do
    end do

    start = 0
    do while (nv - start >= 2)
      n = nv - start
      m = n + 1
      allocate(b(0:m-1, 0:m-1), bs(0:m-1, 0:m-1), rhs(0:m-1))
      b = 0.0_dp
      do i = 0, n - 1
        do j = 0, n - 1
          b(i, j) = gram(start + i, start + j)
        end do
        b(i, n) = 1.0_dp
        b(n, i) = 1.0_dp
      end do
      b(n, n) = 0.0_dp
      rhs = 0.0_dp
      rhs(n) = 1.0_dp

      ok = all(b == b)                       ! NaN screen
      do i = 0, m - 1
        do j = 0, m - 1
          if (abs(b(i, j)) > huge(1.0_dp)) ok = .false.
        end do
      end do

      if (ok) then
        ! 2-norm condition number from the singular values, matching
        ! numpy.linalg.cond's default.
        bs = b
        allocate(sv(0:m-1))
        call dgesvd('N', 'N', m, m, bs, m, sv, bs, m, bs, m, qwork, -1, info)
        if (info == 0) then
          lwork = int(qwork(1))
          allocate(work(lwork))
          bs = b
          call dgesvd('N', 'N', m, m, bs, m, sv, bs, m, bs, m, work, lwork, info)
          deallocate(work)
        end if
        if (info /= 0) then
          ok = .false.
        else if (sv(m-1) <= 0.0_dp) then
          ok = .false.
        else
          cond = sv(0) / sv(m-1)
          if (.not. (cond <= condmax)) ok = .false.
        end if
        deallocate(sv)
      end if

      if (ok) then
        allocate(ipiv(m))
        bs = b
        call dgesv(m, 1, bs, m, ipiv, rhs, m, info)
        deallocate(ipiv)
        if (info /= 0) then
          ok = .false.
        else
          do i = 0, n - 1
            if (rhs(i) /= rhs(i) .or. abs(rhs(i)) > huge(1.0_dp)) ok = .false.
          end do
        end if
      end if

      if (ok) then
        do i = 0, n - 1
          coef(i) = rhs(i)
        end do
        nused(0) = int(n, c_int32_t)
        deallocate(b, bs, rhs, gram)
        return
      end if

      deallocate(b, bs, rhs)
      start = start + 1
    end do

    deallocate(gram)
  end subroutine casscf_diis_coeffs

end module casscf_ah_mod
