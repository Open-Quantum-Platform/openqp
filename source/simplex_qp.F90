module simplex_qp
  use precision, only: dp
  use eigen, only: diag_symm_full
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  integer, parameter, public :: SIMPLEX_QP_EXACT_MAX_N = 12
  integer, parameter, public :: SIMPLEX_QP_SUCCESS = 0
  integer, parameter, public :: SIMPLEX_QP_APPROXIMATE = 1
  integer, parameter, public :: SIMPLEX_QP_INVALID_INPUT = 2
  integer, parameter, public :: SIMPLEX_QP_NONFINITE_INPUT = 3
  integer, parameter, public :: SIMPLEX_QP_EIGEN_FAILURE = 4

  public :: solve_simplex_qp, simplex_qp_objective

contains

  subroutine solve_simplex_qp(h, g, x, value, status, forbidden, nforbidden, &
                              forbid_vertices_before, preferred, repeat_tolerance)
    real(dp), intent(in) :: h(:, :), g(:)
    real(dp), intent(out) :: x(:), value
    integer, intent(out) :: status
    real(dp), optional, intent(in) :: forbidden(:, :)
    integer, optional, intent(in) :: nforbidden, forbid_vertices_before, preferred
    real(dp), optional, intent(in) :: repeat_tolerance
    real(dp), allocatable :: hs(:, :), gs(:), best(:), best_any(:), candidate(:)
    real(dp) :: input_scale, best_value, best_any_value, repeat_tol
    integer :: n, pref, nforbid, old_vertices
    logical :: have_best, have_best_any, eigen_failed

    n = size(g)
    x = 0.0_dp
    value = huge(1.0_dp)
    status = SIMPLEX_QP_INVALID_INPUT
    if (n < 1 .or. size(x) < n) return

    pref = n
    if (present(preferred)) pref = max(1, min(n, preferred))
    x(pref) = 1.0_dp
    if (size(h, 1) < n .or. size(h, 2) < n) return

    if (.not. all(ieee_is_finite(h(1:n, 1:n))) .or. &
        .not. all(ieee_is_finite(g))) then
      status = SIMPLEX_QP_NONFINITE_INPUT
      return
    end if

    repeat_tol = 1.0e-4_dp
    if (present(repeat_tolerance)) repeat_tol = max(0.0_dp, repeat_tolerance)
    nforbid = 0
    if (present(forbidden) .and. present(nforbidden)) then
      if (size(forbidden, 1) >= n) &
        nforbid = max(0, min(nforbidden, size(forbidden, 2)))
    end if
    old_vertices = 0
    if (present(forbid_vertices_before)) &
      old_vertices = max(0, min(n, forbid_vertices_before))

    allocate(hs(n, n), gs(n), best(n), best_any(n), candidate(n))
    input_scale = max(maxval(abs(h(1:n, 1:n))), maxval(abs(g)))
    if (input_scale > 0.0_dp) then
      hs = 0.5_dp*(h(1:n, 1:n)/input_scale + &
                   transpose(h(1:n, 1:n))/input_scale)
      ! Subtracting a constant from every linear coefficient leaves a simplex
      ! objective unchanged up to an additive constant and improves scaling.
      gs = g/input_scale - g(pref)/input_scale
    else
      hs = 0.0_dp
      gs = 0.0_dp
    end if

    have_best = .false.
    have_best_any = .false.
    eigen_failed = .false.
    best_value = huge(1.0_dp)
    best_any_value = huge(1.0_dp)

    if (n == 1 .or. &
        (maxval(abs(hs)) == 0.0_dp .and. maxval(abs(gs)) == 0.0_dp)) then
      candidate = 0.0_dp
      candidate(pref) = 1.0_dp
      call consider(candidate)
      status = SIMPLEX_QP_SUCCESS
    else if (n <= SIMPLEX_QP_EXACT_MAX_N) then
      call enumerate_faces()
      status = SIMPLEX_QP_SUCCESS
      if (eigen_failed) then
        call projected_fallback()
        status = SIMPLEX_QP_EIGEN_FAILURE
      end if
    else
      call projected_fallback()
      status = SIMPLEX_QP_APPROXIMATE
    end if

    ! If the unconstrained minimizer was rejected as a repeated solution, probe
    ! deterministic points on the segment toward the preferred (latest) state.
    if (have_best_any) call consider_repeat_boundary()

    x = 0.0_dp
    if (have_best) then
      x(1:n) = best
    else
      x(pref) = 1.0_dp
    end if
    value = simplex_qp_objective(h(1:n, 1:n), g, x(1:n))

  contains

    subroutine enumerate_faces()
      integer(int64) :: mask, last_mask
      integer :: i, j, k, m, info
      integer, allocatable :: idx(:)
      real(dp), allocatable :: hf(:, :), z(:, :), r(:, :), eval(:), &
                               grad0(:), rhs(:), spectral(:), y(:), cf(:)
      real(dp) :: curv_tol, compat_tol, feas_tol, denom
      logical :: compatible

      allocate(idx(n), hf(n, n), z(n, n - 1), r(n - 1, n - 1), &
               eval(n - 1), grad0(n), rhs(n - 1), spectral(n - 1), &
               y(n - 1), cf(n))
      last_mask = 2_int64**n - 1_int64
      do mask = 1_int64, last_mask
        k = 0
        do i = 1, n
          if (btest(mask, i - 1)) then
            k = k + 1
            idx(k) = i
          end if
        end do

        candidate = 0.0_dp
        if (k == 1) then
          candidate(idx(1)) = 1.0_dp
          call consider(candidate)
          cycle
        end if

        hf(1:k, 1:k) = hs(idx(1:k), idx(1:k))
        grad0(1:k) = matmul(hf(1:k, 1:k), &
                            spread(1.0_dp/real(k, dp), 1, k)) + gs(idx(1:k))

        ! Helmert basis for the tangent space sum(c)=0 on this face.
        z(1:k, 1:k - 1) = 0.0_dp
        do j = 1, k - 1
          denom = sqrt(real(j*(j + 1), dp))
          z(1:j, j) = 1.0_dp/denom
          z(j + 1, j) = -real(j, dp)/denom
        end do

        m = k - 1
        r(1:m, 1:m) = matmul(transpose(z(1:k, 1:m)), &
                              matmul(hf(1:k, 1:k), z(1:k, 1:m)))
        r(1:m, 1:m) = 0.5_dp*(r(1:m, 1:m) + transpose(r(1:m, 1:m)))
        call diag_symm_full(1, m, r, n - 1, eval, info)
        if (info /= 0) then
          eigen_failed = .true.
          cycle
        end if

        curv_tol = 512.0_dp*epsilon(1.0_dp)* &
                   max(1.0_dp, maxval(abs(eval(1:m))))
        if (minval(eval(1:m)) < -curv_tol) cycle

        rhs(1:m) = matmul(transpose(z(1:k, 1:m)), grad0(1:k))
        spectral(1:m) = matmul(transpose(r(1:m, 1:m)), rhs(1:m))
        compat_tol = 2048.0_dp*epsilon(1.0_dp)* &
                     max(1.0_dp, maxval(abs(rhs(1:m))))
        compatible = .true.
        do i = 1, m
          if (abs(eval(i)) > curv_tol) then
            spectral(i) = -spectral(i)/eval(i)
          else if (abs(spectral(i)) <= compat_tol) then
            spectral(i) = 0.0_dp
          else
            compatible = .false.
            exit
          end if
        end do
        if (.not. compatible) cycle

        y(1:m) = matmul(r(1:m, 1:m), spectral(1:m))
        cf(1:k) = 1.0_dp/real(k, dp) + matmul(z(1:k, 1:m), y(1:m))
        feas_tol = 4096.0_dp*epsilon(1.0_dp)* &
                   max(1.0_dp, maxval(abs(cf(1:k))))
        if (minval(cf(1:k)) < -feas_tol) cycle
        candidate(idx(1:k)) = max(cf(1:k), 0.0_dp)
        call consider(candidate)
      end do
    end subroutine enumerate_faces

    subroutine projected_fallback()
      integer :: i, start, iter, nstarts, best_linear
      real(dp), allocatable :: current(:), trial(:), gradient(:)
      real(dp) :: lipschitz, step, qcurrent, qtrial

      allocate(current(n), trial(n), gradient(n))
      candidate = 1.0_dp/real(n, dp)
      call consider(candidate)
      do i = 1, n
        candidate = 0.0_dp
        candidate(i) = 1.0_dp
        call consider(candidate)
      end do

      lipschitz = maxval(sum(abs(hs), dim=2))
      best_linear = minloc(gs, dim=1)
      nstarts = min(n, 32) + 3
      do start = 1, nstarts
        current = 0.0_dp
        select case (start)
        case (1)
          current = 1.0_dp/real(n, dp)
        case (2)
          current(pref) = 1.0_dp
        case (3)
          current(best_linear) = 1.0_dp
        case default
          current(start - 3) = 1.0_dp
        end select

        if (lipschitz <= 64.0_dp*epsilon(1.0_dp)) then
          call consider(current)
          cycle
        end if
        step = 1.0_dp/lipschitz
        qcurrent = scaled_objective(current)
        do iter = 1, 500
          gradient = matmul(hs, current) + gs
          trial = current - step*gradient
          call project_simplex(trial)
          qtrial = scaled_objective(trial)
          if (qtrial > qcurrent + &
              256.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(qcurrent))) then
            step = 0.5_dp*step
            if (step <= tiny(1.0_dp)) exit
            cycle
          end if
          if (sqrt(sum((trial - current)**2)) <= &
              128.0_dp*epsilon(1.0_dp)* &
              max(1.0_dp, sqrt(sum(current**2)))) then
            current = trial
            exit
          end if
          current = trial
          qcurrent = qtrial
        end do
        call consider(current)
      end do
    end subroutine projected_fallback

    subroutine project_simplex(v)
      real(dp), intent(inout) :: v(:)
      real(dp), allocatable :: sorted(:)
      real(dp) :: cumulative, theta, key
      integer :: i, j

      allocate(sorted(size(v)))
      sorted = v
      do i = 2, size(sorted)
        key = sorted(i)
        j = i - 1
        do while (j >= 1)
          if (sorted(j) >= key) exit
          sorted(j + 1) = sorted(j)
          j = j - 1
        end do
        sorted(j + 1) = key
      end do
      cumulative = 0.0_dp
      theta = sorted(1) - 1.0_dp
      do i = 1, size(sorted)
        cumulative = cumulative + sorted(i)
        if (sorted(i) - (cumulative - 1.0_dp)/real(i, dp) > 0.0_dp) &
          theta = (cumulative - 1.0_dp)/real(i, dp)
      end do
      v = max(v - theta, 0.0_dp)
      if (sum(v) > 0.0_dp) v = v/sum(v)
    end subroutine project_simplex

    subroutine consider(raw)
      real(dp), intent(in) :: raw(:)
      real(dp) :: normalized(n), q, total

      if (.not. all(ieee_is_finite(raw(1:n))) .or. &
          minval(raw(1:n)) < -1.0e-9_dp) return
      normalized = max(raw(1:n), 0.0_dp)
      total = sum(normalized)
      if (.not. ieee_is_finite(total) .or. total <= tiny(1.0_dp)) return
      normalized = normalized/total
      q = scaled_objective(normalized)
      if (.not. ieee_is_finite(q)) return

      ! Keep the branches separate: Fortran does not require short-circuit
      ! evaluation, and each incumbent is undefined until its first candidate.
      if (.not. have_best_any) then
        best_any = normalized
        best_any_value = q
        have_best_any = .true.
      else if (better(q, normalized, best_any_value, best_any)) then
        best_any = normalized
        best_any_value = q
      end if

      if (candidate_allowed(normalized)) then
        if (.not. have_best) then
          best = normalized
          best_value = q
          have_best = .true.
        else if (better(q, normalized, best_value, best)) then
          best = normalized
          best_value = q
        end if
      end if
    end subroutine consider

    logical function candidate_allowed(c)
      real(dp), intent(in) :: c(:)
      integer :: j

      candidate_allowed = .false.
      if (old_vertices > 0) then
        if (any(c(1:old_vertices) >= 1.0_dp - repeat_tol)) return
      end if
      if (nforbid > 0) then
        do j = 1, nforbid
          if (.not. all(ieee_is_finite(forbidden(1:n, j)))) cycle
          if (sqrt(sum((c - forbidden(1:n, j))**2)) < repeat_tol) return
        end do
      end if
      candidate_allowed = .true.
    end function candidate_allowed

    logical function better(q, c, incumbent_q, incumbent)
      real(dp), intent(in) :: q, c(:), incumbent_q, incumbent(:)
      real(dp) :: tie_tol
      integer :: i

      ! Equal objectives favor the requested (normally latest) state first,
      ! then successively newer/higher-index states. This makes every
      ! semidefinite face reproducible without privileging old SCF history.
      tie_tol = 1024.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp, abs(q), abs(incumbent_q))
      if (q < incumbent_q - tie_tol) then
        better = .true.
        return
      end if
      if (abs(q - incumbent_q) > tie_tol) then
        better = .false.
        return
      end if
      if (c(pref) > incumbent(pref) + tie_tol) then
        better = .true.
        return
      else if (c(pref) < incumbent(pref) - tie_tol) then
        better = .false.
        return
      end if
      do i = n, 1, -1
        if (c(i) > incumbent(i) + tie_tol) then
          better = .true.
          return
        else if (c(i) < incumbent(i) - tie_tol) then
          better = .false.
          return
        end if
      end do
      better = .false.
    end function better

    subroutine consider_repeat_boundary()
      real(dp), parameter :: alpha(18) = [ &
        2.0e-4_dp, 5.0e-4_dp, 1.0e-3_dp, 2.0e-3_dp, 5.0e-3_dp, &
        1.0e-2_dp, 2.0e-2_dp, 5.0e-2_dp, 1.0e-1_dp, 2.0e-1_dp, &
        3.0e-1_dp, 4.0e-1_dp, 5.0e-1_dp, 6.0e-1_dp, 7.0e-1_dp, &
        8.0e-1_dp, 9.0e-1_dp, 1.0_dp]
      real(dp) :: vertex(n)
      integer :: i

      vertex = 0.0_dp
      vertex(pref) = 1.0_dp
      do i = 1, size(alpha)
        candidate = (1.0_dp - alpha(i))*best_any + alpha(i)*vertex
        call consider(candidate)
      end do
    end subroutine consider_repeat_boundary

    real(dp) function scaled_objective(c)
      real(dp), intent(in) :: c(:)

      scaled_objective = 0.5_dp*dot_product(c, matmul(hs, c)) + &
                         dot_product(gs, c)
    end function scaled_objective

  end subroutine solve_simplex_qp

  real(dp) function simplex_qp_objective(h, g, x) result(value)
    real(dp), intent(in) :: h(:, :), g(:), x(:)

    value = 0.5_dp*dot_product(x, matmul(h, x)) + dot_product(g, x)
  end function simplex_qp_objective

end module simplex_qp
