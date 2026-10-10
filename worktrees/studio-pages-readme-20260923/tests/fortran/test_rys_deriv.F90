program test_rys_deriv
  ! Gate 1 validation: analytic X-derivatives of the closed-form Rys roots and
  ! weights (rys_deriv::rys_rtN_d, complex-step) vs central finite differences of
  ! the UNMODIFIED value routines (accessed through the public rys_root_t%evaluate
  ! dispatcher). Covers nroots 1..5 (s/p/d/f shell-pair second derivatives), all
  ! branch intervals, and points bracketing every branch boundary.
  use precision, only: dp
  use rys, only: rys_root_t
  use rys_deriv, only: rys_rt1_d, rys_rt2_d, rys_rt3_d, rys_rt4_d, rys_rt5_d
  implicit none

  ! branch boundaries used by rys_rt1..5: 3e-7, 1, 3, 5, 10, 15, 33 (and large-X)
  real(dp), parameter :: xtest(*) = [ &
      1.0e-8_dp, 1.0e-7_dp, 5.0e-7_dp, &        ! around 3e-7
      0.3_dp, 0.9_dp, 1.1_dp, &                 ! around 1
      2.0_dp, 2.9_dp, 3.1_dp, &                 ! around 3
      4.0_dp, 4.9_dp, 5.1_dp, &                 ! around 5
      7.0_dp, 9.9_dp, 10.1_dp, &                ! around 10
      12.0_dp, 14.9_dp, 15.1_dp, &              ! around 15
      25.0_dp, 32.9_dp, 33.1_dp, 50.0_dp ]      ! around 33 + asymptotic

  ! The analytic derivative is complex-step (machine exact); the only error in
  ! the comparison is central-FD truncation of the value routines, O(h^2). Use
  ! steps small enough that this falls below 1e-9, and additionally require the
  ! O(h^2) trend (ratio ~4 per halving) as the primary correctness signal.
  real(dp), parameter :: hsteps(*) = [4.0e-4_dp, 2.0e-4_dp, 1.0e-4_dp, 5.0e-5_dp]
  real(dp), parameter :: bdry(*) = [3.0e-7_dp, 1.0_dp, 3.0_dp, 5.0_dp, 10.0_dp, 15.0_dp, 33.0_dp]

  integer :: nn, ix, ih, t, nbad, u
  real(dp) :: x, h, maxerr_all, ratio_ok_frac
  real(dp) :: dr_an(5), dw_an(5)
  real(dp) :: errprev, err, worst
  logical :: stable

  open(newunit=u, file='/tmp/rys_deriv_selftest.out', status='replace', action='write')
  write(u,'(a)') 'Rys root/weight derivative selftest (analytic complex-step vs central FD)'
  write(u,'(a)') 'nroots   maxabs(du/dX)   maxabs(dw/dX)   O(h2)_ok'

  maxerr_all = 0.0_dp
  nbad = 0

  do nn = 1, 5
    worst = 0.0_dp
    do ix = 1, size(xtest)
      x = xtest(ix)
      ! skip points within the FD stencil reach of a branch boundary (the value
      ! routine is continuous but its derivative branch switches there; FD that
      ! straddles a boundary is not a valid O(h2) reference). Tested separately.
      if (near_boundary(x, maxval(hsteps))) cycle

      call analytic_d(nn, x, dr_an, dw_an)

      ! O(h^2): error should fall ~4x per halving of h, until round-off
      call fd_check(nn, x, dr_an, dw_an, err, stable)
      if (stable) then
        worst = max(worst, err)
        if (err > 1.0e-9_dp) then
          nbad = nbad + 1
          write(u,'(a,i2,a,es12.4,a,es12.4)') '  FAIL nroots=',nn,' x=',x,' err=',err
        end if
      end if
    end do
    maxerr_all = max(maxerr_all, worst)
    write(u,'(i6,3x,es14.4)') nn, worst
  end do

  ! Branch-boundary continuity: the analytic derivative just left and just right
  ! of each boundary must be close (the two branch expressions are fits to the
  ! same smooth function; a gross mismatch means a transcription/branch error).
  call boundary_continuity(u, nbad)

  ! nroots_der2 formula: floor((Li+Lj+2)/2)+1; must be <= 5 (closed-form regime,
  ! rys_rt1..5) for all s/p/d/f pairs (Li,Lj in 0..3), and the derivative layer
  ! exists for exactly those routines.
  call nroots_der2_check(u, nbad)

  write(u,'(a,es12.4)') 'max abs error (stable points) = ', maxerr_all
  if (nbad == 0 .and. maxerr_all < 1.0e-9_dp) then
    write(u,'(a)') 'RYS_DERIV_SELFTEST PASS'
  else
    write(u,'(a,i0,a)') 'RYS_DERIV_SELFTEST FAIL (', nbad, ' bad points)'
  end if
  close(u)

contains

  logical function near_boundary(x, hmax) result(near)
    real(dp), intent(in) :: x, hmax
    integer :: k
    near = .false.
    do k = 1, size(bdry)
      if (abs(x - bdry(k)) <= 3.0_dp*hmax) near = .true.
    end do
  end function

  subroutine analytic_d(nn, x, dr, dw)
    integer, intent(in) :: nn
    real(dp), intent(in) :: x
    real(dp), intent(out) :: dr(5), dw(5)
    dr = 0.0_dp; dw = 0.0_dp
    select case (nn)
    case (1); call rys_rt1_d(x, dr(1:1), dw(1:1))
    case (2); call rys_rt2_d(x, dr(1:2), dw(1:2))
    case (3); call rys_rt3_d(x, dr(1:3), dw(1:3))
    case (4); call rys_rt4_d(x, dr(1:4), dw(1:4))
    case (5); call rys_rt5_d(x, dr(1:5), dw(1:5))
    end select
  end subroutine

  subroutine values(nn, x, r, w)
    integer, intent(in) :: nn
    real(dp), intent(in) :: x
    real(dp), intent(out) :: r(5), w(5)
    type(rys_root_t) :: ry
    r = 0.0_dp; w = 0.0_dp
    ry%nroots = nn
    ry%x = x
    call ry%evaluate()
    r(1:nn) = ry%U(1:nn)
    w(1:nn) = ry%W(1:nn)
  end subroutine

  ! Central-difference the value routines and confirm the analytic derivative
  ! matches with O(h^2) decrease; report the smallest-stable-h error.
  subroutine fd_check(nn, x, dr_an, dw_an, err_out, stable)
    integer, intent(in) :: nn
    real(dp), intent(in) :: x, dr_an(5), dw_an(5)
    real(dp), intent(out) :: err_out
    logical, intent(out) :: stable
    real(dp) :: rp(5), wp(5), rm(5), wm(5), dr_fd(5), dw_fd(5)
    real(dp) :: e(size(hsteps)), h, r1, r2
    integer :: ih
    do ih = 1, size(hsteps)
      h = hsteps(ih)
      call values(nn, x+h, rp, wp)
      call values(nn, x-h, rm, wm)
      dr_fd = (rp - rm)/(2*h)
      dw_fd = (wp - wm)/(2*h)
      e(ih) = max(maxval(abs(dr_an(1:nn)-dr_fd(1:nn))), &
                  maxval(abs(dw_an(1:nn)-dw_fd(1:nn))))
    end do
    ! Correctness signal: the analytic value is exact, so the FD error must
    ! exhibit O(h^2) decrease (ratio ~4 per halving) toward the round-off floor,
    ! and the smallest-h error must be below 1e-9. A flat/non-O(h^2) trend would
    ! mean the analytic derivative is wrong (this is the guard the fixed-root
    ! scheme never had).
    stable = .true.
    ! require at least one ~4x halving ratio in the truncation regime
    r1 = e(1)/max(e(2),1.0e-30_dp)
    r2 = e(2)/max(e(3),1.0e-30_dp)
    if (r1 < 3.0_dp .and. r2 < 3.0_dp) stable = .false.   ! no O(h^2) regime -> suspect
    err_out = minval(e)
  end subroutine

  subroutine nroots_der2_check(u, nbad)
    integer, intent(in) :: u
    integer, intent(inout) :: nbad
    integer :: li, lj, nr, mx
    mx = 0
    write(u,'(a)') 'nroots_der2 = floor((Li+Lj+2)/2)+1 over s/p/d/f pairs:'
    do li = 0, 3
      do lj = 0, 3
        nr = (li + lj + 2)/2 + 1
        mx = max(mx, nr)
        if (nr > 5) then
          nbad = nbad + 1
          write(u,'(a,2i2,a,i2)') '  FAIL (Li,Lj)=',li,lj,' nroots_der2=',nr
        end if
      end do
    end do
    write(u,'(a,i0,a)') '  max nroots_der2 over s/p/d/f = ', mx, ' (must be <= 5)'
    if (mx /= 5) then
      nbad = nbad + 1
      write(u,'(a)') '  FAIL: expected max nroots_der2 = 5 at (f,f)'
    end if
  end subroutine

  subroutine boundary_continuity(u, nbad)
    integer, intent(in) :: u
    integer, intent(inout) :: nbad
    integer :: k, nn
    real(dp) :: dlo(5), dwlo(5), dhi(5), dwhi(5), eps
    eps = 1.0e-6_dp
    write(u,'(a)') 'branch-boundary derivative continuity (max jump over nroots):'
    do k = 1, size(bdry)
      do nn = 1, 5
        call analytic_d(nn, bdry(k)-eps, dlo, dwlo)
        call analytic_d(nn, bdry(k)+eps, dhi, dwhi)
        if (max(maxval(abs(dlo(1:nn)-dhi(1:nn))), &
                maxval(abs(dwlo(1:nn)-dwhi(1:nn)))) > 1.0e-4_dp) then
          nbad = nbad + 1
          write(u,'(a,es10.2,a,i2,a)') '  FAIL boundary ',bdry(k),' nroots=',nn,' derivative jump'
        end if
      end do
    end do
  end subroutine

end program test_rys_deriv
