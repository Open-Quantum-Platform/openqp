! The augmented-Hessian (default) TRAH micro-solver must deliver a Newton step
! at every gradient size, so the macro iteration stays quadratic down to a
! tight tolerance.  The model is E = 1/2 x.H x with a dense Hessian whose
! diagonal is a poor preconditioner (spectrum 0.01..2, eigenvectors mixing all
! coordinates), as the DFT orbital Hessian is.  A Davidson stopped by a fixed
! absolute residual returns the diagonal-preconditioned gradient step once
! |g| is below that residual; that step contracts by about 1 - 0.01/diag(H)
! per iteration, and the run stalls far above the tolerance.
module trah_ah_model
 use precision, only: dp
 use trah_core_mod
 implicit none
 integer, parameter :: nm = 24
 type, extends(trah_provider_t) :: quad
  real(dp) :: h(nm,nm) = 0, x(nm) = 0
  integer :: steps = 0
 contains
  procedure :: grad_hdiag => q_gh
  procedure :: hess_vec => q_hv
  procedure :: trial_energy => q_te
  procedure :: apply_step => q_ap
 end type
contains
 subroutine q_gh(this, g, hdiag, e, ierr)
  class(quad), intent(inout) :: this
  real(dp), intent(out) :: g(:), hdiag(:), e
  integer, intent(out) :: ierr
  integer :: i
  g = matmul(this%h, this%x)
  e = 0.5_dp*dot_product(this%x, g)
  do i = 1, nm
   hdiag(i) = this%h(i,i)
  end do
  ierr = 0
 end subroutine
 subroutine q_hv(this, v, hv, ierr)
  class(quad), intent(inout) :: this
  real(dp), intent(in) :: v(:)
  real(dp), intent(out) :: hv(:)
  integer, intent(out) :: ierr
  hv = matmul(this%h, v)
  ierr = 0
 end subroutine
 subroutine q_te(this, p, e, ierr)
  class(quad), intent(inout) :: this
  real(dp), intent(in) :: p(:)
  real(dp), intent(out) :: e
  integer, intent(out) :: ierr
  e = 0.5_dp*dot_product(this%x + p, matmul(this%h, this%x + p))
  ierr = 0
 end subroutine
 subroutine q_ap(this, p, ierr)
  class(quad), intent(inout) :: this
  real(dp), intent(in) :: p(:)
  integer, intent(out) :: ierr
  this%x = this%x + p
  this%steps = this%steps + 1
  ierr = 0
 end subroutine
end module

program check_trah_ah_newton
 use trah_ah_model
 implicit none
 type(quad) :: m
 type(trah_params_t) :: par
 type(trah_result_t) :: res
 real(dp) :: q(nm,nm), w1(nm), w2(nm), lam(nm), g(nm), t, gconv
 integer :: i, j
 ! Q = (I - 2 w1 w1^T)(I - 2 w2 w2^T): a deterministic orthogonal matrix that
 ! mixes every coordinate.
 do i = 1, nm
  w1(i) = sin(1.3_dp*i) + 0.2_dp
  w2(i) = cos(0.7_dp*i*i) - 0.1_dp
 end do
 w1 = w1/norm2(w1); w2 = w2/norm2(w2)
 q = 0
 do i = 1, nm
  q(i,i) = 1
 end do
 q = q - 2*spread(w1, 2, nm)*spread(w1, 1, nm)
 q = matmul(q, reflector(w2))
 do i = 1, nm
  lam(i) = 0.01_dp*(200.0_dp)**(real(i-1, dp)/real(nm-1, dp))
 end do
 m%h = matmul(q, matmul(diag(lam), transpose(q)))
 ! Start with an RMS gradient of 2e-6: the regime the DFT runs reached
 ! quadratically before stalling.
 do j = 1, nm
  m%x(j) = sin(2.1_dp*j) + 0.3_dp*cos(5.0_dp*j)
 end do
 g = matmul(m%h, m%x)
 t = sqrt(dot_product(g, g)/nm)
 m%x = m%x*(2.0e-6_dp/t)
 m%nparam = nm
 m%hess_vec_is_linear = .true.
 par%sub_solver = 1; par%nrtv = 1; par%nmic = 32
 par%conv_tol = 1.0e-11_dp; par%nmac = 30; par%verbose = .false.
 call trah_run(m, par, res)
 g = matmul(m%h, m%x)
 if (.not. res%converged .or. res%ierr /= 0) error stop 'AH TRAH did not converge'
 if (sqrt(dot_product(g, g)/nm) >= par%conv_tol) error stop 'AH TRAH gradient above tolerance'
 ! Quadratic: a few Newton steps, not a long diagonal-preconditioned descent.
 if (m%steps > 6) error stop 'AH TRAH took too many steps (not Newton)'
 gconv = sqrt(dot_product(g, g)/nm)
 ! The |g|-relative stop must never be looser than the old absolute 1e-6.
 ! With conv_tol = 1e-5 the 0.1*conv_tol*sqrt(n) floor is 4.9e-6; at
 ! |g| ~ 1e-3 the old solver stopped at 1e-6, so the step must too.  For the
 ! normalised Ritz vector [1; p]/|[1; p]| the residual is
 ! (g.p - theta, g + (H - theta) p)/|[1; p]| with theta = g.p.
 block
  type(trah_params_t) :: loose
  real(dp) :: p(nm), hp(nm), hd(nm), theta, res, pred
  integer :: used, ierr, k, ntrial
  loose = par
  loose%conv_tol = 1.0e-5_dp
  ntrial = 0
  do k = 1, 6
   do j = 1, nm
    m%x(j) = sin(0.37_dp*j*k) + 0.5_dp*cos(1.9_dp*j + k)
   end do
   g = matmul(m%h, m%x)
   g = g*(1.0e-3_dp/norm2(g))
   do i = 1, nm
    hd(i) = m%h(i,i)
   end do
   call trah_micro_step(m, loose, g, hd, 10.0_dp, nm, p, pred, used, ierr)
   if (ierr /= 0) error stop 'loose-tolerance micro-step failed'
   hp = matmul(m%h, p)
   theta = dot_product(g, p)
   res = sqrt(dot_product(g + hp - theta*p, g + hp - theta*p))/sqrt(1.0_dp + dot_product(p, p))
   if (res > 1.0e-6_dp*(1.0_dp + 1.0e-6_dp)) then
    print '(a,i0,a,es10.3,a,i0)', 'trial ', k, ': residual ', res, ' after products ', used
    error stop 'Davidson stopped above the old 1e-6 residual'
   end if
   ntrial = ntrial + 1
  end do
 end block
 print '(a,i0,a,es10.3)', 'PASS: augmented-Hessian Newton steps=', m%steps, &
       ' |g|rms=', gconv
contains
 function reflector(w) result(r)
  real(dp), intent(in) :: w(:)
  real(dp) :: r(size(w), size(w))
  integer :: k
  r = -2*spread(w, 2, size(w))*spread(w, 1, size(w))
  do k = 1, size(w)
   r(k,k) = r(k,k) + 1
  end do
 end function
 function diag(d) result(r)
  real(dp), intent(in) :: d(:)
  real(dp) :: r(size(d), size(d))
  integer :: k
  r = 0
  do k = 1, size(d)
   r(k,k) = d(k)
  end do
 end function
end program
