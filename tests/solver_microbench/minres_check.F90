program minres_check
  use precision, only: dp
  use iso_c_binding, only: c_null_ptr
  use bench_common
  use minres_mod
  use pcg_mod
  implicit none

  integer, parameter :: n = 5000
  integer, parameter :: mxit = 20000
  real(dp), parameter :: tol = 1.0e-9_dp
  integer, target :: dummy
  real(dp), allocatable :: b(:), xm(:), xp(:), tmp(:)
  real(dp) :: relm, relp, sol_diff, errm, errp
  integer :: k, it_m, it_p

  allocate(b(n), xm(n), xp(n), tmp(n))
  do k = 1, n
    b(k) = 1.0_dp + 0.5_dp*sin(real(k,dp))
  end do

  ! ================= TEST 1: SPD system (CG valid) =========================
  a_diag = 2.0_dp; a_offd = 0.49_dp; op_shift = 0.0_dp; precond_identity = .false.
  inject_mode = 0

  xm = b
  call minres_optimize(xm, matvec, precond, dummy, mxit, tol=tol, err=errm, iters=it_m)
  relm = true_resid(b, xm)

  call solve_pcg(b, xp, it_p, relp)

  tmp = xm - xp
  sol_diff = maxval(abs(tmp))
  write(*,'(A)') '===== TEST 1: SPD system (diag-dominant, Jacobi precond) ====='
  write(*,'(A,I6,A,ES11.4)') ' MINRES iters=', it_m, '  ||b-Ax||/||b||=', relm
  write(*,'(A,I6,A,ES11.4)') ' PCG    iters=', it_p, '  ||b-Ax||/||b||=', relp
  write(*,'(A,ES11.4)')      ' max|x_minres - x_pcg| = ', sol_diff
  write(*,'(A,L1)')          ' PASS (both converged, solutions agree): ', &
       (relm < 1.0e-7_dp .and. relp < 1.0e-7_dp .and. sol_diff < 1.0e-6_dp)
  write(*,*)

  ! ================= TEST 2: INDEFINITE system (CG breaks) =================
  ! A = tridiagonal Laplacian (2,-1) minus 2*I  -> eigenvalues in (-2,2).
  a_diag = 2.0_dp; a_offd = 1.0_dp; op_shift = 2.0_dp; precond_identity = .true.
  inject_mode = 0

  xm = b
  call minres_optimize(xm, matvec, precond, dummy, mxit, tol=tol, err=errm, iters=it_m)
  relm = true_resid(b, xm)

  call solve_pcg(b, xp, it_p, relp)

  write(*,'(A)') '===== TEST 2: INDEFINITE system (Laplacian - 2I), M = I ====='
  write(*,'(A,I6,A,ES11.4,A,L1)') ' MINRES iters=', it_m, '  ||b-Ax||/||b||=', relm, &
       '   converged=', (relm < 1.0e-6_dp)
  write(*,'(A,I6,A,ES11.4,A,L1)') ' PCG    iters=', it_p, '  ||b-Ax||/||b||=', relp, &
       '   converged=', (relp < 1.0e-6_dp)
  write(*,'(A)') ' (CG on an indefinite operator is not guaranteed to converge and,'
  write(*,'(A)') '  when it does, needs far more iterations; MINRES is the safe choice.)'
  write(*,'(A,L1)') ' PASS (MINRES solves indefinite system robustly): ', (relm < 1.0e-6_dp)
  write(*,*)

  ! ========== TEST 3: textbook CG breakdown, A = diag(1,-1) ================
  ! p0^T A p0 = 0 on the very first step -> CG divides by ~0; MINRES is fine.
  block
    real(dp) :: b2(2), x2(2), xm2(2), rel2m, rel2p
    integer  :: itm2, itp2
    indef_diag = .true.; precond_identity = .true.; inject_mode = 0
    b2 = [1.0_dp, 1.0_dp]
    xm2 = b2
    call minres_optimize(xm2, matvec, precond, dummy, 100, tol=1.0e-12_dp, iters=itm2)
    rel2m = true_resid(b2, xm2)
    call solve_pcg(b2, x2, itp2, rel2p)
    write(*,'(A)') '===== TEST 3: A = diag(+1,-1), b = (1,1) -> CG breaks at step 1 ====='
    write(*,'(A,I4,A,ES11.4,A,L1)') ' MINRES iters=', itm2, '  ||b-Ax||/||b||=', rel2m, &
         '   solved=', (rel2m < 1.0e-8_dp)
    write(*,'(A,I4,A,ES11.4,A,L1)') ' PCG    iters=', itp2, '  ||b-Ax||/||b||=', rel2p, &
         '   solved=', (rel2p < 1.0e-8_dp)
    write(*,'(A,L1)') ' PASS (MINRES solves it, CG does NOT): ', &
         (rel2m < 1.0e-8_dp .and. .not. (rel2p < 1.0e-8_dp))
  end block

contains

  real(dp) function true_resid(b, x)
    real(dp), intent(in) :: b(:), x(:)
    real(dp) :: ax(size(b)), r(size(b))
    integer :: save_mode
    save_mode = inject_mode; inject_mode = 0
    call matvec(ax, x, c_null_ptr)
    inject_mode = save_mode
    r = b - ax
    true_resid = sqrt(sum(r*r)) / sqrt(sum(b*b))
  end function true_resid

  ! Drive PCG with init/step so we can read iters and the true residual even
  ! when it fails to converge on the indefinite system.
  subroutine solve_pcg(b, x, iters, rel)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: iters
    real(dp), intent(out) :: rel
    type(pcg_t) :: s
    integer :: kk
    call s%init(b=b, update=matvec, precond=precond, dat=dummy, tol=tol)
    iters = 0
    if (s%errcode == PCG_OK) then
      do kk = 1, mxit
        call s%step()
        if (s%errcode /= PCG_OK) exit
        iters = iters + 1
      end do
    end if
    x = s%x
    call s%clean()
    rel = true_resid(b, x)
  end subroutine solve_pcg

end program minres_check
