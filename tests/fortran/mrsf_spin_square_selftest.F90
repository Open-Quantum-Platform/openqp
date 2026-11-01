! Pin the MRSF/UMRSF <S^2> of response vectors.
!
! ROHF limit (alpha and beta MOs identical): the UMRSF expression umrsfssqu must
! reduce to the ROHF result mrsfssq, singlet 2*w_CV and triplet 2, for random
! vectors that populate every configuration class, including the G/D open-open
! slots whose normalization depends on the alpha/beta overlaps.
module mrsf_spin_square_selftest_m
  use iso_c_binding, only: c_double, c_int
  use precision, only: dp
  use tdhf_mrsf_lib, only: umrsfssqu, umrsfssq_quintet, mrsfssq
  implicit none
  public :: mrsf_spin_square_selftest, mrsf_spin_square_evaluate
contains
  subroutine mrsf_spin_square_selftest(error, checks) bind(C)
    real(c_double), intent(out) :: error
    integer(c_int), intent(out) :: checks
    integer, parameter :: sizes(3) = [7, 9, 6], occ(3) = [4, 5, 3]
    real(dp), allocatable :: mo_a(:,:), mo_b(:,:), smat(:), wrk1(:,:), wrk2(:), x(:,:,:)
    real(dp) :: ss_u, ss_r, wcv, th
    integer :: shape, n, na, nb, nv, i, j, k, mult
    logical :: singlet

    error = 0; checks = 0
    do shape = 1, 3
      n = sizes(shape); na = occ(shape); nb = na-2; nv = n-nb
      allocate(mo_a(n,n), mo_b(n,n), smat(n*(n+1)/2), wrk1(n,n), wrk2(n), x(na,nv,1))
      smat = 0
      k = 0
      do j = 1, n
        do i = 1, j
          k = k+1
          if (i == j) smat(k) = 1
        end do
      end do
      mo_a = 0
      do i = 1, n
        mo_a(i,i) = 1
      end do

      ! ROHF limit, every configuration class populated
      mo_b = mo_a
      do mult = 1, 3, 2
        singlet = mult == 1
        do j = 1, nv
          do i = 1, na
            x(i,j,1) = sin(real(5*i+11*j+3*shape+mult, dp))
          end do
        end do
        x(nb+2,2,1) = 0                       ! O2->O2 slot is not an MRSF amplitude
        if (.not. singlet) then
          x(nb+2,1,1) = 0; x(nb+1,2,1) = 0    ! G and D exist only for singlets
        end if
        x = x/sqrt(sum(x**2))
        wcv = sum(x(1:nb,3:nv,1)**2)
        call umrsfssqu(ss_u, mo_a, mo_b, smat, wrk1, wrk2, n, n*(n+1)/2, na*nv, 1, &
                       n, x, na, nb, singlet, .not. singlet)
        call mrsfssq(ss_r, reshape(x, [na*nv]), na, nb, n, singlet)
        if (singlet) then
          error = max(error, abs(ss_r-2*wcv), abs(ss_u-2*wcv))
        else
          error = max(error, abs(ss_r-2), abs(ss_u-2))
        end if
        checks = checks+2
      end do

      ! A beta core/added-virtual rotation leaves their occupied wedge unchanged.
      ! Both target parities must give exactly 2, not merely agree with each other.
      mo_b = mo_a
      th = 0.6_dp
      mo_b(:,[nb,n]) = matmul(mo_b(:,[nb,n]), &
        reshape([cos(th),sin(th),-sin(th),cos(th)],[2,2]))
      x = 0; x(1,nv,1) = 1
      call umrsfssqu(ss_u, mo_a, mo_b, smat, wrk1, wrk2, n, n*(n+1)/2, na*nv, 1, &
                     n, x, na, nb, .true., .false.)
      call umrsfssqu(ss_r, mo_a, mo_b, smat, wrk1, wrk2, n, n*(n+1)/2, na*nv, 1, &
                     n, x, na, nb, .false., .true.)
      error = max(error, abs(ss_u-2), abs(ss_r-2))
      checks = checks+1

      deallocate(mo_a, mo_b, smat, wrk1, wrk2, x)
    end do
  end subroutine mrsf_spin_square_selftest
  ! Test interface to the same AO-overlap/MO/packed-vector route used after Davidson.
  subroutine mrsf_spin_square_evaluate(n_c, nc_c, mult_c, mo_a, mo_b, smat, x, ss) bind(C)
    integer(c_int), value, intent(in) :: n_c, nc_c, mult_c
    real(c_double), intent(in) :: mo_a(n_c,n_c), mo_b(n_c,n_c), smat(n_c*(n_c+1)/2)
    real(c_double), intent(in) :: x(*)
    real(c_double), intent(out) :: ss
    real(dp), allocatable :: wrk1(:,:), wrk2(:)
    integer :: n, nc
    n = int(n_c); nc = int(nc_c)
    allocate(wrk1(n,n),wrk2(n))
    if (mult_c == 5) then
      call umrsfssq_quintet(ss,mo_a,mo_b,smat,wrk1,n,n,x,nc+2,nc,1)
    else
      call umrsfssqu(ss,mo_a,mo_b,smat,wrk1,wrk2,n,n*(n+1)/2,(nc+2)*(n-nc),1, &
                    n,x,nc+2,nc,mult_c==1,mult_c==3)
    end if
    deallocate(wrk1,wrk2)
  end subroutine mrsf_spin_square_evaluate
end module mrsf_spin_square_selftest_m
