! Pin the MRSF/UMRSF <S^2> of response vectors.
!
! ROHF limit (alpha and beta MOs identical): the UMRSF expression umrsfssqu must
! reduce to the ROHF result mrsfssq, singlet 2*w_CV and triplet 2, for random
! vectors that populate every configuration class, including the G/D open-open
! slots whose normalization depends on the alpha/beta overlaps.
! Unequal alpha/beta MOs: a vector with only core->virtual amplitudes must have
! the same <S^2> for both multiplicities, because E+ and E- of a C->V pair differ
! by two alpha<->beta exchanges and S^2 couples them only through one.
module mrsf_spin_square_selftest_m
  use iso_c_binding, only: c_double, c_int
  use precision, only: dp
  use tdhf_mrsf_lib, only: umrsfssqu, mrsfssq
  implicit none
  public :: mrsf_spin_square_selftest
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

      ! Unequal alpha/beta MOs: C->V-only vector, singlet == triplet
      mo_b = mo_a
      do i = 1, n-1
        th = 0.1_dp*real(i+shape, dp)
        mo_b(:,[i,i+1]) = matmul(mo_b(:,[i,i+1]), &
          reshape([cos(th), sin(th), -sin(th), cos(th)], [2,2]))
      end do
      x = 0
      do j = 3, nv
        do i = 1, nb
          x(i,j,1) = cos(real(7*i+3*j+shape, dp))
        end do
      end do
      x = x/sqrt(sum(x**2))
      call umrsfssqu(ss_u, mo_a, mo_b, smat, wrk1, wrk2, n, n*(n+1)/2, na*nv, 1, &
                     n, x, na, nb, .true., .false.)
      call umrsfssqu(ss_r, mo_a, mo_b, smat, wrk1, wrk2, n, n*(n+1)/2, na*nv, 1, &
                     n, x, na, nb, .false., .true.)
      ! the alpha/beta overlaps must be non-trivial (a single core orbital
      ! gives exactly 2 for any overlaps, so test the overlaps themselves)
      if (minval([(abs(dot_product(mo_a(:,k), mo_b(:,k))), k = 1, na)]) > 0.99_dp) &
        error = max(error, 1.0_dp)
      error = max(error, abs(ss_u-ss_r))
      checks = checks+1

      deallocate(mo_a, mo_b, smat, wrk1, wrk2, x)
    end do
  end subroutine mrsf_spin_square_selftest
end module mrsf_spin_square_selftest_m
