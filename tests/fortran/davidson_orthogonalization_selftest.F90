module davidson_orthogonalization_selftest_m
  use iso_c_binding, only: c_double, c_int
  use precision, only: dp
  use tdhf_lib, only: rpanewb
  use tdhf_sf_lib, only: sfresvec, sfqvec
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
  implicit none
  private
  public :: davidson_orthogonalization_selftest
contains
  subroutine davidson_orthogonalization_selftest(error, checks, failures) bind(C)
    real(c_double), intent(out) :: error
    integer(c_int), intent(out) :: checks, failures
    real(dp), parameter :: scales(4) = [1.0e-300_dp, 1.0e-12_dp, 1.0_dp, 1.0e300_dp]
    real(dp) :: basis(6,6), candidates(6,4), expected(6,6), gram(6,6), local_error
    integer :: mode, scale_index, kind, nvec, novec, status, wanted, j
    logical :: tda, accepted

    error=0; checks=0; failures=0
    expected=0
    expected(1,1)=1
    expected(2:3,2)=[2.0_dp,3.0_dp]/sqrt(13.0_dp)
    expected(2:4,3)=[3.0_dp,-2.0_dp,5.0_dp]/sqrt(38.0_dp)
    expected(5,4)=1; expected(6,5)=1
    do mode=1,2
      tda=mode==1
      wanted=merge(3,5,tda)
      do scale_index=1,size(scales)
        basis=0; basis(1,1)=1
        candidates=0
        candidates(1:3,1)=[1.0_dp,2.0_dp,3.0_dp]
        candidates(1:4,2)=[-4.0_dp,3.0_dp,-2.0_dp,5.0_dp]
        candidates(1,3)=7; candidates(5,3)=2
        candidates(1,4)=-3; candidates(6,4)=4
        candidates=candidates*scales(scale_index)
        nvec=1; novec=-1
        call rpanewb(2,basis,candidates,novec,nvec,status,tda)
        checks=checks+1
        if (status/=0 .or. nvec/=wanted .or. novec/=1) failures=failures+1
        if (.not.all(ieee_is_finite(basis))) then
          failures=failures+1
        else
          local_error=maxval(abs(basis(:,1:wanted)-expected(:,1:wanted)))
          gram=matmul(transpose(basis),basis)
          do j=1,wanted
            gram(j,j)=gram(j,j)-1
          end do
          local_error=max(local_error,maxval(abs(gram)))
          error=max(error,local_error)
          if (local_error>1.0e-12_dp) failures=failures+1
        end if
      end do
      do scale_index=1,size(scales)
        do kind=0,5
          basis=0; basis(1,1)=1; candidates=0
          select case(kind)
          case(0) ! Exactly zero update.
          case(1)
            candidates(2,1)=ieee_value(0.0_dp,ieee_quiet_nan)
          case(2)
            candidates(2,1)=ieee_value(0.0_dp,ieee_positive_inf)
          case(3) ! Exactly dependent, at every amplitude.
            candidates(1,1)=scales(scale_index)
          case(4) ! Relative remainder below the dependence cutoff.
            candidates(1:2,1)=scales(scale_index)*[1.0_dp,1.0e-12_dp]
          case(5) ! Small but resolvable independent direction.
            candidates(1:2,1)=scales(scale_index)*[1.0_dp,1.0e-6_dp]
          end select
          candidates(:,2)=candidates(:,1)
          nvec=1; novec=-1
          call rpanewb(1,basis,candidates,novec,nvec,status,tda)
          accepted=kind==5
          checks=checks+1
          if (nvec/=merge(2,1,accepted) .or. status/=merge(0,3,accepted) .or. novec/=1) &
            failures=failures+1
          if (.not.all(ieee_is_finite(basis))) then
            failures=failures+1
          else
            ! Original basis is always retained; accepted near-dependent
            ! candidates have the unique orthogonal direction +e2.
            local_error=abs(basis(1,1)-1)+sum(abs(basis(2:,1)))
            if (accepted) then
              local_error=max(local_error,abs(basis(2,2)-1),abs(basis(1,2)),maxval(abs(basis(3:,2))))
            end if
            error=max(error,local_error)
            if (local_error>1.0e-12_dp) failures=failures+1
          end if
        end do
      end do
    end do
    call sf_converged_root_cases(error, checks, failures)
  end subroutine davidson_orthogonalization_selftest

  ! SF/MRSF chain sfresvec -> sfqvec -> rpanewb with two tracked roots in the
  ! subspace span{e1,e2} of a signed 5x5 model matrix.  Root 1 couples to e4
  ! only through c1; root 2 is unconverged through couplings to e3 and e5.
  ! The drivers mark roots with ||r||^2 <= cnvtol as optional candidates.
  subroutine sf_converged_root_cases(error, checks, failures)
    real(c_double), intent(inout) :: error
    integer(c_int), intent(inout) :: checks, failures
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: q2(5), e4(5)

    ! Root 2's preconditioned correction (-0.3 e3)/(2-3) + (0.2 e5)/(2-5) and
    ! root 1's c1 e4/(1-4), both normalized.
    q2 = 0.0_dp; q2(3) = 0.3_dp; q2(5) = -0.2_dp/3.0_dp
    q2 = q2/norm2(q2)
    e4 = 0.0_dp; e4(4) = -1.0_dp

    ! Roundoff-sized residual on root 1, one free slot: the unconverged root
    ! gets the slot.
    call run_case(1.0e-14_dp, 3, .true., 0, 3, q2, q2)
    ! Same data without the mask: sfresvec returns a zero correction for a
    ! residual at its own rounding floor, so the noise no longer takes the
    ! slot (before that floor it did, and the solver stopped with ick=2).
    call run_case(1.0e-14_dp, 3, .false., 0, 3, q2, q2)
    ! Roundoff-sized residual on root 1 with two free slots: its direction is
    ! rounding noise and must not enter the subspace (before the floor e4 was
    ! appended as a full unit vector, making the subspace depend on the
    ! summation order of the sigma build).
    call run_case(1.0e-14_dp, 4, .true., 0, 3, q2, q2)
    ! ||r1||^2 = 2.5e-11 <= tol with two free slots: the converged root still
    ! refines the subspace, after the required candidate.
    call run_case(5.0e-6_dp, 4, .true., 0, 4, q2, e4)
    ! ||r1||^2 = 4e-10 > tol: both roots are required; lack of space for one
    ! of them is still reported as ick=2.
    call run_case(2.0e-5_dp, 3, .true., 2, 3, e4, e4)

  contains

    subroutine run_case(c1, mxvec, use_mask, wanted_status, wanted_nvec, col3, col4)
      real(dp), intent(in) :: c1, col3(5), col4(5)
      integer, intent(in) :: mxvec, wanted_status, wanted_nvec
      logical, intent(in) :: use_mask
      real(dp) :: amat(5,5), basis(5,4), sigma(5,2), vec(2,2), eigv(2), xm(5), &
                  q(5,2), rnorm(2), local_error
      integer :: i, nvec, novec, status

      amat = 0.0_dp
      do i = 1, 5
        amat(i,i) = real(i, dp)
        xm(i) = real(i, dp)
      end do
      amat(2,3) = -0.3_dp; amat(3,2) = -0.3_dp
      amat(2,5) =  0.2_dp; amat(5,2) =  0.2_dp
      amat(1,4) = c1; amat(4,1) = c1
      basis = 0.0_dp; basis(1,1) = 1.0_dp; basis(2,2) = 1.0_dp
      sigma = matmul(amat, basis(:,1:2))
      vec = 0.0_dp; vec(1,1) = 1.0_dp; vec(2,2) = 1.0_dp
      eigv = [1.0_dp, 2.0_dp]
      call sfresvec(q, basis(:,1:mxvec), sigma, vec, eigv, 2, rnorm, 2)
      call sfqvec(q, xm, eigv, 2)
      nvec = 2; novec = -1
      if (use_mask) then
        call rpanewb(2, basis(:,1:mxvec), q, novec, nvec, status, .true., &
                     optional_candidate=rnorm<=tol)
      else
        call rpanewb(2, basis(:,1:mxvec), q, novec, nvec, status, .true.)
      end if

      checks = checks+1
      if (status /= wanted_status .or. nvec /= wanted_nvec .or. novec /= 2) &
        failures = failures+1
      ! The floor zeroes only the correction: rnorm stays the computed
      ! ||r||^2 = c1^2 of root 1 and 0.13 of root 2.
      if (abs(rnorm(1) - c1**2) > 1.0e-12_dp*c1**2 .or. abs(rnorm(2) - 0.13_dp) > 1.0e-14_dp) &
        failures = failures+1
      ! Original basis retained, appended columns are the expected unique
      ! directions with their signs.
      local_error = abs(basis(1,1)-1) + abs(basis(2,2)-1) &
                  + sum(abs(basis(2:,1))) + abs(basis(1,2)) + sum(abs(basis(3:,2)))
      local_error = max(local_error, maxval(abs(basis(:,3)-col3)))
      if (wanted_nvec == 4) local_error = max(local_error, maxval(abs(basis(:,4)-col4)))
      error = max(error, local_error)
      if (local_error > 1.0e-12_dp) failures = failures+1
    end subroutine run_case
  end subroutine sf_converged_root_cases
end module davidson_orthogonalization_selftest_m
