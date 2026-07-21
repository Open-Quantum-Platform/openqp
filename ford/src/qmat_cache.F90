!> @brief Cached canonical orthogonalizer Q = S^(-1/2)
!
!> @details Q is needed by every initial guess and by the SCF setup;
!>          without caching it is computed (an O(nbf^3) eigendecomposition)
!>          at least twice per run. The cache lives in the tagarray under
!>          OQP_QMAT and is invalidated by int1e whenever the overlap
!>          matrix is recomputed (new geometry or basis set).
module qmat_cache

  use precision, only: dp

  implicit none

  private
  public get_qmat_cached

contains

!> @brief Return Q = S^(-1/2), reusing the cached copy if available
!
!> @param[in]  infos  OQP run information (holds the tagarray)
!> @param[in]  smat   packed overlap matrix
!> @param[out] qmat   canonical orthogonalizer, (nbf x nbf);
!>                    columns beyond the rank are zero
!> @param[in]  nbf    number of basis functions
!> @param[out] qrnk   optional, rank of Q (number of linearly
!>                    independent basis functions)
 subroutine get_qmat_cached(infos, smat, qmat, nbf, qrnk)
    use types, only: information
    use oqp_tagarray_driver
    use mathlib, only: matrix_invsqrt
    use, intrinsic :: iso_c_binding, only: c_int32_t

    implicit none

    type(information), intent(inout) :: infos
    real(kind=dp), intent(in) :: smat(*)
    real(kind=dp), intent(out) :: qmat(nbf,*)
    integer, intent(in) :: nbf
    integer, intent(out), optional :: qrnk

    character(len=*), parameter :: tags_qmat(1) = (/ character(len=80) :: OQP_QMAT /)
    real(kind=dp), contiguous, pointer :: q_st(:,:)
    integer(c_int32_t) :: tag_id
    integer :: j

    if (infos%dat%contains(tags_qmat, tag_id)) then
      call tagarray_get_data(infos%dat, OQP_QMAT, q_st)
      if (size(q_st,1) == nbf .and. size(q_st,2) == nbf) then
        qmat(:,1:nbf) = q_st
        if (present(qrnk)) then
!         matrix_invsqrt zeroes the columns beyond the rank
          do j = nbf, 1, -1
            if (any(qmat(:,j) /= 0.0_dp)) exit
          end do
          qrnk = j
        end if
        return
      end if
!     dimension mismatch (e.g. stale record): recompute below
    end if

    call matrix_invsqrt(smat, qmat, nbf, qrnk=qrnk)

    call infos%dat%alloc_or_die(OQP_QMAT, (/ nbf, nbf /), q_st, description=OQP_QMAT_comment)
    q_st = qmat(:,1:nbf)

 end subroutine get_qmat_cached

end module qmat_cache
