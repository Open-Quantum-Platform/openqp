module tdhf_mrsf_hessian_operator_mod

  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_compute_t
  use tdhf_mrsf_hessian_space_mod, only: mrsf_physical_dimensions, &
    build_mrsf_packed_transform
  use tdhf_mrsf_sigma_mod, only: apply_mrsf_tda_sigma

  implicit none

  private
  public :: build_mrsf_tda_matrix

contains

!###############################################################################

  subroutine build_mrsf_tda_matrix(infos,int2_driver,mo_a,mo_b,fa,fb, &
                                    matrix,asymmetry,status,column_batch)
    ! Dense diagnostic construction of the physical spin-adapted MRSF-TDA
    ! matrix.  Each column is obtained from the same seven-density production
    ! sigma action used by Davidson.  The routine is intended for small-system
    ! columnwise verification and the dense reference Hessian response, not as
    ! the production eigensolver.

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),fa(:,:),fb(:,:)
    real(kind=dp), intent(out) :: matrix(:,:),asymmetry
    integer, intent(out) :: status
    integer, intent(in), optional :: column_batch

    real(kind=dp), allocatable :: packed_basis(:,:),packed_sigma(:,:)
    integer :: alloc_status,batch,first,last,ncol
    integer :: nbf,nocca,noccb,packed,physical

    status = 0
    asymmetry = 0.0_dp
    matrix = 0.0_dp
    nbf = size(mo_a,1)
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    if (infos%basis%nbf /= nbf .or. &
        any(shape(mo_a) /= [nbf,nbf]) .or. &
        any(shape(mo_b) /= [nbf,nbf]) .or. &
        any(shape(fa) /= [nbf,nbf]) .or. &
        any(shape(fb) /= [nbf,nbf])) then
      status = -13
      return
    end if
    call mrsf_physical_dimensions(nbf,nocca,noccb,infos%tddft%mult, &
                                   packed,physical,status)
    if (status /= 0) return
    if (any(shape(matrix) /= [physical,physical])) then
      status = -10
      return
    end if
    ! A dense response matrix is intrinsically O(physical**2), but the integral
    ! consumer also replicates every simultaneous sigma column over its worker
    ! threads.  Default to one column at a time so the contraction workspace is
    ! bounded independently of the response-space dimension.  Small diagnostic
    ! callers may opt into a larger explicit batch after sizing their memory.
    batch = 1
    if (present(column_batch)) batch = column_batch
    if (batch <= 0) then
      status = -11
      return
    end if
    batch = min(batch,physical)

    allocate(packed_basis(packed,physical),stat=alloc_status)
    if (alloc_status /= 0) then
      status = -12
      return
    end if
    allocate(packed_sigma(packed,batch),stat=alloc_status)
    if (alloc_status /= 0) then
      deallocate(packed_basis)
      status = -12
      return
    end if
    call build_mrsf_packed_transform(nbf,nocca,noccb,infos%tddft%mult, &
                                     packed_basis,status)
    if (status /= 0) then
      deallocate(packed_basis,packed_sigma)
      return
    end if
    do first=1,physical,batch
      last = min(first+batch-1,physical)
      ncol = last-first+1
      call apply_mrsf_tda_sigma(infos,int2_driver,mo_a,mo_b,fa,fb, &
                                packed_basis(:,first:last), &
                                packed_sigma(:,1:ncol),status)
      if (status /= 0) exit
      matrix(:,first:last) = matmul(transpose(packed_basis), &
                                    packed_sigma(:,1:ncol))
    end do
    if (status == 0) asymmetry = maxval(abs(matrix-transpose(matrix)))
    deallocate(packed_basis,packed_sigma)
  end subroutine build_mrsf_tda_matrix

end module tdhf_mrsf_hessian_operator_mod
