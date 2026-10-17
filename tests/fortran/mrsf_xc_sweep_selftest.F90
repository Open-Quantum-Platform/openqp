module mrsf_xc_sweep_selftest_m
  use iso_c_binding, only: c_double
  use precision, only: dp
  use c_interop, only: oqp_handle_t, oqp_handle_get_info
  use types, only: information
  use basis_tools, only: basis_set
  use mod_dft, only: dft_initialize, dftclean
  use mod_dft_molgrid, only: dft_grid_t
  use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
  use mod_dft_gridint_response_cache, only: response_cache_t
  use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_DM_B, OQP_td_p
  use mathlib, only: unpack_matrix
  implicit none
  public :: mrsf_xc_sweep_selftest
contains
  ! Use the production relaxed density from the preceding MRSF gradient.
  subroutine mrsf_xc_sweep_selftest(c_handle,error,moving_norm) bind(C)
    type(oqp_handle_t) :: c_handle
    real(c_double), intent(out) :: error,moving_norm
    type(information), pointer :: infos
    type(basis_set) :: basis
    type(dft_grid_t) :: grid
    type(response_cache_t), target :: cache
    real(dp), contiguous, pointer :: da(:),db(:),tp(:,:)
    real(dp), allocatable :: d(:,:,:),p(:,:,:),d0(:,:,:),p0(:,:,:)
    real(dp), allocatable :: fixed(:,:),moving(:,:),fused(:,:)
    integer :: n,nat
    infos=>oqp_handle_get_info(c_handle)
    basis=infos%basis;basis%atoms=>infos%atoms
    n=basis%nbf;nat=infos%mol_prop%natom
    allocate(d(n,n,2),p(n,n,2),d0(n,n,2),p0(n,n,2),fixed(3,nat),moving(3,nat),fused(3,nat))
    call tagarray_get_data(infos%dat,OQP_DM_A,da)
    call tagarray_get_data(infos%dat,OQP_DM_B,db)
    call tagarray_get_data(infos%dat,OQP_td_p,tp)
    call unpack_matrix(da,d0(:,:,1));call unpack_matrix(db,d0(:,:,2))
    call unpack_matrix(tp(:,1),p0(:,:,1));call unpack_matrix(tp(:,2),p0(:,:,2))
    call dft_initialize(infos,basis,grid)
    d=d0;p=p0;fixed=0
    call utddft_xc_gradient(basis,grid,fixed,d(:,:,1),d(:,:,2),p(:,:,1:1),p(:,:,2:2), &
      nMtx=1,threshold=0.0_dp,infos=infos,include_ground_state=.false.,cache=cache)
    ! Match the original sequence, including its density normalization round trip.
    moving=0
    call utddft_xc_gradient(basis,grid,moving,d(:,:,1),d(:,:,2),p(:,:,1:1),p(:,:,2:2), &
      nMtx=1,threshold=0.0_dp,infos=infos,include_ground_state=.false., &
      include_weight_derivative=.true.,weight_derivative_only=.true.,cache=cache)
    d=d0;p=p0;fused=0.125_dp
    call utddft_xc_gradient(basis,grid,fused,d(:,:,1),d(:,:,2),p(:,:,1:1),p(:,:,2:2), &
      nMtx=1,threshold=0.0_dp,infos=infos,include_ground_state=.false., &
      include_weight_derivative=.true.,cache=cache)
    error=maxval(abs(fused-0.125_dp-fixed-moving))
    moving_norm=maxval(abs(moving))
    call cache%free()
    call dftclean(infos)
  end subroutine
end module
