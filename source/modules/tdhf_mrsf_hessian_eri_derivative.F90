module tdhf_mrsf_hessian_eri_derivative_mod

  use precision, only: dp
  use types, only: information
  use basis_tools, only: basis_set
  use fock_deriv_mod, only: fock_deriv_matrix_mrsf_scaled_batch
  use tdhf_mrsf_conventions_mod, only: mrsf_raw_spc_multiplier

  implicit none

  private
  public :: build_mrsf_explicit_eri_derivative

contains

!###############################################################################

  subroutine build_mrsf_explicit_eri_derivative(infos,density,response_scale, &
      spc_coco,spc_ovov,spc_coov,multiplicity,derivative_fock,status)
    ! Explicit nuclear derivative of the two-electron action on the seven
    ! spin-adapted MRSF AO densities.  Channels 1:4 contain Coulomb and
    ! exchange, channels 5:7 exchange only.  The special-channel scales are
    ! applied directly rather than through spc/HFscale, so pure semilocal DFT
    ! with a user-specified SPC coefficient is finite and well defined.

    type(information), target, intent(inout) :: infos
    real(kind=dp), target, intent(in) :: density(:,:,:)
    real(kind=dp), intent(in) :: response_scale,spc_coco,spc_ovov,spc_coov
    integer, intent(in) :: multiplicity
    real(kind=dp), intent(out) :: derivative_fock(:,:,:,:)
    integer, intent(out) :: status

    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: matrix5(:,:,:,:,:), &
      channel_scale(:),coulomb_scale(:)
    real(kind=dp) :: spin_multiplier
    integer :: atom,cart,channel,natom,nbf,ncoord

    basis=>infos%basis
    basis%atoms=>infos%atoms
    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    ncoord=3*natom
    status=0
    derivative_fock=0.0_dp
    spin_multiplier=real(mrsf_raw_spc_multiplier(multiplicity),dp)
    if(nbf<=0 .or. natom<=0 .or. spin_multiplier==0.0_dp .or. &
       any(shape(density)/=[7,nbf,nbf]) .or. &
       any(shape(derivative_fock)/=[7,nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    allocate(matrix5(nbf,nbf,3,natom,7),channel_scale(7), &
      coulomb_scale(7),source=0.0_dp)
    do channel=1,7
      select case(channel)
      case(1:4)
        channel_scale(channel)=spc_coov
        coulomb_scale(channel)=channel_scale(channel)
      case(5)
        channel_scale(channel)=spc_ovov
      case(6)
        channel_scale(channel)=spc_coco
      case(7)
        channel_scale(channel)=response_scale
      end select
    end do
    call fock_deriv_matrix_mrsf_scaled_batch(infos,basis, &
      transpose_density_channels(density),coulomb_scale,channel_scale,matrix5)
    do channel=1,7
      if(channel<=6) matrix5(:,:,:,:,channel)= &
        spin_multiplier*matrix5(:,:,:,:,channel)
      do atom=1,natom
        do cart=1,3
          derivative_fock(channel,:,:,3*(atom-1)+cart)= &
            matrix5(:,:,cart,atom,channel)
        end do
      end do
    end do
    deallocate(matrix5,channel_scale,coulomb_scale)

  contains

    function transpose_density_channels(channels) result(batch)
      real(kind=dp), intent(in) :: channels(:,:,:)
      real(kind=dp) :: batch(size(channels,2),size(channels,3),size(channels,1))
      integer :: channel_index
      do channel_index=1,size(channels,1)
        batch(:,:,channel_index)=channels(channel_index,:,:)
      end do
    end function transpose_density_channels
  end subroutine build_mrsf_explicit_eri_derivative

end module tdhf_mrsf_hessian_eri_derivative_mod
