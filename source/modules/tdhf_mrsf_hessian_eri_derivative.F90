module tdhf_mrsf_hessian_eri_derivative_mod

  use precision, only: dp
  use types, only: information
  use basis_tools, only: basis_set
  use fock_deriv_mod, only: fock_deriv_matrix_general_scaled
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
    real(kind=dp), allocatable :: matrix4(:,:,:,:)
    real(kind=dp) :: channel_scale,coulomb_scale,spin_multiplier
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
    allocate(matrix4(nbf,nbf,3,natom))
    do channel=1,7
      select case(channel)
      case(1:4)
        channel_scale=spc_coov
        coulomb_scale=channel_scale
      case(5)
        channel_scale=spc_ovov
        coulomb_scale=0.0_dp
      case(6)
        channel_scale=spc_coco
        coulomb_scale=0.0_dp
      case(7)
        channel_scale=response_scale
        coulomb_scale=0.0_dp
      end select
      call fock_deriv_matrix_general_scaled(infos,basis,density(channel,:,:), &
        coulomb_scale,channel_scale,matrix4)
      if(channel<=6) matrix4=spin_multiplier*matrix4
      do atom=1,natom
        do cart=1,3
          derivative_fock(channel,:,:,3*(atom-1)+cart)=matrix4(:,:,cart,atom)
        end do
      end do
    end do
    deallocate(matrix4)
  end subroutine build_mrsf_explicit_eri_derivative

end module tdhf_mrsf_hessian_eri_derivative_mod
