module tdhf_mrsf_density_contraction_mod

  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_compute_t
  use tdhf_mrsf_lib, only: int2_mrsf_data_t
  use tdhf_mrsf_conventions_mod, only: mrsf_raw_spc_multiplier

  implicit none

  private
  public :: contract_mrsf_seven_density_batch

contains

!###############################################################################

  subroutine contract_mrsf_seven_density_batch(infos,int2_driver,density, &
      response_scale,spc_coco,spc_ovov,spc_coov,contracted,status)
    ! Contract the seven spin-adapted AO response densities with independent
    ! channel scales.  Four integral batches avoid every spc/HFscale division:
    ! CO--OV (1:4), OV--OV (5), CO--CO (6), and ordinary A0 (7).  This remains
    ! finite for pure semilocal DFT with explicitly requested SPC coefficients.

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), target, intent(in) :: density(:,:,:,:)
    real(kind=dp), intent(in) :: response_scale,spc_coco,spc_ovov,spc_coov
    real(kind=dp), intent(out) :: contracted(:,:,:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable, target :: group_density(:,:,:,:)
    real(kind=dp), allocatable :: envelope(:,:,:)
    integer :: allocation_status,nbf,nvec,spin_multiplier

    status=0
    contracted=0.0_dp
    nvec=size(density,1)
    nbf=infos%basis%nbf
    spin_multiplier=mrsf_raw_spc_multiplier(infos%tddft%mult)
    if(nvec<=0 .or. nbf<=0 .or. spin_multiplier==0 .or. &
       any(shape(density)/=[nvec,7,nbf,nbf]) .or. &
       any(shape(contracted)/=[nvec,7,nbf,nbf])) then
      status=-1
      return
    end if
    allocate(group_density(nvec,7,nbf,nbf),envelope(nvec,nbf,nbf), &
      source=0.0_dp,stat=allocation_status)
    if(allocation_status/=0) then
      status=-2
      return
    end if

    call run_group(1,4,spc_coov)
    if(status==0) call run_group(5,5,spc_ovov)
    if(status==0) call run_group(6,6,spc_coco)
    if(status==0) call run_group(7,7,response_scale)
    if(status==0 .and. spin_multiplier<0) &
      contracted(:,1:6,:,:)=-contracted(:,1:6,:,:)
    deallocate(group_density,envelope)

  contains

    subroutine run_group(first_channel,last_channel,scale)
      integer, intent(in) :: first_channel,last_channel
      real(kind=dp), intent(in) :: scale
      type(int2_mrsf_data_t), allocatable, target :: data
      real(kind=dp), pointer :: result(:,:,:,:)
      integer :: channel,local_status

      if(status/=0 .or. abs(scale)<=epsilon(1.0_dp)) return
      group_density=0.0_dp
      group_density(:,first_channel:last_channel,:,:)= &
        density(:,first_channel:last_channel,:,:)
      if(last_channel<7) then
        envelope=0.0_dp
        do channel=first_channel,last_channel
          envelope=max(envelope,abs(density(:,channel,:,:)))
        end do
        ! int2_mrsf_data_t screens shell pairs using its final density channel.
        ! A positive envelope is used only for screening; its channel-7 output
        ! is discarded below and cannot enter the requested group.
        group_density(:,7,:,:)=envelope
      end if
      allocate(data,source=int2_mrsf_data_t(d3=group_density, &
        tamm_dancoff=.true.,scale_exchange=scale,scale_coulomb=scale), &
        stat=local_status)
      if(local_status/=0) then
        status=-3
        return
      end if
      call int2_driver%run(data,cam=infos%control%hamilton==20 .and. &
        infos%dft%cam_flag,alpha=infos%tddft%cam_alpha, &
        alpha_coulomb=infos%tddft%cam_alpha,beta=infos%tddft%cam_beta, &
        beta_coulomb=infos%tddft%cam_beta,mu=infos%tddft%cam_mu)
      if (.not.allocated(data%f3)) then
        nullify(data%d3)
        deallocate(data)
        status=-4
        return
      end if
      if (size(data%f3,1)/=nvec .or. size(data%f3,2)/=7 .or. &
          size(data%f3,3)/=nbf .or. size(data%f3,4)/=nbf .or. &
          size(data%f3,5)<1) then
        call data%clean()
        deallocate(data)
        status=-4
        return
      end if
      result=>data%f3(:,:,:,:,1)
      contracted(:,first_channel:last_channel,:,:)= &
        result(:,first_channel:last_channel,:,:)
      nullify(result)
      call data%clean()
      deallocate(data)
    end subroutine run_group

  end subroutine contract_mrsf_seven_density_batch

end module tdhf_mrsf_density_contraction_mod
