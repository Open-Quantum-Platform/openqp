module tdhf_mrsf_sigma_mod

  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_compute_t
  use tdhf_lib, only: iatogen
  use tdhf_mrsf_lib, only: int2_mrsf_data_t, mrsfcbc, mrsfmntoia, mrsfesum
  use tdhf_mrsf_density_contraction_mod, only: &
    contract_mrsf_seven_density_batch
  use tdhf_mrsf_conventions_mod, only: mrsf_raw_spc_multiplier
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public :: apply_mrsf_tda_sigma
  public :: prepare_mrsf_response_scaling
  public :: mrsf_sigma_status_message

contains

!###############################################################################

  subroutine prepare_mrsf_response_scaling(infos,response_scale,status)
    ! Resolve the TDDFT response and SPC sentinels in one place.  Both the
    ! Davidson driver and standalone Hessian response must use this routine;
    ! otherwise a fresh information object can expose the -1 input sentinel as
    ! a physical exact-exchange scale.

    type(information), intent(inout) :: infos
    real(kind=dp), intent(out) :: response_scale
    integer, intent(out) :: status

    logical :: dft

    status = 0
    dft = infos%control%hamilton == 20

    if (dft) then
      if (infos%tddft%hfscale == -1.0_dp) &
        infos%tddft%hfscale = infos%dft%hfscale
      if (infos%dft%cam_flag) then
        if (infos%tddft%cam_alpha == -1.0_dp) &
          infos%tddft%cam_alpha = infos%dft%cam_alpha
        if (infos%tddft%cam_beta == -1.0_dp) &
          infos%tddft%cam_beta = infos%dft%cam_beta
        if (infos%tddft%cam_mu == -1.0_dp) &
          infos%tddft%cam_mu = infos%dft%cam_mu
        infos%tddft%hfscale = infos%tddft%cam_alpha
      end if
      response_scale = infos%tddft%hfscale
    else
      response_scale = 1.0_dp
      infos%tddft%hfscale = response_scale
    end if

    if (infos%tddft%spc_coco == -1.0_dp) &
      infos%tddft%spc_coco = response_scale
    if (infos%tddft%spc_ovov == -1.0_dp) &
      infos%tddft%spc_ovov = response_scale
    if (infos%tddft%spc_coov == -1.0_dp) &
      infos%tddft%spc_coov = response_scale

    if (.not.ieee_is_finite(response_scale) .or. &
        .not.ieee_is_finite(infos%tddft%spc_coco) .or. &
        .not.ieee_is_finite(infos%tddft%spc_ovov) .or. &
        .not.ieee_is_finite(infos%tddft%spc_coov)) then
      status = -4
      return
    end if
    if (dft .and. infos%dft%cam_flag) then
      if (.not.ieee_is_finite(infos%tddft%cam_alpha) .or. &
          .not.ieee_is_finite(infos%tddft%cam_beta) .or. &
          .not.ieee_is_finite(infos%tddft%cam_mu)) status = -4
    end if
  end subroutine prepare_mrsf_response_scaling

!###############################################################################

  pure function mrsf_sigma_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=160) :: message

    select case(status)
    case(-1)
      message = 'Native MRSF sigma supports only spin-adapted singlet/triplet MRSF.'
    case(-2)
      message = 'Native MRSF sigma received inconsistent basis, electron, or array dimensions.'
    case(-3)
      message = 'Zero-HFscale SPC fallback does not support a range-separated CAM response.'
    case(-4)
      message = 'Native MRSF sigma received a non-finite response, CAM, or SPC scale.'
    case(-5)
      message = 'Cannot allocate memory for the native MRSF sigma application.'
    case(-6)
      message = 'Native MRSF sigma requires an initialized two-electron integral driver.'
    case(-7)
      message = 'Native MRSF sigma integral contraction did not return seven channel images.'
    case(-8)
      message = 'Independent-channel MRSF density contraction failed.'
    case default
      message = 'Native MRSF sigma application failed.'
    end select
  end function mrsf_sigma_status_message

!###############################################################################

  subroutine apply_mrsf_tda_sigma(infos,int2_driver,mo_a,mo_b,fa,fb, &
                                   trial,sigma,status)
    ! Apply the production spin-adapted MRSF-TDA response operator to a batch
    ! of packed physical trial vectors.  This is the single native path shared
    ! by Davidson and the Hessian response:
    !
    !   packed X -> mrsfcbc seven AO densities -> integral contraction
    !            -> mrsfmntoia -> mrsfesum -> packed A X.
    !
    ! Channels 1:6 form K_raw=-C_SPC^physical and use the named spin
    ! conversion.  Channel 7 (ball) belongs to A0 and is never spin-flipped.

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),fa(:,:),fb(:,:)
    real(kind=dp), intent(in) :: trial(:,:)
    real(kind=dp), intent(out) :: sigma(:,:)
    integer, intent(out) :: status

    type(int2_mrsf_data_t), target :: int2_data
    real(kind=dp), allocatable, target :: density(:,:,:,:)
    real(kind=dp), allocatable, target :: grouped_contracted(:,:,:,:)
    real(kind=dp), pointer :: contracted(:,:,:,:)
    real(kind=dp), allocatable :: work(:,:)
    real(kind=dp) :: response_scale
    integer :: alloc_status,ivec,mrst,nbf,nocca,noccb,nvec,xvec_dim
    logical :: dft,grouped_channels,nonzero_spc

    status = 0
    sigma = 0.0_dp
    nbf = size(mo_a,1)
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    mrst = infos%tddft%mult
    nvec = size(trial,2)
    xvec_dim = nocca*(nbf-noccb)
    dft = infos%control%hamilton == 20

    if (infos%tddft%umrsf .or. (mrst /= 1 .and. mrst /= 3)) then
      status = -1
      return
    end if
    if (nbf <= 0 .or. infos%basis%nbf /= nbf .or. &
        noccb < 0 .or. nocca < 2 .or. nocca > nbf .or. &
        nocca-noccb /= 2 .or. nvec <= 0 .or. &
        any(shape(mo_a) /= [nbf,nbf]) .or. &
        any(shape(mo_b) /= [nbf,nbf]) .or. &
        any(shape(fa) /= [nbf,nbf]) .or. &
        any(shape(fb) /= [nbf,nbf]) .or. &
        size(trial,1) /= xvec_dim .or. &
        any(shape(sigma) /= [xvec_dim,nvec])) then
      status = -2
      return
    end if

    if (.not.associated(int2_driver%basis)) then
      status = -6
      return
    end if
    if (int2_driver%basis%nbf /= nbf .or. &
        .not.allocated(int2_driver%schwarz_ints_regular)) then
      status = -6
      return
    end if

    call prepare_mrsf_response_scaling(infos,response_scale,status)
    if (status /= 0) return
    nonzero_spc = abs(infos%tddft%spc_coco) > epsilon(1.0_dp) .or. &
                  abs(infos%tddft%spc_ovov) > epsilon(1.0_dp) .or. &
                  abs(infos%tddft%spc_coov) > epsilon(1.0_dp)
    grouped_channels = abs(response_scale) <= epsilon(1.0_dp) .and. nonzero_spc
    if (grouped_channels .and. dft .and. infos%dft%cam_flag) then
      ! In a CAM run int2_run_cam replaces the consumer scales with alpha/beta,
      ! and its second pass updates only channel 7.  The four-group full-range
      ! fallback is therefore valid for pure semilocal DFT, not zero-alpha CAM.
      status = -3
      return
    end if

    allocate(density(nvec,7,nbf,nbf),source=0.0_dp,stat=alloc_status)
    if (alloc_status /= 0) then
      status = -5
      return
    end if
    allocate(work(nbf,nbf),source=0.0_dp,stat=alloc_status)
    if (alloc_status /= 0) then
      deallocate(density)
      status = -5
      return
    end if
    do ivec=1,nvec
      call iatogen(trial(:,ivec),work,nocca,noccb)
      call mrsfcbc(infos,mo_a,mo_b,work,density(ivec,:,:,:))
    end do

    if (grouped_channels) then
      allocate(grouped_contracted(nvec,7,nbf,nbf),source=0.0_dp, &
               stat=alloc_status)
      if (alloc_status /= 0) then
        deallocate(density,work)
        status = -5
        return
      end if
      call contract_mrsf_seven_density_batch( &
        infos,int2_driver,density,response_scale, &
        infos%tddft%spc_coco,infos%tddft%spc_ovov, &
        infos%tddft%spc_coov,grouped_contracted,status)
      if (status /= 0) then
        deallocate(grouped_contracted,density,work)
        status = -8
        return
      end if
      contracted => grouped_contracted
    else
      int2_data = int2_mrsf_data_t( &
        d3=density,tamm_dancoff=.true.,scale_exchange=response_scale, &
        scale_coulomb=response_scale)
      call int2_driver%run(int2_data, &
        cam=dft.and.infos%dft%cam_flag, &
        alpha=infos%tddft%cam_alpha, &
        alpha_coulomb=infos%tddft%cam_alpha, &
        beta=infos%tddft%cam_beta, &
        beta_coulomb=infos%tddft%cam_beta, &
        mu=infos%tddft%cam_mu)
      if (.not.allocated(int2_data%f3)) then
        nullify(int2_data%d3)
        deallocate(density,work)
        status = -7
        return
      end if
      if (size(int2_data%f3,2) /= 7) then
        call int2_data%clean()
        deallocate(density,work)
        status = -7
        return
      end if
      contracted => int2_data%f3(:,:,:,:,1)

      if (mrsf_raw_spc_multiplier(mrst)<0) &
        contracted(:,1:6,:,:) = -contracted(:,1:6,:,:)
      if (abs(response_scale) > epsilon(1.0_dp)) then
        contracted(:,6,:,:) = contracted(:,6,:,:) * &
          infos%tddft%spc_coco/response_scale
        contracted(:,5,:,:) = contracted(:,5,:,:) * &
          infos%tddft%spc_ovov/response_scale
        contracted(:,1:4,:,:) = contracted(:,1:4,:,:) * &
          infos%tddft%spc_coov/response_scale
      end if
    end if

    do ivec=1,nvec
      call mrsfmntoia(infos,contracted(ivec,:,:,:),sigma,mo_a,mo_b,ivec)
      call iatogen(trial(:,ivec),work,nocca,noccb)
      call mrsfesum(infos,work,fa,fb,sigma,ivec)
    end do
    nullify(contracted)
    if (grouped_channels) then
      deallocate(grouped_contracted)
    else
      call int2_data%clean()
    end if
    deallocate(density,work)
  end subroutine apply_mrsf_tda_sigma

end module tdhf_mrsf_sigma_mod
