module tdhf_mrsf_hessian_intermediates_mod

  use precision, only: dp

  implicit none

  private

  type, public :: mrsf_hessian_z_intermediates_t
    real(kind=dp), allocatable :: mo_energy(:),fa(:,:),fb(:,:)
    real(kind=dp), allocatable :: hxa(:,:),hxb(:,:),ab1_a(:,:),ab1_b(:,:)
    real(kind=dp), allocatable :: z_rhs_hxa(:,:),z_rhs_hxb(:,:)
    real(kind=dp), allocatable :: z_ab1_a(:,:),z_ab1_b(:,:)
    real(kind=dp), allocatable :: tij(:,:),tab(:,:)
    real(kind=dp), allocatable :: dmo_energy(:,:),dfa(:,:,:),dfb(:,:,:)
    real(kind=dp), allocatable :: dhxa(:,:,:),dhxb(:,:,:)
    real(kind=dp), allocatable :: dz_rhs_hxa(:,:,:),dz_rhs_hxb(:,:,:)
    real(kind=dp), allocatable :: dab1_a(:,:,:),dab1_b(:,:,:)
    real(kind=dp), allocatable :: dz_ab1_a(:,:,:),dz_ab1_b(:,:,:)
    real(kind=dp), allocatable :: dtij(:,:,:),dtab(:,:,:)
    real(kind=dp), allocatable :: dseven(:,:,:,:),dabxc(:,:,:)
    real(kind=dp), allocatable :: dta_packed(:,:),dtb_packed(:,:)
    real(kind=dp), allocatable :: ppija(:,:),ppijb(:,:)
    real(kind=dp), allocatable :: dppija(:,:,:),dppijb(:,:,:)
  contains
    procedure :: clean => mrsf_hessian_z_intermediates_clean
  end type mrsf_hessian_z_intermediates_t

  public :: build_mrsf_hf_z_intermediates
  public :: build_mrsf_hf_w_intermediates

contains

!###############################################################################

  subroutine build_mrsf_hf_z_intermediates(infos,int2_driver,mo_a,mo_b, &
      dmo_a,dmo_b,fock_a_ao,fock_b_ao,dfock_a_ao,dfock_b_ao,x,dx,z, &
      data,status)
    ! Build the baseline and first nuclear derivatives entering the MRSF
    ! orbital-adjoint equation.  This is Hiroya Nakata's TDHF/TDDFT Hessian
    ! differentiation strategy applied to the physical seven-density,
    ! two-SOMO spin-adapted MRSF response.  No configuration-state expansion
    ! or displaced nuclear geometry is introduced.

    use types, only: information
    use int2_compute, only: int2_compute_t
    use oqp_tagarray_driver, only: tagarray_get_data,data_has_tags, &
      OQP_E_MO_A,OQP_td_t,OQP_td_mrsf_density,OQP_td_mrsf_hxa, &
      OQP_td_mrsf_hxb
    use mathlib, only: unpack_matrix
    use tdhf_lib, only: iatogen
    use tdhf_sf_lib, only: sfrogen
    use tdhf_mrsf_lib, only: mrsfxvec,mrsfsp
    use tdhf_mrsf_sigma_mod, only: prepare_mrsf_response_scaling
    use tdhf_mrsf_density_contraction_mod, only: &
      contract_mrsf_seven_density_batch
    use tdhf_mrsf_hessian_eri_derivative_mod, only: &
      build_mrsf_explicit_eri_derivative
    use tdhf_mrsf_hessian_unrelaxed_density_mod, only: &
      build_mrsf_unrelaxed_density_derivatives, &
      build_mrsf_mo_difference_density_derivatives
    use tdhf_mrsf_hessian_mo_response_mod, only: &
      build_mrsf_mo_fock_derivatives
    use messages, only: WITH_ABORT

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: fock_a_ao(:,:),fock_b_ao(:,:)
    real(kind=dp), intent(in) :: dfock_a_ao(:,:,:),dfock_b_ao(:,:,:)
    real(kind=dp), intent(in) :: x(:),dx(:,:),z(:)
    type(mrsf_hessian_z_intermediates_t), intent(inout) :: data
    integer, intent(out) :: status

    character(len=*), parameter :: required(*)=[character(len=80) :: &
      OQP_E_MO_A,OQP_td_t,OQP_td_mrsf_density,OQP_td_mrsf_hxa, &
      OQP_td_mrsf_hxb]
    real(kind=dp), contiguous, pointer :: stored_energy(:),td_t(:,:), &
      seven(:,:,:),stored_hxa(:,:),stored_hxb(:,:)
    real(kind=dp), allocatable, target :: base_batch(:,:,:,:),derivative_batch(:,:,:,:)
    real(kind=dp), allocatable :: base_fock(:,:,:,:),response_fock(:,:,:,:), &
      explicit_fock(:,:,:,:),total_fock(:,:,:,:),ta(:,:),tb(:,:), &
      dta(:,:,:),dtb(:,:,:),f_ta(:,:),f_tb(:,:),df_ta(:,:,:),df_tb(:,:,:), &
      f_za(:,:),f_zb(:,:),df_za(:,:,:),df_zb(:,:,:),ava(:,:),avb(:,:), &
      dava(:,:,:),davb(:,:,:),x_expanded(:),dx_expanded(:),x_matrix(:,:), &
      dx_matrix(:,:),g7(:,:),dg7(:,:,:),tmp_a(:,:),tmp_b(:,:), &
      beta_energy_derivative(:,:)
    real(kind=dp) :: response_scale,baseline_error_a,baseline_error_b
    integer :: coordinate,nbf,nbf2,nocca,noccb,nvira,nvirb,ncoord, &
      packed,lzdim,local_status

    call data%clean()
    nbf=infos%basis%nbf
    nbf2=nbf*(nbf+1)/2
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    nvira=nbf-nocca
    nvirb=nbf-noccb
    ncoord=size(dmo_a,3)
    packed=nocca*nvirb
    lzdim=noccb*(nocca-noccb+nvira)+(nocca-noccb)*nvira
    status=0
    if(infos%control%hamilton==20) then
      status=-2
      return
    end if
    if(nbf<=0 .or. ncoord<=0 .or. nocca-noccb/=2 .or. &
       size(x)/=packed .or. any(shape(dx)/=[packed,ncoord]) .or. &
       size(z)/=lzdim .or. any(shape(mo_a)/=[nbf,nbf]) .or. &
       any(shape(mo_b)/=[nbf,nbf]) .or. &
       any(shape(dmo_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dmo_b)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dfock_a_ao)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dfock_b_ao)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    call prepare_mrsf_response_scaling(infos,response_scale,local_status)
    if(local_status/=0) then
      status=-3
      return
    end if
    call data_has_tags(infos%dat,required, &
      'tdhf_mrsf_hessian_intermediates_mod', &
      'build_mrsf_hf_z_intermediates',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_E_MO_A,stored_energy)
    call tagarray_get_data(infos%dat,OQP_td_t,td_t)
    call tagarray_get_data(infos%dat,OQP_td_mrsf_density,seven)
    call tagarray_get_data(infos%dat,OQP_td_mrsf_hxa,stored_hxa)
    call tagarray_get_data(infos%dat,OQP_td_mrsf_hxb,stored_hxb)

    allocate(data%mo_energy(nbf),data%fa(nbf,nbf),data%fb(nbf,nbf), &
      data%hxa(nbf,nocca),data%hxb(nbf,nbf), &
      data%z_rhs_hxa(nbf,nocca),data%z_rhs_hxb(nbf,nbf), &
      data%ab1_a(nocca,nvira),data%ab1_b(noccb,nvirb), &
      data%z_ab1_a(nocca,nvira),data%z_ab1_b(noccb,nvirb), &
      data%tij(nocca,nocca),data%tab(nvirb,nvirb), &
      data%dmo_energy(nbf,ncoord),data%dfa(nbf,nbf,ncoord), &
      data%dfb(nbf,nbf,ncoord),data%dhxa(nbf,nocca,ncoord), &
      data%dhxb(nbf,nbf,ncoord),data%dab1_a(nocca,nvira,ncoord), &
      data%dz_rhs_hxa(nbf,nocca,ncoord), &
      data%dz_rhs_hxb(nbf,nbf,ncoord), &
      data%dab1_b(noccb,nvirb,ncoord), &
      data%dz_ab1_a(nocca,nvira,ncoord), &
      data%dz_ab1_b(noccb,nvirb,ncoord), &
      data%dtij(nocca,nocca,ncoord),data%dtab(nvirb,nvirb,ncoord), &
      data%dseven(7,nbf,nbf,ncoord),data%dabxc(nbf,nbf,ncoord), &
      data%dta_packed(nbf2,ncoord),data%dtb_packed(nbf2,ncoord), &
      source=0.0_dp)
    allocate(beta_energy_derivative(nbf,ncoord),source=0.0_dp)
    data%mo_energy=stored_energy
    data%fa=matmul(transpose(mo_a),matmul(fock_a_ao,mo_a))
    data%fb=matmul(transpose(mo_b),matmul(fock_b_ao,mo_b))
    call build_mrsf_mo_fock_derivatives(mo_a,dmo_a,fock_a_ao, &
      dfock_a_ao,data%dfa,data%dmo_energy,local_status)
    if(local_status==0) call build_mrsf_mo_fock_derivatives(mo_b,dmo_b, &
      fock_b_ao,dfock_b_ao,data%dfb,beta_energy_derivative,local_status)
    if(local_status==0) then
      ! The orbital energies entering sfrolhs are eigenvalues of the
      ! Guest--Saunders effective ROHF Fock operator.  At zero level shift its
      ! diagonal is the equal alpha/beta average in each C, O, and V sector;
      ! differentiating only the physical alpha Fock gives the wrong CO, CV,
      ! and especially OV energy-gap response.
      data%dmo_energy=0.5_dp*(data%dmo_energy+beta_energy_derivative)
    end if
    deallocate(beta_energy_derivative)
    if(local_status/=0) then
      status=-4
      call data%clean()
      return
    end if

    allocate(ta(nbf,nbf),tb(nbf,nbf),dta(nbf,nbf,ncoord), &
      dtb(nbf,nbf,ncoord),source=0.0_dp)
    call unpack_matrix(td_t(:,1),ta)
    call unpack_matrix(td_t(:,2),tb)
    call build_mrsf_unrelaxed_density_derivatives(infos,mo_a,mo_b,dmo_a, &
      dmo_b,x,dx,data%dseven,data%dabxc,data%dta_packed, &
      data%dtb_packed,local_status)
    if(local_status==0) call build_mrsf_mo_difference_density_derivatives( &
      infos,x,dx,data%tij,data%tab,data%dtij,data%dtab,local_status)
    if(local_status/=0) then
      status=-5
      call data%clean()
      deallocate(ta,tb,dta,dtb)
      return
    end if
    do coordinate=1,ncoord
      call unpack_matrix(data%dta_packed(:,coordinate),dta(:,:,coordinate))
      call unpack_matrix(data%dtb_packed(:,coordinate),dtb(:,:,coordinate))
    end do
    allocate(f_ta(nbf,nbf),f_tb(nbf,nbf),df_ta(nbf,nbf,ncoord), &
      df_tb(nbf,nbf,ncoord),source=0.0_dp)
    call build_hf_response_fock(infos,ta,tb,dta,dtb,f_ta,f_tb,df_ta, &
      df_tb,local_status)
    if(local_status/=0) then
      status=-6
      call data%clean()
      deallocate(ta,tb,dta,dtb,f_ta,f_tb,df_ta,df_tb)
      return
    end if
    call extract_mo_response_blocks(mo_a,dmo_a,f_ta,df_ta,nocca, &
      data%ab1_a,data%dab1_a)
    call extract_mo_response_blocks(mo_b,dmo_b,f_tb,df_tb,noccb, &
      data%ab1_b,data%dab1_b)

    allocate(ava(nbf,nbf),avb(nbf,nbf),dava(nbf,nbf,ncoord), &
      davb(nbf,nbf,ncoord),f_za(nbf,nbf),f_zb(nbf,nbf), &
      df_za(nbf,nbf,ncoord),df_zb(nbf,nbf,ncoord),source=0.0_dp)
    call sfrogen(ava,avb,z,nocca,noccb)
    call transform_mo_density_derivative(mo_a,dmo_a,ava,dava)
    call transform_mo_density_derivative(mo_b,dmo_b,avb,davb)
    ta=matmul(mo_a,matmul(ava,transpose(mo_a)))
    tb=matmul(mo_b,matmul(avb,transpose(mo_b)))
    call build_hf_response_fock(infos,ta,tb,dava,davb,f_za,f_zb,df_za, &
      df_zb,local_status)
    if(local_status/=0) then
      status=-7
      call data%clean()
      deallocate(ta,tb,dta,dtb,f_ta,f_tb,df_ta,df_tb,ava,avb,dava,davb, &
        f_za,f_zb,df_za,df_zb)
      return
    end if
    call extract_mo_response_blocks(mo_a,dmo_a,f_za,df_za,nocca, &
      data%z_ab1_a,data%dz_ab1_a)
    call extract_mo_response_blocks(mo_b,dmo_b,f_zb,df_zb,noccb, &
      data%z_ab1_b,data%dz_ab1_b)

    allocate(base_batch(1,7,nbf,nbf),base_fock(1,7,nbf,nbf), &
      derivative_batch(ncoord,7,nbf,nbf), &
      response_fock(ncoord,7,nbf,nbf), &
      explicit_fock(7,nbf,nbf,ncoord),total_fock(7,nbf,nbf,ncoord), &
      source=0.0_dp)
    base_batch(1,:,:,:)=seven
    do coordinate=1,ncoord
      derivative_batch(coordinate,:,:,:)=data%dseven(:,:,:,coordinate)
    end do
    call contract_mrsf_seven_density_batch(infos,int2_driver,base_batch, &
      response_scale,infos%tddft%spc_coco,infos%tddft%spc_ovov, &
      infos%tddft%spc_coov,base_fock,local_status)
    if(local_status==0) call contract_mrsf_seven_density_batch(infos, &
      int2_driver,derivative_batch,response_scale,infos%tddft%spc_coco, &
      infos%tddft%spc_ovov,infos%tddft%spc_coov,response_fock,local_status)
    if(local_status==0) call build_mrsf_explicit_eri_derivative(infos,seven, &
      response_scale,infos%tddft%spc_coco,infos%tddft%spc_ovov, &
      infos%tddft%spc_coov,infos%tddft%mult,explicit_fock,local_status)
    if(local_status/=0) then
      status=-8
      call data%clean()
      call cleanup_local()
      return
    end if
    do coordinate=1,ncoord
      total_fock(:,:,:,coordinate)=explicit_fock(:,:,:,coordinate)+ &
        response_fock(coordinate,:,:,:)
    end do

    allocate(x_expanded(packed),dx_expanded(packed),x_matrix(nbf,nbf), &
      dx_matrix(nbf,nbf),g7(nbf,nbf),dg7(nbf,nbf,ncoord), &
      tmp_a(nbf,nocca),tmp_b(nbf,nbf),source=0.0_dp)
    call mrsfxvec(infos,x,x_expanded)
    call iatogen(x_expanded,x_matrix,nocca,noccb)
    g7=matmul(transpose(mo_a),matmul(base_fock(1,7,:,:),mo_a))
    data%hxa=2.0_dp*matmul(g7,transpose(x_matrix(1:nocca,:)))
    data%hxb=2.0_dp*matmul(transpose(g7(1:nocca,:)),x_matrix(1:nocca,:))
    call mrsfsp(data%hxa,data%hxb,mo_a,mo_a,x_matrix,base_fock(1,:,:,:), &
      nocca,noccb)
    do coordinate=1,ncoord
      dg7(:,:,coordinate)= &
        matmul(transpose(dmo_a(:,:,coordinate)), &
          matmul(base_fock(1,7,:,:),mo_a))+ &
        matmul(transpose(mo_a), &
          matmul(total_fock(7,:,:,coordinate),mo_a))+ &
        matmul(transpose(mo_a), &
          matmul(base_fock(1,7,:,:),dmo_a(:,:,coordinate)))
      call mrsfxvec(infos,dx(:,coordinate),dx_expanded)
      call iatogen(dx_expanded,dx_matrix,nocca,noccb)
      data%dhxa(:,:,coordinate)=2.0_dp*( &
        matmul(dg7(:,:,coordinate),transpose(x_matrix(1:nocca,:)))+ &
        matmul(g7,transpose(dx_matrix(1:nocca,:))))
      data%dhxb(:,:,coordinate)=2.0_dp*( &
        matmul(transpose(dg7(1:nocca,:,coordinate)),x_matrix(1:nocca,:))+ &
        matmul(transpose(g7(1:nocca,:)),dx_matrix(1:nocca,:)))

      tmp_a=0.0_dp
      tmp_b=0.0_dp
      call mrsfsp(tmp_a,tmp_b,dmo_a(:,:,coordinate),mo_a,x_matrix, &
        base_fock(1,:,:,:),nocca,noccb)
      data%dhxa(:,:,coordinate)=data%dhxa(:,:,coordinate)+tmp_a
      data%dhxb(:,:,coordinate)=data%dhxb(:,:,coordinate)+tmp_b
      tmp_a=0.0_dp
      tmp_b=0.0_dp
      call mrsfsp(tmp_a,tmp_b,mo_a,dmo_a(:,:,coordinate),x_matrix, &
        base_fock(1,:,:,:),nocca,noccb)
      data%dhxa(:,:,coordinate)=data%dhxa(:,:,coordinate)+tmp_a
      data%dhxb(:,:,coordinate)=data%dhxb(:,:,coordinate)+tmp_b
      tmp_a=0.0_dp
      tmp_b=0.0_dp
      call mrsfsp(tmp_a,tmp_b,mo_a,mo_a,dx_matrix,base_fock(1,:,:,:), &
        nocca,noccb)
      data%dhxa(:,:,coordinate)=data%dhxa(:,:,coordinate)+tmp_a
      data%dhxb(:,:,coordinate)=data%dhxb(:,:,coordinate)+tmp_b
      tmp_a=0.0_dp
      tmp_b=0.0_dp
      call mrsfsp(tmp_a,tmp_b,mo_a,mo_a,x_matrix, &
        total_fock(:,:,:,coordinate),nocca,noccb)
      data%dhxa(:,:,coordinate)=data%dhxa(:,:,coordinate)+tmp_a
      data%dhxb(:,:,coordinate)=data%dhxb(:,:,coordinate)+tmp_b
    end do
    ! sfrorhs mutates its Hx arguments by adding 2*F*T.  The stored gradient
    ! Hx values therefore belong to W, whereas differentiating the Z-vector
    ! right-hand side must start from the pre-sfrorhs values and let sfrorhs
    ! add 2*F*T exactly once.  Preserve both definitions explicitly.
    data%z_rhs_hxa=data%hxa
    data%z_rhs_hxb=data%hxb
    data%dz_rhs_hxa=data%dhxa
    data%dz_rhs_hxb=data%dhxb
    baseline_error_a=maxval(abs(data%z_rhs_hxa-(stored_hxa-2.0_dp* &
      matmul(data%fa(:,1:nocca),data%tij))))
    baseline_error_b=maxval(abs(data%z_rhs_hxb(:,noccb+1:nbf)-( &
      stored_hxb(:,noccb+1:nbf)-2.0_dp* &
      matmul(data%fb(:,noccb+1:nbf),data%tab))))
    if(max(baseline_error_a,baseline_error_b)>1.0e-8_dp) then
      status=-9
      call data%clean()
      call cleanup_local()
      return
    end if
    data%hxa=stored_hxa
    data%hxb=stored_hxb
    do coordinate=1,ncoord
      data%dhxa(:,:,coordinate)=data%dz_rhs_hxa(:,:,coordinate)+ &
        2.0_dp*matmul(data%dfa(:,1:nocca,coordinate),data%tij)+ &
        2.0_dp*matmul(data%fa(:,1:nocca),data%dtij(:,:,coordinate))
      data%dhxb(:,:,coordinate)=data%dz_rhs_hxb(:,:,coordinate)
      data%dhxb(:,noccb+1:nbf,coordinate)= &
        data%dhxb(:,noccb+1:nbf,coordinate)+ &
        2.0_dp*matmul(data%dfb(:,noccb+1:nbf,coordinate),data%tab)+ &
        2.0_dp*matmul(data%fb(:,noccb+1:nbf),data%dtab(:,:,coordinate))
    end do
    call cleanup_local()

  contains

    subroutine cleanup_local()
      if(allocated(ta)) deallocate(ta,tb,dta,dtb)
      if(allocated(f_ta)) deallocate(f_ta,f_tb,df_ta,df_tb)
      if(allocated(ava)) deallocate(ava,avb,dava,davb,f_za,f_zb,df_za,df_zb)
      if(allocated(base_batch)) deallocate(base_batch,derivative_batch, &
        base_fock,response_fock,explicit_fock,total_fock)
      if(allocated(x_expanded)) deallocate(x_expanded,dx_expanded,x_matrix, &
        dx_matrix,g7,dg7,tmp_a,tmp_b)
    end subroutine cleanup_local

  end subroutine build_mrsf_hf_z_intermediates

!###############################################################################

  subroutine build_mrsf_hf_w_intermediates(infos,mo_a,mo_b,dmo_a,dmo_b, &
      tij,tab,dtij,dtab,z,dz,data,drelaxed_a,drelaxed_b,status)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data,data_has_tags,OQP_td_p, &
      OQP_td_mrsf_ppija,OQP_td_mrsf_ppijb
    use mathlib, only: unpack_matrix
    use tdhf_mrsf_hessian_relaxed_density_mod, only: &
      build_mrsf_relaxed_density_derivatives
    use messages, only: WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: tij(:,:),tab(:,:),dtij(:,:,:),dtab(:,:,:)
    real(kind=dp), intent(in) :: z(:),dz(:,:)
    type(mrsf_hessian_z_intermediates_t), intent(inout) :: data
    real(kind=dp), intent(out) :: drelaxed_a(:,:,:),drelaxed_b(:,:,:)
    integer, intent(out) :: status

    real(kind=dp), contiguous, pointer :: td_p(:,:),stored_ppija(:,:), &
      stored_ppijb(:,:)
    real(kind=dp), allocatable :: pa(:,:),pb(:,:),fpa(:,:),fpb(:,:), &
      dfpa(:,:,:),dfpb(:,:,:),dpa_fock(:,:,:),dpb_fock(:,:,:), &
      full_a(:,:),full_b(:,:),dfull_a(:,:,:),dfull_b(:,:,:)
    real(kind=dp) :: baseline_error_a,baseline_error_b
    integer :: nbf,nocca,noccb,ncoord,local_status

    nbf=infos%basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    ncoord=size(dmo_a,3)
    status=0
    if(infos%control%hamilton==20) then
      status=-2
      return
    end if
    if(any(shape(drelaxed_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(drelaxed_b)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    call build_mrsf_relaxed_density_derivatives(mo_a,mo_b,dmo_a,dmo_b, &
      tij,tab,dtij,dtab,z,dz,nocca,noccb,drelaxed_a,drelaxed_b, &
      local_status)
    if(local_status/=0) then
      status=-3
      return
    end if
    call data_has_tags(infos%dat,[character(len=80) :: OQP_td_p, &
      OQP_td_mrsf_ppija,OQP_td_mrsf_ppijb], &
      'tdhf_mrsf_hessian_intermediates_mod', &
      'build_mrsf_hf_w_intermediates',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_td_p,td_p)
    call tagarray_get_data(infos%dat,OQP_td_mrsf_ppija,stored_ppija)
    call tagarray_get_data(infos%dat,OQP_td_mrsf_ppijb,stored_ppijb)
    allocate(pa(nbf,nbf),pb(nbf,nbf),fpa(nbf,nbf),fpb(nbf,nbf), &
      dfpa(nbf,nbf,ncoord),dfpb(nbf,nbf,ncoord), &
      dpa_fock(nbf,nbf,ncoord),dpb_fock(nbf,nbf,ncoord), &
      full_a(nbf,nbf),full_b(nbf,nbf),dfull_a(nbf,nbf,ncoord), &
      dfull_b(nbf,nbf,ncoord),source=0.0_dp)
    call unpack_matrix(td_p(:,1),pa)
    call unpack_matrix(td_p(:,2),pb)
    dpa_fock=drelaxed_a
    dpb_fock=drelaxed_b
    call build_hf_response_fock(infos,pa,pb,dpa_fock,dpb_fock,fpa,fpb, &
      dfpa,dfpb,local_status)
    if(local_status/=0) then
      status=-4
      deallocate(pa,pb,fpa,fpb,dfpa,dfpb,dpa_fock,dpb_fock,full_a,full_b, &
        dfull_a,dfull_b)
      return
    end if
    call full_mo_response(mo_a,dmo_a,fpa,dfpa,full_a,dfull_a)
    call full_mo_response(mo_b,dmo_b,fpb,dfpb,full_b,dfull_b)
    allocate(data%ppija(nocca,nocca),data%ppijb(noccb,noccb), &
      data%dppija(nocca,nocca,ncoord), &
      data%dppijb(noccb,noccb,ncoord),source=0.0_dp)
    data%ppija=full_a(1:nocca,1:nocca)
    data%ppijb=full_b(1:noccb,1:noccb)
    baseline_error_a=maxval(abs(data%ppija-stored_ppija))
    baseline_error_b=maxval(abs(data%ppijb-stored_ppijb))
    if(max(baseline_error_a,baseline_error_b)>1.0e-8_dp) then
      status=-5
      deallocate(pa,pb,fpa,fpb,dfpa,dfpb,dpa_fock,dpb_fock,full_a,full_b, &
        dfull_a,dfull_b)
      return
    end if
    data%ppija=stored_ppija
    data%ppijb=stored_ppijb
    data%dppija=dfull_a(1:nocca,1:nocca,:)
    data%dppijb=dfull_b(1:noccb,1:noccb,:)
    deallocate(pa,pb,fpa,fpb,dfpa,dfpb,dpa_fock,dpb_fock,full_a,full_b, &
      dfull_a,dfull_b)
  end subroutine build_mrsf_hf_w_intermediates

!###############################################################################

  subroutine build_hf_response_fock(infos,pa,pb,dpa,dpb,fa,fb,dfa,dfb, &
      status)
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use fock_deriv_mod, only: fock_deriv_matrix_os

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: pa(:,:),pb(:,:),dpa(:,:,:),dpb(:,:,:)
    real(kind=dp), intent(out) :: fa(:,:),fb(:,:),dfa(:,:,:),dfb(:,:,:)
    integer, intent(out) :: status

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(int2_tdgrd_data_t) :: int2_data
    real(kind=dp), allocatable, target :: density(:,:,:)
    real(kind=dp), allocatable :: explicit_a(:,:,:,:),explicit_b(:,:,:,:), &
      ptot(:,:),sym_a(:,:),sym_b(:,:)
    integer :: atom,cart,coordinate,natom,nbf,ncoord

    basis=>infos%basis
    basis%atoms=>infos%atoms
    nbf=basis%nbf
    natom=size(infos%atoms%xyz,2)
    ncoord=3*natom
    status=0
    fa=0.0_dp
    fb=0.0_dp
    dfa=0.0_dp
    dfb=0.0_dp
    if(any(shape(pa)/=[nbf,nbf]) .or. any(shape(pb)/=[nbf,nbf]) .or. &
       any(shape(dpa)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dpb)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    allocate(density(nbf,nbf,2), &
      explicit_a(nbf,nbf,3,natom),explicit_b(nbf,nbf,3,natom), &
      ptot(nbf,nbf),sym_a(nbf,nbf),sym_b(nbf,nbf),source=0.0_dp)
    call int2_driver%init(basis,infos)
    call int2_driver%set_screening()
    density(:,:,1)=pa
    density(:,:,2)=pb
    int2_data=int2_tdgrd_data_t(d2=density,int_apb=.true., &
      int_amb=.false.,tamm_dancoff=.false.,scale_exchange=1.0_dp)
    call int2_driver%run(int2_data)
    fa=int2_data%apb(:,:,1,1)
    fb=int2_data%apb(:,:,2,1)
    call int2_data%clean()
    sym_a=0.5_dp*(pa+transpose(pa))
    sym_b=0.5_dp*(pb+transpose(pb))
    ptot=sym_a+sym_b
    call fock_deriv_matrix_os(infos,basis,ptot,sym_a,1.0_dp,explicit_a)
    call fock_deriv_matrix_os(infos,basis,ptot,sym_b,1.0_dp,explicit_b)
    coordinate=0
    do atom=1,natom
      do cart=1,3
        coordinate=coordinate+1
        density(:,:,1)=dpa(:,:,coordinate)
        density(:,:,2)=dpb(:,:,coordinate)
        int2_data=int2_tdgrd_data_t(d2=density,int_apb=.true., &
          int_amb=.false.,tamm_dancoff=.false.,scale_exchange=1.0_dp)
        call int2_driver%run(int2_data)
        ! int2_tdgrd_data_t returns the spin A+B response to D+D^T.  Its
        ! explicit nuclear derivative must follow the same convention; the
        ! ordinary open-shell Fock derivative is therefore doubled here.
        dfa(:,:,coordinate)=2.0_dp*explicit_a(:,:,cart,atom)+ &
          int2_data%apb(:,:,1,1)
        dfb(:,:,coordinate)=2.0_dp*explicit_b(:,:,cart,atom)+ &
          int2_data%apb(:,:,2,1)
        call int2_data%clean()
      end do
    end do
    call int2_driver%clean()
    deallocate(density,explicit_a,explicit_b,ptot,sym_a,sym_b)
  end subroutine build_hf_response_fock

!###############################################################################

  subroutine transform_mo_density_derivative(mo,dmo,density_mo, &
      density_ao_derivative)
    real(kind=dp), intent(in) :: mo(:,:),dmo(:,:,:),density_mo(:,:)
    real(kind=dp), intent(out) :: density_ao_derivative(:,:,:)
    integer :: coordinate
    do coordinate=1,size(dmo,3)
      density_ao_derivative(:,:,coordinate)= &
        matmul(dmo(:,:,coordinate),matmul(density_mo,transpose(mo)))+ &
        matmul(mo,matmul(density_mo,transpose(dmo(:,:,coordinate))))
    end do
  end subroutine transform_mo_density_derivative

!###############################################################################

  subroutine full_mo_response(mo,dmo,fock,dfock,fock_mo,dfock_mo)
    real(kind=dp), intent(in) :: mo(:,:),dmo(:,:,:),fock(:,:),dfock(:,:,:)
    real(kind=dp), intent(out) :: fock_mo(:,:),dfock_mo(:,:,:)
    integer :: coordinate
    fock_mo=matmul(transpose(mo),matmul(fock,mo))
    do coordinate=1,size(dmo,3)
      dfock_mo(:,:,coordinate)= &
        matmul(transpose(dmo(:,:,coordinate)),matmul(fock,mo))+ &
        matmul(transpose(mo),matmul(dfock(:,:,coordinate),mo))+ &
        matmul(transpose(mo),matmul(fock,dmo(:,:,coordinate)))
    end do
  end subroutine full_mo_response

!###############################################################################

  subroutine extract_mo_response_blocks(mo,dmo,fock,dfock,nocc,block,dblock)
    real(kind=dp), intent(in) :: mo(:,:),dmo(:,:,:),fock(:,:),dfock(:,:,:)
    integer, intent(in) :: nocc
    real(kind=dp), intent(out) :: block(:,:),dblock(:,:,:)
    real(kind=dp), allocatable :: full(:,:),dfull(:,:,:)
    integer :: nbf,ncoord
    nbf=size(mo,1)
    ncoord=size(dmo,3)
    allocate(full(nbf,nbf),dfull(nbf,nbf,ncoord))
    call full_mo_response(mo,dmo,fock,dfock,full,dfull)
    block=full(1:nocc,nocc+1:nbf)
    dblock=dfull(1:nocc,nocc+1:nbf,:)
    deallocate(full,dfull)
  end subroutine extract_mo_response_blocks

!###############################################################################

  subroutine mrsf_hessian_z_intermediates_clean(this)
    class(mrsf_hessian_z_intermediates_t), intent(inout) :: this
    if(allocated(this%mo_energy)) deallocate(this%mo_energy)
    if(allocated(this%fa)) deallocate(this%fa)
    if(allocated(this%fb)) deallocate(this%fb)
    if(allocated(this%hxa)) deallocate(this%hxa)
    if(allocated(this%hxb)) deallocate(this%hxb)
    if(allocated(this%z_rhs_hxa)) deallocate(this%z_rhs_hxa)
    if(allocated(this%z_rhs_hxb)) deallocate(this%z_rhs_hxb)
    if(allocated(this%ab1_a)) deallocate(this%ab1_a)
    if(allocated(this%ab1_b)) deallocate(this%ab1_b)
    if(allocated(this%z_ab1_a)) deallocate(this%z_ab1_a)
    if(allocated(this%z_ab1_b)) deallocate(this%z_ab1_b)
    if(allocated(this%tij)) deallocate(this%tij)
    if(allocated(this%tab)) deallocate(this%tab)
    if(allocated(this%dmo_energy)) deallocate(this%dmo_energy)
    if(allocated(this%dfa)) deallocate(this%dfa)
    if(allocated(this%dfb)) deallocate(this%dfb)
    if(allocated(this%dhxa)) deallocate(this%dhxa)
    if(allocated(this%dhxb)) deallocate(this%dhxb)
    if(allocated(this%dz_rhs_hxa)) deallocate(this%dz_rhs_hxa)
    if(allocated(this%dz_rhs_hxb)) deallocate(this%dz_rhs_hxb)
    if(allocated(this%dab1_a)) deallocate(this%dab1_a)
    if(allocated(this%dab1_b)) deallocate(this%dab1_b)
    if(allocated(this%dz_ab1_a)) deallocate(this%dz_ab1_a)
    if(allocated(this%dz_ab1_b)) deallocate(this%dz_ab1_b)
    if(allocated(this%dtij)) deallocate(this%dtij)
    if(allocated(this%dtab)) deallocate(this%dtab)
    if(allocated(this%dseven)) deallocate(this%dseven)
    if(allocated(this%dabxc)) deallocate(this%dabxc)
    if(allocated(this%dta_packed)) deallocate(this%dta_packed)
    if(allocated(this%dtb_packed)) deallocate(this%dtb_packed)
    if(allocated(this%ppija)) deallocate(this%ppija)
    if(allocated(this%ppijb)) deallocate(this%ppijb)
    if(allocated(this%dppija)) deallocate(this%dppija)
    if(allocated(this%dppijb)) deallocate(this%dppijb)
  end subroutine mrsf_hessian_z_intermediates_clean

end module tdhf_mrsf_hessian_intermediates_mod
