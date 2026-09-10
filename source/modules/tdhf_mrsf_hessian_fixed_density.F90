module tdhf_mrsf_hessian_fixed_density_mod

  use precision, only: dp

  implicit none

  private
  public :: build_mrsf_fixed_density_hessian

contains

!###############################################################################

  subroutine build_mrsf_fixed_density_hessian(infos,h_fixed,status)
    ! Fixed relaxed-density part of the total MRSF Cartesian Hessian.  This is
    ! the direct second derivative of the nuclear, one-electron, Pulay, ECP,
    ! and two-electron integrals while the converged state-specific densities
    ! are held fixed.  Semilocal XC quadrature derivatives and all first-order
    ! density/Z/amplitude response rows are separate terms.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data,data_has_tags,OQP_DM_A, &
      OQP_DM_B,OQP_WAO,OQP_td_p,OQP_td_mrsf_density
    use mathlib, only: unpack_matrix
    use grd1, only: eijden,hess_nn,hess_ee_overlap,hess_ee_kinetic,hess_en
    use grd2, only: grd2_hess_driver
    use ecp_tool, only: add_ecphess
    use tdhf_mrsf_gradient_mod, only: grd2_mrsf_compute_data_t
    use tdhf_mrsf_sigma_mod, only: prepare_mrsf_response_scaling
    use messages, only: WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: h_fixed(:,:)
    integer, intent(out) :: status

    character(len=*), parameter :: required(*)=[character(len=80) :: &
      OQP_DM_A,OQP_DM_B,OQP_WAO,OQP_td_p,OQP_td_mrsf_density]
    type(basis_set), pointer :: basis
    type(grd2_mrsf_compute_data_t) :: gcomp
    real(kind=dp), contiguous, pointer :: dmat_a(:),dmat_b(:),wao(:), &
      td_p(:,:),td_mrsf_density(:,:,:)
    real(kind=dp), allocatable, target :: d(:,:,:),p(:,:,:),spc(:,:,:)
    real(kind=dp), allocatable :: overlap_density(:),one_density(:)
    real(kind=dp) :: response_scale,reference_scale
    integer :: nbf,nbf2,natom,ncoord,local_status

    basis=>infos%basis
    basis%atoms=>infos%atoms
    nbf=basis%nbf
    nbf2=nbf*(nbf+1)/2
    natom=size(basis%atoms%xyz,2)
    ncoord=3*natom
    status=0
    h_fixed=0.0_dp
    if(nbf<=0 .or. natom<=0 .or. any(shape(h_fixed)/=[ncoord,ncoord]) .or. &
       infos%mol_prop%nelec_a-infos%mol_prop%nelec_b/=2 .or. &
       infos%tddft%umrsf .or. &
       (infos%tddft%mult/=1 .and. infos%tddft%mult/=3)) then
      status=-1
      return
    end if
    call prepare_mrsf_response_scaling(infos,response_scale,local_status)
    if(local_status/=0) then
      status=-2
      return
    end if
    reference_scale=1.0_dp
    if(infos%control%hamilton==20) reference_scale=infos%dft%hfscale
    call data_has_tags(infos%dat,required, &
      'tdhf_mrsf_hessian_fixed_density_mod', &
      'build_mrsf_fixed_density_hessian',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_DM_A,dmat_a)
    call tagarray_get_data(infos%dat,OQP_DM_B,dmat_b)
    call tagarray_get_data(infos%dat,OQP_WAO,wao)
    call tagarray_get_data(infos%dat,OQP_td_p,td_p)
    call tagarray_get_data(infos%dat,OQP_td_mrsf_density,td_mrsf_density)
    allocate(d(nbf,nbf,2),p(nbf,nbf,2),spc(7,nbf,nbf), &
      overlap_density(nbf2),one_density(nbf2),source=0.0_dp)
    call unpack_matrix(dmat_a,d(:,:,1))
    call unpack_matrix(dmat_b,d(:,:,2))
    call unpack_matrix(td_p(:,1),p(:,:,1))
    call unpack_matrix(td_p(:,2),p(:,:,2))
    spc=td_mrsf_density

    call hess_nn(basis%atoms,basis%ecp_zn_num,h_fixed)
    call eijden(overlap_density,nbf,infos)
    overlap_density=overlap_density+2.0_dp*wao
    one_density=dmat_a+dmat_b+td_p(:,1)+td_p(:,2)
    call hess_ee_overlap(basis,overlap_density,h_fixed)
    call hess_ee_kinetic(basis,one_density,h_fixed)
    call hess_en(basis,basis%atoms%xyz, &
      basis%atoms%zn-basis%ecp_zn_num,one_density,h_fixed)
    call add_ecphess(basis,basis%atoms%xyz,one_density,h_fixed)

    gcomp=grd2_mrsf_compute_data_t(d2=d,p2=p,spc2=spc,nbf=nbf, &
      hfscale=reference_scale,hfscale2=response_scale, &
      spcscale=[infos%tddft%spc_coco,infos%tddft%spc_ovov, &
                infos%tddft%spc_coov],mrst=infos%tddft%mult)
    call gcomp%init()
    call gcomp%build_cart(basis)
    call grd2_hess_driver(infos,basis,h_fixed,gcomp, &
      cam=infos%control%hamilton==20.and.infos%dft%cam_flag, &
      alpha=infos%tddft%cam_alpha,beta=infos%tddft%cam_beta, &
      mu=infos%tddft%cam_mu)
    call gcomp%clean()
    h_fixed=0.5_dp*(h_fixed+transpose(h_fixed))
    deallocate(d,p,spc,overlap_density,one_density)
  end subroutine build_mrsf_fixed_density_hessian

end module tdhf_mrsf_hessian_fixed_density_mod
