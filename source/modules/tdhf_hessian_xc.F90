module tdhf_hessian_xc_mod
  use precision, only: dp
  implicit none
  private
  public :: build_tdhf_xc_fixed_hessian, add_tdhf_xc_response_rows
contains
  subroutine build_tdhf_xc_fixed_hessian(infos,hxc)
    ! Excitation-only XC skeleton: d/dR [gxc(D,P,X)-gxc(D,0,0)].
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, OQP_TD_P, OQP_TD_XPY
    use mathlib, only: unpack_matrix, symmetrize_matrix, orthogonal_transform
    use tdhf_lib, only: iatogen
    use mod_dft, only: dft_initialize, dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use messages, only: show_message, WITH_ABORT
    type(information),target,intent(inout)::infos
    real(dp),intent(out)::hxc(:,:)
    type(basis_set),pointer::basis
    type(dft_grid_t)::grid
    real(dp),contiguous,pointer::dpk(:),c(:,:),ppk(:,:),xpy(:,:)
    real(dp),allocatable,target::d(:,:),p(:,:,:),x(:,:,:)
    real(dp),allocatable::wrk(:,:),tmp(:,:)
    integer::nbf,nocc,natom,ncart
    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=basis%nbf; nocc=infos%mol_prop%nocc; natom=size(basis%atoms%xyz,2); ncart=3*natom
    if(any(shape(hxc)/=[ncart,ncart])) call show_message('TDDFT XC Hessian has wrong shape.',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_DM_A,dpk); call tagarray_get_data(infos%dat,OQP_VEC_MO_A,c)
    call tagarray_get_data(infos%dat,OQP_TD_P,ppk); call tagarray_get_data(infos%dat,OQP_TD_XPY,xpy)
    allocate(d(nbf,nbf),p(nbf,nbf,1),x(nbf,nbf,1),wrk(nbf,nbf),tmp(nbf,nbf))
    call unpack_matrix(dpk,d); call unpack_matrix(ppk(:,1),p(:,:,1))
    call iatogen(xpy(:,infos%tddft%target_state),wrk,nocc,nocc)
    call symmetrize_matrix(wrk,nbf); wrk=0.5_dp*wrk
    call orthogonal_transform('t',nbf,c,wrk,x(:,:,1),tmp); hxc=0.0_dp
    if (infos%functional%needGrd) then
      block
        use mod_dft_gridint_tdgga_hessian_driver, only: tddft_gga_fixed_hessian
        call dft_initialize(infos,basis,grid)
        call tddft_gga_fixed_hessian(basis,grid,d,p(:,:,1),x(:,:,1),hxc,infos)
        call dftclean(infos)
      end block
      deallocate(d,p,x,wrk,tmp)
      return
    end if
    block
      use mod_dft_gridint_tdlda_driver, only: tddft_lda_fixed_hessian
      call dft_initialize(infos,basis,grid)
      call tddft_lda_fixed_hessian(basis,grid,d,p(:,:,1),x(:,:,1),hxc,infos)
      call dftclean(infos)
    end block
    deallocate(d,p,x,wrk,tmp)
  end subroutine build_tdhf_xc_fixed_hessian

  subroutine add_tdhf_xc_response_rows(infos,dground,dprel,dxpy,rows)
    ! GAMESS DXR by density polarization of the analytic XC gradient.
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, OQP_TD_P, OQP_TD_XPY
    use mathlib, only: unpack_matrix, symmetrize_matrix, orthogonal_transform
    use tdhf_lib, only: iatogen
    use mod_dft, only: dft_initialize, dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use messages, only: show_message, WITH_ABORT
    type(information),target,intent(inout)::infos
    real(dp),intent(in)::dground(:,:,:),dprel(:,:,:),dxpy(:,:,:)
    real(dp),intent(inout)::rows(:,:)
    type(basis_set),pointer::basis
    type(dft_grid_t)::grid
    real(dp),contiguous,pointer::dpk(:),c(:,:),ppk(:,:),xpy(:,:)
    real(dp),allocatable,target::d0(:,:),p0(:,:,:),x0(:,:,:)
    real(dp),allocatable::wrk(:,:),tmp(:,:)
    integer::nbf,nocc,natom,ncart
    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=basis%nbf; nocc=infos%mol_prop%nocc; natom=size(basis%atoms%xyz,2); ncart=3*natom
    if(any(shape(rows)/=[ncart,ncart])) call show_message('TDDFT XC rows have wrong shape.',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_DM_A,dpk); call tagarray_get_data(infos%dat,OQP_VEC_MO_A,c)
    call tagarray_get_data(infos%dat,OQP_TD_P,ppk); call tagarray_get_data(infos%dat,OQP_TD_XPY,xpy)
    allocate(d0(nbf,nbf),p0(nbf,nbf,1),x0(nbf,nbf,1),wrk(nbf,nbf),tmp(nbf,nbf))
    call unpack_matrix(dpk,d0); call unpack_matrix(ppk(:,1),p0(:,:,1))
    call iatogen(xpy(:,infos%tddft%target_state),wrk,nocc,nocc); call symmetrize_matrix(wrk,nbf); wrk=0.5_dp*wrk
    call orthogonal_transform('t',nbf,c,wrk,x0(:,:,1),tmp)
    if (infos%functional%needGrd) then
      block
        use mod_dft_gridint_tdgga_response_driver, only: tddft_gga_response_rows
        real(dp), allocatable :: direct_rows(:,:)
        allocate(direct_rows(ncart,ncart))
        call dft_initialize(infos,basis,grid)
        call tddft_gga_response_rows(basis,grid,d0,p0(:,:,1),x0(:,:,1), &
          dground,dxpy,direct_rows,infos)
        call dftclean(infos)
        rows=rows+direct_rows
      end block
      deallocate(d0,p0,x0,wrk,tmp)
      return
    end if
    block
      use mod_dft_gridint_tdlda_driver, only: tddft_lda_response_rows
      real(dp), allocatable :: direct_rows(:,:)
      allocate(direct_rows(ncart,ncart))
      call dft_initialize(infos,basis,grid)
      call tddft_lda_response_rows(basis,grid,d0,p0(:,:,1),x0(:,:,1), &
        dground,dxpy,direct_rows,infos)
      call dftclean(infos)
      rows=rows+direct_rows
    end block
    deallocate(d0,p0,x0,wrk,tmp)
  end subroutine add_tdhf_xc_response_rows
end module tdhf_hessian_xc_mod
