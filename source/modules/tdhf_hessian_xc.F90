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
    real(dp),allocatable,target::d(:,:),p(:,:,:),x(:,:,:),z(:,:,:)
    real(dp),allocatable::wrk(:,:),tmp(:,:),gp(:,:),gm(:,:),g0p(:,:),g0m(:,:),hgrad(:,:),hhalf(:,:)
    real(dp),parameter::step=1.0e-3_dp
    integer::k,cc,aa,nbf,nocc,natom,ncart
    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=basis%nbf; nocc=infos%mol_prop%nocc; natom=size(basis%atoms%xyz,2); ncart=3*natom
    if(any(shape(hxc)/=[ncart,ncart])) call show_message('TDDFT XC Hessian has wrong shape.',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_DM_A,dpk); call tagarray_get_data(infos%dat,OQP_VEC_MO_A,c)
    call tagarray_get_data(infos%dat,OQP_TD_P,ppk); call tagarray_get_data(infos%dat,OQP_TD_XPY,xpy)
    allocate(d(nbf,nbf),p(nbf,nbf,1),x(nbf,nbf,1),z(nbf,nbf,1),wrk(nbf,nbf),tmp(nbf,nbf))
    allocate(gp(3,natom),gm(3,natom),g0p(3,natom),g0m(3,natom),source=0.0_dp)
    call unpack_matrix(dpk,d); call unpack_matrix(ppk(:,1),p(:,:,1))
    call iatogen(xpy(:,infos%tddft%target_state),wrk,nocc,nocc)
    call symmetrize_matrix(wrk,nbf); wrk=0.5_dp*wrk
    call orthogonal_transform('t',nbf,c,wrk,x(:,:,1),tmp); z=0.0_dp; hxc=0.0_dp
    if (infos%functional%needGrd) then
      block
        use mod_dft_gridint_tdgga_hessian_driver, only: tddft_gga_fixed_hessian
        call dft_initialize(infos,basis,grid)
        call tddft_gga_fixed_hessian(basis,grid,d,p(:,:,1),x(:,:,1),hxc,infos)
        call dftclean(infos)
      end block
      deallocate(d,p,x,z,wrk,tmp,gp,gm,g0p,g0m)
      return
    end if
    do k=1,ncart
      cc=mod(k-1,3)+1; aa=(k-1)/3+1
      basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)+step; call basis%init_shell_centers()
      call dft_initialize(infos,basis,grid)
      call xc_gradient(basis,grid,d,p,x,gp,infos); call xc_gradient(basis,grid,d,z,z,g0p,infos)
      call dftclean(infos)
      basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)-2.0_dp*step; call basis%init_shell_centers()
      call dft_initialize(infos,basis,grid)
      call xc_gradient(basis,grid,d,p,x,gm,infos); call xc_gradient(basis,grid,d,z,z,g0m,infos)
      call dftclean(infos)
      basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)+step; call basis%init_shell_centers()
      hxc(:,k)=reshape(((gp-g0p)-(gm-g0m))/(2.0_dp*step),[ncart])
    end do
    allocate(hgrad(ncart,ncart),source=hxc); allocate(hhalf(ncart,ncart))
    call build_scalar_fixed_hessian(infos,basis,d,p(:,:,1),x(:,:,1),1.0e-2_dp,hxc)
    call build_scalar_fixed_hessian(infos,basis,d,p(:,:,1),x(:,:,1),5.0e-3_dp,hhalf)
    block
      use io_constants, only: iw
      open(unit=iw,file=infos%log_filename,position='append')
      write(iw,'(A,1P,E12.4)') 'TDDFT fixed XC scalar-vs-gradient FD maximum: ',maxval(abs(hxc-hgrad))
      write(iw,'(A,1P,E12.4)') 'TDDFT fixed XC density-step stability maximum: ',maxval(abs(hhalf-hxc))
      close(iw)
    end block
    ! dmatd_density_blk receives the alpha density in the restricted path.
    ! The OpenQP transition-density contraction already carries the paired
    ! occupied-virtual normalization; apply only the remaining spin factor.
    hxc=2.0_dp*(4.0_dp*hhalf-hxc)/3.0_dp
    deallocate(d,p,x,z,wrk,tmp,gp,gm,g0p,g0m,hgrad,hhalf)
  end subroutine build_tdhf_xc_fixed_hessian

  subroutine build_scalar_fixed_hessian(infos,basis,d,p,x,density_step,h)
    use types, only: information
    use basis_tools, only: basis_set
    type(information),target,intent(inout)::infos
    type(basis_set),intent(inout)::basis
    real(dp),intent(in)::d(:,:),p(:,:),x(:,:)
    real(dp),intent(in)::density_step
    real(dp),intent(out)::h(:,:)
    real(dp),parameter::step=1.0e-3_dp
    real(dp)::q0,qpp,qpm,qmp,qmm,qp,qm
    integer::k,l,ck,ak,cl,al,ncart
    ncart=size(h,1); h=0.0_dp
    q0=xc_lagrangian_value(infos,basis,d,p,x,density_step)
    do k=1,ncart
      ck=mod(k-1,3)+1; ak=(k-1)/3+1
      basis%atoms%xyz(ck,ak)=basis%atoms%xyz(ck,ak)+step; call basis%init_shell_centers()
      qp=xc_lagrangian_value(infos,basis,d,p,x,density_step)
      basis%atoms%xyz(ck,ak)=basis%atoms%xyz(ck,ak)-2.0_dp*step; call basis%init_shell_centers()
      qm=xc_lagrangian_value(infos,basis,d,p,x,density_step)
      basis%atoms%xyz(ck,ak)=basis%atoms%xyz(ck,ak)+step; call basis%init_shell_centers()
      h(k,k)=(qp-2.0_dp*q0+qm)/(step*step)
      do l=k+1,ncart
        cl=mod(l-1,3)+1; al=(l-1)/3+1
        basis%atoms%xyz(ck,ak)=basis%atoms%xyz(ck,ak)+step
        basis%atoms%xyz(cl,al)=basis%atoms%xyz(cl,al)+step; call basis%init_shell_centers()
        qpp=xc_lagrangian_value(infos,basis,d,p,x,density_step)
        basis%atoms%xyz(cl,al)=basis%atoms%xyz(cl,al)-2.0_dp*step; call basis%init_shell_centers()
        qpm=xc_lagrangian_value(infos,basis,d,p,x,density_step)
        basis%atoms%xyz(ck,ak)=basis%atoms%xyz(ck,ak)-2.0_dp*step
        basis%atoms%xyz(cl,al)=basis%atoms%xyz(cl,al)+2.0_dp*step; call basis%init_shell_centers()
        qmp=xc_lagrangian_value(infos,basis,d,p,x,density_step)
        basis%atoms%xyz(cl,al)=basis%atoms%xyz(cl,al)-2.0_dp*step; call basis%init_shell_centers()
        qmm=xc_lagrangian_value(infos,basis,d,p,x,density_step)
        basis%atoms%xyz(ck,ak)=basis%atoms%xyz(ck,ak)+step
        basis%atoms%xyz(cl,al)=basis%atoms%xyz(cl,al)+step; call basis%init_shell_centers()
        h(k,l)=(qpp-qpm-qmp+qmm)/(4.0_dp*step*step); h(l,k)=h(k,l)
      end do
    end do
  end subroutine build_scalar_fixed_hessian

  function xc_lagrangian_value(infos,basis,d,p,x,dens_step) result(q)
    use types, only: information
    use basis_tools, only: basis_set
    use mod_dft, only: dft_initialize,dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_energy, only: dmatd_density_blk
    type(information),target,intent(inout)::infos
    type(basis_set),intent(inout)::basis
    real(dp),intent(in)::d(:,:),p(:,:),x(:,:)
    real(dp),intent(in)::dens_step
    real(dp)::q,ep,em,ex0,exp,exm,ne,ke
    type(dft_grid_t)::grid
    real(dp),allocatable::fa(:),fb(:,:),den(:,:)
    integer::nbf,nbf2,nang
    nbf=size(d,1); nbf2=nbf*(nbf+1)/2; nang=maxval(basis%am)+2
    allocate(fa(nbf2),fb(nbf2,1),den(nbf,nbf)); call dft_initialize(infos,basis,grid)
    den=d+dens_step*p; call xc_energy(den,ep)
    den=d-dens_step*p; call xc_energy(den,em)
    den=d; call xc_energy(den,ex0)
    den=d+dens_step*x; call xc_energy(den,exp)
    den=d-dens_step*x; call xc_energy(den,exm)
    q=(ep-em)/(2.0_dp*dens_step) &
      +2.0_dp*(exp-2.0_dp*ex0+exm)/(dens_step*dens_step)
    call dftclean(infos); deallocate(fa,fb,den)
  contains
    subroutine xc_energy(dena,e)
      real(dp),intent(in),target::dena(:,:)
      real(dp),intent(out)::e
      call dmatd_density_blk(basis,grid,dena,dena,fa,fb,e,ne,ke,nang,nbf, &
        infos%dft%grid_density_cutoff,.false.,infos)
    end subroutine xc_energy
  end function xc_lagrangian_value

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
    real(dp),allocatable,target::d0(:,:),p0(:,:,:),x0(:,:,:),z(:,:,:),dp1(:,:),dm1(:,:),pp(:,:,:),pm(:,:,:),xp(:,:,:),xm(:,:,:)
    real(dp),allocatable::wrk(:,:),tmp(:,:),gp(:,:),gm(:,:),g0p(:,:),g0m(:,:)
    real(dp),parameter::step=1.0e-3_dp
    integer::k,nbf,nocc,natom,ncart
    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=basis%nbf; nocc=infos%mol_prop%nocc; natom=size(basis%atoms%xyz,2); ncart=3*natom
    if(any(shape(rows)/=[ncart,ncart])) call show_message('TDDFT XC rows have wrong shape.',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_DM_A,dpk); call tagarray_get_data(infos%dat,OQP_VEC_MO_A,c)
    call tagarray_get_data(infos%dat,OQP_TD_P,ppk); call tagarray_get_data(infos%dat,OQP_TD_XPY,xpy)
    allocate(d0(nbf,nbf),p0(nbf,nbf,1),x0(nbf,nbf,1),z(nbf,nbf,1),dp1(nbf,nbf),dm1(nbf,nbf))
    allocate(pp(nbf,nbf,1),pm(nbf,nbf,1),xp(nbf,nbf,1),xm(nbf,nbf,1),wrk(nbf,nbf),tmp(nbf,nbf))
    allocate(gp(3,natom),gm(3,natom),g0p(3,natom),g0m(3,natom),source=0.0_dp)
    call unpack_matrix(dpk,d0); call unpack_matrix(ppk(:,1),p0(:,:,1))
    call iatogen(xpy(:,infos%tddft%target_state),wrk,nocc,nocc); call symmetrize_matrix(wrk,nbf); wrk=0.5_dp*wrk
    call orthogonal_transform('t',nbf,c,wrk,x0(:,:,1),tmp); z=0.0_dp
    call dft_initialize(infos,basis,grid)
    do k=1,ncart
      dp1=d0+step*dground(:,:,k); dm1=d0-step*dground(:,:,k)
      ! GAMESS TDHXR1 keeps Peff frozen.  Its response row contains dD_K
      ! and dPv_K, but no dPeff_K; the latter belongs to the ordinary D1G/D2G
      ! response rows and would be counted twice here.
      pp=p0; pm=p0
      xp=x0+step*dxpy(:,:,k:k); xm=x0-step*dxpy(:,:,k:k)
      call xc_gradient(basis,grid,dp1,pp,xp,gp,infos); call xc_gradient(basis,grid,dm1,pm,xm,gm,infos)
      call xc_gradient(basis,grid,dp1,z,z,g0p,infos); call xc_gradient(basis,grid,dm1,z,z,g0m,infos)
      rows(:,k)=rows(:,k)+reshape(((gp-g0p)-(gm-g0m))/(2.0_dp*step),[ncart])
    end do
    call dftclean(infos)
    deallocate(d0,p0,x0,z,dp1,dm1,pp,pm,xp,xm,wrk,tmp,gp,gm,g0p,g0m)
  end subroutine add_tdhf_xc_response_rows

  subroutine xc_gradient(basis,grid,d,p,x,g,infos)
    use types, only: information
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: tddft_xc_gradient
    type(basis_set),intent(inout)::basis
    type(dft_grid_t),target,intent(in)::grid
    real(dp),contiguous,target,intent(inout)::d(:,:)
    real(dp),target,intent(inout)::p(:,:,:),x(:,:,:)
    real(dp),intent(out)::g(:,:)
    type(information),target,intent(in)::infos
    call tddft_xc_gradient(basis,grid,g,d,p,x,1,1.0e-14_dp,infos,.true.)
  end subroutine xc_gradient
end module tdhf_hessian_xc_mod
