subroutine response_cache_storage_selftest(failed) bind(C)
  use iso_c_binding, only:c_int
  use precision, only:fp,i8b
  use mod_dft_gridint_response_cache
  implicit none
  integer(c_int),intent(out)::failed
  type(response_cache_t)::cache
  type(response_slice_t)::sample
  real(fp)::ao(16),weights(2),kernel(32,2),stats(5)
  integer::ids(2),i
  integer(i8b)::meta
  failed=0; ao=3; weights=4; kernel=5; stats=6; ids=[2,5]
  meta=3_i8b*int(storage_size(sample)/8,i8b)
  cache%max_bytes=meta+512+200
  call cache%prepare(1_i8b,3,2)
  call cache%store_ao(1,.false.,2,2,.false.,4,ids,ao,weights)
  call cache%store_ao(2,.false.,2,2,.false.,4,ids,ao,weights)
  call cache%store_ao(3,.true.,0,0,.true.,4,ids,ao,weights)
  call cache%store_kernel(1,.false.,kernel,stats)
  call cache%store_kernel(2,.false.,kernel,stats)
  if(.not.cache%slices(1)%ao_ready.or.cache%slices(2)%ao_ready) failed=failed+1
  if(.not.cache%slices(3)%ao_ready.or..not.cache%slices(3)%ao_skip) failed=failed+1
  if(.not.cache%slices(1)%kernel_ready.or.cache%slices(2)%kernel_ready) failed=failed+1
  if(cache%bytes>cache%max_bytes) failed=failed+1
  if(any(cache%slices(1)%ao/=ao).or.any(cache%slices(1)%kernel/=kernel)) failed=failed+1
  call cache%prepare(1_i8b,3,2)
  if(.not.cache%slices(1)%ao_ready) failed=failed+1
  call cache%prepare(2_i8b,3,2)
  if(cache%slices(1)%ao_ready.or.cache%slices(1)%kernel_ready) failed=failed+1
  cache%max_bytes=1
  call cache%prepare(2_i8b,3,2)
  if(cache%enabled.or.allocated(cache%slices).or.cache%bytes/=0) failed=failed+1
  call cache%free();call cache%free()
end subroutine

subroutine response_cache_compare(c_handle,error,failed) bind(C)
  use iso_c_binding,only:c_double,c_int
  use precision,only:fp,i8b
  use c_interop,only:oqp_handle_t,oqp_handle_get_info
  use types,only:information
  use basis_tools,only:basis_set
  use atomic_structure_m,only:atomic_structure
  use mod_dft,only:dft_initialize,dftclean
  use mod_dft_molgrid,only:dft_grid_t
  use mod_dft_gridint_fxc,only:utddft_fxc
  use mod_dft_gridint_response_cache,only:response_cache_t
  implicit none
  type(oqp_handle_t)::c_handle
  real(c_double),intent(out)::error
  integer(c_int),intent(out)::failed
  type(information),pointer::infos
  type(basis_set)::basis
  type(atomic_structure),target::atoms
  type(dft_grid_t)::grid
  type(response_cache_t),target::cache
  real(fp),allocatable::wa(:,:),wb(:,:),fa(:,:,:),fb(:,:,:),ra(:,:,:),rb(:,:,:)
  real(fp),allocatable,target::da(:,:,:),db(:,:,:)
  real(fp)::threshold
  integer::n,i,j,step,mode,point(2)
  integer(i8b)::limits(3)=[0_i8b,1048576_i8b,67108864_i8b]
  infos=>oqp_handle_get_info(c_handle)
  atoms=infos%atoms
  basis=infos%basis; basis%atoms=>atoms
  n=basis%nbf
  allocate(wa(n,n),wb(n,n),fa(n,n,1),fb(n,n,1),ra(n,n,1),rb(n,n,1),da(n,n,1),db(n,n,1))
  error=0;failed=0
  do mode=1,3
    atoms=infos%atoms
    basis=infos%basis;basis%atoms=>atoms
    call dft_initialize(infos,basis,grid)
    wa=0;wb=0
    do i=1,n
      wa(i,i)=0.4_fp;wb(i,i)=0.3_fp
    end do
    cache%max_bytes=limits(mode)
    threshold=1.0e-15_fp
    do step=1,8
      if(step==3) wa(1,1)=0.5_fp
      if(step==4) basis%cc(1)=basis%cc(1)*1.01_fp
      if(step==5) then
        point=maxloc(abs(grid%totWts))
        grid%totWts(point(1),point(2))=grid%totWts(point(1),point(2))*1.001_fp
      end if
      if(step==6) threshold=1.0e-12_fp
      if(step==7) then
        call dftclean(infos)
        call dft_initialize(infos,basis,grid)
      end if
      if(step==8) atoms%xyz(1,1)=atoms%xyz(1,1)+0.01_fp
      do j=1,n
        do i=1,n
          da(i,j,1)=0.01_fp*cos(real(i+j+step,fp))
          db(i,j,1)=0.02_fp*sin(real(i+j+2*step,fp))
        end do
      end do
      fa=0;fb=0;ra=0;rb=0
      call utddft_fxc(basis,grid,.false.,wa,wb,fa,fb,da,db,1,threshold,infos,cache=cache)
      call utddft_fxc(basis,grid,.false.,wa,wb,ra,rb,da,db,1,threshold,infos)
      error=max(error,maxval(abs(fa-ra)),maxval(abs(fb-rb)))
      if(cache%bytes>limits(mode)) failed=failed+1
      if(mode==2.and.step==2) then
        if(cache%kernel_hits==0.or.cache%kernel_hits>=grid%nSlices) failed=failed+1
      end if
      if(mode==3.and.step==2) then
        if(cache%ao_hits==0.or.cache%kernel_hits==0) failed=failed+1
      end if
      if(step>=3.and.cache%kernel_hits/=0) failed=failed+1
    end do
    call cache%free()
    call dftclean(infos)
  end do
  if(error>1.0e-10_fp) failed=failed+1
end subroutine

subroutine gradient_cache_compare(c_handle,error,failed) bind(C)
  use iso_c_binding,only:c_double,c_int
  use precision,only:fp,i8b
  use c_interop,only:oqp_handle_t,oqp_handle_get_info
  use types,only:information
  use basis_tools,only:basis_set
  use mod_dft,only:dft_initialize,dftclean
  use mod_dft_molgrid,only:dft_grid_t
  use mod_dft_gridint_tdxc_grad,only:utddft_xc_gradient,tddft_xc_gradient
  use mod_dft_gridint_grad,only:derexc_blk
  use mod_dft_gridint_response_cache,only:response_cache_t
  implicit none
  type(oqp_handle_t)::c_handle
  real(c_double),intent(out)::error
  integer(c_int),intent(out)::failed
  type(information),pointer::infos
  type(basis_set)::basis
  type(dft_grid_t)::grid
  type(response_cache_t),target::cache
  real(fp),allocatable,target::wa(:,:),wb(:,:),pa(:,:,:),pb(:,:,:),xa(:,:,:),xb(:,:,:)
  real(fp),allocatable::g(:,:),ref(:,:)
  real(fp)::ne,ke
  integer::n,nat,i,j,step,mode,kind,pass
  integer(i8b)::limits(3)=[0_i8b,1048576_i8b,67108864_i8b]
  infos=>oqp_handle_get_info(c_handle)
  basis=infos%basis;basis%atoms=>infos%atoms
  n=basis%nbf;nat=infos%mol_prop%natom
  allocate(wa(n,n),wb(n,n),pa(n,n,1),pb(n,n,1),xa(n,n,1),xb(n,n,1),g(3,nat),ref(3,nat))
  call dft_initialize(infos,basis,grid)
  error=0;failed=0
  do kind=1,4
    do mode=1,3
      call cache%free()
      cache%max_bytes=limits(mode)
      do step=1,3
        do pass=1,2
          wa=0;wb=0
          do i=1,n
            wa(i,i)=0.4_fp;wb(i,i)=0.3_fp
          end do
          if(step==3) wa(1,1)=0.5_fp
          do j=1,n
            do i=1,n
              pa(i,j,1)=0.01_fp*cos(real(i+j+step,fp))
              pb(i,j,1)=0.02_fp*sin(real(i+j+step,fp))
            end do
          end do
          xa=pa*0.2_fp;xb=pb*0.3_fp
          g=0;ne=0;ke=0
          ! First pass uses the cache; second is the uncached comparator.
          if(pass==1) then
            select case(kind)
            case(1)
              call utddft_xc_gradient(basis,grid,g,wa,wb,pa,pb,nMtx=1,threshold=1e-15_fp, &
                infos=infos,include_weight_derivative=.true.,cache=cache)
            case(2)
              call utddft_xc_gradient(basis,grid,g,wa,wb,pa,pb,xa,xb,1,1e-15_fp,infos,cache=cache)
            case(3)
              call tddft_xc_gradient(basis,grid,g,wa,pa,xa,1,1e-15_fp,infos,cache=cache)
            case(4)
              call derexc_blk(basis,grid,wa,wb,g,ne,ke,basis%mxam+2,n,1e-15_fp,.true.,infos,cache=cache)
            end select
            ref=g
            if(cache%bytes>limits(mode).or.cache%kernel_bytes/=0) failed=failed+1
            if(mode>1.and.step>1.and.cache%ao_hits==0) failed=failed+1
          else
            select case(kind)
            case(1)
              call utddft_xc_gradient(basis,grid,g,wa,wb,pa,pb,nMtx=1,threshold=1e-15_fp, &
                infos=infos,include_weight_derivative=.true.)
            case(2)
              call utddft_xc_gradient(basis,grid,g,wa,wb,pa,pb,xa,xb,1,1e-15_fp,infos)
            case(3)
              call tddft_xc_gradient(basis,grid,g,wa,pa,xa,1,1e-15_fp,infos)
            case(4)
              call derexc_blk(basis,grid,wa,wb,g,ne,ke,basis%mxam+2,n,1e-15_fp,.true.,infos)
            end select
            error=max(error,maxval(abs(g-ref)))
          end if
        end do
      end do
    end do
  end do
  call cache%free()
  call dftclean(infos)
  if(error>1e-10_fp) failed=failed+1
end subroutine
