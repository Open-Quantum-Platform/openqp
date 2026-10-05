! Compare signed Fock/response matrices, not just their final eigenvalues.
! The legacy serializer is the reference for this representation-only change.
module int2_shell_block_selftest_m
  use iso_c_binding, only: c_int,c_double
  use precision, only: dp
  use basis_tools, only: basis_set
  use int2_compute
  use tdhf_lib, only: int2_td_data_t
  use tdhf_mrsf_lib, only: int2_mrsf_data_t,int2_umrsf_data_t,mrsf_set_fp32
  implicit none
contains
  subroutine int2_shell_block_selftest(error, checks, failures) bind(C)
    real(c_double), intent(out) :: error
    integer(c_int), intent(out) :: checks, failures
    type(basis_set) :: basis
    class(int2_compute_data_t), allocatable :: old, new
    type(eri_data_t) :: eri
    type(int2_storage_t) :: legacy_buf, block_buf
    real(dp), allocatable,target :: den(:,:),td(:,:,:),mr(:,:,:,:),raw(:,:,:,:)
    integer :: kind,mode,trial,nthreads,pass,a,b,c,d,flip,i,j,k,l,nold,nnew,nbf,slot
    integer :: q(4),r(4),dims(4)
    integer, parameter :: perms(4,8)=reshape([1,2,3,4,2,1,3,4,1,2,4,3,2,1,4,3, &
                                            3,4,1,2,4,3,1,2,3,4,2,1,4,3,2,1],[4,8])
    real(dp) :: delta,cut
    error=0.0_dp; checks=0; failures=0
    basis%nshell=4
    basis%naos=[1,3,5,7]
    basis%ao_offset=[1,2,5,10]
    basis%nbf=sum(basis%naos);nbf=basis%nbf
    allocate(den(nbf*(nbf+1)/2,4),td(nbf,nbf,3),mr(3,11,nbf,nbf))
    do i=1,size(den,1)
      do j=1,4
        den(i,j)=sin(real(7*i+3*j,dp))
      end do
    end do
    do i=1,nbf
      do j=1,nbf
        do k=1,3
          td(i,j,k)=sin(real(7*i+3*j+k,dp))
          do l=1,11
            mr(k,l,i,j)=cos(real(3*i+7*j+5*k+l,dp))
          end do
        end do
      end do
    end do
    call mrsf_set_fp32(0)
    ! Include all spin channels, CAM's two passes, TDA, and full A+B/A-B.
    do kind=1,5
      do trial=1,7
        ! Four forced-layout cases, then the actual consumer defaults with
        ! one, two and four thread images. Do not override the latter: the
        ! default mixture of packed and direct quartets must agree as well.
        mode=min(trial,5)
        nthreads=2
        if (trial==5) nthreads=1
        if (trial==7) nthreads=4
        call mrsf_set_fp32(merge(1,0,mode==3.or.mode==4))
        select case(kind)
        case(1)
          allocate(int2_rhf_data_t::old,new)
        case(2)
          allocate(int2_urohf_data_t::old,new)
        case(3)
          allocate(int2_td_data_t::old,new)
        case(4)
          allocate(int2_mrsf_data_t::old,new)
        case(5)
          allocate(int2_umrsf_data_t::old,new)
        end select
        call attach(old);call attach(new)
        old%scale_coulomb=0.7_dp;new%scale_coulomb=0.7_dp
        old%scale_exchange=-0.3_dp;new%scale_exchange=-0.3_dp
        old%num_passes=2;new%num_passes=2
        call legacy_buf%init(17);call block_buf%init(17)
        do pass=1,2
          old%cur_pass=pass;new%cur_pass=pass
          call old%parallel_start(basis,nthreads);call new%parallel_start(basis,nthreads)
          old%shell_blocks=.false.
          if (mode<=4) then
            new%shell_blocks=.true.
            new%shell_block_min=merge(256,0,mode==2)
          end if
          nold=0;nnew=0
          do a=1,4
            do b=1,a
              do c=1,a
                do d=1,merge(b,c,a==c)
                  q=[a,b,c,d]
                  do flip=1,8
                    r=q(perms(:,flip));dims=basis%naos(r)
                    eri%ids=q;eri%flips=perms(:,flip);eri%nbf=dims
                    eri%weight=real(mode,dp);eri%weighted_cutoff=modulo(mode,2)==0
                    allocate(raw(dims(4),dims(3),dims(2),dims(1)),source=0.0_dp)
                    do i=1,dims(1)
                      do j=1,dims(2)
                        do k=1,dims(3)
                          do l=1,dims(4)
                            raw(l,k,j,i)=sin(real(3*i+7*j+5*k+11*l,dp))*0.01_dp
                          end do
                        end do
                      end do
                    end do
                    eri%pints=>raw
                    slot=1+modulo(flip,nthreads)
                    if (legacy_buf%ncur>0) call old%update(legacy_buf)
                    legacy_buf%thread_id=slot;block_buf%thread_id=slot
                    cut=merge(0.005_dp,1e-14_dp,modulo(mode,2)==0)
                    call int2_compute_data_t_storeints(old,basis,eri,legacy_buf,cut,nold)
                    call old%update(legacy_buf)
                    call int2_compute_data_t_storeints(new,basis,eri,block_buf,cut,nnew)
                    call new%update(block_buf)
                    nullify(eri%pints)
                    deallocate(raw)
                    call compare(delta)
                    error=max(error,delta)
                    checks=checks+1
                    if (delta>2e-10_dp) failures=failures+1
                  end do
                end do
              end do
            end do
          end do
          if (nold/=nnew) failures=failures+1
        end do
        call old%clean();call new%clean()
        call old%clean();call new%clean()
        call legacy_buf%clean();call block_buf%clean()
        deallocate(old,new)
      end do
    end do
    call mrsf_set_fp32(0)
  contains
    subroutine attach(x)
      class(int2_compute_data_t),intent(inout)::x
      select type(x)
      type is(int2_rhf_data_t)
        x%d=>den
      type is(int2_urohf_data_t)
        x%d=>den
      type is(int2_td_data_t)
        x%d2=>td
        x%tamm_dancoff=modulo(mode,2)==1
        x%tamm_dancoff_coulomb=.true.
        x%int_apb=.true.;x%int_amb=.true.
      type is(int2_mrsf_data_t)
        x%d3=>mr(:,1:7,:,:)
      type is(int2_umrsf_data_t)
        x%d3=>mr
      end select
    end subroutine
    subroutine compare(delta)
      real(dp),intent(out)::delta
      select type(old)
      type is(int2_rhf_data_t)
        select type(new)
        type is(int2_rhf_data_t)
          delta=maxval(abs(old%f-new%f))
        end select
      type is(int2_urohf_data_t)
        select type(new)
        type is(int2_urohf_data_t)
          delta=maxval(abs(old%f-new%f))
        end select
      type is(int2_td_data_t)
        select type(new)
        type is(int2_td_data_t)
          delta=max(maxval(abs(old%apb-new%apb)),maxval(abs(old%amb-new%amb)))
        end select
      type is(int2_mrsf_data_t)
        select type(new)
        type is(int2_mrsf_data_t)
          delta=maxval(abs(old%f3-new%f3))
          if (allocated(old%f3s)) delta=max(delta,real(maxval(abs(old%f3s-new%f3s)),dp))
        end select
      type is(int2_umrsf_data_t)
        select type(new)
        type is(int2_umrsf_data_t)
          delta=maxval(abs(old%f3-new%f3))
        end select
      end select
    end subroutine
  end subroutine
end module

! Exercise the compact inactive-image contract used by response workspace reuse.
! The packed reference retains full images, independent of that optimization.
subroutine int2_td_shell_images_selftest(error, checks, failures) bind(C)
  use iso_c_binding, only: c_double,c_int
  use precision, only: dp
  use basis_tools, only: basis_set
  use int2_compute, only: int2_storage_t,int2_shell_block_t
  use tdhf_lib, only: int2_td_data_t
  implicit none
  real(c_double), intent(out) :: error
  integer(c_int), intent(out) :: checks,failures
  type(basis_set) :: basis
  type(int2_td_data_t) :: old,new
  type(int2_storage_t) :: buf(3)
  type(int2_shell_block_t) :: block
  real(dp), target :: density(4,4,2), value(1,1,1,1)
  real(dp) :: delta
  integer :: flags,ns,t,i,j,v
  logical :: ap,am,tda
  error=0;checks=0;failures=0
  basis%nbf=4;basis%nshell=4
  basis%ao_offset=[1,2,3,4];basis%naos=[1,1,1,1]
  do v=1,2
    do j=1,4
      do i=1,4
        density(i,j,v)=sin(real(3*i+7*j+v,dp))
      end do
    end do
  end do
  old%d2=>density;new%d2=>density
  block%shells=[4,3,2,1];block%dims=1;block%offsets=[3,2,1,0]
  value=0.25_dp;block%values=>value
  do ns=1,3,2
    do flags=0,7
      ap=btest(flags,0);am=btest(flags,1);tda=btest(flags,2)
      old%int_apb=ap;new%int_apb=ap
      old%int_amb=am;new%int_amb=am
      old%tamm_dancoff=tda;new%tamm_dancoff=tda
      old%tamm_dancoff_coulomb=.true.;new%tamm_dancoff_coulomb=.true.
      old%cur_pass=1;new%cur_pass=1
      call old%parallel_start(basis,ns);call new%parallel_start(basis,ns)
      if (.not.ap) then
        deallocate(new%apb)
        allocate(new%apb(4,4,2,1),source=0.0_dp)
      end if
      if (.not.(am.or.tda)) then
        deallocate(new%amb)
        allocate(new%amb(4,4,2,1),source=0.0_dp)
      end if
      do t=1,ns
        call buf(t)%init(1)
        buf(t)%thread_id=t;buf(t)%ncur=1
        buf(t)%ids(:,1)=[4,3,2,1];buf(t)%ints(1)=value(1,1,1,1)
      end do
      !$omp parallel do num_threads(ns) default(shared) private(t) schedule(static)
      do t=1,ns
        call old%update(buf(t))
        call new%consume_shell(block,t)
      end do
      !$omp end parallel do
      delta=max(maxval(abs(sum(old%apb,dim=4)-sum(new%apb,dim=4))), &
                maxval(abs(sum(old%amb,dim=4)-sum(new%amb,dim=4))))
      error=max(error,delta);checks=checks+1
      if (delta>2e-12_dp) failures=failures+1
      if (.not.ap .and. any(new%apb/=0)) failures=failures+1
      if (.not.(am.or.tda) .and. any(new%amb/=0)) failures=failures+1
      do t=1,ns
        call buf(t)%clean()
      end do
    end do
  end do
  call old%clean();call new%clean()
  call old%clean();call new%clean()
end subroutine
