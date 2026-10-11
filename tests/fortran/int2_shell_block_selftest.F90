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
    integer :: q(4),r(4),dims(4),layout_status,expected_min
    character(16) :: layout
    logical :: expected_blocks
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
    call get_environment_variable('OQP_INT2_LAYOUT',layout,status=layout_status)
    if (layout_status/=0) layout=''
    call mrsf_set_fp32(0)
    ! Include all spin channels, CAM's two passes, TDA, and full A+B/A-B.
    do kind=1,5
      do trial=1,10
        ! Four forced-layout cases, then production layout resolution with
        ! one, two and four thread images at both screening cutoffs.
        mode=trial
        nthreads=2
        if (trial>=5) then
          mode=5+modulo(trial-5,2)
          nthreads=2**((trial-5)/2)
        end if
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
          else
            call int2_resolve_layout(new)
            ! Assert the final configuration for every consumer, including
            ! RHF/UHF/TD's huge-cutoff normalization and both driver overrides.
            expected_min=merge(256,huge(0),kind>=4)
            expected_blocks=kind>=4
            select case(trim(layout))
            case('shell')
              expected_min=0
              expected_blocks=.true.
            case('legacy')
              expected_blocks=.false.
            end select
            if (new%shell_blocks.neqv.expected_blocks) failures=failures+1
            if (new%shell_block_min/=expected_min) failures=failures+1
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
          if (allocated(old%f3_private)) then
            if (size(old%f3_private)>0) delta=max(delta,maxval(abs(old%f3_private-new%f3_private)))
          end if
          if (allocated(old%f3s)) delta=max(delta,real(maxval(abs(old%f3s-new%f3s)),dp))
          if (allocated(old%f3_ex_sp)) delta=max(delta,real(maxval(abs(old%f3_ex_sp-new%f3_ex_sp)),dp))
        end select
      type is(int2_umrsf_data_t)
        select type(new)
        type is(int2_umrsf_data_t)
          delta=maxval(abs(old%f3-new%f3))
          if (allocated(old%f3_private)) then
            if (size(old%f3_private)>0) delta=max(delta,maxval(abs(old%f3_private-new%f3_private)))
          end if
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

! Independent signed-channel oracle for the shared MRSF/UMRSF accumulator.
! Dyadic, nonsymmetric inputs preserve exact FP32 arithmetic as well as FP64;
! every rank contributes a distinct nonzero factor, so a missing allreduce fails.
subroutine int2_mrsf_shared_selftest(comm, use_mpi, ranks, threads, error, checks, failures) bind(C)
  use iso_c_binding, only: c_int,c_double,c_bool
  use precision, only: dp,sp
  use parallel, only: PARALLEL_INT
  use basis_tools, only: basis_set
  use int2_compute, only: int2_storage_t,int2_shell_block_t
  use tdhf_mrsf_lib, only: int2_mrsf_data_t,int2_umrsf_data_t,mrsf_set_fp32
  implicit none
  integer(c_int), value :: comm,ranks,threads
  logical(c_bool), value :: use_mpi
  real(c_double), intent(out) :: error
  integer(c_int), intent(out) :: checks,failures
  class(int2_mrsf_data_t), allocatable :: data
  type(basis_set) :: basis
  type(int2_storage_t), allocatable :: buf(:)
  type(int2_shell_block_t), allocatable :: shell(:)
  real(dp), allocatable, target :: density(:,:,:,:),values(:,:,:,:,:)
  real(dp), allocatable :: reference(:,:,:,:)
  real(sp), allocatable :: reference_sp(:,:,:,:)
  real(sp) :: xs
  integer :: kind,layout,route,trial,pass,t,ncoul,nmatrix,nbf,nf,flag,v,c,i,j,k,l
  integer :: ao_sizes(3),nf_sizes(3),flags(3),target(2,8),read_at(2,8),a
  real(dp) :: val,xval,cval,global_factor,delta
  ao_sizes=[4,8,3]; nf_sizes=[2,3,1]; flags=[1,0,1]
  error=0;checks=0;failures=0
  allocate(buf(threads),shell(threads),values(1,1,1,1,threads))
  do kind=1,2
    if(kind==1)then
      allocate(int2_mrsf_data_t::data)
    else
      allocate(int2_umrsf_data_t::data)
    end if
    call data%pe%init(int(comm,PARALLEL_INT),use_mpi)
    if(data%pe%size/=ranks) failures=failures+1
    if(data%pe%use_mpi.neqv.use_mpi) failures=failures+1
    ncoul=merge(4,8,kind==1);nmatrix=ncoul+3
    do layout=1,3
    data%shared_exchange=layout==2
    do route=1,2
      do trial=1,3
        nbf=ao_sizes(trial);nf=nf_sizes(trial);flag=flags(trial)
        call mrsf_set_fp32(flag)
        basis%nbf=nbf;basis%nshell=nbf
        basis%naos=spread(1,1,nbf);basis%ao_offset=[(i,i=1,nbf)]
        allocate(density(nf,nmatrix,nbf,nbf),reference(nf,nmatrix,nbf,nbf),source=0.0_dp)
        do j=1,nbf
          do i=1,nbf
            do c=1,nmatrix
              do v=1,nf
                density(v,c,i,j)=real(modulo(7*i-3*j+5*c+v,31)-15,dp)/16.0_dp
              end do
            end do
          end do
        end do
        data%d3=>density
        data%num_passes=2
        data%scale_coulomb=0.5_dp;data%scale_exchange=-0.25_dp
        do pass=1,2
          data%cur_pass=pass
          call data%parallel_start(basis,int(threads))
          if(size(data%f3,2)/=nmatrix) failures=failures+1
          if(data%shared_exchange) then
            if(size(data%f3,5)/=1) failures=failures+1
            if(size(data%f3_private,2)/=ncoul.or.size(data%f3_private,5)/=threads-1) failures=failures+1
            ! Exactly ncoul*T+3 DP images: no second full result at finalization.
            if(size(data%f3)+size(data%f3_private)/=nf*(ncoul*threads+3)*nbf*nbf) failures=failures+1
          else
            if(size(data%f3,5)/=threads) failures=failures+1
            if(allocated(data%f3_private)) failures=failures+1
            if(size(data%f3)/=nf*nmatrix*threads*nbf*nbf) failures=failures+1
          end if
          if(pass==1) then
            if(allocated(data%f3s).neqv.(kind==1.and.flag==1)) failures=failures+1
            if(allocated(data%f3_ex_sp).neqv.(kind==1.and.flag==1.and.data%shared_exchange)) failures=failures+1
          end if
          ! All distinct AO indices when available; small repeated-index case too.
          i=nbf;j=merge(3,2,nbf>=4);k=2;l=1
          data%scale_exchange=merge(-0.25_dp,-0.5_dp,pass==1)
          do t=1,threads
            call buf(t)%init(1)
            buf(t)%thread_id=t;buf(t)%ncur=1
            buf(t)%ids(:,1)=[i,j,k,l]
            val=real((data%pe%rank+1)*t*(3-2*pass),dp)/8.0_dp
            buf(t)%ints(1)=val
            shell(t)%shells=[i,j,k,l];shell(t)%dims=1
            shell(t)%offsets=[i-1,j-1,k-1,l-1]
            values(1,1,1,1,t)=val;shell(t)%values=>values(:,:,:,:,t)
          end do
          !$omp parallel do num_threads(threads) default(shared) private(t) schedule(static)
          do t=1,threads
            if(route==1) then
              call data%update(buf(t))
            else
              call data%consume_shell(shell(t),t)
            end if
          end do
          !$omp end parallel do
          ! Oracle: full old image, independently assembled with all eight ERI
          ! permutations rather than the production's folded Coulomb slab.
          global_factor=real(ranks*(ranks+1)/2,dp)
          do t=1,threads
            val=global_factor*real(t*(3-2*pass),dp)/8.0_dp
            xval=val*data%scale_exchange;cval=val*data%scale_coulomb
            target=reshape([i,k,k,i,i,l,l,i,j,k,k,j,j,l,l,j],[2,8])
            read_at=reshape([j,l,l,j,j,k,k,j,i,l,l,i,i,k,k,i],[2,8])
            do c=1,nmatrix
              if(pass==2.and.c/=nmatrix) cycle
              if(pass==1.and.c<=ncoul) then
                reference(:,c,i,j)=reference(:,c,i,j)+cval*(density(:,c,k,l)+density(:,c,l,k))
                reference(:,c,j,i)=reference(:,c,j,i)+cval*(density(:,c,k,l)+density(:,c,l,k))
                reference(:,c,k,l)=reference(:,c,k,l)+cval*(density(:,c,i,j)+density(:,c,j,i))
                reference(:,c,l,k)=reference(:,c,l,k)+cval*(density(:,c,i,j)+density(:,c,j,i))
              end if
              do a=1,8
                if(kind==2.and.c>=9.and.c<=10) then
                  reference(:,c,target(1,a),target(2,a))= &
                    reference(:,c,target(1,a),target(2,a))-xval*density(:,c,read_at(2,a),read_at(1,a))
                else
                  reference(:,c,target(1,a),target(2,a))= &
                    reference(:,c,target(1,a),target(2,a))-xval*density(:,c,read_at(1,a),read_at(2,a))
                end if
              end do
            end do
            call buf(t)%clean()
            nullify(shell(t)%values)
          end do
          if(kind==1.and.flag==1.and.pass==1) then
            if(any(data%f3(:,ncoul+1:nmatrix,:,:,1)/=0.0_dp)) failures=failures+1
            if(data%shared_exchange) then
              if(.not.any(data%f3_ex_sp/=0.0_sp)) failures=failures+1
            else
              if(.not.any(data%f3s(:,ncoul+1:nmatrix,:,:,:)/=0.0_sp)) failures=failures+1
            end if
          end if
          call data%parallel_stop()
          checks=checks+1
          if(pass==1) then
            if(allocated(data%f3_private).neqv.data%shared_exchange) failures=failures+1
          else
            if(size(data%f3,2)/=nmatrix.or.size(data%f3,5)/=merge(1,int(threads),data%shared_exchange)) failures=failures+1
            if(allocated(data%f3_private).or.allocated(data%f3s).or.allocated(data%f3_ex_sp)) failures=failures+1
            if(.not.any(reference(:,nmatrix,:,:)/=0.0_dp)) failures=failures+1
            delta=maxval(abs(data%f3(:,:,:,:,1)-reference))
            error=max(error,delta)
            if(delta/=0.0_dp) failures=failures+1
          end if
        end do
        ! Deliberately re-enter without clean: next start must resize/reset all
        ! owned images and clear a previous FP32 mode. d3 remains borrowed.
        nullify(data%d3)
        deallocate(density,reference)
      end do
    end do
    if(kind==1) then
      nbf=4;nf=1;nmatrix=7
      basis%nbf=nbf;basis%nshell=nbf
      basis%naos=[1,1,1,1];basis%ao_offset=[1,2,3,4]
      allocate(density(nf,nmatrix,nbf,nbf),reference_sp(nf,3,nbf,nbf))
      do j=1,nbf
        do i=1,nbf
          do c=1,nmatrix
            density(1,c,i,j)=sin(real(7*i-3*j+5*c,dp))
          end do
        end do
      end do
      data%d3=>density;data%num_passes=1;data%cur_pass=1
      data%scale_exchange=-0.3_dp;data%scale_coulomb=0.7_dp
      call mrsf_set_fp32(1)
      i=4;j=3;k=2;l=1
      target=reshape([i,k,k,i,i,l,l,i,j,k,k,j,j,l,l,j],[2,8])
      read_at=reshape([j,l,l,j,j,k,k,j,i,l,l,i,i,k,k,i],[2,8])
      do route=1,2
        call data%parallel_start(basis,int(threads))
        call buf(1)%init(1);buf(1)%thread_id=1;buf(1)%ncur=1
        buf(1)%ids(:,1)=[i,j,k,l];buf(1)%ints(1)=0.13_dp
        shell(1)%shells=[i,j,k,l];shell(1)%dims=1
        shell(1)%offsets=[i-1,j-1,k-1,l-1]
        values(1,1,1,1,1)=0.13_dp;shell(1)%values=>values(:,:,:,:,1)
        reference_sp=0.0_sp
        xs=real(buf(1)%ints(1)*data%scale_exchange,sp)
        do c=5,7
          do a=1,8
            reference_sp(:,c-4,target(1,a),target(2,a))= &
              reference_sp(:,c-4,target(1,a),target(2,a))-xs*real(density(:,c,read_at(1,a),read_at(2,a)),sp)
          end do
        end do
        if(route==1) then
          call data%update(buf(1))
        else
          call data%consume_shell(shell(1),1)
        end if
        if(any(data%f3(:,ncoul+1:nmatrix,:,:,1)/=0.0_dp)) failures=failures+1
        if(data%shared_exchange) then
          delta=real(maxval(abs(data%f3_ex_sp-reference_sp)),dp)
        else
          delta=real(maxval(abs(data%f3s(:,5:7,:,:,1)-reference_sp)),dp)
        end if
        error=max(error,delta);checks=checks+1
        if(delta/=0.0_dp) failures=failures+1
        call data%parallel_stop()
        call buf(1)%clean();nullify(shell(1)%values)
      end do
      nullify(data%d3);deallocate(density,reference_sp)
    end if
    end do ! storage layouts, including transition without clean
    call data%clean()
    ! Leave a CAM calculation after its first pass, while the borrowed density
    ! is still live. Cleanup owns only the accumulator images; it must neither
    ! overwrite nor release the caller's density. Re-enter with private storage
    ! and FP32 disabled after repeated cleanup of the shared FP32 case.
    nbf=2;nf=1;nmatrix=ncoul+3
    basis%nbf=nbf;basis%nshell=nbf
    basis%naos=[1,1];basis%ao_offset=[1,2]
    allocate(density(nf,nmatrix,nbf,nbf),reference(nf,nmatrix,nbf,nbf))
    do j=1,nbf
      do i=1,nbf
        do c=1,nmatrix
          density(1,c,i,j)=real(7*i-3*j+5*c,dp)/16.0_dp
        end do
      end do
    end do
    reference=density
    do route=1,2
      data%d3=>density
      data%shared_exchange=route==1
      data%num_passes=2;data%cur_pass=1
      call mrsf_set_fp32(merge(1,0,route==1))
      call data%parallel_start(basis,int(threads))
      if(.not.associated(data%d3,density)) failures=failures+1
      if(allocated(data%f3s).neqv.(kind==1.and.route==1)) failures=failures+1
      call data%parallel_stop() ! non-final pass preserves owned images
      if(.not.allocated(data%f3)) failures=failures+1
      if(any(density/=reference)) failures=failures+1
      call data%clean();call data%clean()
      if(associated(data%d3)) failures=failures+1
      if(allocated(data%f3).or.allocated(data%f3_private).or.allocated(data%f3_ex_sp).or. &
         allocated(data%f3s).or.allocated(data%ds).or.allocated(data%ds_sp).or. &
         allocated(data%d3_sp).or.allocated(data%dsh)) failures=failures+1
      if(.not.allocated(density)) then
        failures=failures+1
      else
        if(any(density/=reference)) failures=failures+1
      end if
    end do
    deallocate(density,reference)
    deallocate(data)
  end do
  call mrsf_set_fp32(0)
end subroutine int2_mrsf_shared_selftest
