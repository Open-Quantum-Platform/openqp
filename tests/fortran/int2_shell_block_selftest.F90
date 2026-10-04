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
    integer :: kind,mode,pass,a,b,c,d,flip,i,j,k,l,nold,nnew,nbf,slot
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
      do mode=1,4
        call mrsf_set_fp32(merge(1,0,mode>2))
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
          call old%parallel_start(basis,2);call new%parallel_start(basis,2)
          old%shell_blocks=.false.;new%shell_blocks=.true.
          new%shell_block_min=merge(256,0,mode==2)
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
                    slot=1+modulo(flip,2)
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
