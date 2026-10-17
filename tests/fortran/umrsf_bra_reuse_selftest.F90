! Pin reused images against independent production bra-density construction
! and the original eleven-channel contraction, with unequal spin orbitals.
module umrsf_bra_reuse_selftest_m
  use iso_c_binding, only: c_double, c_int
  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_storage_t
  use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc, mrsf_set_fp32
  use tdhf_umrsf_gradient_mod, only: umrsf_bra_density, umrsf_bra_image
  implicit none
  public :: umrsf_bra_reuse_selftest
contains
  subroutine umrsf_bra_reuse_selftest(error, checks) bind(C)
    real(c_double), intent(out) :: error
    integer(c_int), intent(out) :: checks
    type(information) :: inf
    type(int2_umrsf_data_t) :: ket, bra
    type(int2_storage_t) :: buf
    real(dp), allocatable, target :: d(:,:,:,:), b(:,:,:,:)
    real(dp), allocatable :: va(:,:), vb(:,:), x(:,:), f(:,:,:), out(:,:)
    integer, parameter :: sizes(4)=[6,4,4,2], occ(4)=[4,2,4,2]
    integer :: shape,n,na,spin,workers,v,i,j,k,l,pq,rs,q,c,same
    real(dp) :: value, sign
    error=0; checks=0
    call mrsf_set_fp32(0)
    call buf%init(1)
    inf%tddft%debug_mode=.false.
    do shape=1,4
      n=sizes(shape);na=occ(shape)
      inf%basis%nbf=n;inf%basis%nshell=n
      inf%basis%naos=[(1,i=1,n)];inf%basis%ao_offset=[(i,i=1,n)]
      inf%mol_prop%nelec_a=na;inf%mol_prop%nelec_b=na-2
      allocate(d(3,11,n,n),b(3,11,n,n),va(n,n),vb(n,n),x(n,n),f(n,n,11),out(n,n))
      do same=0,1
        do j=1,n
          do i=1,n
            va(i,j)=sin(real(3*i+7*j+i*j,dp))/sqrt(real(n,dp))
            vb(i,j)=cos(real(5*i+2*j+i*j,dp))/sqrt(real(n,dp))
          end do
        end do
        if(same==1)vb=va
        do spin=1,3,2
          inf%tddft%mult=spin
          d=0;b=0
          do v=1,3
            x=0
            do j=na-1,n
              do i=1,na
                x(i,j)=sin(real(11*i+5*j+13*v,dp))
              end do
            end do
            call umrsfcbc(inf,va,vb,x,d(v,:,:,:))
            call umrsf_bra_density(inf,va,vb,x,b(v,:,:,:))
          end do
          do workers=1,4,3
            ket%d3=>d;bra%d3=>b
            ket%scale_exchange=0.6_dp;bra%scale_exchange=0.6_dp
            ket%scale_coulomb=0.6_dp;bra%scale_coulomb=0.6_dp
            ket%tamm_dancoff=.true.;bra%tamm_dancoff=.true.
            ket%cur_pass=1;bra%cur_pass=1;ket%num_passes=1;bra%num_passes=1
            call ket%parallel_start(inf%basis,workers)
            call bra%parallel_start(inf%basis,workers)
            q=0
            do i=1,n
              do j=1,i
                pq=i*(i-1)/2+j
                do k=1,i
                  do l=1,k
                    rs=k*(k-1)/2+l
                    if(rs>pq)cycle
                    value=sin(real(7*pq+13*rs,dp))*0.01_dp
                    if(i==j)value=value*0.5_dp
                    if(k==l)value=value*0.5_dp
                    if(pq==rs)value=value*0.5_dp
                    q=q+1;buf%thread_id=1+modulo(q,workers)
                    buf%ids(:,1)=[i,j,k,l];buf%ints(1)=value
                    buf%ncur=1;call ket%update(buf)
                    buf%ncur=1;call bra%update(buf)
                  end do
                end do
              end do
            end do
            call ket%parallel_stop();call bra%parallel_stop()
            do v=1,3
              do c=1,11
                sign=1.0_dp
                if(spin==3.and.c<=10)sign=-1.0_dp
                f(:,:,c)=sign*ket%f3(v,c,:,:,1)
              end do
              do c=1,11
                sign=1.0_dp
                if(spin==3.and.c<=10)sign=-1.0_dp
                call umrsf_bra_image(f,c,out)
                error=max(error,maxval(abs(out-sign*bra%f3(v,c,:,:,1))))
                checks=checks+1
              end do
            end do
            call ket%clean();call bra%clean()
          end do
        end do
      end do
      deallocate(d,b,va,vb,x,f,out)
    end do
    call buf%clean()
  end subroutine
end module
