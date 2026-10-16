module mrsf_pair_selftest_m
  use iso_c_binding,only:c_double,c_int
  use precision,only:dp
  use types,only:information
  use basis_tools,only:basis_set
  use int2_compute,only:int2_storage_t
  use tdhf_mrsf_lib
  implicit none
  public :: mrsf_pair_selftest
contains
  subroutine mrsf_pair_selftest(error,checks,failures) bind(C)
    real(c_double),intent(out)::error
    integer(c_int),intent(out)::checks,failures
    type(information),target::inf
    type(basis_set),pointer::basis
    type(int2_mrsf_data_t)::old,new,stat
    type(mrsf_pair_t)::pair
    type(int2_storage_t)::buf
    real(dp),allocatable,target::den(:,:,:,:),sd(:,:,:,:)
    real(dp),allocatable::c(:,:),x(:,:),aold(:,:),anew(:,:),work(:,:),b(:,:,:),image(:,:,:)
    integer::shape,n,na,nc,dim,nf,spin,spc,workers,mode,i,j,k,l,v,idx,q,pq,rs,reuse
    real(dp)::value,delta,scales(6)
    integer,parameter::sizes(5)=[6,8,4,4,2],occup(5)=[4,5,2,4,2]
    error=0;checks=0;failures=0
    call mrsf_set_fp32(0)
    basis=>inf%basis
    inf%tddft%debug_mode=.false.;inf%tddft%umrsf=.false.;inf%control%mrsf_fp32=0
    inf%tddft%hfscale=0.6_dp
    call buf%init(1)
    do shape=1,5
      n=sizes(shape);na=occup(shape);nc=na-2;dim=na*(n-nc)
      basis%nbf=n;basis%nshell=n
      basis%naos=[(1,i=1,n)];basis%ao_offset=[(i,i=1,n)]
      inf%mol_prop%nelec_a=na;inf%mol_prop%nelec_b=nc
      allocate(c(n,n),work(n,n),b(n,n,3),sd(3,7,n,n),image(7,n,n))
      do j=1,n
        do i=1,n
          c(i,j)=sin(real(3*i+7*j+i*j,dp))/sqrt(real(n,dp))
        end do
      end do
      if(shape==1)then
        inf%tddft%mult=1;inf%tddft%pair_cache=.false.
        if(mrsf_pair_allowed(inf,c,c,.false.,.false.))failures=failures+1
        inf%tddft%pair_cache=.true.
        if(.not.mrsf_pair_allowed(inf,c,c,.false.,.false.))failures=failures+1
        if(mrsf_pair_allowed(inf,c,c,.false.,.false.,0_8))failures=failures+1
        if(mrsf_pair_allowed(inf,c,c,.true.,.false.))failures=failures+1
        if(mrsf_pair_allowed(inf,c,c,.false.,.true.))failures=failures+1
        work=c;work(1,1)=work(1,1)+0.1_dp
        if(mrsf_pair_allowed(inf,c,work,.false.,.false.))failures=failures+1
        inf%tddft%umrsf=.true.
        if(mrsf_pair_allowed(inf,c,c,.false.,.false.))failures=failures+1
        inf%tddft%umrsf=.false.;inf%control%mrsf_fp32=1
        if(mrsf_pair_allowed(inf,c,c,.false.,.false.))failures=failures+1
        inf%control%mrsf_fp32=0;inf%tddft%mult=5
        if(mrsf_pair_allowed(inf,c,c,.false.,.false.))failures=failures+1
        inf%tddft%mult=1;inf%tddft%hfscale=0
        if(mrsf_pair_allowed(inf,c,c,.false.,.false.))failures=failures+1
        inf%tddft%hfscale=0.6_dp;inf%mol_prop%nelec_b=nc+1
        if(mrsf_pair_allowed(inf,c,c,.false.,.false.))failures=failures+1
        inf%mol_prop%nelec_b=nc
      end if
      sd=0
      do j=1,n
        do i=1,n
          sd(1,7,i,j)=c(i,na-1)*c(j,na-1)
          sd(2,7,i,j)=c(i,na)*c(j,na)
          sd(3,7,i,j)=c(i,na-1)*c(j,na)
        end do
      end do
      do spin=1,3,2
        inf%tddft%mult=spin
        do spc=1,2
          scales=merge(1.0_dp,-1.0_dp,spin==1)
          if(spc==2)scales=scales*[0.7_dp,0.7_dp,0.7_dp,0.7_dp,1.3_dp,-0.4_dp]
          inf%tddft%spc_coov=0.6_dp*abs(scales(1))
          inf%tddft%spc_ovov=0.6_dp*abs(scales(5))
          inf%tddft%spc_coco=0.6_dp*merge(1.0_dp,-0.4_dp,spc==1)
          do workers=1,4,3
            do mode=1,2
              nf=merge(3,dim+1,mode==1)
              allocate(x(dim,nf),aold(dim,nf),anew(dim,nf),den(nf,7,n,n))
              do v=1,nf
                do i=1,dim
                  x(i,v)=sin(real(11*i+5*v,dp))
                end do
              end do
              if(mode==2)then
                x(:,1:dim)=0
                do i=1,dim
                  x(i,i)=1
                end do
              end if
              den=0
              do v=1,nf
                work=0
                work(1:na,nc+1:n)=reshape(x(:,v),[na,n-nc])
                call mrsfcbc(inf,c,c,work,den(v,:,:,:))
              end do
              old%d3=>den;new%d3=>den;stat%d3=>sd
              old%scale_exchange=0.6_dp;old%scale_coulomb=0.6_dp
              new%scale_exchange=0.6_dp;new%scale_coulomb=0.6_dp
              stat%scale_exchange=0.6_dp;stat%scale_coulomb=0.6_dp
              new%general_only=.true.;stat%general_only=.true.
              old%cur_pass=1;new%cur_pass=1;stat%cur_pass=1
              old%num_passes=1;new%num_passes=1;stat%num_passes=1
              call old%parallel_start(basis,workers)
              call new%parallel_start(basis,workers)
              call stat%parallel_start(basis,workers)
              if(size(new%f3,2)/=1.or.size(stat%f3,2)/=1)failures=failures+1
              if(size(new%ds)/=0.or.size(stat%ds)/=0)failures=failures+1
              if(maxval(abs(old%dsh-new%dsh))>0)failures=failures+1
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
                      buf%ncur=1;call old%update(buf)
                      buf%ncur=1;call new%update(buf)
                      buf%ncur=1;call stat%update(buf)
                    end do
                  end do
                end do
              end do
              call old%parallel_stop();call new%parallel_stop();call stat%parallel_stop()
              do k=1,3
                b(:,:,k)=-matmul(transpose(c),matmul(stat%f3(k,1,:,:,1),c))
              end do
              call mrsf_pair_matrix(pair,inf,b)
              do v=1,nf
                image=old%f3(v,:,:,:,1)
                do k=1,6
                  image(k,:,:)=image(k,:,:)*scales(k)
                end do
                call mrsfmntoia(inf,image,aold,c,c,v)
                image=0.0_dp
                image(7,:,:)=new%f3(v,1,:,:,1)
                call mrsfmntoia(inf,image,anew,c,c,v)
                block
                  real(dp) :: reference(size(anew,1))
                  reference=anew(:,v)
                  call mrsf_general_mntoia(inf,new%f3(v,1,:,:,1),anew,c,c,v)
                  error=max(error,maxval(abs(reference-anew(:,v))))
                  if(maxval(abs(reference-anew(:,v)))>2e-10_dp)failures=failures+1
                end block
              end do
              call mrsf_pair_apply(pair,x,anew)
              delta=maxval(abs(aold-anew))
              error=max(error,delta);checks=checks+1
              if(delta>2e-10_dp)failures=failures+1
              ! Re-enter without clean; retained output must be reinitialized.
              call new%parallel_start(basis,workers)
              call new%parallel_stop()
              if(maxval(abs(new%f3))>0)failures=failures+1
              call old%clean();call new%clean();call stat%clean()
              call old%clean();call new%clean();call stat%clean()
              deallocate(x,aold,anew,den)
            end do
          end do
        end do
      end do
      deallocate(c,work,b,sd,image)
    end do
    call buf%clean()
  end subroutine
end module
