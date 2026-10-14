! Exercise the production context, both storage modes, and repeated ownership.
program casscf_anhess_selftest
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use casscf_anhess_mod
  implicit none
  integer, parameter :: i8=c_int64_t, dp=c_double
  type(cas_anhess_ctx_t) :: ctx
  integer :: na, n, nd, mode, pass, i, j, p, q, r, s, a, b, pos, expected
  integer :: threads, roots(0:0)
  integer(i8) :: status, tiny(0:0)
  integer(i8), allocatable :: dets(:)
  integer(c_int32_t), allocatable :: pairs(:)
  real(dp), allocatable :: h1(:), eri(:), ci(:), hess(:), ref(:), first(:)
  real(dp) :: weights(0:0), factor, err, scale
  character(len=16) :: arg

  call get_command_argument(1, arg)
  read(arg, *) threads
  weights=1.0_dp
  roots=0
  tiny=0_i8
  do na=5,6
    n=na+1
    nd=merge(100,400,na==5)
    expected=merge(0,1,na==5)
    allocate(dets(0:nd-1), pairs(0:2*na-1), h1(0:n*n-1), &
             eri(0:n**4-1), ci(0:nd-1), hess(0:na*na-1), &
             ref(0:na*na-1), first(0:na*na-1))
    pos=0
    do a=0,2**na-1
      if (popcnt(a)/=3) cycle
      do b=0,2**na-1
        if (popcnt(b)/=na-3) cycle
        dets(pos)=ior(int(a,i8),ishft(int(b,i8),na))
        ci(pos)=sin(real(pos+1,dp))
        pos=pos+1
      end do
    end do
    if (pos/=nd) error stop 'determinant extent'
    ci=ci/sqrt(sum(ci*ci))
    h1=0.0_dp
    do i=0,na-1
      pairs(2*i:2*i+1)=[int(na,c_int32_t),int(i,c_int32_t)]
    end do
    do p=0,n-1
      do q=0,n-1
        h1(p*n+q)=0.01_dp*sin(real((p+1)*(q+1),dp))
        if (p==q) h1(p*n+q)=h1(p*n+q)+real(2**p,dp)
        do r=0,n-1
          do s=0,n-1
            factor=0.001_dp*cos(real((p+1)*(q+1),dp))*cos(real((r+1)*(s+1),dp))
            eri(((p*n+q)*n+r)*n+s)=factor
          end do
        end do
      end do
    end do
    do mode=0,2
      do pass=1,2
        status=cas_anhess_init(ctx,1,0,1,0,0,1,1,threads,1_i8,tiny)
        if (status/=0_i8) error stop 'small context init'
        if (mode==2) then
          status=cas_anhess_init(ctx,n,0,na,3,na-3,na,1,threads,int(nd,i8),dets)
        else
          status=cas_anhess_init(ctx,n,0,na,3,na-3,na,1,threads,int(nd,i8),dets,mode==1)
        end if
        if (status/=0_i8) error stop 'target context init'
        if (allocated(ctx%stack) .eqv. ctx%use_mf) error stop 'stack ownership'
        if (mode==2 .and. (merge(1,0,ctx%use_mf)/=expected)) error stop 'automatic mode'
        status=cas_anhess_build(ctx,pairs,weights,roots,1,dets,h1,eri,ci,hess)
        if (status/=0_i8) error stop 'Hessian build'
        if (ctx%use_mf) then
          if (any(ctx%mfperm/=ctx%sperm-1_i8)) error stop 'permutation origin'
          if (any(ctx%mfperm<0_i8) .or. any(ctx%mfperm>=int(nd,i8))) error stop 'permutation bounds'
        end if
        if (pass==1) then
          first=hess
        else
          if (any(first/=hess)) error stop 'repeat changed Hessian'
        end if
      end do
      if (mode==0) ref=hess
      err=maxval(abs(hess-ref))
      scale=max(1.0_dp,maxval(abs(ref)))
      if (err>1.0e-12_dp*scale) error stop 'direct vs dense Hessian'
      print *, 'nact=',na,' ndet=',nd,' threads=',threads,' mode=',mode,' maxdiff=',err
    end do
    status=cas_anhess_init(ctx,33,0,32,0,0,1,1,threads,1_i8,tiny)
    if (status/=CAS_HESS_ERR_INPUT .or. allocated(ctx%wmat)) error stop 'invalid init cleanup'
    call cas_anhess_free(ctx)
    call cas_anhess_free(ctx)
    deallocate(dets,pairs,h1,eri,ci,hess,ref,first)
  end do
end program casscf_anhess_selftest
