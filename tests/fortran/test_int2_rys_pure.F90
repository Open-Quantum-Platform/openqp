program test_int2_rys_pure
  ! Compare the contracted spherical path against independently evaluated
  ! Cartesian ERIs followed by the established cart2sph transformation.
  use precision, only: dp
  use basis_tools, only: basis_set
  use atomic_structure_m, only: atomic_structure
  use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
  use int2e_rys, only: int2_rys_data_t, int2_rys_compute, int2_rys_reduce_pure
  use constants, only: num_cart_bf, shells_pnrm2, HARMONIC_ACTIVE
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(basis_set) :: basis
  type(atomic_structure), target :: atoms
  type(int2_pair_storage) :: pairs
  type(int2_cutoffs_t) :: cuts
  type(int2_rys_data_t) :: gd
  integer, parameter :: classes(4,8)=reshape([ &
    0,0,3,3, 3,3,3,3, 0,1,2,3, 4,4,4,4, &
    1,3,0,2, 0,0,5,6, 2,0,0,0, 0,1,0,0 ],[4,8])
  integer, parameter :: masks(4)=[0,1,9,15], contractions(3)=[1,3,6]
  integer :: c,nc,s,j,k,st,mode,np,nbf(4),nbf_ref(4),mu,pass,mask,ngsaved,ncase
  real(dp), allocatable :: ints(:), ref(:)
  real(dp) :: err,worst,mu2
  logical :: zero

  allocate(ints(15**4),ref(15**4))
  allocate(atoms%xyz(3,4))
  atoms%xyz(:,1)=[0.0_dp,0.1_dp,0.0_dp]
  atoms%xyz(:,2)=[0.8_dp,0.2_dp,0.3_dp]
  atoms%xyz(:,3)=[-0.3_dp,0.9_dp,0.2_dp]
  atoms%xyz(:,4)=[0.4_dp,-0.6_dp,0.7_dp]
  basis%atoms=>atoms
  basis%nshell=4
  allocate(basis%am(4),basis%harmonic(4),basis%origin(4),basis%ncontr(4),basis%g_offset(4))
  allocate(basis%shell_centers(4,3))
  basis%shell_centers=transpose(atoms%xyz)
  basis%origin=[1,2,3,4]
  call cuts%set(1e-14_dp,1e-15_dp,1e-15_dp,40.0_dp)
  worst=0
  ncase=0
  do pass=1,2
    ! Reinitialize the same object after cleanup, including projection caches.
    call gd%init(6,cuts,st)
    if(st/=0) error stop 'allocation failed'
    HARMONIC_ACTIVE=pass==1
    do nc=1,3
      np=contractions(nc)
      basis%nprim=4*np
      basis%ncontr=np
      allocate(basis%ex(4*np),basis%cc(4*np))
      do s=1,4
        basis%g_offset(s)=1+(s-1)*np
        do j=1,np
          k=(s-1)*np+j
          basis%ex(k)=0.15_dp+0.2_dp*j+0.04_dp*s
          basis%cc(k)=(-1.0_dp)**j/np
        end do
      end do
      call pairs%alloc(basis,cuts)
      do c=1,size(classes,2)
        ! Alternating large/small and differently ordered quartets exercises
        ! buffer growth, reuse, dimension changes and shell permutation.
        basis%am=classes(:,c)
        call pairs%compute(basis,cuts)
        do mask=1,size(masks)
          do s=1,4
            basis%harmonic(s)=merge(1,0,btest(masks(mask),s-1))
          end do
          do mu=0,1
            mu2=0.16_dp
            do mode=0,1
              call gd%set_ids(basis,[1,2,3,4])
              if(mu==0) then
                call int2_rys_compute(ints,gd,pairs,zero,basis=basis,direct_pure=mode==1)
              else
                call int2_rys_compute(ints,gd,pairs,zero,mu2=mu2,basis=basis,direct_pure=mode==1)
              end if
              if(zero) error stop 'unexpected screened quartet'
              nbf=gd%nbf
              if(.not.gd%direct_pure) then
                call normalize_cart()
                call int2_rys_reduce_pure(basis,gd,ints,nbf)
              end if
              if(.not.all(ieee_is_finite(ints(1:product(nbf))))) error stop 'nonfinite ERI'
              if(mode==0) then
                nbf_ref=nbf
                ref(1:product(nbf))=ints(1:product(nbf))
              else
                if(any(nbf/=nbf_ref)) error stop 'spherical dimension mismatch'
                err=maxval(abs(ref(1:product(nbf))-ints(1:product(nbf))))
                worst=max(worst,err)
                if(err>2e-11_dp*max(1.0_dp,maxval(abs(ref(1:product(nbf)))))) then
                  print *, 'FAIL',pass,np,c,mask,mu,err
                  error stop 'spherical ERI mismatch'
                end if
                ncase=ncase+1
              end if
            end do
          end do
        end do
        ! An entirely screened pair must not expose a previous contraction.
        ngsaved=pairs%ppid(1,2)
        pairs%ppid(1,2)=0
        call gd%set_ids(basis,[1,2,3,4])
        call int2_rys_compute(ints,gd,pairs,zero,basis=basis,direct_pure=.true.)
        if(.not.zero) error stop 'screened quartet not reported'
        pairs%ppid(1,2)=ngsaved
      end do
      call pairs%clean()
      deallocate(basis%ex,basis%cc)
    end do
    call gd%clean()
    if(allocated(gd%proj_cache)) error stop 'projection cache not released'
    call gd%clean()
  end do
  HARMONIC_ACTIVE=.true.
  print *, 'RYS_PURE_PASS',ncase,'max absolute error',worst
contains
  subroutine normalize_cart()
    integer :: i,j,k,l,n
    n=0
    do i=1,gd%nbf(1)
      do j=1,gd%nbf(2)
        do k=1,gd%nbf(3)
          do l=1,gd%nbf(4)
            n=n+1
            ints(n)=ints(n)*shells_pnrm2(i,gd%am(1))*shells_pnrm2(j,gd%am(2))* &
              shells_pnrm2(k,gd%am(3))*shells_pnrm2(l,gd%am(4))
          end do
        end do
      end do
    end do
  end subroutine normalize_cart
end program test_int2_rys_pure
