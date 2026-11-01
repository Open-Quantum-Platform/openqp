subroutine test_int2_rotaxis_pure() bind(C, name="oqp_test_int2_rotaxis_pure")
  ! Compare fused rotation/projection with Cartesian rotation followed by
  ! the established spherical transformation, for all s/p/d shell orders.
  use precision, only: dp
  use basis_tools, only: basis_set
  use atomic_structure_m, only: atomic_structure
  use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
  use int2e_rys, only: int2_rys_data_t, int2_rys_compute
  use int2e_rotaxis, only: genr22, genr22_pure, genr22_reduce_pure
  use constants, only: num_cart_bf, HARMONIC_ACTIVE
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(basis_set) :: basis
  type(atomic_structure), target :: atoms
  type(int2_pair_storage) :: pairs
  type(int2_cutoffs_t) :: cuts
  type(int2_rys_data_t) :: gd
  logical :: zero
  integer, parameter :: contractions(3)=[1,3,6]
  integer :: geom,nc,np,s,j,k,c,mask,mu,nbf(4),nbf_ref(4),flips(4),flips_ref(4),ncase,nsaved,stat
  real(dp) :: ints(1296), ref(1296), ref_rys(1296), err, worst, worst_rys

  allocate(atoms%xyz(3,4))
  basis%atoms=>atoms
  basis%nshell=4
  allocate(basis%am(4),basis%harmonic(4),basis%origin(4),basis%ncontr(4),basis%g_offset(4))
  allocate(basis%shell_centers(4,3))
  basis%origin=[1,2,3,4]
  HARMONIC_ACTIVE=.true.
  call cuts%set(1e-14_dp,1e-15_dp,1e-15_dp,40.0_dp)
  call gd%init(1,cuts,stat)
  if(stat/=0) error stop 'Rys reference allocation failed'
  worst_rys=0
  worst=0
  ncase=0
  do geom=1,6
    select case(geom)
    case(1)
      atoms%xyz(:,1)=[0.0_dp,0.1_dp,0.0_dp]
      atoms%xyz(:,2)=[0.8_dp,0.2_dp,0.3_dp]
      atoms%xyz(:,3)=[-0.3_dp,0.9_dp,0.2_dp]
      atoms%xyz(:,4)=[0.4_dp,-0.6_dp,0.7_dp]
    case(2)
      ! Coincident centers exercise identity and indeterminate axes.
      atoms%xyz=0
    case(3)
      ! Collinear centers exercise the parallel-axis construction.
      atoms%xyz=0
      atoms%xyz(3,:)=[0.0_dp,0.8_dp,1.9_dp,2.4_dp]
    case(4)
      ! Well-separated shell pairs exercise the asymptotic Boys evaluation.
      atoms%xyz(:,1)=[0.0_dp,0.1_dp,0.0_dp]
      atoms%xyz(:,2)=[0.8_dp,0.2_dp,0.3_dp]
      atoms%xyz(:,3)=[70.0_dp,0.9_dp,0.2_dp]
      atoms%xyz(:,4)=[70.4_dp,-0.6_dp,0.7_dp]
    case(5)
      atoms%xyz=0
      atoms%xyz(3,:)=[0.0_dp,0.8_dp,1.9_dp,2.4_dp]
      atoms%xyz(1,3)=1e-11_dp
      atoms%xyz(2,4)=-2e-11_dp
    case(6)
      ! A translated nonplanar quartet checks product-center differences.
      atoms%xyz(:,1)=[10.0_dp,-19.9_dp,30.0_dp]
      atoms%xyz(:,2)=[10.8_dp,-19.8_dp,30.3_dp]
      atoms%xyz(:,3)=[9.7_dp,-19.1_dp,30.2_dp]
      atoms%xyz(:,4)=[10.4_dp,-20.6_dp,30.7_dp]
    end select
    basis%shell_centers=transpose(atoms%xyz)
    do nc=1,size(contractions)
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
      do c=0,80
        do s=1,4
          basis%am(s)=mod(c/3**(s-1),3)
        end do
        call pairs%compute(basis,cuts)
        do mask=0,15
          do s=1,4
            basis%harmonic(s)=merge(1,0,btest(mask,s-1))
          end do
          do mu=0,1
            if(mu==0) then
              call genr22(basis,pairs,ref,[1,2,3,4],flips_ref,cuts)
              call genr22_pure(basis,pairs,ints,[1,2,3,4],flips,cuts,nbf)
            else
              call genr22(basis,pairs,ref,[1,2,3,4],flips_ref,cuts,0.16_dp)
              call genr22_pure(basis,pairs,ints,[1,2,3,4],flips,cuts,nbf,0.16_dp)
            end if
            nbf_ref=num_cart_bf(basis%am(flips_ref))
            call genr22_reduce_pure(basis,[1,2,3,4],flips_ref,ref,nbf_ref)
            if(any(flips/=flips_ref).or.any(nbf/=nbf_ref)) error stop 'shell order or dimension mismatch'
            if(.not.all(ieee_is_finite(ints(:product(nbf)))) ) error stop 'nonfinite ERI'
            err=maxval(abs(ints(:product(nbf))-ref(:product(nbf))))
            worst=max(worst,err)
            if(err>2e-11_dp*max(1.0_dp,maxval(abs(ref(:product(nbf)))))) then
              print *, 'FAIL',geom,np,c,mask,mu,err
              error stop 'rotated-axis spherical ERI mismatch'
            end if
            if (sum(basis%am)==0 .or. (np==1 .and. sum(basis%am)==1)) then
              ! The unrotated s/p formulas must also agree with independent Rys
              ! quadrature, including components close to the axis cutoff.
              call gd%set_ids(basis,[1,2,3,4])
              if(mu==0) then
                call int2_rys_compute(ref_rys,gd,pairs,zero)
              else
                call int2_rys_compute(ref_rys,gd,pairs,zero,mu2=0.16_dp)
              end if
              if(zero) error stop 'unexpected screened Rys reference'
              err=maxval(abs(ref(:product(nbf_ref))-ref_rys(:product(nbf_ref))))
              if(.not.all(ieee_is_finite(ref_rys(:product(nbf_ref))))) &
                error stop 'nonfinite Rys reference'
              if(geom==5.and.sum(basis%am)==1) then
                if(maxval(abs(ref(:2)-ref_rys(:2)))>1e-20_dp) &
                  error stop 'small transverse s/p components lost'
              end if
              worst_rys=max(worst_rys,err)
              ! Both paths retain their existing, different Boys approximations.
              if(err>2e-11_dp*max(1.0_dp,maxval(abs(ref_rys(:product(nbf_ref)))))) then
                print *, 'RYS_FAIL',geom,np,c,mask,mu,err,ref(:product(nbf_ref)),ref_rys(:product(nbf_ref))
                error stop 'unrotated s/p integral disagrees with Rys quadrature'
              end if
            end if
            ncase=ncase+1
          end do
        end do
        if (c == 0) then
          ! An empty bra pair must overwrite the scalar output with zero.
          nsaved=pairs%ppid(1,2)
          pairs%ppid(1,2)=0
          ints=123.0_dp
          call genr22(basis,pairs,ints,[1,2,3,4],flips,cuts)
          if(ints(1)/=0.0_dp) error stop 'screened ssss output not cleared'
          pairs%ppid(1,2)=nsaved
        end if
      end do
      call pairs%clean()
      deallocate(basis%ex,basis%cc)
    end do
  end do
  call gd%clean()
  print *, 'SP_RYS_REFERENCE',worst_rys
  print *, 'ROTAXIS_PURE_PASS',ncase,'max absolute error',worst
end subroutine test_int2_rotaxis_pure
