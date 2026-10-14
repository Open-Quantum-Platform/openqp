! Compiled by tests/test_linalg_wrapper_exports.py from production routines.
program matrix_views
  use precision, only: dp
  use types, only: information
  use matrix_contracts
  implicit none
  type(information) :: info
  integer :: k,n,i,j,c,nocc,nvir,mode,spin,ip
  integer, parameter :: sizes(4)=[5,7,3,5]
  real(dp), parameter :: sentinel=-9182.25_dp
  real(dp), allocatable :: a(:,:),b(:,:),x(:,:),ab(:,:),bb(:,:),xb(:,:),saved(:,:)
  real(dp), allocatable :: f(:,:,:),fb(:,:,:),p(:,:),pb(:,:),v(:,:),vv(:,:),ov(:,:),ev(:)
  real(dp), allocatable :: xp(:,:),xm(:,:),z(:,:),w(:,:),wb(:,:),expected(:,:),density(:,:),packed(:),ao(:,:)
  do k=1,size(sizes)
    n=sizes(k)
    allocate(a(n,n),b(n,n),x(n,n),ab(2*n,n),bb(2*n,n),xb(2*n,n),saved(2*n,n), &
             f(11,n,n),fb(11,2*n,n),p(3*(n-1),1),pb(3*(n-1),1),ov(n,1),ev(n))
    do j=1,n
      do i=1,n
        a(i,j)=sin(real(3*i+7*j,dp))+0.1_dp*i*j
        b(i,j)=cos(real(5*i-2*j,dp))+0.03_dp*i*j
        x(i,j)=sin(real(11*i-j,dp))-0.09_dp*i*j
        do c=1,11
          f(c,i,j)=cos(real(c+7*i-3*j,dp))+0.02_dp*c*i*j
        end do
      end do
      ev(j)=0.1_dp*j-0.7_dp
    end do
    ab=sentinel;bb=sentinel;xb=sentinel
    ab(1:2*n:2,:)=a;bb(1:2*n:2,:)=b;xb(1:2*n:2,:)=x;saved=xb
    info%basis%nbf=n; info%mol_prop%nelec_a=3; info%mol_prop%nelec_b=1
    info%tddft%debug_mode=.false.
    call ov_contract(b,x,ov,n,3,3)
    call check(ov(:,1),matmul(b(:,4:),x(3,4:)),'OV MATMUL')
    call ov_contract(bb(1:2*n:2,:),xb(1:2*n:2,:),ov,n,3,3)
    call check(ov(:,1),matmul(b(:,4:),x(3,4:)),'strided OV MATMUL')
    do spin=1,3,2
      info%tddft%mult=spin
      do mode=7,11,4
        ! Density channels remain interleaved and must retain packing/scatter.
        f=0.0_dp;fb=sentinel;fb(:,1:2*n:2,:)=0.0_dp
        if(mode==7) then
          call mrsfcbc(info,a,b,x,f(:7,:,:))
          call mrsfcbc(info,ab(1:2*n:2,:),bb(1:2*n:2,:),xb(1:2*n:2,:),fb(:7,1:2*n:2,:))
        else
          call umrsfcbc(info,a,b,x,f)
          call umrsfcbc(info,ab(1:2*n:2,:),bb(1:2*n:2,:),xb(1:2*n:2,:),fb(:,1:2*n:2,:))
        end if
        call check(reshape(f(:mode,:,:),[mode*n*n]),reshape(fb(:mode,1:2*n:2,:),[mode*n*n]),'density stride')
        if(any(fb(:,2:2*n:2,:)/=sentinel)) error stop 'density scatter clobbered sentinel'
        ! Nonsymmetric signed AO channels exercise every AO->MO correction.
        do j=1,n
          do i=1,n
            do c=1,11
              f(c,i,j)=sin(real(c+3*i-5*j,dp))+0.011_dp*c*i*j
            end do
          end do
        end do
        fb=sentinel;fb(:,1:2*n:2,:)=f
        if(mode==7) then
          call mrsfmntoia(info,f(:7,:,:),p,a,b,1)
          call mrsfmntoia(info,fb(:7,1:2*n:2,:),pb,ab(1:2*n:2,:),bb(1:2*n:2,:),1)
        else
          call umrsfmntoia(info,f,p,a,b,1)
          call umrsfmntoia(info,fb(:,1:2*n:2,:),pb,ab(1:2*n:2,:),bb(1:2*n:2,:),1)
        end if
        call check(p(:,1),pb(:,1),'AO to MO stride')
      end do
      p=0.13_dp;pb=p
      call mrsfesum(info,x,a,b,p,1)
      call mrsfesum(info,xb(1:2*n:2,:),ab(1:2*n:2,:),bb(1:2*n:2,:),pb,1)
      call check(p(:,1),pb(:,1),'response sum stride')
    end do
    if(any(xb/=saved)) error stop 'borrowed response input changed'
    if(any(ab(1:2*n:2,:)/=a).or.any(bb(1:2*n:2,:)/=b)) error stop 'borrowed MO input changed'
    if(any(ab(2:2*n:2,:)/=sentinel).or.any(bb(2:2*n:2,:)/=sentinel)) error stop 'MO sentinel changed'
    nocc=2
    if(n==3) nocc=n ! Empty virtual space must never form a row offset past the storage.
    nvir=n-nocc
    allocate(xp(nocc,max(1,nvir)),xm(nocc,max(1,nvir)),z(nocc,max(1,nvir)), &
             w(n,n),wb(2*n,n),expected(n,n),density(n,n),packed(n*(n+1)/2),ao(n,n))
    do j=1,max(1,nvir)
      do i=1,nocc
        xp(i,j)=sin(real(i+3*j,dp));xm(i,j)=cos(real(2*i-j,dp));z(i,j)=0.07_dp*i-0.03_dp*j
      end do
    end do
    do mode=1,2
      density=0.0_dp
      density(:nocc,:nocc)=-matmul(xp(:,:nvir),transpose(xp(:,:nvir)))
      density(nocc+1:,nocc+1:)=matmul(transpose(xp(:,:nvir)),xp(:,:nvir))
      if(mode==2) then
        density(:nocc,:nocc)=0.5_dp*(density(:nocc,:nocc)-matmul(xm(:,:nvir),transpose(xm(:,:nvir))))
        density(nocc+1:,nocc+1:)=0.5_dp*(density(nocc+1:,nocc+1:)+matmul(transpose(xm(:,:nvir)),xm(:,:nvir)))
      end if
      ao=matmul(matmul(a,density),transpose(a))
      call tdhf_unrelaxed_density(xm,xp,a,packed,nocc,mode==1)
      ip=0
      do j=1,n
        do i=1,j
          ip=ip+1
          if(abs(packed(ip)-ao(i,j))>1e-11_dp) error stop 'unrelaxed density MATMUL mismatch'
        end do
      end do
    end do
    expected=0.0_dp
    expected(:nocc,:nocc)=0.37_dp*(matmul(xp(:,:nvir),transpose(xm(:,:nvir)))+ &
      matmul(xm(:,:nvir),transpose(xp(:,:nvir))))+a(:nocc,:nocc)
    do j=1,nvir
      do i=1,nocc
        do c=1,nocc
          expected(i,c)=expected(i,c)-ev(nocc+j)*(xp(i,j)*xp(c,j)+xm(i,j)*xm(c,j))
        end do
      end do
    end do
    do j=1,nvir
      do i=1,nvir
        do c=1,nocc
          expected(nocc+i,nocc+j)=expected(nocc+i,nocc+j)+ &
            0.37_dp*(xp(c,i)*xm(c,j)+xm(c,i)*xp(c,j))+ev(c)*(xp(c,i)*xp(c,j)+xm(c,i)*xm(c,j))
        end do
      end do
      do i=1,nocc
        expected(i,nocc+j)=dot_product(b(:nocc,i),xp(:,j))+dot_product(x(:nocc,i),xm(:,j))+z(i,j)*ev(i)
      end do
    end do
    do j=1,n
      expected(j+1:,j)=expected(j,j+1:)
    end do
    call compute_w_mo(w,xp,xm,z,a,b,x,ev,0.37_dp,nocc)
    call check(reshape(w,[n*n]),reshape(expected,[n*n]),'W MATMUL reference')
    wb=sentinel
    call compute_w_mo(wb(1:2*n:2,:),xp,xm,z,a,b,x,ev,0.37_dp,nocc)
    call check(reshape(wb(1:2*n:2,:),[n*n]),reshape(expected,[n*n]),'strided W MATMUL reference')
    if(any(wb(2:2*n:2,:)/=sentinel)) error stop 'W output clobbered sentinel'
    deallocate(a,b,x,ab,bb,xb,saved,f,fb,p,pb,ov,ev,xp,xm,z,w,wb,expected,density,packed,ao)
  end do
  print *, 'MATRIX_STORAGE_CONTRACT_PASS'
contains
  subroutine check(actual,reference,label)
    real(dp),intent(in) :: actual(:),reference(:)
    character(*),intent(in) :: label
    if(size(actual)/=size(reference)) error stop 'size mismatch'
    if(maxval(abs(actual-reference))>1e-11_dp) then
      print *,label,maxval(abs(actual-reference))
      error stop 'matrix contract mismatch'
    end if
  end subroutine
end program
