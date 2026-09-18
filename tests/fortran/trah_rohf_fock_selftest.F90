module precision
  implicit none
  integer, parameter :: dp=kind(1.0d0)
end module

module fixture
  use precision
  implicit none
  type :: trah_converger
    integer :: nbf=3,nocc_a=2,nocc_b=1,nvir_a=1,nvir_b=2,scf_type=3
    real(dp) :: fock_ao(6,2),mo_a(3,3),mo_b(3,3)
  end type
contains
  subroutine unpack_matrix(p,a)
    real(dp),intent(in) :: p(:)
    real(dp),intent(out) :: a(:,:)
    integer :: i,j,k
    k=0
    do i=1,size(a,1)
      do j=1,i
        k=k+1;a(i,j)=p(k);a(j,i)=p(k)
      end do
    end do
  end subroutine
! PRODUCTION_ROUTINES
end module

program check
  use fixture
  implicit none
  type(trah_converger) :: c
  real(dp) :: fa(3,3),fb(3,3),ka(3,3),kb(3,3),k(3,3),x(3),ha(1,2),hb(2,1)
  real(dp) :: ta(3,3),tb(3,3),original(3,3),corrected(3,3),fd(3,3),gp(3),gm(3)
  real(dp),parameter :: h=1.0e-5_dp
  integer :: i,j,l,t
  ! Stationary ROHF reference with nonzero individual spin-Fock CV blocks.
  fa=reshape([-1._dp,.2_dp,.3_dp,.2_dp,-.5_dp,0._dp,.3_dp,0._dp,.7_dp],[3,3])
  fb=reshape([-.8_dp,0._dp,-.3_dp,0._dp,.1_dp,.4_dp,-.3_dp,.4_dp,.9_dp],[3,3])
  c%mo_a=0
  do i=1,3
    c%mo_a(i,i)=1
  end do
  c%mo_b=c%mo_a
  t=0
  do i=1,3
    do j=1,i
      t=t+1;c%fock_ao(t,:)=[fa(i,j),fb(i,j)]
    end do
  end do
  do l=1,3
    x=0;x(l)=1
    call skew_sym_k(c,x,k,2)
    ka=k;ka(1:2,1:2)=0
    kb=k;kb(2:3,2:3)=0
    ta=matmul(fa,ka)-matmul(ka,fa)
    tb=matmul(fb,kb)-matmul(kb,fb)
    ha=ta(3:3,1:2);hb=tb(2:3,1:1)
    original(:,l)=2*[hb(1,1),ha(1,1)+hb(2,1),ha(1,2)]
    call rohf_missing_fock_terms(c,x,ha,hb)
    corrected(:,l)=2*[hb(1,1),ha(1,1)+hb(2,1),ha(1,2)]
    call rotated_gradient(l,h,gp)
    call rotated_gradient(l,-h,gm)
    fd(:,l)=(gp-gm)/(2*h)
  end do
  if (maxval(abs(corrected-fd))>1.e-8_dp) error stop 'corrected curvature disagrees with finite differences'
  if (maxval(abs(original-fd))<.5_dp) error stop 'fixture did not detect the omitted CO-OV coupling'
  if (maxval(abs(corrected-transpose(corrected)))>1.e-12_dp) error stop 'stationary Hessian is not symmetric'
  print *, 'PASS: ROHF Fock curvature finite differences'
contains
  subroutine rotated_gradient(l,angle,g)
    integer,intent(in) :: l
    real(dp),intent(in) :: angle
    real(dp),intent(out) :: g(3)
    real(dp) :: u(3,3),a(3,3),b(3,3)
    integer :: p,q,n
    ! Independent exact plane rotation for CO, CV, and OV, respectively.
    select case(l)
    case(1);p=1;q=2
    case(2);p=1;q=3
    case(3);p=2;q=3
    end select
    u=0
    do n=1,3
      u(n,n)=1
    end do
    u(p,p)=cos(angle);u(q,q)=cos(angle)
    u(q,p)=sin(angle);u(p,q)=-sin(angle)
    a=matmul(transpose(u),matmul(fa,u));b=matmul(transpose(u),matmul(fb,u))
    g=2*[b(2,1),a(3,1)+b(3,1),a(3,2)]
  end subroutine
end program
