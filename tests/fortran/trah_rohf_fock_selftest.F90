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
  real(dp) :: ei(3),ej(3)
  real(dp),parameter :: h=1.0e-4_dp
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
    ej=0;ej(l)=h
    do i=1,3
      ei=0;ei(i)=h
      fd(i,l)=(energy(ei+ej)-energy(ei-ej)-energy(-ei+ej)+energy(-ei-ej))/(4*h*h)
    end do
  end do
  if (maxval(abs(corrected-fd))>2.e-6_dp) error stop 'corrected curvature disagrees with finite differences'
  if (maxval(abs(original-fd))<.1_dp) error stop 'fixture did not detect the omitted CO-OV coupling'
  if (maxval(abs(corrected-transpose(corrected)))>1.e-12_dp) error stop 'Hessian is not symmetric'
  print *, 'PASS: ROHF Fock curvature finite differences'
contains
  function energy(x) result(e)
    real(dp),intent(in) :: x(3)
    real(dp) :: e,u(3,3),term(3,3),rot(3,3),a(3,3),b(3,3)
    integer :: n
    ! Independent exponential-coordinate energy, not the moving gradient.
    rot=0
    rot(2,1)=x(1);rot(3,1)=x(2);rot(3,2)=x(3)
    rot=rot-transpose(rot)
    u=0
    do n=1,3
      u(n,n)=1
    end do
    term=u
    do n=1,8
      term=matmul(term,rot)/n
      u=u+term
    end do
    a=matmul(transpose(u),matmul(fa,u));b=matmul(transpose(u),matmul(fb,u))
    e=a(1,1)+a(2,2)+b(1,1)
  end function
end program
