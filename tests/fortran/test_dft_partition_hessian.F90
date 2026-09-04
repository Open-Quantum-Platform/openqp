program test_dft_partition_hessian
  use precision, only: fp
  use mod_dft_partfunc, only: PTYPE_SSF, PTYPE_BECKE4
  use mod_dft_partition_hessian, only: partition_weight_nuclear_derivatives, &
    partition_weight_nuclear_first_derivatives
  implicit none

  integer, parameter :: nat=3, nq=3*nat
  real(fp), parameter :: h=2.0e-4_fp
  real(fp) :: xyz0(3,nat), point0(3), shift(nat,nat)
  real(fp) :: w(nat), dw(3,nat,nat), d2w(3,nat,3,nat,nat)
  real(fp) :: w_first(nat),dw_first(3,nat,nat)
  logical :: dummy(nat)
  integer :: ptype, status

  xyz0 = reshape([0.0_fp,0.0_fp,0.0_fp, 1.4_fp,0.2_fp,-0.1_fp, &
                  -0.3_fp,1.2_fp,0.4_fp],[3,nat])
  point0 = xyz0(:,2)+[0.31_fp,-0.27_fp,0.19_fp]
  shift = 0.0_fp
  shift(1,2)=0.08_fp; shift(1,3)=-0.05_fp; shift(2,3)=0.04_fp
  dummy = .false.

  do ptype = PTYPE_SSF, PTYPE_BECKE4, PTYPE_BECKE4-PTYPE_SSF
    call partition_weight_nuclear_derivatives(xyz0,point0,2,dummy,ptype,.true., &
      shift,w,dw,d2w,status)
    if (status /= 0) error stop 'analytic partition derivative failed'
    call partition_weight_nuclear_first_derivatives(xyz0,point0,2,dummy,ptype, &
      .true.,shift,w_first,dw_first,status)
    if(status/=0 .or. maxval(abs(w_first-w))>2.0e-14_fp .or. &
       maxval(abs(dw_first-dw))>2.0e-13_fp) &
      error stop 'first-order partition path differs from full derivative'
    call check_finite_difference(ptype,xyz0,point0,shift,dw,d2w)
    if (abs(sum(w)-1.0_fp) > 2.0e-13_fp) error stop 'weights do not normalize'
    if (maxval(abs(sum(dw,dim=3))) > 2.0e-12_fp) &
      error stop 'first derivatives do not normalize'
    if (maxval(abs(sum(d2w,dim=5))) > 2.0e-11_fp) &
      error stop 'second derivatives do not normalize'
    if (maxval(abs(sum(dw,dim=2))) > 2.0e-12_fp) &
      error stop 'first derivative translation sum rule failed'
  end do
  print '(a)', 'dft partition Hessian selftest passed'

contains

  subroutine displaced_weights(ptype,xyz0,point0,shift,q1,s1,q2,s2,wout)
    integer, intent(in) :: ptype,q1,s1,q2,s2
    real(fp), intent(in) :: xyz0(3,nat),point0(3),shift(nat,nat)
    real(fp), intent(out) :: wout(nat)
    real(fp) :: xyz(3,nat),point(3),dw0(3,nat,nat),d20(3,nat,3,nat,nat)
    integer :: atom,comp,st
    xyz=xyz0
    atom=(q1-1)/3+1; comp=mod(q1-1,3)+1
    xyz(comp,atom)=xyz(comp,atom)+real(s1,fp)*h
    if (q2 > 0) then
      atom=(q2-1)/3+1; comp=mod(q2-1,3)+1
      xyz(comp,atom)=xyz(comp,atom)+real(s2,fp)*h
    end if
    point=point0+(xyz(:,2)-xyz0(:,2))
    call partition_weight_nuclear_derivatives(xyz,point,2,dummy,ptype,.true., &
      shift,wout,dw0,d20,st)
    if (st /= 0) error stop 'displaced partition evaluation failed'
  end subroutine displaced_weights

  subroutine check_finite_difference(ptype,xyz,point,shift,dw0,d20)
    integer, intent(in) :: ptype
    real(fp), intent(in) :: xyz(3,nat),point(3),shift(nat,nat)
    real(fp), intent(in) :: dw0(3,nat,nat),d20(3,nat,3,nat,nat)
    real(fp) :: wp(nat),wm(nat),wpp(nat),wpm(nat),wmp(nat),wmm(nat)
    real(fp) :: fd1,fd2,max1,max2
    integer :: q,r,a,b,c,ra
    max1=0.0_fp; max2=0.0_fp
    do q=1,nq
      a=(q-1)/3+1; b=mod(q-1,3)+1
      call displaced_weights(ptype,xyz,point,shift,q,1,0,0,wp)
      call displaced_weights(ptype,xyz,point,shift,q,-1,0,0,wm)
      do c=1,nat
        fd1=(wp(c)-wm(c))/(2.0_fp*h)
        max1=max(max1,abs(fd1-dw0(b,a,c)))
      end do
      do r=1,nq
        ra=int(ceiling(real(r,fp)/3.0_fp))
        call displaced_weights(ptype,xyz,point,shift,q,1,r,1,wpp)
        call displaced_weights(ptype,xyz,point,shift,q,1,r,-1,wpm)
        call displaced_weights(ptype,xyz,point,shift,q,-1,r,1,wmp)
        call displaced_weights(ptype,xyz,point,shift,q,-1,r,-1,wmm)
        do c=1,nat
          fd2=(wpp(c)-wpm(c)-wmp(c)+wmm(c))/(4.0_fp*h*h)
          max2=max(max2,abs(fd2-d20(b,a,mod(r-1,3)+1,ra,c)))
        end do
      end do
    end do
    if (max1 > 2.0e-7_fp) error stop 'partition first derivative FD mismatch'
    if (max2 > 3.0e-5_fp) error stop 'partition second derivative FD mismatch'
  end subroutine check_finite_difference
end program test_dft_partition_hessian
