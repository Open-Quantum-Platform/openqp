program test_tdhf_mrsf_hessian_components

  use precision, only: dp
  use tdhf_mrsf_hessian_components_mod, only: mrsf_hessian_is_applicable, &
    assemble_mrsf_cartesian_hessian,project_mrsf_rigid_translations

  implicit none

  real(kind=dp) :: fixed(6,6),xc(6,6),rows(6,6),hessian(6,6),asymmetry
  integer :: atom,cart,index,status

  if(.not.mrsf_hessian_is_applicable(3,3,4,2,1,.false.,.true., &
      .false.,.false.,.false.,1)) error stop 'valid MRSF scope rejected'
  if(mrsf_hessian_is_applicable(3,3,4,2,1,.true.,.true., &
      .false.,.false.,.false.,1)) error stop 'UMRSF scope was not rejected'
  if(mrsf_hessian_is_applicable(3,3,4,2,1,.false.,.true., &
      .false.,.false.,.true.,1)) error stop 'CAM scope was not rejected'

  fixed=0.0_dp
  xc=0.0_dp
  rows=0.0_dp
  do index=1,6
    fixed(index,index)=1.0_dp+0.1_dp*real(index,dp)
  end do
  rows(1,2)=0.7_dp
  rows(2,1)=-0.1_dp
  call assemble_mrsf_cartesian_hessian(fixed,xc,rows,hessian,asymmetry,status)
  if(status/=0) error stop 'MRSF Hessian assembly failed'
  if(abs(asymmetry-0.8_dp)>1.0e-14_dp) &
    error stop 'MRSF row asymmetry diagnostic is wrong'
  if(abs(hessian(1,2)-0.3_dp)>1.0e-14_dp .or. &
     abs(hessian(2,1)-0.3_dp)>1.0e-14_dp) &
    error stop 'MRSF response rows were not symmetrized once'

  call project_mrsf_rigid_translations(hessian,2,status)
  if(status/=0) error stop 'translation projection failed'
  do cart=1,3
    do index=1,6
      if(abs(sum([(hessian(3*(atom-1)+cart,index),atom=1,2)]))>1.0e-13_dp) &
        error stop 'projected MRSF Hessian has a row translation component'
      if(abs(sum([(hessian(index,3*(atom-1)+cart),atom=1,2)]))>1.0e-13_dp) &
        error stop 'projected MRSF Hessian has a column translation component'
    end do
  end do

end program test_tdhf_mrsf_hessian_components
