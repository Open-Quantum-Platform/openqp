program test_mrsf_xc_fock_total_derivative

  use precision, only: fp
  use mod_dft_gridint_mrsf_xc_fock_deriv_point, only: &
    lda_spin_fock_point_derivative,gga_spin_fock_point_derivative, &
    moving_ao_pair_derivative,moving_density_derivative

  implicit none

  integer, parameter :: nc=6
  real(fp), parameter :: h=1.0e-30_fp
  real(fp) :: qscale,weight,dweight(nc),pair,dpair(nc)
  real(fp) :: vr(2),dvr(2,nc),vs(3),dvs(3,nc),result(2,nc)
  real(fp) :: expected(2,nc),grad_rho(3,2),dgrad_rho(3,2,nc)
  real(fp) :: grad_pair(3),dgrad_pair(3,nc),rho(2),drho(2,nc)
  real(fp) :: aov(2),aog1(2,3),aog2(2,6),density(2,2)
  real(fp) :: density_drho(nc),density_dgrad(3,nc)
  integer :: ao_atom(2),coordinate,spin,status
  complex(fp) :: z2(2),g2(2),value
  complex(fp) :: z5(5),g5(5),grc(3,2),gpc(3),coefficient(3)
  complex(fp) :: weightc,pairc

  qscale=1.7_fp
  weight=0.43_fp
  dweight=[0.03_fp,-0.02_fp,0.01_fp,-0.04_fp,0.015_fp,0.005_fp]
  pair=0.61_fp
  dpair=[0.07_fp,-0.05_fp,0.02_fp,0.03_fp,-0.01_fp,-0.06_fp]
  rho=[0.8_fp,0.35_fp]
  drho(1,:)=[0.04_fp,-0.02_fp,0.03_fp,-0.01_fp,0.05_fp,-0.04_fp]
  drho(2,:)=[-0.03_fp,0.06_fp,-0.02_fp,0.04_fp,-0.01_fp,0.02_fp]

  ! LDA oracle: complex-step differentiation of a local polynomial potential.
  ! This is a point-variable derivative, not a displaced nuclear geometry.
  z2=cmplx(rho,0.0_fp,fp)
  g2=lda_polynomial_gradient(z2)
  vr=real(g2,fp)
  do coordinate=1,nc
    z2=cmplx(rho,h*drho(:,coordinate),fp)
    g2=lda_polynomial_gradient(z2)
    dvr(:,coordinate)=aimag(g2)/h
    weightc=cmplx(weight,h*dweight(coordinate),fp)
    pairc=cmplx(pair,h*dpair(coordinate),fp)
    do spin=1,2
      value=qscale*weightc*g2(spin)*pairc
      expected(spin,coordinate)=aimag(value)/h
    end do
  end do
  call lda_spin_fock_point_derivative(qscale,weight,dweight,vr,dvr,pair, &
    dpair,result,status)
  if(status/=0) error stop 'LDA point derivative rejected valid input'
  if(maxval(abs(result-expected))>2.0e-13_fp) &
    error stop 'LDA point derivative differs from complex-step oracle'

  grad_rho(:,1)=[0.25_fp,-0.31_fp,0.17_fp]
  grad_rho(:,2)=[-0.12_fp,0.22_fp,0.28_fp]
  dgrad_rho=reshape([ &
    0.02_fp,-0.01_fp,0.03_fp, -0.04_fp,0.01_fp,0.02_fp, &
   -0.03_fp,0.05_fp,0.01_fp,  0.02_fp,-0.02_fp,0.04_fp, &
    0.01_fp,0.02_fp,-0.04_fp, -0.01_fp,0.03_fp,-0.02_fp, &
    0.04_fp,-0.02_fp,0.01_fp,  0.03_fp,0.01_fp,-0.03_fp, &
   -0.02_fp,0.01_fp,0.05_fp,  0.01_fp,-0.04_fp,0.02_fp, &
    0.03_fp,0.04_fp,-0.01_fp, -0.02_fp,0.02_fp,0.01_fp],[3,2,nc])
  grad_pair=[0.19_fp,-0.27_fp,0.33_fp]
  dgrad_pair=reshape([ &
    0.02_fp,-0.01_fp,0.04_fp, -0.03_fp,0.05_fp,0.01_fp, &
    0.01_fp,0.03_fp,-0.02_fp,  0.04_fp,-0.01_fp,0.02_fp, &
   -0.02_fp,0.04_fp,0.03_fp,  0.01_fp,0.02_fp,-0.05_fp],[3,nc])

  ! GGA oracle: the same complex-step construction differentiates a
  ! five-variable spin polynomial and the complete AO-gradient integrand.
  z5=cmplx(gga_variables(rho,grad_rho),0.0_fp,fp)
  g5=gga_polynomial_gradient(z5)
  vr=real(g5(1:2),fp)
  vs=real(g5(3:5),fp)
  do coordinate=1,nc
    grc=cmplx(grad_rho,h*dgrad_rho(:,:,coordinate),fp)
    z5(1:2)=cmplx(rho,h*drho(:,coordinate),fp)
    z5(3)=sum(grc(:,1)*grc(:,1))
    z5(4)=sum(grc(:,2)*grc(:,2))
    z5(5)=sum(grc(:,1)*grc(:,2))
    g5=gga_polynomial_gradient(z5)
    dvr(:,coordinate)=aimag(g5(1:2))/h
    dvs(:,coordinate)=aimag(g5(3:5))/h
    weightc=cmplx(weight,h*dweight(coordinate),fp)
    pairc=cmplx(pair,h*dpair(coordinate),fp)
    gpc=cmplx(grad_pair,h*dgrad_pair(:,coordinate),fp)
    coefficient=2.0_fp*g5(3)*grc(:,1)+g5(5)*grc(:,2)
    value=qscale*weightc*(g5(1)*pairc+sum(coefficient*gpc))
    expected(1,coordinate)=aimag(value)/h
    coefficient=2.0_fp*g5(4)*grc(:,2)+g5(5)*grc(:,1)
    value=qscale*weightc*(g5(2)*pairc+sum(coefficient*gpc))
    expected(2,coordinate)=aimag(value)/h
  end do
  call gga_spin_fock_point_derivative(qscale,weight,dweight,vr,vs,dvr,dvs, &
    grad_rho,dgrad_rho,pair,grad_pair,dpair,dgrad_pair,result,status)
  if(status/=0) error stop 'GGA point derivative rejected valid input'
  if(maxval(abs(result-expected))>4.0e-13_fp) &
    error stop 'GGA point derivative differs from complex-step oracle'

  ! Moving-AO and density derivatives cancel exactly under a rigid
  ! translation of all atoms for every Cartesian component.
  aov=[0.73_fp,-0.41_fp]
  aog1=reshape([0.12_fp,-0.23_fp,0.31_fp, &
                -0.17_fp,0.29_fp,0.08_fp],[2,3])
  aog2=reshape([0.11_fp,-0.07_fp,0.04_fp,0.03_fp,-0.02_fp,0.06_fp, &
                -0.09_fp,0.05_fp,0.13_fp,-0.04_fp,0.07_fp,-0.01_fp],[2,6])
  ao_atom=[1,2]
  call moving_ao_pair_derivative(1,ao_atom(1),ao_atom(2),aov(1),aov(2), &
    aog1(1,:),aog1(2,:),aog2(1,:),aog2(2,:),pair,grad_pair,dpair, &
    dgrad_pair,status)
  if(status/=0) error stop 'moving AO-pair derivative failed'
  do coordinate=1,3
    if(abs(dpair(coordinate)+dpair(coordinate+3))>1.0e-14_fp .or. &
       maxval(abs(dgrad_pair(:,coordinate)+ &
       dgrad_pair(:,coordinate+3)))>1.0e-14_fp) &
      error stop 'moving AO-pair derivative violates translation invariance'
  end do
  density=reshape([0.8_fp,-0.2_fp,-0.2_fp,0.45_fp],[2,2])
  call moving_density_derivative(density,ao_atom,1,aov,aog1,aog2,.true., &
    density_drho,density_dgrad,status)
  if(status/=0) error stop 'moving density derivative failed'
  do coordinate=1,3
    if(abs(density_drho(coordinate)+density_drho(coordinate+3))>1.0e-14_fp &
       .or. maxval(abs(density_dgrad(:,coordinate)+ &
       density_dgrad(:,coordinate+3)))>1.0e-14_fp) &
      error stop 'moving density derivative violates translation invariance'
  end do

  ! Translation of the complete GGA point contribution, including a
  ! normalized partition derivative and coefficient response, also vanishes.
  dweight(1:3)=[0.02_fp,-0.03_fp,0.01_fp]
  dweight(4:6)=-dweight(1:3)
  drho(1,:)=density_drho
  drho(2,:)=0.6_fp*density_drho
  dgrad_rho(:,1,:)=density_dgrad
  dgrad_rho(:,2,:)=0.6_fp*density_dgrad
  do coordinate=1,nc
    dvr(1,coordinate)=0.3_fp*drho(1,coordinate) &
      +0.1_fp*drho(2,coordinate)
    dvr(2,coordinate)=0.1_fp*drho(1,coordinate) &
      +0.4_fp*drho(2,coordinate)
    dvs(1,coordinate)=0.05_fp*drho(1,coordinate)
    dvs(2,coordinate)=-0.03_fp*drho(2,coordinate)
    dvs(3,coordinate)=0.02_fp*sum(dgrad_rho(:,:,coordinate))
  end do
  call gga_spin_fock_point_derivative(qscale,weight,dweight,vr,vs,dvr,dvs, &
    grad_rho,dgrad_rho,pair,grad_pair,dpair,dgrad_pair,result,status)
  if(status/=0) error stop 'translation test GGA derivative failed'
  do coordinate=1,3
    if(maxval(abs(result(:,coordinate)+result(:,coordinate+3)))>1.0e-13_fp) &
      error stop 'complete GGA point derivative violates translation invariance'
  end do

contains

  pure function lda_polynomial_gradient(z) result(g)
    complex(fp), intent(in) :: z(2)
    complex(fp) :: g(2)
    g(1)=0.7_fp*z(1)**2+0.13_fp*z(2)+0.09_fp*z(1)*z(2)
    g(2)=0.45_fp*z(2)**2+0.13_fp*z(1)+0.045_fp*z(1)**2
  end function lda_polynomial_gradient

  pure function gga_variables(local_rho,local_gradient) result(z)
    real(fp), intent(in) :: local_rho(2),local_gradient(3,2)
    real(fp) :: z(5)
    z(1:2)=local_rho
    z(3)=dot_product(local_gradient(:,1),local_gradient(:,1))
    z(4)=dot_product(local_gradient(:,2),local_gradient(:,2))
    z(5)=dot_product(local_gradient(:,1),local_gradient(:,2))
  end function gga_variables

  pure function gga_polynomial_gradient(z) result(g)
    complex(fp), intent(in) :: z(5)
    complex(fp) :: g(5)
    g(1)=0.31_fp*z(1)**2+0.07_fp*z(2)+0.05_fp*z(3)+0.03_fp*z(5)
    g(2)=0.27_fp*z(2)**2+0.07_fp*z(1)+0.04_fp*z(4)-0.02_fp*z(5)
    g(3)=0.05_fp*z(1)+0.11_fp*z(3)+0.025_fp*z(5)
    g(4)=0.04_fp*z(2)+0.09_fp*z(4)-0.015_fp*z(5)
    g(5)=0.03_fp*z(1)-0.02_fp*z(2)+0.025_fp*z(3) &
      -0.015_fp*z(4)+0.08_fp*z(5)
  end function gga_polynomial_gradient

end program test_mrsf_xc_fock_total_derivative
