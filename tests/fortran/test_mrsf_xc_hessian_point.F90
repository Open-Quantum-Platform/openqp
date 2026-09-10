program test_mrsf_xc_hessian_point

  use precision, only: fp
  use mod_dft_gridint_mrsf_xc_hessian_point, only: &
    mrsf_xc_weighted_fixed_hessian,mrsf_xc_weighted_response_row, &
    mrsf_xc_kernel_fock_coefficients

  implicit none

  integer, parameter :: nv=5,nat=2,nc=3*nat
  real(fp), parameter :: step=1.0e-30_fp,scale=1.37_fp
  real(fp) :: a(nv),b(nv,nv),t(nv,nv,nv),g(nv),p(nv)
  real(fp) :: dg(nv,nc),dp(nv,nc),d2g(nv,nc,nc),d2p(nv,nc,nc)
  real(fp) :: response_g(nv),response_p(nv),dresponse_g(nv,nc)
  real(fp) :: dresponse_p(nv,nc),first(nv),second(nv,nv),third(nv,nv,nv)
  real(fp) :: weight,dweight(nc),d2weight(nc,nc),hessian(nc,nc),row(nc)
  real(fp) :: hessian_block(3,3,nat,nat),row_block(3,nat)
  real(fp) :: d2g_block(nv,3,3,nat,nat),d2p_block(nv,3,3,nat,nat)
  real(fp) :: d2weight_block(3,nat,3,nat)
  real(fp) :: expected_hessian(nc,nc),expected_row(nc),exc
  real(fp) :: grad_ground(3,2),grad_probe(3,2)
  real(fp) :: dgrad_ground(3,2,nc),dgrad_probe(3,2,nc)
  real(fp) :: v_r(2),coefficient(3,2),dv_r(2,nc),dcoefficient(3,2,nc)
  real(fp) :: expected_v_r(2),expected_coefficient(3,2)
  real(fp) :: expected_dv_r(2,nc),expected_dcoefficient(3,2,nc)
  complex(fp) :: coordinates(nc),direction
  complex(fp) :: complex_v_r(2),complex_coefficient(3,2)
  integer :: i,j,k,status,atom_a,atom_b,cart_a,cart_b

  call initialize_polynomial(a,b,t)
  g=[0.61_fp,0.27_fp,0.13_fp,0.08_fp,-0.04_fp]
  p=[0.07_fp,-0.03_fp,0.02_fp,0.05_fp,-0.01_fp]
  weight=0.43_fp
  do j=1,nc
    dweight(j)=0.006_fp*real(2*j-5,fp)
    do i=1,nv
      dg(i,j)=0.004_fp*real(3*i-2*j,fp)
      dp(i,j)=0.003_fp*real(i+j-4,fp)
      dresponse_g(i,j)=0.002_fp*real(2*i+j-5,fp)
      dresponse_p(i,j)=-0.0015_fp*real(i-2*j+1,fp)
    end do
  end do
  do k=1,nc
    do j=1,nc
      d2weight(j,k)=0.0007_fp*real(j+k-4,fp)
      do i=1,nv
        d2g(i,j,k)=0.0005_fp*real(i+j+k-5,fp)
        d2p(i,j,k)=-0.0004_fp*real(2*i-j-k,fp)
      end do
    end do
  end do
  response_g=[0.015_fp,-0.012_fp,0.009_fp,0.006_fp,-0.004_fp]
  response_p=[-0.008_fp,0.011_fp,-0.005_fp,0.007_fp,0.003_fp]
  do atom_b=1,nat
    do cart_b=1,3
      k=3*(atom_b-1)+cart_b
      do atom_a=1,nat
        do cart_a=1,3
          j=3*(atom_a-1)+cart_a
          d2g_block(:,cart_a,cart_b,atom_a,atom_b)=d2g(:,j,k)
          d2p_block(:,cart_a,cart_b,atom_a,atom_b)=d2p(:,j,k)
          d2weight_block(cart_a,atom_a,cart_b,atom_b)=d2weight(j,k)
        end do
      end do
    end do
  end do

  call polynomial_derivatives(cmplx(g,0.0_fp,fp),a,b,t,exc,first,second)
  third=t
  hessian_block=0.0_fp
  call mrsf_xc_weighted_fixed_hessian(exc,first,second,third,p, &
    reshape(dg,[nv,3,nat]),reshape(dp,[nv,3,nat]), &
    d2g_block,d2p_block,scale,weight,reshape(dweight,[3,nat]), &
    d2weight_block, &
    hessian_block,status)
  if(status/=0) error stop 'fixed XC Hessian rejected valid input'
  do atom_b=1,nat
    do cart_b=1,3
      k=3*(atom_b-1)+cart_b
      do atom_a=1,nat
        do cart_a=1,3
          j=3*(atom_a-1)+cart_a
          hessian(j,k)=hessian_block(cart_a,cart_b,atom_a,atom_b)
        end do
      end do
    end do
  end do

  do k=1,nc
    coordinates=cmplx(0.0_fp,0.0_fp,fp)
    coordinates(k)=cmplx(0.0_fp,step,fp)
    expected_hessian(:,k)=aimag(nuclear_gradient(coordinates,a,b,t,g,p, &
      dg,dp,d2g,d2p,weight,dweight,d2weight))/step
  end do
  if(maxval(abs(hessian-expected_hessian))>2.0e-12_fp) &
    error stop 'fixed XC Hessian differs from complex-step oracle'

  row_block=0.0_fp
  call mrsf_xc_weighted_response_row(first,second,third,p, &
    reshape(dg,[nv,3,nat]),reshape(dp,[nv,3,nat]),response_g,response_p, &
    reshape(dresponse_g,[nv,3,nat]), &
    reshape(dresponse_p,[nv,3,nat]),scale,weight, &
    reshape(dweight,[3,nat]),row_block,status)
  if(status/=0) error stop 'XC response row rejected valid input'
  row=reshape(row_block,[nc])
  direction=cmplx(0.0_fp,step,fp)
  expected_row=aimag(response_gradient(direction,a,b,t,g,p,dg,dp, &
    response_g,response_p,dresponse_g,dresponse_p,weight,dweight))/step
  if(maxval(abs(row-expected_row))>2.0e-12_fp) &
    error stop 'XC response row differs from complex-step oracle'

  grad_ground=reshape([0.21_fp,-0.17_fp,0.09_fp, &
    -0.13_fp,0.19_fp,0.07_fp],[3,2])
  grad_probe=reshape([0.05_fp,-0.03_fp,0.02_fp, &
    0.04_fp,0.01_fp,-0.06_fp],[3,2])
  do k=1,nc
    do j=1,2
      do i=1,3
        dgrad_ground(i,j,k)=0.002_fp*real(i+2*j-k,fp)
        dgrad_probe(i,j,k)=-0.001_fp*real(2*i-j+k,fp)
      end do
    end do
  end do
  call mrsf_xc_kernel_fock_coefficients(first,second,third,p,dg,dp,.true., &
    grad_ground,grad_probe,dgrad_ground,dgrad_probe,v_r,coefficient,dv_r, &
    dcoefficient,status)
  if(status/=0) error stop 'XC-kernel coefficients rejected valid input'
  call complex_kernel_coefficients(cmplx(g,0.0_fp,fp), &
    cmplx(p,0.0_fp,fp),cmplx(grad_ground,0.0_fp,fp), &
    cmplx(grad_probe,0.0_fp,fp),a,b,t,complex_v_r,complex_coefficient)
  expected_v_r=real(complex_v_r,fp)
  expected_coefficient=real(complex_coefficient,fp)
  do k=1,nc
    direction=cmplx(0.0_fp,step,fp)
    call complex_kernel_coefficients(cmplx(g,0.0_fp,fp)+direction*dg(:,k), &
      cmplx(p,0.0_fp,fp)+direction*dp(:,k), &
      cmplx(grad_ground,0.0_fp,fp)+direction*dgrad_ground(:,:,k), &
      cmplx(grad_probe,0.0_fp,fp)+direction*dgrad_probe(:,:,k), &
      a,b,t,complex_v_r,complex_coefficient)
    expected_dv_r(:,k)=aimag(complex_v_r)/step
    expected_dcoefficient(:,:,k)=aimag(complex_coefficient)/step
  end do
  if(maxval(abs(v_r-expected_v_r))>2.0e-13_fp .or. &
     maxval(abs(coefficient-expected_coefficient))>2.0e-13_fp .or. &
     maxval(abs(dv_r-expected_dv_r))>2.0e-12_fp .or. &
     maxval(abs(dcoefficient-expected_dcoefficient))>2.0e-12_fp) &
    error stop 'XC-kernel AO coefficients differ from complex-step oracle'

contains

  subroutine initialize_polynomial(linear,quadratic,cubic)
    real(fp), intent(out) :: linear(nv),quadratic(nv,nv),cubic(nv,nv,nv)
    integer :: ii,jj,kk
    linear=[0.21_fp,-0.17_fp,0.09_fp,0.04_fp,-0.03_fp]
    do jj=1,nv
      do ii=1,nv
        quadratic(ii,jj)=0.011_fp*real(ii+jj,fp)
        do kk=1,nv
          cubic(ii,jj,kk)=0.0013_fp*real(ii+jj+kk,fp)
        end do
      end do
    end do
  end subroutine initialize_polynomial

  subroutine polynomial_derivatives(z,linear,quadratic,cubic,value,grad,hess)
    complex(fp), intent(in) :: z(nv)
    real(fp), intent(in) :: linear(nv),quadratic(nv,nv),cubic(nv,nv,nv)
    real(fp), intent(out) :: value,grad(nv),hess(nv,nv)
    complex(fp) :: complex_value,complex_grad(nv),complex_hess(nv,nv)
    call complex_polynomial_derivatives(z,linear,quadratic,cubic, &
      complex_value,complex_grad,complex_hess)
    value=real(complex_value,fp)
    grad=real(complex_grad,fp)
    hess=real(complex_hess,fp)
  end subroutine polynomial_derivatives

  subroutine complex_polynomial_derivatives(z,linear,quadratic,cubic, &
      value,grad,hess)
    complex(fp), intent(in) :: z(nv)
    real(fp), intent(in) :: linear(nv),quadratic(nv,nv),cubic(nv,nv,nv)
    complex(fp), intent(out) :: value,grad(nv),hess(nv,nv)
    integer :: ii,jj,kk
    value=cmplx(0.19_fp,0.0_fp,fp)
    grad=cmplx(linear,0.0_fp,fp)
    hess=cmplx(quadratic,0.0_fp,fp)
    do ii=1,nv
      value=value+linear(ii)*z(ii)
      do jj=1,nv
        value=value+0.5_fp*quadratic(ii,jj)*z(ii)*z(jj)
        grad(ii)=grad(ii)+quadratic(ii,jj)*z(jj)
        do kk=1,nv
          value=value+cubic(ii,jj,kk)*z(ii)*z(jj)*z(kk)/6.0_fp
          grad(ii)=grad(ii)+0.5_fp*cubic(ii,jj,kk)*z(jj)*z(kk)
          hess(ii,jj)=hess(ii,jj)+cubic(ii,jj,kk)*z(kk)
        end do
      end do
    end do
  end subroutine complex_polynomial_derivatives

  function nuclear_gradient(x,linear,quadratic,cubic,g0,p0,g1,p1,g2,p2, &
      w0,w1,w2) result(result)
    complex(fp), intent(in) :: x(nc)
    real(fp), intent(in) :: linear(nv),quadratic(nv,nv),cubic(nv,nv,nv)
    real(fp), intent(in) :: g0(nv),p0(nv),g1(nv,nc),p1(nv,nc)
    real(fp), intent(in) :: g2(nv,nc,nc),p2(nv,nc,nc)
    real(fp), intent(in) :: w0,w1(nc),w2(nc,nc)
    complex(fp) :: result(nc),gz(nv),pz(nv),gx(nv,nc),px(nv,nc),w,wx(nc)
    complex(fp) :: value,grad(nv),hess(nv,nv),q,qx
    integer :: ii,jj
    gz=cmplx(g0,0.0_fp,fp);pz=cmplx(p0,0.0_fp,fp)
    gx=cmplx(g1,0.0_fp,fp);px=cmplx(p1,0.0_fp,fp)
    w=cmplx(w0,0.0_fp,fp);wx=cmplx(w1,0.0_fp,fp)
    do ii=1,nc
      gz=gz+g1(:,ii)*x(ii)
      pz=pz+p1(:,ii)*x(ii)
      w=w+w1(ii)*x(ii)
      do jj=1,nc
        gz=gz+0.5_fp*g2(:,ii,jj)*x(ii)*x(jj)
        pz=pz+0.5_fp*p2(:,ii,jj)*x(ii)*x(jj)
        w=w+0.5_fp*w2(ii,jj)*x(ii)*x(jj)
        gx(:,ii)=gx(:,ii)+g2(:,ii,jj)*x(jj)
        px(:,ii)=px(:,ii)+p2(:,ii,jj)*x(jj)
        wx(ii)=wx(ii)+w2(ii,jj)*x(jj)
      end do
    end do
    call complex_polynomial_derivatives(gz,linear,quadratic,cubic, &
      value,grad,hess)
    q=value+sum(grad*pz)
    do ii=1,nc
      qx=sum((grad+matmul(hess,pz))*gx(:,ii))+sum(grad*px(:,ii))
      result(ii)=scale*(wx(ii)*q+w*qx)
    end do
  end function nuclear_gradient

  function response_gradient(s,linear,quadratic,cubic,g0,p0,g1,p1, &
      rg0,rp0,rg1,rp1,w0,w1) result(result)
    complex(fp), intent(in) :: s
    real(fp), intent(in) :: linear(nv),quadratic(nv,nv),cubic(nv,nv,nv)
    real(fp), intent(in) :: g0(nv),p0(nv),g1(nv,nc),p1(nv,nc)
    real(fp), intent(in) :: rg0(nv),rp0(nv),rg1(nv,nc),rp1(nv,nc)
    real(fp), intent(in) :: w0,w1(nc)
    complex(fp) :: result(nc),gz(nv),pz(nv),gx(nv,nc),px(nv,nc)
    complex(fp) :: value,grad(nv),hess(nv,nv),q,qx
    integer :: ii
    gz=cmplx(g0,0.0_fp,fp)+s*rg0
    pz=cmplx(p0,0.0_fp,fp)+s*rp0
    gx=cmplx(g1,0.0_fp,fp)+s*rg1
    px=cmplx(p1,0.0_fp,fp)+s*rp1
    call complex_polynomial_derivatives(gz,linear,quadratic,cubic, &
      value,grad,hess)
    q=value+sum(grad*pz)
    do ii=1,nc
      qx=sum((grad+matmul(hess,pz))*gx(:,ii))+sum(grad*px(:,ii))
      result(ii)=scale*(w1(ii)*q+w0*qx)
    end do
  end function response_gradient

  subroutine complex_kernel_coefficients(z,pz,grz,pgrz,linear,quadratic, &
      cubic,vr,coefficient_out)
    complex(fp), intent(in) :: z(nv),pz(nv),grz(3,2),pgrz(3,2)
    real(fp), intent(in) :: linear(nv),quadratic(nv,nv),cubic(nv,nv,nv)
    complex(fp), intent(out) :: vr(2),coefficient_out(3,2)
    complex(fp) :: value,grad(nv),hess(nv,nv),kernel(nv)
    call complex_polynomial_derivatives(z,linear,quadratic,cubic, &
      value,grad,hess)
    kernel=matmul(hess,pz)
    vr=kernel(1:2)
    coefficient_out(:,1)=2.0_fp*kernel(3)*grz(:,1)+kernel(5)*grz(:,2)+ &
      2.0_fp*grad(3)*pgrz(:,1)+grad(5)*pgrz(:,2)
    coefficient_out(:,2)=2.0_fp*kernel(4)*grz(:,2)+kernel(5)*grz(:,1)+ &
      2.0_fp*grad(4)*pgrz(:,2)+grad(5)*pgrz(:,1)
  end subroutine complex_kernel_coefficients

end program test_mrsf_xc_hessian_point
