program test_dft_gga_nuclear_point

  use precision, only: fp
  use mod_dft_gga_nuclear_point, only: gga_density_nuclear_point

  implicit none

  integer, parameter :: nao = 3, nat = 2
  real(fp), parameter :: h1 = 1.0e-5_fp, h2 = 1.0e-5_fp
  integer :: ao_atom(nao), a, b, atom_a, atom_b
  real(fp) :: density(nao,nao), centers(3,nat), point(3)
  real(fp) :: aov(nao), g1(nao,3), g2(nao,6), g3(nao,10)
  real(fp) :: drho(3,nat), dgrho(3,3,nat)
  real(fp) :: d2rho(3,3,nat,nat), d2grho(3,3,3,nat,nat)
  real(fp) :: cp(3,nat), cm(3,nat), rp, rm, gp(3), gm(3)
  real(fp) :: drp(3,nat), drm(3,nat), dgp(3,3,nat), dgm(3,3,nat)
  real(fp) :: err1, err2

  ao_atom = [1, 1, 2]
  centers(:,1) = [-0.35_fp, 0.11_fp, -0.08_fp]
  centers(:,2) = [ 0.42_fp,-0.16_fp,  0.21_fp]
  point = [0.13_fp, -0.24_fp, 0.37_fp]
  density = reshape([1.1_fp,0.2_fp,-0.1_fp, 0.2_fp,0.7_fp,0.3_fp, &
                    -0.1_fp,0.3_fp,0.9_fp], [nao,nao])

  call ao_data(centers, aov, g1, g2, g3)
  call gga_density_nuclear_point(density, ao_atom, aov, g1, g2, g3, &
                                 drho, dgrho, d2rho, d2grho)

  err1 = 0.0_fp
  do atom_a = 1, nat
    do a = 1, 3
      cp = centers; cm = centers
      cp(a,atom_a) = cp(a,atom_a) + h1
      cm(a,atom_a) = cm(a,atom_a) - h1
      call rho_data(cp, rp, gp)
      call rho_data(cm, rm, gm)
      err1 = max(err1, abs(drho(a,atom_a)-(rp-rm)/(2*h1)))
      err1 = max(err1, maxval(abs(dgrho(:,a,atom_a)-(gp-gm)/(2*h1))))
    end do
  end do

  err2 = 0.0_fp
  do atom_b = 1, nat
    do atom_a = 1, nat
      do b = 1, 3
        do a = 1, 3
          cp=centers; cm=centers
          cp(b,atom_b)=cp(b,atom_b)+h2
          cm(b,atom_b)=cm(b,atom_b)-h2
          call rho_derivatives(cp,drp,dgp)
          call rho_derivatives(cm,drm,dgm)
          err2=max(err2,abs(d2rho(a,b,atom_a,atom_b) &
            -(drp(a,atom_a)-drm(a,atom_a))/(2*h2)))
          err2=max(err2,maxval(abs(d2grho(:,a,b,atom_a,atom_b) &
            -(dgp(:,a,atom_a)-dgm(:,a,atom_a))/(2*h2))))
        end do
      end do
    end do
  end do

  print '(a,es12.4)', 'GGA point first-derivative max error  = ', err1
  print '(a,es12.4)', 'GGA point second-derivative max error = ', err2
  if (err1 > 2.0e-9_fp .or. err2 > 2.0e-6_fp) &
    error stop 'GGA nuclear point derivatives failed finite-difference check'

contains

  subroutine rho_data(c, rho, grad)
    real(fp), intent(in) :: c(3,nat)
    real(fp), intent(out) :: rho, grad(3)
    real(fp) :: v(nao), dg(nao,3), dh(nao,6), dt(nao,10)
    call ao_data(c,v,dg,dh,dt)
    rho = dot_product(v,matmul(density,v))
    grad = 2.0_fp*matmul(transpose(dg),matmul(density,v))
  end subroutine rho_data

  subroutine rho_derivatives(c, d_rho, d_grad)
    real(fp), intent(in) :: c(3,nat)
    real(fp), intent(out) :: d_rho(3,nat), d_grad(3,3,nat)
    real(fp) :: v(nao), dg(nao,3), dh(nao,6), dt(nao,10)
    real(fp) :: d2r(3,3,nat,nat), d2g(3,3,3,nat,nat)
    call ao_data(c,v,dg,dh,dt)
    call gga_density_nuclear_point(density,ao_atom,v,dg,dh,dt, &
                                   d_rho,d_grad,d2r,d2g)
  end subroutine rho_derivatives

  subroutine ao_data(c, v, dg, dh, dt)
    real(fp), intent(in) :: c(3,nat)
    real(fp), intent(out) :: v(nao), dg(nao,3), dh(nao,6), dt(nao,10)
    integer :: i
    real(fp) :: x, y, z, alpha, e
    do i=1,nao
      alpha = 0.45_fp + 0.17_fp*i
      x=point(1)-c(1,ao_atom(i)); y=point(2)-c(2,ao_atom(i))
      z=point(3)-c(3,ao_atom(i)); e=exp(-alpha*(x*x+y*y+z*z))
      v(i)=e
      dg(i,:)=[-2*alpha*x,-2*alpha*y,-2*alpha*z]*e
      dh(i,:)=[4*alpha**2*x*x-2*alpha,4*alpha**2*y*y-2*alpha, &
        4*alpha**2*z*z-2*alpha,4*alpha**2*x*y,4*alpha**2*y*z, &
        4*alpha**2*x*z]*e
      dt(i,:)=[12*alpha**2*x-8*alpha**3*x**3, &
        12*alpha**2*y-8*alpha**3*y**3,12*alpha**2*z-8*alpha**3*z**3, &
        (4*alpha**2*y-8*alpha**3*x*x*y), &
        (4*alpha**2*z-8*alpha**3*x*x*z), &
        (4*alpha**2*x-8*alpha**3*y*y*x), &
        (4*alpha**2*z-8*alpha**3*y*y*z), &
        (4*alpha**2*x-8*alpha**3*z*z*x), &
        (4*alpha**2*y-8*alpha**3*z*z*y),-8*alpha**3*x*y*z]*e
    end do
  end subroutine ao_data

end program test_dft_gga_nuclear_point
