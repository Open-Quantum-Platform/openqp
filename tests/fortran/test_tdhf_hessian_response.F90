program test_tdhf_hessian_response

  use precision, only: dp
  use tdhf_hessian_response_mod, only: solve_tdhf_amplitude_response, &
    solve_tdhf_z_response, complete_rhf_orbital_response, &
    assemble_tdhf_sigma_derivative, solve_mrsf_tda_response_dense, &
    assemble_mrsf_tda_eigenvalue_hessian

  implicit none

  real(kind=dp) :: amb(1,1), apb(1,1), u0(1), v0(1)
  real(kind=dp) :: dambu(1,1), dapbv(1,1), du(1,1), dv(1,1)
  real(kind=dp) :: domega(1), residual
  real(kind=dp) :: ohess(2,2), drhs(2,1), dhz(2,1), dz(2,1)
  real(kind=dp) :: mo(3,3), sx(3,3), umat(3,3), dmo(3,1), dpao(3,3), uov(2)
  real(kind=dp) :: umat3(3,3,1), epsx(3,1), gmat(3,3)
  real(kind=dp) :: deri(2,1), inner(2,1), dsigma(2,1), ztest(2)
  real(kind=dp) :: mrsf_a(3,3), mrsf_x(3), mrsf_dax(3,2)
  real(kind=dp) :: mrsf_dx(3,2), mrsf_domega(2)
  real(kind=dp) :: mrsf_d2a(2,2),mrsf_hess(2,2),mrsf_asym
  integer :: status

  amb(1,1) = 2.0_dp
  apb(1,1) = 8.0_dp
  u0(1) = 1.0_dp
  v0(1) = 0.5_dp
  dambu(1,1) = 0.2_dp*u0(1)
  dapbv(1,1) = -0.1_dp*v0(1)

  call solve_tdhf_amplitude_response(amb, apb, 4.0_dp, u0, v0, &
                                      dambu, dapbv, du, dv, domega, &
                                      residual, status)

  if (status /= 0) error stop 'projected TDHF response did not converge'
  if (abs(domega(1) - 0.175_dp) > 1.0e-13_dp) &
    error stop 'excitation-energy derivative is incorrect'
  if (abs(du(1,1) + 0.028125_dp) > 1.0e-13_dp) &
    error stop 'd(X+Y) is incorrect'
  if (abs(dv(1,1) - 0.0140625_dp) > 1.0e-13_dp) &
    error stop 'd(X-Y) is incorrect'
  if (abs(v0(1)*du(1,1) + u0(1)*dv(1,1)) > 1.0e-13_dp) &
    error stop 'TDHF response violates differentiated normalization'
  if (residual > 1.0e-13_dp) error stop 'TDHF response residual is too large'

  ohess = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2,2])
  drhs(:,1) = [6.0_dp, 7.0_dp]
  dhz(:,1) = [1.0_dp, 1.0_dp]
  call solve_tdhf_z_response(ohess, drhs, dhz, dz, residual, status)
  if (status /= 0) error stop 'derivative Z-vector response did not converge'
  if (maxval(abs(dz(:,1) - [9.0_dp/11.0_dp, 19.0_dp/11.0_dp])) > &
      1.0e-13_dp) error stop 'derivative Z-vector is incorrect'
  if (residual > 1.0e-13_dp) &
    error stop 'derivative Z-vector residual is too large'

  mo = 0.0_dp
  mo(1,1) = 1.0_dp; mo(2,2) = 1.0_dp; mo(3,3) = 1.0_dp
  sx = reshape([0.2_dp, 0.1_dp, -0.2_dp, &
                0.1_dp, 0.4_dp,  0.3_dp, &
               -0.2_dp, 0.3_dp, -0.6_dp], [3,3])
  uov = [0.7_dp, -0.4_dp]
  call complete_rhf_orbital_response(mo, 1, sx, uov, umat, dmo, dpao, status)
  if (status /= 0) error stop 'orbital-response completion rejected valid data'
  if (maxval(abs(umat + transpose(umat) + sx)) > 1.0e-14_dp) &
    error stop 'orbital response violates differentiated orthonormality'
  if (maxval(abs(dmo(:,1) - [-0.1_dp, 0.7_dp, -0.4_dp])) > 1.0e-14_dp) &
    error stop 'occupied MO derivative is incorrect'
  if (maxval(abs(dpao - reshape([-0.4_dp, 1.4_dp, -0.8_dp, &
                                  1.4_dp, 0.0_dp,  0.0_dp, &
                                 -0.8_dp, 0.0_dp,  0.0_dp], [3,3]))) > 1.0e-14_dp) &
    error stop 'closed-shell AO density derivative is incorrect'

  umat3(:,:,1) = umat
  epsx(:,1) = [0.1_dp, 0.4_dp, -0.2_dp]
  gmat = reshape([1.0_dp, 0.2_dp, 0.3_dp, &
                  0.4_dp, 1.5_dp, 0.6_dp, &
                  0.7_dp, 0.8_dp, 2.0_dp], [3,3])
  deri(:,1) = [0.05_dp, -0.07_dp]
  inner(:,1) = [0.02_dp, 0.09_dp]
  ztest = [0.5_dp, -0.25_dp]
  call assemble_tdhf_sigma_derivative(ztest, 1, umat3, epsx, gmat, &
                                       deri, inner, dsigma, status)
  if (status /= 0) error stop 'sigma-derivative assembly rejected valid data'
  if (maxval(abs(dsigma(:,1) - [0.005_dp, 0.535_dp])) > 1.0e-14_dp) &
    error stop 'sigma-derivative primitive contraction is incorrect'

  mrsf_a = 0.0_dp
  mrsf_a(1,1) = 1.0_dp
  mrsf_a(2,2) = 2.0_dp
  mrsf_a(3,3) = 4.0_dp
  mrsf_x = [0.0_dp,1.0_dp,0.0_dp]
  ! This is (dA)X for a symmetric derivative with dA_12=0.3,
  ! dA_22=0.2, and dA_32=-0.5.
  mrsf_dax(:,1) = [0.3_dp,0.2_dp,-0.5_dp]
  mrsf_dax(:,2) = [0.4_dp,-0.1_dp,0.2_dp]
  call solve_mrsf_tda_response_dense(mrsf_a,2.0_dp,mrsf_x,mrsf_dax, &
    mrsf_dx,mrsf_domega,residual,status)
  if (status /= 0) error stop 'dense MRSF-TDA response rejected valid data'
  if (abs(mrsf_domega(1)-0.2_dp) > 1.0e-14_dp) &
    error stop 'MRSF excitation-energy derivative is incorrect'
  if (maxval(abs(mrsf_dx(:,1)-[0.3_dp,0.0_dp,0.25_dp])) > 1.0e-13_dp) &
    error stop 'MRSF amplitude derivative is incorrect'
  if (maxval(abs(mrsf_dx(:,2)-[0.4_dp,0.0_dp,-0.1_dp])) > 1.0e-13_dp) &
    error stop 'second MRSF amplitude derivative is incorrect'
  if (abs(dot_product(mrsf_x,mrsf_dx(:,1))) > 1.0e-14_dp) &
    error stop 'MRSF amplitude response violates its gauge condition'
  if (residual > 1.0e-13_dp) &
    error stop 'MRSF amplitude response residual is too large'

  mrsf_d2a = reshape([0.7_dp,0.11_dp,0.11_dp,-0.4_dp],[2,2])
  call assemble_mrsf_tda_eigenvalue_hessian(mrsf_x,mrsf_dax,mrsf_dx, &
    mrsf_d2a,mrsf_hess,mrsf_asym,status)
  if (status /= 0) error stop 'MRSF eigenvalue Hessian assembly failed'
  if (maxval(abs(mrsf_hess-reshape([0.63_dp,0.45_dp,0.45_dp,-0.12_dp], &
      [2,2]))) > 1.0e-13_dp) error stop 'MRSF eigenvalue Hessian is incorrect'
  if (mrsf_asym > 1.0e-13_dp) &
    error stop 'MRSF eigenvalue Hessian response is asymmetric'

end program test_tdhf_hessian_response
