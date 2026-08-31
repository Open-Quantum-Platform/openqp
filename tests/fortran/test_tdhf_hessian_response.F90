program test_tdhf_hessian_response

  use precision, only: dp
  use tdhf_hessian_response_mod, only: solve_tdhf_amplitude_response, &
    solve_tdhf_z_response, complete_rhf_orbital_response

  implicit none

  real(kind=dp) :: amb(1,1), apb(1,1), u0(1), v0(1)
  real(kind=dp) :: dambu(1,1), dapbv(1,1), du(1,1), dv(1,1)
  real(kind=dp) :: domega(1), residual
  real(kind=dp) :: ohess(2,2), drhs(2,1), dhz(2,1), dz(2,1)
  real(kind=dp) :: mo(3,3), sx(3,3), umat(3,3), dmo(3,1), dpao(3,3), uov(2)
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

end program test_tdhf_hessian_response
