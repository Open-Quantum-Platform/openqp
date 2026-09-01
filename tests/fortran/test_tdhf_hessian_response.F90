program test_tdhf_hessian_response

  use precision, only: dp
  use tdhf_hessian_response_mod, only: solve_tdhf_amplitude_response, &
    solve_tdhf_z_response, complete_rhf_orbital_response, &
    assemble_tdhf_sigma_derivative, solve_mrsf_tda_response_dense, &
    solve_mrsf_tda_response_matrix_free, &
    solve_mrsf_tda_cluster_response_matrix_free, &
    solve_mrsf_z_response_matrix_free, &
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
  real(kind=dp) :: cluster_a(5,5),cluster_e(2),cluster_x(5,2)
  real(kind=dp) :: cluster_dax(5,2),cluster_z(5,2),cluster_w(2,2)
  real(kind=dp) :: expected_z(5,2),base_z(5,2),base_w(2,2)
  real(kind=dp) :: rotated_x(5,2),rotated_dax(5,2),rotated_z(5,2)
  real(kind=dp) :: rotated_w(2,2),rotation(2,2),projector_deriv(5,5)
  real(kind=dp) :: rotated_projector_deriv(5,5),near_limit_z(5,2)
  real(kind=dp) :: isolated_x(3,1),isolated_dax(3,1),isolated_z(3,1)
  real(kind=dp) :: isolated_w(1,1),angle,c,s
  real(kind=dp) :: zorb(3,3),zrhs(3,2),zdhz(3,2),zdz(3,2),zexpected(3,2)
  integer :: i,j,status

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

  ! MRSF uses the same differentiated Lagrange-multiplier equation, but the
  ! ROHF orbital Hessian can be indefinite.  The matrix-free GMRES path must
  ! therefore converge without a positive-definite assumption.
  zorb=reshape([0.0_dp,2.0_dp,0.0_dp, &
                2.0_dp,0.0_dp,1.0_dp, &
                0.0_dp,1.0_dp,-1.0_dp],[3,3])
  zexpected(:,1)=[0.4_dp,-0.2_dp,0.7_dp]
  zexpected(:,2)=[-0.3_dp,0.8_dp,-0.5_dp]
  zdhz(:,1)=[0.1_dp,-0.4_dp,0.2_dp]
  zdhz(:,2)=[-0.2_dp,0.3_dp,0.6_dp]
  zrhs=matmul(zorb,zexpected)+zdhz
  call solve_mrsf_z_response_matrix_free(apply_test_z_orbital_hessian,zrhs, &
    zdhz,zdz,residual,status,tol=1.0e-13_dp,maxit=40,restart=3)
  if(status/=0) error stop 'matrix-free MRSF derivative Z solve failed'
  if(maxval(abs(zdz-zexpected))>1.0e-11_dp) &
    error stop 'matrix-free MRSF derivative Z is incorrect'
  if(residual>1.0e-11_dp) &
    error stop 'matrix-free MRSF derivative Z residual is too large'

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

  mrsf_dx=0.0_dp; mrsf_domega=0.0_dp; residual=0.0_dp
  call solve_mrsf_tda_response_matrix_free(apply_test_mrsf_operator,2.0_dp, &
    mrsf_x,mrsf_dax,mrsf_dx,mrsf_domega,residual,status,tol=1.0e-13_dp, &
    maxit=30,restart=3)
  if (status /= 0) error stop 'matrix-free MRSF-TDA response did not converge'
  if (maxval(abs(mrsf_dx(:,1)-[0.3_dp,0.0_dp,0.25_dp])) > 1.0e-12_dp) &
    error stop 'matrix-free MRSF response is incorrect'
  if (maxval(abs(mrsf_dx(:,2)-[0.4_dp,0.0_dp,-0.1_dp])) > 1.0e-12_dp) &
    error stop 'matrix-free second MRSF response is incorrect'
  if (residual > 1.0e-11_dp) &
    error stop 'matrix-free MRSF response residual is too large'

  ! Force the iterative path to stop before convergence.  The bounded dense
  ! projected fallback must recover the same spin-adapted response without
  ! introducing a determinant representation.
  mrsf_dx=0.0_dp; mrsf_domega=0.0_dp; residual=0.0_dp
  call solve_mrsf_tda_response_matrix_free(apply_test_mrsf_operator,2.0_dp, &
    mrsf_x,mrsf_dax,mrsf_dx,mrsf_domega,residual,status,tol=1.0e-13_dp, &
    maxit=1,restart=1)
  if(status/=0) error stop 'dense projected MRSF fallback failed'
  if(maxval(abs(mrsf_dx(:,1)-[0.3_dp,0.0_dp,0.25_dp]))>1.0e-12_dp .or. &
     maxval(abs(mrsf_dx(:,2)-[0.4_dp,0.0_dp,-0.1_dp]))>1.0e-12_dp) &
    error stop 'dense projected MRSF fallback is incorrect'
  if(residual>1.0e-11_dp) &
    error stop 'dense projected MRSF fallback residual is too large'

  mrsf_d2a = reshape([0.7_dp,0.11_dp,0.11_dp,-0.4_dp],[2,2])
  call assemble_mrsf_tda_eigenvalue_hessian(mrsf_x,mrsf_dax,mrsf_dx, &
    mrsf_d2a,mrsf_hess,mrsf_asym,status)
  if (status /= 0) error stop 'MRSF eigenvalue Hessian assembly failed'
  if (maxval(abs(mrsf_hess-reshape([0.63_dp,0.45_dp,0.45_dp,-0.12_dp], &
      [2,2]))) > 1.0e-13_dp) error stop 'MRSF eigenvalue Hessian is incorrect'
  if (mrsf_asym > 1.0e-13_dp) &
    error stop 'MRSF eigenvalue Hessian response is asymmetric'

  ! Exact two-root cluster in a five-dimensional symmetric physical space.
  ! The first two basis vectors span the root cluster; no determinant or
  ! Fock-space representation enters this test.
  cluster_a=0.0_dp
  cluster_a(1,1)=1.0_dp
  cluster_a(2,2)=1.0001_dp
  cluster_a(3,3)=3.0_dp
  cluster_a(4,4)=4.0_dp
  cluster_a(5,5)=6.0_dp
  cluster_e=[1.0_dp,1.0001_dp]
  cluster_x=0.0_dp
  cluster_x(1,1)=1.0_dp
  cluster_x(2,2)=1.0_dp
  cluster_dax(:,1)=[0.2_dp,-0.3_dp,0.4_dp,-0.6_dp,0.5_dp]
  cluster_dax(:,2)=[-0.3_dp,-0.1_dp,0.7_dp,0.2_dp,-0.8_dp]
  call solve_mrsf_tda_cluster_response_matrix_free( &
    apply_test_cluster_operator,cluster_e,cluster_x,cluster_dax,cluster_z, &
    cluster_w,residual,status,tol=1.0e-13_dp,maxit=100,restart=10)
  if(status/=0) error stop 'near-degenerate MRSF cluster response failed'
  expected_z=0.0_dp
  do j=1,2
    do i=3,5
      expected_z(i,j)=-cluster_dax(i,j)/(cluster_a(i,i)-cluster_e(j))
    end do
  end do
  if(maxval(abs(cluster_z-expected_z))>1.0e-11_dp) &
    error stop 'MRSF cluster complement response is incorrect'
  if(maxval(abs(cluster_w-cluster_dax(1:2,:)))>1.0e-13_dp) &
    error stop 'MRSF cluster effective derivative is incorrect'
  if(maxval(abs(matmul(transpose(cluster_x),cluster_z)))>1.0e-13_dp) &
    error stop 'MRSF cluster response violates the subspace gauge'
  if(residual>1.0e-11_dp) &
    error stop 'MRSF cluster response residual is too large'

  ! A nontrivial rotation of unequal clustered roots makes X^T A X dense.
  ! The response and effective derivative must transform covariantly, while
  ! the derivative of the cluster projector remains invariant.
  angle=0.37_dp
  c=cos(angle); s=sin(angle)
  rotation=reshape([c,s,-s,c],[2,2])
  base_z=cluster_z
  base_w=cluster_w
  projector_deriv=matmul(base_z,transpose(cluster_x))+ &
    matmul(cluster_x,transpose(base_z))
  rotated_x=matmul(cluster_x,rotation)
  rotated_dax=matmul(cluster_dax,rotation)
  call solve_mrsf_tda_cluster_response_matrix_free( &
    apply_test_cluster_operator,cluster_e,rotated_x,rotated_dax,rotated_z, &
    rotated_w,residual,status,tol=1.0e-13_dp,maxit=100,restart=10)
  if(status/=0) error stop 'rotated MRSF cluster response failed'
  if(maxval(abs(rotated_z-matmul(base_z,rotation)))>1.0e-11_dp) &
    error stop 'MRSF cluster response is not basis covariant'
  if(maxval(abs(rotated_w-matmul(transpose(rotation), &
      matmul(base_w,rotation))))>1.0e-12_dp) &
    error stop 'MRSF effective derivative is not basis covariant'
  rotated_projector_deriv=matmul(rotated_z,transpose(rotated_x))+ &
    matmul(rotated_x,transpose(rotated_z))
  if(maxval(abs(rotated_projector_deriv-projector_deriv))>1.0e-11_dp) &
    error stop 'MRSF cluster projector derivative is basis dependent'

  ! The cluster equation remains regular as its two internal eigenvalues
  ! coalesce, because no inverse internal root gap is formed.
  cluster_a(2,2)=1.0_dp+1.0e-10_dp
  cluster_e=[1.0_dp,1.0_dp+1.0e-10_dp]
  call solve_mrsf_tda_cluster_response_matrix_free( &
    apply_test_cluster_operator,cluster_e,cluster_x,cluster_dax,near_limit_z, &
    cluster_w,residual,status,tol=1.0e-13_dp,maxit=100,restart=10)
  if(status/=0) error stop 'near-coalescent MRSF cluster response failed'
  cluster_a(2,2)=1.0_dp
  cluster_e=[1.0_dp,1.0_dp]
  call solve_mrsf_tda_cluster_response_matrix_free( &
    apply_test_cluster_operator,cluster_e,cluster_x,cluster_dax,cluster_z, &
    cluster_w,residual,status,tol=1.0e-13_dp,maxit=100,restart=10)
  if(status/=0) error stop 'degenerate MRSF cluster response failed'
  if(maxval(abs(cluster_z-near_limit_z))>1.0e-9_dp) &
    error stop 'MRSF cluster response is discontinuous at degeneracy'

  ! A one-vector cluster must reduce to the established isolated-root GMRES
  ! response and its scalar Hellmann--Feynman derivative.
  isolated_x(:,1)=mrsf_x
  isolated_dax(:,1)=mrsf_dax(:,1)
  call solve_mrsf_tda_cluster_response_matrix_free( &
    apply_test_mrsf_operator,[2.0_dp],isolated_x,isolated_dax,isolated_z, &
    isolated_w,residual,status,tol=1.0e-13_dp,maxit=30,restart=3)
  if(status/=0) error stop 'one-root MRSF cluster response failed'
  if(maxval(abs(isolated_z(:,1)-mrsf_dx(:,1)))>1.0e-12_dp .or. &
     abs(isolated_w(1,1)-mrsf_domega(1))>1.0e-13_dp) &
    error stop 'one-root cluster response differs from isolated response'

  ! Fail closed when the supplied cluster basis or energy spectrum does not
  ! describe the invariant operator subspace.
  rotated_x=cluster_x
  rotated_x(1,1)=2.0_dp
  call solve_mrsf_tda_cluster_response_matrix_free( &
    apply_test_cluster_operator,[1.0_dp,1.0_dp],rotated_x,cluster_dax, &
    rotated_z,rotated_w,residual,status)
  if(status/=-2) error stop 'nonorthonormal MRSF cluster was not rejected'
  call solve_mrsf_tda_cluster_response_matrix_free( &
    apply_test_cluster_operator,[1.0_dp,1.2_dp],cluster_x,cluster_dax, &
    rotated_z,rotated_w,residual,status)
  if(status/=-5) error stop 'inconsistent MRSF cluster spectrum was not rejected'

contains

  subroutine apply_test_mrsf_operator(vector,result,operator_status)
    real(kind=dp), intent(in) :: vector(:)
    real(kind=dp), intent(out) :: result(:)
    integer, intent(out) :: operator_status
    if(size(vector)/=3 .or. size(result)/=3) then
      operator_status=-1
      result=0.0_dp
      return
    end if
    result=matmul(mrsf_a,vector)
    operator_status=0
  end subroutine apply_test_mrsf_operator

  subroutine apply_test_cluster_operator(vector,result,operator_status)
    real(kind=dp), intent(in) :: vector(:)
    real(kind=dp), intent(out) :: result(:)
    integer, intent(out) :: operator_status
    if(size(vector)/=5 .or. size(result)/=5) then
      operator_status=-1
      result=0.0_dp
      return
    end if
    result=matmul(cluster_a,vector)
    operator_status=0
  end subroutine apply_test_cluster_operator

  subroutine apply_test_z_orbital_hessian(vector,result,operator_status)
    real(kind=dp), intent(in) :: vector(:)
    real(kind=dp), intent(out) :: result(:)
    integer, intent(out) :: operator_status
    if(size(vector)/=3 .or. size(result)/=3) then
      operator_status=-1
      result=0.0_dp
      return
    end if
    result=matmul(zorb,vector)
    operator_status=0
  end subroutine apply_test_z_orbital_hessian

end program test_tdhf_hessian_response
