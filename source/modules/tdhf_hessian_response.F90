module tdhf_hessian_response_mod

  use precision, only: dp

  implicit none

  private
  public :: solve_tdhf_amplitude_response
  public :: solve_tdhf_z_response
  public :: complete_rhf_orbital_response

contains

!###############################################################################

  subroutine complete_rhf_orbital_response(mo, nocc, sx_mo, uov, umat, dmo_occ, dp_ao, status)
    ! Complete the occupied-virtual CPHF amplitudes into a full first-order
    ! MO connection.  The occupied-occupied and virtual-virtual rotations use
    ! the symmetric gauge, U_pq=-S^x_pq/2.  The transposed occupied-virtual
    ! block follows from differentiated orthonormality,
    !
    !                 U + U^T + S^x = 0.
    !
    ! The returned AO density is the spin-summed closed-shell response.

    real(kind=dp), intent(in) :: mo(:,:), sx_mo(:,:), uov(:)
    integer, intent(in) :: nocc
    real(kind=dp), intent(out) :: umat(:,:), dmo_occ(:,:), dp_ao(:,:)
    integer, intent(out) :: status

    integer :: a, i, ia, nbf, nvir

    nbf = size(mo, 1)
    nvir = nbf - nocc
    status = 0
    umat = 0.0_dp
    dmo_occ = 0.0_dp
    dp_ao = 0.0_dp

    if (nbf <= 0 .or. size(mo,2) /= nbf .or. nocc <= 0 .or. nocc >= nbf .or. &
        any(shape(sx_mo) /= [nbf,nbf]) .or. size(uov) /= nocc*nvir .or. &
        any(shape(umat) /= [nbf,nbf]) .or. &
        any(shape(dmo_occ) /= [nbf,nocc]) .or. &
        any(shape(dp_ao) /= [nbf,nbf])) then
      status = -1
      return
    end if

    umat(1:nocc,1:nocc) = -0.5_dp*sx_mo(1:nocc,1:nocc)
    umat(nocc+1:nbf,nocc+1:nbf) = &
      -0.5_dp*sx_mo(nocc+1:nbf,nocc+1:nbf)
    ia = 0
    do a = 1, nvir
      do i = 1, nocc
        ia = ia + 1
        umat(nocc+a,i) = uov(ia)
        umat(i,nocc+a) = -sx_mo(i,nocc+a) - uov(ia)
      end do
    end do

    dmo_occ = matmul(mo, umat(:,1:nocc))
    dp_ao = 2.0_dp*(matmul(dmo_occ, transpose(mo(:,1:nocc))) + &
                    matmul(mo(:,1:nocc), transpose(dmo_occ)))
  end subroutine complete_rhf_orbital_response

!###############################################################################

  subroutine solve_tdhf_amplitude_response(amb, apb, omega, u0, v0, &
                                            dambu, dapbv, du, dv, domega, &
                                            residual_max, status, tol, maxit)
    ! Solve the differentiated full-response TDHF/TDDFT equations
    !
    !   (A-B) dU - omega dV = domega V - d(A-B) U,
    !   (A+B) dV - omega dU = domega U - d(A+B) V.
    !
    ! The coupled matrix has the null vector (U,V).  Each right-hand side
    ! is projected perpendicular to that vector.  After the projected
    ! conjugate-gradient solution, the remaining null-vector component is
    ! chosen so that d(U^T V)=0.  This is the same gauge condition used by
    ! the GAMESS TDHBLD implementation.

    real(kind=dp), intent(in) :: amb(:,:), apb(:,:)
    real(kind=dp), intent(in) :: omega, u0(:), v0(:)
    real(kind=dp), intent(in) :: dambu(:,:), dapbv(:,:)
    real(kind=dp), intent(out) :: du(:,:), dv(:,:), domega(:)
    real(kind=dp), intent(out) :: residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    real(kind=dp), allocatable :: b(:), x(:), r(:), p(:), q(:), nullvec(:)
    real(kind=dp) :: alpha, beta, bnorm, gauge_denominator, gauge_value
    real(kind=dp) :: n2, pq, rr, rr_new, solve_tol, uvnorm
    integer :: i, iteration, ncoord, nexc, niter

    nexc = size(u0)
    ncoord = size(dambu, 2)
    status = 0
    residual_max = 0.0_dp
    domega = 0.0_dp
    du = 0.0_dp
    dv = 0.0_dp

    if (nexc <= 0 .or. size(v0) /= nexc .or. &
        any(shape(amb) /= [nexc, nexc]) .or. &
        any(shape(apb) /= [nexc, nexc]) .or. &
        size(dambu, 1) /= nexc .or. &
        any(shape(dapbv) /= [nexc, ncoord]) .or. &
        any(shape(du) /= [nexc, ncoord]) .or. &
        any(shape(dv) /= [nexc, ncoord]) .or. &
        size(domega) /= ncoord) then
      status = -1
      return
    end if

    solve_tol = 1.0e-12_dp
    if (present(tol)) solve_tol = tol
    niter = 400
    if (present(maxit)) niter = maxit
    if (solve_tol <= 0.0_dp .or. niter <= 0) then
      status = -2
      return
    end if

    uvnorm = dot_product(u0, v0)
    gauge_denominator = 2.0_dp*uvnorm
    if (abs(gauge_denominator) <= tiny(1.0_dp)) then
      status = -3
      return
    end if

    allocate(b(2*nexc), x(2*nexc), r(2*nexc), p(2*nexc), q(2*nexc), &
             nullvec(2*nexc))
    nullvec(1:nexc) = u0
    nullvec(nexc+1:2*nexc) = v0
    n2 = dot_product(nullvec, nullvec)

    do i = 1, ncoord
      domega(i) = (dot_product(u0, dambu(:,i)) + &
                   dot_product(v0, dapbv(:,i)))/(2.0_dp*uvnorm)
      b(1:nexc) = domega(i)*v0 - dambu(:,i)
      b(nexc+1:2*nexc) = domega(i)*u0 - dapbv(:,i)
      call project_from_null(b, nullvec, n2)

      bnorm = sqrt(dot_product(b, b))
      if (bnorm <= tiny(1.0_dp)) cycle

      x = 0.0_dp
      r = b
      p = r
      rr = dot_product(r, r)

      do iteration = 1, niter
        call apply_coupled_matrix(amb, apb, omega, p, q, nexc)
        call project_from_null(q, nullvec, n2)
        pq = dot_product(p, q)
        if (pq <= tiny(1.0_dp)) then
          status = i
          exit
        end if
        alpha = rr/pq
        x = x + alpha*p
        r = r - alpha*q
        rr_new = dot_product(r, r)
        if (sqrt(rr_new) <= solve_tol*bnorm) then
          rr = rr_new
          exit
        end if
        beta = rr_new/rr
        p = r + beta*p
        rr = rr_new
      end do
      if (status /= 0) exit
      if (sqrt(rr) > solve_tol*bnorm) then
        status = i
        exit
      end if

      ! Enforce V^T dU + U^T dV = 0 without changing the equations.
      gauge_value = dot_product(v0, x(1:nexc)) + &
                    dot_product(u0, x(nexc+1:2*nexc))
      x = x - (gauge_value/gauge_denominator)*nullvec
      du(:,i) = x(1:nexc)
      dv(:,i) = x(nexc+1:2*nexc)

      call apply_coupled_matrix(amb, apb, omega, x, q, nexc)
      residual_max = max(residual_max, maxval(abs(q-b)))
    end do

    deallocate(b, x, r, p, q, nullvec)

  end subroutine solve_tdhf_amplitude_response

!###############################################################################

  subroutine solve_tdhf_z_response(orbital_hessian, rhs_derivative, &
                                    operator_derivative_z, dz, residual_max, &
                                    status, tol, maxit)
    ! Differentiate H Z = R for every Cartesian perturbation:
    !
    !   H dZ = dR - (dH) Z.

    real(kind=dp), intent(in) :: orbital_hessian(:,:)
    real(kind=dp), intent(in) :: rhs_derivative(:,:), operator_derivative_z(:,:)
    real(kind=dp), intent(out) :: dz(:,:)
    real(kind=dp), intent(out) :: residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    real(kind=dp), allocatable :: b(:), p(:), q(:), r(:), x(:)
    real(kind=dp) :: alpha, beta, bnorm, pq, rr, rr_new, solve_tol
    integer :: i, iteration, ncoord, nexc, niter

    nexc = size(orbital_hessian, 1)
    ncoord = size(rhs_derivative, 2)
    status = 0
    residual_max = 0.0_dp
    dz = 0.0_dp

    if (nexc <= 0 .or. size(orbital_hessian, 2) /= nexc .or. &
        size(rhs_derivative, 1) /= nexc .or. &
        any(shape(operator_derivative_z) /= [nexc, ncoord]) .or. &
        any(shape(dz) /= [nexc, ncoord])) then
      status = -1
      return
    end if

    solve_tol = 1.0e-12_dp
    if (present(tol)) solve_tol = tol
    niter = 400
    if (present(maxit)) niter = maxit
    if (solve_tol <= 0.0_dp .or. niter <= 0) then
      status = -2
      return
    end if

    allocate(b(nexc), p(nexc), q(nexc), r(nexc), x(nexc))
    do i = 1, ncoord
      b = rhs_derivative(:,i) - operator_derivative_z(:,i)
      bnorm = sqrt(dot_product(b, b))
      if (bnorm <= tiny(1.0_dp)) cycle
      x = 0.0_dp
      r = b
      p = r
      rr = dot_product(r, r)

      do iteration = 1, niter
        q = matmul(orbital_hessian, p)
        pq = dot_product(p, q)
        if (pq <= tiny(1.0_dp)) then
          status = i
          exit
        end if
        alpha = rr/pq
        x = x + alpha*p
        r = r - alpha*q
        rr_new = dot_product(r, r)
        if (sqrt(rr_new) <= solve_tol*bnorm) then
          rr = rr_new
          exit
        end if
        beta = rr_new/rr
        p = r + beta*p
        rr = rr_new
      end do
      if (status /= 0) exit
      if (sqrt(rr) > solve_tol*bnorm) then
        status = i
        exit
      end if
      dz(:,i) = x
      residual_max = max(residual_max, &
        maxval(abs(matmul(orbital_hessian, x)-b)))
    end do
    deallocate(b, p, q, r, x)

  end subroutine solve_tdhf_z_response

!###############################################################################

  subroutine apply_coupled_matrix(amb, apb, omega, x, y, nexc)
    real(kind=dp), intent(in) :: amb(:,:), apb(:,:), omega, x(:)
    real(kind=dp), intent(out) :: y(:)
    integer, intent(in) :: nexc

    y(1:nexc) = matmul(amb, x(1:nexc)) - omega*x(nexc+1:2*nexc)
    y(nexc+1:2*nexc) = matmul(apb, x(nexc+1:2*nexc)) - &
                       omega*x(1:nexc)
  end subroutine apply_coupled_matrix

!###############################################################################

  subroutine project_from_null(vector, nullvec, norm2)
    real(kind=dp), intent(inout) :: vector(:)
    real(kind=dp), intent(in) :: nullvec(:), norm2

    vector = vector - nullvec*(dot_product(nullvec, vector)/norm2)
  end subroutine project_from_null

end module tdhf_hessian_response_mod
