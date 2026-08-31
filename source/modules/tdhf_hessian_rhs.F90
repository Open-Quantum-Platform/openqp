module tdhf_hessian_rhs_mod

  use precision, only: dp

  implicit none

  private
  public :: build_tdhf_amplitude_derivative_actions

contains

!###############################################################################

  subroutine build_tdhf_amplitude_derivative_actions(infos, umat, eps_deriv, &
                                                       u0, v0, dambu, dapbv)
    ! Analytic nuclear derivatives d(A-B)U and d(A+B)V for closed-shell TDHF.
    ! The four terms are the derivative-ERI contraction, outer MO rotation,
    ! response to the differentiated transition density, and orbital-energy
    ! derivative, following GAMESS TDHSGC.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use fock_deriv_mod, only: fock_deriv_matrix_general
    use tdhf_response_operator_mod, only: apply_tdhf_ao_operators
    use tdhf_hessian_response_mod, only: assemble_tdhf_sigma_derivative
    use messages, only: show_message, WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: umat(:,:,:), eps_deriv(:,:), u0(:), v0(:)
    real(kind=dp), intent(out) :: dambu(:,:), dapbv(:,:)

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: mo(:,:)
    real(kind=dp), allocatable, target :: pu(:,:), pv(:,:), put(:,:), pvt(:,:)
    real(kind=dp), allocatable :: fpu(:,:,:,:), fput(:,:,:,:), fpv(:,:,:,:), fpvt(:,:,:,:)
    real(kind=dp), allocatable :: dp_u(:,:,:), dp_v(:,:,:), densities(:,:,:)
    real(kind=dp), allocatable :: amb_ao(:,:,:), apb_ao(:,:,:)
    real(kind=dp), allocatable :: gminus(:,:), gplus(:,:), deri_m(:,:), deri_p(:,:)
    real(kind=dp), allocatable :: inner_m(:,:), inner_p(:,:), work(:,:), dmo(:,:,:)
    integer :: c, i, k, nbf, ncart, ncoord, nocc, nvir, nexc, status

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    nexc = nocc*nvir
    ncart = 3*size(basis%atoms%xyz,2)
    ncoord = size(umat,3)
    if (infos%control%hamilton >= 20) then
      call show_message('TDDFT amplitude derivatives require the XC skeleton derivative.', &
                        WITH_ABORT)
    end if
    if (ncoord /= ncart .or. any(shape(umat) /= [nbf,nbf,ncart]) .or. &
        any(shape(eps_deriv) /= [nbf,ncart]) .or. size(u0) /= nexc .or. &
        size(v0) /= nexc .or. any(shape(dambu) /= [nexc,ncart]) .or. &
        any(shape(dapbv) /= [nexc,ncart])) then
      call show_message('TDHF amplitude-derivative arrays have incompatible shapes.', WITH_ABORT)
    end if
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)

    allocate(pu(nbf,nbf), pv(nbf,nbf), put(nbf,nbf), pvt(nbf,nbf), &
      fpu(nbf,nbf,3,ncart/3), fput(nbf,nbf,3,ncart/3), &
      fpv(nbf,nbf,3,ncart/3), fpvt(nbf,nbf,3,ncart/3))
    call transition_density(mo, nocc, u0, pu)
    call transition_density(mo, nocc, v0, pv)
    put = transpose(pu)
    pvt = transpose(pv)
    call fock_deriv_matrix_general(infos, basis, pu, 1.0_dp, fpu)
    call fock_deriv_matrix_general(infos, basis, put, 1.0_dp, fput)
    call fock_deriv_matrix_general(infos, basis, pv, 1.0_dp, fpv)
    call fock_deriv_matrix_general(infos, basis, pvt, 1.0_dp, fpvt)

    allocate(dp_u(nbf,nbf,ncart), dp_v(nbf,nbf,ncart), dmo(nbf,nbf,ncart))
    do k = 1, ncart
      dmo(:,:,k) = matmul(mo, umat(:,:,k))
      call transition_density_derivative(mo, dmo(:,:,k), nocc, u0, dp_u(:,:,k))
      call transition_density_derivative(mo, dmo(:,:,k), nocc, v0, dp_v(:,:,k))
    end do

    allocate(densities(nbf,nbf,2*ncart+2), amb_ao(nbf,nbf,2*ncart+2), &
             apb_ao(nbf,nbf,2*ncart+2))
    densities(:,:,1) = pu
    densities(:,:,2) = pv
    densities(:,:,3:ncart+2) = dp_u
    densities(:,:,ncart+3:2*ncart+2) = dp_v
    call apply_tdhf_ao_operators(infos, densities, amb_ao, apb_ao)

    allocate(gminus(nbf,nbf), gplus(nbf,nbf), deri_m(nexc,ncart), &
      deri_p(nexc,ncart), inner_m(nexc,ncart), inner_p(nexc,ncart), &
      work(nbf,nbf), source=0.0_dp)
    call ao_to_mo(amb_ao(:,:,1), mo, gminus, work)
    call ao_to_mo(apb_ao(:,:,2), mo, gplus, work)
    do k = 1, ncart
      c = mod(k-1,3)+1
      i = (k-1)/3+1
      ! OpenQP stores the transition density as C_occ Z C_vir^T, the
      ! transpose of the GAMESS TDHPTD convention; hence the A-B sign below.
      work = 2.0_dp*(fpu(:,:,c,i)-fput(:,:,c,i))
      call ao_to_ov(work, mo, nocc, deri_m(:,k))
      work = 2.0_dp*(fpv(:,:,c,i)+fpvt(:,:,c,i))
      call ao_to_ov(work, mo, nocc, deri_p(:,k))
      call ao_to_ov(amb_ao(:,:,k+2), mo, nocc, inner_m(:,k))
      call ao_to_ov(apb_ao(:,:,ncart+k+2), mo, nocc, inner_p(:,k))
    end do
    call assemble_tdhf_sigma_derivative(u0, nocc, umat, eps_deriv, gminus, &
      deri_m, inner_m, dambu, status)
    if (status /= 0) call show_message('Failed to assemble d(A-B)U.', WITH_ABORT)
    call assemble_tdhf_sigma_derivative(v0, nocc, umat, eps_deriv, gplus, &
      deri_p, inner_p, dapbv, status)
    if (status /= 0) call show_message('Failed to assemble d(A+B)V.', WITH_ABORT)

    deallocate(pu, pv, put, pvt, fpu, fput, fpv, fpvt, dp_u, dp_v, &
      densities, amb_ao, apb_ao, gminus, gplus, deri_m, deri_p, inner_m, &
      inner_p, work, dmo)
  contains
    subroutine transition_density(coeff, no, z, p)
      real(kind=dp), intent(in) :: coeff(:,:), z(:)
      integer, intent(in) :: no
      real(kind=dp), intent(out) :: p(:,:)
      real(kind=dp) :: zmat(no,size(coeff,1)-no)
      zmat = reshape(z, shape(zmat))
      p = matmul(matmul(coeff(:,1:no), zmat), transpose(coeff(:,no+1:)))
    end subroutine transition_density

    subroutine transition_density_derivative(coeff, dcoeff, no, z, dpmat)
      real(kind=dp), intent(in) :: coeff(:,:), dcoeff(:,:), z(:)
      integer, intent(in) :: no
      real(kind=dp), intent(out) :: dpmat(:,:)
      real(kind=dp) :: zmat(no,size(coeff,1)-no)
      zmat = reshape(z, shape(zmat))
      dpmat = matmul(matmul(dcoeff(:,1:no), zmat), transpose(coeff(:,no+1:))) &
            + matmul(matmul(coeff(:,1:no), zmat), transpose(dcoeff(:,no+1:)))
    end subroutine transition_density_derivative

    subroutine ao_to_mo(ao, coeff, mo_mat, tmp)
      real(kind=dp), intent(in) :: ao(:,:), coeff(:,:)
      real(kind=dp), intent(out) :: mo_mat(:,:), tmp(:,:)
      tmp = matmul(ao, coeff)
      mo_mat = matmul(transpose(coeff), tmp)
    end subroutine ao_to_mo

    subroutine ao_to_ov(ao, coeff, no, ov)
      real(kind=dp), intent(in) :: ao(:,:), coeff(:,:)
      integer, intent(in) :: no
      real(kind=dp), intent(out) :: ov(:)
      real(kind=dp) :: tmp(size(coeff,1),size(coeff,1)-no)
      real(kind=dp) :: block(no,size(coeff,1)-no)
      tmp = matmul(ao, coeff(:,no+1:))
      block = matmul(transpose(coeff(:,1:no)), tmp)
      ov = reshape(block, [size(ov)])
    end subroutine ao_to_ov
  end subroutine build_tdhf_amplitude_derivative_actions

end module tdhf_hessian_rhs_mod
