module tdhf_hessian_rhs_mod

  use precision, only: dp

  implicit none

  private
  public :: build_tdhf_amplitude_derivative_actions
  logical, parameter :: enable_tddft_kxc_ground_response = .true.

contains

!###############################################################################

  subroutine build_tdhf_amplitude_derivative_actions(infos, umat, eps_deriv, dground, &
                                                       u0, v0, dambu, dapbv)
    ! Analytic nuclear derivatives d(A-B)U and d(A+B)V for closed-shell TDHF.
    ! The four terms are the derivative-ERI contraction, outer MO rotation,
    ! response to the differentiated transition density, and orbital-energy
    ! derivative, following GAMESS TDHSGC.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use tdhf_response_operator_mod, only: apply_tdhf_ao_operators
    use tdhf_hessian_z_rhs_mod, only: explicit_channel_derivative_matrix, BLK_OV
    use tdhf_hessian_response_mod, only: assemble_tdhf_sigma_derivative
    use messages, only: show_message, WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: umat(:,:,:), eps_deriv(:,:), dground(:,:,:), u0(:), v0(:)
    real(kind=dp), intent(out) :: dambu(:,:), dapbv(:,:)

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: mo(:,:)
    real(kind=dp), allocatable, target :: pu(:,:), pv(:,:), put(:,:), pvt(:,:)
    real(kind=dp), allocatable :: dp_u(:,:,:), dp_v(:,:,:), densities(:,:,:)
    real(kind=dp), allocatable :: amb_ao(:,:,:), apb_ao(:,:,:)
    real(kind=dp), allocatable :: gminus(:,:), gplus(:,:), deri_m(:,:), deri_p(:,:)
    real(kind=dp), allocatable :: inner_m(:,:), inner_p(:,:), work(:,:), dmo(:,:,:)
    real(kind=dp), allocatable :: deri_full_m(:,:,:),deri_full_p(:,:,:)
    integer :: c, k, nbf, ncart, ncoord, nocc, nvir, nexc, status

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    nexc = nocc*nvir
    ncart = 3*size(basis%atoms%xyz,2)
    ncoord = size(umat,3)
    if (ncoord /= ncart .or. any(shape(umat) /= [nbf,nbf,ncart]) .or. &
        any(shape(dground) /= [nbf,nbf,ncart]) .or. &
        any(shape(eps_deriv) /= [nbf,ncart]) .or. size(u0) /= nexc .or. &
        size(v0) /= nexc .or. any(shape(dambu) /= [nexc,ncart]) .or. &
        any(shape(dapbv) /= [nexc,ncart])) then
      call show_message('TDHF amplitude-derivative arrays have incompatible shapes.', WITH_ABORT)
    end if
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)

    allocate(pu(nbf,nbf), pv(nbf,nbf), put(nbf,nbf), pvt(nbf,nbf))
    call transition_density(mo, nocc, u0, pu)
    call transition_density(mo, nocc, v0, pv)
    put = transpose(pu)
    pvt = transpose(pv)

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
    allocate(deri_full_m(nbf,nbf,ncart),deri_full_p(nbf,nbf,ncart))
    ! Only the occupied-virtual block is read (deri_m/deri_p below).
    call explicit_channel_derivative_matrix(infos,mo,pu,-1,deri_full_m,blocks=BLK_OV)
    call explicit_channel_derivative_matrix(infos,mo,pv,+1,deri_full_p,blocks=BLK_OV)
    if (infos%control%hamilton == 20 .and. enable_tddft_kxc_ground_response) then
      block
        use mod_dft, only: dft_initialize,dftclean
        use mod_dft_molgrid, only: dft_grid_t
        use mod_dft_gridint_gxc, only: tddft_gxc
        type(dft_grid_t) :: grid
        real(dp),allocatable,target :: dx(:,:,:)
        real(dp),allocatable :: fx(:,:,:),pvs(:,:),kxc_mo(:,:),kxc_work(:,:)
        integer :: kk
        allocate(dx(nbf,nbf,3*ncart),fx(nbf,nbf,3*ncart),pvs(nbf,nbf), &
          kxc_mo(nbf,nbf),kxc_work(nbf,nbf),source=0.0_dp)
        pvs=0.5_dp*(pv+transpose(pv))
        do kk=1,ncart
          ! dground is spin summed, whereas the restricted grid consumer
          ! accepts a one-spin density.  The quadratic Gxc polarization must
          ! therefore use dground/2 at this API boundary.
          dx(:,:,3*kk-2)=pvs+0.5_dp*dground(:,:,kk)
          dx(:,:,3*kk-1)=pvs
          dx(:,:,3*kk)=0.5_dp*dground(:,:,kk)
        end do
        call dft_initialize(infos,basis,grid)
        call tddft_gxc(basis,grid,.true.,mo,fx,dx,3*ncart,1.0e-14_dp,infos)
        call dftclean(infos)
        do kk=1,ncart
          fx(:,:,3*kk-2)=fx(:,:,3*kk-2)-fx(:,:,3*kk-1)-fx(:,:,3*kk)
          call ao_to_mo(fx(:,:,3*kk-2),mo,kxc_mo,kxc_work)
          deri_full_p(:,:,kk)=deri_full_p(:,:,kk)+kxc_mo
        end do
        deallocate(dx,fx,pvs,kxc_mo,kxc_work)
      end block
    end if
    do k = 1, ncart
      deri_m(:,k)=reshape(deri_full_m(1:nocc,nocc+1:,k),[nexc])
      deri_p(:,k)=reshape(deri_full_p(1:nocc,nocc+1:,k),[nexc])
      call ao_to_ov(amb_ao(:,:,k+2), mo, nocc, inner_m(:,k))
      call ao_to_ov(apb_ao(:,:,ncart+k+2), mo, nocc, inner_p(:,k))
    end do
    call assemble_tdhf_sigma_derivative(u0, nocc, umat, eps_deriv, gminus, &
      deri_m, inner_m, dambu, status)
    if (status /= 0) call show_message('Failed to assemble d(A-B)U.', WITH_ABORT)
    call assemble_tdhf_sigma_derivative(v0, nocc, umat, eps_deriv, gplus, &
      deri_p, inner_p, dapbv, status)
    if (status /= 0) call show_message('Failed to assemble d(A+B)V.', WITH_ABORT)

    ! d(A-B)U and d(A+B)V are always the assembled fixed-AO-frame actions.
    ! A former HF branch replaced them by moving-frame derivative-gradient
    ! polarizations whenever every U was symmetric (H2/STO-3G, where parity
    ! makes U_ov vanish); it only fitted the matching shortcut in
    ! tdhf_hessian, and both together gave 0.20 instead of 1.159 Eh/bohr^2
    ! for the H2 excited-state zz element.
    deallocate(pu, pv, put, pvt, dp_u, dp_v, &
      densities, amb_ao, apb_ao, gminus, gplus, deri_m, deri_p, inner_m, &
      inner_p, work, dmo, deri_full_m, deri_full_p)
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
