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
    use tdhf_hessian_z_rhs_mod, only: explicit_channel_derivative_matrix
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
    real(kind=dp), allocatable :: orbital_m(:,:),orbital_p(:,:)
    integer :: c, k, nbf, ncart, ncoord, nocc, nvir, nexc, status
    logical :: zero_orbital_connection

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
    call explicit_channel_derivative_matrix(infos,mo,pu,-1,deri_full_m)
    call explicit_channel_derivative_matrix(infos,mo,pv,+1,deri_full_p)
    if (infos%control%hamilton == 20 .and. enable_tddft_kxc_ground_response) then
      block
        use mod_dft, only: dft_initialize,dftclean
        use mod_dft_molgrid, only: dft_grid_t
        use mod_dft_gridint_gxc, only: tddft_gxc
        type(dft_grid_t) :: grid
        real(dp),allocatable,target :: dx(:,:,:)
        real(dp),allocatable :: fx(:,:,:),pvs(:,:)
        integer :: kk
        allocate(dx(nbf,nbf,3*ncart),fx(nbf,nbf,3*ncart),pvs(nbf,nbf),source=0.0_dp)
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
          fx(:,:,3*kk-2)=0.5_dp*(fx(:,:,3*kk-2)-fx(:,:,3*kk-1)-fx(:,:,3*kk))
          deri_full_p(:,:,kk)=deri_full_p(:,:,kk)+fx(:,:,3*kk-2)
        end do
        deallocate(dx,fx,pvs)
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

    allocate(orbital_m(nexc,ncart),orbital_p(nexc,ncart))
    call orbital_channel_derivative(u0,-1,orbital_m)
    call orbital_channel_derivative(v0,+1,orbital_p)
    ! Never suppress the KS orbital/density connection.  The DFT kxc and XC
    ! response rows require dD even when the antisymmetric part of U vanishes.
    zero_orbital_connection = infos%control%hamilton < 20
    do k=1,ncart
      if (maxval(abs(umat(:,:,k)-transpose(umat(:,:,k)))) > 1.0e-10_dp) &
        zero_orbital_connection = .false.
    end do
    do k=1,ncart
      ! The derivative-gradient polarization is already expressed in the
      ! moving AO frame used by the analytic gradient.  It therefore contains
      ! the metric/orbital-energy connection terms which would otherwise be
      ! added explicitly in a fixed AO frame.  Adding C*U and eps^R here would
      ! count those contributions twice.
      if (zero_orbital_connection) then
        dambu(:,k)=deri_m(:,k)+orbital_m(:,k)
        dapbv(:,k)=deri_p(:,k)+orbital_p(:,k)
      end if
    end do

    deallocate(pu, pv, put, pvt, dp_u, dp_v, &
      densities, amb_ao, apb_ao, gminus, gplus, deri_m, deri_p, inner_m, &
      inner_p, work, dmo, deri_full_m, deri_full_p, orbital_m, orbital_p)
  contains
    subroutine orbital_channel_derivative(z,channel,result)
      real(dp),intent(in)::z(:)
      integer,intent(in)::channel
      real(dp),intent(out)::result(:,:)
      real(dp),allocatable::ct(:,:),dc(:,:),dens(:,:,:),am(:,:,:),ap(:,:,:),vals(:,:),tmp(:,:),momat(:,:)
      real(dp),parameter::steps(4)=[0.5_dp,-0.5_dp,1.0_dp,-1.0_dp]
      integer::kk,s
      allocate(ct(nbf,nbf),dc(nbf,nbf),dens(nbf,nbf,4),am(nbf,nbf,4),ap(nbf,nbf,4), &
        vals(nexc,4),tmp(nbf,nbf),momat(nbf,nbf))
      do kk=1,ncart
        ! Only the anti-Hermitian MO rotation is an independent response in
        ! the moving-AO frame.  The symmetric metric connection is already
        ! contained in the explicit derivative-gradient polarization.
        dc=0.5_dp*matmul(mo,umat(:,:,kk)-transpose(umat(:,:,kk)))
        do s=1,4
          ct=mo+steps(s)*dc
          call transition_density(ct,nocc,z,dens(:,:,s))
        end do
        call apply_tdhf_ao_operators(infos,dens,am,ap)
        do s=1,4
          ct=mo+steps(s)*dc
          if(channel<0) then
            call ao_to_mo(am(:,:,s),ct,momat,tmp)
          else
            call ao_to_mo(ap(:,:,s),ct,momat,tmp)
          end if
          vals(:,s)=reshape(momat(1:nocc,nocc+1:),[nexc])
        end do
        result(:,kk)=(4.0_dp*(vals(:,1)-vals(:,2))-(vals(:,3)-vals(:,4))/2.0_dp)/3.0_dp
      end do
      deallocate(ct,dc,dens,am,ap,vals,tmp,momat)
    end subroutine orbital_channel_derivative

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
