module tdhf_hessian_orbital_mod

  use precision, only: dp

  implicit none

  private
  public :: build_tdhf_ground_orbital_response

contains

!###############################################################################

  subroutine build_tdhf_ground_orbital_response(infos, sx_mo, umat, eps_deriv, dp_ao)
    ! Closed-shell HF orbital response to all Cartesian nuclear perturbations.
    ! This supplies the S^K, U^K, epsilon^K, and P^K quantities entering the
    ! differentiated TD response operators.  DFT requires the additional
    ! moving-grid XC derivative and is rejected until that term is supplied.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, &
      OQP_VEC_MO_A, OQP_E_MO_A
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix
    use ecp_tool, only: ecp_deriv_ints
    use fock_deriv_mod, only: fock_deriv_matrix
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve
    use tdhf_hessian_response_mod, only: complete_rhf_orbital_response
    use messages, only: show_message, WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: sx_mo(:,:,:), umat(:,:,:)
    real(kind=dp), intent(out) :: eps_deriv(:,:), dp_ao(:,:,:)

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dmat(:), mo(:,:), eps(:)
    real(kind=dp), allocatable, target :: pfull(:,:)
    real(kind=dp), allocatable :: ds(:,:,:,:), dt(:,:,:,:), dv(:,:,:,:), decp(:,:,:,:)
    real(kind=dp), allocatable :: dg(:,:,:,:), f0ao(:,:), f0mo(:,:), gd0mo(:,:)
    real(kind=dp), allocatable :: d0(:,:), d0pack(:,:), gpack(:,:), gfull(:,:)
    real(kind=dp), allocatable :: probe1(:,:), bvec(:,:), uvec(:,:)
    real(kind=dp), allocatable :: dmo_occ(:,:)
    integer :: a, c, i, ia, k, mu, nbf, nbf2, natom, ncart, nocc, nvir, status
    logical :: converged

    if (infos%control%scftype /= 1 .or. infos%mol_prop%mult /= 1) then
      call show_message('TD Hessian orbital response requires a closed-shell restricted reference.', &
                        WITH_ABORT)
    end if
    if (infos%control%hamilton >= 20) then
      call show_message('TDDFT orbital response requires the moving-grid XC Fock derivative.', &
                        WITH_ABORT)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz,2)
    ncart = 3*natom
    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    if (any(shape(sx_mo) /= [nbf,nbf,ncart]) .or. &
        any(shape(umat) /= [nbf,nbf,ncart]) .or. &
        any(shape(eps_deriv) /= [nbf,ncart]) .or. &
        any(shape(dp_ao) /= [nbf,nbf,ncart])) then
      call show_message('TD Hessian orbital-response output has the wrong shape.', WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, eps)
    allocate(pfull(nbf,nbf)); call unpack_matrix(dmat, pfull)
    allocate(ds(nbf,nbf,3,natom), dt(nbf,nbf,3,natom), &
             dv(nbf,nbf,3,natom), decp(nbf,nbf,3,natom), &
             dg(nbf,nbf,3,natom))
    call der_overlap_matrix(basis, ds)
    call der_kinetic_matrix(basis, dt)
    call der_nucattr_matrix(basis, basis%atoms%xyz, &
      basis%atoms%zn-basis%ecp_zn_num, dv)
    do k = 1, natom
      do c = 1, 3
        do i = 1, nbf
          do mu = 1, nbf
            ds(mu,i,c,k) = ds(mu,i,c,k)*basis%bfnrm(mu)*basis%bfnrm(i)
            dt(mu,i,c,k) = dt(mu,i,c,k)*basis%bfnrm(mu)*basis%bfnrm(i)
            dv(mu,i,c,k) = dv(mu,i,c,k)*basis%bfnrm(mu)*basis%bfnrm(i)
          end do
        end do
      end do
    end do
    call ecp_deriv_ints(basis, basis%atoms%xyz, decp)
    dv = dv + decp
    call fock_deriv_matrix(infos, basis, pfull, 1.0_dp, dg)

    allocate(f0ao(nbf,nbf), f0mo(nbf,nbf), gd0mo(nbf,nbf), &
             d0(nbf,nbf), d0pack(nbf2,1), gpack(nbf2,1), gfull(nbf,nbf), &
             probe1(nbf,nbf), bvec(nocc*nvir,ncart), &
             uvec(nocc*nvir,ncart), dmo_occ(nbf,nocc), source=0.0_dp)
    ia = 0
    do k = 1, natom
      do c = 1, 3
        ia = ia + 1
        call ao_to_mo(ds(:,:,c,k), mo, sx_mo(:,:,ia), probe1)
        f0ao = dt(:,:,c,k) + dv(:,:,c,k) + dg(:,:,c,k)
        call ao_to_mo(f0ao, mo, f0mo, probe1)
        d0 = -2.0_dp*matmul(matmul(mo(:,1:nocc), sx_mo(1:nocc,1:nocc,ia)), &
                            transpose(mo(:,1:nocc)))
        call pack_matrix(d0, d0pack(:,1))
        gpack = 0.0_dp
        call fock_jk(basis, d=d0pack, f=gpack, scale_exch=1.0_dp, infos=infos)
        call unpack_symmetric(gpack(:,1), gfull, nbf)
        call ao_to_mo(gfull, mo, gd0mo, probe1)
        do a = 1, nvir
          do i = 1, nocc
            bvec((a-1)*nocc+i,ia) = -f0mo(i,nocc+a) &
              + eps(i)*sx_mo(i,nocc+a,ia) - gd0mo(i,nocc+a)
          end do
        end do
      end do
    end do

    call cphf_solve(infos, ncart, bvec, uvec, converged=converged)
    if (.not. converged) then
      call show_message('Ground-state CPHF did not converge for the TD Hessian.', WITH_ABORT)
    end if
    do ia = 1, ncart
      call complete_rhf_orbital_response(mo, nocc, sx_mo(:,:,ia), uvec(:,ia), &
        umat(:,:,ia), dmo_occ, dp_ao(:,:,ia), status)
      if (status /= 0) call show_message('Failed to complete the TD Hessian orbital response.', &
                                         WITH_ABORT)
      c = mod(ia-1,3)+1
      k = (ia-1)/3+1
      f0ao = dt(:,:,c,k) + dv(:,:,c,k) + dg(:,:,c,k)
      call pack_matrix(dp_ao(:,:,ia), d0pack(:,1))
      gpack = 0.0_dp
      call fock_jk(basis, d=d0pack, f=gpack, scale_exch=1.0_dp, infos=infos)
      call unpack_symmetric(gpack(:,1), gfull, nbf)
      call ao_to_mo(f0ao+gfull, mo, f0mo, probe1)
      eps_deriv(:,ia) = [(f0mo(i,i)-eps(i)*sx_mo(i,i,ia), i=1,nbf)]
    end do

    deallocate(pfull, ds, dt, dv, decp, dg, f0ao, f0mo, gd0mo, d0, &
      d0pack, gpack, gfull, probe1, bvec, uvec, dmo_occ)
  contains
    subroutine ao_to_mo(ao, coeff, mo_mat, work)
      real(kind=dp), intent(in) :: ao(:,:), coeff(:,:)
      real(kind=dp), intent(out) :: mo_mat(:,:), work(:,:)
      work = matmul(ao, coeff)
      mo_mat = matmul(transpose(coeff), work)
    end subroutine ao_to_mo

    subroutine unpack_symmetric(gpk, gfu, n)
      real(kind=dp), intent(in) :: gpk(:)
      real(kind=dp), intent(out) :: gfu(:,:)
      integer, intent(in) :: n
      integer :: ii, ij, jj
      ij = 0
      do ii = 1, n
        do jj = 1, ii
          ij = ij + 1
          gfu(ii,jj) = gpk(ij)
          gfu(jj,ii) = gpk(ij)
        end do
      end do
    end subroutine unpack_symmetric
  end subroutine build_tdhf_ground_orbital_response

end module tdhf_hessian_orbital_mod
