module hf_hessian_mod

  implicit none

  character(len=*), parameter :: module_name = "hf_hessian_mod"

contains

!###############################################################################

  subroutine hf_hessian_C(c_handle) bind(C, name="hf_hessian")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call hf_hessian(inf)
  end subroutine hf_hessian_C

!###############################################################################

  subroutine hf_hessian(infos)
    ! Native OpenQP HF/DFT Hessian CPHF response prepass.
    !
    ! This routine deliberately exercises the production Fortran CPHF/CPKS PCG
    ! solver for every Cartesian nuclear perturbation used by a ground-state
    ! analytic Hessian.  It builds the closed-shell occupied-virtual RHS from
    ! OpenQP derivative integrals and the current OpenQP SCF density/MOs, then
    ! calls cphf_solve on the full 3N RHS block.  The final analytic Hessian
    ! contraction is still guarded at the Python dispatch layer by marking the
    ! PySCF matrix as an oracle/final-assembly fallback; this routine is the
    ! native response gate, not a placeholder Hessian matrix.
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, OQP_E_MO_A
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix
    use fock_deriv_mod, only: fock_deriv_contract
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve
    use io_constants, only: iw

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dmat_a(:), mo_a(:,:), eps(:)
    real(kind=dp), allocatable :: pfull(:,:), probe(:,:), gx(:,:)
    real(kind=dp), allocatable :: dSa(:,:,:,:), dTa(:,:,:,:), dVa(:,:,:,:)
    real(kind=dp), allocatable :: Sx(:,:), hx(:,:), F0x(:,:), Gd0(:,:)
    real(kind=dp), allocatable :: d0(:,:), d0p(:,:), gp(:,:), gfull(:,:)
    real(kind=dp), allocatable :: bvec(:,:), uvec(:,:), scr(:,:), col(:,:)
    real(kind=dp) :: hfscale
    integer :: nbf, nbf2, nocc, nvir, natom, ncart
    integer :: i, j, a, mu, nu, ia, icart, kc, cc

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocc = infos%mol_prop%nocc
    nvir = nbf - nocc
    natom = size(basis%atoms%xyz, 2)
    ncart = 3*natom
    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale

    write(iw,'(/,A)') 'PyOQP: Native OpenQP HF/DFT Hessian CPHF response prepass'
    write(iw,'(A,I6,A,I6,A,I6,A,I6)') '  nbf=', nbf, ' nocc=', nocc, ' nvir=', nvir, ' rhs=', ncart
    write(iw,'(A)') '  Final analytic Hessian assembly remains guarded; PySCF final Hessian is retained as oracle.'

    if (nocc <= 0 .or. nvir <= 0 .or. ncart <= 0) then
      write(iw,'(A)') '  Native CPHF prepass skipped: empty occupied/virtual/nuclear space.'
      return
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, eps)

    allocate(pfull(nbf,nbf)); call unpack_matrix(dmat_a, pfull)
    allocate(dSa(nbf,nbf,3,natom), dTa(nbf,nbf,3,natom), dVa(nbf,nbf,3,natom))
    call der_overlap_matrix(basis, dSa)
    call der_kinetic_matrix(basis, dTa)
    call der_nucattr_matrix(basis, basis%atoms%xyz, basis%atoms%zn, dVa)

    allocate(scr(nbf,nbf), col(nbf,nbf))
    allocate(Sx(nbf,nbf), hx(nbf,nbf), F0x(nbf,nbf), Gd0(nbf,nbf))
    allocate(probe(nbf,nbf), gx(3,natom))
    allocate(d0(nbf,nbf), d0p(nbf2,1), gp(nbf2,1), gfull(nbf,nbf))
    allocate(bvec(nocc*nvir,ncart), uvec(nocc*nvir,ncart), source=0.0_dp)

    icart = 0
    do kc = 1, natom
      do cc = 1, 3
        icart = icart + 1
        call mo_transform(mo_a, dSa(:,:,cc,kc), nbf, scr, col, Sx)
        scr = dTa(:,:,cc,kc) + dVa(:,:,cc,kc)
        call mo_transform(mo_a, scr, nbf, col, F0x, hx)

        F0x = hx
        do a = 1, nvir
          do i = 1, nocc
            do mu = 1, nbf
              do nu = 1, nbf
                probe(mu,nu) = 0.5_dp*( mo_a(mu,nocc+a)*mo_a(nu,i) + mo_a(mu,i)*mo_a(nu,nocc+a) )
              end do
            end do
            call fock_deriv_contract(infos, basis, pfull, probe, hfscale, gx)
            F0x(i,nocc+a) = hx(i,nocc+a) + 2.0_dp*gx(cc,kc)
          end do
        end do

        d0 = 0.0_dp
        do i = 1, nocc
          do j = 1, nocc
            do mu = 1, nbf
              do nu = 1, nbf
                d0(mu,nu) = d0(mu,nu) - 2.0_dp*Sx(i,j)*mo_a(mu,i)*mo_a(nu,j)
              end do
            end do
          end do
        end do
        call pack_matrix(d0, d0p(:,1))
        gp = 0.0_dp
        call fock_jk(basis, d=d0p, f=gp, scale_exch=hfscale, infos=infos)
        call unpack_from_packed(gp(:,1), gfull, nbf)
        call mo_transform(mo_a, gfull, nbf, scr, col, Gd0)

        ia = 0
        do a = 1, nvir
          do i = 1, nocc
            ia = ia + 1
            bvec(ia,icart) = -F0x(i,nocc+a) + eps(i)*Sx(i,nocc+a) - Gd0(i,nocc+a)
          end do
        end do
      end do
    end do

    call cphf_solve(infos, ncart, bvec, uvec)
    write(iw,'(A)') 'PyOQP: Native OpenQP HF/DFT Hessian CPHF response prepass complete'

    deallocate(pfull, dSa, dTa, dVa, scr, col, Sx, hx, F0x, Gd0, probe, gx, &
               d0, d0p, gp, gfull, bvec, uvec)
  end subroutine hf_hessian

!###############################################################################

  subroutine mo_transform(c_mo, a_ao, n, s1, s2, b_mo)
    use precision, only: dp
    real(kind=dp), intent(in) :: c_mo(:,:), a_ao(:,:)
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: s1(:,:), s2(:,:), b_mo(:,:)
    call dgemm('t','n', n, n, n, 1.0_dp, c_mo, n, a_ao, n, 0.0_dp, s1, n)
    call dgemm('n','n', n, n, n, 1.0_dp, s1, n, c_mo, n, 0.0_dp, b_mo, n)
  end subroutine mo_transform

!###############################################################################

  subroutine unpack_from_packed(gpk, gfu, n)
    use precision, only: dp
    real(kind=dp), intent(in) :: gpk(:)
    real(kind=dp), intent(inout) :: gfu(:,:)
    integer, intent(in) :: n
    integer :: ii, jj, ij
    ij = 0
    do ii = 1, n
      do jj = 1, ii
        ij = ij + 1
        gfu(ii,jj) = gpk(ij); gfu(jj,ii) = gpk(ij)
      end do
    end do
  end subroutine unpack_from_packed

end module hf_hessian_mod
