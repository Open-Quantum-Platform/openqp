module cphf_dpdx_selftest_mod
!> @brief End-to-end native closed-shell (RHF) CPHF response chain
!>   F^x -> B^x -> U^x -> dpx/dx, computing the analytic relaxed alpha-density
!>   derivative for one nuclear coordinate (atom 1, z). A Python harness compares
!>   it to a central finite difference of the converged SCF alpha density at
!>   displaced geometries (the reference).
!>
!> @warning WORK IN PROGRESS, NOT YET VALIDATED. The assembled dP/dx currently
!>   disagrees with the finite-difference reference by a CONSTANT (h-independent)
!>   amount -- a flat plateau, i.e. a systematic formula error in this assembly,
!>   not finite-difference noise (the underlying ingredients F^x, the M!=P
!>   contraction, dS/dx, dT/dx, dV/dx, and cphf_solve are each independently
!>   validated). The error is non-uniform (~1.6-2.2x where the analytic is large),
!>   pointing to a convention bug in this end-to-end assembly: most likely the
!>   CPHF RHS sign/scale relative to cphf_solve's operator, the closed-shell
!>   density factor (alpha vs total), or the occ-occ overlap-response term. This
!>   harness is committed to record the chain and the open issue; it is NOT wired
!>   into any production path and hf_hessian stays guarded. Next: isolate by
!>   validating U_ai against a finite difference of the (phase-aligned) MO
!>   coefficients before the dC -> dP assembly.
!>
!>   Conventions (closed shell, OQP stores the alpha density P = sum_occ C C^T):
!>     S^x_pq  = (C^T dS/dx C)_pq
!>     h^x_pq  = (C^T (dT/dx+dV/dx) C)_pq
!>     G^x_pq  = sum_uv C_up C_vq dG[P]_uv/dx, from fock_deriv_contract:
!>               G^x_pq = 2 fock_deriv_contract(sym(C_p C_q^T), P)_x
!>     Skeleton derivative Fock (MO):  F0^x = h^x + G^x.
!>   Occupied-occupied orbital response is fixed by the overlap derivative:
!>     U_ij = -1/2 S^x_ij.
!>   The "B0" density from that fixed occ-occ rotation,
!>     D0_uv = sum_ij(occ) U_ij ( C_ui C_vj + C_uj C_vi ) ... (built directly),
!>   contributes a response Fock G[D0] that enters the occ-vir RHS. The CPHF
!>   equation for the occ-vir amplitudes U_ai (solved by cphf_solve, whose
!>   operator is (eps_a-eps_i) + response):
!>     (eps_a-eps_i) U_ai + sum_bj A_ai,bj U_bj
!>        = -F0^x_ai + eps_i S^x_ai - G[D0]^MO_ai
!>   i.e. cphf_solve is called with the right-hand side
!>     B_ai = -F0^x_ai + eps_i S^x_ai - G[D0]^MO_ai.
!>   Then
!>     dC_ui = sum_a C_u,a+nocc U_ai + sum_j C_uj U_ji,  U_ji = -1/2 S^x_ji
!>     dP_uv = sum_i ( dC_ui C_vi + C_ui dC_vi ).

  implicit none
  character(len=*), parameter :: module_name = "cphf_dpdx_selftest_mod"

contains

  subroutine cphf_dpdx_selftest_C(c_handle) bind(C, name="cphf_dpdx_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call cphf_dpdx_selftest(inf)
  end subroutine cphf_dpdx_selftest_C

  subroutine cphf_dpdx_selftest(infos)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, OQP_E_MO_A
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix
    use fock_deriv_mod, only: fock_deriv_contract
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve
    type(information), target, intent(inout) :: infos

    integer, parameter :: KC = 1, CC = 3   ! atom 1, z coordinate
    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dmat_a(:), mo_a(:,:), eps(:)
    real(kind=dp), allocatable :: pfull(:,:), probe(:,:), gx(:,:)
    real(kind=dp), allocatable :: dSa(:,:,:,:), dTa(:,:,:,:), dVa(:,:,:,:)
    real(kind=dp), allocatable :: Sx(:,:), hx(:,:), F0x(:,:), Gd0(:,:)
    real(kind=dp), allocatable :: d0(:,:), d0p(:,:), gp(:,:), gfull(:,:)
    real(kind=dp), allocatable :: bvec(:,:), uvec(:,:), scr(:,:), col(:,:)
    real(kind=dp), allocatable :: dcx(:,:), dpx(:,:)
    real(kind=dp) :: hfscale
    integer :: nbf, nbf2, nocc, nvir, natom, i, j, a, mu, nu, ia, u

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocc = infos%mol_prop%nocc
    nvir = nbf - nocc
    natom = size(basis%atoms%xyz, 2)
    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, eps)

    allocate(pfull(nbf,nbf)); call unpack_matrix(dmat_a, pfull)

    ! derivative integral matrices for coordinate (KC,CC)
    allocate(dSa(nbf,nbf,3,natom), dTa(nbf,nbf,3,natom), dVa(nbf,nbf,3,natom))
    call der_overlap_matrix(basis, dSa)
    call der_kinetic_matrix(basis, dTa)
    call der_nucattr_matrix(basis, basis%atoms%xyz, basis%atoms%zn, dVa)

    allocate(scr(nbf,nbf), col(nbf,nbf))
    allocate(Sx(nbf,nbf), hx(nbf,nbf), F0x(nbf,nbf), source=0.0_dp)

    call mo_transform(mo_a, dSa(:,:,CC,KC), nbf, scr, col, Sx)
    scr = dTa(:,:,CC,KC) + dVa(:,:,CC,KC)
    call mo_transform(mo_a, scr, nbf, col, F0x, hx)   ! hx = C^T(dT+dV)C

    ! skeleton 2e derivative G^x (MO occ-vir block)
    allocate(probe(nbf,nbf), gx(3,natom))
    F0x = hx
    do a = 1, nvir
      do i = 1, nocc
        do mu = 1, nbf
          do nu = 1, nbf
            probe(mu,nu) = 0.5_dp*( mo_a(mu,nocc+a)*mo_a(nu,i) + mo_a(mu,i)*mo_a(nu,nocc+a) )
          end do
        end do
        call fock_deriv_contract(infos, basis, pfull, probe, hfscale, gx)
        F0x(i,nocc+a) = hx(i,nocc+a) + 2.0_dp*gx(CC,KC)
      end do
    end do

    ! occ-occ overlap response density D0_uv = sum_ij U_ij C_ui C_vj, U_ij=-1/2 Sx_ij
    allocate(d0(nbf,nbf), source=0.0_dp)
    do i = 1, nocc
      do j = 1, nocc
        do mu = 1, nbf
          do nu = 1, nbf
            d0(mu,nu) = d0(mu,nu) - 0.5_dp*Sx(i,j)*mo_a(mu,i)*mo_a(nu,j)
          end do
        end do
      end do
    end do
    ! response Fock G[D0] (AO) then MO
    allocate(d0p(nbf2,1), gp(nbf2,1), gfull(nbf,nbf), Gd0(nbf,nbf), source=0.0_dp)
    call pack_matrix(d0, d0p(:,1))
    gp = 0.0_dp
    call fock_jk(basis, d=d0p, f=gp, scale_exch=hfscale, infos=infos)
    call unpack_from_packed(gp(:,1), gfull, nbf)
    call mo_transform(mo_a, gfull, nbf, scr, col, Gd0)

    ! CPHF RHS (occ-vir): B_ai = -F0x_ai + eps_i Sx_ai - Gd0_ai
    allocate(bvec(nocc*nvir,1), uvec(nocc*nvir,1), source=0.0_dp)
    ia = 0
    do a = 1, nvir
      do i = 1, nocc
        ia = ia + 1
        bvec(ia,1) = -F0x(i,nocc+a) + eps(i)*Sx(i,nocc+a) - Gd0(i,nocc+a)
      end do
    end do

    call cphf_solve(infos, 1, bvec, uvec)

    ! dC_ui = sum_a C_u,a+nocc U_ai - 1/2 sum_j C_uj Sx_ji
    allocate(dcx(nbf,nocc), source=0.0_dp)
    ia = 0
    do a = 1, nvir
      do i = 1, nocc
        ia = ia + 1
        do mu = 1, nbf
          dcx(mu,i) = dcx(mu,i) + mo_a(mu,nocc+a)*uvec(ia,1)
        end do
      end do
    end do
    do i = 1, nocc
      do j = 1, nocc
        do mu = 1, nbf
          dcx(mu,i) = dcx(mu,i) - 0.5_dp*mo_a(mu,j)*Sx(j,i)
        end do
      end do
    end do

    ! relaxed alpha-density derivative dP_uv = sum_i (dC_ui C_vi + C_ui dC_vi)
    allocate(dpx(nbf,nbf), source=0.0_dp)
    do i = 1, nocc
      do mu = 1, nbf
        do nu = 1, nbf
          dpx(mu,nu) = dpx(mu,nu) + dcx(mu,i)*mo_a(nu,i) + mo_a(mu,i)*dcx(nu,i)
        end do
      end do
    end do

    open(newunit=u, file='/tmp/cphf_dpdx_analytic.out', status='replace', action='write')
    write(u,'(i6)') nbf
    do i = 1, nbf
      do j = 1, nbf
        write(u,'(2i5,es20.10)') i, j, dpx(i,j)
      end do
    end do
    close(u)

    deallocate(pfull, dSa, dTa, dVa, scr, col, Sx, hx, F0x, probe, gx, &
               d0, d0p, gp, gfull, Gd0, bvec, uvec, dcx, dpx)
  end subroutine cphf_dpdx_selftest

  subroutine mo_transform(c_mo, a_ao, n, s1, s2, b_mo)
    use precision, only: dp
    real(kind=dp), intent(in) :: c_mo(:,:), a_ao(:,:)
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: s1(:,:), s2(:,:), b_mo(:,:)
    call dgemm('t','n', n, n, n, 1.0_dp, c_mo, n, a_ao, n, 0.0_dp, s1, n)
    call dgemm('n','n', n, n, n, 1.0_dp, s1, n, c_mo, n, 0.0_dp, b_mo, n)
  end subroutine mo_transform

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

end module cphf_dpdx_selftest_mod
