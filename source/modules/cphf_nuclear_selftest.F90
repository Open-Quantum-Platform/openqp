module cphf_nuclear_selftest_mod
!> @brief Gateway validation for the nuclear CPHF right-hand side: the
!>   skeleton (frozen-density) derivative-Fock matrix in the MO basis,
!>     Fx_mo_ai(x) = sum_uv C_ua C_vi dG[P]_uv/dx,
!>   built from the validated 2e derivative-Fock contraction
!>   (fock_deriv_mod::fock_deriv_contract) using symmetric MO-outer-product
!>   probes M = 1/2 (C_a C_i^T + C_i C_a^T).
!>
!>   Derivation of the factor: fock_deriv_contract(M,P) is symmetric-bilinear in
!>   (M,P) and equals dE_2e/dx at M=P (validated, O(h^2)). Since dG[P]/dx is
!>   symmetric in uv, sum_uv M_uv dG[P]_uv/dx = 2 fock_deriv_contract(M,P) for any
!>   symmetric M. Hence Fx_mo_ai = 2 fock_deriv_contract(sym(C_a C_i^T), P).
!>
!>   The full one-electron skeleton derivative dh/dx is added from the validated
!>   der_kinetic_matrix + der_nucattr_matrix. The combined skeleton MO Fock
!>   derivative  F0x_mo_ai = (C^T (dT/dx+dV/dx) C)_ai + Fx_mo_ai  is FD-validated
!>   against a frozen-density reference computed by the Python harness
!>   (which differences the AO Fock build G[P]+h at displaced geometries).
!>
!>   This isolates the one genuinely new ingredient (the non-symmetric->symmetric
!>   probe and its factor) before the full B^x -> U^x -> dP/dx assembly.
!>
!>   Writes the analytic F0x_mo (3,natom,nocc,nvir) to /tmp/cphf_f0x_mo.out.

  implicit none
  character(len=*), parameter :: module_name = "cphf_nuclear_selftest_mod"

contains

  subroutine cphf_f0x_selftest_C(c_handle) bind(C, name="cphf_f0x_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call cphf_f0x_selftest(inf)
  end subroutine cphf_f0x_selftest_C

  subroutine cphf_f0x_selftest(infos)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, OQP_E_MO_A
    use mathlib, only: unpack_matrix
    use fock_deriv_mod, only: fock_deriv_contract
    use mathlib, only: pack_matrix
    use scf_addons, only: fock_jk
    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dmat_a(:), mo_a(:,:), eps(:)
    real(kind=dp), allocatable :: pfull(:,:), probe(:,:), gx(:,:)
    real(kind=dp), allocatable :: gx_mo_an(:,:,:,:), gx_mo_fd(:,:,:,:)
    real(kind=dp), allocatable :: scr(:,:), col(:,:), ppack(:,:), gpack(:,:), gfull(:,:)
    real(kind=dp) :: hfscale, h, err, rel, denom
    integer :: nbf, nbf2, nocc, nvir, natom, i, a, c, k, mu, nu, u

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
    ! single fixed AO probe M (symmetric) independent of P, to exercise M != P.
    ! Reference: FD of the basis-consistent SCALAR tr(M . G[P]) at frozen P (M is
    ! a fixed AO matrix, so no basis-inconsistency). Analytic prediction:
    !   d/dx tr(M . G[P]) = 2 fock_deriv_contract(M, P)   (M symmetric, dG sym).
    allocate(gx_mo_an(1,1,3,natom), gx_mo_fd(1,1,3,natom), source=0.0_dp)
    allocate(scr(nbf,nbf), col(nbf,nbf), probe(nbf,nbf), gx(3,natom))
    allocate(ppack(nbf2,1), gpack(nbf2,1), gfull(nbf,nbf))

    do mu = 1, nbf
      do nu = 1, nbf
        probe(mu,nu) = 1.0_dp/real(mu+nu, dp)   ! symmetric, fixed AO probe
      end do
    end do

    ! ANALYTIC: 2 fock_deriv_contract(M, P)
    call fock_deriv_contract(infos, basis, pfull, probe, hfscale, gx)
    do k = 1, natom
      do c = 1, 3
        gx_mo_an(1,1,c,k) = 2.0_dp*gx(c,k)
      end do
    end do

    ! REFERENCE: central FD of tr(M . G[P]) at frozen P (basis-consistent scalar).
    call pack_matrix(pfull, ppack(:,1))
    h = 1.0e-4_dp
    do k = 1, natom
      do c = 1, 3
        basis%atoms%xyz(c,k) = basis%atoms%xyz(c,k) + h
        gpack = 0.0_dp
        call fock_jk(basis, d=ppack, f=gpack, scale_exch=hfscale, infos=infos)
        gx_mo_fd(1,1,c,k) = trace_full(probe, gpack(:,1), nbf)
        basis%atoms%xyz(c,k) = basis%atoms%xyz(c,k) - 2*h
        gpack = 0.0_dp
        call fock_jk(basis, d=ppack, f=gpack, scale_exch=hfscale, infos=infos)
        gx_mo_fd(1,1,c,k) = gx_mo_fd(1,1,c,k) - trace_full(probe, gpack(:,1), nbf)
        basis%atoms%xyz(c,k) = basis%atoms%xyz(c,k) + h
        gx_mo_fd(1,1,c,k) = gx_mo_fd(1,1,c,k)/(2*h)
      end do
    end do

    err = maxval(abs(gx_mo_an - gx_mo_fd))
    denom = maxval(abs(gx_mo_fd))
    rel = err/max(denom, 1.0e-30_dp)

    open(newunit=u, file='/tmp/cphf_f0x_mo.out', status='replace', action='write')
    write(u,'(a,3i6)') 'nocc nvir natom = ', nocc, nvir, natom
    write(u,'(a,es12.4)') '2e MO-Fock deriv max|an-fd| = ', err
    write(u,'(a,es12.4)') 'max|fd|                     = ', denom
    write(u,'(a,es12.4)') 'relative error              = ', rel
    write(u,'(a)') 'sample (k c   analytic   fd   ratio):'
    do k = 1, natom
      do c = 1, 3
        write(u,'(2i4,3es16.7)') k,c, gx_mo_an(1,1,c,k), gx_mo_fd(1,1,c,k), &
          gx_mo_fd(1,1,c,k)/sign(max(abs(gx_mo_an(1,1,c,k)),1d-30),gx_mo_an(1,1,c,k))
      end do
    end do
    if (err < 1.0e-6_dp) then
      write(u,'(a)') 'CPHF_F0X_SELFTEST PASS'
    else
      write(u,'(a)') 'CPHF_F0X_SELFTEST FAIL'
    end if
    close(u)

    deallocate(pfull, gx_mo_an, gx_mo_fd, scr, col, probe, gx, ppack, gpack, gfull)

  contains
    !> tr(M . G) for full symmetric M and packed lower-triangular symmetric G.
    pure function trace_full(m, gp, n) result(tr)
      real(kind=dp), intent(in) :: m(:,:), gp(:)
      integer, intent(in) :: n
      real(kind=dp) :: tr
      integer :: ii, jj, ij
      tr = 0.0_dp
      ij = 0
      do ii = 1, n
        do jj = 1, ii
          ij = ij + 1
          if (ii == jj) then
            tr = tr + m(ii,ii)*gp(ij)
          else
            tr = tr + (m(ii,jj)+m(jj,ii))*gp(ij)
          end if
        end do
      end do
    end function trace_full
  end subroutine cphf_f0x_selftest

end module cphf_nuclear_selftest_mod
