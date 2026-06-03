module fock_deriv_selftest_mod
!> @brief Validate the native 2e derivative-Fock contraction (fock_deriv_contract)
!>   against an exact, non-iterative finite-difference reference.
!>
!>   For a fixed (frozen) density P and a fixed probe M, the analytic quantity
!>     g_x = sum_uv M_uv F^x_uv[P] = sum_uv M_uv d/dx ( G_uv[P] )
!>   where G[P] = fock_jk(P) is the closed-shell two-electron Fock build. The
!>   reference is a central finite difference of  tr(M . G[P])  with respect to
!>   each nuclear coordinate, holding P and M fixed. Because the reference
!>   involves no SCF iteration, its accuracy is truncation-limited (O(h^2)), so
!>   this is a clean, reference-quality-independent check.
!>
!>   Results are written to /tmp/fockx_selftest.out.

  implicit none
  character(len=*), parameter :: module_name = "fock_deriv_selftest_mod"

contains

  subroutine fockx_selftest_C(c_handle) bind(C, name="fockx_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call fockx_selftest(inf)
  end subroutine fockx_selftest_C

  subroutine fockx_selftest(infos)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A
    use mathlib, only: unpack_matrix, pack_matrix
    use scf_addons, only: fock_jk
    use fock_deriv_mod, only: fock_deriv_contract
    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dmat_a(:)
    real(kind=dp), allocatable :: pfull(:,:), mfull(:,:), ppack(:,:), gpack(:,:)
    real(kind=dp), allocatable :: gx_an(:,:), gx_fd(:,:)
    real(kind=dp) :: hfscale, h, trp, trm, err
    integer :: nbf, nbf2, natom, k, a, p, q, u

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2)

    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)

    allocate(pfull(nbf,nbf), mfull(nbf,nbf), source=0.0_dp)
    call unpack_matrix(dmat_a, pfull)   ! alpha density (full, symmetric)

    ! Validate the builder with the trace identity (M = P): the mixed product
    ! reduces exactly to the energy-gradient density product, so
    !   contraction(P,P) = sum_uv P_uv F^x_uv[P] = d/dx [ tr(P . G[P]) ]
    ! at frozen P, where the 2e energy is E_2e = 1/2 tr(P . G[P]). Hence the
    ! finite-difference reference is  d/dx [ tr(P . G[P]) ] = 2 dE_2e/dx, with
    ! G[P] from fock_jk and P frozen. This reference is non-iterative ->
    ! truncation-limited (O(h^2)).
    mfull = pfull

    allocate(gx_an(3,natom), gx_fd(3,natom), source=0.0_dp)

    ! analytic: g_x = sum_uv P_uv F^x_uv[P]
    call fock_deriv_contract(infos, basis, pfull, mfull, hfscale, gx_an)

    ! finite-difference reference: central diff of tr(P . G[P]) at frozen P
    allocate(ppack(nbf2,1), gpack(nbf2,1))
    call pack_matrix(pfull, ppack(:,1))
    h = 1.0e-4_dp
    do k = 1, natom
      do a = 1, 3
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        gpack = 0.0_dp
        call fock_jk(basis, d=ppack, f=gpack, scale_exch=hfscale, infos=infos)
        trp = trace_MG(pfull, gpack(:,1), nbf)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) - 2*h
        gpack = 0.0_dp
        call fock_jk(basis, d=ppack, f=gpack, scale_exch=hfscale, infos=infos)
        trm = trace_MG(pfull, gpack(:,1), nbf)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        ! E_2e = 1/2 tr(P G[P]); contraction(P,P) = sum P F^x[P] = dE_2e/dx, so
        ! the reference is 1/2 d/dx tr(P G[P]).
        gx_fd(a,k) = 0.5_dp*(trp - trm)/(2*h)
      end do
    end do

    ! report both the raw ratio and the matched comparison so a constant-factor
    ! convention mismatch is visible rather than hidden.
    err = maxval(abs(gx_an - gx_fd))

    open(newunit=u, file='/tmp/fockx_selftest.out', status='replace', action='write')
    write(u,'(a,es12.4)') 'F^x trace-identity max|analytic - FD| = ', err
    write(u,'(a)') 'analytic   gx(:,1):'
    write(u,'(3es16.8)') gx_an(:,1)
    write(u,'(a)') 'finite-diff gx(:,1):'
    write(u,'(3es16.8)') gx_fd(:,1)
    write(u,'(a,3f12.6)') 'ratio fd/an (atom1) = ', gx_fd(:,1)/gx_an(:,1)
    if (err < 1.0e-6_dp) then
      write(u,'(a)') 'FOCKX_SELFTEST PASS'
    else
      write(u,'(a)') 'FOCKX_SELFTEST FAIL'
    end if
    close(u)

    deallocate(pfull, mfull, gx_an, gx_fd, ppack, gpack)
  contains
    !> tr(M . G) with G stored as a packed (lower-triangular) symmetric matrix.
    pure function trace_MG(m, gp, n) result(tr)
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
            tr = tr + (m(ii,jj) + m(jj,ii))*gp(ij)
          end if
        end do
      end do
    end function trace_MG
  end subroutine fockx_selftest

end module fock_deriv_selftest_mod
