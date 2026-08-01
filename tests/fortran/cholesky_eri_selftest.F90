module cholesky_eri_selftest_mod
!> @brief Validate the Cholesky factorisation of the two-electron integrals
!>   against the integrals themselves.
!>
!>   The decomposition claims
!>
!>     (mu nu | lambda sigma) = sum_J L(P,J) L(Q,J)
!>
!>   to within the requested tolerance.  The only check worth making is the
!>   direct one: build the AO integrals with the shared engine, factorise them,
!>   then reconstruct every element and compare.  A rank or a residual norm
!>   would not catch a wrong pivot or a mis-indexed column; a reconstruction
!>   does.
!>
!>   Results are written to /tmp/cholesky_eri_selftest.out.

  implicit none
  character(len=*), parameter :: module_name = "cholesky_eri_selftest_mod"

contains

  subroutine cholesky_eri_selftest_C(c_handle) bind(C, name="cholesky_eri_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call cholesky_eri_selftest(inf)
  end subroutine cholesky_eri_selftest_C

  subroutine cholesky_eri_selftest(infos)

    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use cc_ao2mo, only: cc_eri_collect_t, cc_packed_length
    use cholesky_eri, only: cholesky_eri_decompose, cholesky_eri_max_vectors, &
                            cholesky_packed_index
    use cholesky_direct, only: cholesky_direct_decompose

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(cc_eri_collect_t), target :: eri_data

    real(dp), allocatable, target :: gao(:)
    real(dp), allocatable :: lvec(:,:), dvec(:,:)
    integer :: nbf, npair, nchol, maxchol, p, q, iu
    integer :: dchol, npass
    real(dp) :: derr, dmax_d, dapprox
    logical :: dtrunc
    real(dp) :: tol, err, exact, approx, dmax, drms
    integer(8) :: nel
    logical :: truncated

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    npair = nbf*(nbf+1)/2

    allocate(gao(cc_packed_length(nbf)))
    gao = 0.0_dp

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    eri_data%g => gao
    eri_data%nbf = nbf
    eri_data%npair = nbf*(nbf+1)/2
    call int2_driver%run(eri_data)
    call eri_data%clean()
    call int2_driver%clean()

    tol = 1.0e-8_dp
    maxchol = cholesky_eri_max_vectors(nbf, 20)
    allocate(lvec(npair, maxchol))

    call cholesky_eri_decompose(nbf, gao, tol, maxchol, lvec, nchol, err, &
                                truncated)

    ! Reconstruct and compare against the integrals we started from.
    dmax = 0.0_dp
    drms = 0.0_dp
    nel = 0
    do p = 1, npair
      do q = 1, p
        exact = gao(cholesky_packed_index(p, q))
        approx = dot_product(lvec(p,1:nchol), lvec(q,1:nchol))
        dmax = max(dmax, abs(exact-approx))
        drms = drms + (exact-approx)**2
        nel = nel + 1
      end do
    end do
    drms = sqrt(drms/real(nel, dp))

    ! Integral-direct factorisation of the same integrals -- it never builds
    ! the packed store, and is checked against it all the same.
    allocate(dvec(npair, maxchol))
    call cholesky_direct_decompose(basis, infos, tol, maxchol, dvec, dchol, &
                                   derr, dtrunc, npass)
    dmax_d = 0.0_dp
    do p = 1, npair
      do q = 1, p
        exact = gao(cholesky_packed_index(p, q))
        dapprox = dot_product(dvec(p,1:dchol), dvec(q,1:dchol))
        dmax_d = max(dmax_d, abs(exact-dapprox))
      end do
    end do

    open(newunit=iu, file='/tmp/cholesky_eri_selftest.out', status='replace')
    write(iu,'(A,I0)')    'nbf          = ', nbf
    write(iu,'(A,I0)')    'npair        = ', npair
    write(iu,'(A,I0)')    'nchol        = ', nchol
    write(iu,'(A,L1)')    'truncated    = ', truncated
    write(iu,'(A,ES14.6)')'tol          = ', tol
    write(iu,'(A,ES14.6)')'residual_max = ', err
    write(iu,'(A,ES14.6)')'recon_max    = ', dmax
    write(iu,'(A,ES14.6)')'recon_rms    = ', drms
    write(iu,'(A,F10.4)') 'nchol_per_nbf= ', real(nchol,dp)/real(nbf,dp)
    write(iu,'(A,I0)')    'direct_nchol = ', dchol
    write(iu,'(A,I0)')    'direct_npass = ', npass
    write(iu,'(A,L1)')    'direct_trunc = ', dtrunc
    write(iu,'(A,ES14.6)')'direct_recon_max = ', dmax_d
    close(iu)

    deallocate(gao, lvec, dvec)

  end subroutine cholesky_eri_selftest

end module cholesky_eri_selftest_mod
