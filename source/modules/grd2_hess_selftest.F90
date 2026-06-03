module grd2_hess_selftest_mod
!> @brief Diagnostic self-test for the analytic two-electron (ERI) second
!>   derivative skeleton. It checks the analytic per-quartet Hessian driver
!>   grd2::grd2_hess_driver against a central finite difference of the
!>   production analytic 2e gradient grd2::grd2_driver, contracted with the
!>   same fixed converged density. The identity d/dR [ sum P P dX/dR ] =
!>   sum P P d2X/dR2 holds for a fixed density matrix P, so no SCF response is
!>   involved; this isolates the ERI second-derivative integral contraction.
!>
!>   Validation harness only: it does not compute or return a physical Hessian
!>   and is independent of the hf_hessian kernel.

  implicit none

  character(len=*), parameter :: module_name = "grd2_hess_selftest_mod"

contains

!###############################################################################

  subroutine grd2_hess_selftest_C(c_handle) bind(C, name="grd2_hess_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call grd2_hess_selftest(inf)
  end subroutine grd2_hess_selftest_C

!###############################################################################

  subroutine grd2_hess_selftest(infos)
    use types, only: information
    use precision, only: dp
    use io_constants, only: iw
    use grd2, only: grd2_driver, grd2_hess_driver, grd2_compute_data_t
    use hf_gradient_mod, only: grd2_rhf_compute_data_t
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A

    implicit none

    type(information), target, intent(inout) :: infos

    class(grd2_compute_data_t), allocatable :: gcomp
    real(dp), contiguous, pointer :: dmat_a(:)
    real(dp), allocatable :: hess_an(:,:), hess_fd(:,:), de_p(:,:), de_m(:,:)
    integer :: natom, n3, k, a, ka, p
    real(dp) :: h, err, sym, errel

    associate(basis => infos%basis)

    basis%atoms => infos%atoms
    natom = size(basis%atoms%xyz, 2)
    n3 = 3*natom

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)

    ! Fixed closed-shell RHF density-fetch object (same one the gradient uses)
    gcomp = grd2_rhf_compute_data_t( da = dmat_a, hfscale = 1.0_dp, nbf = basis%nbf )
    call gcomp%init()

    allocate(hess_an(n3,n3), hess_fd(n3,n3), source=0.0_dp)
    allocate(de_p(3,natom), de_m(3,natom))

    ! Analytic 2e second-derivative skeleton
    call grd2_hess_driver(infos, basis, hess_an, gcomp)

    ! Central finite difference of the analytic 2e gradient (density fixed)
    h = 1.0e-4_dp
    do k = 1, natom
      do a = 1, 3
        ka = 3*(k-1) + a

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        call basis%init_shell_centers()
        de_p = 0.0_dp; call grd2_driver(infos, basis, de_p, gcomp)

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) - 2*h
        call basis%init_shell_centers()
        de_m = 0.0_dp; call grd2_driver(infos, basis, de_m, gcomp)

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        call basis%init_shell_centers()

        hess_fd(:,ka) = reshape((de_p - de_m)/(2*h), [n3])
      end do
    end do

    call gcomp%clean()

    err = maxval(abs(hess_an - hess_fd))
    sym = maxval(abs(hess_an - transpose(hess_an)))
    errel = err / max(maxval(abs(hess_fd)), 1.0e-12_dp)

    block
      integer :: u
      open(newunit=u, file='/tmp/grd2_hess_selftest.out', status='replace', action='write')
      write(u,'(a,es12.4)') '2e hess max|an-fd|   = ', err
      write(u,'(a,es12.4)') '2e hess rel error    = ', errel
      write(u,'(a,es12.4)') '2e hess asymmetry    = ', sym
      write(u,'(a,es12.4)') 'max|analytic|        = ', maxval(abs(hess_an))
      write(u,'(a,es12.4)') 'max|finite-diff|     = ', maxval(abs(hess_fd))
      if (errel < 1.0e-4_dp) then
        write(u,'(a)') 'GRD2_HESS_SELFTEST PASS'
      else
        write(u,'(a)') 'GRD2_HESS_SELFTEST FAIL'
      end if
      if (n3 <= 12) then
        write(u,'(a)') '--- analytic ---'
        do p = 1, n3; write(u,'(9es13.5)') hess_an(p,1:n3); end do
        write(u,'(a)') '--- finite difference ---'
        do p = 1, n3; write(u,'(9es13.5)') hess_fd(p,1:n3); end do
        write(u,'(a)') '--- analytic - finite difference ---'
        do p = 1, n3; write(u,'(9es13.5)') hess_an(p,1:n3) - hess_fd(p,1:n3); end do
      end if
      close(u)
    end block

    deallocate(hess_an, hess_fd, de_p, de_m)

    end associate
  end subroutine grd2_hess_selftest

end module grd2_hess_selftest_mod
