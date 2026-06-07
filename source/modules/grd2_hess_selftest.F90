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

!###############################################################################

  subroutine hess_skel_selftest_C(c_handle) bind(C, name="hess_skel_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call hess_skel_selftest(inf)
  end subroutine hess_skel_selftest_C

!###############################################################################

  subroutine hess_skel_selftest(infos)
    ! Compare the analytic fixed-density Hessian SKELETON (1e + 2e, excluding
    ! nuclear repulsion) against a central finite difference of the full
    ! frozen-density energy gradient (the same 1e + 2e pieces). With the
    ! Lagrangian (W) and density (P) matrices held fixed, FD of the gradient
    ! must reproduce the analytic second-derivative skeleton. This isolates the
    ! skeleton assembly (signs, W) from the CPHF response term.
    use types, only: information
    use precision, only: dp
    use grd1, only: eijden, grad_ee_overlap, grad_ee_kinetic, &
                    grad_en_hellman_feynman, grad_en_pulay, &
                    hess_ee_overlap, hess_ee_kinetic, hess_en
    use grd2, only: grd2_driver, grd2_hess_driver, grd2_compute_data_t
    use hf_gradient_mod, only: grd2_rhf_compute_data_t
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A

    implicit none

    type(information), target, intent(inout) :: infos

    class(grd2_compute_data_t), allocatable :: gcomp
    real(dp), contiguous, pointer :: dmat_a(:)
    real(dp), allocatable :: wlag(:), pden(:), hcc(:,:)
    real(dp), allocatable :: han(:,:), hfd(:,:), gp(:,:), gm(:,:), de2(:,:)
    real(dp), allocatable :: zneff(:)
    integer :: nbf, nbf2, natom, n3, k, a, ka, p
    real(dp) :: h, err

    associate(basis => infos%basis)

    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2)
    n3 = 3*natom

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)

    allocate(wlag(nbf2), pden(nbf2), hcc(n3,n3), source=0.0_dp)
    allocate(han(n3,n3), hfd(n3,n3), source=0.0_dp)
    allocate(gp(3,natom), gm(3,natom), de2(3,natom))
    allocate(zneff(natom))
    zneff = basis%atoms%zn - basis%ecp_zn_num

    ! frozen Lagrangian (W) and density (P)
    call eijden(wlag, nbf, infos)
    pden = dmat_a

    gcomp = grd2_rhf_compute_data_t( da = dmat_a, hfscale = 1.0_dp, nbf = nbf )
    call gcomp%init()

    ! analytic skeleton (1e + 2e, no nuclear repulsion)
    call hess_ee_overlap(basis, wlag, han)
    call hess_ee_kinetic(basis, pden, han)
    call hess_en(basis, basis%atoms%xyz, zneff, pden, han, hess_cc=hcc)
    call grd2_hess_driver(infos, basis, han, gcomp)

    ! finite difference of the frozen-density gradient (same pieces)
    h = 1.0e-4_dp
    do k = 1, natom
      do a = 1, 3
        ka = 3*(k-1) + a

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        call basis%init_shell_centers()
        gp = 0.0_dp
        call grad_ee_overlap(basis, wlag, gp)
        call grad_ee_kinetic(basis, pden, gp)
        call grad_en_hellman_feynman(basis, basis%atoms%xyz, zneff, pden, gp)
        call grad_en_pulay(basis, basis%atoms%xyz, zneff, pden, gp)
        de2 = 0.0_dp; call grd2_driver(infos, basis, de2, gcomp); gp = gp + de2

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) - 2*h
        call basis%init_shell_centers()
        gm = 0.0_dp
        call grad_ee_overlap(basis, wlag, gm)
        call grad_ee_kinetic(basis, pden, gm)
        call grad_en_hellman_feynman(basis, basis%atoms%xyz, zneff, pden, gm)
        call grad_en_pulay(basis, basis%atoms%xyz, zneff, pden, gm)
        de2 = 0.0_dp; call grd2_driver(infos, basis, de2, gcomp); gm = gm + de2

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        call basis%init_shell_centers()

        hfd(:,ka) = reshape((gp - gm)/(2*h), [n3])
      end do
    end do

    call gcomp%clean()

    err = maxval(abs(han - hfd)) / max(maxval(abs(hfd)), 1.0e-12_dp)

    block
      integer :: u
      open(newunit=u, file='/tmp/hess_skel_selftest.out', status='replace', action='write')
      write(u,'(a,es12.4)') 'skeleton max|an-fd|  = ', maxval(abs(han - hfd))
      write(u,'(a,es12.4)') 'skeleton rel error   = ', err
      write(u,'(a,es12.4)') 'max|analytic|        = ', maxval(abs(han))
      if (err < 1.0e-4_dp) then
        write(u,'(a)') 'HESS_SKEL_SELFTEST PASS'
      else
        write(u,'(a)') 'HESS_SKEL_SELFTEST FAIL'
      end if
      if (n3 <= 12) then
        write(u,'(a)') '--- analytic skeleton ---'
        do p = 1, n3; write(u,'(9es13.5)') han(p,1:n3); end do
        write(u,'(a)') '--- finite-difference skeleton ---'
        do p = 1, n3; write(u,'(9es13.5)') hfd(p,1:n3); end do
        write(u,'(a)') '--- analytic - FD ---'
        do p = 1, n3; write(u,'(9es13.5)') han(p,1:n3) - hfd(p,1:n3); end do
      end if
      close(u)
    end block

    deallocate(wlag, pden, hcc, han, hfd, gp, gm, de2, zneff)

    end associate
  end subroutine hess_skel_selftest

!###############################################################################

  subroutine hess_skel_open_selftest_C(c_handle) bind(C, name="hess_skel_open_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call hess_skel_open_selftest(inf)
  end subroutine hess_skel_open_selftest_C

!###############################################################################

  subroutine hess_skel_open_selftest(infos)
    ! Open-shell (UHF / ROHF) analog of hess_skel_selftest.  It compares the
    ! analytic fixed-density Hessian SKELETON (1e + 2e, excluding nuclear
    ! repulsion) against a central finite difference of the open-shell
    ! frozen-density energy gradient.
    !
    ! Open-shell specifics (vs the closed-shell hess_skel_selftest):
    !   * the total AO density is P = P_alpha + P_beta (OQP_DM_A + OQP_DM_B);
    !   * the energy-weighted (Lagrangian) density W comes from eijden, whose
    !     scftype>=2 branch already builds W = -(Fa.Pa + Fb.Pb) for U/ROHF;
    !   * the two-electron contraction uses grd2_uhf_compute_data_t (Coulomb
    !     from the total density, exchange from each spin density), exactly the
    !     object the open-shell production gradient (hf_2e_grad) uses.
    !
    ! With W and the spin densities held fixed, FD of the open-shell gradient
    ! reproduces the analytic skeleton; this validates the skeleton assembly
    ! (signs, W, total-vs-spin density routing) for both UHF and ROHF before the
    ! open-shell CPHF response term is added.  Works for scftype = 2 (UHF) and
    ! scftype = 3 (ROHF); aborts (cleanly) if invoked on a closed-shell run.
    use types, only: information
    use precision, only: dp
    use messages, only: show_message, WITH_ABORT
    use grd1, only: eijden, grad_ee_overlap, grad_ee_kinetic, &
                    grad_en_hellman_feynman, grad_en_pulay, &
                    hess_ee_overlap, hess_ee_kinetic, hess_en
    use grd2, only: grd2_driver, grd2_hess_driver, grd2_compute_data_t
    use hf_gradient_mod, only: grd2_uhf_compute_data_t
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_DM_B

    implicit none

    type(information), target, intent(inout) :: infos

    class(grd2_compute_data_t), allocatable :: gcomp
    real(dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    real(dp), allocatable :: wlag(:), pden(:), hcc(:,:)
    real(dp), allocatable :: han(:,:), hfd(:,:), gp(:,:), gm(:,:), de2(:,:)
    real(dp), allocatable :: zneff(:)
    integer :: nbf, nbf2, natom, n3, k, a, ka, p
    real(dp) :: h, err

    if (infos%control%scftype < 2) then
      call show_message('hess_skel_open_selftest requires an open-shell '// &
        '(UHF/ROHF) reference; scftype<2 detected.', WITH_ABORT)
    end if

    associate(basis => infos%basis)

    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2)
    n3 = 3*natom

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)

    allocate(wlag(nbf2), pden(nbf2), hcc(n3,n3), source=0.0_dp)
    allocate(han(n3,n3), hfd(n3,n3), source=0.0_dp)
    allocate(gp(3,natom), gm(3,natom), de2(3,natom))
    allocate(zneff(natom))
    zneff = basis%atoms%zn - basis%ecp_zn_num

    ! frozen Lagrangian (W, open-shell) and total density (P = Pa + Pb)
    call eijden(wlag, nbf, infos)
    pden = dmat_a + dmat_b

    ! open-shell 2e density-fetch object (same one the open-shell gradient uses)
    gcomp = grd2_uhf_compute_data_t( da = dmat_a, db = dmat_b, &
                                     hfscale = 1.0_dp, nbf = nbf )
    call gcomp%init()

    ! analytic skeleton (1e + 2e, no nuclear repulsion)
    call hess_ee_overlap(basis, wlag, han)
    call hess_ee_kinetic(basis, pden, han)
    call hess_en(basis, basis%atoms%xyz, zneff, pden, han, hess_cc=hcc)
    call grd2_hess_driver(infos, basis, han, gcomp)

    ! finite difference of the frozen-density open-shell gradient (same pieces)
    h = 1.0e-4_dp
    do k = 1, natom
      do a = 1, 3
        ka = 3*(k-1) + a

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        call basis%init_shell_centers()
        gp = 0.0_dp
        call grad_ee_overlap(basis, wlag, gp)
        call grad_ee_kinetic(basis, pden, gp)
        call grad_en_hellman_feynman(basis, basis%atoms%xyz, zneff, pden, gp)
        call grad_en_pulay(basis, basis%atoms%xyz, zneff, pden, gp)
        de2 = 0.0_dp; call grd2_driver(infos, basis, de2, gcomp); gp = gp + de2

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) - 2*h
        call basis%init_shell_centers()
        gm = 0.0_dp
        call grad_ee_overlap(basis, wlag, gm)
        call grad_ee_kinetic(basis, pden, gm)
        call grad_en_hellman_feynman(basis, basis%atoms%xyz, zneff, pden, gm)
        call grad_en_pulay(basis, basis%atoms%xyz, zneff, pden, gm)
        de2 = 0.0_dp; call grd2_driver(infos, basis, de2, gcomp); gm = gm + de2

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        call basis%init_shell_centers()

        hfd(:,ka) = reshape((gp - gm)/(2*h), [n3])
      end do
    end do

    call gcomp%clean()

    err = maxval(abs(han - hfd)) / max(maxval(abs(hfd)), 1.0e-12_dp)

    block
      integer :: u
      open(newunit=u, file='/tmp/hess_skel_open_selftest.out', status='replace', action='write')
      write(u,'(a,i0)')     'scftype              = ', infos%control%scftype
      write(u,'(a,es12.4)') 'skeleton max|an-fd|  = ', maxval(abs(han - hfd))
      write(u,'(a,es12.4)') 'skeleton rel error   = ', err
      write(u,'(a,es12.4)') 'max|analytic|        = ', maxval(abs(han))
      if (err < 1.0e-4_dp) then
        write(u,'(a)') 'HESS_SKEL_OPEN_SELFTEST PASS'
      else
        write(u,'(a)') 'HESS_SKEL_OPEN_SELFTEST FAIL'
      end if
      if (n3 <= 12) then
        write(u,'(a)') '--- analytic skeleton ---'
        do p = 1, n3; write(u,'(9es13.5)') han(p,1:n3); end do
        write(u,'(a)') '--- finite-difference skeleton ---'
        do p = 1, n3; write(u,'(9es13.5)') hfd(p,1:n3); end do
        write(u,'(a)') '--- analytic - FD ---'
        do p = 1, n3; write(u,'(9es13.5)') han(p,1:n3) - hfd(p,1:n3); end do
      end if
      close(u)
    end block

    deallocate(wlag, pden, hcc, han, hfd, gp, gm, de2, zneff)

    end associate
  end subroutine hess_skel_open_selftest

end module grd2_hess_selftest_mod
