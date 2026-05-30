module hess1_selftest_mod
!> @brief Diagnostic self-test for the native 1e (overlap + kinetic) Hessian
!>   contributions. It checks the analytic second-derivative drivers
!>   (grd1::hess_ee_overlap, grd1::hess_ee_kinetic) against a central finite
!>   difference of the production gradient routines (grd1::grad_ee_overlap,
!>   grd1::grad_ee_kinetic), contracted with the same fixed matrix. The identity
!>   d/dR [ sum M*dX/dR ] = sum M*d2X/dR2 holds for any fixed symmetric M, so a
!>   deterministic placeholder matrix is used and no SCF density is required.
!>
!>   This is a validation harness only; it does not compute or return a physical
!>   Hessian and is independent of the guarded hf_hessian kernel.

  implicit none

  character(len=*), parameter :: module_name = "hess1_selftest_mod"

contains

!###############################################################################

  subroutine hess1_selftest_C(c_handle) bind(C, name="hess1_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call hess1_selftest(inf)
  end subroutine hess1_selftest_C

!###############################################################################

  subroutine hess1_selftest(infos)
    use types, only: information
    use precision, only: dp
    use io_constants, only: iw
    use grd1, only: grad_ee_overlap, grad_ee_kinetic, &
                    hess_ee_overlap, hess_ee_kinetic, &
                    grad_en_hellman_feynman, grad_en_pulay, hess_en, &
                    der_overlap_matrix
    use basis_tools, only: bas_norm_matrix
    use mathlib, only: unpack_matrix

    implicit none

    type(information), target, intent(inout) :: infos

    integer :: nbf, nbf_tri, natom, n3, p, k, a, ka
    real(dp), allocatable :: m_packed(:)
    real(dp), allocatable :: hess_o_an(:,:), hess_k_an(:,:), hess_v_an(:,:)
    real(dp), allocatable :: hess_o_fd(:,:), hess_k_fd(:,:), hess_v_fd(:,:)
    real(dp), allocatable :: gp(:,:), gm(:,:), gp2(:,:), gm2(:,:)
    real(dp) :: h, err_o, err_k, err_v, sym_o, sym_k, sym_v, err_ds

    associate(basis => infos%basis)

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2)
    n3 = 3*natom

    ! deterministic, bounded fixed matrix (packed lower triangle)
    allocate(m_packed(nbf_tri))
    do p = 1, nbf_tri
      m_packed(p) = 1.0_dp / real(p+1, dp)
    end do

    allocate(hess_o_an(n3,n3), hess_k_an(n3,n3), hess_v_an(n3,n3), source=0.0_dp)
    allocate(hess_o_fd(n3,n3), hess_k_fd(n3,n3), hess_v_fd(n3,n3), source=0.0_dp)
    allocate(gp(3,natom), gm(3,natom), gp2(3,natom), gm2(3,natom))

    ! analytic second-derivative contributions
    call hess_ee_overlap(basis, m_packed, hess_o_an)
    call hess_ee_kinetic(basis, m_packed, hess_k_an)
    call hess_en(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, hess_v_an)

    ! central finite difference of the analytic gradient routines
    h = 1.0e-4_dp
    do k = 1, natom
      do a = 1, 3
        ka = 3*(k-1) + a

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        gp = 0.0_dp; call grad_ee_overlap(basis, m_packed, gp)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) - 2*h
        gm = 0.0_dp; call grad_ee_overlap(basis, m_packed, gm)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        hess_o_fd(:,ka) = reshape((gp-gm)/(2*h), [n3])

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        gp = 0.0_dp; call grad_ee_kinetic(basis, m_packed, gp)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) - 2*h
        gm = 0.0_dp; call grad_ee_kinetic(basis, m_packed, gm)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        hess_k_fd(:,ka) = reshape((gp-gm)/(2*h), [n3])

        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        gp = 0.0_dp; gp2 = 0.0_dp
        call grad_en_hellman_feynman(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, gp)
        call grad_en_pulay(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, gp2)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) - 2*h
        gm = 0.0_dp; gm2 = 0.0_dp
        call grad_en_hellman_feynman(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, gm)
        call grad_en_pulay(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, gm2)
        basis%atoms%xyz(a,k) = basis%atoms%xyz(a,k) + h
        hess_v_fd(:,ka) = reshape(((gp+gp2)-(gm+gm2))/(2*h), [n3])
      end do
    end do

    err_o = maxval(abs(hess_o_an - hess_o_fd))
    err_k = maxval(abs(hess_k_an - hess_k_fd))
    err_v = maxval(abs(hess_v_an - hess_v_fd))
    sym_o = maxval(abs(hess_o_an - transpose(hess_o_an)))
    sym_k = maxval(abs(hess_k_an - transpose(hess_k_an)))
    sym_v = maxval(abs(hess_v_an - transpose(hess_v_an)))

    ! --- CPHF RHS building block: dS/dR matrix vs the validated overlap gradient.
    !   sum_uv (bfnrm_u bfnrm_v M_uv) dS(u,v,c,A) must equal grad_ee_overlap(M).
    block
      real(dp), allocatable :: dSmat(:,:,:,:), mnorm(:,:), g_an(:,:), g_ref(:,:)
      integer :: cc, kk, mu, nu
      allocate(dSmat(nbf,nbf,3,natom), mnorm(nbf,nbf), g_an(3,natom), g_ref(3,natom))
      call der_overlap_matrix(basis, dSmat)
      call unpack_matrix(m_packed, mnorm)
      call bas_norm_matrix(mnorm, basis%bfnrm, nbf)
      g_an = 0.0_dp
      do kk = 1, natom
        do cc = 1, 3
          do mu = 1, nbf
            do nu = 1, nbf
              g_an(cc,kk) = g_an(cc,kk) + mnorm(mu,nu)*dSmat(mu,nu,cc,kk)
            end do
          end do
        end do
      end do
      g_ref = 0.0_dp
      call grad_ee_overlap(basis, m_packed, g_ref)
      err_ds = maxval(abs(g_an - g_ref))
      deallocate(dSmat, mnorm, g_an, g_ref)
    end block

    ! Separated finite-difference diagnostic (small systems only): FD the basis
    ! (Pulay) and charge (Hellmann-Feynman) gradient pieces against the basis and
    ! charge positions INDEPENDENTLY, since grad_en_* take coord separately from
    ! basis%atoms%xyz. This isolates d2E/dbasis2, d2E/dbasis dcharge, d2E/dcharge2.
    if (n3 <= 12) then
      block
        real(dp) :: fd_bb(n3,n3), fd_bc(n3,n3), fd_cb(n3,n3), fd_cc(n3,n3)
        real(dp), allocatable :: coord0(:,:), cmov(:,:)
        integer :: kk, aa, kak
        allocate(coord0(3,natom), cmov(3,natom))
        coord0 = basis%atoms%xyz   ! fixed reference positions
        fd_bb = 0; fd_bc = 0; fd_cb = 0; fd_cc = 0
        do kk = 1, natom
          do aa = 1, 3
            kak = 3*(kk-1) + aa
            ! ---- move BASIS, hold charge=coord0 fixed ----
            basis%atoms%xyz(aa,kk) = coord0(aa,kk) + h
            gp = 0; gp2 = 0
            call grad_en_pulay(basis, coord0, basis%atoms%zn, m_packed, gp2)
            call grad_en_hellman_feynman(basis, coord0, basis%atoms%zn, m_packed, gp)
            basis%atoms%xyz(aa,kk) = coord0(aa,kk) - h
            gm = 0; gm2 = 0
            call grad_en_pulay(basis, coord0, basis%atoms%zn, m_packed, gm2)
            call grad_en_hellman_feynman(basis, coord0, basis%atoms%zn, m_packed, gm)
            basis%atoms%xyz(aa,kk) = coord0(aa,kk)
            fd_bb(:,kak) = reshape((gp2-gm2)/(2*h), [n3])   ! d(Pulay)/d(basis)
            fd_cb(:,kak) = reshape((gp -gm )/(2*h), [n3])   ! d(HF)/d(basis)
            ! ---- move CHARGE, hold basis fixed ----
            cmov = coord0; cmov(aa,kk) = coord0(aa,kk) + h
            gp = 0; gp2 = 0
            call grad_en_pulay(basis, cmov, basis%atoms%zn, m_packed, gp2)
            call grad_en_hellman_feynman(basis, cmov, basis%atoms%zn, m_packed, gp)
            cmov(aa,kk) = coord0(aa,kk) - h
            gm = 0; gm2 = 0
            call grad_en_pulay(basis, cmov, basis%atoms%zn, m_packed, gm2)
            call grad_en_hellman_feynman(basis, cmov, basis%atoms%zn, m_packed, gm)
            fd_bc(:,kak) = reshape((gp2-gm2)/(2*h), [n3])   ! d(Pulay)/d(charge)
            fd_cc(:,kak) = reshape((gp -gm )/(2*h), [n3])   ! d(HF)/d(charge)
          end do
        end do
        block
          integer :: u2
          open(newunit=u2, file='/tmp/hess1_sep.out', status='replace', action='write')
          write(u2,'(a,es12.4)') 'sum(fd_bb+fd_bc+fd_cb+fd_cc) vs total fd = ', &
            maxval(abs((fd_bb+fd_bc+fd_cb+fd_cc) - hess_v_fd))
          write(u2,'(a)') '--- fd_bb (d2E/dbasis dbasis) ---'
          do p=1,n3; write(u2,'(12es12.4)') fd_bb(p,1:n3); end do
          write(u2,'(a)') '--- fd_cb (d2E/dcharge dbasis) ---'
          do p=1,n3; write(u2,'(12es12.4)') fd_cb(p,1:n3); end do
          write(u2,'(a)') '--- fd_cc (d2E/dcharge dcharge) ---'
          do p=1,n3; write(u2,'(12es12.4)') fd_cc(p,1:n3); end do
          write(u2,'(a)') '--- hess_en analytic (total) ---'
          do p=1,n3; write(u2,'(12es12.4)') hess_v_an(p,1:n3); end do
          close(u2)
        end block
        deallocate(coord0, cmov)
      end block
    end if

    block
      integer :: u
      open(newunit=u, file='/tmp/hess1_selftest.out', status='replace', action='write')
      write(u,'(a,es12.4)') 'overlap max|an-fd| = ', err_o
      write(u,'(a,es12.4)') 'kinetic max|an-fd| = ', err_k
      write(u,'(a,es12.4)') 'overlap asymmetry  = ', sym_o
      write(u,'(a,es12.4)') 'kinetic asymmetry  = ', sym_k
      ! Nuclear-attraction (hess_en) is WIP: reported but NOT part of PASS/FAIL.
      write(u,'(a,es12.4)') 'nucattr max|an-fd| (WIP) = ', err_v
      write(u,'(a,es12.4)') 'nucattr asymmetry  (WIP) = ', sym_v
      write(u,'(a,es12.4)') 'dS/dR matrix vs grad     = ', err_ds
      if (max(err_o, err_k) < 1.0e-6_dp) then
        write(u,'(a)') 'HESS1E_SELFTEST PASS'
      else
        write(u,'(a)') 'HESS1E_SELFTEST FAIL'
      end if
      close(u)
    end block

    deallocate(m_packed, hess_o_an, hess_k_an, hess_v_an, &
               hess_o_fd, hess_k_fd, hess_v_fd, gp, gm, gp2, gm2)
    end associate
  end subroutine hess1_selftest

end module hess1_selftest_mod
