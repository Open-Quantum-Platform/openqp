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
                    der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix
    use basis_tools, only: bas_norm_matrix
    use mathlib, only: unpack_matrix, orthogonal_transform_sym
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_SM, OQP_VEC_MO_A, OQP_E_MO_A

    implicit none

    type(information), target, intent(inout) :: infos

    integer :: nbf, nbf_tri, natom, n3, p, k, a, ka
    real(dp), allocatable :: m_packed(:)
    real(dp), allocatable :: hess_o_an(:,:), hess_k_an(:,:), hess_v_an(:,:), hess_v_an_b9(:,:)
    real(dp), allocatable :: hess_o_fd(:,:), hess_k_fd(:,:), hess_v_fd(:,:)
    real(dp), allocatable :: gp(:,:), gm(:,:), gp2(:,:), gm2(:,:)
    real(dp) :: h, err_o, err_k, err_v, sym_o, sym_k, sym_v, err_ds, err_dt, err_dv, err_ortho

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

    allocate(hess_o_an(n3,n3), hess_k_an(n3,n3), hess_v_an(n3,n3), hess_v_an_b9(n3,n3), source=0.0_dp)
    allocate(hess_o_fd(n3,n3), hess_k_fd(n3,n3), hess_v_fd(n3,n3), source=0.0_dp)
    allocate(gp(3,natom), gm(3,natom), gp2(3,natom), gm2(3,natom))

    ! analytic second-derivative contributions
    call hess_ee_overlap(basis, m_packed, hess_o_an)
    call hess_ee_kinetic(basis, m_packed, hess_k_an)
    call hess_en(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, hess_v_an, &
                 hess_cc=hess_v_an_b9)

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

      ! same check for the kinetic derivative matrix vs grad_ee_kinetic
      call der_kinetic_matrix(basis, dSmat)
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
      call grad_ee_kinetic(basis, m_packed, g_ref)
      err_dt = maxval(abs(g_an - g_ref))

      ! nuclear-attraction derivative matrix vs grad_en_pulay + grad_en_HF
      call der_nucattr_matrix(basis, basis%atoms%xyz, basis%atoms%zn, dSmat)
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
      call grad_en_pulay(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, g_ref)
      call grad_en_hellman_feynman(basis, basis%atoms%xyz, basis%atoms%zn, m_packed, g_ref)
      err_dv = maxval(abs(g_an - g_ref))
      deallocate(dSmat, mnorm, g_an, g_ref)
    end block

    ! CPHF prerequisite: MO data access + transform path. Verify C^T S C = I
    ! (MO orthonormality) using the production overlap and MO coefficients.
    err_ortho = -1.0_dp
    block
      real(dp), contiguous, pointer :: smat(:), mo_a(:,:)
      real(dp), allocatable :: smo(:), ident(:)
      integer :: p2, q2, ij
      call tagarray_get_data(infos%dat, OQP_SM, smat)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
      allocate(smo(nbf*(nbf+1)/2), ident(nbf*(nbf+1)/2))
      call orthogonal_transform_sym(nbf, nbf, smat, mo_a, nbf, smo)
      ident = 0.0_dp
      ij = 0
      do p2 = 1, nbf
        do q2 = 1, p2
          ij = ij + 1
          if (p2 == q2) ident(ij) = 1.0_dp
        end do
      end do
      err_ortho = maxval(abs(smo - ident))
      deallocate(smo, ident)
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
          write(u2,'(a)') '--- hess_en FD (total) ---'
          do p=1,n3; write(u2,'(12es12.4)') hess_v_fd(p,1:n3); end do
          write(u2,'(a)') '--- analytic - FD ---'
          do p=1,n3; write(u2,'(12es12.4)') hess_v_an(p,1:n3)-hess_v_fd(p,1:n3); end do
          write(u2,'(a)') '--- analytic block9 (p_CC only) ---'
          do p=1,n3; write(u2,'(12es12.4)') hess_v_an_b9(p,1:n3); end do
          write(u2,'(a)') '--- block9 - fd_cc ---'
          do p=1,n3; write(u2,'(12es12.4)') hess_v_an_b9(p,1:n3)-fd_cc(p,1:n3); end do
          close(u2)
        end block
        deallocate(coord0, cmov)
      end block
    end if

    ! ---- Gate 2 PRIMARY oracle data ------------------------------------------
    ! Emit per-nucleus uncontracted basis-basis second-derivative blocks
    !   PAA(a,b,mu,nu) = d2/dA_a dA_b <mu| 1/|r-C| |nu>   (bra-bra)
    !   PAB(a,b,mu,nu) = d2/dA_a dB_b <mu| 1/|r-C| |nu>   (bra-ket mixed)
    ! for each nucleus C, chargeless (znuc=1) and bfnrm-normalized to match the
    ! stored overlap OQP_SM. The Python oracle (tests/test_hess_nuc_oracle.py)
    ! applies -Z_C and the explicit AO permutation/scale map, then compares
    ! element-wise to PySCF int1e_ipiprinv / int1e_iprinvip evaluated with
    ! with_rinv_at_nucleus(C). This is the primary correctness gate for the
    ! native Rys nuclear-attraction Hessian integrals (AM-shift formulation).
    if (nbf <= 64) then
      block
        use mod_1e_primitives, only: comp_coulomb_der2_blocks
        use mod_shell_tools, only: shell_t, shpair_t
        use constants, only: CART_X, CART_Y, CART_Z
        type(shell_t) :: eshi, eshj
        type(shpair_t) :: ecab
        real(dp), allocatable :: PAA(:,:,:,:), PAB(:,:,:,:), blkAA(:,:,:,:), blkAB(:,:,:,:)
        real(dp), allocatable :: Sfull(:,:)
        real(dp), contiguous, pointer :: smat2(:)
        integer, allocatable :: ao_atom(:), ao_lx(:), ao_ly(:), ao_lz(:)
        integer :: ic2, ii2, jj2, i2, j2, o_i, o_j, a2, b2, u3

        allocate(ao_atom(nbf), ao_lx(nbf), ao_ly(nbf), ao_lz(nbf))
        do ii2 = 1, basis%nshell
          call eshi%fetch_by_id(basis, ii2)
          do i2 = 1, eshi%nao
            o_i = basis%ao_offset(ii2) + i2 - 1
            ao_atom(o_i) = eshi%atid - 1            ! 0-based to match PySCF
            ao_lx(o_i) = CART_X(i2, eshi%ang)
            ao_ly(o_i) = CART_Y(i2, eshi%ang)
            ao_lz(o_i) = CART_Z(i2, eshi%ang)
          end do
        end do

        allocate(Sfull(nbf,nbf))
        call tagarray_get_data(infos%dat, OQP_SM, smat2)
        call unpack_matrix(smat2, Sfull)

        allocate(PAA(3,3,nbf,nbf), PAB(3,3,nbf,nbf))
        call ecab%alloc(basis)
        open(newunit=u3, file='/tmp/hess_nuc_blocks.txt', status='replace', action='write')
        write(u3,'(2i6)') nbf, natom
        do o_i = 1, nbf
          write(u3,'(4i5)') ao_atom(o_i), ao_lx(o_i), ao_ly(o_i), ao_lz(o_i)
        end do
        ! atoms: nuclear charge (integer) and coordinates in bohr
        do ic2 = 1, natom
          write(u3,'(i5,3es24.16)') nint(basis%atoms%zn(ic2)), basis%atoms%xyz(1:3,ic2)
        end do
        do i2 = 1, nbf
          write(u3,'(*(es24.16))') Sfull(i2,1:nbf)
        end do
        do ic2 = 1, natom
          PAA = 0.0_dp; PAB = 0.0_dp
          do ii2 = 1, basis%nshell
            call eshi%fetch_by_id(basis, ii2)
            do jj2 = 1, basis%nshell
              call eshj%fetch_by_id(basis, jj2)
              call ecab%shell_pair(basis, eshi, eshj)
              if (ecab%numpairs == 0) cycle
              allocate(blkAA(3,3,ecab%inao,ecab%jnao), blkAB(3,3,ecab%inao,ecab%jnao))
              call comp_coulomb_der2_blocks(ecab, basis%atoms%xyz(:,ic2), 1.0_dp, blkAA, blkAB)
              do i2 = 1, ecab%inao
                o_i = basis%ao_offset(ii2) + i2 - 1
                do j2 = 1, ecab%jnao
                  o_j = basis%ao_offset(jj2) + j2 - 1
                  PAA(:,:,o_i,o_j) = blkAA(:,:,i2,j2)*basis%bfnrm(o_i)*basis%bfnrm(o_j)
                  PAB(:,:,o_i,o_j) = blkAB(:,:,i2,j2)*basis%bfnrm(o_i)*basis%bfnrm(o_j)
                end do
              end do
              deallocate(blkAA, blkAB)
            end do
          end do
          do a2 = 1, 3
            do b2 = 1, 3
              do i2 = 1, nbf
                write(u3,'(*(es24.16))') PAA(a2,b2,i2,1:nbf)
              end do
            end do
          end do
          do a2 = 1, 3
            do b2 = 1, 3
              do i2 = 1, nbf
                write(u3,'(*(es24.16))') PAB(a2,b2,i2,1:nbf)
              end do
            end do
          end do
        end do
        close(u3)
        deallocate(PAA, PAB, Sfull, ao_atom, ao_lx, ao_ly, ao_lz)
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
      write(u,'(a,es12.4)') 'dT/dR matrix vs grad     = ', err_dt
      write(u,'(a,es12.4)') 'dV/dR matrix vs grad     = ', err_dv
      write(u,'(a,es12.4)') 'MO orthonormality C^TSC-I = ', err_ortho
      if (max(err_o, err_k) < 1.0e-6_dp) then
        write(u,'(a)') 'HESS1E_SELFTEST PASS'
      else
        write(u,'(a)') 'HESS1E_SELFTEST FAIL'
      end if
      close(u)
    end block

    deallocate(m_packed, hess_o_an, hess_k_an, hess_v_an, hess_v_an_b9, &
               hess_o_fd, hess_k_fd, hess_v_fd, gp, gm, gp2, gm2)
    end associate
  end subroutine hess1_selftest

end module hess1_selftest_mod
