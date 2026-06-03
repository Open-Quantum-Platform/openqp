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
    ! calls cphf_solve on the full 3N RHS block and stores the native Hessian
    ! matrix in OQP::hf_hessian for the Python frequency driver.
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, OQP_E_MO_A, &
      OQP_hf_hessian, TA_TYPE_REAL64
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix, hess_nn
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
    real(kind=dp), allocatable :: bvec(:,:), uvec(:,:), scr(:,:), col(:,:), hess_native(:,:)
    real(kind=dp), contiguous, pointer :: hess_store(:,:)
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
    write(iw,'(A)') '  Storing native OpenQP HF/DFT analytic Hessian matrix in OQP::hf_hessian.'

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

    ! der_* matrices are returned in the UNNORMALIZED basis; bring them into the
    ! same normalized (bfnrm) convention as the MO coefficients / density so the
    ! CPHF RHS and the response contractions are correct for d/f functions
    ! (bfnrm /= 1). Invisible for s/p-only bases (e.g. STO-3G).
    block
      integer :: kc2, cc2, mu2, nu2
      do kc2 = 1, natom
        do cc2 = 1, 3
          do nu2 = 1, nbf
            do mu2 = 1, nbf
              dSa(mu2,nu2,cc2,kc2) = dSa(mu2,nu2,cc2,kc2)*basis%bfnrm(mu2)*basis%bfnrm(nu2)
              dTa(mu2,nu2,cc2,kc2) = dTa(mu2,nu2,cc2,kc2)*basis%bfnrm(mu2)*basis%bfnrm(nu2)
              dVa(mu2,nu2,cc2,kc2) = dVa(mu2,nu2,cc2,kc2)*basis%bfnrm(mu2)*basis%bfnrm(nu2)
            end do
          end do
        end do
      end do
    end block

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

        ! --- XC contribution to the CPKS right-hand side (DFT only) -------------
        ! The A-matrix (cphf_apbx) includes the XC kernel fxc, so the perturbation
        ! RHS must carry BOTH XC pieces or the relaxed response dPx is wrong (the
        ! HF response, ~2x too large for DFT):
        !   (i)  skeleton  dVxc/dR  (fixed orbitals, basis+grid move)  -> in F0x
        !   (ii) fxc[d0], d0 = reorthonormalization density            -> in Gd0
        ! Both enter B_ai with the SAME (minus) sign as the other Fock terms, so
        ! they are captured together by ONE central FD of the XC Fock matrix
        ! (dftexcor) along the combined path: geometry R +/- h AND occupied MOs
        ! reorthonormalized by dmo_i = -1/2 sum_j C_j S^x_ji.  dftexcor handles all
        ! density/spin scale factors internally, so no manual convention factors.
        if (infos%control%hamilton == 20) then
          block
            use dft, only: dft_initialize, dftclean, dftexcor
            use mod_dft_molgrid, only: dft_grid_t
            type(dft_grid_t) :: mgr
            real(dp), allocatable :: dmoR(:,:), mop(:,:), frp(:), frm(:), dVxcR(:,:), hxcR(:,:)
            real(dp) :: hxr, telr, tknr, exr
            integer :: ir, jr
            allocate(dmoR(nbf,nocc), mop(nbf,nbf), frp(nbf2), frm(nbf2), dVxcR(nbf,nbf), hxcR(nbf,nbf))
            hxr = 1.0d-3
            dmoR = 0.0_dp
            do ir = 1, nocc
              do jr = 1, nocc
                dmoR(:,ir) = dmoR(:,ir) - 0.5_dp*mo_a(:,jr)*Sx(jr,ir)
              end do
            end do
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) + hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mop = mo_a; mop(:,1:nocc) = mo_a(:,1:nocc) + hxr*dmoR
            frp = 0.0_dp
            call dftexcor(basis, mgr, 1, frp, frp, mop, mop, nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) - 2*hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mop = mo_a; mop(:,1:nocc) = mo_a(:,1:nocc) - hxr*dmoR
            frm = 0.0_dp
            call dftexcor(basis, mgr, 1, frm, frm, mop, mop, nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) + hxr
            call basis%init_shell_centers()
            call unpack_from_packed((frp - frm)/(2*hxr), dVxcR, nbf)
            call mo_transform(mo_a, dVxcR, nbf, scr, col, hxcR)
            do a = 1, nvir
              do i = 1, nocc
                F0x(i,nocc+a) = F0x(i,nocc+a) + hxcR(i,nocc+a)
              end do
            end do
            deallocate(dmoR, mop, frp, frm, dVxcR, hxcR)
          end block
        end if

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

    ! ===== CPHF orbital-relaxation response (PySCF hessian/rhf.hess_elec) =====
    ! H^resp_xy = 4 Tr[F^x dm1^y] - 4 Tr[S^x (eps.dm1^y)] - 2 Tr[s1oo^x mo_e1^y]
    !   dm1^y_pq      = sum_k dC^y_pk C_qk                       (one-sided)
    !   mo_e1^y_kl    = (h^y + G[P]^y + G[dP^y])^MO_kl - 1/2 (eps_k+eps_l) s1oo^y_kl
    ! F^x = h^x + G[P]^x; dC^y from the validated CPHF amplitudes U^y. The first
    ! two terms equal Tr[dP^y F^x] and the eps-weighted overlap term; the third
    ! is the FULL occ-occ energy-weighted term (the off-diagonal part is what a
    ! diagonal dε approximation misses). 2e traces use fock_deriv_contract
    ! (=1/2 Tr[M G[P]^x]) and fock_jk (G[dP^y]).
    allocate(hess_native(ncart,ncart), source=0.0_dp)
    block
      real(dp), allocatable :: sflat(:,:,:), hflat(:,:,:)
      real(dp), allocatable :: dCx(:,:,:), dPx(:,:,:), Gdp(:,:,:)
      real(dp), allocatable :: s1oo(:,:,:), hMOoo(:,:,:), GdpMOoo(:,:,:), moe1a(:,:,:)
      real(dp), allocatable :: Mi(:,:), gxy(:,:), A2(:,:), tGP(:,:), hresp(:,:)
      real(dp), allocatable :: s1(:,:), s2(:,:), bMO(:,:), dpp(:,:), gpp(:,:), gfl(:,:)
      real(dp), allocatable :: cocc(:,:), tmpno(:,:)
      real(dp) :: a1v, a3v, t3a, dcsx
      integer :: x, yy, ii, jj, kk, ll, aa, ia2, mu2, nu2, ccx, kcx

      allocate(sflat(nbf,nbf,ncart), hflat(nbf,nbf,ncart))
      do x = 1, ncart
        ccx = mod(x-1,3)+1; kcx = (x-1)/3+1
        sflat(:,:,x) = dSa(:,:,ccx,kcx)
        hflat(:,:,x) = dTa(:,:,ccx,kcx) + dVa(:,:,ccx,kcx)
      end do
      allocate(cocc(nbf,nocc)); cocc = mo_a(:,1:nocc)

      ! occ-occ MO blocks of S^x and h^x
      allocate(s1oo(nocc,nocc,ncart), hMOoo(nocc,nocc,ncart), source=0.0_dp)
      allocate(s1(nbf,nbf), s2(nbf,nbf), bMO(nbf,nbf), tmpno(nbf,nocc))
      do x = 1, ncart
        call dgemm('n','n',nbf,nocc,nbf,1.0_dp,sflat(:,:,x),nbf,cocc,nbf,0.0_dp,tmpno,nbf)
        call dgemm('t','n',nocc,nocc,nbf,1.0_dp,cocc,nbf,tmpno,nbf,0.0_dp,s1oo(:,:,x),nocc)
        call dgemm('n','n',nbf,nocc,nbf,1.0_dp,hflat(:,:,x),nbf,cocc,nbf,0.0_dp,tmpno,nbf)
        call dgemm('t','n',nocc,nocc,nbf,1.0_dp,cocc,nbf,tmpno,nbf,0.0_dp,hMOoo(:,:,x),nocc)
      end do

      ! relaxed orbital derivative dC^y, density dP^y (total), response Fock G[dP^y]
      allocate(dCx(nbf,nocc,ncart), dPx(nbf,nbf,ncart), Gdp(nbf,nbf,ncart), source=0.0_dp)
      allocate(GdpMOoo(nocc,nocc,ncart), source=0.0_dp)
      allocate(dpp(nbf2,1), gpp(nbf2,1), gfl(nbf,nbf))
      do yy = 1, ncart
        ia2 = 0
        do aa = 1, nvir
          do ii = 1, nocc
            ia2 = ia2 + 1
            dCx(:,ii,yy) = dCx(:,ii,yy) + mo_a(:,nocc+aa)*uvec(ia2,yy)
          end do
        end do
        do ii = 1, nocc
          do jj = 1, nocc
            dCx(:,ii,yy) = dCx(:,ii,yy) - 0.5_dp*mo_a(:,jj)*s1oo(jj,ii,yy)
          end do
        end do
        do ii = 1, nocc
          do mu2 = 1, nbf
            do nu2 = 1, nbf
              dPx(mu2,nu2,yy) = dPx(mu2,nu2,yy) &
                + 2.0_dp*(dCx(mu2,ii,yy)*mo_a(nu2,ii) + mo_a(mu2,ii)*dCx(nu2,ii,yy))
            end do
          end do
        end do
        call pack_matrix(dPx(:,:,yy), dpp(:,1))
        gpp = 0.0_dp
        call fock_jk(basis, d=dpp, f=gpp, scale_exch=hfscale, infos=infos)
        call unpack_from_packed(gpp(:,1), gfl, nbf); Gdp(:,:,yy) = gfl
        call dgemm('n','n',nbf,nocc,nbf,1.0_dp,gfl,nbf,cocc,nbf,0.0_dp,tmpno,nbf)
        call dgemm('t','n',nocc,nocc,nbf,1.0_dp,cocc,nbf,tmpno,nbf,0.0_dp,GdpMOoo(:,:,yy),nocc)
      end do

      ! mo_e1 without the G[P]^y part (added via Mi trick in term3)
      allocate(moe1a(nocc,nocc,ncart))
      do yy = 1, ncart
        do ll = 1, nocc
          do kk = 1, nocc
            moe1a(kk,ll,yy) = hMOoo(kk,ll,yy) + GdpMOoo(kk,ll,yy) &
                            - 0.5_dp*(eps(kk)+eps(ll))*s1oo(kk,ll,yy)
          end do
        end do
      end do

      ! 2e traces: A2(x,y)=Tr[dP^y G[P]^x]; tGP(x,y)=Tr[M^x G[P]^y]
      ! with M^x = sum_kl s1oo^x_kl C_k C_l^T
      allocate(gxy(3,natom), A2(ncart,ncart), tGP(ncart,ncart), Mi(nbf,nbf), source=0.0_dp)
      do yy = 1, ncart
        gxy = 0.0_dp
        call fock_deriv_contract(infos, basis, pfull, dPx(:,:,yy), hfscale, gxy)
        A2(:,yy) = 2.0_dp*reshape(gxy, [ncart])
      end do
      do x = 1, ncart
        call dgemm('n','n',nbf,nocc,nocc,1.0_dp,cocc,nbf,s1oo(:,:,x),nocc,0.0_dp,tmpno,nbf)
        call dgemm('n','t',nbf,nbf,nocc,1.0_dp,tmpno,nbf,cocc,nbf,0.0_dp,Mi,nbf)
        gxy = 0.0_dp
        call fock_deriv_contract(infos, basis, pfull, Mi, hfscale, gxy)
        tGP(x,:) = 2.0_dp*reshape(gxy, [ncart])
      end do

      ! assemble response  hresp(x,y) = 4Tr[F^x dm1^y]-4Tr[S^x eps.dm1^y]-2Tr[s1oo^x mo_e1^y]
      !   = (Tr[dP^y h^x] + A2) - 4 A3 - 2 (sum_kl s1oo^x_kl moe1a^y_kl) - 2 tGP
      allocate(hresp(ncart,ncart), source=0.0_dp)
      do x = 1, ncart
        do yy = 1, ncart
          a1v = sum(dPx(:,:,yy)*hflat(:,:,x))
          a3v = 0.0_dp
          do ii = 1, nocc
            dcsx = 0.0_dp
            do mu2 = 1, nbf
              do nu2 = 1, nbf
                dcsx = dcsx + dCx(mu2,ii,yy)*sflat(mu2,nu2,x)*mo_a(nu2,ii)
              end do
            end do
            a3v = a3v + eps(ii)*dcsx
          end do
          t3a = 0.0_dp
          do ll = 1, nocc
            do kk = 1, nocc
              t3a = t3a + s1oo(kk,ll,x)*moe1a(kk,ll,yy)
            end do
          end do
          hresp(x,yy) = (a1v + A2(x,yy)) - 4.0_dp*a3v - 2.0_dp*t3a - 2.0_dp*tGP(x,yy)
        end do
      end do

      hess_native = 0.5_dp*(hresp + transpose(hresp))

      ! --- DFT exchange-correlation second-derivative contribution -----------
      ! The XC part of the Hessian is obtained by central finite differencing the
      ! analytic XC nuclear gradient (derexc_blk) over geometry while displacing
      ! the density by the analytic relaxed density derivative dP^y. This adds
      ! both the XC skeleton (d2Exc/dR2 at fixed density) and the XC response
      ! (through dP^y) in one shot, with no re-SCF. The HF-exchange fraction is
      ! already in the Coulomb/exchange terms above (hfscale); derexc_blk
      ! supplies the remaining DFT exchange-correlation functional.
      if (infos%control%hamilton == 20) then
        block
          use dft, only: dft_initialize, dftclean, dftexcor
          use mod_dft_gridint_grad, only: derexc_blk
          use mod_dft_molgrid, only: dft_grid_t
          type(dft_grid_t) :: mg
          real(dp), allocatable :: dap(:,:), dedp(:,:), dedm(:,:)
          real(dp), allocatable :: mop(:,:), frp(:), frm(:), dFxc(:,:), dFoo(:,:)
          real(dp), allocatable :: tmpn(:,:), dHse(:,:), dHt3(:,:)
          real(dp) :: hx, tele, tkin, eexc
          integer :: yy2, ccy, kcy, nang, x2, kk2, ll2
          hx = 1.0d-3; nang = maxval(basis%am) + 2
          ! XC contribution split exactly as PySCF hessian/rks but realised through
          ! the OpenQP moving-grid XC machinery so it stays consistent with the
          ! OpenQP numerical Hessian:
          !   dHse : skeleton + density-response (term1).  Central FD of the analytic
          !          XC gradient (derexc) along the relaxed path R+lambda, P+lambda*dP.
          !          This is the genuine total derivative d/dR[g_XC(R,P(R))] of the
          !          OpenQP XC gradient, so the moving-grid weight derivatives are
          !          handled identically to the SCF/numerical-gradient convention.
          !   dHt3 : -2 Tr[s1oo^x (vxc^y+fxc[dP^y])_oo], the XC part of the
          !          energy-weighted (mo_e1) term, from the FD of the XC Fock
          !          matrix (dftexcor) along the same relaxed orbital path.
          allocate(dap(nbf,nbf), dedp(3,natom), dedm(3,natom))
          allocate(mop(nbf,nbf), frp(nbf2), frm(nbf2), dFxc(nbf,nbf), dFoo(nocc,nocc))
          allocate(tmpn(nbf,nocc), dHse(ncart,ncart), dHt3(ncart,ncart))
          dHt3 = 0.0_dp
          ! warm-up to flush any stale grid state left by the CPHF solver
          call dft_initialize(infos, basis, mg); call dftclean(infos)
          do yy2 = 1, ncart
            ccy = mod(yy2-1,3)+1; kcy = (yy2-1)/3+1
            basis%atoms%xyz(ccy,kcy) = basis%atoms%xyz(ccy,kcy) + hx
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mg)
            dap = pfull + hx*dPx(:,:,yy2); dedp = 0.0_dp            ! skeleton + density response
            call derexc_blk(basis, mg, dap, dap, dedp, tele, tkin, nang, nbf, &
                            infos%dft%grid_density_cutoff, .false., infos)
            mop = mo_a; mop(:,1:nocc) = mo_a(:,1:nocc) + hx*dCx(:,:,yy2)
            call dftexcor(basis, mg, 1, frp, frp, mop, mop, nbf, nbf2, eexc, tele, tkin, infos)
            call dftclean(infos)
            basis%atoms%xyz(ccy,kcy) = basis%atoms%xyz(ccy,kcy) - 2*hx
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mg)
            dap = pfull - hx*dPx(:,:,yy2); dedm = 0.0_dp
            call derexc_blk(basis, mg, dap, dap, dedm, tele, tkin, nang, nbf, &
                            infos%dft%grid_density_cutoff, .false., infos)
            mop = mo_a; mop(:,1:nocc) = mo_a(:,1:nocc) - hx*dCx(:,:,yy2)
            call dftexcor(basis, mg, 1, frm, frm, mop, mop, nbf, nbf2, eexc, tele, tkin, infos)
            call dftclean(infos)
            basis%atoms%xyz(ccy,kcy) = basis%atoms%xyz(ccy,kcy) + hx
            call basis%init_shell_centers()
            dHse(:,yy2) = reshape((dedp - dedm)/(2*hx), [ncart])
            ! term3: -2 s1oo^x (vxc^y + fxc[dP^y])_oo
            call unpack_from_packed((frp - frm)/(2*hx), dFxc, nbf)
            call dgemm('n','n',nbf,nocc,nbf,1.0_dp,dFxc,nbf,mo_a,nbf,0.0_dp,tmpn,nbf)
            call dgemm('t','n',nocc,nocc,nbf,1.0_dp,mo_a,nbf,tmpn,nbf,0.0_dp,dFoo,nocc)
            do x2 = 1, ncart
              do ll2 = 1, nocc
                do kk2 = 1, nocc
                  dHt3(x2,yy2) = dHt3(x2,yy2) - 2.0_dp*s1oo(kk2,ll2,x2)*dFoo(kk2,ll2)
                end do
              end do
            end do
          end do
          hess_native = hess_native + 0.5_dp*(dHse + transpose(dHse)) &
                                    + 0.5_dp*(dHt3 + transpose(dHt3))
          deallocate(dap, dedp, dedm, mop, frp, frm, dFxc, dFoo, tmpn, dHse, dHt3)
        end block
      end if

      deallocate(sflat, hflat, dCx, dPx, Gdp, s1oo, hMOoo, GdpMOoo, moe1a, &
                 Mi, gxy, A2, tGP, hresp, s1, s2, bMO, dpp, gpp, gfl, cocc, tmpno)
    end block

    call hess_nn(basis%atoms, basis%ecp_zn_num, hess_native)

    ! --- One-electron + Pulay second-derivative skeleton (fixed density) ------
    ! Mirrors the production HF gradient assembly (hf_1e_grad): the analytic
    ! Hessian skeleton is d/dx of [grad_ee_overlap(W) + grad_ee_kinetic(P)
    ! + grad_en(P)] evaluated at the fixed converged density, i.e. the
    ! second-derivative integral contractions hess_ee_overlap / hess_ee_kinetic
    ! / hess_en.  This is distinct from (and additive to) the CPHF response
    ! term above; the 2e ERI second-derivative skeleton is added separately.
    block
      use grd1, only: eijden, hess_ee_overlap, hess_ee_kinetic, hess_en
      real(kind=dp), allocatable :: wlag(:), pden(:), hcc(:,:)
      allocate(wlag(nbf2), pden(nbf2), hcc(ncart,ncart), source=0.0_dp)
      call eijden(wlag, nbf, infos)                 ! energy-weighted (Lagrangian) density
      pden = dmat_a                                 ! total density (closed-shell RHF)
      call hess_ee_overlap(basis, wlag, hess_native)            ! overlap / Pulay
      call hess_ee_kinetic(basis, pden, hess_native)            ! kinetic
      call hess_en(basis, basis%atoms%xyz, &
                   basis%atoms%zn - basis%ecp_zn_num, pden, hess_native, hess_cc=hcc)
      deallocate(wlag, pden, hcc)
    end block

    ! --- Two-electron (ERI) second-derivative skeleton (fixed density) --------
    ! d^2/dR^2 of the analytic 2e gradient contraction at the converged density,
    ! i.e. sum P P d^2/dR^2 [ (ij|kl) - 1/4 c_x (ik|jl) ]. Validated against a
    ! finite difference of grd2_driver (see grd2_hess_selftest). Additive to the
    ! CPHF response and 1e skeleton above.
    block
      use grd2, only: grd2_hess_driver, grd2_compute_data_t
      use hf_gradient_mod, only: grd2_rhf_compute_data_t
      class(grd2_compute_data_t), allocatable :: gcomp
      gcomp = grd2_rhf_compute_data_t( da = dmat_a, hfscale = hfscale, nbf = nbf )
      call gcomp%init()
      call grd2_hess_driver(infos, basis, hess_native, gcomp)
      call gcomp%clean()
    end block

    call infos%dat%reserve_data(OQP_hf_hessian, TA_TYPE_REAL64, ncart*ncart, (/ ncart, ncart /), &
      comment='Native OpenQP HF/DFT analytic Hessian matrix')
    call tagarray_get_data(infos%dat, OQP_hf_hessian, hess_store)
    hess_store = hess_native
    write(iw,'(A)') 'PyOQP: Native OpenQP HF/DFT Hessian matrix stored'

    deallocate(pfull, dSa, dTa, dVa, scr, col, Sx, hx, F0x, Gd0, probe, gx, &
               d0, d0p, gp, gfull, bvec, uvec, hess_native)
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
