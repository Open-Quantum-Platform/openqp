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
    use messages, only: show_message, WITH_ABORT

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

    ! Unsupported-feature guards (apply to ALL references, RHF/RKS included).
    ! Effective-core-potential (ECP) second derivatives ARE supported: RHF/UHF
    ! contract the ECP skeleton d^2 V_ECP/dR^2 analytically (add_ecphess, libecpint
    ! deriv order 2) plus the ECP core-derivative in the CPHF response; ROHF folds
    ! the ECP gradient (add_ecpder) into its semi-numerical resp_grad.
    ! Range-separated (CAM/LC) functionals are also supported: the 2e derivative
    ! integrals are erfc-attenuation capable, so grd2_hess_driver (skeleton),
    ! grd2_driver (fock_deriv_contract response) and fock_jk (cphf) all run the
    ! long-range Coulomb + short-range erfc-exchange two-pass split when
    ! infos%dft%cam_flag is set.

    ! Open-shell (UHF/ROHF) dispatch.  The body below is the closed-shell
    ! (RHF/RKS) kernel: it reads only the alpha density/MOs (OQP_DM_A, mo_a, eps)
    ! and treats nocc as doubly occupied, so it must never run on an open-shell
    ! SCF.  UHF (scftype==2) -> hf_hessian_uhf, ROHF (scftype==3) -> hf_hessian_rohf
    ! (both HF and DFT, finite-difference validated).
    if (infos%control%scftype == 2) then
      call hf_hessian_uhf(infos)
      return
    else if (infos%control%scftype == 3) then
      call hf_hessian_rohf(infos)
      return
    else if (infos%control%scftype > 3) then
      call show_message('Native analytic Hessian supports RHF/RKS, UHF (HF) '// &
        'and ROHF (HF) references only for this scftype. Use [hess] '// &
        'type=numerical.', WITH_ABORT)
    end if

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

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/,A)') 'PyOQP: Native OpenQP HF/DFT Hessian CPHF response prepass'
    write(iw,'(A,I6,A,I6,A,I6,A,I6)') '  nbf=', nbf, ' nocc=', nocc, ' nvir=', nvir, ' rhs=', ncart
    write(iw,'(A)') '  Storing native OpenQP HF/DFT analytic Hessian matrix in OQP::hf_hessian.'

    if (nocc <= 0 .or. nvir <= 0 .or. ncart <= 0) then
      write(iw,'(A)') '  Native CPHF prepass skipped: empty occupied/virtual/nuclear space.'
      close(iw)
      return
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, eps)

    allocate(pfull(nbf,nbf)); call unpack_matrix(dmat_a, pfull)
    allocate(dSa(nbf,nbf,3,natom), dTa(nbf,nbf,3,natom), dVa(nbf,nbf,3,natom))
    call der_overlap_matrix(basis, dSa)
    call der_kinetic_matrix(basis, dTa)
    call der_nucattr_matrix(basis, basis%atoms%xyz, &
                            basis%atoms%zn - basis%ecp_zn_num, dVa)  ! ECP-screened point charge

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

    ! ECP first-derivative integrals enter the core-Hamiltonian derivative
    ! dHcore/dR (added into dVa, the nuclear-attraction derivative tensor), so the
    ! ECP contributes to the CPHF right-hand side and the orbital-relaxation
    ! response exactly as point-charge nuclear attraction does.  libecpint returns
    ! these already in the OpenQP normalized convention, hence added AFTER the
    ! bfnrm scaling above.  No-op for non-ECP bases.
    block
      use ecp_tool, only: ecp_deriv_ints
      real(kind=dp), allocatable :: dVecp(:,:,:,:)
      allocate(dVecp(nbf,nbf,3,natom))
      call ecp_deriv_ints(basis, basis%atoms%xyz, dVecp)
      dVa = dVa + dVecp
      deallocate(dVecp)
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

    ! ===== CPHF orbital-relaxation response =====
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
          ! XC contribution split into a skeleton+density-response term and an
          ! energy-weighting term, realised through the OpenQP moving-grid XC
          ! machinery so it stays consistent with the OpenQP numerical Hessian:
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
      use ecp_tool, only: add_ecphess
      real(kind=dp), allocatable :: wlag(:), pden(:), hcc(:,:)
      allocate(wlag(nbf2), pden(nbf2), hcc(ncart,ncart), source=0.0_dp)
      call eijden(wlag, nbf, infos)                 ! energy-weighted (Lagrangian) density
      pden = dmat_a                                 ! total density (closed-shell RHF)
      call hess_ee_overlap(basis, wlag, hess_native)            ! overlap / Pulay
      call hess_ee_kinetic(basis, pden, hess_native)            ! kinetic
      call hess_en(basis, basis%atoms%xyz, &
                   basis%atoms%zn - basis%ecp_zn_num, pden, hess_native, hess_cc=hcc)
      call add_ecphess(basis, basis%atoms%xyz, pden, hess_native) ! ECP skeleton (if any)
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
    close(iw)

    deallocate(pfull, dSa, dTa, dVa, scr, col, Sx, hx, F0x, Gd0, probe, gx, &
               d0, d0p, gp, gfull, bvec, uvec, hess_native)
  end subroutine hf_hessian

!###############################################################################

  subroutine hf_hessian_uhf(infos)
    ! Native open-shell (UHF) analytic HF Hessian.
    !
    ! Mirrors the closed-shell hf_hessian response assembly per spin, summed over
    ! s in {alpha, beta} with single (not doubled) occupation factors.  Each spin
    ! uses its own MO set C^s, orbital energies eps^s and density P^s; the
    ! two-electron couplings are open-shell (Coulomb from the total density
    ! P = Pa + Pb, exchange from the spin density P^s):
    !
    !   B^s_ia   = -(h^x_ia + G^{s,x}[P]_ia) + eps^s_i S^x_ia - G^s[d0]_ia ,
    !   d0^s     = -sum_ij S^x,s_ij C^s_i C^s_j^T          (reorthonormalization),
    !   G^s[.]   = J[.^a + .^b] - c_x K[.^s]               (scf_addons::fock_jk),
    !   G^{s,x}[P] via fock_deriv_mod::fock_deriv_contract_os (Coulomb P, exch P^s).
    !
    ! The 3N right-hand sides are solved with cphf_mod::cphf_solve_uhf, and the
    ! orbital-relaxation response is assembled as (per spin, summed):
    !
    !   H^resp_xy = sum_s [ Tr[dP^s,y h^x] + Tr[dP^s,y G^{s,x}[P]] ]
    !             - 2 sum_s sum_i eps^s_i (dC^s,y_i . S^x . C^s_i)
    !             -   sum_s sum_kl s1oo^s,x_kl moe1^s,y_kl
    !             -   sum_s Tr[Mi^s,x G^{s,y}[P]] ,
    !
    !   moe1^s,y_kl = h^x_kl(MO) + G^s[dP^y]_kl(MO) - 1/2(eps^s_k+eps^s_l) s1oo^s,y_kl ,
    !   Mi^s,x      = sum_kl s1oo^s,x_kl C^s_k C^s_l^T .
    !
    ! The fixed-density skeleton (1e total density + open-shell Lagrangian W, 2e
    ! via grd2_uhf_compute_data_t) and the nuclear-repulsion term are added on
    ! top, exactly as in hess_skel_open_selftest.  HF only (the UKS f_xc response
    ! is not finite-difference validated).
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_DM_B, &
      OQP_VEC_MO_A, OQP_VEC_MO_B, OQP_E_MO_A, OQP_E_MO_B, OQP_hf_hessian, TA_TYPE_REAL64
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix, hess_nn
    use fock_deriv_mod, only: fock_deriv_contract_os
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve_uhf
    use io_constants, only: iw

    implicit none

    type(information), target, intent(inout) :: infos

    !> Per-spin work container (alpha/beta have different nocc/nvir).
    type :: uhf_spin_t
      real(dp), allocatable :: mo(:,:)            ! MO coefficients (nbf,nbf)
      real(dp), allocatable :: eps(:)             ! orbital energies (nbf)
      real(dp), allocatable :: p(:,:)             ! spin AO density (nbf,nbf)
      integer :: nocc = 0, nvir = 0, loff = 0     ! occ/vir count, CPHF block offset
      real(dp), allocatable :: s1oo(:,:,:)        ! occ-occ MO of S^x (nocc,nocc,ncart)
      real(dp), allocatable :: hoo(:,:,:)         ! occ-occ MO of h^x
      real(dp), allocatable :: g2e(:,:)           ! G^{s,x}[P]_ia for all coords (nocc*nvir,ncart)
      real(dp), allocatable :: dCx(:,:,:)         ! relaxed dC (nbf,nocc,ncart)
      real(dp), allocatable :: dPx(:,:,:)         ! relaxed spin density derivative
      real(dp), allocatable :: gdpoo(:,:,:)       ! occ-occ MO of G^s[dP^y]
      real(dp), allocatable :: moe1(:,:,:)        ! occ-occ energy-weighted derivative
    end type

    type(basis_set), pointer :: basis
    real(dp), contiguous, pointer :: dma(:), dmb(:), moa(:,:), mob(:,:), epsa(:), epsb(:)
    real(dp), contiguous, pointer :: hess_store(:,:)
    real(dp), allocatable :: ptot(:,:), dSa(:,:,:,:), dTa(:,:,:,:), dVa(:,:,:,:)
    real(dp), allocatable :: sflat(:,:,:), hflat(:,:,:)
    real(dp), allocatable :: bvec(:,:), uvec(:,:), hess_native(:,:)
    real(dp), allocatable :: scr(:,:), tmp(:,:), gx(:,:), probe(:,:)
    real(dp), allocatable :: SxMO(:,:), hxMO(:,:), d0a(:,:), d0b(:,:)
    real(dp), allocatable :: dpck(:,:), fpck(:,:), gfull(:,:)
    real(dp), allocatable :: Gd0(:,:), Mi(:,:)
    real(dp), allocatable :: A2(:,:), tGP(:,:), hresp(:,:)
    type(uhf_spin_t) :: sp(2)
    real(dp) :: hfscale, a1v, a3v, t3a, dcsx
    integer :: nbf, nbf2, natom, ncart, nocca, noccb, nvira, nvirb, la, lb, ltot
    integer :: s, i, j, a, ia, icart, kc, cc, x, yy, kk, ll, mu, nu

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2)
    ncart = 3*natom
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    la = nocca*nvira
    lb = noccb*nvirb
    ltot = la + lb
    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale

    write(iw,'(/,A)') 'PyOQP: Native OpenQP open-shell (UHF) HF Hessian CPHF response prepass'
    write(iw,'(A,I6,A,I6,A,I6,A,I6,A,I6)') '  nbf=', nbf, ' nocca=', nocca, &
      ' noccb=', noccb, ' rhs=', ncart, ' ltot=', ltot
    write(iw,'(A)') '  Storing native OpenQP open-shell HF analytic Hessian in OQP::hf_hessian.'

    if (ncart <= 0 .or. (la <= 0 .and. lb <= 0)) then
      write(iw,'(A)') '  UHF CPHF prepass skipped: empty occupied/virtual/nuclear space.'
      return
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dma)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmb)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, moa)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mob)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, epsa)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, epsb)

    ! per-spin containers
    sp(1)%nocc = nocca; sp(1)%nvir = nvira; sp(1)%loff = 0
    sp(2)%nocc = noccb; sp(2)%nvir = nvirb; sp(2)%loff = la
    allocate(sp(1)%mo(nbf,nbf), sp(1)%eps(nbf), sp(1)%p(nbf,nbf))
    allocate(sp(2)%mo(nbf,nbf), sp(2)%eps(nbf), sp(2)%p(nbf,nbf))
    sp(1)%mo = moa; sp(1)%eps = epsa
    sp(2)%mo = mob; sp(2)%eps = epsb
    call unpack_matrix(dma, sp(1)%p)
    call unpack_matrix(dmb, sp(2)%p)
    allocate(ptot(nbf,nbf)); ptot = sp(1)%p + sp(2)%p

    ! derivative integrals (normalized into the bfnrm convention of the MOs)
    allocate(dSa(nbf,nbf,3,natom), dTa(nbf,nbf,3,natom), dVa(nbf,nbf,3,natom))
    call der_overlap_matrix(basis, dSa)
    call der_kinetic_matrix(basis, dTa)
    call der_nucattr_matrix(basis, basis%atoms%xyz, &
                            basis%atoms%zn - basis%ecp_zn_num, dVa)  ! ECP-screened point charge
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

    ! ECP first-derivative integrals -> core-Hamiltonian derivative dHcore/dR (see
    ! the RHF kernel for the rationale).  Already in the normalized convention, so
    ! added after the bfnrm scaling.  No-op for non-ECP bases.
    block
      use ecp_tool, only: ecp_deriv_ints
      real(dp), allocatable :: dVecp(:,:,:,:)
      allocate(dVecp(nbf,nbf,3,natom))
      call ecp_deriv_ints(basis, basis%atoms%xyz, dVecp)
      dVa = dVa + dVecp
      deallocate(dVecp)
    end block

    ! flat (ncart) AO views of S^x and h^x = (T+V)^x
    allocate(sflat(nbf,nbf,ncart), hflat(nbf,nbf,ncart))
    do x = 1, ncart
      cc = mod(x-1,3)+1; kc = (x-1)/3+1
      sflat(:,:,x) = dSa(:,:,cc,kc)
      hflat(:,:,x) = dTa(:,:,cc,kc) + dVa(:,:,cc,kc)
    end do

    ! occ-occ and occ MO transforms needed by the response assembly
    allocate(scr(nbf,nbf), tmp(nbf,nbf), SxMO(nbf,nbf), hxMO(nbf,nbf))
    do s = 1, 2
      allocate(sp(s)%s1oo(sp(s)%nocc, sp(s)%nocc, ncart), source=0.0_dp)
      allocate(sp(s)%hoo (sp(s)%nocc, sp(s)%nocc, ncart), source=0.0_dp)
      do x = 1, ncart
        call mo_transform(sp(s)%mo, sflat(:,:,x), nbf, scr, tmp, SxMO)
        call mo_transform(sp(s)%mo, hflat(:,:,x), nbf, scr, tmp, hxMO)
        sp(s)%s1oo(:,:,x) = SxMO(1:sp(s)%nocc,1:sp(s)%nocc)
        sp(s)%hoo (:,:,x) = hxMO(1:sp(s)%nocc,1:sp(s)%nocc)
      end do
    end do

    ! ===== CPHF right-hand sides B^s (occ-vir) for all 3N perturbations =====
    allocate(bvec(ltot,ncart), uvec(ltot,ncart), source=0.0_dp)
    allocate(probe(nbf,nbf), gx(3,natom))
    allocate(d0a(nbf,nbf), d0b(nbf,nbf), gfull(nbf,nbf), Gd0(nbf,nbf))
    allocate(dpck(nbf2,2), fpck(nbf2,2))

    ! 2e response-Fock skeleton  G^{s,x}[P]_ia  for ALL 3N coordinates.  The
    ! occ-vir probe C^s_a C^s_i^T is geometry-independent, so a single open-shell
    ! derivative-Fock contraction per occ-vir pair yields every Cartesian
    ! component at once (avoids an ncart-fold redundant grd2 sweep).
    do s = 1, 2
      allocate(sp(s)%g2e(sp(s)%nocc*sp(s)%nvir, ncart), source=0.0_dp)
      do a = 1, sp(s)%nvir
        do i = 1, sp(s)%nocc
          do mu = 1, nbf
            do nu = 1, nbf
              probe(mu,nu) = 0.5_dp*( sp(s)%mo(mu,sp(s)%nocc+a)*sp(s)%mo(nu,i) &
                                    + sp(s)%mo(mu,i)*sp(s)%mo(nu,sp(s)%nocc+a) )
            end do
          end do
          gx = 0.0_dp
          call fock_deriv_contract_os(infos, basis, ptot, sp(s)%p, probe, hfscale, gx)
          ia = (a-1)*sp(s)%nocc + i
          sp(s)%g2e(ia,:) = reshape(gx, [ncart])
        end do
      end do
    end do

    icart = 0
    do kc = 1, natom
      do cc = 1, 3
        icart = icart + 1

        ! reorthonormalization density per spin: d0^s = -sum_ij S^x,s_ij C^s_i C^s_j^T
        d0a = 0.0_dp; d0b = 0.0_dp
        do s = 1, 2
          call mo_transform(sp(s)%mo, dSa(:,:,cc,kc), nbf, scr, tmp, SxMO)
          do i = 1, sp(s)%nocc
            do j = 1, sp(s)%nocc
              do mu = 1, nbf
                do nu = 1, nbf
                  if (s == 1) then
                    d0a(mu,nu) = d0a(mu,nu) - SxMO(i,j)*sp(s)%mo(mu,i)*sp(s)%mo(nu,j)
                  else
                    d0b(mu,nu) = d0b(mu,nu) - SxMO(i,j)*sp(s)%mo(mu,i)*sp(s)%mo(nu,j)
                  end if
                end do
              end do
            end do
          end do
        end do
        call pack_matrix(d0a, dpck(:,1))
        call pack_matrix(d0b, dpck(:,2))
        fpck = 0.0_dp
        call fock_jk(basis, d=dpck, f=fpck, scale_exch=hfscale, infos=infos)

        do s = 1, 2
          call mo_transform(sp(s)%mo, dSa(:,:,cc,kc), nbf, scr, tmp, SxMO)
          call mo_transform(sp(s)%mo, dTa(:,:,cc,kc)+dVa(:,:,cc,kc), nbf, scr, tmp, hxMO)
          call unpack_from_packed(fpck(:,s), gfull, nbf)   ! G^s[d0]
          call mo_transform(sp(s)%mo, gfull, nbf, scr, tmp, Gd0)
          do a = 1, sp(s)%nvir
            do i = 1, sp(s)%nocc
              ia = (a-1)*sp(s)%nocc + i
              bvec(sp(s)%loff+ia,icart) = &
                  -(hxMO(i,sp(s)%nocc+a) + sp(s)%g2e(ia,icart)) &
                  + sp(s)%eps(i)*SxMO(i,sp(s)%nocc+a) &
                  - Gd0(i,sp(s)%nocc+a)
            end do
          end do
        end do

        ! --- XC contribution to the CPKS right-hand side (UKS only) -------------
        ! Mirror the closed-shell RKS RHS XC: one central FD of the spin XC Fock
        ! matrices (open-shell dftexcor) along R +/- h AND occupied MOs reorthonor-
        ! malized by dmoR^s = -1/2 sum_j C^s_j S^x,s_ji captures both the XC
        ! skeleton dVxc/dR and f_xc[d0]; subtract the vir-occ MO blocks from B
        ! (which carries -F0x).
        if (infos%control%hamilton >= 20) then
          block
            use dft, only: dft_initialize, dftclean, dftexcor
            use mod_dft_molgrid, only: dft_grid_t
            type(dft_grid_t) :: mgr
            real(dp), allocatable :: dmoa(:,:), dmob(:,:), mopa(:,:), mopb(:,:)
            real(dp), allocatable :: SxMOa(:,:), SxMOb(:,:)
            real(dp), allocatable :: frap(:), frbp(:), fram(:), frbm(:), dvx(:,:), hxc(:,:)
            real(dp) :: hxr, exr, telr, tknr
            integer :: ir, jr
            allocate(dmoa(nbf,nocca), dmob(nbf,noccb), mopa(nbf,nbf), mopb(nbf,nbf))
            allocate(SxMOa(nbf,nbf), SxMOb(nbf,nbf))
            allocate(frap(nbf2), frbp(nbf2), fram(nbf2), frbm(nbf2), dvx(nbf,nbf), hxc(nbf,nbf))
            hxr = 1.0d-3
            call mo_transform(sp(1)%mo, dSa(:,:,cc,kc), nbf, scr, tmp, SxMOa)
            call mo_transform(sp(2)%mo, dSa(:,:,cc,kc), nbf, scr, tmp, SxMOb)
            dmoa = 0.0_dp
            do ir = 1, nocca
              do jr = 1, nocca
                dmoa(:,ir) = dmoa(:,ir) - 0.5_dp*SxMOa(jr,ir)*sp(1)%mo(:,jr)
              end do
            end do
            dmob = 0.0_dp
            do ir = 1, noccb
              do jr = 1, noccb
                dmob(:,ir) = dmob(:,ir) - 0.5_dp*SxMOb(jr,ir)*sp(2)%mo(:,jr)
              end do
            end do
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) + hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mopa = sp(1)%mo; mopa(:,1:nocca) = sp(1)%mo(:,1:nocca) + hxr*dmoa
            mopb = sp(2)%mo; mopb(:,1:noccb) = sp(2)%mo(:,1:noccb) + hxr*dmob
            frap = 0.0_dp; frbp = 0.0_dp
            call dftexcor(basis, mgr, infos%control%scftype, frap, frbp, mopa, mopb, &
                          nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) - 2*hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mopa = sp(1)%mo; mopa(:,1:nocca) = sp(1)%mo(:,1:nocca) - hxr*dmoa
            mopb = sp(2)%mo; mopb(:,1:noccb) = sp(2)%mo(:,1:noccb) - hxr*dmob
            fram = 0.0_dp; frbm = 0.0_dp
            call dftexcor(basis, mgr, infos%control%scftype, fram, frbm, mopa, mopb, &
                          nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) + hxr
            call basis%init_shell_centers()
            call unpack_from_packed((frap - fram)/(2*hxr), dvx, nbf)
            call mo_transform(sp(1)%mo, dvx, nbf, scr, tmp, hxc)
            do a = 1, nvira
              do i = 1, nocca
                ia = (a-1)*nocca + i
                bvec(ia,icart) = bvec(ia,icart) - hxc(i,nocca+a)
              end do
            end do
            call unpack_from_packed((frbp - frbm)/(2*hxr), dvx, nbf)
            call mo_transform(sp(2)%mo, dvx, nbf, scr, tmp, hxc)
            do a = 1, nvirb
              do i = 1, noccb
                ia = (a-1)*noccb + i
                bvec(la+ia,icart) = bvec(la+ia,icart) - hxc(i,noccb+a)
              end do
            end do
            deallocate(dmoa, dmob, mopa, mopb, SxMOa, SxMOb, frap, frbp, fram, frbm, dvx, hxc)
          end block
        end if
      end do
    end do

    call cphf_solve_uhf(infos, ncart, bvec, uvec)

    ! ===== open-shell CPHF orbital-relaxation response =====
    ! relaxed dC^s, spin density derivative dP^s
    do s = 1, 2
      allocate(sp(s)%dCx(nbf, sp(s)%nocc, ncart), source=0.0_dp)
      allocate(sp(s)%dPx(nbf, nbf, ncart), source=0.0_dp)
      do yy = 1, ncart
        do a = 1, sp(s)%nvir
          do i = 1, sp(s)%nocc
            ia = (a-1)*sp(s)%nocc + i
            sp(s)%dCx(:,i,yy) = sp(s)%dCx(:,i,yy) &
              + sp(s)%mo(:,sp(s)%nocc+a)*uvec(sp(s)%loff+ia, yy)
          end do
        end do
        do i = 1, sp(s)%nocc
          do j = 1, sp(s)%nocc
            sp(s)%dCx(:,i,yy) = sp(s)%dCx(:,i,yy) - 0.5_dp*sp(s)%mo(:,j)*sp(s)%s1oo(j,i,yy)
          end do
        end do
        do i = 1, sp(s)%nocc
          do mu = 1, nbf
            do nu = 1, nbf
              sp(s)%dPx(mu,nu,yy) = sp(s)%dPx(mu,nu,yy) &
                + sp(s)%dCx(mu,i,yy)*sp(s)%mo(nu,i) + sp(s)%mo(mu,i)*sp(s)%dCx(nu,i,yy)
            end do
          end do
        end do
      end do
      allocate(sp(s)%gdpoo(sp(s)%nocc, sp(s)%nocc, ncart), source=0.0_dp)
      allocate(sp(s)%moe1 (sp(s)%nocc, sp(s)%nocc, ncart), source=0.0_dp)
    end do

    ! G^s[dP^y] (couples both spins via fock_jk) -> occ-occ MO block
    do yy = 1, ncart
      call pack_matrix(sp(1)%dPx(:,:,yy), dpck(:,1))
      call pack_matrix(sp(2)%dPx(:,:,yy), dpck(:,2))
      fpck = 0.0_dp
      call fock_jk(basis, d=dpck, f=fpck, scale_exch=hfscale, infos=infos)
      do s = 1, 2
        call unpack_from_packed(fpck(:,s), gfull, nbf)
        call mo_transform(sp(s)%mo, gfull, nbf, scr, tmp, hxMO)
        sp(s)%gdpoo(:,:,yy) = hxMO(1:sp(s)%nocc,1:sp(s)%nocc)
      end do
    end do

    ! energy-weighted derivative occ-occ block (without the G[P]^y piece, which
    ! is folded into tGP via the Mi^x probe below)
    do s = 1, 2
      do yy = 1, ncart
        do ll = 1, sp(s)%nocc
          do kk = 1, sp(s)%nocc
            sp(s)%moe1(kk,ll,yy) = sp(s)%hoo(kk,ll,yy) + sp(s)%gdpoo(kk,ll,yy) &
              - 0.5_dp*(sp(s)%eps(kk)+sp(s)%eps(ll))*sp(s)%s1oo(kk,ll,yy)
          end do
        end do
      end do
    end do

    ! 2e response traces, summed over spin:
    !   A2(x,y)  = sum_s Tr[dP^s,y G^{s,x}[P]]
    !   tGP(x,y) = sum_s Tr[Mi^s,x G^{s,y}[P]],  Mi^s,x = sum_kl s1oo^s,x_kl C^s_k C^s_l^T
    allocate(A2(ncart,ncart), tGP(ncart,ncart), Mi(nbf,nbf), source=0.0_dp)
    do yy = 1, ncart
      do s = 1, 2
        gx = 0.0_dp
        call fock_deriv_contract_os(infos, basis, ptot, sp(s)%p, sp(s)%dPx(:,:,yy), hfscale, gx)
        A2(:,yy) = A2(:,yy) + reshape(gx, [ncart])
      end do
    end do
    do x = 1, ncart
      do s = 1, 2
        Mi = 0.0_dp
        do ll = 1, sp(s)%nocc
          do kk = 1, sp(s)%nocc
            do mu = 1, nbf
              do nu = 1, nbf
                Mi(mu,nu) = Mi(mu,nu) + sp(s)%s1oo(kk,ll,x)*sp(s)%mo(mu,kk)*sp(s)%mo(nu,ll)
              end do
            end do
          end do
        end do
        gx = 0.0_dp
        call fock_deriv_contract_os(infos, basis, ptot, sp(s)%p, Mi, hfscale, gx)
        tGP(x,:) = tGP(x,:) + reshape(gx, [ncart])
      end do
    end do

    ! assemble  H^resp_xy
    allocate(hresp(ncart,ncart), source=0.0_dp)
    do x = 1, ncart
      do yy = 1, ncart
        a1v = 0.0_dp; a3v = 0.0_dp; t3a = 0.0_dp
        do s = 1, 2
          a1v = a1v + sum(sp(s)%dPx(:,:,yy)*hflat(:,:,x))
          do i = 1, sp(s)%nocc
            dcsx = 0.0_dp
            do mu = 1, nbf
              do nu = 1, nbf
                dcsx = dcsx + sp(s)%dCx(mu,i,yy)*sflat(mu,nu,x)*sp(s)%mo(nu,i)
              end do
            end do
            a3v = a3v + sp(s)%eps(i)*dcsx
          end do
          do ll = 1, sp(s)%nocc
            do kk = 1, sp(s)%nocc
              t3a = t3a + sp(s)%s1oo(kk,ll,x)*sp(s)%moe1(kk,ll,yy)
            end do
          end do
        end do
        hresp(x,yy) = (a1v + A2(x,yy)) - 2.0_dp*a3v - t3a - tGP(x,yy)
      end do
    end do

    allocate(hess_native(ncart,ncart))
    hess_native = 0.5_dp*(hresp + transpose(hresp))

    ! --- DFT (UKS) exchange-correlation second-derivative contribution --------
    ! Open-shell analog of the closed-shell RKS XC block:
    !   dHse : XC skeleton + density-response, central FD of the analytic
    !          open-shell XC gradient (derexc_blk) along R +/- h, P^s +/- h dP^s.
    !   dHt3 : -2 sum_s sum_kl s1oo^s,x_kl (vxc^s,y + fxc[dP^s,y])_kl, the XC part
    !          of the energy-weighted term, from the FD of the spin XC Fock
    !          (dftexcor) along R, C^s +/- h dC^s.
    ! The HF-exchange fraction is already in the Coulomb/exchange terms (hfscale).
    if (infos%control%hamilton >= 20) then
      block
        use dft, only: dft_initialize, dftclean, dftexcor
        use mod_dft_gridint_grad, only: derexc_blk
        use mod_dft_molgrid, only: dft_grid_t
        type(dft_grid_t) :: mg
        real(dp), allocatable :: dapa(:,:), dapb(:,:), dedp(:,:), dedm(:,:)
        real(dp), allocatable :: mopa(:,:), mopb(:,:), frap(:), frbp(:), fram(:), frbm(:)
        real(dp), allocatable :: dFxc(:,:), tmpn(:,:), dHse(:,:), dHt3(:,:), dFoo(:,:)
        real(dp) :: hx, exr, telr, tknr
        integer :: yy2, ccy, kcy, nang, x2, kk2, ll2, ss
        hx = 1.0d-3; nang = maxval(basis%am) + 2
        allocate(dapa(nbf,nbf), dapb(nbf,nbf), dedp(3,natom), dedm(3,natom))
        allocate(mopa(nbf,nbf), mopb(nbf,nbf), frap(nbf2), frbp(nbf2), fram(nbf2), frbm(nbf2))
        allocate(dFxc(nbf,nbf), dHse(ncart,ncart), dHt3(ncart,ncart))
        dHse = 0.0_dp; dHt3 = 0.0_dp
        call dft_initialize(infos, basis, mg); call dftclean(infos)   ! warm-up
        do yy2 = 1, ncart
          ccy = mod(yy2-1,3)+1; kcy = (yy2-1)/3+1
          basis%atoms%xyz(ccy,kcy) = basis%atoms%xyz(ccy,kcy) + hx
          call basis%init_shell_centers()
          call dft_initialize(infos, basis, mg)
          dapa = sp(1)%p + hx*sp(1)%dPx(:,:,yy2); dapb = sp(2)%p + hx*sp(2)%dPx(:,:,yy2)
          dedp = 0.0_dp
          call derexc_blk(basis, mg, dapa, dapb, dedp, telr, tknr, nang, nbf, &
                          infos%dft%grid_density_cutoff, .true., infos)
          mopa = sp(1)%mo; mopa(:,1:nocca) = sp(1)%mo(:,1:nocca) + hx*sp(1)%dCx(:,:,yy2)
          mopb = sp(2)%mo; mopb(:,1:noccb) = sp(2)%mo(:,1:noccb) + hx*sp(2)%dCx(:,:,yy2)
          frap = 0.0_dp; frbp = 0.0_dp
          call dftexcor(basis, mg, infos%control%scftype, frap, frbp, mopa, mopb, &
                        nbf, nbf2, exr, telr, tknr, infos)
          call dftclean(infos)
          basis%atoms%xyz(ccy,kcy) = basis%atoms%xyz(ccy,kcy) - 2*hx
          call basis%init_shell_centers()
          call dft_initialize(infos, basis, mg)
          dapa = sp(1)%p - hx*sp(1)%dPx(:,:,yy2); dapb = sp(2)%p - hx*sp(2)%dPx(:,:,yy2)
          dedm = 0.0_dp
          call derexc_blk(basis, mg, dapa, dapb, dedm, telr, tknr, nang, nbf, &
                          infos%dft%grid_density_cutoff, .true., infos)
          mopa = sp(1)%mo; mopa(:,1:nocca) = sp(1)%mo(:,1:nocca) - hx*sp(1)%dCx(:,:,yy2)
          mopb = sp(2)%mo; mopb(:,1:noccb) = sp(2)%mo(:,1:noccb) - hx*sp(2)%dCx(:,:,yy2)
          fram = 0.0_dp; frbm = 0.0_dp
          call dftexcor(basis, mg, infos%control%scftype, fram, frbm, mopa, mopb, &
                        nbf, nbf2, exr, telr, tknr, infos)
          call dftclean(infos)
          basis%atoms%xyz(ccy,kcy) = basis%atoms%xyz(ccy,kcy) + hx
          call basis%init_shell_centers()
          dHse(:,yy2) = reshape((dedp - dedm)/(2*hx), [ncart])
          ! dHt3: -2 sum_s s1oo^s,x (dVxc^s,y + fxc[dP^s,y])_oo
          do ss = 1, 2
            if (ss == 1) then
              call unpack_from_packed((frap - fram)/(2*hx), dFxc, nbf)
            else
              call unpack_from_packed((frbp - frbm)/(2*hx), dFxc, nbf)
            end if
            allocate(tmpn(nbf,sp(ss)%nocc), dFoo(sp(ss)%nocc,sp(ss)%nocc))
            call dgemm('n','n', nbf, sp(ss)%nocc, nbf, 1.0_dp, dFxc, nbf, sp(ss)%mo, nbf, 0.0_dp, tmpn, nbf)
            call dgemm('t','n', sp(ss)%nocc, sp(ss)%nocc, nbf, 1.0_dp, sp(ss)%mo, nbf, tmpn, nbf, 0.0_dp, dFoo, sp(ss)%nocc)
            do x2 = 1, ncart
              do ll2 = 1, sp(ss)%nocc
                do kk2 = 1, sp(ss)%nocc
                  dHt3(x2,yy2) = dHt3(x2,yy2) - 1.0_dp*sp(ss)%s1oo(kk2,ll2,x2)*dFoo(kk2,ll2)
                end do
              end do
            end do
            deallocate(tmpn, dFoo)
          end do
        end do
        hess_native = hess_native + 0.5_dp*(dHse + transpose(dHse)) &
                                  + 0.5_dp*(dHt3 + transpose(dHt3))
        deallocate(dapa, dapb, dedp, dedm, mopa, mopb, frap, frbp, fram, frbm, dFxc, dHse, dHt3)
      end block
    end if

    ! nuclear repulsion
    call hess_nn(basis%atoms, basis%ecp_zn_num, hess_native)

    ! --- one-electron + Pulay second-derivative skeleton (fixed density) ------
    block
      use grd1, only: eijden, hess_ee_overlap, hess_ee_kinetic, hess_en
      use ecp_tool, only: add_ecphess
      real(dp), allocatable :: wlag(:), pden(:), hcc(:,:)
      allocate(wlag(nbf2), pden(nbf2), hcc(ncart,ncart), source=0.0_dp)
      call eijden(wlag, nbf, infos)                 ! open-shell Lagrangian W
      pden = dma + dmb                              ! total density (Pa + Pb)
      call hess_ee_overlap(basis, wlag, hess_native)
      call hess_ee_kinetic(basis, pden, hess_native)
      call hess_en(basis, basis%atoms%xyz, &
                   basis%atoms%zn - basis%ecp_zn_num, pden, hess_native, hess_cc=hcc)
      call add_ecphess(basis, basis%atoms%xyz, pden, hess_native) ! ECP skeleton (if any)
      deallocate(wlag, pden, hcc)
    end block

    ! --- two-electron (ERI) second-derivative skeleton (fixed density) --------
    block
      use grd2, only: grd2_hess_driver, grd2_compute_data_t
      use hf_gradient_mod, only: grd2_uhf_compute_data_t
      class(grd2_compute_data_t), allocatable :: gcomp
      gcomp = grd2_uhf_compute_data_t( da = dma, db = dmb, hfscale = hfscale, nbf = nbf )
      call gcomp%init()
      call grd2_hess_driver(infos, basis, hess_native, gcomp)
      call gcomp%clean()
    end block

    call infos%dat%reserve_data(OQP_hf_hessian, TA_TYPE_REAL64, ncart*ncart, (/ ncart, ncart /), &
      comment='Native OpenQP open-shell (UHF) HF analytic Hessian matrix')
    call tagarray_get_data(infos%dat, OQP_hf_hessian, hess_store)
    hess_store = hess_native
    write(iw,'(A)') 'PyOQP: Native OpenQP open-shell (UHF) HF Hessian matrix stored'

    deallocate(ptot, dSa, dTa, dVa, sflat, hflat, bvec, uvec, scr, tmp, SxMO, hxMO, &
               probe, gx, d0a, d0b, gfull, Gd0, dpck, fpck, &
               A2, tGP, Mi, hresp, hess_native)
  end subroutine hf_hessian_uhf

!###############################################################################

  subroutine hf_hessian_rohf(infos)
    ! Native open-shell (ROHF) analytic HF Hessian (HF only).
    !
    ! ROHF uses a SINGLE MO set with a docc/socc/virt partition, so the orbital
    ! response is solved over the ROHF rotation space (cphf_solve_rohf) rather
    ! than the UHF spin blocks.  The ROHF energy has the same functional form as
    ! UHF in terms of (Pa, Pb), so the Hessian decomposes identically into
    !   H = E_nn'' + skeleton(1e total density + open-shell W, 2e via grd2_uhf)
    !       + response(orbital relaxation),
    ! where the skeleton + nuclear repulsion are exactly hess_skel_open.
    !
    ! The orbital-relaxation response is evaluated SEMI-NUMERICALLY, reusing the
    ! validated analytic open-shell gradient: with the relaxed orbital derivative
    ! dC^b (from the ROHF CPHF amplitudes) the response is the central finite
    ! difference, AT FIXED GEOMETRY, of the density/Lagrangian-dependent gradient
    ! along the orbital path C_occ +/- h dC^b:
    !   H^resp(:,b) = [ g(C + h dC^b) - g(C - h dC^b) ] / 2h ,
    !   g(C') = grad_ee_overlap(W') + grad_ee_kinetic(P') + grad_en(P')
    !           + grad_2e(Pa', Pb') ,  W' = -(Pa' Fa' Pa' + Pb' Fb' Pb') ,
    ! with Fa'/Fb' rebuilt from the perturbed densities (Hcore + fock_jk).  This
    ! captures BOTH the relaxed-density and the energy-weighted (W) response
    ! through the gradient's own W build (eijden convention), so no ROHF-specific
    ! Lagrangian-derivative algebra is required.  The CPHF right-hand side is the
    ! non-canonical Pulay form (orbital energies replaced by the full Fock occ-occ
    ! blocks), reducing to the validated UHF RHS in the canonical limit.
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_DM_B, &
      OQP_VEC_MO_A, OQP_FOCK_A, OQP_FOCK_B, OQP_Hcore, OQP_hf_hessian, TA_TYPE_REAL64
    use mathlib, only: unpack_matrix, pack_matrix, orthogonal_transform_sym
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix, hess_nn, &
      grad_ee_overlap, grad_ee_kinetic, grad_en_hellman_feynman, grad_en_pulay
    use grd2, only: grd2_driver, grd2_compute_data_t
    use hf_gradient_mod, only: grd2_uhf_compute_data_t
    use fock_deriv_mod, only: fock_deriv_contract_os
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve_rohf, rohf_pack_trial, rohf_unpack_trial
    use io_constants, only: iw

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(dp), contiguous, pointer :: dma(:), dmb(:), mo(:,:), focka(:), fockb(:), hcore(:)
    real(dp), contiguous, pointer :: hess_store(:,:)
    real(dp), allocatable :: pa(:,:), pb(:,:), ptot(:,:)
    real(dp), allocatable :: dSa(:,:,:,:), dTa(:,:,:,:), dVa(:,:,:,:)
    real(dp), allocatable :: faMO(:,:), fbMO(:,:)
    real(dp), allocatable :: scr(:,:), tmp(:,:), SxMO(:,:), hxMO(:,:), probe(:,:)
    real(dp), allocatable :: ga2e(:,:,:), gb2e(:,:,:)
    real(dp), allocatable :: d0a(:,:), d0b(:,:), dpck(:,:), fpck(:,:), gfull(:,:), Gd0(:,:)
    real(dp), allocatable :: ba(:,:), bb(:,:), bvec(:,:), uvec(:,:)
    real(dp), allocatable :: xa(:,:), xb(:,:), dCa(:,:), dCb(:,:), gp(:,:), gm(:,:)
    real(dp), allocatable :: zneff(:), hess_native(:,:), hresp(:,:)
    real(dp), allocatable :: faop(:), fbop(:)
    integer, allocatable :: iecp_atom(:)
    real(dp) :: hfscale, hstep, gx(3, size(infos%atoms%xyz,2))
    integer :: nbf, nbf2, natom, ncart, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: i, j, a, icart, kc, cc, x, mu, nu, ie, nec

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2)
    ncart = 3*natom
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira
    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale
    hstep = 1.0d-3

    write(iw,'(/,A)') 'PyOQP: Native OpenQP open-shell (ROHF) HF Hessian CPHF response prepass'
    write(iw,'(A,I6,A,I6,A,I6,A,I6,A,I6)') '  nbf=', nbf, ' nocca=', nocca, &
      ' noccb=', noccb, ' rhs=', ncart, ' rotdim=', ltot
    write(iw,'(A)') '  Storing native OpenQP open-shell (ROHF) HF analytic Hessian in OQP::hf_hessian.'

    if (ncart <= 0 .or. ltot <= 0) then
      write(iw,'(A)') '  ROHF CPHF prepass skipped: empty rotation/nuclear space.'
      return
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dma)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmb)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, focka)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fockb)
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore)

    allocate(pa(nbf,nbf), pb(nbf,nbf), ptot(nbf,nbf))
    call unpack_matrix(dma, pa); call unpack_matrix(dmb, pb); ptot = pa + pb
    allocate(zneff(natom)); zneff = basis%atoms%zn - basis%ecp_zn_num

    ! Map each atom to its ECP-centre index in ecp_coord (which is sized
    ! 3*num_ecps, i.e. one (x,y,z) triple per ECP centre, NOT per atom).  The
    ! semi-numerical resp_grad displaces atoms one Cartesian at a time and must
    ! move the matching ECP centre in lockstep; iecp_atom(kc)=0 means atom kc
    ! carries no ECP (its centre must not be touched).
    allocate(iecp_atom(natom)); iecp_atom = 0
    if (basis%ecp_params%is_ecp) then
      nec = size(basis%ecp_params%n_expo)
      do ie = 1, nec
        do i = 1, natom
          if (all(abs(basis%ecp_params%ecp_coord(3*(ie-1)+1:3*ie) &
                      - basis%atoms%xyz(:,i)) < 1.0e-6_dp)) then
            iecp_atom(i) = ie
            exit
          end if
        end do
      end do
    end if

    ! derivative integrals (normalized into the bfnrm/MO convention)
    allocate(dSa(nbf,nbf,3,natom), dTa(nbf,nbf,3,natom), dVa(nbf,nbf,3,natom))
    call der_overlap_matrix(basis, dSa)
    call der_kinetic_matrix(basis, dTa)
    call der_nucattr_matrix(basis, basis%atoms%xyz, &
                            basis%atoms%zn - basis%ecp_zn_num, dVa)  ! ECP-screened point charge
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

    ! ECP first-derivative integrals -> core-Hamiltonian derivative dHcore/dR,
    ! feeding the non-canonical CPHF RHS (hxMO below).  The ECP skeleton + response
    ! is then completed by add_ecpder inside resp_grad (semi-numerical).  Already
    ! normalized, so added after the bfnrm scaling.  No-op for non-ECP bases.
    block
      use ecp_tool, only: ecp_deriv_ints
      real(dp), allocatable :: dVecp(:,:,:,:)
      allocate(dVecp(nbf,nbf,3,natom))
      call ecp_deriv_ints(basis, basis%atoms%xyz, dVecp)
      dVa = dVa + dVecp
      deallocate(dVecp)
    end block

    ! occ-occ Fock blocks (MO) of the converged spin Fock matrices (non-canonical)
    allocate(scr(nbf,nbf), tmp(nbf,nbf), SxMO(nbf,nbf), hxMO(nbf,nbf))
    allocate(faMO(nbf,nbf), fbMO(nbf,nbf))
    call unpack_matrix(focka, scr); call mo_transform(mo, scr, nbf, tmp, hxMO, faMO)
    call unpack_matrix(fockb, scr); call mo_transform(mo, scr, nbf, tmp, hxMO, fbMO)

    ! 2e response-Fock skeleton  G^{s,x}[P]_ai  for all coordinates (per spin)
    allocate(ga2e(nvira,nocca,ncart), gb2e(nvirb,noccb,ncart), source=0.0_dp)
    allocate(probe(nbf,nbf))
    do a = 1, nvira
      do i = 1, nocca
        do mu = 1, nbf
          do nu = 1, nbf
            probe(mu,nu) = 0.5_dp*( mo(mu,nocca+a)*mo(nu,i) + mo(mu,i)*mo(nu,nocca+a) )
          end do
        end do
        gx = 0.0_dp
        call fock_deriv_contract_os(infos, basis, ptot, pa, probe, hfscale, gx)
        ga2e(a,i,:) = reshape(gx, [ncart])
      end do
    end do
    do a = 1, nvirb
      do i = 1, noccb
        do mu = 1, nbf
          do nu = 1, nbf
            probe(mu,nu) = 0.5_dp*( mo(mu,noccb+a)*mo(nu,i) + mo(mu,i)*mo(nu,noccb+a) )
          end do
        end do
        gx = 0.0_dp
        call fock_deriv_contract_os(infos, basis, ptot, pb, probe, hfscale, gx)
        gb2e(a,i,:) = reshape(gx, [ncart])
      end do
    end do

    ! ===== CPHF right-hand sides (non-canonical Pulay form), packed =====
    allocate(d0a(nbf,nbf), d0b(nbf,nbf), gfull(nbf,nbf), Gd0(nbf,nbf))
    allocate(dpck(nbf2,2), fpck(nbf2,2))
    allocate(ba(nvira,nocca), bb(nvirb,noccb))
    allocate(bvec(ltot,ncart), uvec(ltot,ncart), source=0.0_dp)
    icart = 0
    do kc = 1, natom
      do cc = 1, 3
        icart = icart + 1
        call mo_transform(mo, dSa(:,:,cc,kc), nbf, scr, tmp, SxMO)
        call mo_transform(mo, dTa(:,:,cc,kc)+dVa(:,:,cc,kc), nbf, scr, tmp, hxMO)

        ! reorthonormalization densities d0^s = -sum_ij S^x_ij C_i C_j (per spin occ)
        d0a = 0.0_dp; d0b = 0.0_dp
        do i = 1, nocca
          do j = 1, nocca
            do mu = 1, nbf
              do nu = 1, nbf
                d0a(mu,nu) = d0a(mu,nu) - SxMO(i,j)*mo(mu,i)*mo(nu,j)
              end do
            end do
          end do
        end do
        do i = 1, noccb
          do j = 1, noccb
            do mu = 1, nbf
              do nu = 1, nbf
                d0b(mu,nu) = d0b(mu,nu) - SxMO(i,j)*mo(mu,i)*mo(nu,j)
              end do
            end do
          end do
        end do
        call pack_matrix(d0a, dpck(:,1)); call pack_matrix(d0b, dpck(:,2))
        fpck = 0.0_dp
        call fock_jk(basis, d=dpck, f=fpck, scale_exch=hfscale, infos=infos)

        ! Non-canonical Pulay RHS.  The reorthonormalization Fock-coupling is the
        ! occupied-projected anticommutator of S^x and the spin Fock:
        !   B^s_ai = -(h^x + G2e + G[d0])_ai
        !            + sum_{j in occ} ( S^x_aj F^s_ji + F^s_aj S^x_ji ) .
        ! The first sum is the usual eps_i S^x_ai in the canonical (diagonal-Fock)
        ! limit; the second vanishes there (F^s_aj is a vir-occ Fock element) and
        ! supplies the non-canonical correction needed for the socc rotations.
        call unpack_from_packed(fpck(:,1), gfull, nbf)
        call mo_transform(mo, gfull, nbf, scr, tmp, Gd0)
        do i = 1, nocca
          do a = 1, nvira
            ba(a,i) = -(hxMO(i,nocca+a) + ga2e(a,i,icart) + Gd0(i,nocca+a)) &
                    + dot_product(SxMO(nocca+a,1:nocca), faMO(1:nocca,i)) &
                    + dot_product(faMO(nocca+a,1:nocca), SxMO(1:nocca,i))
          end do
        end do
        ! beta block
        call unpack_from_packed(fpck(:,2), gfull, nbf)
        call mo_transform(mo, gfull, nbf, scr, tmp, Gd0)
        do i = 1, noccb
          do a = 1, nvirb
            bb(a,i) = -(hxMO(i,noccb+a) + gb2e(a,i,icart) + Gd0(i,noccb+a)) &
                    + dot_product(SxMO(noccb+a,1:noccb), fbMO(1:noccb,i)) &
                    + dot_product(fbMO(noccb+a,1:noccb), SxMO(1:noccb,i))
          end do
        end do

        ! --- XC contribution to the CPKS right-hand side (ROKS only) -----------
        ! Central FD of the spin XC Fock (open-shell dftexcor) along R +/- h AND
        ! occupied MOs reorthonormalized by dmoR^s = -1/2 sum_j C_j S^x_ji; the XC
        ! skeleton dVxc/dR + f_xc[d0], subtracted from B (which carries -F0x).
        if (infos%control%hamilton >= 20) then
          block
            use dft, only: dft_initialize, dftclean, dftexcor
            use mod_dft_molgrid, only: dft_grid_t
            type(dft_grid_t) :: mgr
            real(dp), allocatable :: dmoa(:,:), dmob(:,:), mopa(:,:), mopb(:,:)
            real(dp), allocatable :: frap(:), frbp(:), fram(:), frbm(:), dvx(:,:), hxc(:,:)
            real(dp) :: hxr, exr, telr, tknr
            integer :: ir, jr
            allocate(dmoa(nbf,nocca), dmob(nbf,noccb), mopa(nbf,nbf), mopb(nbf,nbf))
            allocate(frap(nbf2), frbp(nbf2), fram(nbf2), frbm(nbf2), dvx(nbf,nbf), hxc(nbf,nbf))
            hxr = 1.0d-3
            dmoa = 0.0_dp
            do ir = 1, nocca
              do jr = 1, nocca
                dmoa(:,ir) = dmoa(:,ir) - 0.5_dp*SxMO(jr,ir)*mo(:,jr)
              end do
            end do
            dmob = 0.0_dp
            do ir = 1, noccb
              do jr = 1, noccb
                dmob(:,ir) = dmob(:,ir) - 0.5_dp*SxMO(jr,ir)*mo(:,jr)
              end do
            end do
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) + hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mopa = mo; mopa(:,1:nocca) = mo(:,1:nocca) + hxr*dmoa
            mopb = mo; mopb(:,1:noccb) = mo(:,1:noccb) + hxr*dmob
            frap = 0.0_dp; frbp = 0.0_dp
            call dftexcor(basis, mgr, infos%control%scftype, frap, frbp, mopa, mopb, &
                          nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) - 2*hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mopa = mo; mopa(:,1:nocca) = mo(:,1:nocca) - hxr*dmoa
            mopb = mo; mopb(:,1:noccb) = mo(:,1:noccb) - hxr*dmob
            fram = 0.0_dp; frbm = 0.0_dp
            call dftexcor(basis, mgr, infos%control%scftype, fram, frbm, mopa, mopb, &
                          nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) + hxr
            call basis%init_shell_centers()
            call unpack_from_packed((frap - fram)/(2*hxr), dvx, nbf)
            call mo_transform(mo, dvx, nbf, scr, tmp, hxc)
            do i = 1, nocca
              do a = 1, nvira
                ba(a,i) = ba(a,i) - hxc(i,nocca+a)
              end do
            end do
            call unpack_from_packed((frbp - frbm)/(2*hxr), dvx, nbf)
            call mo_transform(mo, dvx, nbf, scr, tmp, hxc)
            do i = 1, noccb
              do a = 1, nvirb
                bb(a,i) = bb(a,i) - hxc(i,noccb+a)
              end do
            end do
            deallocate(dmoa, dmob, mopa, mopb, frap, frbp, fram, frbm, dvx, hxc)
          end block
        end if

        call rohf_pack_trial(bvec(:,icart), ba, bb, nbf, nocca, noccb)
      end do
    end do

    call cphf_solve_rohf(infos, ncart, bvec, uvec)

    ! ===== semi-numerical orbital-relaxation response =====
    ! Build the relaxed alpha/beta orbital derivatives independently, UHF-style:
    !   dCa_i = sum_a C^{vir_a}_a xa(a,i) - 1/2 sum_{j in docc+socc} S^x_ji C_j
    !   dCb_i = sum_a C^{vir_b}_a xb(a,i) - 1/2 sum_{j in docc}      S^x_ji C_j
    ! The socc-docc rotation lives in xb (socc is beta-virtual), so it relaxes Pb
    ! and leaves Pa invariant (it is an alpha occ-occ rotation) -- exactly the
    ! ROHF physics, with no socc-docc cross term needed in dCa.
    allocate(xa(nvira,nocca), xb(nvirb,noccb), dCa(nbf,nocca), dCb(nbf,noccb))
    allocate(gp(3,natom), gm(3,natom), hresp(ncart,ncart), source=0.0_dp)
    allocate(faop(nbf2), fbop(nbf2))
    if (infos%control%hamilton >= 20) then   ! flush stale grid state from the CPHF solver
      block
        use dft, only: dft_initialize, dftclean
        use mod_dft_molgrid, only: dft_grid_t
        type(dft_grid_t) :: mgw
        call dft_initialize(infos, basis, mgw); call dftclean(infos)
      end block
    end if
    do x = 1, ncart
      cc = mod(x-1,3)+1; kc = (x-1)/3+1
      call rohf_unpack_trial(uvec(:,x), xa, xb, nbf, nocca, noccb)
      call mo_transform(mo, dSa(:,:,cc,kc), nbf, scr, tmp, SxMO)

      dCa = 0.0_dp
      do i = 1, nocca
        do a = 1, nvira
          dCa(:,i) = dCa(:,i) + mo(:,nocca+a)*xa(a,i)
        end do
        do j = 1, nocca
          dCa(:,i) = dCa(:,i) - 0.5_dp*SxMO(j,i)*mo(:,j)
        end do
      end do
      dCb = 0.0_dp
      do i = 1, noccb
        do a = 1, nvirb
          dCb(:,i) = dCb(:,i) + mo(:,noccb+a)*xb(a,i)
        end do
        do j = 1, noccb
          dCb(:,i) = dCb(:,i) - 0.5_dp*SxMO(j,i)*mo(:,j)
        end do
      end do

      call resp_grad( 1.0_dp, gp)
      call resp_grad(-1.0_dp, gm)
      hresp(:,x) = reshape((gp - gm)/(2.0_dp*hstep), [ncart])
    end do

    ! The central difference of the ELECTRONIC gradient over geometry AND the
    ! relaxed orbital path already contains the full electronic Hessian (skeleton
    ! + orbital-relaxation response); only the (orbital-independent) nuclear
    ! repulsion second derivative is added analytically.
    allocate(hess_native(ncart,ncart))
    hess_native = 0.5_dp*(hresp + transpose(hresp))
    call hess_nn(basis%atoms, basis%ecp_zn_num, hess_native)

    call infos%dat%reserve_data(OQP_hf_hessian, TA_TYPE_REAL64, ncart*ncart, (/ ncart, ncart /), &
      comment='Native OpenQP open-shell (ROHF) HF analytic Hessian matrix')
    call tagarray_get_data(infos%dat, OQP_hf_hessian, hess_store)
    hess_store = hess_native
    write(iw,'(A)') 'PyOQP: Native OpenQP open-shell (ROHF) HF Hessian matrix stored'

    deallocate(pa, pb, ptot, dSa, dTa, dVa, faMO, fbMO, scr, tmp, &
               SxMO, hxMO, probe, ga2e, gb2e, d0a, d0b, dpck, fpck, gfull, Gd0, &
               ba, bb, bvec, uvec, xa, xb, dCa, dCb, gp, gm, hresp, zneff, &
               hess_native, faop, fbop, iecp_atom)

  contains

    !> Electronic gradient (1e + 2e + Pulay-W; NO nuclear repulsion) at the
    !> geometry displaced by sgn*hstep in coordinate (cc,kc) AND the alpha-occ
    !> MOs displaced by sgn*hstep*dC (host-associated cc,kc,dC,hstep).  Central
    !> differencing over sgn therefore captures the electronic skeleton AND the
    !> orbital-relaxation response together: the one-electron Hamiltonian, all
    !> gradient integrals, the densities Pa'/Pb' and the energy-weighted density
    !> W' = -(Pa' Fa' Pa' + Pb' Fb' Pb') (eijden convention, Fock rebuilt as
    !> Hcore' + fock_jk) are all evaluated at the displaced point, so no
    !> ROHF-specific Lagrangian-derivative algebra is needed.
    subroutine resp_grad(sgn, gout)
      use int1, only: omp_hst
      real(dp), intent(in) :: sgn
      real(dp), intent(out) :: gout(:,:)
      real(dp), allocatable :: cocc(:,:), pap(:,:), pbp(:,:)
      real(dp), allocatable, target :: paP_tri(:), pbP_tri(:)
      real(dp), allocatable :: ptP_tri(:), wlag(:), ta(:), hc(:), sm(:), tm(:)
      real(dp) :: tol
      integer :: ii, ij
      class(grd2_compute_data_t), allocatable :: gc

      allocate(cocc(nbf,nocca), pap(nbf,nbf), pbp(nbf,nbf))
      allocate(paP_tri(nbf2), pbP_tri(nbf2), ptP_tri(nbf2), wlag(nbf2), ta(nbf2))
      allocate(hc(nbf2), sm(nbf2), tm(nbf2))

      ! displace geometry and rebuild the one-electron Hamiltonian there.  The ECP
      ! center (ecp_coord) is a separate array from atoms%xyz, so it must be moved
      ! in lockstep or the displaced add_ecpint/add_ecpder would see the basis and
      ! the ECP at mismatched centers (catastrophic for the ECP atom).
      basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) + sgn*hstep
      if (iecp_atom(kc) > 0) &
        basis%ecp_params%ecp_coord(3*(iecp_atom(kc)-1)+cc) = &
          basis%ecp_params%ecp_coord(3*(iecp_atom(kc)-1)+cc) + sgn*hstep
      call basis%init_shell_centers()
      tol = log(10.0d0)*20.0_dp
      call omp_hst(basis, basis%atoms%xyz, basis%atoms%zn - basis%ecp_zn_num, &
                   hc, sm, tm, logtol=tol, comm=infos%mpiinfo%comm, usempi=infos%mpiinfo%usempi)
      ! NB: the ECP one-electron potential is deliberately NOT added to hc here.
      ! The full ECP gradient (operator + basis-centre/Pulay derivatives) is the
      ! analytic add_ecpder below; folding the ECP into the spin Fock used to build
      ! the energy-weighted density W' would double-count its Pulay contribution
      ! (verified: doing so gives ~1.5e-2 vs the numerical Hessian, omitting it
      ! gives ~2e-5).

      ! relaxed orbitals -> perturbed densities and spin Fock matrices
      cocc(:,1:nocca) = mo(:,1:nocca) + sgn*hstep*dCa
      call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, cocc, nbf, cocc, nbf, 0.0_dp, pap, nbf)
      cocc(:,1:noccb) = mo(:,1:noccb) + sgn*hstep*dCb
      call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, cocc, nbf, cocc, nbf, 0.0_dp, pbp, nbf)
      call pack_matrix(pap, paP_tri); call pack_matrix(pbp, pbP_tri)
      ptP_tri = paP_tri + pbP_tri
      dpck(:,1) = paP_tri; dpck(:,2) = pbP_tri
      fpck = 0.0_dp
      call fock_jk(basis, d=dpck, f=fpck, scale_exch=hfscale, infos=infos)
      faop = hc + fpck(:,1); fbop = hc + fpck(:,2)

      gout = 0.0_dp
      ! DFT (ROKS): add the XC potential to the spin Focks (so W' is the full KS
      ! energy-weighted density) and the explicit open-shell XC gradient to gout.
      ! Both are evaluated at the displaced geometry with the relaxed orbitals, so
      ! the geometry+orbital FD gives the full KS Hessian (Pulay/W XC + explicit XC)
      ! with no separate analytic XC term.
      if (infos%control%hamilton >= 20) then
        block
          use dft, only: dft_initialize, dftclean, dftexcor
          use mod_dft_gridint_grad, only: derexc_blk
          use mod_dft_molgrid, only: dft_grid_t
          type(dft_grid_t) :: mg
          real(dp), allocatable :: mopa(:,:), mopb(:,:), fra(:), frb(:), dedft(:,:)
          real(dp) :: exr, telr, tknr
          integer :: nang
          allocate(mopa(nbf,nbf), mopb(nbf,nbf), fra(nbf2), frb(nbf2), dedft(3,natom))
          nang = maxval(basis%am) + 2
          call dft_initialize(infos, basis, mg)
          mopa = mo; mopa(:,1:nocca) = mo(:,1:nocca) + sgn*hstep*dCa
          mopb = mo; mopb(:,1:noccb) = mo(:,1:noccb) + sgn*hstep*dCb
          fra = 0.0_dp; frb = 0.0_dp
          call dftexcor(basis, mg, infos%control%scftype, fra, frb, mopa, mopb, &
                        nbf, nbf2, exr, telr, tknr, infos)
          faop = faop + fra; fbop = fbop + frb
          dedft = 0.0_dp
          call derexc_blk(basis, mg, pap, pbp, dedft, telr, tknr, nang, nbf, &
                          infos%dft%grid_density_cutoff, .true., infos)
          call dftclean(infos)
          gout = gout + dedft
          deallocate(mopa, mopb, fra, frb, dedft)
        end block
      end if

      call orthogonal_transform_sym(nbf, nbf, faop, pap, nbf, ta)
      call orthogonal_transform_sym(nbf, nbf, fbop, pbp, nbf, wlag)
      wlag = -wlag - ta
      ij = 0
      do ii = 1, nbf
        ij = ij + ii
        wlag(ij) = 0.5_dp*wlag(ij)
      end do

      call grad_ee_overlap(basis, wlag, gout)
      call grad_ee_kinetic(basis, ptP_tri, gout)
      call grad_en_hellman_feynman(basis, basis%atoms%xyz, zneff, ptP_tri, gout)
      call grad_en_pulay(basis, basis%atoms%xyz, zneff, ptP_tri, gout)
      ! ECP gradient at the displaced geometry/density: central FD over the
      ! geometry+orbital path then yields BOTH the ECP skeleton second derivative
      ! and the ECP orbital-relaxation response.  No-op for non-ECP bases.
      block
        use ecp_tool, only: add_ecpder
        call add_ecpder(basis, basis%atoms%xyz, ptP_tri, gout)
      end block
      gc = grd2_uhf_compute_data_t( da = paP_tri, db = pbP_tri, hfscale = hfscale, nbf = nbf )
      call gc%init()
      call grd2_driver(infos, basis, gout, gc)
      call gc%clean()

      ! restore geometry (and the ECP center moved above)
      basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) - sgn*hstep
      if (iecp_atom(kc) > 0) &
        basis%ecp_params%ecp_coord(3*(iecp_atom(kc)-1)+cc) = &
          basis%ecp_params%ecp_coord(3*(iecp_atom(kc)-1)+cc) - sgn*hstep
      call basis%init_shell_centers()

      deallocate(cocc, pap, pbp, paP_tri, pbP_tri, ptP_tri, wlag, ta, hc, sm, tm)
    end subroutine resp_grad

  end subroutine hf_hessian_rohf

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
