module hf_hessian_mod

  implicit none

  character(len=*), parameter :: module_name = "hf_hessian_mod"

  ! g2e(ia,x) = 1/2 Tr[probe_ia G^x[P]] (fock_deriv_contract convention) in
  ! terms of the channel operator of eri_derivative_operator_mo.
  real(kind=8), parameter :: G2E_OPERATOR_SCALE = 0.125d0
  ! Open shell: Tr[M (J^x[Ptot] - c_x K^x[P^s])] (fock_deriv_contract_os) in
  ! terms of the Coulomb-only and exchange-only channel operators.
  real(kind=8), parameter :: OS_J_SCALE = 0.25d0, OS_K_SCALE = 0.5d0

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
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, OQP_E_MO_A, &
      OQP_hf_hessian, TA_TYPE_REAL64
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix, hess_nn
    use fock_deriv_mod, only: fock_deriv_contract
    use tdhf_hessian_z_rhs_mod, only: eri_derivative_operator_mo
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve
    use io_constants, only: iw
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dmat_a(:), mo_a(:,:), eps(:)
    real(kind=dp), allocatable :: pfull(:,:), probe(:,:), gx(:,:), g2e(:,:), gop(:,:,:)
    ! XC response-assembly pieces produced on the shared displaced grids of the
    ! IR/Raman XC terms (DFT with intensities only)
    real(kind=dp), allocatable :: xc_dhse(:,:), xc_dfxc(:,:,:)
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
    ! contract the ECP skeleton d^2 V_ECP/dR^2 analytically (add_ecphess, ecp_raw_ints
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
    if (infos%control%verbose >= 2) then
      write(iw,'(/,A)') 'PyOQP: Native OpenQP HF/DFT Hessian CPHF response prepass'
      write(iw,'(A,I6,A,I6,A,I6,A,I6)') '  nbf=', nbf, ' nocc=', nocc, ' nvir=', nvir, ' rhs=', ncart
      write(iw,'(A)') '  Storing native OpenQP HF/DFT analytic Hessian matrix in OQP::hf_hessian.'
    end if

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
    call hess_tick('start')
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
    ! response exactly as point-charge nuclear attraction does.  ecp_tool returns
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

    ! 2e response-Fock skeleton  G^x[P]_ia  for ALL 3N coordinates.  The occ-vir
    ! probe is geometry-independent and one derivative-Fock contraction returns
    ! every Cartesian component, so evaluate it once per occ-vir pair here
    ! instead of once per pair AND per coordinate (an ncart-fold redundant grd2
    ! sweep).  Same scheme as hf_hessian_uhf.
    call hess_tick('1e derivative integrals')
    allocate(g2e(nocc*nvir,ncart), source=0.0_dp)
    if (.not. infos%dft%cam_flag .and. .not. infos%mpiinfo%usempi) then
      ! One blocked derivative-ERI traversal assembles C^T G^x[P] C for all
      ! 3N coordinates at once, instead of one full traversal per occ-vir
      ! pair (nocc*nvir traversals).  Range-separated functionals (two
      ! attenuated passes) and MPI keep the per-pair contraction below.
      ! gop keeps the AO operator O^x (fock_deriv_contract(P,M) =
      ! G2E_OPERATOR_SCALE*Tr[M O^x]) for the 2e traces after the CPHF solve.
      block
        real(kind=dp), allocatable :: eye(:,:), t(:,:), gov(:,:)
        allocate(eye(nbf,nbf), source=0.0_dp)
        do mu = 1, nbf
          eye(mu,mu) = 1.0_dp
        end do
        allocate(gop(nbf,nbf,ncart), t(nbf,nvir), gov(nocc,nvir))
        call eri_derivative_operator_mo(infos, eye, pfull, 1, hfscale, gop)
        do icart = 1, ncart
          call dgemm('n','n',nbf,nvir,nbf,1.0_dp,gop(:,:,icart),nbf,mo_a(:,nocc+1:),nbf,0.0_dp,t,nbf)
          call dgemm('t','n',nocc,nvir,nbf,G2E_OPERATOR_SCALE,mo_a(:,1:nocc),nbf,t,nbf,0.0_dp,gov,nocc)
          do a = 1, nvir
            do i = 1, nocc
              g2e((a-1)*nocc+i,icart) = gov(i,a)
            end do
          end do
        end do
      end block
      call check_g2e_operator(infos, basis, mo_a, pfull, hfscale, nocc, g2e)
    else
      do a = 1, nvir
        do i = 1, nocc
          do mu = 1, nbf
            do nu = 1, nbf
              probe(mu,nu) = 0.5_dp*( mo_a(mu,nocc+a)*mo_a(nu,i) + mo_a(mu,i)*mo_a(nu,nocc+a) )
            end do
          end do
          call fock_deriv_contract(infos, basis, pfull, probe, hfscale, gx)
          g2e((a-1)*nocc+i,:) = reshape(gx, [ncart])
        end do
      end do
    end if

    call hess_tick('2e response skeleton g2e')
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
            F0x(i,nocc+a) = hx(i,nocc+a) + 2.0_dp*g2e((a-1)*nocc+i,icart)
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
            use mod_dft, only: dft_initialize, dftclean, dftexcor
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

    call hess_tick('CPHF right-hand sides')
    call cphf_solve(infos, ncart, bvec, uvec)
    call hess_tick('CPHF solve')

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

      ! analytic dipole derivatives (IR intensities) from the relaxed dP^y
      if (hf_hess_properties_wanted(infos)) then
      block
        real(dp), allocatable :: dipf(:,:,:), dmu(:,:)
        call hf_dipder_init(infos, pfull, dipf, dmu)
        do yy = 1, ncart
          call hf_dipder_add_response(dipf, dPx(:,:,yy), dmu(:,yy))
        end do
        call hf_dipder_store(infos, dmu)
      end block

      ! analytic polarizability derivatives (Raman activities); the ECP enters
      ! only through h^x (ecp_deriv_ints is folded into dVa above)
      block
        use oqp_tagarray_driver, only: OQP_hf_polarizability_derivatives
        real(dp), allocatable :: dpol(:,:,:)
        real(dp), contiguous, pointer :: pstore(:,:,:)
        allocate(dpol(3,3,ncart))
        if (infos%control%hamilton == 20) then
          allocate(xc_dhse(ncart,ncart), xc_dfxc(nbf,nbf,ncart))
          call hf_polder_rhf(infos, mo_a, eps, pfull, sflat, hflat, uvec, dPx, &
                             nocc, nvir, hfscale, dpol, xc_dhse, xc_dfxc)
        else
          call hf_polder_rhf(infos, mo_a, eps, pfull, sflat, hflat, uvec, dPx, &
                             nocc, nvir, hfscale, dpol)
        end if
        call infos%dat%alloc_or_die(OQP_hf_polarizability_derivatives, (/ 3, 3, ncart /), pstore, &
          description='Analytic nuclear derivatives of the static polarizability (a.u.), (3,3,3N)')
        pstore = dpol
        deallocate(dpol)
      end block
      end if

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
      call hess_tick('dP/dR, response Fock, IR/Raman')
      allocate(gxy(3,natom), A2(ncart,ncart), tGP(ncart,ncart), Mi(nbf,nbf), source=0.0_dp)
      if (allocated(gop)) then
        ! Both traces from the stored AO operator: two GEMMs, no ERI pass.
        block
          real(kind=dp), allocatable :: mall(:,:,:)
          allocate(mall(nbf,nbf,ncart))
          do x = 1, ncart
            call dgemm('n','n',nbf,nocc,nocc,1.0_dp,cocc,nbf,s1oo(:,:,x),nocc,0.0_dp,tmpno,nbf)
            call dgemm('n','t',nbf,nbf,nocc,1.0_dp,tmpno,nbf,cocc,nbf,0.0_dp,mall(:,:,x),nbf)
          end do
          call dgemm('t','n',ncart,ncart,nbf*nbf,2.0_dp*G2E_OPERATOR_SCALE,gop,nbf*nbf, &
                     dPx,nbf*nbf,0.0_dp,A2,ncart)
          call dgemm('t','n',ncart,ncart,nbf*nbf,2.0_dp*G2E_OPERATOR_SCALE,mall,nbf*nbf, &
                     gop,nbf*nbf,0.0_dp,tGP,ncart)
          call check_trace_operator(infos, basis, pfull, hfscale, dPx(:,:,1), A2(:,1), 'A2')
          call check_trace_operator(infos, basis, pfull, hfscale, mall(:,:,1), tGP(1,:), 'tGP')
        end block
        deallocate(gop)
      else
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
      end if

      call hess_tick('2e traces A2/tGP')
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
          use mod_dft, only: dft_initialize, dftclean, dftexcor
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
          if (allocated(xc_dhse)) then
            ! already evaluated on the IR/Raman displaced grids
            dHse = xc_dhse
            do yy2 = 1, ncart
              call dgemm('n','n',nbf,nocc,nbf,1.0_dp,xc_dfxc(:,:,yy2),nbf,mo_a,nbf,0.0_dp,tmpn,nbf)
              call dgemm('t','n',nocc,nocc,nbf,1.0_dp,mo_a,nbf,tmpn,nbf,0.0_dp,dFoo,nocc)
              do x2 = 1, ncart
                do ll2 = 1, nocc
                  do kk2 = 1, nocc
                    dHt3(x2,yy2) = dHt3(x2,yy2) - 2.0_dp*s1oo(kk2,ll2,x2)*dFoo(kk2,ll2)
                  end do
                end do
              end do
            end do
            deallocate(xc_dhse, xc_dfxc)
          else
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
          end if
          hess_native = hess_native + 0.5_dp*(dHse + transpose(dHse)) &
                                    + 0.5_dp*(dHt3 + transpose(dHt3))
          deallocate(dap, dedp, dedm, mop, frp, frm, dFxc, dFoo, tmpn, dHse, dHt3)
        end block
      end if

      deallocate(sflat, hflat, dCx, dPx, Gdp, s1oo, hMOoo, GdpMOoo, moe1a, &
                 Mi, gxy, A2, tGP, hresp, s1, s2, bMO, dpp, gpp, gfl, cocc, tmpno)
    end block

    call hess_tick('response assembly + XC')
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
      type(grd2_rhf_compute_data_t) :: gcomp
      gcomp = grd2_rhf_compute_data_t( da = dmat_a, hfscale = hfscale, nbf = nbf )
      call gcomp%init()
      call gcomp%build_cart(basis)
      call grd2_hess_driver(infos, basis, hess_native, gcomp)
      call gcomp%clean()
    end block

    call hess_tick('skeleton second derivatives')
    call infos%dat%alloc_or_die(OQP_hf_hessian, (/ ncart, ncart /), hess_store, &
      description='Native OpenQP HF/DFT analytic Hessian matrix')
    hess_store = hess_native
    if (infos%control%verbose >= 2) then
      write(iw,'(A)') 'PyOQP: Native OpenQP HF/DFT Hessian matrix stored'
    end if
    close(iw)

    deallocate(pfull, dSa, dTa, dVa, scr, col, Sx, hx, F0x, Gd0, probe, gx, g2e, &
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
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, OQP_DM_A, OQP_DM_B, &
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
    real(dp), allocatable :: gsao(:,:,:,:)      ! AO G^{s,x} (operator path only)
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

    if (infos%control%verbose >= 2) then
      write(iw,'(/,A)') 'PyOQP: Native OpenQP open-shell (UHF) HF Hessian CPHF response prepass'
      write(iw,'(A,I6,A,I6,A,I6,A,I6,A,I6)') '  nbf=', nbf, ' nocca=', nocca, &
        ' noccb=', noccb, ' rhs=', ncart, ' ltot=', ltot
      write(iw,'(A)') '  Storing native OpenQP open-shell HF analytic Hessian in OQP::hf_hessian.'
    end if

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
    if (os_operator_path(infos)) then
      ! Three blocked derivative-ERI traversals (J^x[Ptot], K^x[Pa], K^x[Pb])
      ! instead of one per occ-vir pair and spin.
      call os_g2e_from_operators(infos, basis, ptot, sp(1)%p, sp(2)%p, hfscale, &
        sp(1)%mo, sp(1)%nocc, sp(2)%mo, sp(2)%nocc, sp(1)%g2e, sp(2)%g2e, gsao)
    else
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
    end if

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
            use mod_dft, only: dft_initialize, dftclean, dftexcor
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
            call dftexcor(basis, mgr, int(infos%control%scftype), frap, frbp, mopa, mopb, &
                          nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) - 2*hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mopa = sp(1)%mo; mopa(:,1:nocca) = sp(1)%mo(:,1:nocca) - hxr*dmoa
            mopb = sp(2)%mo; mopb(:,1:noccb) = sp(2)%mo(:,1:noccb) - hxr*dmob
            fram = 0.0_dp; frbm = 0.0_dp
            call dftexcor(basis, mgr, int(infos%control%scftype), fram, frbm, mopa, mopb, &
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
    if (allocated(gsao)) then
      ! Both traces from the stored AO operators: GEMMs, no ERI pass.
      block
        real(dp), allocatable :: mall(:,:,:), t(:,:)
        allocate(mall(nbf,nbf,ncart))
        do s = 1, 2
          allocate(t(nbf,sp(s)%nocc))
          do x = 1, ncart
            call dgemm('n','n',nbf,sp(s)%nocc,sp(s)%nocc,1.0_dp,sp(s)%mo,nbf, &
                       sp(s)%s1oo(:,:,x),sp(s)%nocc,0.0_dp,t,nbf)
            call dgemm('n','t',nbf,nbf,sp(s)%nocc,1.0_dp,t,nbf,sp(s)%mo,nbf, &
                       0.0_dp,mall(:,:,x),nbf)
          end do
          deallocate(t)
          call dgemm('t','n',ncart,ncart,nbf*nbf,1.0_dp,gsao(:,:,:,s),nbf*nbf, &
                     sp(s)%dPx,nbf*nbf,1.0_dp,A2,ncart)
          call dgemm('t','n',ncart,ncart,nbf*nbf,1.0_dp,mall,nbf*nbf, &
                     gsao(:,:,:,s),nbf*nbf,1.0_dp,tGP,ncart)
        end do
      end block
      deallocate(gsao)
    else
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
    end if

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
        use mod_dft, only: dft_initialize, dftclean, dftexcor
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
          call dftexcor(basis, mg, int(infos%control%scftype), frap, frbp, mopa, mopb, &
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
          call dftexcor(basis, mg, int(infos%control%scftype), fram, frbm, mopa, mopb, &
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
      type(grd2_uhf_compute_data_t) :: gcomp
      gcomp = grd2_uhf_compute_data_t( da = dma, db = dmb, hfscale = hfscale, nbf = nbf )
      call gcomp%init()
      call gcomp%build_cart(basis)
      call grd2_hess_driver(infos, basis, hess_native, gcomp)
      call gcomp%clean()
    end block

    ! analytic dipole derivatives (IR intensities) from the relaxed dP^y
    if (hf_hess_properties_wanted(infos)) then
    block
      real(dp), allocatable :: dipf(:,:,:), dmu(:,:)
      integer :: yd
      call hf_dipder_init(infos, ptot, dipf, dmu)
      do yd = 1, ncart
        call hf_dipder_add_response(dipf, sp(1)%dPx(:,:,yd) + sp(2)%dPx(:,:,yd), dmu(:,yd))
      end do
      call hf_dipder_store(infos, dmu)
    end block

    ! analytic polarizability derivatives (Raman activities); the ECP enters
    ! only through h^x (ecp_deriv_ints is folded into dVa above)
    block
      use oqp_tagarray_driver, only: OQP_hf_polarizability_derivatives
      real(dp), allocatable :: dpol(:,:,:)
      real(dp), contiguous, pointer :: pstore(:,:,:)
      allocate(dpol(3,3,ncart))
      call hf_polder_uhf(infos, sp(1)%mo, sp(2)%mo, sp(1)%eps, sp(2)%eps, sp(1)%p, sp(2)%p, &
                         sp(1)%nocc, sp(2)%nocc, sflat, hflat, uvec, sp(1)%dPx, sp(2)%dPx, &
                         hfscale, dpol)
      call infos%dat%alloc_or_die(OQP_hf_polarizability_derivatives, (/ 3, 3, ncart /), pstore, &
        description='Analytic nuclear derivatives of the static polarizability (a.u.), (3,3,3N)')
      pstore = dpol
      deallocate(dpol)
    end block
    end if

    call infos%dat%alloc_or_die(OQP_hf_hessian, (/ ncart, ncart /), hess_store, &
      description='Native OpenQP open-shell (UHF) HF analytic Hessian matrix')
    hess_store = hess_native
    if (infos%control%verbose >= 2) then
      write(iw,'(A)') 'PyOQP: Native OpenQP open-shell (UHF) HF Hessian matrix stored'
    end if

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
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, OQP_DM_A, OQP_DM_B, &
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
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(dp), contiguous, pointer :: dma(:), dmb(:), mo(:,:), focka(:), fockb(:), hcore(:)
    real(dp), contiguous, pointer :: hess_store(:,:)
    real(dp), allocatable :: pa(:,:), pb(:,:), ptot(:,:)
    real(dp), allocatable :: dipf(:,:,:), dmu(:,:), dptx(:,:)
    logical :: want_props
    real(dp), allocatable :: dSa(:,:,:,:), dTa(:,:,:,:), dVa(:,:,:,:)
    real(dp), allocatable :: faMO(:,:), fbMO(:,:)
    real(dp), allocatable :: scr(:,:), tmp(:,:), SxMO(:,:), hxMO(:,:), probe(:,:)
    real(dp), allocatable :: ga2e(:,:,:), gb2e(:,:,:)
    real(dp), allocatable :: d0a(:,:), d0b(:,:), dpck(:,:), fpck(:,:), gfull(:,:), Gd0(:,:)
    real(dp), allocatable :: ba(:,:), bb(:,:), bvec(:,:), uvec(:,:)
    real(dp), allocatable :: nac_rohf_bvec_hf_jk_pulay(:,:)
    real(dp), allocatable :: xa(:,:), xb(:,:), dCa(:,:), dCb(:,:), gp(:,:), gm(:,:)
    real(dp), allocatable :: zneff(:), hess_native(:,:), hresp(:,:)
    real(dp), allocatable :: faop(:), fbop(:)
    integer, allocatable :: iecp_atom(:)
    real(dp) :: hfscale, hstep, gx(3, size(infos%atoms%xyz,2))
    integer :: nbf, nbf2, natom, ncart, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: i, j, a, icart, kc, cc, x, mu, nu, ie, nec
    integer :: nac_dump_env_length, nac_dump_env_status
    logical :: nac_dump_rohf_response
    character(len=32) :: nac_dump_env

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

    nac_dump_env = ''
    call get_environment_variable('NAC_DUMP_ROHF_RESPONSE', nac_dump_env, &
                                  length=nac_dump_env_length, status=nac_dump_env_status)
    nac_dump_rohf_response = nac_dump_env_status == 0 .and. nac_dump_env_length > 0

    if (infos%control%verbose >= 2) then
      write(iw,'(/,A)') 'PyOQP: Native OpenQP open-shell (ROHF) HF Hessian CPHF response prepass'
      write(iw,'(A,I6,A,I6,A,I6,A,I6,A,I6)') '  nbf=', nbf, ' nocca=', nocca, &
        ' noccb=', noccb, ' rhs=', ncart, ' rotdim=', ltot
      write(iw,'(A)') '  Storing native OpenQP open-shell (ROHF) HF analytic Hessian in OQP::hf_hessian.'
    end if

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
      ! Every ECP centre must have been matched to an atom: an unmapped centre
      ! would stay fixed while its atom is displaced, silently corrupting the
      ! semi-numerical response. Abort loudly instead.
      if (count(iecp_atom > 0) /= nec) then
        call show_message('hf_hessian (ROHF): could not map every ECP centre '// &
          'to an atom (coordinate mismatch > 1e-6 bohr); analytic Hessian '// &
          'would be wrong - use [hess] type=numerical for this system.', WITH_ABORT)
      end if
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
    if (os_operator_path(infos)) then
      ! Three blocked derivative-ERI traversals instead of one per pair/spin.
      block
        real(dp), allocatable :: gva(:,:), gvb(:,:)
        integer :: x
        call os_g2e_from_operators(infos, basis, ptot, pa, pb, hfscale, &
          mo, nocca, mo, noccb, gva, gvb)
        do x = 1, ncart
          do a = 1, nvira
            do i = 1, nocca
              ga2e(a,i,x) = gva((a-1)*nocca+i,x)
            end do
          end do
          do a = 1, nvirb
            do i = 1, noccb
              gb2e(a,i,x) = gvb((a-1)*noccb+i,x)
            end do
          end do
        end do
      end block
    else
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
    end if

    ! ===== CPHF right-hand sides (non-canonical Pulay form), packed =====
    allocate(d0a(nbf,nbf), d0b(nbf,nbf), gfull(nbf,nbf), Gd0(nbf,nbf))
    allocate(dpck(nbf2,2), fpck(nbf2,2))
    allocate(ba(nvira,nocca), bb(nvirb,noccb))
    allocate(bvec(ltot,ncart), uvec(ltot,ncart), source=0.0_dp)
    if (nac_dump_rohf_response) &
      allocate(nac_rohf_bvec_hf_jk_pulay(ltot,ncart), source=0.0_dp)
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

        ! Diagnostic snapshot before the ROKS XC skeleton is added.  The native
        ! Fortran layout is (ROHF rotation, Cartesian coordinate) =
        ! (ltot,ncart).  OQPData performs a C-order reshape of this column-major
        ! buffer, so a non-square Python consumer must recover it with
        ! raw.reshape(ncart,ltot).T (a bare transpose is insufficient).
        if (nac_dump_rohf_response) &
          call rohf_pack_trial(nac_rohf_bvec_hf_jk_pulay(:,icart), ba, bb, &
                               nbf, nocca, noccb)

        ! --- XC contribution to the CPKS right-hand side (ROKS only) -----------
        ! Central FD of the spin XC Fock (open-shell dftexcor) along R +/- h AND
        ! occupied MOs reorthonormalized by dmoR^s = -1/2 sum_j C_j S^x_ji; the XC
        ! skeleton dVxc/dR + f_xc[d0], subtracted from B (which carries -F0x).
        if (infos%control%hamilton >= 20) then
          block
            use mod_dft, only: dft_initialize, dftclean, dftexcor
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
            call dftexcor(basis, mgr, int(infos%control%scftype), frap, frbp, mopa, mopb, &
                          nbf, nbf2, exr, telr, tknr, infos)
            call dftclean(infos)
            basis%atoms%xyz(cc,kc) = basis%atoms%xyz(cc,kc) - 2*hxr
            call basis%init_shell_centers()
            call dft_initialize(infos, basis, mgr)
            mopa = mo; mopa(:,1:nocca) = mo(:,1:nocca) - hxr*dmoa
            mopb = mo; mopb(:,1:noccb) = mo(:,1:noccb) - hxr*dmob
            fram = 0.0_dp; frbm = 0.0_dp
            call dftexcor(basis, mgr, int(infos%control%scftype), fram, frbm, mopa, mopb, &
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

    if (nac_dump_rohf_response) then
      block
        real(dp), contiguous, pointer :: dump_bvec_hf_jk_pulay(:,:)
        real(dp), contiguous, pointer :: dump_bvec_full(:,:), dump_uvec(:,:)

        call infos%dat%erase((/ character(len=80) :: &
          'OQP::nac_rohf_bvec_hf_jk_pulay', &
          'OQP::nac_rohf_bvec_full', &
          'OQP::nac_rohf_uvec' /))
        ! All three records are stored in native Fortran order (ltot,ncart):
        ! first index = packed ds/dv/sv ROHF rotation, second = Cartesian
        ! coordinate.  Python must use raw.reshape(ncart,ltot).T to undo the
        ! C-order view of the column-major storage.
        call tagarray_reserve_data(infos%dat, 'OQP::nac_rohf_bvec_hf_jk_pulay', &
          TA_TYPE_REAL64, ltot*ncart, (/ ltot, ncart /), &
          comment='ROHF HF+JK/Pulay CPHF RHS; Fortran (rotation,Cartesian)')
        call tagarray_reserve_data(infos%dat, 'OQP::nac_rohf_bvec_full', &
          TA_TYPE_REAL64, ltot*ncart, (/ ltot, ncart /), &
          comment='ROHF full HF+JK+XC CPHF RHS; Fortran (rotation,Cartesian)')
        call tagarray_reserve_data(infos%dat, 'OQP::nac_rohf_uvec', &
          TA_TYPE_REAL64, ltot*ncart, (/ ltot, ncart /), &
          comment='ROHF CPHF response U; Fortran (rotation,Cartesian)')
        call tagarray_get_data(infos%dat, 'OQP::nac_rohf_bvec_hf_jk_pulay', &
                               dump_bvec_hf_jk_pulay)
        call tagarray_get_data(infos%dat, 'OQP::nac_rohf_bvec_full', dump_bvec_full)
        call tagarray_get_data(infos%dat, 'OQP::nac_rohf_uvec', dump_uvec)
        dump_bvec_hf_jk_pulay = nac_rohf_bvec_hf_jk_pulay
        dump_bvec_full = bvec
        dump_uvec = uvec
      end block
      deallocate(nac_rohf_bvec_hf_jk_pulay)
    end if

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
    want_props = hf_hess_properties_wanted(infos)
    if (want_props) then
      allocate(dptx(nbf,nbf))
      call hf_dipder_init(infos, ptot, dipf, dmu)
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

      ! relaxed total density derivative -> dipole derivative (IR intensities)
      if (want_props) then
        call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, dCa, nbf, mo, nbf, 0.0_dp, dptx, nbf)
        call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, dCb, nbf, mo, nbf, 1.0_dp, dptx, nbf)
        dptx = dptx + transpose(dptx)
        call hf_dipder_add_response(dipf, dptx, dmu(:,x))
      end if
    end do

    ! analytic polarizability derivatives (Raman activities); the ECP enters
    ! only through h^x (ecp_deriv_ints is folded into dVa above)
    if (want_props) then
    call hf_dipder_store(infos, dmu)
    deallocate(dipf, dmu, dptx)
    block
      use oqp_tagarray_driver, only: OQP_hf_polarizability_derivatives
      real(dp), allocatable :: dpol(:,:,:), fa_ao(:,:), fb_ao(:,:)
      real(dp), contiguous, pointer :: pstore(:,:,:)
      allocate(dpol(3,3,ncart), fa_ao(nbf,nbf), fb_ao(nbf,nbf))
      call unpack_matrix(focka, fa_ao)
      call unpack_matrix(fockb, fb_ao)
      call hf_polder_rohf(infos, mo, fa_ao, fb_ao, pa, pb, nocca, noccb, &
                          dSa, dTa + dVa, uvec, hfscale, dpol)
      call infos%dat%alloc_or_die(OQP_hf_polarizability_derivatives, (/ 3, 3, ncart /), pstore, &
        description='Analytic nuclear derivatives of the static polarizability (a.u.), (3,3,3N)')
      pstore = dpol
      deallocate(dpol, fa_ao, fb_ao)
    end block
    end if

    ! The central difference of the ELECTRONIC gradient over geometry AND the
    ! relaxed orbital path already contains the full electronic Hessian (skeleton
    ! + orbital-relaxation response); only the (orbital-independent) nuclear
    ! repulsion second derivative is added analytically.
    allocate(hess_native(ncart,ncart))
    hess_native = 0.5_dp*(hresp + transpose(hresp))
    call hess_nn(basis%atoms, basis%ecp_zn_num, hess_native)

    call infos%dat%alloc_or_die(OQP_hf_hessian, (/ ncart, ncart /), hess_store, &
      description='Native OpenQP open-shell (ROHF) HF analytic Hessian matrix')
    hess_store = hess_native
    if (infos%control%verbose >= 2) then
      write(iw,'(A)') 'PyOQP: Native OpenQP open-shell (ROHF) HF Hessian matrix stored'
    end if

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
      type(grd2_uhf_compute_data_t) :: gc

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
          use mod_dft, only: dft_initialize, dftclean, dftexcor
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
          call dftexcor(basis, mg, int(infos%control%scftype), fra, frb, mopa, mopb, &
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
      call gc%build_cart(basis)
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

!> @brief Analytic nuclear derivatives of the static dipole polarizability,
!>        closed-shell HF (Raman activities).
!> @details  With the field responses U^a (A U^a = -mu^a) the polarizability is
!>   the stationary value of the Hylleraas functional
!>     alpha_ab = -4 [ mu^a.U^b + mu^b.U^a + U^a.A.U^b ],
!>   so (2n+1 rule) its nuclear derivative needs no second-order response: keep
!>   U^a, U^b fixed and differentiate the MO-basis quantities along the relaxed
!>   orbital path dC = C T^x with T_vo = U^x, T_oo = -S^x_oo/2,
!>   T_vv = -S^x_vv/2, T_ov = -U^x' - S^x_ov:
!>     d alpha_ab/dx = -4 [ mu^a(x).U^b + mu^b(x).U^a + (U^a.A.U^b)(x) ].
!>   With Pt^a = Co U^a Cv' + Cv U^a' Co', Y = U^a'U^b, Z = U^a U^b',
!>   W = Cv Y Cv' - Co Z Co' and G the closed-shell response Fock (fock_jk),
!>     U^a.A.U^b = Tr[F_vv Y] - Tr[F_oo Z] + Tr[Pt^a G[Pt^b]],
!>   whose derivative is
!>     Tr[(h^x + G^x[P] + G[dP^x]) W] - Tr[S^x_vv/2 (eps.Y + Y.eps)]
!>       + Tr[S^x_oo/2 (eps.Z + Z.eps)] + Tr[Pt^a G^x[Pt^b]]
!>       + 2 Tr[T W^a G^b] + 2 Tr[T W^b G^a]   (W^a, G^b in the MO basis).
!>   The derivative-integral traces come from fock_deriv_contract (one
!>   traversal per field pair for G^x[P].W and one for Pt^a.G^x[Pt^b]), so the
!>   whole tensor costs 3 CPHF right-hand sides, 9 response Fock builds and 12
!>   derivative-integral traversals, independent of the number of atoms.
!>   For Kohn-Sham the exchange-correlation pieces -- Tr[dVxc/dx W] and the
!>   derivative of 1/2 Tr[Pt^a f_xc[Pt^b]] -- are taken by central differences
!>   along the relaxed path (R +/- h, C +/- h C T^x) with the moving grid, the
!>   same treatment the Hessian uses for its XC terms: no SCF or CPKS re-solve,
!>   two Vxc and six f_xc builds per coordinate.
  subroutine hf_polder_rhf(infos, mo, eps, pfull, sflat, hflat, uvec, dPx, &
                           nocc, nvir, hfscale, dpol, xc_dhse, xc_dfxc)
    use oqp_linalg
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int1, only: multipole_integrals
    use grd1, only: der_dipole_matrix
    use mathlib, only: unpack_matrix, pack_matrix
    use fock_deriv_mod, only: fock_deriv_contract
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve
    use io_constants, only: iw
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: mo(:,:), eps(:), pfull(:,:), sflat(:,:,:), hflat(:,:,:)
    real(dp), intent(in) :: uvec(:,:), dPx(:,:,:), hfscale
    integer, intent(in) :: nocc, nvir
    real(dp), intent(out) :: dpol(:,:,:)            ! (3,3,ncart)
    !> Optional (DFT): while the XC terms below displace every coordinate on
    !> the moving grid, also return the two XC pieces of the Hessian response
    !> assembly that use the same displaced grids and the same relaxed
    !> occupied orbitals, so the caller does not rebuild them:
    !>   xc_dhse(:,y)   = d/dR_y of the XC gradient along R+l, P+l dP^y
    !>   xc_dfxc(:,:,y) = d/dR_y of the XC Fock matrix along the relaxed path
    real(dp), intent(out), optional :: xc_dhse(:,:), xc_dfxc(:,:,:)

    type(basis_set), pointer :: basis
    real(dp), allocatable :: mints(:,:), dfull(:,:,:), mmo(:,:,:), dD(:,:,:,:,:)
    real(dp), allocatable :: bF(:,:), uF(:,:), ua(:,:,:), xa(:,:,:), pta(:,:,:)
    real(dp), allocatable :: ga(:,:,:), gamo(:,:,:), wsym(:,:,:), gw(:,:,:)
    real(dp), allocatable :: gwx(:,:,:), gppx(:,:,:), yab(:,:,:), zab(:,:,:)
    real(dp), allocatable :: tmat(:,:), smo(:,:), scr(:,:), scr2(:,:), wmo(:,:)
    real(dp), allocatable :: dpk(:,:), fpk(:,:), gx(:,:), ux(:,:)
    real(dp) :: origin(3), t1(3,3), t2, t3, alpha(3,3), hyl, a1, a2
    integer :: nbf, nbf2, natom, ncart, lexc, a, b, q, i, j, k, x, kc, cc, ia, ip

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2)
    ncart = 3*natom
    lexc = nocc*nvir
    origin = 0.0_dp

    ! ---- dipole integrals (normalized), MO dipole, field CPHF -------------
    allocate(mints(nbf2,19), source=0.0_dp)
    call multipole_integrals(basis, mints, origin, 3)
    allocate(dfull(nbf,nbf,3), mmo(nbf,nbf,3), scr(nbf,nbf), scr2(nbf,nbf))
    do a = 1, 3
      call unpack_matrix(mints(:,a), dfull(:,:,a))
      call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,dfull(:,:,a),nbf,0.0_dp,scr,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,mmo(:,:,a),nbf)
    end do
    deallocate(mints)
    allocate(bF(lexc,3), uF(lexc,3), source=0.0_dp)
    do a = 1, 3
      do q = 1, nvir
        do i = 1, nocc
          bF((q-1)*nocc+i,a) = -mmo(i,nocc+q,a)
        end do
      end do
    end do
    call cphf_solve(infos, 3, bF, uF)
    do b = 1, 3
      do a = 1, 3
        alpha(a,b) = 4.0_dp*sum(bF(:,a)*uF(:,b))
      end do
    end do

    ! ---- field-response matrices -----------------------------------------
    allocate(ua(nocc,nvir,3), xa(nbf,nbf,3), pta(nbf,nbf,3), ga(nbf,nbf,3), gamo(nbf,nbf,3))
    allocate(dpk(nbf2,1), fpk(nbf2,1))
    do a = 1, 3
      ua(:,:,a) = reshape(uF(:,a), [nocc,nvir])
      call dgemm('n','n',nbf,nvir,nocc,1.0_dp,mo,nbf,ua(:,:,a),nocc,0.0_dp,scr,nbf)
      call dgemm('n','t',nbf,nbf,nvir,1.0_dp,scr,nbf,mo(:,nocc+1:),nbf,0.0_dp,xa(:,:,a),nbf)
      pta(:,:,a) = xa(:,:,a) + transpose(xa(:,:,a))
      call pack_matrix(pta(:,:,a), dpk(:,1))
      fpk = 0.0_dp
      call fock_jk(basis, d=dpk, f=fpk, scale_exch=hfscale, infos=infos)
      call unpack_from_packed(fpk(:,1), ga(:,:,a), nbf)
      call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,ga(:,:,a),nbf,0.0_dp,scr,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,gamo(:,:,a),nbf)
    end do

    ! Self-check of the conventions: the Hylleraas value equals alpha.
    if (infos%control%verbose >= 2) then
      do b = 1, 3
        do a = 1, 3
          hyl = sum(bF(:,a)*uF(:,b)) + sum(bF(:,b)*uF(:,a))        ! = -(mu^a.U^b + mu^b.U^a)
          a1 = 0.0_dp
          do q = 1, nvir
            do i = 1, nocc
              a1 = a1 + ua(i,q,a)*ua(i,q,b)*(eps(nocc+q) - eps(i))
            end do
          end do
          a2 = sum(pta(:,:,a)*ga(:,:,b))
          write(iw,'(A,2I2,3ES16.8)') '  polder Hylleraas check a b, alpha, -4[..], diff:', a, b, &
            alpha(a,b), 4.0_dp*hyl - 4.0_dp*(a1 + a2), alpha(a,b) - (4.0_dp*hyl - 4.0_dp*(a1 + a2))
        end do
      end do
    end if

    ! ---- pair intermediates: W_ab, G[W_ab], derivative-integral traces ---
    allocate(yab(nvir,nvir,6), zab(nocc,nocc,6), wsym(nbf,nbf,6), gw(nbf,nbf,6))
    allocate(gwx(3,natom,6), gppx(3,natom,6), gx(3,natom))
    ip = 0
    do b = 1, 3
      do a = 1, b
        ip = ip + 1
        call dgemm('t','n',nvir,nvir,nocc,1.0_dp,ua(:,:,a),nocc,ua(:,:,b),nocc,0.0_dp,yab(:,:,ip),nvir)
        call dgemm('n','t',nocc,nocc,nvir,1.0_dp,ua(:,:,a),nocc,ua(:,:,b),nocc,0.0_dp,zab(:,:,ip),nocc)
        call dgemm('n','n',nbf,nvir,nvir,1.0_dp,mo(:,nocc+1:),nbf,yab(:,:,ip),nvir,0.0_dp,scr,nbf)
        call dgemm('n','t',nbf,nbf,nvir,1.0_dp,scr,nbf,mo(:,nocc+1:),nbf,0.0_dp,scr2,nbf)
        call dgemm('n','n',nbf,nocc,nocc,1.0_dp,mo,nbf,zab(:,:,ip),nocc,0.0_dp,scr,nbf)
        call dgemm('n','t',nbf,nbf,nocc,-1.0_dp,scr,nbf,mo,nbf,1.0_dp,scr2,nbf)
        wsym(:,:,ip) = 0.5_dp*(scr2 + transpose(scr2))
        call pack_matrix(wsym(:,:,ip), dpk(:,1))
        fpk = 0.0_dp
        call fock_jk(basis, d=dpk, f=fpk, scale_exch=hfscale, infos=infos)
        call unpack_from_packed(fpk(:,1), gw(:,:,ip), nbf)
        call fock_deriv_contract(infos, basis, pfull, wsym(:,:,ip), hfscale, gx)
        gwx(:,:,ip) = 2.0_dp*gx
        call fock_deriv_contract(infos, basis, pta(:,:,b), pta(:,:,a), hfscale, gx)
        gppx(:,:,ip) = 2.0_dp*gx
      end do
    end do

    ! ---- dipole-integral derivatives -------------------------------------
    allocate(dD(nbf,nbf,3,natom,3))
    call der_dipole_matrix(basis, origin, dD)
    do a = 1, 3
      do k = 1, natom
        do cc = 1, 3
          do j = 1, nbf
            dD(:,j,cc,k,a) = dD(:,j,cc,k,a)*basis%bfnrm(:)*basis%bfnrm(j)
          end do
        end do
      end do
    end do

    ! ---- per nuclear coordinate ------------------------------------------
    allocate(tmat(nbf,nbf), smo(nbf,nbf), wmo(nbf,nbf), ux(nocc,nvir))
    dpol = 0.0_dp
    do x = 1, ncart
      cc = mod(x-1,3) + 1; kc = (x-1)/3 + 1
      call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,sflat(:,:,x),nbf,0.0_dp,scr,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,smo,nbf)
      ux = reshape(uvec(:,x), [nocc,nvir])
      tmat = -0.5_dp*smo
      tmat(nocc+1:,1:nocc) = transpose(ux)
      tmat(1:nocc,nocc+1:) = -ux - smo(1:nocc,nocc+1:)

      ! mu^a(x).U^b = Tr[D_a^x X^b] + sum_ia U^b_ia ([T'M_a]_ia + [M_a T]_ia)
      do a = 1, 3
        call dgemm('t','n',nbf,nbf,nbf,1.0_dp,tmat,nbf,mmo(:,:,a),nbf,0.0_dp,scr,nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mmo(:,:,a),nbf,tmat,nbf,1.0_dp,scr,nbf)
        do b = 1, 3
          t1(a,b) = sum(dD(:,:,cc,kc,a)*xa(:,:,b)) + sum(ua(:,:,b)*scr(1:nocc,nocc+1:))
        end do
      end do

      ip = 0
      do b = 1, 3
        do a = 1, b
          ip = ip + 1
          ! F part: Tr[(h^x + G^x[P] + G[dP^x]) W] and the S^x energy weights
          t2 = sum(hflat(:,:,x)*wsym(:,:,ip)) + gwx(cc,kc,ip) + sum(dPx(:,:,x)*gw(:,:,ip))
          do j = 1, nvir
            do i = 1, nvir
              t2 = t2 - 0.5_dp*smo(nocc+i,nocc+j)*(eps(nocc+i) + eps(nocc+j))*yab(i,j,ip)
            end do
          end do
          do j = 1, nocc
            do i = 1, nocc
              t2 = t2 + 0.5_dp*smo(i,j)*(eps(i) + eps(j))*zab(i,j,ip)
            end do
          end do
          ! two-electron part: Tr[Pt^a G^x[Pt^b]] + 2Tr[T W^a G^b] + 2Tr[T W^b G^a]
          t3 = gppx(cc,kc,ip) + 2.0_dp*tw_g(a, b) + 2.0_dp*tw_g(b, a)
          dpol(a,b,x) = -4.0_dp*(t1(a,b) + t1(b,a) + t2 + t3)
          dpol(b,a,x) = dpol(a,b,x)
        end do
      end do
    end do

    if (infos%control%hamilton == 20) call add_xc_terms()

    deallocate(dfull, mmo, dD, bF, uF, ua, xa, pta, ga, gamo, wsym, gw, gwx, gppx, &
               yab, zab, tmat, smo, scr, scr2, wmo, dpk, fpk, gx, ux)

  contains

    !> Kohn-Sham exchange-correlation contributions by central differences
    !> along the relaxed orbital path with the moving grid.
    subroutine add_xc_terms()
      use mod_dft, only: dft_initialize, dftclean, dftexcor
      use mod_dft_molgrid, only: dft_grid_t
      use mod_dft_gridint_fxc, only: tddft_fxc
      use mod_dft_gridint_grad, only: derexc_blk
      type(dft_grid_t) :: mg
      real(dp), parameter :: hxc = 1.0e-3_dp
      real(dp), allocatable :: mos(:,:), dmo(:,:), fr(:), vxc(:,:,:), fx(:,:,:), dx(:,:,:), pts(:,:,:)
      real(dp), allocatable :: dap(:,:), ded(:,:,:)
      real(dp) :: sab(3,3,2), exr, telr, tknr, sgn, tele, tkin
      integer :: is, xx, kx, cx, ipx, aa, bb, nang
      logical :: share

      share = present(xc_dhse) .and. present(xc_dfxc)
      nang = maxval(basis%am) + 2
      allocate(mos(nbf,nbf), dmo(nbf,nbf), fr(nbf2), vxc(nbf,nbf,2), &
               fx(nbf,nbf,3), dx(nbf,nbf,3), pts(nbf,nbf,3))
      if (share) allocate(dap(nbf,nbf), ded(3,natom,2))
      ! same warm-up as the Hessian response-assembly loop: flush grid state
      ! left by the CPHF solves
      if (share) then
        call dft_initialize(infos, basis, mg); call dftclean(infos)
      end if
      do xx = 1, ncart
        cx = mod(xx-1,3) + 1; kx = (xx-1)/3 + 1
        call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,sflat(:,:,xx),nbf,0.0_dp,scr,nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,smo,nbf)
        ux = reshape(uvec(:,xx), [nocc,nvir])
        tmat = -0.5_dp*smo
        tmat(nocc+1:,1:nocc) = transpose(ux)
        tmat(1:nocc,nocc+1:) = -ux - smo(1:nocc,nocc+1:)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mo,nbf,tmat,nbf,0.0_dp,dmo,nbf)
        do is = 1, 2
          sgn = merge(1.0_dp, -1.0_dp, is == 1)
          basis%atoms%xyz(cx,kx) = basis%atoms%xyz(cx,kx) + sgn*hxc
          call basis%init_shell_centers()
          call dft_initialize(infos, basis, mg)
          if (share) then
            ! skeleton + density response of the XC gradient (Hessian dHse)
            dap = pfull + sgn*hxc*dPx(:,:,xx); ded(:,:,is) = 0.0_dp
            call derexc_blk(basis, mg, dap, dap, ded(:,:,is), tele, tkin, nang, nbf, &
                            infos%dft%grid_density_cutoff, .false., infos)
          end if
          mos = mo + sgn*hxc*dmo
          fr = 0.0_dp
          call dftexcor(basis, mg, 1, fr, fr, mos, mos, nbf, nbf2, exr, telr, tknr, infos)
          call unpack_from_packed(fr, vxc(:,:,is), nbf)
          do bb = 1, 3
            wmo = 0.0_dp
            wmo(1:nocc,nocc+1:) = ua(:,:,bb)
            wmo(nocc+1:,1:nocc) = transpose(ua(:,:,bb))
            call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mos,nbf,wmo,nbf,0.0_dp,scr2,nbf)
            call dgemm('n','t',nbf,nbf,nbf,1.0_dp,scr2,nbf,mos,nbf,0.0_dp,pts(:,:,bb),nbf)
          end do
          dx = pts
          fx = 0.0_dp
          call tddft_fxc(basis=basis, molGrid=mg, isVecs=.true., wf=mos, fx=fx, dx=dx, &
                         nmtx=3, threshold=0.0_dp, infos=infos)
          do bb = 1, 3
            do aa = 1, 3
              sab(aa,bb,is) = sum(pts(:,:,aa)*fx(:,:,bb))
            end do
          end do
          call dftclean(infos)
          basis%atoms%xyz(cx,kx) = basis%atoms%xyz(cx,kx) - sgn*hxc
          call basis%init_shell_centers()
        end do
        vxc(:,:,1) = (vxc(:,:,1) - vxc(:,:,2))/(2.0_dp*hxc)
        if (share) then
          xc_dhse(:,xx) = reshape((ded(:,:,1) - ded(:,:,2))/(2.0_dp*hxc), [ncart])
          xc_dfxc(:,:,xx) = vxc(:,:,1)
        end if
        ipx = 0
        do bb = 1, 3
          do aa = 1, bb
            ipx = ipx + 1
            dpol(aa,bb,xx) = dpol(aa,bb,xx) - 4.0_dp*( sum(vxc(:,:,1)*wsym(:,:,ipx)) &
                + 0.25_dp*(sab(aa,bb,1) + sab(bb,aa,1) - sab(aa,bb,2) - sab(bb,aa,2))/(2.0_dp*hxc) )
            dpol(bb,aa,xx) = dpol(aa,bb,xx)
          end do
        end do
      end do
      deallocate(mos, dmo, fr, vxc, fx, dx, pts)
      if (share) deallocate(dap, ded)
    end subroutine add_xc_terms

    !> Tr[T W^a G^b] with W^a the symmetric MO matrix carrying U^a in its
    !> occ-vir and vir-occ blocks and G^b = C' G[Pt^b] C.
    real(dp) function tw_g(ia_, ib_)
      integer, intent(in) :: ia_, ib_
      wmo = 0.0_dp
      wmo(1:nocc,nocc+1:) = ua(:,:,ia_)
      wmo(nocc+1:,1:nocc) = transpose(ua(:,:,ia_))
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,tmat,nbf,wmo,nbf,0.0_dp,scr2,nbf)
      tw_g = sum(scr2*transpose(gamo(:,:,ib_)))
    end function tw_g
  end subroutine hf_polder_rhf

!###############################################################################

!> @brief Analytic nuclear derivatives of the static polarizability, UHF/UKS.
!> @details Spin-unrestricted counterpart of hf_polder_rhf. With
!>   alpha_ab = -2 sum_s mu^{a,s}.U^{b,s}, the Hylleraas functional is
!>   -2 [ mu^a.U^b + mu^b.U^a + U^a.A.U^b ] with
!>     U^a.A.U^b = sum_s ( Tr[F^s_vv Y^s] - Tr[F^s_oo Z^s] )
!>               + 1/2 sum_s Tr[Pt^{a,s} G^s[Pt^b]],
!>   G^s[D] = J[D^a+D^b] - c K[D^s] (fock_jk, open shell). Each spin follows
!>   its own relaxed orbital path T^{x,s}; the derivative-integral traces use
!>   fock_deriv_contract_os, which returns Tr[M G^{s,x}[P]] directly. The
!>   Kohn-Sham XC pieces are central differences along the relaxed path with
!>   the moving grid, as in the closed-shell routine.
  subroutine hf_polder_uhf(infos, moa, mob, epsa, epsb, pa, pb, nocca, noccb, &
                           sflat, hflat, uvec, dpxa, dpxb, hfscale, dpol)
    use oqp_linalg
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int1, only: multipole_integrals
    use grd1, only: der_dipole_matrix
    use mathlib, only: unpack_matrix, pack_matrix
    use fock_deriv_mod, only: fock_deriv_contract_os
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve_uhf
    use io_constants, only: iw
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: moa(:,:), mob(:,:), epsa(:), epsb(:), pa(:,:), pb(:,:)
    integer, intent(in) :: nocca, noccb
    real(dp), intent(in) :: sflat(:,:,:), hflat(:,:,:), uvec(:,:), dpxa(:,:,:), dpxb(:,:,:)
    real(dp), intent(in) :: hfscale
    real(dp), intent(out) :: dpol(:,:,:)

    type(basis_set), pointer :: basis
    real(dp), allocatable :: mo(:,:,:), eps(:,:), mints(:,:), dfull(:,:,:), mmo(:,:,:,:), dD(:,:,:,:,:)
    real(dp), allocatable :: bF(:,:), uF(:,:), xa(:,:,:,:), pta(:,:,:,:), ga(:,:,:,:), gamo(:,:,:,:)
    real(dp), allocatable :: wsym(:,:,:,:), gw(:,:,:,:), gwx(:,:,:), gppx(:,:,:), ptot(:,:), psp(:,:,:)
    real(dp), allocatable :: tmat(:,:,:), smo(:,:), scr(:,:), scr2(:,:), wmo(:,:), dpk(:,:), fpk(:,:), gx(:,:)
    real(dp) :: origin(3), t1(3,3), t2, t3, alpha(3,3)
    integer :: nbf, nbf2, natom, ncart, s, no(2), nv(2), loff(2), ltot
    !> field responses U^{a,s} as (nocc,nvir,3) per spin
    type :: spin_resp_t
      real(dp), allocatable :: u(:,:,:)
    end type
    type(spin_resp_t) :: us(2)
    integer :: a, b, q, i, j, k, x, kc, cc, ip

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf; nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2); ncart = 3*natom
    no = [nocca, noccb]; nv = nbf - no
    loff = [0, nocca*(nbf-nocca)]; ltot = loff(2) + noccb*(nbf-noccb)
    origin = 0.0_dp
    allocate(mo(nbf,nbf,2), eps(nbf,2), psp(nbf,nbf,2), ptot(nbf,nbf))
    mo(:,:,1) = moa; mo(:,:,2) = mob; eps(:,1) = epsa; eps(:,2) = epsb
    psp(:,:,1) = pa; psp(:,:,2) = pb; ptot = pa + pb

    allocate(mints(nbf2,19), source=0.0_dp)
    call multipole_integrals(basis, mints, origin, 3)
    allocate(dfull(nbf,nbf,3), mmo(nbf,nbf,3,2), scr(nbf,nbf), scr2(nbf,nbf))
    allocate(bF(ltot,3), uF(ltot,3), source=0.0_dp)
    do a = 1, 3
      call unpack_matrix(mints(:,a), dfull(:,:,a))
      do s = 1, 2
        call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo(:,:,s),nbf,dfull(:,:,a),nbf,0.0_dp,scr,nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo(:,:,s),nbf,0.0_dp,mmo(:,:,a,s),nbf)
        do q = 1, nv(s)
          do i = 1, no(s)
            bF(loff(s)+(q-1)*no(s)+i,a) = -mmo(i,no(s)+q,a,s)
          end do
        end do
      end do
    end do
    deallocate(mints)
    call cphf_solve_uhf(infos, 3, bF, uF)
    do s = 1, 2
      allocate(us(s)%u(no(s), nv(s), 3))
      do a = 1, 3
        us(s)%u(:,:,a) = reshape(uF(loff(s)+1:loff(s)+no(s)*nv(s), a), [no(s), nv(s)])
      end do
    end do
    do b = 1, 3
      do a = 1, 3
        alpha(a,b) = 2.0_dp*sum(bF(:,a)*uF(:,b))
      end do
    end do

    allocate(xa(nbf,nbf,3,2), pta(nbf,nbf,3,2), ga(nbf,nbf,3,2), gamo(nbf,nbf,3,2))
    allocate(dpk(nbf2,2), fpk(nbf2,2))
    do a = 1, 3
      do s = 1, 2
        call dgemm('n','n',nbf,nv(s),no(s),1.0_dp,mo(:,:,s),nbf,us(s)%u(:,:,a),no(s),0.0_dp,scr,nbf)
        call dgemm('n','t',nbf,nbf,nv(s),1.0_dp,scr,nbf,mo(:,no(s)+1:,s),nbf,0.0_dp,xa(:,:,a,s),nbf)
        pta(:,:,a,s) = xa(:,:,a,s) + transpose(xa(:,:,a,s))
        call pack_matrix(pta(:,:,a,s), dpk(:,s))
      end do
      fpk = 0.0_dp
      call fock_jk(basis, d=dpk, f=fpk, scale_exch=hfscale, infos=infos)
      do s = 1, 2
        call unpack_from_packed(fpk(:,s), ga(:,:,a,s), nbf)
        call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo(:,:,s),nbf,ga(:,:,a,s),nbf,0.0_dp,scr,nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo(:,:,s),nbf,0.0_dp,gamo(:,:,a,s),nbf)
      end do
    end do

    if (infos%control%verbose >= 2) then
      do b = 1, 3
        do a = 1, 3
          t2 = 0.0_dp
          do s = 1, 2
            do q = 1, nv(s)
              do i = 1, no(s)
                t2 = t2 + uF(loff(s)+(q-1)*no(s)+i,a)*uF(loff(s)+(q-1)*no(s)+i,b) &
                          *(eps(no(s)+q,s) - eps(i,s))
              end do
            end do
            t2 = t2 + 0.5_dp*sum(pta(:,:,a,s)*ga(:,:,b,s))
          end do
          write(iw,'(A,2I2,2ES16.8)') '  polder(UHF) Hylleraas check a b, alpha, diff:', a, b, alpha(a,b), &
            alpha(a,b) - (2.0_dp*(sum(bF(:,a)*uF(:,b)) + sum(bF(:,b)*uF(:,a))) - 2.0_dp*t2)
        end do
      end do
    end if

    ! pair intermediates
    allocate(wsym(nbf,nbf,6,2), gw(nbf,nbf,6,2), gwx(3,natom,6), gppx(3,natom,6), gx(3,natom))
    gwx = 0.0_dp; gppx = 0.0_dp
    ip = 0
    do b = 1, 3
      do a = 1, b
        ip = ip + 1
        do s = 1, 2
          ! W^s = Cv (U^a' U^b) Cv' - Co (U^a U^b') Co'
          scr2 = 0.0_dp
          call vv_oo_w(a, b, s, scr2)
          wsym(:,:,ip,s) = 0.5_dp*(scr2 + transpose(scr2))
          call pack_matrix(wsym(:,:,ip,s), dpk(:,s))
          call fock_deriv_contract_os(infos, basis, ptot, psp(:,:,s), wsym(:,:,ip,s), hfscale, gx)
          gwx(:,:,ip) = gwx(:,:,ip) + gx
          call fock_deriv_contract_os(infos, basis, pta(:,:,b,1) + pta(:,:,b,2), pta(:,:,b,s), &
                                      pta(:,:,a,s), hfscale, gx)
          gppx(:,:,ip) = gppx(:,:,ip) + 0.5_dp*gx
        end do
        fpk = 0.0_dp
        call fock_jk(basis, d=dpk, f=fpk, scale_exch=hfscale, infos=infos)
        do s = 1, 2
          call unpack_from_packed(fpk(:,s), gw(:,:,ip,s), nbf)
        end do
      end do
    end do

    allocate(dD(nbf,nbf,3,natom,3))
    call der_dipole_matrix(basis, origin, dD)
    do a = 1, 3
      do k = 1, natom
        do cc = 1, 3
          do j = 1, nbf
            dD(:,j,cc,k,a) = dD(:,j,cc,k,a)*basis%bfnrm(:)*basis%bfnrm(j)
          end do
        end do
      end do
    end do

    allocate(tmat(nbf,nbf,2), smo(nbf,nbf), wmo(nbf,nbf))
    dpol = 0.0_dp
    do x = 1, ncart
      cc = mod(x-1,3) + 1; kc = (x-1)/3 + 1
      t1 = 0.0_dp
      do s = 1, 2
        call build_t(x, s)
        do a = 1, 3
          call dgemm('t','n',nbf,nbf,nbf,1.0_dp,tmat(:,:,s),nbf,mmo(:,:,a,s),nbf,0.0_dp,scr,nbf)
          call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mmo(:,:,a,s),nbf,tmat(:,:,s),nbf,1.0_dp,scr,nbf)
          do b = 1, 3
            t1(a,b) = t1(a,b) + sum(dD(:,:,cc,kc,a)*xa(:,:,b,s)) &
                      + sum(us(s)%u(:,:,b)*scr(1:no(s),no(s)+1:))
          end do
        end do
      end do
      ip = 0
      do b = 1, 3
        do a = 1, b
          ip = ip + 1
          t2 = gwx(cc,kc,ip) + sum(dpxa(:,:,x)*gw(:,:,ip,1)) + sum(dpxb(:,:,x)*gw(:,:,ip,2))
          t3 = gppx(cc,kc,ip)
          do s = 1, 2
            t2 = t2 + sum(hflat(:,:,x)*wsym(:,:,ip,s)) + sw_terms(a, b, s, x)
            t3 = t3 + tw_g(a, b, s) + tw_g(b, a, s)
          end do
          dpol(a,b,x) = -2.0_dp*(t1(a,b) + t1(b,a) + t2 + t3)
          dpol(b,a,x) = dpol(a,b,x)
        end do
      end do
    end do

    if (infos%control%hamilton == 20) call add_xc_terms_uhf()

    deallocate(mo, eps, psp, ptot, dfull, mmo, dD, bF, uF, xa, pta, ga, gamo, wsym, gw, gwx, gppx, &
               tmat, smo, scr, scr2, wmo, dpk, fpk, gx)

  contains


    !> w = Cv (U^a' U^b) Cv' - Co (U^a U^b') Co' for spin s
    subroutine vv_oo_w(a_, b_, s_, w)
      integer, intent(in) :: a_, b_, s_
      real(dp), intent(out) :: w(:,:)
      real(dp), allocatable :: y(:,:), z(:,:), t(:,:)
      allocate(y(nv(s_),nv(s_)), z(no(s_),no(s_)), t(nbf,max(nv(s_),no(s_))))
      y = matmul(transpose(us(s_)%u(:,:,a_)), us(s_)%u(:,:,b_))
      z = matmul(us(s_)%u(:,:,a_), transpose(us(s_)%u(:,:,b_)))
      call dgemm('n','n',nbf,nv(s_),nv(s_),1.0_dp,mo(:,no(s_)+1:,s_),nbf,y,nv(s_),0.0_dp,t,nbf)
      call dgemm('n','t',nbf,nbf,nv(s_),1.0_dp,t,nbf,mo(:,no(s_)+1:,s_),nbf,0.0_dp,w,nbf)
      call dgemm('n','n',nbf,no(s_),no(s_),1.0_dp,mo(:,:,s_),nbf,z,no(s_),0.0_dp,t,nbf)
      call dgemm('n','t',nbf,nbf,no(s_),-1.0_dp,t,nbf,mo(:,:,s_),nbf,1.0_dp,w,nbf)
    end subroutine vv_oo_w

    !> S^x energy-weight terms of the Fock-block derivative for spin s
    real(dp) function sw_terms(a_, b_, s_, x_)
      integer, intent(in) :: a_, b_, s_, x_
      real(dp), allocatable :: y(:,:), z(:,:), sm(:,:)
      integer :: ii, jj
      allocate(sm(nbf,nbf))
      call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo(:,:,s_),nbf,sflat(:,:,x_),nbf,0.0_dp,scr2,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr2,nbf,mo(:,:,s_),nbf,0.0_dp,sm,nbf)
      y = matmul(transpose(us(s_)%u(:,:,a_)), us(s_)%u(:,:,b_))
      z = matmul(us(s_)%u(:,:,a_), transpose(us(s_)%u(:,:,b_)))
      sw_terms = 0.0_dp
      do jj = 1, nv(s_)
        do ii = 1, nv(s_)
          sw_terms = sw_terms - 0.5_dp*sm(no(s_)+ii,no(s_)+jj)*(eps(no(s_)+ii,s_) + eps(no(s_)+jj,s_))*y(ii,jj)
        end do
      end do
      do jj = 1, no(s_)
        do ii = 1, no(s_)
          sw_terms = sw_terms + 0.5_dp*sm(ii,jj)*(eps(ii,s_) + eps(jj,s_))*z(ii,jj)
        end do
      end do
    end function sw_terms

    !> relaxed orbital-path generator T^{x,s}
    subroutine build_t(x_, s_)
      integer, intent(in) :: x_, s_
      real(dp), allocatable :: ux(:,:)
      call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo(:,:,s_),nbf,sflat(:,:,x_),nbf,0.0_dp,scr2,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr2,nbf,mo(:,:,s_),nbf,0.0_dp,smo,nbf)
      ux = reshape(uvec(loff(s_)+1:loff(s_)+no(s_)*nv(s_), x_), [no(s_), nv(s_)])
      tmat(:,:,s_) = -0.5_dp*smo
      tmat(no(s_)+1:,1:no(s_),s_) = transpose(ux)
      tmat(1:no(s_),no(s_)+1:,s_) = -ux - smo(1:no(s_),no(s_)+1:)
    end subroutine build_t

    !> 1/2 * 2 Tr[T^s W^{a,s} G^{b,s}] (the 1/2 of the spin-resolved 2e term)
    real(dp) function tw_g(a_, b_, s_)
      integer, intent(in) :: a_, b_, s_
      wmo = 0.0_dp
      wmo(1:no(s_),no(s_)+1:) = us(s_)%u(:,:,a_)
      wmo(no(s_)+1:,1:no(s_)) = transpose(us(s_)%u(:,:,a_))
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,tmat(:,:,s_),nbf,wmo,nbf,0.0_dp,scr2,nbf)
      tw_g = sum(scr2*transpose(gamo(:,:,b_,s_)))
    end function tw_g

    subroutine add_xc_terms_uhf()
      use mod_dft, only: dft_initialize, dftclean, dftexcor
      use mod_dft_molgrid, only: dft_grid_t
      use mod_dft_gridint_fxc, only: utddft_fxc
      type(dft_grid_t) :: mg
      real(dp), parameter :: hxc = 1.0e-3_dp
      real(dp), allocatable :: mos(:,:,:), dmo(:,:,:), fra(:), frb(:), vxc(:,:,:,:)
      real(dp), allocatable :: fxa(:,:,:), fxb(:,:,:), dxa(:,:,:), dxb(:,:,:), pts(:,:,:,:)
      real(dp) :: sab(3,3,2), exr, telr, tknr, sgn
      integer :: is, xx, kx, cx, ipx, aa, bb, ss

      allocate(mos(nbf,nbf,2), dmo(nbf,nbf,2), fra(nbf2), frb(nbf2), vxc(nbf,nbf,2,2), &
               fxa(nbf,nbf,3), fxb(nbf,nbf,3), dxa(nbf,nbf,3), dxb(nbf,nbf,3), pts(nbf,nbf,3,2))
      do xx = 1, ncart
        cx = mod(xx-1,3) + 1; kx = (xx-1)/3 + 1
        do ss = 1, 2
          call build_t(xx, ss)
          call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mo(:,:,ss),nbf,tmat(:,:,ss),nbf,0.0_dp,dmo(:,:,ss),nbf)
        end do
        do is = 1, 2
          sgn = merge(1.0_dp, -1.0_dp, is == 1)
          basis%atoms%xyz(cx,kx) = basis%atoms%xyz(cx,kx) + sgn*hxc
          call basis%init_shell_centers()
          call dft_initialize(infos, basis, mg)
          mos = mo + sgn*hxc*dmo
          fra = 0.0_dp; frb = 0.0_dp
          call dftexcor(basis, mg, int(infos%control%scftype), fra, frb, mos(:,:,1), mos(:,:,2), &
                        nbf, nbf2, exr, telr, tknr, infos)
          call unpack_from_packed(fra, vxc(:,:,1,is), nbf)
          call unpack_from_packed(frb, vxc(:,:,2,is), nbf)
          do ss = 1, 2
            do bb = 1, 3
              wmo = 0.0_dp
              wmo(1:no(ss),no(ss)+1:) = us(ss)%u(:,:,bb)
              wmo(no(ss)+1:,1:no(ss)) = transpose(us(ss)%u(:,:,bb))
              call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mos(:,:,ss),nbf,wmo,nbf,0.0_dp,scr2,nbf)
              call dgemm('n','t',nbf,nbf,nbf,1.0_dp,scr2,nbf,mos(:,:,ss),nbf,0.0_dp,pts(:,:,bb,ss),nbf)
            end do
          end do
          dxa = pts(:,:,:,1); dxb = pts(:,:,:,2)
          fxa = 0.0_dp; fxb = 0.0_dp
          call utddft_fxc(basis=basis, molGrid=mg, isVecs=.true., wfa=mos(:,:,1), wfb=mos(:,:,2), &
                          fxa=fxa, fxb=fxb, dxa=dxa, dxb=dxb, nmtx=3, threshold=0.0_dp, infos=infos)
          do bb = 1, 3
            do aa = 1, 3
              sab(aa,bb,is) = sum(pts(:,:,aa,1)*fxa(:,:,bb)) + sum(pts(:,:,aa,2)*fxb(:,:,bb))
            end do
          end do
          call dftclean(infos)
          basis%atoms%xyz(cx,kx) = basis%atoms%xyz(cx,kx) - sgn*hxc
          call basis%init_shell_centers()
        end do
        ipx = 0
        do bb = 1, 3
          do aa = 1, bb
            ipx = ipx + 1
            dpol(aa,bb,xx) = dpol(aa,bb,xx) - 2.0_dp*( &
                sum((vxc(:,:,1,1) - vxc(:,:,1,2))*wsym(:,:,ipx,1))/(2.0_dp*hxc) &
              + sum((vxc(:,:,2,1) - vxc(:,:,2,2))*wsym(:,:,ipx,2))/(2.0_dp*hxc) &
              + 0.25_dp*(sab(aa,bb,1) + sab(bb,aa,1) - sab(aa,bb,2) - sab(bb,aa,2))/(2.0_dp*hxc) )
            dpol(bb,aa,xx) = dpol(aa,bb,xx)
          end do
        end do
      end do
      deallocate(mos, dmo, fra, frb, vxc, fxa, fxb, dxa, dxb, pts)
    end subroutine add_xc_terms_uhf
  end subroutine hf_polder_uhf

!###############################################################################

!> @brief Analytic nuclear derivatives of the static polarizability, ROHF/ROKS.
!> @details The ROHF response lives in the docc/socc/virt rotation space
!>   theta (rohf_pack_trial layout); the CPHF operator of cphf_apbx_rohf is
!>   the UHF-form operator on the embedded spin rotations x^s = unpack(theta),
!>   with the exact commutator [F^s_MO, K]_vo of the non-canonical spin Fock
!>   matrices (K the common antisymmetric generator) instead of
!>   F_vv x - x F_oo. With alpha_ab = -2 mu.theta the Hylleraas functional is
!>     -2 [ mu^a.theta^b + mu^b.theta^a + theta^a.H.theta^b ],
!>     theta^a.H.theta^b = sum_s Tr[F^s_MO M^s_ab] + 1/2 sum_s Tr[Pt^{a,s} G^s[Pt^b]],
!>     M^s_ab = K^b X^{a,s}' - X^{a,s}' K^b,
!>   X^{a,s} carrying x^{a,s} in its vir-occ block. Its derivative along the
!>   common relaxed orbital path dC = C T (T_vo = theta^x including the
!>   socc-docc rotation, symmetric part -S^x/2) uses the same derivative-Fock
!>   contractions as the UHF routine; T^x' F + F T^x replaces the canonical
!>   orbital-energy weights. Kohn-Sham XC pieces are central differences along
!>   the relaxed path with the moving grid.
  subroutine hf_polder_rohf(infos, mo, fa_ao, fb_ao, pa, pb, nocca, noccb, &
                            dsa, dha, uvec, hfscale, dpol)
    use oqp_linalg
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int1, only: multipole_integrals
    use grd1, only: der_dipole_matrix
    use mathlib, only: unpack_matrix, pack_matrix
    use fock_deriv_mod, only: fock_deriv_contract_os
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve_rohf, rohf_pack_trial, rohf_unpack_trial
    use io_constants, only: iw
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: mo(:,:), fa_ao(:,:), fb_ao(:,:), pa(:,:), pb(:,:)
    integer, intent(in) :: nocca, noccb
    real(dp), intent(in) :: dsa(:,:,:,:), dha(:,:,:,:), uvec(:,:), hfscale
    real(dp), intent(out) :: dpol(:,:,:)

    type(basis_set), pointer :: basis
    real(dp), allocatable :: mints(:,:), dfull(:,:,:), mmo(:,:,:), dD(:,:,:,:,:)
    real(dp), allocatable :: bF(:,:), uF(:,:), xa(:,:), xb(:,:)
    real(dp), allocatable :: xmo(:,:,:,:), kmo(:,:,:), xao(:,:,:,:), pta(:,:,:,:)
    real(dp), allocatable :: ga(:,:,:,:), gamo(:,:,:,:), fmo(:,:,:), psp(:,:,:), ptot(:,:)
    real(dp), allocatable :: msym(:,:,:,:), wsym(:,:,:,:), gw(:,:,:,:), gwx(:,:,:), gppx(:,:,:)
    real(dp), allocatable :: tmat(:,:), smo(:,:), scr(:,:), scr2(:,:), wmo(:,:), dps(:,:,:)
    real(dp), allocatable :: dpk(:,:), fpk(:,:), gx(:,:), occ(:,:)
    real(dp) :: origin(3), t1(3,3), t2, t3, alpha(3,3), hyl
    integer :: nbf, nbf2, natom, ncart, ltot, nvira, nvirb, offset, no(2)
    integer :: a, b, s, i, j, k, x, kc, cc, ip, iv

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf; nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz, 2); ncart = 3*natom
    nvira = nbf - nocca; nvirb = nbf - noccb; offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira
    no = [nocca, noccb]
    origin = 0.0_dp
    allocate(psp(nbf,nbf,2), ptot(nbf,nbf), fmo(nbf,nbf,2), occ(nbf,2), source=0.0_dp)
    psp(:,:,1) = pa; psp(:,:,2) = pb; ptot = pa + pb
    occ(1:nocca,1) = 1.0_dp; occ(1:noccb,2) = 1.0_dp
    allocate(scr(nbf,nbf), scr2(nbf,nbf), wmo(nbf,nbf), tmat(nbf,nbf), smo(nbf,nbf))
    call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,fa_ao,nbf,0.0_dp,scr,nbf)
    call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,fmo(:,:,1),nbf)
    call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,fb_ao,nbf,0.0_dp,scr,nbf)
    call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,fmo(:,:,2),nbf)

    ! ---- dipole, field CPHF over the ROHF rotation space ------------------
    allocate(mints(nbf2,19), source=0.0_dp)
    call multipole_integrals(basis, mints, origin, 3)
    allocate(dfull(nbf,nbf,3), mmo(nbf,nbf,3), xa(nvira,nocca), xb(nvirb,noccb))
    allocate(bF(ltot,3), uF(ltot,3), source=0.0_dp)
    do a = 1, 3
      call unpack_matrix(mints(:,a), dfull(:,:,a))
      call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,dfull(:,:,a),nbf,0.0_dp,scr,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,mmo(:,:,a),nbf)
      xa = -mmo(nocca+1:,1:nocca,a)
      xb = -mmo(noccb+1:,1:noccb,a)
      call rohf_pack_trial(bF(:,a), xa, xb, nbf, nocca, noccb)
    end do
    deallocate(mints)
    call cphf_solve_rohf(infos, 3, bF, uF)
    do b = 1, 3
      do a = 1, 3
        alpha(a,b) = 2.0_dp*sum(bF(:,a)*uF(:,b))
      end do
    end do

    ! ---- embedded spin rotations, generators, response Fock ---------------
    allocate(xmo(nbf,nbf,3,2), kmo(nbf,nbf,3), xao(nbf,nbf,3,2), pta(nbf,nbf,3,2), &
             ga(nbf,nbf,3,2), gamo(nbf,nbf,3,2), source=0.0_dp)
    allocate(dpk(nbf2,2), fpk(nbf2,2))
    do a = 1, 3
      call rohf_unpack_trial(uF(:,a), xa, xb, nbf, nocca, noccb)
      xmo(nocca+1:,1:nocca,a,1) = xa
      xmo(noccb+1:,1:noccb,a,2) = xb
      do i = 1, nocca
        do iv = 1, nvira
          kmo(nocca+iv,i,a) = xa(iv,i)
          kmo(i,nocca+iv,a) = -xa(iv,i)
        end do
      end do
      do j = 1, noccb
        do iv = 1, offset
          kmo(noccb+iv,j,a) = kmo(noccb+iv,j,a) + xb(iv,j)
          kmo(j,noccb+iv,a) = kmo(j,noccb+iv,a) - xb(iv,j)
        end do
      end do
      do s = 1, 2
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mo,nbf,xmo(:,:,a,s),nbf,0.0_dp,scr,nbf)
        call dgemm('n','t',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,xao(:,:,a,s),nbf)
        pta(:,:,a,s) = xao(:,:,a,s) + transpose(xao(:,:,a,s))
        call pack_matrix(pta(:,:,a,s), dpk(:,s))
      end do
      fpk = 0.0_dp
      call fock_jk(basis, d=dpk, f=fpk, scale_exch=hfscale, infos=infos)
      do s = 1, 2
        call unpack_from_packed(fpk(:,s), ga(:,:,a,s), nbf)
        call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,ga(:,:,a,s),nbf,0.0_dp,scr,nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,gamo(:,:,a,s),nbf)
      end do
    end do

    ! ---- pair intermediates ------------------------------------------------
    allocate(msym(nbf,nbf,6,2), wsym(nbf,nbf,6,2), gw(nbf,nbf,6,2), gwx(3,natom,6), &
             gppx(3,natom,6), gx(3,natom), source=0.0_dp)
    ip = 0
    do b = 1, 3
      do a = 1, b
        ip = ip + 1
        do s = 1, 2
          ! M^s_ab = K^b X^{a,s}' - X^{a,s}' K^b, symmetrized in a<->b
          call dgemm('n','t',nbf,nbf,nbf,0.5_dp,kmo(:,:,b),nbf,xmo(:,:,a,s),nbf,0.0_dp,scr,nbf)
          call dgemm('t','n',nbf,nbf,nbf,-0.5_dp,xmo(:,:,a,s),nbf,kmo(:,:,b),nbf,1.0_dp,scr,nbf)
          call dgemm('n','t',nbf,nbf,nbf,0.5_dp,kmo(:,:,a),nbf,xmo(:,:,b,s),nbf,1.0_dp,scr,nbf)
          call dgemm('t','n',nbf,nbf,nbf,-0.5_dp,xmo(:,:,b,s),nbf,kmo(:,:,a),nbf,1.0_dp,scr,nbf)
          msym(:,:,ip,s) = scr
          scr2 = 0.5_dp*(scr + transpose(scr))
          call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mo,nbf,scr2,nbf,0.0_dp,scr,nbf)
          call dgemm('n','t',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,wsym(:,:,ip,s),nbf)
          call pack_matrix(wsym(:,:,ip,s), dpk(:,s))
          call fock_deriv_contract_os(infos, basis, ptot, psp(:,:,s), wsym(:,:,ip,s), hfscale, gx)
          gwx(:,:,ip) = gwx(:,:,ip) + gx
          call fock_deriv_contract_os(infos, basis, pta(:,:,b,1) + pta(:,:,b,2), pta(:,:,b,s), &
                                      pta(:,:,a,s), hfscale, gx)
          gppx(:,:,ip) = gppx(:,:,ip) + 0.5_dp*gx
        end do
        fpk = 0.0_dp
        call fock_jk(basis, d=dpk, f=fpk, scale_exch=hfscale, infos=infos)
        do s = 1, 2
          call unpack_from_packed(fpk(:,s), gw(:,:,ip,s), nbf)
        end do
      end do
    end do

    if (infos%control%verbose >= 2) then
      ip = 0
      do b = 1, 3
        do a = 1, b
          ip = ip + 1
          hyl = 0.0_dp
          do s = 1, 2
            hyl = hyl + sum(fmo(:,:,s)*msym(:,:,ip,s)) + 0.5_dp*sum(pta(:,:,a,s)*ga(:,:,b,s))
          end do
          write(iw,'(A,2I2,2ES16.8)') '  polder(ROHF) Hylleraas check a b, alpha, diff:', a, b, alpha(a,b), &
            alpha(a,b) - (2.0_dp*(sum(bF(:,a)*uF(:,b)) + sum(bF(:,b)*uF(:,a))) - 2.0_dp*hyl)
        end do
      end do
    end if

    allocate(dD(nbf,nbf,3,natom,3))
    call der_dipole_matrix(basis, origin, dD)
    do a = 1, 3
      do k = 1, natom
        do cc = 1, 3
          do j = 1, nbf
            dD(:,j,cc,k,a) = dD(:,j,cc,k,a)*basis%bfnrm(:)*basis%bfnrm(j)
          end do
        end do
      end do
    end do

    ! ---- per nuclear coordinate ------------------------------------------
    allocate(dps(nbf,nbf,2))
    dpol = 0.0_dp
    do x = 1, ncart
      cc = mod(x-1,3) + 1; kc = (x-1)/3 + 1
      call build_t(x)
      ! relaxed spin-density derivatives dP^s = C (T O^s + O^s T') C'
      do s = 1, 2
        do j = 1, nbf
          scr2(:,j) = tmat(:,j)*occ(j,s)
        end do
        scr2 = scr2 + transpose(scr2)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mo,nbf,scr2,nbf,0.0_dp,scr,nbf)
        call dgemm('n','t',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf,0.0_dp,dps(:,:,s),nbf)
      end do
      t1 = 0.0_dp
      do a = 1, 3
        call dgemm('t','n',nbf,nbf,nbf,1.0_dp,tmat,nbf,mmo(:,:,a),nbf,0.0_dp,scr,nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mmo(:,:,a),nbf,tmat,nbf,1.0_dp,scr,nbf)
        do b = 1, 3
          do s = 1, 2
            t1(a,b) = t1(a,b) + sum(dD(:,:,cc,kc,a)*xao(:,:,b,s)) + sum(xmo(:,:,b,s)*scr)
          end do
        end do
      end do
      ip = 0
      do b = 1, 3
        do a = 1, b
          ip = ip + 1
          t2 = gwx(cc,kc,ip)
          t3 = gppx(cc,kc,ip)
          do s = 1, 2
            ! (T'F + F T) contracted with M^s
            call dgemm('t','n',nbf,nbf,nbf,1.0_dp,tmat,nbf,fmo(:,:,s),nbf,0.0_dp,scr,nbf)
            call dgemm('n','n',nbf,nbf,nbf,1.0_dp,fmo(:,:,s),nbf,tmat,nbf,1.0_dp,scr,nbf)
            t2 = t2 + sum(dha(:,:,cc,kc)*wsym(:,:,ip,s)) + sum(dps(:,:,s)*gw(:,:,ip,s)) &
                    + sum(scr*msym(:,:,ip,s))
            t3 = t3 + tw_g(a, b, s) + tw_g(b, a, s)
          end do
          dpol(a,b,x) = -2.0_dp*(t1(a,b) + t1(b,a) + t2 + t3)
          dpol(b,a,x) = dpol(a,b,x)
        end do
      end do
    end do

    if (infos%control%hamilton == 20) call add_xc_terms_rohf()

    deallocate(psp, ptot, fmo, occ, scr, scr2, wmo, tmat, smo, dfull, mmo, xa, xb, bF, uF, &
               xmo, kmo, xao, pta, ga, gamo, dpk, fpk, msym, wsym, gw, gwx, gppx, gx, dD, dps)

  contains

    !> common relaxed orbital-path generator T^x of the ROHF orbitals
    subroutine build_t(x_)
      integer, intent(in) :: x_
      integer :: cx_, kx_
      cx_ = mod(x_-1,3) + 1; kx_ = (x_-1)/3 + 1
      call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,dsa(:,:,cx_,kx_),nbf,0.0_dp,scr2,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr2,nbf,mo,nbf,0.0_dp,smo,nbf)
      call rohf_unpack_trial(uvec(:,x_), xa, xb, nbf, nocca, noccb)
      tmat = -0.5_dp*smo
      tmat(nocca+1:,1:nocca) = xa
      tmat(1:nocca,nocca+1:) = -transpose(xa) - smo(1:nocca,nocca+1:)
      if (offset > 0) then
        tmat(noccb+1:nocca,1:noccb) = xb(1:offset,:)
        tmat(1:noccb,noccb+1:nocca) = -transpose(xb(1:offset,:)) - smo(1:noccb,noccb+1:nocca)
      end if
    end subroutine build_t

    !> Tr[T W^{a,s} G^{b,s}] (W^{a,s} = X^{a,s} + X^{a,s}' in the MO basis)
    real(dp) function tw_g(a_, b_, s_)
      integer, intent(in) :: a_, b_, s_
      wmo = xmo(:,:,a_,s_) + transpose(xmo(:,:,a_,s_))
      call dgemm('n','n',nbf,nbf,nbf,1.0_dp,tmat,nbf,wmo,nbf,0.0_dp,scr2,nbf)
      tw_g = sum(scr2*transpose(gamo(:,:,b_,s_)))
    end function tw_g

    subroutine add_xc_terms_rohf()
      use mod_dft, only: dft_initialize, dftclean, dftexcor
      use mod_dft_molgrid, only: dft_grid_t
      use mod_dft_gridint_fxc, only: utddft_fxc
      type(dft_grid_t) :: mg
      real(dp), parameter :: hxc = 1.0e-3_dp
      real(dp), allocatable :: mos(:,:), dmo(:,:), fra(:), frb(:), vxc(:,:,:,:)
      real(dp), allocatable :: fxa(:,:,:), fxb(:,:,:), dxa(:,:,:), dxb(:,:,:), pts(:,:,:,:)
      real(dp) :: sab(3,3,2), exr, telr, tknr, sgn
      integer :: is, xx, kx, cx, ipx, aa, bb, ss

      allocate(mos(nbf,nbf), dmo(nbf,nbf), fra(nbf2), frb(nbf2), vxc(nbf,nbf,2,2), &
               fxa(nbf,nbf,3), fxb(nbf,nbf,3), dxa(nbf,nbf,3), dxb(nbf,nbf,3), pts(nbf,nbf,3,2))
      do xx = 1, ncart
        cx = mod(xx-1,3) + 1; kx = (xx-1)/3 + 1
        call build_t(xx)
        call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mo,nbf,tmat,nbf,0.0_dp,dmo,nbf)
        do is = 1, 2
          sgn = merge(1.0_dp, -1.0_dp, is == 1)
          basis%atoms%xyz(cx,kx) = basis%atoms%xyz(cx,kx) + sgn*hxc
          call basis%init_shell_centers()
          call dft_initialize(infos, basis, mg)
          mos = mo + sgn*hxc*dmo
          fra = 0.0_dp; frb = 0.0_dp
          call dftexcor(basis, mg, int(infos%control%scftype), fra, frb, mos, mos, &
                        nbf, nbf2, exr, telr, tknr, infos)
          call unpack_from_packed(fra, vxc(:,:,1,is), nbf)
          call unpack_from_packed(frb, vxc(:,:,2,is), nbf)
          do ss = 1, 2
            do bb = 1, 3
              wmo = xmo(:,:,bb,ss) + transpose(xmo(:,:,bb,ss))
              call dgemm('n','n',nbf,nbf,nbf,1.0_dp,mos,nbf,wmo,nbf,0.0_dp,scr2,nbf)
              call dgemm('n','t',nbf,nbf,nbf,1.0_dp,scr2,nbf,mos,nbf,0.0_dp,pts(:,:,bb,ss),nbf)
            end do
          end do
          dxa = pts(:,:,:,1); dxb = pts(:,:,:,2)
          fxa = 0.0_dp; fxb = 0.0_dp
          call utddft_fxc(basis=basis, molGrid=mg, isVecs=.true., wfa=mos, wfb=mos, &
                          fxa=fxa, fxb=fxb, dxa=dxa, dxb=dxb, nmtx=3, threshold=0.0_dp, infos=infos)
          do bb = 1, 3
            do aa = 1, 3
              sab(aa,bb,is) = sum(pts(:,:,aa,1)*fxa(:,:,bb)) + sum(pts(:,:,aa,2)*fxb(:,:,bb))
            end do
          end do
          call dftclean(infos)
          basis%atoms%xyz(cx,kx) = basis%atoms%xyz(cx,kx) - sgn*hxc
          call basis%init_shell_centers()
        end do
        ipx = 0
        do bb = 1, 3
          do aa = 1, bb
            ipx = ipx + 1
            dpol(aa,bb,xx) = dpol(aa,bb,xx) - 2.0_dp*( &
                sum((vxc(:,:,1,1) - vxc(:,:,1,2))*wsym(:,:,ipx,1))/(2.0_dp*hxc) &
              + sum((vxc(:,:,2,1) - vxc(:,:,2,2))*wsym(:,:,ipx,2))/(2.0_dp*hxc) &
              + 0.25_dp*(sab(aa,bb,1) + sab(bb,aa,1) - sab(aa,bb,2) - sab(bb,aa,2))/(2.0_dp*hxc) )
            dpol(bb,aa,xx) = dpol(aa,bb,xx)
          end do
        end do
      end do
      deallocate(mos, dmo, fra, frb, vxc, fxa, fxb, dxa, dxb, pts)
    end subroutine add_xc_terms_rohf
  end subroutine hf_polder_rohf

!###############################################################################

!> @brief Whether the caller wants the IR/Raman property derivatives.
!> @details Opt-in: only a nonzero OQP::hess_properties computes them. The
!>   Python ground-state driver sets it (0 for matrix-only native TS/IRC or
!>   analysis=False). Callers that never set it, such as the TD Hessian's
!>   internal ground-state call, skip the field-response solves.
  logical function hf_hess_properties_wanted(infos) result(want)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_hess_properties, ta_ok
    use iso_c_binding, only: c_int64_t
    type(information), target, intent(inout) :: infos
    integer(c_int64_t), contiguous, pointer :: flag(:)
    integer(4) :: status
    want = .false.
    call tagarray_get_data(infos%dat, OQP_hess_properties, flag, status)
    if (status == ta_ok) then
      if (size(flag) > 0) want = flag(1) /= 0
    end if
  end function hf_hess_properties_wanted

!###############################################################################

!> @brief Nuclear-coordinate derivatives of the electric dipole moment from the
!>        relaxed density derivatives of the analytic Hessian.
!> @details  dmu_a/dR_x = Z_eff(A) delta_{a,c(x)}
!>                       - Tr[P dD_a/dR_x] - Tr[dP/dR_x D_a],
!>   with D_a the AO dipole integrals about the origin used by
!>   electric_dipole_au (so the result is the derivative of that dipole, also
!>   for ions) and dP/dR_x the total relaxed density derivative the CPHF step
!>   already built. No further response solve is needed.
!>   hf_dipder_init returns the normalized D_a and dD_a/dR and the static
!>   (nuclear + integral-derivative) part; hf_dipder_add_response subtracts
!>   Tr[dP^x D_a] for one coordinate; hf_dipder_store saves (3,3N) in
!>   OQP::hf_dipole_derivatives.
  subroutine hf_dipder_init(infos, ptot, dfull, dmu)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int1, only: multipole_integrals
    use grd1, only: der_dipole_matrix
    use mathlib, only: unpack_matrix
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: ptot(:,:)
    real(dp), allocatable, intent(out) :: dfull(:,:,:), dmu(:,:)

    type(basis_set), pointer :: basis
    real(dp), allocatable :: mints(:,:), dD(:,:,:,:,:)
    real(dp) :: origin(3), z
    integer :: nbf, natom, a, c, k, x, mu, nu

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    natom = size(basis%atoms%xyz, 2)
    origin = 0.0_dp

    allocate(mints(nbf*(nbf+1)/2,19), source=0.0_dp)
    call multipole_integrals(basis, mints, origin, 3)
    allocate(dfull(nbf,nbf,3))
    do a = 1, 3
      call unpack_matrix(mints(:,a), dfull(:,:,a))
    end do
    deallocate(mints)

    allocate(dD(nbf,nbf,3,natom,3))
    call der_dipole_matrix(basis, origin, dD)
    allocate(dmu(3,3*natom), source=0.0_dp)
    do a = 1, 3
      do k = 1, natom
        z = basis%atoms%zn(k) - basis%ecp_zn_num(k)
        do c = 1, 3
          x = 3*(k-1) + c
          do nu = 1, nbf
            do mu = 1, nbf
              dmu(a,x) = dmu(a,x) - ptot(mu,nu)*dD(mu,nu,c,k,a) &
                                    *basis%bfnrm(mu)*basis%bfnrm(nu)
            end do
          end do
          if (a == c) dmu(a,x) = dmu(a,x) + z
        end do
      end do
    end do
    deallocate(dD)
  end subroutine hf_dipder_init

  subroutine hf_dipder_add_response(dfull, dptot, dmu_x)
    use precision, only: dp
    real(dp), intent(in) :: dfull(:,:,:), dptot(:,:)
    real(dp), intent(inout) :: dmu_x(3)
    integer :: a
    do a = 1, 3
      dmu_x(a) = dmu_x(a) - sum(dptot*dfull(:,:,a))
    end do
  end subroutine hf_dipder_add_response

  subroutine hf_dipder_store(infos, dmu)
    use precision, only: dp
    use types, only: information
    use oqp_tagarray_driver, only: OQP_hf_dipole_derivatives
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: dmu(:,:)
    real(dp), contiguous, pointer :: store(:,:)
    call infos%dat%alloc_or_die(OQP_hf_dipole_derivatives, (/ 3, size(dmu,2) /), store, &
      description='Analytic nuclear derivatives of the electric dipole (a.u.), (3,3N)')
    store = dmu
  end subroutine hf_dipder_store

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


!###############################################################################

!> @brief Opt-in cross-check of the blocked operator path for the RHF 2e
!>        response skeleton: OQP_HESS_G2E_CHECK=1 recomputes a few occ-vir
!>        pairs with the per-pair fock_deriv_contract and prints both.
  subroutine check_g2e_operator(infos, basis, mo_a, pfull, hfscale, nocc, g2e)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use fock_deriv_mod, only: fock_deriv_contract
    use io_constants, only: iw
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: mo_a(:,:), pfull(:,:), g2e(:,:), hfscale
    integer, intent(in) :: nocc
    real(kind=dp), allocatable :: probe(:,:), gx(:,:)
    real(kind=dp), allocatable :: ref(:)
    character(len=8) :: envs
    integer :: st, nbf, nvir, ncart, a, i, k, mu, nu, ia
    real(kind=dp) :: num, den
    call get_environment_variable('OQP_HESS_G2E_CHECK', envs, status=st)
    if (st /= 0) return
    if (trim(adjustl(envs)) /= '1') return
    nbf = size(mo_a,1); nvir = nbf - nocc; ncart = size(g2e,2)
    allocate(probe(nbf,nbf), gx(3,ncart/3), ref(ncart))
    do k = 0, 2
      a = 1 + k*(nvir-1)/2
      i = nocc - k*(nocc-1)/2
      do mu = 1, nbf
        do nu = 1, nbf
          probe(mu,nu) = 0.5_dp*( mo_a(mu,nocc+a)*mo_a(nu,i) + mo_a(mu,i)*mo_a(nu,nocc+a) )
        end do
      end do
      call fock_deriv_contract(infos, basis, pfull, probe, hfscale, gx)
      ref = reshape(gx, [ncart])
      ia = (a-1)*nocc + i
      num = dot_product(ref, g2e(ia,:)); den = dot_product(g2e(ia,:), g2e(ia,:))
      write(iw,'(A,2I5,A,ES12.4,A,ES12.4,A,ES18.10)') '  g2e check i,a=', i, a, &
        '  max|ref|=', maxval(abs(ref)), '  max|ref-op|=', maxval(abs(ref-g2e(ia,:))), &
        '  ref/op=', num/max(den, tiny(1.0_dp))
    end do
  end subroutine check_g2e_operator


!###############################################################################

!> @brief The blocked operator path needs one MPI rank and no attenuated
!>        (range-separated) exchange pass.
  logical function os_operator_path(infos) result(ok)
    use types, only: information
    type(information), intent(in) :: infos
    ok = .not. infos%dft%cam_flag .and. .not. infos%mpiinfo%usempi
  end function os_operator_path

!###############################################################################

!> @brief Open-shell 2e response skeleton g2e_s(ia,x) = Tr[probe_ia G^{s,x}]
!>        with G^{s,x} = J^x[Ptot] - c_x K^x[P^s], from three blocked
!>        derivative-ERI traversals (J on Ptot, K on Pa and on Pb) assembled
!>        as AO operators and transformed to each spin's occ-vir MO block.
!>        Layout of g2ea/g2eb: ((a-1)*nocc+i, x).
  subroutine os_g2e_from_operators(infos, basis, ptot, pa, pb, hfscale, &
                                   moa, nocca, mob, noccb, g2ea, g2eb, gao)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use tdhf_hessian_z_rhs_mod, only: eri_derivative_operator_mo
    use oqp_linalg
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: ptot(:,:), pa(:,:), pb(:,:), hfscale
    real(kind=dp), intent(in) :: moa(:,:), mob(:,:)
    integer, intent(in) :: nocca, noccb
    real(kind=dp), allocatable, intent(out) :: g2ea(:,:), g2eb(:,:)
    !> optional: keep the AO operators G^{s,x} (s = alpha, beta), normalized so
    !> that fock_deriv_contract_os(Ptot, P^s, M) = Tr[M G^{s,x}]
    real(kind=dp), allocatable, intent(out), optional :: gao(:,:,:,:)
    real(kind=dp), allocatable :: eye(:,:), jao(:,:,:), kao(:,:,:)
    integer :: nbf, ncart, k
    nbf = size(moa,1); ncart = 3*size(infos%atoms%xyz,2)
    allocate(eye(nbf,nbf), source=0.0_dp)
    do k = 1, nbf
      eye(k,k) = 1.0_dp
    end do
    allocate(jao(nbf,nbf,ncart), kao(nbf,nbf,ncart))
    ! Coulomb only (exchange scale 0) on the total density
    call eri_derivative_operator_mo(infos, eye, ptot, 1, 0.0_dp, jao)
    jao = OS_J_SCALE*jao
    if (present(gao)) allocate(gao(nbf,nbf,ncart,2))
    call spin_block(pa, moa, nocca, g2ea, 1)
    call spin_block(pb, mob, noccb, g2eb, 2)
    call check_os_g2e_operator(infos, basis, ptot, pa, moa, nocca, hfscale, g2ea)
    deallocate(eye, jao, kao)
  contains
    subroutine spin_block(ps, mo, nocc, g2e, ispin)
      real(kind=dp), intent(in) :: ps(:,:), mo(:,:)
      integer, intent(in) :: nocc, ispin
      real(kind=dp), allocatable, intent(out) :: g2e(:,:)
      real(kind=dp), allocatable :: g(:,:), t(:,:), gov(:,:)
      integer :: nvir, x, a, i
      nvir = nbf - nocc
      allocate(g2e(nocc*nvir,ncart), g(nbf,nbf), t(nbf,nvir), gov(nocc,nvir))
      if (hfscale /= 0.0_dp) then
        ! exchange only (Coulomb scale 0) on the spin density
        call eri_derivative_operator_mo(infos, eye, ps, 1, hfscale, kao, coulscale=0.0_dp)
      else
        kao = 0.0_dp
      end if
      do x = 1, ncart
        g = jao(:,:,x) + OS_K_SCALE*kao(:,:,x)
        if (present(gao)) gao(:,:,x,ispin) = g
        call dgemm('n','n',nbf,nvir,nbf,1.0_dp,g,nbf,mo(:,nocc+1:),nbf,0.0_dp,t,nbf)
        call dgemm('t','n',nocc,nvir,nbf,1.0_dp,mo(:,1:nocc),nbf,t,nbf,0.0_dp,gov,nocc)
        do a = 1, nvir
          do i = 1, nocc
            g2e((a-1)*nocc+i,x) = gov(i,a)
          end do
        end do
      end do
    end subroutine spin_block
  end subroutine os_g2e_from_operators

!###############################################################################

!> @brief Opt-in cross-check (OQP_HESS_G2E_CHECK=1) of the open-shell operator
!>        path against fock_deriv_contract_os for three alpha occ-vir pairs.
  subroutine check_os_g2e_operator(infos, basis, ptot, pa, mo, nocc, hfscale, g2e)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use fock_deriv_mod, only: fock_deriv_contract_os
    use io_constants, only: iw
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: ptot(:,:), pa(:,:), mo(:,:), hfscale, g2e(:,:)
    integer, intent(in) :: nocc
    real(kind=dp), allocatable :: probe(:,:), gx(:,:), ref(:)
    character(len=8) :: envs
    integer :: st, nbf, nvir, ncart, a, i, k, mu, nu, ia
    real(kind=dp) :: num, den
    call get_environment_variable('OQP_HESS_G2E_CHECK', envs, status=st)
    if (st /= 0) return
    if (trim(adjustl(envs)) /= '1') return
    nbf = size(mo,1); nvir = nbf - nocc; ncart = size(g2e,2)
    allocate(probe(nbf,nbf), gx(3,ncart/3), ref(ncart))
    do k = 0, 2
      a = 1 + k*(nvir-1)/2
      i = nocc - k*(nocc-1)/2
      do mu = 1, nbf
        do nu = 1, nbf
          probe(mu,nu) = 0.5_dp*( mo(mu,nocc+a)*mo(nu,i) + mo(mu,i)*mo(nu,nocc+a) )
        end do
      end do
      call fock_deriv_contract_os(infos, basis, ptot, pa, probe, hfscale, gx)
      ref = reshape(gx, [ncart])
      ia = (a-1)*nocc + i
      num = dot_product(ref, g2e(ia,:)); den = dot_product(g2e(ia,:), g2e(ia,:))
      write(iw,'(A,2I5,A,ES12.4,A,ES12.4,A,ES18.10)') '  os g2e check i,a=', i, a, &
        '  max|ref|=', maxval(abs(ref)), '  max|ref-op|=', maxval(abs(ref-g2e(ia,:))), &
        '  ref/op=', num/max(den, tiny(1.0_dp))
    end do
    ! separate Coulomb / exchange ratios for calibration (one pair)
    call os_part_ratio('J', ptot, 0.0_dp*pa, 1.0_dp, 0.0_dp)
    call os_part_ratio('K', 0.0_dp*ptot, pa, 0.0_dp, hfscale)
  contains
    subroutine os_part_ratio(tag, pc, px, cs, hs)
      use tdhf_hessian_z_rhs_mod, only: eri_derivative_operator_mo
      character(len=*), intent(in) :: tag
      real(kind=dp), intent(in) :: pc(:,:), px(:,:), cs, hs
      real(kind=dp), allocatable :: op(:,:,:), v(:)
      if (hs == 0.0_dp .and. cs == 0.0_dp) return
      a = 1; i = nocc
      do mu = 1, nbf
        do nu = 1, nbf
          probe(mu,nu) = 0.5_dp*( mo(mu,nocc+a)*mo(nu,i) + mo(mu,i)*mo(nu,nocc+a) )
        end do
      end do
      call fock_deriv_contract_os(infos, basis, pc, px, probe, hfscale, gx)
      ref = reshape(gx, [ncart])
      allocate(op(nbf,nbf,ncart), v(ncart))
      if (cs /= 0.0_dp) then
        call eri_derivative_operator_mo(infos, mo, pc, 1, 0.0_dp, op)
      else
        call eri_derivative_operator_mo(infos, mo, px, 1, hs, op, coulscale=0.0_dp)
      end if
      v = op(i,nocc+a,:)
      write(iw,'(A,A,A,ES18.10,A,ES12.4)') '  os g2e part ', tag, '  ref/op=', &
        dot_product(ref,v)/max(dot_product(v,v), tiny(1.0_dp)), '  max|ref|=', maxval(abs(ref))
    end subroutine os_part_ratio
  end subroutine check_os_g2e_operator


!###############################################################################

!> @brief Opt-in cross-check (OQP_HESS_G2E_CHECK=1) of one 2e-trace column
!>        built from the stored operator against fock_deriv_contract.
  subroutine check_trace_operator(infos, basis, pfull, hfscale, m, col, tag)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use fock_deriv_mod, only: fock_deriv_contract
    use io_constants, only: iw
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: pfull(:,:), hfscale, m(:,:), col(:)
    character(len=*), intent(in) :: tag
    real(kind=dp), allocatable :: gx(:,:), ref(:)
    character(len=8) :: envs
    integer :: st
    call get_environment_variable('OQP_HESS_G2E_CHECK', envs, status=st)
    if (st /= 0) return
    if (trim(adjustl(envs)) /= '1') return
    allocate(gx(3,size(col)/3), ref(size(col)))
    call fock_deriv_contract(infos, basis, pfull, m, hfscale, gx)
    ref = 2.0_dp*reshape(gx, [size(col)])
    write(iw,'(A,A,A,ES12.4,A,ES12.4)') '  trace check ', tag, '  max|ref|=', maxval(abs(ref)), &
      '  max|ref-op|=', maxval(abs(ref-col))
  end subroutine check_trace_operator


!###############################################################################

!> @brief Opt-in wall-clock section timer (OQP_HESS_TIMERS=1): prints the time
!>        since the previous call; label 'start' only resets the clock.
  subroutine hess_tick(label, logfile)
    use io_constants, only: iw
    character(len=*), intent(in) :: label
    !> log to append to when unit iw is not open (callers outside hf_hessian)
    character(len=*), intent(in), optional :: logfile
    integer(8), save :: last = -1_8
    integer(8) :: now, rate
    character(len=8) :: envs
    integer :: st
    logical :: opened, mine
    call get_environment_variable('OQP_HESS_TIMERS', envs, status=st)
    if (st /= 0) return
    if (trim(adjustl(envs)) /= '1') return
    call system_clock(now, rate)
    if (label /= 'start' .and. last >= 0_8) then
      inquire(unit=iw, opened=opened)
      mine = .not. opened .and. present(logfile)
      if (mine) open(unit=iw, file=logfile, position='append')
      if (opened .or. mine) write(iw,'(A,A40,F10.3,A)') '  hess timer: ', label, &
        real(now-last,8)/real(rate,8), ' s'
      if (mine) close(iw)
    end if
    last = now
  end subroutine hess_tick

end module hf_hessian_mod
