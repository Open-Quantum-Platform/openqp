!> @file qmrsf_dk_paper_native.F90
!> @brief Paper-faithful QMRSF-DK implementation using only OpenQP data.
!>
!> Quintet mixed-reference spin-flip TDDFT with a dressed exchange-correlation
!> kernel.  The module obtains the converged ROKS/ROHF orbitals, spin Fock
!> matrices and spin-resolved v_xc directly from OpenQP, transforms the
!> four-SOMO active integrals with qmrsf_ao2mo_mod, and constructs the spin-pure
!> CAS(4,4) response blocks with the analytic 20-singlet / 15-triplet /
!> 1-quintet transformation.
module qmrsf_dk_paper_native_mod
  use precision, only: dp
  implicit none
  private
  public :: qmrsf_dk_paper_native

  integer, parameter :: NACT = 4
  integer, parameter :: NSO = 8
  integer, parameter :: NDET = 36
  integer, parameter :: NSING = 20
  integer, parameter :: NTRIP = 15

  !> Seam conventions for the exchange tensor that carries c_H; see
  !> @ref native_seam_mask.  The active pair whose orientation the value-based
  !> seam is sensitive to is the middle one of the four energy-ordered SOMOs.
  integer, parameter :: SEAM_NATIVE = 0
  integer, parameter :: SEAM_S3R    = 1
  integer, parameter :: SEAM_HAAR   = 2
  integer, parameter :: NSEAM = 1
  !> Value that both [dftgrid] and [tdhf] use for an unspecified CAM parameter.
  real(dp), parameter :: CAM_UNSET = -1.0_dp
  character(len=6), parameter :: SEAM_NAME(0:NSEAM-1) = (/ 'native' /)
  !> An adjacent orbital pair is treated as a doublet when it is degenerate to
  !> within PAIR_GAP_THRESHOLD, which is loose enough to cover the
  !> distance-degenerate frontier pairs of a separated dimer and tight enough to
  !> exclude genuinely split orbitals.  When the reference carries no such
  !> degeneracy the frontier pair is used instead, so that the construction is
  !> continuous as a degeneracy is approached and does not switch between
  !> subspaces along a potential energy surface.
  real(dp), parameter :: PAIR_GAP_THRESHOLD = 1.0e-2_dp

contains

  subroutine qmrsf_dk_paper_native(infos)
    use types, only: information
    use io_constants, only: iw
    use printing, only: print_module_info
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A, OQP_E_MO_A
    use qmrsf_ao2mo_mod, only: qmrsf_active_integrals
    use eigen, only: diag_symm_full
    use oqp_linalg

    type(information), target, intent(inout) :: infos

    real(dp), contiguous, pointer :: mo_a(:,:)
    real(dp), allocatable :: h_act(:,:), eri_act(:,:,:,:), cact(:,:)
    real(dp), allocatable :: eri_lr(:,:,:,:)
    real(dp) :: va_act(NACT,NACT), vb_act(NACT,NACT)
    real(dp) :: h_eff(NACT,NACT), kcore_act(NACT,NACT), kcore_lr(NACT,NACT)
    real(dp) :: hzero(NACT,NACT), eri_k(NACT,NACT,NACT,NACT)
    real(dp) :: eri_klr(NACT,NACT,NACT,NACT)
    real(dp) :: h1klr(NSO,NSO), gklr(NSO,NSO,NSO,NSO), hklrdet(NDET,NDET)
    real(dp) :: a_ref, b_ref, a_h, b_h, mu_ref, mu_h
    logical  :: is_cam
    real(dp) :: h1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp) :: h1k(NSO,NSO), gk(NSO,NSO,NSO,NSO)
    real(dp) :: hdet(NDET,NDET), hkdet(NDET,NDET), ucsf(NDET,NDET), hcsf(NDET,NDET)
    real(dp) :: hbare(NDET,NDET), work36(NDET,NDET)
    real(dp) :: asing(NSING,NSING), atrip(NTRIP,NTRIP), aquint(1,1)
    real(dp) :: asing_keep(NSING,NSING), atrip_keep(NTRIP,NTRIP)
    real(dp) :: asing_all(NSING,NSING,0:NSEAM-1), atrip_all(NTRIP,NTRIP,0:NSEAM-1)
    real(dp) :: esing(NSING), etrip(NTRIP), equint(1)
    real(dp) :: esing_all(NSING,0:NSEAM-1), etrip_all(NTRIP,0:NSEAM-1)
    real(dp) :: equint_all(0:NSEAM-1)
    real(dp) :: ecore, c_h, c_ref, orth_err, cross_err, vxc_max
    real(dp) :: eref, thresh, qvec(1,1)
    real(dp) :: diag36v(NDET), dvcsf(NDET)
    real(dp) :: pair_thc, pair_aniso
    real(dp) :: eps_act(NACT)
    real(dp), contiguous, pointer :: e_mo(:)
    integer :: npair, pairs(2,2)
    character(len=8) :: csf_label(NDET)
    integer :: dets(4,NDET), act(NACT)
    integer :: nbf, ncore, i, j, ierr, iseam
    logical :: is_dft

    open(unit=iw, file=infos%log_filename, position='append')
    call print_module_info('QMRSF-DK', &
         'Quintet-reference MRSF-TDDFT with a dressed exchange-correlation kernel')

    nbf = int(infos%basis%nbf)
    ncore = int(infos%mol_prop%nelec_B)
    if (int(infos%mol_prop%nelec_A) - ncore /= NACT) then
      write(iw,'(/,5x,a)') 'QMRSF-DK requires a high-spin quintet reference '// &
           '(N_alpha - N_beta = 4).'
      call flush(iw)
      close(iw)
      return
    end if

    do i = 1, NACT
      act(i) = ncore + i
    end do

    is_dft = (infos%control%hamilton == 20)
    is_cam = is_dft .and. infos%dft%cam_flag
    c_ref = 1.0_dp
    if (is_dft) c_ref = infos%dft%hfscale
    c_h = infos%tddft%hfscale
    if (c_h < 0.0_dp) c_h = c_ref

    ! A range-separated functional treats alpha + beta*erf(mu r) of the
    ! interaction by exact exchange, so its exchange operator is
    ! alpha*K + beta*K_lr rather than a single scaled K.  The reference values
    ! come from [dftgrid]; the response may override them through [tdhf], with
    ! -1 meaning "inherit", as in the MRSF response.
    a_ref = c_ref; b_ref = 0.0_dp; a_h = c_h; b_h = 0.0_dp
    mu_ref = 0.0_dp; mu_h = 0.0_dp
    if (is_cam) then
      a_ref  = infos%dft%cam_alpha
      b_ref  = infos%dft%cam_beta
      mu_ref = infos%dft%cam_mu
      ! CAM_UNSET is the "not specified" sentinel of both input sections.  A
      ! genuine beta may be negative -- the DTCAM presets set beta = -0.20 for
      ! the reference and -0.10 for the response -- so only the exact sentinel
      ! may be treated as absent.  The two range parameters are likewise
      ! independent: those presets use mu = 0.33 for the reference and 0.30 for
      ! the response.
      a_h  = infos%tddft%cam_alpha; if (a_h  == CAM_UNSET) a_h  = a_ref
      b_h  = infos%tddft%cam_beta;  if (b_h  == CAM_UNSET) b_h  = b_ref
      mu_h = infos%tddft%cam_mu;    if (mu_h == CAM_UNSET) mu_h = mu_ref
      if (a_ref == CAM_UNSET .or. b_ref == CAM_UNSET .or. &
          mu_ref <= 0.0_dp .or. mu_h <= 0.0_dp) then
        write(iw,'(/,5x,a)') 'QMRSF-DK: the range-separation parameters of this '// &
             'reference are incomplete.'
        write(iw,'(5x,a,3f10.6)') 'reference alpha, beta, mu = ',a_ref,b_ref,mu_ref
        write(iw,'(5x,a,3f10.6)') 'kernel    alpha, beta, mu = ',a_h,b_h,mu_h
        call flush(iw)
        close(iw)
        return
      end if
    end if

    allocate(h_act(NACT,NACT), eri_act(NACT,NACT,NACT,NACT))
    allocate(eri_lr(NACT,NACT,NACT,NACT))
    if (is_cam) then
      call qmrsf_active_integrals(infos, NACT, act, ncore, h_act, eri_act, ecore, &
                                  kcore_act, mu_h, eri_lr, kcore_lr, mu_ref)
    else
      call qmrsf_active_integrals(infos, NACT, act, ncore, h_act, eri_act, ecore, kcore_act)
      eri_lr = 0.0_dp
      kcore_lr = 0.0_dp
    end if

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    allocate(cact(nbf,NACT))
    do i = 1, NACT
      cact(:,i) = mo_a(:,act(i))
    end do

    va_act = 0.0_dp
    vb_act = 0.0_dp
    if (is_dft) call native_active_vxc(infos, mo_a, cact, va_act, vb_act)
    ! After stripping the spin-averaged semilocal v_xc from F^(0), the paper's
    ! one-electron operator is h_act+(1-c_ref)K_core.  This direct expression is
    ! independent of OpenQP's stored effective-ROHF Fock representation.
    h_eff = h_act + (1.0_dp-a_ref)*kcore_act - b_ref*kcore_lr

    ! Paper seam: scale exchange-VALUE integrals, not merely the exchange role
    ! in each antisymmetrized contraction.  H_K is built after zeroing every
    ! Hartree density pair (ii|**) and (**|jj), exactly as _exchange_eri in the
    ! published reference implementation; H_DK = H - (1-c_H) H_K.
    !
    ! That mask is not invariant to the orientation of a degenerate active pair,
    ! so the two covariant alternatives of native_seam_mask are evaluated in the
    ! same run and reported side by side.  Only the seam differs between them;
    ! the reference orbitals, h_eff, the active integrals, the v_xc diagonal, the
    ! spin adaptation and every block structure are shared.
    call native_build_spinorb(h_eff,eri_act,1.0_dp,h1,g)
    call native_gen_dets(dets)
    call native_build_h(dets,h1,g,hbare)
    call native_build_ucsf(dets,ucsf,orth_err,csf_label)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, e_mo)
    do i = 1, NACT
      eps_act(i) = e_mo(act(i))
    end do
    call native_find_pairs(eps_act,npair,pairs)
    call native_pair_doublet(eri_act,pairs(1,1),pairs(2,1),pair_thc,pair_aniso)
    hzero = 0.0_dp

    do iseam = 0, NSEAM-1
      call native_seam_mask(eri_act,iseam,npair,pairs,eri_k)
      call native_build_spinorb(hzero,eri_k,1.0_dp,h1k,gk)
      call native_build_h(dets,h1k,gk,hkdet)
      ! Remove the exchange the functional already provides: K - K^HF, which is
      ! (1-c_h) K for a global hybrid and (1-alpha) K - beta K^lr for a
      ! range-separated one.  The seam mask is applied to both ranges, so the
      ! covariant conventions act on the full exchange operator.
      hdet = hbare - (1.0_dp-a_h)*hkdet
      if (is_cam) then
        call native_seam_mask(eri_lr,iseam,npair,pairs,eri_klr)
        call native_build_spinorb(hzero,eri_klr,1.0_dp,h1klr,gklr)
        call native_build_h(dets,h1klr,gklr,hklrdet)
        hdet = hdet + b_h*hklrdet
      end if
      ! The open-shell diagonals carry v_xc inside the full Kohn-Sham orbital
      ! energies and the two reference corrections; only the closed-shell
      ! (dressed-kernel) rows receive the separate occupation-based shift
      ! Delta_d^{v_xc}, applied after the spin adaptation.
      call native_relative_diagonals(dets,h_eff,eri_act,a_h,va_act,vb_act,eps_act,diag36v)
      do i = 1, NDET
        hdet(i,i) = diag36v(i)
      end do
      call dgemm('N','T',NDET,NDET,NDET,1.0_dp,hdet,NDET,ucsf,NDET, &
                 0.0_dp,work36,NDET)
      call dgemm('N','N',NDET,NDET,NDET,1.0_dp,ucsf,NDET,work36,NDET, &
                 0.0_dp,hcsf,NDET)

      if (iseam == SEAM_NATIVE) then
        cross_err = 0.0_dp
        do i = 1, NDET
          do j = 1, NDET
            if ((i <= NSING .and. j > NSING) .or. &
                (i > NSING .and. i <= NSING+NTRIP .and. &
                 (j <= NSING .or. j > NSING+NTRIP)) .or. &
                (i > NSING+NTRIP .and. j <= NSING+NTRIP)) then
              cross_err = max(cross_err,abs(hcsf(i,j)))
            end if
          end do
        end do
      end if

      asing = hcsf(1:NSING,1:NSING)
      atrip = hcsf(NSING+1:NSING+NTRIP,NSING+1:NSING+NTRIP)
      aquint = hcsf(NDET:NDET,NDET:NDET)
      if (iseam == SEAM_NATIVE) then
        ! main-text convention: open-shell rows already carry their v_xc via
        ! the full-KS eps and corrections; add Delta_d^{v_xc} to 0OS rows only.
        call native_vxc_csf(dets,ucsf,va_act,vb_act,dvcsf)
        do i = 1, NSING
          if (index(csf_label(i),'s') == 0) asing(i,i) = asing(i,i) + dvcsf(i)
        end do
      end if
      call diag_symm_full(0,NSING,asing,NSING,esing,ierr)
      if (ierr /= 0) then
        write(iw,'(5x,a,i0)') 'QMRSF-DK: diagonalization of the singlet manifold failed, info = ',ierr
        close(iw)
        return
      end if
      call diag_symm_full(0,NTRIP,atrip,NTRIP,etrip,ierr)
      if (ierr /= 0) then
        write(iw,'(5x,a,i0)') 'QMRSF-DK: diagonalization of the triplet manifold failed, info = ',ierr
        close(iw)
        return
      end if
      esing_all(:,iseam) = esing
      etrip_all(:,iseam) = etrip
      equint_all(iseam)  = aquint(1,1)
      asing_all(:,:,iseam) = asing
      atrip_all(:,:,iseam) = atrip
      if (iseam == SEAM_NATIVE) then
        asing_keep = asing
        atrip_keep = atrip
      end if
    end do

    esing = esing_all(:,SEAM_NATIVE)
    etrip = etrip_all(:,SEAM_NATIVE)
    asing = asing_keep
    atrip = atrip_keep
    equint(1) = equint_all(SEAM_NATIVE)
    qvec(1,1) = 1.0_dp

    vxc_max = max(maxval(abs(va_act)),maxval(abs(vb_act)))
    write(iw,'(/,5x,a)') 'QMRSF-DK response space'
    write(iw,'(5x,a)')   '-----------------------'
    write(iw,'(5x,a)')   'The four singly occupied orbitals of the quintet reference span a'
    write(iw,'(5x,a)')   'CAS(4,4).  Its 36 M_s = 0 determinants are spin adapted into 20'
    write(iw,'(5x,a)')   'singlet, 15 triplet and 1 quintet configuration state functions,'
    write(iw,'(5x,a)')   'which are reported below as three separate manifolds.'
    write(iw,'(/,5x,a,i0)')     'Basis functions                        = ',nbf
    write(iw,'(5x,a,i0)')       'Doubly occupied (inactive) orbitals    = ',ncore
    write(iw,'(5x,a,4(i0,1x))') 'Active (singly occupied) orbitals      = ',act
    write(iw,'(5x,a,f18.10)')   'Quintet ROKS reference energy E_ROKS   = ', &
                                infos%mol_energy%energy
    write(iw,'(5x,a,f10.6)')    'Exact exchange in the reference        = ',a_ref
    write(iw,'(5x,a,f10.6)')    'Exact exchange in the dressed kernel   = ',a_h
    if (is_cam) then
      write(iw,'(5x,a)')        'Range-separated dressed kernel: the exact-exchange operator'
      write(iw,'(5x,a)')        'is alpha*K + beta*K(erf(mu r)/r), applied to both the'
      write(iw,'(5x,a)')        'reference core exchange and the seam.'
      write(iw,'(5x,a,3f10.6)') 'reference alpha, beta, mu              = ',a_ref,b_ref,mu_ref
      write(iw,'(5x,a,3f10.6)') 'kernel    alpha, beta, mu              = ',a_h,b_h,mu_h
    end if
    write(iw,'(/,5x,a)')        'Numerical checks'
    write(iw,'(5x,a,es12.4)')   'Spin-adaptation orthonormality error   = ',orth_err
    write(iw,'(5x,a,es12.4)')   'Largest neglected inter-multiplicity coupling = ',cross_err
    write(iw,'(5x,a,es12.4)')   'Largest active-space v_xc element      = ',vxc_max
    write(iw,'(/,5x,a)')        'Spin-resolved v_xc on the active orbitals (Hartree)'
    write(iw,'(7x,a)')          'orbital        alpha             beta'
    do i = 1, NACT
      write(iw,'(9x,i0,4x,2f16.10)') act(i),va_act(i,i),vb_act(i,i)
    end do

    eref = infos%mol_energy%energy
    thresh = infos%control%conf_print_threshold

    call native_print_block(infos, 'singlet', 'S', 0.0_dp, NSING, esing, asing, &
                            csf_label(1:NSING), eref, esing(1), thresh)
    call native_print_block(infos, 'triplet', 'T', 2.0_dp, NTRIP, etrip, atrip, &
                            csf_label(NSING+1:NSING+NTRIP), eref, esing(1), thresh)
    call native_print_block(infos, 'quintet', 'Q', 6.0_dp, 1, equint, qvec, &
                            csf_label(NDET:NDET), eref, esing(1), thresh)

    call native_write_dump(a_h,a_ref,equint(1),esing,etrip,orth_err,cross_err, &
                           esing_all,etrip_all,equint_all,pair_aniso)
    write(iw,'(/,5x,a)') 'Machine-readable QMRSF-DK results written to qmrsf_dk_full_live.dat'

    deallocate(h_act,eri_act,eri_lr,cact)
    call flush(iw)
    close(iw)
  end subroutine qmrsf_dk_paper_native


  subroutine native_active_vxc(infos,mo,cact,va,vb)
    use types, only: information
    use mod_dft, only: dft_initialize, dftclean, dftexcor
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: unpack_matrix
    use oqp_linalg
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: mo(:,:), cact(:,:)
    real(dp), intent(out) :: va(NACT,NACT), vb(NACT,NACT)
    type(dft_grid_t) :: molgrid
    real(dp), allocatable :: fap(:), fbp(:), fsq(:,:), tmp(:,:), ca(:,:), cb(:,:)
    real(dp) :: eexc, tele, tkin
    integer :: nbf, ntri, isc

    nbf = size(mo,1)
    ntri = nbf*(nbf+1)/2
    isc = max(2,int(infos%control%scftype))
    allocate(fap(ntri),fbp(ntri),fsq(nbf,nbf),tmp(nbf,NACT),ca(nbf,nbf),cb(nbf,nbf))
    ca = mo
    cb = mo
    call dft_initialize(infos,infos%basis,molgrid)
    call dftexcor(infos%basis,molgrid,isc,fap,fbp,ca,cb,nbf,ntri, &
                  eexc,tele,tkin,infos)
    call unpack_matrix(fap,fsq,nbf,'U')
    call dgemm('N','N',nbf,NACT,nbf,1.0_dp,fsq,nbf,cact,nbf, &
               0.0_dp,tmp,nbf)
    call dgemm('T','N',NACT,NACT,nbf,1.0_dp,cact,nbf,tmp,nbf, &
               0.0_dp,va,NACT)
    call unpack_matrix(fbp,fsq,nbf,'U')
    call dgemm('N','N',nbf,NACT,nbf,1.0_dp,fsq,nbf,cact,nbf, &
               0.0_dp,tmp,nbf)
    call dgemm('T','N',NACT,NACT,nbf,1.0_dp,cact,nbf,tmp,nbf, &
               0.0_dp,vb,NACT)
    call dftclean(infos)
    deallocate(fap,fbp,fsq,tmp,ca,cb)
  end subroutine native_active_vxc


  subroutine native_build_spinorb(h,eri,c_h,h1,g)
    real(dp), intent(in) :: h(NACT,NACT), eri(NACT,NACT,NACT,NACT), c_h
    real(dp), intent(out) :: h1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: p,q,r,s,spat(NSO),spin(NSO),i
    real(dp) :: direct,exchange
    do i=1,NSO
      if (i <= NACT) then
        spat(i)=i; spin(i)=0
      else
        spat(i)=i-NACT; spin(i)=1
      end if
    end do
    h1=0.0_dp
    do p=1,NSO; do q=1,NSO
      if (spin(p)==spin(q)) h1(p,q)=h(spat(p),spat(q))
    end do; end do
    g=0.0_dp
    do p=1,NSO; do q=1,NSO; do r=1,NSO; do s=1,NSO
      direct=0.0_dp; exchange=0.0_dp
      if (spin(p)==spin(r) .and. spin(q)==spin(s)) &
        direct=eri(spat(p),spat(r),spat(q),spat(s))
      if (spin(p)==spin(s) .and. spin(q)==spin(r)) &
        exchange=eri(spat(p),spat(s),spat(q),spat(r))
      g(p,q,r,s)=direct-c_h*exchange
    end do; end do; end do; end do
  end subroutine native_build_spinorb


  subroutine native_exchange_eri(eri,eri_k)
    real(dp), intent(in) :: eri(NACT,NACT,NACT,NACT)
    real(dp), intent(out) :: eri_k(NACT,NACT,NACT,NACT)
    integer :: p,q,r,s
    eri_k = eri
    do p=1,NACT; do q=1,NACT; do r=1,NACT; do s=1,NACT
      if (p==q .or. r==s) eri_k(p,q,r,s)=0.0_dp
    end do; end do; end do; end do
  end subroutine native_exchange_eri


!> @brief Rotate the active two-electron tensor inside one orbital plane.
!> @detail Returns g'_{pqrs} = R_{pa} R_{qb} R_{rc} R_{sd} g_{abcd} for the
!>         SO(2) rotation R of the (ia,ib) plane by @p theta, i.e. the tensor
!>         as seen from the orbital frame a' = cos(theta) a + sin(theta) b,
!>         b' = -sin(theta) a + cos(theta) b.
  subroutine native_rotate_pair(g,ia,ib,theta,gout)
    real(dp), intent(in) :: g(NACT,NACT,NACT,NACT), theta
    integer, intent(in) :: ia,ib
    real(dp), intent(out) :: gout(NACT,NACT,NACT,NACT)
    real(dp) :: r(NACT,NACT), t1(NACT,NACT,NACT,NACT), t2(NACT,NACT,NACT,NACT)
    real(dp) :: c,s
    integer :: p,q,rr,ss,a
    c = cos(theta); s = sin(theta)
    r = 0.0_dp
    do p=1,NACT; r(p,p)=1.0_dp; end do
    r(ia,ia)= c; r(ia,ib)= s
    r(ib,ia)=-s; r(ib,ib)= c
    t1 = 0.0_dp
    do p=1,NACT; do a=1,NACT
      if (r(p,a)==0.0_dp) cycle
      t1(p,:,:,:) = t1(p,:,:,:) + r(p,a)*g(a,:,:,:)
    end do; end do
    t2 = 0.0_dp
    do q=1,NACT; do a=1,NACT
      if (r(q,a)==0.0_dp) cycle
      t2(:,q,:,:) = t2(:,q,:,:) + r(q,a)*t1(:,a,:,:)
    end do; end do
    t1 = 0.0_dp
    do rr=1,NACT; do a=1,NACT
      if (r(rr,a)==0.0_dp) cycle
      t1(:,:,rr,:) = t1(:,:,rr,:) + r(rr,a)*t2(:,:,a,:)
    end do; end do
    gout = 0.0_dp
    do ss=1,NACT; do a=1,NACT
      if (r(ss,a)==0.0_dp) cycle
      gout(:,:,:,ss) = gout(:,:,:,ss) + r(ss,a)*t1(:,:,:,a)
    end do; end do
  end subroutine native_rotate_pair


!> @brief Gauge-invariant anisotropy of the degenerate pair (ia,ib).
!> @detail The pair charge distributions split into the SO(2)-invariant
!>         rho_s = (a^2+b^2)/sqrt2 and the Lambda doublet
!>         (d,x) = ((a^2-b^2)/sqrt2, sqrt2 ab), which rotates by 2*theta.  The
!>         doublet interaction matrix V = [[V_dd,V_dx],[V_dx,V_xx]] with
!>         V_dd = ((aa|aa)+(bb|bb))/2 - (aa|bb), V_xx = 2(ab|ba),
!>         V_dx = (aa|ab) - (bb|ab) obeys V(theta) = R(2 theta) V R(2 theta)^T,
!>         so its eigenvalue splitting is the frame-independent measure of the
!>         anisotropy that makes the value-based seam non-covariant, and
!>         @p thc is the frame angle that diagonalizes it.
  subroutine native_pair_doublet(eri,ia,ib,thc,aniso)
    real(dp), intent(in) :: eri(NACT,NACT,NACT,NACT)
    integer, intent(in) :: ia,ib
    real(dp), intent(out) :: thc,aniso
    real(dp) :: vdd,vxx,vdx
    vdd = 0.5_dp*(eri(ia,ia,ia,ia)+eri(ib,ib,ib,ib)) - eri(ia,ia,ib,ib)
    vxx = 2.0_dp*eri(ia,ib,ib,ia)
    vdx = eri(ia,ia,ia,ib) - eri(ib,ib,ia,ib)
    thc = 0.25_dp*atan2(2.0_dp*vdx, vdd-vxx)
    aniso = sqrt((vdd-vxx)**2 + 4.0_dp*vdx**2)
  end subroutine native_pair_doublet


!> @brief Build the seam exchange tensor under the requested convention.
!> @detail The value-based mask of @ref native_exchange_eri is not covariant
!>         under a rotation of a degenerate active pair: it sends the doublet
!>         component d to the Hartree side and x to the exchange side, so every
!>         dressed element containing K_ab or J_ab moves by
!>         (1/2) sin^2(2 theta) Lambda.  Both alternatives below remove that
!>         freedom by treating the doublet isotropically; each reduces to the
!>         value-based mask when the pair carries no anisotropy, and both are
!>         inert at c_H = 1, where the seam is multiplied by (1 - c_H).
!>           mode 0  SEAM_NATIVE  the published value-based mask
!>           mode 1  SEAM_S3R     mean of the mask evaluated in the two
!>                                symmetry frames of the doublet (the frame
!>                                that diagonalizes V and its 45-degree
!>                                partner); doublet channel weight 1/2
!>           mode 2  SEAM_HAAR    mean of the mask over the full SO(2) orbit
!>                                (Reynolds projection); the traceless part of
!>                                the doublet-doublet channel carries 1/4
!> @brief Split a charge-distribution matrix into the four seam channels.
!> @detail For the pair plane (ia,ib) the symmetric distribution matrices
!>         decompose orthogonally into
!>           1 H   the SO(2)-invariant Hartree span: the diagonal outside the
!>                 pair plus the pair trace rho_s = (a^2+b^2)/sqrt2,
!>           2 Dd  the doublet component d = (a^2-b^2)/sqrt2,
!>           3 Dx  the doublet component x = sqrt2 ab,
!>           4 R   every remaining off-diagonal distribution.
!>         (Dd,Dx) is the Lambda doublet: it rotates by twice the orbital angle,
!>         while H and R are unmixed by the pair rotation.
  pure subroutine native_pair_channel(d,npair,pairs,ich,p)
    real(dp), intent(in) :: d(NACT,NACT)
    integer, intent(in) :: npair,pairs(2,2),ich
    real(dp), intent(out) :: p(NACT,NACT)
    real(dp) :: h(NACT,NACT),db(NACT,NACT),tr,df
    integer :: t,k,ia,ib
    logical :: inpair(NACT)
    h = 0.0_dp; db = 0.0_dp; inpair = .false.
    do k = 1, npair
      ia = pairs(1,k); ib = pairs(2,k)
      inpair(ia) = .true.; inpair(ib) = .true.
      tr = 0.5_dp*(d(ia,ia)+d(ib,ib))
      df = 0.5_dp*(d(ia,ia)-d(ib,ib))
      h(ia,ia) = tr;  h(ib,ib) = tr
      db(ia,ia) = df; db(ib,ib) = -df
      db(ia,ib) = d(ia,ib); db(ib,ia) = d(ib,ia)
    end do
    do t = 1, NACT
      if (.not.inpair(t)) h(t,t) = d(t,t)
    end do
    select case (ich)
    case (1); p = h                 ! SO(2)-invariant Hartree span
    case (2); p = db                ! every Lambda doublet
    case default; p = d - h - db    ! remaining off-diagonal distributions
    end select
  end subroutine native_pair_channel


!> @brief Identify the degenerate-pair candidates among the four active orbitals.
!> @detail Every adjacent pair that is degenerate to within PAIR_GAP_THRESHOLD is
!>         treated, up to the two disjoint pairs that a doubly degenerate
!>         frontier set such as the e1g/e2u pi orbitals of benzene or the
!>         distance-degenerate combinations of a separated dimer requires.  A
!>         reference without such a degeneracy has no orientation freedom to
!>         remove; the frontier pair is then used, which leaves the construction
!>         continuous as a degeneracy is approached and keeps the same subspace
!>         along a potential energy surface.
  pure subroutine native_find_pairs(eps,npair,pairs)
    real(dp), intent(in) :: eps(NACT)
    integer, intent(out) :: npair, pairs(2,2)
    real(dp) :: gap(NACT-1)
    integer :: k, kbest, ksecond
    do k = 1, NACT-1
      gap(k) = abs(eps(k+1)-eps(k))
    end do
    npair = 0
    pairs = 0
    kbest = minloc(gap,1)
    if (gap(kbest) < PAIR_GAP_THRESHOLD) then
      npair = 1
      pairs(:,1) = (/ kbest, kbest+1 /)
      ksecond = 0
      do k = 1, NACT-1
        if (abs(k-kbest) < 2) cycle
        if (gap(k) >= PAIR_GAP_THRESHOLD) cycle
        if (ksecond == 0) then
          ksecond = k
        else if (gap(k) < gap(ksecond)) then
          ksecond = k
        end if
      end do
      if (ksecond /= 0) then
        npair = 2
        pairs(:,2) = (/ ksecond, ksecond+1 /)
      end if
    else
      npair = 1
      pairs(:,1) = (/ NACT/2, NACT/2+1 /)
    end if
  end subroutine native_find_pairs


!> @brief Channel decomposition of the active two-electron tensor.
!> @detail ch(:,:,:,:,i,j) is the part of the tensor whose first
!>         charge-distribution index pair lies in channel i and whose second
!>         lies in channel j; the sixteen parts sum back to the tensor.
  subroutine native_seam_channels(g,npair,pairs,ch)
    real(dp), intent(in) :: g(NACT,NACT,NACT,NACT)
    integer, intent(in) :: npair,pairs(2,2)
    real(dp), intent(out) :: ch(NACT,NACT,NACT,NACT,3,3)
    real(dp) :: p1(NACT,NACT,NACT,NACT,3), m(NACT,NACT), p(NACT,NACT)
    integer :: i,j,q,r,s,t
    do i = 1, 3
      do r = 1, NACT; do s = 1, NACT
        m = g(:,:,r,s)
        call native_pair_channel(m,npair,pairs,i,p)
        p1(:,:,r,s,i) = p
      end do; end do
    end do
    do i = 1, 3; do j = 1, 3
      do t = 1, NACT; do q = 1, NACT
        m = p1(t,q,:,:,i)
        call native_pair_channel(m,npair,pairs,j,p)
        ch(t,q,:,:,i,j) = p
      end do; end do
    end do; end do
  end subroutine native_seam_channels


!> @brief Build the seam exchange tensor under the requested convention.
!> @detail The value-based mask of @ref native_exchange_eri keeps exactly the
!>         off-diagonal charge distributions, i.e. the channels (Dx,R) on both
!>         index pairs, and discards the doublet component d together with the
!>         Hartree span.  Splitting one doublet that way is what breaks the
!>         invariance: under a rotation of the pair d and x mix, so every
!>         dressed element containing K_ab or J_ab moves by
!>         (1/2) sin^2(2 theta) Lambda, while the antisymmetrized J - K does
!>         not move at all.  The alternatives below weight the doublet
!>         isotropically and are therefore orientation independent by
!>         construction, with no reference frame entering their definition:
!>           mode 0  SEAM_NATIVE  the published mask, (Dx,R) x (Dx,R)
!>           mode 1  SEAM_S3R     every channel touching the doublet at 1/2
!>           mode 2  SEAM_HAAR    the same, except that the doublet-doublet
!>                                block keeps 1/2 of its trace and 1/4 of its
!>                                traceless part -- the Reynolds projection,
!>                                i.e. the mean of the published mask over the
!>                                whole SO(2) orbit
!>         Both reduce to the published mask when the pair carries no doublet
!>         content, and both are inert at c_H = 1 because the seam enters
!>         multiplied by (1 - c_H).
  subroutine native_seam_mask(eri,mode,npair,pairs,eri_k)
    real(dp), intent(in) :: eri(NACT,NACT,NACT,NACT)
    integer, intent(in) :: mode,npair,pairs(2,2)
    real(dp), intent(out) :: eri_k(NACT,NACT,NACT,NACT)
    integer, parameter :: ICH_H = 1, ICH_D = 2, ICH_R = 3
    integer, parameter :: NQUAD = 64
    real(dp) :: ch(NACT,NACT,NACT,NACT,3,3)
    real(dp) :: g1(NACT,NACT,NACT,NACT), g2(NACT,NACT,NACT,NACT)
    real(dp) :: m1(NACT,NACT,NACT,NACT), b1(NACT,NACT,NACT,NACT)
    real(dp) :: acc(NACT,NACT,NACT,NACT)
    real(dp) :: tha,thb,pi,wgt
    integer :: ka,kb,nb

    select case (mode)

    case (SEAM_NATIVE)
      call native_exchange_eri(eri,eri_k)

    case (SEAM_S3R)
      ! Keep the off-diagonal distributions as the published mask does, but give
      ! every channel that touches a Lambda doublet the uniform weight 1/2, so
      ! that d and x are treated alike and no orientation survives.  For a single
      ! pair this is the mean of the published mask taken in the two symmetry
      ! frames of the doublet.
      call native_seam_channels(eri,npair,pairs,ch)
      eri_k = ch(:,:,:,:,ICH_R,ICH_R) &
            + 0.5_dp*(ch(:,:,:,:,ICH_D,ICH_R) + ch(:,:,:,:,ICH_R,ICH_D) &
                    + ch(:,:,:,:,ICH_D,ICH_D))

    case (SEAM_HAAR)
      ! Reynolds projection: the mean of the published mask over the whole orbit
      ! of the group that leaves the reference invariant, which is one SO(2) per
      ! degenerate pair.  Evaluated by quadrature over the product group, so no
      ! channel algebra and no reference frame enter.
      pi = acos(-1.0_dp)
      nb = 1
      if (npair >= 2) nb = NQUAD
      wgt = 1.0_dp/real(NQUAD*nb,dp)
      acc = 0.0_dp
      do ka = 0, NQUAD-1
        tha = real(ka,dp)*(0.5_dp*pi)/real(NQUAD,dp)
        call native_rotate_pair(eri,pairs(1,1),pairs(2,1),tha,g1)
        do kb = 0, nb-1
          if (npair >= 2) then
            thb = real(kb,dp)*(0.5_dp*pi)/real(nb,dp)
            call native_rotate_pair(g1,pairs(1,2),pairs(2,2),thb,g2)
          else
            thb = 0.0_dp
            g2 = g1
          end if
          call native_exchange_eri(g2,m1)
          if (npair >= 2) then
            call native_rotate_pair(m1,pairs(1,2),pairs(2,2),-thb,b1)
          else
            b1 = m1
          end if
          call native_rotate_pair(b1,pairs(1,1),pairs(2,1),-tha,g2)
          acc = acc + g2
        end do
      end do
      eri_k = wgt*acc

    case default
      call native_exchange_eri(eri,eri_k)
    end select
  end subroutine native_seam_mask


  subroutine native_gen_dets(dets)
    integer, intent(out) :: dets(4,NDET)
    integer :: a1,a2,b1,b2,k,t(4),i,j,z
    k=0
    do a1=1,NACT-1; do a2=a1+1,NACT
      do b1=1,NACT-1; do b2=b1+1,NACT
        k=k+1; t=(/a1,a2,b1+NACT,b2+NACT/)
        do i=1,3; do j=i+1,4
          if (t(j)<t(i)) then; z=t(i); t(i)=t(j); t(j)=z; end if
        end do; end do
        dets(:,k)=t
      end do; end do
    end do; end do
  end subroutine native_gen_dets


  pure logical function native_inset(x,d)
    integer, intent(in) :: x,d(4)
    native_inset=any(d==x)
  end function native_inset


  real(dp) function native_melem(d1,d2,h1,g)
    integer, intent(in) :: d1(4),d2(4)
    real(dp), intent(in) :: h1(NSO,NSO),g(NSO,NSO,NSO,NSO)
    integer :: holes(4),parts(4),common(4),nh,np,nc,occ(4),nocc
    integer :: i,k,idx,pp,hh,qc
    real(dp) :: sgn,val
    nh=0; np=0; nc=0
    do i=1,4
      if (.not.native_inset(d2(i),d1)) then; nh=nh+1; holes(nh)=d2(i); end if
      if (.not.native_inset(d1(i),d2)) then; np=np+1; parts(np)=d1(i); end if
      if (native_inset(d1(i),d2)) then; nc=nc+1; common(nc)=d1(i); end if
    end do
    if (nh>2) then; native_melem=0.0_dp; return; end if
    occ=d2; nocc=4; sgn=1.0_dp
    do k=1,nh
      idx=0
      do i=1,nocc
        if (occ(i)==holes(k)) then; idx=i; exit; end if
      end do
      if (mod(idx-1,2)==1) sgn=-sgn
      do i=idx,nocc-1; occ(i)=occ(i+1); end do
      nocc=nocc-1
    end do
    do k=np,1,-1
      idx=1
      do i=1,nocc
        if (occ(i)<parts(k)) idx=idx+1
      end do
      if (mod(idx-1,2)==1) sgn=-sgn
      do i=nocc,idx,-1; occ(i+1)=occ(i); end do
      occ(idx)=parts(k); nocc=nocc+1
    end do
    if (nh==0) then
      val=0.0_dp
      do i=1,4; val=val+h1(d1(i),d1(i)); end do
      do i=1,3; do k=i+1,4
        val=val+g(d1(i),d1(k),d1(i),d1(k))
      end do; end do
      native_melem=val
    else if (nh==1) then
      pp=parts(1); hh=holes(1); val=h1(pp,hh)
      do i=1,nc; qc=common(i); val=val+g(pp,qc,hh,qc); end do
      native_melem=sgn*val
    else
      native_melem=sgn*g(parts(1),parts(2),holes(1),holes(2))
    end if
  end function native_melem


  subroutine native_build_h(dets,h1,g,hmat)
    integer, intent(in) :: dets(4,NDET)
    real(dp), intent(in) :: h1(NSO,NSO),g(NSO,NSO,NSO,NSO)
    real(dp), intent(out) :: hmat(NDET,NDET)
    integer :: i,j
    do i=1,NDET
      do j=i,NDET
        hmat(i,j)=native_melem(dets(:,i),dets(:,j),h1,g)
        hmat(j,i)=hmat(i,j)
      end do
    end do
  end subroutine native_build_h


  subroutine native_add_vxc_diag(dets,va,vb,hmat)
    integer, intent(in) :: dets(4,NDET)
    real(dp), intent(in) :: va(NACT,NACT),vb(NACT,NACT)
    real(dp), intent(inout) :: hmat(NDET,NDET)
    integer :: d,p,na,nb
    do d=1,NDET
      do p=1,NACT
        na=merge(1,0,any(dets(:,d)==p))
        nb=merge(1,0,any(dets(:,d)==p+NACT))
        hmat(d,d)=hmat(d,d)+real(na-1,dp)*va(p,p)+real(nb,dp)*vb(p,p)
      end do
    end do
  end subroutine native_add_vxc_diag





  !> Delta^orb for a 2OS flip (q -> p): the exact, path-averaged difference
  !> between the generating-reference eps and the canonical ROKS eps,
  !>   Delta^orb_(pq) = (eps_p - eps_q)[Xi^{+-1}] - (eps_p - eps_q)[ROKS]
  !>                  = (c_h/2)(K_pp + K_qq) + c_h K_pq .
  !> All h and J cancel (both determinants occupy the same four orbitals);
  !> only exchange survives.  Identity: 1/2(two generating-path eps-diffs)
  !> == (eps_p - eps_q)[ROKS] + Delta^orb_(pq), exact at any h_eff.
  pure real(dp) function native_dorb_2os(eri,c_h,va,vb,p,q)
    real(dp), intent(in) :: eri(NACT,NACT,NACT,NACT), c_h
    real(dp), intent(in) :: va(NACT,NACT), vb(NACT,NACT)
    integer, intent(in) :: p, q
    ! exchange part + the v_xc part of the same eps difference:
    ! the generating-reference eps carry v^beta_p / v^alpha_q, the canonical
    ! ones 1/2(v^a+v^b); the difference is 1/2(v^b-v^a) on each flip orbital.
    native_dorb_2os = 0.5_dp*c_h*(eri(p,p,p,p) + eri(q,q,q,q)) &
                    + c_h*eri(p,q,q,p) &
                    + 0.5_dp*((vb(p,p)-va(p,p)) + (vb(q,q)-va(q,q)))
  end function native_dorb_2os

  !> Delta^orb for a 4OS self-flip (j -> j) in the reference whose same-pair
  !> partner of j is jp and opposite pair is (m1,m2):
  !>   Delta^orb[j] = (eps_j^to - eps_j^from)[Xi] - (eps_j - eps_j)[ROKS]
  !>               = c_h [ K_j,m1 + K_j,m2 + (jj|jj) - K_j,jp ] ,
  !> the canonical part being identically zero for a self-flip.
  pure real(dp) function native_dorb_selfflip(eri,c_h,va,vb,j,jp,m1,m2)
    real(dp), intent(in) :: eri(NACT,NACT,NACT,NACT), c_h
    real(dp), intent(in) :: va(NACT,NACT), vb(NACT,NACT)
    integer, intent(in) :: j, jp, m1, m2
    ! exchange part + the v_xc spin splitting of level j
    ! (canonical part of a self-flip is identically zero).
    native_dorb_selfflip = c_h*( eri(j,m1,m1,j) + eri(j,m2,m2,j) &
                               + eri(j,j,j,j) - eri(j,jp,jp,j) ) &
                         + (vb(j,j)-va(j,j))
  end function native_dorb_selfflip

  !> Reference-energy alignment Delta^SR (main-text sign), for the
  !> reference whose minority spin sits on orbital m: same-spin exchange part
  !> + the v_xc content of the reference energy difference (the minority
  !> orbital's spin splitting, occupation-based).
  pure real(dp) function native_delta_sr(eri,c_h,va,vb,m)
    real(dp), intent(in) :: eri(NACT,NACT,NACT,NACT), c_h
    real(dp), intent(in) :: va(NACT,NACT), vb(NACT,NACT)
    integer, intent(in) :: m
    integer :: j
    native_delta_sr = vb(m,m) - va(m,m)
    do j = 1, NACT
      if (j /= m) native_delta_sr = native_delta_sr + c_h*eri(m,j,j,m)
    end do
  end function native_delta_sr

  !> E_uuuu = <Xi^{(+2)}|H_c|Xi^{(+2)}>: all four actives alpha.
  pure real(dp) function native_e_uuuu(h_eff,eri,c_h)
    real(dp), intent(in) :: h_eff(NACT,NACT), eri(NACT,NACT,NACT,NACT), c_h
    integer :: i, j
    native_e_uuuu = 0.0_dp
    do i = 1, NACT
      native_e_uuuu = native_e_uuuu + h_eff(i,i)
      do j = i+1, NACT
        native_e_uuuu = native_e_uuuu + eri(i,i,j,j) - c_h*eri(i,j,j,i)
      end do
    end do
  end function native_e_uuuu

  !> 4OS common origin, relative to the stationary quintet: four self-flip paths
  !> averaged + the 1/4-averaged Delta^SR, evaluated on the
  !> aabb representative with its own pair exchanges added back.
  real(dp) function native_fouros_common(h_eff,eri,c_h,va,vb)
    real(dp), intent(in) :: h_eff(NACT,NACT), eri(NACT,NACT,NACT,NACT), c_h
    real(dp), intent(in) :: va(NACT,NACT), vb(NACT,NACT)
    integer :: b, a, bo, ao, m
    integer :: aocc3(3), bocc1(1), aocc1(1), bocc3(3)
    real(dp) :: acc, dsr
    ! canonical spelling: for a self-flip the canonical ROKS eps difference is
    ! identically zero (eps_j - eps_j), so each path contributes
    !   0 + Delta^orb[path] - c_h(jj|jj)
    ! with Delta^orb[path] the exact generating-vs-canonical conversion; this
    ! is the same number, path by path, as the former eps_ref differences.
    acc = 0.0_dp
    do b = 3, 4                       ! M_s=+1 refs: flip the beta-pair member
      bo = 7 - b
      acc = acc + native_dorb_selfflip(eri,c_h,va,vb,b,bo,1,2) - c_h*eri(b,b,b,b)
    end do
    do a = 1, 2                       ! M_s=-1 refs (spin mirror)
      ao = 3 - a
      acc = acc + native_dorb_selfflip(eri,c_h,va,vb,a,ao,3,4) - c_h*eri(a,a,a,a)
    end do
    dsr = 0.0_dp
    do m = 1, NACT
      dsr = dsr + native_delta_sr(eri,c_h,va,vb,m)
    end do
    native_fouros_common = 0.25_dp*acc + 0.25_dp*dsr &
                         + c_h*(eri(1,2,2,1) + eri(3,4,4,3))
  end function native_fouros_common

  !> Determinant diagonals relative to the stationary quintet, typed as the paper:
  !>   2OS: orbital Hessian, two generating paths, 1/2-averaged Delta^SR.
  !>        kernel = -(qq|pp) + (1-c_h)(qq|qq), equivalently -c_HF(ii|aa) plus
  !>        the determinant-alignment term Delta^DA = (1-c_HF)[(ii|ii)-(ii|aa)].
  !>        (qq|pp) is a Hartree integral and is not scaled by c_HF; the
  !>        (1-c_h)(qq|qq) piece removes the hole orbital's residual hybrid
  !>        self-interaction.  Both vanish at c_h = 1.
  !>   4OS: D_common - c_h(K_pairA + K_pairB)
  !>   0OS: w_ij = 2h+2h+J+J+4J - 2 c_h K - E_uuuu
  subroutine native_relative_diagonals(dets,h_eff,eri,c_h,va,vb,eps_mo,diag)
    integer, intent(in) :: dets(4,NDET)
    real(dp), intent(in) :: h_eff(NACT,NACT), eri(NACT,NACT,NACT,NACT), c_h
    real(dp), intent(in) :: va(NACT,NACT), vb(NACT,NACT)
    !> Converged ROKS quintet MO energies (OQP_E_MO_A).
    !>
    !> eps_i and eps_a of Eq.(master-open) ARE the ROKS orbital energies.  They
    !> are read from the SCF and never rebuilt, and c_HF -- which scales the
    !> EXCHANGE -- must not be applied to them: the rest of the matrix is scaled
    !> to match the MOs, not the MOs to match the kernel.
    real(dp), intent(in) :: eps_mo(NACT)
    real(dp), intent(out) :: diag(NDET)
    character(len=4) :: cfg
    integer :: i, k, p, q, r, ssng, u1, u2, d1, d2, n2
    integer :: aoccA(3), boccA(1), aoccB(1), boccB(3)
    real(dp) :: D4, EU, a1, a2, kern, dsr
    D4 = native_fouros_common(h_eff,eri,c_h,va,vb)
    EU = native_e_uuuu(h_eff,eri,c_h)
    do i = 1, NDET
      call native_det_config(dets(:,i),cfg)
      n2 = 0
      do k = 1, NACT
        if (cfg(k:k) == '2') n2 = n2 + 1
      end do
      if (n2 == 2) then                              ! ---- 0OS (DSF)
        p = 0; q = 0
        do k = 1, NACT
          if (cfg(k:k) == '2') then
            if (p == 0) then; p = k; else; q = k; end if
          end if
        end do
        diag(i) = 2.0_dp*h_eff(p,p) + 2.0_dp*h_eff(q,q) &
                + eri(p,p,p,p) + eri(q,q,q,q) + 4.0_dp*eri(p,p,q,q) &
                - 2.0_dp*c_h*eri(p,q,q,p) - EU
      else if (n2 == 1) then                         ! ---- 2OS
        p = index(cfg,'2'); q = index(cfg,'0')
        r = 0; ssng = 0
        do k = 1, NACT
          if (cfg(k:k) == 'u' .or. cfg(k:k) == 'd') then
            if (r == 0) then; r = k; else; ssng = k; end if
          end if
        end do
        ! canonical ROKS quintet orbital energies + the exact orbital
        ! correction Delta^orb (identical to the two-generating-path average
        ! 1/2(a1+a2) of the Xi^{+-1} eps differences -- pure-exchange identity)
        a1 = eps_mo(p) - eps_mo(q)                         ! (eps_a - eps_i)_ROKS, read from the SCF
        a2 = native_dorb_2os(eri,c_h,va,vb,p,q)            ! Delta^orb_(pq) incl v_xc part
        kern = -eri(q,q,p,p) + (1.0_dp-c_h)*eri(q,q,q,q)
        dsr = 0.5_dp*(native_delta_sr(eri,c_h,va,vb,r) + native_delta_sr(eri,c_h,va,vb,ssng))
        diag(i) = a1 + a2 + kern + dsr
      else                                           ! ---- 4OS
        u1 = 0; u2 = 0; d1 = 0; d2 = 0
        do k = 1, NACT
          if (cfg(k:k) == 'u') then
            if (u1 == 0) then; u1 = k; else; u2 = k; end if
          else if (cfg(k:k) == 'd') then
            if (d1 == 0) then; d1 = k; else; d2 = k; end if
          end if
        end do
        diag(i) = D4 - c_h*(eri(u1,u2,u2,u1) + eri(d1,d2,d2,d1))
      end if
    end do
  end subroutine native_relative_diagonals

  !> Per-CSF v_xc occupation rule, applied after spin
  !> adaptation.  Every CSF in this construction spans determinants of one
  !> spatial configuration, so the rule is exact per CSF and generates no
  !> cross-multiplicity coupling.
  subroutine native_vxc_csf(dets,ucsf,va,vb,dv)
    integer, intent(in) :: dets(4,NDET)
    real(dp), intent(in) :: ucsf(NDET,NDET), va(NACT,NACT), vb(NACT,NACT)
    real(dp), intent(out) :: dv(NDET)
    character(len=4) :: cfg
    integer :: m, d, dbest, k
    real(dp) :: best
    do m = 1, NDET
      best = 0.0_dp; dbest = 1
      do d = 1, NDET
        if (abs(ucsf(m,d)) > best) then
          best = abs(ucsf(m,d)); dbest = d
        end if
      end do
      call native_det_config(dets(:,dbest),cfg)
      dv(m) = 0.0_dp
      do k = 1, NACT
        if (cfg(k:k) == '2') then
          dv(m) = dv(m) + vb(k,k)
        else if (cfg(k:k) == '0') then
          dv(m) = dv(m) - va(k,k)
        else
          dv(m) = dv(m) + 0.5_dp*(vb(k,k) - va(k,k))
        end if
      end do
    end do
  end subroutine native_vxc_csf

  subroutine native_det_config(det,cfg)
    integer, intent(in) :: det(4)
    character(len=4), intent(out) :: cfg
    integer :: p
    logical :: a,b
    cfg='0000'
    do p=1,NACT
      a=any(det==p); b=any(det==p+NACT)
      if (a.and.b) then; cfg(p:p)='2'
      else if (a) then; cfg(p:p)='u'
      else if (b) then; cfg(p:p)='d'
      end if
    end do
  end subroutine native_det_config


  integer function native_os_class(cfg)
    character(len=4), intent(in) :: cfg
    integer :: p
    native_os_class=0
    do p=1,NACT
      if (cfg(p:p)=='u' .or. cfg(p:p)=='d') native_os_class=native_os_class+1
    end do
  end function native_os_class


  integer function native_fermi_phase(cfg)
    character(len=4), intent(in) :: cfg
    integer :: seq(4),n,p,i,j
    n=0
    do p=1,NACT
      select case(cfg(p:p))
      case('2'); n=n+1; seq(n)=p; n=n+1; seq(n)=p+NACT
      case('u'); n=n+1; seq(n)=p
      case('d'); n=n+1; seq(n)=p+NACT
      end select
    end do
    native_fermi_phase=1
    do i=1,n-1; do j=i+1,n
      if (seq(i)>seq(j)) native_fermi_phase=-native_fermi_phase
    end do; end do
  end function native_fermi_phase


  subroutine native_build_ucsf(dets,u,orth_err,label)
    integer, intent(in) :: dets(4,NDET)
    real(dp), intent(out) :: u(NDET,NDET),orth_err
    character(len=8), intent(out), optional :: label(NDET)
    character(len=8) :: lbl(NDET)
    character(len=4) :: cfg(NDET),sw,fo_cfg(6)
    logical :: seen(NDET)
    integer :: phase(NDET),i,j,p,rs,rt,k,fo(6)
    real(dp) :: c4(6,6),gram(NDET,NDET),sq2,sq3,sq6

    sq2=sqrt(2.0_dp); sq3=sqrt(3.0_dp); sq6=sqrt(6.0_dp)
    u=0.0_dp; seen=.false.; rs=0; rt=NSING
    do i=1,NDET
      call native_det_config(dets(:,i),cfg(i))
      phase(i)=native_fermi_phase(cfg(i))
    end do

    lbl=' '
    do i=1,NDET
      if (native_os_class(cfg(i))==0) then
        rs=rs+1; u(rs,i)=real(phase(i),dp)
        lbl(rs)=cfg(i)
      end if
    end do

    do i=1,NDET
      if (native_os_class(cfg(i))/=2 .or. seen(i)) cycle
      sw=cfg(i)
      do p=1,NACT
        if (sw(p:p)=='u') then; sw(p:p)='d'
        else if (sw(p:p)=='d') then; sw(p:p)='u'
        end if
      end do
      j=0
      do k=1,NDET
        if (cfg(k)==sw) then; j=k; exit; end if
      end do
      if (j==0) error stop 'QMRSF-DK native: missing 2OS spin partner'
      seen(i)=.true.; seen(j)=.true.
      rs=rs+1
      u(rs,i)= real(phase(i),dp)/sq2
      u(rs,j)=-real(phase(j),dp)/sq2
      lbl(rs)=native_open_label(cfg(i))
      rt=rt+1
      u(rt,i)=real(phase(i),dp)/sq2
      u(rt,j)=real(phase(j),dp)/sq2
      lbl(rt)=native_open_label(cfg(i))
    end do

    fo_cfg=(/'uudd','udud','uddu','duud','dudu','dduu'/)
    do i=1,6
      fo(i)=0
      do j=1,NDET
        if (cfg(j)==fo_cfg(i)) then; fo(i)=j; exit; end if
      end do
      if (fo(i)==0) error stop 'QMRSF-DK native: missing 4OS determinant'
    end do
    c4=0.0_dp
    c4(:,1)=(/0.0_dp,0.5_dp,-0.5_dp,-0.5_dp,0.5_dp,0.0_dp/)
    c4(:,2)=(/1.0_dp/sq3,-0.5_dp/sq3,-0.5_dp/sq3,-0.5_dp/sq3,-0.5_dp/sq3,1.0_dp/sq3/)
    c4(:,3)=(/0.0_dp,0.5_dp,0.5_dp,-0.5_dp,-0.5_dp,0.0_dp/)
    c4(:,4)=(/1.0_dp/sq3,-0.5_dp/sq3,0.5_dp/sq3,-0.5_dp/sq3,0.5_dp/sq3,-1.0_dp/sq3/)
    c4(:,5)=(/1.0_dp/sq6,1.0_dp/sq6,-1.0_dp/sq6,1.0_dp/sq6,-1.0_dp/sq6,-1.0_dp/sq6/)
    c4(:,6)=(/1.0_dp/sq6,1.0_dp/sq6,1.0_dp/sq6,1.0_dp/sq6,1.0_dp/sq6,1.0_dp/sq6/)
    do k=1,2
      rs=rs+1
      do i=1,6; u(rs,fo(i))=c4(i,k)*real(phase(fo(i)),dp); end do
      write(lbl(rs),'(a,i1,a)') 'ssss[S',k,']'
    end do
    do k=3,5
      rt=rt+1
      do i=1,6; u(rt,fo(i))=c4(i,k)*real(phase(fo(i)),dp); end do
      write(lbl(rt),'(a,i1,a)') 'ssss[T',k-2,']'
    end do
    do i=1,6; u(NDET,fo(i))=c4(i,6)*real(phase(fo(i)),dp); end do
    lbl(NDET)='ssss[Q]'

    if (rs/=NSING .or. rt/=NSING+NTRIP) &
      error stop 'QMRSF-DK native: incorrect spin-adapted block dimensions'
    call dgemm('N','T',NDET,NDET,NDET,1.0_dp,u,NDET,u,NDET, &
               0.0_dp,gram,NDET)
    do i=1,NDET; gram(i,i)=gram(i,i)-1.0_dp; end do
    orth_err=maxval(abs(gram))
    if (present(label)) label=lbl
  end subroutine native_build_ucsf


!> @brief Print one spin manifold of the QMRSF-DK spectrum.
!> @detail A single quintet reference yields the singlet, triplet and quintet
!>         manifolds of the CAS(4,4) at once.  Each manifold is reported with
!>         the same layout as the TD-DFT summary: total energies, excitation
!>         energies measured from the QMRSF-DK ground state, and the leading
!>         spin-adapted configurations of every root.
!> @param[in]  infos      calculation metadata
!> @param[in]  manifold   name of the spin manifold
!> @param[in]  tag        one-letter state tag (S, T or Q)
!> @param[in]  ssq        exact <S^2> of the manifold
!> @param[in]  nstates    number of roots in the manifold
!> @param[in]  eig        response eigenvalues of the manifold
!> @param[in]  vec        eigenvectors, one column per root
!> @param[in]  label      configuration label of every basis CSF
!> @param[in]  eref       E_ROKS(quintet) - A_quintet, the energy-law offset
!> @param[in]  e0         lowest singlet eigenvalue, the zero of excitation
!> @param[in]  thresh     print threshold on the configuration coefficients
  subroutine native_print_block(infos, manifold, tag, ssq, nstates, eig, vec, &
                                label, eref, e0, thresh)
    use types, only: information
    use io_constants, only: iw
    use physical_constants, only: UNITS_EV

    type(information), intent(in) :: infos
    character(len=*), intent(in) :: manifold, tag
    real(dp), intent(in) :: ssq, eref, e0, thresh
    integer, intent(in) :: nstates
    real(dp), intent(in) :: eig(:), vec(:,:)
    character(len=8), intent(in) :: label(:)

    integer :: istate, icsf
    real(dp) :: etot, dele, coeff, weight

    write(iw,'(/,1x,78("^"))')
    write(iw,'(/,8x,3a,/)') 'Summary of the QMRSF-DK ', manifold, ' manifold'
    write(iw,'(A12,A20,A18,A16,A12)') 'State', 'E(Hartree)', 'omega(Hartree)', &
                                      'dE(eV)', '<S^2>'

    do istate = 1, nstates
      etot = eref + eig(istate)
      dele = (eig(istate) - e0)/UNITS_EV
      write(iw,'(3X,A,G0,t13,F20.10,F18.10,F16.6,F12.4)') tag, istate-1, etot, &
                                                          eig(istate), dele, ssq
    end do

    write(iw,'(/,8x,3a,/)') 'Leading configurations of the ', manifold, ' states'

    do istate = 1, nstates
      dele = (eig(istate) - e0)/UNITS_EV
      write(iw,'(/,1x,"State ",a,i0,2x,"Energy =",f12.6,1x,"eV")') trim(tag), istate-1, dele
      write(iw,'(15x,"<S^2> =",1x,f9.4)') ssq
      write(iw,'(8x,"CSF",4x,"Coeff",5x,"Weight(%)",4x,"Configuration")')
      write(iw,'(8x,3("-"),2x,9("-"),2x,9("-"),4x,13("-"))')
      do icsf = 1, nstates
        coeff = vec(icsf,istate)
        if (abs(coeff) <= thresh) cycle
        weight = 100.0_dp*coeff*coeff
        write(iw,'(7x,i4,1x,f9.6,2x,f9.2,6x,a)') icsf, coeff, weight, trim(label(icsf))
      end do
    end do

    write(iw,'(1x,78("=")/)')

  end subroutine native_print_block


!> @brief Occupation label of an open-shell configuration state function.
!> @detail Singly occupied orbitals are singlet- or triplet-coupled in pairs,
!>         so the individual alpha/beta labels carry no meaning; both are
!>         printed as 's', in the notation used for the CAS(4,4) manifold.
  pure function native_open_label(cfg) result(lbl)
    character(len=4), intent(in) :: cfg
    character(len=8) :: lbl
    integer :: p
    lbl=' '
    do p=1,NACT
      if (cfg(p:p)=='u' .or. cfg(p:p)=='d') then
        lbl(p:p)='s'
      else
        lbl(p:p)=cfg(p:p)
      end if
    end do
  end function native_open_label








  subroutine native_write_dump(c_h,c_ref,aq,es,et,orth_err,cross_err, &
                               esing_all,etrip_all,equint_all,aniso)
    real(dp), intent(in) :: c_h,c_ref,aq,es(NSING),et(NTRIP),orth_err,cross_err
    real(dp), intent(in) :: esing_all(NSING,0:NSEAM-1), etrip_all(NTRIP,0:NSEAM-1)
    real(dp), intent(in) :: equint_all(0:NSEAM-1), aniso
    integer :: u,m
    open(newunit=u,file='qmrsf_dk_full_live.dat',status='replace',action='write')
    write(u,'(a,4(1x,i0))') 'QMRSF_DK_NATIVE_V1',NACT,NSING,NTRIP,1
    write(u,'(2es24.16)') c_h,c_ref
    write(u,'(es24.16)') aq
    write(u,'(*(es24.16,1x))') es
    write(u,'(*(es24.16,1x))') et
    write(u,'(es24.16)') aq
    write(u,'(2es24.16)') orth_err,cross_err
    write(u,'(a,1x,i0)') 'QMRSF_DK_SEAMS_V1',NSEAM
    write(u,'(es24.16)') aniso
    do m = 0, NSEAM-1
      write(u,'(a)') trim(SEAM_NAME(m))
      write(u,'(es24.16)') equint_all(m)
      write(u,'(*(es24.16,1x))') esing_all(:,m)
      write(u,'(*(es24.16,1x))') etrip_all(:,m)
    end do
    close(u)
  end subroutine native_write_dump



end module qmrsf_dk_paper_native_mod
