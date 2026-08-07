!> @file qmrsf_dk_paper_native.F90
!> @brief Paper-faithful QMRSF-DK implementation using only OpenQP data.
!>
!> This module is the native counterpart of
!> spikes/qmrsf_dk_original/qmrsf_dk_paper_direct.py.  It obtains the converged
!> ROKS/ROHF orbitals, spin Fock matrices, and spin-resolved v_xc directly from
!> OpenQP, transforms the four-SOMO active integrals with qmrsf_ao2mo_mod, and
!> constructs the spin-pure CAS(4,4) response blocks with the analytic
!> 20-singlet/15-triplet/1-quintet transformation.  PySCF is not imported or
!> required at run time.
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
  integer, parameter :: NSEAM = 3
  integer, parameter :: PAIR_A = 2, PAIR_B = 3
  character(len=6), parameter :: SEAM_NAME(0:NSEAM-1) = &
       (/ 'native', 'S3R   ', 'Haar  ' /)

contains

  subroutine qmrsf_dk_paper_native(infos)
    use types, only: information
    use io_constants, only: iw
    use printing, only: print_module_info
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use qmrsf_ao2mo_mod, only: qmrsf_active_integrals
    use eigen, only: diag_symm_full
    use oqp_linalg

    type(information), target, intent(inout) :: infos

    real(dp), contiguous, pointer :: mo_a(:,:)
    real(dp), allocatable :: h_act(:,:), eri_act(:,:,:,:), cact(:,:)
    real(dp) :: va_act(NACT,NACT), vb_act(NACT,NACT)
    real(dp) :: h_eff(NACT,NACT), kcore_act(NACT,NACT)
    real(dp) :: hzero(NACT,NACT), eri_k(NACT,NACT,NACT,NACT)
    real(dp) :: h1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp) :: h1k(NSO,NSO), gk(NSO,NSO,NSO,NSO)
    real(dp) :: hdet(NDET,NDET), hkdet(NDET,NDET), ucsf(NDET,NDET), hcsf(NDET,NDET)
    real(dp) :: hbare(NDET,NDET), work36(NDET,NDET)
    real(dp) :: asing(NSING,NSING), atrip(NTRIP,NTRIP), aquint(1,1)
    real(dp) :: asing_keep(NSING,NSING), atrip_keep(NTRIP,NTRIP)
    real(dp) :: esing(NSING), etrip(NTRIP), equint(1)
    real(dp) :: esing_all(NSING,0:NSEAM-1), etrip_all(NTRIP,0:NSEAM-1)
    real(dp) :: equint_all(0:NSEAM-1)
    real(dp) :: ecore, c_h, c_ref, orth_err, cross_err, vxc_max
    real(dp) :: eref, thresh, aquint_diag, qvec(1,1)
    real(dp) :: pair_thc, pair_aniso, gauge_dev(0:NSEAM-1)
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
    c_ref = 1.0_dp
    if (is_dft) c_ref = infos%dft%hfscale
    c_h = infos%tddft%hfscale
    if (c_h < 0.0_dp) c_h = c_ref

    allocate(h_act(NACT,NACT), eri_act(NACT,NACT,NACT,NACT))
    call qmrsf_active_integrals(infos, NACT, act, ncore, h_act, eri_act, ecore, kcore_act)

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
    h_eff = h_act + (1.0_dp-c_ref)*kcore_act

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
    call native_pair_doublet(eri_act,PAIR_A,PAIR_B,pair_thc,pair_aniso)
    call native_dump_active(h_eff,eri_act,va_act,vb_act,c_h)
    hzero = 0.0_dp

    do iseam = 0, NSEAM-1
      call native_seam_mask(eri_act,iseam,PAIR_A,PAIR_B,eri_k)
      call native_build_spinorb(hzero,eri_k,1.0_dp,h1k,gk)
      call native_build_h(dets,h1k,gk,hkdet)
      hdet = hbare - (1.0_dp-c_h)*hkdet
      call native_add_vxc_diag(dets,va_act,vb_act,hdet)
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
    aquint_diag = equint_all(SEAM_NATIVE)
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
    write(iw,'(5x,a,f10.6)')    'Exact exchange in the reference        = ',c_ref
    write(iw,'(5x,a,f10.6)')    'Exact exchange in the dressed kernel   = ',c_h
    write(iw,'(/,5x,a)')        'Numerical checks'
    write(iw,'(5x,a,es12.4)')   'Spin-adaptation orthonormality error   = ',orth_err
    write(iw,'(5x,a,es12.4)')   'Largest neglected inter-multiplicity coupling = ',cross_err
    write(iw,'(5x,a,es12.4)')   'Largest active-space v_xc element      = ',vxc_max
    write(iw,'(/,5x,a)')        'Spin-resolved v_xc on the active orbitals (Hartree)'
    write(iw,'(7x,a)')          'orbital        alpha             beta'
    do i = 1, NACT
      write(iw,'(9x,i0,4x,2f16.10)') act(i),va_act(i,i),vb_act(i,i)
    end do
    write(iw,'(/,5x,a)') 'State energies are obtained as'
    write(iw,'(5x,a)')   '   E(state) = E_SCF(quintet reference) + omega(state),'
    write(iw,'(5x,a)')   'where omega is the response eigenvalue measured from the quintet'
    write(iw,'(5x,a)')   'configuration state function.'

    eref = infos%mol_energy%energy - aquint_diag
    thresh = infos%control%conf_print_threshold

    call native_print_block(infos, 'singlet', 'S', 0.0_dp, NSING, esing, asing, &
                            csf_label(1:NSING), eref, esing(1), thresh)
    call native_print_block(infos, 'triplet', 'T', 2.0_dp, NTRIP, etrip, atrip, &
                            csf_label(NSING+1:NSING+NTRIP), eref, esing(1), thresh)
    call native_print_block(infos, 'quintet', 'Q', 6.0_dp, 1, equint, qvec, &
                            csf_label(NDET:NDET), eref, esing(1), thresh)

    call native_seam_covariance(h_eff,eri_act,va_act,vb_act,c_h,dets,ucsf, &
                                esing_all,gauge_dev)
    call native_print_seam_compare(act,pair_thc,pair_aniso,esing_all,equint_all, &
                                   infos%mol_energy%energy,gauge_dev)

    call native_write_dump(c_h,c_ref,equint(1),esing,etrip,orth_err,cross_err, &
                           esing_all,etrip_all,equint_all,pair_aniso)
    write(iw,'(/,5x,a)') 'Machine-readable QMRSF-DK results written to qmrsf_dk_full_live.dat'

    deallocate(h_act,eri_act,cact)
    call flush(iw)
    close(iw)
  end subroutine qmrsf_dk_paper_native


  subroutine native_active_vxc(infos,mo,cact,va,vb)
    use types, only: information
    use dft, only: dft_initialize, dftclean, dftexcor
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
  pure subroutine native_pair_channel(d,ia,ib,ich,p)
    real(dp), intent(in) :: d(NACT,NACT)
    integer, intent(in) :: ia,ib,ich
    real(dp), intent(out) :: p(NACT,NACT)
    real(dp) :: h(NACT,NACT),dd(NACT,NACT),dx(NACT,NACT),tr,df
    integer :: t
    h = 0.0_dp; dd = 0.0_dp; dx = 0.0_dp
    do t = 1, NACT
      if (t /= ia .and. t /= ib) h(t,t) = d(t,t)
    end do
    tr = 0.5_dp*(d(ia,ia)+d(ib,ib))
    df = 0.5_dp*(d(ia,ia)-d(ib,ib))
    h(ia,ia) = tr;  h(ib,ib) = tr
    dd(ia,ia) = df; dd(ib,ib) = -df
    dx(ia,ib) = d(ia,ib); dx(ib,ia) = d(ib,ia)
    select case (ich)
    case (1); p = h
    case (2); p = dd
    case (3); p = dx
    case default; p = d - h - dd - dx
    end select
  end subroutine native_pair_channel


!> @brief Channel decomposition of the active two-electron tensor.
!> @detail ch(:,:,:,:,i,j) is the part of the tensor whose first
!>         charge-distribution index pair lies in channel i and whose second
!>         lies in channel j; the sixteen parts sum back to the tensor.
  subroutine native_seam_channels(g,ia,ib,ch)
    real(dp), intent(in) :: g(NACT,NACT,NACT,NACT)
    integer, intent(in) :: ia,ib
    real(dp), intent(out) :: ch(NACT,NACT,NACT,NACT,4,4)
    real(dp) :: p1(NACT,NACT,NACT,NACT,4), m(NACT,NACT), p(NACT,NACT)
    integer :: i,j,q,r,s,t
    do i = 1, 4
      do r = 1, NACT; do s = 1, NACT
        m = g(:,:,r,s)
        call native_pair_channel(m,ia,ib,i,p)
        p1(:,:,r,s,i) = p
      end do; end do
    end do
    do i = 1, 4; do j = 1, 4
      do t = 1, NACT; do q = 1, NACT
        m = p1(t,q,:,:,i)
        call native_pair_channel(m,ia,ib,j,p)
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
  subroutine native_seam_mask(eri,mode,ia,ib,eri_k)
    real(dp), intent(in) :: eri(NACT,NACT,NACT,NACT)
    integer, intent(in) :: mode,ia,ib
    real(dp), intent(out) :: eri_k(NACT,NACT,NACT,NACT)
    integer, parameter :: ICH_DD = 2, ICH_DX = 3, ICH_R = 4
    real(dp) :: ch(NACT,NACT,NACT,NACT,4,4)
    real(dp) :: ed(NACT,NACT), ex(NACT,NACT)
    real(dp) :: gdd,gxx,gdx,gxd,ndd,nxx,ndx,nxd,s2
    integer :: i,p,q,r,s

    if (mode == SEAM_NATIVE) then
      call native_exchange_eri(eri,eri_k)
      return
    end if

    call native_seam_channels(eri,ia,ib,ch)
    eri_k = ch(:,:,:,:,ICH_R,ICH_R)
    do i = ICH_DD, ICH_DX
      eri_k = eri_k + 0.5_dp*(ch(:,:,:,:,i,ICH_R) + ch(:,:,:,:,ICH_R,i))
    end do

    if (mode == SEAM_S3R) then
      eri_k = eri_k + 0.5_dp*(ch(:,:,:,:,ICH_DD,ICH_DD) + ch(:,:,:,:,ICH_DD,ICH_DX) &
                            + ch(:,:,:,:,ICH_DX,ICH_DD) + ch(:,:,:,:,ICH_DX,ICH_DX))
      return
    end if

    ! Reynolds projection of the doublet-doublet block.  In the orthonormal
    ! doublet directions e_d = (E_aa - E_bb)/sqrt2 and e_x = (E_ab + E_ba)/sqrt2
    ! the block is a 2x2 matrix G; averaging the mask over the SO(2) orbit sends
    ! G to (G + G^T + I tr G)/8, i.e. keeps half of the trace and a quarter of
    ! the traceless part.
    s2 = sqrt(2.0_dp)
    ed = 0.0_dp; ex = 0.0_dp
    ed(ia,ia) =  1.0_dp/s2; ed(ib,ib) = -1.0_dp/s2
    ex(ia,ib) =  1.0_dp/s2; ex(ib,ia) =  1.0_dp/s2
    gdd = 0.0_dp; gxx = 0.0_dp; gdx = 0.0_dp; gxd = 0.0_dp
    do p = 1, NACT; do q = 1, NACT; do r = 1, NACT; do s = 1, NACT
      gdd = gdd + ed(p,q)*ed(r,s)*eri(p,q,r,s)
      gxx = gxx + ex(p,q)*ex(r,s)*eri(p,q,r,s)
      gdx = gdx + ed(p,q)*ex(r,s)*eri(p,q,r,s)
      gxd = gxd + ex(p,q)*ed(r,s)*eri(p,q,r,s)
    end do; end do; end do; end do
    ndd = (3.0_dp*gdd + gxx)/8.0_dp
    nxx = (gdd + 3.0_dp*gxx)/8.0_dp
    ndx = 0.25_dp*gdx
    nxd = 0.25_dp*gxd
    do p = 1, NACT; do q = 1, NACT; do r = 1, NACT; do s = 1, NACT
      eri_k(p,q,r,s) = eri_k(p,q,r,s) &
                     + ndd*ed(p,q)*ed(r,s) + nxx*ex(p,q)*ex(r,s) &
                     + ndx*ed(p,q)*ex(r,s) + nxd*ex(p,q)*ed(r,s)
    end do; end do; end do; end do
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
    write(iw,'(A12,A20,A16,A12)') 'State', 'E(Hartree)', 'dE(eV)', '<S^2>'

    do istate = 1, nstates
      etot = eref + eig(istate)
      dele = (eig(istate) - e0)/UNITS_EV
      write(iw,'(3X,A,G0,t13,F20.10,F16.6,F12.4)') tag, istate-1, etot, dele, ssq
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


!> @brief Gauge-covariance gate for the seam conventions.
!> @detail Rotates the degenerate active pair by a fixed test angle, rebuilds the
!>         singlet spectrum from the rotated one- and two-electron integrals, and
!>         returns the largest excitation-energy change per seam.  The active
!>         orbitals are only re-expressed, so the physical spectrum cannot depend
!>         on the angle: the covariant seams must return zero to machine
!>         precision, while the value-based seam exposes the gauge freedom.  The
!>         v_xc diagonal is left alone because at a degenerate pair its block is
!>         proportional to the identity and therefore already rotation invariant.
  subroutine native_seam_covariance(h_eff,eri,va,vb,c_h,dets,ucsf,esing_ref,dev)
    use eigen, only: diag_symm_full
    use physical_constants, only: UNITS_EV
    use oqp_linalg
    real(dp), intent(in) :: h_eff(NACT,NACT), eri(NACT,NACT,NACT,NACT)
    real(dp), intent(in) :: va(NACT,NACT), vb(NACT,NACT), c_h
    integer, intent(in) :: dets(4,NDET)
    real(dp), intent(in) :: ucsf(NDET,NDET), esing_ref(NSING,0:NSEAM-1)
    real(dp), intent(out) :: dev(0:NSEAM-1)

    real(dp), parameter :: THETA_TEST = 0.6457718232379019_dp   ! 37 degrees
    real(dp) :: hrot(NACT,NACT), grot(NACT,NACT,NACT,NACT)
    real(dp) :: eri_k(NACT,NACT,NACT,NACT), hzero(NACT,NACT)
    real(dp) :: h1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp) :: h1k(NSO,NSO), gk(NSO,NSO,NSO,NSO)
    real(dp) :: hbare(NDET,NDET), hkdet(NDET,NDET), hdet(NDET,NDET)
    real(dp) :: work36(NDET,NDET), hcsf(NDET,NDET)
    real(dp) :: asing(NSING,NSING), esing(NSING)
    real(dp) :: r(NACT,NACT), c, s
    integer :: m, k, p, a, ierr

    c = cos(THETA_TEST); s = sin(THETA_TEST)
    r = 0.0_dp
    do p = 1, NACT; r(p,p) = 1.0_dp; end do
    r(PAIR_A,PAIR_A)= c; r(PAIR_A,PAIR_B)= s
    r(PAIR_B,PAIR_A)=-s; r(PAIR_B,PAIR_B)= c
    hrot = matmul(r,matmul(h_eff,transpose(r)))
    call native_rotate_pair(eri,PAIR_A,PAIR_B,THETA_TEST,grot)

    call native_build_spinorb(hrot,grot,1.0_dp,h1,g)
    call native_build_h(dets,h1,g,hbare)
    hzero = 0.0_dp
    dev = 0.0_dp
    do m = 0, NSEAM-1
      call native_seam_mask(grot,m,PAIR_A,PAIR_B,eri_k)
      call native_build_spinorb(hzero,eri_k,1.0_dp,h1k,gk)
      call native_build_h(dets,h1k,gk,hkdet)
      hdet = hbare - (1.0_dp-c_h)*hkdet
      call native_add_vxc_diag(dets,va,vb,hdet)
      call dgemm('N','T',NDET,NDET,NDET,1.0_dp,hdet,NDET,ucsf,NDET, &
                 0.0_dp,work36,NDET)
      call dgemm('N','N',NDET,NDET,NDET,1.0_dp,ucsf,NDET,work36,NDET, &
                 0.0_dp,hcsf,NDET)
      asing = hcsf(1:NSING,1:NSING)
      call diag_symm_full(0,NSING,asing,NSING,esing,ierr)
      if (ierr /= 0) then
        dev(m) = -1.0_dp
        cycle
      end if
      do k = 2, NSING
        dev(m) = max(dev(m), abs((esing(k)-esing(1)) &
                                -(esing_ref(k,m)-esing_ref(1,m)))/UNITS_EV)
      end do
    end do
  end subroutine native_seam_covariance


!> @brief Report the singlet spectrum under every seam convention.
!> @detail The value-based seam leaves the orientation of a degenerate active
!>         pair free, so its excitation energies are only defined up to that
!>         freedom; the two covariant seams remove it.  Their spread is the
!>         residual ambiguity of the convention itself.
  subroutine native_print_seam_compare(act,thc,aniso,esing_all,equint_all,escf,dev)
    use io_constants, only: iw
    use physical_constants, only: UNITS_EV
    integer, intent(in) :: act(NACT)
    real(dp), intent(in) :: thc,aniso,escf
    real(dp), intent(in) :: esing_all(NSING,0:NSEAM-1), equint_all(0:NSEAM-1)
    real(dp), intent(in) :: dev(0:NSEAM-1)
    integer :: k,m
    real(dp) :: de(0:NSEAM-1), spread

    write(iw,'(/,1x,78("^"))')
    write(iw,'(/,8x,a,/)') 'Seam convention comparison'
    write(iw,'(5x,a,i0,a,i0)') 'Degenerate-pair candidate (active orbitals) = ', &
         act(PAIR_A),', ',act(PAIR_B)
    write(iw,'(5x,a,f12.6)') 'Doublet anisotropy |Lambda| (Hartree)      = ',aniso
    write(iw,'(5x,a,f12.6)') 'Symmetry-frame angle of the pair (degrees) = ', &
         thc*180.0_dp/acos(-1.0_dp)
    write(iw,'(/,5x,a)') 'The value-based seam sends the doublet component d to the Hartree'
    write(iw,'(5x,a)')   'side and x to the exchange side, so its dressed elements move by'
    write(iw,'(5x,a)')   '(1/2) sin^2(2 theta) Lambda when the pair is rotated.  S3R and Haar'
    write(iw,'(5x,a)')   'weight the doublet isotropically and are therefore orientation'
    write(iw,'(5x,a)')   'independent; they coincide with the value-based seam at Lambda = 0'
    write(iw,'(5x,a)')   'and are inert at c_H = 1.'

    write(iw,'(/,5x,a)') 'Gauge-covariance gate: largest change of any singlet'
    write(iw,'(5x,a)')   'excitation energy when the pair is rotated by 37 degrees (eV).'
    write(iw,'(5x,a)')   'Meaningful only at an exact degeneracy, where the rotation is a'
    write(iw,'(5x,a)')   'symmetry of the reference; elsewhere the diagonal v_xc patch is'
    write(iw,'(5x,a)')   'itself orientation dependent and every convention moves.'
    do m = 0, NSEAM-1
      write(iw,'(7x,a8,es14.4)') SEAM_NAME(m), dev(m)
    end do

    ! Each convention has its own quintet configuration state function, so the
    ! energy law E(state) = E_SCF + omega(state) must subtract that convention's
    ! own A_quintet, not the value-based one.
    write(iw,'(/,5x,a)') 'Ground-state total energy (Hartree)'
    do m = 0, NSEAM-1
      write(iw,'(7x,a8,f22.10)') SEAM_NAME(m), escf + esing_all(1,m) - equint_all(m)
    end do

    write(iw,'(/,5x,a)') 'Singlet excitation energies (eV)'
    write(iw,'(7x,a,3(a10),a12)') 'State ', &
         (SEAM_NAME(m), m=0,NSEAM-1), 'S3R-Haar'
    do k = 2, min(NSING,7)
      do m = 0, NSEAM-1
        de(m) = (esing_all(k,m)-esing_all(1,m))/UNITS_EV
      end do
      spread = de(SEAM_S3R)-de(SEAM_HAAR)
      write(iw,'(7x,a1,i0,3x,3f10.4,f12.4)') 'S',k-1,(de(m),m=0,NSEAM-1),spread
    end do
    write(iw,'(1x,78("=")/)')
  end subroutine native_print_seam_compare


!> @brief Dump the active-space inputs of the seam so that an independent
!>        implementation can be applied to exactly the same integrals.
  subroutine native_dump_active(h_eff,eri,va,vb,c_h)
    real(dp), intent(in) :: h_eff(NACT,NACT), eri(NACT,NACT,NACT,NACT)
    real(dp), intent(in) :: va(NACT,NACT), vb(NACT,NACT), c_h
    integer :: u,p,q,r,s
    open(newunit=u,file='qmrsf_dk_active.dat',status='replace',action='write')
    write(u,'(a,1x,i0)') 'QMRSF_DK_ACTIVE_V1',NACT
    write(u,'(es24.16)') c_h
    do p = 1, NACT; write(u,'(*(es24.16,1x))') h_eff(p,:); end do
    do p = 1, NACT; write(u,'(*(es24.16,1x))') va(p,:); end do
    do p = 1, NACT; write(u,'(*(es24.16,1x))') vb(p,:); end do
    do p = 1, NACT; do q = 1, NACT; do r = 1, NACT; do s = 1, NACT
      write(u,'(es24.16)') eri(p,q,r,s)
    end do; end do; end do; end do
    close(u)
  end subroutine native_dump_active


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
