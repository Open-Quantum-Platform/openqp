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

contains

  subroutine qmrsf_dk_paper_native(infos)
    use types, only: information
    use io_constants, only: iw
    use printing, only: print_module_info
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use qmrsf_ao2mo_mod, only: qmrsf_active_integrals
    use eigen, only: diag_symm_full

    type(information), target, intent(inout) :: infos

    real(dp), contiguous, pointer :: mo_a(:,:)
    real(dp), allocatable :: h_act(:,:), eri_act(:,:,:,:), cact(:,:)
    real(dp) :: va_act(NACT,NACT), vb_act(NACT,NACT)
    real(dp) :: h_eff(NACT,NACT), kcore_act(NACT,NACT)
    real(dp) :: hzero(NACT,NACT), eri_k(NACT,NACT,NACT,NACT)
    real(dp) :: h1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp) :: h1k(NSO,NSO), gk(NSO,NSO,NSO,NSO)
    real(dp) :: hdet(NDET,NDET), hkdet(NDET,NDET), ucsf(NDET,NDET), hcsf(NDET,NDET)
    real(dp) :: asing(NSING,NSING), atrip(NTRIP,NTRIP), aquint(1,1)
    real(dp) :: esing(NSING), etrip(NTRIP), equint(1)
    real(dp) :: ecore, c_h, c_ref, orth_err, cross_err, vxc_max
    integer :: dets(4,NDET), act(NACT)
    integer :: nbf, ncore, i, j, ierr
    logical :: is_dft

    open(unit=iw, file=infos%log_filename, position='append')
    call print_module_info('QMRSF_DK_NATIVE', &
         'Paper spin-adapted QMRSF-DK from OpenQP orbitals, Fock, v_xc, and ERIs')

    nbf = int(infos%basis%nbf)
    ncore = int(infos%mol_prop%nelec_B)
    if (int(infos%mol_prop%nelec_A) - ncore /= NACT) then
      write(iw,'(/,5x,a)') 'QMRSF-DK native requires a high-spin quintet reference '// &
           '(N_alpha-N_beta=4).'
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
    call native_build_spinorb(h_eff,eri_act,1.0_dp,h1,g)
    call native_gen_dets(dets)
    call native_build_h(dets,h1,g,hdet)
    call native_exchange_eri(eri_act,eri_k)
    hzero = 0.0_dp
    call native_build_spinorb(hzero,eri_k,1.0_dp,h1k,gk)
    call native_build_h(dets,h1k,gk,hkdet)
    hdet = hdet - (1.0_dp-c_h)*hkdet
    call native_add_vxc_diag(dets,va_act,vb_act,hdet)
    call native_build_ucsf(dets,ucsf,orth_err)
    hcsf = matmul(ucsf,matmul(hdet,transpose(ucsf)))

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

    asing = hcsf(1:NSING,1:NSING)
    atrip = hcsf(NSING+1:NSING+NTRIP,NSING+1:NSING+NTRIP)
    aquint = hcsf(NDET:NDET,NDET:NDET)
    call diag_symm_full(0,NSING,asing,NSING,esing,ierr)
    if (ierr /= 0) then
      write(iw,'(5x,a,i0)') 'QMRSF-DK native singlet diagonalization failed, info=',ierr
      close(iw)
      return
    end if
    call diag_symm_full(0,NTRIP,atrip,NTRIP,etrip,ierr)
    if (ierr /= 0) then
      write(iw,'(5x,a,i0)') 'QMRSF-DK native triplet diagonalization failed, info=',ierr
      close(iw)
      return
    end if
    equint(1) = aquint(1,1)

    vxc_max = max(maxval(abs(va_act)),maxval(abs(vb_act)))
    write(iw,'(/,5x,a)') '================ QMRSF-DK native paper matrix ================'
    write(iw,'(5x,a,i0)')       'basis functions                 = ',nbf
    write(iw,'(5x,a,i0)')       'inactive doubly occupied MOs    = ',ncore
    write(iw,'(5x,a,4(i0,1x))') 'active SOMO indices             = ',act
    write(iw,'(5x,a,f10.6)')    'reference exact exchange c_ref  = ',c_ref
    write(iw,'(5x,a,f10.6)')    'response exact exchange c_H     = ',c_h
    write(iw,'(5x,a,es12.4)')   'hard-coded U orthonormal error  = ',orth_err
    write(iw,'(5x,a,es12.4)')   'discarded cross-spin block max = ',cross_err
    write(iw,'(5x,a,es12.4)')   'max active spin v_xc            = ',vxc_max
    do i = 1, NACT
      write(iw,'(5x,a,i0,a,2f16.10)') 'active v_xc diagonal, orbital ',i, &
           ' (alpha,beta) = ',va_act(i,i),vb_act(i,i)
    end do
    write(iw,'(5x,a)') 'Energy law: E_I = E_ROKS(quintet) + lambda_I - A_quintet'
    write(iw,'(/,5x,a)') 'Lowest singlet excitation energies relative to S0 (eV):'
    do i = 1, min(10,NSING)
      write(iw,'(7x,a,i2,a,f12.6)') 'S',i-1,'  ',(esing(i)-esing(1))*27.211386245988_dp
    end do
    write(iw,'(/,5x,a)') 'Lowest triplet energies relative to singlet S0 (eV):'
    do i = 1, min(6,NTRIP)
      write(iw,'(7x,a,i2,a,f12.6)') 'T',i-1,'  ',(etrip(i)-esing(1))*27.211386245988_dp
    end do

    call native_write_dump(c_h,c_ref,equint(1),esing,etrip,orth_err,cross_err)
    write(iw,'(/,5x,a)') 'QMRSF-DK native results written to qmrsf_dk_full_live.dat.'

    deallocate(h_act,eri_act,cact)
    call flush(iw)
    close(iw)
  end subroutine qmrsf_dk_paper_native


  subroutine native_active_vxc(infos,mo,cact,va,vb)
    use types, only: information
    use dft, only: dft_initialize, dftclean, dftexcor
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: unpack_matrix
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
    tmp = matmul(fsq,cact)
    va = matmul(transpose(cact),tmp)
    call unpack_matrix(fbp,fsq,nbf,'U')
    tmp = matmul(fsq,cact)
    vb = matmul(transpose(cact),tmp)
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


  subroutine native_build_ucsf(dets,u,orth_err)
    integer, intent(in) :: dets(4,NDET)
    real(dp), intent(out) :: u(NDET,NDET),orth_err
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

    do i=1,NDET
      if (native_os_class(cfg(i))==0) then
        rs=rs+1; u(rs,i)=real(phase(i),dp)
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
      rt=rt+1
      u(rt,i)=real(phase(i),dp)/sq2
      u(rt,j)=real(phase(j),dp)/sq2
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
    end do
    do k=3,5
      rt=rt+1
      do i=1,6; u(rt,fo(i))=c4(i,k)*real(phase(fo(i)),dp); end do
    end do
    do i=1,6; u(NDET,fo(i))=c4(i,6)*real(phase(fo(i)),dp); end do

    if (rs/=NSING .or. rt/=NSING+NTRIP) &
      error stop 'QMRSF-DK native: incorrect spin-adapted block dimensions'
    gram=matmul(u,transpose(u))
    do i=1,NDET; gram(i,i)=gram(i,i)-1.0_dp; end do
    orth_err=maxval(abs(gram))
  end subroutine native_build_ucsf


  subroutine native_write_dump(c_h,c_ref,aq,es,et,orth_err,cross_err)
    real(dp), intent(in) :: c_h,c_ref,aq,es(NSING),et(NTRIP),orth_err,cross_err
    integer :: u
    open(newunit=u,file='qmrsf_dk_full_live.dat',status='replace',action='write')
    write(u,'(a,4(1x,i0))') 'QMRSF_DK_NATIVE_V1',NACT,NSING,NTRIP,1
    write(u,'(2es24.16)') c_h,c_ref
    write(u,'(es24.16)') aq
    write(u,'(*(es24.16,1x))') es
    write(u,'(*(es24.16,1x))') et
    write(u,'(es24.16)') aq
    write(u,'(2es24.16)') orth_err,cross_err
    close(u)
  end subroutine native_write_dump

end module qmrsf_dk_paper_native_mod
