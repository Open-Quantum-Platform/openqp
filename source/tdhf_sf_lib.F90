module tdhf_sf_lib

    use oqp_linalg

contains

  subroutine sfroesum(fazzfb,pmo,noca,nocb,ivec)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(in), dimension(:,:) :: fazzfb
    real(kind=dp), intent(inout), dimension(:,:) :: pmo
    integer, intent(in) :: noca, nocb, ivec

    integer :: i, ij, j, nbf

    nbf = ubound(fazzfb, 1)

    ij = 0
    do j = nocb+1, nbf
      do i = 1, noca
        ij = ij+1
        pmo(ij,ivec) = pmo(ij,ivec)+fazzfb(i,j)
      end do
    end do
  end subroutine sfroesum

  subroutine sfresvec(q,a,b,vec,eigv,nvec,rnorm,ndsr)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: q
    real(kind=dp), intent(in), dimension(:,:) :: a, b
    real(kind=dp), intent(inout), dimension(:,:) :: vec
    real(kind=dp), intent(in), dimension(:) :: eigv
    integer, intent(in) :: nvec
    real(kind=dp), intent(out), dimension(:) :: rnorm
    integer, intent(in) :: ndsr

    integer :: ist, xvec_dim

    xvec_dim = ubound(q, 1)

    call dgemm('n','n',xvec_dim,ndsr,nvec, &
               1.0_dp,b,xvec_dim, &
                      vec,nvec, &
               0.0_dp,q,xvec_dim)

    do ist = 1, ndsr
      vec(:,ist) = -vec(:,ist)*eigv(ist)
    end do

    call dgemm('n','n',xvec_dim,ndsr,nvec, &
               1.0_dp,a,xvec_dim, &
                      vec,nvec, &
               1.0_dp,q,xvec_dim)

    do ist = 1, ndsr
      rnorm(ist) = dot_product(q(:,ist),q(:,ist))
    end do

  end subroutine sfresvec

  subroutine sfqvec(q,xm,eigv,ndsr)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: q
    real(kind=dp), intent(in), dimension(:) :: xm, eigv
    integer, intent(in) :: ndsr

    integer :: ii, ist, xvec_dim
    real(kind=dp) :: sign, val1, val2

    xvec_dim = ubound(xm, 1)

    do ist = 1, ndsr
      do ii = 1, xvec_dim
        val1 = eigv(ist)-xm(ii)
        val2 = abs(val1)
        if( val2<1.0D-12 )then
          val1 = 1.0D-05
        else if( val2<1.0D-05 )then
          sign = val2/val1
          val1 = sign*1.0D-05
        end if
        q(ii,ist) = q(ii,ist)/val1
      end do
    end do
  end subroutine sfqvec

  subroutine sfesum(eiga,eigb,pmo,z,noca,nocb,ivec)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(in) :: eiga(:), eigb(:)
    real(kind=dp), intent(inout) :: pmo(:,:)
    real(kind=dp), intent(in) :: z(:,:)
    integer, intent(in) :: noca, nocb, ivec

    integer :: i, ij, j, nbf

    nbf = ubound(eiga, 1)

!   ----- add (ea-ei)*zai -----
    ij = 0
    do j=nocb+1,nbf
      do i=1,noca
        ij = ij+1
        pmo(ij,ivec) = pmo(ij,ivec)+(eigb(j)-eiga(i))*z(ij,ivec)
      end do
    end do
  end subroutine sfesum

  subroutine trfrmb(bvec,vec,nvec,ndsr)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: bvec
    real(kind=dp), intent(in), dimension(:,:) :: vec
    integer, intent(in) :: nvec, ndsr

    real(kind=dp), allocatable, dimension(:,:) :: scr
    integer :: xvec_dim

    xvec_dim = ubound(bvec, 1)
    allocate(scr(xvec_dim,ndsr), &
             source=0.0_dp)

    scr = bvec
    ! Get a new Bvec
    call dgemm('n', 'n', xvec_dim, ndsr, nvec, &
               1.0_dp, scr, xvec_dim, &
                       vec, nvec,&
               0.0_dp, bvec, xvec_dim)
  ! do ii = 1, xvec_dim
  !   do jj = 1, ndsr
  !   bvec(ii,jj) = 0.0_dp
  !     do kk = 1, nvec
  !       bvec(ii,jj) = bvec(ii,jj)+scr(ii,kk)*vec(kk,jj)
  !     end do
  !   end do
  ! end do

    deallocate(scr)
  end subroutine trfrmb

  subroutine sfdmat(bvec,abxc,mo_a,ta,tb, &
                    noca,nocb)
    use precision, only: dp
    use tdhf_lib, only: iatogen
    use mathlib, only: pack_matrix
    use mathlib, only: orthogonal_transform

    implicit none

    real(kind=dp), intent(in), dimension(:) :: bvec
    real(kind=dp), intent(in), dimension(:,:) :: mo_a
    real(kind=dp), intent(inout), dimension(:,:) :: abxc
    real(kind=dp), intent(out), dimension(:) :: ta, tb
    integer, intent(in) :: noca, nocb

    integer :: nvirb, nbf, nbf_tri, xvec_dim
    real(kind=dp), allocatable, dimension(:,:) :: scr1, scr2

    nbf = ubound(mo_a,1)
    nbf_tri = ubound(ta,1)
    xvec_dim = ubound(bvec,1)
    allocate(scr1(nbf,nbf), &
             scr2(nbf,nbf), &
             source=0.0_dp)

  ! MO(I+,A-) -> AO(M,N)
    nvirb = nbf-nocb

    call iatogen(bvec,scr1,noca,nocb)
    call orthogonal_transform('t', nbf, mo_a, scr1, abxc, scr2)

  ! Unrelaxed difference density matrix -----

  ! OCC(Alpha)-OCC(Alpha)
    call dgemm('n','t',noca,noca,nvirb, &
              -1.0_dp,bvec,noca, &
                      bvec,noca, &
               0.0_dp,scr1,noca)

  ! MO(I+,J+) -> AO(M,N)
    call dgemm('n','n',nbf,noca,noca, &
               1.0_dp,mo_a,nbf, &
                      scr1,noca, &
               0.0_dp,scr2,nbf)
    call dgemm('n','t',nbf,nbf,noca, &
               1.0_dp,scr2,nbf, &
                      mo_a,nbf, &
               0.0_dp,scr1,nbf)
    call pack_matrix(scr1,ta)

    call dgemm('t','n',nvirb,nvirb,noca, &
               1.0_dp,bvec,noca, &
                      bvec,noca, &
               0.0_dp,scr1,nvirb)

  ! MO(A-,B-) -> AO(M,N)
    call dgemm('n','n',nbf,nvirb,nvirb, &
               1.0_dp,mo_a(:,nocb+1:),nbf, &
                      scr1,nvirb, &
               0.0_dp,scr2,nbf)
    call dgemm('n','t',nbf,nbf,nvirb, &
               1.0_dp,scr2,nbf, &
                      mo_a(:,nocb+1:),nbf, &
               0.0_dp,scr1,nbf)
    call pack_matrix(scr1,tb)

    deallocate(scr1,scr2)
  end subroutine sfdmat

  subroutine get_transitions(trans, noca, nocb, nbf)

    implicit none

    integer, intent(out), dimension(:,:) :: trans
    integer, intent(in) :: noca, nocb, nbf

    integer :: ij, i, j
    ij = 0
    do j = nocb+1, nbf
       do i = 1, noca
        ij = ij+1
        trans(ij,1) = i
        trans(ij,2) = j
      end do
    end do

  end subroutine get_transitions

  subroutine print_results(infos, bvec_mo, excitation_energy, &
                           trans, dip, spin_square, nstates)
    use precision, only: dp
    use types, only: information
    use physical_constants, only: toev => ev2htree

    implicit none

    type(information), target, intent(in) :: infos
    real(kind=dp), intent(in), dimension(:,:) :: bvec_mo
    real(kind=dp), intent(in), dimension(:) :: excitation_energy
    integer, intent(in), dimension(:,:) :: trans
    real(kind=dp), intent(in), dimension(:,:,:) :: dip
    real(kind=dp), intent(in), dimension(:) :: spin_square
    integer, intent(in) :: nstates

    integer :: istat, jstat, ij, i, j, nocca, noccb, xvec_dim, ndeex
    real(kind=dp) :: ydum, xdum, threshold, ROHF_energy, energ, f

    threshold = infos%control%conf_print_threshold
    xvec_dim = ubound(bvec_mo, 1)

    ROHF_energy = infos%mol_energy%energy
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B

    do istat=1,nstates
      ydum = toev*excitation_energy(istat)
      write(*,'(/,1x,"State #",I4,2X,"Energy =",F12.6,1X,"eV")') istat, ydum
!     write(*,'(3x,"Symmetry of state =",4x,a)') '?a?a?'
      write(*,'(15x,"<S^2> =",1x,f9.4)') spin_square(istat)
      write(*,'(8x,"DRF",4x,"Coeff",8x,"OCC",7x,"VIR")')
      write(*,'(8x,3("-"),2x,8("-"),5x,6("-"),4x,6("-"))')
      do ij=1,xvec_dim
        i = trans(ij,1)
        j = trans(ij,2)
        xdum = bvec_mo(ij,istat)
        if (abs(xdum)>threshold) then
          write(*,'(7x,i4,1x,f9.6,6x,i4,2x,"->",2x,i4,2x)') ij,xdum,i,j
        end if
      end do
    end do

    write(*,'(/5x, "Summary table",/)')
    write(*,'(1x, "State", 6x, "Energy", 7x,"Excitation", 3x, "Excitation(eV)", &
             &2x, "<S^2>", 9x, "Transition dipole moment, a.u.",&
             &8x, "Oscillator")')
    write(*,'(11x, "Hartree", 11x, "eV", 10x, "rel. GS" &
             &18x, "X", 10x, "Y", 10x, "Z", 8x,"Abs.", 6x, "strength")')

    ndeex = 0
    do istat = 1, nstates
      if(excitation_energy(istat)<0.0_dp) ndeex = ndeex+1
    end do

    ! De-excitation
    do istat = 1, ndeex
      energ = excitation_energy(istat)-excitation_energy(1)
      f = 2.0d0 / 3.0d0 * (energ) * sum(dip(:,1,istat)**2)
      write(*,'(x, i3, 1x, f17.10, 2f13.6, 6x, &
               &f5.3, 4(1x,e10.3),2x,e10.3)') &
           istat, ROHF_energy+excitation_energy(istat), toev*excitation_energy(istat), &
           toev*energ, spin_square(istat), dip(1:3,1,istat), sqrt(sum(dip(:,1,istat)**2)), f
    end do

    ! Reference ROHF state
    write(*,'(1x, i3, 1x, f17.10, 2f13.6, 8x,&
            &"(ROHF/UHF Reference state)")') 0, ROHF_energy, 0.0_dp, -excitation_energy(1)*toev

    ! Excitation
    do istat=ndeex+1,nstates
      energ = excitation_energy(istat)-excitation_energy(1)
      f = 2.0d0 / 3.0d0 * (energ) * sum(dip(:,1,istat)**2)
      write(*,'(x, i3, 1x, f17.10, 2f13.6, 6x, &
               &f5.3, 4(1x,e10.3),2x,e10.3)') &
           istat, ROHF_energy+excitation_energy(istat), toev*excitation_energy(istat), &
           toev*energ, spin_square(istat), dip(1:3,1,istat), sqrt(sum(dip(:,1,istat)**2)), f
    end do
    write(*,*)
    write(*,"(2x,'Transition',3x,'Excitation',9x,'Transition dipole, a.u.',19x,'Oscillator',&
          &/18x,'eV',14x,'x',10x,'y',10x,'z',9x,'Abs.',7x,'strength')")
    do istat=1, nstates
       do jstat=istat+1, nstates
          energ = excitation_energy(jstat)-excitation_energy(istat)
          f = 2.0d0 / 3.0d0 * (energ) * sum(dip(:,istat,jstat)**2)
    write(*,"(3x,i0,1x,'->',1x,i0,t11,3x,f11.6,3x,3e11.3,1x,e11.3,2x,e11.3)") &
             istat,jstat,toev*energ,dip(1:3,istat,jstat), sqrt(sum(dip(:,istat,jstat)**2)), f
       enddo
    enddo
    write(*,*)

  end subroutine print_results

  subroutine sfrorhs(rhs,xhxa,xhxb,hpta,hptb,Tij,Tab,Fa,Fb, &
                     noca,nocb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:) :: rhs
    real(kind=dp), intent(inout), dimension(:,:) :: xhxa, xhxb
    real(kind=dp), intent(in), dimension(:,:) :: hpta
    real(kind=dp), intent(in), dimension(:,:) :: hptb
    real(kind=dp), intent(in), dimension(:,:) :: tij
    real(kind=dp), intent(in), dimension(:,:) :: tab
    real(kind=dp), intent(in), dimension(:,:) :: fa, fb
    integer, intent(in) :: noca, nocb

    real(kind=dp), allocatable, dimension(:,:) :: wrk
    integer :: nbf, i, j, ij, k, nconf

    nbf = ubound(fa, 1)

    allocate(wrk(nbf,nbf), &
             source=0.0_dp)

  ! HPTA --> AB1_MO(1)
  ! HPTB --> AB1_MO(2)
  ! TA   --> TIJ
  ! TB   --> TAB

  ! Alpha
  ! XHXA+= 2*FA(P+,I+)*TA(I+,J+)
    call dgemm('n', 'n', nbf, noca, noca, &
               2.0_dp, fa, nbf, &
                       tij, noca, &
               1.0_dp, xhxa, nbf)

  ! Beta
  ! XHXB+= 2*FB(P-,A-)*TB(A-,B-)
    do j = nocb+1, nbf
      do i = nocb+1, nbf
        wrk(i,j) = tab(i-nocb,j-nocb)
      end do
    end do
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, fb, nbf, &
                       wrk, nbf, &
               1.0_dp, xhxb, nbf)

  ! doc-socc
    ij = 0
    do i = nocb+1, noca
      do j = 1, nocb
        ij = ij+1
        rhs(ij) = hptb(j,i-nocb)+xhxa(i,j)-xhxa(j,i)-xhxb(j,i)
      end do
    end do

  ! doc-virt
    do k = noca+1, nbf
      do j = 1, nocb
        ij = ij+1
        rhs(ij) = hpta(j,k-noca)+hptb(j,k-nocb)+xhxa(k,j)-xhxb(j,k)
      end do
    end do

  ! soc-virt
    do k = noca+1, nbf
      do i = nocb+1, noca
        ij = ij+1
        rhs(ij) = hpta(i,k-noca)+xhxa(k,i)+xhxb(k,i)-xhxb(i,k)
      end do
    end do

  ! Multiplied by -1 i.e., RHS of Z-vector eq. -----
    nconf = ij
    rhs(1:nconf) = -rhs(1:nconf)

  end subroutine sfrorhs

  subroutine sfromcal(xm,xminv,energy,fa,fb,noca,nocb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:) :: xm
    real(kind=dp), intent(out), dimension(:) :: xminv
    real(kind=dp), intent(in), dimension(:) :: energy
    real(kind=dp), intent(in), dimension(:,:) :: fa
    real(kind=dp), intent(in), dimension(:,:) :: fb
    integer, intent(in) :: noca, nocb

    integer :: ij, i, j, k, nbf, nsoc, lzdim, nvira

    nbf = ubound(fa, 1)

    nvira = nbf-noca
    nsoc = noca-nocb
    lzdim = nocb*(nsoc+nvira)+nsoc*nvira

  ! doc-socc
    ij = 0
    do i = nocb+1, noca
      do j = 1, nocb
        ij = ij+1
        xm(ij) = (fb(i,i)-fb(j,j))*0.5_dp
      end do
    end do

  ! DOC-VIRT
    do k = noca+1, nbf
      do j = 1, nocb
        ij = ij+1
        xm(ij) = energy(k)-energy(j)
      end do
    end do

  ! SOCC-VIRT
    do k = noca+1, nbf
      do i = nocb+1, noca
        ij = ij+1
        xm(ij) = (fa(k,k)-fa(i,i))*0.5_dp
      end do
    end do

    do j = 1, lzdim
      xminv(j)=1.0_dp/xm(j)
    end do

  end subroutine sfromcal

  subroutine sfrogen(ava,avb,pv,noca,nocb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: ava
    real(kind=dp), intent(out), dimension(:,:) :: avb
    real(kind=dp), intent(in), dimension(:) :: pv
    integer, intent(in) :: noca, nocb

    integer :: ij, i, j, k, nbf

    nbf = ubound(ava, 1)

    ava = 0.0_dp
    avb = 0.0_dp

  ! doc-socc
    ij = 0
    do i = nocb+1, noca
      do j = 1, nocb
        ij = ij+1
        avb(j,i) = pv(ij)
      end do
    end do

  ! doc-virt
    do k = noca+1, nbf
      do j = 1, nocb
        ij = ij+1
        ava(j,k) = pv(ij)
        avb(j,k) = pv(ij)
      end do
    end do

  ! socc-virt
    do k = noca+1, nbf
      do i = nocb+1, noca
        ij = ij+1
        ava(i,k) = pv(ij)
      end do
    end do

  end subroutine sfrogen

  subroutine sfrolhs(pmo, z, e, fa, fb, hpza, hpzb,  &
                     noca, nocb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:) :: pmo
    real(kind=dp), intent(in), dimension(:) :: z
    real(kind=dp), intent(in), dimension(:) :: e
    real(kind=dp), intent(in), dimension(:,:) :: fa, fb
    real(kind=dp), intent(in), dimension(:,:) :: hpza
    real(kind=dp), intent(in), dimension(:,:) :: hpzb
    integer, intent(in) :: noca, nocb

    real(kind=dp), allocatable, dimension(:,:) :: ztmp
    real(kind=dp), allocatable, dimension(:,:) :: wrk
    integer :: ij, i, k, j, nbf, lr1, lr2

    nbf = ubound(fa, 1)
    lr1 = nocb+1
    lr2 = noca
    allocate(ztmp(nbf,nbf), &
             wrk(nbf,nbf), &
             source=0.0_dp)

    ij = 0
    do i = nocb+1, noca
      do j = 1, nocb
        ij = ij+1
        ztmp(j,i) = z(ij)
      end do
    end do

    do k = noca+1, nbf
      do j = 1, nocb
        ij = ij+1
        ztmp(j,k) = z(ij)
      end do
    end do

    do k = noca+1, nbf
      do i = nocb+1, noca
        ij = ij+1
        ztmp(i,k) = z(ij)
      end do
    end do

  ! doc-socc
    do j = 1, nocb
      wrk(j,1) = wrk(j,1)+hpzb(j,1) &
                         -fa(lr1,lr1)*ztmp(j,lr1) &
                         -fa(lr2,lr1)*ztmp(j,lr2)
      wrk(j,2) = wrk(j,2)+hpzb(j,2) &
                         -fa(lr2,lr2)*ztmp(j,lr2) &
                         -fa(lr1,lr2)*ztmp(j,lr1)
    end do

    do j = 1, nocb
      do k = 1, nocb
        wrk(j,1) = wrk(j,1)+fa(k,j)*ztmp(k,lr1)
        wrk(j,2) = wrk(j,2)+fa(k,j)*ztmp(k,lr2)
      end do
    end do

    do j = 1, nocb
      do k = 1, nbf-noca
        wrk(j,1) = wrk(j,1)+fb(noca+k,lr1)*ztmp(j,noca+k) &
                           +fb(noca+k,j)*ztmp(lr1,noca+k)
        wrk(j,2) = wrk(j,2)+fb(noca+k,j)*ztmp(lr2,noca+k) &
                           +fb(noca+k,lr2)*ztmp(j,noca+k)
      end do
    end do

    ij = 0
    wrk = wrk*0.5_dp
    do i = 1, 2
      do j = 1, nocb
        ij = ij+1
        pmo(ij) = (e(nocb+i)-e(j))*z(ij)+wrk(j,i)
      end do
    end do

  ! doc-virt
    wrk = 0.0_dp
    do k = 1, nbf-noca
      do j = 1, nocb
        wrk(j,k) = wrk(j,k)+hpza(j,k) &
                           +hpzb(j,noca-nocb+k) &
                           +fb(lr1,noca+k)*ztmp(j,lr1) &
                           +fb(lr2,noca+k)*ztmp(j,lr2) &
                           -fa(lr1,j)*ztmp(lr1,noca+k) &
                           -fa(lr2,j)*ztmp(lr2,noca+k)
      end do
    end do

    wrk = wrk*0.5_dp
    do k = 1, nbf-noca
      do j = 1, nocb
        ij = ij+1
        pmo(ij) =(e(noca+k)-e(j))*z(ij)+wrk(j,k)
      end do
    end do

  ! socc-virt
    wrk = 0.0_dp
    do k = 1, nbf-noca
      wrk(k,1) = wrk(k,1)+hpza(lr1,k) &
                         +fb(lr1,lr1)*ztmp(lr1,noca+k) &
                         +fb(lr2,lr1)*ztmp(lr2,noca+k)
      wrk(k,2) = wrk(k,2)+hpza(lr2,k) &
                         +fb(lr1,lr2)*ztmp(lr1,noca+k) &
                         +fb(lr2,lr2)*ztmp(lr2,noca+k)
    end do

    do k = 1, nbf-noca
      do j = 1, nocb
        wrk(k,1) = wrk(k,1)-fa(j,noca+k)*ztmp(j,lr1) &
                           -fa(j,lr1)*ztmp(j,noca+k)
        wrk(k,2) = wrk(k,2)-fa(j,noca+k)*ztmp(j,lr2) &
                           -fa(j,lr2)*ztmp(j,noca+k)
      end do
    end do

    do k = 1, nbf-noca
      do j = 1, nbf-noca
        wrk(k,1) = wrk(k,1)-fb(noca+j,noca+k)*ztmp(lr1,noca+j)
        wrk(k,2) = wrk(k,2)-fb(noca+j,noca+k)*ztmp(lr2,noca+j)
      end do
    end do

    wrk = wrk*0.5_dp
    do k = 1, nbf-noca
      do i = 1, noca-nocb
        ij = ij+1
        pmo(ij) = (e(noca+k)-e(nocb+i))*z(ij)+wrk(k,i)
      end do
    end do

    deallocate(ztmp, wrk)
  end subroutine sfrolhs

  subroutine pcgrbpini(r, pk, error, d, xm_in, a_pk)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:) :: r
    real(kind=dp), intent(out), dimension(:) :: pk
    real(kind=dp), intent(out) :: error
    real(kind=dp), intent(in), dimension(:) :: d
    real(kind=dp), intent(in), dimension(:) :: xm_in
    real(kind=dp), intent(in), dimension(:) :: a_pk

    real(kind=dp) :: beta

  ! R ini and R norm(error)
    r = d-a_pk
    error = dot_product(r, r)

  ! Beta ini
    beta = 1.0_dp/dot_product(r**2, xm_in)

  ! pk ini
    pk = beta*xm_in*r

  end subroutine pcgrbpini

  subroutine pcgb(pk,r,xm_in)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out) :: pk(:)
    real(kind=dp), intent(in) :: r(:)
    real(kind=dp), intent(in) :: xm_in(:)

    real(kind=dp) :: beta

    beta = 1.0_dp/sum(r*r*xm_in)

  ! pk ini
    pk = pk + beta*xm_in*r

  end subroutine pcgb

  subroutine sfropcal(pa, pb, ta, tb, z, &
                      noca, nocb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: pa, pb
    real(kind=dp), intent(in), dimension(:,:) :: ta, tb
    real(kind=dp), intent(in), dimension(:) :: z
    integer, intent(in) :: noca, nocb

    integer :: i, j, k, ij, nbf


    nbf = ubound(pa, 1)

  ! Alpha
    pa = 0.0_dp
    do j = 1, noca
      do i = 1, noca
        pa(i,j) = ta(i,j)
      end do
    end do

    pb = 0.0_dp
    do j = nocb+1, nbf
      do i = nocb+1, nbf
        pb(i,j) = tb(i-nocb,j-nocb)
      end do
    end do

  ! add Z contribution

  ! DOC-SOCC
    ij = 0
    do i = nocb+1, noca
      do j = 1, nocb
        ij = ij+1
        pb(j,i) = pb(j,i)+z(ij)*0.5_dp
      end do
    end do

  ! DOC-VIRT
    do k = noca+1, nbf
      do j = 1, nocb
        ij = ij+1
        pa(j,k) = pa(j,k)+z(ij)*0.5_dp
        pb(j,k) = pb(j,k)+z(ij)*0.5_dp
      end do
    end do

  ! SOCC-VIRT
    do k = noca+1, nbf
      do i = nocb+1, noca
        ij = ij+1
        pa(i,k) = pa(i,k)+z(ij)*0.5_dp
      end do
    end do

  end subroutine sfropcal

  subroutine sfrowcal(wmo, target_energy, mo_energy_a, fa, fb, bvec, xk, &
                      xhxa, xhxb, hppija, hppijb, noca, nocb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: wmo
    real(kind=dp), intent(in) :: target_energy
    real(kind=dp), intent(in), dimension(:) :: mo_energy_a
    real(kind=dp), intent(in), dimension(:,:) :: fa, fb
    real(kind=dp), intent(in), dimension(:) :: bvec
    real(kind=dp), intent(in), dimension(:) :: xk
    real(kind=dp), intent(in), dimension(:,:) :: xhxa, xhxb
    real(kind=dp), intent(in), dimension(:,:) :: hppija, hppijb
    integer, intent(in) :: noca, nocb

    real(kind=dp), allocatable, dimension(:,:) :: wrk, wrk1, wrk2
    integer :: i, a, k, x, y, j, b, ij, nbf, nvirb, lr1, lr2

    nbf = ubound(fa, 1)
    lr1 = nocb+1
    lr2 = noca
    nvirb = nbf-nocb

    allocate(wrk(nbf,nbf), &
             wrk1(nbf,nbf), &
             wrk2(nbf,nbf), &
             source=0.0_dp)

!   ----- COPY xk -----
    ij = 0
    do i = nocb+1, noca
      do j = 1, nocb
        ij = ij+1
        wrk1(j,i) = xk(ij)
      end do
    end do

    do i = noca+1, nbf
      do j = 1, nocb
        ij = ij+1
        wrk1(j,i) = xk(ij)
      end do
    end do

    do k = noca+1, nbf
      do i = nocb+1, noca
        ij = ij+1
        wrk1(i,k) = xk(ij)
      end do
    end do

! ! W_ix
    wrk = 0.0_dp
    do x = 1, nocb
      do k = 1, nocb
        wrk(x,1) = wrk(x,1)-wrk1(k,lr1)*fa(k,x)
        wrk(x,2) = wrk(x,2)-wrk1(k,lr2)*fa(k,x)
      end do
    end do

    do x = 1, nocb
      do k = 1, nbf-noca
        wrk(x,1) = wrk(x,1)+wrk1(lr1,noca+k)*fa(noca+k,x)
        wrk(x,2) = wrk(x,2)+wrk1(lr2,noca+k)*fa(noca+k,x)
      end do
    end do

    wmo(1:nocb,lr1:lr2) = wrk(1:nocb,1:2)*0.5_dp &
                        + xhxa(1:nocb,lr1:lr2) &
                        + xhxb(1:nocb,lr1:lr2) &
                        + hppija(1:nocb,lr1:lr2)
    wmo(1:nocb,lr1) = wmo(1:nocb,lr1) &
                    + mo_energy_a(1:nocb)*wrk1(1:nocb,lr1)
    wmo(1:nocb,lr2) = wmo(1:nocb,lr2) &
                    + mo_energy_a(1:nocb)*wrk1(1:nocb,lr2)

!   ----- W_IA -----
    wrk = 0.0_dp
    do i = 1, nocb
      do a = 1, nbf-noca
        wrk(i,a) = wrk(i,a)+fa(lr1,i)*wrk1(lr1,noca+a)
        wrk(i,a) = wrk(i,a)+fa(lr2,i)*wrk1(lr2,noca+a)
      end do
    end do

    wmo(1:nocb,noca+1:nbf) = wrk(1:nocb,1:nbf-noca)*0.5_dp &
                          + xhxb(1:nocb,noca+1:nbf)

    do a = noca+1, nbf
      wmo(1:nocb,a) = wmo(1:nocb,a) &
                    + mo_energy_a(1:nocb)*wrk1(1:nocb,a)
    end do

!   ----- W_XA -----
    wrk = 0.0_dp
    do a = 1, nbf-noca
      do k = 1, nocb
        wrk(1,a) = wrk(1,a)+fa(k,lr1)*wrk1(k,noca+a)
        wrk(2,a) = wrk(2,a)+fa(k,lr2)*wrk1(k,noca+a)
      end do
    end do

    do a = 1, nbf-noca
      wrk(1,a) = wrk(1,a)-fb(lr1,lr1)*wrk1(lr1,noca+a) &
                         -fb(lr2,lr1)*wrk1(lr2,noca+a)
      wrk(2,a) = wrk(2,a)-fb(lr1,lr2)*wrk1(lr1,noca+a) &
                         -fb(lr2,lr2)*wrk1(lr2,noca+a)
    end do

    wmo(lr1:lr2,noca+1:nbf) = wrk(1:2,1:nbf-noca)*0.5_dp &
                           + xhxb(lr1:lr2,noca+1:nbf)

    do a = noca+1, nbf
      wmo(lr1:lr2,a) = wmo(lr1:lr2,a) &
                     + mo_energy_a(lr1:lr2)*wrk1(lr1:lr2,a)
    end do

!   Alpha intermediate
    wrk = - fb
    do i = nocb+1, nbf
        wrk(i,i) = wrk(i,i)+target_energy
    end do

    wrk1(1:nbf-nocb,1:nbf-nocb) = wrk(nocb+1:nbf,nocb+1:nbf)*2.0_dp

    call dgemm('n', 'n', noca, nvirb, nvirb, &
               1.0_dp, bvec, noca, &
                       wrk1, nbf, &
               0.0_dp, wrk2, noca)
    call dgemm('n', 't', noca, noca, nvirb, &
               1.0_dp, wrk2, noca, &
                        bvec, noca, &
               0.0_dp, wrk1, nbf)

!   beta intermediate
    wrk = fa
    do i = 1, noca
        wrk(i,i) = wrk(i,i)+target_energy
    end do
    wrk(1:noca,1:noca) = wrk(1:noca,1:noca)*2.0_dp

    call dgemm('n', 'n', noca, nvirb, noca,  &
               1.0_dp, wrk, nbf, &
                       bvec, noca, &
               0.0_dp, wrk2, noca)
    call dgemm('t', 'n', nvirb, nvirb, noca, &
               1.0_dp, bvec, noca, &
                       wrk2, noca, &
               0.0_dp, wrk, nbf)

  ! W_ij
    do i = 1, nocb
      do j = 1, i
        wmo(i,j) = hppija(i,j)+hppijb(i,j)+wrk1(i,j)
      end do
    end do

  ! W_xy
    do x = nocb+1, noca
      do y = nocb+1, x
        wmo(x,y) = hppija(x,y)+wrk1(x,y)+wrk(x-nocb,y-nocb)
      end do
    end do

  ! W_ab
    do a = noca+1, nbf
      do b = noca+1, a
        wmo(a,b) = wrk(a-nocb,b-nocb)
      end do
    end do

  ! Scale diagonal elements
    do i = 1, nbf
      wmo(i,i) = wmo(i,i)*0.5_dp
    end do

    wmo = -wmo
    deallocate(wrk, wrk1, wrk2)

  end subroutine sfrowcal

  function get_spin_square(dmat_a,dmat_b,ta,tb,abxc,Smat,nocb) result(s2)
  ! dmat_a / dmat_b -- alpha/beta density of the excited state
  ! ta / tb -- alpha/beta difference density matrix
    use precision, only : dp
    use mathlib, only: symmetrize_matrix, traceprod_sym_packed
    use mathlib, only: pack_matrix, unpack_matrix
    use messages, only: show_message, with_abort

    implicit none

    real(kind=dp), intent(in), dimension(:) :: &
      dmat_a, dmat_b, ta, tb
    real(kind=dp), intent(in), dimension(:,:) :: abxc
    real(kind=dp), intent(in), dimension(:) :: smat
    integer, intent(in) :: nocb
    real(kind=dp) :: s2

    real(kind=dp), allocatable :: scr1(:), dmat_t(:), &
      dmat_t_sq(:,:), smat_sq(:,:), tmp1(:,:), tmp2(:,:)
    integer :: nbf, nbf_tri, ok
    real(kind=dp) :: dum1, dum2, dum3, dum4

    nbf = ubound(abxc, 1)
    nbf_tri = ubound(dmat_a, 1)

    allocate(scr1(nbf_tri), &
             dmat_t(nbf_tri), &
             dmat_t_sq(nbf,nbf), &
             smat_sq(nbf,nbf), &
             tmp1(nbf,nbf), &
             tmp2(nbf,nbf), &
             source=0.0_dp, stat=ok)
    if (ok/=0) call show_message('Cannot allocate memory in qet_spin_square',with_abort)

   ! Calculate spin expectation values
     dum1 = nocb+1

   ! Symmetric matrix scr1 = Smat*Dmat_a*Smat
     dmat_t = dmat_a + ta
     call unpack_matrix(dmat_t,dmat_t_sq)
     call unpack_matrix(smat,smat_sq)
     call dgemm('n', 'n', nbf, nbf, nbf, &
                1.0_dp, smat_sq, nbf, &
                        dmat_t_sq, nbf, &
                0.0_dp, tmp1, nbf)
     call dgemm('n', 'n', nbf, nbf, nbf, &
                1.0_dp, tmp1, nbf, &
                        smat_sq, nbf, &
                0.0_dp, tmp2, nbf)
     call pack_matrix(tmp2,scr1)
   ! -tr[ Dmat_b*Smat*Dmat_a*Smat ]
     dmat_t = dmat_b + tb
     dum2 = -traceprod_sym_packed(dmat_t,scr1,nbf)

   ! Symmetric matrix scr1 = Smat*Ta*Smat
     call unpack_matrix(Ta,tmp1)
     call dgemm('n', 'n', nbf, nbf, nbf, &
                1.0_dp, smat_sq, nbf, &
                        tmp1, nbf, &
                0.0_dp, tmp2, nbf)
     call dgemm('n', 'n', nbf, nbf, nbf, &
                1.0_dp, tmp2, nbf, &
                        smat_sq, nbf, &
                0.0_dp, tmp1, nbf)
     call pack_matrix(tmp1,scr1)
   ! -tr[ Tb*Smat*Ta*Smat ])
     dum3 =-traceprod_sym_packed(tb,scr1,nbf)

   ! +tr[ abxc*Smat ]
     tmp1 = abxc
     call symmetrize_matrix(tmp1, nbf)
     call pack_matrix(tmp1, scr1)
     dum4 = traceprod_sym_packed(scr1, smat, nbf)/2.0_dp

     s2 = dum1 + dum2 - dum3 + dum4**2

 end function get_spin_square

  subroutine get_transition_density(trden, bvec_mo, nbf, nocca, noccb, &
                                    nstates)
  ! compute transition density between ground state and excited states
    use precision, only: dp
    use tdhf_lib, only: iatogen
    use messages, only: show_message, with_abort

    implicit none

    real(kind=dp), intent(out), dimension(:,:,:,:) :: trden
    real(kind=dp), intent(in), dimension(:,:) :: bvec_mo
    integer, intent(in) :: nbf, nocca, noccb, nstates

    real(kind=dp), allocatable :: tmp(:,:)
    integer :: jst, ok

    allocate(tmp(nbf,nbf), &
             source=0.0d0, stat=ok)

    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    ! Compute transition dipole between the ground state and all excited
    do jst = 1, nstates
      ! Compute transition density
      ! unpack X
      call iatogen(bvec_mo(:,jst), trden(:,:,1,jst), nocca, noccb)
    end do

  end subroutine get_transition_density

  subroutine get_transition_dipole(basis, dip, mo_a, trden, nstates)
    use precision, only: dp
    use int1
!   use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use mathlib, only: orthogonal_transform, symmetrize_matrix, traceprod_sym_packed
    use mathlib, only: pack_matrix, unpack_matrix

    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: trden(:,:,:,:), mo_a(:,:)
    real(kind=dp), intent(out) :: dip(:,:,:)
    integer, intent(in) :: nstates

    real(kind=dp) :: center_of_mass(3)
    real(kind=dp), allocatable :: mints(:,:), trden_ao(:,:)
    real(kind=dp), allocatable, target :: tmp(:,:)
    real(kind=dp), pointer :: tmp2(:)
    integer :: nbf, nbf2, ok
    integer :: ist, jst

    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    allocate(mints(nbf2,3), &
             trden_ao(nbf,nbf), &
             tmp(nbf,nbf), &
             source=0.0_dp, stat=ok)

    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    ! Compute dipole integrals at the center of mass
    center_of_mass = basis%atoms%center(weight='mass')

    call multipole_integrals(basis, mints, center_of_mass, 1)

    do ist = 1, nstates
      do jst = 1, nstates
        if (ist==jst) cycle

        ! Convert transition density from MO to AO basis
        call orthogonal_transform('t', nbf, mo_a, trden(:,:,ist,jst), trden_ao, tmp)

        tmp2(1:nbf2) => tmp
        call symmetrize_matrix(trden_ao, nbf)
        call pack_matrix(trden_ao, tmp2)

        ! Compute dipole moment:
        ! D_i = Tr(T * dipole_ints_i), i = x, y, z
        dip(1,ist,jst) = -traceprod_sym_packed(tmp2, mints(:,1), nbf)*0.5_dp
        dip(2,ist,jst) = -traceprod_sym_packed(tmp2, mints(:,2), nbf)*0.5_dp
        dip(3,ist,jst) = -traceprod_sym_packed(tmp2, mints(:,3), nbf)*0.5_dp
      end do
    end do

  end subroutine
end module tdhf_sf_lib
