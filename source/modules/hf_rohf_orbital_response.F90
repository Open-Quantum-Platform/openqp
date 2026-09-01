module hf_rohf_orbital_response_mod

  use precision, only: dp

  implicit none

  private

  type, public :: rohf_nuclear_response_t
    integer :: nbf = 0
    integer :: nocca = 0
    integer :: noccb = 0
    integer :: ncart = 0
    real(dp), allocatable :: dmo_alpha(:,:,:)
    real(dp), allocatable :: dmo_beta(:,:,:)
    real(dp), allocatable :: dmo_common(:,:,:)
    real(dp), allocatable :: ds_ao(:,:,:)
    real(dp), allocatable :: dhcore_ao(:,:,:)
  contains
    procedure :: clean => rohf_nuclear_response_clean
  end type rohf_nuclear_response_t

  public :: build_rohf_nuclear_response
  public :: complete_rohf_orbital_connection
  public :: rohf_response_status_message

contains

!###############################################################################

  subroutine build_rohf_nuclear_response(infos,response,status, &
      require_two_somo)
    ! Construct the ROHF/ROKS nuclear orbital response for all 3N Cartesian
    ! perturbations.  The CPHF unknowns retain the common ROHF docc/SOMO/virtual
    ! rotation space, while the returned alpha and beta orbital connections use
    ! their distinct occupied/virtual boundaries.
    !
    ! Method-development provenance: the analytical TDHF Hessian implementation
    ! initiated by Hiroya Nakata is the methodological starting point for this
    ! reusable OpenQP response construction.  The equations below preserve the
    ! established OpenQP ROHF/ROKS CPHF right-hand side exactly.
    !
    ! The returned dS and dHcore are explicit AO integral derivatives.  A total
    ! dF/dR is deliberately not returned: its production construction additionally
    ! requires a full open-shell derivative-ERI matrix and, for ROKS, the moving-
    ! grid XC contribution.  The available scalar-probe contraction would require
    ! O(nbf**2) complete derivative-integral evaluations per spin and coordinate.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_DM_B, &
      OQP_VEC_MO_A, OQP_FOCK_A, OQP_FOCK_B
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, &
      der_nucattr_matrix
    use fock_deriv_mod, only: fock_deriv_contract_os
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve_rohf, rohf_pack_trial, &
      rohf_unpack_trial

    type(information), target, intent(inout) :: infos
    type(rohf_nuclear_response_t), intent(inout) :: response
    integer, intent(out) :: status
    logical, intent(in), optional :: require_two_somo

    type(basis_set), pointer :: basis
    real(dp), contiguous, pointer :: dma(:),dmb(:),mo(:,:)
    real(dp), contiguous, pointer :: focka(:),fockb(:)
    real(dp), allocatable :: pa(:,:),pb(:,:),ptot(:,:)
    real(dp), allocatable :: dsa(:,:,:,:),dta(:,:,:,:),dva(:,:,:,:)
    real(dp), allocatable :: famo(:,:),fbmo(:,:),scr(:,:),tmp(:,:)
    real(dp), allocatable :: sxmo(:,:),hxmo(:,:),probe(:,:)
    real(dp), allocatable :: ga2e(:,:,:),gb2e(:,:,:)
    real(dp), allocatable :: d0a(:,:),d0b(:,:),dpck(:,:),fpck(:,:)
    real(dp), allocatable :: gfull(:,:),gd0(:,:),ba(:,:),bb(:,:)
    real(dp), allocatable :: bvec(:,:),uvec(:,:),xa(:,:),xb(:,:)
    real(dp), allocatable :: connection(:,:),connection_alpha(:,:), &
      connection_common(:,:)
    real(dp), allocatable :: dvecp(:,:,:,:)
    real(dp) :: hfscale,gx(3,size(infos%atoms%xyz,2))
    integer :: nbf,nbf2,natom,ncart,nocca,noccb,nvira,nvirb
    integer :: offset,ltot,alloc_status
    integer :: a,i,j,mu,nu,kc,cc,x
    logical :: two_somo_only,dft

    call response%clean()
    status=0
    two_somo_only=.true.
    if(present(require_two_somo)) two_somo_only=require_two_somo

    basis=>infos%basis
    basis%atoms=>infos%atoms
    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    ncart=3*natom
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    nvira=nbf-nocca
    nvirb=nbf-noccb
    offset=nocca-noccb
    ltot=noccb*(offset+nvira)+offset*nvira
    nbf2=nbf*(nbf+1)/2
    dft=infos%control%hamilton==20

    if(infos%control%scftype/=3 .or. offset<=0) then
      status=-1
      return
    end if
    if(infos%control%hamilton>20) then
      status=-7
      return
    end if
    if(nbf<=0 .or. natom<=0 .or. noccb<0 .or. nocca>nbf .or. &
       nvira<=0 .or. nvirb<=0 .or. ltot<=0) then
      status=-2
      return
    end if
    if(dft .and. infos%dft%cam_flag) then
      status=-3
      return
    end if
    if(two_somo_only .and. offset/=2) then
      status=-4
      return
    end if

    hfscale=1.0_dp
    if(dft) hfscale=infos%dft%hfscale

    call tagarray_get_data(infos%dat,OQP_DM_A,dma)
    call tagarray_get_data(infos%dat,OQP_DM_B,dmb)
    call tagarray_get_data(infos%dat,OQP_VEC_MO_A,mo)
    call tagarray_get_data(infos%dat,OQP_FOCK_A,focka)
    call tagarray_get_data(infos%dat,OQP_FOCK_B,fockb)

    response%nbf=nbf
    response%nocca=nocca
    response%noccb=noccb
    response%ncart=ncart
    allocate(response%dmo_alpha(nbf,nbf,ncart), &
             response%dmo_beta(nbf,nbf,ncart), &
             response%dmo_common(nbf,nbf,ncart), &
             response%ds_ao(nbf,nbf,ncart), &
             response%dhcore_ao(nbf,nbf,ncart),source=0.0_dp, &
             stat=alloc_status)
    if(alloc_status/=0) then
      call response%clean()
      status=-5
      return
    end if

    allocate(pa(nbf,nbf),pb(nbf,nbf),ptot(nbf,nbf), &
             dsa(nbf,nbf,3,natom),dta(nbf,nbf,3,natom), &
             dva(nbf,nbf,3,natom),dvecp(nbf,nbf,3,natom), &
             famo(nbf,nbf),fbmo(nbf,nbf),scr(nbf,nbf),tmp(nbf,nbf), &
             sxmo(nbf,nbf),hxmo(nbf,nbf),probe(nbf,nbf), &
             stat=alloc_status)
    if(alloc_status/=0) then
      call response%clean()
      status=-5
      return
    end if
    call unpack_matrix(dma,pa)
    call unpack_matrix(dmb,pb)
    ptot=pa+pb

    call der_overlap_matrix(basis,dsa)
    call der_kinetic_matrix(basis,dta)
    call der_nucattr_matrix(basis,basis%atoms%xyz, &
      basis%atoms%zn-basis%ecp_zn_num,dva)
    do kc=1,natom
      do cc=1,3
        do nu=1,nbf
          do mu=1,nbf
            dsa(mu,nu,cc,kc)=dsa(mu,nu,cc,kc)* &
              basis%bfnrm(mu)*basis%bfnrm(nu)
            dta(mu,nu,cc,kc)=dta(mu,nu,cc,kc)* &
              basis%bfnrm(mu)*basis%bfnrm(nu)
            dva(mu,nu,cc,kc)=dva(mu,nu,cc,kc)* &
              basis%bfnrm(mu)*basis%bfnrm(nu)
          end do
        end do
      end do
    end do
    block
      use ecp_tool, only: ecp_deriv_ints
      call ecp_deriv_ints(basis,basis%atoms%xyz,dvecp)
    end block
    dva=dva+dvecp
    do x=1,ncart
      cc=mod(x-1,3)+1
      kc=(x-1)/3+1
      response%ds_ao(:,:,x)=dsa(:,:,cc,kc)
      response%dhcore_ao(:,:,x)=dta(:,:,cc,kc)+dva(:,:,cc,kc)
    end do

    call unpack_matrix(focka,scr)
    call transform_ao_to_mo(mo,scr,tmp,hxmo,famo)
    call unpack_matrix(fockb,scr)
    call transform_ao_to_mo(mo,scr,tmp,hxmo,fbmo)

    allocate(ga2e(nvira,nocca,ncart),gb2e(nvirb,noccb,ncart), &
             source=0.0_dp,stat=alloc_status)
    if(alloc_status/=0) then
      call response%clean()
      status=-5
      return
    end if
    do a=1,nvira
      do i=1,nocca
        do nu=1,nbf
          do mu=1,nbf
            probe(mu,nu)=0.5_dp*(mo(mu,nocca+a)*mo(nu,i)+ &
              mo(mu,i)*mo(nu,nocca+a))
          end do
        end do
        gx=0.0_dp
        call fock_deriv_contract_os(infos,basis,ptot,pa,probe, &
          hfscale,gx)
        ga2e(a,i,:)=reshape(gx,[ncart])
      end do
    end do
    do a=1,nvirb
      do i=1,noccb
        do nu=1,nbf
          do mu=1,nbf
            probe(mu,nu)=0.5_dp*(mo(mu,noccb+a)*mo(nu,i)+ &
              mo(mu,i)*mo(nu,noccb+a))
          end do
        end do
        gx=0.0_dp
        call fock_deriv_contract_os(infos,basis,ptot,pb,probe, &
          hfscale,gx)
        gb2e(a,i,:)=reshape(gx,[ncart])
      end do
    end do

    allocate(d0a(nbf,nbf),d0b(nbf,nbf),dpck(nbf2,2),fpck(nbf2,2), &
             gfull(nbf,nbf),gd0(nbf,nbf),ba(nvira,nocca), &
             bb(nvirb,noccb),bvec(ltot,ncart),uvec(ltot,ncart), &
             xa(nvira,nocca),xb(nvirb,noccb),connection(nbf,nbf), &
             connection_alpha(nbf,nbf),connection_common(nbf,nbf), &
             source=0.0_dp,stat=alloc_status)
    if(alloc_status/=0) then
      call response%clean()
      status=-5
      return
    end if

    do x=1,ncart
      call transform_ao_to_mo(mo,response%ds_ao(:,:,x),scr,tmp,sxmo)
      call transform_ao_to_mo(mo,response%dhcore_ao(:,:,x),scr,tmp,hxmo)

      d0a=0.0_dp
      d0b=0.0_dp
      do i=1,nocca
        do j=1,nocca
          do nu=1,nbf
            do mu=1,nbf
              d0a(mu,nu)=d0a(mu,nu)-sxmo(i,j)*mo(mu,i)*mo(nu,j)
            end do
          end do
        end do
      end do
      do i=1,noccb
        do j=1,noccb
          do nu=1,nbf
            do mu=1,nbf
              d0b(mu,nu)=d0b(mu,nu)-sxmo(i,j)*mo(mu,i)*mo(nu,j)
            end do
          end do
        end do
      end do
      call pack_matrix(d0a,dpck(:,1))
      call pack_matrix(d0b,dpck(:,2))
      fpck=0.0_dp
      call fock_jk(basis,d=dpck,f=fpck,scale_exch=hfscale,infos=infos)

      call unpack_matrix(fpck(:,1),gfull)
      call transform_ao_to_mo(mo,gfull,scr,tmp,gd0)
      do i=1,nocca
        do a=1,nvira
          ba(a,i)=-(hxmo(i,nocca+a)+ga2e(a,i,x)+ &
            gd0(i,nocca+a)) &
            +dot_product(sxmo(nocca+a,1:nocca),famo(1:nocca,i)) &
            +dot_product(famo(nocca+a,1:nocca),sxmo(1:nocca,i))
        end do
      end do
      call unpack_matrix(fpck(:,2),gfull)
      call transform_ao_to_mo(mo,gfull,scr,tmp,gd0)
      do i=1,noccb
        do a=1,nvirb
          bb(a,i)=-(hxmo(i,noccb+a)+gb2e(a,i,x)+ &
            gd0(i,noccb+a)) &
            +dot_product(sxmo(noccb+a,1:noccb),fbmo(1:noccb,i)) &
            +dot_product(fbmo(noccb+a,1:noccb),sxmo(1:noccb,i))
        end do
      end do

      if(dft) call add_semilocal_roks_rhs(infos,basis,mo,sxmo,x,ba,bb, &
        scr,tmp,status)
      if(status/=0) then
        call response%clean()
        return
      end if
      call rohf_pack_trial(bvec(:,x),ba,bb,nbf,nocca,noccb)
    end do

    ! Hessian rows differentiate this response once more.  The ordinary
    ! gradient tolerance leaves a few microhartree/bohr**2 of antisymmetric
    ! noise in the unsymmetrized matrix, so solve the nuclear response to the
    ! tighter accuracy required by the second derivative.
    call cphf_solve_rohf(infos,ncart,bvec,uvec,tol=1.0d-13)

    if(dft) then
      block
        use mod_dft, only: dft_initialize,dftclean
        use mod_dft_molgrid, only: dft_grid_t
        type(dft_grid_t) :: grid
        call dft_initialize(infos,basis,grid)
        call dftclean(infos)
      end block
    end if

    do x=1,ncart
      call rohf_unpack_trial(uvec(:,x),xa,xb,nbf,nocca,noccb)
      call transform_ao_to_mo(mo,response%ds_ao(:,:,x),scr,tmp,sxmo)
      call complete_rohf_orbital_connection(mo,sxmo,nocca,xa, &
        response%dmo_alpha(:,:,x),connection,status)
      if(status/=0) then
        call response%clean()
        return
      end if
      connection_alpha=connection
      call complete_rohf_orbital_connection(mo,sxmo,noccb,xb, &
        response%dmo_beta(:,:,x),connection,status)
      if(status/=0) then
        call response%clean()
        return
      end if
      ! A restricted open-shell calculation stores one common spatial-MO
      ! coefficient matrix.  The alpha completion above supplies its O-V and
      ! C-V rotations, whereas the beta completion supplies C-O and the same
      ! C-V rotations.  Each spin-specific completion is sufficient for its
      ! occupied density, but neither alone is the derivative of the complete
      ! common ROHF orbital set required by an MRSF response operator.
      connection_common=connection_alpha
      connection_common(1:noccb,noccb+1:nocca)= &
        connection(1:noccb,noccb+1:nocca)
      connection_common(noccb+1:nocca,1:noccb)= &
        connection(noccb+1:nocca,1:noccb)
      response%dmo_common(:,:,x)=matmul(mo,connection_common)
    end do
  end subroutine build_rohf_nuclear_response

!###############################################################################

  subroutine add_semilocal_roks_rhs(infos,basis,mo,sxmo,x,ba,bb, &
      scr,tmp,status)
    ! Preserve the existing ROKS CPKS perturbation: a central difference of the
    ! spin XC Fock at R+/-h with the occupied orbitals reorthonormalized by dS.
    use types, only: information
    use basis_tools, only: basis_set
    use mathlib, only: unpack_matrix
    use mod_dft, only: dft_initialize,dftclean,dftexcor
    use mod_dft_molgrid, only: dft_grid_t

    type(information), target, intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(dp), intent(in) :: mo(:,:),sxmo(:,:)
    integer, intent(in) :: x
    real(dp), intent(inout) :: ba(:,:),bb(:,:)
    real(dp), intent(inout) :: scr(:,:),tmp(:,:)
    integer, intent(out) :: status

    type(dft_grid_t) :: grid
    real(dp), allocatable :: dmoa(:,:),dmob(:,:),mopa(:,:),mopb(:,:)
    real(dp), allocatable :: frap(:),frbp(:),fram(:),frbm(:)
    real(dp), allocatable :: dvx(:,:),hxc(:,:)
    real(dp) :: step,exr,telr,tknr
    integer :: nbf,nbf2,nocca,noccb,nvira,nvirb,natom
    integer :: cc,kc,i,j,a,alloc_status

    status=0
    nbf=size(mo,1)
    nbf2=nbf*(nbf+1)/2
    natom=size(basis%atoms%xyz,2)
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    nvira=nbf-nocca
    nvirb=nbf-noccb
    cc=mod(x-1,3)+1
    kc=(x-1)/3+1
    if(x<1 .or. x>3*natom .or. infos%dft%cam_flag) then
      status=-3
      return
    end if
    allocate(dmoa(nbf,nocca),dmob(nbf,noccb),mopa(nbf,nbf), &
             mopb(nbf,nbf),frap(nbf2),frbp(nbf2),fram(nbf2), &
             frbm(nbf2),dvx(nbf,nbf),hxc(nbf,nbf), &
             stat=alloc_status)
    if(alloc_status/=0) then
      status=-5
      return
    end if
    dmoa=0.0_dp
    do i=1,nocca
      do j=1,nocca
        dmoa(:,i)=dmoa(:,i)-0.5_dp*sxmo(j,i)*mo(:,j)
      end do
    end do
    dmob=0.0_dp
    do i=1,noccb
      do j=1,noccb
        dmob(:,i)=dmob(:,i)-0.5_dp*sxmo(j,i)*mo(:,j)
      end do
    end do

    step=1.0d-3
    basis%atoms%xyz(cc,kc)=basis%atoms%xyz(cc,kc)+step
    call basis%init_shell_centers()
    call dft_initialize(infos,basis,grid)
    mopa=mo
    mopb=mo
    mopa(:,1:nocca)=mo(:,1:nocca)+step*dmoa
    mopb(:,1:noccb)=mo(:,1:noccb)+step*dmob
    frap=0.0_dp
    frbp=0.0_dp
    call dftexcor(basis,grid,int(infos%control%scftype),frap,frbp, &
      mopa,mopb,nbf,nbf2,exr,telr,tknr,infos)
    call dftclean(infos)

    basis%atoms%xyz(cc,kc)=basis%atoms%xyz(cc,kc)-2.0_dp*step
    call basis%init_shell_centers()
    call dft_initialize(infos,basis,grid)
    mopa=mo
    mopb=mo
    mopa(:,1:nocca)=mo(:,1:nocca)-step*dmoa
    mopb(:,1:noccb)=mo(:,1:noccb)-step*dmob
    fram=0.0_dp
    frbm=0.0_dp
    call dftexcor(basis,grid,int(infos%control%scftype),fram,frbm, &
      mopa,mopb,nbf,nbf2,exr,telr,tknr,infos)
    call dftclean(infos)

    basis%atoms%xyz(cc,kc)=basis%atoms%xyz(cc,kc)+step
    call basis%init_shell_centers()
    call unpack_matrix((frap-fram)/(2.0_dp*step),dvx)
    call transform_ao_to_mo(mo,dvx,scr,tmp,hxc)
    do i=1,nocca
      do a=1,nvira
        ba(a,i)=ba(a,i)-hxc(i,nocca+a)
      end do
    end do
    call unpack_matrix((frbp-frbm)/(2.0_dp*step),dvx)
    call transform_ao_to_mo(mo,dvx,scr,tmp,hxc)
    do i=1,noccb
      do a=1,nvirb
        bb(a,i)=bb(a,i)-hxc(i,noccb+a)
      end do
    end do
  end subroutine add_semilocal_roks_rhs

!###############################################################################

  subroutine complete_rohf_orbital_connection(mo,sx_mo,nocc,xvo,dmo, &
      connection,status)
    ! Complete the occupied-virtual CPHF response in the symmetric metric gauge.
    ! U=C^T S dC obeys U+U^T=-S^x_MO.  The solved U_vo block is kept exactly;
    ! U_ov follows from the metric identity, while U_oo and U_vv are -S^x/2.
    real(dp), intent(in) :: mo(:,:),sx_mo(:,:),xvo(:,:)
    integer, intent(in) :: nocc
    real(dp), intent(out) :: dmo(:,:),connection(:,:)
    integer, intent(out) :: status

    integer :: nbf,nvir

    status=0
    nbf=size(mo,1)
    nvir=nbf-nocc
    if(nbf<=0 .or. nocc<=0 .or. nvir<=0 .or. &
       any(shape(mo)/=[nbf,nbf]) .or. &
       any(shape(sx_mo)/=[nbf,nbf]) .or. &
       any(shape(xvo)/=[nvir,nocc]) .or. &
       any(shape(dmo)/=[nbf,nbf]) .or. &
       any(shape(connection)/=[nbf,nbf])) then
      status=-6
      return
    end if
    connection=-0.5_dp*sx_mo
    connection(nocc+1:nbf,1:nocc)=xvo
    connection(1:nocc,nocc+1:nbf)= &
      -sx_mo(1:nocc,nocc+1:nbf)-transpose(xvo)
    dmo=matmul(mo,connection)
  end subroutine complete_rohf_orbital_connection

!###############################################################################

  pure function rohf_response_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=160) :: message

    select case(status)
    case(-1)
      message='ROHF nuclear response requires an open-shell restricted reference.'
    case(-2)
      message='ROHF nuclear response received inconsistent orbital or nuclear dimensions.'
    case(-3)
      message='ROHF nuclear response does not yet support CAM/range-separated CPKS.'
    case(-4)
      message='Two-SOMO ROHF nuclear response requires nelec_alpha-nelec_beta=2.'
    case(-5)
      message='Cannot allocate memory for the ROHF nuclear orbital response.'
    case(-6)
      message='Cannot complete the full ROHF orbital connection from the CPHF response.'
    case(-7)
      message='ROHF nuclear response does not support this electronic Hamiltonian.'
    case default
      message='ROHF nuclear orbital response failed.'
    end select
  end function rohf_response_status_message

!###############################################################################

  subroutine rohf_nuclear_response_clean(this)
    class(rohf_nuclear_response_t), intent(inout) :: this

    if(allocated(this%dmo_alpha)) deallocate(this%dmo_alpha)
    if(allocated(this%dmo_beta)) deallocate(this%dmo_beta)
    if(allocated(this%dmo_common)) deallocate(this%dmo_common)
    if(allocated(this%ds_ao)) deallocate(this%ds_ao)
    if(allocated(this%dhcore_ao)) deallocate(this%dhcore_ao)
    this%nbf=0
    this%nocca=0
    this%noccb=0
    this%ncart=0
  end subroutine rohf_nuclear_response_clean

!###############################################################################

  subroutine transform_ao_to_mo(mo,ao,scr,tmp,mo_matrix)
    real(dp), intent(in) :: mo(:,:),ao(:,:)
    real(dp), intent(inout) :: scr(:,:),tmp(:,:),mo_matrix(:,:)
    integer :: nbf

    nbf=size(mo,1)
    call dgemm('t','n',nbf,nbf,nbf,1.0_dp,mo,nbf,ao,nbf, &
      0.0_dp,scr,nbf)
    call dgemm('n','n',nbf,nbf,nbf,1.0_dp,scr,nbf,mo,nbf, &
      0.0_dp,mo_matrix,nbf)
  end subroutine transform_ao_to_mo

end module hf_rohf_orbital_response_mod
