module hf_rohf_orbital_response_mod

  use precision, only: dp

  implicit none

  private

  integer, parameter :: derivative_probe_block_size=64

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
    use mathlib, only: unpack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, &
      der_nucattr_matrix
    use fock_deriv_mod, only: fock_deriv_contract_os_batch
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use cphf_mod, only: cphf_solve_rohf, rohf_pack_trial, &
      rohf_unpack_trial
    use mrsf_xc_fock_total_derivative_mod, only: &
      mrsf_xc_fock_total_derivative

    type(information), target, intent(inout) :: infos
    type(rohf_nuclear_response_t), intent(inout) :: response
    integer, intent(out) :: status
    logical, intent(in), optional :: require_two_somo

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(int2_tdgrd_data_t) :: int2_data
    real(dp), contiguous, pointer :: dma(:),dmb(:),mo(:,:)
    real(dp), contiguous, pointer :: focka(:),fockb(:)
    real(dp), allocatable :: pa(:,:),pb(:,:),ptot(:,:)
    real(dp), allocatable :: dsa(:,:,:,:),dta(:,:,:,:),dva(:,:,:,:)
    real(dp), allocatable :: famo(:,:),fbmo(:,:),scr(:,:),tmp(:,:)
    real(dp), allocatable :: sxmo(:,:),hxmo(:,:)
    real(dp), allocatable :: ga2e(:,:,:),gb2e(:,:,:)
    real(dp), allocatable :: d0a(:,:),d0b(:,:)
    real(dp), allocatable :: d0a_all(:,:,:),d0b_all(:,:,:), &
      dvxc_a(:,:,:),dvxc_b(:,:,:),sxmo_all(:,:,:),hxmo_all(:,:,:)
    real(dp), allocatable, target :: response_density(:,:,:)
    real(dp), allocatable, target :: probe_batch(:,:,:), &
      pcoul_batch(:,:,:),pexch_batch(:,:,:)
    real(dp), allocatable :: gx_batch(:,:,:)
    integer, allocatable :: probe_spin(:),probe_virtual(:),probe_occupied(:)
    real(dp), allocatable :: gfull(:,:),gd0(:,:),ba(:,:),bb(:,:)
    real(dp), allocatable :: bvec(:,:),uvec(:,:),xa(:,:),xb(:,:)
    real(dp), allocatable :: connection(:,:),connection_alpha(:,:), &
      connection_common(:,:)
    real(dp), allocatable :: dvecp(:,:,:,:)
    real(dp) :: hfscale
    integer :: nbf,natom,ncart,nocca,noccb,nvira,nvirb
    integer :: offset,ltot,alloc_status
    integer :: a,i,j,mu,nu,kc,cc,x,nprobe,probe_index,first,nactive,slot
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
             sxmo(nbf,nbf),hxmo(nbf,nbf), &
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
    nprobe=nvira*nocca+nvirb*noccb
    allocate(probe_batch(nbf,nbf,derivative_probe_block_size), &
      pcoul_batch(nbf,nbf,derivative_probe_block_size), &
      pexch_batch(nbf,nbf,derivative_probe_block_size), &
      gx_batch(3,natom,derivative_probe_block_size), &
      probe_spin(nprobe),probe_virtual(nprobe),probe_occupied(nprobe), &
      stat=alloc_status)
    if(alloc_status/=0) then
      call response%clean()
      status=-5
      return
    end if
    probe_index=0
    do a=1,nvira
      do i=1,nocca
        probe_index=probe_index+1
        probe_spin(probe_index)=1
        probe_virtual(probe_index)=a
        probe_occupied(probe_index)=i
      end do
    end do
    do a=1,nvirb
      do i=1,noccb
        probe_index=probe_index+1
        probe_spin(probe_index)=2
        probe_virtual(probe_index)=a
        probe_occupied(probe_index)=i
      end do
    end do
    do first=1,nprobe,derivative_probe_block_size
      nactive=min(derivative_probe_block_size,nprobe-first+1)
      probe_batch=0.0_dp
      do slot=1,nactive
        probe_index=first+slot-1
        a=probe_virtual(probe_index)
        i=probe_occupied(probe_index)
        pcoul_batch(:,:,slot)=ptot
        if(probe_spin(probe_index)==1) then
          pexch_batch(:,:,slot)=pa
          do nu=1,nbf
            do mu=1,nbf
              probe_batch(mu,nu,slot)=0.5_dp*( &
                mo(mu,nocca+a)*mo(nu,i)+mo(mu,i)*mo(nu,nocca+a))
            end do
          end do
        else
          pexch_batch(:,:,slot)=pb
          do nu=1,nbf
            do mu=1,nbf
              probe_batch(mu,nu,slot)=0.5_dp*( &
                mo(mu,noccb+a)*mo(nu,i)+mo(mu,i)*mo(nu,noccb+a))
            end do
          end do
        end if
      end do
      call fock_deriv_contract_os_batch(infos,basis, &
        pcoul_batch(:,:,1:nactive),pexch_batch(:,:,1:nactive), &
        probe_batch(:,:,1:nactive),hfscale,gx_batch(:,:,1:nactive))
      do slot=1,nactive
        probe_index=first+slot-1
        a=probe_virtual(probe_index)
        i=probe_occupied(probe_index)
        if(probe_spin(probe_index)==1) then
          ga2e(a,i,:)=reshape(gx_batch(:,:,slot),[ncart])
        else
          gb2e(a,i,:)=reshape(gx_batch(:,:,slot),[ncart])
        end if
      end do
    end do

    allocate(d0a(nbf,nbf),d0b(nbf,nbf), &
             gfull(nbf,nbf),gd0(nbf,nbf),ba(nvira,nocca), &
             bb(nvirb,noccb),bvec(ltot,ncart),uvec(ltot,ncart), &
             xa(nvira,nocca),xb(nvirb,noccb),connection(nbf,nbf), &
             connection_alpha(nbf,nbf),connection_common(nbf,nbf), &
             d0a_all(nbf,nbf,ncart),d0b_all(nbf,nbf,ncart), &
             dvxc_a(nbf,nbf,ncart),dvxc_b(nbf,nbf,ncart), &
             sxmo_all(nbf,nbf,ncart),hxmo_all(nbf,nbf,ncart), &
             response_density(nbf,nbf,2*ncart), &
             source=0.0_dp,stat=alloc_status)
    if(alloc_status/=0) then
      call response%clean()
      status=-5
      return
    end if

    ! Form the overlap-metric density response for every perturbation first.
    ! The semilocal ROKS contribution is then one exact analytic all-coordinate
    ! grid contraction, replacing the historical +/- displaced-grid loop.
    do x=1,ncart
      call transform_ao_to_mo(mo,response%ds_ao(:,:,x),scr,tmp,sxmo)
      call transform_ao_to_mo(mo,response%dhcore_ao(:,:,x),scr,tmp,hxmo)
      sxmo_all(:,:,x)=sxmo
      hxmo_all(:,:,x)=hxmo

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
      d0a_all(:,:,x)=d0a
      d0b_all(:,:,x)=d0b
      response_density(:,:,2*x-1)=d0a
      response_density(:,:,2*x)=d0b
    end do

    ! Apply the exact open-shell Coulomb/exchange operator to every metric
    ! density derivative during one shell-quartet traversal.  This is the
    ! multi-right-hand-side analogue of the former coordinate-wise fock_jk
    ! loop; no density fitting, integral approximation, or changed screening
    ! threshold is introduced.
    call int2_driver%init(basis,infos)
    call int2_driver%set_screening()
    int2_data=int2_tdgrd_data_t(d2=response_density,int_apb=.true., &
      int_amb=.false.,tamm_dancoff=.false.,scale_exchange=hfscale)
    call int2_driver%run(int2_data)
    if(dft) then
      block
        use mod_dft, only: dft_initialize,dftclean
        use mod_dft_molgrid, only: dft_grid_t
        type(dft_grid_t) :: grid
        call dft_initialize(infos,basis,grid)
        call mrsf_xc_fock_total_derivative(basis,grid,pa,pb,d0a_all, &
          d0b_all,dvxc_a,dvxc_b,infos,status,threshold=0.0_dp)
        call dftclean(infos)
      end block
      if(status/=0) then
        call response%clean()
        return
      end if
    end if

    do x=1,ncart
      sxmo=sxmo_all(:,:,x)
      hxmo=hxmo_all(:,:,x)
      d0a=d0a_all(:,:,x)
      d0b=d0b_all(:,:,x)
      ! int2_tdgrd_data_t acts on D+D^T.  The metric density derivatives
      ! above are symmetric, so divide the returned A+B image by two to
      ! recover the ordinary open-shell Fock response G[D^x].
      gfull=0.5_dp*int2_data%apb(:,:,2*x-1,1)
      call transform_ao_to_mo(mo,gfull,scr,tmp,gd0)
      do i=1,nocca
        do a=1,nvira
          ba(a,i)=-(hxmo(i,nocca+a)+ga2e(a,i,x)+ &
            gd0(i,nocca+a)) &
            +dot_product(sxmo(nocca+a,1:nocca),famo(1:nocca,i)) &
            +dot_product(famo(nocca+a,1:nocca),sxmo(1:nocca,i))
        end do
      end do
      gfull=0.5_dp*int2_data%apb(:,:,2*x,1)
      call transform_ao_to_mo(mo,gfull,scr,tmp,gd0)
      do i=1,noccb
        do a=1,nvirb
          bb(a,i)=-(hxmo(i,noccb+a)+gb2e(a,i,x)+ &
            gd0(i,noccb+a)) &
            +dot_product(sxmo(noccb+a,1:noccb),fbmo(1:noccb,i)) &
            +dot_product(fbmo(noccb+a,1:noccb),sxmo(1:noccb,i))
        end do
      end do

      if(dft) then
        call transform_ao_to_mo(mo,dvxc_a(:,:,x),scr,tmp,gd0)
        do i=1,nocca
          do a=1,nvira
            ba(a,i)=ba(a,i)-gd0(i,nocca+a)
          end do
        end do
        call transform_ao_to_mo(mo,dvxc_b(:,:,x),scr,tmp,gd0)
        do i=1,noccb
          do a=1,nvirb
            bb(a,i)=bb(a,i)-gd0(i,noccb+a)
          end do
        end do
      end if
      call rohf_pack_trial(bvec(:,x),ba,bb,nbf,nocca,noccb)
    end do
    call int2_data%clean()
    call int2_driver%clean()

    ! Hessian rows differentiate this response once more.  The ordinary
    ! gradient tolerance leaves a few microhartree/bohr**2 of antisymmetric
    ! noise in the unsymmetrized matrix, so solve the nuclear response to the
    ! tighter accuracy required by the second derivative.
    call cphf_solve_rohf(infos,ncart,bvec,uvec,tol=1.0d-13,status=status)
    if(status/=0) then
      call response%clean()
      return
    end if

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
    if(nbf<=0 .or. nocc<0 .or. nvir<=0 .or. &
       any(shape(mo)/=[nbf,nbf]) .or. &
       any(shape(sx_mo)/=[nbf,nbf]) .or. &
       any(shape(xvo)/=[nvir,nocc]) .or. &
       any(shape(dmo)/=[nbf,nbf]) .or. &
       any(shape(connection)/=[nbf,nbf])) then
      status=-6
      return
    end if
    connection=-0.5_dp*sx_mo
    if(nocc>0) then
      connection(nocc+1:nbf,1:nocc)=xvo
      connection(1:nocc,nocc+1:nbf)= &
        -sx_mo(1:nocc,nocc+1:nbf)-transpose(xvo)
    end if
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
