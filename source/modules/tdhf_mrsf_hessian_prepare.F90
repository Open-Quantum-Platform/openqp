module tdhf_mrsf_hessian_prepare_mod

  use precision, only: dp

  implicit none

  private

  type, public :: mrsf_hessian_first_response_t
    integer :: nbf=0
    integer :: ncoord=0
    integer :: packed=0
    real(kind=dp) :: omega=0.0_dp
    real(kind=dp) :: residual=0.0_dp
    real(kind=dp), allocatable :: mo_a(:,:),mo_b(:,:),fock_a(:,:),fock_b(:,:)
    real(kind=dp), allocatable :: x(:),dmo_a(:,:,:),dmo_b(:,:,:), &
      dmo_common(:,:,:),ds(:,:,:)
    real(kind=dp), allocatable :: dhcore(:,:,:),dpa(:,:,:),dpb(:,:,:)
    real(kind=dp), allocatable :: dfock_a(:,:,:),dfock_b(:,:,:)
    real(kind=dp), allocatable :: dax(:,:),dx(:,:),domega(:)
  contains
    procedure :: clean => mrsf_hessian_first_response_clean
  end type mrsf_hessian_first_response_t

  public :: prepare_mrsf_hessian_first_response

contains

!###############################################################################

  subroutine prepare_mrsf_hessian_first_response(infos,response,status)
    ! Build the complete first nuclear response needed before differentiating
    ! the stationary MRSF gradient.  Hiroya Nakata's analytical TDHF/TDDFT
    ! Hessian is the methodological starting point.  The state response here
    ! remains entirely in the physical two-SOMO spin-adapted representation.

    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use oqp_tagarray_driver, only: tagarray_get_data,data_has_tags, &
      OQP_VEC_MO_A,OQP_VEC_MO_B,OQP_FOCK_A,OQP_FOCK_B,OQP_td_bvec_mo, &
      OQP_td_energies
    use mathlib, only: unpack_matrix
    use hf_rohf_orbital_response_mod, only: rohf_nuclear_response_t, &
      build_rohf_nuclear_response
    use tdhf_mrsf_hessian_fock_response_mod, only: &
      build_orbital_density_derivatives
    use tdhf_mrsf_hessian_state_response_mod, only: &
      solve_mrsf_first_nuclear_response
    use mrsf_xc_fock_total_derivative_mod, only: &
      mrsf_xc_fock_total_derivative
    use messages, only: WITH_ABORT
    use io_constants, only: iw

    type(information), target, intent(inout) :: infos
    type(mrsf_hessian_first_response_t), intent(inout) :: response
    integer, intent(out) :: status

    character(len=*), parameter :: required(*)=[character(len=80) :: &
      OQP_VEC_MO_A,OQP_VEC_MO_B,OQP_FOCK_A,OQP_FOCK_B,OQP_td_bvec_mo, &
      OQP_td_energies]
    type(basis_set), pointer :: basis
    type(rohf_nuclear_response_t) :: orbital_response
    type(int2_compute_t) :: int2_driver
    real(kind=dp), contiguous, pointer :: stored_mo_a(:,:),stored_mo_b(:,:), &
      packed_fock_a(:),packed_fock_b(:),state_vectors(:,:),energies(:)
    real(kind=dp), allocatable :: pa(:,:),pb(:,:),dvxc_a(:,:,:),dvxc_b(:,:,:)
    integer :: nbf,nocca,noccb,ncoord,packed,target,local_status, &
      clock_start,clock_now,clock_rate
    real(kind=dp) :: time_orbital,time_xc_derivative,time_state
    logical :: dft

    call response%clean()
    basis=>infos%basis
    basis%atoms=>infos%atoms
    nbf=basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    ncoord=3*size(infos%atoms%xyz,2)
    packed=nocca*(nbf-noccb)
    target=min(infos%tddft%target_state,infos%tddft%nstate)
    dft=infos%control%hamilton==20
    status=0
    if(nbf<=0 .or. ncoord<=0 .or. packed<=0 .or. target<=0 .or. &
       infos%control%scftype/=3 .or. nocca-noccb/=2 .or. &
       infos%tddft%umrsf .or. &
       (infos%tddft%mult/=1 .and. infos%tddft%mult/=3)) then
      status=-1
      return
    end if
    if(infos%control%hamilton>20 .or. &
       (dft .and. (infos%dft%cam_flag .or. infos%functional%needTau .or. &
        infos%functional%needLapl))) then
      status=-2
      return
    end if

    call data_has_tags(infos%dat,required,'tdhf_mrsf_hessian_prepare_mod', &
      'prepare_mrsf_hessian_first_response',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_VEC_MO_A,stored_mo_a)
    call tagarray_get_data(infos%dat,OQP_VEC_MO_B,stored_mo_b)
    call tagarray_get_data(infos%dat,OQP_FOCK_A,packed_fock_a)
    call tagarray_get_data(infos%dat,OQP_FOCK_B,packed_fock_b)
    call tagarray_get_data(infos%dat,OQP_td_bvec_mo,state_vectors)
    call tagarray_get_data(infos%dat,OQP_td_energies,energies)
    if(any(shape(stored_mo_a)/=[nbf,nbf]) .or. &
       any(shape(stored_mo_b)/=[nbf,nbf]) .or. &
       size(state_vectors,1)/=packed .or. target>size(state_vectors,2) .or. &
       target>size(energies)) then
      status=-3
      return
    end if

    response%nbf=nbf
    response%ncoord=ncoord
    response%packed=packed
    response%omega=energies(target)
    allocate(response%mo_a(nbf,nbf),response%mo_b(nbf,nbf), &
      response%fock_a(nbf,nbf),response%fock_b(nbf,nbf),response%x(packed), &
      response%dmo_a(nbf,nbf,ncoord),response%dmo_b(nbf,nbf,ncoord), &
      response%dmo_common(nbf,nbf,ncoord), &
      response%ds(nbf,nbf,ncoord),response%dhcore(nbf,nbf,ncoord), &
      response%dpa(nbf,nbf,ncoord),response%dpb(nbf,nbf,ncoord), &
      response%dfock_a(nbf,nbf,ncoord), &
      response%dfock_b(nbf,nbf,ncoord),response%dax(packed,ncoord), &
      response%dx(packed,ncoord),response%domega(ncoord),source=0.0_dp)
    response%mo_a=stored_mo_a
    response%mo_b=stored_mo_b
    response%x=state_vectors(:,target)
    call unpack_matrix(packed_fock_a,response%fock_a)
    call unpack_matrix(packed_fock_b,response%fock_b)

    time_orbital=0.0_dp;time_xc_derivative=0.0_dp;time_state=0.0_dp
    call system_clock(clock_start,clock_rate)
    call build_rohf_nuclear_response(infos,orbital_response,local_status, &
      require_two_somo=.true.)
    call system_clock(clock_now)
    time_orbital=real(clock_now-clock_start,dp)/real(clock_rate,dp)
    if(local_status/=0) then
      status=-4
      call response%clean()
      return
    end if
    response%dmo_a=orbital_response%dmo_alpha
    response%dmo_b=orbital_response%dmo_beta
    response%dmo_common=orbital_response%dmo_common
    response%ds=orbital_response%ds_ao
    response%dhcore=orbital_response%dhcore_ao
    call orbital_response%clean()

    allocate(pa(nbf,nbf),pb(nbf,nbf),source=0.0_dp)
    pa=matmul(response%mo_a(:,1:nocca), &
      transpose(response%mo_a(:,1:nocca)))
    pb=matmul(response%mo_b(:,1:noccb), &
      transpose(response%mo_b(:,1:noccb)))
    call build_orbital_density_derivatives(response%mo_a,response%dmo_a, &
      nocca,response%dpa,local_status)
    if(local_status==0) call build_orbital_density_derivatives( &
      response%mo_b,response%dmo_b,noccb,response%dpb,local_status)

    if(dft .and. local_status==0) then
      call system_clock(clock_start)
      block
        use mod_dft, only: dft_initialize,dftclean
        use mod_dft_molgrid, only: dft_grid_t
        type(dft_grid_t) :: grid
        allocate(dvxc_a(nbf,nbf,ncoord),dvxc_b(nbf,nbf,ncoord), &
          source=0.0_dp)
        call dft_initialize(infos,basis,grid)
        call mrsf_xc_fock_total_derivative(basis,grid,pa,pb,response%dpa, &
          response%dpb,dvxc_a,dvxc_b,infos,local_status,threshold=0.0_dp)
        call dftclean(infos)
      end block
      call system_clock(clock_now)
      time_xc_derivative=real(clock_now-clock_start,dp)/real(clock_rate,dp)
    end if

    if(local_status==0) then
      call system_clock(clock_start)
      call int2_driver%init(basis,infos)
      call int2_driver%set_screening()
      if(dft) then
        call solve_mrsf_first_nuclear_response(infos,int2_driver, &
          response%mo_a,response%mo_b,response%dmo_a,response%dmo_b, &
          response%dmo_common, &
          response%fock_a,response%fock_b,response%dhcore,response%omega, &
          response%x,response%dax,response%dx,response%domega, &
          response%residual,local_status,dvxc_a,dvxc_b, &
          response%dpa,response%dpb,response%dfock_a,response%dfock_b)
      else
        call solve_mrsf_first_nuclear_response(infos,int2_driver, &
          response%mo_a,response%mo_b,response%dmo_a,response%dmo_b, &
          response%dmo_common, &
          response%fock_a,response%fock_b,response%dhcore,response%omega, &
          response%x,response%dax,response%dx,response%domega, &
          response%residual,local_status,dpa_out=response%dpa, &
          dpb_out=response%dpb,dfock_a_out=response%dfock_a, &
          dfock_b_out=response%dfock_b)
      end if
      call int2_driver%clean()
      call system_clock(clock_now)
      time_state=real(clock_now-clock_start,dp)/real(clock_rate,dp)
    end if
    if(allocated(dvxc_a)) deallocate(dvxc_a,dvxc_b)
    deallocate(pa,pb)
    if(local_status/=0) then
      status=local_status
      call response%clean()
    else
      write(iw,'(A,3(F10.3,1X))') &
        'MRSF first-response substage wall times (s): orbital, XC-Fock, state = ', &
        time_orbital,time_xc_derivative,time_state
    end if
  end subroutine prepare_mrsf_hessian_first_response

!###############################################################################

  subroutine mrsf_hessian_first_response_clean(this)
    class(mrsf_hessian_first_response_t), intent(inout) :: this

    if(allocated(this%mo_a)) deallocate(this%mo_a)
    if(allocated(this%mo_b)) deallocate(this%mo_b)
    if(allocated(this%fock_a)) deallocate(this%fock_a)
    if(allocated(this%fock_b)) deallocate(this%fock_b)
    if(allocated(this%x)) deallocate(this%x)
    if(allocated(this%dmo_a)) deallocate(this%dmo_a)
    if(allocated(this%dmo_b)) deallocate(this%dmo_b)
    if(allocated(this%dmo_common)) deallocate(this%dmo_common)
    if(allocated(this%ds)) deallocate(this%ds)
    if(allocated(this%dhcore)) deallocate(this%dhcore)
    if(allocated(this%dpa)) deallocate(this%dpa)
    if(allocated(this%dpb)) deallocate(this%dpb)
    if(allocated(this%dfock_a)) deallocate(this%dfock_a)
    if(allocated(this%dfock_b)) deallocate(this%dfock_b)
    if(allocated(this%dax)) deallocate(this%dax)
    if(allocated(this%dx)) deallocate(this%dx)
    if(allocated(this%domega)) deallocate(this%domega)
    this%nbf=0
    this%ncoord=0
    this%packed=0
    this%omega=0.0_dp
    this%residual=0.0_dp
  end subroutine mrsf_hessian_first_response_clean

end module tdhf_mrsf_hessian_prepare_mod
