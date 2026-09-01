module tdhf_mrsf_hessian_mod

  implicit none

  private
  public :: tdhf_mrsf_hessian,tdhf_mrsf_hessian_C

contains

!###############################################################################

  subroutine tdhf_mrsf_hessian_C(c_handle) bind(C,name='tdhf_mrsf_hessian')
    use c_interop, only: oqp_handle_t,oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: infos
    infos=>oqp_handle_get_info(c_handle)
    call tdhf_mrsf_hessian(infos)
  end subroutine tdhf_mrsf_hessian_C

!###############################################################################

  subroutine tdhf_mrsf_hessian(infos)
    ! Native analytical Hessian of the founding two-SOMO MRSF-TDHF method.
    ! Hiroya Nakata's TDHF/TDDFT analytical Hessian is the methodological
    ! starting point; all state quantities below remain in the physical
    ! spin-adapted seven-density representation.

    use precision, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use mod_dft_molgrid, only: dft_grid_t
    use oqp_tagarray_driver, only: tagarray_get_data,data_has_tags, &
      OQP_DM_A,OQP_DM_B,OQP_td_p,OQP_td_mrsf_density,OQP_td_z,OQP_WAO, &
      OQP_tdhf_hessian
    use mathlib, only: unpack_matrix
    use tdhf_mrsf_hessian_components_mod, only: &
      mrsf_hessian_is_applicable,assemble_mrsf_cartesian_hessian, &
      project_mrsf_rigid_translations
    use tdhf_mrsf_hessian_prepare_mod, only: &
      mrsf_hessian_first_response_t,prepare_mrsf_hessian_first_response
    use tdhf_mrsf_hessian_intermediates_mod, only: &
      mrsf_hessian_z_intermediates_t,build_mrsf_hf_z_intermediates, &
      build_mrsf_hf_w_intermediates
    use tdhf_mrsf_hessian_z_response_mod, only: &
      solve_mrsf_z_response_from_mo_derivatives
    use tdhf_mrsf_hessian_w_response_mod, only: &
      build_mrsf_w_ao,build_mrsf_w_ao_derivative
    use tdhf_mrsf_hessian_fixed_density_mod, only: &
      build_mrsf_fixed_density_hessian
    use tdhf_mrsf_hessian_rows_mod, only: build_tdhf_mrsf_response_rows
    use mod_dft_gridint_mrsf_xc_hessian, only: &
      mrsf_tddft_xc_hessian_rows
    use tdhf_mrsf_z_vector_mod, only: apply_z_operator
    use parallel, only: par_env_t
    use mod_dft, only: dft_initialize,dftclean
    use messages, only: show_message,WITH_ABORT
    use io_constants, only: iw
!$  use omp_lib, only: omp_get_max_threads,omp_set_num_threads

    type(information), target, intent(inout) :: infos

    character(len=*), parameter :: required(*)=[character(len=80) :: &
      OQP_DM_A,OQP_DM_B,OQP_td_p,OQP_td_mrsf_density,OQP_td_z,OQP_WAO]
    type(basis_set), pointer :: basis
    type(int2_compute_t), target :: int2_driver
    type(dft_grid_t), target :: unused_grid
    type(mrsf_hessian_first_response_t) :: response
    type(mrsf_hessian_z_intermediates_t) :: intermediate
    type(par_env_t) :: pe
    real(kind=dp), contiguous, pointer :: dm_a(:),dm_b(:),td_p(:,:), &
      seven(:,:,:),z(:),wao_packed(:),hstore(:,:)
    real(kind=dp), allocatable :: drhs(:,:),operator_derivative_z(:,:),dz(:,:), &
      drelaxed_a(:,:,:),drelaxed_b(:,:,:),dw(:,:,:),w_baseline(:,:), &
      w_stored(:,:),reference_spin(:,:,:),reference_fock(:,:,:), &
      relaxed_spin(:,:,:),dreference_spin(:,:,:,:), &
      dreference_fock(:,:,:,:),drelaxed_spin(:,:,:,:),hfixed(:,:), &
      hxc(:,:),xc_rows(:,:),rows(:,:),rows_one(:,:),rows_two(:,:), &
      rows_xc(:,:),htotal(:,:)
    real(kind=dp) :: amplitude_residual,z_residual,row_asymmetry,w_error, &
      orbital_exchange_scale
    integer :: nbf,nbf2,natom,ncoord,nocca,noccb,lzdim,status,local_status
    integer :: omp_saved_threads
    character(len=160) :: error_message
    logical :: is_dft

    omp_saved_threads=1
!$  omp_saved_threads=omp_get_max_threads()
!$  call omp_set_num_threads(1)
    call pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    nbf=infos%basis%nbf
    nbf2=nbf*(nbf+1)/2
    natom=size(infos%atoms%xyz,2)
    ncoord=3*natom
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    is_dft=infos%control%hamilton==20
    orbital_exchange_scale=1.0_dp
    if(is_dft) orbital_exchange_scale=infos%dft%hfscale
    lzdim=noccb*(nocca-noccb+nbf-nocca)+(nocca-noccb)*(nbf-nocca)
    if(.not.mrsf_hessian_is_applicable(infos%control%scftype, &
       infos%mol_prop%mult,nocca,noccb,infos%tddft%mult, &
       logical(infos%tddft%umrsf,kind=kind(.false.)), &
       infos%control%hamilton==20, &
       logical(infos%functional%needTau,kind=kind(.false.)), &
       logical(infos%functional%needLapl,kind=kind(.false.)), &
       logical(infos%dft%cam_flag,kind=kind(.false.)),int(pe%size))) then
      call show_message('Analytic MRSF Hessian requires a serial two-SOMO '// &
        'triplet ROHF reference and a singlet or triplet MRSF target.',WITH_ABORT)
    end if
    call data_has_tags(infos%dat,required,'tdhf_mrsf_hessian_mod', &
      'tdhf_mrsf_hessian',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_DM_A,dm_a)
    call tagarray_get_data(infos%dat,OQP_DM_B,dm_b)
    call tagarray_get_data(infos%dat,OQP_td_p,td_p)
    call tagarray_get_data(infos%dat,OQP_td_mrsf_density,seven)
    call tagarray_get_data(infos%dat,OQP_td_z,z)
    call tagarray_get_data(infos%dat,OQP_WAO,wao_packed)
    if(size(z)/=lzdim) call show_message( &
      'Stored MRSF Z-vector has the wrong dimension.',WITH_ABORT)

    call prepare_mrsf_hessian_first_response(infos,response,status)
    if(status/=0) then
      write(error_message,'(A,I0,A)') &
        'MRSF first nuclear response failed with status ',status,'.'
      call show_message(trim(error_message),WITH_ABORT)
    end if
    amplitude_residual=response%residual
    if(.not.ieee_is_finite(amplitude_residual) .or. &
       amplitude_residual>1.0e-6_dp .or. &
       any(.not.ieee_is_finite(response%dx)) .or. &
       any(.not.ieee_is_finite(response%dax))) call show_message( &
      'MRSF amplitude response contains a non-finite or unconverged result.', &
      WITH_ABORT)
    allocate(reference_spin(nbf,nbf,2),reference_fock(nbf,nbf,2), &
      dreference_spin(nbf,nbf,2,ncoord), &
      dreference_fock(nbf,nbf,2,ncoord),source=0.0_dp)
    call unpack_matrix(dm_a,reference_spin(:,:,1))
    call unpack_matrix(dm_b,reference_spin(:,:,2))
    reference_fock(:,:,1)=response%fock_a
    reference_fock(:,:,2)=response%fock_b
    dreference_spin(:,:,1,:)=response%dpa
    dreference_spin(:,:,2,:)=response%dpb
    dreference_fock(:,:,1,:)=response%dfock_a
    dreference_fock(:,:,2,:)=response%dfock_b
    basis=>infos%basis
    basis%atoms=>infos%atoms
    call int2_driver%init(basis,infos)
    call int2_driver%set_screening()
    call build_mrsf_hf_z_intermediates(infos,int2_driver,response%mo_a, &
      response%mo_b,response%dmo_common,response%dmo_common,response%fock_a, &
      response%fock_b,response%dfock_a,response%dfock_b,response%x, &
      response%dx,z,reference_spin,dreference_spin,intermediate,status)
    if(status/=0) call show_message( &
      'MRSF Z-response intermediates could not be built.',WITH_ABORT)

    allocate(drhs(lzdim,ncoord),operator_derivative_z(lzdim,ncoord), &
      dz(lzdim,ncoord),source=0.0_dp)
    if(is_dft) call dft_initialize(infos,basis,unused_grid)
    call solve_mrsf_z_response_from_mo_derivatives(orbital_hessian_action, &
      infos%tddft%mult,logical(infos%tddft%umrsf,kind=kind(.false.)), &
      .false.,nocca,noccb,intermediate%mo_energy,intermediate%fa, &
      intermediate%fb,intermediate%z_rhs_hxa,intermediate%z_rhs_hxb, &
      intermediate%ab1_a,intermediate%ab1_b,intermediate%z_ab1_a, &
      intermediate%z_ab1_b,intermediate%tij,intermediate%tab,z, &
      intermediate%dmo_energy,intermediate%dfa,intermediate%dfb, &
      intermediate%dz_rhs_hxa,intermediate%dz_rhs_hxb, &
      intermediate%dab1_a, &
      intermediate%dab1_b,intermediate%dz_ab1_a, &
      intermediate%dz_ab1_b,intermediate%dtij,intermediate%dtab,drhs, &
      operator_derivative_z,dz,z_residual,status, &
      ! The differentiated equation is driven by ROHF CPHF quantities whose
      ! production convergence is 1e-9.  Requiring a matrix-free solve below
      ! 1e-12 is not meaningful and can reject an otherwise converged response.
      tol=max(1.0e-10_dp,infos%tddft%zvconv), &
      maxit=max(200,int(infos%control%maxit_zv)), &
      restart=lzdim)
    if(is_dft) call dftclean(infos)
    if(status/=0) call show_message( &
      'MRSF differentiated Z-vector did not converge.',WITH_ABORT)
    if(.not.ieee_is_finite(z_residual) .or. z_residual>1.0e-6_dp .or. &
       any(.not.ieee_is_finite(dz))) call show_message( &
      'MRSF differentiated Z-vector contains a non-finite or unconverged '// &
      'result.',WITH_ABORT)

    allocate(drelaxed_a(nbf,nbf,ncoord),drelaxed_b(nbf,nbf,ncoord), &
      source=0.0_dp)
    call build_mrsf_hf_w_intermediates(infos,response%mo_a,response%mo_b, &
      response%dmo_common,response%dmo_common,reference_spin, &
      dreference_spin,intermediate%tij, &
      intermediate%tab, &
      intermediate%dtij,intermediate%dtab,z,dz,intermediate,drelaxed_a, &
      drelaxed_b,status)
    if(status/=0) then
      write(error_message,'(A,I0,A)') &
        'MRSF relaxed-density response failed with status ',status,'.'
      call show_message(trim(error_message),WITH_ABORT)
    end if
    allocate(dw(nbf,nbf,ncoord),w_baseline(nbf,nbf),w_stored(nbf,nbf), &
      source=0.0_dp)
    call build_mrsf_w_ao(response%mo_a,intermediate%mo_energy, &
      intermediate%fa,intermediate%fb,z,intermediate%hxa, &
      intermediate%hxb,intermediate%ppija,intermediate%ppijb,nocca,noccb, &
      .true.,3,infos%tddft%mult,.false.,.false.,w_baseline,status)
    if(status/=0) call show_message('MRSF baseline W reconstruction failed.', &
      WITH_ABORT)
    call unpack_matrix(wao_packed,w_stored)
    w_error=maxval(abs(w_baseline-w_stored))
    if(w_error>1.0e-8_dp) then
      write(error_message,'(A,1P,E15.7,A)') &
        'MRSF baseline W reconstruction differs from the gradient by ', &
        w_error,'.'
      call show_message(trim(error_message),WITH_ABORT)
    end if
    call build_mrsf_w_ao_derivative(response%mo_a,response%dmo_common, &
      intermediate%mo_energy,intermediate%dmo_energy,intermediate%fa, &
      intermediate%dfa,intermediate%fb,intermediate%dfb,z,dz, &
      intermediate%hxa,intermediate%dhxa,intermediate%hxb, &
      intermediate%dhxb,intermediate%ppija,intermediate%dppija, &
      intermediate%ppijb,intermediate%dppijb,nocca,noccb,.true.,3, &
      infos%tddft%mult,.false.,.false.,dw,status)
    if(status/=0) call show_message('MRSF W response failed.',WITH_ABORT)

    allocate(relaxed_spin(nbf,nbf,2), &
      drelaxed_spin(nbf,nbf,2,ncoord),source=0.0_dp)
    call unpack_matrix(td_p(:,1),relaxed_spin(:,:,1))
    call unpack_matrix(td_p(:,2),relaxed_spin(:,:,2))
    drelaxed_spin(:,:,1,:)=drelaxed_a
    drelaxed_spin(:,:,2,:)=drelaxed_b
    allocate(hfixed(ncoord,ncoord),hxc(ncoord,ncoord),xc_rows(ncoord,ncoord), &
      rows(ncoord,ncoord),rows_one(ncoord,ncoord),rows_two(ncoord,ncoord), &
      rows_xc(ncoord,ncoord),htotal(ncoord,ncoord),source=0.0_dp)
    call build_mrsf_fixed_density_hessian(infos,hfixed,status)
    if(status/=0) call show_message('MRSF fixed-density Hessian failed.', &
      WITH_ABORT)
    if(infos%control%hamilton==20) then
      block
        use mod_dft, only: dft_initialize,dftclean
        use mod_dft_molgrid, only: dft_grid_t
        type(dft_grid_t) :: xc_grid
        call dft_initialize(infos,basis,xc_grid)
        call mrsf_tddft_xc_hessian_rows(basis,xc_grid,reference_spin, &
          relaxed_spin,dreference_spin,drelaxed_spin,hxc,xc_rows,infos,status)
        call dftclean(infos)
      end block
      if(status/=0) call show_message( &
        'MRSF-TDDFT semilocal XC Hessian integration failed.',WITH_ABORT)
    end if
    call build_tdhf_mrsf_response_rows(infos,reference_spin,reference_fock, &
      relaxed_spin,seven,dreference_spin,dreference_fock,drelaxed_spin,dw, &
      intermediate%dseven, &
      intermediate%dabxc,xc_rows,.true.,rows,rows_one,rows_two,rows_xc, &
      status)
    if(status/=0) then
      write(error_message,'(A,I0,A)') &
        'MRSF Cartesian response rows failed with status ',status,'.'
      call show_message(trim(error_message),WITH_ABORT)
    end if
    call assemble_mrsf_cartesian_hessian(hfixed,hxc,rows,htotal, &
      row_asymmetry,status)
    if(status/=0) call show_message('MRSF Hessian assembly failed.',WITH_ABORT)
    call project_mrsf_rigid_translations(htotal,natom,status)
    if(status/=0) call show_message('MRSF Hessian projection failed.',WITH_ABORT)
    call infos%dat%alloc_or_die(OQP_tdhf_hessian,(/ncoord,ncoord/),hstore, &
      description='Native OpenQP analytic spin-adapted MRSF Hessian')
    hstore=htotal
    open(unit=iw,file=infos%log_filename,position='append')
    write(iw,'(/,A,1P,E12.4)') &
      'MRSF Hessian maximum amplitude-response residual: ',amplitude_residual
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum Z-response residual: ',z_residual
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian baseline W reconstruction error: ',w_error
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum amplitude derivative: ',maxval(abs(response%dx))
    write(iw,'(A,1P,100(E16.8,1X))') &
      'MRSF Hessian excitation-energy derivatives: ',response%domega
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum differentiated Z element: ',maxval(abs(dz))
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum relaxed-density derivative: ', &
      max(maxval(abs(drelaxed_a)),maxval(abs(drelaxed_b)))
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum energy-weighted-density derivative: ', &
      maxval(abs(dw))
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum seven-density derivative: ', &
      maxval(abs(intermediate%dseven))
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian unsymmetrized response-row asymmetry: ',row_asymmetry
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum fixed-density element: ',maxval(abs(hfixed))
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum one-electron response-row element: ', &
      maxval(abs(rows_one))
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian maximum two-electron response-row element: ', &
      maxval(abs(rows_two))
    if(infos%control%hamilton==20) then
      write(iw,'(A,1P,E12.4)') &
        'MRSF Hessian maximum fixed XC element: ',maxval(abs(hxc))
      write(iw,'(A,1P,E12.4)') &
        'MRSF Hessian maximum XC response-row element: ',maxval(abs(rows_xc))
    end if
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian one-electron response-row asymmetry: ', &
      maxval(abs(rows_one-transpose(rows_one)))
    write(iw,'(A,1P,E12.4)') &
      'MRSF Hessian two-electron response-row asymmetry: ', &
      maxval(abs(rows_two-transpose(rows_two)))
    close(iw)
    call int2_driver%clean()
    call response%clean()
    call intermediate%clean()
!$  call omp_set_num_threads(omp_saved_threads)

  contains

    subroutine orbital_hessian_action(vector,result,callback_status)
      real(kind=dp), intent(in) :: vector(:)
      real(kind=dp), intent(out) :: result(:)
      integer, intent(out) :: callback_status
      call apply_z_operator(vector,result,infos,basis,unused_grid,int2_driver, &
        nocca,noccb,nbf,response%mo_a,response%mo_b, &
        intermediate%mo_energy,intermediate%fa,intermediate%fb, &
        orbital_exchange_scale,is_dft)
      callback_status=0
      if(any(.not.ieee_is_finite(result))) callback_status=-1
    end subroutine orbital_hessian_action

  end subroutine tdhf_mrsf_hessian

end module tdhf_mrsf_hessian_mod
