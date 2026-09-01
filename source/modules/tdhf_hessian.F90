module tdhf_hessian_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_hessian_mod"

contains

!###############################################################################

  subroutine tdhf_hessian_C(c_handle) bind(C, name="tdhf_hessian")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call tdhf_hessian(inf)
  end subroutine tdhf_hessian_C

!###############################################################################

  subroutine tdhf_hessian(infos)
    use precision, only: dp
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_TD_XPY, OQP_TD_XMY, &
      OQP_TD_Z, OQP_TD_ENERGIES, OQP_hf_hessian, OQP_tdhf_hessian
    use tdhf_hessian_components_mod, only: assemble_tdhf_cartesian_hessian, &
      tdhf_hessian_is_applicable,tdhf_hessian_functional_is_verified, &
      tdhf_hessian_target_is_supported,tdhf_hessian_lowest_root_is_isolated
    use tdhf_hessian_fixed_density_mod, only: build_tdhf_fixed_density_hessian
    use tdhf_response_operator_mod, only: build_tdhf_response_matrices
    use tdhf_hessian_orbital_mod, only: build_tdhf_ground_orbital_response
    use tdhf_hessian_rhs_mod, only: build_tdhf_amplitude_derivative_actions
    use tdhf_hessian_response_mod, only: solve_tdhf_amplitude_response, solve_tdhf_z_response
    use tdhf_hessian_z_rhs_mod, only: build_tdhf_z_rhs_derivative
    use tdhf_hessian_density_mod, only: build_tdhf_relaxed_density_derivatives
    use tdhf_hessian_rows_mod, only: build_tdhf_response_rows_hf
    use tdhf_hessian_xc_mod, only: build_tdhf_xc_fixed_hessian, &
      add_tdhf_xc_response_rows
    use hf_hessian_mod, only: hf_hessian
    use io_constants, only: iw
    use parallel, only: par_env_t
    use messages, only: show_message, WITH_ABORT
    use, intrinsic :: iso_c_binding, only: c_null_char
!$  use omp_lib, only: omp_get_max_threads, omp_set_num_threads

    implicit none

    type(information), target, intent(inout) :: infos
    real(dp),contiguous,pointer::xpy(:,:),xmy(:,:),zstore(:),energies(:),hstore(:,:),hground(:,:)
    real(dp),allocatable::u0(:),v0(:),z0(:),amb(:,:),apb(:,:),sx(:,:,:),umat(:,:,:),umat_ao(:,:,:)
    real(dp),allocatable::deps(:,:),deps_ao(:,:),dground(:,:,:),dxc_skeleton(:,:,:),dambu(:,:),dapbv(:,:)
    real(dp),allocatable::du(:,:),dv(:,:),domega(:)
    real(dp),allocatable::drhs(:,:),dapbz(:,:),dummy(:,:),dz(:,:),dprel(:,:,:),dw(:,:,:),dua(:,:,:),dva(:,:,:)
    real(dp),allocatable::hfixed(:,:),hground_fixed(:,:),hxc(:,:)
    real(dp),allocatable::rows(:,:),rowsxc(:,:),rows_one(:,:),rows_two(:,:),htotal(:,:)
    real(dp)::amp_res,z_res,asym,omega,projection_mean
    logical::zero_orbital_connection
    logical::verified_functional
    character(len=size(infos%dft%xc_functional_name))::functional_name
    integer::nbf,nocc,nvir,nexc,ncart,natom,target,status,cart,atom,omp_saved_threads,i
    type(par_env_t) :: pe

    omp_saved_threads=1
!$  omp_saved_threads=omp_get_max_threads()
!$  call omp_set_num_threads(1)
    call pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    verified_functional=.true.
    if(infos%control%hamilton==20) then
      ! XC_functional_name is a fixed-size C buffer and may occupy all twenty
      ! bytes without a terminating NUL.  Decode it within its declared bound
      ! instead of letting the unbounded C-string helper read adjacent fields.
      functional_name=' '
      do i=1,size(infos%dft%xc_functional_name)
        if(infos%dft%xc_functional_name(i)==c_null_char) exit
        functional_name(i:i)=infos%dft%xc_functional_name(i)
      end do
      verified_functional=tdhf_hessian_functional_is_verified(functional_name)
    end if
    if (.not.tdhf_hessian_is_applicable(infos%control%scftype,infos%tddft%mult, &
        logical(infos%tddft%tda,kind=kind(.false.)), &
        infos%control%hamilton==20, &
        verified_functional, &
        logical(infos%functional%needGrd,kind=kind(.false.)), &
        logical(infos%functional%needTau,kind=kind(.false.)), &
        logical(infos%dft%cam_flag,kind=kind(.false.)), int(pe%size))) then
      call show_message('Analytic TD Hessian requires serial closed-shell restricted '// &
        'singlet full response and supports TDHF or restricted LDA/GGA TDDFT.', &
        WITH_ABORT)
    end if
    nbf=infos%basis%nbf; nocc=infos%mol_prop%nocc; nvir=nbf-nocc; nexc=nocc*nvir
    natom=size(infos%atoms%xyz,2); ncart=3*natom; target=infos%tddft%target_state
    if (.not.tdhf_hessian_target_is_supported(target)) then
      call show_message('Analytic TD Hessian currently supports only the lowest '// &
                        'excited root (target state 1).', WITH_ABORT)
    end if
    call tagarray_get_data(infos%dat,OQP_TD_XPY,xpy); call tagarray_get_data(infos%dat,OQP_TD_XMY,xmy)
    call tagarray_get_data(infos%dat,OQP_TD_Z,zstore); call tagarray_get_data(infos%dat,OQP_TD_ENERGIES,energies)
    if(.not.tdhf_hessian_lowest_root_is_isolated(energies,1.0e-10_dp)) then
      call show_message('Analytic TD Hessian requires at least two computed excited '// &
                        'roots and an isolated lowest root; use tdhf.nstate>=2 or a '// &
                        'numerical Hessian.',WITH_ABORT)
    end if
    ! The coupled solver stores its first block in the A-B channel and its
    ! second block in the A+B channel.  OpenQP tags these as X-Y and X+Y,
    ! respectively.
    allocate(u0(nexc),v0(nexc),z0(nexc)); u0=xmy(:,target); v0=xpy(:,target); z0=zstore; omega=energies(target)
    allocate(amb(nexc,nexc),apb(nexc,nexc)); call build_tdhf_response_matrices(infos,amb,apb)
    allocate(sx(nbf,nbf,ncart),umat(nbf,nbf,ncart),deps(nbf,ncart), &
      dground(nbf,nbf,ncart),dxc_skeleton(nbf,nbf,ncart))
    call build_tdhf_ground_orbital_response(infos,sx,umat,deps,dground,dxc_skeleton)
    allocate(dambu(nexc,ncart),dapbv(nexc,ncart)); call build_tdhf_amplitude_derivative_actions( &
      infos,umat,deps,dground,u0,v0,dambu,dapbv)
    allocate(du(nexc,ncart),dv(nexc,ncart),domega(ncart))
    call solve_tdhf_amplitude_response(amb,apb,omega,u0,v0,dambu,dapbv,du,dv,domega,amp_res,status)
    if(status/=0) call show_message('TDHF amplitude response did not converge.',WITH_ABORT)
    allocate(drhs(nexc,ncart)); call build_tdhf_z_rhs_derivative(infos,umat,v0,u0,dv,du,drhs)
    allocate(dapbz(nexc,ncart),dummy(nexc,ncart)); call build_tdhf_amplitude_derivative_actions( &
      infos,umat,deps,dground,u0*0.0_dp,z0,dummy,dapbz)
    allocate(dz(nexc,ncart)); call solve_tdhf_z_response(apb,drhs,dapbz,dz,z_res,status)
    if(status/=0) call show_message('TDHF derivative Z-vector did not converge.',WITH_ABORT)
    allocate(dprel(nbf,nbf,ncart),dw(nbf,nbf,ncart),dua(nbf,nbf,ncart),dva(nbf,nbf,ncart))
    allocate(umat_ao(nbf,nbf,ncart),deps_ao(nbf,ncart)); umat_ao=umat; deps_ao=deps
    ! This historical minimal-basis shortcut is not valid for KS response:
    ! GAMESS propagates U, orbital-energy, and dD terms through every DFT row.
    zero_orbital_connection=infos%control%hamilton<20
    do status=1,ncart
      if(maxval(abs(umat(:,:,status)-transpose(umat(:,:,status))))>1.0e-10_dp) &
        zero_orbital_connection=.false.
    end do
    if(zero_orbital_connection) then
      ! In the minimal symmetric limit the displaced canonical orbitals and
      ! their eigenvalues are constant in the AO representation.  Suppress
      ! the spurious canonical metric response in every AO density derivative;
      ! the complete native RHF Hessian already supplies the ground response.
      umat_ao=0.0_dp
      deps_ao=0.0_dp
      dground=0.0_dp
    end if
    call build_tdhf_relaxed_density_derivatives(infos,umat_ao,deps_ao,omega,domega,v0,u0,dv,du,z0,dz,dprel,dw,dua,dva)
    allocate(hfixed(ncart,ncart),hground_fixed(ncart,ncart),hxc(ncart,ncart), &
      rows(ncart,ncart),rowsxc(ncart,ncart),rows_one(ncart,ncart), &
      rows_two(ncart,ncart),htotal(ncart,ncart)); hxc=0.0_dp
    call build_tdhf_fixed_density_hessian(infos,hfixed,hground_fixed)
    if(infos%control%hamilton==20) call build_tdhf_xc_fixed_hessian(infos,hxc)
    call build_tdhf_response_rows_hf(infos,umat,deps,dground,dprel,dw,dua,dva, &
      dxc_skeleton,rows,rows_one,rows_two)
    rowsxc=0.0_dp
    if(infos%control%hamilton==20) call add_tdhf_xc_response_rows(infos,dground,dprel,dua,rowsxc)
    rows=rows+rowsxc
    ! Replace the ground-state fixed-density skeleton by the complete native
    ! RHF Hessian.  The remaining fixed and response terms are strictly the
    ! excitation-energy contribution, avoiding both omission and double count
    ! of the ground-state CPHF response.
    call hf_hessian(infos)
    call tagarray_get_data(infos%dat,OQP_hf_hessian,hground)
    if (.not.zero_orbital_connection) hfixed=hfixed-hground_fixed+hground
    call assemble_tdhf_cartesian_hessian(hfixed,hxc,rows,htotal,asym,status)
    if(status/=0) call show_message('TDHF Hessian assembly failed.',WITH_ABORT)
    ! Apply the translational projector on both Cartesian indices.  The exact
    ! molecular Hessian annihilates all three rigid translations; enforcing
    ! this acoustic sum rule removes only accumulated atom-centred quadrature
    ! and response-solver residuals while preserving symmetry.
    do cart=1,3
      do status=1,ncart
        projection_mean=0.0_dp
        do atom=1,natom
          projection_mean=projection_mean+htotal(3*(atom-1)+cart,status)
        end do
        projection_mean=projection_mean/real(natom,dp)
        do atom=1,natom
          htotal(3*(atom-1)+cart,status)=htotal(3*(atom-1)+cart,status)-projection_mean
        end do
      end do
      do status=1,ncart
        projection_mean=0.0_dp
        do atom=1,natom
          projection_mean=projection_mean+htotal(status,3*(atom-1)+cart)
        end do
        projection_mean=projection_mean/real(natom,dp)
        do atom=1,natom
          htotal(status,3*(atom-1)+cart)=htotal(status,3*(atom-1)+cart)-projection_mean
        end do
      end do
    end do
    call infos%dat%alloc_or_die(OQP_tdhf_hessian,(/ncart,ncart/),hstore, &
      description='Native OpenQP analytic TDHF Hessian matrix')
    hstore=htotal
    open(unit=iw,file=infos%log_filename,position='append')
    write(iw,'(/,A,1P,E12.4)') 'TDHF Hessian maximum amplitude-response residual: ',amp_res
    write(iw,'(A,1P,E12.4)') 'TDHF Hessian maximum Z-response residual: ',z_res
    write(iw,'(A,1P,E12.4)') 'TDHF Hessian unsymmetrized response-row asymmetry: ',asym
    close(iw)
!$  call omp_set_num_threads(omp_saved_threads)
  end subroutine tdhf_hessian

end module tdhf_hessian_mod
