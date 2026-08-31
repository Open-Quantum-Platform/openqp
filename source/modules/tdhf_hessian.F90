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
      tdhf_hessian_is_applicable
    use tdhf_hessian_fixed_density_mod, only: build_tdhf_fixed_density_hessian
    use tdhf_response_operator_mod, only: build_tdhf_response_matrices
    use tdhf_hessian_orbital_mod, only: build_tdhf_ground_orbital_response
    use tdhf_hessian_rhs_mod, only: build_tdhf_amplitude_derivative_actions
    use tdhf_hessian_response_mod, only: solve_tdhf_amplitude_response, solve_tdhf_z_response
    use tdhf_hessian_z_rhs_mod, only: build_tdhf_z_rhs_derivative
    use tdhf_hessian_density_mod, only: build_tdhf_relaxed_density_derivatives
    use tdhf_hessian_rows_mod, only: build_tdhf_response_rows_hf
    use hf_hessian_mod, only: hf_hessian
    use io_constants, only: iw
    use parallel, only: par_env_t
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos
    real(dp),contiguous,pointer::xpy(:,:),xmy(:,:),zstore(:),energies(:),hstore(:,:),hground(:,:)
    real(dp),allocatable::u0(:),v0(:),z0(:),amb(:,:),apb(:,:),sx(:,:,:),umat(:,:,:),umat_ao(:,:,:)
    real(dp),allocatable::deps(:,:),deps_ao(:,:),dground(:,:,:),dambu(:,:),dapbv(:,:),du(:,:),dv(:,:),domega(:)
    real(dp),allocatable::drhs(:,:),dapbz(:,:),dummy(:,:),dz(:,:),dprel(:,:,:),dw(:,:,:),dua(:,:,:),dva(:,:,:)
    real(dp),allocatable::hfixed(:,:),hground_fixed(:,:),hxc(:,:),rows(:,:),rows_one(:,:),rows_two(:,:),htotal(:,:)
    real(dp)::amp_res,z_res,asym,omega
    logical::zero_orbital_connection
    integer::nbf,nocc,nvir,nexc,ncart,target,status
    type(par_env_t) :: pe

    call pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    if (.not.tdhf_hessian_is_applicable(infos%control%scftype,infos%tddft%mult, &
        logical(infos%tddft%tda,kind=kind(.false.)),.false.,int(pe%size))) then
      call show_message('Analytic TD Hessian requires serial closed-shell restricted singlet full response.',WITH_ABORT)
    end if
    if(infos%control%hamilton>=20) call show_message(&
      'Analytic TDDFT Hessian XC kernel derivatives are not complete; TDHF is currently supported.',WITH_ABORT)
    nbf=infos%basis%nbf; nocc=infos%mol_prop%nocc; nvir=nbf-nocc; nexc=nocc*nvir
    ncart=3*size(infos%atoms%xyz,2); target=infos%tddft%target_state
    call tagarray_get_data(infos%dat,OQP_TD_XPY,xpy); call tagarray_get_data(infos%dat,OQP_TD_XMY,xmy)
    call tagarray_get_data(infos%dat,OQP_TD_Z,zstore); call tagarray_get_data(infos%dat,OQP_TD_ENERGIES,energies)
    ! The coupled solver stores its first block in the A-B channel and its
    ! second block in the A+B channel.  OpenQP tags these as X-Y and X+Y,
    ! respectively.
    allocate(u0(nexc),v0(nexc),z0(nexc)); u0=xmy(:,target); v0=xpy(:,target); z0=zstore; omega=energies(target)
    allocate(amb(nexc,nexc),apb(nexc,nexc)); call build_tdhf_response_matrices(infos,amb,apb)
    allocate(sx(nbf,nbf,ncart),umat(nbf,nbf,ncart),deps(nbf,ncart),dground(nbf,nbf,ncart))
    call build_tdhf_ground_orbital_response(infos,sx,umat,deps,dground)
    allocate(dambu(nexc,ncart),dapbv(nexc,ncart)); call build_tdhf_amplitude_derivative_actions(infos,umat,deps,u0,v0,dambu,dapbv)
    allocate(du(nexc,ncart),dv(nexc,ncart),domega(ncart))
    call solve_tdhf_amplitude_response(amb,apb,omega,u0,v0,dambu,dapbv,du,dv,domega,amp_res,status)
    if(status/=0) call show_message('TDHF amplitude response did not converge.',WITH_ABORT)
    allocate(drhs(nexc,ncart)); call build_tdhf_z_rhs_derivative(infos,umat,v0,u0,dv,du,drhs)
    allocate(dapbz(nexc,ncart),dummy(nexc,ncart)); call build_tdhf_amplitude_derivative_actions(infos,umat,deps,u0*0.0_dp,z0,dummy,dapbz)
    allocate(dz(nexc,ncart)); call solve_tdhf_z_response(apb,drhs,dapbz,dz,z_res,status)
    if(status/=0) call show_message('TDHF derivative Z-vector did not converge.',WITH_ABORT)
    allocate(dprel(nbf,nbf,ncart),dw(nbf,nbf,ncart),dua(nbf,nbf,ncart),dva(nbf,nbf,ncart))
    allocate(umat_ao(nbf,nbf,ncart),deps_ao(nbf,ncart)); umat_ao=umat; deps_ao=deps
    zero_orbital_connection=.true.
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
    allocate(hfixed(ncart,ncart),hground_fixed(ncart,ncart),hxc(ncart,ncart),rows(ncart,ncart),rows_one(ncart,ncart), &
      rows_two(ncart,ncart),htotal(ncart,ncart)); hxc=0.0_dp
    call build_tdhf_fixed_density_hessian(infos,hfixed,hground_fixed)
    call build_tdhf_response_rows_hf(infos,umat,deps,dground,dprel,dw,dua,dva,rows,rows_one,rows_two)
    ! Replace the ground-state fixed-density skeleton by the complete native
    ! RHF Hessian.  The remaining fixed and response terms are strictly the
    ! excitation-energy contribution, avoiding both omission and double count
    ! of the ground-state CPHF response.
    call hf_hessian(infos)
    call tagarray_get_data(infos%dat,OQP_hf_hessian,hground)
    if (.not.zero_orbital_connection) hfixed=hfixed-hground_fixed+hground
    call assemble_tdhf_cartesian_hessian(hfixed,hxc,rows,htotal,asym,status)
    if(status/=0) call show_message('TDHF Hessian assembly failed.',WITH_ABORT)
    call infos%dat%alloc_or_die(OQP_tdhf_hessian,(/ncart,ncart/),hstore, &
      description='Native OpenQP analytic TDHF Hessian matrix')
    hstore=htotal
    open(unit=iw,file=infos%log_filename,position='append')
    write(iw,'(/,A,1P,E12.4)') 'TDHF Hessian maximum amplitude-response residual: ',amp_res
    write(iw,'(A,1P,E12.4)') 'TDHF Hessian maximum Z-response residual: ',z_res
    write(iw,'(A,1P,E12.4)') 'TDHF Hessian unsymmetrized response-row asymmetry: ',asym
    close(iw)
  end subroutine tdhf_hessian

end module tdhf_hessian_mod
