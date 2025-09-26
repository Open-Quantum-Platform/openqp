!==============================================================================
! Module: trah_opentrustregion_interface
! A simple interface to call OpenTrustRegion solver from OpenQP without needing
! opentrustregion_solver_class
!==============================================================================
module otr_interface

  use, intrinsic :: iso_c_binding, only: c_bool
  use opentrustregion, only: solver, update_orbs_type,&
          obj_func_type, hess_x_type, logger_type, rp, ip
  use mathlib, only: unpack_matrix
  use scf_converger, only: trah_converger, scf_conv_trah_result, scf_conv_result
  use scf_addons,only: calc_fock,compute_energy,calc_fock2
  use precision, only: dp
  use types, only:information
  use mod_dft_molgrid, only: dft_grid_t
  use basis_tools,      only: basis_set
  use guess, only:get_ab_initio_density
  use scf_addons, only: scf_energy_t

  implicit none

  ! Module-level state for callbacks
  class(information), pointer :: infos       ! OpenQP information object
  type(dft_grid_t), pointer :: molgrid
  type(trah_converger), pointer :: conv
  type(scf_energy_t), pointer :: energy

  real(dp), allocatable :: work1(:,:), work2(:,:)


contains

  subroutine init_trah_solver(infos_in, molgrid_in, conv_in, energy_in)
    class(information), intent(inout), target :: infos_in
    type(dft_grid_t), intent(in), target :: molgrid_in
    class(trah_converger), intent(inout), target :: conv_in
    class(scf_energy_t), intent(inout), target :: energy_in
    type(basis_set), pointer :: basis
    ! Initialize module state
    infos   => infos_in
    molgrid => molgrid_in
    conv => conv_in
    energy => energy_in

    basis => infos%basis
    allocate(work1(conv%nbf,conv%nbf), work2(conv%nbf,conv%nbf))

  end subroutine init_trah_solver

  subroutine run_trah_solver(res)
    procedure(update_orbs_type), pointer :: p_update
    procedure(obj_func_type),   pointer :: p_obj
    procedure(logger_type),     pointer :: p_log
    class(scf_conv_result), intent(inout) :: res
    logical(kind=4) :: error, stability, line_search, davidson,&
                         jacobi_davidson, prefer_jacobi_davidson
    integer(ip) :: n_random_trial_vectors, n_micro,&
                        n_param, max_iter, verbose
    real(dp) :: start_trust_radius, global_red_factor,&
                          local_red_factor, conv_tol

    n_param = conv%n_param
    max_iter = int(infos%control%maxit, kind=ip)
    conv_tol = real(infos%control%conv, kind=rp)
    verbose  = int(3, kind=ip)
    stability              = (infos%control%trh_stab .eqv. .true._c_bool)
    line_search            = (infos%control%trh_ls   .eqv. .true._c_bool)
    davidson               = (infos%control%trh_dav  .eqv. .true._c_bool)
    jacobi_davidson        = (infos%control%trh_jd   .eqv. .true._c_bool)
    prefer_jacobi_davidson = (infos%control%trh_pjd  .eqv. .true._c_bool)

    n_random_trial_vectors = int(infos%control%trh_nrtv, kind=ip)
    start_trust_radius     = real(infos%control%trh_r0, kind=ip)
    n_micro                = int(infos%control%trh_nmic, kind=ip)
    global_red_factor      = real(infos%control%trh_gred, kind=ip)
    local_red_factor       = real(infos%control%trh_lred, kind=ip)

    ! Bind callbacks
    p_update => update_orbs
    p_obj    => obj_func
    p_log    => logger

    call solver(p_update, p_obj, n_param, error=error, &
                stability=stability, line_search=line_search, davidson=davidson, &
                jacobi_davidson=jacobi_davidson, &
                prefer_jacobi_davidson=prefer_jacobi_davidson, &
                conv_tol=conv_tol, n_random_trial_vectors=n_random_trial_vectors, &
                start_trust_radius=start_trust_radius, n_macro=max_iter, &
                n_micro=n_micro, global_red_factor=global_red_factor, &
                local_red_factor=local_red_factor, verbose=verbose)


    conv%dat%buffer(conv%dat%slot)%mo_a = conv%mo_a
    conv%dat%buffer(conv%dat%slot)%focks = conv%fock_ao
    if (infos%control%scftype>1) then
      conv%dat%buffer(conv%dat%slot)%mo_b = conv%mo_b
    end if

    select type (res)
    class is (scf_conv_trah_result)
      res%etot = conv%etot
    end select

    res%error = 0
    if (error) then
      write(*,*) 'OpenTrustRegion solver failed.'
      res%error = 4
    end if

    if (allocated(work1)) deallocate(work1)
    if (allocated(work2)) deallocate(work2)

  end subroutine run_trah_solver

  subroutine update_orbs(kappa, func, grad, h_diag, hess_x_funptr)
    real(dp), intent(in)                     :: kappa(:)
    real(dp), intent(out)                    :: func
    real(dp), intent(out)                    :: grad(:), h_diag(:)
    procedure(hess_x_type), pointer, intent(out) :: hess_x_funptr
    type(basis_set), pointer :: basis
    integer :: nschwz

    basis => infos%basis
    ! Rotate orbitals
    select case (infos%control%scftype)
    case (1)
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, conv%mo_a)
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,1), conv%mo_a,infos,basis)
      call calc_fock2(basis, infos, molgrid, conv%fock_ao, energy, conv%mo_a, conv%dens, nschwz=nschwz)
!      call calc_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens)
      call conv%calc_g_h(grad, h_diag)
    case (2)
      call conv%rotate_orbs(kappa(1:conv%nocc_a*conv%nvir_a), conv%nbf, conv%nocc_a, conv%mo_a)
      call conv%rotate_orbs(kappa(conv%nocc_a*conv%nvir_a+1:), conv%nbf, conv%nocc_b, conv%mo_b)
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,2), conv%mo_b,infos,basis)
      call calc_fock2(basis, infos, molgrid, conv%fock_ao, energy, conv%mo_a, conv%dens, conv%mo_b, nschwz=nschwz)

!      call calc_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens, conv%mo_b)
      call conv%calc_g_h(grad, h_diag)
    case (3)
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, conv%mo_a)
      conv%mo_b = conv%mo_a
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,2), conv%mo_b,infos,basis)
      call calc_fock2(basis, infos, molgrid, conv%fock_ao, energy, conv%mo_a, conv%dens, conv%mo_b, nschwz=nschwz)
!      call calc_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens, conv%mo_b)
      call conv%calc_g_h(grad, h_diag)
    end select
    func = energy%etot!compute_energy(infos)
    conv%etot = func
    hess_x_funptr => hess_x_cb
    h_diag = 2.0_dp * h_diag
    grad = 2.0_dp * grad
  end subroutine update_orbs


  function hess_x_cb(x) result(hx)
    real(dp), intent(in) :: x(:)
    real(dp)             :: hx(size(x))
    call conv%calc_h_op(infos, x, hx)
    hx = 2.0_dp * hx
  end function hess_x_cb


  function obj_func(kappa) result(val)
    real(dp), intent(in) :: kappa(:)
    real(dp)             :: val
    type(basis_set), pointer :: basis
    basis => infos%basis
    select case(infos%control%scftype)
    case (1)
      work1 = conv%mo_a
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, work1)
      call get_ab_initio_density(conv%dens(:,1), work1, conv%dens(:,1), work1,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, work1, conv%dens)
      val = compute_energy(infos)
    case (2)
      work1 = conv%mo_a
      work2 = conv%mo_b
      call conv%rotate_orbs(kappa(1:conv%nvir_a*conv%nocc_a), conv%nbf, conv%nocc_a, work1)
      call conv%rotate_orbs(kappa(conv%nvir_a*conv%nocc_a+1:), conv%nbf, conv%nocc_b, work2)
      call get_ab_initio_density(conv%dens(:,1), work1, conv%dens(:,2), work2,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, work1, conv%dens, work2)
      val = compute_energy(infos)
    case (3)
      work1 = conv%mo_a
      work2 = conv%mo_b
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, work1)
      work2 = work1
      call get_ab_initio_density(conv%dens(:,1), work1, conv%dens(:,2), work2,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, work1, conv%dens, work2)
      val = compute_energy(infos)
    end select
  end function obj_func


  subroutine logger(message)
    use io_constants, only: IW
    implicit none
    character(*), intent(in) :: message
    write(IW, "(A)") trim(message)
  end subroutine

end module otr_interface
