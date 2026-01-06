!> @brief Thin interface between OpenQP and OpenTrustRegion.
!> @detail Provides callback glue so OpenTrustRegion’s generic trust-region
!>         optimizer can drive OpenQP’s TRAH SCF updates without requiring
!>         an external solver class. Exposes:
!>           - init_trah_solver: bind OpenQP state to module pointers
!>           - run_trah_solver : configure and invoke the OTR solver
!>           - update_orbs     : objective/gradient/Hdiag callback
!>           - hess_x_cb       : Hessian–vector product callback
!>           - obj_func        : energy-only evaluation (trial move)
!>           - logger          : forwards OTR log lines to OpenQP I/O
!> @author Mohsen Mazaherifar
!> @date August 2025
module otr_interface

  use, intrinsic :: iso_c_binding, only: c_bool
  use opentrustregion, only: solver, update_orbs_type,&
          obj_func_type, hess_x_type, logger_type, rp, ip, &
          solver_settings_type, stability_settings_type, &
          default_solver_settings
  use mathlib, only: unpack_matrix
  use scf_converger, only: trah_converger, scf_conv_trah_result, scf_conv_result
  use scf_addons,only: compute_energy,calc_fock
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
  integer :: iter_otr
  real(dp) :: grad_norm

  real(dp), allocatable :: work1(:,:), work2(:,:)


contains
  !> @brief Initialize the OTR–OpenQP bridge and working buffers.
  !> @detail Stores references to OpenQP objects (infos, molgrid, TRAH converger,
  !>         energy accumulator), allocates temporary work arrays, and zeros the
  !>         incremental Fock/Density buffers used for ΔD updates.
  !> @param[inout] infos_in   OpenQP information/control object (target).
  !> @param[in]    molgrid_in DFT molecular grid (target).
  !> @param[inout] conv_in    TRAH converger (provides MO/D/Fock buffers).
  !> @param[inout] energy_in  SCF energy structure to be updated.
  !> @author Mohsen Mazaherifar
  !> @date August 2025
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
    iter_otr = 0

    basis => infos%basis
    allocate(work1(conv%nbf,conv%nbf), work2(conv%nbf,conv%nbf))

    conv%f_old = 0.0_dp
    conv%d_old = 0.0_dp

  end subroutine init_trah_solver

  !> @brief Configure and run the OpenTrustRegion driver.
  !> @detail Wires the required callbacks (`update_orbs`, `obj_func`, `logger`),
  !>         maps OpenQP control flags to OTR options (stability, line-search,
  !>         Davidson/Jacobi–Davidson, trust-radius settings), executes the solve,
  !>         and returns iteration/error status in `res`.
  !> @param[inout] res  Output SCF converger result (TRAH-specific fields filled).
  !> @note Updates the active OpenQP buffers (MO/Fock) upon return.
  !> @author Mohsen Mazaherifar
  !> @date August 2025
  subroutine run_trah_solver(res)
    procedure(update_orbs_type), pointer :: p_update
    procedure(obj_func_type),   pointer :: p_obj
    procedure(logger_type),     pointer :: p_log
    class(scf_conv_result), intent(inout) :: res
    type(solver_settings_type) :: settings
    logical(kind=4) :: stability, line_search, davidson,&
                         jacobi_davidson, prefer_jacobi_davidson
    integer(ip) :: error, n_random_trial_vectors, n_micro,&
                        n_param, max_iter, verbose
    real(dp) :: start_trust_radius, global_red_factor,&
                          local_red_factor, conv_tol

    n_param = conv%n_param
    max_iter = int(infos%control%maxit, kind=ip)
    conv_tol = real(infos%control%conv, kind=rp)
    verbose  = int(3, kind=ip)
    settings = default_solver_settings
    settings%conv_tol = conv_tol
    settings%n_random_trial_vectors = int(infos%control%trh_nrtv, kind=ip)
    settings%jacobi_davidson_start = 30
    settings%seed = 42 
    settings%verbose = verbose
    settings%stability = (infos%control%trh_stab .eqv. .true._c_bool)
    settings%line_search = (infos%control%trh_ls   .eqv. .true._c_bool)
    settings%start_trust_radius = real(infos%control%trh_r0, kind=ip)
    settings%global_red_factor = real(infos%control%trh_gred, kind=ip)
    settings%local_red_factor = real(infos%control%trh_lred, kind=ip)
    settings%n_macro = max_iter
    settings%n_micro = int(infos%control%trh_nmic, kind=ip)
!    settings = solver_settings_type(precond = null(), conv_check = null(), logger = null(), &
!                             stability = .false., line_search = .false., &
!                             initialized = .true., conv_tol = 1e-5_rp, &
!                             start_trust_radius = 0.4_rp, global_red_factor = 1e-3_rp, &
!                             local_red_factor = 1e-4_rp, n_random_trial_vectors = 1, &
!                             n_macro = 150, n_micro = 50, jacobi_davidson_start = 30, &
!                             seed = 42, verbose = 3, subsystem_solver = "davidson")

    settings%logger => logger 

    ! Bind callbacks
    p_update => update_orbs
    p_obj    => obj_func

    call solver(p_update, p_obj, n_param, error, settings)

    conv%dat%buffer(conv%dat%slot)%mo_a = conv%mo_a
    conv%dat%buffer(conv%dat%slot)%focks = conv%fock_ao
    if (infos%control%scftype>1) then
      conv%dat%buffer(conv%dat%slot)%mo_b = conv%mo_b
    end if

    select type (res)
    class is (scf_conv_trah_result)
      res%iter = iter_otr
    end select

    if (error /= 0) then
      write(*,*) 'OpenTrustRegion solver failed.'
      res%ierr = 4
      select type (res)
      class is (scf_conv_trah_result)
        res%iter = max_iter
      end select
    else
      if(grad_norm>conv_tol) then
        write(*,*) 'Trust radius too small. Convergence criterion&
        is not fulfilled but calculation should be converged up to floating&
        point precision.'
        res%error = min(conv_tol*0.99,grad_norm)
      else
        res%error = grad_norm
      end if
    endif


    if (allocated(work1)) deallocate(work1)
    if (allocated(work2)) deallocate(work2)

  end subroutine run_trah_solver

  !> @brief Objective/gradient/Hessian-diagonal callback used by OTR.
  !> @detail Applies orbital rotations `kappa` to (α[,β]) MOs, rebuilds densities,
  !>         constructs Fock via `calc_fock` (using incremental ΔD/ΔF when available),
  !>         then forms orbital-rotation gradient and Hessian diagonal with
  !>         `conv%calc_g_h`. Also binds the Hessian–vector product callback.
  !> @param[in]   kappa    Packed rotation vector(s).
  !> @param[out]  func     Objective value (total electronic energy).
  !> @param[out]  grad     Objective gradient in rotation coordinates.
  !> @param[out]  h_diag   Diagonal of approximate Hessian in rotation space.
  !> @param[out]  hess_x_funptr Pointer to Hessian–vector product routine.
  !> @author Mohsen Mazaherifar
  !> @date August 2025
  subroutine update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
    real(dp), intent(in), target :: kappa(:)
    real(dp), intent(out) :: func
    real(dp), intent(out), target :: grad(:), h_diag(:)
    procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
    integer(ip), intent(out) :: error
    type(basis_set), pointer :: basis
    integer :: nschwz

    basis => infos%basis
    iter_otr = iter_otr + 1
    ! Rotate orbitals
    select case (infos%control%scftype)
    case (1)
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, conv%mo_a)
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,1), conv%mo_a,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, energy, conv%mo_a, conv%dens, conv%mo_b, nschwz, conv%f_old, conv%d_old)
      call conv%calc_g_h(grad, h_diag)
    case (2)
      call conv%rotate_orbs(kappa(1:conv%nocc_a*conv%nvir_a), conv%nbf, conv%nocc_a, conv%mo_a)
      call conv%rotate_orbs(kappa(conv%nocc_a*conv%nvir_a+1:), conv%nbf, conv%nocc_b, conv%mo_b)
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,2), conv%mo_b,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, energy, conv%mo_a, conv%dens, conv%mo_b, nschwz, conv%f_old, conv%d_old)
      call conv%calc_g_h(grad, h_diag)
    case (3)
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, conv%mo_a)
      conv%mo_b = conv%mo_a
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,2), conv%mo_b,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, energy, conv%mo_a, conv%dens, conv%mo_b, nschwz, conv%f_old, conv%d_old)
      call conv%calc_g_h(grad, h_diag)
    end select
    grad_norm = sqrt(dot_product(grad, grad)/conv%n_param)
    func = compute_energy(energy)
    conv%etot = func
    hess_x_funptr => hess_x_cb
    h_diag = 2.0_dp * h_diag
    grad = 2.0_dp * grad
  end subroutine update_orbs

  !> @brief hess_x_cb.
  !> @author Mohsen Mazaherifar
  !> @date August 2025
  subroutine hess_x_cb(x, hx, error)
    real(dp), intent(in), target :: x(:)
    real(dp), intent(out), target :: hx(:)
    integer(ip), intent(out) :: error
    call conv%calc_h_op(infos, x, hx)
    hx = 2.0_dp * hx
  end subroutine hess_x_cb

  !> @brief Energy-only objective for a trial move (no gradient).
  !> @detail Rotates temporary copies of the MOs according to `kappa`, rebuilds
  !>         densities, recomputes Fock and energies, and returns the total energy.
  !>         Used by line-search/auxiliary steps in OTR.
  !> @param[in]  kappa  Packed rotation vector(s).
  !> @return     val    Total electronic energy at the trial point.
  !> @author Mohsen Mazaherifar
  !> @date August 2025
  function obj_func(kappa, error) result(val)
    real(dp), intent(in), target :: kappa(:)
    integer(ip), intent(out) :: error
    real(dp)             :: val
    type(basis_set), pointer :: basis
    integer :: nschwz
    basis => infos%basis
    select case(infos%control%scftype)
    case (1)
      work1 = conv%mo_a
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, work1)
      call get_ab_initio_density(conv%dens(:,1), work1, conv%dens(:,1), work1,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, energy, work1, conv%dens, work1, nschwz, conv%f_old, conv%d_old)
    case (2)
      work1 = conv%mo_a
      work2 = conv%mo_b
      call conv%rotate_orbs(kappa(1:conv%nvir_a*conv%nocc_a), conv%nbf, conv%nocc_a, work1)
      call conv%rotate_orbs(kappa(conv%nvir_a*conv%nocc_a+1:), conv%nbf, conv%nocc_b, work2)
      call get_ab_initio_density(conv%dens(:,1), work1, conv%dens(:,2), work2,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, energy, work1, conv%dens, work2, nschwz, conv%f_old, conv%d_old)
    case (3)
      work1 = conv%mo_a
      work2 = conv%mo_b
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, work1)
      work2 = work1
      call get_ab_initio_density(conv%dens(:,1), work1, conv%dens(:,2), work2,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, energy, work1, conv%dens, work2, nschwz, conv%f_old, conv%d_old) 
    end select 
    val = compute_energy(energy)
  end function obj_func

  subroutine logger(message)
    use io_constants, only: IW
    implicit none
    character(*), intent(in) :: message
    write(IW, "(A)") trim(message)
  end subroutine

end module otr_interface
