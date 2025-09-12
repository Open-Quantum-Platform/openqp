!==============================================================================
! Module: trah_opentrustregion_interface
! A simple interface to call OpenTrustRegion solver from OpenQP without needing
! opentrustregion_solver_class
!==============================================================================
module otr_interface

  use, intrinsic :: iso_fortran_env, only: int32
  use opentrustregion, only: solver, update_orbs_type, obj_func_type, hess_x_type, logger_type
  use mathlib, only: unpack_matrix
  use scf_converger, only: trah_converger
  use scf_addons,only: calc_fock,compute_energy
  use precision, only: dp
  use types, only:information
  use mod_dft_molgrid, only: dft_grid_t
  use basis_tools,      only: basis_set
  use guess, only:get_ab_initio_density
  implicit none

  ! Module-level state for callbacks
  class(information), pointer :: infos       ! OpenQP information object
  type(dft_grid_t), pointer :: molgrid
  type(trah_converger), pointer :: conv
  real(dp), allocatable :: work1(:,:), work2(:,:)
  integer(int32)            :: n_param      ! number of parameters = nocc*nvir
  integer(int32)            :: max_iter     ! macro iteration limit
  real(dp)                  :: conv_tol     ! convergence tolerance
  integer(int32)            :: verbose      ! verbosity level

contains

  subroutine init_trah_solver(infos_in, molgrid_in, conv_in)
    class(information), intent(inout), target :: infos_in
    type(dft_grid_t), intent(in), target :: molgrid_in
    class(trah_converger), intent(inout), target :: conv_in
    type(basis_set), pointer :: basis
!    character(len=*), intent(in)          :: print_level
    ! Initialize module state
    infos   => infos_in
    molgrid => molgrid_in
    conv => conv_in

    basis => infos%basis
    basis => infos%basis
    n_param = conv%n_param
    max_iter = int(infos%control%maxit, kind=int32)
    conv_tol = infos%control%conv
    verbose  = int(3, kind=int32)
    allocate(work1(conv%nbf,conv%nbf), work2(conv%nbf,conv%nbf))
!    select case (infos%control%scftype)
!    case (1)
!      call get_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens)
!    case (2)
     conv%mo_b = conv%mo_a
     call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,2), conv%mo_b,infos,basis)
     call calc_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens, conv%mo_b)

!    case (3)
!      self%mo_a = self%mo_b
!      call get_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens, conv%mo_b)
!    end select

  end subroutine init_trah_solver

  subroutine run_trah_solver()
    logical(kind=4) :: error
    procedure(update_orbs_type), pointer :: p_update
    procedure(obj_func_type),   pointer :: p_obj
    procedure(logger_type),     pointer :: p_log

    ! Bind callbacks
    p_update => update_orbs
    p_obj    => obj_func
    p_log    => logger

    ! Call OpenTrustRegion solver
    call solver(p_update, p_obj, n_param, error, &
!                conv_tol=conv_tol, n_macro=max_iter, &
                verbose=verbose)! , logger=p_log)
    if (error) then
      write(*,*) 'OpenTrustRegion solver failed.'
    end if
  end subroutine run_trah_solver

  subroutine update_orbs(kappa, func, grad, h_diag, hess_x_funptr)
    real(dp), intent(in)                     :: kappa(:)
    real(dp), intent(out)                    :: func
    real(dp), intent(out)                    :: grad(:), h_diag(:)
    procedure(hess_x_type), pointer, intent(out) :: hess_x_funptr
    type(basis_set), pointer :: basis

    basis => infos%basis
    ! Rotate orbitals
    select case (infos%control%scftype)
    case (1)
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, conv%mo_a)
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,1), conv%mo_a,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens)
      call conv%calc_g_h(grad, h_diag)
    case (2)
      call conv%rotate_orbs(kappa(1:conv%nocc_a*conv%nvir_a), conv%nbf, conv%nocc_a, conv%mo_a)
      call conv%rotate_orbs(kappa(conv%nocc_a*conv%nvir_a+1:), conv%nbf, conv%nocc_b, conv%mo_b)
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,2), conv%mo_b,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens, conv%mo_b)
      call conv%calc_g_h(grad, h_diag)
    case (3)
      call conv%rotate_orbs(kappa, conv%nbf, conv%nocc_a, conv%mo_a)
      conv%mo_b = conv%mo_a
      call get_ab_initio_density(conv%dens(:,1), conv%mo_a, conv%dens(:,2), conv%mo_b,infos,basis)
      call calc_fock(basis, infos, molgrid, conv%fock_ao, conv%mo_a, conv%dens, conv%mo_b)
      call conv%calc_g_h(grad, h_diag)
    end select
    func = compute_energy(infos)
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
