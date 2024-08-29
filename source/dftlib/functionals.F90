!> @brief  MODULE functionals
!> @brief  The part of libxc driver
!> @detail This module save information about DFT functional,
!>         which will work
!> @author Igor S. Gerasimov
!> @date   July, 2019
!>         Adding ability for Ground State and TD calculations
module functionals
  use xc_f03_lib_m
  use iso_c_binding, only: C_INT, C_SIZE_T
  use precision, only: fp
  implicit none
  private
  ! Code of errors if LibXC can not perform calculation of some derivatives
  integer, parameter :: ENERGY_ERROR   = 1, & !< energy calculation error
                        FIRST_ERROR    = 2, & !< first derivatives calculation error
                        SECOND_ERROR   = 3, & !< second derivatives calculation error
                        THIRD_ERROR    = 4, & !< third derivatives calculation error
                        MGGA_3RD_ERROR = 5    !< TD-DFT gradient is not allowed for meta-GGA functionals
  type functional_t
    type(xc_f03_func_t),      dimension(:), allocatable, private :: functionals_list
    type(xc_f03_func_info_t), dimension(:), allocatable, private :: functionals_info
    real(kind=fp),            dimension(:), allocatable, private :: coefficients

    logical :: needgrd  = .false. !< toggles calculation of density gradient
    logical :: needtau  = .false. !< toggles calculation of tau (\sum dot_product(\nabla \phi, \nabla \phi))
    logical :: needlapl = .false. !< toggles calculation of lapl (\nabla^2 \rho)
  contains
    procedure :: add_functional, can_calculate, destroy
    procedure :: calc_evxc, calc_evfxc, calc_xc
  end type functional_t
  public functional_t
contains
  !> @brief  Add functional into internal array of functionals
  !> @author Igor S. Gerasimov
  !> @date   July,  2019 --Initial release--
  !> @date   March, 2021 Add optional hfex, alpha, beta, omega parameters
  !> @date   July,  2021 Using messages module
  !> @params func_id             - (in)            internal key of functional in libxc
  !> @params coeff               - (in)            coefficient before functional for parametric schemes like B3LYP of PBE0
  !> @params external_parameters - (in, optional)  setting external parameters for some LibXC functionals
  !> @params hfex                - (out, optional) returns % of HF of added functional (should be used only for hybrid functionals)
  !> @params alpha               - (out, optional) returns short-range % of HF of added range-saparated functional (should be used only for range-saparated hybrid functionals)
  !> @params beta                - (out, optional) returns long-range % of HF of added range-saparated functional (should be used only for range-saparated hybrid functionals)
  !> @params omega               - (out, optional) returns error function parameter of added range-saparated functional (should be used only for range-saparated hybrid functionals)
  subroutine add_functional(this, func_id, coeff, external_parameters, hfex, alpha, beta, omega)
    use messages, only: show_message, WITH_ABORT
    class(functional_t), intent(inout)                     :: this
    integer(C_INT),    intent(in)                          :: func_id
    real(kind=fp),     intent(in)                          :: coeff
    real(kind=fp),     dimension(:), intent(in),  optional :: external_parameters
    real(kind=fp),                   intent(out), optional :: hfex, alpha, beta, omega
    ! internal variables
    type(xc_f03_func_t),      dimension(:), allocatable :: tmp_functionals
    type(xc_f03_func_info_t), dimension(:), allocatable :: tmp_functionals_info
    real(kind=fp),            dimension(:), allocatable :: tmp_coefficients
    type(xc_f03_func_reference_t)                       :: xc_ref
    type(xc_f03_func_t)                                 :: xc_func
    integer(C_INT)                                      :: refnum
    integer                                             :: refnumt
    if(.not.allocated(this%functionals_list)) then
      allocate(this%functionals_list(0))
      allocate(this%functionals_info(0))
      allocate(this%coefficients    (0))
    end if
    allocate(tmp_functionals     (size(this%functionals_list)+1))
    allocate(tmp_functionals_info(size(this%functionals_info)+1))
    allocate(tmp_coefficients    (size(this%coefficients    )+1))
    tmp_functionals     (1:size(this%functionals_list)) = this%functionals_list(1:size(this%functionals_list))
    tmp_functionals_info(1:size(this%functionals_info)) = this%functionals_info(1:size(this%functionals_info))
    tmp_coefficients    (1:size(this%coefficients    )) = this%coefficients    (1:size(this%coefficients    ))
    call move_alloc(tmp_functionals     , this%functionals_list)
    call move_alloc(tmp_functionals_info, this%functionals_info)
    call move_alloc(tmp_coefficients    , this%coefficients    )
    call xc_f03_func_init(xc_func, func_id, XC_POLARIZED)
    select case(xc_f03_func_info_get_kind(xc_f03_func_get_info(xc_func)))
      case(XC_EXCHANGE)
        call show_message("(A,ES16.8E2,A)", "The " // trim(xc_f03_func_info_get_name(xc_f03_func_get_info(xc_func))) // &
                   " exchange functional will be used with a coefficient ", coeff, ".")
      case(XC_CORRELATION)
        call show_message("(A,ES16.8E2,A)", "The " // trim(xc_f03_func_info_get_name(xc_f03_func_get_info(xc_func))) // &
          " correlation functional will be used with a coefficient ", coeff, ".")
      case(XC_EXCHANGE_CORRELATION)
        call show_message("(A,ES16.8E2,A)", "The " // trim(xc_f03_func_info_get_name(xc_f03_func_get_info(xc_func))) // &
          " exchange-correlation functional will be used with a coefficient ", coeff, ".")
      case(XC_KINETIC)
        call show_message("(A,ES16.8E2,A)", "The " // trim(xc_f03_func_info_get_name(xc_f03_func_get_info(xc_func))) // &
          " kinetic functional will be used with a coefficient ", coeff, ".")
    end select
    select case (xc_f03_func_info_get_family(xc_f03_func_get_info(xc_func)))
      case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        this%needgrd = .true.
      case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
        this%needtau = .true.
      case(XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
      case default
    end select
    ! checking, that functional can be used.
    if(IAND(xc_f03_func_info_get_flags(xc_f03_func_get_info(xc_func)), XC_FLAGS_DEVELOPMENT) .ne. 0) then
      call show_message("The behavior of this functional can be changed in the next versions of LibXC.")
    end if
    if(IAND(xc_f03_func_info_get_flags(xc_f03_func_get_info(xc_func)), XC_FLAGS_NEEDS_LAPLACIAN) .ne. 0) then
      call show_message("This functional requires laplacian, but the calculation of laplacian is not implemented" // &
          " in the current version of OQP.", WITH_ABORT)
    end if
    if(IAND(xc_f03_func_info_get_flags(xc_f03_func_get_info(xc_func)), XC_FLAGS_VV10) .ne. 0) then
      call show_message("This functional uses VV10 correlation, but the calculation of VV10 correlation is not" // &
          " implemented in the current version of OQP.", WITH_ABORT)
    end if
    ! Then showing referencies
    refnum = 0_C_INT
    call show_message("The functional has been described in the following articles:")
    do refnumt = 1, XC_MAX_REFERENCES
      xc_ref = xc_f03_func_info_get_references(xc_f03_func_get_info(xc_func),refnum)
      call show_message("(A,I1,A)", "[", refnumt, "] " // trim(xc_f03_func_reference_get_ref(xc_ref)) // "; DOI: " // &
        trim(xc_f03_func_reference_get_doi(xc_ref)))
      if(refnum .lt. 0) exit
    end do
    if(present(external_parameters)) then
      call xc_f03_func_set_ext_params(xc_func, external_parameters)
    end if
    if(present(hfex)) then
      hfex = xc_f03_hyb_exx_coef(xc_func)
    end if
    if(present(alpha).and.present(beta).and.present(omega)) then
      call xc_f03_hyb_cam_coef(xc_func, omega, alpha, beta)
      ! LibXC has alpha as fraction of full exchange and beta - additional fraction for short-range
      ! In the same time, OQP has alpha as a fraction for short-range, and beta - additional fraction for long-range
      alpha = alpha + beta
      beta = -beta
    else if(present(alpha).or.present(beta).or.present(omega)) then
      call show_message("Check this range-separated functional: it has incorrect calling of add_functional routine" // &
          " (see functionals.src and libxc.src)", WITH_ABORT)
    end if
    this%functionals_list(size(this%functionals_list)) = xc_func
    this%functionals_info(size(this%functionals_info)) = xc_f03_func_get_info(xc_func)
    this%coefficients    (size(this%coefficients    )) = coeff
  end subroutine add_functional
  !> @brief  Checking that at least one functional is selected
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 --Initial release--
  logical function can_calculate(this) result(OK)
    class(functional_t), intent(inout) :: this
    OK = (allocated(this%functionals_list) .and. size(this%functionals_list) .gt. 0)
  end function can_calculate
  !> @brief  Destroy internal variables
  !> @author Igor S. Gerasimov
  !> @date   Dec, 2020 --Initial release--
  !> @date   Mar, 2021 Add destruction of functionals_list pointers
  subroutine destroy(this)
    class(functional_t), intent(inout) :: this
    integer :: i
    if(allocated(this%functionals_list)) then
      do i = 1, size(this%functionals_list)
        call xc_f03_func_end(this%functionals_list(i))
      end do
    end if
    if(allocated(this%functionals_list)) deallocate(this%functionals_list)
    if(allocated(this%functionals_info)) deallocate(this%functionals_info)
    if(allocated(this%coefficients    )) deallocate(this%coefficients    )
  end subroutine destroy
  !> @brief  Perform DFT energy and its first derivatives calculation
  !> @author Igor S. Gerasimov
  !> @date   July, 2019
  !> @params NPoints  - (in)  number of points
  !> @params rho      - (in)  density at points
  !> @params sigma    - (in)  gradient of density at points
  !> @params tau      - (in)  local kinetic energy at points
  !> @params lapl     - (in)  laplacian of density at points
  !> @params energy   - (out) energy at points
  !> @params dedrho   - (out) first derivative energy by density
  !> @params dedsigma - (out) first derivative energy by normed gradient
  !> @params dedtau   - (out) first derivative energy by kinetic energy
  !> @params dedlapl  - (out) first derivative energy by laplacian
  subroutine calc_evxc(this, NPoints, rho, sigma, tau, lapl, energy, dedrho, dedsigma, dedtau, dedlapl)
    class(functional_t), intent(inout) :: this
    integer, intent(in) :: NPoints
    real(kind=fp), dimension(*), intent(in)  :: rho, sigma, tau, lapl
    real(kind=fp), dimension(*), intent(out) :: energy
    real(kind=fp), dimension(*), intent(out) :: dedrho, dedsigma, dedtau, dedlapl
    ! iterators
    integer :: i
    ! temporary arrays for containing energy and derivatives
    real(kind=fp), dimension(1*NPoints) :: tmp_energy
    real(kind=fp), dimension(2*NPoints) :: tmp_dEdrho, tmp_dEdtau
    real(kind=fp), dimension(3*NPoints) :: tmp_dEdsigma
    real(kind=fp), dimension(2*NPoints) :: tmp_dedlapl
    ! count of point for LibXC
    integer(C_SIZE_T)                   :: libxc_int
    ! coefficient of functional
    real(kind=fp) :: coefficient
    ! saving total density of each point
    real(kind=fp), dimension(NPoints)   :: rhosum
    ! Build array for multiplying energy at each point
    rhosum = rho(1:2*npoints:2) + rho(2:2*npoints:2)
    ! Convert npoints from integer(8) to integer(4), that was used by LibXC
    libxc_int = int(npoints)
    ! returned values
    energy   (1:1*npoints) = 0.0_fp
    dedrho   (1:2*npoints) = 0.0_fp
    dedsigma (1:3*npoints) = 0.0_fp
    dedtau   (1:2*npoints) = 0.0_fp
    dedlapl  (1:2*npoints) = 0.0_fp

    if (.not.allocated(this%functionals_list)) return

    ! The sum by arrays must be from LDAs to GGAs to meta-GGAs for avoiding errors
    ! when first meta-GGA functional is calculated, and then LDA functional also is calculated,
    ! which leads to double counting of some derivatives of mGGA
    do i = 1, size(this%functionals_list)
      coefficient = this%coefficients(i)
      select case (xc_f03_func_info_get_family(this%functionals_info(i)))
        case(XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
          call xc_f03_lda_exc_vxc(this%functionals_list(i), libxc_int, rho, tmp_energy, tmp_dedrho)
          energy   (1:1*npoints) = energy   (1:1*npoints) + tmp_energy   * coefficient
          dedrho   (1:2*npoints) = dedrho   (1:2*npoints) + tmp_dedrho   * coefficient
        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc_vxc(this%functionals_list(i), libxc_int, rho, sigma, &
            tmp_energy, tmp_dedrho, tmp_dedsigma)
          energy   (1:1*npoints) = energy   (1:1*npoints) + tmp_energy   * coefficient
          dedrho   (1:2*npoints) = dedrho   (1:2*npoints) + tmp_dedrho   * coefficient
          dedsigma (1:3*npoints) = dedsigma (1:3*npoints) + tmp_dedsigma * coefficient
        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc_vxc(this%functionals_list(i), libxc_int, rho, sigma, lapl, tau, &
            tmp_energy, tmp_dedrho, tmp_dedsigma, tmp_dedlapl, tmp_dedtau)
          energy   (1:1*npoints) = energy   (1:1*npoints) + tmp_energy   * coefficient
          dedrho   (1:2*npoints) = dedrho   (1:2*npoints) + tmp_dedrho   * coefficient
          dedsigma (1:3*npoints) = dedsigma (1:3*npoints) + tmp_dedsigma * coefficient
          dedtau   (1:2*npoints) = dedtau   (1:2*npoints) + tmp_dedtau   * coefficient
          dedlapl  (1:2*npoints) = dedlapl  (1:2*npoints) + tmp_dedlapl  * coefficient
        case default
          call write_error(FIRST_ERROR, this%functionals_info(i))
      end select
    end do
    ! LibXC returns density of energy per particle
    energy(1:npoints) = energy(1:npoints) * rhosum
  end subroutine calc_evxc
  !> @brief  Perform DFT energy and its first and second derivatives calculation
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019
  !> @params NPoints     - (in)  number of points
  !> @params rho         - (in)  density at point
  !> @params sigma       - (in)  gradient of density at point
  !> @params tau         - (in)  local kinetic energy at point
  !> @params lapl        - (in)  laplacian of density at point
  !> @params energy      - (out) energy at point
  !> @params dedrho      - (out) first derivatives energy by density
  !> @params dedsigma    - (out) first derivatives energy by normed gradient
  !> @params dedtau      - (out) first derivatives energy by kinetic energy
  !> @params dedlapl     - (out) first derivatives energy by laplacian
  !> @params v2rho2      - (out) second derivatives by density and density
  !> @params v2sigma2    - (out) second derivatives by normed gradient and normed gradient
  !> @params v2tau2      - (out) second derivatives by kinetic energy and kinetic energy
  !> @params v2lapl2     - (out) second derivatives by laplacian and laplacian
  !> @params v2rhosigma  - (out) second derivatives by density and normed gradient
  !> @params v2rhotau    - (out) second derivatives by density and kinetic energy
  !> @params v2rholapl   - (out) second derivatives by density and laplacian
  !> @params v2sigmatau  - (out) second derivatives by normed gradient and kinetic energy
  !> @params v2sigmalapl - (out) second derivatives by normed gradient and laplacian
  !> @params v2lapltau   - (out) second derivatives by laplacian and kinetic energy
  subroutine calc_evfxc(this, NPoints, &
                        rho, sigma, tau, lapl, energy, dedrho, dedsigma, dedtau, dedlapl, &
                        v2rho2, v2sigma2, v2tau2, v2lapl2, &
                        v2rhosigma, v2rhotau, v2rholapl, v2sigmatau, v2sigmalapl, v2lapltau)
    class(functional_t), intent(inout) :: this
    integer, intent(in) :: NPoints
    real(kind=fp), dimension(*), intent(in)  :: rho, sigma, tau, lapl
    real(kind=fp), dimension(*), intent(out) :: energy
    real(kind=fp), dimension(*), intent(out) :: dedrho, dedsigma, dedtau, dedlapl
    real(kind=fp), dimension(*), intent(out) :: v2rho2, v2sigma2, v2tau2, v2lapl2
    real(kind=fp), dimension(*), intent(out) :: v2rhosigma, v2rhotau, v2rholapl, v2sigmatau, v2sigmalapl, v2lapltau
    ! internal variables
    integer :: i
    real(kind=fp), dimension(1*NPoints) :: tmp_energy
    real(kind=fp), dimension(2*NPoints) :: tmp_dedrho, tmp_dedtau, tmp_dedlapl
    real(kind=fp), dimension(3*NPoints) :: tmp_dedsigma, tmp_v2rho2, tmp_v2tau2, tmp_v2lapl2
    real(kind=fp), dimension(4*NPoints) :: tmp_v2rhotau, tmp_v2rholapl, tmp_v2lapltau
    real(kind=fp), dimension(6*NPoints) :: tmp_v2sigma2, tmp_v2rhosigma, tmp_v2sigmatau, tmp_v2sigmalapl
    ! count of point for LibXC
    integer(C_SIZE_T)                   :: libxc_int
    ! coefficient of functional
    real(kind=fp) :: coefficient
    ! saving total density of each point
    real(kind=fp), dimension(NPoints)   :: rhosum
    ! Build array for multiplying energy at each point
    rhosum = rho(1:2*npoints:2) + rho(2:2*npoints:2)
    ! Convert npoints from integer(8) to integer(4), that was used by LibXC
    libxc_int = int(npoints)
    energy         (1:1*npoints)  = 0.0_fp
    dedrho         (1:2*npoints)  = 0.0_fp
    dedsigma       (1:3*npoints)  = 0.0_fp
    dedlapl        (1:2*npoints)  = 0.0_fp
    dedtau         (1:2*npoints)  = 0.0_fp
    v2rho2         (1:3*npoints)  = 0.0_fp
    v2rhosigma     (1:6*npoints)  = 0.0_fp
    v2rholapl      (1:4*npoints)  = 0.0_fp
    v2rhotau       (1:4*npoints)  = 0.0_fp
    v2sigma2       (1:6*npoints)  = 0.0_fp
    v2sigmalapl    (1:6*npoints)  = 0.0_fp
    v2sigmatau     (1:6*npoints)  = 0.0_fp
    v2lapl2        (1:3*npoints)  = 0.0_fp
    v2lapltau      (1:4*npoints)  = 0.0_fp
    v2tau2         (1:3*npoints)  = 0.0_fp

    if (.not.allocated(this%functionals_list)) return

    ! The sum by arrays must be in "select case" block for avoiding errors
    ! when firstly mGGA functional is calculated, and then LDA functional also is calculated,
    ! that leads to double counting of some derivatives of mGGA
    do i = 1, size(this%functionals_list)
      coefficient = this%coefficients(i)
      select case (xc_f03_func_info_get_family(this%functionals_info(i)))
        case(XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
          call xc_f03_lda_exc_vxc_fxc(this%functionals_list(i), libxc_int, rho, &
            tmp_energy, tmp_dedrho, tmp_v2rho2)
          energy        (1:1*npoints) = energy (1:1*npoints)         + tmp_energy         * coefficient
          dedrho        (1:2*npoints) = dedrho (1:2*npoints)         + tmp_dedrho         * coefficient
          v2rho2        (1:3*npoints) = v2rho2 (1:3*npoints)         + tmp_v2rho2         * coefficient
        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc_vxc_fxc(this%functionals_list(i), libxc_int, rho, sigma, &
            tmp_energy, tmp_dedrho, tmp_dedsigma, tmp_v2rho2, tmp_v2rhosigma, tmp_v2sigma2)
          energy        (1: 1*npoints) = energy      (1: 1*npoints)   + tmp_energy         * coefficient
          dedrho        (1: 2*npoints) = dedrho      (1: 2*npoints)   + tmp_dedrho         * coefficient
          dedsigma      (1: 3*npoints) = dedsigma    (1: 3*npoints)   + tmp_dedsigma       * coefficient
          v2rho2        (1: 3*npoints) = v2rho2      (1: 3*npoints)   + tmp_v2rho2         * coefficient
          v2rhosigma    (1: 6*npoints) = v2rhosigma  (1: 6*npoints)   + tmp_v2rhosigma     * coefficient
          v2sigma2      (1: 6*npoints) = v2sigma2    (1: 6*npoints)   + tmp_v2sigma2       * coefficient
        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc_vxc_fxc(this%functionals_list(i), libxc_int, rho, sigma, lapl, tau, &
            tmp_energy, tmp_dedrho, tmp_dedsigma, tmp_dedlapl, tmp_dedtau, &
            tmp_v2rho2, tmp_v2rhosigma, tmp_v2rholapl, tmp_v2rhotau, &
            tmp_v2sigma2, tmp_v2sigmalapl, tmp_v2sigmatau, tmp_v2lapl2, tmp_v2lapltau, tmp_v2tau2)
          energy        (1:1*npoints)   = energy        (1:1*npoints)   + tmp_energy         * coefficient
          dedrho        (1:2*npoints)   = dedrho        (1:2*npoints)   + tmp_dedrho         * coefficient
          dedsigma      (1:3*npoints)   = dedsigma      (1:3*npoints)   + tmp_dedsigma       * coefficient
          dedlapl       (1:2*npoints)   = dedlapl       (1:2*npoints)   + tmp_dedlapl        * coefficient
          dedtau        (1:2*npoints)   = dedtau        (1:2*npoints)   + tmp_dedtau         * coefficient
          v2rho2        (1:3*npoints)   = v2rho2        (1:3*npoints)   + tmp_v2rho2         * coefficient
          v2rhosigma    (1:6*npoints)   = v2rhosigma    (1:6*npoints)   + tmp_v2rhosigma     * coefficient
          v2rholapl     (1:4*npoints)   = v2rholapl     (1:4*npoints)   + tmp_v2rholapl      * coefficient
          v2rhotau      (1:4*npoints)   = v2rhotau      (1:4*npoints)   + tmp_v2rhotau       * coefficient
          v2sigma2      (1:6*npoints)   = v2sigma2      (1:6*npoints)   + tmp_v2sigma2       * coefficient
          v2sigmalapl   (1:6*npoints)   = v2sigmalapl   (1:6*npoints)   + tmp_v2sigmalapl    * coefficient
          v2sigmatau    (1:6*npoints)   = v2sigmatau    (1:6*npoints)   + tmp_v2sigmatau     * coefficient
          v2lapl2       (1:3*npoints)   = v2lapl2       (1:3*npoints)   + tmp_v2lapl2        * coefficient
          v2lapltau     (1:4*npoints)   = v2lapltau     (1:4*npoints)   + tmp_v2lapltau      * coefficient
          v2tau2        (1:3*npoints)   = v2tau2        (1:3*npoints)   + tmp_v2tau2         * coefficient
        case default
          call write_error(SECOND_ERROR, this%functionals_info(i))
      end select
    end do
    ! LibXC returns density of energy per particle
    energy(1:npoints) = energy(1:npoints) * rhosum
  end subroutine calc_evfxc
  !> @brief  Perform DFT energy and its first, second and third derivatives calculation
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019
  !> @params npoints        - (in)  number of points
  !> @params rho            - (in)  density at point
  !> @params sigma          - (in)  gradient of density at point
  !> @params tau            - (in)  local kinetic energy at point
  !> @params lapl           - (in)  laplacian of density at point
  !> @params energy         - (out) energy at point
  !> @params dedrho         - (out) first  derivatives energy by density
  !> @params dedsigma       - (out) first  derivatives energy by normed gradient
  !> @params dedlapl        - (out) first  derivatives energy by laplacian
  !> @params dedtau         - (out) first  derivatives energy by kinetic energy
  !> @params v2rho2         - (out) second derivatives by density         and density
  !> @params v2rhosigma     - (out) second derivatives by density         and normed gradient
  !> @params v2rholapl      - (out) second derivatives by density         and laplacian
  !> @params v2rhotau       - (out) second derivatives by density         and kinetic energy
  !> @params v2sigma2       - (out) second derivatives by normed gradient and normed gradient
  !> @params v2sigmalapl    - (out) second derivatives by normed gradient and laplacian
  !> @params v2sigmatau     - (out) second derivatives by normed gradient and kinetic energy
  !> @params v2lapl2        - (out) second derivatives by laplacian       and laplacian
  !> @params v2lapltau      - (out) second derivatives by laplacian       and kinetic energy
  !> @params v2tau2         - (out) second derivatives by kinetic energy  and kinetic energy
  !> @params v3rho3         - (out) third  derivatives by density         and density         and density
  !> @params v3rho2sigma    - (out) third  derivatives by density         and density         and normed gradient
  !> @params v3rho2lapl     - (out) third  derivatives by density         and density         and laplacian
  !> @params v3rho2tau      - (out) third  derivatives by density         and density         and kinetic energy
  !> @params v3rhosigma2    - (out) third  derivatives by density         and normed gradient and normed gradient
  !> @params v3rhosigmalapl - (out) third  derivatives by density         and normed gradient and laplacian
  !> @params v3rhosigmatau  - (out) third  derivatives by density         and normed gradient and density
  !> @params v3rholapl2     - (out) third  derivatives by density         and laplacian       and laplacian
  !> @params v3rholapltau   - (out) third  derivatives by density         and laplacian       and density
  !> @params v3rhotau2      - (out) third  derivatives by density         and kinetic energy  and kinetic energy
  !> @params v3sigma3       - (out) third  derivatives by normed gradient and normed gradient and normed gradient
  !> @params v3sigma2lapl   - (out) third  derivatives by normed gradient and normed gradient and laplacian
  !> @params v3sigma2tau    - (out) third  derivatives by normed gradient and normed gradient and kinetic energy
  !> @params v3sigmalapl2   - (out) third  derivatives by normed gradient and laplacian       and laplacian
  !> @params v3sigmalapltau - (out) third  derivatives by normed gradient and laplacian       and kinetic energy
  !> @params v3sigmatau2    - (out) third  derivatives by normed gradient and kinetic energy  and kinetic energy
  !> @params v3lapl3        - (out) third  derivatives by laplacian       and laplacian       and laplacian
  !> @params v3lapl2tau     - (out) third  derivatives by laplacian       and laplacian       and kinetic energy
  !> @params v3lapltau2     - (out) third  derivatives by laplacian       and kinetic energy  and kinetic energy
  !> @params v3tau3         - (out) third  derivatives by kinetic energy  and kinetic energy  and kinetic energy
  subroutine calc_xc(this, npoints, &
                     rho, sigma, tau, lapl, energy, dedrho, dedsigma, dedlapl, dedtau, &
                     v2rho2, v2rhosigma, v2rholapl, v2rhotau, &
                     v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2, &
                     v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
                     v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
                     v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
                     v3lapltau2, v3tau3)
    class(functional_t), intent(inout) :: this
    integer, intent(in) :: NPoints
    real(kind=fp), dimension(*), intent(in)  :: rho, sigma, tau, lapl
    real(kind=fp), dimension(*), intent(out) :: energy
    real(kind=fp), dimension(*), intent(out) :: dedrho, dedsigma, dedlapl, dedtau
    real(kind=fp), dimension(*), intent(out) :: v2rho2, v2rhosigma, v2rholapl, v2rhotau
    real(kind=fp), dimension(*), intent(out) :: v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2
    real(kind=fp), dimension(*), intent(out) :: v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl
    real(kind=fp), dimension(*), intent(out) :: v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl
    real(kind=fp), dimension(*), intent(out) :: v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau
    real(kind=fp), dimension(*), intent(out) :: v3lapltau2, v3tau3
    ! internal variables
    integer :: i
    real(kind=fp), dimension(1*NPoints)  :: tmp_energy
    real(kind=fp), dimension(2*NPoints)  :: tmp_dedrho, tmp_dedtau, tmp_dedlapl
    real(kind=fp), dimension(3*NPoints)  :: tmp_dedsigma, tmp_v2rho2, tmp_v2tau2, tmp_v2lapl2
    real(kind=fp), dimension(4*NPoints)  :: tmp_v2rhotau, tmp_v2rholapl, tmp_v2lapltau
    real(kind=fp), dimension(6*NPoints)  :: tmp_v2sigma2, tmp_v2rhosigma, tmp_v2sigmatau, tmp_v2sigmalapl
    !for third derivatives
    real(kind=fp), dimension(4*NPoints)  :: tmp_v3rho3, tmp_v3tau3, tmp_v3lapl3
    real(kind=fp), dimension(6*NPoints)  :: tmp_v3rho2tau, tmp_v3rhotau2, tmp_v3rho2lapl, tmp_v3rholapl2, tmp_v3lapl2tau, &
                                            tmp_v3lapltau2
    real(kind=fp), dimension(8*NPoints)  :: tmp_v3rholapltau
    real(kind=fp), dimension(9*NPoints)  :: tmp_v3rho2sigma, tmp_v3sigmatau2, tmp_v3sigmalapl2
    real(kind=fp), dimension(10*NPoints) :: tmp_v3sigma3
    real(kind=fp), dimension(12*NPoints) :: tmp_v3rhosigma2, tmp_v3rhosigmatau, tmp_v3sigma2tau, tmp_v3rhosigmalapl, &
                                            tmp_v3sigma2lapl, tmp_v3sigmalapltau
    ! count of point for LibXC
    integer(C_SIZE_T)                    :: libxc_int
    ! coefficient of functional
    real(kind=fp) :: coefficient
    ! saving total density of each point
    real(kind=fp), dimension(NPoints)    :: rhosum
    ! Build array for multiplying energy at each point
    rhosum = rho(1:2*npoints:2) + rho(2:2*npoints:2)
    ! Convert npoints from integer(8) to integer(4), that was used by LibXC
    libxc_int = int(npoints)
    energy         (1:1*npoints)  = 0.0_fp
    dedrho         (1:2*npoints)  = 0.0_fp
    dedsigma       (1:3*npoints)  = 0.0_fp
    dedlapl        (1:2*npoints)  = 0.0_fp
    dedtau         (1:2*npoints)  = 0.0_fp
    v2rho2         (1:3*npoints)  = 0.0_fp
    v2rhosigma     (1:6*npoints)  = 0.0_fp
    v2rholapl      (1:4*npoints)  = 0.0_fp
    v2rhotau       (1:4*npoints)  = 0.0_fp
    v2sigma2       (1:6*npoints)  = 0.0_fp
    v2sigmalapl    (1:6*npoints)  = 0.0_fp
    v2sigmatau     (1:6*npoints)  = 0.0_fp
    v2lapl2        (1:3*npoints)  = 0.0_fp
    v2lapltau      (1:4*npoints)  = 0.0_fp
    v2tau2         (1:3*npoints)  = 0.0_fp
    v3rho3         (1:4*npoints)  = 0.0_fp
    v3rho2sigma    (1:9*npoints)  = 0.0_fp
    v3rho2lapl     (1:6*npoints)  = 0.0_fp
    v3rho2tau      (1:6*npoints)  = 0.0_fp
    v3rhosigma2    (1:12*npoints) = 0.0_fp
    v3rhosigmalapl (1:12*npoints) = 0.0_fp
    v3rhosigmatau  (1:12*npoints) = 0.0_fp
    v3rholapl2     (1:6*npoints)  = 0.0_fp
    v3rholapltau   (1:8*npoints)  = 0.0_fp
    v3rhotau2      (1:6*npoints)  = 0.0_fp
    v3sigma3       (1:10*npoints) = 0.0_fp
    v3sigma2lapl   (1:12*npoints) = 0.0_fp
    v3sigma2tau    (1:12*npoints) = 0.0_fp
    v3sigmalapl2   (1:9*npoints)  = 0.0_fp
    v3sigmalapltau (1:12*npoints) = 0.0_fp
    v3sigmatau2    (1:9*npoints)  = 0.0_fp
    v3lapl3        (1:4*npoints)  = 0.0_fp
    v3lapl2tau     (1:6*npoints)  = 0.0_fp
    v3lapltau2     (1:6*npoints)  = 0.0_fp
    v3tau3         (1:4*npoints)  = 0.0_fp

    if (.not.allocated(this%functionals_list)) return

    ! The sum by arrays must be in "select case" block for avoiding errors
    ! when firstly mGGA functional is calculated, and then LDA functional also is calculated,
    ! that leads to double counting of some derivatives of mGGA
    do i = 1, size(this%functionals_list)
      coefficient = this%coefficients(i)
      select case (xc_f03_func_info_get_family(this%functionals_info(i)))
        case(XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
          call xc_f03_lda_exc_vxc_fxc_kxc(this%functionals_list(i), libxc_int, rho, &
            tmp_energy, tmp_dedrho, tmp_v2rho2, tmp_v3rho3)
          energy        (1:1*npoints) = energy (1:1*npoints)         + tmp_energy         * coefficient
          dedrho        (1:2*npoints) = dedrho (1:2*npoints)         + tmp_dedrho         * coefficient
          v2rho2        (1:3*npoints) = v2rho2 (1:3*npoints)         + tmp_v2rho2         * coefficient
          v3rho3        (1:4*npoints) = v3rho3 (1:4*npoints)         + tmp_v3rho3         * coefficient
        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc_vxc_fxc_kxc(this%functionals_list(i), libxc_int, rho, sigma,  &
            tmp_energy, tmp_dedrho, tmp_dedsigma, tmp_v2rho2, tmp_v2rhosigma, tmp_v2sigma2, &
            tmp_v3rho3, tmp_v3rho2sigma, tmp_v3rhosigma2, tmp_v3sigma3)
          energy        (1: 1*npoints) = energy      (1: 1*npoints)   + tmp_energy         * coefficient
          dedrho        (1: 2*npoints) = dedrho      (1: 2*npoints)   + tmp_dedrho         * coefficient
          dedsigma      (1: 3*npoints) = dedsigma    (1: 3*npoints)   + tmp_dedsigma       * coefficient
          v2rho2        (1: 3*npoints) = v2rho2      (1: 3*npoints)   + tmp_v2rho2         * coefficient
          v2rhosigma    (1: 6*npoints) = v2rhosigma  (1: 6*npoints)   + tmp_v2rhosigma     * coefficient
          v2sigma2      (1: 6*npoints) = v2sigma2    (1: 6*npoints)   + tmp_v2sigma2       * coefficient
          v3rho3        (1: 4*npoints) = v3rho3      (1: 4*npoints)   + tmp_v3rho3         * coefficient
          v3rho2sigma   (1: 9*npoints) = v3rho2sigma (1: 9*npoints)   + tmp_v3rho2sigma    * coefficient
          v3rhosigma2   (1:12*npoints) = v3rhosigma2 (1:12*npoints)   + tmp_v3rhosigma2    * coefficient
          v3sigma3      (1:10*npoints) = v3sigma3    (1:10*npoints)   + tmp_v3sigma3       * coefficient
        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc_vxc_fxc_kxc(this%functionals_list(i), libxc_int, rho, sigma, lapl, tau, &
            tmp_energy, tmp_dedrho, tmp_dedsigma, tmp_dedlapl, tmp_dedtau, tmp_v2rho2, tmp_v2rhosigma, tmp_v2rholapl, &
            tmp_v2rhotau, tmp_v2sigma2, tmp_v2sigmalapl, tmp_v2sigmatau, tmp_v2lapl2, tmp_v2lapltau, tmp_v2tau2,   &
            tmp_v3rho3, tmp_v3rho2sigma, tmp_v3rho2lapl, tmp_v3rho2tau, tmp_v3rhosigma2, tmp_v3rhosigmalapl,       &
            tmp_v3rhosigmatau, tmp_v3rholapl2, tmp_v3rholapltau, tmp_v3rhotau2, tmp_v3sigma3, tmp_v3sigma2lapl,    &
            tmp_v3sigma2tau, tmp_v3sigmalapl2, tmp_v3sigmalapltau, tmp_v3sigmatau2, tmp_v3lapl3, tmp_v3lapl2tau,   &
            tmp_v3lapltau2, tmp_v3tau3)
          energy        (1:1*npoints)   = energy        (1:1*npoints)   + tmp_energy         * coefficient
          dedrho        (1:2*npoints)   = dedrho        (1:2*npoints)   + tmp_dedrho         * coefficient
          dedsigma      (1:3*npoints)   = dedsigma      (1:3*npoints)   + tmp_dedsigma       * coefficient
          dedlapl       (1:2*npoints)   = dedlapl       (1:2*npoints)   + tmp_dedlapl        * coefficient
          dedtau        (1:2*npoints)   = dedtau        (1:2*npoints)   + tmp_dedtau         * coefficient
          v2rho2        (1:3*npoints)   = v2rho2        (1:3*npoints)   + tmp_v2rho2         * coefficient
          v2rhosigma    (1:6*npoints)   = v2rhosigma    (1:6*npoints)   + tmp_v2rhosigma     * coefficient
          v2rholapl     (1:4*npoints)   = v2rholapl     (1:4*npoints)   + tmp_v2rholapl      * coefficient
          v2rhotau      (1:4*npoints)   = v2rhotau      (1:4*npoints)   + tmp_v2rhotau       * coefficient
          v2sigma2      (1:6*npoints)   = v2sigma2      (1:6*npoints)   + tmp_v2sigma2       * coefficient
          v2sigmalapl   (1:6*npoints)   = v2sigmalapl   (1:6*npoints)   + tmp_v2sigmalapl    * coefficient
          v2sigmatau    (1:6*npoints)   = v2sigmatau    (1:6*npoints)   + tmp_v2sigmatau     * coefficient
          v2lapl2       (1:3*npoints)   = v2lapl2       (1:3*npoints)   + tmp_v2lapl2        * coefficient
          v2lapltau     (1:4*npoints)   = v2lapltau     (1:4*npoints)   + tmp_v2lapltau      * coefficient
          v2tau2        (1:3*npoints)   = v2tau2        (1:3*npoints)   + tmp_v2tau2         * coefficient
          v3rho3        (1:4*npoints)   = v3rho3        (1:4*npoints)   + tmp_v3rho3         * coefficient
          v3rho2sigma   (1:9*npoints)   = v3rho2sigma   (1:9*npoints)   + tmp_v3rho2sigma    * coefficient
          v3rho2lapl    (1:6*npoints)   = v3rho2lapl    (1:6*npoints)   + tmp_v3rho2lapl     * coefficient
          v3rho2tau     (1:6*npoints)   = v3rho2tau     (1:6*npoints)   + tmp_v3rho2tau      * coefficient
          v3rhosigma2   (1:12*npoints)  = v3rhosigma2   (1:12*npoints)  + tmp_v3rhosigma2    * coefficient
          v3rhosigmalapl(1:12*npoints)  = v3rhosigmalapl(1:12*npoints)  + tmp_v3rhosigmalapl * coefficient
          v3rhosigmatau (1:12*npoints)  = v3rhosigmatau (1:12*npoints)  + tmp_v3rhosigmatau  * coefficient
          v3rholapl2    (1:6*npoints)   = v3rholapl2    (1:6*npoints)   + tmp_v3rholapl2     * coefficient
          v3rholapltau  (1:8*npoints)   = v3rholapltau  (1:8*npoints)   + tmp_v3rholapltau   * coefficient
          v3rhotau2     (1:6*npoints)   = v3rhotau2     (1:6*npoints)   + tmp_v3rhotau2      * coefficient
          v3sigma3      (1:10*npoints)  = v3sigma3      (1:10*npoints)  + tmp_v3sigma3       * coefficient
          v3sigma2lapl  (1:12*npoints)  = v3sigma2lapl  (1:12*npoints)  + tmp_v3sigma2lapl   * coefficient
          v3sigma2tau   (1:12*npoints)  = v3sigma2tau   (1:12*npoints)  + tmp_v3sigma2tau    * coefficient
          v3sigmalapl2  (1:9*npoints)   = v3sigmalapl2  (1:9*npoints)   + tmp_v3sigmalapl2   * coefficient
          v3sigmalapltau(1:12*npoints)  = v3sigmalapltau(1:12*npoints)  + tmp_v3sigmalapltau * coefficient
          v3sigmatau2   (1:9*npoints)   = v3sigmatau2   (1:9*npoints)   + tmp_v3sigmatau2    * coefficient
          v3lapl3       (1:4*npoints)   = v3lapl3       (1:4*npoints)   + tmp_v3lapl3        * coefficient
          v3lapl2tau    (1:6*npoints)   = v3lapl2tau    (1:6*npoints)   + tmp_v3lapl2tau     * coefficient
          v3lapltau2    (1:6*npoints)   = v3lapltau2    (1:6*npoints)   + tmp_v3lapltau2     * coefficient
          v3tau3        (1:4*npoints)   = v3tau3        (1:4*npoints)   + tmp_v3tau3         * coefficient
        case default
          call write_error(THIRD_ERROR, this%functionals_info(i))
      end select
    end do
    ! LibXC returns density of energy per particle
    energy(1:npoints) = energy(1:npoints) * rhosum
  end subroutine calc_xc
  !> @brief  Internal procedure for writing error
  !> @author Igor S. Gerasimov
  !> @date   Oct, 2019 - Initial release -
  !> @date   Jul, 2021 Using messages module
  !> @params error_code      - (in)  code of error
  !> @params functional_info - (in)  info about functional for getting of name
  subroutine write_error(error_code, functional_info)
    use messages, only: WITH_ABORT, show_message
    integer,                  intent(in) :: error_code
    type(xc_f03_func_info_t), intent(in) :: functional_info
    character(len=:), allocatable :: functional_name
    character(len=:), allocatable :: error_line
    select case(xc_f03_func_info_get_kind(functional_info))
      case(XC_EXCHANGE)
        functional_name = trim(xc_f03_func_info_get_name(functional_info)) // " exchange functional."
      case(XC_CORRELATION)
        functional_name = trim(xc_f03_func_info_get_name(functional_info)) // " correlation functional."
      case(XC_EXCHANGE_CORRELATION)
        functional_name = trim(xc_f03_func_info_get_name(functional_info)) // " exchange-correlation functional."
      case(XC_KINETIC)
        functional_name = trim(xc_f03_func_info_get_name(functional_info)) // " kinetic functional."
      case default
        functional_name = "Unnamed functional"
      end select
    select case(error_code)
      case(ENERGY_ERROR)
        error_line = "Something went wrong while the energy was tried to calculate using " // functional_name
      case(FIRST_ERROR)
        error_line = "Something went wrong while the first derivatives were tried to calculate using " // functional_name
      case(SECOND_ERROR)
        error_line = "Something went wrong while the second derivatives were tried to calculate using " // functional_name
      case(THIRD_ERROR)
        error_line = "Something went wrong while the third derivatives were tried to calculate using " // functional_name
      case(MGGA_3RD_ERROR)
        error_line = "TD-DFT third derivatives do not support for meta-GGA functionals like " // functional_name
      case default
        error_line = "Explore OQP error"
    end select
    call show_message(error_line)
    call show_message("Abort was produced by LibXC interface...", WITH_ABORT)
  end subroutine write_error
end module functionals
