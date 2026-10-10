module gpu_xc_response_backend
  !! Minimal Fortran-side ABI stub for the experimental CUDA TDHF/TDDFT
  !! XC-response branch.  The routines intentionally remain opt-in behind
  !! OQP_CUDA_ENABLE so CPU-only builds keep the normal response path.
  use iso_c_binding, only: c_double, c_int
  implicit none
  private

  public :: gpu_xc_response_enabled
  public :: gpu_xc_response_describe
  public :: gpu_xc_response_contract
  public :: gpu_xc_response_status_message
  public :: gpu_xc_response_preflight_status
  public :: GPU_XC_RESPONSE_STATUS_INVALID_INPUT
  public :: GPU_XC_RESPONSE_STATUS_OVERFLOW

  integer(c_int), parameter :: GPU_XC_RESPONSE_STATUS_DISABLED = 1_c_int
  integer(c_int), parameter :: GPU_XC_RESPONSE_STATUS_INVALID_INPUT = 2_c_int
  integer(c_int), parameter :: GPU_XC_RESPONSE_STATUS_OVERFLOW = 3_c_int

#ifdef OQP_CUDA_ENABLE
  interface
    integer(c_int) function oqp_gpu_xc_response_contract(nbasis, nstate, density, kernel, response) bind(C, name="oqp_gpu_xc_response_contract")
      import :: c_double, c_int
      integer(c_int), value :: nbasis
      integer(c_int), value :: nstate
      real(c_double), intent(in) :: density(*)
      real(c_double), intent(in) :: kernel(*)
      real(c_double), intent(inout) :: response(*)
    end function oqp_gpu_xc_response_contract
  end interface
#else
  ! Keep the C ABI symbol name visible in source-level tests even when CUDA is
  ! disabled for ordinary CPU-only builds: oqp_gpu_xc_response_contract.
#endif

contains

  logical function gpu_xc_response_enabled()
#ifdef OQP_CUDA_ENABLE
    gpu_xc_response_enabled = .true.
#else
    gpu_xc_response_enabled = .false.
#endif
  end function gpu_xc_response_enabled

  subroutine gpu_xc_response_describe(message)
    character(len=*), intent(out) :: message
#ifdef OQP_CUDA_ENABLE
    message = "CUDA XC-response backend enabled"
#else
    message = "CUDA XC-response backend disabled; using CPU fallback"
#endif
  end subroutine gpu_xc_response_describe

  subroutine gpu_xc_response_status_message(status, message)
    integer(c_int), intent(in) :: status
    character(len=*), intent(out) :: message

    select case (status)
    case (0_c_int)
      message = "CUDA XC-response contract completed"
    case (GPU_XC_RESPONSE_STATUS_DISABLED)
      message = "CUDA XC-response backend disabled; using CPU fallback"
    case (GPU_XC_RESPONSE_STATUS_INVALID_INPUT)
      message = "CUDA XC-response invalid input dimensions or null buffer"
    case (GPU_XC_RESPONSE_STATUS_OVERFLOW)
      message = "CUDA XC-response dimension overflow before launch"
    case default
      write(message, '(A,I0)') "CUDA XC-response runtime status ", status
    end select
  end subroutine gpu_xc_response_status_message

  integer(c_int) function gpu_xc_response_preflight_status(nbasis, nstate) result(status)
    integer(c_int), value :: nbasis
    integer(c_int), value :: nstate

    if (nbasis <= 0_c_int .or. nstate <= 0_c_int) then
      status = GPU_XC_RESPONSE_STATUS_INVALID_INPUT
    else if (nbasis > huge(0_c_int) / nstate) then
      status = GPU_XC_RESPONSE_STATUS_OVERFLOW
#ifdef OQP_CUDA_ENABLE
    else
      status = 0_c_int
#else
    else
      status = GPU_XC_RESPONSE_STATUS_DISABLED
#endif
    end if
  end function gpu_xc_response_preflight_status

  integer(c_int) function gpu_xc_response_contract(nbasis, nstate, density, kernel, response) result(status)
    integer(c_int), value :: nbasis
    integer(c_int), value :: nstate
    real(c_double), intent(in) :: density(*)
    real(c_double), intent(in) :: kernel(*)
    real(c_double), intent(inout) :: response(*)
    status = gpu_xc_response_preflight_status(nbasis, nstate)
    if (status /= 0_c_int) return
#ifdef OQP_CUDA_ENABLE
    status = oqp_gpu_xc_response_contract(nbasis, nstate, density, kernel, response)
#else
    status = GPU_XC_RESPONSE_STATUS_DISABLED
#endif
  end function gpu_xc_response_contract

end module gpu_xc_response_backend
