module gpu_backend

  use iso_c_binding, only: c_bool, c_char, c_double, c_int, c_null_char, c_ptr, c_loc

  implicit none

  logical :: gpu_metc_requested = .false.
  integer(c_int) :: gpu_device_id = 0_c_int

  interface
#ifdef OQP_CUDA_ENABLE
    function oqp_gpu_metc_contract(ids, ints, ncur, f3, d3, nf, nmatrix, nbf, &
                                   cur_pass, scale_exchange, scale_coulomb, is_umrsf) &
                                   bind(C, name="oqp_gpu_metc_contract") result(ierr)
      import :: c_bool, c_double, c_int, c_ptr
      type(c_ptr), value :: ids
      type(c_ptr), value :: ints
      integer(c_int), value :: ncur
      type(c_ptr), value :: f3
      type(c_ptr), value :: d3
      integer(c_int), value :: nf
      integer(c_int), value :: nmatrix
      integer(c_int), value :: nbf
      integer(c_int), value :: cur_pass
      real(c_double), value :: scale_exchange
      real(c_double), value :: scale_coulomb
      logical(c_bool), value :: is_umrsf
      integer(c_int) :: ierr
    end function oqp_gpu_metc_contract
#endif
  end interface

contains

  function gpu_backend_available() result(available)
    logical :: available
#ifdef OQP_CUDA_ENABLE
    available = .true.
#else
    available = .false.
#endif
  end function gpu_backend_available

  function gpu_backend_metc_enabled() result(enabled)
    logical :: enabled
    character(len=16) :: env_value
    integer :: status

    enabled = .false.
#ifdef OQP_CUDA_ENABLE
    enabled = gpu_metc_requested
    call get_environment_variable("OQP_GPU_METC", env_value, status=status)
    if (status == 0) then
      select case (trim(adjustl(env_value)))
      case ("1", "true", "TRUE", "on", "ON", "yes", "YES")
        enabled = .true.
      case ("0", "false", "FALSE", "off", "OFF", "no", "NO")
        enabled = .false.
      end select
    end if
#endif
  end function gpu_backend_metc_enabled

  subroutine gpu_backend_describe(buffer, buffer_len) bind(C, name="oqp_gpu_backend_describe")
    character(kind=c_char), intent(out) :: buffer(*)
    integer(c_int), value, intent(in) :: buffer_len

#ifdef OQP_CUDA_ENABLE
    call copy_c_string("cuda", buffer, buffer_len)
#else
    call copy_c_string("cpu-fallback", buffer, buffer_len)
#endif
  end subroutine gpu_backend_describe

  subroutine gpu_backend_configure(enabled, device) bind(C, name="oqp_gpu_backend_configure")
    logical(c_bool), value, intent(in) :: enabled
    integer(c_int), value, intent(in) :: device

    if (enabled .and. device < 0_c_int) error stop "GPU device id must be non-negative"
    gpu_metc_requested = enabled
    gpu_device_id = device
  end subroutine gpu_backend_configure

  subroutine gpu_backend_metc_contract(ids, ints, ncur, f3, d3, nf, nmatrix, nbf, &
                                       cur_pass, scale_exchange, scale_coulomb, &
                                       is_umrsf, ierr)
    integer, intent(in) :: ids(:,:)
    real(c_double), intent(in), target, contiguous :: ints(:)
    integer, intent(in) :: ncur, nf, nmatrix, nbf, cur_pass
    real(c_double), intent(inout), target, contiguous :: f3(:,:,:,:)
    real(c_double), intent(in), target, contiguous :: d3(:,:,:,:)
    real(c_double), intent(in) :: scale_exchange, scale_coulomb
    logical, intent(in) :: is_umrsf
    integer, intent(out) :: ierr

    integer(c_int), allocatable, target :: ids32(:,:)
    integer :: n
    logical(c_bool) :: c_is_umrsf

    ierr = 1
#ifdef OQP_CUDA_ENABLE
    if (ncur <= 0) then
      ierr = 0
      return
    end if

    allocate(ids32(4, ncur))
    do n = 1, ncur
      ids32(1:4, n) = int(ids(1:4, n), c_int)
    end do
    c_is_umrsf = is_umrsf

    ierr = int(oqp_gpu_metc_contract(c_loc(ids32(1,1)), c_loc(ints(1)), int(ncur, c_int), &
                                     c_loc(f3(1,1,1,1)), c_loc(d3(1,1,1,1)), &
                                     int(nf, c_int), int(nmatrix, c_int), int(nbf, c_int), &
                                     int(cur_pass, c_int), scale_exchange, scale_coulomb, c_is_umrsf))
    deallocate(ids32)
#endif
  end subroutine gpu_backend_metc_contract

  subroutine copy_c_string(text, buffer, buffer_len)
    character(len=*), intent(in) :: text
    character(kind=c_char), intent(out) :: buffer(*)
    integer(c_int), value, intent(in) :: buffer_len
    integer :: i, ncopy

    if (buffer_len <= 0) return

    ncopy = min(len_trim(text), int(buffer_len) - 1)
    do i = 1, ncopy
      buffer(i) = text(i:i)
    end do
    buffer(ncopy + 1) = c_null_char
  end subroutine copy_c_string

end module gpu_backend
