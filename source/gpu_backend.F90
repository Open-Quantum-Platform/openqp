module gpu_backend

  use iso_c_binding, only: c_bool, c_char, c_int, c_null_char

  implicit none

contains

  function gpu_backend_available() result(available)
    logical :: available
#ifdef OQP_CUDA_ENABLE
    available = .true.
#else
    available = .false.
#endif
  end function gpu_backend_available

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

    ! The first METC GPU milestone only records a stable ABI boundary.  The
    ! CUDA implementation will attach device allocation and cuBLAS handles here.
    if (enabled .and. device < 0_c_int) error stop "GPU device id must be non-negative"
  end subroutine gpu_backend_configure

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
