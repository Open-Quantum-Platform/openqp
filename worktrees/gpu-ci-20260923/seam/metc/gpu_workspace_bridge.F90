module gpu_workspace_bridge
  !! Source-side workspace manifest / ABI bridge for the unified GPU workspace
  !! manager (METC-C0).
  !!
  !! This module is the Fortran/CUDA-side counterpart of the Python
  !! ``GpuWorkspaceManager`` (``pyoqp/oqp/utils/gpu_workspace.py``).  It declares
  !! the shared ABI vocabulary -- target namespace codes, residency classes, the
  !! f3 accumulator policy, the float64 element size, and a manifest-row derived
  !! type -- plus validation-only verb stubs.  Its sole purpose at this stage is
  !! to prove that the source side and the Python control plane agree on buffer
  !! identity, byte sizes, offsets, residency, and f3 policy *before* the METC-C1
  !! residency refit borrows resident pointers.
  !!
  !! IMPORTANT (METC-C0 boundary):
  !!   * No CUDA kernels are launched here.
  !!   * No device memory is allocated, copied, or freed here.
  !!   * The verb stubs are validation-only; METC-C1 will back them with real
  !!     workspace residency (d3 -> RESIDENT_INPUT, per-thread f3 -> RESIDENT_ACCUM).
  !!
  !! The integer parameter values below are a hard contract: the parity test
  !! ``tests/test_gpu_workspace_bridge.py`` asserts they equal the Python
  !! ``TARGET_ABI_CODES`` / ``RESIDENCY_ABI_CODES`` / ``F3_POLICY_ABI_CODES``.

  use iso_c_binding, only: c_int, c_int64_t, c_ptr

  implicit none
  private

  ! --- ABI schema + element size (match Python WORKSPACE_* constants) ---------
  integer(c_int), parameter, public :: GPU_WS_MANIFEST_SCHEMA = 1_c_int
  integer(c_int), parameter, public :: GPU_WS_ELEMENT_BYTES   = 8_c_int

  ! --- target namespace codes (match Python TARGET_ABI_CODES) -----------------
  integer(c_int), parameter, public :: GPU_WS_TARGET_METC        = 0_c_int
  integer(c_int), parameter, public :: GPU_WS_TARGET_XC_RESPONSE = 1_c_int

  ! --- residency class codes (match Python RESIDENCY_ABI_CODES) ---------------
  integer(c_int), parameter, public :: GPU_WS_RESIDENCY_HOST_ONLY       = 0_c_int
  integer(c_int), parameter, public :: GPU_WS_RESIDENCY_DEVICE_RESIDENT = 1_c_int
  integer(c_int), parameter, public :: GPU_WS_RESIDENCY_MIRRORED        = 2_c_int
  integer(c_int), parameter, public :: GPU_WS_RESIDENCY_BORROWED        = 3_c_int

  ! --- f3 accumulator policy codes (match Python F3_POLICY_ABI_CODES) ---------
  ! PER_THREAD is the default; SINGLE_ATOMIC must be selected explicitly.  The
  ! bridge never assumes a single global f3 accumulator.
  integer(c_int), parameter, public :: GPU_WS_F3_PER_THREAD    = 0_c_int
  integer(c_int), parameter, public :: GPU_WS_F3_SINGLE_ATOMIC = 1_c_int

  !> One manifest row describing a planned workspace buffer.  Mirrors a single
  !> row of the Python ``workspace_manifest_rows`` output.
  type, public :: gpu_ws_buffer_t
    character(len=32)  :: name          = ''
    integer(c_int)     :: target        = GPU_WS_TARGET_METC
    integer(c_int)     :: residency     = GPU_WS_RESIDENCY_MIRRORED
    integer(c_int64_t) :: nbytes        = 0_c_int64_t
    integer(c_int64_t) :: offset        = 0_c_int64_t
    integer(c_int)     :: element_bytes = GPU_WS_ELEMENT_BYTES
  end type gpu_ws_buffer_t

  public :: gpu_ws_target_is_valid
  public :: gpu_ws_residency_is_valid
  public :: gpu_ws_f3_policy_is_valid
  public :: gpu_ws_validate_table
  public :: gpu_ws_table_total_bytes

  ! C-ABI workspace runtime (implemented in source/gpu_workspace_runtime.c).
  ! METC-C1a backs these verbs with a real contiguous arena allocator + pointer
  ! lookup by byte offset.  Under OQP_CUDA_ENABLE the arena is a CUDA device
  ! allocation; otherwise it is a host allocation (test-safe fallback).
  interface
    !> Reserve a contiguous arena of `total_bytes` for `target`.
    !> Returns a positive handle id on success, or 0 on failure.
    function oqp_gpu_ws_acquire(target, total_bytes) result(handle) &
        bind(C, name="oqp_gpu_ws_acquire")
      import :: c_int, c_int64_t
      integer(c_int), value :: target
      integer(c_int64_t), value :: total_bytes
      integer(c_int) :: handle
    end function oqp_gpu_ws_acquire

    !> Confirm [offset, offset+nbytes) lies within the handle's arena.
    !> Returns 0 on success, nonzero error code otherwise.
    function oqp_gpu_ws_validate(handle, offset, nbytes) result(ierr) &
        bind(C, name="oqp_gpu_ws_validate")
      import :: c_int, c_int64_t
      integer(c_int), value :: handle
      integer(c_int64_t), value :: offset
      integer(c_int64_t), value :: nbytes
      integer(c_int) :: ierr
    end function oqp_gpu_ws_validate

    !> Return base_ptr + offset for the handle's arena, or C_NULL_PTR on an
    !> invalid handle or an out-of-range offset.
    function oqp_gpu_ws_ptr(handle, offset) result(p) &
        bind(C, name="oqp_gpu_ws_ptr")
      import :: c_int, c_int64_t, c_ptr
      integer(c_int), value :: handle
      integer(c_int64_t), value :: offset
      type(c_ptr) :: p
    end function oqp_gpu_ws_ptr

    !> Return the handle's arena size in bytes, or -1 for an invalid handle.
    function oqp_gpu_ws_total_bytes(handle) result(nbytes) &
        bind(C, name="oqp_gpu_ws_total_bytes")
      import :: c_int, c_int64_t
      integer(c_int), value :: handle
      integer(c_int64_t) :: nbytes
    end function oqp_gpu_ws_total_bytes

    !> Free the handle's arena.  Double release is safe (returns nonzero).
    function oqp_gpu_ws_release(handle) result(ierr) &
        bind(C, name="oqp_gpu_ws_release")
      import :: c_int
      integer(c_int), value :: handle
      integer(c_int) :: ierr
    end function oqp_gpu_ws_release
  end interface

  public :: oqp_gpu_ws_acquire
  public :: oqp_gpu_ws_validate
  public :: oqp_gpu_ws_ptr
  public :: oqp_gpu_ws_total_bytes
  public :: oqp_gpu_ws_release

contains

  pure logical function gpu_ws_target_is_valid(code) result(ok)
    integer(c_int), intent(in) :: code
    ok = (code == GPU_WS_TARGET_METC) .or. (code == GPU_WS_TARGET_XC_RESPONSE)
  end function gpu_ws_target_is_valid

  pure logical function gpu_ws_residency_is_valid(code) result(ok)
    integer(c_int), intent(in) :: code
    ok = (code == GPU_WS_RESIDENCY_HOST_ONLY) .or. &
         (code == GPU_WS_RESIDENCY_DEVICE_RESIDENT) .or. &
         (code == GPU_WS_RESIDENCY_MIRRORED) .or. &
         (code == GPU_WS_RESIDENCY_BORROWED)
  end function gpu_ws_residency_is_valid

  pure logical function gpu_ws_f3_policy_is_valid(code) result(ok)
    integer(c_int), intent(in) :: code
    ok = (code == GPU_WS_F3_PER_THREAD) .or. (code == GPU_WS_F3_SINGLE_ATOMIC)
  end function gpu_ws_f3_policy_is_valid

  !> Validate a workspace manifest table.  Checks that every row has a known
  !> target and residency class, a positive byte size, and a contiguous,
  !> non-overlapping offset.  Returns 0 on success, otherwise the 1-based index
  !> of the first offending row.  No device memory is touched.
  function gpu_ws_validate_table(buffers, n) result(ierr)
    type(gpu_ws_buffer_t), intent(in) :: buffers(:)
    integer, intent(in) :: n
    integer :: ierr
    integer :: i
    integer(c_int64_t) :: expected_offset

    ierr = 0
    expected_offset = 0_c_int64_t
    do i = 1, n
      if (.not. gpu_ws_target_is_valid(buffers(i)%target)) then
        ierr = i
        return
      end if
      if (.not. gpu_ws_residency_is_valid(buffers(i)%residency)) then
        ierr = i
        return
      end if
      if (buffers(i)%nbytes <= 0_c_int64_t) then
        ierr = i
        return
      end if
      if (buffers(i)%offset /= expected_offset) then
        ierr = i
        return
      end if
      expected_offset = expected_offset + buffers(i)%nbytes
    end do
  end function gpu_ws_validate_table

  pure function gpu_ws_table_total_bytes(buffers, n) result(total)
    type(gpu_ws_buffer_t), intent(in) :: buffers(:)
    integer, intent(in) :: n
    integer(c_int64_t) :: total
    integer :: i
    total = 0_c_int64_t
    do i = 1, n
      total = total + buffers(i)%nbytes
    end do
  end function gpu_ws_table_total_bytes

end module gpu_workspace_bridge
