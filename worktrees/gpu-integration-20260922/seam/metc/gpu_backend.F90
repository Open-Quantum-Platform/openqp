module gpu_backend

  use iso_c_binding, only: c_bool, c_char, c_double, c_int, c_null_char, c_ptr, c_loc
  use iso_fortran_env, only: error_unit

  implicit none

  logical :: gpu_metc_requested = .false.
  integer(c_int) :: gpu_device_id = 0_c_int

  ! Active METC residency session (0 = none) and the pass whose d3 is resident.
  ! Module-global: one METC GPU session at a time, driven across Davidson passes.
  integer(c_int), private :: metc_session = 0_c_int
  integer, private :: metc_pass_uploaded = -1

  ! ------------------------------------------------------------------------
  ! METC session state machine (all-or-nothing correctness guard).
  !
  !   GPU_METC_INACTIVE : no device accumulation has happened yet.  The GPU may
  !                       be attempted, and a failure can safely fall back to a
  !                       FULL host accumulation (host f3 is pristine).
  !   GPU_METC_ACTIVE   : a session is running and the resident device f3 holds
  !                       (partial) sums.  The host f3 is NOT being written, so a
  !                       later flush failure can no longer fall back per-flush
  !                       without mixing partial device f3 into host f3.
  !   GPU_METC_FAILED   : the session failed while ACTIVE (or, under strict mode,
  !                       failed to start).  Sticky and fatal: the caller MUST
  !                       abort rather than silently substitute / mix a host f3.
  !
  ! Integral data is streamed flush-by-flush and discarded after each flush, so
  ! once we are ACTIVE there is no way to replay earlier flushes on the host --
  ! a clean full host rebuild is impossible mid-accumulation.  The only safe
  ! response to an ACTIVE failure is therefore to abort, never to fall back.
  ! ------------------------------------------------------------------------
  integer, parameter, public :: GPU_METC_INACTIVE = 0
  integer, parameter, public :: GPU_METC_ACTIVE   = 1
  integer, parameter, public :: GPU_METC_FAILED   = 2
  integer, private :: metc_state = GPU_METC_INACTIVE

  ! Sticky "do not attempt the GPU again during THIS accumulation" flag, set
  ! after a soft (non-fatal) begin failure so we never start a session late and
  ! mix a partial device f3 with the host f3 the rest of the flushes built.
  logical, private :: metc_soft_disabled = .false.

  ! Strict mode: when .true., even a begin-time failure (GPU never started) is
  ! fatal instead of falling back to the CPU.  Loaded once from the
  ! OQP_GPU_METC_STRICT environment variable so the Fortran and Python layers
  ! agree without a dedicated binding.
  logical, private :: metc_strict = .false.
  logical, private :: metc_strict_read = .false.

  ! Non-zero ierr returned when the GPU path is intentionally skipped (poisoned
  ! or soft-disabled) so the hot path is not entered.
  integer, parameter, private :: METC_IERR_SKIP = 1

  ! C workspace/METC residency runtime (source/gpu_workspace_runtime.c).  These
  ! are always present (host-fallback or CUDA); the Fortran side only calls them
  ! when gpu_backend_metc_enabled() is true.
  interface
    function oqp_gpu_metc_session_begin(nf, nmatrix, nbf, nthreads, max_ncur) &
        result(handle) bind(C, name="oqp_gpu_metc_session_begin")
      import :: c_int
      integer(c_int), value :: nf, nmatrix, nbf, nthreads, max_ncur
      integer(c_int) :: handle
    end function oqp_gpu_metc_session_begin

    function oqp_gpu_metc_session_zero_f3(session) result(ierr) &
        bind(C, name="oqp_gpu_metc_session_zero_f3")
      import :: c_int
      integer(c_int), value :: session
      integer(c_int) :: ierr
    end function oqp_gpu_metc_session_zero_f3

    function oqp_gpu_metc_session_upload_d3(session, d3) result(ierr) &
        bind(C, name="oqp_gpu_metc_session_upload_d3")
      import :: c_int, c_ptr
      integer(c_int), value :: session
      type(c_ptr), value :: d3
      integer(c_int) :: ierr
    end function oqp_gpu_metc_session_upload_d3

    function oqp_gpu_metc_session_contract(session, thread, ids, ints, ncur, &
        cur_pass, scale_exchange, scale_coulomb, is_umrsf) result(ierr) &
        bind(C, name="oqp_gpu_metc_session_contract")
      import :: c_int, c_ptr, c_double
      integer(c_int), value :: session, thread, ncur, cur_pass, is_umrsf
      type(c_ptr), value :: ids, ints
      real(c_double), value :: scale_exchange, scale_coulomb
      integer(c_int) :: ierr
    end function oqp_gpu_metc_session_contract

    function oqp_gpu_metc_session_download_f3(session, f3) result(ierr) &
        bind(C, name="oqp_gpu_metc_session_download_f3")
      import :: c_int, c_ptr
      integer(c_int), value :: session
      type(c_ptr), value :: f3
      integer(c_int) :: ierr
    end function oqp_gpu_metc_session_download_f3

    function oqp_gpu_metc_session_end(session) result(ierr) &
        bind(C, name="oqp_gpu_metc_session_end")
      import :: c_int
      integer(c_int), value :: session
      integer(c_int) :: ierr
    end function oqp_gpu_metc_session_end
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

  function gpu_backend_metc_active() result(active)
    logical :: active
    active = (metc_session /= 0_c_int)
  end function gpu_backend_metc_active

  !> Current METC session state (GPU_METC_INACTIVE/ACTIVE/FAILED).
  function gpu_backend_metc_state() result(st)
    integer :: st
    st = metc_state
  end function gpu_backend_metc_state

  !> True once a session has failed while ACTIVE (or failed to start under
  !> strict mode).  The caller MUST treat this as fatal: a partial device f3
  !> cannot be safely combined with, or replaced by, a host f3.
  function gpu_backend_metc_failed() result(failed)
    logical :: failed
    failed = (metc_state == GPU_METC_FAILED)
  end function gpu_backend_metc_failed

  !> True when strict (no-CPU-fallback) mode is requested.
  function gpu_backend_metc_strict() result(strict)
    logical :: strict
    call ensure_strict_loaded()
    strict = metc_strict
  end function gpu_backend_metc_strict

  !> Explicit strict-mode override (primarily for tests); also marks the
  !> env-derived default as already loaded.
  subroutine gpu_backend_metc_set_strict(flag)
    logical, intent(in) :: flag
    metc_strict = flag
    metc_strict_read = .true.
  end subroutine gpu_backend_metc_set_strict

  !> Load strict mode from OQP_GPU_METC_STRICT exactly once.
  subroutine ensure_strict_loaded()
    character(len=16) :: env_value
    integer :: status
    if (metc_strict_read) return
    metc_strict_read = .true.
    metc_strict = .false.
    call get_environment_variable("OQP_GPU_STRICT", env_value, status=status)
    if (status /= 0) then
      call get_environment_variable("OQP_GPU_METC_STRICT", env_value, status=status)
    end if
    if (status == 0) then
      select case (trim(adjustl(env_value)))
      case ("1", "true", "TRUE", "on", "ON", "yes", "YES")
        metc_strict = .true.
      end select
    end if
  end subroutine ensure_strict_loaded

  !> Tear down any live session and free the device arena.  MUST be called with
  !> the oqp_gpu_metc_session critical region held (it mutates module session
  !> handles without its own locking).
  subroutine gpu_backend_metc_teardown_locked()
    integer :: rc
    if (metc_session /= 0_c_int) then
      rc = int(oqp_gpu_metc_session_end(metc_session))
    end if
    metc_session = 0_c_int
    metc_pass_uploaded = -1
  end subroutine gpu_backend_metc_teardown_locked

  !> Reset the METC state machine for a fresh accumulation scope.  Frees any
  !> stray session and clears the sticky soft-disable/poison flags so the next
  !> accumulation may attempt the GPU again.  Call at the start (cur_pass == 1)
  !> of a new build, NEVER mid-accumulation.
  subroutine gpu_backend_metc_reset()
    !$omp critical (oqp_gpu_metc_session)
    call gpu_backend_metc_teardown_locked()
    metc_state = GPU_METC_INACTIVE
    metc_soft_disabled = .false.
    !$omp end critical (oqp_gpu_metc_session)
  end subroutine gpu_backend_metc_reset

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

  !> Begin a resident METC session: acquire the arena, then zero the per-thread
  !> f3 accumulator once at the start of the Davidson accumulation scope.
  subroutine gpu_backend_metc_begin(nf, nmatrix, nbf, nthreads, max_ncur, ierr)
    integer, intent(in) :: nf, nmatrix, nbf, nthreads, max_ncur
    integer, intent(out) :: ierr

    ierr = 1
    metc_session = oqp_gpu_metc_session_begin(int(nf, c_int), int(nmatrix, c_int), &
                     int(nbf, c_int), int(nthreads, c_int), int(max_ncur, c_int))
    metc_pass_uploaded = -1
    if (metc_session /= 0_c_int) ierr = 0
  end subroutine gpu_backend_metc_begin

  subroutine gpu_backend_metc_zero_f3(ierr)
    integer, intent(out) :: ierr
    ierr = int(oqp_gpu_metc_session_zero_f3(metc_session))
  end subroutine gpu_backend_metc_zero_f3

  !> Upload/refresh the shared d3 input tensor for the current pass.  Idempotent
  !> within a pass: re-uploads only when the pass changes.
  subroutine gpu_backend_metc_upload_d3(d3, cur_pass, ierr)
    real(c_double), intent(in), target, contiguous :: d3(:,:,:,:)
    integer, intent(in) :: cur_pass
    integer, intent(out) :: ierr

    ierr = 0
    if (metc_pass_uploaded == cur_pass) return
    ierr = int(oqp_gpu_metc_session_upload_d3(metc_session, c_loc(d3(1,1,1,1))))
    if (ierr == 0) metc_pass_uploaded = cur_pass
  end subroutine gpu_backend_metc_upload_d3

  !> Resident hot-path contraction for one flush on one thread.  Copies this
  !> flush's ids/ints into the thread's resident scratch and launches the kernel
  !> against borrowed d3/f3 pointers.  No cudaMalloc/cudaFree.
  subroutine gpu_backend_metc_flush(thread, ids, ints, ncur, cur_pass, &
                                    scale_exchange, scale_coulomb, is_umrsf, ierr)
    integer, intent(in) :: thread, ncur, cur_pass
    integer(2), intent(in) :: ids(:,:)
    real(c_double), intent(in), target, contiguous :: ints(:)
    real(c_double), intent(in) :: scale_exchange, scale_coulomb
    logical, intent(in) :: is_umrsf
    integer, intent(out) :: ierr

    integer(c_int), allocatable, target :: ids32(:,:)
    integer :: n, umrsf_flag

    if (ncur <= 0) then
      ierr = 0
      return
    end if

    ! int16 buffer indices -> int32 for the C/CUDA ABI (host-side conversion).
    allocate(ids32(4, ncur))
    do n = 1, ncur
      ids32(1:4, n) = int(ids(1:4, n), c_int)
    end do
    umrsf_flag = 0
    if (is_umrsf) umrsf_flag = 1

    ierr = int(oqp_gpu_metc_session_contract(metc_session, int(thread, c_int), &
                 c_loc(ids32(1,1)), c_loc(ints(1)), int(ncur, c_int), &
                 int(cur_pass, c_int), scale_exchange, scale_coulomb, &
                 int(umrsf_flag, c_int)))
    deallocate(ids32)
  end subroutine gpu_backend_metc_flush

  !> One-call resident contraction attempt used by the MRSF/UMRSF update path.
  !> Lazily begins the session (single-threaded via OpenMP critical) and uploads
  !> d3 once per pass, then performs the per-thread flush outside the critical
  !> (each thread writes only its own f3/ids/ints slice).  Returns ierr/=0 to let
  !> the caller fall back to the CPU contraction.
  subroutine gpu_backend_metc_try_contract(nf, nmatrix, nbf, nthreads, cur_pass, &
                                           max_ncur, d3, thread, ids, ints, ncur, &
                                           scale_exchange, scale_coulomb, is_umrsf, ierr)
    integer, intent(in) :: nf, nmatrix, nbf, nthreads, cur_pass, max_ncur, thread, ncur
    real(c_double), intent(in), target, contiguous :: d3(:,:,:,:)
    integer(2), intent(in) :: ids(:,:)
    real(c_double), intent(in), target, contiguous :: ints(:)
    real(c_double), intent(in) :: scale_exchange, scale_coulomb
    logical, intent(in) :: is_umrsf
    integer, intent(out) :: ierr

    ierr = 0
    call ensure_strict_loaded()

    !$omp critical (oqp_gpu_metc_session)
    if (metc_state == GPU_METC_FAILED .or. metc_soft_disabled) then
      ! Poisoned, or the GPU was already abandoned for this accumulation.  Never
      ! start a late session: that would mix host-built and device-built f3.
      ierr = METC_IERR_SKIP
    else if (metc_session == 0_c_int) then
      ! First use this accumulation: try to start the resident session.
      call gpu_backend_metc_begin(nf, nmatrix, nbf, nthreads, max_ncur, ierr)
      if (ierr == 0) call gpu_backend_metc_zero_f3(ierr)
      if (ierr == 0) then
        ! Session live and device f3 zeroed: we are now committed to the GPU.
        metc_state = GPU_METC_ACTIVE
      else
        ! Begin failed before ANY device accumulation: nothing is resident, so
        ! the host can safely own the WHOLE accumulation -- unless the user
        ! demanded the GPU (strict), in which case this is fatal.
        call gpu_backend_metc_teardown_locked()
        if (metc_strict) then
          metc_state = GPU_METC_FAILED
        else
          metc_state = GPU_METC_INACTIVE
          metc_soft_disabled = .true.
        end if
      end if
    end if

    if (ierr == 0 .and. metc_state == GPU_METC_ACTIVE) then
      ! Refresh d3 for this pass.  A failure here is post-activation: device f3
      ! is already the live accumulator, so poison rather than fall back.
      call gpu_backend_metc_upload_d3(d3, cur_pass, ierr)
      if (ierr /= 0) metc_state = GPU_METC_FAILED
    end if
    !$omp end critical (oqp_gpu_metc_session)

    if (metc_state == GPU_METC_FAILED) then
      ! Strict begin failure or a failed d3 upload after activation.  Either way
      ! there is no safe host fallback that preserves correctness, so abort.
      call gpu_backend_metc_fatal()
      return
    end if
    if (ierr /= 0 .or. metc_state /= GPU_METC_ACTIVE) then
      ! Soft fallback: the GPU never became active for this accumulation, so the
      ! host path can safely own the entire build.  Signal the caller.
      ierr = METC_IERR_SKIP
      return
    end if

    call gpu_backend_metc_flush(thread, ids, ints, ncur, cur_pass, &
                                scale_exchange, scale_coulomb, is_umrsf, ierr)
    if (ierr /= 0) then
      ! A flush failed while the session was ACTIVE: the resident device f3 may
      ! hold partial sums.  Mixing them with a host fallback would corrupt the
      ! Fock matrix, and the streamed integrals for earlier flushes are already
      ! gone, so a clean host rebuild is impossible.  The only correct response
      ! is to abort (all-or-nothing contract).
      metc_state = GPU_METC_FAILED
      call gpu_backend_metc_fatal()
    end if
  end subroutine gpu_backend_metc_try_contract

  !> Fatal handler for an unrecoverable METC GPU failure.  Invoked only once a
  !> session is committed (ACTIVE) or when strict mode forbids a CPU fallback.
  !> Tears down the device session and stops execution rather than returning a
  !> mixed / silently incorrect Fock matrix.
  subroutine gpu_backend_metc_fatal()
    !$omp critical (oqp_gpu_metc_session)
    call gpu_backend_metc_teardown_locked()
    !$omp end critical (oqp_gpu_metc_session)
    write(error_unit, '(A)') "FATAL(gpu_backend): GPU METC contraction failed after the "// &
      "resident session became active (or strict mode is set)."
    write(error_unit, '(A)') "  Partial device results cannot be safely combined with a "// &
      "host fallback; aborting instead of returning an incorrect Fock matrix."
    write(error_unit, '(A)') "  Re-run with the GPU METC path disabled (unset "// &
      "OQP_GPU_METC) to use the CPU path."
    flush(error_unit)
    error stop 1
  end subroutine gpu_backend_metc_fatal

  !> Download the resident per-thread f3 accumulator into the host f3 array at
  !> the end of accumulation, then release the session.  The host f3 keeps its
  !> (nfocks,nmatrix,nbf,nbf,nthreads) layout so the existing thread reduction
  !> still applies.
  subroutine gpu_backend_metc_finalize(f3, ierr)
    real(c_double), intent(inout), target, contiguous :: f3(:,:,:,:,:)
    integer, intent(out) :: ierr
    integer :: rc

    ierr = 0

    !$omp critical (oqp_gpu_metc_session)
    if (metc_state /= GPU_METC_ACTIVE) then
      ! Defensive: never download a partial/poisoned device f3 into the host
      ! array.  (The caller is expected to have already aborted on FAILED.)
      call gpu_backend_metc_teardown_locked()
      metc_state = GPU_METC_INACTIVE
      metc_soft_disabled = .false.
    else
      ierr = int(oqp_gpu_metc_session_download_f3(metc_session, c_loc(f3(1,1,1,1,1))))
      rc = int(oqp_gpu_metc_session_end(metc_session))
      metc_session = 0_c_int
      metc_pass_uploaded = -1
      metc_state = GPU_METC_INACTIVE
      metc_soft_disabled = .false.
      if (ierr == 0) ierr = rc
    end if
    !$omp end critical (oqp_gpu_metc_session)
  end subroutine gpu_backend_metc_finalize

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
