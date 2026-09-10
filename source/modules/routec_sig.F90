!> @brief OQP_ROUTEC_SIG — device-resident MRSF sigma-session bridge.
!>
!> Replaces the ENTIRE per-vector sigma triple
!>   mrsfcbc (6a) -> int2_driver%run (6b) -> mrsfmntoia+mrsfesum (6c)
!> with a single device call: given the MO trial vectors bvec_mo(:,ist:iend)
!> it returns amo(:,ist:iend) = (A-B).X in the MO frame. The Davidson subspace
!> solver is untouched.
!>
!> Two ways to reach the GPU library (identical ABI, identical behaviour):
!>
!>   * COMPILE-TIME LINK (build OpenQP with -DOQP_GPU_LINKED and link
!>     libopenqp_gpu; the CMake option OPENQP_WITH_GPU auto-downloads and does
!>     this).  routec_sig_init/set_scale/iter/free are resolved by the linker;
!>     the session is on by default and can be turned off at runtime with
!>     OQP_ROUTEC_SIG=0 / off / none.
!>
!>   * RUNTIME dlopen (default when not linked).  Inert unless $OQP_ROUTEC_SIG
!>     names a loadable dylib exporting the mandatory symbols.
!>
!> When inert the driver runs the byte-for-byte native path (the referee).
module routec_sig
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: routec_sig_available, routec_sig_begin, routec_sig_apply, &
            routec_sig_end

  abstract interface
    !> int routec_sig_init(const int* nbf, const double* mo_a, const double* mo_b,
    !>     const double* fmo_a, const double* fmo_b, const int* nocca,
    !>     const int* noccb, const int* kind)
    function routec_sig_init_i(nbf, mo_a, mo_b, fmo_a, fmo_b, &
                               nocca, noccb, kind) result(ierr) bind(C)
      import :: c_int, c_double
      integer(c_int), intent(in) :: nbf, nocca, noccb, kind
      real(c_double), intent(in) :: mo_a(*), mo_b(*), fmo_a(*), fmo_b(*)
      integer(c_int) :: ierr
    end function
    !> void routec_sig_set_scale(const double* s)
    subroutine routec_sig_setscale_i(s) bind(C)
      import :: c_double
      real(c_double), intent(in) :: s
    end subroutine
    !> void routec_sig_iter(const double* bvec_mo, const int* nv_new,
    !>     double* sigma_mo, int* info)
    subroutine routec_sig_iter_i(bvec_mo, nv_new, sigma_mo, info) bind(C)
      import :: c_double, c_int
      real(c_double), intent(in)  :: bvec_mo(*)
      integer(c_int), intent(in)  :: nv_new
      real(c_double), intent(out) :: sigma_mo(*)
      integer(c_int), intent(out) :: info
    end subroutine
    !> void routec_sig_free(void)
    subroutine routec_sig_free_i() bind(C)
    end subroutine
  end interface

#ifdef OQP_GPU_LINKED
  ! Compile-time link: the symbols are resolved by the linker.
  interface
    function routec_sig_init(nbf, mo_a, mo_b, fmo_a, fmo_b, nocca, noccb, kind) &
             result(ierr) bind(C, name="routec_sig_init")
      import :: c_int, c_double
      integer(c_int), intent(in) :: nbf, nocca, noccb, kind
      real(c_double), intent(in) :: mo_a(*), mo_b(*), fmo_a(*), fmo_b(*)
      integer(c_int) :: ierr
    end function
    subroutine routec_sig_set_scale(s) bind(C, name="routec_sig_set_scale")
      import :: c_double
      real(c_double), intent(in) :: s
    end subroutine
    subroutine routec_sig_iter(bvec_mo, nv_new, sigma_mo, info) &
               bind(C, name="routec_sig_iter")
      import :: c_double, c_int
      real(c_double), intent(in)  :: bvec_mo(*)
      integer(c_int), intent(in)  :: nv_new
      real(c_double), intent(out) :: sigma_mo(*)
      integer(c_int), intent(out) :: info
    end subroutine
    subroutine routec_sig_free() bind(C, name="routec_sig_free")
    end subroutine
  end interface
#else
  ! Runtime dynamic-loader path.
#ifdef OQP_DL_WIN32
  ! Windows has no libdl.  kernel32's LoadLibraryA/GetProcAddress are the
  ! equivalents (and are always available to the link), so bind those and wrap
  ! them below in c_dlopen/c_dlsym so every call site stays platform-neutral.
  ! LoadLibraryA has no counterpart to dlopen's mode, which is ignored.
  interface
    function c_loadlibrary(file) bind(C, name="LoadLibraryA")
      import :: c_ptr, c_char
      character(kind=c_char), intent(in) :: file(*)
      type(c_ptr) :: c_loadlibrary
    end function
    function c_getprocaddress(handle, name) bind(C, name="GetProcAddress")
      import :: c_ptr, c_char, c_funptr
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: name(*)
      type(c_funptr) :: c_getprocaddress
    end function
  end interface
#else
  interface
    function c_dlopen(file, mode) bind(C, name="dlopen")
      import :: c_ptr, c_char, c_int
      character(kind=c_char), intent(in) :: file(*)
      integer(c_int), value :: mode
      type(c_ptr) :: c_dlopen
    end function
    function c_dlsym(handle, name) bind(C, name="dlsym")
      import :: c_ptr, c_char, c_funptr
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: name(*)
      type(c_funptr) :: c_dlsym
    end function
  end interface
#endif
  integer(c_int), parameter :: RTLD_NOW = 2_c_int
#endif

  logical :: tried = .false.
  logical :: noted = .false.
  procedure(routec_sig_init_i),     pointer :: ext_init  => null()
  procedure(routec_sig_setscale_i), pointer :: ext_scale => null()
  procedure(routec_sig_iter_i),     pointer :: ext_iter  => null()
  procedure(routec_sig_free_i),     pointer :: ext_free  => null()

contains

  !> Resolve the sigma-session entry points, once.  In the linked build the
  !> pointers bind straight to the C symbols (subject to an OQP_ROUTEC_SIG
  !> runtime opt-out); otherwise dlopen the dylib named by OQP_ROUTEC_SIG.
  subroutine routec_sig_load()
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=4096) :: path
    integer :: plen, stat
#ifndef OQP_GPU_LINKED
    type(c_ptr) :: h
    type(c_funptr) :: fp
#endif
    tried = .true.
#ifdef OQP_GPU_LINKED
    ! Linked in: on by default, opt out with OQP_ROUTEC_SIG=0/off/none.
    call get_environment_variable("OQP_ROUTEC_SIG", path, plen, stat)
    if (stat == 0 .and. plen > 0) then
      select case (trim(adjustl(path)))
      case ("0", "off", "OFF", "none", "NONE", "no", "NO")
        return
      end select
    end if
    ext_init  => routec_sig_init
    ext_scale => routec_sig_set_scale
    ext_iter  => routec_sig_iter
    ext_free  => routec_sig_free
#else
    call get_environment_variable("OQP_ROUTEC_SIG", path, plen, stat)
    if (stat /= 0 .or. plen == 0) return
    h = c_dlopen(trim(path)//c_null_char, RTLD_NOW)
    if (.not. c_associated(h)) then
      write(error_unit, '(a)') &
        ' routec_sig: OQP_ROUTEC_SIG set but dlopen failed: '//trim(path)
      return
    end if
    fp = c_dlsym(h, "routec_sig_init"//c_null_char)
    if (.not. c_associated(fp)) then
      write(error_unit, '(a)') &
        ' routec_sig: routec_sig_init symbol not found in '//trim(path)
      return
    end if
    call c_f_procpointer(fp, ext_init)
    fp = c_dlsym(h, "routec_sig_iter"//c_null_char)
    if (.not. c_associated(fp)) then
      write(error_unit, '(a)') &
        ' routec_sig: routec_sig_iter symbol not found in '//trim(path)
      ext_init => null()
      return
    end if
    call c_f_procpointer(fp, ext_iter)
    ! set_scale is MANDATORY: without it the device sigma would silently run at
    ! the library's default exchange scale, which is wrong for scaled/hybrid
    ! DFT (the MRSF gate also enables non-CAM DFT). Treat its absence as an
    ! unavailable session so the driver keeps the native path.
    fp = c_dlsym(h, "routec_sig_set_scale"//c_null_char)
    if (.not. c_associated(fp)) then
      write(error_unit, '(a)') &
        ' routec_sig: routec_sig_set_scale symbol not found in '//trim(path)
      ext_init => null()
      ext_iter => null()
      return
    end if
    call c_f_procpointer(fp, ext_scale)
    ! optional
    fp = c_dlsym(h, "routec_sig_free"//c_null_char)
    if (c_associated(fp)) call c_f_procpointer(fp, ext_free)
#endif
  end subroutine routec_sig_load

  !> .true. iff the sigma-session (init+iter+set_scale) is resolved.
  logical function routec_sig_available()
    if (.not. tried) call routec_sig_load()
    routec_sig_available = associated(ext_init) .and. associated(ext_iter) &
                     .and. associated(ext_scale)
  end function routec_sig_available

  !> Begin a Davidson sigma-session: upload Ca,Cb (MO coeff), Fa,Fb (MO Fock,
  !> the unpacked square C^T F_AO C), set kind from mrst, set HF scale.
  !> Returns 0 on success (engine ready), nonzero otherwise.
  integer function routec_sig_begin(nbf, mo_a, mo_b, fa, fb, &
                                    nocca, noccb, mrst, scale) result(ierr)
    integer, intent(in) :: nbf, nocca, noccb, mrst
    real(c_double), intent(in) :: mo_a(*), mo_b(*), fa(*), fb(*)
    real(kind=c_double), intent(in) :: scale
    integer(c_int) :: kind_c
    ierr = -1
    if (.not. routec_sig_available()) return
    kind_c = merge(3_c_int, 1_c_int, mrst == 3)   ! kind 1=singlet, 3=triplet
    ierr = int(ext_init(int(nbf, c_int), mo_a, mo_b, fa, fb, &
                        int(nocca, c_int), int(noccb, c_int), kind_c))
    if (ierr == 0 .and. associated(ext_scale)) &
        call ext_scale(real(scale, c_double))
  end function routec_sig_begin

  !> One Davidson iteration: bvec_mo(ntrial,nv_new) IN, sigma_mo(ntrial,nv_new)
  !> OUT = (A-B).X for the nv_new new columns. Returns .true. on success.
  logical function routec_sig_apply(bvec_mo, nv_new, sigma_mo) result(ok)
    use, intrinsic :: iso_fortran_env, only: error_unit
    real(c_double), intent(in)  :: bvec_mo(*)
    integer, intent(in)         :: nv_new
    real(c_double), intent(out) :: sigma_mo(*)
    integer(c_int) :: info
    info = 1_c_int
    call ext_iter(bvec_mo, int(nv_new, c_int), sigma_mo, info)
    ok = (info == 0_c_int)
    if (ok .and. .not. noted) then
      noted = .true.
      write(error_unit, '(a)') &
        ' [routec_sig] MRSF Davidson diverted to GPU sigma-session'
    end if
  end function routec_sig_apply

  subroutine routec_sig_end()
    if (associated(ext_free)) call ext_free()
  end subroutine routec_sig_end

#ifdef OQP_DL_WIN32
  !> dlopen/dlsym shims over the Win32 loader (see the interface block above).
  function c_dlopen(file, mode) result(handle)
    character(kind=c_char, len=*), intent(in) :: file
    integer(c_int), intent(in) :: mode
    type(c_ptr) :: handle
    integer(c_int) :: ignored
    ignored = mode
    handle = c_loadlibrary(file)
  end function c_dlopen

  function c_dlsym(handle, name) result(fp)
    type(c_ptr), intent(in) :: handle
    character(kind=c_char, len=*), intent(in) :: name
    type(c_funptr) :: fp
    fp = c_getprocaddress(handle, name)
  end function c_dlsym
#endif

end module routec_sig
