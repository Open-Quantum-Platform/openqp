!> @brief openqp-gpu bridge for the SCF Fock (J/K), the DFT exchange-correlation
!>        potential, and the two-electron nuclear gradient (ground state + MRSF).
!>
!> Companion to routec_sig (the MRSF-TDDFT sigma-session). Where routec_sig
!> diverts the excited-state response, this module diverts the ground-state
!> two-electron work: the SCF Fock build, the semilocal Vxc, and the analytic
!> 2e gradient. All four entry points live in one library (libopenqp_gpu).
!>
!> Two ways to reach the library, selected at compile time (as in routec_sig):
!>
!>   * COMPILE-TIME LINK (build OpenQP with -DOQP_GPU_LINKED and link
!>     libopenqp_gpu; the CMake option OPENQP_WITH_GPU does this).  Each
!>     capability is on by default and can be turned off per run with its env
!>     var set to 0/off/none: OQP_ROUTEC_LIB (Fock J/K), OQP_ROUTEC_XC_LIB
!>     (Vxc), OQP_ROUTEC_GRAD_LIB (gradient).
!>
!>   * RUNTIME dlopen (default when not linked).  Each capability is inert
!>     unless its env var names a loadable dylib exporting the symbol.
!>
!> Every wrapper returns .true. only when the device call succeeds (info==0);
!> on absence or any decline the caller runs the native path unchanged. So the
!> seam is safe to ship inert and degrades to stock OpenQP.
module routec_bridge
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: routec_try_fock_jk, routec_try_vxc, &
            routec_try_grad2, routec_try_grad2_mrsf

  abstract interface
    !> J/K Fock build.  d,f: packed lower-triangular (nbf_tri, nfocks); f is
    !> returned ready to use (0.5 scaling + diagonal doubling done device-side).
    subroutine routec_fock_jk_i(d, f, nbf, nfocks, se, sc, info) bind(C)
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*)
      real(c_double), intent(inout) :: f(*)
      integer(c_int), intent(in)    :: nbf, nfocks
      real(c_double), intent(in)    :: se, sc
      integer(c_int), intent(out)   :: info
    end subroutine
    !> Semilocal Vxc.  d: packed total density; fxc: packed Vxc out; eexc: XC
    !> energy; totele: grid electron-count diagnostic.
    subroutine routec_vxc_i(d, fxc, nbf, nfocks, eexc, totele, info) bind(C)
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*)
      real(c_double), intent(inout) :: fxc(*)
      integer(c_int), intent(in)    :: nbf, nfocks
      real(c_double), intent(inout) :: eexc, totele
      integer(c_int), intent(out)   :: info
    end subroutine
    !> Ground-state 2e gradient.  d: packed total density; xyz: (3,natm) Bohr;
    !> de: (3,natm) receives dE2e/dx (Coulomb + hfscale*exchange).
    subroutine routec_grad2_i(d, xyz, de, nbf, natm, hfscale, coulscale, &
                              info) bind(C)
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*), xyz(*)
      real(c_double), intent(inout) :: de(*)
      integer(c_int), intent(in)    :: nbf, natm
      real(c_double), intent(in)    :: hfscale, coulscale
      integer(c_int), intent(out)   :: info
    end subroutine
    !> MRSF 2e gradient.  d,p: (nbf,nbf,2) RAW alpha/beta SCF and relaxed
    !> difference densities (before the grd2_mrsf init sum/diff); spc:
    !> (7,nbf,nbf) MRSF response densities; de: (3,natm).
    subroutine routec_grad2_mrsf_i(d, p, spc, xyz, de, nbf, natm, hfscale, &
                                   hfscale2, coulscale, spcscale, mrst, &
                                   info) bind(C)
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*), p(*), spc(*), xyz(*)
      real(c_double), intent(inout) :: de(*)
      integer(c_int), intent(in)    :: nbf, natm, mrst
      real(c_double), intent(in)    :: hfscale, hfscale2, coulscale
      real(c_double), intent(in)    :: spcscale(*)
      integer(c_int), intent(out)   :: info
    end subroutine
  end interface

#ifdef OQP_GPU_LINKED
  ! Compile-time link: the symbols are resolved by the linker.
  interface
    subroutine routec_fock_jk(d, f, nbf, nfocks, se, sc, info) &
               bind(C, name="routec_fock_jk")
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*)
      real(c_double), intent(inout) :: f(*)
      integer(c_int), intent(in)    :: nbf, nfocks
      real(c_double), intent(in)    :: se, sc
      integer(c_int), intent(out)   :: info
    end subroutine
    subroutine routec_vxc(d, fxc, nbf, nfocks, eexc, totele, info) &
               bind(C, name="routec_vxc")
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*)
      real(c_double), intent(inout) :: fxc(*)
      integer(c_int), intent(in)    :: nbf, nfocks
      real(c_double), intent(inout) :: eexc, totele
      integer(c_int), intent(out)   :: info
    end subroutine
    subroutine routec_grad2(d, xyz, de, nbf, natm, hfscale, coulscale, info) &
               bind(C, name="routec_grad2")
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*), xyz(*)
      real(c_double), intent(inout) :: de(*)
      integer(c_int), intent(in)    :: nbf, natm
      real(c_double), intent(in)    :: hfscale, coulscale
      integer(c_int), intent(out)   :: info
    end subroutine
    subroutine routec_grad2_mrsf(d, p, spc, xyz, de, nbf, natm, hfscale, &
                                 hfscale2, coulscale, spcscale, mrst, info) &
               bind(C, name="routec_grad2_mrsf")
      import :: c_double, c_int
      real(c_double), intent(in)    :: d(*), p(*), spc(*), xyz(*)
      real(c_double), intent(inout) :: de(*)
      integer(c_int), intent(in)    :: nbf, natm, mrst
      real(c_double), intent(in)    :: hfscale, hfscale2, coulscale
      real(c_double), intent(in)    :: spcscale(*)
      integer(c_int), intent(out)   :: info
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

  logical :: tried_jk = .false., tried_xc = .false., tried_grad = .false.
  logical :: noted_jk = .false., noted_xc = .false., noted_grad = .false.
  procedure(routec_fock_jk_i),    pointer :: ext_jk         => null()
  procedure(routec_vxc_i),        pointer :: ext_vxc        => null()
  procedure(routec_grad2_i),      pointer :: ext_grad2      => null()
  procedure(routec_grad2_mrsf_i), pointer :: ext_grad2_mrsf => null()

contains

  !> Runtime opt-out check: returns .true. iff $var is set to 0/off/no/none.
  logical function disabled(var)
    character(len=*), intent(in) :: var
    character(len=32) :: v
    integer :: n, s
    disabled = .false.
    call get_environment_variable(var, v, n, s)
    if (s /= 0 .or. n == 0) return
    select case (trim(adjustl(v)))
    case ("0", "off", "OFF", "no", "NO", "none", "NONE"); disabled = .true.
    end select
  end function disabled

#ifndef OQP_GPU_LINKED
  !> dlopen the dylib named by env var `var` and dlsym `sym` into `fp`.
  subroutine dlbind(var, sym, fp)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: var, sym
    type(c_funptr), intent(out)  :: fp
    type(c_ptr) :: h
    character(len=4096) :: path
    integer :: plen, stat
    fp = c_null_funptr
    call get_environment_variable(var, path, plen, stat)
    if (stat /= 0 .or. plen == 0) return
    h = c_dlopen(trim(path)//c_null_char, RTLD_NOW)
    if (.not. c_associated(h)) then
      write(error_unit,'(a)') ' routec_bridge: '//var//' set but dlopen failed: '//trim(path)
      return
    end if
    fp = c_dlsym(h, trim(sym)//c_null_char)
    if (.not. c_associated(fp)) &
      write(error_unit,'(a)') ' routec_bridge: '//trim(sym)//' not found in '//trim(path)
  end subroutine dlbind
#endif

  subroutine load_jk()
    tried_jk = .true.
    if (disabled("OQP_ROUTEC_LIB")) return
#ifdef OQP_GPU_LINKED
    ext_jk => routec_fock_jk
#else
    block
      type(c_funptr) :: fp
      call dlbind("OQP_ROUTEC_LIB", "routec_fock_jk", fp)
      if (c_associated(fp)) call c_f_procpointer(fp, ext_jk)
    end block
#endif
  end subroutine load_jk

  subroutine load_xc()
    tried_xc = .true.
    if (disabled("OQP_ROUTEC_XC_LIB")) return
#ifdef OQP_GPU_LINKED
    ext_vxc => routec_vxc
#else
    block
      type(c_funptr) :: fp
      call dlbind("OQP_ROUTEC_XC_LIB", "routec_vxc", fp)
      if (c_associated(fp)) call c_f_procpointer(fp, ext_vxc)
    end block
#endif
  end subroutine load_xc

  subroutine load_grad()
    tried_grad = .true.
    if (disabled("OQP_ROUTEC_GRAD_LIB")) return
#ifdef OQP_GPU_LINKED
    ext_grad2      => routec_grad2
    ext_grad2_mrsf => routec_grad2_mrsf
#else
    block
      type(c_funptr) :: fp
      call dlbind("OQP_ROUTEC_GRAD_LIB", "routec_grad2", fp)
      if (c_associated(fp)) call c_f_procpointer(fp, ext_grad2)
      call dlbind("OQP_ROUTEC_GRAD_LIB", "routec_grad2_mrsf", fp)
      if (c_associated(fp)) call c_f_procpointer(fp, ext_grad2_mrsf)
    end block
#endif
  end subroutine load_grad

  !> SCF J/K Fock build.  Returns .true. and fills f on device success.
  function routec_try_fock_jk(d, f, nbf, nfocks, se, sc) result(done)
    use, intrinsic :: iso_fortran_env, only: error_unit
    real(c_double), intent(in)    :: d(*)
    real(c_double), intent(inout) :: f(*)
    integer, intent(in)           :: nbf, nfocks
    real(c_double), intent(in)    :: se, sc
    logical :: done
    integer(c_int) :: info
    done = .false.
    if (.not. tried_jk) call load_jk()
    if (.not. associated(ext_jk)) return
    info = 1_c_int
    call ext_jk(d, f, int(nbf,c_int), int(nfocks,c_int), se, sc, info)
    done = (info == 0_c_int)
    if (done .and. .not. noted_jk) then
      noted_jk = .true.
      write(error_unit,'(a)') ' [routec_bridge] SCF Fock J/K on GPU'
    end if
  end function routec_try_fock_jk

  !> Semilocal Vxc build.  Returns .true. and fills fxc/eexc/totele on success.
  function routec_try_vxc(d, fxc, nbf, nfocks, eexc, totele) result(done)
    use, intrinsic :: iso_fortran_env, only: error_unit
    real(c_double), intent(in)    :: d(*)
    real(c_double), intent(inout) :: fxc(*)
    integer, intent(in)           :: nbf, nfocks
    real(c_double), intent(inout) :: eexc, totele
    logical :: done
    integer(c_int) :: info
    done = .false.
    if (.not. tried_xc) call load_xc()
    if (.not. associated(ext_vxc)) return
    info = 1_c_int
    call ext_vxc(d, fxc, int(nbf,c_int), int(nfocks,c_int), eexc, totele, info)
    done = (info == 0_c_int)
    if (done .and. .not. noted_xc) then
      noted_xc = .true.
      write(error_unit,'(a)') ' [routec_bridge] DFT Vxc on GPU'
    end if
  end function routec_try_vxc

  !> Ground-state 2e gradient.  Returns .true. and fills de(3,natm) on success.
  function routec_try_grad2(d, xyz, de, nbf, natm, hfscale, coulscale) result(done)
    use, intrinsic :: iso_fortran_env, only: error_unit
    real(c_double), intent(in)    :: d(*), xyz(*)
    real(c_double), intent(inout) :: de(*)
    integer, intent(in)           :: nbf, natm
    real(c_double), intent(in)    :: hfscale, coulscale
    logical :: done
    integer(c_int) :: info
    done = .false.
    if (.not. tried_grad) call load_grad()
    if (.not. associated(ext_grad2)) return
    info = 1_c_int
    call ext_grad2(d, xyz, de, int(nbf,c_int), int(natm,c_int), hfscale, coulscale, info)
    done = (info == 0_c_int)
    if (done .and. .not. noted_grad) then
      noted_grad = .true.
      write(error_unit,'(a)') ' [routec_bridge] 2e nuclear gradient on GPU'
    end if
  end function routec_try_grad2

  !> MRSF 2e gradient.  Densities d,p,spc are RAW (before grd2_mrsf init).
  function routec_try_grad2_mrsf(d, p, spc, xyz, de, nbf, natm, hfscale, &
                                 hfscale2, coulscale, spcscale, mrst) result(done)
    real(c_double), intent(in)    :: d(*), p(*), spc(*), xyz(*)
    real(c_double), intent(inout) :: de(*)
    integer, intent(in)           :: nbf, natm, mrst
    real(c_double), intent(in)    :: hfscale, hfscale2, coulscale
    real(c_double), intent(in)    :: spcscale(3)
    logical :: done
    integer(c_int) :: info
    done = .false.
    if (.not. tried_grad) call load_grad()
    if (.not. associated(ext_grad2_mrsf)) return
    info = 1_c_int
    call ext_grad2_mrsf(d, p, spc, xyz, de, int(nbf,c_int), int(natm,c_int), &
                        hfscale, hfscale2, coulscale, spcscale, int(mrst,c_int), info)
    done = (info == 0_c_int)
  end function routec_try_grad2_mrsf

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

end module routec_bridge
