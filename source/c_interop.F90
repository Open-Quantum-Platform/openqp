module c_interop

  use types, only: information
  use messages, only: without_abort
!  use types, only: oqp_handle_t
  use iso_c_binding, only: c_int, c_ptr, c_loc, c_f_pointer, c_associated, c_null_ptr, c_double, c_int32_t, c_int64_t, c_char

  implicit none

  private

  public oqp_handle_t

  public oqp_init
  public oqp_handle_refresh_ptr
  public oqp_handle_get_info

  interface oqp_handle_get_info
    module procedure oqp_handle_get_info_f
    module procedure oqp_handle_get_info_c
  end interface oqp_handle_get_info


  type, bind(C) :: oqp_handle_t
    type(c_ptr) :: inf
    type(c_ptr) :: xyz
    type(c_ptr) :: qn
    type(c_ptr) :: mass
    type(c_ptr) :: grad
    type(c_ptr) :: mol_prop
    type(c_ptr) :: mol_energy
    type(c_ptr) :: dft
    type(c_ptr) :: tddft
    type(c_ptr) :: control
    type(c_ptr) :: mpiinfo
    type(c_ptr) :: elshell
  end type

  ! Export buffers for oqp_get_basis. The public C API is fixed at int64_t, but
  ! basis_set%origin/am/ncontr are default Fortran integers, whose width follows
  ! the build (4 bytes in LP64 builds, e.g. native macOS Accelerate). Returning
  ! c_loc() of those internal arrays directly makes the Python side read pairs
  ! of 32-bit values as single 64-bit integers => corrupted basis metadata
  ! (e.g. centers [0,1] read back as [8589934592, -1]). Convert into these
  ! int64 buffers and export their addresses instead. Contents stay valid until
  ! the next oqp_get_basis call; callers (pyoqp get_basis) copy immediately.
  integer(c_int64_t), allocatable, target :: basis_am_i64(:)
  integer(c_int64_t), allocatable, target :: basis_origin_i64(:)
  integer(c_int64_t), allocatable, target :: basis_ncontr_i64(:)
  integer(c_int64_t), allocatable, target :: basis_spherical_i64(:)

contains

!--------------------------------------------------------------------------------

  function oqp_init() bind(C, name='oqp_init') result(res)

    implicit none
    type(c_ptr) :: res
    type(oqp_handle_t), pointer :: c_handle
    type(information), pointer :: inf
    integer :: ok

    res = c_null_ptr

    allocate(inf, stat=ok)
    if (ok /= 0) return

    allocate(c_handle, stat=ok)
    if (ok /= 0) return

    c_handle%inf = c_loc(inf)
    call oqp_handle_refresh_ptr(c_handle)

    res = c_loc(c_handle)

    call inf%dat%new("OQP")

  end function oqp_init

!--------------------------------------------------------------------------------

  function oqp_clean(c_handle) bind(C, name='oqp_clean') result(ok)

    implicit none
    integer(c_int) :: ok
    type(c_ptr), value :: c_handle

    type(oqp_handle_t), pointer :: f_handle
    type(information), pointer :: inf

    call c_f_pointer(c_handle, f_handle)
    call c_f_pointer(f_handle%inf, inf)

    call inf%dat%delete()

    deallocate(inf, stat=ok)

    if (ok/=0) return

    deallocate(f_handle, stat=ok)

  end function oqp_clean

!--------------------------------------------------------------------------------

  subroutine oqp_handle_refresh_ptr(c_handle)

    implicit none
    type(oqp_handle_t), intent(inout) :: c_handle
    type(information), pointer :: inf

    call c_f_pointer(c_handle%inf, inf)

    c_handle%mol_prop    = c_loc(inf%mol_prop)
    c_handle%mol_energy  = c_loc(inf%mol_energy)
    c_handle%dft         = c_loc(inf%dft)
    c_handle%control     = c_loc(inf%control)
    c_handle%tddft       = c_loc(inf%tddft)
    c_handle%mpiinfo     = c_loc(inf%mpiinfo)
    c_handle%elshell     = c_loc(inf%elshell)
    if (allocated(inf%atoms%xyz)) then
        c_handle%xyz  = c_loc(inf%atoms%xyz)
        c_handle%qn   = c_loc(inf%atoms%zn)
        c_handle%mass = c_loc(inf%atoms%mass)
    end if
    if (allocated(inf%atoms%grad)) then
        c_handle%grad = c_loc(inf%atoms%grad)
    end if

  end subroutine oqp_handle_refresh_ptr

!--------------------------------------------------------------------------------

  function oqp_handle_get_info_f(f_handle) result(res)

    implicit none
    type(oqp_handle_t), target :: f_handle
    type(information), pointer :: res

    call c_f_pointer(f_handle%inf, res)

  end function oqp_handle_get_info_f

!--------------------------------------------------------------------------------

  function oqp_handle_get_info_c(c_handle) result(res)

    implicit none
    type(c_ptr) :: c_handle
    type(information), pointer :: res
    type(oqp_handle_t), pointer :: f_handle

    call c_f_pointer(c_handle, f_handle)
    call c_f_pointer(f_handle%inf, res)

  end function oqp_handle_get_info_c

!--------------------------------------------------------------------------------

  function oqp_set_atoms(c_handle, natoms, x, y, z, q, mass) bind(C, name='oqp_set_atoms') result(ok)

    implicit none
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t), value :: natoms
    real(c_double) :: x(*), y(*), z(*), q(*)
    real(c_double), optional :: mass(*)
    integer(c_int) :: ok

    type(information), pointer :: inf

    ok = 10
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    ok = inf%set_atoms_arr(natoms, x, y, z, q, mass)
    if (ok/=0) return

    call oqp_handle_refresh_ptr(c_handle)

  end function oqp_set_atoms

!--------------------------------------------------------------------------------

  function oqp_get_atoms(c_handle, xyz) result(ok)
    type(oqp_handle_t) :: c_handle
    integer(c_int) :: ok
    real(c_double) :: xyz(3,*)

    integer :: nat
    type(information), pointer :: inf

    ok = 10
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    ok = 1
    nat = ubound(inf%atoms%xyz,2)
    xyz(:,1:nat) = inf%atoms%xyz

  end function

!--------------------------------------------------------------------------------

  function oqp_get_natom(c_handle) result(n) bind(C, name='oqp_get_natom')
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t) :: n

    type(information), pointer :: inf

    n = -1
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    n = ubound(inf%atoms%xyz,2)

  end function

!--------------------------------------------------------------------------------

  function oqp_get_nbf(c_handle) result(n) bind(C, name='oqp_get_nbf')
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t) :: n

    type(information), pointer :: inf

    n = -1
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    n = inf%basis%nbf

  end function

!--------------------------------------------------------------------------------

  function oqp_get_basis(c_handle, nsh, nprim, nbf, am, at, cdeg, ex, cc) result(ret) bind(C, name='oqp_get_basis')
    use basis_tools, only: basis_set
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t) :: nsh, nprim, nbf
    integer(c_int64_t) :: ret
    type(c_ptr), intent(out) :: am, at, cdeg, ex, cc

    type(information), pointer :: inf
    type(basis_set), pointer :: bas

    ret = -1
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    bas => inf%basis

    nbf = bas%nbf
    nprim = bas%nprim
    nsh = bas%nshell

    if (nbf <= 0) return

#define ADDRESSOF(a,b) if(allocated(a))then;b=c_loc(a);else;return;endif
    ADDRESSOF(bas%ex,    ex)
    ADDRESSOF(bas%cc,    cc)
#undef ADDRESSOF

    ! Integer arrays: do NOT export c_loc() of the internal default-integer
    ! arrays -- their width follows the build (4 bytes in LP64 builds) while
    ! the C API promises int64_t. Convert into the module-level int64 export
    ! buffers and hand out those addresses (see declarations above).
    if (.not. allocated(bas%am))     return
    if (.not. allocated(bas%origin)) return
    if (.not. allocated(bas%ncontr)) return
    basis_am_i64     = int(bas%am,     c_int64_t)
    basis_origin_i64 = int(bas%origin, c_int64_t)
    basis_ncontr_i64 = int(bas%ncontr, c_int64_t)
    am   = c_loc(basis_am_i64)
    at   = c_loc(basis_origin_i64)
    cdeg = c_loc(basis_ncontr_i64)

    ret = 0

  end function

!--------------------------------------------------------------------------------

!> @brief Export the per-shell "stored as spherical harmonics" flag.
!> @detail This is the EFFECTIVE flag, deliberately not a raw copy of
!>   `bas%harmonic`. OpenQP transforms a shell to the solid-harmonic set only
!>   when the run is harmonic, the shell is tagged pure, AND l >= 2 -- see
!>   c2s_ncomp in source/integrals/cart2sph.F90, which keeps s and p Cartesian
!>   (x, y, z) even in a spherical basis.
!>
!>   Consumers need that exact combination to know the AO layout, and deriving
!>   it outside this file has already gone wrong twice: tagging every shell
!>   pure applies the spherical component permutation to p shells and produces
!>   confidently WRONG signs, while tagging none pure (the previous state of
!>   the symmetry staging code) miscounts the AOs and silently disables the
!>   reduction on every spherical basis. Export the answer, not the ingredients.
!>
!> @warning LIFETIME AND THREADING. The returned address is that of the
!>   module-level `basis_spherical_i64`, which is SHARED across handles and
!>   reallocated on every call -- the same contract as `oqp_get_basis`'s
!>   `basis_am_i64` / `basis_origin_i64` / `basis_ncontr_i64` above, kept
!>   deliberately so the C API stays uniform. The caller must COPY the data
!>   before the next basis query and must not call this concurrently from two
!>   threads. pyoqp honours this: `oqpdata.py` copies into a fresh numpy array
!>   at the point of the call. A per-handle buffer would be the fix if
!>   concurrent C callers ever need supporting.
!> @param[out] nsh        number of shells
!> @param[out] spherical  address of an nsh-long int64 array, 1 = spherical
!> @return 0 on success, -1 if the handle or the basis is not usable
  function oqp_get_basis_spherical(c_handle, nsh, spherical) result(ret) &
      bind(C, name='oqp_get_basis_spherical')
    use basis_tools, only: basis_set
    use constants, only: HARMONIC_ACTIVE
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t) :: nsh
    type(c_ptr), intent(out) :: spherical
    integer(c_int64_t) :: ret

    type(information), pointer :: inf
    type(basis_set), pointer :: bas

    ret = -1
    nsh = 0
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    bas => inf%basis
    nsh = bas%nshell
    if (nsh <= 0) return
    if (.not. allocated(bas%harmonic)) return
    if (.not. allocated(bas%am)) return

    basis_spherical_i64 = merge(1_c_int64_t, 0_c_int64_t, &
        HARMONIC_ACTIVE .and. &
        bas%harmonic(1:bas%nshell) == 1 .and. &
        bas%am(1:bas%nshell) >= 2)
    spherical = c_loc(basis_spherical_i64)

    ret = 0

  end function

!--------------------------------------------------------------------------------

!> @brief Get calculation results from OQP handle
!> @param[in]    c_handle[in]  OQP handle
!> @param[in]    code[in]      Request string
!> @param[in]    v[out]        Pointer to data
!> @return       positive value:  success, returns size of the data
!>               negative values: -1 - handle not initialized;
!>                                -2 - data not available
!>                                -3 - unknown request code
  function oqp_get(c_handle, code, type_id, ndims, dims, v) result(n) bind(C, name='oqp_get')
    use strings, only: c_f_char
    use oqp_tagarray_driver, only: tagarray_get_cptr
    use tagarray_defines

    type(oqp_handle_t) :: c_handle
    integer(c_int64_t) :: n
    character(kind=c_char) :: code(*)
    type(c_ptr), intent(out) :: v
    integer(c_int32_t) :: type_id
    integer(c_int32_t) :: ndims
    integer(c_int64_t) :: dims(TA_MAX_DIMENSIONS_LENGTH)

    type(information), pointer :: inf
    character(:), allocatable :: code_str

    n = -1
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    code_str = trim(adjustl(c_f_char(code)))

    n = tagarray_get_cptr(inf%dat, code_str, v, type_id, ndims, dims)
    if (.not.c_associated(v))  n = -2


  end function

!--------------------------------------------------------------------------------

!> @brief Allocate storage in OQP handle
!> @param[in]    c_handle[in]  OQP handle
!> @param[in]    tag[in]       Tag to store data at
!> @param[in]    v[out]        Pointer to the data
!> @return       positive value:  success, returns size of the data
!>               negative values: -1 - handle not initialized;
!>                                -2 - data not available
  function oqp_alloc(c_handle, tag, type_id, ndims, dims, v) result(n) bind(C, name='oqp_alloc')
    use strings, only: c_f_char
    use oqp_tagarray_driver, only: tagarray_get_cptr
    use tagarray_defines

    type(oqp_handle_t) :: c_handle
    integer(c_int64_t) :: n
    character(kind=c_char) :: tag(*)
    type(c_ptr), intent(out) :: v
    integer(c_int32_t) :: type_id
    integer(c_int32_t) :: ndims
    integer(c_int64_t) :: dims(TA_MAX_DIMENSIONS_LENGTH)

    type(information), pointer :: inf
    character(:), allocatable :: tag_str
    integer(c_int64_t) :: data_size
    integer(c_int32_t) :: status_

    n = -1
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    tag_str = trim(adjustl(c_f_char(tag)))

    ! 1. allocate the memory in container (override, if it already exists)
    status_ = inf%dat%create(tag_str, type_id, dims(:ndims), override=.true.)
    ! 2. Get the pointer to the freshly allocated data
    n = tagarray_get_cptr(inf%dat, tag_str, v, type_id, ndims, dims, data_size)
    if (.not.c_associated(v))  n = -2

  end function

!--------------------------------------------------------------------------------

!> @brief Clean an entry in OQP handle
!> @param[in]    c_handle[in]  OQP handle
!> @param[in]    tag[in]       Data tag
!> @return       positive value:  success, returns size of the data
!>               negative values: -1 - handle not initialized;
!>                                -2 - tag not found
!>                                -3 - error removing data
  function oqp_del(c_handle, tag) result(n) bind(C, name='oqp_del')
    use strings, only: c_f_char

    type(oqp_handle_t) :: c_handle
    integer(c_int64_t) :: n
    character(kind=c_char) :: tag(*)

    type(information), pointer :: inf
    character(:), allocatable :: tag_str

    n = -1
    if (.not.c_associated(c_handle%inf)) return
    call c_f_pointer(c_handle%inf, inf)

    tag_str = trim(adjustl(c_f_char(tag)))

    n = -2
    ! A missing record is an ordinary query result for the C/Python deletion
    ! interface.  In particular, symmetry invalidation deliberately attempts
    ! to erase every optional staging record.  Do not route that expected
    ! condition through data_has_tags(), which writes a diagnostic to the
    ! default Fortran unit and creates an otherwise spurious fort.6 file.
    if (.not. inf%dat%contains(tag_str)) return

    n = -3
    call inf%dat%erase([tag_str])

    n = 0

  end function

!--------------------------------------------------------------------------------
end module c_interop
