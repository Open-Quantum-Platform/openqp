!> Shared helpers for the TDHF/SF/MRSF z-vector (CPHF/CPKS) solvers.
!>
!> These routines were previously duplicated, one copy per response module
!> (`tdhf_z_vector`, `tdhf_sf_z_vector`, `tdhf_mrsf_z_vector`). They are
!> collected here so the three z-vector drivers share a single, tested
!> implementation.
module zvector_common

  use precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public :: sanitize_zvector_preconditioner
  public :: zv_opts_t, zv_read_opts, zv_prog_tau, zv_warm_get, zv_warm_put
  public :: ZV_SLOT_SF, ZV_SLOT_TDHF

  !> Shared performance opt-ins for the (RHF/SF) z-vector CG solvers, mirroring
  !> the MRSF implementation. Read once per solve from env OQP_<PREFIX>_ZV_*.
  type :: zv_opts_t
    logical :: warm_on = .false.   !< warm-start across geometry/MD steps (opt-in)
    logical :: prog_on = .false.   !< progressive (iteration-dependent) screening (opt-in)
    logical :: diag_on = .false.   !< Jacobi cold-start guess x0 = M^-1 rhs (opt-in)
    logical :: timers  = .false.   !< per-section profiler
    real(dp) :: conv_user = -1.0_dp   !< override cnvtol (<0 => use input)
    real(dp) :: prog_k    = 1.0e-2_dp !< tau = prog_k * ||r||
    real(dp) :: prog_cap  = 1.0e-6_dp !< loosest tau
    real(dp) :: prog_pin  = 1.0e-6_dp !< pin tight once ||r||^2 < this
  end type zv_opts_t

  ! Warm-start stores, one per method (persist across steps in one process).
  integer, parameter :: ZV_SLOT_SF = 1, ZV_SLOT_TDHF = 2, ZV_NSLOT = 4
  type :: zv_store_t
    real(dp), allocatable :: vec(:,:)
    logical,  allocatable :: has(:)
    integer :: lzdim = 0
  end type zv_store_t
  type(zv_store_t), save :: zv_stores(ZV_NSLOT)

contains

  logical function zv_is_false(s)
    character(len=*), intent(in) :: s
    character :: c
    c = s(1:1)
    zv_is_false = (c=='0' .or. c=='n' .or. c=='N' .or. c=='f' .or. c=='F')
  end function zv_is_false

  !> Read OQP_<prefix>_ZV_* opt-ins (warm-start/progressive/Jacobi default ON).
  subroutine zv_read_opts(opts, prefix)
    type(zv_opts_t), intent(out) :: opts
    character(len=*), intent(in) :: prefix
    character(len=64) :: e_
    integer :: ios
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_WARMSTART', e_)
    if (len_trim(e_) > 0) opts%warm_on = .not. zv_is_false(e_)
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_PROG', e_)
    if (len_trim(e_) > 0) opts%prog_on = .not. zv_is_false(e_)
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_DIAGGUESS', e_)
    if (len_trim(e_) > 0) opts%diag_on = .not. zv_is_false(e_)
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_TIMERS', e_)
    opts%timers = len_trim(e_) > 0
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_CONV', e_)
    if (len_trim(e_) > 0) then
      read(e_,*,iostat=ios) opts%conv_user
      if (ios /= 0) opts%conv_user = -1.0_dp
    end if
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_PROG_CAP', e_)
    if (len_trim(e_) > 0) read(e_,*,iostat=ios) opts%prog_cap
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_PROG_K', e_)
    if (len_trim(e_) > 0) read(e_,*,iostat=ios) opts%prog_k
    call get_environment_variable('OQP_'//trim(prefix)//'_ZV_PROG_PIN', e_)
    if (len_trim(e_) > 0) read(e_,*,iostat=ios) opts%prog_pin
  end subroutine zv_read_opts

  !> Progressive screening threshold for a CG step (tight once pinned).
  function zv_prog_tau(opts, error, tight) result(tau)
    type(zv_opts_t), intent(in) :: opts
    real(dp), intent(in) :: error, tight
    real(dp) :: tau, resid
    if (.not. ieee_is_finite(error) .or. error < opts%prog_pin) then
      tau = tight; return
    end if
    resid = sqrt(max(error, 0.0_dp))
    tau = opts%prog_k * resid
    if (tau < tight) tau = tight
    if (tau > opts%prog_cap) tau = opts%prog_cap
  end function zv_prog_tau

  !> Cache a converged z-vector for warm-starting the next step (per method slot).
  subroutine zv_warm_put(slot, vec, lzdim, state, nstate)
    integer, intent(in) :: slot, lzdim, state, nstate
    real(dp), intent(in) :: vec(:)
    integer :: ncol
    if (slot < 1 .or. slot > ZV_NSLOT) return
    if (any(.not. ieee_is_finite(vec))) return
    ncol = max(nstate, state)
    if (zv_stores(slot)%lzdim /= lzdim .or. .not. allocated(zv_stores(slot)%vec)) then
      if (allocated(zv_stores(slot)%vec)) deallocate(zv_stores(slot)%vec)
      if (allocated(zv_stores(slot)%has)) deallocate(zv_stores(slot)%has)
      allocate(zv_stores(slot)%vec(lzdim, ncol), source=0.0_dp)
      allocate(zv_stores(slot)%has(ncol), source=.false.)
      zv_stores(slot)%lzdim = lzdim
    else if (size(zv_stores(slot)%has) < state) then
      return
    end if
    zv_stores(slot)%vec(:,state) = vec
    zv_stores(slot)%has(state) = .true.
  end subroutine zv_warm_put

  !> Fetch a warm-start guess; returns .true. only on a finite, matching hit.
  function zv_warm_get(slot, vec, lzdim, state) result(used)
    integer, intent(in) :: slot, lzdim, state
    real(dp), intent(out) :: vec(:)
    logical :: used
    used = .false.
    vec = 0.0_dp
    if (slot < 1 .or. slot > ZV_NSLOT) return
    if (.not. allocated(zv_stores(slot)%vec) .or. .not. allocated(zv_stores(slot)%has)) return
    if (zv_stores(slot)%lzdim /= lzdim) return
    if (state < 1 .or. state > size(zv_stores(slot)%has)) return
    if (.not. zv_stores(slot)%has(state)) return
    if (any(.not. ieee_is_finite(zv_stores(slot)%vec(:,state)))) return
    vec = zv_stores(slot)%vec(:,state)
    used = .true.
  end function zv_warm_get

!> @brief Build a finite diagonal (Jacobi) preconditioner `xminv = 1/xm`.
!>
!> Each denominator that is non-finite or smaller in magnitude than `floor`
!> is replaced by `+/-floor` (sign preserved), so the preconditioner can never
!> introduce a NaN/Inf or an overflow. The number of regularized entries is
!> reported on `log_unit`, tagged with `tag` (e.g. "RHF", "SF", "MRSF").
!>
!> `floor` is supplied by the caller because the response modules use slightly
!> different thresholds (1e-12 for RHF/SF, 1e-14 for MRSF).
  subroutine sanitize_zvector_preconditioner(xm, xminv, log_unit, floor, tag)
    real(kind=dp),    intent(in)  :: xm(:)
    real(kind=dp),    intent(out) :: xminv(:)
    integer,          intent(in)  :: log_unit
    real(kind=dp),    intent(in)  :: floor
    character(len=*), intent(in), optional :: tag

    integer :: i, regularized
    real(kind=dp) :: denom
    character(len=16) :: prefix

    prefix = ''
    if (present(tag)) prefix = trim(tag)//' '

    regularized = 0
    do i = 1, size(xm)
      denom = xm(i)
      if (.not. ieee_is_finite(denom) .or. abs(denom) < floor) then
        if (ieee_is_finite(denom) .and. denom < 0.0_dp) then
          denom = -floor
        else
          denom =  floor
        end if
        regularized = regularized + 1
      end if
      xminv(i) = 1.0_dp / denom
    end do

    if (regularized > 0) then
      write(log_unit,'(1x,A,"z-vector preconditioner regularized ",I0," denominator(s)")') &
            trim(prefix), regularized
      call flush(log_unit)
    end if

  end subroutine sanitize_zvector_preconditioner

end module zvector_common
