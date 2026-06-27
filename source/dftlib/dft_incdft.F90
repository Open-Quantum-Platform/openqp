!> @brief Incremental DFT (IncDFT): reuse the exchange-correlation Fock/energy
!>        across SCF iterations where the density is no longer changing.
!>
!> @details OpenQP already builds the J/K (HF Coulomb/exchange) part of the Fock
!>   matrix incrementally from the density change (scf.F90 dold/fold), but the XC
!>   matrix is rebuilt from the full density on every iteration. As the SCF
!>   converges the density change dP sparsifies and the XC contribution stops
!>   changing, so rebuilding it is wasted quadrature.
!>
!>   XC is NONLINEAR in the density, so the reused matrix is an APPROXIMATION
!>   that must be refreshed. This module reuses the most recent full XC build
!>   (V_xc[D_ref], E_xc[D_ref]) only inside a controlled "late-SCF" window of the
!>   DIIS error, with a periodic forced full rebuild (drift hygiene, mirroring
!>   the J/K incremental reset cadence) and -- crucially -- a return to FULL XC
!>   builds once the error drops below `incdft_stop`, so the converged density is
!>   the true fixed point of the exact Fock and the converged energy/gradient are
!>   unchanged from the non-incremental baseline.
!>
!>   Opt-in via env OQP_XC_INCDFT (default off). Complementary to the Phi cache
!>   (Opt 1): the Phi cache removes the geometry-only collocation cost from every
!>   XC build, while IncDFT skips the density-driven XC work entirely on reused
!>   iterations.
!>
!> @author Claude (Anthropic), 2026
module mod_dft_incdft

  use precision, only: fp
  implicit none

  private

!> @brief Stored reference XC contribution from the last full build
  type, public :: incdft_t
    logical :: valid  = .false.
    integer :: ntri   = 0          !< packed matrix length nbf*(nbf+1)/2
    integer :: nf     = 0          !< number of spin blocks
    real(fp), allocatable :: vxc(:,:)   !< packed XC Fock (ntri, nf)
    real(fp) :: eexc   = 0.0_fp
    real(fp) :: totele = 0.0_fp
    real(fp) :: totkin = 0.0_fp
    integer  :: reuse_run = 0      !< consecutive reuses since the last full build
    integer  :: n_full  = 0        !< full XC builds this SCF (diagnostic)
    integer  :: n_reuse = 0        !< reused XC builds this SCF (diagnostic)
  contains
    procedure :: free
  end type

  type(incdft_t), public, save :: g_xc_ref

  ! Reuse-window controls (env-overridable; see incdft_load_controls)
  real(fp), public :: incdft_start    = 3.0e-2_fp  !< start reusing below this DIIS error
  real(fp), public :: incdft_stop     = 1.0e-4_fp  !< stop reusing below this (force exact final)
  integer,  public :: incdft_refresh  = 5          !< forced full rebuild every N reuses
  integer,  public :: incdft_min_iter = 2          !< never reuse before this iteration

  public :: incdft_env_enabled
  public :: incdft_should_reuse
  public :: incdft_reset
  public :: incdft_store

contains

!> @brief Is IncDFT enabled via the environment? (cached after first read)
  function incdft_env_enabled() result(en)
    logical :: en
    logical, save :: known = .false.
    logical, save :: value = .false.
    character(len=32) :: s
    integer :: st, ln
    if (.not. known) then
      call get_environment_variable('OQP_XC_INCDFT', s, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
        s = adjustl(s)
        value = (trim(s) == '1' .or. trim(s) == 'true' .or. trim(s) == 'TRUE' &
                 .or. trim(s) == 'on' .or. trim(s) == 'ON'  .or. trim(s) == 'yes' &
                 .or. trim(s) == 'YES' .or. trim(s) == 'T'  .or. trim(s) == 't')
      end if
      call incdft_load_controls()
      known = .true.
    end if
    en = value
  end function

!> @brief Optional env overrides for the reuse window (advanced/tuning).
  subroutine incdft_load_controls()
    character(len=32) :: s
    integer :: st, ln
    real(fp) :: v
    integer :: iv
    call get_environment_variable('OQP_XC_INCDFT_START', s, length=ln, status=st)
    if (st == 0 .and. ln > 0) then
      read(s, *, iostat=st) v
      if (st == 0 .and. v > 0.0_fp) incdft_start = v
    end if
    call get_environment_variable('OQP_XC_INCDFT_STOP', s, length=ln, status=st)
    if (st == 0 .and. ln > 0) then
      read(s, *, iostat=st) v
      if (st == 0 .and. v > 0.0_fp) incdft_stop = v
    end if
    call get_environment_variable('OQP_XC_INCDFT_REFRESH', s, length=ln, status=st)
    if (st == 0 .and. ln > 0) then
      read(s, *, iostat=st) iv
      if (st == 0 .and. iv >= 1) incdft_refresh = iv
    end if
  end subroutine

!> @brief Clear the reference store (call at SCF entry / geometry change).
  subroutine incdft_reset()
    call g_xc_ref%free()
    g_xc_ref%reuse_run = 0
    g_xc_ref%n_full    = 0
    g_xc_ref%n_reuse   = 0
  end subroutine

!> @brief Decide whether this iteration may reuse the stored XC matrix.
!> @param[in] diis_error  current DIIS error (proxy for closeness to convergence)
!> @param[in] iter        SCF iteration index
!> @details Reuse only inside the window [incdft_stop, incdft_start): far from
!>   convergence the XC changes too much; very close to convergence we force full
!>   builds so the fixed point is exact. A periodic forced rebuild bounds drift.
  function incdft_should_reuse(diis_error, iter) result(reuse)
    real(fp), intent(in) :: diis_error
    integer,  intent(in) :: iter
    logical :: reuse
    reuse = .false.
    if (.not. g_xc_ref%valid)             return  ! need a full reference first
    if (iter < incdft_min_iter)           return
    if (diis_error >= incdft_start)        return  ! too far from convergence
    if (diis_error <  incdft_stop)         return  ! too close -> exact final builds
    if (g_xc_ref%reuse_run >= incdft_refresh) return  ! periodic forced rebuild
    reuse = .true.
  end function

!> @brief Store the result of a full XC build as the new reference.
  subroutine incdft_store(vxc, eexc, totele, totkin)
    real(fp), intent(in) :: vxc(:,:)
    real(fp), intent(in) :: eexc, totele, totkin
    if (allocated(g_xc_ref%vxc)) then
      if (size(g_xc_ref%vxc,1) /= size(vxc,1) .or. size(g_xc_ref%vxc,2) /= size(vxc,2)) &
        deallocate(g_xc_ref%vxc)
    end if
    if (.not. allocated(g_xc_ref%vxc)) allocate(g_xc_ref%vxc(size(vxc,1), size(vxc,2)))
    g_xc_ref%vxc    = vxc
    g_xc_ref%ntri   = size(vxc,1)
    g_xc_ref%nf     = size(vxc,2)
    g_xc_ref%eexc   = eexc
    g_xc_ref%totele = totele
    g_xc_ref%totkin = totkin
    g_xc_ref%valid  = .true.
    g_xc_ref%reuse_run = 0
    g_xc_ref%n_full = g_xc_ref%n_full + 1
  end subroutine

  subroutine free(self)
    class(incdft_t), intent(inout) :: self
    if (allocated(self%vxc)) deallocate(self%vxc)
    self%valid = .false.
    self%ntri  = 0
    self%nf    = 0
  end subroutine

end module mod_dft_incdft
