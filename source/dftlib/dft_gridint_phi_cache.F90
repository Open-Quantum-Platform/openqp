!> @brief Cross-iteration cache of the collocation matrix Phi (AO values and
!>        derivatives on the quadrature grid) used by the DFT XC build.
!>
!> @details The collocation block Phi[mu,i] = phi_mu(r_i) (and its grid
!>   derivatives), together with the grid weights and the per-slice
!>   significant-AO pruning metadata, depends ONLY on geometry + basis + grid
!>   + integration thresholds -- NOT on the electron density. It is therefore
!>   identical on every SCF iteration, yet the native grid loop recomputes it
!>   (compAOs / pruneAOs) on every Fock build. This module stores the
!>   post-pruning Phi block per grid slice once per geometry and replays it on
!>   subsequent iterations, mirroring the incremental reuse OpenQP already does
!>   for the J/K (HF) part of the Fock matrix.
!>
!>   The cache is OPT-IN (env var OQP_XC_PHI_CACHE) and is only populated by the
!>   repeated SCF energy/Fock build (mod_dft_gridint_energy::dmatd_blk); one-shot
!>   consumers (gradients, response) never set the opt-in flag. Validity is keyed
!>   by a geometry hash + derivative order + DFT threshold, so any change
!>   transparently rebuilds the cache. The replayed Phi is bit-for-bit identical
!>   to the recomputed Phi, so the converged energy/gradient are unchanged.
!>
!> @author Claude (Anthropic), 2026
module mod_dft_gridint_phi_cache

  use precision, only: fp, i8b
  implicit none

  private

!> @brief Cached collocation data for a single grid slice
  type :: phi_slice_t
    logical :: stored   = .false.  !< this slice was populated during the build pass
    logical :: skip     = .false.  !< slice contributes nothing (empty or fully pruned)
    integer :: numPts   = 0        !< number of significant grid points in the slice
    integer :: numAOs_p = 0        !< effective leading AO dimension (pruned or full)
    logical :: skip_p   = .true.   !< .true. if AO pruning was skipped (dense slice)
    integer,  allocatable :: indices_p(:)  !< significant AO indices (numAOs_p)
    real(fp), allocatable :: ao(:)         !< compacted Phi block (numAOs_p*numPts*numAOVecs)
    real(fp), allocatable :: wts(:)        !< grid weights (numPts)
  end type

!> @brief Geometry-keyed cache of the collocation matrix across SCF iterations
  type, public :: phi_cache_t
    logical :: active = .false.   !< this run participates in caching
    logical :: replay = .false.   !< .true. => load cached Phi; .false. => build it
    logical :: ready  = .false.   !< cache is fully populated and valid
    integer :: nSlices   = 0
    integer :: numAOVecs = 0
    integer :: numAOs    = 0
    ! validity signature
    integer(i8b) :: sig_geom   = 0_i8b
    integer      :: sig_vecs   = 0
    integer      :: sig_natoms = 0
    real(fp)     :: sig_thr    = -1.0_fp
    integer(i8b) :: nbytes     = 0_i8b
    type(phi_slice_t), allocatable :: slices(:)
  contains
    procedure :: begin_run
    procedure :: finish_run
    procedure :: store
    procedure :: get_meta
    procedure :: get_bulk
    procedure :: free
  end type

  !> Module-global singleton: persists across run_xc calls / SCF iterations
  type(phi_cache_t), public, save :: g_phi_cache

  public :: phi_cache_geom_hash
  public :: phi_cache_env_enabled

contains

!> @brief Deterministic 64-bit hash of the molecular geometry (FNV-1a over the
!>        raw IEEE bits of every coordinate). Any coordinate change invalidates.
  pure function phi_cache_geom_hash(xyz) result(h)
    real(fp), intent(in) :: xyz(:,:)
    integer(i8b) :: h
    integer(i8b) :: bits
    integer :: i, j
    h = -3750763034362895579_i8b      ! FNV-1a 64-bit offset basis (1469598103934665603 as signed)
    do j = 1, size(xyz, 2)
      do i = 1, size(xyz, 1)
        bits = transfer(real(xyz(i,j), fp), 1_i8b)
        h = ieor(h, bits)
        h = h * 1099511628211_i8b     ! FNV-1a 64-bit prime (modular wraparound is fine)
      end do
    end do
  end function

!> @brief Is the Phi cache enabled via the environment? (cached after first read)
  function phi_cache_env_enabled() result(en)
    logical :: en
    logical, save :: known = .false.
    logical, save :: value = .false.
    character(len=32) :: s
    integer :: st, ln
    if (.not. known) then
      call get_environment_variable('OQP_XC_PHI_CACHE', s, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
        s = adjustl(s)
        value = (trim(s) == '1' .or. trim(s) == 'true' .or. trim(s) == 'TRUE' &
                 .or. trim(s) == 'on' .or. trim(s) == 'ON'  .or. trim(s) == 'yes' &
                 .or. trim(s) == 'YES' .or. trim(s) == 'T'  .or. trim(s) == 't')
      end if
      known = .true.
    end if
    en = value
  end function

!> @brief Decide build/replay mode for the current run and (re)allocate as needed.
!> @param[in] enable      master on/off (env AND opt-in caller)
!> @param[in] nSlices     number of grid slices
!> @param[in] numAOs      number of basis functions
!> @param[in] numAOVecs   AO vectors per point (1=val, 4=+grad, 10=+2nd der)
!> @param[in] natoms      number of atoms
!> @param[in] geom_hash   geometry hash (see phi_cache_geom_hash)
!> @param[in] thr         DFT integration threshold (affects nonzero-point set)
  subroutine begin_run(self, enable, nSlices, numAOs, numAOVecs, natoms, geom_hash, thr)
    class(phi_cache_t), intent(inout) :: self
    logical,      intent(in) :: enable
    integer,      intent(in) :: nSlices, numAOs, numAOVecs, natoms
    integer(i8b), intent(in) :: geom_hash
    real(fp),     intent(in) :: thr
    logical :: match

    self%active = enable
    self%replay = .false.
    if (.not. enable) return

    match = self%ready &
      .and. self%nSlices   == nSlices  &
      .and. self%numAOs    == numAOs   &
      .and. self%numAOVecs == numAOVecs &
      .and. self%sig_natoms == natoms  &
      .and. self%sig_geom  == geom_hash &
      .and. self%sig_thr   == thr

    if (match) then
      self%replay = .true.
      return
    end if

    ! (Re)build: drop any stale cache and prepare fresh per-slice storage.
    call self%free()
    self%nSlices   = nSlices
    self%numAOs    = numAOs
    self%numAOVecs = numAOVecs
    self%sig_natoms = natoms
    self%sig_geom  = geom_hash
    self%sig_thr   = thr
    allocate(self%slices(nSlices))
    self%ready  = .false.
    self%replay = .false.
  end subroutine

!> @brief Finalize a build pass: mark the cache ready and tally its footprint.
  subroutine finish_run(self)
    class(phi_cache_t), intent(inout) :: self
    integer :: i
    if (.not. self%active) return
    if (self%replay) return
    self%nbytes = 0_i8b
    do i = 1, self%nSlices
      if (allocated(self%slices(i)%ao))  self%nbytes = self%nbytes + 8_i8b*size(self%slices(i)%ao, kind=i8b)
      if (allocated(self%slices(i)%wts)) self%nbytes = self%nbytes + 8_i8b*size(self%slices(i)%wts, kind=i8b)
    end do
    self%ready = .true.
  end subroutine

!> @brief Store one slice's post-pruning collocation block (build pass).
!> @param[in] iSlice     slice index
!> @param[in] skip       whether the slice contributes nothing
!> @param[in] numPts     significant points in the slice
!> @param[in] numAOs_p   effective leading AO dimension
!> @param[in] skip_p     whether pruning was skipped
!> @param[in] indices_p  significant AO indices (length >= numAOs_p)
!> @param[in] ao         flat AO memory (xce%aoMem_); first numAOs_p*numPts*numAOVecs used
!> @param[in] wts        grid weights (length >= numPts)
  subroutine store(self, iSlice, skip, numPts, numAOs_p, skip_p, indices_p, ao, wts)
    class(phi_cache_t), intent(inout) :: self
    integer,  intent(in) :: iSlice, numPts, numAOs_p
    logical,  intent(in) :: skip, skip_p
    integer,  intent(in) :: indices_p(:)
    real(fp), intent(in) :: ao(:)
    real(fp), intent(in) :: wts(:)
    integer :: n

    associate (s => self%slices(iSlice))
      s%stored = .true.
      s%skip   = skip
      if (skip) return
      s%numPts   = numPts
      s%numAOs_p = numAOs_p
      s%skip_p   = skip_p
      n = numAOs_p*numPts*self%numAOVecs
      ! Significant-AO indices are only consumed by the pruned (sparse) scatter
      ! path; the dense path adds the full block, so skip storing them there.
      if (.not. skip_p) s%indices_p = indices_p(1:numAOs_p)
      s%ao  = ao(1:n)
      s%wts = wts(1:numPts)
    end associate
  end subroutine

!> @brief Read a slice's scalar metadata (replay pass, phase 1).
  subroutine get_meta(self, iSlice, skip, numPts, numAOs_p, skip_p)
    class(phi_cache_t), intent(in) :: self
    integer, intent(in)  :: iSlice
    logical, intent(out) :: skip, skip_p
    integer, intent(out) :: numPts, numAOs_p
    associate (s => self%slices(iSlice))
      if (.not. s%stored) then       ! defensive: treat unpopulated slice as empty
        skip = .true.; numPts = 0; numAOs_p = 0; skip_p = .true.
        return
      end if
      skip     = s%skip
      numPts   = s%numPts
      numAOs_p = s%numAOs_p
      skip_p   = s%skip_p
    end associate
  end subroutine

!> @brief Copy a slice's bulk data back into the engine arrays (replay, phase 2).
!> @param[out] indices_p  receives significant AO indices (first numAOs_p)
!> @param[out] ao         receives the Phi block into xce%aoMem_ (first n)
!> @param[out] wts        receives grid weights into xce%xyzw(:,4) (first numPts)
  subroutine get_bulk(self, iSlice, indices_p, ao, wts)
    class(phi_cache_t), intent(in) :: self
    integer,  intent(in)  :: iSlice
    integer,  intent(out) :: indices_p(:)
    real(fp), intent(out) :: ao(:)
    real(fp), intent(out) :: wts(:)
    integer :: n, np
    associate (s => self%slices(iSlice))
      np = s%numAOs_p
      n  = np*s%numPts*self%numAOVecs
      if (.not. s%skip_p) indices_p(1:np) = s%indices_p(1:np)
      ao(1:n)             = s%ao(1:n)
      wts(1:s%numPts)     = s%wts(1:s%numPts)
    end associate
  end subroutine

!> @brief Release all cached storage and reset to the empty/invalid state.
  subroutine free(self)
    class(phi_cache_t), intent(inout) :: self
    if (allocated(self%slices)) deallocate(self%slices)
    self%ready   = .false.
    self%replay  = .false.
    self%nSlices = 0
    self%nbytes  = 0_i8b
  end subroutine

end module mod_dft_gridint_phi_cache
