!> @brief Grid ingredient export: n, |grad n|^2, tau and xi^alpha for a list of
!>   fractional orders at every grid point, for feature extraction / NN-DFT
!>   training and for comparing ingredients.
!> @details One run_xc pass with the ordinary tau collects the coordinates,
!>   weights, spin densities, sigma and tau; one further pass per alpha collects
!>   xi^alpha (the engine's tau slot).  Every pass visits the same grid points in
!>   the same per-slice order, so the results are addressed by (slice, point).
!>   Results are stored in the tagarray container:
!>     OQP::XI_GRID_ALPHA  (nalpha)           orders
!>     OQP::XI_GRID_P      (nalpha)           inner integer orders used
!>     OQP::XI_GRID_XYZW   (4*npts, point-major) x, y, z (bohr), weight
!>     OQP::XI_GRID_RHO    (2, npts)          n_alpha, n_beta
!>     OQP::XI_GRID_SIGMA  (3, npts)          sigma_aa, sigma_ab, sigma_bb
!>     OQP::XI_GRID_TAU    (2, npts)          tau_alpha, tau_beta (1/2 convention)
!>     OQP::XI_GRID_XI     (2*nalpha, npts)   xi^alpha_k alpha/beta, k = 1..nalpha
!>   Points of slices whose density is below the engine threshold are absent.
module mod_dft_gridint_xi_dump
  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t, xc_options_t, run_xc
  implicit none
  private
  public :: xi_grid_ingredients
  public :: oqp_xi_grid_ingredients

  integer, parameter :: NBASE = 11   !< xyzw(4) + rho(2) + sigma(3) + tau(2)

  type, extends(xc_consumer_t) :: xc_consumer_xi_dump_t
    integer :: nfeat = 0
    integer :: stride = 0                   !< maxSlicePts
    integer :: nslots = 0                   !< nSlices*maxSlicePts
    integer :: first = 1, last = 1          !< feature columns written by this pass
    logical :: base = .true.                !< base pass (xyzw/rho/sigma/tau) or xi pass
    integer :: k = 0                        !< xi pass index
    real(kind=fp), allocatable :: buf(:,:)  !< (nfeat, nslots)
    real(kind=fp), allocatable :: filled(:) !< (nslots) 1 where a point was written
  contains
    procedure :: parallel_start => dump_parallel_start
    procedure :: parallel_stop => dump_parallel_stop
    procedure :: update => dump_update
    procedure :: postUpdate => dump_post_update
    procedure :: clean => dump_clean
  end type

contains

  subroutine dump_parallel_start(self, xce, nthreads)
    class(xc_consumer_xi_dump_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    ! buffers are allocated once by the driver; nothing per thread (each slice
    ! owns a disjoint block of buf)
  end subroutine

  subroutine dump_parallel_stop(self)
    class(xc_consumer_xi_dump_t), intent(inout) :: self
    ! MPI: every slot is written by exactly one rank; zeros elsewhere
    call self%pe%allreduce(self%buf(self%first:self%last, :), (self%last - self%first + 1)*self%nslots)
    if (self%base) call self%pe%allreduce(self%filled, self%nslots)
  end subroutine

  subroutine dump_update(self, xce, mythread)
    class(xc_consumer_xi_dump_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    integer :: i, j0
    j0 = (xce%currSlice - 1)*self%stride
    if (self%base) then
      do i = 1, xce%numPts
        self%buf(1:4, j0 + i) = xce%xyzw(i, 1:4)
        self%buf(5:6, j0 + i) = xce%XCLib%rho(1:2, i)
        self%buf(7:9, j0 + i) = xce%XCLib%sig(1:3, i)
        self%buf(10:11, j0 + i) = xce%XCLib%tau(1:2, i)
        self%filled(j0 + i) = 1.0_fp
      end do
    else
      do i = 1, xce%numPts
        self%buf(self%first:self%first + 1, j0 + i) = xce%XCLib%tau(1:2, i)
      end do
    end if
  end subroutine

  subroutine dump_post_update(self, xce, mythread)
    class(xc_consumer_xi_dump_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
  end subroutine

  subroutine dump_clean(self)
    class(xc_consumer_xi_dump_t), intent(inout) :: self
    if (allocated(self%buf)) deallocate(self%buf)
    if (allocated(self%filled)) deallocate(self%filled)
  end subroutine

  !> C entry point: alphas(nalpha), ps(nalpha) (-1 = auto, 0, 1), scale (0/1),
  !> cutoff (bohr, 0 = exact).
  subroutine oqp_xi_grid_ingredients(c_handle, nalpha, alphas, ps, scale, cutoff) &
      bind(C, name="oqp_xi_grid_ingredients")
    use iso_c_binding, only: c_int64_t, c_double
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t), value :: nalpha
    real(c_double), intent(in) :: alphas(*)
    integer(c_int64_t), intent(in) :: ps(*)
    integer(c_int64_t), value :: scale
    real(c_double), value :: cutoff
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call xi_grid_ingredients(inf, int(nalpha), alphas(1:nalpha), int(ps(1:nalpha)), int(scale), real(cutoff, fp))
  end subroutine oqp_xi_grid_ingredients

  subroutine xi_grid_ingredients(infos, nalpha, alphas, ps, scale, cutoff)
    use types, only: information
    use mod_dft, only: dft_initialize
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: unpack_matrix
    use oqp_tagarray_driver, only: tagarray_get_data, tagarray_reserve_data, OQP_DM_A, OQP_DM_B, ta_type_real64
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nalpha
    real(kind=fp), intent(in) :: alphas(nalpha)
    integer, intent(in) :: ps(nalpha)
    integer, intent(in) :: scale
    real(kind=fp), intent(in) :: cutoff

    type(dft_grid_t), target :: molgrid
    type(xc_consumer_xi_dump_t) :: dat
    type(xc_options_t) :: xc_opts
    real(kind=fp), contiguous, pointer :: dm_a(:), dm_b(:), out(:)
    real(kind=fp), allocatable, target :: da(:,:), db(:,:)
    integer :: nbf, k, j, npts, nfeat, nslots, i, nang
    logical :: urohf

    associate(basis => infos%basis)
      nbf = basis%nbf
      urohf = infos%control%scftype /= 1   ! 1 = RHF (scf_addons::scf_rhf)
      nang = maxval(basis%am) + 1 + 1

      call dft_initialize(infos, basis, molgrid, need_functional=.false.)

      ! converged densities (packed, normalised like dmatd_density_blk)
      call tagarray_get_data(infos%dat, OQP_DM_A, dm_a)
      allocate(da(nbf, nbf), source=0.0_fp)
      call unpack_matrix(dm_a, da, nbf, "U")
      allocate(db(nbf, nbf), source=0.0_fp)
      if (urohf) then
        call tagarray_get_data(infos%dat, OQP_DM_B, dm_b)
        call unpack_matrix(dm_b, db, nbf, "U")
      else
        db = da
      end if
      do j = 1, nbf
        da(:, j) = da(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
        db(:, j) = db(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
      end do

      nfeat = NBASE + 2*nalpha
      nslots = molgrid%nSlices*molgrid%maxSlicePts
      dat%nfeat = nfeat
      dat%stride = molgrid%maxSlicePts
      dat%nslots = nslots
      allocate(dat%buf(nfeat, nslots), source=0.0_fp)
      allocate(dat%filled(nslots), source=0.0_fp)
      call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

      xc_opts%isGGA = .true.
      xc_opts%needTau = .true.
      xc_opts%functional => infos%functional
      xc_opts%hasBeta = urohf
      xc_opts%isWFVecs = .false.
      xc_opts%numAOs = nbf
      xc_opts%maxPts = molgrid%maxSlicePts
      xc_opts%limPts = molgrid%maxNRadTimesNAng
      xc_opts%numAtoms = infos%mol_prop%natom
      xc_opts%maxAngMom = nang
      xc_opts%nDer = 0
      xc_opts%numOccAlpha = infos%mol_prop%nelec_A
      xc_opts%numOccBeta = infos%mol_prop%nelec_B
      xc_opts%wfAlpha => da
      if (urohf) xc_opts%wfBeta => db
      xc_opts%molGrid => molgrid
      xc_opts%dft_threshold = infos%dft%grid_density_cutoff
      xc_opts%ao_threshold = infos%dft%grid_ao_threshold
      xc_opts%ao_sparsity_ratio = 0.0_fp
      xc_opts%use_phi_cache = .false.

      ! pass 0: ordinary ingredients
      dat%base = .true.; dat%first = 1; dat%last = NBASE
      xc_opts%xi_mode = 0
      call run_xc(xc_opts, dat, basis)

      ! passes 1..nalpha: xi^alpha in the tau slot
      do k = 1, nalpha
        dat%base = .false.; dat%k = k
        dat%first = NBASE + 2*k - 1; dat%last = NBASE + 2*k
        xc_opts%xi_mode = 1
        xc_opts%xi_alpha = alphas(k)
        xc_opts%xi_p = ps(k)
        xc_opts%xi_scale = scale
        xc_opts%xi_cutoff = cutoff
        call run_xc(xc_opts, dat, basis)
      end do

      ! compact and store
      npts = count(dat%filled > 0.5_fp)
      call infos%dat%erase((/ character(len=80) :: 'OQP::XI_GRID_ALPHA', 'OQP::XI_GRID_P', &
           'OQP::XI_GRID_XYZW', 'OQP::XI_GRID_RHO', 'OQP::XI_GRID_SIGMA', 'OQP::XI_GRID_TAU', 'OQP::XI_GRID_XI' /))
      call tagarray_reserve_data(infos%dat, 'OQP::XI_GRID_ALPHA', ta_type_real64, nalpha, (/ nalpha /))
      call tagarray_reserve_data(infos%dat, 'OQP::XI_GRID_P', ta_type_real64, nalpha, (/ nalpha /))
      call tagarray_reserve_data(infos%dat, 'OQP::XI_GRID_XYZW', ta_type_real64, 4*npts, (/ 4*npts /))
      call tagarray_reserve_data(infos%dat, 'OQP::XI_GRID_RHO', ta_type_real64, 2*npts, (/ 2*npts /))
      call tagarray_reserve_data(infos%dat, 'OQP::XI_GRID_SIGMA', ta_type_real64, 3*npts, (/ 3*npts /))
      call tagarray_reserve_data(infos%dat, 'OQP::XI_GRID_TAU', ta_type_real64, 2*npts, (/ 2*npts /))
      call tagarray_reserve_data(infos%dat, 'OQP::XI_GRID_XI', ta_type_real64, 2*nalpha*npts, (/ 2*nalpha*npts /))

      call tagarray_get_data(infos%dat, 'OQP::XI_GRID_ALPHA', out); out(1:nalpha) = alphas
      call tagarray_get_data(infos%dat, 'OQP::XI_GRID_P', out)
      do k = 1, nalpha
        if (ps(k) >= 0) then
          out(k) = real(ps(k), fp)
        else
          out(k) = real(max(0, ceiling(alphas(k))), fp)
        end if
      end do
      call store('OQP::XI_GRID_XYZW', 1, 4)
      call store('OQP::XI_GRID_RHO', 5, 2)
      call store('OQP::XI_GRID_SIGMA', 7, 3)
      call store('OQP::XI_GRID_TAU', 10, 2)
      call store('OQP::XI_GRID_XI', NBASE + 1, 2*nalpha)

      call dat%clean()
      deallocate(da, db)
    end associate

  contains

    subroutine store(tag, f0, nf)
      character(len=*), intent(in) :: tag
      integer, intent(in) :: f0, nf
      real(kind=fp), contiguous, pointer :: o(:)
      integer :: m, n
      call tagarray_get_data(infos%dat, tag, o)
      n = 0
      do m = 1, dat%nslots
        if (dat%filled(m) > 0.5_fp) then
          o(n*nf + 1:n*nf + nf) = dat%buf(f0:f0 + nf - 1, m)
          n = n + 1
        end if
      end do
    end subroutine store

  end subroutine xi_grid_ingredients

end module mod_dft_gridint_xi_dump
