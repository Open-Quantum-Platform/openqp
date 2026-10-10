!> @brief Grid integration of the superposition-of-atomic-potentials (SAP)
!>        one-electron operator.
!>
!> Builds the SAP potential matrix
!>     V_{mu,nu} = sum_g w_g phi_mu(r_g) V_SAP(r_g) phi_nu(r_g),
!>     V_SAP(r) = sum_A V_A(|r - R_A|) = sum_A -Z_eff^A(|r - R_A|)/|r - R_A|,
!> on the existing DFT molecular grid, reusing the AO-on-grid evaluation and
!> AO-pruning machinery of mod_dft_gridint (run_grid_aos). The construction
!> mirrors the Kohn-Sham matrix accumulation in mod_dft_gridint_energy.
!>
!> Reference: S. Lehtola, "Assessment of Initial Guesses for Self-Consistent
!> Field Calculations. Superposition of Atomic Potentials: Simple yet
!> Efficient", J. Chem. Theory Comput. 15, 1593 (2019).
module mod_dft_gridint_sap

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t
  use sap_lut, only: sap_table_t
  use oqp_linalg

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_sap_t
    real(kind=fp), allocatable :: va2(:,:)     !< full SAP matrix accumulator
    real(kind=fp), allocatable :: vmat_(:,:)   !< per-thread pruned block
    real(kind=fp), allocatable :: tmp_(:,:)    !< per-thread scratch
    ! SAP inputs (read-only during the grid loop)
    real(kind=fp), pointer :: atxyz(:,:) => null()  !< atom coordinates (3,nat), bohr
    integer, allocatable :: atz(:)                  !< nuclear charge per atom
    type(sap_table_t), pointer :: sap => null()     !< radial Z_eff table
  contains
    procedure :: parallel_start
    procedure :: parallel_stop
    procedure :: resetOrbPointers
    procedure :: update
    procedure :: postUpdate
    procedure :: clean
  end type

!-------------------------------------------------------------------------------

  private
  public sap_potential_matrix

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

  subroutine parallel_start(self, xce, nthreads)
    implicit none
    class(xc_consumer_sap_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    ! Note: do NOT call self%clean() here -- it would wipe the SAP inputs
    ! (atxyz/atz/sap) set by the driver before the parallel region.
    if (allocated(self%va2)) deallocate(self%va2)
    if (allocated(self%vmat_)) deallocate(self%vmat_)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
    allocate( self%va2(xce%numAOs*xce%numAOs, nthreads) &
            , self%vmat_(xce%numAOs*xce%numAOs, nthreads) &
            , self%tmp_(xce%numAOs*xce%maxPts, nthreads) &
            , source=0.0d0)
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_sap_t), intent(inout) :: self
    if (ubound(self%va2,2) /= 1) then
      self%va2(:,lbound(self%va2,2)) = sum(self%va2, dim=2)
    end if
    call self%pe%allreduce(self%va2(:,1), size(self%va2(:,1)))
  end subroutine

!-------------------------------------------------------------------------------

  subroutine clean(self)
    implicit none
    class(xc_consumer_sap_t), intent(inout) :: self
    if (allocated(self%va2)) deallocate(self%va2)
    if (allocated(self%vmat_)) deallocate(self%vmat_)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
    if (allocated(self%atz)) deallocate(self%atz)
  end subroutine

!-------------------------------------------------------------------------------

  subroutine resetOrbPointers(self, xce, vmat, tmp, va, myThread)
    class(xc_consumer_sap_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer, optional :: vmat(:,:)
    real(kind=fp), intent(out), pointer, optional :: tmp(:,:)
    real(kind=fp), intent(out), pointer, optional :: va(:,:)
    integer, intent(in) :: myThread

    associate ( numAOs   => xce%numAOs &
              , numAOs_p => xce%numAOs_p &
              , numPts   => xce%numPts )
      if (present(vmat)) &
        vmat(1:numAOs_p, 1:numAOs_p) => self%vmat_(1:numAOs_p*numAOs_p, myThread)
      if (present(tmp)) &
        tmp(1:numAOs_p, 1:numPts) => self%tmp_(1:numAOs_p*numPts, myThread)
      if (present(va)) &
        va(1:numAOs, 1:numAOs) => self%va2(1:numAOs*numAOs, myThread)
    end associate
  end subroutine

!-------------------------------------------------------------------------------

  subroutine update(self, xce, mythread)
    implicit none
    class(xc_consumer_sap_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    integer :: i, iat
    real(kind=fp) :: vsap, dist, dx, dy, dz
    real(kind=fp), pointer :: vmat(:,:)
    real(kind=fp), pointer :: tmp(:,:)

    call self%resetOrbPointers(xce, vmat=vmat, tmp=tmp, myThread=myThread)

    associate ( aoV    => xce%aoV &
              , xyzw   => xce%xyzw &
              , wts    => xce%wts &
              , numAOs => xce%numAOs_p &
              , numPts => xce%numPts )

      do i = 1, numPts
        vsap = 0.0_fp
        do iat = 1, size(self%atz)
          dx = xyzw(i,1) - self%atxyz(1,iat)
          dy = xyzw(i,2) - self%atxyz(2,iat)
          dz = xyzw(i,3) - self%atxyz(3,iat)
          dist = sqrt(dx*dx + dy*dy + dz*dz)
          vsap = vsap + self%sap%potential(self%atz(iat), dist)
        end do
        ! factor 0.5: dsyr2k forms aoV*tmp^T + tmp*aoV^T
        tmp(:, i) = 0.5_fp * wts(i) * vsap * aoV(:, i)
      end do

      call dsyr2k('U', 'N', numAOs, numPts, 1.0_fp, &
                          aoV, numAOs, &
                          tmp, numAOs, &
                  0.0_fp, vmat, numAOs)

    end associate
  end subroutine

!-------------------------------------------------------------------------------

  subroutine postUpdate(self, xce, mythread)
    implicit none
    class(xc_consumer_sap_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(kind=fp), pointer :: vmat(:,:)
    real(kind=fp), pointer :: va(:,:)

    call self%resetOrbPointers(xce, vmat=vmat, va=va, myThread=myThread)

    associate( numAOs  => xce%numAOs_p &
             , indices => xce%indices_p )
      if (xce%skip_p) then
        va = va + vmat
      else
        va(indices(1:numAOs), indices(1:numAOs)) = &
          va(indices(1:numAOs), indices(1:numAOs)) + vmat
      end if
    end associate
  end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute the SAP potential matrix (packed lower-triangular).
!> @param[in]    basis          basis set
!> @param[in]    molGrid        molecular grid
!> @param[out]   vsap_tri       packed (lower-triangular) SAP matrix, nbf*(nbf+1)/2
!> @param[in]    nbf            number of basis functions
!> @param[in]    sap            radial Z_eff table
!> @param[in]    infos          calculation information
  subroutine sap_potential_matrix(basis, molGrid, vsap_tri, nbf, sap, infos)
    use mod_dft_molgrid, only: dft_grid_t
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_grid_aos
    use types, only: information

    implicit none

    type(dft_grid_t), target, intent(in) :: molGrid
    type(information), target, intent(in) :: infos
    type(basis_set) :: basis
    integer, intent(in) :: nbf
    real(kind=fp), intent(out) :: vsap_tri(*)
    type(sap_table_t), target, intent(in) :: sap

    type(xc_consumer_sap_t) :: dat
    type(xc_options_t) :: xc_opts
    integer :: i0, i, maxl, nang, nat

    nat = infos%mol_prop%natom

    maxl = maxval(basis%am)
    nang = maxl + 1 + 1

    xc_opts%isGGA = .false.
    xc_opts%needTau = .false.
    xc_opts%hasBeta = .false.
    xc_opts%isWFVecs = .false.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = nat
    xc_opts%maxAngMom = nang
    xc_opts%nDer = 0
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%molGrid => molGrid
    xc_opts%dft_threshold = infos%dft%grid_density_cutoff
    xc_opts%ao_threshold = infos%dft%grid_ao_threshold
    ! Disable AO-space pruning: the engine's pruned-AO path compresses the
    ! wavefunction (wfAlpha), which the SAP guess does not provide. Keeping all
    ! AOs (skip_p always true) avoids that path at a small performance cost.
    xc_opts%ao_sparsity_ratio = 0.0_fp

    ! SAP inputs
    dat%atxyz => infos%atoms%xyz
    dat%sap => sap
    allocate(dat%atz(nat))
    do i = 1, nat
      dat%atz(i) = nint(infos%atoms%zn(i))
    end do

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_grid_aos(xc_opts, dat, basis)

    ! Pack accumulated (raw-normalized) matrix into lower triangle, applying
    ! the basis-function normalization, exactly as in dmatd_blk.
    i0 = 0
    do i = 1, nbf
      vsap_tri(i0+1:i0+i) = dat%va2((i-1)*nbf+1:(i-1)*nbf+i, 1) &
                            * basis%bfnrm(i) * basis%bfnrm(1:i)
      i0 = i0 + i
    end do

    call dat%clean()

  end subroutine

!-------------------------------------------------------------------------------

end module mod_dft_gridint_sap
