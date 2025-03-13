module mod_dft_gridint

  use precision, only: fp
!  use params, only: dft_wt_der
  use basis_tools, only: basis_set
  use io_constants, only: iw
  use mod_dft_xc_libxc, only: xc_libxc_t
  use mod_dft_molgrid, only: dft_grid_t
  use functionals, only: functional_t
  use oqp_linalg
  use parallel, only: par_env_t
  implicit none

!###############################################################################

  integer, parameter, public :: &
      OQP_FUNTYP_LDA  = 0, &
      OQP_FUNTYP_GGA  = 1, &
      OQP_FUNTYP_MGGA = 2

  integer, parameter, public :: &
      X__ = 1, Y__ = 2, Z__ = 3

  integer, parameter, public :: &
      XX_ = 1, YY_ = 2, ZZ_ = 3, &
      XY_ = 4, YZ_ = 5, XZ_ = 6

  integer, parameter, public :: &
      dXX = 1, dXY = 2, dXZ = 3, &
      dYX = 4, dYY = 5, dYZ = 6, &
      dZX = 7, dZY = 8, dZZ = 9

  integer, parameter, public :: &
      XXX = 1, YYY = 2, ZZZ = 3, &
      XXY = 4, XXZ = 5, YYX = 6, YYZ = 7, &
      ZZX = 8, ZZY = 9, XYZ = 10

! Convert 'square' XYZ indices to triangular,
! needed in hessian code
  integer, parameter, public :: &
      SQ_TO_TR(3,3) = reshape( [ &
                      XX_, XY_, XZ_ &
                    , XY_, YY_, YZ_ &
                    , XZ_, YZ_, ZZ_ &
                    ], shape(SQ_TO_TR) )

!###############################################################################

!> @brief Basic type to consume XC values on a grid
  type, abstract :: xc_consumer_t
    real(kind=fp) :: E_xc
    real(kind=fp) :: E_exch
    real(kind=fp) :: E_corr
    real(kind=fp) :: N_elec
    real(kind=fp) :: E_kin
    real(kind=fp) :: G_total(3)
    type(par_env_t) :: pe
  contains
    procedure(xc_consumer_parallel_start), deferred, pass :: parallel_start
    procedure(xc_consumer_parallel_stop), deferred, pass :: parallel_stop
    procedure(xc_consumer_update), deferred, pass :: update
    procedure(xc_consumer_postUpdate), deferred, pass :: postUpdate
    procedure(xc_consumer_clean), deferred, pass :: clean
  end type

!###############################################################################

!> @brief Interface structure to set up XC engine options
  type :: xc_options_t
    logical :: isGGA = .false.
    logical :: needTau = .false.
    logical :: hasBeta = .false.
    !< .T./.F.  - wfA and wfB are MO vectors/densities
    logical :: isWFVecs = .true.

    integer :: numAOs = 0
    integer :: maxPts = 0
    integer :: limPts = 0
    integer :: numAtoms = 0
    integer :: maxAngMom = 0
    integer :: nDer = 0
    integer :: nXCDer = 0
    integer :: numAOVecs = 0
    integer :: numTmpVec = 0
    integer :: numOccAlpha = 0
    integer :: numOccBeta = 0

    real(kind=fp) :: dft_threshold = 0.0d0
    real(kind=fp) :: ao_threshold = 0.0d0
    real(kind=fp) :: ao_sparsity_ratio = 0.0d0

    !< alpha spin wavefunction
    real(KIND=fp), contiguous, pointer :: wfAlpha(:, :) => null()
    !< beta spin wavefunction
    real(KIND=fp), contiguous, pointer :: wfBeta(:, :) => null()

    !< Molecular grid data
    type(dft_grid_t), pointer :: molGrid => null()
    type(functional_t), pointer :: functional

  end type

!###############################################################################

!> @brief Main class which knows how to compute XC functional values, AO and MO
!>   values and gradients on a grid
!> @details It is complemented with xc_consumer_t class to use calculation results
  type :: xc_engine_t
!    private
    real(KIND=fp), allocatable :: xyzw(:, :)
    real(KIND=fp), allocatable :: aoMem_(:) !< AO memory
    real(KIND=fp), allocatable :: moMemA_(:) !< MO memory (alpha)
    real(KIND=fp), allocatable :: moMemB_(:) !< MO memory (beta)
    real(KIND=fp), allocatable :: tmpWfAlpha(:) !< tmp alpha spin wavefunction
    real(KIND=fp), allocatable :: tmpWfBeta(:) !< tmp beta spin wavefunction
    integer, allocatable :: indices_p(:) !< AO significant indices

    real(KIND=fp), contiguous, pointer :: &
        aoMem(:, :, :) => null() & !< AO memory
      , moMemA(:, :, :) => null() & !< MO memory (alpha)
      , moMemB(:, :, :) => null() & !< MO memory (beta)
      , aoV(:, :) => null() & !< AO values
      , moVA(:, :) => null() & !< MO values (alpha)
      , moVB(:, :) => null() & !< MO values (beta)
      , aoG1(:, :, :) => null() & !< AO gradient
      , aoG2(:, :, :) => null() & !< AO 2nd der.
      , moG1A(:, :, :) => null() & !< MO gradient (alpha)
      , moG2A(:, :, :) => null() & !< MO 2nd der. (alpha)
      , moG1B(:, :, :) => null() & !< MO gradient (beta)
      , moG2B(:, :, :) => null() & !< MO 2nd der. (beta)
      , wts(:) => null() & !< weights
      , wfAlpha(:, :) => null() & !< alpha spin wavefunction
      , wfBeta(:, :) => null() & !< beta spin wavefunction
      , wfAlpha_p(:, :) => null() & !< pruned alpha spin wavefunction
      , wfBeta_p(:, :) => null()   !< pruned beta spin wavefunction

    logical :: isGGA = .false.
    logical :: needTau = .false.
    logical :: hasBeta = .false.
    logical :: isWFVecs = .true.   !< .TRUE.  - wfA and wfB are MO vectors
                                   !< .FALSE. - wfA and wfB are densities
    integer :: numAOs = 0 !< number of AOs
    integer :: numAOs_p = 0 !< number of pruned AOs
    logical :: skip_p = .true. !< skip if no pruned numAOs
    integer :: numPts = 0
    integer :: numAtoms = 0
    integer :: maxPts = 0
    integer :: maxAngMom = 0
    integer :: nAODer = 0
    integer :: nXCDer = 1
    integer :: funTyp = 0   !< 0 - LDA, 1 - GGA, 2 - MGGA

    integer :: numAOVecs = 0
    integer :: numTmpVec = 0

    integer :: numOccAlpha = 0
    integer :: numOccBeta = 0

    real(kind=fp) :: threshold  = 1.0d-15
    real(kind=fp) :: ao_threshold  = 1.0d-15
    real(kind=fp) :: ao_sparsity_ratio = 0.90d+0 !< Cut off if more than 90% pruned AOs
                                                 !< (skip_p becomes False).

    type(xc_libxc_t), allocatable :: XCLib

    integer       :: dbgLevel   = 0
    real(kind=fp) :: N_elec     = 0.0
    real(kind=fp) :: E_kin      = 0.0
    real(kind=fp) :: G_total(3) = 0.0

    procedure(compute_density), pointer, pass :: compRho => null()
    procedure(compute_density_grad), pointer, pass :: compDRho => null()
    procedure(compute_density_tau), pointer, pass :: compTau => null()

  contains
    procedure :: init
    procedure :: echo => echoVars
    procedure :: getStats
    procedure :: resetPointers
    procedure :: resetOrbPointers
    procedure :: resetXCPointers
    procedure :: compAOs
    procedure :: pruneAOs
    procedure :: resetPrunedPointers
    procedure :: compMOs
    procedure :: compRMOs
    procedure :: compRMOGs

    generic :: compRRho  => compRRho_ab, compRRho_a
    generic :: compRDRho => compRDRho_ab, compRDRho_a
    generic :: compRTau => compRTau_ab, compRTau_a

    procedure :: compRhoAll

    procedure :: compXC

    procedure, private :: compRRho_ab
    procedure, private :: compRDRho_ab
    procedure, private :: compRTau_ab
    procedure, private :: compRRho_a
    procedure, private :: compRDRho_a
    procedure, private :: compRTau_a

  end type xc_engine_t

!###############################################################################

  abstract interface
!> @brief Initialization of the data for XC consumer
!> @note This class should handle multithreaded runs by its own means
!> @note This subroutine is executed inside the parallel region by master thread only
    subroutine xc_consumer_parallel_start(self, xce, nthreads)
      import :: xc_consumer_t, xc_engine_t
      implicit none
      class(xc_consumer_t), target, intent(inout) :: self
      class(xc_engine_t), intent(in) :: xce
      integer, intent(in) :: nthreads
    end subroutine
!> @brief Finalization of data in parallel run
!> @note This subroutine is executed outside of the parallel region
    subroutine xc_consumer_parallel_stop(self)
      import :: xc_consumer_t
      implicit none
      class(xc_consumer_t), intent(inout) :: self
    end subroutine
!> @brief Release resources of xc_consumer_t
!> @note This subroutine is executed outside of the parallel region
    subroutine xc_consumer_clean(self)
      import :: xc_consumer_t
      implicit none
      class(xc_consumer_t), intent(inout) :: self
    end subroutine
!> @brief Main subroutine to consume XC functional values provided by xc_engine_t
!> @note This subroutine is executed inside the parallel region by every thread
    subroutine xc_consumer_update(self, xce, mythread)
      import :: xc_consumer_t, xc_engine_t
      implicit none
      class(xc_consumer_t), intent(inout) :: self
      class(xc_engine_t), intent(in) :: xce
      integer :: mythread
    end subroutine
!> @note This subroutine is executed inside the parallel region by every thread
    subroutine xc_consumer_postUpdate(self, xce, mythread)
      import :: xc_consumer_t, xc_engine_t
      implicit none
      class(xc_consumer_t), intent(inout) :: self
      class(xc_engine_t), intent(in) :: xce
      integer :: mythread
    end subroutine
  end interface

!###############################################################################

  abstract interface
    subroutine compute_density(self, rho)
        import
        class(xc_engine_t) :: self
        real(kind=fp), intent(out) :: rho(:,:)
    end subroutine
    subroutine compute_density_grad(self, drho, sigma)
        import
        class(xc_engine_t) :: self
        real(kind=fp), intent(out) :: drho(:,:), sigma(:,:)
    end subroutine
    subroutine compute_density_tau(self, tau)
        import
        class(xc_engine_t) :: self
        real(kind=fp), intent(out) :: tau(:,:)
    end subroutine
  end interface

!###############################################################################

  private
  public xc_engine_t
  public xc_consumer_t
  public xc_options_t
  public run_xc
  public mo_tran_symm_
  public mo_tran_gemm_
  public xc_der1
  public xc_der2_contr
  public xc_der3_contr

  public compAtGradRho
  public compAtGradDRho
  public compAtGradTau

contains

!###############################################################################
!###############################################################################

 subroutine mo_tran_symm_(numAOs, nVecs, nPts, A, B, C)
    integer, intent(in) :: numAOs, nVecs, nPts
    real(kind=fp), intent(in) :: A(*), B(*)
    real(kind=fp), intent(inout) :: C(*)
    call dsymm('L', 'U', &
        numAOs, nVecs*nPts, &
        1.0_fp, A, numAOs, &
                B, numAOs, &
        0.0_fp, C, numAOs)
 end subroutine

!###############################################################################

 subroutine mo_tran_gemm_(numMOs, numAOs, nVecs, nPts, numAOs_active, A, B, C)
    integer, intent(in) :: numMOs, numAOs, nVecs, nPts, numAOs_active
    real(kind=fp), intent(in) :: A(*), B(*)
    real(kind=fp), intent(inout) :: C(*)
    call dgemm('T', 'N', &
        numMOs, nVecs*nPts, numAOs, &
        1.0_fp, a, numAOs, &
                b, numAOs, &
        0.0_fp, c, numAOs_active)
 end subroutine

!###############################################################################

!> @brief Scale 2d array along 1st dimension by a given
!>  vector of weights
 subroutine scale_2d(array, weights)
    real(kind=fp), intent(inout) :: array(:,:)
    real(kind=fp), intent(in) :: weights(:)
    integer :: i
    do i = lbound(array,2), ubound(array,2)
        array(:,i) = array(:,i) * weights
    end do
 end subroutine

!###############################################################################

!> @brief Print parameters of the xc_engine_t instance
!> @author Vladimir Mironov
  subroutine echoVars(self)
    class(xc_engine_t) :: self
    real(kind=fp) :: exc, ex, ec
    call self%XCLib%getEnergy(exc, ex, ec)

    write (*, *) 'isGGA=', self%isGGA
    write (*, *) 'needTau =', self%needTau
    write (*, *) 'hasBeta =', self%hasBeta
    write (*, *) 'numAOs  =', self%numAOs
    write (*, *) 'numPts  =', self%numPts
    write (*, *) 'nAODer  =', self%nAODer
    write (*, *) 'nXCDer  =', self%nXCDer
    write (*, *) 'numOccA =', self%numOccAlpha

    write (*, *) 'N_elec  =', self%N_elec
    write (*, *) 'E_kin   =', self%E_kin
    write (*, *) 'G_total =', self%G_total
    write (*, *) 'E_xc    =', exc
  end subroutine

!###############################################################################

!> @brief Get debug statistics
!> @author Vladimir Mironov
 subroutine getStats(self, E_xc, E_exch, E_corr, N_elec, E_kin, G_total)
    class(xc_engine_t) :: self
    real(kind=fp), optional, intent(out) :: &
        E_xc, E_exch, E_corr, N_elec, E_kin, G_total(3)
    real(kind=fp) :: exc, ex, ec

    call self%XCLib%getEnergy(exc, ex, ec)

    if (present(E_xc   )) E_xc       = exc
    if (present(E_exch )) E_exch     = ex
    if (present(E_corr )) E_corr     = ec
    if (present(N_elec )) N_elec     = self%N_elec
    if (present(E_kin  )) E_kin      = self%E_kin
    if (present(G_total)) G_total(3) = self%G_total(3)

 end subroutine

!###############################################################################

!> @brief Initialize xc_engine_t instance
!> @param[in] numAOs    number of atomic orbitals in a basis
!> @param[in] nAt       number of atoms in a system
!> @param[in] maxAngMom maximum angular momentum of basis functions and their derivatives
!> @param[in] maxPts    maximum known number of non-zero points in a slice
!> @param[in] limPts    maximum possible number of points (i.e. max(nRad*nAng)) in a slice
!> @param[in] nDer      degree of energy derivative needed
!> @param[in] hasBeta   .TRUE. if open-shell calculation
!> @param[in] isGGA     .TRUE. if GGA/metaGGA functional
!> @param[in] needTau   .TRUE. if metaGGA functional
!> @param[in] vec_or_dens .TRUE./.FALSE. - wavefunction is MO vectors/density
!> @param[in] nOccAlpha number of occupied orbitals, alpha spin
!> @param[in] nOccBeta  number of occupied orbitals, beta spin
!> @param[in] wfAlpha   wavefunction, alpha spin
!> @param[in] wfBeta    wavefunction, beta spin
!> @author Vladimir Mironov
  subroutine init(self, xco)
    implicit none
    class(xc_engine_t), target, intent(inout) :: self
    type(xc_options_t), target, intent(in) :: xco
!   Will be possibly needed to use LibXC:
    integer, parameter :: nAOVecs(0:3) = [1, 4, 10, 20]

    logical :: reqSigma

    self%funTyp = OQP_FUNTYP_LDA
    if (xco%isGGA) self%funTyp = OQP_FUNTYP_GGA
    if (xco%needTau) self%funTyp = OQP_FUNTYP_MGGA

    self%nAODer = xco%nDer
    if (self%funTyp /= OQP_FUNTYP_LDA) self%nAODer = self%nAODer + 1

    self%nXCDer = max(1, xco%nXCDer) ! at least 1st derivative

!   Find out the amount of memory needed
    if (self%nAODer<0 .or. self%nAODer>3) then
      write (*, *) 'Invalid grad level in xc_engine_t % INIT'
      stop
    end if

    self%numAOVecs = nAOVecs(self%nAODer)

    self%numTmpVec = 1
    if (xco%needTau) self%numTmpVec = 4

!   Allocate memory for XC calculations
    allocate ( &
      self%aoMem_(xco%numAOs*self%numAOVecs*xco%maxPts), &
      self%moMemA_(xco%numAOs*self%numAOVecs*xco%maxPts), &
      self%tmpWfAlpha(xco%numAOs*xco%numAOs), &
      self%tmpWfBeta(xco%numAOs*xco%numAOs), &
!     Allocate memory for grid points storage
      self%xyzw(xco%limPts, 4) &
    )
    if (xco%hasBeta) allocate(self%moMemB_(xco%numAOs*self%numAOVecs*xco%maxPts))
!   Allocate memory for AO significant indicis
    allocate (self%indices_p(xco%numAOs)) !< AO significant indices

    self%maxAngMom = xco%maxAngMom+self%nAODer

!   Set up other runtime options
    self%numAOs = xco%numAOs
    self%numAtoms = xco%numAtoms
    self%maxPts = xco%maxPts

    self%isGGA = xco%isGGA
    self%needTau = xco%needTau
    self%hasBeta = xco%hasBeta

!   Manage density/MO vectors
    self%isWFVecs = xco%isWFVecs
    self%wfAlpha => xco%wfAlpha

    self%numOccAlpha = xco%numOccAlpha

    if (self%isWFVecs) then
      self%compRho => compRhoMO
      self%compDRho => compDRhoMO
      self%compTau => compTauMO
    else
      self%compRho => compRhoAO
      self%compDRho => compDRhoAO
      self%compTau => compTauAO
    end if

    if (self%hasBeta) then
      self%numOccBeta = xco%numOccBeta
      self%wfBeta => xco%wfBeta
    else
      self%numOccBeta = xco%numOccAlpha
      self%wfBeta => xco%wfAlpha
    end if

    self%threshold = xco%dft_threshold
    self%ao_threshold  = xco%ao_threshold
    self%ao_sparsity_ratio = xco%ao_sparsity_ratio

!   Initialize XC library
    allocate(self%XCLib)
    reqSigma = self%funTyp /= OQP_FUNTYP_LDA
    call self%XCLib%init(reqSigma, self%needTau, .false., self%hasBeta, self%maxPts, self%nXCDer)

  end subroutine

!> @brief Adjust internal memory storage for a given
!>  number of grid points
!> @param[in] numPts    number of grid points
!> @author Vladimir Mironov
 subroutine resetPointers(self, numPts)
    class(xc_engine_t) :: self
    integer, intent(in) :: numPts

    self%numPts = numPts

    call self%resetOrbPointers
    call self%resetXCPointers
    call self%XCLib%setPts(numPts)

 end subroutine

!###############################################################################

!> @brief Adjust XC memory storage for a given
!>  number of grid points
!> @author Vladimir Mironov
 subroutine resetXCPointers(self)
   class(xc_engine_t), target :: self

   associate( numPts    => self%numPts &
        )

     self%wts(1:numPts) => self%xyzw(1:numPts,4)

   end associate

 end subroutine

!###############################################################################

!> @brief Adjust internal AO/MO memory storage for a given
!>  number of grid points
!> @author Vladimir Mironov
  subroutine resetOrbPointers(self)
    class(xc_engine_t), target :: self

    associate(  numAOs    => self%numAOs &
              , numPts    => self%numPts &
              , numAOVecs => self%numAOVecs &
        )

      self%aoMem(1:numAOs, 1:numPts, 1:numAOVecs) => self%aoMem_(1:)
      self%moMemA(1:numAOs, 1:numPts, 1:numAOVecs) => self%moMemA_(1:)

      if (self%hasBeta) then
        self%moMemB(1:numAOs, 1:numPts, 1:numAOVecs) => self%moMemB_(1:)
      end if

    end associate

    select case (self%nAODer)
    case (0)
      self%aoV => self%aoMem(:, :, 1)
      self%moVA => self%moMemA(:, :, 1)

      if (self%hasBeta) then
        self%moVB => self%moMemB(:, :, 1)
      end if
    case (1)
      self%aoV => self%aoMem(:, :, 1)
      self%aoG1 => self%aoMem(:, :, 2:4)

      self%moVA => self%moMemA(:, :, 1)
      self%moG1A => self%moMemA(:, :, 2:4)

      if (self%hasBeta) then
        self%moVB => self%moMemB(:, :, 1)
        self%moG1B => self%moMemB(:, :, 2:4)
      end if
    case (2)
      self%aoV => self%aoMem(:, :, 1)
      self%aoG1 => self%aoMem(:, :, 2:4)
      self%aoG2 => self%aoMem(:, :, 5:10)

      self%moVA => self%moMemA(:, :, 1)
      self%moG1A => self%moMemA(:, :, 2:4)
      self%moG2A => self%moMemA(:, :, 5:10)

      if (self%hasBeta) then
        self%moVB => self%moMemB(:, :, 1)
        self%moG1B => self%moMemB(:, :, 2:4)
        self%moG2B => self%moMemB(:, :, 5:10)
      end if
    end select

  end subroutine

!###############################################################################

!> @brief Compute atomic orbital values/gradient/hessian in a grid point
!> @param[in]  iPtIn      index of the point in self%xyzw array
!> @param[in]  iPtOut     index of the point in AO/MO arrays
!> @param[out] nnz        number of non-zero AOs in the point
!> @author Vladimir Mironov
  subroutine compAOs(self, basis, nDer, xyz)
    class(xc_engine_t) :: self
    type(basis_set),intent(in) :: basis
    integer, intent(in) :: nDer
    real(kind=fp), intent(in) :: xyz(:,:)

    integer :: nnz, ipt
    real(kind=fp) :: ptxyz(3)

    select case (nDer)
    case (0)
      do iPt = 1, ubound(xyz,1)
        ptxyz = xyz(iPt,1:3)

        call basis%aoval(ptxyz, nnz, &
                     self%aoV(:, iPt))
      end do
    case (1)
      do iPt = 1, ubound(xyz,1)
        ptxyz = xyz(iPt,1:3)

        call basis%aoval(ptxyz, nnz, &
                      self%aoV(:, iPt), &
                      self%aoG1(:, iPt, X__), &
                      self%aoG1(:, iPt, Y__), &
                      self%aoG1(:, iPt, Z__))
      end do
    case (2)
      do iPt = 1, ubound(xyz,1)
        ptxyz = xyz(iPt,1:3)

        call basis%aoval(ptxyz, nnz, &
                       self%aoV(:, iPt), &
                       self%aoG1(:, iPt, X__), &
                       self%aoG1(:, iPt, Y__), &
                       self%aoG1(:, iPt, Z__), &
                       self%aoG2(:, iPt, XX_), &
                       self%aoG2(:, iPt, YY_), &
                       self%aoG2(:, iPt, ZZ_), &
                       self%aoG2(:, iPt, XY_), &
                       self%aoG2(:, iPt, YZ_), &
                       self%aoG2(:, iPt, XZ_))

      end do
    case default
      write (*,'("Invalid grad level=",I2," in xc_engine_t % COMPAOS")') nDer
      stop
    end select

  end subroutine

  subroutine pruneAOs(self, skip)
    class(xc_engine_t), target :: self

    logical :: skip
    integer :: i, numAOs_p

    numAOs_p = 0
    ! Save significant indices
    do i = 1, self%numAOs
      if (maxval(abs(self%aoV(i, :))) > self%ao_threshold) then
        numAOs_p = numAOs_p + 1
        self%indices_p(numAOs_p) = i
      end if
    end do

    ! Cycle if all AOs are pruned
    skip = numAOs_p == 0
    if (skip) return

    ! Check if the number of runed AOs is less
    ! than the prune cutoff (approximately 90%); if so, then
    ! grid pruning should be skipped.
    self%skip_p = real(numAOs_p) / real(self%numAOs) > self%ao_sparsity_ratio

    if (self%skip_p) then
      ! Set the full number of AOs since we skip pruning AOs
      self%numAOs_p = self%numAOs

      self%wfAlpha_p => self%wfAlpha
      if (self%hasbeta) &
        self%wfBeta_p => self%wfBeta

    else
      ! Set the number of pruned AOs
      self%numAOs_p = numAOs_p

      call self%ResetPrunedPointers

    end if

  end subroutine
!###############################################################################

!> @brief Adjust XC memory storage for a given
!>  number of pruned grid points
!> @author Konstantin Komarov
  subroutine resetPrunedPointers(self)
    class(xc_engine_t), target :: self

    real(kind=fp), pointer, dimension(:,:,:) :: reorderable_data

    associate(  numAOs    => self%numAOs &
              , numAOs_p  => self%numAOs_p &
              , numPts    => self%numPts &
              , hasBeta   => self%hasBeta &
              , isWFVecs  => self%isWFVecs &
              , numAOVecs => self%numAOVecs &
              , numTmpVec => self%numTmpVec &
              , indices   => self%indices_p &
      )
      ! Setup the pointer for reorderable data
      reorderable_data => self%aoMem(:, :, :)

      ! Update pointers with pruned AOs
      self%aoMem(1:numAOs_p, 1:numPts, 1:numAOVecs) => self%aoMem_(1:)

      if (isWFVecs) then

        ! Set pointer for pruned
        self%wfAlpha_p(1:numAOs_p, 1:numAOs) => self%tmpWfAlpha(1:numAOs_p*numAOs)
        ! Compress array
        self%wfAlpha_p(:numAOs_p,:) = self%wfAlpha(indices(:numAOs_p),:)

        if (hasBeta) then
          ! Set pointer for pruned
          self%wfBeta_p(1:numAOs_p, 1:numAOs) => self%tmpWfBeta(1:numAOs_p*numAOs)
          ! Compress array
          self%wfBeta_p(:numAOs_p, :) = self%wfBeta(indices(:numAOs_p),:)
        end if

      else

        ! Set pointer for pruned
        self%moMemA(1:numAOs_p, 1:numPts, 1:numAOVecs) => self%moMemA_(1:)
        self%wfAlpha_p(1:numAOs_p, 1:numAOs_p) => self%tmpWfAlpha(1:numAOs_p*numAOs_p)
        ! Compress array
        self%wfAlpha_p(:numAOs_p, :numAOs_p) = self%wfAlpha(indices(:numAOs_p), indices(:numAOs_p))

        if (hasBeta) then
          ! Set pointer for pruned
          self%moMemB(1:numAOs_p, 1:numPts, 1:numAOVecs) => self%moMemB_(1:)
          self%wfBeta_p(1:numAOs_p, 1:numAOs_p) => self%tmpWfBeta(1:numAOs_p*numAOs_p)
          ! Compress array
          self%wfBeta_p(:numAOs_p, :numAOs_p) = self%wfBeta(indices(:numAOs_p), indices(:numAOs_p))
        end if

      end if

      select case (self%nAODer)
      case (0)
        ! Compress array
        self%aoMem(1:numAOs_p, :, 1:1) = reorderable_data(indices(1:numAOs_p), :, 1:1)

        ! Update pointers for pruned
        self%aoV => self%aoMem(:, :, 1)
        self%moVA => self%moMemA(:, :, 1)
        if (hasBeta) &
          self%moVB => self%moMemB(:, :, 1)

      case (1)
        ! Compress array
        self%aoMem(1:numAOs_p, :, 1:4) = reorderable_data(indices(1:numAOs_p), :, 1:4)

        ! Update pointers for pruned
        self%aoV => self%aoMem(:, :, 1)
        self%aoG1 => self%aoMem(:, :, 2:4)
        self%moVA => self%moMemA(:, :, 1)
        self%moG1A => self%moMemA(:, :, 2:4)
        if (hasBeta) then
          self%moVB => self%moMemB(:, :, 1)
          self%moG1B => self%moMemB(:, :, 2:4)
        end if

      case (2)
        ! Compress array
        self%aoMem(1:numAOs_p, :, 1:10) = reorderable_data(indices(1:numAOs_p), :, 1:10)

        ! Update pointers for pruned
        self%aoV => self%aoMem(:, :, 1)
        self%aoG1 => self%aoMem(:, :, 2:4)
        self%aoG2 => self%aoMem(:, :, 5:10)

        self%moVA => self%moMemA(:, :, 1)
        self%moG1A => self%moMemA(:, :, 2:4)
        self%moG2A => self%moMemA(:, :, 5:10)
        if (hasBeta) then
          self%moVB => self%moMemB(:, :, 1)
          self%moG1B => self%moMemB(:, :, 2:4)
          self%moG2B => self%moMemB(:, :, 5:10)
        end if
      end select

    end associate

  end subroutine

!###############################################################################

!> @brief Transform AOs to "MOs"
!> @details Multiply AO vector to the MO coefficient matrix or density matrix.
!>  True MOs are only obtained in the former case.
!> @author Vladimir Mironov
 subroutine compMOs(self)
    class(xc_engine_t) :: self
    integer :: nVecs, nPts

    nVecs = min(self%numAOVecs, 4) ! Don't transform second derivatives
    nPts  = ubound(self%aoMem,2)

    associate( nAlpha   => self%numOccAlpha &
             , nBeta    => self%numOccBeta  &
             , isWFVecs => self%isWFVecs    &
             , numAOs   => self%numAOs      &
             , numAOs_p => self%numAOs_p &
             , hasBeta  => self%hasBeta     &
      )

      if (isWFVecs) then
        call mo_tran_gemm_(nAlpha, numAOs_p, nVecs, nPts, numAOs, self%wfAlpha_p, self%aoMem, self%moMemA)
      else
        call mo_tran_symm_(numAOs_p, nVecs, nPts, self%wfAlpha_p, self%aoMem, self%moMemA)
      end if

      if (.not. hasBeta) return

      if (isWFVecs) then
        call mo_tran_gemm_(nBeta, numAOs_p, nVecs, nPts, numAOs, self%wfBeta_p, self%aoMem, self%moMemB)
      else
        call mo_tran_symm_(numAOs_p, nVecs, nPts, self%wfBeta_p, self%aoMem, self%moMemB)
      end if

    end associate
 end subroutine

!###############################################################################

 subroutine compRhoAll(self, skip)
    class(xc_engine_t) :: self
    logical, intent(out) :: skip
    real(kind=fp) :: rhoab

    skip = .false.

    ! electronic density
    call self%compRho(self%XCLib%rho)

    rhoab = dot_product(self%wts, sum(self%XCLib%rho, dim=1))
    if (rhoab < 1.0d-12) then
      skip = .true.
      return
    end if

    self%N_elec = self%N_elec + rhoab

    if (self%funTyp /= OQP_FUNTYP_LDA) then
      ! electronic density 1st derivative
      CALL self%compDRho(self%XCLib%drho, self%XCLib%sig)

      ! The total electron density gradient
      if (self%dbgLevel > 1) then
        self%G_total(1) = self%G_total(1) &
                        + dot_product(self%wts, self%XCLib%sig(1,:))
        self%G_total(2) = self%G_total(2) &
                        + dot_product(self%wts, self%XCLib%sig(2,:))
        self%G_total(3) = self%G_total(3) &
                        + dot_product(self%wts, self%XCLib%sig(3,:))
      end if
    end if

    if (self%funTyp == OQP_FUNTYP_MGGA) then
      ! electronic density 2nd derivative
      call self%compTau(self%XCLib%tau)

      if (self%dbgLevel > 1) then
          self%E_kin = self%E_kin &
                 + dot_product(self%wts, sum(self%XCLib%tau, dim=1))
      end if
    end if

 end subroutine


 subroutine compRMOs(xce, da, mo)
    class(xc_engine_t) :: xce
    real(kind=fp) :: da(:,:,:)
    real(kind=fp) :: mo(:,:,:)
    integer :: nPts, nMtx, i

    nMtx = ubound(da, 3)
    nPts = xce%numPts

    do i = 1, nMtx
      call mo_tran_symm_( &
              xce%numAOs_p, 1, nPts, da(:,:,i), xce%aoV, mo(:,:,i))
    end do

 end subroutine

 subroutine compRMOGs(xce, da, moG1)
    class(xc_engine_t) :: xce
    real(kind=fp) :: da(:,:,:)
    real(kind=fp) :: moG1(:,:,:,:)
    integer :: nPts, nMtx, i, j

    nMtx = ubound(da, 3)
    nPts = xce%numPts

    do i = 1, nMtx
      do j = 1, 3
        call mo_tran_symm_(&
                xce%numAOs_p, 1, nPts, da(:,:,i), &
                xce%aoG1(:,:,j), moG1(:,:,j,i))
      end do
    end do

 end subroutine

!> @brief Compute electronic density in a grid point, density-driven calculation
!> @param[in]  xce      XC engine, parameters
!> @param[in]  mo       "Molecular orbitals"
!> @param[out] rho      electronic density
!> @author Vladimir Mironov
 subroutine compRRho_ab(xce, mo, rho)

    class(xc_engine_t) :: xce
    real(kind=fp), intent(out) :: rho(:,:,:)
    real(kind=fp), intent(in) :: mo(:,:,:,:)
    integer :: i, j, k, nMtx, nSpin

    nMtx = ubound(mo,3)
    nSpin = ubound(mo,4)

    do k = 1, nSpin
      do j = 1, nMtx
        do i = 1, xce%numPts
          rho(k,i,j) = dot_product(xce%aoV(:,i), mo(:,i,j,k))
        end do
      end do
    end do

 end subroutine

!> @brief Compute electronic density gradient in a grid point, density-driven calculation
!> @param[in]  xce      XC engine, parameters
!> @param[in]  mo       "Molecular orbitals"
!> @param[out] drho     density directional derivative (along X, Y, and Z axes)
!> @param[out] drrho    dRho/d[x,y,z] vector and its dot product with `drho`
!> @author Vladimir Mironov
 subroutine compRDRho_ab(xce, mo, drrho)

    class(xc_engine_t) :: xce
    real(kind=fp), intent(in) :: mo(:,:,:,:)
    real(kind=fp), intent(out) :: drrho(:,:,:,:)

    integer :: i, j, k, nMtx, nSpin

    nMtx = ubound(mo,3)
    nSpin = ubound(mo,4)

    do k = 1, nSpin
      do j = 1, nMtx
        do i = 1, xce%numPts
          drrho(1,k,i,j) = 2*dot_product(xce%aoG1(:,i,1), mo(:,i,j,k))
          drrho(2,k,i,j) = 2*dot_product(xce%aoG1(:,i,2), mo(:,i,j,k))
          drrho(3,k,i,j) = 2*dot_product(xce%aoG1(:,i,3), mo(:,i,j,k))
        end do
      end do
    end do

 end subroutine

!> @brief Compute tau: (MO)' times (AO)'
!> @param[in]  xce      XC engine, parameters
!> @param[in]  mog1    MO directional derivatives
!> @param[out] rtau     kinetic energy density
!> @author Vladimir Mironov
 subroutine compRTau_ab(xce, moG1, rtau)

    class(xc_engine_t) :: xce
    real(kind=fp), intent(out) :: rtau(:,:,:)
    real(kind=fp), intent(in) :: moG1(:,:,:,:,:)
    integer :: i, j, k, nSpin, nMtx

    nSpin = ubound(moG1,5)
    nMtx = ubound(moG1, 4)

    do k = 1, nSpin
      do j = 1, nMtx
        do i = 1, xce%numPts
          rtau(k,i,j) = 0.5*sum(xce%aoG1(:,i,1:3)*moG1(:,i,1:3,j,k))
        end do
      end do
    end do

 end subroutine


!> @brief Compute electronic density in a grid point, density-driven calculation
!> @param[in]  xce      XC engine, parameters
!> @param[in]  mo       "Molecular orbitals"
!> @param[out] rho      electronic density
!> @author Vladimir Mironov
 subroutine compRRho_a(xce, mo, rho)

    class(xc_engine_t) :: xce
    real(kind=fp), intent(out) :: rho(:,:)
    real(kind=fp), intent(in) :: mo(:,:,:)
    integer :: j, i, nMtx

    nMtx = ubound(mo,3)

    do j = 1, nMtx
      do i = 1, xce%numPts
        rho(i,j) = dot_product(xce%aoV(:,i), mo(:,i,j))
      end do
    end do

 end subroutine

!> @brief Compute electronic density gradient in a grid point, density-driven calculation
!> @param[in]  xce      XC engine, parameters
!> @param[in]  mo       "Molecular orbitals"
!> @param[out] drho     density directional derivative (along X, Y, and Z axes)
!> @param[out] drrho    dRho/d[x,y,z] vector and its dot product with `drho`
!> @author Vladimir Mironov
 subroutine compRDRho_a(xce, mo, drrho)

    class(xc_engine_t) :: xce
    real(kind=fp), intent(in) :: mo(:,:,:)
    real(kind=fp), intent(out) :: drrho(:,:,:)

    integer :: i, j, nMtx

    nMtx = ubound(mo,3)

    do j = 1, nMtx
      do i = 1, xce%numPts
        drrho(1,i,j) = 2*dot_product(xce%aoG1(:,i,1), mo(:,i,j))
        drrho(2,i,j) = 2*dot_product(xce%aoG1(:,i,2), mo(:,i,j))
        drrho(3,i,j) = 2*dot_product(xce%aoG1(:,i,3), mo(:,i,j))
      end do
    end do

 end subroutine

!> @brief Compute tau: (MO)' times (AO)'
!> @param[in]  xce     XC engine, parameters
!> @param[in]  mog1    MO directional derivatives
!> @param[out] tau     kinetic energy density
!> @author Vladimir Mironov
 subroutine compRTau_a(xce, moG1, tau)

    class(xc_engine_t) :: xce
    real(kind=fp), intent(out) :: tau(:,:)
    real(kind=fp), intent(in) :: moG1(:,:,:,:)
    integer :: i, j

    do j = 1, ubound(moG1,4)
      do i = 1, xce%numPts
        tau(i,j) = 0.5*sum(xce%aoG1(:,i,1:3)*moG1(:,i,1:3,j))
      end do
    end do

 end subroutine
!###############################################################################

!> @brief Compute XC contribution to the energy
!> @param[inout]  bfGrad     array of gradient contributinos per basis function
!> @param[inout]  exec       XC energy integral
!> @param[inout]  ecorl      correlation energy integral
!> @param[inout]  totele     density integral == number of electrons
!> @param[inout]  totkin     kinetic energy integral
!> @param[inout]  togradxyz  density gradient integral
!> @author Vladimir Mironov
  subroutine compXC(self, functional, skip)
    class(xc_engine_t) :: self
    type(functional_t) :: functional
    logical :: skip

!   Compute MOs
    call self%compMOs

    call self%compRhoAll(skip)

    if (skip) return

    call self%XCLib%compute(functional, self%wts)

  end subroutine

!###############################################################################

!> @brief Compute electronic density in a grid point, AO-driven calculation
!> @param[out] rho     electronic density, alpha-spin
!> @author Vladimir Mironov
  subroutine compRhoAO(self, rho)

    class(xc_engine_t) :: self
    real(kind=fp), intent(out) :: rho(:,:)
    integer :: i


    if (self%hasBeta) then
      do i = 1, self%numPts
        rho(1,i) = dot_product(self%aoV(:,i), self%moVA(:,i))
        rho(2,i) = dot_product(self%aoV(:,i), self%moVB(:,i))
      end do
    else
      do i = 1, self%numPts
        rho(1,i) = 0.5_fp*dot_product(self%aoV(:,i), self%moVA(:,i))
        rho(2,i) = rho(1,i)
      end do
    end if

  end subroutine

!###############################################################################

!> @brief Compute electronic density in a grid point, MO-driven calculation
!> @param[out] rho     electronic density, alpha-spin
!> @author Vladimir Mironov
  subroutine compRhoMO(self, rho)

    class(xc_engine_t) :: self
    real(kind=fp), intent(out) :: rho(:,:)
    integer :: i, noa, nob

    noa = self%numOccAlpha
    nob = self%numOccBeta

    if (self%hasBeta) then
      do i = 1, self%numPts
        rho(1,i) = dot_product(self%moVA(:noa,i), self%moVA(:noa,i))
        rho(2,i) = dot_product(self%moVB(:nob,i), self%moVB(:nob,i))
      end do
    else
      do i = 1, self%numPts
        rho(1,i) = dot_product(self%moVA(:noa,i), self%moVA(:noa,i))
        rho(2,i) = rho(1,i)
      end do
    end if

  end subroutine

!###############################################################################

!> @brief Compute electronic density gradient in a grid point, AO-driven calculation
!> @param[out] drhoa   dRho, alpha-spin
!> @param[out] drhob   dRho, beta-spin
!> @author Vladimir Mironov
  subroutine compDRhoAO(self, drho, sigma)

    class(xc_engine_t) :: self
    real(kind=fp), intent(out) :: drho(:,:), sigma(:,:)
    integer :: i, j
    real(kind=fp) :: drhoa(3), drhob(3)

    do i = 1, self%numPts
      if (self%hasBeta) then
        do j = 1, 3
          drhoa(j) = 2*dot_product(self%aoG1(:,i,j), self%moVA(:,i))
          drhob(j) = 2*dot_product(self%aoG1(:,i,j), self%moVB(:,i))
        end do
      else
        do j = 1, 3
          drhoa(j) = dot_product(self%aoG1(:,i,j), self%moVA(:,i))
        end do
        drhob = drhoa
      end if
      drho(1:3,i) = drhoa
      drho(4:6,i) = drhob
      sigma(1,i) = dot_product(drhoa, drhoa)
      sigma(2,i) = dot_product(drhoa, drhob)
      sigma(3,i) = dot_product(drhob, drhob)
    end do

  end subroutine

!###############################################################################

!> @brief Compute electronic density gradient in a grid point, MO-driven calculation
!> @param[out] drhoa   dRho, alpha-spin
!> @param[out] drhob   dRho, beta-spin
!> @author Vladimir Mironov
  subroutine compDRhoMO(self, drho, sigma)

    class(xc_engine_t) :: self
    real(kind=fp), intent(out) :: drho(:,:), sigma(:,:)
    integer :: i, j, noa, nob
    real(kind=fp) :: drhoa(3), drhob(3)

    noa = self%numOccAlpha
    nob = self%numOccBeta

    do i = 1, self%numPts
      if (self%hasBeta) then
        do j = 1, 3
          drhoa(j) = 2*dot_product(self%moG1A(:noa,i,j), self%moVA(:noa,i))
          drhob(j) = 2*dot_product(self%moG1B(:nob,i,j), self%moVB(:nob,i))
        end do
      else
        do j = 1, 3
          drhoa(j) = 2*dot_product(self%moG1A(:noa,i,j), self%moVA(:noa,i))
        end do
        drhob = drhoa
      end if
      drho(1:3,i) = drhoa
      drho(4:6,i) = drhob
      sigma(1,i) = dot_product(drhoa, drhoa)
      sigma(2,i) = dot_product(drhoa, drhob)
      sigma(3,i) = dot_product(drhob, drhob)
    end do

  end subroutine

!###############################################################################

!> @brief Compute electronic density 2nd derivatives in a grid point, AO-driven calculation
!> @param[out] taua   d^2(Rho) alpha-spin
!> @param[out] taub   d^2(Rho), beta-spin
!> @author Vladimir Mironov
  subroutine compTauAO(self, tau)

    class(xc_engine_t) :: self
    real(kind=fp), intent(out) :: tau(:,:)
    real(kind=fp) :: taua(3), taub(3)
    integer :: i, j

    do i = 1, self%numPts
      if (self%hasBeta) then
        do j = 1, 3
          taua(j) = dot_product(self%aoG1(:,i,j), self%moG1A(:,i,j))
          taub(j) = dot_product(self%aoG1(:,i,j), self%moG1B(:,i,j))
        end do
      else
        do j = 1, 3
          taua(j) = 0.5_fp*dot_product(self%aoG1(:,i,j), self%moG1A(:,i,j))
        end do
        taub = taua
      end if

      tau(1,i) = 0.5*sum(taua)
      tau(2,i) = 0.5*sum(taub)

    end do

  end subroutine

!###############################################################################

!> @brief Compute electronic density 2nd derivatives in a grid point, AO-driven calculation
!> @param[out] taua   d^2(Rho) alpha-spin
!> @param[out] taub   d^2(Rho), beta-spin
!> @author Vladimir Mironov
  subroutine compTauMO(self, tau)

    class(xc_engine_t) :: self
    real(kind=fp), intent(out) :: tau(:,:)
    real(kind=fp) :: taua(3), taub(3)
    integer :: i, j, noa, nob

    noa = self%numOccAlpha
    nob = self%numOccBeta

    do i = 1, self%numPts
      if (self%hasBeta) then
        do j = 1, 3
          taua(j) = dot_product(self%moG1A(:noa,i,j), self%moG1A(:noa,i,j))
          taub(j) = dot_product(self%moG1B(:nob,i,j), self%moG1B(:nob,i,j))
        end do
      else
        do j = 1, 3
          taua(j) = dot_product(self%moG1A(:noa,i,j), self%moG1A(:noa,i,j))
        end do
        taub = taua
      end if

      tau(1,i) = 0.5*sum(taua)
      tau(2,i) = 0.5*sum(taub)

    end do

  end subroutine

!###############################################################################

!> @brief Compute XC contributions to the gradient from a grid point, LDA part
!> @param[in]    iPt      index of a grid point
!> @param[inout] bfGrad   array of gradient contributions per basis function
!> @param[in]    fgrad    XC gradient
!> @param[in]    moV      MO-like orbital values
!> @author Vladimir Mironov
  subroutine compAtGradRho(bfGrad, fgrad, moV, aoG1, nPts)
    real(kind=fp), contiguous, intent(in) :: moV(:,:), aoG1(:,:,:)
    real(kind=fp), intent(inout) :: bfGrad(:,:)
    real(kind=fp), intent(in) :: fGrad(:)
    integer, intent(in) :: nPts

    integer :: i, j
    real(kind=fp) :: f

    do i = 1, nPts
        f = fGrad(i)*2.0_fp
        do j = 1, 3
          bfGrad(:,j) = bfGrad(:,j) + aoG1(:,i,j)*moV(:,i)*f
        end do
    end do

  end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute XC contributions to the gradient from a grid point, GGA part
!> @param[inout] bfGrad   array of gradient contributions per basis function
!> @param[in]    fgrad    XC gradient
!> @param[in]    moV      MO-like orbital values
!> @param[in]    moG1     MO-like orbital gradients
!> @author Vladimir Mironov
  subroutine compAtGradDRho(bfGrad, fGrad, moV, moG1, aoG1, aoG2, nPts)
    real(kind=fp), intent(in) :: fGrad(:,:)
    real(kind=fp), intent(inout) :: bfGrad(:,:)
    real(kind=fp), contiguous, intent(in) :: moV(:,:), moG1(:,:,:)
    real(kind=fp), contiguous, intent(in) :: aoG1(:,:,:), aoG2(:,:,:)
    integer, intent(in) :: npts

    integer :: i, jt, j1, j2
    real(kind=fp) :: f(3)

    do i = 1, nPts
      f = fGrad(i,:)*2.0_fp
      do j1 = 1, 3 ! X__ Y__ Z__
          do j2 = 1, 3 ! X__ Y__ Z__
              jt = SQ_TO_TR(j1,j2)
              bfGrad(:,j1) = bfGrad(:,j1) + f(j2)* &
                  (  aoG2(:,i,jt)  *  moV(:,i) &
                   + aoG1(:,i,j1)  * moG1(:,i,j2) &
                  )
          end do
      end do
    end do

 end subroutine


!> @brief Compute XC contributions to the gradient from a grid point, mGGA part
!> @param[in]    iPt      index of a grid point
!> @param[inout] bfGrad   array of gradient contributinos per basis function
!> @param[in]    dedta    XC energy, mGGA contribution, alpha-spin
!> @param[in]    dedtb    XC energy, mGGA contribution, beta-spin
!> @author Vladimir Mironov
 subroutine compAtGradTau(bfGrad, fgrad, moG1, aoG2, npts)
    real(kind=fp), intent(in) :: fgrad(:)
    real(kind=fp), intent(inout) :: bfGrad(:,:)
    real(kind=fp), contiguous, intent(in) :: moG1(:,:,:)
    real(kind=fp), contiguous, intent(in) :: aoG2(:,:,:)
    integer, intent(in) :: npts

    integer :: i, j1, j2, jt
    real(kind=fp) :: f

    do i = 1, npts
        f = fgrad(i)
        do j1 = 1, 3 ! x__ y__ z__
            do j2 = 1, 3 ! x__ y__ z__
                jt = sq_to_tr(j1,j2)
                bfgrad(:,j1) = bfgrad(:,j1) + f*aoG2(:,i,jt)*moG1(:,i,j2)
            end do
        end do
    end do

 end subroutine


!> @brief Get first derivative of the XC functional
!> @param[in] xce    XC engine
!> @param[in] beta   Whether to return spin-polarized quantities
!> @param[out] d_r   dE_xc / d_rho (alpha, beta)
!> @param[out] d_s   dE_xc / d_sigma (alpha-alpha, beta-beta, alpha-beta)
!> @param[out] d_t   dE_xc / d_tau (alpha, beta)
  subroutine xc_der1(xce, beta, ipt, &
                     d_r, d_s, d_t)
    class(xc_engine_t) :: xce
    logical :: beta
    real(kind=fp), intent(out) :: d_r(2), d_s(3), d_t(2)
    integer, intent(in) :: ipt

    associate (xc => xce%XCLib, ids => xce%XCLib%ids, i => ipt)
      if (beta) then
        d_r(1) = xc%d1dr(ids%ra, i)
        d_r(2) = xc%d1dr(ids%rb, i)

        d_s(1) = xc%d1ds(ids%ga, i)
        d_s(2) = xc%d1ds(ids%gb, i)
        d_s(3) = xc%d1ds(ids%gc, i)

        d_t(1) = xc%d1dt(ids%ta, i)
        d_t(2) = xc%d1dt(ids%tb, i)
      else
        d_r      = xc%d1dr(ids%ra, i)
        d_s(1:2) = xc%d1ds(ids%ga, i)
        d_s(3)   = xc%d1ds(ids%gc, i)
        d_t      = xc%d1dt(ids%ta, i)
      end if
    end associate
  end subroutine

!###############################################################################

!> @brief Get second derivative of the XC functional contracted with response densities
!> @param[in] xce    XC engine
!> @param[in] beta   Whether to return spin-polarized quantities
!> @param[in] d_r    \delta_rho (alpha, beta)
!> @param[in] d_s    \delta_sigma (alpha-alpha, beta-beta, alpha-beta)
!> @param[in] d_t    \delta_tau (alpha, beta)
!> @param[out] f_r   \sum_i d2E_xc / (d_rho   * d_zeta_i) (alpha, beta)
!> @param[out] f_s   \sum_i d2E_xc / (d_sigma * d_zeta_i) (alpha-alpha, beta-beta, alpha-beta)
!> @param[out] f_t   \sum_i d2E_xc / (d_tau   * d_zeta_i) (alpha, beta)
  subroutine xc_der2_contr(xce, beta, ipt, &
                  dr, ds, dt, &
                  f_r, f_s, f_t)
    class(xc_engine_t) :: xce
    real(kind=fp), intent(in) :: dr(2), ds(3), dt(2)
    real(kind=fp), intent(out) :: f_r(2), f_s(3), f_t(2)
    integer, intent(in) :: ipt
    logical, intent(in) :: beta

    if (beta) then
        call xc_ab_der2_contr(xce, ipt, &
                  dr, ds, dt, &
                  f_r, f_s, f_t)
    else
        call xc_a_der2_contr(xce, ipt, &
                  dr(1), ds(1), dt(1), &
                  f_r, f_s, f_t)
    end if

  end subroutine

!###############################################################################

!> @brief Get second derivative of the XC functional contracted with response densities,
!>  spin-polarized version
!> @param[in] xce    XC engine
!> @param[in] d_r    \delta_rho (alpha, beta)
!> @param[in] d_s    \delta_sigma (alpha-alpha, beta-beta, alpha-beta)
!> @param[in] d_t    \delta_tau (alpha, beta)
!> @param[out] f_r   \sum_i d2E_xc / (d_rho   * d_zeta_i) (alpha, beta)
!> @param[out] f_s   \sum_i d2E_xc / (d_sigma * d_zeta_i) (alpha-alpha, beta-beta, alpha-beta)
!> @param[out] f_t   \sum_i d2E_xc / (d_tau   * d_zeta_i) (alpha, beta)
  subroutine xc_ab_der2_contr(xce, ipt, &
                  dr, ds, dt, &
                  f_r, f_s, f_t)
    class(xc_engine_t) :: xce
    real(kind=fp), intent(in) :: dr(2), ds(3), dt(2)
    real(kind=fp), intent(out) :: f_r(2), f_s(3), f_t(2)
    integer, intent(in) :: ipt

    real(kind=fp) :: cr_r(2)
    real(kind=fp) :: cr_s(2), cs_r(3), cs_s(3)
    real(kind=fp) :: cr_t(2), cs_t(3), ct_r(2), ct_s(2), ct_t(2)

    associate (xc => xce%XCLib, ids => xce%XCLib%ids, i => ipt)

      f_r = 0
      f_s = 0
      f_t = 0

      cr_r(1) = xc%d2r2(ids%rara, i) * dr(1) &
              + xc%d2r2(ids%rarb, i) * dr(2)

      cr_r(2) = xc%d2r2(ids%rarb, i) * dr(1) &
              + xc%d2r2(ids%rbrb, i) * dr(2)

      f_r = f_r + cr_r

      if (xce%funTyp /= OQP_FUNTYP_LDA) then

        cr_s(1) =  xc%d2rs(ids%raga, i) * ds(1) &
                +  xc%d2rs(ids%ragb, i) * ds(2) &
                +  xc%d2rs(ids%ragc, i) * ds(3)

        cr_s(2) =  xc%d2rs(ids%rbga, i) * ds(1) &
                +  xc%d2rs(ids%rbgb, i) * ds(2) &
                +  xc%d2rs(ids%rbgc, i) * ds(3)

        f_r = f_r + cr_s


        cs_r(1) = xc%d2rs(ids%raga, i) * dr(1) &
                + xc%d2rs(ids%rbga, i) * dr(2)

        cs_r(2) = xc%d2rs(ids%ragb, i) * dr(1) &
                + xc%d2rs(ids%rbgb, i) * dr(2)

        cs_r(3) = xc%d2rs(ids%ragc, i) * dr(1) &
                + xc%d2rs(ids%rbgc, i) * dr(2)

        f_s = f_s + cs_r

        cs_s(1) = xc%d2s2(ids%gaga, i) * ds(1) &
                + xc%d2s2(ids%gagb, i) * ds(2) &
                + xc%d2s2(ids%gagc, i) * ds(3)

        cs_s(2) = xc%d2s2(ids%gagb, i) * ds(1) &
                + xc%d2s2(ids%gbgb, i) * ds(2) &
                + xc%d2s2(ids%gbgc, i) * ds(3)

        cs_s(3) = xc%d2s2(ids%gagc, i) * ds(1) &
                + xc%d2s2(ids%gbgc, i) * ds(2) &
                + xc%d2s2(ids%gcgc, i) * ds(3)

        f_s = f_s + cs_s

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          cr_t(1) =  xc%d2rt(ids%rata, i) * dt(1) &
                  +  xc%d2rt(ids%ratb, i) * dt(2)

          cr_t(2) =  xc%d2rt(ids%rbta, i) * dt(1) &
                  +  xc%d2rt(ids%rbtb, i) * dt(2)

          f_r = f_r + cr_t

          cs_t(1) = xc%d2st(ids%gata, i) * dt(1) &
                  + xc%d2st(ids%gatb, i) * dt(2)

          cs_t(2) = xc%d2st(ids%gbta, i) * dt(1) &
                  + xc%d2st(ids%gbtb, i) * dt(2)

          cs_t(3) = xc%d2st(ids%gcta, i) * dt(1) &
                  + xc%d2st(ids%gctb, i) * dt(2)

          f_s = f_s + cs_t

          ct_r(1) = xc%d2rt(ids%rata, i) * dr(1) &
                  + xc%d2rt(ids%rbta, i) * dr(2)

          ct_r(2) = xc%d2rt(ids%ratb, i) * dr(1) &
                  + xc%d2rt(ids%rbtb, i) * dr(2)

          f_t = f_t + ct_r

          ct_s(1) = xc%d2st(ids%gata, i) * ds(1) &
                  + xc%d2st(ids%gbta, i) * ds(2) &
                  + xc%d2st(ids%gcta, i) * ds(3)

          ct_s(2) = xc%d2st(ids%gatb, i) * ds(1) &
                  + xc%d2st(ids%gbtb, i) * ds(2) &
                  + xc%d2st(ids%gctb, i) * ds(3)


          f_t = f_t + ct_s

          ct_t(1) = xc%d2t2(ids%tata, i) * dt(1) &
                  + xc%d2t2(ids%tatb, i) * dt(2)

          ct_t(2) = xc%d2t2(ids%tatb, i) * dt(1) &
                  + xc%d2t2(ids%tbtb, i) * dt(2)

          f_t = f_t + ct_t

        end if
      end if

    end associate

  end subroutine

!###############################################################################

!> @brief Get second derivative of the XC functional contracted with response densities,
!>   not spin-polarized version
!> @param[in] xce    XC engine
!> @param[in] d_r    \delta_rho (alpha, beta)
!> @param[in] d_s    \delta_sigma (alpha-alpha, beta-beta, alpha-beta)
!> @param[in] d_t    \delta_tau (alpha, beta)
!> @param[out] f_r   \sum_i d2E_xc / (d_rho   * d_zeta_i) (alpha, beta)
!> @param[out] f_s   \sum_i d2E_xc / (d_sigma * d_zeta_i) (alpha-alpha, beta-beta, alpha-beta)
!> @param[out] f_t   \sum_i d2E_xc / (d_tau   * d_zeta_i) (alpha, beta)
  subroutine xc_a_der2_contr(xce, ipt, &
                  dr, ds, dt, &
                  f_r, f_s, f_t)
    class(xc_engine_t) :: xce
    real(kind=fp), intent(in) :: dr, ds, dt
    real(kind=fp), intent(out) :: f_r(2), f_s(3), f_t(2)
    integer, intent(in) :: ipt

    real(kind=fp) :: cr_r(2)
    real(kind=fp) :: cr_s(2), cs_r(3), cs_s(3)
    real(kind=fp) :: cr_t(2), cs_t(3), ct_r(2), ct_s(2), ct_t(2)

    associate (xc => xce%XCLib, ids => xce%XCLib%ids, i => ipt)

      f_r = 0
      f_s = 0
      f_t = 0

      cr_r(1:2) = xc%d2r2(ids%rara, i) * dr &
                + xc%d2r2(ids%rarb, i) * dr

      f_r = f_r + cr_r

      if (xce%funTyp /= OQP_FUNTYP_LDA) then

        cr_s    =  xc%d2rs(ids%raga, i) * ds &
                +  xc%d2rs(ids%ragb, i) * ds &
                +  xc%d2rs(ids%ragc, i) * ds

        f_r = f_r + cr_s


        cs_r(1) = xc%d2rs(ids%raga, i) * dr &
                + xc%d2rs(ids%rbga, i) * dr

        cs_r(2) = cs_r(1)


        cs_r(3) = xc%d2rs(ids%ragc, i) * dr &
                + xc%d2rs(ids%rbgc, i) * dr

        f_s = f_s + cs_r

        cs_s(1) = xc%d2s2(ids%gaga, i) * ds &
                + xc%d2s2(ids%gagb, i) * ds &
                + xc%d2s2(ids%gagc, i) * ds

        cs_s(2) = cs_s(1)

        cs_s(3) = xc%d2s2(ids%gagc, i) * ds &
                + xc%d2s2(ids%gbgc, i) * ds &
                + xc%d2s2(ids%gcgc, i) * ds

        f_s = f_s + cs_s

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          cr_t    =  xc%d2rt(ids%rata, i) * dt &
                  +  xc%d2rt(ids%ratb, i) * dt

          f_r = f_r + cr_t

          cs_t(1) = xc%d2st(ids%gata, i) * dt &
                  + xc%d2st(ids%gatb, i) * dt

          cs_t(2) = cs_t(1)

          cs_t(3) = xc%d2st(ids%gcta, i) * dt &
                  + xc%d2st(ids%gctb, i) * dt

          f_s = f_s + cs_t

          ct_r = xc%d2rt(ids%rata, i) * dr &
               + xc%d2rt(ids%rbta, i) * dr

          f_t = f_t + ct_r

          ct_s = xc%d2st(ids%gata, i) * ds &
               + xc%d2st(ids%gbta, i) * ds &
               + xc%d2st(ids%gcta, i) * ds

          f_t = f_t + ct_s

          ct_t = xc%d2t2(ids%tata, i) * dt &
               + xc%d2t2(ids%tatb, i) * dt

          f_t = f_t + ct_t

        end if
      end if

    end associate

  end subroutine

!###############################################################################

!> @brief Get third derivative of the XC functional contracted with response densities,
!>  spin-polarized version
!> @param[in] xce    XC engine
!> @param[in] d_r    \delta_rho (alpha, beta)
!> @param[in] d_s    \delta_sigma (alpha-alpha, beta-beta, alpha-beta)
!> @param[in] d_t    \delta_tau (alpha, beta)
!> @param[in] ss     (\nabla\rho(T)*\nabla\rho(T)) (alpha-alpha, beta-beta, alpha-beta)
!> @param[out] g_r   \sum_i,j d2E_xc / (d_rho   * d_zeta_i*d_zeta_j) (alpha, beta)
!> @param[out] g_s   \sum_i,j d2E_xc / (d_sigma * d_zeta_i*d_zeta_j) (alpha-alpha, beta-beta, alpha-beta)
!> @param[out] g_t   \sum_i,j d2E_xc / (d_tau   * d_zeta_i*d_zeta_j) (alpha, beta)
  subroutine xc_der3_contr(xce, ipt, &
                  dr, ds, dt, &
                  ss, &
                  f_s, &
                  g_r, g_s, g_t)
    class(xc_engine_t) :: xce
    real(kind=fp), intent(in) :: dr(2), ds(3), dt(2)
    real(kind=fp), intent(in) :: ss(3)
    real(kind=fp), intent(out) :: f_s(3)
    real(kind=fp), intent(out) :: g_r(2), g_s(3), g_t(2)
    integer, intent(in) :: ipt

    real(kind=fp) :: cr(2), cs(3), ct(2)

    associate (xc => xce%XCLib, ids => xce%XCLib%ids, i => ipt)

      f_s = 0
      g_r = 0
      g_s = 0
      g_t = 0

      cr(1) = xc%d3r3(ids%rarara, i) * dr(1)*dr(1) &
            + xc%d3r3(ids%rararb, i) * dr(1)*dr(2) &
            + xc%d3r3(ids%rararb, i) * dr(2)*dr(1) &
            + xc%d3r3(ids%rarbrb, i) * dr(2)*dr(2)

      cr(2) = xc%d3r3(ids%rararb, i) * dr(1)*dr(1) &
            + xc%d3r3(ids%rarbrb, i) * dr(1)*dr(2) &
            + xc%d3r3(ids%rarbrb, i) * dr(2)*dr(1) &
            + xc%d3r3(ids%rbrbrb, i) * dr(2)*dr(2)

      g_r = g_r + cr

      if (xce%funTyp /= OQP_FUNTYP_LDA) then

        cs(1) = xc%d2rs(ids%raga, i) * dr(1) &
              + xc%d2rs(ids%rbga, i) * dr(2)

        cs(2) = xc%d2rs(ids%ragb, i) * dr(1) &
              + xc%d2rs(ids%rbgb, i) * dr(2)

        cs(3) = xc%d2rs(ids%ragc, i) * dr(1) &
              + xc%d2rs(ids%rbgc, i) * dr(2)

        f_s = f_s + cs

        cs(1) = xc%d2s2(ids%gaga, i) * ds(1) &
              + xc%d2s2(ids%gagb, i) * ds(2) &
              + xc%d2s2(ids%gagc, i) * ds(3)

        cs(2) = xc%d2s2(ids%gagb, i) * ds(1) &
              + xc%d2s2(ids%gbgb, i) * ds(2) &
              + xc%d2s2(ids%gbgc, i) * ds(3)

        cs(3) = xc%d2s2(ids%gagc, i) * ds(1) &
              + xc%d2s2(ids%gbgc, i) * ds(2) &
              + xc%d2s2(ids%gcgc, i) * ds(3)

        f_s = f_s + cs

        cr(1) = xc%d2rs(ids%raga, i) * ss(1) &
              + xc%d2rs(ids%ragb, i) * ss(2) &
              + xc%d2rs(ids%ragc, i) * ss(3)

        cr(2) = xc%d2rs(ids%rbga, i) * ss(1) &
              + xc%d2rs(ids%rbgb, i) * ss(2) &
              + xc%d2rs(ids%rbgc, i) * ss(3)

        g_r = g_r + cr

        cs(1) = xc%d2s2(ids%gaga, i) * ss(1) &
              + xc%d2s2(ids%gagb, i) * ss(2) &
              + xc%d2s2(ids%gagc, i) * ss(3)

        cs(2) = xc%d2s2(ids%gagb, i) * ss(1) &
              + xc%d2s2(ids%gbgb, i) * ss(2) &
              + xc%d2s2(ids%gbgc, i) * ss(3)

        cs(3) = xc%d2s2(ids%gagc, i) * ss(1) &
              + xc%d2s2(ids%gbgc, i) * ss(2) &
              + xc%d2s2(ids%gcgc, i) * ss(3)

        g_s = g_s + cs


        cr(1) =  xc%d3r2s(ids%raraga, i) * dr(1) * ds(1) &
              +  xc%d3r2s(ids%rarbga, i) * dr(2) * ds(1) &

              +  xc%d3r2s(ids%raragb, i) * dr(1) * ds(2) &
              +  xc%d3r2s(ids%rarbgb, i) * dr(2) * ds(2) &

              +  xc%d3r2s(ids%raragc, i) * dr(1) * ds(3) &
              +  xc%d3r2s(ids%rarbgc, i) * dr(2) * ds(3)

        cr(2) =  xc%d3r2s(ids%rarbga, i) * dr(1) * ds(1) &
              +  xc%d3r2s(ids%rbrbga, i) * dr(2) * ds(1) &

              +  xc%d3r2s(ids%rarbgb, i) * dr(1) * ds(2) &
              +  xc%d3r2s(ids%rbrbgb, i) * dr(2) * ds(2) &

              +  xc%d3r2s(ids%rarbgc, i) * dr(1) * ds(3) &
              +  xc%d3r2s(ids%rbrbgc, i) * dr(2) * ds(3)

        g_r = g_r + 2*cr

        cr(1) =  xc%d3rs2(ids%ragaga, i) * ds(1) * ds(1) &
              +  xc%d3rs2(ids%ragagb, i) * ds(1) * ds(2) &
              +  xc%d3rs2(ids%ragagc, i) * ds(1) * ds(3) &

              +  xc%d3rs2(ids%ragagb, i) * ds(2) * ds(1) &
              +  xc%d3rs2(ids%ragbgb, i) * ds(2) * ds(2) &
              +  xc%d3rs2(ids%ragbgc, i) * ds(2) * ds(3) &

              +  xc%d3rs2(ids%ragagc, i) * ds(3) * ds(1) &
              +  xc%d3rs2(ids%ragbgc, i) * ds(3) * ds(2) &
              +  xc%d3rs2(ids%ragcgc, i) * ds(3) * ds(3)


        cr(2) =  xc%d3rs2(ids%rbgaga, i) * ds(1) * ds(1) &
              +  xc%d3rs2(ids%rbgagb, i) * ds(1) * ds(2) &
              +  xc%d3rs2(ids%rbgagc, i) * ds(1) * ds(3) &

              +  xc%d3rs2(ids%rbgagb, i) * ds(2) * ds(1) &
              +  xc%d3rs2(ids%rbgbgb, i) * ds(2) * ds(2) &
              +  xc%d3rs2(ids%rbgbgc, i) * ds(2) * ds(3) &

              +  xc%d3rs2(ids%rbgagc, i) * ds(3) * ds(1) &
              +  xc%d3rs2(ids%rbgbgc, i) * ds(3) * ds(2) &
              +  xc%d3rs2(ids%rbgcgc, i) * ds(3) * ds(3)

        g_r = g_r + cr


        cs(1) = xc%d3r2s(ids%raraga, i) * dr(1) * dr(1) &
              + xc%d3r2s(ids%rarbga, i) * dr(1) * dr(2) &
              + xc%d3r2s(ids%rarbga, i) * dr(2) * dr(1) &
              + xc%d3r2s(ids%rbrbga, i) * dr(2) * dr(2)

        cs(2) = xc%d3r2s(ids%raragb, i) * dr(1) * dr(1) &
              + xc%d3r2s(ids%rarbgb, i) * dr(1) * dr(2) &
              + xc%d3r2s(ids%rarbgb, i) * dr(2) * dr(1) &
              + xc%d3r2s(ids%rbrbgb, i) * dr(2) * dr(2)

        cs(3) = xc%d3r2s(ids%raragc, i) * dr(1) * dr(1) &
              + xc%d3r2s(ids%rarbgc, i) * dr(1) * dr(2) &
              + xc%d3r2s(ids%rarbgc, i) * dr(2) * dr(1) &
              + xc%d3r2s(ids%rbrbgc, i) * dr(2) * dr(2)

        g_s = g_s + cs

        cs(1) = xc%d3rs2(ids%ragaga, i) * dr(1) * ds(1) &
              + xc%d3rs2(ids%ragagb, i) * dr(1) * ds(2) &
              + xc%d3rs2(ids%ragagc, i) * dr(1) * ds(3) &

              + xc%d3rs2(ids%rbgaga, i) * dr(2) * ds(1) &
              + xc%d3rs2(ids%rbgagb, i) * dr(2) * ds(2) &
              + xc%d3rs2(ids%rbgagc, i) * dr(2) * ds(3)

        cs(2) = xc%d3rs2(ids%ragagb, i) * dr(1) * ds(1) &
              + xc%d3rs2(ids%ragbgb, i) * dr(1) * ds(2) &
              + xc%d3rs2(ids%ragbgc, i) * dr(1) * ds(3) &

              + xc%d3rs2(ids%rbgagb, i) * dr(2) * ds(1) &
              + xc%d3rs2(ids%rbgbgb, i) * dr(2) * ds(2) &
              + xc%d3rs2(ids%rbgbgc, i) * dr(2) * ds(3)

        cs(3) = xc%d3rs2(ids%ragagc, i) * dr(1) * ds(1) &
              + xc%d3rs2(ids%ragbgc, i) * dr(1) * ds(2) &
              + xc%d3rs2(ids%ragcgc, i) * dr(1) * ds(3) &

              + xc%d3rs2(ids%rbgagc, i) * dr(2) * ds(1) &
              + xc%d3rs2(ids%rbgbgc, i) * dr(2) * ds(2) &
              + xc%d3rs2(ids%rbgcgc, i) * dr(2) * ds(3)

        g_s = g_s + 2*cs

        cs(1) = xc%d3s3(ids%gagaga, i) * ds(1) * ds(1) &
              + xc%d3s3(ids%gagagb, i) * ds(1) * ds(2) &
              + xc%d3s3(ids%gagagc, i) * ds(1) * ds(3) &

              + xc%d3s3(ids%gagagb, i) * ds(2) * ds(1) &
              + xc%d3s3(ids%gagbgb, i) * ds(2) * ds(2) &
              + xc%d3s3(ids%gagbgc, i) * ds(2) * ds(3) &

              + xc%d3s3(ids%gagagc, i) * ds(3) * ds(1) &
              + xc%d3s3(ids%gagbgc, i) * ds(3) * ds(2) &
              + xc%d3s3(ids%gagcgc, i) * ds(3) * ds(3)

        cs(2) = xc%d3s3(ids%gagagb, i) * ds(1) * ds(1) &
              + xc%d3s3(ids%gagbgb, i) * ds(1) * ds(2) &
              + xc%d3s3(ids%gagbgc, i) * ds(1) * ds(3) &

              + xc%d3s3(ids%gagbgb, i) * ds(2) * ds(1) &
              + xc%d3s3(ids%gbgbgb, i) * ds(2) * ds(2) &
              + xc%d3s3(ids%gbgbgc, i) * ds(2) * ds(3) &

              + xc%d3s3(ids%gagbgc, i) * ds(3) * ds(1) &
              + xc%d3s3(ids%gbgbgc, i) * ds(3) * ds(2) &
              + xc%d3s3(ids%gbgcgc, i) * ds(3) * ds(3)

        cs(3) = xc%d3s3(ids%gagagc, i) * ds(1) * ds(1) &
              + xc%d3s3(ids%gagbgc, i) * ds(1) * ds(2) &
              + xc%d3s3(ids%gagcgc, i) * ds(1) * ds(3) &

              + xc%d3s3(ids%gagbgc, i) * ds(2) * ds(1) &
              + xc%d3s3(ids%gbgbgc, i) * ds(2) * ds(2) &
              + xc%d3s3(ids%gbgcgc, i) * ds(2) * ds(3) &

              + xc%d3s3(ids%gagcgc, i) * ds(3) * ds(1) &
              + xc%d3s3(ids%gbgcgc, i) * ds(3) * ds(2) &
              + xc%d3s3(ids%gcgcgc, i) * ds(3) * ds(3)

        g_s = g_s + cs

        if (xce%funTyp == OQP_FUNTYP_MGGA) then

          cs(1) = xc%d2st(ids%gata, i) * dt(1) &
                + xc%d2st(ids%gatb, i) * dt(2)

          cs(2) = xc%d2st(ids%gbta, i) * dt(1) &
                + xc%d2st(ids%gbtb, i) * dt(2)

          cs(3) = xc%d2st(ids%gcta, i) * dt(1) &
                + xc%d2st(ids%gctb, i) * dt(2)

          f_s = f_s + cs

          ct(1) = xc%d2st(ids%gata, i) * ss(1) &
                + xc%d2st(ids%gbta, i) * ss(2) &
                + xc%d2st(ids%gcta, i) * ss(3)

          ct(2) = xc%d2st(ids%gatb, i) * ss(1) &
                + xc%d2st(ids%gbtb, i) * ss(2) &
                + xc%d2st(ids%gctb, i) * ss(3)

          g_t = g_t + ct


          cr(1) = xc%d3r2t(ids%rarata, i) * dr(1)*dt(1) &
                + xc%d3r2t(ids%rarbta, i) * dr(2)*dt(1) &
                + xc%d3r2t(ids%raratb, i) * dr(1)*dt(2) &
                + xc%d3r2t(ids%rarbtb, i) * dr(2)*dt(2)

          cr(2) = xc%d3r2t(ids%rarbta, i) * dr(1)*dt(1) &
                + xc%d3r2t(ids%rbrbta, i) * dr(2)*dt(1) &
                + xc%d3r2t(ids%rarbtb, i) * dr(1)*dt(2) &
                + xc%d3r2t(ids%rbrbtb, i) * dr(2)*dt(2)

          g_r = g_r + 2*cr

          cr(1) = xc%d3rst(ids%ragata, i) * ds(1)*dt(1) &
                + xc%d3rst(ids%ragbta, i) * ds(2)*dt(1) &
                + xc%d3rst(ids%ragcta, i) * ds(3)*dt(1) &

                + xc%d3rst(ids%ragatb, i) * ds(1)*dt(2) &
                + xc%d3rst(ids%ragbtb, i) * ds(2)*dt(2) &
                + xc%d3rst(ids%ragctb, i) * ds(3)*dt(2)


          cr(2) = xc%d3rst(ids%rbgata, i) * ds(1)*dt(1) &
                + xc%d3rst(ids%rbgbta, i) * ds(2)*dt(1) &
                + xc%d3rst(ids%rbgcta, i) * ds(3)*dt(1) &

                + xc%d3rst(ids%rbgatb, i) * ds(1)*dt(2) &
                + xc%d3rst(ids%rbgbtb, i) * ds(2)*dt(2) &
                + xc%d3rst(ids%rbgctb, i) * ds(3)*dt(2)

          g_r = g_r + 2*cr

          cr(1) = xc%d3rt2(ids%ratata, i) * dt(1)*dt(1) &
                + xc%d3rt2(ids%ratatb, i) * dt(2)*dt(1) &
                + xc%d3rt2(ids%ratatb, i) * dt(1)*dt(2) &
                + xc%d3rt2(ids%ratbtb, i) * dt(2)*dt(2)

          cr(2) = xc%d3rt2(ids%ratatb, i) * dt(1)*dt(1) &
                + xc%d3rt2(ids%rbtatb, i) * dt(2)*dt(1) &
                + xc%d3rt2(ids%ratbtb, i) * dt(1)*dt(2) &
                + xc%d3rt2(ids%rbtbtb, i) * dt(2)*dt(2)

          g_r = g_r + cr

          cs(1) = xc%d3rst(ids%ragata, i) * dr(1) * dt(1) &
                + xc%d3rst(ids%rbgata, i) * dr(2) * dt(1) &
                + xc%d3rst(ids%ragatb, i) * dr(1) * dt(2) &
                + xc%d3rst(ids%rbgatb, i) * dr(2) * dt(2)

          cs(2) = xc%d3rst(ids%ragbta, i) * dr(1) * dt(1) &
                + xc%d3rst(ids%rbgbta, i) * dr(2) * dt(1) &
                + xc%d3rst(ids%ragbtb, i) * dr(1) * dt(2) &
                + xc%d3rst(ids%rbgbtb, i) * dr(2) * dt(2)

          cs(3) = xc%d3rst(ids%ragcta, i) * dr(1) * dt(1) &
                + xc%d3rst(ids%rbgcta, i) * dr(2) * dt(1) &
                + xc%d3rst(ids%ragctb, i) * dr(1) * dt(2) &
                + xc%d3rst(ids%rbgctb, i) * dr(2) * dt(2)

          g_s = g_s + 2*cs

          cs(1) = xc%d3s2t(ids%gagata, i) * ds(1) * dt(1) &
                + xc%d3s2t(ids%gagbta, i) * ds(2) * dt(1) &
                + xc%d3s2t(ids%gagcta, i) * ds(3) * dt(1) &

                + xc%d3s2t(ids%gagatb, i) * ds(1) * dt(2) &
                + xc%d3s2t(ids%gagbtb, i) * ds(2) * dt(2) &
                + xc%d3s2t(ids%gagctb, i) * ds(3) * dt(2)

          cs(2) = xc%d3s2t(ids%gagbta, i) * ds(1) * dt(1) &
                + xc%d3s2t(ids%gbgbta, i) * ds(2) * dt(1) &
                + xc%d3s2t(ids%gbgcta, i) * ds(3) * dt(1) &

                + xc%d3s2t(ids%gagbtb, i) * ds(1) * dt(2) &
                + xc%d3s2t(ids%gbgbtb, i) * ds(2) * dt(2) &
                + xc%d3s2t(ids%gbgctb, i) * ds(3) * dt(2)

          cs(3) = xc%d3s2t(ids%gagcta, i) * ds(1) * dt(1) &
                + xc%d3s2t(ids%gbgcta, i) * ds(2) * dt(1) &
                + xc%d3s2t(ids%gcgcta, i) * ds(3) * dt(1) &

                + xc%d3s2t(ids%gagctb, i) * ds(1) * dt(2) &
                + xc%d3s2t(ids%gbgctb, i) * ds(2) * dt(2) &
                + xc%d3s2t(ids%gcgctb, i) * ds(3) * dt(2)

          g_s = g_s + 2*cs

          cs(1) = xc%d3st2(ids%gatata, i) * dt(1) * dt(1) &
                + xc%d3st2(ids%gatatb, i) * dt(2) * dt(1) &
                + xc%d3st2(ids%gatatb, i) * dt(1) * dt(2) &
                + xc%d3st2(ids%gatbtb, i) * dt(2) * dt(2)

          cs(2) = xc%d3st2(ids%gbtata, i) * dt(1) * dt(1) &
                + xc%d3st2(ids%gbtatb, i) * dt(2) * dt(1) &
                + xc%d3st2(ids%gbtatb, i) * dt(1) * dt(2) &
                + xc%d3st2(ids%gbtbtb, i) * dt(2) * dt(2)

          cs(3) = xc%d3st2(ids%gctata, i) * dt(1) * dt(1) &
                + xc%d3st2(ids%gctatb, i) * dt(2) * dt(1) &
                + xc%d3st2(ids%gctatb, i) * dt(1) * dt(2) &
                + xc%d3st2(ids%gctbtb, i) * dt(2) * dt(2)

          g_s = g_s + cs


          ct(1) = xc%d3r2t(ids%rarata, i) * dr(1)*dr(1) &
                + xc%d3r2t(ids%rarbta, i) * dr(1)*dr(2) &
                + xc%d3r2t(ids%rarbta, i) * dr(2)*dr(1) &
                + xc%d3r2t(ids%rbrbta, i) * dr(2)*dr(2)

          ct(2) = xc%d3r2t(ids%raratb, i) * dr(1)*dr(1) &
                + xc%d3r2t(ids%rarbtb, i) * dr(1)*dr(2) &
                + xc%d3r2t(ids%rarbtb, i) * dr(2)*dr(1) &
                + xc%d3r2t(ids%rbrbtb, i) * dr(2)*dr(2)

          g_t = g_t + ct


          ct(1) =  xc%d3rst(ids%ragata, i) * dr(1) * ds(1) &
                +  xc%d3rst(ids%rbgata, i) * dr(2) * ds(1) &

                +  xc%d3rst(ids%ragbta, i) * dr(1) * ds(2) &
                +  xc%d3rst(ids%rbgbta, i) * dr(2) * ds(2) &

                +  xc%d3rst(ids%ragcta, i) * dr(1) * ds(3) &
                +  xc%d3rst(ids%rbgcta, i) * dr(2) * ds(3)

          ct(2) =  xc%d3rst(ids%ragatb, i) * dr(1) * ds(1) &
                +  xc%d3rst(ids%rbgatb, i) * dr(2) * ds(1) &

                +  xc%d3rst(ids%ragbtb, i) * dr(1) * ds(2) &
                +  xc%d3rst(ids%rbgbtb, i) * dr(2) * ds(2) &

                +  xc%d3rst(ids%ragctb, i) * dr(1) * ds(3) &
                +  xc%d3rst(ids%rbgctb, i) * dr(2) * ds(3)

          g_t = g_t + 2*ct

          ct(1) =  xc%d3s2t(ids%gagata, i) * ds(1) * ds(1) &
                +  xc%d3s2t(ids%gagbta, i) * ds(1) * ds(2) &
                +  xc%d3s2t(ids%gagcta, i) * ds(1) * ds(3) &

                +  xc%d3s2t(ids%gagbta, i) * ds(2) * ds(1) &
                +  xc%d3s2t(ids%gbgbta, i) * ds(2) * ds(2) &
                +  xc%d3s2t(ids%gbgcta, i) * ds(2) * ds(3) &

                +  xc%d3s2t(ids%gagcta, i) * ds(3) * ds(1) &
                +  xc%d3s2t(ids%gbgcta, i) * ds(3) * ds(2) &
                +  xc%d3s2t(ids%gcgcta, i) * ds(3) * ds(3)


          ct(2) =  xc%d3s2t(ids%gagatb, i) * ds(1) * ds(1) &
                +  xc%d3s2t(ids%gagbtb, i) * ds(1) * ds(2) &
                +  xc%d3s2t(ids%gagctb, i) * ds(1) * ds(3) &

                +  xc%d3s2t(ids%gagbtb, i) * ds(2) * ds(1) &
                +  xc%d3s2t(ids%gbgbtb, i) * ds(2) * ds(2) &
                +  xc%d3s2t(ids%gbgctb, i) * ds(2) * ds(3) &

                +  xc%d3s2t(ids%gagctb, i) * ds(3) * ds(1) &
                +  xc%d3s2t(ids%gbgctb, i) * ds(3) * ds(2) &
                +  xc%d3s2t(ids%gcgctb, i) * ds(3) * ds(3)

          g_t = g_t + ct

          ct(1) = xc%d3rt2(ids%ratata, i) * dr(1)*dt(1) &
                + xc%d3rt2(ids%ratatb, i) * dr(1)*dt(2) &
                + xc%d3rt2(ids%rbtata, i) * dr(2)*dt(1) &
                + xc%d3rt2(ids%rbtatb, i) * dr(2)*dt(2)

          ct(2) = xc%d3rt2(ids%ratatb, i) * dr(1)*dt(1) &
                + xc%d3rt2(ids%ratbtb, i) * dr(1)*dt(2) &
                + xc%d3rt2(ids%rbtatb, i) * dr(2)*dt(1) &
                + xc%d3rt2(ids%rbtbtb, i) * dr(2)*dt(2)

          g_t = g_t + 2*ct

          ct(1) = xc%d3st2(ids%gatata, i) * ds(1)*dt(1) &
                + xc%d3st2(ids%gbtata, i) * ds(2)*dt(1) &
                + xc%d3st2(ids%gctata, i) * ds(3)*dt(1) &

                + xc%d3st2(ids%gatatb, i) * ds(1)*dt(2) &
                + xc%d3st2(ids%gbtatb, i) * ds(2)*dt(2) &
                + xc%d3st2(ids%gctatb, i) * ds(3)*dt(2)

          ct(2) = xc%d3st2(ids%gatatb, i) * ds(1)*dt(1) &
                + xc%d3st2(ids%gbtatb, i) * ds(2)*dt(1) &
                + xc%d3st2(ids%gctatb, i) * ds(3)*dt(1) &

                + xc%d3st2(ids%gatbtb, i) * ds(1)*dt(2) &
                + xc%d3st2(ids%gbtbtb, i) * ds(2)*dt(2) &
                + xc%d3st2(ids%gctbtb, i) * ds(3)*dt(2)

          g_t = g_t + 2*ct

          ct(1) = xc%d3t3(ids%tatata, i) * dt(1)*dt(1) &
                + xc%d3t3(ids%tatatb, i) * dt(1)*dt(2) &
                + xc%d3t3(ids%tatatb, i) * dt(2)*dt(1) &
                + xc%d3t3(ids%tatbtb, i) * dt(2)*dt(2)

          ct(2) = xc%d3t3(ids%tatatb, i) * dt(1)*dt(1) &
                + xc%d3t3(ids%tatbtb, i) * dt(1)*dt(2) &
                + xc%d3t3(ids%tatbtb, i) * dt(2)*dt(1) &
                + xc%d3t3(ids%tbtbtb, i) * dt(2)*dt(2)

          g_t = g_t + ct

        end if
      end if

    end associate

  end subroutine


!###############################################################################

  subroutine run_xc(xc_opts, xc_dat, basis)
    use basis_tools, only: basis_set
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num

    implicit none
    class(xc_consumer_t), intent(inout) :: xc_dat
    type(xc_options_t), intent(in) :: xc_opts
    type(basis_set), intent(in) :: basis

    type(xc_engine_t), allocatable :: xce

    logical :: skip


    real(KIND=fp) :: exc, totgradxyz(3), totele, totkin
    real(KIND=fp) :: dftthr, wcutoff

    integer :: next
    integer :: iSlice
    integer :: iAtom
    integer :: npt

    integer :: i
    integer :: numNzPts
    logical :: done

    integer :: myjob
    integer :: iChunk, chunkSize
    integer :: myThread, numThreads



    next = -1
    myjob = -1

    npt = xc_opts%molGrid%nMolPts
!   Set cut-offs for the weight WCUTOFF
!   WCUTOFF is a cell volume and we set it to a fixed value.
!   Most cells have large volume (about 97% have volume >1e-14)
    dftthr = 1.0d-04/npt
    if (xc_opts%dft_threshold > 0.0d0) dftthr = xc_opts%dft_threshold

    wcutoff = 1.0d-15
    if (dftthr > 1.1d-15) wcutoff = 1.0d-08/npt

    exc = 0
    totele = 0
    totgradxyz = 0
    totkin = 0

!$omp parallel &
!$omp   private(iChunk, chunkSize, numThreads, done) &
!$omp   private(iSlice, numNzPts, xce) &
!$omp   private(myThread), &
!$omp   private(i), &
!$omp   private(iAtom), &
!$omp   private(skip), &
!$omp   reduction(+:exc, totele, totgradxyz, totkin)

    numThreads = 1
    myThread = 1
!$  numThreads = omp_get_num_threads()
!$  myThread = omp_get_thread_num()+1
    chunkSize = max(1, xc_opts%molGrid%nSlices/(xc_dat%pe%size*4))
    if (chunkSize/numThreads > 40) then
      chunkSize = 40*numThreads
    end if
    allocate (xce)
    call xce%init(xc_opts)

!$omp master
    call xc_dat%parallel_start(xce, numThreads)
!$omp end master
!$omp barrier

    done = .false.

    do iChunk = 1, xc_opts%molGrid%nSlices, chunkSize
      if (mod(iChunk/chunkSize, xc_dat%pe%size) /= xc_dat%pe%rank) cycle

!$omp do schedule(dynamic)
      slc: do iSlice = iChunk, min(xc_opts%molGrid%nSlices, iChunk-1+chunkSize)

        call xc_opts%molgrid%getSliceNonZero(wcutoff, iSlice, xce%xyzw, numNzPts)
        if (numNzPts==0) CYCLE

        iAtom = xc_opts%molGrid%idOrigin(iSlice)

        call xce%resetPointers(numNzPts)

        do i = 1, numNzPts
          xce%xyzw(i,:3) = &
            xce%xyzw(i,:3) + basis%atoms%xyz(:3,iAtom)
        end do

        call xce%compAOs(basis, xce%nAODer, xce%xyzw(:numNzPts,:3))

        call xce%pruneAOs(skip)

        IF (skip) CYCLE

        call xce%compXC(xc_opts%functional, skip)

        IF (skip) CYCLE

        call xc_dat%update(xce, myThread)

        call xc_dat%postUpdate(xce, myThread)

      end do slc
!$omp end do nowait
    end do

    call xce%getStats( &
            E_xc=exc, &
            N_elec=totele, &
            E_kin=totkin, &
            G_total=totgradxyz)

    deallocate (xce)
!$omp end parallel

    call xc_dat%parallel_stop()
    call xc_dat%pe%allreduce(exc, 1)
    call xc_dat%pe%allreduce(totele, 1)
    call xc_dat%pe%allreduce(totkin, 1)
    call xc_dat%pe%allreduce(totgradxyz, 1)

    xc_dat%E_xc = exc
    xc_dat%N_elec = totele
    xc_dat%E_kin = totkin
    xc_dat%G_total = totgradxyz

  end subroutine

end module mod_dft_gridint
