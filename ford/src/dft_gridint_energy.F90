module mod_dft_gridint_energy

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t, OQP_FUNTYP_LDA, OQP_FUNTYP_MGGA
  use oqp_linalg

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_ks_t
    real(kind=fp), allocatable :: fa2(:,:)
    real(kind=fp), allocatable :: fb2(:,:)
    real(kind=fp), allocatable :: focks_(:,:)
    real(kind=fp), allocatable :: tmp_(:,:)
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
  public dmatd_blk

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

  subroutine parallel_start(self, xce, nthreads)
    implicit none
    class(xc_consumer_ks_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: nSpin
    call self%clean()
    nSpin = 1
    if (xce%hasBeta) nSpin = 2
    allocate( self%fa2(xce%numAOs*xce%numAOs,nthreads) &
            , self%focks_(xce%numAOs*xce%numAOs*nSpin,nthreads) &
            , self%tmp_(xce%numAOs*xce%maxPts*xce%numTmpVec,nthreads) &
            , source=0.0d0)
    if (xce%hasBeta) then
      allocate(self%fb2(xce%numAOs*xce%numAOs,nthreads), source=0.0d0)
    end if
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_ks_t), intent(inout) :: self
    if (ubound(self%fa2,2) /= 1 ) then
      self%fa2(:,lbound(self%fa2,2)) = sum(self%fa2, dim=size(shape(self%fa2)))
    end if
    call self%pe%allreduce(self%fa2(:,1), &
              size(self%fa2(:,1)))
    if (allocated(self%fb2)) then
      if (ubound(self%fa2,2) /= 1 ) then
        self%fb2(:,lbound(self%fb2,2)) = sum(self%fb2, dim=size(shape(self%fb2)))
      end if
      call self%pe%allreduce(self%fb2(:,1), &
                size(self%fb2(:,1)))
    end if
  end subroutine

!-------------------------------------------------------------------------------

  subroutine clean(self)
    implicit none
    class(xc_consumer_ks_t), intent(inout) :: self
    if (allocated(self%fa2)) deallocate(self%fa2)
    if (allocated(self%fb2)) deallocate(self%fb2)
    if (allocated(self%focks_)) deallocate(self%focks_)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
  end subroutine

!-------------------------------------------------------------------------------
!> @brief Adjust internal memory storage for a given
!>  number of pruned grid points
!> @author Konstantin Komarov
 subroutine resetOrbPointers(self, xce, focks, tmp, fock_a, fock_b, myThread)
    class(xc_consumer_ks_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer :: focks(:,:,:)
    real(kind=fp), intent(out), pointer, optional :: tmp(:,:,:)
    real(kind=fp), intent(out), pointer, optional :: fock_a(:,:)
    real(kind=fp), intent(out), pointer, optional :: fock_b(:,:)
    integer, intent(in) :: myThread

    integer :: nSpin

    associate ( numAOs   => xce%numAOs &
              , numAOs_p => xce%numAOs_p &
              , numPts   => xce%numPts &
              , TmpVec => xce%numTmpVec &
              , hasBeta => xce%hasBeta &
      )
      nSpin = 1
      if (hasBeta) nSpin = 2

      focks(1:numAOs_p,1:numAOs_p, 1:nSpin) => self%focks_(1:numAOs_p*numAOs_p*nSpin, myThread)

      if (present(tmp)) &
        tmp(1:numAOs_p, 1:numPts, 1:TmpVec) => self%tmp_(1:numAOs_p*numPts*TmpVec, myThread)

      if (present(fock_a)) fock_a(1:numAOs,1:numAOs) => self%fa2(1:numAOs*numAOs, myThread)
      if (present(fock_b)) fock_b(1:numAOs,1:numAOs) => self%fb2(1:numAOs*numAOs, myThread)

    end associate

 end subroutine

!-------------------------------------------------------------------------------

  subroutine update(self, xce, mythread)
    implicit none
    class(xc_consumer_ks_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    integer :: i, j
    real(kind=fp) :: c(3)
    real(kind=fp), pointer :: focks(:,:,:)
    real(kind=fp), pointer :: tmp(:,:,:)

    call self%resetOrbPointers(xce, focks=focks, tmp=tmp, myThread=myThread)

    associate ( d1dr => xce%XCLib%d1dr &
              , d1ds => xce%XCLib%d1ds &
              , d1dt => xce%XCLib%d1dt &
              , ra   => xce%XCLib%ids%ra &
              , rb   => xce%XCLib%ids%rb &
              , ga   => xce%XCLib%ids%ga &
              , gb   => xce%XCLib%ids%gb &
              , gc   => xce%XCLib%ids%gc &
              , ta   => xce%XCLib%ids%ta &
              , tb   => xce%XCLib%ids%tb &
              , drho => xce%XCLib%drho &
              , aoV    => xce%aoV &
              , aoG1   => xce%aoG1 &
              , numAOs => xce%numAOs_p &  ! number of pruned AOs
              , numPts => xce%numPts &
      )

      ! LDA case
      do i = 1, numPts
        tmp(:, i, 1) = 0.5*d1dr(ra, i)*aoV(:, i)
      end do

      ! GGA case
      if (xce%funTyp /= OQP_FUNTYP_LDA) then
        do i = 1, numPts
          c = 2*d1ds(ga,i)*drho(1:3,i) + d1ds(gc,i)*drho(4:6,i)
          tmp(:,i,1) = tmp(:,i,1) &
            +c(1)*aoG1(:, i, 1) &
            +c(2)*aoG1(:, i, 2) &
            +c(3)*aoG1(:, i, 3)
        end do
      end if

      call dsyr2k('U', 'N', numAOs, numPts, 1.0_fp, &
                          aoV, numAOs, &
                          tmp(:,:,1), numAOs, &
                  0.0_fp, focks(:,:,1), numAOs)

      ! metaGGA case
      if (xce%funTyp == OQP_FUNTYP_MGGA) then
        do i = 1, numPts
          tmp(:,i,2:4) = d1dt(ta,i)*aoG1(:,i,1:3)
        end do

        do j = 1, 3
          call dsyr2k('U', 'N', numAOs, numPts, 0.25_fp, &
                              aoG1(:,:,j), numAOs, &
                              tmp(:,:,j+1), numAOs, &
                      1.0_fp, focks(:,:,1), numAOs)
        end do
      end if

      if (xce%hasBeta) then

        ! LDA case
        do i = 1, numPts
          tmp(:, i, 1) = 0.5*d1dr(rb, i)*aoV(:, i)
        end do

        ! GGA case
        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          do i = 1, numPts
            c = 2*d1ds(gb,i)*drho(4:6,i) + d1ds(gc,i)*drho(1:3,i)
            tmp(:,i,1) = tmp(:,i,1) &
              +c(1)*aoG1(:, i, 1) &
              +c(2)*aoG1(:, i, 2) &
              +c(3)*aoG1(:, i, 3)
          end do
        end if

        call dsyr2k('U', 'N', numAOs, numPts, 1.0_fp, &
                           aoV, numAOs, &
                           tmp(:,:,1), numAOs, &
                   0.0_fp, focks(:,:,2), numAOs)

        ! metaGGA case
        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          do i = 1, numPts
            tmp(:,i,2:4) = d1dt(tb,i)*aoG1(:,i,1:3)
          end do

          do j = 1, 3
            call dsyr2k('U', 'N', numAOs, numPts, 0.25_fp, &
                                aoG1(:,:,j), numAOs, &
                                tmp(:,:,j+1), numAOs, &
                        1.0_fp, focks(:,:,2), numAOs)
          end do
        end if
      end if

    end associate

  end subroutine

  subroutine postUpdate(self, xce, mythread)
    implicit none
    class(xc_consumer_ks_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(kind=fp), pointer :: focks(:,:,:)
    real(kind=fp), pointer :: fock_a(:,:)
    real(kind=fp), pointer :: fock_b(:,:)

    call self%resetOrbPointers(xce, focks=focks, fock_a=fock_a, &
                               fock_b=fock_b, myThread=myThread)

    associate(  numAOs  => xce%numAOs_p &  ! number of pruned AOs
              , indices => xce%indices_p &
      )

      if (xce%skip_p) then

         fock_a = fock_a + focks(:,:,1)
         if (xce%hasBeta) &
            fock_b = fock_b + focks(:,:,2)

      else

         fock_a(indices(1:numAOs), &
                indices(1:numAOs)) = fock_a(indices(1:numAOs), &
                                            indices(1:numAOs)) &
                                   + focks(:,:,1)
         if (xce%hasBeta) &
            fock_b(indices(1:numAOs), &
                   indices(1:numAOs)) = fock_b(indices(1:numAOs), &
                                               indices(1:numAOs)) &
                                      + focks(:,:,2)

      end if

    end associate

  end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute grid XC contribution to the Kohn-Sham matrix
!> @param[in]    coeffa     MO coefficients, alpha-spin
!> @param[in]    coeffb     MO coefficients, beta-spin
!> @param[inout] fa         KS matrix, alpha-spin
!> @param[inout] fb         KS matrix, beta-spin
!> @param[out]   exc        XC energy
!> @param[out]   totele     electronic denisty integral
!> @param[out]   totkin     kinetic energy integral
!> @param[in]    mxAngMom   max. needed ang. mom. value (incl. derivatives)
!> @param[in]    nbf         basis set size
!> @param[in]    urohf      .TRUE. if open-shell calculation
!> @author Vladimir Mironov
  subroutine dmatd_blk(basis, molGrid, coeffa, coeffb, fa, fb, &
                       exc, totele, totkin, &
                       mxAngMom, nbf, dft_threshold, urohf, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use mod_dft_molgrid, only: dft_grid_t
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information

    implicit none

    type(dft_grid_t), target, intent(in) :: molGrid
    type(information), target, intent(in) :: infos

    type(basis_set) :: basis
    logical, intent(IN) :: urohf
    integer, intent(IN) :: mxAngMom, nbf
    real(kind=fp), intent(inout) :: exc, totele, totkin
    real(kind=fp), target, intent(inout) :: coeffa(nbf, *), coeffb(nbf, *)
    real(kind=fp), intent(inout) :: fa(*), fb(*)
    real(kind=fp), intent(in) :: dft_threshold

    type(xc_consumer_ks_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: i0, i, j

    do j = 1, nbf
      coeffa(:, j) = coeffa(:, j)*basis%bfnrm(:)
    end do
    if (urohf) then
      do j = 1, nbf
        coeffb(:, j) = coeffb(:, j)*basis%bfnrm(:)
      end do
    end if

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%functional => infos%functional
    xc_opts%hasBeta = urohf
    xc_opts%isWFVecs = .true.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = mxAngMom
    xc_opts%nDer = 0
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => coeffa(:,1:nbf)
    xc_opts%wfBeta => coeffb(:,1:nbf)
    xc_opts%molGrid => molGrid
    xc_opts%dft_threshold = dft_threshold
    xc_opts%ao_threshold = infos%dft%grid_ao_threshold
    xc_opts%ao_sparsity_ratio = infos%dft%grid_ao_sparsity_ratio
    ! skip ao_prune_grid if it is pruned grid (SG1)
    if(infos%dft%grid_pruned) xc_opts%ao_sparsity_ratio = 0.0_fp

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_xc(xc_opts, dat, basis)

    exc = dat%E_xc
    totele = dat%N_elec
    totkin = dat%E_kin

!   Normalize KS matrices
    do j = 1, nbf
      coeffa(:, j) = coeffa(:, j)/basis%bfnrm(:)
    end do

    i0 = 0
    do i = 1, nbf
      fa(i0+1:i0+i) = fa(i0+1:i0+i)+ &
                      dat%fa2((i-1)*nbf+1:(i-1)*nbf+i,1) &
                      *basis%bfnrm(i) &
                      *basis%bfnrm(1:i)
      i0 = i0+i
    end do

    if (urohf) then
      do j = 1, nbf
        coeffb(:, j) = coeffb(:, j)/basis%bfnrm(:)
      end do

      i0 = 0
      do i = 1, nbf
        fb(i0+1:i0+i) = fb(i0+1:i0+i)+ &
                        dat%fb2((i-1)*nbf+1:(i-1)*nbf+i,1) &
                        *basis%bfnrm(i) &
                        *basis%bfnrm(1:i)
        i0 = i0+i
      end do
    end if

    call dat%clean()

  end subroutine

!-------------------------------------------------------------------------------

end module mod_dft_gridint_energy
