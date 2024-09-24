module mod_dft_gridint_grad

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t
  use mod_dft_gridint, only: OQP_FUNTYP_LDA, OQP_FUNTYP_MGGA
  use mod_dft_gridint, only: compAtGradRho, compAtGradDRho, compAtGradTau

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_grad_t
    real(kind=fp), allocatable :: bfgrad(:,:,:)
    real(kind=fp), allocatable :: tmp_(:,:)
    real(kind=fp), allocatable :: d1dsx(:,:,:) !< Temporary storage for dE/d\sigma
  contains
    procedure :: parallel_start
    procedure :: parallel_stop
    procedure :: resetGradPointers
    procedure :: update
    procedure :: postUpdate
    procedure :: clean
  end type

!-------------------------------------------------------------------------------

  private
  public derexc_blk

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

  subroutine parallel_start(self, xce, nthreads)
    implicit none
    class(xc_consumer_grad_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    call self%clean()
    allocate( self%bfGrad(xce%numAOs, 3, nthreads) &
            , self%d1dsx(xce%maxPts, 3, nthreads) &
            , self%tmp_(xce%numAOs*3, nthreads) &
            , source=0.0d0)
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_grad_t), intent(inout) :: self
    if (ubound(self%bfGrad,3) /= 1) then
      self%bfGrad(:,:,lbound(self%bfGrad,3)) = sum(self%bfGrad, dim=3)
    end if

    call self%pe%allreduce(self%bfGrad(:,:,1), &
              size(self%bfGrad(:,:,1)))
  end subroutine

!-------------------------------------------------------------------------------

  subroutine clean(self)
    implicit none
    class(xc_consumer_grad_t), intent(inout) :: self
    if (allocated(self%bfGrad)) deallocate(self%bfGrad)
    if (allocated(self%d1dsx)) deallocate(self%d1dsx)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
  end subroutine

!-------------------------------------------------------------------------------
!> @brief Adjust internal memory storage for a given
!>  number of pruned grid points
!> @author Konstantin Komarov
 subroutine resetGradPointers(self, xce, tmp, myThread)
    class(xc_consumer_grad_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer :: tmp(:,:)
    integer, intent(in) :: myThread

!   pruned AOs or no pruned AOs
    associate ( numAOs => xce%numAOs_p &  ! number of pruned AOs
      )
      tmp(1:numAOs, 1:3) => self%tmp_(1:numAOs*3, myThread)
    end associate

 end subroutine

!-------------------------------------------------------------------------------

 subroutine update(self, xce, mythread)

    class(xc_consumer_grad_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread, i
    real(kind=fp), pointer :: tmpGrad(:,:)

    call self%resetGradPointers(xce, tmpGrad,  myThread)

    associate ( bfgrad  => self%bfgrad(:,:,mythread) &
              , d1dsx   => self%d1dsx(:,:,mythread) &
              , aoG1    => xce%aoG1 &
              , aoG2    => xce%aoG2 &
              , moVA    => xce%moVA &
              , moVB    => xce%moVB &
              , moG1A   => xce%moG1A &
              , moG1B   => xce%moG1B &
              , hasBeta => xce%hasBeta &
              , numPts  => xce%numPts &
              , xc      => xce%XCLib &
              , drho    => xce%xclib%drho  &
              , ids     => xce%XCLib%ids &
              , d1ds    => xce%XCLib%d1ds &
              , d1dr    => xce%XCLib%d1dr &
              , d1dt    => xce%XCLib%d1dt &
      )
      tmpGrad = 0.0d0

!     LDA gradient
      call compAtGradRho(tmpGrad, d1dr(1,:), moVA, aoG1, numPts)

!     GGA gradient
      if (xce%funTyp /= OQP_FUNTYP_LDA) then
          do i = 1, numPts
            d1dsx(i,1:3) = 2*d1ds(ids%ga,i)*drho(1:3,i)+d1ds(ids%gc,i)*drho(4:6,i)
          end do

          call compAtGradDRho(tmpGrad, d1dsx, moVA, moG1A, aoG1, aoG2, numPts)
      end if

!     Meta-GGA gradient
      if (xce%funTyp == OQP_FUNTYP_MGGA) &
        call compAtGradTau(tmpGrad, d1dt(1,:), moG1A, aoG2, numPts)


      if (hasBeta) then
        call compAtGradRho(tmpGrad, d1dr(2,:), moVB, aoG1, numPts)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then
            do i = 1, numPts
              d1dsx(i,1:3) = 2*d1ds(ids%gb,i)*drho(4:6,i)+d1ds(ids%gc,i)*drho(1:3,i)
            end do
            call compAtGradDRho(tmpGrad, d1dsx, moVB, moG1B, aoG1, aoG2, numPts)
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) &
          call compAtGradTau(tmpGrad, d1dt(2,:), moG1B, aoG2, numPts)
      end if

    end associate
 end subroutine

 subroutine postUpdate(self, xce, mythread)

    class(xc_consumer_grad_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(kind=fp), pointer :: tmpGrad(:,:)

    call self%resetGradPointers(xce, tmpGrad,  myThread)

    associate ( numAOs  => xce%numAOs_p &  ! number of pruned AOs
              , indices => xce%indices_p &
      )

      if (xce%skip_p) then
        self%bfGrad(:,:,myThread) = self%bfGrad(:,:,myThread) + tmpGrad
      else
        self%bfGrad(indices(1:numAOs), :, mythread) = &
          self%bfGrad(indices(1:numAOs), :, mythread) + tmpGrad(1:numAOs, :)
      end if

   end associate

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute grid XC contribution to the nuclear gradient
!> @note  Weight derivatives are not applied here. The gradient seems
!>  to be good enough even using fairly poor grids. However, I do
!>  not recomment to use it for numerical Hessian calculation until
!>  weight derivatives are implemented
!> @param[in]    da        density matrix, alpha-spin
!> @param[in]    db        density matrix, beta-spin
!> @param[inout] dedft     nuclear gradient
!> @param[out]   totele    electronic denisty integral
!> @param[out]   totkin    kinetic energy integral
!> @param[in]    mxAngMom  max. needed ang. mom. value (incl. derivatives)
!> @param[in]    nbf        basis set size
!> @param[in]    isGGA     .TRUE. if GGA/mGGA functional used
!> @param[in]    urohf     .TRUE. if open-shell calculation
!> @author Vladimir Mironov
  subroutine derexc_blk(basis, molGrid, da, db, dedft, &
                        totele, totkin, &
                        mxAngMom, nbf, dft_threshold, urohf, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid

    type(basis_set) :: basis
    logical, intent(IN) :: urohf
    integer, intent(IN) :: mxAngMom, nbf
    real(KIND=fp), intent(INOUT) :: totele, totkin
    real(KIND=fp), intent(INOUT) :: da(nbf, *), db(nbf, *), dedft(:, :)
    real(kind=fp), intent(in) :: dft_threshold

    type(xc_consumer_grad_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: j

    integer :: nat

    real(KIND=fp), target, allocatable :: da2(:, :), db2(:, :)


    nat = infos%mol_prop%natom

    allocate (da2(nbf, nbf))
    do j = 1, nbf
      da2(:, j) = da(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
    end do
    if (urohf) then
      allocate (db2(nbf, nbf))
      do j = 1, nbf
        db2(:, j) = db(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
      end do
    end if

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%functional => infos%functional
    xc_opts%hasBeta = urohf
    xc_opts%isWFVecs = .false.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = mxAngMom
    xc_opts%nDer = 1
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => da2
    xc_opts%wfBeta => db2
    xc_opts%dft_threshold = dft_threshold
    xc_opts%molGrid => molGrid

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_xc(xc_opts, dat, basis)

    totele = dat%N_elec
    totkin = dat%E_kin

    do j = 1, basis%nshell
      associate (atom => basis%origin(j), &
                 offset => basis%ao_offset(j), &
                 naos => basis%naos(j))
        dedft(1, atom) = dedft(1, atom)-sum(dat%bfGrad(offset:offset+naos-1, 1, 1))
        dedft(2, atom) = dedft(2, atom)-sum(dat%bfGrad(offset:offset+naos-1, 2, 1))
        dedft(3, atom) = dedft(3, atom)-sum(dat%bfGrad(offset:offset+naos-1, 3, 1))
      end associate
    end do

    deallocate (da2)
    if (urohf) deallocate (db2)

    call dat%clean()

  end subroutine

!-------------------------------------------------------------------------------

end module mod_dft_gridint_grad
