module mod_dft_gridint_fxc

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t
  use mod_dft_gridint, only: X__, Y__, Z__
  use mod_dft_gridint, only: OQP_FUNTYP_LDA, OQP_FUNTYP_GGA, OQP_FUNTYP_MGGA
  use oqp_linalg

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_tde_t
    integer :: nMtx = 1
    real(kind=fp), pointer :: da(:,:,:) => null()
    real(kind=fp), pointer :: db(:,:,:) => null()
    real(kind=fp), allocatable :: focks(:,:,:,:,:)
    real(kind=fp), allocatable :: mo(:,:,:,:,:)
    real(kind=fp), allocatable :: rRho(:,:,:,:)
    real(kind=fp), allocatable :: drRho(:,:,:,:,:)
    real(kind=fp), allocatable :: rTau(:,:,:,:)
!   Temporary storage
    real(kind=fp), allocatable :: focks_(:,:)
    real(kind=fp), allocatable :: tmpMO_(:,:)
    real(kind=fp), allocatable :: tmpDensity_(:,:,:)
    real(kind=fp), allocatable :: moG1_(:,:)
    real(kind=fp), allocatable :: tmp_(:,:)
  contains
    procedure :: parallel_start
    procedure :: parallel_stop
    procedure :: update
    procedure :: postUpdate
    procedure :: clean
    procedure :: resetXCPointers
    procedure :: computeRAll
    procedure :: resetOrbPointers
    procedure :: RUpdate
    procedure :: UUpdate
  end type

!-------------------------------------------------------------------------------

  private
  public tddft_fxc
  public utddft_fxc
  public xc_consumer_tde_t

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

  subroutine parallel_start(self, xce, nThreads)
    implicit none
    class(xc_consumer_tde_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nThreads
    integer :: nSpin
    call self%clean()

    nSpin = 1
    if (xce%hasBeta) nSpin = 2

    allocate( &
        self%focks(xce%numAOs, xce%numAOs, self%nMtx, nSpin, nThreads) &
      , self%mo(xce%numAOs, xce%maxPts, self%nMtx, nSpin, nThreads) &
      , self%rRho(nSpin, xce%maxPts, self%nMtx, nThreads) &
      , self%drRho(4, nSpin, xce%maxPts, self%nMtx, nThreads) &
!       Temporary storage
      , self%focks_(xce%numAOs * xce%numAOs * self%nMtx * nSpin, nThreads) &
      , self%tmpMO_(xce%numAOs * xce%maxPts * self%nMtx * nSpin, nThreads) &
      , self%tmpDensity_(xce%numAOs * xce%numAOs * self%nMtx, nSpin, nThreads) &
      , self%tmp_(xce%numAOs * xce%maxPts * xce%numTmpVec, nthreads) &
      , source=0.0d0)

    if (xce%funTyp == OQP_FUNTYP_MGGA) then
        allocate( &
            self%rTau(nSpin, xce%maxPts, self%nMtx, nThreads) &
!       Temporary storage
          , self%moG1_(xce%numAOs*xce%maxPts*3*self%nMtx, nThreads) &
          , source=0.0d0)
    end if
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_tde_t), intent(inout) :: self

    if (ubound(self%focks,5) /= 1) then
      self%focks(:,:,:,:,lbound(self%focks,5)) = sum(self%focks, dim=5)
    end if
    call self%pe%allreduce(self%focks(:,:,:,:,1), &
              size(self%focks(:,:,:,:,1)))
  end subroutine

!-------------------------------------------------------------------------------

  subroutine clean(self)
    implicit none
    class(xc_consumer_tde_t), intent(inout) :: self
    if (allocated(self%focks)) deallocate(self%focks)
    if (allocated(self%mo)) deallocate(self%mo)
    if (allocated(self%rRho)) deallocate(self%rRho)
    if (allocated(self%drRho)) deallocate(self%drRho)
    if (allocated(self%rTau)) deallocate(self%rTau)
!       Temporary storage
    if (allocated(self%focks_)) deallocate(self%focks_)
    if (allocated(self%tmpMO_)) deallocate(self%tmpMO_)
    if (allocated(self%tmpDensity_)) deallocate(self%tmpDensity_)
    if (allocated(self%moG1_)) deallocate(self%moG1_)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
  end subroutine

!-------------------------------------------------------------------------------

 subroutine update(self, xce, myThread)

    class(xc_consumer_tde_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: myThread

    call self%computeRAll(xce, myThread)

    if (xce%hasBeta) then
      call self%UUpdate(xce, myThread)
    else
      call self%RUpdate(xce, myThread)
    end if
 end subroutine

 subroutine postUpdate(self, xce, myThread)

    class(xc_consumer_tde_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: myThread

    real(kind=fp), pointer :: focks(:,:,:,:)

    call self%resetOrbPointers(xce, focks=focks, myThread=myThread)

    associate(  numAOs  => xce%numAOs_p &  ! number of pruned AOs
              , indices => xce%indices_p &
      )
      if (xce%skip_p) then

         self%focks(:,:,:,:,myThread) = self%focks(:,:,:,:,myThread) + focks

      else

         self%focks(indices(1:numAOs), indices(1:numAOs), :, :, myThread) &
            = self%focks(indices(1:numAOs), indices(1:numAOs), :, :, myThread) &
            + focks

      end if

    end associate

 end subroutine

!-------------------------------------------------------------------------------
!> @brief Adjust internal memory storage for a given
!>  number of pruned grid points
!> @author Konstantin Komarov
 subroutine resetXCPointers(self, xce, da, db, mo, moG1, myThread)
    class(xc_consumer_tde_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer :: da(:,:,:)
    real(kind=fp), intent(out), pointer :: db(:,:,:)
    real(kind=fp), intent(out), pointer :: mo(:,:,:,:)
    real(kind=fp), intent(out), pointer :: moG1(:,:,:,:)
    integer, intent(in) :: myThread
    integer :: nSpin

    if (xce%skip_p) then

!     no pruned AOs
      da => self%da
      if (xce%hasBeta) db => self%db
      mo => self%mo(:,:,:,:,myThread)

    else

!     pruned AOs
      nSpin = 1
      if (xce%hasBeta) nSpin = 2

      associate ( indices => xce%indices_p &
                , numAOs  => xce%numAOs_p & ! number of pruned AOs
                , numPts  => xce%numPts &
                , nMtx    => self%nMtx)
        ! Set pointer for Dens A
        da(1:numAOs, 1:numAOs, 1:nMtx) => &
              self%tmpDensity_(1:numAOs*numAOs*nMtx, 1, myThread)
        ! Compress Dens A
        da(1:numAOs, 1:numAOs,:) = &
              self%da(indices(1:numAOs), indices(1:numAOs), :)
        ! Set pointer for MOs
        da(1:numAOs, 1:numAOs, 1:nMtx) => &
              self%tmpDensity_(1:numAOs*numAOs*nMtx, 1, myThread)
        mo(1:numAOs, 1:numPts, 1:nMtx, 1:nSpin) => &
              self%tmpMO_(1:numAOs*numPts*nMtx*nSpin, myThread)

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          ! Set pointer for tmp mGGA
          moG1(1:numAOs, 1:numPts, 1:3, 1:nMtx) => &
              self%moG1_(1:numAOs*numPts*3*nMtx, myThread)
        end if

        if (xce%hasBeta) then
            ! Set pointer for Dens B
            db(1:numAOs, 1:numAOs, 1:nMtx) => &
                  self%tmpDensity_(1:numAOs*numAOs*nMtx, 2, myThread)
            ! Compress Dens B
            db(1:numAOs, 1:numAOs, :) = &
                  self%db(indices(1:numAOs), indices(1:numAOs), :)
        end if

      end associate

    end if
    ! In any case set pointer for moG1
    if (xce%funTyp == OQP_FUNTYP_MGGA) then
      ! Set pointer for tmp mGGA
      moG1(1:xce%numAOs_p, 1:xce%numPts, 1:3, 1:self%nMtx) => &
          self%moG1_(1:xce%numAOs_p*xce%numPts*3*self%nMtx, myThread)
    end if


 end subroutine

 subroutine resetOrbPointers(self, xce, focks, tmp, myThread)
    class(xc_consumer_tde_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer :: focks(:,:,:,:)
    real(kind=fp), intent(out), pointer, optional :: tmp(:,:,:)
    integer, intent(in) :: myThread
    integer :: nSpin

    nSpin = 1
    if (xce%hasBeta) nSpin = 2
!   pruned AOs or no pruned AOs
    associate ( numAOs => xce%numAOs_p &  ! number of pruned AOs
      )
      focks(1:numAOs, 1:numAOs, 1:self%nMtx, 1:nSpin) => &
         self%focks_(1:numAOs * numAOs * self%nMtx * nSpin, myThread)

      if (present(tmp)) &
        tmp(1:numAOs, 1:xce%numPts, 1:xce%numTmpVec) => &
          self%tmp_(1:numAOs * xce%numPts * xce%numTmpVec, myThread)
    end associate

 end subroutine

!> @brief Compute required MO-like values as well as \rho, \sigma, and \tau
 subroutine computeRAll(self, xce, myThread)

    class(xc_consumer_tde_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: myThread
    real(kind=fp), pointer :: da(:,:,:)
    real(kind=fp), pointer :: db(:,:,:)
    real(kind=fp), pointer :: mo(:,:,:,:)
    real(kind=fp), pointer :: moG1(:,:,:,:)

    call self%resetXCPointers( &
             xce, da, db, mo, moG1, myThread)

    associate ( hasBeta => xce%hasBeta &
              , numAOs  => xce%numAOs_p &  ! number of pruned AOs
              , rrho    => self%rrho(:,:,:,mythread) &
              , drrho   => self%drrho(:,:,:,:,mythread) &
              , rtau    => self%rtau(:,:,:,mythread) &
              , indices => xce%indices_p &
      )

      call xce%compRMOs(da, mo(:,:,:,1))
      if (hasBeta) then
        call xce%compRMOs(db, mo(:,:,:,2))
      end if

      if (.not. xce%skip_p) self%mo(indices(1:numAOs),:,:,:,myThread) = mo(1:numAOs,:,:,:)

      call xce%compRRho(mo, rRho)
      if (xce%funTyp /= OQP_FUNTYP_LDA) then
        call xce%compRDRho(mo, drRho)
      end if

      if (xce%funTyp == OQP_FUNTYP_MGGA) then
        call xce%compRMOGs(da, moG1)
        call compRTau(xce, moG1, rTau, 1)
        if (xce%hasBeta) then
          call xce%compRMOGs(db, moG1)
          call compRTau(xce, moG1, rTau, 2)
        end if
      end if
    end associate

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute Tau: (MO)' times (AO)'
!> @param[in]  xce      XC engine, parameters
!> @param[in]  moG1     MO directional derivatives
!> @param[out] rTau     kinetic energy density
!> @author Vladimir Mironov
 subroutine compRTau(xce, moG1, rTau, nSpin)

    class(xc_engine_t) :: xce
    real(kind=fp), intent(out) :: rTau(:,:,:)
    real(kind=fp), intent(in) :: moG1(:,:,:,:)
    integer :: i, j, nSpin, nMtx

    nMtx = ubound(moG1, 4)

    do j = 1, nMtx
      do i = 1, xce%numPts
        rTau(nSpin,i,j) = 0.5*sum(xce%aoG1(:,i,1:3)*moG1(:,i,1:3,j))
      end do
    end do

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Add XC derivative contribution to the Kohn-Sham-like matrices
!> @param[inout] fa   Kohn-Sham matrices
!> @author Vladimir Mironov
 subroutine UUpdate(dat, xce, myThread)
    use mod_dft_gridint, only: xc_der1, xc_der2_contr
    class(xc_engine_t) :: xce
    class(xc_consumer_tde_t) :: dat
    integer :: myThread

    integer :: i, j, k
    real(kind=fp) :: cs(3)
    real(kind=fp) :: d_r(2), d_s(3), d_t(2)
    real(kind=fp) :: f_r(2), f_s(3), f_t(2)
    real(kind=fp) :: rhoab(2), tauab(2), Sigma(3), dsaa, dsab, dsba, dsbb
    real(kind=fp), pointer :: focks(:,:,:,:)
    real(kind=fp), pointer :: tmp(:,:,:)

    call dat%resetOrbPointers(xce, focks, tmp, myThread)

    associate ( aoV    => xce%aoV &
              , aoG1   => xce%aoG1 &
              , rRho   => dat%rRho(:,:,:,myThread)  &
              , drRho  => dat%drRho(:,:,:,:,myThread)  &
              , rTau   => dat%rTau(:,:,:,myThread)  &
              , dRho   => xce%xclib%dRho  &
              , numAOs => xce%numAOs_p &  ! number of pruned AOs
              , numPts => xce%numPts &
              , nMtx   => dat%nMtx &
              , xc     => xce%XCLib &
              , ids    => xce%XCLib%ids &
      )
      do j = 1, nMtx

        do i = 1, numPts

          rhoab = rRho(1:2,i,j)

          Sigma = 0
          tauab = 0

          if (xce%funTyp /= OQP_FUNTYP_LDA) then
            dsaa = dot_product(drRho(:,1,i,j), dRho(1:3,i))
            dsab = dot_product(drRho(:,1,i,j), dRho(4:6,i))
            dsbb = dot_product(drRho(:,2,i,j), dRho(4:6,i))
            dsba = dot_product(drRho(:,2,i,j), dRho(1:3,i))
            Sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]
          end if

          if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rTau(1:2,i,j)

          call xc_der1(xce, .true., i, d_r, d_s, d_t)
          call xc_der2_contr(xce, .true., i, &
                  rhoab, Sigma, tauab, &
                  f_r, f_s, f_t)

          if (maxval(abs([dsaa,dsbb,dsab,dsba])) < xce%threshold) then
            f_s = 0
          end if
          if (maxval(abs(tauab)) < xce%threshold) then
            f_t = 0
          end if

          tmp(:,i,1) = 0.5_fp*f_r(1)*aoV(:,i)

          if (xce%funTyp /= OQP_FUNTYP_LDA) then
            cs =    2*f_s(1) * dRho(1:3,i) &
                  +   f_s(3) * dRho(4:6,i) &
                  + 2*d_s(1) * drRho(:,1,i,j) &
                  +   d_s(3) * drRho(:,2,i,j)
            tmp(:,i,1) = tmp(:,i,1) &
                   +   cs(X__)*aoG1(:,i,X__) &
                   +   cs(Y__)*aoG1(:,i,Y__) &
                   +   cs(Z__)*aoG1(:,i,Z__)
          end if

          if (xce%funTyp == OQP_FUNTYP_MGGA) then
            tmp(:,i,2:4) = f_t(1)*aoG1(:,i,X__:Z__)
          end if

        end do

        call dsyr2k('u', 'n', numAOs, numPts, 1.0_fp, &
                             aoV, numAOs, &
                             tmp(:,:,1), numAOs, &
                     0.0_fp, focks(:,:,j,1), numAOs)

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          do k = 1, 3
              call dsyr2k('u', 'n', numAOs, numPts, 0.25_fp, &
                                   aoG1(:,:,k), numAOs, &
                                   tmp(:,:,k+1), numAOs, &
                           1.0_fp, focks(:,:,j,1), numAOs)
          end do
        end if

        do i = 1, numPts

          rhoab = rRho(1:2,i,j)

          Sigma = 0
          tauab = 0

          if (xce%funTyp /= OQP_FUNTYP_LDA) then
            dsaa = dot_product(drRho(:,1,i,j), dRho(1:3,i))
            dsab = dot_product(drRho(:,1,i,j), dRho(4:6,i))
            dsbb = dot_product(drRho(:,2,i,j), dRho(4:6,i))
            dsba = dot_product(drRho(:,2,i,j), dRho(1:3,i))
            Sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]
          end if

          if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rTau(1:2,i,j)

          call xc_der1(xce, .true., i, d_r, d_s, d_t)
          call xc_der2_contr(xce, .true., i, &
                  rhoab, Sigma, tauab, &
                  f_r, f_s, f_t)

          if (maxval(abs([dsaa,dsbb,dsab,dsba])) < xce%threshold) then
            f_s = 0
          end if
          if (maxval(abs(tauab)) < xce%threshold) then
            f_t = 0
          end if

          tmp(:,i,1) = 0.5_fp*f_r(2)*aoV(:,i)

          if (xce%funTyp /= OQP_FUNTYP_LDA) then
            cs =    2*f_s(2) * dRho(4:6,i) &
                  +   f_s(3) * dRho(1:3,i) &
                  + 2*d_s(2) * drRho(:,2,i,j) &
                  +   d_s(3) * drRho(:,1,i,j)
            tmp(:,i,1) = tmp(:,i,1) &
                   +   cs(X__)*aoG1(:,i,X__) &
                   +   cs(Y__)*aoG1(:,i,Y__) &
                   +   cs(Z__)*aoG1(:,i,Z__)
          end if

          if (xce%funTyp == OQP_FUNTYP_MGGA) then
            tmp(:,i,2:4) = f_t(2)*aoG1(:,i,X__:Z__)
          end if

        end do

        call dsyr2k('u', 'n', numAOs, numPts, 1.0_fp, &
                             aoV, numAOs, &
                             tmp(:,:,1), numAOs, &
                     0.0_fp, focks(:,:,j,2), numAOs)

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          do k = 1, 3
              call dsyr2k('u', 'n', numAOs, numPts, 0.25_fp, &
                                   aoG1(:,:,k), numAOs, &
                                   tmp(:,:,k+1), numAOs, &
                           1.0_fp, focks(:,:,j,2), numAOs)
          end do
        end if

      end do

    end associate

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Add XC derivative contribution to the Kohn-Sham-like matrices
!> @param[inout] fa   Kohn-Sham matrices
!> @author Vladimir Mironov
 subroutine RUpdate(dat, xce, myThread)
    use mod_dft_gridint, only: xc_der1, xc_der2_contr
    class(xc_engine_t) :: xce
    class(xc_consumer_tde_t) :: dat
    integer :: myThread

    integer :: i, j, k
    real(kind=fp) :: c3(3)
    real(kind=fp) :: d_r(2), d_s(3), d_t(2)
    real(kind=fp) :: f_r(2), f_s(3), f_t(2)
    real(kind=fp) :: rhoab(2), tauab(2), Sigma(3)
    real(kind=fp), pointer :: focks(:,:,:,:)
    real(kind=fp), pointer :: tmp(:,:,:)

    call dat%resetOrbPointers(xce, focks, tmp, myThread)

    associate ( aoV    => xce%aoV &
              , aoG1   => xce%aoG1 &
              , rRho   => dat%rRho(:,:,:,myThread)  &
              , drRho  => dat%drRho(:,:,:,:,myThread)  &
              , rTau   => dat%rTau(:,:,:,myThread)  &
              , dRho   => xce%xclib%dRho  &
              , numAOs => xce%numAOs_p &  ! number of pruned AOs
              , numPts => xce%numPts &
              , nMtx   => dat%nMtx &
      )
      do j = 1, nMtx

        do i = 1, numPts
            rhoab = rRho(1,i,j)
            Sigma = 0
            tauab = 0
            if (xce%funTyp /= OQP_FUNTYP_LDA) Sigma = 2*sum(drRho(1:3,1,i,j)*dRho(1:3,i))
            if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rTau(1,i,j)

            call xc_der1(xce, .false., i, d_r, d_s, d_t)
            call xc_der2_contr(xce, .false., i, &
                    rhoab, Sigma, tauab, &
                    f_r, f_s, f_t)
            ! LDA
            tmp(:,i,1) = 0.5_fp*f_r(1)*aoV(:,i)
            if (xce%funTyp /= OQP_FUNTYP_LDA) then
              ! GGA
              c3 = (2*f_s(1)+f_s(3)) * dRho(1:3,i) &
                 + (2*d_s(1)+d_s(3)) * drRho(1:3,1,i,j)

              tmp(:,i,1) = tmp(:,i,1) &
                   +   c3(X__)*aoG1(:,i,X__) &
                   +   c3(Y__)*aoG1(:,i,Y__) &
                   +   c3(Z__)*aoG1(:,i,Z__)
            end if

            if (xce%funTyp == OQP_FUNTYP_MGGA) then
              ! mGGA
              tmp(:,i,2) = f_t(1)*aoG1(:,i,X__)
              tmp(:,i,3) = f_t(1)*aoG1(:,i,Y__)
              tmp(:,i,4) = f_t(1)*aoG1(:,i,Z__)
            end if
        end do

        call dsyr2k('U', 'N', numAOs, numPts, 1.0_fp, &
                             aoV, numAOs, &
                             tmp(:,:,1), numAOs, &
                     0.0_fp, focks(:,:,j,1), numAOs)

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          do k = 1, 3
            call dsyr2k('U', 'N', numAOs, numPts, 0.25_fp, &
                                 aoG1(:,:,k), numAOs, &
                                 tmp(:,:,k+1), numAOs, &
                         1.0_fp, focks(:,:,j,1), numAOs)
          end do
        end if

      end do

    end associate

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute derivative XC contribution to the TD-DFT KS-like matrices
!> @param[in]    basis     basis set
!> @param[in]    isVecs    .true. if orbitals are provided instead of density matrix
!> @param[in]    wf        density matrix/orbitals
!> @param[inout] fx        fock-like matrices
!> @param[inout] dx        densities
!> @param[in]    nMtx      number of density/Fock-like matrices
!> @param[in]    threshold tolerance
!> @param[in]    isGGA     .TRUE. if GGA/mGGA functional used
!> @param[in]    infos     OQP metadata
!> @author Vladimir Mironov
  subroutine utddft_fxc(basis, molGrid, isVecs, &
                  wfa, wfb, &
                  fxa, fxb, &
                  dxa, dxb, &
                  nMtx, threshold, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: triangular_to_full

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid

    type(basis_set) :: basis
    logical, intent(in) :: isVecs
    integer, intent(in) :: nMtx
    real(kind=fp), intent(in) :: wfa(:,:), wfb(:,:)
    real(kind=fp), intent(inout), target :: dxa(:,:,:), dxb(:,:,:)
    real(kind=fp), intent(inout) :: fxa(:,:,:), fxb(:,:,:)
    real(kind=fp), intent(in) :: threshold

    type(xc_consumer_tde_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: i, j, nbf

    real(kind=fp), allocatable, target :: d2a(:,:), d2b(:,:)

    nbf = ubound(wfa,1)

    ! Scale w.f. by B.F. norms
    allocate(d2a(nbf,nbf), d2b(nbf,nbf))
    if (isVecs) then
      do i = 1, nbf
        d2a(:,i) = wfa(:,i) * basis%bfnrm(:)
        d2b(:,i) = wfb(:,i) * basis%bfnrm(:)
      end do
    else
      do i = 1, nbf
        d2a(:,i) = wfa(:,i) &
                * basis%bfnrm(i) &
                * basis%bfnrm(:)
        d2b(:,i) = wfb(:,i) &
                * basis%bfnrm(i) &
                * basis%bfnrm(:)
      end do
    end if

    ! Scale densities by B.F. norms
    do j = 1, nMtx
      do i = 1, nbf
        dxa(:,i,j) = dxa(:,i,j) &
                  * basis%bfnrm(i) &
                  * basis%bfnrm(:)
        dxb(:,i,j) = dxb(:,i,j) &
                  * basis%bfnrm(i) &
                  * basis%bfnrm(:)
      end do
    end do

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%functional => infos%functional
    xc_opts%hasBeta = .true.
    xc_opts%isWFVecs = isVecs
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = basis%mxam
    xc_opts%nDer = 0
    xc_opts%nXCDer = 2
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => d2a
    xc_opts%wfBeta  => d2b
    xc_opts%dft_threshold = threshold
    xc_opts%molGrid => molGrid

    dat%da => dxa
    dat%db => dxb
    dat%nMtx = nMtx

    call run_xc(xc_opts, dat, basis)

    deallocate (d2a, d2b)

    do j = 1, nMtx
      call triangular_to_full(dat%focks(:,:,j,1,1), nbf, 'u')
      do i = 1, nbf
        fxa(:,i,j) = fxa(:,i,j) &
                    +   dat%focks(:,i,j,1,1) &
                      * basis%bfnrm(i) &
                      * basis%bfnrm(:)
      end do
    end do

    do j = 1, nMtx
      call triangular_to_full(dat%focks(:,:,j,2,1), nbf, 'u')
      do i = 1, nbf
        fxb(:,i,j) = fxb(:,i,j) &
                    +   dat%focks(:,i,j,2,1) &
                      * basis%bfnrm(i) &
                      * basis%bfnrm(:)
      end do
    end do

    ! Scale densities back
    do j = 1, nMtx
      do i = 1, nbf
        dxa(:,i,j) = dxa(:,i,j) &
                  / basis%bfnrm(i) &
                  / basis%bfnrm(:)
        dxb(:,i,j) = dxb(:,i,j) &
                  / basis%bfnrm(i) &
                  / basis%bfnrm(:)
      end do
    end do

    call dat%clean()

  end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute derivative XC contribution to the TD-DFT KS-like matrices
!> @param[in]    basis     basis set
!> @param[in]    isVecs    .true. if orbitals are provided instead of density matrix
!> @param[in]    wf        density matrix/orbitals
!> @param[inout] fx        fock-like matrices
!> @param[inout] dx        densities
!> @param[in]    nMtx      number of density/Fock-like matrices
!> @param[in]    threshold tolerance
!> @param[in]    infos     OQP metadata
!> @author Vladimir Mironov
  subroutine tddft_fxc(basis, molGrid, isVecs, wf, fx, dx, &
                       nMtx, threshold, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: triangular_to_full

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid

    type(basis_set) :: basis
    logical, intent(in) :: isVecs
    integer, intent(in) :: nMtx
    real(kind=fp), intent(in) :: wf(:,:)
    real(kind=fp), intent(inout), target :: dx(:,:,:)
    real(kind=fp), intent(inout) :: fx(:,:,:)
    real(kind=fp), intent(in) :: threshold

    type(xc_consumer_tde_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: i, j, nbf

    real(kind=fp), allocatable, target :: d2(:,:)

    nbf = ubound(wf,1)

    ! Scale w.f. by B.F. norms
    allocate(d2(nbf,nbf))
    if (isVecs) then
      do i = 1, nbf
        d2(:,i) = wf(:,i) * basis%bfnrm(:)
      end do
    else
      do i = 1, nbf
        d2(:,i) = wf(:,i) &
                * basis%bfnrm(i) &
                * basis%bfnrm(:)
      end do
    end if

    ! Scale densities by B.F. norms
    do j = 1, nMtx
      do i = 1, nbf
        dx(:,i,j) = dx(:,i,j) &
                  * basis%bfnrm(i) &
                  * basis%bfnrm(:)
      end do
    end do

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%functional => infos%functional
    xc_opts%hasBeta = .false.
    xc_opts%isWFVecs = isVecs
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = basis%mxam
    xc_opts%nDer = 0
    xc_opts%nXCDer = 2
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => d2
    xc_opts%dft_threshold = threshold
    xc_opts%molGrid => molGrid

    dat%da => dx
    dat%nMtx = nMtx

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_xc(xc_opts, dat, basis)

    deallocate (d2)

    do j = 1, nMtx
      call triangular_to_full(dat%focks(:,:,j,1,1), nbf, 'u')
      do i = 1, nbf
        fx(:,i,j) = fx(:,i,j) &
                    +   dat%focks(:,i,j,1,1) &
                      * basis%bfnrm(i) &
                      * basis%bfnrm(:)
      end do
    end do

    ! Scale densities back
    do j = 1, nMtx
      do i = 1, nbf
        dx(:,i,j) = dx(:,i,j) &
                  / basis%bfnrm(i) &
                  / basis%bfnrm(:)
      end do
    end do

    call dat%clean()

  end subroutine

!-------------------------------------------------------------------------------

end module mod_dft_gridint_fxc
