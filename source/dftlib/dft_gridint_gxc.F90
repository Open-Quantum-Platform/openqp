module mod_dft_gridint_gxc

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t
  use mod_dft_gridint, only: X__, Y__, Z__
  use mod_dft_gridint, only: OQP_FUNTYP_LDA, OQP_FUNTYP_GGA, OQP_FUNTYP_MGGA
  use mod_dft_gridint_fxc, only: xc_consumer_tde_t
  use oqp_linalg

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_tde_t) :: xc_consumer_gxc_t
  contains
    procedure :: RUpdate => GxcRUpdate
    procedure :: UUpdate => GxcUUpdate
  end type

!-------------------------------------------------------------------------------

  private
  public tddft_gxc

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

!> @brief Compute XC 3rd derivative contractions needed for TD-DFT gradient, namely:
!>  \sum_{mns,kls'} (G_xc)_{mns,kls',pqs''} * (X+Y)_{mns} * (X+Y)_{kls'}
!> @details spin-polarized version
!> @param[inout] fa   Kohn-Sham matrices
!> @author Vladimir Mironov
 subroutine GxcUUpdate(dat, xce, mythread)
    use mod_dft_gridint, only: xc_der2_contr, xc_der3_contr
    class(xc_engine_t) :: xce
    class(xc_consumer_gxc_t) :: dat
    integer :: mythread

    integer :: i, j, k
    real(kind=fp) :: c3(3)
    real(kind=fp) :: f_s(3)
    real(kind=fp) :: g_r(2), g_s(3), g_t(2)
    real(kind=fp) :: rhoab(2), tauab(2), sigma(3), ssigma(3), dsaa, dsab, dsba, dsbb
    real(kind=fp), pointer :: focks(:,:,:,:)
    real(kind=fp), pointer :: tmp(:,:,:)

    call dat%resetOrbPointers(xce, focks, tmp, myThread)

    associate ( aoV    => xce%aoV &
              , aoG1   => xce%aoG1 &
              , rrho   => dat%rrho(:,:,:,mythread)  &
              , drrho  => dat%drrho(:,:,:,:,mythread)  &
              , rtau   => dat%rtau(:,:,:,mythread)  &
              , drho   => xce%xclib%drho  &
              , numAOs => xce%numAOs_p &  ! number of pruned AOs
              , numPts => xce%numPts &
              , nMtx   => dat%nMtx &
      )
      do j = 1, nMtx

        do i = 1, numPts
            rhoab = rrho(1:2,i,j)
            sigma = 0
            ssigma = 0
            tauab = 0
            if (xce%funTyp /= OQP_FUNTYP_LDA) then
              dsaa = dot_product(drrho(:,1,i,j), drho(1:3,i))
              dsab = dot_product(drrho(:,1,i,j), drho(4:6,i))
              dsbb = dot_product(drrho(:,2,i,j), drho(4:6,i))
              dsba = dot_product(drrho(:,2,i,j), drho(1:3,i))
              sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]

              dsaa = dot_product(drrho(:,1,i,j), drrho(:,1,i,j))
              dsab = dot_product(drrho(:,1,i,j), drrho(:,2,i,j))
              dsbb = dot_product(drrho(:,2,i,j), drrho(:,2,i,j))
              dsba = dot_product(drrho(:,2,i,j), drrho(:,1,i,j))
              sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]
            end if
            if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rtau(1:2,i,j)

            call xc_der3_contr(xce, i, &
                    rhoab, sigma, tauab, &
                    ssigma, &
                    f_s, &
                    g_r, g_s, g_t)

            ! LDA
            tmp(:,i,1) = 0.5_fp*g_r(1)*aoV(:,i)
            if (xce%funTyp /= OQP_FUNTYP_LDA) then
              ! GGA
              c3 =    2*g_s(1) * drho(1:3,i) &
                 +      g_s(3) * drho(4:6,i) &
                 +  2*2*f_s(1) * drrho(1:3,1,i,j) &
                 +    2*f_s(3) * drrho(1:3,2,i,j)

              tmp(:,i,1) = tmp(:,i,1) &
                   +   c3(X__)*aoG1(:,i,X__) &
                   +   c3(Y__)*aoG1(:,i,Y__) &
                   +   c3(Z__)*aoG1(:,i,Z__)
            end if

            if (xce%funTyp == OQP_FUNTYP_MGGA) then
              ! mGGA
              tmp(:,i,2:4) = g_t(1)*aoG1(:,i,X__:Z__)
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

        do i = 1, numPts
            rhoab = rrho(1:2,i,j)
            sigma = 0
            ssigma = 0
            tauab = 0
            if (xce%funTyp /= OQP_FUNTYP_LDA) then
              dsaa = dot_product(drrho(:,1,i,j), drho(1:3,i))
              dsab = dot_product(drrho(:,1,i,j), drho(4:6,i))
              dsbb = dot_product(drrho(:,2,i,j), drho(4:6,i))
              dsba = dot_product(drrho(:,2,i,j), drho(1:3,i))
              sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]

              dsaa = dot_product(drrho(:,1,i,j), drrho(:,1,i,j))
              dsab = dot_product(drrho(:,1,i,j), drrho(:,2,i,j))
              dsbb = dot_product(drrho(:,2,i,j), drrho(:,2,i,j))
              dsba = dot_product(drrho(:,2,i,j), drrho(:,1,i,j))
              sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]
            end if
            if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rtau(1:2,i,j)

            call xc_der3_contr(xce, i, &
                    rhoab, sigma, tauab, &
                    ssigma, &
                    f_s, &
                    g_r, g_s, g_t)

            ! LDA
            tmp(:,i,1) = 0.5_fp*g_r(2)*aoV(:,i)
            if (xce%funTyp /= OQP_FUNTYP_LDA) then
              ! GGA
              c3 =    2*g_s(2) * drho(4:6,i) &
                 +      g_s(3) * drho(1:3,i) &
                 +  2*2*f_s(2) * drrho(1:3,2,i,j) &
                 +    2*f_s(3) * drrho(1:3,1,i,j)

              tmp(:,i,1) = tmp(:,i,1) &
                   +   c3(X__)*aoG1(:,i,X__) &
                   +   c3(Y__)*aoG1(:,i,Y__) &
                   +   c3(Z__)*aoG1(:,i,Z__)
            end if

            if (xce%funTyp == OQP_FUNTYP_MGGA) then
              ! mGGA
              tmp(:,i,2:4) = g_t(2)*aoG1(:,i,X__:Z__)
            end if
        end do

        call dsyr2k('U', 'N', numAOs, numPts, 1.0_fp, &
                             aoV, numAOs, &
                             tmp(:,:,1), numAOs, &
                     0.0_fp, focks(:,:,j,2), numAOs)

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          do k = 1, 3
            call dsyr2k('U', 'N', numAOs, numPts, 0.25_fp, &
                                 aoG1(:,:,k), numAOs, &
                                 tmp(:,:,k+1), numAOs, &
                         1.0_fp, focks(:,:,j,2), numAOs)
          end do
        end if

      end do

    end associate

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute XC 3rd derivative contractions needed for TD-DFT gradient, namely:
!>  \sum_{mns,kls'} (G_xc)_{mns,kls',pqs''} * (X+Y)_{mns} * (X+Y)_{kls'}
!> @details not spin-polarized version
!> @param[inout] fa   Kohn-Sham matrices
!> @author Vladimir Mironov
 subroutine GxcRUpdate(dat, xce, mythread)
    use mod_dft_gridint, only: xc_der2_contr, xc_der3_contr
    class(xc_engine_t) :: xce
    class(xc_consumer_gxc_t) :: dat
    integer :: mythread

    integer :: i, j, k
    real(kind=fp) :: c3(3)
    real(kind=fp) :: f_s(3)
    real(kind=fp) :: g_r(2), g_s(3), g_t(2)
    real(kind=fp) :: rhoab(2), tauab(2), sigma(3), ssigma(3)
    real(kind=fp), pointer :: focks(:,:,:,:)
    real(kind=fp), pointer :: tmp(:,:,:)

    call dat%resetOrbPointers(xce, focks, tmp, myThread)

    associate ( aoV    => xce%aoV &
              , aoG1   => xce%aoG1 &
              , rrho   => dat%rrho(:,:,:,mythread)  &
              , drrho  => dat%drrho(:,:,:,:,mythread)  &
              , rtau   => dat%rtau(:,:,:,mythread)  &
              , drho   => xce%xclib%drho  &
              , numAOs => xce%numAOs_p &  ! number of pruned AOs
              , numPts => xce%numPts &
              , nMtx   => dat%nMtx &
      )
      do j = 1, nMtx

        do i = 1, numPts
            rhoab = rrho(1,i,j)
            sigma = 0
            ssigma = 0
            tauab = 0
            if (xce%funTyp /= OQP_FUNTYP_LDA) then
              sigma = 2*sum(drrho(1:3,1,i,j)*drho(1:3,i))
              ssigma = 2*sum(drrho(1:3,1,i,j)*drrho(1:3,1,i,j))
            end if
            if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rtau(1,i,j)

            call xc_der3_contr(xce, i, &
                    rhoab, sigma, tauab, &
                    ssigma, &
                    f_s, &
                    g_r, g_s, g_t)

            ! LDA
            tmp(:,i,1) = 0.5_fp*g_r(1)*aoV(:,i)
            if (xce%funTyp /= OQP_FUNTYP_LDA) then
              ! GGA
              c3 =   (2*g_s(1)+g_s(3)) * drho(1:3,i) &
                 + 2*(2*f_s(1)+f_s(3)) * drrho(1:3,1,i,j)

              tmp(:,i,1) = tmp(:,i,1) &
                   +   c3(X__)*aoG1(:,i,X__) &
                   +   c3(Y__)*aoG1(:,i,Y__) &
                   +   c3(Z__)*aoG1(:,i,Z__)
            end if

            if (xce%funTyp == OQP_FUNTYP_MGGA) then
              ! mGGA
              tmp(:,i,2) = g_t(1)*aoG1(:,i,X__)
              tmp(:,i,3) = g_t(1)*aoG1(:,i,Y__)
              tmp(:,i,4) = g_t(1)*aoG1(:,i,Z__)
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
!> @param[in]    infos     OQP metadata
!> @author Vladimir Mironov
  subroutine tddft_gxc(basis, molGrid, isVecs, wf, fx, dx, &
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

    type(xc_consumer_gxc_t) :: dat
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
    xc_opts%nXCDer = 3
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => d2
    xc_opts%dft_threshold = threshold
    xc_opts%molGrid => molGrid

    dat%da => dx
    dat%nMtx = nMtx

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

end module mod_dft_gridint_gxc
