module mod_dft_fuzzycell

  use precision, only: fp

  use basis_tools, only: basis_set
  use mod_dft_partfunc, only: partition_function
  use mod_dft_molgrid, only: dft_grid_t

  implicit none

  real(KIND=fp), parameter :: HUGEFP = huge(1.0_fp)
  real(KIND=fp), parameter :: PI = 3.141592653589793238463_fp
  real(KIND=fp), parameter :: FOUR_PI = 4.0_fp*PI

  private
  public prune_basis

  public dft_fc_blk

contains

!-------------------------------------------------------------------------------

!> @brief Find shells and primitives which are significant in a given set
!>  of 3D coordinates
!  TODO: move to basis set related source file
!> @author Vladimir Mironov
  subroutine prune_basis(inBas, xyzv, xyzat, nSh, nPrim, nBf, &
                         outSh, outShNG, outPrim, atoms)
    use constants, only: num_cart_bf
    use atomic_structure_m, only: atomic_structure
    type(basis_set), intent(IN) :: inBas
    real(KIND=fp), intent(IN) :: xyzv(:, :), xyzat(:)
    integer, intent(OUT) :: nSh, nPrim, nBf
    integer, contiguous, intent(OUT) :: outSh(:), outShNG(:), outPrim(:)
    type(atomic_structure), intent(in) :: atoms

!    integer :: nat, ich, mul, num, nqmt, ne, na, nb, ian
!    real(KIND=fp) :: zan, c

    integer :: ish, ig, nCur

    nSh = 0
    nPrim = 0
    nBf = 0
    do ish = 1, inBas%nshell
      associate (ncontr => inBas%ncontr(ish), &
                 g0 => inBas%g_offset(ish), &
                 am => inBas%am(ish), &
                 xyz => atoms%xyz(1:3,inBas%origin(ish))-xyzat(1:3))

        nCur = 0
        do ig = g0, g0+ncontr-1
          if (bfnz(xyz, xyzv, ig)) then
            nPrim = nPrim+1
            nCur = nCur+1
            outPrim(nPrim) = ig
          end if
        end do
        if (nCur > 0) then
          nSh = nSh+1
          nBf = nBf+num_cart_bf(am)
          outSh(nSh) = ish
          outShNG(nSh) = nCur
        end if
      end associate
    end do

  contains

    logical function bfnz(xyz, xyzv, iPrim)
      real(KIND=fp), intent(IN) :: xyz(:), xyzv(:, :)
      integer, intent(IN) :: iPrim
      integer :: i
      bfnz = .false.
      do i = 1, ubound(xyzv, 1)
        bfnz = sum((xyz(1:3)-xyzv(i, 1:3))**2) < inBas%prim_mx_dist2(iPrim)
        if (bfnz) exit
      end do
    end function

  end subroutine

!-------------------------------------------------------------------------------

!> @brief Assemble numerical atomic DFT grids to a molecular grid
!> @param[in]    atmxvec  array of atomic X coordinates
!> @param[in]    atmyvec  array of atomic Y coordinates
!> @param[in]    atmzvec  array of atomic Z coordinates
!> @param[in]    rij      interatomic distances
!> @param[in]    nat      number of atoms
!> @param[in]    curAt    index of current atom
!> @param[in]    rad      effective (e.g. Bragg-Slater) radius of current atom
!> @param[inout] wtab     normalized cell function values for LRD
!> @param[in]    aij      surface shifting factors for Becke's method
!> @author Vladimir Mironov
  subroutine dft_fc_blk(molGrid, dft_partfun, atmxyz, at_mx_dist2, rij, nat, wtab, aij)

!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set

    type(dft_grid_t), intent(inout) :: molGrid
    integer, intent(IN) :: dft_partfun
    integer, intent(IN) :: nat
    real(KIND=fp), intent(IN) :: rij(nat, nat), atmxyz(:,:), at_mx_dist2(:)
    real(KIND=fp), allocatable, intent(INOUT) :: wtab(:,:,:)
    real(KIND=fp), intent(IN), contiguous, optional :: aij(:,:)

    real(KIND=fp), allocatable :: rijInv(:, :), ri(:), wtintr(:)

    integer :: me, nProc
    integer :: iChunk, iSlice
    integer :: maxSize, maxPts, nMolPts

    type(partition_function) :: partfunc

!   Set up selected partition function
    call partfunc%set(dft_partfun)

!   Initialize temporary space for atomic cell functions
!   and point-to-atom distances
    allocate (wtintr(nat), ri(nat))

!   Compute inverse interatomic distances
    allocate (rijInv(nat, nat))
    where (rij /= 0.0_fp)
      rijInv = 1.0_fp/rij
    elsewhere
      rijInv = 0.0_fp
    end where

!   Clear point counters
    maxSize = 0
    maxPts = 0
    nMolPts = 0

    nProc = 1
    me = 0

!   Compute total weights for grid points
!   MPI - static LB, OpenMP - dynamic LB
!$omp parallel do &
!$omp   private(iSlice, iChunk, ri, wtintr) &
!$omp   reduction(max:maxPts, maxSize) &
!$omp   reduction(+:nMolPts) &
!$omp   schedule(dynamic)
    do iChunk = 1, molGrid%nSlices, nProc
      iSlice = iChunk+me

      if (iSlice <= molGrid%nSlices) then
!         Apply selected algorithm on the current slice
        call do_bfcSlice(molGrid, iSlice, partfunc, nat, &
                         atmxyz, at_mx_dist2, &
                         ri, rij, rijInv, wtintr, wtab, aij)

!         Count points:
        nMolPts = nMolPts+molGrid%nTotPts(iSlice)
        maxPts = max(maxPts, molGrid%nTotPts(iSlice))
        maxSize = max(maxSize, molGrid%nAngPts(iSlice) &
                      *molGrid%nRadPts(iSlice))
      end if

    end do
!$omp end parallel do

!   Finalize grid info
    molGrid%maxNRadTimesNAng = maxSize
    molGrid%maxSlicePts = maxPts
    molGrid%nMolPts = nMolPts

  end subroutine

!> @brief Compute total weights for points in a slice
!> @param[in]    iSlice    index of current slice
!> @param[in]    partfunc  partition function
!> @param[in]    nat       number of atoms
!> @param[in]    atmxvec   array of atomic X coordinates
!> @param[in]    atmyvec   array of atomic Y coordinates
!> @param[in]    atmzvec   array of atomic Z coordinates
!> @param[inout] ri        tmp array to store point to atoms distances
!> @param[in]    rij       interatomic distances
!> @param[in]    rijInv    inverse interatomic distances
!> @param[inout] wtintr    tmp array to store cell function values
!> @param[inout] wtab      normalized cell function values for LRD
!> @param[in]    aij       surface shifting factors for Becke's method
!> @author Vladimir Mironov
  subroutine do_bfcSlice(molGrid, iSlice, partfunc, nAt, &
                         atmxyz, at_mx_dist2, &
                         ri, rij, rijInv, wtintr, wtab, aij)
    use mod_grid_storage, only: grid_3d_t

    type(dft_grid_t), intent(inout) :: molGrid
    integer, intent(IN) :: iSlice, nAt
    type(partition_function), intent(IN) :: partfunc
    real(KIND=fp), contiguous :: atmxyz(:,:), at_mx_dist2(:), &
                                 ri(:), rij(:, :), rijInv(:, :), wtintr(:)
    real(KIND=fp), allocatable :: wtab(:,:,:)
    real(KIND=fp), intent(IN), contiguous, optional :: aij(:, :)

    type(grid_3d_t), pointer :: curGrid
    integer :: iAng, iRad, iPt, iAtm
    real(KIND=fp) :: radwt, r1, ptxyz(3), wtAngRad
    real(KIND=fp) :: wtnrm

    logical :: lrd_flag

    lrd_flag = allocated(wtab)

    associate ( &
      dummyAtom => molGrid%dummyAtom, &
      rInner => molGrid%rInner, &
      iAngStart => molGrid%iAngStart(iSlice), &
      iRadStart => molGrid%iRadStart(iSlice), &
      nAngPts => molGrid%nAngPts(iSlice), &
      nRadPts => molGrid%nRadPts(iSlice), &
      wtStart => molGrid%wtStart(iSlice)-1, &
      totWts => molGrid%totWts, &
      ntp => molGrid%nTotPts(iSlice), &
      isInner => molGrid%isInner(iSlice), &
      curAt => molGrid%idOrigin(iSlice), &
      rad => molGrid%rAtm(iSlice))

      ntp = 0

      isInner = 0
      if (molGrid%rad_pts(iRadStart+nRadPts-1)*rad < rInner(curAt)) then
!           Weights of the whole slice are unchanged
        isInner = 1
        ntp = nAngPts*nRadPts
!           Nothing left to do for inner slice
        return
      end if

      curGrid => molGrid%spherical_grids%getbyid(molGrid%idAng(iSlice))
      associate ( &
        xAng => curGrid%x(iAngStart:iAngStart+nAngPts-1), &
        yAng => curGrid%y(iAngStart:iAngStart+nAngPts-1), &
        zAng => curGrid%z(iAngStart:iAngStart+nAngPts-1), &
        wAng => curGrid%w(iAngStart:iAngStart+nAngPts-1))

        do iAng = 1, nAngPts
          rloop: do iRad = 1, nRadPts

            r1 = rad*molGrid%rad_pts(iRadStart+iRad-1)
            radWt = rad*rad*rad*molGrid%rad_wts(iRadStart+iRad-1)
            wtAngRad = FOUR_PI*radWt*wAng(iAng)

            iPt = (iAng-1)*nRadPts+iRad

!               Quick check for inner quadrature points
            if (r1 < rInner(curAt)) then
!                   Point weight is unchanged
              totWts(wtStart+iPt, curAt) = wtAngRad
              ntp = ntp+1
!                   Next point
              cycle rloop
            end if

            ptxyz(1) = r1*xAng(iAng) + atmxyz(1,curAt)
            ptxyz(2) = r1*yAng(iAng) + atmxyz(2,curAt)
            ptxyz(3) = r1*zAng(iAng) + atmxyz(3,curAt)

            do iAtm = 1, nAt
              if (dummyAtom(iAtm)) cycle
              ri(iAtm) = norm2(ptxyz-atmxyz(:,iAtm))

!                   Check if the point belongs to any other atom
              if (iAtm /= curAt .and. &
                  (r1-ri(iAtm) > rij(iAtm, curAt)* &
                   0.5d0*(1.0d0+partfunc%limit))) then
!                       This and all next points along the same direction
!                       belong to another atom. Switch to the next angular
!                       point. It is never happen when Becke's function is
!                       selected.
                exit rloop
              end if
            end do

!           Pre-screen small density.
!           Other points may be significant. Skip to next
!           iteration of radial loop.
!           Note, weight is set non-zero to discriminate between
!           points, belonging to the space of another atom.
            if (ALL(ri*ri>at_mx_dist2 .or. dummyAtom)) then
                totWts(wtStart+iPt,curAt) = tiny(1.0_fp)
                cycle rloop
            end if


!               If screening fails, follow regular BFC procedure and check
!               all atom pairs
            call do_bfc(totWts(wtStart+iPt, curAt), wtnrm, wtintr, wtAngRad, &
                          ri, rijInv, nAt, curAt, dummyAtom, partfunc, aij)

            if (lrd_flag) wtab(1:nat, curAt, iPt) = wtintr(1:nat)*wtnrm

!               Same as before, if weight is zero, the current and
!               all next points along the same direction belong
!               to another atom. Switch to the next angular point.
!               It is never happen when Becke's function is selected.
            if (totWts(wtStart+iPt, curAt) == 0.0d0) then
                exit rloop
            end if

            ntp = ntp+1

          end do rloop
        end do
      end associate

    end associate

  end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute total weight of a grid point
!> @param[out]   wt        total weight
!> @param[out]   wtNorm    sum of cell function values
!> @param[out]   cells     array of cell function values
!> @param[in]    wtAngRad  unmodified weight
!> @param[inout] ri        tmp array to store point to atoms distances
!> @param[in]    dummyAtom array to indicate which atoms are "dummy"
!> @param[in]    rijInv    inverse interatomic distances
!> @param[in]    numAt     number of atoms
!> @param[in]    curAtom   current atom
!> @param[in]    partfunc  partition function
!> @param[in]    aij       surface shifting factors for Becke's method
!> @author Vladimir Mironov
  subroutine do_bfc(wt, wtNorm, cells, wtAngRad, &
                    ri, rijInv, numAt, curAtom, dummyAtom, partfunc, aij)
    logical, contiguous, intent(IN) :: dummyAtom(:)
    real(KIND=fp), contiguous, intent(IN) :: ri(:), rijInv(:, :)
    integer, intent(IN) :: numAt, curAtom
    real(KIND=fp), intent(IN) :: wtAngRad
    type(partition_function), intent(IN) :: partfunc
    real(KIND=fp), intent(OUT) :: wt, wtNorm
    real(KIND=fp), contiguous, intent(OUT) :: cells(:)
    real(KIND=fp), contiguous, intent(IN), optional :: aij(:, :)

    integer :: i, j
    real(KIND=fp) :: f, mu

    where (dummyAtom)
      cells = 0.0d0
    elsewhere
      cells = 1.0d0
    end where

    do i = 2, numAt
      if (dummyAtom(i)) cycle
      do j = 1, i-1
        if (dummyAtom(j)) cycle
        mu = (ri(i)-ri(j))*rijInv(j, i)
        if (present(aij)) mu = mu+aij(j, i)*(1.0d0-mu*mu)
        f = partfunc%eval(mu)
        cells(i) = cells(i)*abs(f)
        cells(j) = cells(j)*abs(1.0d0-f)
      end do
    end do

    wtNorm = 1.0d0/sum(cells(1:numAt))

    wt = cells(curAtom)*wtNorm*wtAngRad

  end subroutine

!-------------------------------------------------------------------------------
end module mod_dft_fuzzycell
