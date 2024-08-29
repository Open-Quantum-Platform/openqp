module mod_dft_molgrid

  use precision, only: fp, dp
  use mod_grid_storage, only: list_grid_3d_t, grid_3d_t

  use lebedev, only: lebedev_get_grid
  use bragg_slater_radii, only: BRSL_NUM_ELEMENTS

  implicit none

  integer, parameter :: DEFAULT_NSLICES_PER_ATOM = 100
  integer, parameter :: MAXGRID = 10
  real(kind=fp), parameter :: GROWTH_FACTOR = 1.5
  integer, parameter :: ALLOC_PAD = 4

  integer, parameter :: MAXDEPTH = 2
  real(KIND=fp), parameter :: HUGEFP = huge(1.0_fp)

  type, extends(list_grid_3d_t) :: sorted_grid_t
!< Array to store vertices for triangulation of sphere, for internal use
    real(KIND=fp), allocatable :: triangles(:, :, :, :)
  contains
    procedure :: add_grid => get_sorted_lebedev_pts
    procedure :: init => initSortedListGrid
  end type


!> @brief Type to store molecular grid information
  type :: dft_grid_t
    !< current number of slices
    integer :: nSlices = 0
    !< maximum number of slices
    integer :: maxSlices = 0
    !< maximum number of grid points per atom
    integer :: maxAtomPts = 0
    !< maximum number of nonzero grid points per slice
    integer :: maxSlicePts = 0
    !< maximum possible number of grid points per slice
    integer :: maxNRadTimesNAng = 0
    !< total number of nonzero grid points
    integer :: nMolPts = 0

    !< spherical atomic grids used in this molecular grid
    type(sorted_grid_t) :: spherical_grids

    ! Every slice is a spacially localized subgrid of molecular grid:
    ! slice = (few radial points)x(few angular points)
    ! The following arrays contains data for each slice

    !< index of angular grid
    integer, allocatable :: idAng(:)
    !< index of the first point in angular grid array
    integer, allocatable :: iAngStart(:)
    !< number of angular grid points
    integer, allocatable :: nAngPts(:)
    !< index of the first point in radial grid array
    integer, allocatable :: iRadStart(:)
    !< number of radial grid points
    integer, allocatable :: nRadPts(:)
    !< number of grid points with nonzero weight
    integer, allocatable :: nTotPts(:)
    !< index of atom which owns the slice
    integer, allocatable :: idOrigin(:)
    !< type of chunk -- used for different processing
    integer, allocatable :: chunkType(:)
    !< index of the first point in weights array
    integer, allocatable :: wtStart(:)
    !< isInner == 1 means that the slice is 'inner' -- i.e. its weight is not modified and is equal to wtRad*wtAng
    integer, allocatable :: isInner(:)

    ! The following arrays contains information about atoms
    !< effective radii of atoms
    real(KIND=fp), allocatable :: rAtm(:)
    !< max radii for 'inner' points
    real(KIND=fp), allocatable :: rInner(:)
    !< .TRUE. for dummy atoms
    logical, allocatable :: dummyAtom(:)

    ! The following arrays contains information about the radial grid
    !< Radial grid points
    real(KIND=fp), allocatable :: rad_pts(:)
    !< Radial grid weights
    real(KIND=fp), allocatable :: rad_wts(:)

    !< array to store grid point weights
    real(KIND=fp), allocatable :: totWts(:, :)
    !< number of nonzero grid points per atom
    integer, allocatable, private :: wt_top(:)

  contains
    procedure, pass :: getSliceData => getSliceData
    procedure, pass :: getSliceNonZero => getSliceNonZero
    procedure, pass :: exportGrid => exportGrid
    procedure, pass :: setSlice => setSlice
    procedure, pass :: reset => reset_dft_grid_t
    procedure, pass :: compress => compress_dft_grid_t
    procedure, pass :: extend => extend_dft_grid_t
    procedure, pass :: add_atomic_grid => add_atomic_grid
    procedure, pass :: add_slices => add_slices
    procedure, pass :: find_neighbours => find_neighbours
  end type dft_grid_t

  private
  public dft_grid_t

contains

!*******************************************************************************
! Legacy Fortran wrappers
!-------------------------------------------------------------------------------

!> @brief Initialize grid list, possibly specifying its size
!> @param[in]  iSize   desired list size
!> @author Vladimir Mironov
  subroutine initSortedListGrid(list, iSize)
    class(sorted_grid_t) :: list
    integer, optional, intent(IN) :: iSize
    integer :: ntrngoct, ntrng, i

    ntrngoct = 4**MAXDEPTH
    ntrng = 8*ntrngoct

    if (.not.allocated(list%triangles)) then
      allocate (list%triangles(3, 3, ntrng, MAXDEPTH+1))
      do i = 1, MAXDEPTH+1
        call triangulate_sphere(list%triangles(:, :, :, i), i-1)
      end do
    end if

    call list%list_grid_3d_t%init(iSize)

  end subroutine

!> @brief Clear grid list
!> @author Vladimir Mironov
  subroutine clearSortedListGrid(list)
    class(sorted_grid_t) :: list
    if (allocated(list%triangles)) deallocate(list%triangles)
    call list%list_grid_3d_t%clear()
  end subroutine

!> @brief Find nearest neghbours of every atom and the corresponding distance
!>  The results are stored into `grid` dataset
!> @details This is needed for screening inner points in SSF variant of
!>  Becke's fuzzy cell method
!> @param[in]    rij       array of interatomic distances
!> @param[in]    nAt       number of atoms
!> @author Vladimir Mironov
  subroutine find_neighbours(grid, rij, partFunType)

    use mod_dft_partfunc, only: partition_function

    class(dft_grid_t) :: grid
    integer, intent(IN) :: partFunType
    real(KIND=fp), intent(IN) :: rij(:,:)
    type(partition_function) :: partFunc
    integer :: nAt, i, j
    real(KIND=fp) :: distNN

    call partFunc%set(partFunType)

    nAt = ubound(rij,1)

    do j = 1, nAt
!     Find distance to the nearest neghbour of current atom
      distNN = HUGEFP
      do i = 1, nAt
        if (grid%dummyAtom(i) .or. i == j) cycle
        distNN = min(distNN, rij(i, j))
      end do
      grid%rInner(j) = 0.5*distNN*(1.0-partFunc%limit)
    end do

  end subroutine

!> @brief Get spherical Lebedev grid of nReq number of points
!> @details If the grid is not yet computed, then compute it,
!>  sort according to pre-defined triangulation scheme and save
!>  the result for later use. If the same grid was already computed,
!>  fetch it from the in-memory storage
!> @author Vladimir Mironov
  subroutine get_sorted_lebedev_pts(grids, nreq)
    class(sorted_grid_t), intent(inout) :: grids
    integer, intent(IN) :: nreq
    real(KIND=fp), allocatable :: xyz(:,:), w(:)
    integer, allocatable :: icnt(:, :)
    type(grid_3d_t), pointer :: curGrid

    allocate (icnt(8*4**MAXDEPTH, MAXDEPTH+2))

    curGrid => grids%get(nreq)
    if (associated(curGrid)) return

    allocate(xyz(nreq,3), w(nreq), source=0.0d0)
    call lebedev_get_grid(nreq, xyz, w)
    call tesselate_layer(xyz(:,1), xyz(:,2), xyz(:,3), w, nreq, icnt, grids%triangles)
    curGrid => grids%set( &
               grid_3d_t(x=xyz(:,1), &
                         y=xyz(:,2), &
                         z=xyz(:,3), &
                         w=w, &
                         npts=nreq, &
                         izones=int(icnt,2)) &
    )

  end subroutine

!> @brief Append slice data of an atom to `grid` arrays
!> @param[in]   idAtom      index of current atom
!> @param[in]   layers      meta-data for atom grid layers: radial grid, angular grid, and angular grid separation depth
!> @param[in]   nLayers     number of atom grid layers
!> @param[in]   rAtm        effective radius of current atom
!> @author Vladimir Mironov
  subroutine add_slices(grid, idAtm, layers, nLayers, rAtm)
    class(dft_grid_t) :: grid
    integer, intent(IN) :: idAtm, layers(:, :), nLayers
    real(KIND=fp), intent(IN) :: rAtm
    integer :: i, j, ilen, idrad0, idrad1, idang0, idang1, nAngTr, nPtSlc, depth
    type(grid_3d_t), pointer :: curAng

    idrad1 = 0

    do i = 1, nLayers
      associate ( &
        nCur => grid%nSlices, &
        maxSlices => grid%maxSlices, &
        nextRad => layers(1, i), &
        depth0 => layers(3, i), &
        wtTopAtm => grid%wt_top(idAtm))

        curAng => grid%spherical_grids%get(layers(2, i))

        idrad0 = idrad1+1
        idrad1 = nextRad
        depth = depth0
        if (grid%rad_pts(idrad0) < 2.0 .and. depth == -1) then
          if (curAng%npts*(idrad1-idrad0+1) >= 10*32) then
            depth = 1
          else if (curAng%npts*(idrad1-idrad0+1) >= 10*8) then
            depth = 0
          end if
        end if

        select case (depth)
        case (-1)

          nPtSlc = (idrad1-idrad0+1)*curAng%nPts

          if (nPtSlc == 0) cycle
          if (nCur == maxSlices) call grid%extend

          nCur = nCur+1

          call grid%setSlice(nCur, &
                                idrad0, idrad1-idrad0+1, &
                                curAng%idGrid, 1, curAng%nPts, &
                                wtTopAtm+1, &
                                idAtm, rAtm, 0)

          wtTopAtm = wtTopAtm+nPtSlc

        case (0:3)
          ilen = 4**depth
          idang1 = 0
          do j = 1, 8*ilen
            nAngTr = curAng%izones(j, depth+2)
            if (nAngTr == 0) cycle
            if (nCur == maxSlices) call grid%extend

            nCur = nCur+1

            idang0 = idang1+1
            idang1 = idang0+nAngTr-1

            nPtSlc = nAngTr*(idrad1-idrad0+1)

            call grid%setSlice(nCur, &
                                  idrad0, idrad1-idrad0+1, &
                                  curAng%idGrid, idAng0, nAngTr, &
                                  wtTopAtm+1, &
                                  idAtm, rAtm, 1)

            wtTopAtm = wtTopAtm+nPtSlc
          end do

        case DEFAULT
          nPtSlc = (idRad1-idRad0+1)*curAng%nPts
          if (nPtSlc == 0) cycle
          if (nCur == maxSlices) call grid%extend

          nCur = nCur+1

          call grid%setSlice(nCur, &
                                idrad0, idrad1-idrad0+1, &
                                curAng%idGrid, 1, curAng%nPts, &
                                wtTopAtm+1, &
                                idAtm, rAtm, 2)

          wtTopAtm = wtTopAtm+nPtSlc

        end select

      end associate
    end do

  end subroutine

!> @brief Split atomic 3D grid on spacially localized clusters ("slices")
!>  and appended the data to `grid` dataset
!> @details First, decompose atomic grid on layers, depending on the pruning
!>  scheme and the distance from the nuclei. It is done in this procedure.
!>  Then, split layers on slices and append their data to the `grid`
!>  dataset. It is done by calling subroutine `addSlices`.
!> @param[in]   idAtm       index of current atom
!> @param[in]   nGrids      number of grids in pruning scheme
!> @param[in]   rAtm        effective radius of current atom
!> @param[in]   pruneRads   radii of the pruning scheme
!> @author Vladimir Mironov
  subroutine add_atomic_grid(grid, atomic_grid)
    use mod_grid_storage, only: atomic_grid_t
    class(dft_grid_t) :: grid
    type(atomic_grid_t) :: atomic_grid
    real(KIND=fp), parameter :: EPS = tiny(0.0_fp)

    integer, parameter :: depths(5) = [-1, 0, 0, 1, -2]
    real(KIND=fp), parameter :: stepSz(5) = [1.0, 1.5, 4.0, 8.0, 1.0e9]
    real(KIND=fp) :: limits(5)

    integer :: i, nLayers, nGrids, iSpl, iRMin, iRMax, iRNext, nRad!, k
    integer :: idAtm
    real(KIND=fp) :: rAtm, rMin, rMax, rNext

    integer :: layers(3, MAXGRID*5)

    idAtm = atomic_grid%idAtm
    nRad = ubound(atomic_grid%rad_pts, 1)
    nGrids = ubound(atomic_grid%sph_npts, 1)
    rAtm = atomic_grid%rAtm

    limits = [1.0, 2.0, 8.0, 12.0, 20.0]
    if (grid%rInner(idAtm) /= 0.0_fp) then
      limits(1) = grid%rInner(idAtm)+EPS
    end if

    nLayers = 0
    iRMax = 0
    iRNext = 0
    do i = 1, nGrids
      iRMin = iRMax+1
      iRMax = findl(atomic_grid%sph_radii(i), atomic_grid%rad_pts, iRMin)

      if (nGrids == 1) iRMax = nRad

      rMin = rAtm*atomic_grid%rad_pts(iRMin)
      rMax = rAtm*atomic_grid%rad_pts(iRMax)

      iSpl = findl(rAtm*rMin, limits, 1)

      rNext = rMin

      do
        rNext = min(rMax, rNext+stepSz(iSpl))
        iRNext = min(iRMax, findl(rNext/rAtm, atomic_grid%rad_pts, iRNext+1))

        nLayers = nLayers+1
        layers(1, nLayers) = iRNext
        layers(2, nLayers) = atomic_grid%sph_npts(i)
        layers(3, nLayers) = depths(iSpl)

        if (rNext >= limits(iSpl)) iSpl = iSpl+1
        if (iRNext == iRMax) exit
      end do

    end do

    call grid%add_slices(idAtm, layers, nLayers, rAtm)

  contains

!> @brief Find the location of the largest element in
!>  sorted (in ascending order) real-valued array,
!>  which is smaller than the specified value. Return last array index
!>  if `ALL(X>=V)`
!> @note Linear search is used due to problem specifics:
!>  tiny arrays and sequential access
!> @param[in]   x       value to found
!> @param[in]   v       sorted real array
!> @param[in]   hint    the index from which the lookup is started
!> @author Vladimir Mironov
    function findl(x, v, hint) result(id)
      integer :: id
      integer :: hint
      real(KIND=fp) :: x, v(:)
      do id = max(lbound(v, 1), hint), ubound(v, 1)
        if (x < v(id)) return
      end do
      id = ubound(v, 1)
    end function

  end subroutine

!> @brief Split spherical grid on localized clusters of points
!> @param[inout] xs     grid point X coordinates
!> @param[inout] ys     grid point Y coordinates
!> @param[inout] zs     grid point Z coordinates
!> @param[inout] wt     grid point weights
!> @param[in]    npts   number of points
!> @param[out]   icnts  sizes of clusters
!> @author Vladimir Mironov
  subroutine tesselate_layer(xs, ys, zs, wt, npts, icnts, triangles)
    real(KIND=fp), intent(INOUT) :: xs(:), ys(:), zs(:)
    real(KIND=fp), intent(INOUT) :: wt(:)
    real(KIND=fp), target, intent(in) :: triangles(:, :, :, :)
    integer, intent(IN) :: npts
    integer, intent(OUT) :: icnts(:, -1:)
    real(KIND=fp), allocatable :: tmp(:, :, :)
    integer :: i, j, itrng, itrngoct, ntrng, ntrngoct, m1, m2, m3, i0, i1, ilen
    real(KIND=fp), pointer :: trng(:, :, :)

    ntrngoct = 4**MAXDEPTH
    ntrng = 8*ntrngoct

    trng => triangles(:, :, :, MAXDEPTH+1)

    allocate (tmp(size(xs, dim=1), 4, ntrng))

    icnts(1:ntrng, :) = 0

    associate (icnt => icnts(:, MAXDEPTH))
      do i = 1, npts
        associate (x => xs(i), y => ys(i), z => zs(i), w => wt(i))
!        Find the octant of this point
          m1 = int(sign(0.5_fp, x)+0.5)
          m2 = int(sign(0.5_fp, y)+0.5)
          m3 = int(sign(0.5_fp, z)+0.5)
          itrngoct = 8-(m1+2*m2+4*m3)
          do itrng = ntrngoct*(itrngoct-1)+1, ntrngoct*itrngoct
!           If point is on the boundary between two cells - first
!           one gets it
            if (is_vector_inside_shape([x, y, z], trng(:, :, itrng))) exit
          end do
          icnt(itrng) = icnt(itrng)+1
          tmp(icnt(itrng), 1, itrng) = x
          tmp(icnt(itrng), 2, itrng) = y
          tmp(icnt(itrng), 3, itrng) = z
          tmp(icnt(itrng), 4, itrng) = w
        end associate
      end do

      i0 = 0
      do i = 1, ntrng
        i1 = i0+icnt(i)
        xs(i0+1:i1) = tmp(1:icnt(i), 1, i)
        ys(i0+1:i1) = tmp(1:icnt(i), 2, i)
        zs(i0+1:i1) = tmp(1:icnt(i), 3, i)
        wt(i0+1:i1) = tmp(1:icnt(i), 4, i)
        i0 = i1
      end do

    end associate

    do i = MAXDEPTH-1, 0, -1
      ilen = 4**i
      do j = 1, 8*ilen
        icnts(j, i) = sum(icnts((j-1)*4+1:j*4, i+1))
      end do
    end do

    icnts(1, -1) = npts

    deallocate (tmp)

  end subroutine

!> @brief Compute triangulation map of the sphere
!> @details This subroutine returns inscribed polyhedron with triangular
!>  faces and octahedral symmetry. Polyhedron faces are computed by recursive
!>  subdivision of the octahedron faces on four spherical
!>  equilateral triangles.
!> @param[in]   depth   depth of recursion
!> @param[out]  vect    vertices of triangles
!> @author Vladimir Mironov
  subroutine triangulate_sphere(vect, depth)
    integer, intent(IN) :: depth
    real(KIND=fp), intent(OUT) :: vect(3, 3, *)
    real(KIND=fp), dimension(3) :: v1, v2, v3
    real(KIND=fp), parameter :: xf(8) = [1, -1, 1, -1, 1, -1, 1, -1]
    real(KIND=fp), parameter :: yf(8) = [1, 1, -1, -1, 1, 1, -1, -1]
    real(KIND=fp), parameter :: zf(8) = [1, 1, 1, 1, -1, -1, -1, -1]
    integer :: ip
    integer :: i, ip0, ip1

    ip = 0

!   Fill in first octant
    v1 = [1.0, 0.0, 0.0]
    v2 = [0.0, 1.0, 0.0]
    v3 = [0.0, 0.0, 1.0]
    call subdivide_triangle(v1, v2, v3, vect, ip, depth)

!   Fill in other octants - use octaherdal symmetry of the polyhedron
    do i = 2, 8
      ip0 = (i-1)*ip+1
      ip1 = i*ip
      vect(1, 1:3, ip0:ip1) = xf(i)*vect(1, 1:3, 1:ip)
      vect(2, 1:3, ip0:ip1) = yf(i)*vect(2, 1:3, 1:ip)
      vect(3, 1:3, ip0:ip1) = zf(i)*vect(3, 1:3, 1:ip)
    end do

  end subroutine

!> @brief Recursively subdivide spherical triangle on four equal triangles
!> @param[in]       v1      vertex of initial triangle
!> @param[in]       v2      vertex of initial triangle
!> @param[in]       v3      vertex of initial triangle
!> @param[inout]    res     stack of triangles; dim: (xyz,vertices,triangles)
!> @param[in]       ip      index of the last triangle on stack
!> @param[in]       depth   current depth of recursion
!> @author Vladimir Mironov
  recursive subroutine subdivide_triangle(v1, v2, v3, res, ip, depth)
    integer, intent(IN) :: depth
    integer, intent(INOUT) :: ip
    real(KIND=8), intent(IN) :: v1(3), v2(3), v3(3)
    real(KIND=8), intent(INOUT) :: res(3, 3, *)
    real(KIND=8) :: c1(3), c2(3), c3(3)

    if (depth < 0 .or. depth > MAXDEPTH) then
      write (*, *) 'DEPTH=', depth, ' IS .GT. MAXDEPTH=', MAXDEPTH
    end if

    if (depth == 0) then
!     add lowest level triangle to the result array
      ip = ip+1
      res(:, 1, ip) = v1
      res(:, 2, ip) = v2
      res(:, 3, ip) = v3
      return
    end if

!   find centers of arcs between input vectors
    c1 = v1+v2
    c2 = v2+v3
    c3 = v1+v3
    c1 = c1/sqrt(sum(c1*c1))
    c2 = c2/sqrt(sum(c2*c2))
    c3 = c3/sqrt(sum(c3*c3))

!   subdivide the resulting four triangles
    call subdivide_triangle(v1, c1, c3, res, ip, depth-1)
    call subdivide_triangle(c1, v2, c2, res, ip, depth-1)
    call subdivide_triangle(c3, c2, v3, res, ip, depth-1)
    call subdivide_triangle(c1, c2, c3, res, ip, depth-1)

  end subroutine

!> @brief Clean up `grid` dataset and reallocate arrays if needed
!> @param[in]   maxSlices   guess for maximum number of slices
!> @param[in]   nAt         number of atoms in a system
!> @param[in]   maxPtPerAt  maximum number of points per atom
!> @param[in]   nRad        size of the radial grid
!> @author Vladimir Mironov
  subroutine reset_dft_grid_t(grid, nAt, maxPtPerAt, nRad)
    class(dft_grid_t), intent(INOUT) :: grid
    integer, intent(IN) :: nAt, maxPtPerAt, nRad
    integer :: maxSlices

    maxSlices = DEFAULT_NSLICES_PER_ATOM * nAt

    if (grid%maxSlices < maxSlices) then
      grid%maxSlices = maxSlices
      if (allocated(grid%idAng)) then
        deallocate (grid%idAng, grid%iAngStart, grid%nAngPts, &
                    grid%iRadStart, grid%nRadPts, grid%nTotPts, &
                    grid%idOrigin, grid%chunkType, grid%rAtm, &
                    grid%wtStart, grid%isInner)
      end if

      allocate (grid%idAng(maxSlices), source=0)
      allocate (grid%iAngStart(maxSlices), source=0)
      allocate (grid%nAngPts(maxSlices), source=0)

      allocate (grid%iRadStart(maxSlices), source=0)
      allocate (grid%nRadPts(maxSlices), source=0)

      allocate (grid%nTotPts(maxSlices), source=0)

      allocate (grid%idOrigin(maxSlices), source=0)
      allocate (grid%chunkType(maxSlices), source=0)
      allocate (grid%rAtm(maxSlices), source=0.0_fp)
      allocate (grid%wtStart(maxSlices), source=0)
      allocate (grid%isInner(maxSlices), source=0)

    end if

!   Initialize spherical grids storage
    call grid%spherical_grids%init()

!   Initialize radial grids storage
    if (allocated(grid%rad_pts)) deallocate(grid%rad_pts)
    if (allocated(grid%rad_wts)) deallocate(grid%rad_wts)
    allocate(grid%rad_pts(nRad), source=0.0_fp)
    allocate(grid%rad_wts(nRad), source=0.0_fp)

    if (allocated(grid%totWts)) deallocate (grid%totWts)
    if (allocated(grid%wt_top)) deallocate (grid%wt_top)
    allocate (grid%totWts(maxPtPerAt, nAt), source=0.0_fp)
    allocate (grid%wt_top(nAt), source=0)

    if (allocated(grid%rInner)) deallocate (grid%rInner)
    if (allocated(grid%dummyAtom)) deallocate (grid%dummyAtom)
    allocate (grid%rInner(nAt), source=0.0_fp)
    allocate (grid%dummyAtom(nAt), source=.false.)

    grid%maxAtomPts = maxPtPerAt
    grid%maxSlicePts = 0
    grid%maxNRadTimesNAng = 0
    grid%nSlices = 0
    grid%nMolPts = 0

  end subroutine reset_dft_grid_t

!> @brief Compress the grid for eaach atom
!> @note This is a very basic implementation and it will be changed in future
!> @author Vladimir Mironov
  subroutine compress_dft_grid_t(grid)
    class(dft_grid_t), intent(INOUT) :: grid

    integer :: i
    integer :: iCur, nPts
    integer, allocatable :: iCurWt(:)

    allocate (iCurWt(ubound(grid%totWts, 2)))

    iCur = 0
    iCurWt = 1
    do i = 1, grid%nSlices
      if (grid%nTotPts(i) > 0) then
        iCur = iCur+1
        grid%idAng(iCur) = grid%idAng(i)
        grid%iAngStart(iCur) = grid%iAngStart(i)
        grid%nAngPts(iCur) = grid%nAngPts(i)
        grid%iRadStart(iCur) = grid%iRadStart(i)
        grid%nRadPts(iCur) = grid%nRadPts(i)
        grid%idOrigin(iCur) = grid%idOrigin(i)
        grid%chunkType(iCur) = grid%chunkType(i)
        grid%rAtm(iCur) = grid%rAtm(i)
        grid%isInner(iCur) = grid%isInner(i)

        grid%nTotPts(iCur) = grid%nTotPts(i)

        associate (totWts => grid%totWts, &
                   iAt => grid%idOrigin(i), wt0 => grid%wtStart(i), &
                   nAng => grid%nAngPts(i), nRad => grid%nRadPts(i))

          nPts = nAng*nRad

          totWts(iCurWt(iAt):iCurWt(iAt)+nPts-1, iAt) = &
            totWts(wt0:wt0+nPts-1, iAt)

          grid%wtStart(iCur) = iCurWt(iAt)
          iCurWt(iAt) = iCurWt(iAt)+nPts
        end associate
      end if
    end do

    grid%nSlices = iCur

    deallocate (iCurWt)

  end subroutine compress_dft_grid_t

!> @brief Extend arrays in molecular grid type
!> @author Vladimir Mironov
  subroutine extend_dft_grid_t(grid)
    class(dft_grid_t), intent(INOUT) :: grid
    integer :: maxSlices

    maxSlices = int( (grid%maxSlices+1)*GROWTH_FACTOR )
    maxSlices = ( (maxSlices - 1) / ALLOC_PAD + 1 ) * ALLOC_PAD

    call reallocate_int(grid%idAng, maxSlices)
    call reallocate_int(grid%iAngStart, maxSlices)
    call reallocate_int(grid%nAngPts, maxSlices)
    call reallocate_int(grid%iRadStart, maxSlices)
    call reallocate_int(grid%nRadPts, maxSlices)
    call reallocate_int(grid%nTotPts, maxSlices)
    call reallocate_int(grid%idOrigin, maxSlices)
    call reallocate_int(grid%chunkType, maxSlices)
    call reallocate_int(grid%wtStart, maxSlices)
    call reallocate_int(grid%isInner, maxSlices)

    call reallocate_real(grid%rAtm, maxSlices)

    grid%maxSlices = maxSlices

  contains
!> @brief Extend allocatable array of integers to a new size
!> @param[inout]  v         the array
!> @param[in]     newsz     new size of the array
!> @author Vladimir Mironov
    subroutine reallocate_int(v, newsz)
      integer, allocatable, intent(INOUT) :: v(:)
      integer, intent(IN) :: newsz
      integer, allocatable :: nv(:)
      if (allocated(v)) then
        allocate (nv(1:newsz), source=v)
        call move_alloc(from=nv, to=v)
      end if
    end subroutine reallocate_int

!> @brief Extend allocatable array of integers to a new size
!> @param[inout]  v         the array
!> @param[in]     newsz     new size of the array
!> @author Vladimir Mironov
    subroutine reallocate_real(v, newsz)
      real(kind=fp), allocatable, intent(INOUT) :: v(:)
      integer, intent(IN) :: newsz
      real(kind=fp), allocatable :: nv(:)
      if (allocated(v)) then
        allocate (nv(1:newsz), source=v)
        call move_alloc(from=nv, to=v)
      end if
    end subroutine reallocate_real

  end subroutine extend_dft_grid_t

!> @brief Get coordinates and weight for quadrature points, which
!>   belong to a slice
!> @param[in]     iSlice    index of a slice
!> @param[out]    xyzw      coordinates and weights
!> @author Vladimir Mironov
  subroutine getSliceData(grid, iSlice, xyzw)
    class(dft_grid_t), intent(IN) :: grid
    integer, intent(IN) :: iSlice
    real(KIND=fp), contiguous, intent(OUT) :: xyzw(:, :)

    real(KIND=fp), parameter :: FOUR_PI = 4.0_fp*3.141592653589793238463_fp

    type(grid_3d_t), pointer :: curGrid
    integer :: iAng, iRad, iPt
    real(KIND=fp) :: r1

    curGrid => grid%spherical_grids%getbyid(grid%idAng(iSlice))
    associate ( &
      rad => grid%rAtm(iSlice), &
      isInner => grid%isInner(iSlice), &
      wtStart => grid%wtStart(iSlice)-1, &
      curAt => grid%idOrigin(iSlice), &
      iAngStart => grid%iAngStart(iSlice), &
      iRadStart => grid%iRadStart(iSlice), &
      nAngPts => grid%nAngPts(iSlice), &
      nRadPts => grid%nRadPts(iSlice))

      associate ( &
        xAng => curGrid%x(iAngStart:iAngStart+nAngPts-1), &
        yAng => curGrid%y(iAngStart:iAngStart+nAngPts-1), &
        zAng => curGrid%z(iAngStart:iAngStart+nAngPts-1), &
        wAng => curGrid%w(iAngStart:iAngStart+nAngPts-1))

        do iAng = 1, nAngPts
        do iRad = 1, nRadPts
          r1 = rad*grid%rad_pts(iRadStart+iRad-1)
          iPt = (iAng-1)*nRadPts+iRad
          xyzw(iPt, 1) = r1*xAng(iAng)
          xyzw(iPt, 2) = r1*yAng(iAng)
          xyzw(iPt, 3) = r1*zAng(iAng)
          if (isInner == 0) then
            xyzw(iPt, 4) = grid%totWts(wtStart+iPt, curAt)
          else
            xyzw(iPt, 4) = FOUR_PI*rad*rad*rad* &
                           grid%rad_wts(iRadStart+iRad-1)*wAng(iAng)
          end if
        end do
        end do

      end associate
    end associate

  end subroutine

!> @brief Export grid for use in legacy TD-DFT code
!> @param[out]   xyz       coordinates of grid points
!> @param[out]   w         weights of grid points
!> @param[out]   kcp       atoms, which the points belongs to
!> @param[out]   npts      number of nonzero point
!> @param[in]    cutoff    cutoff to skip small weights
!> @author Vladimir Mironov
  subroutine exportGrid(grid, xyz, w, kcp, c, npts, cutoff)
    class(dft_grid_t), intent(IN) :: grid
    real(KIND=fp), intent(OUT) :: xyz(3, *), w(*)
    real(KIND=fp), intent(in) :: c(3, *)
    integer, intent(OUT) :: kcp(*)
    integer, intent(OUT) :: npts
    real(KIND=fp), intent(IN) :: cutoff

    integer :: iSlice, iPt, iPtOld, nPt, iAt

    real(KIND=fp), allocatable :: xyzw(:, :)

    iPt = 0

!$omp parallel private(xyzw, iSlice, nPt, iPtOld, iAt)

    allocate (xyzw(grid%maxSlicePts, 4))

!$omp do schedule(dynamic)
    do iSlice = 1, grid%nSlices

      call grid%getSliceNonZero(cutoff, iSlice, xyzw, nPt)

      if (nPt == 0) cycle

      !$omp atomic capture
      iPtOld = iPt
      iPt = iPt+nPt
      !$omp end atomic

      iAt = grid%idOrigin(iSlice)

      xyz(1, iPtOld+1:iPtOld+nPt) = c(1, iAt)+xyzw(1:nPt, 1)
      xyz(2, iPtOld+1:iPtOld+nPt) = c(2, iAt)+xyzw(1:nPt, 2)
      xyz(3, iPtOld+1:iPtOld+nPt) = c(3, iAt)+xyzw(1:nPt, 3)
      w(iPtOld+1:iPtOld+nPt) = xyzw(1:nPt, 4)
      kcp(iPtOld+1:iPtOld+nPt) = iAt

    end do
!$omp end do
!$omp end parallel

    npts = iPt

  end subroutine

!> @brief Get grid points from a slice, which weights are larger
!>  than a cutoff.
!> @param[in]    cutoff    cutoff to skip small weights
!> @param[in]    iSlice    index of current slice
!> @param[out]   xyzw      coordinates and weights of grid points
!> @param[out]   nPt       number of nonzero point for slice
!> @author Vladimir Mironov
  subroutine getSliceNonZero(grid, cutoff, iSlice, xyzw, nPt)
    use constants, only: pi
    class(dft_grid_t), intent(IN) :: grid
    real(KIND=fp), intent(IN) :: cutoff
    integer, intent(IN) :: iSlice
    real(KIND=fp), contiguous, intent(OUT) :: xyzw(:, :)
    integer, intent(OUT) :: nPt

    real(KIND=fp), parameter :: FOUR_PI = 4.0_fp*pi

    type(grid_3d_t), pointer :: curGrid
    integer :: iAng, iRad, iPt
    real(KIND=fp) :: r1, wtCur

    nPt = 0
    curGrid => grid%spherical_grids%getbyid(grid%idAng(iSlice))

    associate ( &
      rad => grid%rAtm(iSlice), &
      isInner => grid%isInner(iSlice), &
      wtStart => grid%wtStart(iSlice)-1, &
      totWts => grid%totWts, &
      curAt => grid%idOrigin(iSlice), &
      iAngStart => grid%iAngStart(iSlice), &
      iRadStart => grid%iRadStart(iSlice), &
      nAngPts => grid%nAngPts(iSlice), &
      nRadPts => grid%nRadPts(iSlice))

      associate ( &
        xAng => curGrid%x(iAngStart:iAngStart+nAngPts-1), &
        yAng => curGrid%y(iAngStart:iAngStart+nAngPts-1), &
        zAng => curGrid%z(iAngStart:iAngStart+nAngPts-1), &
        wAng => curGrid%w(iAngStart:iAngStart+nAngPts-1))

        do iAng = 1, nAngPts
        do iRad = 1, nRadPts
          iPt = (iAng-1)*nRadPts+iRad
          if (isInner == 0) then
            wtCur = totWts(wtStart+iPt, curAt)
          else
            wtCur = FOUR_PI*rad*rad*rad* &
                    grid%rad_wts(iRadStart+iRad-1)*wAng(iAng)
          end if
          if (wtCur == 0.0_fp) exit
          if (wtCur < cutoff) cycle
          nPt = nPt+1
          r1 = rad*grid%rad_pts(iRadStart+iRad-1)
          xyzw(nPt, 1) = r1*xAng(iAng)
          xyzw(nPt, 2) = r1*yAng(iAng)
          xyzw(nPt, 3) = r1*zAng(iAng)
          xyzw(nPt, 4) = wtCur
        end do
        end do
      end associate
    end associate

  end subroutine

!> @brief Set slice data
!> @param[in]    iSlice     index of current slice
!> @param[in]    iRadStart  index of the first point in radial grid array
!> @param[in]    nRadPts    number of radial grid points
!> @param[in]    idAng      index of angular grid
!> @param[in]    idAngStart index of the first point in angular grid array
!> @param[in]    nAngPts    number of angular grid points
!> @param[in]    wtStart    index of the first point in weights array
!> @param[in]    idAtm      index of atom which owns the slice
!> @param[in]    rAtm       effective radius of current atom
!> @param[in]    chunkType  type of chunk -- used for different processing
!> @author Vladimir Mironov
  subroutine setSlice(grid, iSlice, iRadStart, nRadPts, &
                      idAng, iAngStart, nAngPts, &
                      wtStart, idAtm, rAtm, chunkType)
    class(dft_grid_t), intent(INOUT) :: grid

    integer, intent(IN) :: iSlice
    integer, intent(IN) :: iRadStart, nRadPts
    integer, intent(IN) :: idAng, iAngStart, nAngPts
    integer, intent(IN) :: wtStart
    integer, intent(IN) :: idAtm, chunkType
    real(KIND=fp), intent(IN) :: rAtm

    grid%idAng(iSlice) = idAng

    grid%iRadStart(iSlice) = iRadStart
    grid%iAngStart(iSlice) = iAngStart

    grid%nRadPts(iSlice) = nRadPts
    grid%nAngPts(iSlice) = nAngPts
    grid%nTotPts(iSlice) = 0

    grid%idOrigin(iSlice) = idAtm
    grid%rAtm(iSlice) = rAtm

    grid%chunkType(iSlice) = chunkType
    grid%wtStart(iSlice) = wtStart

  end subroutine setSlice

!> @brief Compute average vector of two 3D vectors
!> @param[in]    vectors    input 3D vectors
!> @author Vladimir Mironov
  pure function vector_average(vectors) result(avg)
    real(KIND=fp), intent(IN) :: vectors(:, :)
    real(KIND=fp) :: avg(3)

    real(KIND=fp) :: norm

    avg = sum(vectors(1:3, :), dim=2)

    norm = 1.0/sqrt(sum(avg**2))
    avg = avg*norm

  end function vector_average

!> @brief Compute cross product of two 3D vectors
!> @param[in]    a    first vector
!> @param[in]    b    second vector
!> @author Vladimir Mironov
  pure function cross_product(a, b) result(c)
    real(KIND=fp), intent(IN) :: a(:), b(:)
    real(KIND=fp) :: c(3)
    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)
  end function

!> @brief Check if test vector crosses spherical polygon
!> @param[in]    test  test vector
!> @param[in]    vecs  vectors, defining the polygon
!> @author Vladimir Mironov
  pure function is_vector_inside_shape(test, vecs) result(inside)
    real(KIND=fp), intent(IN) :: test(:), vecs(:, :)
    logical :: inside

    integer :: i
    real(KIND=fp) :: sgn, sgnold, vectmp(3)

!     Test vector must at least has same direction, as average vector
!     vectmp(:) = sum(vecs(:,:),dim=1)
!     if (sum(test*vectmp)<0.0) then
!         inside = .false.
!         return
!     end if

    inside = .true.
    vectmp = cross_product(vecs(:, ubound(vecs, 2)), vecs(:, lbound(vecs, 2)))
    sgnold = sign(1.0_fp, sum(test(1:3)*vectmp(1:3)))

    do i = lbound(vecs, 2), ubound(vecs, 2)-1
      vectmp = cross_product(vecs(:, i), vecs(:, i+1))
      sgn = sign(1.0_fp, sum(test(1:3)*vectmp(1:3)))
      if (sgn * sgnold < 0.0d0) then ! different signs
        inside = .false.
        return
      end if
      sgnold = sgn
    end do

  end function

end module mod_dft_molgrid
