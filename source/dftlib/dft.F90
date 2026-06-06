module dft
! A Module for grid based DFT
  use messages,  only: show_message, WITH_ABORT
  use precision, only: dp
  use io_constants, only: iw
  use basis_tools, only: basis_set
  use mod_dft_molgrid, only: dft_grid_t

  implicit none

  character(len=*), parameter :: module_name = "dft"

  private
  public dft_initialize
  public dftclean
  public dftexcor
  public dftder

!> @brief Pruned-grid specification
!> @details A pruned grid is defined per atom type by up to `ngrids`
!>   radial regions: region i of type t covers radii (in units of the
!>   atomic radius) up to `radii(i,t)` and uses a `nang(i,t)`-point
!>   Lebedev sphere.  `rad_id` maps each atom to its type.
!>   Alternatively (SG-2/SG-3), regions are given as counts of
!>   consecutive radial shells: when `nradPerRegion(i,t) > 0`, region i
!>   of type t spans the next `nradPerRegion(i,t)` shells of the radial
!>   grid and `radii` is ignored for that type (regions with a zero
!>   count are unused).
!>   Several radial grids may coexist: `radial_id` maps each atom to
!>   one of `nrad_types` radial grids.  Radial type 1 is always the
!>   standard unit-radius grid (scaled by the Bragg-Slater radius);
!>   types >= 2 are element-specific grids in absolute bohr: DE2
!>   (`de2_alpha`/`de2_rmax` give alpha and the outermost node) or,
!>   when `me_rscale` is allocated and positive, MultiExp with
!>   `rad_npts` nodes and scaling radius `me_rscale` (SG-0).
!>   If `nang_override` is allocated and non-zero for an atom type,
!>   that type is unpruned: a single `nang_override(t)`-point Lebedev
!>   sphere is used at ALL radii (heavy-atom fallback).
  type dft_grid_pruned_t
    integer :: nrad = 0
    integer :: ngrids = 1
    integer, allocatable :: nang(:,:)    !< (region, atom type)
    real(kind=dp), allocatable :: radii(:,:) !< (region, atom type)
    integer, allocatable :: rad_id(:)    !< atom -> atom type
    integer, allocatable :: nradPerRegion(:,:) !< (region, atom type); 0 = unused
    integer, allocatable :: nang_override(:) !< per atom type; 0 = no override
    integer :: nrad_types = 1            !< number of radial grids
    integer, allocatable :: radial_id(:) !< atom -> radial grid type
    real(kind=dp), allocatable :: de2_alpha(:) !< DE2 alpha of radial type
    real(kind=dp), allocatable :: de2_rmax(:)  !< DE2 outermost node, bohr
    integer, allocatable :: rad_npts(:)  !< nodes of radial type (0: global nrad)
    real(kind=dp), allocatable :: me_rscale(:) !< MultiExp R of radial type (0: DE2)
  end type

!  SG1 region boundaries (in units of the atomic radius) and Lebedev
!  orders, from P.M.W. Gill, B.G. Johnson, J.A. Pople,
!  Chem. Phys. Lett. 209 (1993) 506: rows are H-He, Li-Ne, Na-Ar.
!  SG1 is only defined up to Ar; heavier atoms (row 4) fall back to
!  the unpruned 194-point grid at all radii.  Row 4 is fully
!  overridden via nang_override (set in dft_set_options), so its
!  boundaries are never used.
  real(kind=dp), parameter :: sg1rads(5,4) = reshape(&
        [0.2500d0, 0.500d0, 1.0d00, 4.50d0, 9999999.9d0,   &
         0.1667d0, 0.500d0, 0.90d0, 3.50d0, 9999999.9d0,   &
         0.1000d0, 0.400d0, 0.80d0, 2.5d0,  9999999.9d0,   &
         9999999.9d0, 9999999.9d0, 9999999.9d0, 9999999.9d0, 9999999.9d0], &
         shape(sg1rads))
  integer, parameter :: sg1atoms(4) =  [2,  10,  18,  137]
  integer, parameter :: sg1grids(5) =  [6, 38, 86, 194, 86]

!  SG-2 / SG-3 pruned grids: S. Dasgupta, J.M. Herbert,
!  J. Comput. Chem. 38, 869 (2017).  Radial grid: Mitani
!  double-exponential (DE2), M. Mitani, Theor. Chem. Acc. 130,
!  645 (2011), with element-specific alpha and Nr = 75 (SG-2) or
!  Nr = 99 (SG-3).  The first/last radial nodes are pinned to
!  r = 1e-7 bohr and the element-specific R_max below (values from
!  NVIDIA cuEST's SG-2/SG-3 implementation, CUDALibrarySamples,
!  cuest_molecular_grid.py; R_max is 10x the EML scaling radius).
!  The pruning sectors are counts of consecutive radial shells
!  (ascending radius), each integrated on the given Lebedev sphere.
!  Defined for Z in {1, 3-9, 11-17}; other elements fall back to the
!  unpruned 302/590-point grid on the standard radial grid.
  integer, parameter :: SG_NELEM = 15
  integer, parameter :: sg_elem_z(SG_NELEM) = &
        [1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17]
  real(kind=dp), parameter :: SG_DE2_RMIN = 1.0d-7
  real(kind=dp), parameter :: sg_de2_rmax(SG_NELEM) = &
        [15.0d0, 38.7d0, 26.5d0, 22.0d0, 17.1d0, 14.1d0, 12.3d0, &
         10.8d0, 42.1d0, 32.5d0, 34.3d0, 27.5d0, 23.2d0, 20.6d0, &
         18.4d0]

  integer, parameter :: SG2_NRAD = 75
  integer, parameter :: SG2_MAXSEC = 5
  real(kind=dp), parameter :: sg2_alpha(SG_NELEM) = &
        [2.6d0, 3.2d0, 2.4d0, 2.4d0, 2.2d0, 2.2d0, 2.2d0, 2.2d0, &
         3.2d0, 2.4d0, 2.5d0, 2.3d0, 2.5d0, 2.5d0, 2.5d0]
!  number of radial shells per sector
  integer, parameter :: sg2_cnt(SG2_MAXSEC, SG_NELEM) = reshape([ &
        35, 12, 16, 7, 5,   & ! H
        35, 12, 17, 7, 4,   & ! Li
        35, 12, 17, 7, 4,   & ! Be
        35, 12, 17, 7, 4,   & ! B
        35, 12, 17, 7, 4,   & ! C
        35, 12, 17, 7, 4,   & ! N
        30, 14, 18, 8, 5,   & ! O
        26, 16, 19, 8, 6,   & ! F
        35, 12, 17, 7, 4,   & ! Na
        35, 12, 17, 7, 4,   & ! Mg
        32, 15, 17, 7, 4,   & ! Al
        32, 15, 17, 7, 4,   & ! Si
        30, 14, 17, 7, 7,   & ! P
        30, 14, 17, 7, 7,   & ! S
        26, 16, 19, 8, 6],  & ! Cl
        shape(sg2_cnt))
!  Lebedev order of each sector
  integer, parameter :: sg2_leb(SG2_MAXSEC, SG_NELEM) = reshape([ &
        6, 110, 302,  86, 26,   & ! H
        6, 110, 302,  86, 50,   & ! Li
        6, 110, 302,  86, 50,   & ! Be
        6, 110, 302, 146, 26,   & ! B
        6, 110, 302, 146, 26,   & ! C
        6, 110, 302,  86, 26,   & ! N
        6, 110, 302, 146, 50,   & ! O
        6, 110, 302, 110, 50,   & ! F
        6, 110, 302,  86, 50,   & ! Na
        6, 110, 302,  86, 50,   & ! Mg
        6, 110, 302, 146, 86,   & ! Al
        6, 110, 302, 146, 50,   & ! Si
        6, 110, 302, 146, 38,   & ! P
        6, 110, 302, 146, 38,   & ! S
        6, 110, 302, 110, 50],  & ! Cl
        shape(sg2_leb))

  integer, parameter :: SG3_NRAD = 99
  integer, parameter :: SG3_MAXSEC = 9
  real(kind=dp), parameter :: sg3_alpha(SG_NELEM) = &
        [2.7d0, 3.0d0, 2.4d0, 2.4d0, 2.4d0, 2.4d0, 2.6d0, 2.1d0, &
         3.2d0, 2.6d0, 2.6d0, 2.8d0, 2.4d0, 2.4d0, 2.6d0]
  integer, parameter :: sg3_nsec(SG_NELEM) = &
        [5, 5, 7, 6, 7, 5, 9, 7, 5, 5, 7, 6, 8, 8, 7]
  integer, parameter :: sg3_cnt(SG3_MAXSEC, SG_NELEM) = reshape([ &
        45, 16, 21, 10,  7,  0,  0,  0,  0,   & ! H
        46, 16, 22,  9,  6,  0,  0,  0,  0,   & ! Li
        42,  6, 14, 22,  3,  6,  6,  0,  0,   & ! Be
        42,  6, 14, 22,  9,  6,  0,  0,  0,   & ! B
        46, 16, 22,  1,  2,  6,  6,  0,  0,   & ! C
        40, 18, 24, 11,  6,  0,  0,  0,  0,   & ! N
        40, 14,  2,  2, 24,  1,  1,  8,  7,   & ! O
        35, 17,  4, 25,  2,  8,  8,  0,  0,   & ! F
        46, 16, 22,  9,  6,  0,  0,  0,  0,   & ! Na
        48, 15, 20,  7,  9,  0,  0,  0,  0,   & ! Mg
        42,  6, 14, 22,  3,  6,  6,  0,  0,   & ! Al
        42,  6, 14, 22,  9,  6,  0,  0,  0,   & ! Si
        35,  1, 18,  4, 25,  2,  8,  6,  0,   & ! P
        35,  1, 18,  4, 25,  2,  8,  6,  0,   & ! S
        35, 17,  4, 25,  2,  8,  8,  0,  0],  & ! Cl
        shape(sg3_cnt))
  integer, parameter :: sg3_leb(SG3_MAXSEC, SG_NELEM) = reshape([ &
        6, 110, 590, 194,  50,   0,   0,   0,  0,   & ! H
        6, 110, 590, 146,  50,   0,   0,   0,  0,   & ! Li
        6,  86, 110, 590, 194, 146,  50,   0,  0,   & ! Be
        6,  86, 110, 590, 194,  50,   0,   0,  0,   & ! B
        6, 146, 590, 302, 194, 146,  86,   0,  0,   & ! C
        6, 110, 590, 146,  50,   0,   0,   0,  0,   & ! N
        6, 110, 194, 302, 590, 302, 194, 146, 50,   & ! O
        6, 110, 194, 590, 194, 110,  50,   0,  0,   & ! F
        6, 110, 590, 146,  50,   0,   0,   0,  0,   & ! Na
        6, 110, 590, 146,  50,   0,   0,   0,  0,   & ! Mg
        6,  86, 110, 590, 194, 146,  50,   0,  0,   & ! Al
        6,  86, 110, 590, 194,  50,   0,   0,  0,   & ! Si
        6,  86, 110, 194, 590, 194, 146,  50,  0,   & ! P
        6,  86, 110, 194, 590, 194, 146,  50,  0,   & ! S
        6, 110, 194, 590, 194, 110,  50,   0,  0],  & ! Cl
        shape(sg3_leb))

!  SG-0 pruned grid: S.-H. Chien, P.M.W. Gill, J. Comput. Chem. 27,
!  730 (2006), Table 1 (counts/orders as listed in Psi4's
!  cubature.cc), applied in ASCENDING radial order like the
!  SG-2/SG-3 sectors (validated here: H2O/NH3/CH4/thymine BHHLYP
!  energies agree with dense grids to ~1e-4, while the reversed
!  order fails by up to 1e-2).  Radial grid: MultiExp (Gauss
!  quadrature on (0,1) for the weight ln^2 x; P.M.W. Gill, S.-H.
!  Chien, J. Comput. Chem. 24, 732 (2003)) with Nr = 23 (Z = 1,
!  3-9) or 26 (Z = 11-17) and an element-specific scaling radius R:
!  r_i = -R ln(x_i), w_i = R^3 omega_i / x_i (incl. r^2 Jacobian).
!  Two spherical-rule substitutions are applied: the original
!  18-point rule (Abramowitz & Stegun, not a Lebedev grid) is not
!  available here and the 74-point Lebedev rule carries a negative
!  weight, which the weight-screening machinery here (positive
!  cutoffs in getSliceNonZero etc.) cannot represent; both are
!  replaced by the next safe Lebedev order, 26 and 86 respectively.
!  Defined for Z in {1, 3-9, 11-17}; other elements (He, Ne,
!  Z >= 18) fall back to the SG1 scheme on the standard radial grid.
  integer, parameter :: SG0_MAXSEC = 15
  integer, parameter :: sg0_nrad(SG_NELEM) = &
        [23, 23, 23, 23, 23, 23, 23, 23, 26, 26, 26, 26, 26, 26, 26]
  real(kind=dp), parameter :: sg0_rscale(SG_NELEM) = &
        [1.30d0, 1.95d0, 2.20d0, 1.45d0, 1.20d0, 1.10d0, 1.10d0, &
         1.20d0, 2.30d0, 2.20d0, 2.10d0, 1.30d0, 1.30d0, 1.10d0, &
         1.45d0]
  integer, parameter :: sg0_nsec(SG_NELEM) = &
        [11, 11, 12, 7, 13, 10, 11, 10, 8, 13, 15, 11, 11, 12, 12]
  integer, parameter :: sg0_cnt(SG0_MAXSEC, SG_NELEM) = reshape([ &
        6, 3, 1, 1, 1, 1, 6, 1, 1, 1, 1, 0, 0, 0, 0,   & ! H
        6, 3, 1, 1, 1, 1, 6, 1, 1, 1, 1, 0, 0, 0, 0,   & ! Li
        4, 2, 1, 2, 1, 1, 2, 5, 1, 1, 1, 2, 0, 0, 0,   & ! Be
        4, 4, 3, 3, 6, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0,   & ! B
        6, 2, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 0, 0,   & ! C
        6, 3, 1, 2, 2, 1, 2, 3, 1, 2, 0, 0, 0, 0, 0,   & ! N
        5, 1, 2, 1, 4, 1, 5, 1, 1, 1, 1, 0, 0, 0, 0,   & ! O
        4, 2, 4, 2, 2, 2, 2, 3, 1, 1, 0, 0, 0, 0, 0,   & ! F
        6, 2, 3, 1, 2, 8, 2, 2, 0, 0, 0, 0, 0, 0, 0,   & ! Na
        5, 2, 2, 2, 2, 1, 2, 4, 1, 1, 2, 1, 1, 0, 0,   & ! Mg
        6, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1,   & ! Al
        5, 4, 4, 3, 1, 2, 1, 3, 1, 1, 1, 0, 0, 0, 0,   & ! Si
        5, 4, 4, 3, 1, 2, 1, 3, 1, 1, 1, 0, 0, 0, 0,   & ! P
        4, 1, 8, 2, 1, 2, 1, 3, 1, 1, 1, 1, 0, 0, 0,   & ! S
        4, 7, 2, 2, 1, 1, 2, 3, 1, 1, 1, 1, 0, 0, 0],  & ! Cl
        shape(sg0_cnt))
!  Lebedev order of each sector (18 -> 26 and 74 -> 86 substitutions
!  applied, see above)
  integer, parameter :: sg0_leb(SG0_MAXSEC, SG_NELEM) = reshape([ &
        6, 26, 26,  38,  86, 110, 146,  86,  50,  38, 26,  0,  0,  0, 0,  & ! H
        6, 26, 26,  38,  86, 110, 146,  86,  50,  38, 26,  0,  0,  0, 0,  & ! Li
        6, 26, 26,  38,  86,  86, 110, 146,  50,  38, 26,  6,  0,  0, 0,  & ! Be
        6, 26, 38,  86, 146,  38,   6,   0,   0,   0,  0,  0,  0,  0, 0,  & ! B
        6, 26, 26,  38,  50,  86, 110, 146, 170, 146, 86, 38, 26,  0, 0,  & ! C
        6, 26, 26,  38,  86, 110, 170, 146,  86,  50,  0,  0,  0,  0, 0,  & ! N
        6, 26, 26,  38,  50,  86, 110,  86,  50,  38,  6,  0,  0,  0, 0,  & ! O
        6, 38, 50,  86, 110, 146, 110,  86,  50,   6,  0,  0,  0,  0, 0,  & ! F
        6, 26, 26,  38,  50, 110,  86,   6,   0,   0,  0,  0,  0,  0, 0,  & ! Na
        6, 26, 26,  38,  50,  86, 110, 146, 110,  86, 38, 26,  6,  0, 0,  & ! Mg
        6, 26, 26,  38,  50,  86,  86, 146, 170, 110, 86, 86, 26, 26, 6,  & ! Al
        6, 26, 38,  50,  86, 110, 146, 170,  86,  50,  6,  0,  0,  0, 0,  & ! Si
        6, 26, 38,  50,  86, 110, 146, 170,  86,  50,  6,  0,  0,  0, 0,  & ! P
        6, 26, 26,  38,  50,  86, 110, 170, 146, 110, 50,  6,  0,  0, 0,  & ! S
        6, 26, 26,  38,  50,  86, 110, 170, 146, 110, 86,  6,  0,  0, 0], & ! Cl
        shape(sg0_leb))

!  Gauss nodes/weights on (0,1) for the weight ln^2 x (moments
!  m_k = 2/(k+1)^3), used by the SG-0 MultiExp radial grid.  Generated
!  at 200-digit precision via Golub-Welsch (scripts/
!  sg0_multiexp_nodes.py); moments reproduced to ~1e-194.  Sorted by
!  ascending radius r = -ln x (descending x).
  real(kind=dp), parameter :: me23_x(23) = [ &
        0.98868121412417917d0, 0.96979055758591744d0, 0.94296495577173701d0, &
        0.90865042304171716d0, 0.86743287798275914d0, 0.82001826057974256d0, &
        0.76721868628476330d0, 0.70993783394310550d0, 0.64915496391179920d0, &
        0.58590767062942112d0, 0.52127362052638902d0, 0.45635157223202704d0, &
        0.39224199352412352d0, 0.33002759271630692d0, 0.27075407329850232d0, &
        0.21541139722601381d0, 0.16491579716600480d0, 0.12009269427713703d0, &
        0.081660512821455929d0, 0.050215014094683999d0, 0.026212787562513862d0, &
        0.0099491128468611853d0, 0.0015058924745840717d0]
  real(kind=dp), parameter :: me23_w(23) = [ &
        1.9205788879728201d-6, 2.1565953939261300d-5, 0.00010572868056731378d0, &
        0.00034755788975770020d0, 0.00089889253062980984d0, 0.0019785635068497310d0, &
        0.0038757714610574311d0, 0.0069477499303844705d0, 0.011611019656801442d0, &
        0.018325604666401505d0, 0.027571576765107484d0, 0.039817159733692048d0, &
        0.055477247414684803d0, 0.074860364849398577d0, 0.098100418872550932d0, &
        0.12506617470764026d0, 0.15523427291530433d0, 0.18749582392684797d0, &
        0.21982849560935217d0, 0.24866202472840894d0, 0.26742812676397013d0, &
        0.26236396365964760d0, 0.19397997519811811d0]
  real(kind=dp), parameter :: me26_x(26) = [ &
        0.99104255389177476d0, 0.97606029628496512d0, 0.95471297844715660d0, &
        0.92727872393893840d0, 0.89412727813918990d0, 0.85570731112497829d0, &
        0.81253906250257684d0, 0.76520686385149628d0, 0.71435096447486728d0, &
        0.66065863649024736d0, 0.60485464446878273d0, 0.54769119749545300d0, &
        0.48993751506621983d0, 0.43236914470191774d0, 0.37575717144832967d0, &
        0.32085745807052075d0, 0.26840004931968460d0, 0.21907886291165968d0, &
        0.17354177117918734d0, 0.13238114512422221d0, 0.096124873966380717d0, &
        0.065227756094948038d0, 0.040062890303615981d0, 0.020911970701014455d0, &
        0.0079508350834211655d0, 0.0012118959531442052d0]
  real(kind=dp), parameter :: me26_w(26) = [ &
        9.5039183893267110d-7, 1.0687185195322948d-5, 5.2504259615681300d-5, &
        0.00017307202304493877d0, 0.00044916603235275513d0, 0.00099280471150256854d0, &
        0.0019544295871422983d0, 0.0035237959667716519d0, 0.0059282775235989670d0, &
        0.0094283196275203881d0, 0.014309796296644065d0, 0.020873023457481862d0, &
        0.029418139625084388d0, 0.040226455838504916d0, 0.053537150852983894d0, &
        0.069518256800206744d0, 0.088230076743000778d0, 0.10957766175447101d0, &
        0.13324603342450455d0, 0.15860582735470473d0, 0.18456387321242483d0, &
        0.20930163443099709d0, 0.22975862025927898d0, 0.24044079379163951d0, &
        0.22996905094874790d0, 0.16590959790074127d0]

  type :: saved_HF_info !< keeps HF exchange from input
    logical :: alpha = .false.
    logical :: beta = .false.
    logical :: mu = .false.
    logical :: hfscale = .false.
    logical :: do = .false.
    real(kind=dp) :: saved_alpha = -1.0_dp
    real(kind=dp) :: saved_beta = -1.0_dp
    real(kind=dp) :: saved_mu = -1.0_dp
    real(kind=dp) :: saved_hfscale = -1.0_dp

  contains

    procedure :: save_HF => save_dft_HF_exchange_from_input
    procedure :: update_HF => update_dft_HF_exchange_from_input
  end type

contains

  subroutine save_dft_HF_exchange_from_input(this, infos)
    use types, only: information
    implicit none
    class(saved_HF_info), intent(inout) :: this
    type(information), intent(inout) :: infos

    if (infos%dft%cam_flag) then
      if (infos%dft%cam_alpha /= -1.0_dp) then
        this%saved_alpha = infos%dft%cam_alpha
        this%alpha = .true.
      end if
      if (infos%dft%cam_beta /= -1.0_dp) then
        this%saved_beta = infos%dft%cam_beta
        this%beta = .true.
      end if
      if (infos%dft%cam_mu /= -1.0_dp) then
        this%saved_mu = infos%dft%cam_mu
        this%mu = .true.
      end if
      if (this%alpha.or.this%beta.or.this%mu) &
          this%do = .true.
    else
      if (infos%dft%hfscale /= -1.0_dp) then
        this%saved_hfscale = infos%dft%hfscale
        this%hfscale = .true.
        this%do = .true.
      end if
    end if


  end subroutine save_dft_HF_exchange_from_input

  subroutine update_dft_HF_exchange_from_input(this, infos)
    use types, only: information
    implicit none
    class(saved_HF_info), intent(inout) :: this
    type(information), intent(inout) :: infos

    real(kind=dp) :: scale
    character(len=80), parameter :: format = &
          '(11x,a,":",t22,"|", t24, e12.5, t37, "-|>", t41, e12.5, t54, "|")'

    if (infos%dft%cam_flag) then
      write(*,'(2x,a)') "CAM-B3LYP with tuned Hartree-Fock exchange from the input."
      write(*, '(5x,"CAM parametres: |   It was     |   It become    |")')
      if (this%alpha) then
         scale = this%saved_alpha
      else
         scale =  infos%dft%cam_alpha
      end if
      write(*, fmt=format) "Alpha", 0.19_dp, scale
      if (this%alpha) infos%dft%cam_alpha = this%saved_alpha

      if (this%beta) then
         scale = this%saved_beta
      else
         scale =  infos%dft%cam_beta
      end if
      write(*, fmt=format) "Beta", 0.46_dp, scale
      if (this%beta) infos%dft%cam_beta = this%saved_beta

      if (this%mu) then
         scale = this%saved_mu
      else
         scale =  infos%dft%cam_mu
      end if
      write(*, fmt=format) "mu", 0.33_dp, scale
      if (this%mu) infos%dft%cam_mu = this%saved_mu
    else
      write(*,'(2x,a)') "Tuned Hartree-Fock exchange from the input."
      write(*, '(10x,"Exact HF exchange:")')
      if (this%hfscale) then
         scale = this%saved_hfscale
      else
         scale =  infos%dft%hfscale
      end if
      write(*, fmt=format) "HF scale", infos%dft%hfscale, scale
      if (this%hfscale) infos%dft%hfscale = this%saved_hfscale
      write(*, '(2x,a)') "Please cite the following works when using this option:"
      write(*,fmt='(3a)') "[1] W. Park, A. Lashkaripour, K. Komarov, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., ??, ?? (2024); ", &
            "DOI: 10.1021/acs.jctc.4c00640"
      write(*,fmt='(3a)') "[2] K. Komarov, W. Park, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., 19, 7671-7684 (2023); ", &
            "DOI: 10.1021/acs.jctc.3c00884"
    end if
    write(*,*)

  end subroutine update_dft_HF_exchange_from_input

  subroutine dft_initialize(infos, basis, molGrid, orbitals_cutoff, verbose, need_functional)
    use basis_tools, only: basis_set
    use types, only: information

    implicit none

    type(basis_set), intent(inout) :: basis
    type(information), intent(inout) :: infos
    type(dft_grid_t), intent(inout) :: molGrid
    real(kind=dp), optional :: orbitals_cutoff
    logical, optional :: verbose
    logical, optional :: need_functional

    real(kind=dp) :: logtol
    type(dft_grid_pruned_t) :: pruned

!   Setup sreening parameters
    logtol = -log(1.0e-10_dp)
    if (present(orbitals_cutoff)) logtol = -log(orbitals_cutoff)
    call basis%set_screening(logtol)

!   Set grid DFT options
    call dft_set_options(infos, pruned, need_functional)

!   Initialize grid
    call dft_prepare_grid(infos, basis, molGrid, pruned, verbose)

  end subroutine

!>  @brief Calculates atomic distances
  subroutine get_atomic_distances(xyz, rij)

      implicit none
      real(kind=dp), intent(in) :: xyz(:,:)

      real(kind=dp), intent(out) :: rij(:,:)

      integer :: i, j

      do i = 1, ubound(xyz,2)
         rij(i,i) = 0.0d0
         do j = 1, i-1
            rij(i,j) = norm2(xyz(:,i)-xyz(:,j))
            rij(j,i) = rij(i,j)
        end do
      end do
  end subroutine

  subroutine emovlp(nrad,rads,wts,lmn,zeta,bragg,s)
    use constants, only: pi
    implicit none
    integer, intent(in) :: nrad, lmn
    real(kind=dp), intent(in) :: zeta, bragg
    real(kind=dp), intent(in) :: rads(:), wts(:)
    real(kind=dp), intent(out) :: s
    integer :: idf, i
    real(kind=dp) :: gnorm, r, w, gto

    idf = 1
    do i = 0, lmn
     idf = idf*(2*i+1)
    end do
    gnorm = zeta**(2*lmn+3) * 2**(4*lmn+7)
    gnorm = gnorm/(pi*idf**2)
    gnorm = gnorm**(0.25d+00)

    s = 0
    do i = 1, nrad
      r = bragg*rads(i)
      w = (bragg**3)*wts(i)
      gto = gnorm*r**lmn*exp(-zeta*r*r)
      s = s + w*(gto*gto)
    end do
  end subroutine

  subroutine dftclean(infos)
    use types, only: information
    use libxc, only: libxc_destroy
    type(information), intent(inout) :: infos
    call libxc_destroy(infos%functional)
  end subroutine

  subroutine dft_set_options(infos, pruned, need_functional)
    use iso_c_binding, only: c_null_char
    use messages, only: show_message, WITH_ABORT
    use strings, only: c_f_char
    use types, only: information
    use libxc, only: libxc_input

    implicit none

    type(information), intent(inout) :: infos
    type(dft_grid_pruned_t), intent(inout) :: pruned
    logical, optional, intent(in) :: need_functional
    type(saved_HF_info) :: saved_hf
    logical :: need_func

    integer :: iatm, nrad
    character(len=20) :: xc_func_name
    integer :: nat, i, slen, ntyps
    character(:), allocatable :: pruned_name
    logical :: is_sg3
    integer :: z, ie, nsec, maxsec, nang_fallback
    integer :: zmap(SG_NELEM)

    need_func = .true.
    if (present(need_functional)) need_func = need_functional

!   Default radial/angular grid is 96/302 for LDA/GGA.
    nrad     = infos%dft%grid_rad_size
    nat = ubound(infos%atoms%zn, 1)
    allocate(pruned%rad_id(nat), source=1)

    xc_func_name = c_f_char(infos%dft%xc_functional_name)

    if (.not. infos%dft%grid_pruned) then
      pruned%ngrids = 1
      allocate(pruned%nang(1,1), pruned%radii(1,1))
      pruned%nang(1,1) = infos%dft%grid_ang_size
      pruned%radii(1,1) = 1.0d+30

      write(iw,'(/5X,"Lebedev grid-based DFT options"/&
                &5X,30("-")/&
                &5X,"XC functional: ",A/&
                &5X,"NRAD  =",I8,5X,"NLEB  =",I8/&
                &5X,"THRESH=",1P,E12.2)') &
                  trim(xc_func_name), &
                  nrad, pruned%nang(1,1), &
                  infos%dft%grid_density_cutoff

    else

!     Set parameters for pruned grids
      slen = ubound(infos%dft%grid_pruned_name, 1)
      allocate(character(len=slen) :: pruned_name)
      do i = 1, slen
        if (infos%dft%grid_pruned_name(i) == c_null_char) exit
        pruned_name(i:i) = infos%dft%grid_pruned_name(i)
      end do

      select case (trim(pruned_name))
      case ("SG1")
        pruned%ngrids = 5
        ntyps = 4
        allocate(pruned%nang(pruned%ngrids, ntyps), &
                 pruned%radii(pruned%ngrids, ntyps))
        pruned%radii = sg1rads
        do i = 1, ntyps
          pruned%nang(:,i) = sg1grids
        end do
        do iatm = 1, nat
          ! sg1atoms are inclusive upper bounds of the period:
          ! H-He (Z<=2), Li-Ne (Z<=10), Na-Ar (Z<=18), heavier
          do i = 1, ntyps
            if (int(infos%atoms%zn(iatm))<=sg1atoms(i)) exit
          end do
          pruned%rad_id(iatm) = min(i, ntyps)
        end do
        ! SG1 is undefined above Ar: heavy atoms (type 4) are truly
        ! unpruned, i.e. a single 194-point sphere at all radii
        allocate(pruned%nang_override(ntyps), source=0)
        pruned%nang_override(4) = 194

        write(iw,'(/5X,"Standard Grid 1 (SG1)"/&
                  &5X,21("-")/&
                  &5X,"XC functional: ",A/&
                  &5X,"THRESH=",1P,E12.2)') &
                    trim(xc_func_name), &
                    infos%dft%grid_density_cutoff

      case ("SG2", "SG3")
        is_sg3 = trim(pruned_name) == "SG3"
        if (is_sg3) then
          pruned%nrad = SG3_NRAD
          maxsec = SG3_MAXSEC
          nang_fallback = 590
        else
          pruned%nrad = SG2_NRAD
          maxsec = SG2_MAXSEC
          nang_fallback = 302
        end if

!       One atom type per supported element present in the system;
!       type 1 is the fallback (He, Ne, Ar, Z > 18, dummy atoms):
!       unpruned nang_fallback-point grid on the standard radial grid.
        zmap = 0
        ntyps = 1
        do iatm = 1, nat
          z = int(abs(infos%atoms%zn(iatm))+1.0d-5)
          ie = 0
          if (abs(abs(infos%atoms%zn(iatm))-z) <= 1.0d-5 .and. &
              z >= 1 .and. z <= 17) then
            do i = 1, SG_NELEM
              if (sg_elem_z(i) == z) then
                ie = i
                exit
              end if
            end do
          end if
          if (ie > 0) then
            if (zmap(ie) == 0) then
              ntyps = ntyps+1
              zmap(ie) = ntyps
            end if
            pruned%rad_id(iatm) = zmap(ie)
          else
            pruned%rad_id(iatm) = 1
          end if
        end do

        pruned%ngrids = maxsec
        pruned%nrad_types = ntyps
        allocate(pruned%nang(maxsec, ntyps), source=0)
        allocate(pruned%nradPerRegion(maxsec, ntyps), source=0)
        allocate(pruned%radii(maxsec, ntyps), source=1.0d30)
        allocate(pruned%de2_alpha(ntyps), source=0.0_dp)
        allocate(pruned%de2_rmax(ntyps), source=0.0_dp)
        pruned%radial_id = pruned%rad_id

!       Fallback type: single unpruned region, standard radial grid
        pruned%nang(1,1) = nang_fallback

!       Element types: index-based sectors on the per-element DE2 grid
        do ie = 1, SG_NELEM
          i = zmap(ie)
          if (i == 0) cycle
          if (is_sg3) then
            nsec = sg3_nsec(ie)
            pruned%nang(1:nsec, i) = sg3_leb(1:nsec, ie)
            pruned%nradPerRegion(1:nsec, i) = sg3_cnt(1:nsec, ie)
            pruned%de2_alpha(i) = sg3_alpha(ie)
          else
            nsec = SG2_MAXSEC
            pruned%nang(1:nsec, i) = sg2_leb(1:nsec, ie)
            pruned%nradPerRegion(1:nsec, i) = sg2_cnt(1:nsec, ie)
            pruned%de2_alpha(i) = sg2_alpha(ie)
          end if
          pruned%de2_rmax(i) = sg_de2_rmax(ie)
        end do

        write(iw,'(/5X,"Standard Grid ",A," (",A,") of Dasgupta and Herbert"/&
                  &5X,40("-")/&
                  &5X,"XC functional: ",A/&
                  &5X,"NRAD  =",I8,"   (Mitani DE2 radial grid)"/&
                  &5X,"THRESH=",1P,E12.2)') &
                    pruned_name(3:3), trim(pruned_name), &
                    trim(xc_func_name), pruned%nrad, &
                    infos%dft%grid_density_cutoff

      case ("SG0")
        maxsec = SG0_MAXSEC
!       The radial storage must fit both the MultiExp grids (up to 26
!       nodes) and the standard radial grid of the fallback atoms
        pruned%nrad = max(nrad, 26)

!       One atom type per supported element present in the system;
!       types 1-4 are the SG1 fallback (He, Ne, Z >= 18, non-integer
!       nuclear charges), typed by period as in the SG1 case; element
!       types are 5, 6, ...  Radial types: 1 is the standard grid
!       (fallback); element type 4+k uses MultiExp radial column 1+k.
        zmap = 0
        ntyps = 4
        allocate(pruned%radial_id(nat), source=1)
        do iatm = 1, nat
          z = int(abs(infos%atoms%zn(iatm))+1.0d-5)
          ie = 0
          if (abs(abs(infos%atoms%zn(iatm))-z) <= 1.0d-5 .and. &
              z >= 1 .and. z <= 17) then
            do i = 1, SG_NELEM
              if (sg_elem_z(i) == z) then
                ie = i
                exit
              end if
            end do
          end if
          if (ie > 0) then
            if (zmap(ie) == 0) then
              ntyps = ntyps+1
              zmap(ie) = ntyps
            end if
            pruned%rad_id(iatm) = zmap(ie)
            pruned%radial_id(iatm) = zmap(ie)-3
          else
            ! SG1 fallback type by period (see the SG1 case)
            do i = 1, 4
              if (int(infos%atoms%zn(iatm))<=sg1atoms(i)) exit
            end do
            pruned%rad_id(iatm) = min(i, 4)
            pruned%radial_id(iatm) = 1
          end if
        end do

        pruned%ngrids = maxsec
        pruned%nrad_types = ntyps-3
        allocate(pruned%nang(maxsec, ntyps), source=0)
        allocate(pruned%nradPerRegion(maxsec, ntyps), source=0)
        allocate(pruned%radii(maxsec, ntyps), source=1.0d30)
        allocate(pruned%nang_override(ntyps), source=0)
        allocate(pruned%de2_alpha(pruned%nrad_types), source=0.0_dp)
        allocate(pruned%de2_rmax(pruned%nrad_types), source=0.0_dp)
        allocate(pruned%rad_npts(pruned%nrad_types), source=0)
        allocate(pruned%me_rscale(pruned%nrad_types), source=0.0_dp)

!       Fallback types 1-4: the SG1 scheme (radius-based regions) on
!       the standard radial grid; heavy atoms (type 4) are truly
!       unpruned, i.e. a single 194-point sphere at all radii
        do i = 1, 4
          pruned%nang(1:5, i) = sg1grids
          pruned%radii(1:5, i) = sg1rads(:, i)
        end do
        pruned%nang_override(4) = 194

!       Element types: index-based sectors on the per-element MultiExp
!       radial grid
        do ie = 1, SG_NELEM
          i = zmap(ie)
          if (i == 0) cycle
          nsec = sg0_nsec(ie)
          pruned%nang(1:nsec, i) = sg0_leb(1:nsec, ie)
          pruned%nradPerRegion(1:nsec, i) = sg0_cnt(1:nsec, ie)
          pruned%rad_npts(i-3) = sg0_nrad(ie)
          pruned%me_rscale(i-3) = sg0_rscale(ie)
        end do

        write(iw,'(/5X,"Standard Grid 0 (SG0) of Chien and Gill"/&
                  &5X,39("-")/&
                  &5X,"XC functional: ",A/&
                  &5X,"NRAD  =   23/26   (MultiExp radial grid)"/&
                  &5X,"THRESH=",1P,E12.2)') &
                    trim(xc_func_name), &
                    infos%dft%grid_density_cutoff

      case default
        call show_message('Unknown pruned grid name', WITH_ABORT)
      end select

    end if


!   New we set the DFT XC functionals...
    if (trim(xc_func_name) /= "") then
      ! save HFscale, or cam_alpha,beta,mu from input
      call saved_HF%save_HF(infos)

      call libxc_input(functional_name=trim(xc_func_name), &
                       dft_params=infos%dft, &
                       tddft_params=infos%tddft, &
                       functional=infos%functional)

      ! update HFscale, or cam_alpha,beta,mu from input
      if(saved_HF%do) call saved_HF%update_HF(infos)
    else if (need_func) then
      call show_message('Please, specify functional in the input file', WITH_ABORT)
    end if

  end subroutine

  subroutine dft_prepare_grid(infos, basis, molGrid, pruned, verbose)
      use basis_tools, only: basis_set
      use dft_radial_grid_types, only: get_radial_grid
      use mod_dft_fuzzycell, only: dft_fc_blk
      use mod_grid_storage, only: atomic_grid_t
      use bragg_slater_radii, only: set_bragg_slater, &
          BRSL_NUM_ELEMENTS, &
          BRSL_TYPE_GILL, &
          BRSL_TYPE_TA
      use types, only: information

      implicit none

      type(basis_set), intent(in) :: basis
      type(information), intent(in) :: infos
      type(dft_grid_t), intent(inout) :: molGrid
      type(dft_grid_pruned_t), intent(in) :: pruned
      logical, optional :: verbose

      integer :: i, igrid, nat, iat, maxpt_per_atom, nrad
      integer :: bstype
      integer :: grid_id
      integer :: max_ang_pts
      integer :: ngr, rtid, nrad_at, override
      real(kind=dp) :: dftthr0
      real(KIND=dp) :: brsl_radii(BRSL_NUM_ELEMENTS)
      logical :: verbose_

      real(kind=dp), allocatable :: txyz(:), twght(:)
      real(kind=dp), allocatable :: wtab(:,:,:)
      real(kind=dp), allocatable :: rij(:,:), aij(:,:)
      real(kind=dp), allocatable :: bsrad(:)

      type(atomic_grid_t) :: atomic_grid

      verbose_ = .false.
      if (present(verbose)) verbose_ = verbose

      nat = ubound(infos%atoms%zn, 1)
      max_ang_pts = maxval(pruned%nang)
      if (allocated(pruned%nang_override)) &
        max_ang_pts = max(max_ang_pts, maxval(pruned%nang_override))
      nrad = infos%dft%grid_rad_size
      ! A pruned grid may prescribe its own radial grid size
      if (pruned%nrad > 0) nrad = pruned%nrad
      maxpt_per_atom = nrad*max_ang_pts

      allocate(&
        txyz(max_ang_pts*3), &
        twght(max_ang_pts), &
        rij(nat,nat), &
        aij(nat,nat), &
        bsrad(nat), &
        source=0.0d0)

!     Init storage for the grid
      call molGrid%reset(nat, maxpt_per_atom, nRad, pruned%nrad_types)

!     Print out DFT info
      if (verbose_) then
        dftthr0=1.0d-03/(maxpt_per_atom*nat)
        if(dftthr0.lt.1.1d-15) then
          write(iw,'(5x, "All DFT thresholds are turned off.")')
        else
          write(iw,'(5x, "DFT Threshold         =",e10.3)') dftthr0
        end if
      end if

!     Set up Bragg-Slater radii for atoms
      select case(infos%dft%rad_grid_type)
      case (3)
        bstype = BRSL_TYPE_TA
      case default
        bstype = BRSL_TYPE_GILL
      end select
      call set_bragg_slater(brsl_radii, bstype)

      do i = 1, nat
        bsrad(i) = bragg_slater_radius(brsl_radii, infos%atoms%zn(i))
      end do

!     Set up radial grid (the standard grid is radial type 1)
      call get_radial_grid(molGrid%rad_pts(:,1), molGrid%rad_wts(:,1), &
              nrad, infos%dft%rad_grid_type)

!     Element-specific radial grids, absolute radii.
!     MultiExp (SG-0): per-element node count and scaling radius;
!     DE2 (SG-2/SG-3): the innermost/outermost nodes are pinned to
!     SG_DE2_RMIN and the element-specific R_max (cuEST convention).
!     Unused trailing rows of a MultiExp column stay zero and are
!     never referenced (per-atom grids are sliced to the per-type
!     node count below).
      do i = 2, pruned%nrad_types
        if (allocated(pruned%me_rscale)) then
          nrad_at = pruned%rad_npts(i)
          call multiexp_radial_grid(nrad_at, pruned%me_rscale(i), &
                  molGrid%rad_pts(1:nrad_at,i), molGrid%rad_wts(1:nrad_at,i))
        else
          call de2_radial_grid(nrad, pruned%de2_alpha(i), &
                  SG_DE2_RMIN, pruned%de2_rmax(i), &
                  molGrid%rad_pts(:,i), molGrid%rad_wts(:,i))
        end if
      end do

!     Per-atom radial grid types
      if (allocated(pruned%radial_id)) &
        molGrid%radTypeId(1:nat) = pruned%radial_id(1:nat)

!     Pre-compute atomic distances
      call get_atomic_distances(infos%atoms%xyz, rij)

!     Tag non-real atoms present in the system
      molGrid%dummyAtom(:nat) = bsrad(:nat) == 0.0_dp

!     Find nearest neighbours for all atoms
      call molGrid%find_neighbours(rij, partFunType=infos%dft%dft_partfun)

!     Compute atomic grids for each atom
      do iat = 1, nat
        grid_id = pruned%rad_id(iat)
        rtid = molGrid%radTypeId(iat)

        override = 0
        if (allocated(pruned%nang_override)) &
          override = pruned%nang_override(grid_id)

        if (override > 0) then
!         Unpruned atom type: a single angular grid at all radii
!         (add_atomic_grid extends a single region to all radial shells)
          if (allocated(atomic_grid%sph_nrad)) &
            deallocate(atomic_grid%sph_nrad)
          atomic_grid%sph_npts = [override]
          atomic_grid%sph_radii = [9999999.9d0]
          call molGrid%spherical_grids%add_grid(override)
        else
!         Number of regions used by this atom type and the pruning mode:
!         index-based sectors (nradPerRegion > 0) vs radius-based regions
          ngr = pruned%ngrids
          if (allocated(pruned%nradPerRegion)) then
            if (any(pruned%nradPerRegion(:, grid_id) > 0)) then
              ngr = count(pruned%nradPerRegion(:, grid_id) > 0)
              atomic_grid%sph_nrad = pruned%nradPerRegion(1:ngr, grid_id)
            else
              ngr = count(pruned%nang(:, grid_id) > 0)
              if (allocated(atomic_grid%sph_nrad)) &
                deallocate(atomic_grid%sph_nrad)
            end if
          end if

          atomic_grid%sph_npts = pruned%nang(1:ngr, grid_id)
          atomic_grid%sph_radii = pruned%radii(1:ngr, grid_id)

!         Set angular Lebedev grid(s)
          do igrid = 1, ngr
!           get the unit lebedev sphere (no-op if already stored)
            call molGrid%spherical_grids%add_grid(pruned%nang(igrid, grid_id))
          end do
        end if
        atomic_grid%idAtm = iat

!       Set the radial grid of the atom.  Element-specific grids
!       (types >= 2) store absolute radii: use a unit effective radius;
!       the standard grid (type 1) is scaled by the Bragg-Slater
!       radius.  The per-atom grid is sliced to the per-type node
!       count when one is prescribed (MultiExp grids of SG-0).
        nrad_at = nrad
        if (allocated(pruned%rad_npts)) then
          if (pruned%rad_npts(rtid) > 0) nrad_at = pruned%rad_npts(rtid)
        end if
        atomic_grid%rad_pts = molGrid%rad_pts(1:nrad_at, rtid)
        atomic_grid%rad_wts = molGrid%rad_wts(1:nrad_at, rtid)
        if (rtid > 1) then
          atomic_grid%rAtm = 1.0_dp
        else
          atomic_grid%rAtm = bsrad(iat)
        end if

        call molGrid%add_atomic_grid(atomic_grid)
      end do

!     Assemble molecular grid from atomic grids

!     Do Becke's fuzzy cell
      select case (infos%dft%dft_bfc_algo)
      case(0)
!       SSF algorithm:
!       various partitioning functions, no surface shifting
        call dft_fc_blk(molgrid, infos%dft%dft_partfun, &
                infos%atoms%xyz,basis%at_mx_dist2,rij,nat,wtab)
      case (1)
!       Precompute surface shifting parameters
        call setaij(aij, nat, bsrad)
!       Becke's algorithm:
!       4th deg. Becke's polynomial and surface shifting
        call dft_fc_blk(molgrid, infos%dft%dft_partfun, &
                infos%atoms%xyz,basis%at_mx_dist2,rij,nat,wtab,aij)

      end select

      call molGrid%compress

      if (verbose_) then
        write(iw,'(5X,"Molecular grid: ",I0," points in ",I0," slices")') &
              sum(molGrid%nTotPts(1:molGrid%nSlices)), molGrid%nSlices
      end if

  end subroutine

!> @brief Mitani double-exponential (DE2) radial quadrature
!> @details M. Mitani, Theor. Chem. Acc. 130, 645 (2011);
!>   M. Mitani, Y. Yoshioka, Theor. Chem. Acc. 131, 1169 (2012).
!>   Nodes and weights (the weights include the r^2 Jacobian):
!>     x_i = x_start + (i-1)*h,  i = 1..nr
!>     r_i = exp(alpha*x_i - exp(-x_i))                        [bohr]
!>     w_i = h*(alpha + exp(-x_i))*exp(3*alpha*x_i - 3*exp(-x_i))
!>   x_start and x_end are pinned to the innermost/outermost radial
!>   nodes:  alpha*x - exp(-x) = ln(rmin) resp. ln(rmax), solved by
!>   Newton iteration (the left-hand side is strictly increasing),
!>   and h = (x_end - x_start)/(nr - 1).  This matches the SG-2/SG-3
!>   convention of NVIDIA cuEST (CUDALibrarySamples).
!> @param[in]   nr     number of radial points
!> @param[in]   alpha  DE2 alpha parameter (element-specific)
!> @param[in]   rmin   innermost radial node, bohr
!> @param[in]   rmax   outermost radial node, bohr
!> @param[out]  r      radial nodes, absolute bohr
!> @param[out]  w      radial weights including the r^2 Jacobian
  pure subroutine de2_radial_grid(nr, alpha, rmin, rmax, r, w)
    implicit none
    integer, intent(in) :: nr
    real(kind=dp), intent(in) :: alpha, rmin, rmax
    real(kind=dp), intent(out) :: r(:), w(:)

    integer :: i
    real(kind=dp) :: h, x, x_start, x_end

    x_start = de2_solve_x(alpha, log(rmin), -2.3_dp)
    x_end = de2_solve_x(alpha, log(rmax), 1.25_dp)
    h = (x_end-x_start)/(nr-1)

    do i = 1, nr
      x = x_start+(i-1)*h
      r(i) = exp(alpha*x-exp(-x))
      w(i) = h*(alpha+exp(-x))*exp(3.0_dp*alpha*x-3.0_dp*exp(-x))
    end do

  end subroutine de2_radial_grid

!> @brief Solve alpha*x - exp(-x) = lnr for x by Newton iteration
!> @details The left-hand side is strictly increasing in x for
!>   alpha > 0, so the root is unique.
  pure function de2_solve_x(alpha, lnr, x0) result(x)
    implicit none
    real(kind=dp), intent(in) :: alpha, lnr, x0
    real(kind=dp) :: x

    integer :: iter
    real(kind=dp) :: f, xnew

    x = x0
    do iter = 1, 100
      f = alpha*x-exp(-x)-lnr
      xnew = x-f/(alpha+exp(-x))
      if (abs(xnew-x) < 1.0d-14) then
        x = xnew
        exit
      end if
      x = xnew
    end do

  end function de2_solve_x

!> @brief MultiExp radial quadrature (SG-0)
!> @details P.M.W. Gill, S.-H. Chien, J. Comput. Chem. 24, 732 (2003).
!>   Gauss quadrature on (0,1) for the weight function ln^2 x
!>   (moments m_k = 2/(k+1)^3), mapped to (0,inf) by r = -R ln x:
!>     r_i = -R ln(x_i)                                       [bohr]
!>     w_i = R^3 omega_i / x_i
!>   The weights include the r^2 Jacobian, matching the DE2
!>   convention (the consumer multiplies by rAtm^3 = 1).  Nodes and
!>   weights are tabulated for n = 23 and 26 (the SG-0 sizes), sorted
!>   by ascending radius.
!> @param[in]   nr      number of radial points (23 or 26)
!> @param[in]   rscale  element-specific scaling radius R, bohr
!> @param[out]  r       radial nodes, absolute bohr
!> @param[out]  w       radial weights including the r^2 Jacobian
  subroutine multiexp_radial_grid(nr, rscale, r, w)
    implicit none
    integer, intent(in) :: nr
    real(kind=dp), intent(in) :: rscale
    real(kind=dp), intent(out) :: r(:), w(:)

    integer :: i

    select case (nr)
    case (23)
      do i = 1, nr
        r(i) = -rscale*log(me23_x(i))
        w(i) = rscale**3*me23_w(i)/me23_x(i)
      end do
    case (26)
      do i = 1, nr
        r(i) = -rscale*log(me26_x(i))
        w(i) = rscale**3*me26_w(i)/me26_x(i)
      end do
    case default
      call show_message('MultiExp radial grid: unsupported size', &
              WITH_ABORT)
    end select

  end subroutine multiexp_radial_grid

  subroutine dftexcor(basis,molGrid,iscftyp,fa,fb,coeffa,coeffb,nbf,nbf_tri,eexc,totele,totkin, infos)
    use basis_tools, only: basis_set
    use mod_dft_gridint_energy, only: dmatd_blk
    use types, only: information

    implicit none
    type(information), intent(in) :: infos
    type(basis_set) :: basis
    type(dft_grid_t), intent(in) :: molGrid
    real(kind=dp), intent(inout) :: fa(*),fb(*),coeffa(*),coeffb(*)
    integer, intent(in) :: iscftyp, nbf, nbf_tri
    real(kind=dp), intent(out) :: eexc, totele, totkin
    integer :: nang, maxl
    logical :: urohf

    urohf = iscftyp/=1

    fa(1:nbf_tri) = 0.0d0
    if(iscftyp>=2) fb(1:nbf_tri) = 0.0d0

    maxl = maxval(basis%am)
    nang = maxl+1+1

    totele = 0.0d0
    totkin = 0.0d0
    eexc   = 0.0d0
    call dmatd_blk(basis, molGrid, coeffa,coeffb,fa,fb, &
                     eexc,totele,totkin, &
                     nang,nbf,infos%dft%grid_density_cutoff,urohf, infos)

  end subroutine

!> @brief Analytical DFT gradient
  subroutine dftder(basis, infos, molGrid)
    use mathlib, only: unpack_matrix
    use mod_dft_gridint_grad, only: derexc_blk
    use types, only: information
    use oqp_tagarray_driver

    implicit none

    character(len=*), parameter :: subroutine_name = "dftder"

    type(basis_set) :: basis
    type(information), intent(inout) :: infos
    type(dft_grid_t), intent(inout) :: molGrid

    integer :: iscftype, num, nat, nbf, nang, nder, maxl
    integer :: iok
    real(kind=dp) :: totele, totkin
    logical :: urohf
    real(kind=dp), allocatable :: tda(:,:), tdb(:,:), dedft(:,:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    integer(4) :: status

    iscftype = infos%control%scftype

    num = basis%nbf
    urohf = iscftype/=1
    nat = infos%mol_prop%natom
    nbf = num
    maxl = maxval(basis%am)
    nang = maxl+1+1

    nder = 1

    if (.not.allocated(tda)) then
      allocate(&
        tda(nbf,nbf), &
        dedft(3,nat), &
        stat=iok)
      if (iok/=0) call show_message('Cannot allocate memory',WITH_ABORT)
    end if

    if (urohf) then
      iok=0
      if(.not.allocated(tdb)) allocate(tdb(nbf,nbf), stat=iok)
      if(iok/=0) call show_message('Cannot allocate memory',WITH_ABORT)
    end if

!   RHF
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
    call check_status(status, module_name, subroutine_name, OQP_DM_A)
    call unpack_matrix(dmat_a,tda,nbf,'U')
!   UHF/ROHF
    if (urohf) then
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b, status)
      call check_status(status, module_name, subroutine_name, OQP_DM_B)
      call unpack_matrix(dmat_b,tdb,nbf,'U')
    end if

    dedft = 0
    totele = 0
    call derexc_blk(basis,molGrid,tda,tdb,dedft, &
                    totele,totkin, &
                    nang,nbf,infos%dft%grid_density_cutoff,urohf, infos)

    infos%atoms%grad(:,:nat) = infos%atoms%grad(:,:nat)+dedft(:,:nat)

  end subroutine

!> @brief Calculate surface shifting parameters
!> @author Vladimir Mironov
!> @date  : Jan, 2019
!> @param[out]  aij   surface shifting parameters
!> @param[in]   nat   number of atoms
  subroutine setaij(aij, nat, bsrad)

    implicit none

    real(kind=dp), intent(out) :: aij(nat,*)
    integer, intent(in) :: nat
    real(kind=dp), intent(in) :: bsrad(:)
    integer :: iatm, jatm
    real(kind=dp) :: radi, radj, chi, chi2

    do iatm = 1, nat
      aij(iatm,iatm) = 0.0d0
      radi = bsrad(iatm)
      if (radi<0.001) then
        aij(1,iatm) = -1
        cycle
      end if

      do jatm = 1, nat
        if (iatm==jatm) cycle

        radj = bsrad(jatm)
        if (radj<0.001) then
          aij(jatm,iatm) = 1
          cycle
        end if
        chi = radi/radj
        chi2 = (chi-1)/(chi+1)
        aij(jatm,iatm) = chi2/(chi2*chi2-1)
        aij(jatm,iatm) = min(aij(jatm,iatm),  0.5)
        aij(jatm,iatm) = max(aij(jatm,iatm), -0.5)
      end do
    end do
  end subroutine

  pure function bragg_slater_radius(element_radii, nuclear_charge) result(radius)
    use physical_constants, only: angstrom_to_bohr
    implicit none
    real(kind=dp), intent(in) :: element_radii(:)
    real(kind=dp), intent(in) :: nuclear_charge
    real(kind=dp) :: radius
    integer :: atomic_number
    real(kind=dp), parameter :: tolerance = 1e-5_dp

    radius = 0.0_dp

    atomic_number = int(abs(nuclear_charge) + tolerance)
    if (abs(abs(nuclear_charge) - atomic_number) <= tolerance) then
      radius = element_radii(atomic_number)
    end if

    radius = radius * angstrom_to_bohr

  end function bragg_slater_radius

end module dft
