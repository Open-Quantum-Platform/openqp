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
!>   types >= 2 are element-specific DE2 grids in absolute bohr
!>   (`de2_alpha`/`rad_typ_z` give alpha and the element).
  type dft_grid_pruned_t
    integer :: nrad = 0
    integer :: ngrids = 1
    integer, allocatable :: nang(:,:)    !< (region, atom type)
    real(kind=dp), allocatable :: radii(:,:) !< (region, atom type)
    integer, allocatable :: rad_id(:)    !< atom -> atom type
    integer, allocatable :: nradPerRegion(:,:) !< (region, atom type); 0 = unused
    integer :: nrad_types = 1            !< number of radial grids
    integer, allocatable :: radial_id(:) !< atom -> radial grid type
    real(kind=dp), allocatable :: de2_alpha(:) !< DE2 alpha of radial type
    integer, allocatable :: rad_typ_z(:) !< element of radial type (0: standard)
  end type

!  SG1 region boundaries (in units of the atomic radius) and Lebedev
!  orders, from P.M.W. Gill, B.G. Johnson, J.A. Pople,
!  Chem. Phys. Lett. 209 (1993) 506: rows are H-He, Li-Ne, Na-Ar.
!  SG1 is only defined up to Ar; heavier atoms (row 4) fall back to
!  the unpruned 194-point grid at all radii.
  real(kind=dp), parameter :: sg1rads(5,4) = reshape(&
        [0.2500d0, 0.500d0, 1.0d00, 4.50d0, 9999999.9d0,   &
         0.1667d0, 0.500d0, 0.90d0, 3.50d0, 9999999.9d0,   &
         0.1000d0, 0.400d0, 0.80d0, 2.5d0,  9999999.9d0,   &
         0.0d0,    0.0d0,   0.0d0,  9999999.9d0, 9999999.9d0], &
         shape(sg1rads))
  integer, parameter :: sg1atoms(4) =  [2,  10,  18,  137]
  integer, parameter :: sg1grids(5) =  [6, 38, 86, 194, 86]

!  SG-2 / SG-3 pruned grids: S. Dasgupta, J.M. Herbert,
!  J. Comput. Chem. 38, 869 (2017).  Radial grid: Mitani
!  double-exponential (DE2), M. Mitani, Theor. Chem. Acc. 130,
!  645 (2011), with element-specific alpha and Nr = 75 (SG-2) or
!  Nr = 99 (SG-3); R_max = 20 * Bragg-Slater radius.
!  The pruning sectors are counts of consecutive radial shells
!  (ascending radius), each integrated on the given Lebedev sphere.
!  Defined for Z in {1, 3-9, 11-17}; other elements fall back to the
!  unpruned 302/590-point grid on the standard radial grid.
  integer, parameter :: SG_NELEM = 15
  integer, parameter :: sg_elem_z(SG_NELEM) = &
        [1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17]

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
        allocate(pruned%rad_typ_z(ntyps), source=0)
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
          pruned%rad_typ_z(i) = sg_elem_z(ie)
        end do

        write(iw,'(/5X,"Standard Grid ",A," (",A,") of Dasgupta and Herbert"/&
                  &5X,40("-")/&
                  &5X,"XC functional: ",A/&
                  &5X,"NRAD  =",I8,"   (Mitani DE2 radial grid)"/&
                  &5X,"THRESH=",1P,E12.2)') &
                    pruned_name(3:3), trim(pruned_name), &
                    trim(xc_func_name), pruned%nrad, &
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
      integer :: ngr, rtid
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

!     Element-specific DE2 radial grids (SG-2/SG-3), absolute radii;
!     R_max = 20 * Bragg-Slater radius (NWChem convention)
      do i = 2, pruned%nrad_types
        call de2_radial_grid(nrad, pruned%de2_alpha(i), &
                20.0_dp*bragg_slater_radius(brsl_radii, &
                        real(pruned%rad_typ_z(i), kind=dp)), &
                molGrid%rad_pts(:,i), molGrid%rad_wts(:,i))
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

!       Number of regions used by this atom type and the pruning mode:
!       index-based sectors (nradPerRegion > 0) vs radius-based regions
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
        atomic_grid%idAtm = iat

!       Set angular Lebedev grid(s)
        do igrid = 1, ngr
!         get the unit lebedev sphere (no-op if already stored)
          call molGrid%spherical_grids%add_grid(pruned%nang(igrid, grid_id))
        end do

!       Set the radial grid of the atom.  DE2 radial grids (types >= 2)
!       store absolute radii: use a unit effective radius; the standard
!       grid (type 1) is scaled by the Bragg-Slater radius.
        rtid = molGrid%radTypeId(iat)
        atomic_grid%rad_pts = molGrid%rad_pts(:, rtid)
        atomic_grid%rad_wts = molGrid%rad_wts(:, rtid)
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
!>     x_i = i*h,  i = n_lo..n_hi
!>     r_i = exp(alpha*x_i - exp(-x_i))                        [bohr]
!>     w_i = h*(alpha + exp(-x_i))*exp(3*alpha*x_i - 3*exp(-x_i))
!>   The mesh size h is chosen so that the outermost node equals rmax:
!>   alpha*(n_hi*h) - exp(-n_hi*h) = ln(rmax), solved by Newton
!>   iteration (the left-hand side is strictly increasing in h).
!> @param[in]   nr     number of radial points
!> @param[in]   alpha  DE2 alpha parameter (element-specific)
!> @param[in]   rmax   outermost radial node, bohr
!> @param[out]  r      radial nodes, absolute bohr
!> @param[out]  w      radial weights including the r^2 Jacobian
  pure subroutine de2_radial_grid(nr, alpha, rmax, r, w)
    implicit none
    integer, intent(in) :: nr
    real(kind=dp), intent(in) :: alpha, rmax
    real(kind=dp), intent(out) :: r(:), w(:)

    integer :: i, n_lo, n_hi, iter
    real(kind=dp) :: h, x, f, fp, lnr

    if (mod(nr,2) == 1) then
      n_lo = -(nr-1)/2
    else
      n_lo = -nr/2
    end if
    n_hi = nr+n_lo-1

    lnr = log(rmax)
    h = lnr/(alpha*n_hi)
    do iter = 1, 64
      f = alpha*n_hi*h-exp(-n_hi*h)-lnr
      if (abs(f) < 1.0d-14*max(1.0_dp, abs(lnr))) exit
      fp = n_hi*(alpha+exp(-n_hi*h))
      h = h-f/fp
    end do

    do i = n_lo, n_hi
      x = i*h
      r(i-n_lo+1) = exp(alpha*x-exp(-x))
      w(i-n_lo+1) = h*(alpha+exp(-x))*exp(3.0_dp*alpha*x-3.0_dp*exp(-x))
    end do

  end subroutine de2_radial_grid

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
