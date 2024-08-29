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

  integer, parameter :: max_number_of_grid_based_DFT_grids = 10
  integer, parameter :: max_number_of_atom_grid_types_for_DFT_grids = 10

  type dft_grid_pruned_t
    integer :: nrad, nang(max_number_of_grid_based_DFT_grids)
    real(kind=dp) :: radii(max_number_of_grid_based_DFT_grids, &
                               max_number_of_atom_grid_types_for_DFT_grids)
    integer :: ngrids = 1
    integer, allocatable :: rad_id(:)
  end type

  real(kind=dp), parameter :: sg1rads(5,4) = reshape(&
        [0.2500d0, 0.500d0, 1.0d00, 4.50d0, 9999999.9d0,   &
         0.1667d0, 0.500d0, 0.90d0, 3.50d0, 9999999.9d0,   &
         0.1000d0, 0.400d0, 0.80d0, 2.5d0,  9999999.9d0,   &
         1.0d-30, 999999.9d0, 999999.9d0, 999999.9d0, 999999.9d0], &
         shape(sg1rads))
  integer, parameter :: sg1atoms(4) =  [2,  10,  18,  137]
  integer, parameter :: sg1grids(5) =  [6, 38, 86, 194, 86]

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

  subroutine dft_initialize(infos, basis, molGrid, orbitals_cutoff, verbose)
    use basis_tools, only: basis_set
    use types, only: information

    implicit none

    type(basis_set), intent(inout) :: basis
    type(information), intent(inout) :: infos
    type(dft_grid_t), intent(inout) :: molGrid
    real(kind=dp), optional :: orbitals_cutoff
    logical, optional :: verbose

    real(kind=dp) :: logtol
    type(dft_grid_pruned_t) :: pruned

!   Setup sreening parameters
    logtol = -log(1.0e-10_dp)
    if (present(orbitals_cutoff)) logtol = -log(orbitals_cutoff)
    call basis%set_screening(logtol)

!   Set grid DFT options
    call dft_set_options(infos, pruned)

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

  subroutine dft_set_options(infos, pruned)
    use iso_c_binding, only: c_null_char
    use messages, only: show_message, WITH_ABORT
    use strings, only: c_f_char
    use types, only: information
    use libxc, only: libxc_input

    implicit none

    type(information), intent(inout) :: infos
    type(dft_grid_pruned_t), intent(inout) :: pruned
    type(saved_HF_info) :: saved_hf

    integer :: iatm, nrad
    character(len=20) :: xc_func_name
    integer :: nat, i, slen, ntyps
    character(:), allocatable :: pruned_name

!   Default radial/angular grid is 96/302 for LDA/GGA.
    nrad     = infos%dft%grid_rad_size
    pruned%nang(1)  = infos%dft%grid_ang_size
    nat = ubound(infos%atoms%zn, 1)
    allocate(pruned%rad_id(nat), source=1)

    xc_func_name = c_f_char(infos%dft%xc_functional_name)

    if (.not. infos%dft%grid_pruned) then
      pruned%ngrids = 1
      pruned%radii(1,1) = 1.0d+30

      write(iw,'(/5X,"Lebedev grid-based DFT options"/&
                &5X,30("-")/&
                &5X,"XC functional: ",A/&
                &5X,"NRAD  =",I8,5X,"NLEB  =",I8/&
                &5X,"THRESH=",1P,E12.2)') &
                  trim(xc_func_name), &
                  nrad, pruned%nang(1), &
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
        pruned%radii(1:pruned%ngrids,1:ntyps) = sg1rads
        pruned%nang(:pruned%ngrids) = sg1grids(1:pruned%ngrids)
        do i = pruned%ngrids+1, ubound(pruned%radii, 1)
          pruned%radii(i,1:ntyps) = sg1rads(pruned%ngrids,1:ntyps)
        end do
        do iatm = 1, nat
          do i = 1, ntyps
            if (int(infos%atoms%zn(iatm))<sg1atoms(i)) exit
          end do
          pruned%rad_id(iatm) = i
        end do
      case default
        call show_message('Unknown pruned grid name', WITH_ABORT)
      end select

      write(iw,'(/5X,"Standard Grid 1 (SG1)"/&
                &5X,21("-")/&
                &5X,"XC functional: ",A/&
                &5X,"THRESH=",1P,E12.2)') &
                  trim(xc_func_name), &
                  infos%dft%grid_density_cutoff
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
    else
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
      max_ang_pts = maxval(pruned%nang(:pruned%ngrids))
      maxpt_per_atom = infos%dft%grid_rad_size*max_ang_pts
      nrad = infos%dft%grid_rad_size

      allocate(&
        txyz(max_ang_pts*3), &
        twght(max_ang_pts), &
        rij(nat,nat), &
        aij(nat,nat), &
        bsrad(nat), &
        source=0.0d0)

!     Init storage for the grid
      call molGrid%reset(nat, maxpt_per_atom, nRad)

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

!     Set up radial grid
      call get_radial_grid(molGrid%rad_pts, molGrid%rad_wts, &
              nrad, infos%dft%rad_grid_type)

!     Pre-compute atomic distances
      call get_atomic_distances(infos%atoms%xyz, rij)

!     Tag non-real atoms present in the system
      molGrid%dummyAtom(:nat) = bsrad(:nat) == 0.0_dp

!     Find nearest neighbours for all atoms
      call molGrid%find_neighbours(rij, partFunType=infos%dft%dft_partfun)

!     Set radial grid
      atomic_grid%rad_pts = molGrid%rad_pts
      atomic_grid%rad_wts = molGrid%rad_wts

!     Compute atomic grids for each atom
      do iat = 1, nat
        grid_id = pruned%rad_id(iat)

        atomic_grid%sph_npts = pruned%nang(1:pruned%ngrids)
        atomic_grid%sph_radii = pruned%radii(:pruned%ngrids, grid_id)
        atomic_grid%idAtm = iat

!       Set angular Lebedev grid(s)
        do igrid = 1, pruned%ngrids
!         get the unit lebedev sphere
          call molGrid%spherical_grids%add_grid(pruned%nang(igrid))
        end do

        atomic_grid%rAtm = bsrad(iat)
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

  end subroutine

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
