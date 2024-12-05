!> @brief Module for calculating Atomic Orbital (AO) overlap between different geometries
!>
!> @details This module provides functionality to compute overlap between
!>    basis sets (Atomic Orbitals) of two different geometries. It includes
!>    calculations for AO overlap and then Molecular Orbital (MO) overlap.
!>
!> @note The matrices mol.data["OQP::xyz_old"],
!>                    mol.data["OQP::xyz"],
!>                    mol.data["OQP::VEC_MO_A_old"],
!>                    mol.data["OQP::VEC_MO_A"],
!>                    mol.data["OQP::E_MO_A_old"],
!>                    mol.data["OQP::E_MO_A"]
!>       must be defined in advance before running this program.
!>       Output AO overlap will be written to
!>                    mol.data["OQP::overlap_ao_non_orthogonal"]
!>       Output MO overlap will be written to
!>                    mol.data["OQP::overlap_mo_non_orthogonal"]
!>
!> @data Aug 2024
!>
!> @author Konstantin Komarov
!>
module get_structures_ao_overlap_mod

  implicit none

  character(len=*), parameter :: module_name = "get_structures_ao_overlap_mod"

  public get_structures_ao_overlap

contains

!> @brief C-interoperable wrapper for get_states_overlap
!>
!> @param[in] c_handle   C handle for the information structure
!>
  subroutine get_structures_ao_overlap_C(c_handle) bind(C, name="get_structures_ao_overlap")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call get_structures_ao_overlap(inf)
  end subroutine get_structures_ao_overlap_c

!> @brief Calculate AO overlap between two different geometries
!> @param[in,out] infos Information structure containing all necessary data
!>
!> This subroutine calculates the overlap between basis sets (Atomic Orbitals)
!> of two different geometries. It computes both AO and derived MO overlaps
!> and stores the results in the infos structure.
!>
!> @note Output AO overlap is written to mol.data["OQP::overlap_ao_non_orthogonal"]
!>       Output MO overlap is written to mol.data["OQP::overlap_mo_non_orthogonal"]
!>       Current geometry (xyz) is taken from mol.data["OQP::xyz"]
!>       Current MO coefficients (mo_a) are taken from mol.data["OQP::VEC_MO_A"]
!>       Current MO energies (e_a) are taken from mol.data["OQP::E_MO_A"]
!>       Old geometry (xyz_old) is taken from mol.data["OQP::xyz_old"]
!>       Old MO coefficients (mo_a_old) are taken from mol.data["OQP::VEC_MO_A_old"]
!>       Old MO energies (e_a_old) are taken from mol.data["OQP::E_MO_A_old"]
!>
  subroutine get_structures_ao_overlap(infos)
    use, intrinsic :: iso_c_binding, only: c_int32_t
    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use atomic_structure_m, only: atomic_structure
    use messages, only: show_message, with_abort
    use int1, only: basis_overlap
    use constants, only: tol_int
    use util, only: measure_time

    implicit none

    character(len=*), parameter :: subroutine_name = "get_structures_ao_overlap"

    !> Information structure containing all necessary data
    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    type(basis_set), allocatable :: basis_old
    type(atomic_structure), allocatable, target :: atoms_old
    integer :: i, nbf

    ! Tagarray definitions and data pointers
    integer(c_int32_t) :: stat
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
        OQP_XYZ_old, OQP_VEC_MO_A, OQP_E_MO_A, OQP_VEC_MO_A_old, OQP_E_MO_A_old /)
    character(len=*), parameter :: tags_alloc(*) = (/ character(len=80) :: &
        OQP_overlap_mo, OQP_overlap_ao/)
    real(kind=dp), pointer :: xyz_old(:,:), overlap_ao_out(:,:), overlap_mo_out(:,:), &
        mo_a(:,:), mo_a_old(:,:), e_a(:), e_a_old(:)

    ! Open log file
    open (unit=IW, file=infos%log_filename, position="append")

    ! Initialize basis sets and atomic structures
    allocate(basis_old, source=infos%basis)
    allocate(atoms_old, source=infos%atoms)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf

    ! Allocate and prepare data for output
    call infos%dat%erase(tags_alloc)
    stat = infos%dat%create(OQP_overlap_mo, TA_TYPE_REAL64, &
          (/ nbf, nbf /), description=OQP_overlap_mo_comment)
    stat = infos%dat%create(OQP_overlap_ao, TA_TYPE_REAL64, &
          (/ nbf, nbf /), description=OQP_overlap_ao_comment)
    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_overlap_mo, overlap_mo_out)
    call tagarray_get_data(infos%dat, OQP_overlap_ao, overlap_ao_out)

    ! Load data from python level
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_xyz_old, xyz_old)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, e_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A_old, mo_a_old)
    call tagarray_get_data(infos%dat, OQP_E_MO_A_old, e_a_old)

    ! Apply old geometry to old basis
    atoms_old%xyz = xyz_old
    basis_old%atoms => atoms_old

    call print_geo(basis_old, "Previous geometry")
    call print_geo(basis, " Current geometry")

    overlap_ao_out = 0.0_dp
    overlap_mo_out = 0.0_dp

    ! Calculate AO overlap between old and new basis sets
    call basis_overlap(overlap_ao_out, basis, basis_old, tol=log(10.0_dp)*tol_int)

    ! Normalize AO overlap
    do i = 1, nbf
      overlap_ao_out(:,i) = overlap_ao_out(:,i) * basis%bfnrm(i) * basis_old%bfnrm(:)
    end do

    ! Calculate MO overlap: <old(I)|new(J')> = transpose[Cold(AI)] . <old(A)|new(A')> . Cnew(A'J')
    call mo_overlap(overlap_mo_out, mo_a, mo_a_old, overlap_ao_out, nbf)

    ! Output results: overlap between old and current MOs
    call print_results(overlap_mo_out, e_a, e_a_old, nbf, &
                       infos%mol_prop%nelec_a, iw)

!   Print timings
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    close(iw)

  end subroutine get_structures_ao_overlap

!> @brief Calculate Molecular Orbital (MO) overlap between two geometries
!>
!> This subroutine computes the overlap between Molecular Orbitals of two different
!> geometries, given their MO coefficients and the Atomic Orbital (AO) overlap.
!>
!> @param[out] overlap_mo  Resulting MO overlap matrix
!> @param[in]  mo_a        MO coefficients obtained at the current geometry
!> @param[in]  mo_a_old    MO coefficients obtained at the old geometry
!> @param[in]  overlap_ao  AO overlap matrix between the two geometries
!> @param[in]  nbf         Number of basis functions
!>
!> @note The MO overlap is calculated as:
!>       overlap_mo = transpose(mo_a_old) . overlap_ao . mo_a
!>
  subroutine mo_overlap(overlap_mo, mo_a, mo_a_old, overlap_ao, nbf)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: overlap_mo
    real(kind=dp), intent(in), dimension(:,:) :: mo_a
    real(kind=dp), intent(in), dimension(:,:) :: mo_a_old
    real(kind=dp), intent(in), dimension(:,:) :: overlap_ao
    integer :: nbf

    real(kind=dp), allocatable, dimension(:,:) :: scr
    integer :: i

    allocate(scr(nbf,nbf), source=0.0_dp)

    call dgemm('t', 'n', nbf, nbf, nbf, &
                1.0_dp, mo_a_old, nbf, &
                        overlap_ao, nbf, &
                0.0_dp, scr, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
                1.0_dp, scr, nbf, &
                        mo_a, nbf, &
                0.0_dp, overlap_mo, nbf)
    ! Normalize MO overlap
    do i = 1, nbf
      overlap_mo(:,i) = overlap_mo(:,i) / norm2(overlap_mo(:,i))
    end do

  end subroutine mo_overlap

!> @brief Print results of MO overlap between two geometries
!>
!> This subroutine prints the results of Molecular Orbital (MO) overlap
!> calculations between two geometries, including energy comparisons
!> and overlap values.
!>
!> @param[in] overlap_mo  MO overlap matrix
!> @param[in] e           MO energies of the current geometry (in Hartree)
!> @param[in] e_old       MO energies of the old geometry (in Hartree)
!> @param[in] nbf         Number of basis functions
!> @param[in] iw          Output unit number for writing results
!>
!> @note The subroutine prints a table showing:
!>       - Corresponding overlap between old and new MOs
!>       - Energies of corresponding MOs in eV
!>       - Diagonal and maximum overlap values
!>       - Warnings for low overlap or rearranged orbitals
!>
  subroutine print_results(overlap_mo, e, e_old, nbf, na, iw)
    use precision, only: dp
    use physical_constants, only: ev2htree

    implicit none

    real(kind=dp), intent(in) :: overlap_mo(:,:)
    real(kind=dp), intent(in) :: e(:)
    real(kind=dp), intent(in) :: e_old(:)
    integer, intent(in) :: nbf, na, iw

    integer :: i, loc
    real(kind=dp) :: tmp, tmp2, tmp_abs

    write(iw, fmt="(/16X,39('=')/16X,a/16X,39('='))") &
        'Overlap MOs between geometries computed'

    write(iw,fmt='(/,x,65("-"),/,1x,a,/,x,65("-"))') &
    "Maximum Overlap (MaxO) between MOs_old(A) and MOs(B). Delta = A-B"
    write(iw,fmt='(x,a,1x,a,2x,a,2x,a,2x,a,2x,a)') &
   'A_i <- B(MaxO)', 'A_i, eV','Delta, eV', 'B_i, eV', &
   'B_i A_i Overlap','MaxO'
    do i = 1, nbf
      tmp_abs = maxval(abs(overlap_mo(:nbf,i)))
      loc = maxloc(abs(overlap_mo(:nbf,i)), dim=1)
      tmp = overlap_mo(loc,i)
      tmp2 = overlap_mo(i,i)

      write(iw, advance='no', fmt='(x,i3,3x,i4,3x,f9.3,x,f9.5,x,f9.3,x,2f12.6)') &
        i, loc, e_old(i)*ev2htree,(e_old(i)-e(i))*ev2htree, e(i)*ev2htree, &
        tmp2, tmp

      if (i == na-1) write(iw, advance='no',fmt='(2x,a)') ' HOMO'
      if (i == na) write(iw, advance='no',fmt='(2x,a)') ' LUMO'
      if (i /= loc .and. tmp_abs < 0.9_dp) then
        write(iw, fmt='(2x,a)') ' rearranged, WARNING'
      elseif (i == loc .and. tmp_abs < 0.9_dp) then
        write(iw, fmt='(2x,a)') ' WARNING'
      elseif (i /= loc .and. tmp_abs > 0.9_dp) then
        write(iw, fmt='(2x,a)') ' rearranged'
      else
        write(iw,*)
      end if
    end do
    write(iw,*)

    end subroutine print_results

!> @brief Print geometry information
!>
!> @param[in] basis   The basis class containing geometry information
!> @param[in] text    A descriptive text for the geometry output
!>
  subroutine print_geo(basis, text)
    use io_constants, only: iw
    use basis_tools, only: basis_set
    use physical_constants, only: bohr_to_angstrom

    implicit none

    type(basis_set), intent(in) :: basis
    character(len=*), intent(in) :: text
    integer :: i

    write(iw, fmt="(&
              &/26X,17('=')&
              &/26X,a&
              &/26X,17('=')&
              &/8X,'Atom     Znuc',11X,'X',14X,'Y',14X,'Z'&
              &/6X,62('-'))") text

    do i = 1, size(basis%atoms%zn(:))
       write(iw,'(7x,i4,5x,f4.1,3(x,f15.9))') &
               i, basis%atoms%zn(i), basis%atoms%xyz(1:3,i)*bohr_to_angstrom
    end do

  end subroutine print_geo

end module get_structures_ao_overlap_mod
