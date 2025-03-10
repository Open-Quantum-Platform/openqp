module basis_projection_mod

  implicit none

  character(len=*), parameter :: module_name = "basis_projection_mod"

contains

  subroutine proj_dm_newbas_C(c_handle) bind(C, name="proj_dm_newbas")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call proj_dm_newbas(inf)
  end subroutine proj_dm_newbas_C

  !----------------------------------------------------------------------
  ! Main subroutine to perform the MO projection between basis sets.
  ! Projects molecular orbitals (MO) and density matrices (DM) from a
  ! primary basis set to an alternative (initial) basis set.
  ! Input data:
  !    -  OQP::VEC_MO_A are assumed to be present.
  !       OQP::VEC_MO_A
  !       OQP::DM_A
  !       OQP::VEC_MO_B
  !       OQP::DM_B
  ! Output data:
  !    - The projected density matrices and MOs are written into:
  !       OQP::VEC_MO_A
  !       OQP::DM_A
  !       OQP::VEC_MO_B
  !       OQP::DM_B
  !----------------------------------------------------------------------

  subroutine proj_dm_newbas(infos)
    use precision, only: dp
    use types, only: information
    use io_constants, only: IW
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use mathlib, only: matrix_invsqrt
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use printing, only: print_module_info
    use constants, only: tol_int
    use int1, only: basis_overlap
    use oqp_tagarray_driver
    use iso_c_binding, only: c_char
    use parallel, only: par_env_t
    use guess, only: get_ab_initio_density, corresponding_orbital_projection
    use huckel, only: orthogonalize_orbitals
    use messages, only: show_message, with_abort

    implicit none

    character(len=*), parameter :: subroutine_name = "proj_dm_newbas"

    type(information), target, intent(inout) :: infos
    integer :: i, j, nbf, nbf2, nbf_alt, nbf2_alt, nat, nact, ndoc, nproj, l0

    type(basis_set), pointer :: basis
    type(basis_set), pointer :: alt_basis
    character(len=:), allocatable :: basis_file
    logical :: err
    integer , parameter :: root = 0
    type(par_env_t) :: pe
    real(kind=dp), allocatable :: sco(:,:)
    real(kind=dp), contiguous, pointer :: &
      Smat(:), q(:,:), &
      dmat_a(:), mo_a(:,:), mo_energy_a(:), &
      dmat_b(:), mo_b(:,:), mo_energy_b(:)
    real(kind=dp), contiguous, pointer :: &
      Smat_alt(:), &
      dmat_a_alt(:), mo_a_alt(:,:), mo_energy_a_alt(:), &
      dmat_b_alt(:), mo_b_alt(:,:), mo_energy_b_alt(:)
  ! tagarray
    character(len=*), parameter :: tags_alpha(3) = (/ character(len=80) :: &
      OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_alpha_tmp(3) = (/ character(len=80) :: &
      "OQP::DM_A_tmp", "OQP::E_MO_A_tmp", "OQP::VEC_MO_A_tmp" /)
    character(len=*), parameter :: tags_beta(3) = (/ character(len=80) :: &
      OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)
    character(len=*), parameter :: tags_beta_tmp(3) = (/ character(len=80) :: &
      "OQP::DM_B_tmp", "OQP::E_MO_B_tmp", "OQP::VEC_MO_B_tmp" /)
    character(len=*), parameter :: tags_general(1) = (/ character(len=80) :: &
      OQP_SM /)
    character(len=1,kind=c_char), contiguous, pointer :: basis_filename(:)

  ! Files open
  ! 1. XYZ: Read : Geometric data, ATOMS
  ! 3. LOG: Read Write: Main output file
  ! 
    open (unit=IW, file=infos%log_filename, position="append")
    call print_module_info('Basis Projection', "Projecting MOs and DMs from" &
                              // new_line('') // &
                             "initial to primary basis")

  ! Readings
  ! load basis set
    basis => infos%basis
    alt_basis => infos%alt_basis
    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    alt_basis%atoms => infos%atoms
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 =nbf*(nbf+1)/2
    nat = infos%mol_prop%natom

    nbf_alt = alt_basis%nbf
    nbf2_alt =nbf_alt*(nbf_alt+1)/2

    if (nbf_alt > nbf) then
      call show_message("Warning: The initial basis set ("//trim(adjustl(itoa(nbf_alt)))//" functions) " // &
                  "exceeds the primary basis ("//trim(adjustl(itoa(nbf)))//"). " // &
                  "Please select a smaller initial basis set.", WITH_ABORT)
    end if
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    allocate(sco(nbf_alt,nbf), &
            mo_a(nbf,nbf), &
            mo_b(nbf,nbf), &
            q(nbf,nbf), &
            Dmat_a(nbf2), &
            Dmat_b(nbf2))

    ! allocate alpha
    call infos%dat%reserve_data(OQP_DM_A, TA_TYPE_REAL64, nbf2_alt, comment=OQP_DM_A_comment)
    call infos%dat%reserve_data(OQP_E_MO_A, TA_TYPE_REAL64, nbf_alt, comment=OQP_E_MO_A_comment)
    call infos%dat%reserve_data(OQP_VEC_MO_A, TA_TYPE_REAL64, nbf_alt*nbf_alt, (/ nbf_alt, nbf_alt /), comment=OQP_VEC_MO_A_comment)
    ! load alpha data
    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a_alt)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a_alt)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a_alt)
    call infos%dat%reserve_data(OQP_DM_B, TA_TYPE_REAL64, nbf2_alt, comment=OQP_DM_B_comment)
    call infos%dat%reserve_data(OQP_E_MO_B, TA_TYPE_REAL64, nbf_alt, comment=OQP_E_MO_B_comment)
    call infos%dat%reserve_data(OQP_VEC_MO_B, TA_TYPE_REAL64, nbf_alt*nbf_alt, (/ nbf_alt, nbf_alt /), comment=OQP_VEC_MO_B_comment)
!    ! load beta
    call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b_alt)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b_alt)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b_alt)
    ! clean data
    call infos%dat%remove_records(tags_alpha_tmp)
    call infos%dat%remove_records(tags_beta_tmp)
    ! allocate alpha_tmp
    call infos%dat%reserve_data("OQP::DM_A_tmp", TA_TYPE_REAL64, nbf2, comment=OQP_DM_A_comment)
    call infos%dat%reserve_data("OQP::E_MO_A_tmp", TA_TYPE_REAL64, nbf, comment=OQP_E_MO_A_comment)
    call infos%dat%reserve_data("OQP::VEC_MO_A_tmp", TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_A_comment)
    ! load alpha_tmp data
    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, "OQP::DM_A_tmp", dmat_a)
    call tagarray_get_data(infos%dat, "OQP::E_MO_A_tmp", mo_energy_a)
    call tagarray_get_data(infos%dat, "OQP::VEC_MO_A_tmp", mo_a)
    ! allocate beta_tmp
    call infos%dat%reserve_data("OQP::DM_B_tmp", TA_TYPE_REAL64, nbf2, comment=OQP_DM_B_comment)
    call infos%dat%reserve_data("OQP::E_MO_B_tmp", TA_TYPE_REAL64, nbf, comment=OQP_E_MO_B_comment)
    call infos%dat%reserve_data("OQP::VEC_MO_B_tmp", TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_B_comment)
    ! load beta_tmp
    call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, "OQP::DM_B_tmp", dmat_b)
    call tagarray_get_data(infos%dat, "OQP::E_MO_B_tmp", mo_energy_b)
    call tagarray_get_data(infos%dat, "OQP::VEC_MO_B_tmp", mo_b)
!
    Dmat_b = 0_dp
    Dmat_a = 0_dp
    mo_b = 0_dp
    mo_a = 0_dp
    mo_energy_a = 0_dp
    mo_energy_b = 0_dp

    call basis_overlap(sco, infos%basis, infos%alt_basis, tol=log(10.0d0)*tol_int)
    do i =1, nbf
      sco(:,i) = sco(:,i)*basis%bfnrm(i) * alt_basis%bfnrm(:)
    end do

    if (infos%control%scftype == 1) then
      ndoc = infos%mol_prop%nelec/2
      nact = 0
    else if (infos%control%scftype >= 2) then
      ndoc = infos%mol_prop%nelec_b
      nact = infos%mol_prop%nelec_a-infos%mol_prop%nelec_b
    end if

    nproj = nbf_alt
    l0 = nbf
    call matrix_invsqrt(smat, q, nbf, qrnk=l0)
    mo_a(1:nbf,1:nbf) = q(1:nbf,1:nbf)

    call corresponding_orbital_projection(mo_a_alt, sco, mo_a, ndoc, nact, nproj, nbf, nbf_alt, l0)
    call orthogonalize_orbitals(q, smat, mo_a, nproj, l0, nbf, nbf)

    if (infos%control%scftype >= 2) then
      call corresponding_orbital_projection(mo_b_alt, sco, mo_b, ndoc, nact, nproj, nbf, nbf_alt, l0)
      call orthogonalize_orbitals(q, smat, mo_b, nproj, l0, nbf, nbf)
    else
      Dmat_B = Dmat_A
    end if

  ! Calculate Density Matrix
    if (pe%rank == root) then
  ! RHF
      if (infos%control%scftype == 1) then
        call get_ab_initio_density(Dmat_A, MO_A, infos=infos, basis=basis)
  ! ROHF/UHF
      else
        call get_ab_initio_density(Dmat_A, MO_A, Dmat_B, MO_B, infos, infos%basis)
      endif
    endif
    ! Broadcast MO and density matrices to all processes
    call pe%bcast(MO_A, nbf*nbf)
    if (infos%control%scftype >= 2) then
      call pe%bcast(MO_B, nbf*nbf)
    endif
    ! Broadcast the density matrices to all processes
    if (infos%control%scftype == 1) then
      call pe%bcast(Dmat_A, nbf2)
    else
      call pe%bcast(Dmat_A, nbf2)
      call pe%bcast(Dmat_B, nbf2)
    endif
    call pe%barrier()
    write(IW, '(/,a,/)') "...... Completed Basis Projection Computation ......"
    call measure_time(print_total=1, log_unit=iw)
    close(iw)
  end subroutine proj_dm_newbas

  pure function itoa(num) result(str)
    integer, intent(in) :: num
    character(len=16) :: str
    write(str, '(I0)') num
  end function itoa

end module basis_projection_mod
