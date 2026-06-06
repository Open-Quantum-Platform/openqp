module guess_minao_mod

  implicit none

  character(len=*), parameter :: module_name = "guess_minao_mod"

contains

  subroutine guess_minao_C(c_handle) bind(C, name="guess_minao")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call guess_minao(inf)
  end subroutine guess_minao_C

  subroutine guess_minao(infos)
    use precision, only: dp
    use types, only: information
    use io_constants, only: IW
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use guess, only: get_ab_initio_density
    use mathlib, only: unpack_matrix
    use qmat_cache, only: get_qmat_cached
    use eigen, only: diag_symm_full
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use printing, only: print_module_info
    use parallel, only: par_env_t
    use iso_c_binding, only: c_char
    use constants, only: tol_int
    use int1, only: basis_overlap
    use minao_lut, only: minao_table_t

    implicit none

    character(len=*), parameter :: subroutine_name = "guess_minao"

    type(information), target, intent(inout) :: infos
    integer :: nbf, nbf2, nbf_min, ok, i, j, ish, iat, z, a0, m
    integer :: bar

    type(basis_set), pointer :: basis
    type(basis_set) :: min_basis
    type(minao_table_t) :: minao
    character(len=:), allocatable :: paths, basis_file, data_file
    logical :: err
    integer, parameter :: root = 0
    type(par_env_t) :: pe

    real(kind=dp), allocatable :: dmin(:,:), sco(:,:), qmat(:,:), sfull(:,:)
    real(kind=dp), allocatable :: pmat(:,:), dt(:,:), tmpmn(:,:)
    real(kind=dp), allocatable :: wrk(:,:), occ(:), cno(:,:)
    integer, allocatable :: at_ao0(:), at_nao(:)

    real(kind=dp), contiguous, pointer :: &
      Smat(:), &
      dmat_a(:), mo_a(:,:), mo_energy_a(:), &
      dmat_b(:), mo_b(:,:), mo_energy_b(:)
    character(len=*), parameter :: tags_alpha(3) = (/ character(len=80) :: &
      OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(3) = (/ character(len=80) :: &
      OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)
    character(len=*), parameter :: tags_general(2) = (/ character(len=80) :: &
      OQP_SM, OQP_hbasis_filename /)
    character(len=1,kind=c_char), contiguous, pointer :: cfn(:)

    open (unit=IW, file=infos%log_filename, position="append")
    call print_module_info('Guess_MINAO', &
        'Initial guess using projected atomic minimal-basis densities')

    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

  ! Resolve the two paths ("basisfile|datafile") passed from Python
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_hbasis_filename, cfn)
    allocate(character(ubound(cfn,1)) :: paths)
    do i = 1, ubound(cfn,1)
      paths(i:i) = cfn(i)
    end do
    bar = index(paths, '|')
    if (bar <= 0) call show_message('Guess_MINAO: bad path spec', WITH_ABORT)
    basis_file = paths(1:bar-1)
    data_file = paths(bar+1:len(paths))

  ! Load the minimal reference basis and the atomic-density table
    call min_basis%from_file(basis_file, infos%atoms, err)
    infos%control%basis_set_issue = err
    call pe%bcast(infos%control%basis_set_issue, 1)
    if (err) call show_message('Guess_MINAO: cannot read minimal basis '//trim(basis_file), WITH_ABORT)
    min_basis%atoms => infos%atoms

    call minao%load(data_file, err)
    if (err) call show_message('Guess_MINAO: cannot read MINAO data '//trim(data_file), WITH_ABORT)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nbf_min = min_basis%nbf

  ! Per-atom AO ranges in the minimal basis
    allocate(at_ao0(infos%mol_prop%natom), at_nao(infos%mol_prop%natom), source=0)
    do ish = 1, min_basis%nshell
      iat = min_basis%origin(ish)
      if (at_nao(iat) == 0) at_ao0(iat) = min_basis%ao_offset(ish)
      at_nao(iat) = at_nao(iat) + min_basis%naos(ish)
    end do

  ! Assemble block-diagonal minimal-basis density D_min
    allocate(dmin(nbf_min, nbf_min), source=0.0_dp)
    do iat = 1, infos%mol_prop%natom
      z = nint(infos%atoms%zn(iat))
      if (z < 1) cycle
      if (z > minao%zmax) &
        call show_message('Guess_MINAO: element beyond tabulated range (Z<=36)', WITH_ABORT)
      if (minao%elem(z)%nao /= at_nao(iat)) &
        call show_message('Guess_MINAO: minimal-basis size mismatch for atom', WITH_ABORT)
      a0 = at_ao0(iat)
      m = at_nao(iat)
      dmin(a0:a0+m-1, a0:a0+m-1) = minao%elem(z)%dm
    end do

  ! Cross-overlap between minimal and target basis, bfnrm-scaled (as proj_dm_newbas)
    allocate(sco(nbf_min, nbf))
    call basis_overlap(sco, basis, min_basis, tol=log(10.0d0)*tol_int)
    do i = 1, nbf
      sco(:, i) = sco(:, i) * basis%bfnrm(i) * min_basis%bfnrm(:)
    end do

  ! tagarray records
    call infos%dat%remove_records(tags_alpha)
    call infos%dat%remove_records(tags_beta)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call infos%dat%reserve_data(OQP_DM_A, TA_TYPE_REAL64, nbf2, comment=OQP_DM_A_comment)
    call infos%dat%reserve_data(OQP_E_MO_A, TA_TYPE_REAL64, nbf, comment=OQP_E_MO_A_comment)
    call infos%dat%reserve_data(OQP_VEC_MO_A, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_A_comment)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    if (infos%control%scftype >= 2) then
      call infos%dat%reserve_data(OQP_DM_B, TA_TYPE_REAL64, nbf2, comment=OQP_DM_B_comment)
      call infos%dat%reserve_data(OQP_E_MO_B, TA_TYPE_REAL64, nbf, comment=OQP_E_MO_B_comment)
      call infos%dat%reserve_data(OQP_VEC_MO_B, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_B_comment)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    end if

  ! Canonical orthogonalizer Q from matrix_invsqrt (Q = U L^{-1/2}, so that
  ! Q^T S Q = I and Q Q^T = S^{-1}), plus the full overlap S (unpacked).
    allocate(qmat(nbf,nbf), sfull(nbf,nbf))
    call get_qmat_cached(infos, smat, qmat, nbf)
    call unpack_matrix(smat, sfull)

  ! Projector pmat(target,min) = S^{-1} Scross(target,min) = Q (Q^T Scross),
  ! with Scross(target,min) = sco^T (sco is stored (min,target)). Then form the
  ! projected density D_t = pmat * D_min * pmat^T in the target basis.
    allocate(pmat(nbf, nbf_min), tmpmn(nbf, nbf_min), dt(nbf, nbf), wrk(nbf, nbf))
    ! tmpmn = Q^T * sco^T   (= Q^T Scross)
    call dgemm('t','t', nbf, nbf_min, nbf, 1.0_dp, qmat, nbf, sco, nbf_min, 0.0_dp, tmpmn, nbf)
    ! pmat = Q * tmpmn
    call dgemm('n','n', nbf, nbf_min, nbf, 1.0_dp, qmat, nbf, tmpmn, nbf, 0.0_dp, pmat, nbf)
    ! D_t = pmat * D_min * pmat^T
    call dgemm('n','n', nbf, nbf_min, nbf_min, 1.0_dp, pmat, nbf, dmin, nbf_min, 0.0_dp, tmpmn, nbf)
    call dgemm('n','t', nbf, nbf, nbf_min, 1.0_dp, tmpmn, nbf, pmat, nbf, 0.0_dp, dt, nbf)

  ! Natural orbitals in the orthonormal (Q) basis: diagonalize
  ! M = Q^T (S D_t S) Q = W n W^T. Eigenvalues n are occupations; the AO-basis
  ! natural orbitals are C = Q W and reconstruct D_t exactly.
    ! wrk = S * D_t
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, sfull, nbf, dt, nbf, 0.0_dp, wrk, nbf)
    ! dt <- (S D_t) * S
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, wrk, nbf, sfull, nbf, 0.0_dp, dt, nbf)
    ! wrk = Q^T * (S D_t S)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, qmat, nbf, dt, nbf, 0.0_dp, wrk, nbf)
    ! dt = wrk * Q  = Q^T S D_t S Q
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, wrk, nbf, qmat, nbf, 0.0_dp, dt, nbf)
    allocate(occ(nbf), cno(nbf,nbf))
    call diag_symm_full(1, nbf, dt, nbf, occ, ok)
    if (ok /= 0) call show_message('Guess_MINAO: NO diagonalization failed', WITH_ABORT)
    ! NOs in AO basis: C = Q * W   (W = eigenvectors in dt)
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, qmat, nbf, dt, nbf, 0.0_dp, cno, nbf)

  ! diag_symm_full returns ascending occupations; reverse so highest-occupied first
    do j = 1, nbf
      mo_a(:, j) = cno(:, nbf-j+1)
      mo_energy_a(j) = -occ(nbf-j+1)
    end do
    if (infos%control%scftype >= 2) then
      mo_b = mo_a
      mo_energy_b = mo_energy_a
    end if

  ! Build density from aufbau-occupied natural orbitals
    if (infos%control%scftype == 1) then
      call get_ab_initio_density(Dmat_A, MO_A, infos=infos, basis=basis)
    else
      call get_ab_initio_density(Dmat_A, MO_A, Dmat_B, MO_B, infos, basis)
    end if

    call minao%clean()

    write (IW, '(/1x,a/)') '...... End Of Initial Orbital Guess ......'
    call measure_time(print_total=1, log_unit=iw)
    close(IW)

  end subroutine guess_minao

end module guess_minao_mod
