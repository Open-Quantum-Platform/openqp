module mrsf_nac_interchange_mod

  use precision, only: dp
  use oqp_linalg
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public :: mrsf_nac_rohf_zvector
  public :: mrsf_nac_rohf_zvector_batch
  public :: mrsf_nac_rohf_solve
  public :: mrsf_nac_rohf_pair_overlap
  public :: mrsf_nac_pair_accumulator_init
  public :: mrsf_nac_pair_accumulate
  public :: mrsf_nac_pair_accumulate_antisym
  public :: mrsf_nac_pair_finalize
  public :: mrsf_nac_rohf_hf_adjoint
  public :: mrsf_nac_rohf_hf_adjoint_batch
  public :: mrsf_nac_xc_adjoint
  public :: mrsf_nac_xc_adjoint_batch

  ! Converged unordered-pair ROHF adjoints from the two preceding geometries.
  ! Both saved vectors for a pair are expressed in the most recent saved MO
  ! frame.  At the next geometry the overlap-defined block-Procrustes map
  ! transports them together before either a one-step or linear predictor is
  ! supplied to MINRES.  The converged residual criterion remains unchanged.
  real(kind=dp), allocatable, save :: nac_z_recent(:,:), nac_z_earlier(:,:)
  real(kind=dp), allocatable, save :: nac_z_geometry(:,:,:)
  logical, allocatable, save :: nac_z_have_recent(:), nac_z_have_earlier(:), &
                                nac_z_have_geometry(:)
  integer, allocatable, save :: nac_z_steps_since_exact(:), &
                                nac_z_exact_count(:)
  integer, save :: nac_z_ltot = 0, nac_z_npair = 0, nac_z_natom = 0
  integer, save :: nac_z_nbf = 0, nac_z_nocca = 0, nac_z_noccb = 0

contains

!###############################################################################

  subroutine nac_z_cache_reset()
    if (allocated(nac_z_recent)) deallocate(nac_z_recent)
    if (allocated(nac_z_earlier)) deallocate(nac_z_earlier)
    if (allocated(nac_z_geometry)) deallocate(nac_z_geometry)
    if (allocated(nac_z_have_recent)) deallocate(nac_z_have_recent)
    if (allocated(nac_z_have_earlier)) deallocate(nac_z_have_earlier)
    if (allocated(nac_z_have_geometry)) deallocate(nac_z_have_geometry)
    if (allocated(nac_z_steps_since_exact)) deallocate(nac_z_steps_since_exact)
    if (allocated(nac_z_exact_count)) deallocate(nac_z_exact_count)
    nac_z_ltot = 0; nac_z_npair = 0; nac_z_natom = 0
    nac_z_nbf = 0; nac_z_nocca = 0; nac_z_noccb = 0
  end subroutine nac_z_cache_reset

!###############################################################################

  subroutine nac_z_cache_prepare(ltot, npair, natom, nbf, nocca, noccb)
    integer, intent(in) :: ltot, npair, natom, nbf, nocca, noccb
    logical :: dimensions_match

    dimensions_match = allocated(nac_z_recent) .and. &
      nac_z_ltot == ltot .and. nac_z_npair == npair .and. &
      nac_z_natom == natom .and. nac_z_nbf == nbf .and. &
      nac_z_nocca == nocca .and. nac_z_noccb == noccb
    if (dimensions_match) return

    call nac_z_cache_reset()
    allocate(nac_z_recent(ltot,npair), nac_z_earlier(ltot,npair), &
             nac_z_geometry(3,natom,npair), &
             nac_z_have_recent(npair), nac_z_have_earlier(npair), &
             nac_z_have_geometry(npair), nac_z_steps_since_exact(npair), &
             nac_z_exact_count(npair))
    nac_z_recent = 0.0_dp; nac_z_earlier = 0.0_dp
    nac_z_geometry = 0.0_dp
    nac_z_have_recent = .false.; nac_z_have_earlier = .false.
    nac_z_have_geometry = .false.
    nac_z_steps_since_exact = 0
    nac_z_exact_count = 0
    nac_z_ltot = ltot; nac_z_npair = npair; nac_z_natom = natom
    nac_z_nbf = nbf; nac_z_nocca = nocca; nac_z_noccb = noccb
  end subroutine nac_z_cache_prepare

!###############################################################################

  subroutine nac_z_nearest_orthogonal(overlap, rotation, singular_min, ok)
    real(kind=dp), intent(in) :: overlap(:,:)
    real(kind=dp), intent(out) :: rotation(:,:)
    real(kind=dp), intent(out) :: singular_min
    logical, intent(out) :: ok
    real(kind=dp), allocatable :: a(:,:), u(:,:), vt(:,:), singular(:), work(:)
    real(kind=dp) :: work_query(1)
    integer :: n, info, lwork

    n = size(overlap,1)
    ok = n == size(overlap,2) .and. size(rotation,1) == n .and. &
         size(rotation,2) == n
    if (.not. ok) return
    if (n == 0) then
      singular_min = 1.0_dp
      return
    end if
    allocate(a(n,n), u(n,n), vt(n,n), singular(n))
    a = overlap
    call dgesvd('A','A',n,n,a,n,singular,u,n,vt,n,work_query,-1,info)
    if (info /= 0) then
      ok = .false.; singular_min = 0.0_dp
      deallocate(a,u,vt,singular)
      return
    end if
    lwork = max(1,int(work_query(1)))
    allocate(work(lwork))
    a = overlap
    call dgesvd('A','A',n,n,a,n,singular,u,n,vt,n,work,lwork,info)
    ok = info == 0 .and. all(ieee_is_finite(singular))
    if (ok) then
      singular_min = minval(singular)
      call dgemm('n','n',n,n,n,1.0_dp,u,n,vt,n,0.0_dp,rotation,n)
    else
      singular_min = 0.0_dp
      rotation = 0.0_dp
    end if
    deallocate(a,u,vt,singular,work)
  end subroutine nac_z_nearest_orthogonal

!###############################################################################

  subroutine nac_z_transport_vector(old_vector, new_vector, q_closed, q_open, &
                                    q_virtual, nocca, noccb, nbf)
    real(kind=dp), intent(in) :: old_vector(:)
    real(kind=dp), intent(out) :: new_vector(:)
    real(kind=dp), intent(in) :: q_closed(:,:), q_open(:,:), q_virtual(:,:)
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), allocatable :: co(:,:), cv(:,:), ov(:,:), transformed(:,:)
    integer :: offset, nvira, iv, io, ic, k

    offset = nocca - noccb
    nvira = nbf - nocca
    allocate(co(offset,noccb), cv(nvira,noccb), ov(nvira,offset))
    k = 0
    do io = 1, offset
      do ic = 1, noccb
        k = k + 1; co(io,ic) = old_vector(k)
      end do
    end do
    do iv = 1, nvira
      do ic = 1, noccb
        k = k + 1; cv(iv,ic) = old_vector(k)
      end do
    end do
    do iv = 1, nvira
      do io = 1, offset
        k = k + 1; ov(iv,io) = old_vector(k)
      end do
    end do

    if (offset > 0 .and. noccb > 0) then
      transformed = matmul(q_open, matmul(co, transpose(q_closed)))
      co = transformed
    end if
    if (nvira > 0 .and. noccb > 0) then
      transformed = matmul(q_virtual, matmul(cv, transpose(q_closed)))
      cv = transformed
    end if
    if (nvira > 0 .and. offset > 0) then
      transformed = matmul(q_virtual, matmul(ov, transpose(q_open)))
      ov = transformed
    end if

    k = 0
    do io = 1, offset
      do ic = 1, noccb
        k = k + 1; new_vector(k) = co(io,ic)
      end do
    end do
    do iv = 1, nvira
      do ic = 1, noccb
        k = k + 1; new_vector(k) = cv(iv,ic)
      end do
    end do
    do iv = 1, nvira
      do io = 1, offset
        k = k + 1; new_vector(k) = ov(iv,io)
      end do
    end do
    deallocate(co,cv,ov)
  end subroutine nac_z_transport_vector

!###############################################################################

  subroutine nac_z_transport_matrices(infos, q_closed, q_open, q_virtual, &
                                      singular_min, valid)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_overlap_mo
    use, intrinsic :: iso_c_binding, only: c_int32_t
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: q_closed(:,:), q_open(:,:), q_virtual(:,:)
    real(kind=dp), intent(out) :: singular_min
    logical, intent(out) :: valid
    character(len=80) :: tags_overlap(1), tags_state(1)
    integer(c_int32_t) :: tag_id
    real(kind=dp), contiguous, pointer :: overlap(:,:), state_overlap(:)
    real(kind=dp) :: sc, so, sv, state_min, threshold
    character(len=32) :: env
    integer :: ios, nbf, nocca, noccb
    logical :: okc, oko, okv

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    tags_overlap(1) = OQP_overlap_mo
    valid = infos%dat%contains(tags_overlap, tag_id)
    if (.not. valid) then
      singular_min = 0.0_dp
      return
    end if
    call tagarray_get_data(infos%dat, OQP_overlap_mo, overlap)
    if (size(overlap,1) /= nbf .or. size(overlap,2) /= nbf) then
      valid = .false.; singular_min = 0.0_dp
      return
    end if

    ! overlap stores <old|current>.  Its transpose is the orthogonal map from
    ! the old orbital coordinates into the current coordinates, which is the
    ! direction needed to transport a saved response vector.
    call nac_z_nearest_orthogonal(transpose(overlap(1:noccb,1:noccb)), &
                                  q_closed,sc,okc)
    call nac_z_nearest_orthogonal( &
      transpose(overlap(noccb+1:nocca,noccb+1:nocca)), q_open,so,oko)
    call nac_z_nearest_orthogonal( &
      transpose(overlap(nocca+1:nbf,nocca+1:nbf)), q_virtual,sv,okv)
    singular_min = min(sc,min(so,sv))
    threshold = 0.5_dp
    env = ''
    call get_environment_variable('OQP_MRSF_NAC_ZV_OVERLAP_MIN',env)
    if (len_trim(env) > 0) then
      read(env,*,iostat=ios) threshold
      if (ios /= 0) threshold = 0.5_dp
    end if
    valid = okc .and. oko .and. okv .and. &
            singular_min >= max(0.0_dp,min(1.0_dp,threshold))

    tags_state(1) = 'OQP::state_tracking_overlap'
    if (valid .and. infos%dat%contains(tags_state,tag_id)) then
      call tagarray_get_data(infos%dat,tags_state(1),state_overlap)
      if (size(state_overlap) > 0) then
        state_min = minval(abs(state_overlap))
        valid = state_min >= max(0.0_dp,min(1.0_dp,threshold))
      end if
    end if
  end subroutine nac_z_transport_matrices

!###############################################################################

  subroutine mrsf_nac_rohf_zvector_C(c_handle) &
      bind(C, name="mrsf_nac_rohf_zvector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_rohf_zvector(inf)
  end subroutine mrsf_nac_rohf_zvector_C

!###############################################################################

!> ABI-compatible alias for clients built against the provisional API name.
!> Compatibility clients needing the canonical one-RHS entry should use
!> mrsf_nac_rohf_zvector.  The resident production driver instead calls the
!> internal batched routine; neither C entry solves a 3N block of forward CPHF
!> equations.
  subroutine mrsf_nac_rohf_solve_C(c_handle) &
      bind(C, name="mrsf_nac_rohf_solve")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_rohf_solve(inf)
  end subroutine mrsf_nac_rohf_solve_C

!###############################################################################

  subroutine mrsf_nac_rohf_pair_overlap_C(c_handle) &
      bind(C, name="mrsf_nac_rohf_pair_overlap")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_rohf_pair_overlap(inf)
  end subroutine mrsf_nac_rohf_pair_overlap_C

!###############################################################################

  subroutine mrsf_nac_pair_accumulator_init_C(c_handle) &
      bind(C, name="mrsf_nac_pair_accumulator_init")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_pair_accumulator_init(inf)
  end subroutine mrsf_nac_pair_accumulator_init_C

!###############################################################################

  subroutine mrsf_nac_pair_accumulate_C(c_handle, istate, jstate) &
      bind(C, name="mrsf_nac_pair_accumulate")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use, intrinsic :: iso_c_binding, only: c_int32_t

    type(oqp_handle_t) :: c_handle
    integer(c_int32_t), intent(in), value :: istate, jstate
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_pair_accumulate(inf, int(istate), int(jstate))
  end subroutine mrsf_nac_pair_accumulate_C

!###############################################################################

  subroutine mrsf_nac_pair_finalize_C(c_handle) &
      bind(C, name="mrsf_nac_pair_finalize")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_pair_finalize(inf)
  end subroutine mrsf_nac_pair_finalize_C

!> Solve one native ROHF adjoint Z-vector equation for one ordered state pair.
!> This public compatibility entry retains the one-RHS TagArray contract.  The
!> resident production driver uses mrsf_nac_rohf_zvector_batch to solve the
!> antisymmetric unordered-pair sources in one shared CPHF context.  Neither
!> route solves the 3N forward nuclear CPHF equations; subsequent resident
!> Fortran adjoint contractions evaluate z^T B^R directly.  The operator and
!> coordinates are exactly those used by the nuclear response in
!> hf_hessian_rohf, with no legacy sfrolhs/sfropcal metric.
!>
!> Input : OQP::nac_rohf_rhs       (ltot)
!> Output: OQP::nac_rohf_solution  (ltot)
  subroutine mrsf_nac_rohf_zvector(infos)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, TA_TYPE_REAL64
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_rhs = "OQP::nac_rohf_rhs"
    character(len=*), parameter :: tag_solution = "OQP::nac_rohf_solution"
    type(information), target, intent(inout) :: infos
    real(kind=dp), contiguous, pointer :: rhs_in(:)
    real(kind=dp), pointer :: solution_out(:)
    real(kind=dp), allocatable :: rhs(:,:), solution(:,:)
    integer :: nbf, nocca, noccb, nvira, offset, ltot

    if (infos%control%scftype /= 3) then
      call show_message('mrsf_nac_rohf_zvector requires an ROHF/ROKS reference.', &
                        WITH_ABORT)
    end if

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira

    call tagarray_get_data(infos%dat, tag_rhs, rhs_in)
    if (size(rhs_in) /= ltot) then
      call show_message('OQP::nac_rohf_rhs has the wrong ROHF rotation dimension.', &
                        WITH_ABORT)
    end if

    allocate(rhs(ltot,1), solution(ltot,1))
    rhs(:,1) = rhs_in
    call mrsf_nac_rohf_zvector_batch(infos, rhs, solution)

    call infos%dat%erase((/ character(len=80) :: tag_solution /))
    call tagarray_reserve_data(infos%dat, tag_solution, TA_TYPE_REAL64, ltot, (/ ltot /), &
         comment='one-pair ROHF NAC adjoint Z-vector; no 3N forward CPHF')
    call tagarray_get_data(infos%dat, tag_solution, solution_out)
    solution_out = solution(:,1)

    deallocate(rhs, solution)
  end subroutine mrsf_nac_rohf_zvector

!###############################################################################

!> Solve a batch of native ROHF adjoint equations in one shared CPHF context.
!>
!> The resident driver passes one antisymmetric source for every unordered
!> state pair,
!>
!>   rhs^-_IJ = 1/2 (rhs_IJ - rhs_JI),  I < J.
!>
!> Since every pair has the same symmetric ROHF Hessian H, linearity gives
!> H^-1 rhs^-_IJ = 1/2 (z_IJ - z_JI).  cphf_solve_rohf constructs the Fock,
!> preconditioner and XC context once, then certifies the true unpreconditioned
!> residual of every returned MINRES solution.  The public residual convention
!> is squared, so tol=1e-20 requests ||H z - rhs||_2 <= 1e-10.
  subroutine mrsf_nac_rohf_zvector_batch(infos, rhs, solution, pair_offset, &
      predictor, predictor_available, predictor_accepted, &
      initial_residual_out, final_residual_out, iterations_out)
    use types, only: information
    use io_constants, only: iw
    use cphf_mod, only: cphf_solve_rohf
    use messages, only: show_message, WITH_ABORT, WITHOUT_ABORT

    implicit none

    real(kind=dp), parameter :: fallback_rel_norm = &
      10.0_dp*sqrt(epsilon(1.0_dp))
    real(kind=dp), parameter :: fallback_rel_sq = &
      fallback_rel_norm*fallback_rel_norm
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: rhs(:,:)
    real(kind=dp), intent(out) :: solution(:,:)
    integer, intent(in), optional :: pair_offset
    real(kind=dp), intent(out), optional :: predictor(:,:)
    logical, intent(out), optional :: predictor_available(:), &
                                      predictor_accepted(:)
    real(kind=dp), intent(out), optional :: initial_residual_out(:), &
                                            final_residual_out(:)
    integer, intent(out), optional :: iterations_out(:)
    real(kind=dp), allocatable :: residual(:), initial_residual(:), &
      guess(:,:), recent_current(:), earlier_current(:), q_closed(:,:), &
      q_open(:,:), q_virtual(:,:)
    real(kind=dp) :: rhs_scale_sq, scaled_residual_sq, eta, max_disp, &
      displacement, singular_min, z_relative_error
    logical, allocatable :: converged(:), guess_available(:), guess_accepted(:)
    logical :: log_was_open, predictor_mode, linear_mode, approximate_mode, &
      transport_valid, same_geometry, pass_guess, use_approximation
    integer, allocatable :: iterations(:)
    integer :: nbf, nocca, noccb, nvira, offset, ltot, nrhs, irhs, &
      npair_total, first_pair, absolute_pair, ios, exact_every, warmup_exact
    character(len=32) :: mode, env

    if (infos%control%scftype /= 3) then
      call show_message( &
        'mrsf_nac_rohf_zvector_batch requires an ROHF/ROKS reference.', &
        WITH_ABORT)
    end if

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira
    nrhs = size(rhs,2)
    if (nrhs < 1 .or. size(rhs,1) /= ltot .or. &
        size(solution,1) /= ltot .or. size(solution,2) /= nrhs) then
      call show_message( &
        'Batched ROHF NAC sources/solutions have inconsistent dimensions.', &
        WITH_ABORT)
    end if

    first_pair = 1
    if (present(pair_offset)) first_pair = pair_offset + 1
    npair_total = int(infos%tddft%nstate)* &
                  (int(infos%tddft%nstate)-1)/2
    if (first_pair < 1 .or. first_pair + nrhs - 1 > npair_total) then
      call show_message('Invalid unordered-pair offset for NAC Z-vector batch.', &
                        WITH_ABORT)
    end if

    mode = 'off'
    call get_environment_variable('OQP_MRSF_NAC_ZV_PREDICTOR',mode)
    predictor_mode = trim(mode) == 'transport' .or. &
                     trim(mode) == 'TRANSPORT' .or. &
                     trim(mode) == 'linear' .or. trim(mode) == 'LINEAR'
    approximate_mode = trim(mode) == 'transport_approx' .or. &
      trim(mode) == 'TRANSPORT_APPROX' .or. &
      trim(mode) == 'linear_approx' .or. trim(mode) == 'LINEAR_APPROX'
    predictor_mode = predictor_mode .or. approximate_mode
    linear_mode = trim(mode) == 'linear' .or. trim(mode) == 'LINEAR' .or. &
      trim(mode) == 'linear_approx' .or. trim(mode) == 'LINEAR_APPROX'
    eta = 1.0_dp
    env = ''
    call get_environment_variable('OQP_MRSF_NAC_ZV_ETA',env)
    if (len_trim(env) > 0) then
      read(env,*,iostat=ios) eta
      if (ios /= 0 .or. .not. ieee_is_finite(eta)) eta = 1.0_dp
    end if
    eta = max(0.0_dp,min(2.0_dp,eta))
    max_disp = 0.25_dp
    env = ''
    call get_environment_variable('OQP_MRSF_NAC_ZV_MAX_DISP',env)
    if (len_trim(env) > 0) then
      read(env,*,iostat=ios) max_disp
      if (ios /= 0 .or. .not. ieee_is_finite(max_disp)) max_disp = 0.25_dp
    end if
    max_disp = max(0.0_dp,max_disp)
    exact_every = 0
    env = ''
    call get_environment_variable('OQP_MRSF_NAC_ZV_EXACT_EVERY',env)
    if (len_trim(env) > 0) then
      read(env,*,iostat=ios) exact_every
      if (ios /= 0) exact_every = 0
    end if
    exact_every = max(0,exact_every)
    warmup_exact = 1
    if (linear_mode) warmup_exact = 2
    env = ''
    call get_environment_variable('OQP_MRSF_NAC_ZV_WARMUP_EXACT',env)
    if (len_trim(env) > 0) then
      read(env,*,iostat=ios) warmup_exact
      if (ios /= 0) then
        warmup_exact = 1
        if (linear_mode) warmup_exact = 2
      end if
    end if
    warmup_exact = max(1,warmup_exact)

    allocate(converged(nrhs), residual(nrhs), initial_residual(nrhs), &
             iterations(nrhs), guess(ltot,nrhs), guess_available(nrhs), &
             guess_accepted(nrhs), recent_current(ltot), &
             earlier_current(ltot), q_closed(noccb,noccb), &
             q_open(offset,offset), q_virtual(nvira,nvira))
    guess = 0.0_dp; guess_available = .false.; guess_accepted = .false.
    transport_valid = .false.; singular_min = 0.0_dp
    if (predictor_mode) then
      call nac_z_cache_prepare(ltot,npair_total,int(infos%mol_prop%natom),nbf, &
                               nocca,noccb)
      if (any(nac_z_have_recent(first_pair:first_pair+nrhs-1))) then
        call nac_z_transport_matrices(infos,q_closed,q_open,q_virtual, &
                                      singular_min,transport_valid)
      end if
      do irhs = 1, nrhs
        absolute_pair = first_pair + irhs - 1
        if (.not. nac_z_have_recent(absolute_pair) .or. &
            .not. nac_z_have_geometry(absolute_pair)) cycle
        displacement = sqrt(sum((infos%atoms%xyz - &
          nac_z_geometry(:,:,absolute_pair))**2)/ &
          real(max(1,int(infos%mol_prop%natom)),dp))
        same_geometry = maxval(abs(infos%atoms%xyz - &
          nac_z_geometry(:,:,absolute_pair))) <= 1.0e-12_dp
        if (.not. same_geometry .and. &
            (.not. transport_valid .or. displacement > max_disp)) cycle
        if (same_geometry) then
          recent_current = nac_z_recent(:,absolute_pair)
        else
          call nac_z_transport_vector(nac_z_recent(:,absolute_pair), &
            recent_current,q_closed,q_open,q_virtual,nocca,noccb,nbf)
        end if
        guess(:,irhs) = recent_current
        if (linear_mode .and. nac_z_have_earlier(absolute_pair)) then
          if (same_geometry) then
            earlier_current = nac_z_earlier(:,absolute_pair)
          else
            call nac_z_transport_vector(nac_z_earlier(:,absolute_pair), &
              earlier_current,q_closed,q_open,q_virtual,nocca,noccb,nbf)
          end if
          guess(:,irhs) = recent_current + &
                          eta*(recent_current-earlier_current)
        end if
        guess_available(irhs) = all(ieee_is_finite(guess(:,irhs)))
      end do
    else
      call nac_z_cache_reset()
    end if

    pass_guess = predictor_mode .and. all(guess_available)
    use_approximation = approximate_mode .and. pass_guess
    if (use_approximation) use_approximation = all(nac_z_exact_count( &
      first_pair:first_pair+nrhs-1) >= warmup_exact)
    if (use_approximation .and. exact_every > 0) then
      use_approximation = .not. any(nac_z_steps_since_exact( &
        first_pair:first_pair+nrhs-1) >= exact_every-1)
    end if

    ! The energy driver closes IW before the Python NAC orchestrator runs.
    ! Open the requested log locally so the shared CPHF reporter never creates
    ! a process-global fort.6 file, which would collide across worker jobs.
    inquire(unit=iw, opened=log_was_open)
    if (.not. log_was_open) &
      open(unit=iw, file=infos%log_filename, position='append')
    ! Keep requesting ||H z - rhs||_2 <= 1e-10.  Large response spaces can
    ! reach machine-precision stagnation just above that target.  Accept the
    ! best MINRES iterate only when its residual norm, scaled by
    ! max(1,||rhs||_2), is no larger than 10*sqrt(machine epsilon).  This
    ! provides an absolute floor for small sources and a relative criterion
    ! for large sources while remaining independent of response-space size.
    if (use_approximation) then
      ! Deliberately omit the ROHF response solve.  The transported or linearly
      ! extrapolated adjoint replaces z for this nuclear step; all explicit
      ! current-geometry NAC terms and analytic adjoint contractions remain.
      ! A negative residual sentinel denotes "not evaluated", never success
      ! of the certified equation.  Periodic exact refresh is controlled by
      ! OQP_MRSF_NAC_ZV_EXACT_EVERY; tracking/overlap invalidation also falls
      ! through to the full solve below.
      solution = guess
      converged = .true.
      residual = -1.0_dp
      initial_residual = -1.0_dp
      iterations = 0
      guess_accepted = .true.
    else if (pass_guess) then
      call cphf_solve_rohf(infos, nrhs, rhs, solution, tol=1.0e-20_dp, &
                           maxit=max(int(infos%control%maxit_zv), ltot + 5), &
                           converged=converged, residual=residual, &
                           minres_solver=.true., initial_guess=guess, &
                           initial_residual=initial_residual, &
                           iterations=iterations, &
                           initial_guess_accepted=guess_accepted)
    else
      call cphf_solve_rohf(infos, nrhs, rhs, solution, tol=1.0e-20_dp, &
                           maxit=max(int(infos%control%maxit_zv), ltot + 5), &
                           converged=converged, residual=residual, &
                           minres_solver=.true., &
                           initial_residual=initial_residual, &
                           iterations=iterations, &
                           initial_guess_accepted=guess_accepted)
    end if
    do irhs = 1, nrhs
      rhs_scale_sq = max(1.0_dp, sum(rhs(:,irhs)*rhs(:,irhs)))
      scaled_residual_sq = residual(irhs)/rhs_scale_sq
      if (.not. converged(irhs) .and. &
          scaled_residual_sq > fallback_rel_sq) then
        call show_message( &
          'Native ROHF NAC batched Z-vector RHS ' // &
          trim(integer_to_string(irhs)) // &
          ' failed to converge; squared residual=' // &
          trim(real_to_string(residual(irhs))) // &
          '; scaled squared residual=' // &
          trim(real_to_string(scaled_residual_sq)), WITH_ABORT)
      else if (.not. converged(irhs)) then
        call show_message( &
          'Native ROHF NAC batched Z-vector RHS ' // &
          trim(integer_to_string(irhs)) // &
          ' accepted after MINRES stagnation; squared residual=' // &
          trim(real_to_string(residual(irhs))) // &
          '; scaled squared residual=' // &
          trim(real_to_string(scaled_residual_sq)) // &
          ' (fallback relative norm=10*sqrt(epsilon)).', WITHOUT_ABORT)
      end if
      z_relative_error = 0.0_dp
      if (guess_available(irhs)) z_relative_error = &
        sqrt(sum((guess(:,irhs)-solution(:,irhs))**2))/ &
        max(sqrt(sum(solution(:,irhs)**2)),tiny(1.0_dp))
      write(iw,'(A,1X,I0,1X,A,1X,3(L1,1X),5(ES16.8,1X),I0)') &
        'NAC_Z_PREDICTOR', first_pair+irhs-1, trim(mode), &
        guess_available(irhs), guess_accepted(irhs), use_approximation, &
        eta, singular_min, &
        initial_residual(irhs), residual(irhs), z_relative_error, &
        iterations(irhs)
    end do

    if (predictor_mode) then
      do irhs = 1, nrhs
        absolute_pair = first_pair + irhs - 1
        same_geometry = nac_z_have_geometry(absolute_pair) .and. &
          maxval(abs(infos%atoms%xyz - &
          nac_z_geometry(:,:,absolute_pair))) <= 1.0e-12_dp
        if (nac_z_have_recent(absolute_pair) .and. .not. same_geometry .and. &
            transport_valid) then
          call nac_z_transport_vector(nac_z_recent(:,absolute_pair), &
            nac_z_earlier(:,absolute_pair),q_closed,q_open,q_virtual, &
            nocca,noccb,nbf)
          nac_z_have_earlier(absolute_pair) = .true.
        else if (.not. same_geometry) then
          nac_z_have_earlier(absolute_pair) = .false.
        end if
        nac_z_recent(:,absolute_pair) = solution(:,irhs)
        nac_z_have_recent(absolute_pair) = .true.
        nac_z_geometry(:,:,absolute_pair) = infos%atoms%xyz
        nac_z_have_geometry(absolute_pair) = .true.
        if (use_approximation) then
          nac_z_steps_since_exact(absolute_pair) = &
            nac_z_steps_since_exact(absolute_pair) + 1
        else
          nac_z_steps_since_exact(absolute_pair) = 0
          nac_z_exact_count(absolute_pair) = &
            min(huge(nac_z_exact_count(absolute_pair)), &
                nac_z_exact_count(absolute_pair) + 1)
        end if
      end do
    end if

    if (present(predictor)) then
      if (size(predictor,1) /= ltot .or. size(predictor,2) < nrhs) &
        error stop 'mrsf_nac_rohf_zvector_batch: predictor output shape'
      predictor(:,1:nrhs) = guess
    end if
    if (present(predictor_available)) then
      if (size(predictor_available) < nrhs) error stop &
        'mrsf_nac_rohf_zvector_batch: predictor_available output shape'
      predictor_available(1:nrhs) = guess_available
    end if
    if (present(predictor_accepted)) then
      if (size(predictor_accepted) < nrhs) error stop &
        'mrsf_nac_rohf_zvector_batch: predictor_accepted output shape'
      predictor_accepted(1:nrhs) = guess_accepted
    end if
    if (present(initial_residual_out)) then
      if (size(initial_residual_out) < nrhs) error stop &
        'mrsf_nac_rohf_zvector_batch: initial residual output shape'
      initial_residual_out(1:nrhs) = initial_residual
    end if
    if (present(final_residual_out)) then
      if (size(final_residual_out) < nrhs) error stop &
        'mrsf_nac_rohf_zvector_batch: final residual output shape'
      final_residual_out(1:nrhs) = residual
    end if
    if (present(iterations_out)) then
      if (size(iterations_out) < nrhs) error stop &
        'mrsf_nac_rohf_zvector_batch: iterations output shape'
      iterations_out(1:nrhs) = iterations
    end if
    if (.not. log_was_open) close(iw)
    deallocate(converged,residual,initial_residual,iterations,guess, &
               guess_available,guess_accepted,recent_current, &
               earlier_current,q_closed,q_open,q_virtual)
  contains
    function integer_to_string(value) result(text)
      integer, intent(in) :: value
      character(len=16) :: text
      write(text,'(I0)') value
    end function integer_to_string

    function real_to_string(value) result(text)
      real(kind=dp), intent(in) :: value
      character(len=32) :: text
      write(text,'(ES16.8)') value
    end function real_to_string
  end subroutine mrsf_nac_rohf_zvector_batch

!###############################################################################

!> Fortran ABI/source-compatible alias for the provisional entry-point name.
!> Keep this wrapper free of solver logic: the state-pair Z-vector routine above
!> is the sole implementation.
  subroutine mrsf_nac_rohf_solve(infos)
    use types, only: information
    implicit none
    type(information), target, intent(inout) :: infos
    call mrsf_nac_rohf_zvector(infos)
  end subroutine mrsf_nac_rohf_solve

!###############################################################################

!> Form the native ROHF dual source and the explicit overlap/reorthonormalization
!> contribution for one ordered MRSF state pair.  The matrix inputs are flat
!> Fortran-column-major arrays so the C/Python boundary has no ambiguous 2-D
!> TagArray transpose.
!>
!> Input : OQP::nac_mt_frozen   vec(M_IJ^frozen), length nbf**2
!>         OQP::nac_mt_response vec(M_IJ^response), length nbf**2
!>         OQP::nac_gamma_pair  vec(gamma_IJ), length nbf**2
!> Output: OQP::nac_xmat        vec(sum of the three sources), length nbf**2
!>         OQP::nac_rohf_rhs    E^T(M_IJ-M_IJ^T), length ltot
!>         OQP::nac_pair_vmask   = M_IJ:V^R, (3,natom)
!>         OQP::nac_pair_gsk     = gamma_IJ:Sk^R, (3,natom)
!>         OQP::nac_pair_overlap = their sum, (3,natom)
!>
!> V^R is the dependent MO response fixed by orthonormality.  This routine is
!> the resident-Fortran counterpart of the former Python
!> _symmetric_u_contraction plus gamma:Sk coordinate loop.
!> When metric_only is true, the two direct MRSF source records are neither
!> required nor read and xmat=gamma.  The production driver uses that path for
!> the reverse orientation after scaling gamma by -1/2; the forward orientation
!> supplies the direct source once and gamma scaled by +1/2.
  subroutine mrsf_nac_rohf_pair_overlap(infos, metric_only)
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, OQP_VEC_MO_A, &
      TA_TYPE_REAL64
    use grd1, only: der_overlap_matrix_ket, der_overlap_matrix
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_mt_frozen = "OQP::nac_mt_frozen"
    character(len=*), parameter :: tag_mt_response = "OQP::nac_mt_response"
    character(len=*), parameter :: tag_xmat = "OQP::nac_xmat"
    character(len=*), parameter :: tag_gamma = "OQP::nac_gamma_pair"
    character(len=*), parameter :: tag_rhs = "OQP::nac_rohf_rhs"
    character(len=*), parameter :: tag_out = "OQP::nac_pair_overlap"
    character(len=*), parameter :: tag_vmask = "OQP::nac_pair_vmask"
    character(len=*), parameter :: tag_gsk = "OQP::nac_pair_gsk"
    type(information), target, intent(inout) :: infos
    logical, intent(in), optional :: metric_only
    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: mo(:,:), mt_frozen_flat(:), &
      mt_response_flat(:), gflat(:)
    real(kind=dp), pointer :: rhs(:), xflat_out(:), out(:,:), &
      out_vmask(:,:), out_gsk(:,:)
    real(kind=dp), allocatable :: xmat(:,:), gamma(:,:), dsket(:,:,:,:), &
      dsfull(:,:,:,:), overlap_weight(:,:), overlap_weight_ao(:,:), &
      gamma_ao(:,:), half(:,:)
    real(kind=dp) :: value, gsk, coefficient, norm_product
    integer :: nbf, natom, nocca, noccb, nvira, offset, ltot
    integer :: atom, cart, mu, nu, p, q, hi, lo, k
    logical :: only_metric

    if (infos%control%scftype /= 3) then
      call show_message('mrsf_nac_rohf_pair_overlap requires an ROHF/ROKS reference.', &
                        WITH_ABORT)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    natom = infos%mol_prop%natom
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira

    only_metric = .false.
    if (present(metric_only)) only_metric = metric_only
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, tag_gamma, gflat)
    if (size(gflat) /= nbf*nbf) then
      call show_message('MRSF NAC pair matrices have the wrong dimension.', &
                        WITH_ABORT)
    end if
    if (.not. only_metric) then
      call tagarray_get_data(infos%dat, tag_mt_frozen, mt_frozen_flat)
      call tagarray_get_data(infos%dat, tag_mt_response, mt_response_flat)
      if (size(mt_frozen_flat) /= nbf*nbf .or. &
          size(mt_response_flat) /= nbf*nbf) then
        call show_message('MRSF NAC pair matrices have the wrong dimension.', &
                          WITH_ABORT)
      end if
    end if

    allocate(xmat(nbf,nbf), gamma(nbf,nbf), &
             dsket(nbf,nbf,3,natom), dsfull(nbf,nbf,3,natom), &
             overlap_weight(nbf,nbf), overlap_weight_ao(nbf,nbf), &
             gamma_ao(nbf,nbf), half(nbf,nbf))
    gamma = reshape(gflat, (/ nbf, nbf /))
    if (only_metric) then
      xmat = gamma
    else
      xmat = reshape(mt_frozen_flat, (/ nbf, nbf /)) &
           + reshape(mt_response_flat, (/ nbf, nbf /)) + gamma
    end if

    call infos%dat%erase((/ character(len=80) :: tag_xmat, tag_rhs, &
                                                        tag_out, tag_vmask, &
                                                        tag_gsk /))
    call tagarray_reserve_data(infos%dat, tag_xmat, TA_TYPE_REAL64, nbf*nbf, &
         (/ nbf*nbf /), comment='assembled ordered MRSF orbital source')
    call tagarray_reserve_data(infos%dat, tag_rhs, TA_TYPE_REAL64, ltot, (/ ltot /), &
         comment='native ROHF dual of ordered MRSF orbital source')
    call tagarray_reserve_data(infos%dat, tag_out, TA_TYPE_REAL64, 3*natom, &
         (/ 3, natom /), &
         comment='ordered MRSF overlap and dependent-MO response')
    call tagarray_reserve_data(infos%dat, tag_vmask, TA_TYPE_REAL64, 3*natom, &
         (/ 3, natom /), comment='ordered MRSF dependent-MO response')
    call tagarray_reserve_data(infos%dat, tag_gsk, TA_TYPE_REAL64, 3*natom, &
         (/ 3, natom /), comment='ordered MRSF metric-overlap response')
    ! A TagArray reserve can relocate the backing store for records other than
    ! the one being added.  The MO record itself is unchanged, so reacquire its
    ! pointer after the final reserve before the AO-to-MO contractions below.
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, tag_xmat, xflat_out)
    call tagarray_get_data(infos%dat, tag_rhs, rhs)
    call tagarray_get_data(infos%dat, tag_out, out)
    call tagarray_get_data(infos%dat, tag_vmask, out_vmask)
    call tagarray_get_data(infos%dat, tag_gsk, out_gsk)
    xflat_out = reshape(xmat, (/ nbf*nbf /))

    ! Exact tangent-dual projection E^T(M-M^T), in the same SD/DV/SV order
    ! as rohf_unpack_trial and cphf_solve_rohf.  There is deliberately no
    ! extra 1/2: pack_rohf_dual acts on the full generator derivative.
    k = 0
    do p = noccb + 1, nocca
      do q = 1, noccb
        k = k + 1
        rhs(k) = xmat(p,q) - xmat(q,p)
      end do
    end do
    do p = nocca + 1, nbf
      do q = 1, noccb
        k = k + 1
        rhs(k) = xmat(p,q) - xmat(q,p)
      end do
    end do
    do p = nocca + 1, nbf
      do q = noccb + 1, nocca
        k = k + 1
        rhs(k) = xmat(p,q) - xmat(q,p)
      end do
    end do

    ! The coordinate loop used to transform every derivative overlap matrix
    ! AO->MO (two nbf**3 DGEMMs for each of 3*natom coordinates).  Both outputs
    ! are Frobenius contractions, so transform their pair-dependent weights
    ! once in the reverse direction instead:
    !   <W, C^T S^R C> = <C W C^T, S^R>.
    ! overlap_weight reproduces the original lower-triangular
    ! diagonal/same-space/cross-space dependent-MO contraction exactly.  It is
    ! deliberately not symmetrized, so this identity also holds before using
    ! the mathematical symmetry of the full overlap derivative.
    overlap_weight = 0.0_dp
    do p = 1, nbf
      overlap_weight(p,p) = -0.5_dp*xmat(p,p)
      do q = 1, p - 1
        if (orbital_space(p) == orbital_space(q)) then
          coefficient = -0.5_dp*(xmat(p,q)+xmat(q,p))
        else
          if (orbital_space(p) > orbital_space(q)) then
            hi = p; lo = q
          else
            hi = q; lo = p
          end if
          coefficient = -xmat(lo,hi)
        end if
        overlap_weight(p,q) = coefficient
      end do
    end do
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, &
               overlap_weight, nbf, 0.0_dp, half, nbf)
    call dgemm('n','t', nbf, nbf, nbf, 1.0_dp, half, nbf, mo, nbf, &
               0.0_dp, overlap_weight_ao, nbf)
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, gamma, nbf, &
               0.0_dp, half, nbf)
    call dgemm('n','t', nbf, nbf, nbf, 1.0_dp, half, nbf, mo, nbf, &
               0.0_dp, gamma_ao, nbf)
    do nu = 1, nbf
      do mu = 1, nbf
        norm_product = basis%bfnrm(mu)*basis%bfnrm(nu)
        overlap_weight_ao(mu,nu) = overlap_weight_ao(mu,nu)*norm_product
        gamma_ao(mu,nu) = gamma_ao(mu,nu)*norm_product
      end do
    end do

    call der_overlap_matrix_ket(basis, dsket)
    call der_overlap_matrix(basis, dsfull)
    do atom = 1, natom
      do cart = 1, 3
        value = sum(overlap_weight_ao*dsfull(:,:,cart,atom))
        gsk = sum(gamma_ao*dsket(:,:,cart,atom))
        out_vmask(cart,atom) = value
        out_gsk(cart,atom) = gsk
        out(cart,atom) = value + gsk
      end do
    end do

    deallocate(xmat, gamma, dsket, dsfull, overlap_weight, &
               overlap_weight_ao, &
               gamma_ao, half)

  contains
    pure integer function orbital_space(iorb) result(space)
      integer, intent(in) :: iorb
      if (iorb <= noccb) then
        space = 1
      else if (iorb <= nocca) then
        space = 2
      else
        space = 3
      end if
    end function orbital_space

  end subroutine mrsf_nac_rohf_pair_overlap

!###############################################################################

!> Allocate the resident ordered-pair NAC accumulator and discard stale final
!> tensors from an earlier call on the same molecule handle.
!>
!> Output: OQP::nac_dp_ordered (3*natom,nstate,nstate)
  subroutine mrsf_nac_pair_accumulator_init(infos)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, TA_TYPE_REAL64
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_dp = "OQP::nac_dp_ordered"
    character(len=*), parameter :: tag_dcv = "OQP::nac_dcv"
    character(len=*), parameter :: tag_nacv = "OQP::nac_nacv"
    type(information), target, intent(inout) :: infos
    real(kind=dp), pointer :: dp_ordered(:,:,:)
    integer :: natom, nstate, ncoord

    natom = infos%mol_prop%natom
    nstate = infos%tddft%nstate
    ncoord = 3*natom
    if (natom < 1 .or. nstate < 1) then
      call show_message('Cannot initialize an empty MRSF NAC pair tensor.', &
                        WITH_ABORT)
    end if

    call infos%dat%erase((/ character(len=80) :: tag_dp, tag_dcv, &
                                                        tag_nacv /))
    call tagarray_reserve_data(infos%dat, tag_dp, TA_TYPE_REAL64, &
         ncoord*nstate*nstate, (/ ncoord, nstate, nstate /), &
         comment='resident ordered MRSF NAC Lagrangian vectors')
    call tagarray_get_data(infos%dat, tag_dp, dp_ordered)
    dp_ordered = 0.0_dp
  end subroutine mrsf_nac_pair_accumulator_init

!###############################################################################

!> Add all resident analytic components for one ordered pair.  The selected
!> amplitude block and current-pair coordinate records were produced by the
!> immediately preceding Fortran kernels in the production call sequence.
!>
!> Inputs: OQP::nac_amp                  (3*natom,nstate,nstate)
!>         OQP::nac_esum                 (3,natom)
!>         OQP::nac_rohf_hf_adjoint      (3,natom)
!>         OQP::nac_rohf_xc_adjoint      (3,natom)
!>         OQP::nac_pair_overlap         (3,natom)
!> In/out: OQP::nac_dp_ordered           (3*natom,nstate,nstate)
  subroutine mrsf_nac_pair_accumulate(infos, istate, jstate)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_amp = "OQP::nac_amp"
    character(len=*), parameter :: tag_esum = "OQP::nac_esum"
    character(len=*), parameter :: tag_hf = "OQP::nac_rohf_hf_adjoint"
    character(len=*), parameter :: tag_xc = "OQP::nac_rohf_xc_adjoint"
    character(len=*), parameter :: tag_overlap = "OQP::nac_pair_overlap"
    character(len=*), parameter :: tag_dp = "OQP::nac_dp_ordered"
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: istate, jstate
    real(kind=dp), contiguous, pointer :: amp(:,:,:), esum(:,:), &
      z_hf(:,:), z_xc(:,:), pair_overlap(:,:)
    real(kind=dp), pointer :: dp_ordered(:,:,:)
    real(kind=dp) :: t1, z_response
    integer :: natom, nstate, ncoord, atom, cart, coord

    natom = infos%mol_prop%natom
    nstate = infos%tddft%nstate
    ncoord = 3*natom
    if (istate < 1 .or. istate > nstate .or. &
        jstate < 1 .or. jstate > nstate .or. istate == jstate) then
      call show_message('Invalid ordered pair in MRSF NAC accumulation.', &
                        WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, tag_amp, amp)
    call tagarray_get_data(infos%dat, tag_esum, esum)
    call tagarray_get_data(infos%dat, tag_hf, z_hf)
    call tagarray_get_data(infos%dat, tag_xc, z_xc)
    call tagarray_get_data(infos%dat, tag_overlap, pair_overlap)
    call tagarray_get_data(infos%dat, tag_dp, dp_ordered)
    if (size(amp,1) /= ncoord .or. size(amp,2) /= nstate .or. &
        size(amp,3) /= nstate .or. &
        size(esum,1) /= 3 .or. size(esum,2) /= natom .or. &
        size(z_hf,1) /= 3 .or. size(z_hf,2) /= natom .or. &
        size(z_xc,1) /= 3 .or. size(z_xc,2) /= natom .or. &
        size(pair_overlap,1) /= 3 .or. &
        size(pair_overlap,2) /= natom .or. &
        size(dp_ordered,1) /= ncoord .or. &
        size(dp_ordered,2) /= nstate .or. &
        size(dp_ordered,3) /= nstate) then
      call show_message('MRSF NAC resident pair components have inconsistent dimensions.', &
                        WITH_ABORT)
    end if

    do atom = 1, natom
      do cart = 1, 3
        coord = (atom - 1)*3 + cart
        t1 = amp(coord,istate,jstate) + esum(cart,atom)
        z_response = z_hf(cart,atom) + z_xc(cart,atom)
        dp_ordered(coord,istate,jstate) = t1 + z_response &
          + pair_overlap(cart,atom)
      end do
    end do
  end subroutine mrsf_nac_pair_accumulate

!###############################################################################

!> Store one already-antisymmetrized unordered-pair Lagrangian vector.
!>
!> ``nonz_antisym`` contains one half of the ordered-pair difference for all
!> explicit amplitude, one-electron/Fock and overlap terms.  The current
!> OQP::nac_rohf_z record is the solution of the matching half-difference RHS,
!> so the HF and XC adjoints are already the corresponding half-difference by
!> linearity.  Store +/- the complete vector in the legacy ordered accumulator;
!> mrsf_nac_pair_finalize then reproduces the same canonical antisymmetrization
!> while preserving the existing resident output layout and ABI.
  subroutine mrsf_nac_pair_accumulate_antisym(infos, istate, jstate, &
                                               nonz_antisym)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_hf = "OQP::nac_rohf_hf_adjoint"
    character(len=*), parameter :: tag_xc = "OQP::nac_rohf_xc_adjoint"
    character(len=*), parameter :: tag_dp = "OQP::nac_dp_ordered"
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: istate, jstate
    real(kind=dp), intent(in) :: nonz_antisym(:)
    real(kind=dp), contiguous, pointer :: z_hf(:,:), z_xc(:,:)
    real(kind=dp), pointer :: dp_ordered(:,:,:)
    real(kind=dp) :: value
    integer :: natom, nstate, ncoord, atom, cart, coord

    natom = infos%mol_prop%natom
    nstate = infos%tddft%nstate
    ncoord = 3*natom
    if (istate < 1 .or. jstate > nstate .or. istate >= jstate) then
      call show_message( &
        'Antisymmetric MRSF NAC accumulation requires 1 <= I < J <= nstate.', &
        WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, tag_hf, z_hf)
    call tagarray_get_data(infos%dat, tag_xc, z_xc)
    call tagarray_get_data(infos%dat, tag_dp, dp_ordered)
    if (size(nonz_antisym) /= ncoord .or. &
        size(z_hf,1) /= 3 .or. size(z_hf,2) /= natom .or. &
        size(z_xc,1) /= 3 .or. size(z_xc,2) /= natom .or. &
        size(dp_ordered,1) /= ncoord .or. &
        size(dp_ordered,2) /= nstate .or. &
        size(dp_ordered,3) /= nstate) then
      call show_message( &
        'Antisymmetric MRSF NAC pair components have inconsistent dimensions.', &
        WITH_ABORT)
    end if

    do atom = 1, natom
      do cart = 1, 3
        coord = (atom - 1)*3 + cart
        value = nonz_antisym(coord) + z_hf(cart,atom) + z_xc(cart,atom)
        dp_ordered(coord,istate,jstate) = value
        dp_ordered(coord,jstate,istate) = -value
      end do
    end do
  end subroutine mrsf_nac_pair_accumulate_antisym

!###############################################################################

!> Antisymmetrize the complete ordered-pair Lagrangian tensor and form the
!> canonical gap-scaled coupling in resident Fortran:
!>
!>   d_IJ = 1/2 (dp_IJ - dp_JI)
!>   h_IJ = (Omega_J - Omega_I) d_IJ.
!>
!> Input : OQP::nac_dp_ordered (3*natom,nstate,nstate)
!>         OQP::td_energies    (nstate), MRSF excitation energies
!> Output: OQP::nac_dcv        (3*natom,nstate,nstate)
!>         OQP::nac_nacv       (3*natom,nstate,nstate)
  subroutine mrsf_nac_pair_finalize(infos)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, TA_TYPE_REAL64, &
      OQP_td_energies
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_dp = "OQP::nac_dp_ordered"
    character(len=*), parameter :: tag_dcv = "OQP::nac_dcv"
    character(len=*), parameter :: tag_nacv = "OQP::nac_nacv"
    type(information), target, intent(inout) :: infos
    real(kind=dp), contiguous, pointer :: dp_ordered(:,:,:), energies(:)
    real(kind=dp), pointer :: dcv(:,:,:), nacv(:,:,:)
    real(kind=dp) :: gap
    integer :: natom, nstate, ncoord, istate, jstate, coord

    natom = infos%mol_prop%natom
    nstate = infos%tddft%nstate
    ncoord = 3*natom
    call tagarray_get_data(infos%dat, tag_dp, dp_ordered)
    call tagarray_get_data(infos%dat, OQP_td_energies, energies)
    if (size(dp_ordered,1) /= ncoord .or. &
        size(dp_ordered,2) /= nstate .or. &
        size(dp_ordered,3) /= nstate .or. size(energies) /= nstate) then
      call show_message('Cannot finalize inconsistent MRSF NAC pair tensors.', &
                        WITH_ABORT)
    end if

    call infos%dat%erase((/ character(len=80) :: tag_dcv, &
                                                        tag_nacv /))
    call tagarray_reserve_data(infos%dat, tag_dcv, TA_TYPE_REAL64, &
         ncoord*nstate*nstate, (/ ncoord, nstate, nstate /), &
         comment='antisymmetric MRSF derivative coupling')
    call tagarray_reserve_data(infos%dat, tag_nacv, TA_TYPE_REAL64, &
         ncoord*nstate*nstate, (/ ncoord, nstate, nstate /), &
         comment='symmetric gap-scaled MRSF nonadiabatic coupling')
    ! Output reservations may relocate unrelated TagArray records.  Reacquire
    ! both unchanged inputs before using them to finalize the pair tensors.
    call tagarray_get_data(infos%dat, tag_dp, dp_ordered)
    call tagarray_get_data(infos%dat, OQP_td_energies, energies)
    call tagarray_get_data(infos%dat, tag_dcv, dcv)
    call tagarray_get_data(infos%dat, tag_nacv, nacv)
    dcv = 0.0_dp
    nacv = 0.0_dp

    do jstate = 1, nstate
      do istate = 1, nstate
        if (istate == jstate) cycle
        gap = energies(jstate) - energies(istate)
        do coord = 1, ncoord
          dcv(coord,istate,jstate) = 0.5_dp * &
            (dp_ordered(coord,istate,jstate) &
             - dp_ordered(coord,jstate,istate))
          nacv(coord,istate,jstate) = gap*dcv(coord,istate,jstate)
        end do
      end do
    end do
  end subroutine mrsf_nac_pair_finalize

!###############################################################################

  subroutine mrsf_nac_rohf_hf_adjoint_C(c_handle) &
      bind(C, name="mrsf_nac_rohf_hf_adjoint")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_rohf_hf_adjoint(inf)
  end subroutine mrsf_nac_rohf_hf_adjoint_C

!###############################################################################

!> Contract one native ROHF Z vector with the analytic HF/JK/Pulay part of
!> every nuclear stationarity derivative.  This is the adjoint (Z-vector)
!> evaluation of z^T B^R: it never solves the 3N forward CPHF equations.
!>
!> The algebra is the exact transpose of the non-canonical RHS assembled in
!> hf_hessian_rohf.  In particular, the expensive orbital-by-orbital 2e
!> skeleton is collapsed by response symmetry into one derivative contraction
!> for each spin, and G[D^(R,0)] is collapsed into S^R:G[P_z].
!>
!> Input : OQP::nac_rohf_z              (ltot)
!> Output: OQP::nac_rohf_hf_adjoint     (3,natom)
  subroutine mrsf_nac_rohf_hf_adjoint(infos)
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, OQP_DM_A, OQP_DM_B, &
      OQP_VEC_MO_A, OQP_FOCK_A, OQP_FOCK_B, TA_TYPE_REAL64
    use mathlib, only: unpack_matrix, pack_matrix
    use cphf_mod, only: rohf_unpack_trial
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix
    use fock_deriv_mod, only: fock_deriv_contract_os
    use scf_addons, only: fock_jk
    use ecp_tool, only: ecp_deriv_ints
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_z = "OQP::nac_rohf_z"
    character(len=*), parameter :: tag_out = "OQP::nac_rohf_hf_adjoint"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dma(:), dmb(:), mo(:,:), &
      focka(:), fockb(:), z(:)
    real(kind=dp), pointer :: out(:,:)
    real(kind=dp), allocatable :: pa(:,:), pb(:,:), ptot(:,:)
    real(kind=dp), allocatable :: xa(:,:), xb(:,:), pza(:,:), pzb(:,:)
    real(kind=dp), allocatable :: work(:,:), half(:,:), probe(:,:)
    real(kind=dp), allocatable :: fa(:,:), fb(:,:), famo(:,:), fbmo(:,:)
    real(kind=dp), allocatable :: dmz(:,:), vjkz(:,:), vza(:,:), vzb(:,:)
    real(kind=dp), allocatable :: vzamo(:,:), vzbmo(:,:)
    real(kind=dp), allocatable :: dsa(:,:,:,:), dta(:,:,:,:), dva(:,:,:,:), &
      dvecp(:,:,:,:)
    real(kind=dp), allocatable :: ghf(:,:), gx(:,:), sxmo(:,:), hxmo(:,:)
    real(kind=dp) :: hfscale, value
    integer :: nbf, nbf2, natom, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: atom, cart, i, j, a, mu, nu

    if (infos%control%scftype /= 3) then
      call show_message('mrsf_nac_rohf_hf_adjoint requires an ROHF/ROKS reference.', &
                        WITH_ABORT)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = infos%mol_prop%natom
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira
    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale

    call tagarray_get_data(infos%dat, OQP_DM_A, dma)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmb)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, focka)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fockb)
    call tagarray_get_data(infos%dat, tag_z, z)
    if (size(z) /= ltot) then
      call show_message('OQP::nac_rohf_z has the wrong ROHF rotation dimension.', &
                        WITH_ABORT)
    end if

    allocate(pa(nbf,nbf), pb(nbf,nbf), ptot(nbf,nbf))
    allocate(xa(nvira,nocca), xb(nvirb,noccb), pza(nbf,nbf), pzb(nbf,nbf))
    allocate(work(nbf,nbf), half(nbf,nbf), probe(nbf,nbf))
    allocate(fa(nbf,nbf), fb(nbf,nbf), famo(nbf,nbf), fbmo(nbf,nbf))
    allocate(dmz(nbf2,2), vjkz(nbf2,2), vza(nbf,nbf), vzb(nbf,nbf), &
             vzamo(nbf,nbf), vzbmo(nbf,nbf))
    allocate(dsa(nbf,nbf,3,natom), dta(nbf,nbf,3,natom), &
             dva(nbf,nbf,3,natom), dvecp(nbf,nbf,3,natom))
    allocate(ghf(3,natom), gx(3,natom), sxmo(nbf,nbf), hxmo(nbf,nbf), &
             source=0.0_dp)

    call unpack_matrix(dma, pa)
    call unpack_matrix(dmb, pb)
    ptot = pa + pb
    call unpack_matrix(focka, fa)
    call unpack_matrix(fockb, fb)
    call ao_to_mo(fa, famo)
    call ao_to_mo(fb, fbmo)
    call rohf_unpack_trial(z, xa, xb, nbf, nocca, noccb)

    ! Physical first-order spin densities generated by the Z rotations.
    half = 0.0_dp
    call dgemm('n','n', nbf, nocca, nvira, 1.0_dp, &
               mo(:,nocca+1:nbf), nbf, xa, nvira, 0.0_dp, half, nbf)
    call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, half, nbf, &
               mo(:,1:nocca), nbf, 0.0_dp, work, nbf)
    pza = work + transpose(work)
    if (noccb > 0) then
      half = 0.0_dp
      call dgemm('n','n', nbf, noccb, nvirb, 1.0_dp, &
                 mo(:,noccb+1:nbf), nbf, xb, nvirb, 0.0_dp, half, nbf)
      call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, half, nbf, &
                 mo(:,1:noccb), nbf, 0.0_dp, work, nbf)
      pzb = work + transpose(work)
    else
      pzb = 0.0_dp
    end if

    ! One analytic two-electron derivative contraction per spin replaces the
    ! forward code's (virtual,occupied) probe sweep.  The factor 1/2 follows
    ! because P_z contains both AO triangles while each forward probe is
    ! 1/2(C_a C_i^T + C_i C_a^T).
    probe = 0.5_dp*pza
    gx = 0.0_dp
    call fock_deriv_contract_os(infos, basis, ptot, pa, probe, hfscale, gx)
    ghf = ghf - gx
    probe = 0.5_dp*pzb
    gx = 0.0_dp
    call fock_deriv_contract_os(infos, basis, ptot, pb, probe, hfscale, gx)
    ghf = ghf - gx

    ! Build the JK response to P_z once.  Kernel symmetry turns the forward
    ! -G[D^(R,0)] term into +1/2 S^R_occ:G[P_z]_occ for every coordinate.
    call pack_matrix(pza, dmz(:,1))
    call pack_matrix(pzb, dmz(:,2))
    call fock_jk(basis, d=dmz, f=vjkz, scale_exch=hfscale, infos=infos)
    call unpack_matrix(vjkz(:,1), vza)
    call unpack_matrix(vjkz(:,2), vzb)
    call ao_to_mo(vza, vzamo)
    call ao_to_mo(vzb, vzbmo)

    call der_overlap_matrix(basis, dsa)
    call der_kinetic_matrix(basis, dta)
    call der_nucattr_matrix(basis, basis%atoms%xyz, &
                            basis%atoms%zn - basis%ecp_zn_num, dva)
    call ecp_deriv_ints(basis, basis%atoms%xyz, dvecp)
    do atom = 1, natom
      do cart = 1, 3
        do nu = 1, nbf
          do mu = 1, nbf
            dsa(mu,nu,cart,atom) = dsa(mu,nu,cart,atom) * &
              basis%bfnrm(mu)*basis%bfnrm(nu)
            dta(mu,nu,cart,atom) = dta(mu,nu,cart,atom) * &
              basis%bfnrm(mu)*basis%bfnrm(nu)
            dva(mu,nu,cart,atom) = dva(mu,nu,cart,atom) * &
              basis%bfnrm(mu)*basis%bfnrm(nu)
          end do
        end do
        ! ECP derivative integrals are already in the normalized AO convention.
        dva(:,:,cart,atom) = dva(:,:,cart,atom) + dvecp(:,:,cart,atom)
        call ao_to_mo(dsa(:,:,cart,atom), sxmo)
        call ao_to_mo(dta(:,:,cart,atom)+dva(:,:,cart,atom), hxmo)

        value = 0.0_dp
        do i = 1, nocca
          do a = 1, nvira
            value = value + xa(a,i) * ( &
              -hxmo(i,nocca+a) &
              +dot_product(sxmo(nocca+a,1:nocca), famo(1:nocca,i)) &
              +dot_product(famo(nocca+a,1:nocca), sxmo(1:nocca,i)) )
          end do
        end do
        do i = 1, noccb
          do a = 1, nvirb
            value = value + xb(a,i) * ( &
              -hxmo(i,noccb+a) &
              +dot_product(sxmo(noccb+a,1:noccb), fbmo(1:noccb,i)) &
              +dot_product(fbmo(noccb+a,1:noccb), sxmo(1:noccb,i)) )
          end do
        end do
        do j = 1, nocca
          do i = 1, nocca
            value = value + 0.5_dp*sxmo(i,j)*vzamo(j,i)
          end do
        end do
        do j = 1, noccb
          do i = 1, noccb
            value = value + 0.5_dp*sxmo(i,j)*vzbmo(j,i)
          end do
        end do
        ghf(cart,atom) = ghf(cart,atom) + value
      end do
    end do

    call infos%dat%erase((/ character(len=80) :: tag_out /))
    call tagarray_reserve_data(infos%dat, tag_out, TA_TYPE_REAL64, 3*natom, &
         (/ 3, natom /), comment='native ROHF NAC analytic HF/JK/Pulay adjoint')
    call tagarray_get_data(infos%dat, tag_out, out)
    out = ghf

    deallocate(pa, pb, ptot, xa, xb, pza, pzb, work, half, probe, &
      fa, fb, famo, fbmo, dmz, vjkz, vza, vzb, vzamo, vzbmo, &
      dsa, dta, dva, dvecp, ghf, gx, sxmo, hxmo)
  contains
    subroutine ao_to_mo(ao, transformed)
      real(kind=dp), intent(in) :: ao(:,:)
      real(kind=dp), intent(out) :: transformed(:,:)
      half = 0.0_dp
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, ao, nbf, &
                 0.0_dp, half, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, half, nbf, mo, nbf, &
                 0.0_dp, transformed, nbf)
    end subroutine ao_to_mo
  end subroutine mrsf_nac_rohf_hf_adjoint

!###############################################################################

!> Contract several native ROHF Z vectors with the analytic HF/JK/Pulay
!> stationarity derivative.  The ground densities, ground Fock transforms,
!> one-electron derivative integrals and their AO->MO transforms are common to
!> every right-hand side and are therefore evaluated once per batch.
  subroutine mrsf_nac_rohf_hf_adjoint_batch(infos, z_vectors, ghf_vectors)
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, OQP_DM_A, OQP_DM_B, &
      OQP_VEC_MO_A, OQP_FOCK_A, OQP_FOCK_B
    use mathlib, only: unpack_matrix, pack_matrix
    use cphf_mod, only: rohf_unpack_trial
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix
    use fock_deriv_mod, only: fock_deriv_contract_os_batch
    use scf_addons, only: fock_jk
    use ecp_tool, only: ecp_deriv_ints
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: z_vectors(:,:)
    real(kind=dp), intent(out) :: ghf_vectors(:,:,:)

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dma(:), dmb(:), mo(:,:), &
      focka(:), fockb(:)
    real(kind=dp), allocatable :: pa(:,:), pb(:,:), ptot(:,:)
    real(kind=dp), allocatable :: xa(:,:,:), xb(:,:,:)
    real(kind=dp), allocatable :: pza(:,:,:), pzb(:,:,:)
    real(kind=dp), allocatable :: work(:,:), half(:,:)
    real(kind=dp), allocatable :: fa(:,:), fb(:,:), famo(:,:), fbmo(:,:)
    real(kind=dp), allocatable :: dmz(:,:), vjkz(:,:), vza(:,:), vzb(:,:)
    real(kind=dp), allocatable :: vzamo(:,:,:), vzbmo(:,:,:)
    real(kind=dp), allocatable :: dsa(:,:,:,:), dta(:,:,:,:), dva(:,:,:,:), &
      dvecp(:,:,:,:)
    real(kind=dp), allocatable :: gx_batch(:,:,:), sxmo(:,:), hxmo(:,:)
    real(kind=dp) :: hfscale, value
    integer :: nbf, nbf2, natom, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: nrhs, irhs, atom, cart, i, j, a, mu, nu

    if (infos%control%scftype /= 3) then
      call show_message( &
        'mrsf_nac_rohf_hf_adjoint_batch requires an ROHF/ROKS reference.', &
        WITH_ABORT)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = infos%mol_prop%natom
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira
    nrhs = size(z_vectors,2)
    if (nrhs < 1 .or. size(z_vectors,1) /= ltot .or. &
        size(ghf_vectors,1) /= 3 .or. &
        size(ghf_vectors,2) /= natom .or. &
        size(ghf_vectors,3) /= nrhs) then
      call show_message('Batched HF adjoint dimensions are inconsistent.', &
                        WITH_ABORT)
    end if
    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale

    call tagarray_get_data(infos%dat, OQP_DM_A, dma)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmb)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, focka)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fockb)

    allocate(pa(nbf,nbf), pb(nbf,nbf), ptot(nbf,nbf))
    allocate(xa(nvira,nocca,nrhs), xb(nvirb,noccb,nrhs), &
             pza(nbf,nbf,nrhs), pzb(nbf,nbf,nrhs))
    allocate(work(nbf,nbf), half(nbf,nbf))
    allocate(fa(nbf,nbf), fb(nbf,nbf), famo(nbf,nbf), fbmo(nbf,nbf))
    allocate(dmz(nbf2,2*nrhs), vjkz(nbf2,2*nrhs), &
             vza(nbf,nbf), vzb(nbf,nbf), &
             vzamo(nbf,nbf,nrhs), vzbmo(nbf,nbf,nrhs))
    allocate(dsa(nbf,nbf,3,natom), dta(nbf,nbf,3,natom), &
             dva(nbf,nbf,3,natom), dvecp(nbf,nbf,3,natom))
    allocate(gx_batch(3,natom,nrhs), sxmo(nbf,nbf), hxmo(nbf,nbf), &
             source=0.0_dp)

    call unpack_matrix(dma, pa)
    call unpack_matrix(dmb, pb)
    ptot = pa + pb
    call unpack_matrix(focka, fa)
    call unpack_matrix(fockb, fb)
    call ao_to_mo(fa, famo)
    call ao_to_mo(fb, fbmo)

    do irhs = 1, nrhs
      call rohf_unpack_trial(z_vectors(:,irhs), xa(:,:,irhs), xb(:,:,irhs), &
                             nbf, nocca, noccb)

      half = 0.0_dp
      call dgemm('n','n', nbf, nocca, nvira, 1.0_dp, &
                 mo(:,nocca+1:nbf), nbf, xa(:,:,irhs), nvira, &
                 0.0_dp, half, nbf)
      call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, half, nbf, &
                 mo(:,1:nocca), nbf, 0.0_dp, work, nbf)
      pza(:,:,irhs) = work + transpose(work)

      if (noccb > 0) then
        half = 0.0_dp
        call dgemm('n','n', nbf, noccb, nvirb, 1.0_dp, &
                   mo(:,noccb+1:nbf), nbf, xb(:,:,irhs), nvirb, &
                   0.0_dp, half, nbf)
        call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, half, nbf, &
                   mo(:,1:noccb), nbf, 0.0_dp, work, nbf)
        pzb(:,:,irhs) = work + transpose(work)
      else
        pzb(:,:,irhs) = 0.0_dp
      end if

      call pack_matrix(pza(:,:,irhs), dmz(:,2*irhs-1))
      call pack_matrix(pzb(:,:,irhs), dmz(:,2*irhs))
    end do

    ! The response-symmetry contraction uses one half of each physical spin
    ! density.  All 2*nrhs probes now share one derivative-ERI recurrence.
    pza = 0.5_dp*pza
    pzb = 0.5_dp*pzb
    call fock_deriv_contract_os_batch( &
      infos, basis, ptot, pa, pb, pza, pzb, hfscale, gx_batch)
    ghf_vectors = -gx_batch

    ! fock_jk accepts adjacent alpha/beta response pairs.  A single integral
    ! pass now forms the JK response for every state pair in this batch.
    call fock_jk(basis, d=dmz, f=vjkz, scale_exch=hfscale, infos=infos)
    do irhs = 1, nrhs
      call unpack_matrix(vjkz(:,2*irhs-1), vza)
      call unpack_matrix(vjkz(:,2*irhs), vzb)
      call ao_to_mo(vza, vzamo(:,:,irhs))
      call ao_to_mo(vzb, vzbmo(:,:,irhs))
    end do

    ! These integrals and their normalization are independent of the adjoint
    ! vector.  The old scalar call regenerated them for every state pair.
    call der_overlap_matrix(basis, dsa)
    call der_kinetic_matrix(basis, dta)
    call der_nucattr_matrix(basis, basis%atoms%xyz, &
                            basis%atoms%zn - basis%ecp_zn_num, dva)
    call ecp_deriv_ints(basis, basis%atoms%xyz, dvecp)
    do atom = 1, natom
      do cart = 1, 3
        do nu = 1, nbf
          do mu = 1, nbf
            dsa(mu,nu,cart,atom) = dsa(mu,nu,cart,atom) * &
              basis%bfnrm(mu)*basis%bfnrm(nu)
            dta(mu,nu,cart,atom) = dta(mu,nu,cart,atom) * &
              basis%bfnrm(mu)*basis%bfnrm(nu)
            dva(mu,nu,cart,atom) = dva(mu,nu,cart,atom) * &
              basis%bfnrm(mu)*basis%bfnrm(nu)
          end do
        end do
        dva(:,:,cart,atom) = dva(:,:,cart,atom) + dvecp(:,:,cart,atom)
        call ao_to_mo(dsa(:,:,cart,atom), sxmo)
        call ao_to_mo(dta(:,:,cart,atom)+dva(:,:,cart,atom), hxmo)

        do irhs = 1, nrhs
          value = 0.0_dp
          do i = 1, nocca
            do a = 1, nvira
              value = value + xa(a,i,irhs) * ( &
                -hxmo(i,nocca+a) &
                +dot_product(sxmo(nocca+a,1:nocca), famo(1:nocca,i)) &
                +dot_product(famo(nocca+a,1:nocca), sxmo(1:nocca,i)) )
            end do
          end do
          do i = 1, noccb
            do a = 1, nvirb
              value = value + xb(a,i,irhs) * ( &
                -hxmo(i,noccb+a) &
                +dot_product(sxmo(noccb+a,1:noccb), fbmo(1:noccb,i)) &
                +dot_product(fbmo(noccb+a,1:noccb), sxmo(1:noccb,i)) )
            end do
          end do
          do j = 1, nocca
            do i = 1, nocca
              value = value + 0.5_dp*sxmo(i,j)*vzamo(j,i,irhs)
            end do
          end do
          do j = 1, noccb
            do i = 1, noccb
              value = value + 0.5_dp*sxmo(i,j)*vzbmo(j,i,irhs)
            end do
          end do
          ghf_vectors(cart,atom,irhs) = &
            ghf_vectors(cart,atom,irhs) + value
        end do
      end do
    end do

    deallocate(pa, pb, ptot, xa, xb, pza, pzb, work, half, &
      fa, fb, famo, fbmo, dmz, vjkz, vza, vzb, vzamo, vzbmo, &
      dsa, dta, dva, dvecp, gx_batch, sxmo, hxmo)
  contains
    subroutine ao_to_mo(ao, transformed)
      real(kind=dp), intent(in) :: ao(:,:)
      real(kind=dp), intent(out) :: transformed(:,:)
      half = 0.0_dp
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, ao, nbf, &
                 0.0_dp, half, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, half, nbf, mo, nbf, &
                 0.0_dp, transformed, nbf)
    end subroutine ao_to_mo
  end subroutine mrsf_nac_rohf_hf_adjoint_batch

!###############################################################################

  subroutine mrsf_nac_xc_adjoint_C(c_handle) &
      bind(C, name="mrsf_nac_xc_adjoint")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use io_constants, only: iw
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    logical :: log_was_open

    inf => oqp_handle_get_info(c_handle)
    inquire(unit=iw, opened=log_was_open)
    if (.not. log_was_open) &
      open(unit=iw, file=inf%log_filename, position='append')
    call mrsf_nac_xc_adjoint(inf)
    if (.not. log_was_open) close(iw)
  end subroutine mrsf_nac_xc_adjoint_C

!###############################################################################

!> Contract a native-ROHF adjoint vector with the *analytic* nuclear XC part
!> of the SCF stationarity derivative.
!>
!> Let
!>   r = ( F_beta(sd), F_alpha(vd)+F_beta(vd), F_alpha(vs) )
!> be the half-gradient used by the native ROHF CPHF equations and let z use
!> the matching rohf_pack_trial coordinates.  The physical orbital-rotation
!> densities generated by z obey
!>
!>   delta E_xc = Tr[V_xc^a delta P_a] + Tr[V_xc^b delta P_b]
!>              = 2 z^T r_xc .
!>
!> Mixed-derivative symmetry therefore gives the XC contribution to the NAC
!> adjoint contraction as
!>
!>   z^T b_xc = -z^T (d r_xc/dR) = -1/2 d/dR(delta E_xc).
!>
!> utddft_xc_gradient evaluates the last mixed derivative analytically on the
!> moving atom-centred grid, including the normalized fuzzy-cell weight
!> response.  Its probe-independent ground-state contribution is disabled at
!> the consumer level.
!>
!> Input : OQP::nac_rohf_z              (ltot)
!> Output: OQP::nac_rohf_xc_adjoint     (3,natom)
  subroutine mrsf_nac_xc_adjoint(infos)
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, &
      OQP_DM_A, OQP_DM_B, OQP_VEC_MO_A, TA_TYPE_REAL64
    use mathlib, only: unpack_matrix
    use cphf_mod, only: rohf_unpack_trial
    use dft, only: dft_initialize, dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
    use mod_dft_gridint_fxc, only: utddft_fxc
    use grd1, only: der_overlap_matrix
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: tag_z = "OQP::nac_rohf_z"
    character(len=*), parameter :: tag_out = "OQP::nac_rohf_xc_adjoint"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t) :: molgrid

    real(kind=dp), contiguous, pointer :: dma(:), dmb(:), mo(:,:), z(:)
    real(kind=dp), pointer :: out(:,:)
    real(kind=dp), allocatable :: pa(:,:), pb(:,:), pza(:,:), pzb(:,:)
    real(kind=dp), allocatable :: xa(:,:), xb(:,:), work(:,:), half(:,:)
    real(kind=dp), allocatable :: xcd(:,:,:), xcp(:,:,:), gxc(:,:)
    real(kind=dp), allocatable :: vxcza(:,:), vxczb(:,:), vxcmoa(:,:), vxcmob(:,:)
    real(kind=dp), allocatable :: fxca(:,:,:), fxcb(:,:,:), dsa(:,:,:,:)
    real(kind=dp) :: reorth
    integer :: nbf, natom, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: atom, cart, i, j, mu, nu

    if (infos%control%scftype /= 3) then
      call show_message('mrsf_nac_xc_adjoint requires an ROHF/ROKS reference.', &
                        WITH_ABORT)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    natom = infos%mol_prop%natom
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira

    call tagarray_get_data(infos%dat, OQP_DM_A, dma)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmb)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, tag_z, z)
    if (size(z) /= ltot) then
      call show_message('OQP::nac_rohf_z has the wrong ROHF rotation dimension.', &
                        WITH_ABORT)
    end if

    allocate(pa(nbf,nbf), pb(nbf,nbf), pza(nbf,nbf), pzb(nbf,nbf))
    allocate(xa(nvira,nocca), xb(nvirb,noccb), work(nbf,nbf), half(nbf,nbf))
    allocate(xcd(nbf,nbf,2), xcp(nbf,nbf,2), gxc(3,natom), &
             source=0.0_dp)

    call unpack_matrix(dma, pa)
    call unpack_matrix(dmb, pb)
    call rohf_unpack_trial(z, xa, xb, nbf, nocca, noccb)

    ! delta P_alpha = C_v xa C_o^T + transpose.
    half = 0.0_dp
    call dgemm('n','n', nbf, nocca, nvira, 1.0_dp, &
               mo(:,nocca+1:nbf), nbf, xa, nvira, 0.0_dp, half, nbf)
    call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, half, nbf, &
               mo(:,1:nocca), nbf, 0.0_dp, work, nbf)
    pza = work + transpose(work)

    ! delta P_beta = C_(socc+virt) xb C_docc^T + transpose.  The first
    ! offset rows of xb are the socc-docc block; the remaining rows are vd.
    if (noccb > 0) then
      half = 0.0_dp
      call dgemm('n','n', nbf, noccb, nvirb, 1.0_dp, &
                 mo(:,noccb+1:nbf), nbf, xb, nvirb, 0.0_dp, half, nbf)
      call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, half, nbf, &
                 mo(:,1:noccb), nbf, 0.0_dp, work, nbf)
      pzb = work + transpose(work)
    else
      pzb = 0.0_dp
    end if

    if (infos%control%hamilton == 20) then
      ! This public entry may be called independently of the Z solver.  Make
      ! libxc setup deterministic even if an earlier driver left a live
      ! process-global functional list.
      call dftclean(infos)
      call dft_initialize(infos, basis, molgrid, verbose=.false.)

      xcd(:,:,1) = pa
      xcd(:,:,2) = pb
      xcp(:,:,1) = pza
      xcp(:,:,2) = pzb
      call utddft_xc_gradient(basis=basis, molGrid=molgrid, dedft=gxc, &
           da=xcd(:,:,1), db=xcd(:,:,2), pa=xcp(:,:,1:1), pb=xcp(:,:,2:2), &
           nmtx=1, threshold=0.0_dp, infos=infos, &
           include_ground_state=.false., include_weight_derivative=.true.)

      ! The nuclear CPKS skeleton also contains f_xc[D^(R,0)], where
      ! D^(R,0) is the occupied-space reorthonormalization density
      !
      !   D_s^(R,0) = -C_occ,s S^R_occ,s C_occ,s^T .
      !
      ! utddft_xc_gradient above differentiates at fixed AO density and hence
      ! supplies the explicit moving-grid/basis dV_xc/dR part, but not this
      ! coefficient-response term.  By symmetry of the XC kernel,
      !
      !   -1/2 Tr[P_z f_xc[D^(R,0)]]
      !       = -1/2 Tr[D^(R,0) f_xc[P_z]],
      !
      ! so f_xc[P_z] is built only once per adjoint pair and then contracted
      ! with every overlap derivative.  Evaluating the XC kernel directly
      ! avoids double-counting the JK reorthonormalization term already present
      ! in the analytic HF/JK/Pulay right-hand side.
      allocate(vxcza(nbf,nbf), vxczb(nbf,nbf), &
               vxcmoa(nbf,nbf), vxcmob(nbf,nbf), &
               fxca(nbf,nbf,1), fxcb(nbf,nbf,1), &
               dsa(nbf,nbf,3,natom), source=0.0_dp)
      ! get_response_packed forms JK[P_z] + f_xc[P_z].  The previous
      ! implementation built JK[P_z] a second time solely to subtract it
      ! again.  This adjoint needs only f_xc[P_z], so evaluate that linear XC
      ! kernel directly.  Besides removing two ERI traversals, this avoids a
      ! packed JK subtraction and its cancellation roundoff.
      fxca = 0.0_dp
      fxcb = 0.0_dp
      xcp(:,:,1) = pza
      xcp(:,:,2) = pzb
      call utddft_fxc(basis=basis, molGrid=molgrid, isVecs=.true., &
           wfa=mo, wfb=mo, fxa=fxca, fxb=fxcb, &
           dxa=xcp(:,:,1:1), dxb=xcp(:,:,2:2), &
           nMtx=1, threshold=0.0_dp, infos=infos)
      vxcza = fxca(:,:,1)
      vxczb = fxcb(:,:,1)
      call ao_to_mo(vxcza, vxcmoa)
      call ao_to_mo(vxczb, vxcmob)

      call der_overlap_matrix(basis, dsa)
      do atom = 1, natom
        do cart = 1, 3
          do nu = 1, nbf
            do mu = 1, nbf
              dsa(mu,nu,cart,atom) = dsa(mu,nu,cart,atom) * &
                   basis%bfnrm(mu)*basis%bfnrm(nu)
            end do
          end do
          call ao_to_mo(dsa(:,:,cart,atom), work)
          reorth = 0.0_dp
          do j = 1, nocca
            do i = 1, nocca
              reorth = reorth + work(i,j)*vxcmoa(j,i)
            end do
          end do
          do j = 1, noccb
            do i = 1, noccb
              reorth = reorth + work(i,j)*vxcmob(j,i)
            end do
          end do
          ! Tr[P_z f_xc[D^(R,0)]] = -reorth.
          gxc(cart,atom) = gxc(cart,atom) - reorth
        end do
      end do
      call dftclean(infos)

      gxc = -0.5_dp*gxc
      deallocate(vxcza, vxczb, vxcmoa, vxcmob, fxca, fxcb, dsa)
    else
      gxc = 0.0_dp
    end if

    call infos%dat%erase((/ character(len=80) :: tag_out /))
    call tagarray_reserve_data(infos%dat, tag_out, TA_TYPE_REAL64, 3*natom, &
         (/ 3, natom /), comment='native ROHF NAC analytic XC adjoint')
    call tagarray_get_data(infos%dat, tag_out, out)
    out = gxc

    deallocate(pa, pb, pza, pzb, xa, xb, work, half, xcd, xcp, gxc)
  contains
    subroutine ao_to_mo(ao, transformed)
      real(kind=dp), intent(in) :: ao(:,:)
      real(kind=dp), intent(out) :: transformed(:,:)
      half = 0.0_dp
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, ao, nbf, &
                 0.0_dp, half, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, half, nbf, mo, nbf, &
                 0.0_dp, transformed, nbf)
    end subroutine ao_to_mo
  end subroutine mrsf_nac_xc_adjoint

!###############################################################################

!> Contract several native-ROHF adjoint vectors with the analytic nuclear XC
!> stationarity derivative in one molecular-grid traversal.  All probes share
!> the ground-state density, grid, AO values, functional derivatives, moving-
!> grid partition geometry, and overlap derivatives; the probe dimension is
!> retained through every thread-local accumulation and reduction.
  subroutine mrsf_nac_xc_adjoint_batch(infos, z_vectors, gxc_vectors)
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_reserve_data, tagarray_get_data, &
      OQP_DM_A, OQP_DM_B, OQP_VEC_MO_A
    use mathlib, only: unpack_matrix
    use cphf_mod, only: rohf_unpack_trial
    use dft, only: dft_initialize, dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
    use mod_dft_gridint_fxc, only: utddft_fxc
    use grd1, only: der_overlap_matrix
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: z_vectors(:,:)
    real(kind=dp), intent(out) :: gxc_vectors(:,:,:)

    type(basis_set), pointer :: basis
    type(dft_grid_t) :: molgrid
    real(kind=dp), contiguous, pointer :: dma(:), dmb(:), mo(:,:)
    real(kind=dp), allocatable :: pa(:,:), pb(:,:)
    real(kind=dp), allocatable :: pza(:,:,:), pzb(:,:,:)
    real(kind=dp), allocatable :: xa(:,:), xb(:,:), work(:,:), half(:,:)
    real(kind=dp), allocatable :: xcd(:,:,:), xcpa(:,:,:), xcpb(:,:,:)
    real(kind=dp), allocatable :: fxca(:,:,:), fxcb(:,:,:)
    real(kind=dp), allocatable :: vxcmoa(:,:), vxcmob(:,:)
    real(kind=dp), allocatable :: dsa(:,:,:,:), dsa_occ(:,:,:,:)
    real(kind=dp), allocatable :: gxc_sum(:,:)
    real(kind=dp) :: reorth
    integer :: nbf, natom, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: nrhs, irhs, atom, cart, i, j, mu, nu

    if (infos%control%scftype /= 3) then
      call show_message( &
        'mrsf_nac_xc_adjoint_batch requires an ROHF/ROKS reference.', &
        WITH_ABORT)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    natom = infos%mol_prop%natom
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira
    nrhs = size(z_vectors,2)
    if (nrhs < 1 .or. size(z_vectors,1) /= ltot .or. &
        size(gxc_vectors,1) /= 3 .or. &
        size(gxc_vectors,2) /= natom .or. &
        size(gxc_vectors,3) /= nrhs) then
      call show_message('Batched XC adjoint dimensions are inconsistent.', &
                        WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dma)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmb)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)

    allocate(pa(nbf,nbf), pb(nbf,nbf), &
             pza(nbf,nbf,nrhs), pzb(nbf,nbf,nrhs), &
             xa(nvira,nocca), xb(nvirb,noccb), &
             work(nbf,nbf), half(nbf,nbf), &
             xcd(nbf,nbf,2), xcpa(nbf,nbf,nrhs), &
             xcpb(nbf,nbf,nrhs), gxc_sum(3,natom), source=0.0_dp)
    call unpack_matrix(dma, pa)
    call unpack_matrix(dmb, pb)

    do irhs = 1, nrhs
      call rohf_unpack_trial(z_vectors(:,irhs), xa, xb, &
                             nbf, nocca, noccb)

      half = 0.0_dp
      call dgemm('n','n', nbf, nocca, nvira, 1.0_dp, &
                 mo(:,nocca+1:nbf), nbf, xa, nvira, &
                 0.0_dp, half, nbf)
      call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, half, nbf, &
                 mo(:,1:nocca), nbf, 0.0_dp, work, nbf)
      pza(:,:,irhs) = work + transpose(work)

      if (noccb > 0) then
        half = 0.0_dp
        call dgemm('n','n', nbf, noccb, nvirb, 1.0_dp, &
                   mo(:,noccb+1:nbf), nbf, xb, nvirb, &
                   0.0_dp, half, nbf)
        call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, half, nbf, &
                   mo(:,1:noccb), nbf, 0.0_dp, work, nbf)
        pzb(:,:,irhs) = work + transpose(work)
      else
        pzb(:,:,irhs) = 0.0_dp
      end if
    end do

    gxc_vectors = 0.0_dp
    if (infos%control%hamilton == 20) then
      call dftclean(infos)
      call dft_initialize(infos, basis, molgrid, verbose=.false.)

      xcd(:,:,1) = pa
      xcd(:,:,2) = pb
      xcpa = pza
      xcpb = pzb
      call utddft_xc_gradient(basis=basis, molGrid=molgrid, &
           dedft=gxc_sum, da=xcd(:,:,1), db=xcd(:,:,2), &
           pa=xcpa, pb=xcpb, nMtx=nrhs, threshold=0.0_dp, infos=infos, &
           include_ground_state=.false., &
           include_weight_derivative=.true., dedft_mtx=gxc_vectors)

      ! The XC-only coefficient-response kernel is linear in every P_z.
      ! utddft_fxc already carries an nMtx axis, so evaluate all probes while
      ! sharing the AO values and functional derivatives of each grid slice.
      allocate(fxca(nbf,nbf,nrhs), fxcb(nbf,nbf,nrhs), &
               vxcmoa(nocca,nocca), vxcmob(noccb,noccb), &
               dsa(nbf,nbf,3,natom), &
               dsa_occ(nocca,nocca,3,natom), source=0.0_dp)
      xcpa = pza
      xcpb = pzb
      call utddft_fxc(basis=basis, molGrid=molgrid, isVecs=.true., &
           wfa=mo, wfb=mo, fxa=fxca, fxb=fxcb, dxa=xcpa, dxb=xcpb, &
           nMtx=nrhs, threshold=0.0_dp, infos=infos)

      ! S^R and its AO->MO transform are independent of the state pair.  Only
      ! the occupied block enters the trace, so use rectangular DGEMMs and
      ! retain that block once instead of forming nrhs full MO matrices.
      call der_overlap_matrix(basis, dsa)
      do atom = 1, natom
        do cart = 1, 3
          do nu = 1, nbf
            do mu = 1, nbf
              dsa(mu,nu,cart,atom) = dsa(mu,nu,cart,atom) * &
                   basis%bfnrm(mu)*basis%bfnrm(nu)
            end do
          end do
          call ao_to_mo_occ(dsa(:,:,cart,atom), nocca, &
                            dsa_occ(:,:,cart,atom))
        end do
      end do

      do irhs = 1, nrhs
        call ao_to_mo_occ(fxca(:,:,irhs), nocca, vxcmoa)
        call ao_to_mo_occ(fxcb(:,:,irhs), noccb, vxcmob)
        do atom = 1, natom
          do cart = 1, 3
            reorth = 0.0_dp
            do j = 1, nocca
              do i = 1, nocca
                reorth = reorth + &
                  dsa_occ(i,j,cart,atom)*vxcmoa(j,i)
              end do
            end do
            do j = 1, noccb
              do i = 1, noccb
                reorth = reorth + &
                  dsa_occ(i,j,cart,atom)*vxcmob(j,i)
              end do
            end do
            gxc_vectors(cart,atom,irhs) = &
              gxc_vectors(cart,atom,irhs) - reorth
          end do
        end do
      end do
      gxc_vectors = -0.5_dp*gxc_vectors
      call dftclean(infos)

      deallocate(fxca, fxcb, vxcmoa, vxcmob, dsa, dsa_occ)
    end if

    deallocate(pa, pb, pza, pzb, xa, xb, work, half, xcd, xcpa, xcpb, &
               gxc_sum)
  contains
    subroutine ao_to_mo_occ(ao, nocc, transformed)
      real(kind=dp), intent(in) :: ao(:,:)
      integer, intent(in) :: nocc
      real(kind=dp), intent(out) :: transformed(:,:)
      if (nocc == 0) return
      half = 0.0_dp
      call dgemm('t','n', nocc, nbf, nbf, 1.0_dp, mo, nbf, ao, nbf, &
                 0.0_dp, half, nbf)
      call dgemm('n','n', nocc, nocc, nbf, 1.0_dp, half, nbf, mo, nbf, &
                 0.0_dp, transformed, nocc)
    end subroutine ao_to_mo_occ
  end subroutine mrsf_nac_xc_adjoint_batch

end module mrsf_nac_interchange_mod
