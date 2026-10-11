! Test-only C shim: exercise the actual mrsf_tlf call and the previous minors.
subroutine occupied_cofactors_selftest(n, s, actual, previous, used_lu) bind(C)
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_int, c_double
  use precision, only: dp
  use types, only: information
  use get_state_overlap_mod, only: mrsf_tlf, ov_exact, occupied_cofactors
  implicit none
  integer(c_int64_t), value :: n
  real(c_double), intent(in) :: s(n,n)
  real(c_double), intent(out) :: actual(n,n), previous(n,n)
  integer(c_int), intent(out) :: used_lu
  type(information) :: infos
  real(dp), allocatable :: full_s(:,:), ab(:,:), ia(:,:)
  real(dp) :: direct(n,n), expected
  logical :: batch
  integer :: i, j, nbf, nvir

  nbf = int(n)
  nvir = 0
  if (n >= 2) then
    nbf = int(n) + 2
    nvir = 4
  end if
  allocate(full_s(nbf,nbf), ab(nvir,nvir), ia(n,nvir))
  do j = 1, nbf
    do i = 1, nbf
      full_s(i,j) = 0.003_dp * (2*i-j)
    end do
  end do
  full_s(1:n,1:n) = s
  infos%basis%nbf = nbf
  infos%mol_prop%nelec_a = int(n)
  infos%mol_prop%nelec_b = nbf-nvir
  ! mrsf_tlf may create debug records; this local information owns them.
  call infos%dat%new('occupied_cofactors_selftest')
  if (.not. infos%dat%is_allocated()) error stop 'Cannot allocate selftest container'
  call mrsf_tlf(infos, full_s, actual, ab, ia, 0)
  ! Returned arrays are caller-owned and independent of the debug records.
  call infos%dat%delete()
  call occupied_cofactors(s, direct, batch)
  if (any(actual /= direct)) error stop 'mrsf_tlf did not return occupied cofactors'
  used_lu = 0_c_int
  if (batch) used_lu = 1_c_int
  do j = 1, n
    do i = 1, n
      call ov_exact(previous(i,j), i, j, 0, 0, s, 1, int(n)-1, 1)
    end do
  end do
  do j = 1, nvir
    do i = 1, nvir
      call ov_exact(expected, 0, 0, nbf-nvir+i, nbf-nvir+j, full_s, 1, int(n)-1, 2)
      if (ab(i,j) /= expected) error stop 'LU batch changed s_ab'
    end do
    do i = 1, n
      call ov_exact(expected, i, 0, nbf-nvir+j, 0, full_s, 1, int(n)-1, 3)
      if (ia(i,j) /= expected) error stop 'LU batch changed s_ia'
    end do
  end do
end subroutine occupied_cofactors_selftest

! Test-only C shim: bypass size dispatch while retaining every numerical guard.
subroutine occupied_cofactors_batch_selftest(n, s, actual, previous, used_lu) bind(C)
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_int, c_double
  use precision, only: dp
  use get_state_overlap_mod, only: occupied_cofactors_lu, ov_exact
  implicit none
  integer(c_int64_t), value :: n
  real(c_double), intent(in) :: s(n,n)
  real(c_double), intent(out) :: actual(n,n), previous(n,n)
  integer(c_int), intent(out) :: used_lu
  logical :: batch
  integer :: i, j

  call occupied_cofactors_lu(s, actual, batch)
  used_lu = 0_c_int
  if (batch) used_lu = 1_c_int
  do j = 1, n
    do i = 1, n
      call ov_exact(previous(i,j), i, j, 0, 0, s, 1, int(n)-1, 1)
    end do
  end do
end subroutine occupied_cofactors_batch_selftest

! One information/TagArray owner across grow/shrink and rejected record creation.
subroutine occupied_cofactors_persistent_selftest(nmax, nsteps, sizes, s, actual, previous, checks, status) bind(C)
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_int32_t, c_int, c_double, c_ptr, c_associated
  use precision, only: dp
  use types, only: information
  use tagarray, only: recordinfo_t, TA_CONTAINER_RECORD_EXISTS
  use oqp_tagarray_driver, only: tagarray_get_data, TA_TYPE_REAL64, TA_OK
  use get_state_overlap_mod, only: mrsf_tlf, ov_exact
  implicit none
  integer(c_int64_t), value :: nmax, nsteps
  integer(c_int64_t), intent(in) :: sizes(nsteps)
  real(c_double), intent(in) :: s(nmax,nmax,nsteps)
  real(c_double), intent(out) :: actual(nmax,nmax,nsteps), previous(nmax,nmax,nsteps)
  integer(c_int), intent(out) :: checks(6), status
  type(information) :: infos
  type(c_ptr) :: owner
  real(dp), allocatable :: full_s(:,:), ab(:,:), ia(:,:)
  integer(c_int32_t) :: record_status
  integer :: step, n, nbf, i, j, allocation_status
  logical :: ok

  actual = 0.0_dp
  previous = 0.0_dp
  checks = 0_c_int
  status = 1_c_int
  ! Debug records have nonzero extents; no empty-record behavior is claimed.
  if (nmax < 2 .or. nsteps < 1) goto 900
  if (any(sizes < 2) .or. any(sizes > nmax)) goto 900
  status = 2_c_int
  call infos%dat%new('occupied_cofactors_persistent_selftest')
  if (.not. infos%dat%is_allocated()) goto 900
  owner = infos%dat%container_ptr
  record_status = infos%dat%create('OQP::persistent_owner_probe', TA_TYPE_REAL64, [3_c_int64_t])
  if (record_status /= TA_OK) goto 900
  block
    real(dp), pointer :: probe(:)
    nullify(probe)
    call tagarray_get_data(infos%dat, 'OQP::persistent_owner_probe', probe, record_status)
    if (record_status /= TA_OK) goto 900
    if (.not. associated(probe)) goto 900
    if (size(probe) /= 3) goto 900
    probe = [3.0_dp, -2.0_dp, 5.0_dp]
    nullify(probe)
  end block
  checks(1) = 1_c_int
  do step = 1, int(nsteps)
    n = int(sizes(step))
    nbf = n + 2
    allocate(full_s(nbf,nbf), ab(4,4), ia(n,4), stat=allocation_status)
    if (allocation_status /= 0) goto 900
    do j = 1, nbf
      do i = 1, nbf
        full_s(i,j) = 0.003_dp * (2*i-j)
      end do
    end do
    full_s(1:n,1:n) = s(1:n,1:n,step)
    infos%basis%nbf = nbf
    infos%mol_prop%nelec_a = n
    infos%mol_prop%nelec_b = nbf - 4
    ! Only records are replaced. The owning container stays alive throughout.
    call mrsf_tlf(infos, full_s, actual(1:n,1:n,step), ab, ia, 0)
    status = 3_c_int
    if (.not. c_associated(infos%dat%container_ptr, owner)) goto 900
    call check_record('OQP::persistent_owner_probe', [3.0_dp, -2.0_dp, 5.0_dp], ok)
    if (.not. ok) goto 900
    checks(2) = checks(2) + 1_c_int
    do j = 1, n
      do i = 1, n
        call ov_exact(previous(i,j,step), i, j, 0, 0, full_s, 1, n-1, 1)
      end do
    end do
    status = 4_c_int
    call check_record('OQP::dbg_s_ij', reshape(actual(1:n,1:n,step), [n*n]), ok)
    if (.not. ok) goto 900
    call check_record('OQP::dbg_s_ab', reshape(ab, [16]), ok)
    if (.not. ok) goto 900
    call check_record('OQP::dbg_s_ia', reshape(ia, [4*n]), ok)
    if (.not. ok) goto 900
    call check_record('OQP::dbg_s_mo', reshape(full_s, [nbf*nbf]), ok)
    if (.not. ok) goto 900
    checks(3) = checks(3) + 4_c_int
    ! A real backend error: creating an existing tag without override fails.
    ! Its borrowed view has already been nullified by check_record.
    record_status = infos%dat%create('OQP::dbg_s_ij', TA_TYPE_REAL64, [int(n*n,c_int64_t)])
    status = 5_c_int
    if (record_status /= TA_CONTAINER_RECORD_EXISTS) goto 900
    checks(4) = checks(4) + 1_c_int
    status = 6_c_int
    call check_record('OQP::dbg_s_ij', reshape(actual(1:n,1:n,step), [n*n]), ok)
    if (.not. ok) goto 900
    checks(5) = checks(5) + 1_c_int
    deallocate(full_s, ab, ia)
    status = 2_c_int
  end do
  status = 0_c_int
900 continue
  ! All borrowed pointers were scoped/nullified before any erase/delete.
  call infos%dat%delete()
  if (infos%dat%is_allocated()) status = 7_c_int
  checks(6) = checks(6) + 1_c_int
  call infos%dat%delete()
  if (infos%dat%is_allocated()) status = 7_c_int
  checks(6) = checks(6) + 1_c_int
contains
  subroutine check_record(tag, expected, matches)
    character(len=*), intent(in) :: tag
    real(dp), intent(in) :: expected(:)
    logical, intent(out) :: matches
    type(recordinfo_t) :: record
    real(dp), pointer :: borrowed(:)
    integer(c_int32_t) :: get_status
    matches = .false.
    nullify(borrowed)
    if (.not. infos%dat%contains(tag)) return
    record = infos%dat%get(tag)
    if (record%type_id /= TA_TYPE_REAL64 .or. record%ndims /= 1) return
    if (record%count /= size(expected,kind=c_int64_t)) return
    if (.not. associated(record%dims)) return
    if (record%dims(1) /= size(expected,kind=c_int64_t)) return
    call tagarray_get_data(infos%dat, tag, borrowed, get_status)
    if (get_status == TA_OK) then
      if (associated(borrowed)) matches = all(borrowed == expected)
    end if
    nullify(borrowed)
  end subroutine check_record
end subroutine occupied_cofactors_persistent_selftest
