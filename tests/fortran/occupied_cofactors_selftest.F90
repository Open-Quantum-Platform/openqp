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
  call mrsf_tlf(infos, full_s, actual, ab, ia, 0)
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
