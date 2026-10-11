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
  real(dp) :: ab(0,0), ia(n,0), direct(n,n)
  logical :: batch
  integer :: i, j

  infos%basis%nbf = int(n)
  infos%mol_prop%nelec_a = int(n)
  infos%mol_prop%nelec_b = int(n)
  call mrsf_tlf(infos, s, actual, ab, ia, 0)
  call occupied_cofactors(s, direct, batch)
  if (any(actual /= direct)) error stop 'mrsf_tlf did not return occupied cofactors'
  used_lu = 0_c_int
  if (batch) used_lu = 1_c_int
  do j = 1, n
    do i = 1, n
      call ov_exact(previous(i,j), i, j, 0, 0, s, 1, int(n)-1, 1)
    end do
  end do
end subroutine occupied_cofactors_selftest
