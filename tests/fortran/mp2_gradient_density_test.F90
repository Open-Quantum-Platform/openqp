! Test adapter for the production MP2 amplitude-to-AO density contraction.
module mp2_gradient_density_test_mod
  use iso_c_binding, only: c_int64_t, c_double
  use mp2_gradient_mod, only: l2_to_ao_2rdm
  implicit none
  private
  public :: mp2_gradient_density_test
contains
  subroutine mp2_gradient_density_test(n, no, c, l2, gao) bind(C, name="mp2_gradient_density_test")
    integer(c_int64_t), value, intent(in) :: n, no
    real(c_double), intent(in) :: c(n,n), l2(no,no,n-no,n-no)
    real(c_double), intent(out) :: gao(n,n,n,n)
    call l2_to_ao_2rdm(n, no, c, l2, gao)
  end subroutine mp2_gradient_density_test
end module mp2_gradient_density_test_mod
