module xi_fock_selftest_mod
!> @brief Finite-difference check of the xi^alpha Kohn-Sham matrix.
!>   For the converged SCF density P the grid driver returns E_xc[P] and
!>   F = dE_xc/dP.  The test compares the directional derivative
!>   (E_xc[P+eps*D] - E_xc[P-eps*D])/(2 eps) with sum_ij F_ij D_ij for a fixed
!>   symmetric direction D, once with the ordinary tau (xi_mode = 0, the
!>   calibration of the packed-matrix convention) and once with xi^alpha in the
!>   tau slot (xi_mode = 1, order taken from infos%dft%xi_alpha).  The two
!>   ratios FD/analytic must agree, which proves that the xi^alpha Fock
!>   contribution (V_tau D_i.D_j) is the exact derivative of the energy that the
!>   same code evaluates.  Diagnostic harness only.
!>   Results go to the caller's res(1:8) (no shared file, so concurrent pytest
!>   workers cannot overwrite each other):
!>     res = [FD_tau, analytic_tau, ratio_tau, FD_xi, analytic_xi, ratio_xi,
!>            relative ratio mismatch, 1 (pass) / 0 (fail)]
  implicit none
contains

  subroutine xi_fock_selftest_C(c_handle, res) bind(C, name="xi_fock_selftest")
    use, intrinsic :: iso_c_binding, only: c_double
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    real(c_double), intent(out) :: res(8)
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call xi_fock_selftest(inf, res)
  end subroutine xi_fock_selftest_C

  subroutine xi_fock_selftest(infos, res)
    use types, only: information
    use precision, only: dp
    use mod_dft, only: dft_initialize
    use mod_dft_molgrid, only: dft_grid_t
    use scf_addons, only: calc_dft_xc_density, scf_rhf
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_DM_B

    type(information), target, intent(inout) :: infos
    real(dp), intent(out) :: res(8)

    type(dft_grid_t) :: molgrid
    real(dp), contiguous, pointer :: dm_a(:), dm_b(:)
    real(dp), allocatable :: dmat(:,:), pfxc(:,:), delta(:,:), dpert(:,:)
    real(dp) :: e0, ep, em, ne, ek, eps, fd(0:1), an(0:1), ratio(0:1), w
    integer :: nbf, ntri, nspin, p, i, j, m, iu
    integer(8) :: xi_mode_save
    logical :: urohf, ok

    associate(basis => infos%basis)
      nbf = basis%nbf
      ntri = nbf*(nbf + 1)/2
      urohf = infos%control%scftype /= scf_rhf
      nspin = merge(2, 1, urohf)

      ! grid only: the functional attached by the SCF is reused (no second libxc instance)
      call dft_initialize(infos, basis, molgrid, need_functional=.false.)

      call tagarray_get_data(infos%dat, OQP_DM_A, dm_a)
      allocate(dmat(ntri, nspin), delta(ntri, nspin), dpert(ntri, nspin), pfxc(ntri, nspin))
      dmat(:, 1) = dm_a(1:ntri)
      if (urohf) then
        call tagarray_get_data(infos%dat, OQP_DM_B, dm_b)
        dmat(:, 2) = dm_b(1:ntri)
      end if

      ! deterministic, bounded symmetric direction (packed lower triangle)
      p = 0
      do i = 1, nbf
        do j = 1, i
          p = p + 1
          delta(p, 1) = 0.01_dp*sin(0.37_dp*p + 0.11_dp*i)
          if (nspin == 2) delta(p, 2) = 0.01_dp*cos(0.23_dp*p + 0.05_dp*j)
        end do
      end do

      xi_mode_save = infos%dft%xi_mode
      eps = 1.0e-4_dp

      do m = 0, 1
        infos%dft%xi_mode = m
        call calc_dft_xc_density(infos, basis, molgrid, dmat, pfxc, e0, ne, ek)
        ! sum_ij F_ij D_ij over the full matrix from packed storage
        an(m) = 0.0_dp
        do iu = 1, nspin
          p = 0
          do i = 1, nbf
            do j = 1, i
              p = p + 1
              w = merge(1.0_dp, 2.0_dp, i == j)
              an(m) = an(m) + w*pfxc(p, iu)*delta(p, iu)
            end do
          end do
        end do
        dpert = dmat + eps*delta
        call calc_dft_xc_density(infos, basis, molgrid, dpert, pfxc, ep, ne, ek)
        dpert = dmat - eps*delta
        call calc_dft_xc_density(infos, basis, molgrid, dpert, pfxc, em, ne, ek)
        fd(m) = (ep - em)/(2.0_dp*eps)
        ratio(m) = fd(m)/an(m)
      end do
      infos%dft%xi_mode = xi_mode_save

      ok = abs(ratio(1) - ratio(0)) < 1.0e-6_dp*abs(ratio(0)) .and. abs(an(1)) > 1.0e-12_dp

      res = [fd(0), an(0), ratio(0), fd(1), an(1), ratio(1), &
             abs(ratio(1) - ratio(0))/max(abs(ratio(0)), 1.0e-300_dp), merge(1.0_dp, 0.0_dp, ok)]
    end associate
  end subroutine xi_fock_selftest

end module xi_fock_selftest_mod
