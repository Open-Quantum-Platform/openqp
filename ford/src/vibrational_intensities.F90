module vibrational_intensities_mod

  use iso_c_binding, only: c_double, c_int64_t, c_ptr, c_f_pointer

  implicit none

  private
  public :: vibrational_intensities_native_C

  real(c_double), parameter :: IR_INTENSITY_CONVERSION = 42.255d0

contains

  subroutine vibrational_intensities_native_C(c_handle, nmode, ncoord, modes_ptr, dipole_derivs_ptr, &
      polar_derivs_ptr, ir_ptr, mode_dipoles_ptr, raman_ptr, mode_polars_ptr) &
      bind(C, name="vibrational_intensities_native")
    use c_interop, only: oqp_handle_t
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t), value :: nmode, ncoord
    type(c_ptr), value :: modes_ptr, dipole_derivs_ptr, polar_derivs_ptr
    type(c_ptr), value :: ir_ptr, mode_dipoles_ptr, raman_ptr, mode_polars_ptr

    real(c_double), pointer :: modes(:), dipole_derivs(:), polar_derivs(:)
    real(c_double), pointer :: ir(:), mode_dipoles(:), raman(:), mode_polars(:)
    integer(c_int64_t) :: imode, icoord, a, b, p, mode_offset, dip_offset, polar_offset
    real(c_double) :: alpha_prime(3,3), alpha_bar_prime, gamma2

    ! c_handle is present to match the normal OpenQP C ABI wrapper convention.
    associate(unused => c_handle)
    end associate

    call c_f_pointer(modes_ptr, modes, [nmode*ncoord])
    call c_f_pointer(dipole_derivs_ptr, dipole_derivs, [3_c_int64_t*ncoord])
    call c_f_pointer(polar_derivs_ptr, polar_derivs, [9_c_int64_t*ncoord])
    call c_f_pointer(ir_ptr, ir, [nmode])
    call c_f_pointer(mode_dipoles_ptr, mode_dipoles, [nmode*3_c_int64_t])
    call c_f_pointer(raman_ptr, raman, [nmode])
    call c_f_pointer(mode_polars_ptr, mode_polars, [nmode*9_c_int64_t])

    ir = 0.0_c_double
    raman = 0.0_c_double
    mode_dipoles = 0.0_c_double
    mode_polars = 0.0_c_double

    do imode = 1, nmode
      mode_offset = (imode - 1_c_int64_t) * ncoord

      do p = 1, 3
        do icoord = 1, ncoord
          dip_offset = (p - 1_c_int64_t) * ncoord + icoord
          mode_dipoles((imode - 1_c_int64_t) * 3_c_int64_t + p) = &
              mode_dipoles((imode - 1_c_int64_t) * 3_c_int64_t + p) + &
              dipole_derivs(dip_offset) * modes(mode_offset + icoord)
        end do
      end do
      ir(imode) = IR_INTENSITY_CONVERSION * sum(mode_dipoles((imode - 1_c_int64_t) * 3_c_int64_t + 1: &
                                                             (imode - 1_c_int64_t) * 3_c_int64_t + 3)**2)

      alpha_prime = 0.0_c_double
      do a = 1, 3
        do b = 1, 3
          do icoord = 1, ncoord
            polar_offset = ((a - 1_c_int64_t) * 3_c_int64_t + (b - 1_c_int64_t)) * ncoord + icoord
            alpha_prime(a,b) = alpha_prime(a,b) + polar_derivs(polar_offset) * modes(mode_offset + icoord)
          end do
          mode_polars((imode - 1_c_int64_t) * 9_c_int64_t + (a - 1_c_int64_t) * 3_c_int64_t + b) = &
              alpha_prime(a,b)
        end do
      end do

      alpha_bar_prime = (alpha_prime(1,1) + alpha_prime(2,2) + alpha_prime(3,3)) / 3.0_c_double
      gamma2 = 0.5_c_double * ((alpha_prime(1,1) - alpha_prime(2,2))**2 + &
                               (alpha_prime(2,2) - alpha_prime(3,3))**2 + &
                               (alpha_prime(3,3) - alpha_prime(1,1))**2 + &
                               6.0_c_double * (alpha_prime(1,2)**2 + alpha_prime(1,3)**2 + alpha_prime(2,3)**2))
      raman(imode) = 45.0_c_double * alpha_bar_prime**2 + 7.0_c_double * gamma2
    end do
  end subroutine vibrational_intensities_native_C

end module vibrational_intensities_mod
