!> UMRSF-TDDFT analytic nuclear gradient — clean-room implementation (branch uhf-grad-plan).
!> STAGE: pipeline-wiring STUB. Zeroes the gradient so the dispatch runs end-to-end; the real
!> assembly (reference UHF gradient + relaxed difference density P^Δ + energy-weighted W via grd1,
!> 2-particle density Γ via grd2 UHF driver, + M1 two-reference term) is implemented next per
!> DERIVATIONS/stage1_hf_gradient_spec.md, FD-validated by harness/fdgrad.py. Reuses generic
!> derivative-integral primitives (grd1/grd2); never reads the guarded RO-MRSF gradient.
module tdhf_umrsf_gradient_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_umrsf_gradient_mod"

contains

  subroutine tdhf_umrsf_gradient_C(c_handle) bind(C, name="tdhf_umrsf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use io_constants, only: iw
    use printing, only: print_module_info
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    open(unit=iw, file=inf%log_filename, position="append")
    call print_module_info('UMRSF_TDHF_Gradient','UMRSF-TDDFT Gradient (STUB: pipeline wiring)')
    write(iw,'(/2x,a)') 'UMRSF gradient: STUB (zeroes gradient; analytic assembly not yet implemented)'
    close(iw)
    inf%atoms%grad = 0.0d0
  end subroutine tdhf_umrsf_gradient_C

end module tdhf_umrsf_gradient_mod
