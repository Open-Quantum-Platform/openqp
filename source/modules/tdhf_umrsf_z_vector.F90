!> UMRSF-TDDFT Z-vector (CPKS) — clean-room implementation (branch uhf-grad-plan).
!> STAGE: pipeline-wiring STUB. Sets Z_Vector_converged=.true. so the gradient driver proceeds.
!> The real CPKS solve (RHS from the unrelaxed difference density via umrsfcbc->int2(J/K,spc,HFscale)
!> ->umrsfmntoia; coupled alpha-alpha/beta-beta/alpha-beta; zvconv=1e-10; reuse zvector_common) is
!> implemented next per DERIVATIONS/stage1_hf_gradient_spec.md. Never reads the guarded RO z-vector.
module tdhf_umrsf_z_vector_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_umrsf_z_vector_mod"

contains

  subroutine tdhf_umrsf_z_vector_C(c_handle) bind(C, name="tdhf_umrsf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use io_constants, only: iw
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    open(unit=iw, file=inf%log_filename, position="append")
    write(iw,'(/2x,a)') 'UMRSF Z-vector: STUB (pipeline wiring only; CPKS not yet implemented)'
    close(iw)
    inf%mol_energy%Z_Vector_converged = .true.
  end subroutine tdhf_umrsf_z_vector_C

end module tdhf_umrsf_z_vector_mod
