!> @file mp2_energy.F90
!>
!> @brief Driver for the standalone MP2 ground-state energy method.
!>
!> MP2 is a post-SCF ground-state correlation correction: the SCF reference
!> (RHF/UHF/ROHF) is converged first by the usual PyOQP `reference` step, and
!> this driver adds the second-order Moller-Plesset correlation energy on top,
!> reusing the validated two-electron driver via `mp2_lib`.  It is dispatched
!> from Python as `[tdhf] type = mp2` (a ground-state post-SCF method that
!> reports no excitations).
module mp2_energy_mod

  implicit none

  private
  public :: mp2_energy

  character(len=*), parameter :: module_name = "mp2_energy_mod"

contains

  !> C-bound entry point: `[tdhf] type = mp2` dispatches here.
  subroutine mp2_energy_C(c_handle) bind(C, name="mp2_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call mp2_energy(inf)
  end subroutine mp2_energy_C

  !> MP2 ground-state correlation on the converged SCF reference.
  subroutine mp2_energy(infos)
    use precision, only: dp
    use io_constants, only: iw
    use types, only: information
    use printing, only: print_module_info
    use mp2_lib, only: mp2_correlation

    implicit none

    type(information), target, intent(inout) :: infos
    real(kind=dp) :: e_mp2, e_aa, e_bb, e_ab, e_ref
    logical :: computed

    open(unit=iw, file=infos%log_filename, position="append")
    call print_module_info('MP2_Energy', 'Computing MP2 ground-state correlation')

    e_ref = infos%mol_energy%energy

    call mp2_correlation(infos, e_mp2, e_aa, e_bb, e_ab, computed)

    write(iw,'(/,2X,60("="))')
    write(iw,'(2X,A)') 'MP2  (Moller-Plesset second order, ground state)'
    write(iw,'(2X,60("="))')
    if (computed) then
      write(iw,'(2X,A,F20.10)') 'E(reference, SCF)      = ', e_ref
      write(iw,'(2X,A,F20.10)') 'E(MP2, same-spin aa)   = ', e_aa
      write(iw,'(2X,A,F20.10)') 'E(MP2, same-spin bb)   = ', e_bb
      write(iw,'(2X,A,F20.10)') 'E(MP2, opp-spin  ab)   = ', e_ab
      write(iw,'(2X,A,F20.10)') 'E(MP2, correlation)    = ', e_mp2
      write(iw,'(2X,A,F20.10)') 'E(MP2, total)          = ', e_ref + e_mp2
      ! Report the MP2 total as the molecular energy for downstream consumers.
      infos%mol_energy%energy = e_ref + e_mp2
      infos%mol_energy%etot   = e_ref + e_mp2
    else
      write(iw,'(2X,A)') 'MP2 not computed: the system exceeds the per-MO-pair'
      write(iw,'(2X,A)') 'Coulomb-build size guard (raise OQP_ADC2_MAX_JBUILDS,'
      write(iw,'(2X,A)') 'or use a smaller basis).'
    end if
    write(iw,'(2X,60("="),/)')

    close(iw)

  end subroutine mp2_energy

end module mp2_energy_mod
