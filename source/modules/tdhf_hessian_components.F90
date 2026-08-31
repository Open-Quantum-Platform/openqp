module tdhf_hessian_components_mod

  use precision, only: dp

  implicit none

  private
  public :: assemble_tdhf_cartesian_hessian
  public :: tdhf_hessian_is_applicable

contains

!###############################################################################

  pure subroutine assemble_tdhf_cartesian_hessian(h_fixed, h_xc, response_rows, &
                                                    h_total, row_asymmetry, status)
    ! Assemble the Cartesian Hessian of one TDHF/TDDFT state.
    !
    ! The excitation-energy contribution is written in the stationary-
    ! Lagrangian form
    !
    !   d2 E_I = H_fixed + H_xc + 1/2 (G + transpose(G)),
    !
    ! where H_fixed contains the nuclear, one-electron, Pulay, ECP, and
    ! two-electron second-derivative contractions with the total
    ! state-specific relaxed densities.  H_xc contains the corresponding
    ! exchange-correlation quadrature contribution, and G contains the
    ! first-order response rows.  Because H_fixed already includes the
    ! ground-state contribution, it must not be added separately here.
    !
    ! Keep the directional response rows unsymmetrized until this final sum.
    ! Their antisymmetric part is a useful diagnostic of an incomplete or
    ! unconverged response solution and must not be discarded earlier.

    real(kind=dp), intent(in) :: h_fixed(:,:)
    real(kind=dp), intent(in) :: h_xc(:,:)
    real(kind=dp), intent(in) :: response_rows(:,:)
    real(kind=dp), intent(out) :: h_total(:,:)
    real(kind=dp), intent(out) :: row_asymmetry
    integer, intent(out) :: status

    integer :: ncart

    ncart = size(h_total, 1)
    if (size(h_total, 2) /= ncart .or. &
        any(shape(h_fixed) /= [ncart, ncart]) .or. &
        any(shape(h_xc) /= [ncart, ncart]) .or. &
        any(shape(response_rows) /= [ncart, ncart])) then
      h_total = 0.0_dp
      row_asymmetry = 0.0_dp
      status = -1
      return
    end if

    status = 0

    if (ncart == 0) then
      row_asymmetry = 0.0_dp
    else
      row_asymmetry = maxval(abs(response_rows - transpose(response_rows)))
    end if

    h_total = h_fixed + h_xc &
              + 0.5_dp*(response_rows + transpose(response_rows))

  end subroutine assemble_tdhf_cartesian_hessian

!###############################################################################

  pure logical function tdhf_hessian_is_applicable(scf_type, td_multiplicity, &
                                                    tamm_dancoff, needs_tau, &
                                                    mpi_size) result(applicable)
    ! Applicability of the imported GAMESS formulation at its verified limit.
    ! OpenQP SCF type 1 is the closed-shell restricted reference.  The initial
    ! import is restricted to singlet full-response TDHF/TDDFT, functionals
    ! without kinetic-energy-density dependence, and one MPI rank.

    integer, intent(in) :: scf_type
    integer, intent(in) :: td_multiplicity
    logical, intent(in) :: tamm_dancoff
    logical, intent(in) :: needs_tau
    integer, intent(in) :: mpi_size

    applicable = scf_type == 1 .and. td_multiplicity == 1 .and. &
                 .not. tamm_dancoff .and. .not. needs_tau .and. mpi_size == 1

  end function tdhf_hessian_is_applicable

end module tdhf_hessian_components_mod
