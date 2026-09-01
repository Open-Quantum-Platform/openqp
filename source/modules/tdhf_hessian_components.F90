module tdhf_hessian_components_mod

  use precision, only: dp

  implicit none

  private
  public :: assemble_tdhf_cartesian_hessian
  public :: tdhf_hessian_is_applicable
  public :: tdhf_hessian_functional_is_verified
  public :: tdhf_hessian_target_is_supported
  public :: tdhf_hessian_lowest_root_is_isolated

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
                                                    tamm_dancoff, is_dft, &
                                                    verified_local_density, &
                                                    needs_gradient, needs_tau, &
                                                    range_separated, mpi_size) &
                                                    result(applicable)
    ! Applicability of the imported GAMESS formulation at its verified limit.
    ! OpenQP SCF type 1 is the closed-shell restricted reference.  The initial
    ! import is restricted to singlet full-response TDHF and restricted
    ! LDA/GGA or global-hybrid GGA TDDFT on one MPI rank.  Meta-GGA and
    ! range-separated functionals must fail closed.

    integer, intent(in) :: scf_type
    integer, intent(in) :: td_multiplicity
    logical, intent(in) :: tamm_dancoff
    logical, intent(in) :: is_dft
    logical, intent(in) :: verified_local_density
    logical, intent(in) :: needs_gradient
    logical, intent(in) :: needs_tau
    logical, intent(in) :: range_separated
    integer, intent(in) :: mpi_size

    applicable = scf_type == 1 .and. td_multiplicity == 1 .and. &
                 .not. tamm_dancoff .and. mpi_size == 1 .and. &
                 (.not. is_dft .or. &
                  (verified_local_density .and. &
                   .not. needs_tau .and. &
                   .not. range_separated))

  end function tdhf_hessian_is_applicable

!###############################################################################

  pure logical function tdhf_hessian_functional_is_verified(name) result(verified)
    character(*), intent(in) :: name
    character(len=len(name)) :: upper_name
    integer :: i, code

    upper_name = adjustl(name)
    do i = 1, len_trim(upper_name)
      code = iachar(upper_name(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) &
        upper_name(i:i) = achar(code - iachar('a') + iachar('A'))
    end do
    verified = trim(upper_name) == 'SVWN' .or. &
               trim(upper_name) == 'SVWN5' .or. &
               trim(upper_name) == 'LDA' .or. &
               trim(upper_name) == 'BLYP' .or. &
               trim(upper_name) == 'PBE' .or. &
               trim(upper_name) == 'PBEPBE' .or. &
               trim(upper_name) == 'B3LYP5' .or. &
               trim(upper_name) == 'B3LYPV5'
  end function tdhf_hessian_functional_is_verified

!###############################################################################

  pure logical function tdhf_hessian_target_is_supported(target) result(supported)
    integer, intent(in) :: target
    supported = target == 1
  end function tdhf_hessian_target_is_supported

!###############################################################################

  pure logical function tdhf_hessian_lowest_root_is_isolated(energies,tol) result(isolated)
    real(kind=dp), intent(in) :: energies(:),tol
    isolated = size(energies)>=1 .and. tol>=0.0_dp
    if(isolated .and. size(energies)>1) &
      isolated=abs(energies(2)-energies(1))>tol
  end function tdhf_hessian_lowest_root_is_isolated

end module tdhf_hessian_components_mod
