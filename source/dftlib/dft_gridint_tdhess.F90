! Grid-level entry points needed by the restricted TDDFT Hessian driver.
module mod_dft_gridint_tdhess

  use precision, only: fp
  use basis_tools, only: basis_set
  use mod_dft_molgrid, only: dft_grid_t
  use types, only: information
  use mod_dft_gridint_fxc, only: tddft_fxc
  use mod_dft_gridint_gxc, only: tddft_gxc
  use mod_dft_gridint_tdxc_grad, only: tddft_xc_gradient

  implicit none
  private

  public :: tdxc_fxc_action
  public :: tdxc_kxc_diagonal_action
  public :: tdxc_response_rows

contains

  subroutine tdxc_fxc_action(basis, grid, ground, density, result, threshold, infos)
    type(basis_set), intent(inout) :: basis
    type(dft_grid_t), target, intent(in) :: grid
    real(fp), intent(in) :: ground(:,:)
    real(fp), contiguous, target, intent(inout) :: density(:,:,:)
    real(fp), intent(inout) :: result(:,:,:)
    real(fp), intent(in) :: threshold
    type(information), target, intent(in) :: infos

    call require_restricted_lda(infos, 'tdxc_fxc_action')
    call tddft_fxc(basis, grid, .false., ground, result, density, &
                   size(density,3), threshold, infos)
  end subroutine tdxc_fxc_action

  ! Quadratic third-kernel action Gxc[D,D].  Cross terms Gxc[A,B] are obtained
  ! from three calls with A+B, A and B using lda_tdxc_polarized_kxc.
  subroutine tdxc_kxc_diagonal_action(basis, grid, ground, density, result, &
                                      threshold, infos)
    type(basis_set), intent(inout) :: basis
    type(dft_grid_t), target, intent(in) :: grid
    real(fp), intent(in) :: ground(:,:)
    real(fp), contiguous, target, intent(inout) :: density(:,:,:)
    real(fp), intent(inout) :: result(:,:,:)
    real(fp), intent(in) :: threshold
    type(information), target, intent(in) :: infos

    call require_restricted_lda(infos, 'tdxc_kxc_diagonal_action')
    call tddft_gxc(basis, grid, .false., ground, result, density, &
                   size(density,3), threshold, infos)
  end subroutine tdxc_kxc_diagonal_action

  ! Forms one XC nuclear-gradient row per supplied density matrix.  With xa
  ! present this is the derivative of vxc:Peff + 4 fxc:X^2 and therefore the
  ! fixed-density XC Hessian row.  Replacing pa/xa by their first-order
  ! response densities supplies the GAMESS DXR response row.
  subroutine tdxc_response_rows(basis, grid, ground, peff, xpy, rows, &
                                threshold, infos)
    type(basis_set), intent(inout) :: basis
    type(dft_grid_t), target, intent(in) :: grid
    real(fp), contiguous, target, intent(inout) :: ground(:,:)
    real(fp), target, intent(inout) :: peff(:,:,:)
    real(fp), target, intent(inout) :: xpy(:,:,:)
    real(fp), intent(out) :: rows(:,:)
    real(fp), intent(in) :: threshold
    type(information), target, intent(in) :: infos

    call require_restricted_lda(infos, 'tdxc_response_rows')
    if (size(peff,3) /= size(xpy,3)) error stop &
      'tdxc_response_rows: Peff/X+Y batch mismatch'
    call tddft_xc_gradient(basis, grid, rows, ground, peff, xpy, &
                           size(peff,3), threshold, infos)
  end subroutine tdxc_response_rows

  subroutine require_restricted_lda(infos, caller)
    type(information), intent(in) :: infos
    character(*), intent(in) :: caller
    if (infos%functional%needTau) &
      error stop trim(caller)//': meta-GGA response derivatives are not implemented'
  end subroutine require_restricted_lda

end module mod_dft_gridint_tdhess
