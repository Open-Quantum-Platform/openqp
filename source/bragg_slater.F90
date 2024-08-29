!> @brief Bragg-Slater radii for determining the relative size of the
!>        polyhedra in the polyatomic integration scheme
module bragg_slater_radii
  use precision, only: dp

  private
  public BRSL_NUM_ELEMENTS
  public BRSL_TYPE_GAMESS
  public BRSL_TYPE_GILL
  public BRSL_TYPE_TA
  public set_bragg_slater

  integer, parameter :: BRSL_NUM_ELEMENTS  = 137
  integer, parameter :: BRSL_TYPE_GAMESS = 0
  integer, parameter :: BRSL_TYPE_GILL = 1
  integer, parameter :: BRSL_TYPE_TA = 2

! J.C.Slater, Quantum Theory of Molecules and Solids, Volume 2, Chapter 3
! Except that hydrogen is changed from 0.25 -> bohr radius,
! and missing values such as inert gasses are filled in with
! reasonable looking data (source unknown). Slater's table
! stops at the element americium, the extension is probably
! reasonable for actinides but not all the way to z=137!
  real(kind=dp), parameter :: brsl_values_gamess(BRSL_NUM_ELEMENTS) = [&
    0.52917D+00, 0.31D+00, 1.45D+00, 1.05D+00, 0.85D+00, &
    0.70D+00, 0.65D+00, 0.60D+00, 0.50D+00, 0.38D+00, 1.80D+00, &
    1.50D+00, 1.25D+00, 1.10D+00, 1.00D+00, 1.00D+00, 1.00D+00, &
    0.71D+00, 2.20D+00, 1.80D+00, 1.60D+00, 1.40D+00, 1.35D+00, &
    1.40D+00, 1.40D+00, 1.40D+00, 1.35D+00, 1.35D+00, 1.35D+00, &
    1.35D+00, 1.30D+00, 1.25D+00, 1.15D+00, 1.15D+00, 1.15D+00, &
    0.88D+00, 2.35D+00, 2.00D+00, 1.80D+00, 1.55D+00, 1.45D+00, &
    1.45D+00, 1.35D+00, 1.30D+00, 1.35D+00, 1.40D+00, 1.60D+00, &
    1.55D+00, 1.55D+00, 1.45D+00, 1.45D+00, 1.40D+00, 1.40D+00, &
    1.08D+00, 2.60D+00, 2.15D+00, 1.95D+00, 1.85D+00, 1.85D+00, &
    1.85D+00, 1.85D+00, 1.85D+00, 1.85D+00, 1.80D+00, 1.75D+00, &
    1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
    1.55D+00, 1.45D+00, 1.35D+00, 1.35D+00, 1.30D+00, 1.35D+00, &
    1.35D+00, 1.35D+00, 1.50D+00, 1.90D+00, 1.80D+00, 1.60D+00, &
    1.90D+00, 1.27D+00, 1.20D+00, 2.60D+00, 2.15D+00, 1.95D+00, &
    1.80D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75d0, &
    1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, &
    1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, &
    1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, &
    1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, &
    1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, &
    1.75d0, 1.75d0]

! Gill et al., CPL,209, 506 (1993).
  real(kind=dp), parameter :: brsl_values_gill(BRSL_NUM_ELEMENTS) = [ &
     0.52918D+00, 0.31126D+00, 1.62822D+00, 1.08550D+00, &
     0.81414D+00, 0.65131D+00, 0.54272D+00, 0.46520D+00, &
     0.40704D+00, 0.36185D+00, 2.16481D+00, 1.67109D+00, &
     1.36073D+00, 1.14763D+00, 0.99221D+00, 0.87388D+00, &
     0.78075D+00, 0.70555D+00, 2.20D+00,    1.80D+00, &
     1.60D+00, 1.40D+00, 1.35D+00, 1.40D+00, 1.40D+00, &
     1.40D+00, 1.35D+00, 1.35D+00, 1.35D+00, 1.35D+00, &
     1.30D+00, 1.25D+00, 1.15D+00, 1.15D+00, 1.15D+00, &
     0.88D+00, 2.35D+00, 2.00D+00, 1.80D+00, 1.55D+00, &
     1.45D+00, 1.45D+00, 1.35D+00, 1.30D+00, 1.35D+00, &
     1.40D+00, 1.60D+00, 1.55D+00, 1.55D+00, 1.45D+00, &
     1.45D+00, 1.40D+00, 1.40D+00, 1.08D+00, 2.60D+00, &
     2.15D+00, 1.95D+00, 1.85D+00, 1.85D+00, 1.85D+00, &
     1.85D+00, 1.85D+00, 1.85D+00, 1.80D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.55D+00, 1.45D+00, 1.35D+00, 1.35D+00, &
     1.30D+00, 1.35D+00, 1.35D+00, 1.35D+00, 1.50D+00, &
     1.90D+00, 1.80D+00, 1.60D+00, 1.90D+00, 1.27D+00, &
     1.20D+00, 2.60D+00, 2.15D+00, 1.95D+00, 1.80D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, 1.75D+00, &
     1.75D+00, 1.75D+00]

! O.Treutler, R.Ahlrichs, JCP,102(1), 346 (1995).
  real(kind=dp), parameter :: brsl_values_treutler(BRSL_NUM_ELEMENTS) = [&
     0.8D+00, 0.9D+00, &
     1.8D+00, 1.4D+00, 1.3D+00, 1.1D+00, &
     0.9D+00, 0.9D+00, 0.9D+00, 0.9D+00, &
     1.4D+00, 1.3D+00, 1.3D+00, 1.2D+00, &
     1.1D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.5D+00, 1.4D+00, &
     1.3D+00, 1.2D+00, 1.2D+00, 1.2D+00, 1.2D+00, &
     1.2D+00, 1.2D+00, 1.1D+00, 1.1D+00, 1.1D+00, &
     1.1D+00, 1.0D+00, 0.9D+00, 0.9D+00, 0.9D+00, &
     0.9D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
     1.0D+00]

contains

  subroutine set_bragg_slater(array, bstype)
    real(kind=dp), intent(out) :: array(:)
    integer, intent(in) :: bstype

    select case (bstype)
    case (BRSL_TYPE_GAMESS)
      array = brsl_values_gamess

    case (BRSL_TYPE_TA)
      array = brsl_values_treutler

    case (BRSL_TYPE_GILL)
      array = brsl_values_gill

    case default
      array = brsl_values_gill

    end select
  end subroutine

end module
