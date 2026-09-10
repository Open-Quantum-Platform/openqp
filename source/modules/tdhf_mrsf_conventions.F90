module tdhf_mrsf_conventions_mod

  implicit none

  private
  public :: mrsf_raw_spc_multiplier

contains

!###############################################################################

  pure integer function mrsf_raw_spc_multiplier(multiplicity) result(sign)
    ! The six special AO transition-density channels stored by OpenQP define
    ! K_raw=-C_SPC^physical.  The seventh (ball) channel belongs to A0 and is
    ! not spin-pairing coupling.  Consequently the founding physical equations
    ! A(S)=A0-C_SPC and A(T)=A0+C_SPC require +K_raw for a singlet and -K_raw
    ! for a triplet.
    integer, intent(in) :: multiplicity
    select case(multiplicity)
    case(1)
      sign = 1
    case(3)
      sign = -1
    case default
      sign = 0
    end select
  end function mrsf_raw_spc_multiplier

end module tdhf_mrsf_conventions_mod
