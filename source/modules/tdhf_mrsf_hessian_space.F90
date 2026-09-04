module tdhf_mrsf_hessian_space_mod

  use precision, only: dp

  implicit none

  private
  public :: build_mrsf_packed_transform
  public :: mrsf_physical_dimensions

contains

!###############################################################################

  pure subroutine mrsf_physical_dimensions(nbf,nocca,noccb,multiplicity, &
                                             expanded,physical,status)
    integer, intent(in) :: nbf,nocca,noccb,multiplicity
    integer, intent(out) :: expanded,physical,status

    expanded = 0
    physical = 0
    status = 0
    if (nbf <= 0 .or. noccb < 0 .or. nocca > nbf .or. &
        nocca-noccb /= 2) then
      status = -1
      return
    end if
    expanded = nocca*(nbf-noccb)
    select case(multiplicity)
    case(1)
      physical = expanded-1
    case(3)
      physical = expanded-3
    case default
      status = -2
    end select
    if (physical <= 0) status = -3
  end subroutine mrsf_physical_dimensions

!###############################################################################

  pure subroutine build_mrsf_packed_transform(nbf,nocca,noccb, &
                                               multiplicity,transform,status)
    ! Selection from the independent physical singlet/triplet response space
    ! into the packed Davidson vector consumed by mrsfcbc and returned by
    ! mrsfmntoia.  The packed L coefficient is already the normalized physical
    ! OO amplitude; mrsfcbc performs the L/R spin expansion exactly once.  The
    ! four OO positions are
    !
    !   L=(O1,O1), G=(O2,O1), D=(O1,O2), R=(O2,O2).
    !
    ! Singlet packed space: L, G, and D retained; R absent.
    ! Triplet packed space: L retained; G, D, and R absent.

    integer, intent(in) :: nbf,nocca,noccb,multiplicity
    real(kind=dp), intent(out) :: transform(:,:)
    integer, intent(out) :: status

    integer :: column,expanded,physical,position
    integer :: index_l,index_g,index_d,index_r

    call mrsf_physical_dimensions(nbf,nocca,noccb,multiplicity, &
                                   expanded,physical,status)
    transform = 0.0_dp
    if (status /= 0) return
    if (any(shape(transform) /= [expanded,physical])) then
      status = -4
      return
    end if

    index_l = (nocca-noccb-2)*nocca+nocca-1
    index_g = (nocca-noccb-2)*nocca+nocca
    index_d = (nocca-noccb-1)*nocca+nocca-1
    index_r = (nocca-noccb-1)*nocca+nocca
    if (index_l < 1 .or. index_r > expanded) then
      status = -5
      return
    end if

    column = 0
    do position=1,expanded
      if (multiplicity == 1) then
        if (position == index_r) cycle
        column = column+1
        transform(position,column) = 1.0_dp
      else
        if (position == index_g .or. position == index_d .or. &
            position == index_r) cycle
        column = column+1
        transform(position,column) = 1.0_dp
      end if
    end do
    if (column /= physical) status = -6
  end subroutine build_mrsf_packed_transform

end module tdhf_mrsf_hessian_space_mod
