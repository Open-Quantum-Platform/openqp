
module tdhf_umrsf_lib
  use precision,   only: dp
  implicit none
  private
  public :: umrsfcbc, umrsfmntoia

contains

  ! Thin wrappers delegating to the canonical MRSF library.
  ! These work for both ROHF and UHF provided you pass distinct VA/VB.
  subroutine umrsfcbc(infos, va, vb, bvec, fmrsf)
    use tdhf_mrsf_lib, only: mrsfcbc
    use types,        only: information
    type(information), intent(in)                      :: infos
    real(dp), intent(in),  dimension(:,:)              :: va, vb, bvec
    real(dp), intent(inout), target, dimension(:,:,:)  :: fmrsf
    call mrsfcbc(infos, va, vb, bvec, fmrsf)
  end subroutine umrsfcbc

  subroutine umrsfmntoia(infos, fmrsf, pmo, va, vb, ivec)
    use tdhf_mrsf_lib, only: mrsfmntoia
    use types,        only: information
    type(information), intent(in)                     :: infos
    real(dp), intent(in),  target, dimension(:,:,:)  :: fmrsf
    real(dp), intent(out), dimension(:,:)            :: pmo
    real(dp), intent(in),  dimension(:,:)            :: va, vb
    integer, intent(in)                               :: ivec
    call mrsfmntoia(infos, fmrsf, pmo, va, vb, ivec)
  end subroutine umrsfmntoia

end module tdhf_umrsf_lib
