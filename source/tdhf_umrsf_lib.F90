
module tdhf_umrsf_lib
! Thin U‑MRSF helpers that delegate to the canonical MRSF kernels
! in OpenQP's tdhf_mrsf_lib. These routines already support VA≠VB (UHF).
!
! Public entry points:
!   subroutine umrsfcbc(infos, va, vb, bvec, fmrsf)
!   subroutine umrsfmntoia(infos, fmrsf, pmo, va, vb, ivec)
!
  use precision, only: dp
  implicit none
  private
  public :: umrsfcbc, umrsfmntoia

contains

  subroutine umrsfcbc(infos, va, vb, bvec, fmrsf)
    use tdhf_mrsf_lib, only: mrsfcbc
    use types, only: information
    type(information), intent(in) :: infos
    real(dp), intent(in)  :: va(:,:), vb(:,:), bvec(:,:)
    real(dp), intent(inout), target :: fmrsf(:,:,:)
    call mrsfcbc(infos, va, vb, bvec, fmrsf)
  end subroutine umrsfcbc

  subroutine umrsfmntoia(infos, fmrsf, pmo, va, vb, ivec)
    use tdhf_mrsf_lib, only: mrsfmntoia
    use types, only: information
    type(information), intent(in) :: infos
    real(dp), intent(in),  target :: fmrsf(:,:,:)
    real(dp), intent(out) :: pmo(:,:)
    real(dp), intent(in)  :: va(:,:), vb(:,:)
    integer, intent(in)   :: ivec
    call mrsfmntoia(infos, fmrsf, pmo, va, vb, ivec)
  end subroutine umrsfmntoia

end module tdhf_umrsf_lib
