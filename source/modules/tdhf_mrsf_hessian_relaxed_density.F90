module tdhf_mrsf_hessian_relaxed_density_mod

  use precision, only: dp
  use tdhf_sf_lib, only: sfropcal

  implicit none

  private
  public :: build_mrsf_relaxed_density_derivatives

contains

!###############################################################################

  subroutine build_mrsf_relaxed_density_derivatives(mo_a,mo_b,dmo_a,dmo_b, &
      tij,tab,dtij,dtab,z,dz,nocca,noccb,dpa,dpb,status)
    ! Derivative of the relaxed spin densities used by the MRSF gradient.
    ! sfropcal is linear in the unrelaxed density blocks and the orbital
    ! Lagrange multiplier, while the AO transformation is differentiated
    ! explicitly.  The output follows the symmetrization and one-half scaling
    ! of OQP_td_p in tdhf_mrsf_z_vector.

    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: tij(:,:),tab(:,:),dtij(:,:,:),dtab(:,:,:)
    real(kind=dp), intent(in) :: z(:),dz(:,:)
    integer, intent(in) :: nocca,noccb
    real(kind=dp), intent(out) :: dpa(:,:,:),dpb(:,:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: pmo_a(:,:),pmo_b(:,:),dpmo_a(:,:), &
      dpmo_b(:,:),draw(:,:)
    integer :: coordinate,nbf,ncoord,nvirb,lzdim

    nbf=size(mo_a,1)
    ncoord=size(dmo_a,3)
    nvirb=nbf-noccb
    lzdim=noccb*(nocca-noccb+nbf-nocca)+(nocca-noccb)*(nbf-nocca)
    status=0
    dpa=0.0_dp
    dpb=0.0_dp
    if(nbf<=0 .or. ncoord<=0 .or. nocca-noccb/=2 .or. &
       any(shape(mo_a)/=[nbf,nbf]) .or. any(shape(mo_b)/=[nbf,nbf]) .or. &
       any(shape(dmo_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dmo_b)/=[nbf,nbf,ncoord]) .or. &
       any(shape(tij)/=[nocca,nocca]) .or. any(shape(tab)/=[nvirb,nvirb]) .or. &
       any(shape(dtij)/=[nocca,nocca,ncoord]) .or. &
       any(shape(dtab)/=[nvirb,nvirb,ncoord]) .or. size(z)/=lzdim .or. &
       any(shape(dz)/=[lzdim,ncoord]) .or. &
       any(shape(dpa)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dpb)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if

    allocate(pmo_a(nbf,nbf),pmo_b(nbf,nbf),dpmo_a(nbf,nbf), &
      dpmo_b(nbf,nbf),draw(nbf,nbf),source=0.0_dp)
    call sfropcal(pmo_a,pmo_b,tij,tab,z,nocca,noccb)
    do coordinate=1,ncoord
      call sfropcal(dpmo_a,dpmo_b,dtij(:,:,coordinate), &
        dtab(:,:,coordinate),dz(:,coordinate),nocca,noccb)

      draw=matmul(dmo_a(:,:,coordinate),matmul(pmo_a,transpose(mo_a)))+ &
        matmul(mo_a,matmul(dpmo_a,transpose(mo_a)))+ &
        matmul(mo_a,matmul(pmo_a,transpose(dmo_a(:,:,coordinate))))
      dpa(:,:,coordinate)=0.5_dp*(draw+transpose(draw))

      draw=matmul(dmo_b(:,:,coordinate),matmul(pmo_b,transpose(mo_b)))+ &
        matmul(mo_b,matmul(dpmo_b,transpose(mo_b)))+ &
        matmul(mo_b,matmul(pmo_b,transpose(dmo_b(:,:,coordinate))))
      dpb(:,:,coordinate)=0.5_dp*(draw+transpose(draw))
    end do
    deallocate(pmo_a,pmo_b,dpmo_a,dpmo_b,draw)
  end subroutine build_mrsf_relaxed_density_derivatives

end module tdhf_mrsf_hessian_relaxed_density_mod
