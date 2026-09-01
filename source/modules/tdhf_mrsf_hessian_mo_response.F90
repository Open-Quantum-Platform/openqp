module tdhf_mrsf_hessian_mo_response_mod

  use precision, only: dp
  use tdhf_mrsf_hessian_one_e_action_mod, only: &
    transform_mo_operator_derivative

  implicit none

  private
  public :: build_mrsf_mo_fock_derivatives

contains

!###############################################################################

  subroutine build_mrsf_mo_fock_derivatives(mo,dmo,fock_ao,dfock_ao, &
      dfock_mo,orbital_energy_derivative,status)
    ! Moving-basis total derivative of C^T F C and its diagonal.  Hiroya
    ! Nakata's TDHF/TDDFT Hessian construction is the methodological starting
    ! point; the two-SOMO MRSF extension keeps alpha and beta calls separate.

    real(kind=dp), intent(in) :: mo(:,:),dmo(:,:,:),fock_ao(:,:),dfock_ao(:,:,:)
    real(kind=dp), intent(out) :: dfock_mo(:,:,:), &
      orbital_energy_derivative(:,:)
    integer, intent(out) :: status

    integer :: coordinate,nbf,ncoord,orbital

    nbf=size(mo,1)
    ncoord=size(dmo,3)
    status=0
    dfock_mo=0.0_dp
    orbital_energy_derivative=0.0_dp
    if(nbf<=0 .or. ncoord<=0 .or. any(shape(mo)/=[nbf,nbf]) .or. &
       any(shape(fock_ao)/=[nbf,nbf]) .or. &
       any(shape(dmo)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dfock_ao)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dfock_mo)/=[nbf,nbf,ncoord]) .or. &
       any(shape(orbital_energy_derivative)/=[nbf,ncoord])) then
      status=-1
      return
    end if
    call transform_mo_operator_derivative(mo,dmo,fock_ao,dfock_ao, &
      dfock_mo,status)
    if(status/=0) return
    do coordinate=1,ncoord
      do orbital=1,nbf
        orbital_energy_derivative(orbital,coordinate)= &
          dfock_mo(orbital,orbital,coordinate)
      end do
    end do
  end subroutine build_mrsf_mo_fock_derivatives

end module tdhf_mrsf_hessian_mo_response_mod
