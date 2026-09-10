module tdhf_mrsf_hessian_mo_response_mod

  use precision, only: dp
  use tdhf_mrsf_hessian_one_e_action_mod, only: &
    transform_mo_operator_derivative

  implicit none

  private
  public :: build_mrsf_mo_fock_derivatives
  public :: canonicalize_rohf_common_response

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

!###############################################################################

  subroutine canonicalize_rohf_common_response(mo,dmo,fock_a_ao,fock_b_ao, &
      dfock_a_ao,dfock_b_ao,nocca,noccb,status,gap_tolerance)
    ! Complete the common ROHF orbital derivative in the same canonical gauge
    ! as the stored Guest--Saunders orbitals.  CPHF determines rotations
    ! between the C, O, and V sectors; the antisymmetric rotations within each
    ! sector follow by differentiating the effective-Fock eigenvectors.
    real(kind=dp), intent(in) :: mo(:,:),fock_a_ao(:,:),fock_b_ao(:,:), &
      dfock_a_ao(:,:,:),dfock_b_ao(:,:,:)
    real(kind=dp), intent(inout) :: dmo(:,:,:)
    integer, intent(in) :: nocca,noccb
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: gap_tolerance

    real(kind=dp), allocatable :: fa_mo(:,:),fb_mo(:,:),dfa_mo(:,:,:), &
      dfb_mo(:,:,:),rotation(:,:)
    real(kind=dp) :: gap,tolerance,off_diagonal
    integer :: nbf,ncoord,coordinate,sector,p,q,local_status
    integer :: first(3),last(3)

    nbf=size(mo,1)
    ncoord=size(dmo,3)
    status=0
    tolerance=1.0e-8_dp
    if(present(gap_tolerance)) tolerance=gap_tolerance
    if(nbf<=0 .or. ncoord<=0 .or. noccb<0 .or. noccb>=nocca .or. &
       nocca>=nbf .or. tolerance<=0.0_dp .or. &
       any(shape(mo)/=[nbf,nbf]) .or. &
       any(shape(fock_a_ao)/=[nbf,nbf]) .or. &
       any(shape(fock_b_ao)/=[nbf,nbf]) .or. &
       any(shape(dfock_a_ao)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dfock_b_ao)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if

    allocate(fa_mo(nbf,nbf),fb_mo(nbf,nbf), &
      dfa_mo(nbf,nbf,ncoord),dfb_mo(nbf,nbf,ncoord),rotation(nbf,nbf), &
      source=0.0_dp)
    fa_mo=matmul(transpose(mo),matmul(fock_a_ao,mo))
    fb_mo=matmul(transpose(mo),matmul(fock_b_ao,mo))
    call transform_mo_operator_derivative(mo,dmo,fock_a_ao,dfock_a_ao, &
      dfa_mo,local_status)
    if(local_status==0) call transform_mo_operator_derivative(mo,dmo, &
      fock_b_ao,dfock_b_ao,dfb_mo,local_status)
    if(local_status/=0) then
      status=-2
      deallocate(fa_mo,fb_mo,dfa_mo,dfb_mo,rotation)
      return
    end if

    first=[1,noccb+1,nocca+1]
    last=[noccb,nocca,nbf]
    do coordinate=1,ncoord
      rotation=0.0_dp
      do sector=1,3
        do q=first(sector)+1,last(sector)
          do p=first(sector),q-1
            gap=0.5_dp*((fa_mo(p,p)+fb_mo(p,p))- &
              (fa_mo(q,q)+fb_mo(q,q)))
            off_diagonal=0.25_dp*(dfa_mo(p,q,coordinate)+ &
              dfa_mo(q,p,coordinate)+dfb_mo(p,q,coordinate)+ &
              dfb_mo(q,p,coordinate))
            if(abs(gap)<=tolerance) then
              if(abs(off_diagonal)>10.0_dp*tolerance) then
                status=-3
                deallocate(fa_mo,fb_mo,dfa_mo,dfb_mo,rotation)
                return
              end if
            else
              rotation(p,q)=-off_diagonal/gap
              rotation(q,p)=-rotation(p,q)
            end if
          end do
        end do
      end do
      dmo(:,:,coordinate)=dmo(:,:,coordinate)+matmul(mo,rotation)
    end do
    deallocate(fa_mo,fb_mo,dfa_mo,dfb_mo,rotation)
  end subroutine canonicalize_rohf_common_response

end module tdhf_mrsf_hessian_mo_response_mod
