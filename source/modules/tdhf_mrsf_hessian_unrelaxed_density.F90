module tdhf_mrsf_hessian_unrelaxed_density_mod

  use precision, only: dp
  use types, only: information
  use tdhf_lib, only: iatogen
  use tdhf_sf_lib, only: sfdmat
  use tdhf_mrsf_lib, only: mrsfxvec,mrsfcbc
  use tdhf_mrsf_hessian_orbital_maps_mod, only: &
    build_mrsf_density_orbital_derivative

  implicit none

  private
  public :: build_mrsf_unrelaxed_density_derivatives
  public :: build_mrsf_mo_difference_density_derivatives

contains

!###############################################################################

  subroutine build_mrsf_mo_difference_density_derivatives(infos,x_packed, &
      dx_packed,tij,tab,dtij,dtab,status)
    ! MO-basis unrelaxed difference-density blocks and their analytic nuclear
    ! derivatives.  The spin-adapted physical amplitude is expanded exactly
    ! once, and the homogeneous quadratic maps are differentiated directly.

    type(information), intent(in) :: infos
    real(kind=dp), intent(in) :: x_packed(:),dx_packed(:,:)
    real(kind=dp), intent(out) :: tij(:,:),tab(:,:),dtij(:,:,:),dtab(:,:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: x_expanded(:),dx_expanded(:), &
      x_matrix(:,:),dx_matrix(:,:)
    integer :: coordinate,nbf,nocca,noccb,nvirb,ncoord,packed

    nbf=infos%basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    nvirb=nbf-noccb
    ncoord=size(dx_packed,2)
    packed=nocca*nvirb
    status=0
    tij=0.0_dp
    tab=0.0_dp
    dtij=0.0_dp
    dtab=0.0_dp
    if(nbf<=0 .or. ncoord<=0 .or. nocca-noccb/=2 .or. &
       size(x_packed)/=packed .or. any(shape(dx_packed)/=[packed,ncoord]) .or. &
       any(shape(tij)/=[nocca,nocca]) .or. any(shape(tab)/=[nvirb,nvirb]) .or. &
       any(shape(dtij)/=[nocca,nocca,ncoord]) .or. &
       any(shape(dtab)/=[nvirb,nvirb,ncoord])) then
      status=-1
      return
    end if
    if(infos%tddft%mult/=1 .and. infos%tddft%mult/=3) then
      status=-2
      return
    end if

    allocate(x_expanded(packed),dx_expanded(packed), &
      x_matrix(nocca,nvirb),dx_matrix(nocca,nvirb))
    call mrsfxvec(infos,x_packed,x_expanded)
    x_matrix=reshape(x_expanded,[nocca,nvirb])
    tij=-matmul(x_matrix,transpose(x_matrix))
    tab=matmul(transpose(x_matrix),x_matrix)
    do coordinate=1,ncoord
      call mrsfxvec(infos,dx_packed(:,coordinate),dx_expanded)
      dx_matrix=reshape(dx_expanded,[nocca,nvirb])
      dtij(:,:,coordinate)=-matmul(dx_matrix,transpose(x_matrix)) &
        -matmul(x_matrix,transpose(dx_matrix))
      dtab(:,:,coordinate)=matmul(transpose(dx_matrix),x_matrix) &
        +matmul(transpose(x_matrix),dx_matrix)
    end do
    deallocate(x_expanded,dx_expanded,x_matrix,dx_matrix)
  end subroutine build_mrsf_mo_difference_density_derivatives

!###############################################################################

  subroutine build_mrsf_unrelaxed_density_derivatives(infos,mo_a,mo_b, &
      dmo_a,dmo_b,x_packed,dx_packed,dchannel,dabxc,dta,dtb,status)
    ! Analytic first derivative of every unrelaxed state-density object used by
    ! the MRSF gradient.  The physical singlet/triplet amplitude is expanded by
    ! mrsfxvec before sfdmat, so the response remains spin adapted throughout.
    ! Centered evaluations below are algebraic polarization
    ! identities in C or X separately; they are not displaced-geometry finite
    ! differences.  Separating the variables removes all mixed cubic terms.

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: x_packed(:),dx_packed(:,:)
    real(kind=dp), intent(out) :: dchannel(:,:,:,:),dabxc(:,:,:)
    real(kind=dp), intent(out) :: dta(:,:),dtb(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: x_expanded(:),dx_expanded(:),dx_matrix(:,:), &
      orbital_channel(:,:,:,:),amplitude_channel(:,:,:), &
      plus_a(:,:),minus_a(:,:),plus_abxc(:,:),minus_abxc(:,:),work_abxc(:,:), &
      plus_ta(:),minus_ta(:),plus_tb(:),minus_tb(:),work_ta(:),work_tb(:)
    integer :: coordinate,local_status,nbf,nocca,noccb,ncoord,packed,nbf_tri

    nbf=infos%basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    ncoord=size(dmo_a,3)
    packed=nocca*(nbf-noccb)
    nbf_tri=nbf*(nbf+1)/2
    status=0
    dchannel=0.0_dp
    dabxc=0.0_dp
    dta=0.0_dp
    dtb=0.0_dp
    if(nbf<=0 .or. ncoord<=0 .or. nocca-noccb/=2 .or. &
       size(x_packed)/=packed .or. any(shape(dx_packed)/=[packed,ncoord]) .or. &
       any(shape(mo_a)/=[nbf,nbf]) .or. any(shape(mo_b)/=[nbf,nbf]) .or. &
       any(shape(dmo_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dmo_b)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dchannel)/=[7,nbf,nbf,ncoord]) .or. &
       any(shape(dabxc)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dta)/=[nbf_tri,ncoord]) .or. &
       any(shape(dtb)/=[nbf_tri,ncoord])) then
      status=-1
      return
    end if
    if(infos%tddft%mult/=1 .and. infos%tddft%mult/=3) then
      status=-2
      return
    end if

    allocate(x_expanded(packed),dx_expanded(packed),dx_matrix(nbf,nbf), &
      orbital_channel(7,nbf,nbf,ncoord), &
      amplitude_channel(7,nbf,nbf),plus_a(nbf,nbf),minus_a(nbf,nbf), &
      plus_abxc(nbf,nbf),minus_abxc(nbf,nbf),work_abxc(nbf,nbf), &
      plus_ta(nbf_tri),minus_ta(nbf_tri),plus_tb(nbf_tri), &
      minus_tb(nbf_tri),work_ta(nbf_tri),work_tb(nbf_tri),source=0.0_dp)
    call mrsfxvec(infos,x_packed,x_expanded)
    ! The stationary MRSF gradient constructs all seven channels with the
    ! common ROHF response orbital set (mo_a,mo_a), see
    ! build_mrsf_zvector_rhs.  The alpha and beta semicanonical coefficient
    ! arrays need not have the same nuclear connection, even though they span
    ! the same ROHF spaces.  Differentiating (mo_a,mo_b) here would therefore
    ! be the derivative of the Davidson sigma map, not the derivative of the
    ! density actually contracted by the gradient.
    call build_mrsf_density_orbital_derivative(infos,mo_a,mo_a,dmo_a,dmo_a, &
      x_packed,orbital_channel,local_status)
    if(local_status/=0) then
      status=-3
      deallocate(x_expanded,dx_expanded,dx_matrix,orbital_channel, &
        amplitude_channel,plus_a,minus_a,plus_abxc,minus_abxc,work_abxc, &
        plus_ta,minus_ta,plus_tb,minus_tb,work_ta,work_tb)
      return
    end if

    do coordinate=1,ncoord
      call mrsfxvec(infos,dx_packed(:,coordinate),dx_expanded)
      ! The seven-channel mrsfcbc map consumes the stored physical MRSF
      ! response vector, exactly as build_mrsf_zvector_rhs does.  The expanded
      ! vector is reserved for sfdmat and the relaxed-density maps below.
      call iatogen(dx_packed(:,coordinate),dx_matrix,nocca,noccb)
      amplitude_channel=0.0_dp
      call mrsfcbc(infos,mo_a,mo_a,dx_matrix,amplitude_channel)
      dchannel(:,:,:,coordinate)=orbital_channel(:,:,:,coordinate)+ &
        amplitude_channel

      ! Linear-in-X contribution to d(abxc)/dR.
      work_abxc=0.0_dp
      call sfdmat(dx_expanded,work_abxc,mo_a,work_ta,work_tb,nocca,noccb)
      dabxc(:,:,coordinate)=work_abxc

      ! Quadratic-in-C contribution at fixed X.  The same calls also provide
      ! the orbital part of the two quadratic difference-density derivatives.
      plus_a=mo_a+dmo_a(:,:,coordinate)
      minus_a=mo_a-dmo_a(:,:,coordinate)
      plus_abxc=0.0_dp
      minus_abxc=0.0_dp
      call sfdmat(x_expanded,plus_abxc,plus_a,plus_ta,plus_tb,nocca,noccb)
      call sfdmat(x_expanded,minus_abxc,minus_a,minus_ta,minus_tb,nocca,noccb)
      dabxc(:,:,coordinate)=dabxc(:,:,coordinate)+ &
        0.5_dp*(plus_abxc-minus_abxc)
      dta(:,coordinate)=0.5_dp*(plus_ta-minus_ta)
      dtb(:,coordinate)=0.5_dp*(plus_tb-minus_tb)

      ! Quadratic-in-X contribution at fixed C.  A centered polarization is
      ! exact for this homogeneous quadratic map.
      plus_abxc=0.0_dp
      minus_abxc=0.0_dp
      call sfdmat(x_expanded+dx_expanded,plus_abxc,mo_a,plus_ta,plus_tb, &
        nocca,noccb)
      call sfdmat(x_expanded-dx_expanded,minus_abxc,mo_a,minus_ta,minus_tb, &
        nocca,noccb)
      dta(:,coordinate)=dta(:,coordinate)+0.5_dp*(plus_ta-minus_ta)
      dtb(:,coordinate)=dtb(:,coordinate)+0.5_dp*(plus_tb-minus_tb)

      ! The gradient replaces the raw mrsfcbc ball channel by td_abxc before
      ! every response-integral contraction.  Its derivative must therefore
      ! obey the same alias exactly; retaining the derivative of the discarded
      ! raw ball would both violate the stationary-gradient data model and
      ! double count the ordinary response-exchange contribution downstream.
      dchannel(7,:,:,coordinate)=dabxc(:,:,coordinate)
    end do
    deallocate(x_expanded,dx_expanded,dx_matrix,orbital_channel, &
      amplitude_channel,plus_a,minus_a,plus_abxc,minus_abxc,work_abxc, &
      plus_ta,minus_ta,plus_tb,minus_tb,work_ta,work_tb)
  end subroutine build_mrsf_unrelaxed_density_derivatives

end module tdhf_mrsf_hessian_unrelaxed_density_mod
