module tdhf_mrsf_hessian_amplitude_mod

  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_compute_t
  use tdhf_hessian_response_mod, only: &
    solve_mrsf_tda_response_batch_matrix_free
  use tdhf_mrsf_hessian_space_mod, only: mrsf_physical_dimensions, &
    build_mrsf_packed_transform
  use tdhf_mrsf_sigma_mod, only: apply_mrsf_tda_sigma
  use tdhf_mrsf_lib, only: mrinivec
  use oqp_tagarray_driver, only: tagarray_get_data,OQP_E_MO_A

  implicit none

  private
  public :: solve_mrsf_tda_amplitude_derivatives

contains

!###############################################################################

  subroutine solve_mrsf_tda_amplitude_derivatives(infos,int2_driver,mo_a,mo_b, &
      fa,fb,omega,x_packed,dax_packed,dx_packed,domega,residual_max,status, &
      tolerance,max_iterations,restart)
    ! Solve all nuclear amplitude-response equations in the independent
    ! spin-adapted CO/OV/CV/OO space.  The production packed representation is
    ! used only at the seven-density sigma boundary; forbidden OO slots never
    ! enter the Krylov space.

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),fa(:,:),fb(:,:)
    real(kind=dp), intent(in) :: omega,x_packed(:),dax_packed(:,:)
    real(kind=dp), intent(out) :: dx_packed(:,:),domega(:),residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations,restart

    real(kind=dp), contiguous, pointer :: mo_energy(:)
    real(kind=dp), allocatable :: transform(:,:),x_physical(:), &
      dax_physical(:,:),dx_physical(:,:),packed_diagonal(:), &
      physical_diagonal(:),preconditioner(:),seed_vector(:,:)
    real(kind=dp) :: denominator,solve_tolerance
    integer :: component,nbf,nocca,noccb,packed,physical,ncoord, &
      solve_iterations,space_status

    status=0; residual_max=0.0_dp; dx_packed=0.0_dp; domega=0.0_dp
    nbf=infos%basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    ncoord=size(dax_packed,2)
    call mrsf_physical_dimensions(nbf,nocca,noccb,infos%tddft%mult, &
      packed,physical,space_status)
    if(space_status/=0 .or. ncoord<=0 .or. size(x_packed)/=packed .or. &
       size(dax_packed,1)/=packed .or. &
       any(shape(dx_packed)/=[packed,ncoord]) .or. size(domega)/=ncoord) then
      status=-1
      return
    end if
    if(any(shape(mo_a)/=[nbf,nbf]) .or. any(shape(mo_b)/=[nbf,nbf]) .or. &
       any(shape(fa)/=[nbf,nbf]) .or. any(shape(fb)/=[nbf,nbf])) then
      status=-2
      return
    end if
    allocate(transform(packed,physical),x_physical(physical), &
      dax_physical(physical,ncoord),dx_physical(physical,ncoord), &
      packed_diagonal(packed),physical_diagonal(physical), &
      preconditioner(physical),seed_vector(packed,1))
    call build_mrsf_packed_transform(nbf,nocca,noccb,infos%tddft%mult, &
      transform,space_status)
    if(space_status/=0) then
      status=-3
      deallocate(transform,x_physical,dax_physical,dx_physical, &
        packed_diagonal,physical_diagonal,preconditioner,seed_vector)
      return
    end if
    x_physical=matmul(transpose(transform),x_packed)
    if(maxval(abs(matmul(transform,x_physical)-x_packed))>1.0e-12_dp .or. &
       abs(dot_product(x_physical,x_physical)-1.0_dp)>1.0e-9_dp) then
      status=-4
      deallocate(transform,x_physical,dax_physical,dx_physical, &
        packed_diagonal,physical_diagonal,preconditioner,seed_vector)
      return
    end if
    dax_physical=matmul(transpose(transform),dax_packed)
    call tagarray_get_data(infos%dat,OQP_E_MO_A,mo_energy)
    call mrinivec(infos,mo_energy,mo_energy,seed_vector,packed_diagonal,1, &
      report_symmetry_coverage=.false.)
    physical_diagonal=0.0_dp
    do component=1,physical
      physical_diagonal(component)=sum(packed_diagonal* &
        transform(:,component)*transform(:,component), &
        mask=abs(transform(:,component))>epsilon(1.0_dp))
      denominator=physical_diagonal(component)-omega
      if(.not.(abs(denominator)<huge(denominator))) then
        status=-5
        deallocate(transform,x_physical,dax_physical,dx_physical, &
          packed_diagonal,physical_diagonal,preconditioner,seed_vector)
        return
      end if
      preconditioner(component)=1.0_dp/ &
        sign(max(abs(denominator),1.0e-8_dp),denominator)
    end do
    solve_tolerance=1.0e-10_dp
    if(present(tolerance)) solve_tolerance=tolerance
    solve_iterations=max(200,4*physical)
    if(present(max_iterations)) solve_iterations=max_iterations
    call solve_mrsf_tda_response_batch_matrix_free( &
      apply_physical_sigma_batch,omega,x_physical,dax_physical, &
      dx_physical,domega,residual_max,status,tol=solve_tolerance, &
      maxit=solve_iterations, &
      apply_preconditioner=apply_physical_preconditioner_batch)
    if(status/=0) status=-100+status
    if(status==0) dx_packed=matmul(transform,dx_physical)
    deallocate(transform,x_physical,dax_physical,dx_physical, &
      packed_diagonal,physical_diagonal,preconditioner,seed_vector)

  contains

    subroutine apply_physical_sigma(vector,result,operator_status)
      real(kind=dp), intent(in) :: vector(:)
      real(kind=dp), intent(out) :: result(:)
      integer, intent(out) :: operator_status
      real(kind=dp), allocatable :: packed_vector(:,:),packed_sigma(:,:)

      if(size(vector)/=physical .or. size(result)/=physical) then
        result=0.0_dp
        operator_status=-1
        return
      end if
      allocate(packed_vector(packed,1),packed_sigma(packed,1))
      packed_vector(:,1)=matmul(transform,vector)
      call apply_mrsf_tda_sigma(infos,int2_driver,mo_a,mo_b,fa,fb, &
        packed_vector,packed_sigma,operator_status)
      if(operator_status==0) then
        result=matmul(transpose(transform),packed_sigma(:,1))
      else
        result=0.0_dp
      end if
      deallocate(packed_vector,packed_sigma)
    end subroutine apply_physical_sigma

    subroutine apply_physical_sigma_batch(vectors,results,operator_status)
      real(kind=dp), intent(in) :: vectors(:,:)
      real(kind=dp), intent(out) :: results(:,:)
      integer, intent(out) :: operator_status
      real(kind=dp), allocatable :: packed_vectors(:,:),packed_sigma(:,:)
      integer :: nvec

      nvec=size(vectors,2)
      if(size(vectors,1)/=physical .or. &
         any(shape(results)/=[physical,nvec]) .or. nvec<=0) then
        results=0.0_dp
        operator_status=-1
        return
      end if
      allocate(packed_vectors(packed,nvec),packed_sigma(packed,nvec))
      packed_vectors=matmul(transform,vectors)
      call apply_mrsf_tda_sigma(infos,int2_driver,mo_a,mo_b,fa,fb, &
        packed_vectors,packed_sigma,operator_status)
      if(operator_status==0) then
        results=matmul(transpose(transform),packed_sigma)
      else
        results=0.0_dp
      end if
      deallocate(packed_vectors,packed_sigma)
    end subroutine apply_physical_sigma_batch

    subroutine apply_physical_preconditioner_batch(vectors,results, &
                                                   operator_status)
      real(kind=dp), intent(in) :: vectors(:,:)
      real(kind=dp), intent(out) :: results(:,:)
      integer, intent(out) :: operator_status
      integer :: nvec

      nvec=size(vectors,2)
      operator_status=0
      if(size(vectors,1)/=physical .or. &
         any(shape(results)/=[physical,nvec]) .or. nvec<=0) then
        results=0.0_dp
        operator_status=-1
        return
      end if
      results=vectors*spread(preconditioner,2,nvec)
    end subroutine apply_physical_preconditioner_batch

  end subroutine solve_mrsf_tda_amplitude_derivatives

end module tdhf_mrsf_hessian_amplitude_mod
