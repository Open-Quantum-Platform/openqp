module tdhf_mrsf_hessian_components_mod

  ! Methodological starting point: Hiroya Nakata's analytical TDHF/TDDFT
  ! Hessian implementation, extended here to the two-SOMO MRSF Lagrangian.

  use precision, only: dp

  implicit none

  private
  public :: mrsf_hessian_is_applicable
  public :: assemble_mrsf_cartesian_hessian
  public :: project_mrsf_rigid_translations

contains

!###############################################################################

  pure logical function mrsf_hessian_is_applicable(scf_type,reference_mult, &
      nocca,noccb,target_mult,umrsf,is_dft,needs_tau, &
      needs_laplacian,range_separated,mpi_size) result(applicable)
    integer, intent(in) :: scf_type,reference_mult,nocca,noccb,target_mult
    integer, intent(in) :: mpi_size
    logical, intent(in) :: umrsf,is_dft,needs_tau, &
      needs_laplacian,range_separated

    applicable=scf_type==3 .and. reference_mult==3 .and. &
      nocca-noccb==2 .and. (target_mult==1 .or. target_mult==3) .and. &
      .not.umrsf .and. mpi_size==1 .and. &
      (.not.is_dft .or. (.not.needs_tau .and. &
       .not.needs_laplacian .and. .not.range_separated))
  end function mrsf_hessian_is_applicable

!###############################################################################

  pure subroutine assemble_mrsf_cartesian_hessian(fixed_density,xc_fixed, &
      response_rows,hessian,row_asymmetry,status)
    ! Stationary-Lagrangian assembly.  The raw directional rows remain
    ! unsymmetrized until this final boundary so their antisymmetric part tests
    ! the completeness of amplitude, orbital, and Z-vector response.

    real(kind=dp), intent(in) :: fixed_density(:,:),xc_fixed(:,:), &
      response_rows(:,:)
    real(kind=dp), intent(out) :: hessian(:,:),row_asymmetry
    integer, intent(out) :: status

    integer :: ncoord

    ncoord=size(hessian,1)
    status=0
    hessian=0.0_dp
    row_asymmetry=0.0_dp
    if(ncoord<=0 .or. size(hessian,2)/=ncoord .or. &
       any(shape(fixed_density)/=[ncoord,ncoord]) .or. &
       any(shape(xc_fixed)/=[ncoord,ncoord]) .or. &
       any(shape(response_rows)/=[ncoord,ncoord])) then
      status=-1
      return
    end if
    row_asymmetry=maxval(abs(response_rows-transpose(response_rows)))
    hessian=fixed_density+xc_fixed+ &
      0.5_dp*(response_rows+transpose(response_rows))
  end subroutine assemble_mrsf_cartesian_hessian

!###############################################################################

  pure subroutine project_mrsf_rigid_translations(hessian,natom,status)
    ! Apply the Cartesian translation projector on both indices.  This is a
    ! numerical cleanup after full analytic assembly, not a replacement for a
    ! missing response term: the pre-projection row asymmetry remains reported.

    real(kind=dp), intent(inout) :: hessian(:,:)
    integer, intent(in) :: natom
    integer, intent(out) :: status

    real(kind=dp) :: mean_value
    integer :: atom,cart,index,ncoord

    ncoord=3*natom
    status=0
    if(natom<=0 .or. any(shape(hessian)/=[ncoord,ncoord])) then
      status=-1
      return
    end if
    do cart=1,3
      do index=1,ncoord
        mean_value=0.0_dp
        do atom=1,natom
          mean_value=mean_value+hessian(3*(atom-1)+cart,index)
        end do
        mean_value=mean_value/real(natom,dp)
        do atom=1,natom
          hessian(3*(atom-1)+cart,index)= &
            hessian(3*(atom-1)+cart,index)-mean_value
        end do
      end do
      do index=1,ncoord
        mean_value=0.0_dp
        do atom=1,natom
          mean_value=mean_value+hessian(index,3*(atom-1)+cart)
        end do
        mean_value=mean_value/real(natom,dp)
        do atom=1,natom
          hessian(index,3*(atom-1)+cart)= &
            hessian(index,3*(atom-1)+cart)-mean_value
        end do
      end do
    end do
    hessian=0.5_dp*(hessian+transpose(hessian))
  end subroutine project_mrsf_rigid_translations

end module tdhf_mrsf_hessian_components_mod
