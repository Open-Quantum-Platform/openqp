module tdhf_mrsf_hessian_row_primitives_mod

  use precision, only: dp

  implicit none

  private
  public :: MRSF_ROWS_SUCCESS,MRSF_ROWS_BAD_SHAPE,MRSF_ROWS_BAD_MULTIPLICITY
  public :: MRSF_ROWS_XC_INCOMPLETE,MRSF_ROWS_INCONSISTENT_BALL
  public :: MRSF_ROWS_UNSUPPORTED_XC
  public :: assemble_mrsf_one_e_primitive_rows
  public :: combine_mrsf_response_row_blocks
  public :: check_mrsf_seven_density_alias

  integer,parameter :: MRSF_ROWS_SUCCESS=0
  integer,parameter :: MRSF_ROWS_BAD_SHAPE=-1
  integer,parameter :: MRSF_ROWS_BAD_MULTIPLICITY=-2
  integer,parameter :: MRSF_ROWS_XC_INCOMPLETE=-3
  integer,parameter :: MRSF_ROWS_INCONSISTENT_BALL=-4
  integer,parameter :: MRSF_ROWS_UNSUPPORTED_XC=-5

contains

!###############################################################################

  subroutine assemble_mrsf_one_e_primitive_rows(overlap_derivative, &
      core_derivative,d_relaxed_spin,d_energy_weighted,rows,status)
    ! Algebraic oracle for the response part of the one-electron gradient.
    ! The reference-density response is deliberately absent: its contraction
    ! is part of the separately evaluated ground-state Hessian.

    real(kind=dp),intent(in) :: overlap_derivative(:,:,:)
    real(kind=dp),intent(in) :: core_derivative(:,:,:)
    real(kind=dp),intent(in) :: d_relaxed_spin(:,:,:,:)
    real(kind=dp),intent(in) :: d_energy_weighted(:,:,:)
    real(kind=dp),intent(out) :: rows(:,:)
    integer,intent(out) :: status

    real(kind=dp),allocatable :: density(:,:),operator_matrix(:,:)
    integer :: coordinate,response,nbf,ncoordinate,nresponse

    status=MRSF_ROWS_SUCCESS
    rows=0.0_dp
    nbf=size(overlap_derivative,1)
    ncoordinate=size(overlap_derivative,3)
    nresponse=size(d_energy_weighted,3)
    if(nbf<=0 .or. ncoordinate<=0 .or. nresponse<=0 .or. &
       any(shape(overlap_derivative)/=[nbf,nbf,ncoordinate]) .or. &
       any(shape(core_derivative)/=[nbf,nbf,ncoordinate]) .or. &
       any(shape(d_relaxed_spin)/=[nbf,nbf,2,nresponse]) .or. &
       any(shape(d_energy_weighted)/=[nbf,nbf,nresponse]) .or. &
       any(shape(rows)/=[ncoordinate,nresponse])) then
      status=MRSF_ROWS_BAD_SHAPE
      return
    end if
    allocate(density(nbf,nbf),operator_matrix(nbf,nbf))
    do response=1,nresponse
      density=sum(d_relaxed_spin(:,:,:,response),dim=3)
      density=0.5_dp*(density+transpose(density))
      do coordinate=1,ncoordinate
        operator_matrix=0.5_dp*(core_derivative(:,:,coordinate)+ &
          transpose(core_derivative(:,:,coordinate)))
        rows(coordinate,response)=sum(density*operator_matrix)
        density=0.5_dp*(d_energy_weighted(:,:,response)+ &
          transpose(d_energy_weighted(:,:,response)))
        operator_matrix=0.5_dp*(overlap_derivative(:,:,coordinate)+ &
          transpose(overlap_derivative(:,:,coordinate)))
        ! MRSF stores W with the sign already used by the production
        ! eijden+2W overlap contraction; keep this oracle on that convention.
        rows(coordinate,response)=rows(coordinate,response)+ &
          2.0_dp*sum(density*operator_matrix)
        density=sum(d_relaxed_spin(:,:,:,response),dim=3)
        density=0.5_dp*(density+transpose(density))
      end do
    end do
    deallocate(density,operator_matrix)
  end subroutine assemble_mrsf_one_e_primitive_rows

!###############################################################################

  subroutine check_mrsf_seven_density_alias(d_seven_density,d_td_abxc, &
      tolerance,status)
    ! td_abxc and physical channel 7 are two names for the same spin-adapted
    ! transition-density object.  Treating them as separate terms would count
    ! the response-exchange contribution twice.

    real(kind=dp),intent(in) :: d_seven_density(:,:,:,:)
    real(kind=dp),intent(in) :: d_td_abxc(:,:,:)
    real(kind=dp),intent(in) :: tolerance
    integer,intent(out) :: status

    integer :: nbf,nresponse

    status=MRSF_ROWS_SUCCESS
    nbf=size(d_td_abxc,1)
    nresponse=size(d_td_abxc,3)
    if(nbf<=0 .or. nresponse<=0 .or. &
       any(shape(d_td_abxc)/=[nbf,nbf,nresponse]) .or. &
       any(shape(d_seven_density)/=[7,nbf,nbf,nresponse]) .or. &
       tolerance<0.0_dp) then
      status=MRSF_ROWS_BAD_SHAPE
      return
    end if
    if(maxval(abs(d_seven_density(7,:,:,:)-d_td_abxc))>tolerance) &
      status=MRSF_ROWS_INCONSISTENT_BALL
  end subroutine check_mrsf_seven_density_alias

!###############################################################################

  subroutine combine_mrsf_response_row_blocks(one_e,two_e_response, &
      two_e_reference_mixed,xc,is_dft,xc_complete,rows,rows_one,rows_two, &
      rows_xc,status)
    ! Keep the Cartesian rows unsymmetrized.  Hessian symmetrization is a
    ! downstream diagnostic and must not conceal a missing response term here.

    real(kind=dp),intent(in) :: one_e(:,:),two_e_response(:,:)
    real(kind=dp),intent(in) :: two_e_reference_mixed(:,:),xc(:,:)
    logical,intent(in) :: is_dft,xc_complete
    real(kind=dp),intent(out) :: rows(:,:),rows_one(:,:),rows_two(:,:)
    real(kind=dp),intent(out) :: rows_xc(:,:)
    integer,intent(out) :: status

    integer :: nrow,ncolumn

    status=MRSF_ROWS_SUCCESS
    rows=0.0_dp
    rows_one=0.0_dp
    rows_two=0.0_dp
    rows_xc=0.0_dp
    nrow=size(one_e,1)
    ncolumn=size(one_e,2)
    if(nrow<=0 .or. ncolumn<=0 .or. &
       any(shape(two_e_response)/=[nrow,ncolumn]) .or. &
       any(shape(two_e_reference_mixed)/=[nrow,ncolumn]) .or. &
       any(shape(xc)/=[nrow,ncolumn]) .or. &
       any(shape(rows)/=[nrow,ncolumn]) .or. &
       any(shape(rows_one)/=[nrow,ncolumn]) .or. &
       any(shape(rows_two)/=[nrow,ncolumn]) .or. &
       any(shape(rows_xc)/=[nrow,ncolumn])) then
      status=MRSF_ROWS_BAD_SHAPE
      return
    end if
    if(is_dft .and. .not.xc_complete) then
      status=MRSF_ROWS_XC_INCOMPLETE
      return
    end if
    rows_one=one_e
    rows_two=two_e_response+two_e_reference_mixed
    if(is_dft) rows_xc=xc
    rows=rows_one+rows_two+rows_xc
  end subroutine combine_mrsf_response_row_blocks

end module tdhf_mrsf_hessian_row_primitives_mod

#ifndef OQP_MRSF_ROWS_PRIMITIVES_ONLY
module tdhf_mrsf_hessian_rows_mod

  use precision, only: dp
  use tdhf_mrsf_hessian_row_primitives_mod, only: &
    MRSF_ROWS_SUCCESS,MRSF_ROWS_BAD_SHAPE,MRSF_ROWS_BAD_MULTIPLICITY, &
    MRSF_ROWS_XC_INCOMPLETE,MRSF_ROWS_INCONSISTENT_BALL, &
    MRSF_ROWS_UNSUPPORTED_XC,combine_mrsf_response_row_blocks, &
    check_mrsf_seven_density_alias

  implicit none

  private
  public :: build_tdhf_mrsf_response_rows
  public :: MRSF_ROWS_SUCCESS,MRSF_ROWS_BAD_SHAPE,MRSF_ROWS_BAD_MULTIPLICITY
  public :: MRSF_ROWS_XC_INCOMPLETE,MRSF_ROWS_INCONSISTENT_BALL
  public :: MRSF_ROWS_UNSUPPORTED_XC

contains

!###############################################################################

  subroutine build_tdhf_mrsf_response_rows(infos,reference_spin_density, &
      reference_spin_fock,relaxed_spin_density,seven_density, &
      d_reference_spin_density,d_reference_spin_fock, &
      d_relaxed_spin_density,d_energy_weighted_density,d_seven_density, &
      d_td_abxc,xc_rows,xc_complete,rows,rows_one,rows_two,rows_xc,status)
    ! Cartesian response-row assembly for the founding spin-adapted two-SOMO
    ! MRSF-TDDFT method initiated by Hiroya Nakata.  Every electronic-state
    ! quantity is represented by the two spin densities and the seven physical
    ! MRSF density channels.
    !
    ! This routine differentiates only first-derivative gradient contractions.
    ! Fixed-density nuclear and integral second derivatives remain in the
    ! high-precision ground-state Hessian.  Therefore raw, unsymmetrized rows
    ! are returned for a later completeness and symmetry check.

    use types, only: information
    use basis_tools, only: basis_set
    use mathlib, only: pack_matrix
    use grd1, only: grad_ee_overlap,grad_ee_kinetic, &
      grad_en_hellman_feynman,grad_en_pulay,grad_1e_ecp

    type(information),target,intent(inout) :: infos
    real(kind=dp),intent(in) :: reference_spin_density(:,:,:)
    real(kind=dp),intent(in) :: reference_spin_fock(:,:,:)
    real(kind=dp),intent(in) :: relaxed_spin_density(:,:,:)
    real(kind=dp),intent(in) :: seven_density(:,:,:)
    real(kind=dp),intent(in) :: d_reference_spin_density(:,:,:,:)
    real(kind=dp),intent(in) :: d_reference_spin_fock(:,:,:,:)
    real(kind=dp),intent(in) :: d_relaxed_spin_density(:,:,:,:)
    real(kind=dp),intent(in) :: d_energy_weighted_density(:,:,:)
    real(kind=dp),intent(in) :: d_seven_density(:,:,:,:)
    real(kind=dp),intent(in) :: d_td_abxc(:,:,:)
    real(kind=dp),intent(in) :: xc_rows(:,:)
    logical,intent(in) :: xc_complete
    real(kind=dp),intent(out) :: rows(:,:),rows_one(:,:),rows_two(:,:)
    real(kind=dp),intent(out) :: rows_xc(:,:)
    integer,intent(out) :: status

    type(basis_set),pointer :: basis
    real(kind=dp),allocatable :: packed_density(:),packed_response(:),density(:,:)
    real(kind=dp),allocatable :: gradient_one(:,:),gradient_plus(:,:)
    real(kind=dp),allocatable :: gradient_minus(:,:)
    real(kind=dp),allocatable :: one_bare(:,:),two_response(:,:)
    real(kind=dp),allocatable :: two_reference_mixed(:,:)
    real(kind=dp),allocatable,target :: dplus(:,:,:),dminus(:,:,:)
    real(kind=dp),allocatable,target :: pplus(:,:,:),pminus(:,:,:)
    real(kind=dp),allocatable,target :: splus(:,:,:),sminus(:,:,:)
    real(kind=dp) :: alias_tolerance,alias_scale
    integer :: response,nbf,nbf_tri,natom,ncart,local_status,spin,diagonal,ij
    logical :: is_dft

    status=MRSF_ROWS_SUCCESS
    rows=0.0_dp
    rows_one=0.0_dp
    rows_two=0.0_dp
    rows_xc=0.0_dp
    basis=>infos%basis
    basis%atoms=>infos%atoms
    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    ncart=3*natom
    nbf_tri=nbf*(nbf+1)/2
    is_dft=infos%control%hamilton==20
    if(nbf<=0 .or. ncart<=0 .or. &
       any(shape(reference_spin_density)/=[nbf,nbf,2]) .or. &
       any(shape(reference_spin_fock)/=[nbf,nbf,2]) .or. &
       any(shape(relaxed_spin_density)/=[nbf,nbf,2]) .or. &
       any(shape(seven_density)/=[7,nbf,nbf]) .or. &
       any(shape(d_reference_spin_density)/=[nbf,nbf,2,ncart]) .or. &
       any(shape(d_reference_spin_fock)/=[nbf,nbf,2,ncart]) .or. &
       any(shape(d_relaxed_spin_density)/=[nbf,nbf,2,ncart]) .or. &
       any(shape(d_energy_weighted_density)/=[nbf,nbf,ncart]) .or. &
       any(shape(d_seven_density)/=[7,nbf,nbf,ncart]) .or. &
       any(shape(d_td_abxc)/=[nbf,nbf,ncart]) .or. &
       any(shape(xc_rows)/=[ncart,ncart]) .or. &
       any(shape(rows)/=[ncart,ncart]) .or. &
       any(shape(rows_one)/=[ncart,ncart]) .or. &
       any(shape(rows_two)/=[ncart,ncart]) .or. &
       any(shape(rows_xc)/=[ncart,ncart])) then
      status=MRSF_ROWS_BAD_SHAPE
      return
    end if
    if(infos%tddft%mult/=1 .and. infos%tddft%mult/=3) then
      status=MRSF_ROWS_BAD_MULTIPLICITY
      return
    end if
    if(is_dft .and. (infos%dft%cam_flag .or. &
       infos%functional%needTau .or. infos%functional%needLapl)) then
      status=MRSF_ROWS_UNSUPPORTED_XC
      return
    end if
    if(is_dft .and. .not.xc_complete) then
      status=MRSF_ROWS_XC_INCOMPLETE
      return
    end if
    alias_scale=max(1.0_dp,maxval(abs(d_td_abxc)), &
      maxval(abs(d_seven_density(7,:,:,:))))
    alias_tolerance=128.0_dp*epsilon(1.0_dp)*alias_scale
    call check_mrsf_seven_density_alias(d_seven_density,d_td_abxc, &
      alias_tolerance,local_status)
    if(local_status/=MRSF_ROWS_SUCCESS) then
      status=local_status
      return
    end if

    allocate(packed_density(nbf_tri),packed_response(nbf_tri),density(nbf,nbf), &
      gradient_one(3,natom),gradient_plus(3,natom), &
      gradient_minus(3,natom),one_bare(ncart,ncart), &
      two_response(ncart,ncart), &
      two_reference_mixed(ncart,ncart),dplus(nbf,nbf,2), &
      dminus(nbf,nbf,2),pplus(nbf,nbf,2),pminus(nbf,nbf,2), &
      splus(7,nbf,nbf),sminus(7,nbf,nbf),source=0.0_dp)

    do response=1,ncart
      gradient_one=0.0_dp
      density=0.0_dp
      do spin=1,2
        density=density- &
          matmul(d_reference_spin_density(:,:,spin,response), &
            matmul(reference_spin_fock(:,:,spin), &
              reference_spin_density(:,:,spin))) &
          -matmul(reference_spin_density(:,:,spin), &
            matmul(d_reference_spin_fock(:,:,spin,response), &
              reference_spin_density(:,:,spin))) &
          -matmul(reference_spin_density(:,:,spin), &
            matmul(reference_spin_fock(:,:,spin), &
              d_reference_spin_density(:,:,spin,response)))
      end do
      density=0.5_dp*(density+transpose(density))
      call pack_matrix(density,packed_density)
      ij=0
      do diagonal=1,nbf
        ij=ij+diagonal
        packed_density(ij)=0.5_dp*packed_density(ij)
      end do
      density=0.5_dp*(d_energy_weighted_density(:,:,response)+ &
        transpose(d_energy_weighted_density(:,:,response)))
      call pack_matrix(2.0_dp*density,packed_response)
      packed_density=packed_density+packed_response
      call grad_ee_overlap(basis,packed_density,gradient_one)

      density=sum(d_reference_spin_density(:,:,:,response),dim=3)+ &
        sum(d_relaxed_spin_density(:,:,:,response),dim=3)
      density=0.5_dp*(density+transpose(density))
      call pack_matrix(density,packed_density)
      call grad_en_hellman_feynman(basis,basis%atoms%xyz, &
        basis%atoms%zn-basis%ecp_zn_num,packed_density,gradient_one)
      call grad_ee_kinetic(basis,packed_density,gradient_one)
      call grad_en_pulay(basis,basis%atoms%xyz, &
        basis%atoms%zn-basis%ecp_zn_num,packed_density,gradient_one)
      call grad_1e_ecp(infos,basis,basis%atoms%xyz,packed_density, &
        gradient_one)
      one_bare(:,response)=reshape(gradient_one,[ncart])

      dplus=reference_spin_density+d_reference_spin_density(:,:,:,response)
      dminus=reference_spin_density-d_reference_spin_density(:,:,:,response)
      pplus=relaxed_spin_density+d_relaxed_spin_density(:,:,:,response)
      pminus=relaxed_spin_density-d_relaxed_spin_density(:,:,:,response)
      splus=seven_density+d_seven_density(:,:,:,response)
      sminus=seven_density-d_seven_density(:,:,:,response)
      call mrsf_two_e_gradient_row(infos,basis,dplus,pplus,splus, &
        gradient_plus)
      call mrsf_two_e_gradient_row(infos,basis,dminus,pminus,sminus, &
        gradient_minus)
      two_response(:,response)=reshape( &
        0.5_dp*(gradient_plus-gradient_minus),[ncart])

      two_reference_mixed(:,response)=0.0_dp
    end do

    call combine_mrsf_response_row_blocks(one_bare,two_response, &
      two_reference_mixed,xc_rows,is_dft,xc_complete,rows,rows_one, &
      rows_two,rows_xc,status)
    deallocate(packed_density,packed_response,density,gradient_one,gradient_plus, &
      gradient_minus,one_bare,two_response,two_reference_mixed,dplus, &
      dminus,pplus,pminus,splus,sminus)
  end subroutine build_tdhf_mrsf_response_rows

!###############################################################################

  subroutine mrsf_two_e_gradient_row(infos,basis,density,relaxed,seven, &
      gradient)
    use types, only: information
    use basis_tools, only: basis_set
    use grd2, only: grd2_driver
    use tdhf_mrsf_gradient_mod, only: grd2_mrsf_compute_data_t

    type(information),target,intent(inout) :: infos
    type(basis_set),intent(in) :: basis
    real(kind=dp),target,intent(inout) :: density(:,:,:),relaxed(:,:,:)
    real(kind=dp),target,intent(inout) :: seven(:,:,:)
    real(kind=dp),intent(out) :: gradient(:,:)

    type(grd2_mrsf_compute_data_t) :: contraction
    real(kind=dp) :: reference_exchange,response_exchange

    reference_exchange=1.0_dp
    response_exchange=1.0_dp
    if(infos%control%hamilton==20) then
      reference_exchange=infos%dft%HFscale
      response_exchange=infos%tddft%HFscale
    end if
    gradient=0.0_dp
    contraction=grd2_mrsf_compute_data_t(d2=density,p2=relaxed, &
      spc2=seven,nbf=basis%nbf,hfscale=reference_exchange, &
      hfscale2=response_exchange, &
      spcscale=[infos%tddft%spc_coco,infos%tddft%spc_ovov, &
      infos%tddft%spc_coov],mrst=infos%tddft%mult)
    call contraction%init()
    call contraction%build_cart(basis)
    call grd2_driver(infos,basis,gradient,contraction,petite=.false.)
    call contraction%clean()
  end subroutine mrsf_two_e_gradient_row

end module tdhf_mrsf_hessian_rows_mod
#endif
