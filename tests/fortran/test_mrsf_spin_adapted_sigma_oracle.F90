program test_mrsf_spin_adapted_sigma_oracle

  use precision, only: dp
  use types, only: information
  use production_mrsf_extract, only: mrsfcbc, mrsfmntoia
  use tdhf_mrsf_conventions_mod, only: mrsf_raw_spc_multiplier

  implicit none

  integer, parameter :: nbf=4,npacked=9
  real(kind=dp) :: coefficients(nbf,nbf),trial(nbf,nbf)
  real(kind=dp) :: density(7,nbf,nbf),images(7,nbf,nbf)
  real(kind=dp) :: scaled(7,nbf,nbf),sigma(npacked,1)
  type(information) :: infos
  integer :: multiplicity,nphysical,column,channel,p,q,case_id,row

  coefficients=0.0_dp
  do p=1,nbf
    coefficients(p,p)=1.0_dp
  end do

  infos%basis%nbf=nbf
  infos%mol_prop%nelec_a=3
  infos%mol_prop%nelec_b=1
  infos%tddft%debug_mode=.false.

  do multiplicity=1,3,2
    infos%tddft%mult=multiplicity
    if (multiplicity==1) then
      nphysical=8
    else
      nphysical=6
    end if

    do column=1,nphysical
      call physical_unit_vector(multiplicity,column,trial)
      density=0.0_dp
      call mrsfcbc(infos,coefficients,coefficients,trial,density)

      do channel=1,7
        do p=1,nbf
          do q=1,nbf
            write(*,'(a,6(1x,i0),1x,es26.17e3)') &
              'D',multiplicity,column,channel,p,q,0,density(channel,p,q)
          end do
        end do
      end do

      call diagrammatic_images(density,images)
      do case_id=0,3
        call select_channel_family( &
          multiplicity,case_id,images,scaled)
        sigma=0.0_dp
        call mrsfmntoia( &
          infos,scaled,sigma,coefficients,coefficients,1)
        do row=1,npacked
          write(*,'(a,5(1x,i0),1x,es26.17e3)') &
            'Y',multiplicity,case_id,row,column,0,sigma(row,1)
        end do
      end do
    end do
  end do

contains

  subroutine physical_unit_vector(multiplicity,column,vector)
    integer, intent(in) :: multiplicity,column
    real(kind=dp), intent(out) :: vector(nbf,nbf)

    vector=0.0_dp
    select case(column)
    case(1)
      vector(1,2)=1.0_dp ! C -> O1
    case(2)
      vector(1,3)=1.0_dp ! C -> O2
    case(3)
      vector(1,4)=1.0_dp ! C -> V
    case(4)
      vector(2,4)=1.0_dp ! O1 -> V
    case(5)
      vector(3,4)=1.0_dp ! O2 -> V
    case(6)
      if (multiplicity==1) then
        vector(3,2)=1.0_dp ! G
      else
        vector(2,2)=1.0_dp ! folded (L+R)/sqrt(2)
      end if
    case(7)
      if (multiplicity/=1) error stop 'triplet has no D coordinate'
      vector(2,3)=1.0_dp ! D
    case(8)
      if (multiplicity/=1) error stop 'triplet has six coordinates'
      vector(2,2)=1.0_dp ! folded (L-R)/sqrt(2)
    case default
      error stop 'physical column is out of range'
    end select
  end subroutine physical_unit_vector

  subroutine diagrammatic_images(densities,result)
    real(kind=dp), intent(in) :: densities(7,nbf,nbf)
    real(kind=dp), intent(out) :: result(7,nbf,nbf)
    integer :: c,mu,nu,k,l

    result=0.0_dp
    do c=1,7
      do mu=1,nbf
        do nu=1,nbf
          do k=1,nbf
            do l=1,nbf
              if (c<=4) result(c,mu,nu)=result(c,mu,nu) &
                +eri(mu,nu,k,l)*(densities(c,k,l)+densities(c,l,k))
              result(c,mu,nu)=result(c,mu,nu) &
                -eri(mu,k,nu,l)*densities(c,k,l)
            end do
          end do
        end do
      end do
    end do
  end subroutine diagrammatic_images

  subroutine select_channel_family( &
      multiplicity,case_id,images,selected)
    integer, intent(in) :: multiplicity,case_id
    real(kind=dp), intent(in) :: images(7,nbf,nbf)
    real(kind=dp), intent(out) :: selected(7,nbf,nbf)
    real(kind=dp) :: spin_sign

    selected=0.0_dp
    selected(7,:,:)=images(7,:,:)
    spin_sign=real(mrsf_raw_spc_multiplier(multiplicity),dp)
    select case(case_id)
    case(0) ! A0 only
    case(1) ! CO--CO
      selected(6,:,:)=spin_sign*images(6,:,:)
    case(2) ! OV--OV
      selected(5,:,:)=spin_sign*images(5,:,:)
    case(3) ! CO--OV and its transpose
      selected(1:4,:,:)=spin_sign*images(1:4,:,:)
    case default
      error stop 'unknown channel-family case'
    end select
  end subroutine select_channel_family

  pure real(kind=dp) function eri(p,q,r,s) result(value)
    integer, intent(in) :: p,q,r,s
    integer :: left,right,lo,hi

    left=pair_index(p,q)
    right=pair_index(r,s)
    lo=min(left,right)
    hi=max(left,right)
    value=real(17*lo+31*hi+7*lo*hi,dp)/997.0_dp
  end function eri

  pure integer function pair_index(p,q) result(index)
    integer, intent(in) :: p,q
    integer :: i,j,lo,hi

    lo=min(p,q)
    hi=max(p,q)
    index=0
    do i=1,nbf
      do j=i,nbf
        index=index+1
        if (i==lo .and. j==hi) return
      end do
    end do
    error stop 'orbital index is out of range'
  end function pair_index

end program test_mrsf_spin_adapted_sigma_oracle
