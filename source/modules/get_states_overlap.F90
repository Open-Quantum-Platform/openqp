!> @brief Module for calculating state overlaps
!>        and derivative coupling matrix elements
!>
!> @date Aug 2024
!>
!> @author Konstantin Komarov
!>
module get_state_overlap_mod

  implicit none

  character(len=*), parameter :: module_name = "get_state_overlap_mod"

  public get_states_overlap

contains

!> @brief C-interoperable wrapper for get_states_overlap
!>
!> @param[in] c_handle   C handle for the information structure
!>
  subroutine get_state_overlap_C(c_handle) bind(C, name="get_states_overlap")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call get_states_overlap(inf)
  end subroutine get_state_overlap_c

!> @brief Main subroutine for calculating state overlaps
!>        and derivative coupling matrix elements
!>
!> @param[in,out] infos   Information class containing molecule parameters
!>
  subroutine get_states_overlap(infos)

    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use atomic_structure_m, only: atomic_structure
    use messages, only: show_message, with_abort
    use tdhf_mrsf_lib, only: mrsfxvec
    use util, only: measure_time

    implicit none

    character(len=*), parameter :: subroutine_name = "get_states_overlap"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    integer :: nstates, mrst, xvec_dim, nbf, ok, j
    integer :: noca, nocb, ndtlf

    real(kind=dp), allocatable, target :: bvec(:,:), bvec_old(:,:)

    ! Tagarray
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
        OQP_td_bvec_mo_old, OQP_td_bvec_mo, OQP_overlap_mo /)
    character(len=*), parameter :: tags_alloc(*) = (/ character(len=80) :: &
        OQP_nac /)
    real(kind=dp), contiguous, pointer :: bvec_mo(:,:), bvec_mo_old(:,:), &
          nac_out(:,:), overlap_mo(:,:), td_states_phase(:), td_states_overlap(:,:)

    ! Files open
    open (unit=IW, file=infos%log_filename, position="append")

    ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms
    nstates = infos%tddft%nstate
    mrst = infos%tddft%mult
    nbf = basis%nbf
    xvec_dim = infos%mol_prop%nelec_a*(nbf-infos%mol_prop%nelec_b)

!   ndtlf = 0          less accurate
!   ndtlf = 1 : tlf(1)
!   ndtlf = 2 : tlf(2) most accurate
    ndtlf = infos%tddft%tlf

    ! Allocate data for outputing in python level
    call infos%dat%remove_records(tags_alloc)
    call infos%dat%reserve_data(OQP_td_states_phase, ta_type_real64, &
          nstates, (/ nstates /), comment=OQP_td_states_phase_comment)
    call infos%dat%reserve_data(OQP_td_states_overlap, ta_type_real64, &
          nstates*nstates, (/ nstates, nstates /), comment=OQP_td_states_overlap_comment)
    call infos%dat%reserve_data(OQP_nac, ta_type_real64, &
          nstates*nstates, (/ nstates, nstates /), comment=OQP_nac_comment)

    ! Load data from python level
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_overlap_mo, overlap_mo)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo_old, bvec_mo_old)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_td_states_phase, td_states_phase)
    call tagarray_get_data(infos%dat, OQP_td_states_overlap, td_states_overlap)
    call tagarray_get_data(infos%dat, OQP_nac, nac_out)

    allocate(bvec(xvec_dim,nstates), &
             bvec_old(xvec_dim,nstates), &
             source=0.0_dp, stat=ok)
    if( ok/=0 ) call show_message('Cannot allocate memory',with_abort)

    noca = infos%mol_prop%nelec_a
    nocb = infos%mol_prop%nelec_b

    do j = 1, nstates
      call mrsfxvec(infos, bvec_mo_old(:,j), bvec_old(:,j))
      call mrsfxvec(infos, bvec_mo(:,j), bvec(:,j))
    end do

    call check_states_phase(bvec, bvec_old, td_states_phase)

    call compute_states_overlap( &
          infos, overlap_mo, td_states_overlap, bvec, &
          bvec_old, nbf,noca, nocb, nstates, ndtlf)

    call get_dcv(nac_out, td_states_overlap, nstates)

!   call print_nac(infos, td_states_overlap, nac_out)

!   Print timings
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    close(iw)

  end subroutine get_states_overlap

  subroutine check_states_phase(Bvec, Bvec_old, td_states_phase)
    use precision, only: dp
    use io_constants, only: iw
    implicit none

    real(kind=dp), dimension(:,:) :: Bvec
    real(kind=dp), dimension(:,:) :: Bvec_old
    real(kind=dp), dimension(:) :: td_states_phase

    integer :: i

    ! Get overlap and correct the sign of X amplitude before state tracking
    write(iw, fmt='(/x,a,/x,a,/7x,"State  Overlap")') &
      'Check the sign of X amplitude', &
      'with respect to previous geometry'

    do i = 1, ubound(Bvec, 2)
      td_states_phase(i) = dot_product(Bvec_old(:,i), Bvec(:,i))
     if (td_states_phase(i) < 0.0d0) then
       Bvec(:,i) = -1.0d0*Bvec(:,i)
       td_states_phase(i) = -1.0d0*td_states_phase(i)
     endif
      write(iw, fmt='(6x,i4,x,f12.8)') i, td_states_phase(i)
    end do

  end subroutine

!> @brief Compute overlap integrals between MRSF response states
!>        at different MD time steps
!>
!> @details Fast overlap evaluations using the TLF approximation
!>          introduced in JCTC 15 882 (2019)
!>
!> @author Seunghoon Lee, Konstantin Komarov
!>
!> @param[in] infos     Information structure
!> @param[in] s_mo      Overlap matrix in MO basis
!> @param[out] s_st     State overlap matrix
!> @param[in] mo_a      MO coefficients at current time step
!> @param[in] mo_a_old  MO coefficients at previous time step
!> @param[in] nbf       Number of basis functions
!> @param[in] noca      Number of occupied alpha orbitals
!> @param[in] nocb      Number of occupied beta orbitals
!> @param[in] nstates   Number of states
!> @param[in] ndtlf     TLF approximation order
!>
  subroutine compute_states_overlap(&
              infos, s_mo, s_st, mo_a, mo_a_old, &
              nbf, noca, nocb, nstates, ndtlf)

    use precision, only: dp
    use types, only: information
!$  use omp_lib

    implicit none

    type(information) :: infos
    real(kind=dp), intent(inout), dimension(:,:) :: s_st
    real(kind=dp), dimension(noca,nbf-nocb,*) :: mo_a, mo_a_old
    real(kind=dp), dimension(:,:) :: s_mo
    integer :: nbf, nstates, noca, nocb, ndtlf

    integer :: ni, oi, pi, qi, ri, si, i, &
               ioc, ioc1, ioc2, &
               ivir, joc, jvir, nvirb
    real(kind=dp), parameter :: sqrt2 = 1/sqrt(2.0_dp)
    real(kind=dp), allocatable, dimension(:,:,:) :: &
          alpham, betam, deltam, gammam
    real(kind=dp), allocatable, dimension(:,:) :: &
          s_ij, s_ab, s_ia

    nvirb = nbf-nocb

    allocate(alpham(nstates,noca,nvirb), &
             betam(nstates,noca,nvirb), &
             deltam(nstates,noca,noca), &
             gammam(nstates,noca,noca), &
             s_ij(noca,noca), &
             s_ab(nvirb,nvirb), &
             s_ia(noca,nvirb), &
             source=0.0_dp)

!   get S_ij, S_ab, S_ia
    call mrsf_tlf(infos, s_mo, s_ij, s_ab, s_ia, ndtlf)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(oi, ioc, ivir, jvir, ni, ioc1, ioc2, pi, qi, ri, si)

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
    do oi = 1, nstates
      do ioc = 1, noca
        do ivir = 1, nvirb
          alpham(oi, ioc, ivir) = 0.0_dp
          do jvir = 1, nvirb
            if((ioc > nocb) .and. (jvir <= 2)) cycle
            alpham(oi, ioc, ivir) = alpham(oi, ioc, ivir) &
                                   + mo_a_old(ioc,jvir,oi) * s_ab(jvir,ivir)
          end do
        end do
      end do
    end do
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
    do ni = 1, nstates
      do ioc = 1, noca
        do ivir = 1, nvirb
          betam(ni, ioc, ivir) = 0.0_dp
          do joc = 1, noca
            if((joc > nocb) .and. (ivir <= 2)) cycle
            betam(ni,ioc,ivir) = betam(ni,ioc,ivir) &
                                + mo_a(joc,ivir,ni) * s_ij(ioc,joc)
          end do
        end do
      end do
    end do
!$OMP END DO
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
    do oi = 1, nstates
      do ioc1 = 1, noca
        do ioc2 = 1, noca
          gammam(oi, ioc1, ioc2) = 0.0_dp
          do jvir = 1, nvirb
            if((ioc1 > nocb) .and. (jvir <= 2)) cycle
            gammam(oi,ioc1,ioc2) = gammam(oi,ioc1,ioc2) &
                                  + mo_a_old(ioc1,jvir,oi)*s_ia(ioc2,jvir)
          end do
        end do
      end do
    end do
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
    do ni = 1, nstates
      do ioc1 = 1, noca
        do ioc2 = 1, noca
          deltam(ni, ioc1, ioc2) = 0.0_dp
          do jvir = 1, nvirb
            if((ioc2 > nocb) .and. (jvir <= 2)) cycle
            deltam(ni,ioc1,ioc2) = deltam(ni,ioc1,ioc2) &
                                  + mo_a(ioc2,jvir,ni) * s_ia(ioc1,jvir)
          end do
        end do
      end do
    end do
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3) REDUCTION(+:s_st)
   ! 4 index summation
    do oi = 1, nstates
      do ni = 1, nstates
        do pi = nocb+1,noca
          do qi = nocb+1,noca
            do ri = nocb+1,noca
              do si = nocb+1, noca
                s_st(oi,ni) = s_st(oi,ni) &
                             + mo_a_old(pi,qi-nocb,oi) * s_ij(pi,ri) &
                              * mo_a(ri,si-nocb,ni) * s_ab(qi-nocb,si-nocb)
              end do
            end do
          end do
        end do
      end do
    end do
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3) REDUCTION(+:s_st)
    do oi = 1, nstates
      do ni = 1, nstates
        do pi = nocb+1, noca
          do qi = nocb+1, noca
            do ri = 1, noca
              do si = nocb+1, nbf
                if((ri >= nocb+1) .and. (si <= noca)) cycle
                s_st(oi,ni) = s_st(oi,ni) &
                             + mo_a_old(pi,qi-nocb,oi) * mo_a(ri,si-nocb,ni) &
                              * ( s_ij(pi,ri) * s_ab(qi-nocb,si-nocb) &
                                + s_ia(pi,si-nocb) * s_ia(ri,qi-nocb) ) * sqrt2
              end do
            end do
          end do
        end do
      end do
    end do
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3) REDUCTION(+:s_st)
    do oi = 1, nstates
      do ni = 1, nstates
        do pi = 1, noca
          do qi = nocb+1, nbf
            if((pi >= nocb+1) .and. (qi <= noca)) cycle
            do ri = nocb+1, noca
              do si = nocb+1, noca
                s_st(oi,ni) = s_st(oi,ni) &
                             + mo_a_old(pi,qi-nocb,oi) * mo_a(ri,si-nocb,ni) &
                              * ( s_ij(pi,ri) * s_ab(qi-nocb,si-nocb) &
                                + s_ia(pi,si-nocb) * s_ia(ri,qi-nocb) ) * sqrt2
              end do
            end do
          end do
        end do
      end do
    end do
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3) REDUCTION(+:s_st)
    do oi = 1, nstates
      do ni = 1, nstates
        do ioc = 1, noca
          do ivir = 1, nvirb
            s_st(oi,ni) = s_st(oi,ni) &
                         + alpham(oi,ioc,ivir) * betam(ni,ioc,ivir)
          end do
        end do
      end do
    end do
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3) REDUCTION(+:s_st)
    do oi = 1, nstates
      do ni = 1, nstates
        do ioc1 = 1, noca
          do ioc2 = 1, noca
            s_st(oi,ni) = s_st(oi,ni) &
                         + gammam(oi,ioc1,ioc2) * deltam(ni,ioc1,ioc2)
          end do
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    do i = 1, nstates
      s_st(:,i) = s_st(:,i) / norm2(s_st(:,i))
    end do

  end subroutine compute_states_overlap

!>
!>     @brief Compute overlap integrals between CSFs
!>            of MRSF-TDDFT using TLF approximation
!>
!>     @details TLF approximation for SF- and LR-TDDFT
!>              is introduced in JCTC 15 882 (2019)
!>
!>     @author Seunghoon Lee, Kostantin Komarov
!>
  subroutine mrsf_tlf(infos, s_mo, s_ij, s_ab, s_ia, ndtlf)

    use precision, only: dp
    use types, only: information

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), dimension(:,:) :: s_mo
    real(kind=dp), intent(out), dimension(:,:) :: s_ij
    real(kind=dp), intent(out), dimension(:,:) :: s_ab
    real(kind=dp), intent(out), dimension(:,:) :: s_ia
    integer, intent(in) :: ndtlf

    integer :: nbf, noca, nocb, nvirb, noc
    integer :: i, i1, i2, ia1, ia2, j1, j2
    real(kind=dp) :: precomp, temp1, temp2, tmp, tmp1, tmp2

    nbf = infos%basis%nbf
    noca = infos%mol_prop%nelec_a
    nocb = infos%mol_prop%nelec_b
    nvirb = nbf - nocb

    noc = noca-1
    i1 = 0
    i2 = 0
    ia1 = 0
    ia2 = 0

    select case (ndtlf)
    case(0)
!     alpha determinant
      do i1 = 1, noca
         do i2 = 1, noca
            call ov_exact(temp1, i1, i2, ia1, ia2, s_mo, 1, noc, 1)
            s_ij(i1,i2) = temp1
         end do
      end do

!     1-2 det
      do j1 = 1, nvirb
         ia1 = nocb+j1
         do j2 = 1, nvirb
            ia2 = nocb+j2
            call ov_exact(temp2, i1, i2, ia1, ia2, s_mo, 1, noc, 2)
            s_ab(j1,j2) = temp2
         end do
      end do

    case(1)
      precomp = 1.0_dp
      do i = 1, noca
         precomp = precomp*s_mo(i,i)
      end do
!     alpha determinant
      do i1 = 1, noca
         do i2 = 1, noca
         call tlf_exp(tmp, 11, i1, i2, s_mo, precomp, noca, nbf)
         if (i1/=i2) then
            tmp = -1.0_dp*tmp
         end if
         s_ij(i1,i2) = tmp
         end do
      end do

      precomp = 1.0_dp
      do i = 1, noca-2
         precomp = precomp*s_mo(i,i)
      end do
!     1-2 det
      do j1 = 1, nvirb
         ia1 = nocb+j1
         do j2 = 1, nvirb
            ia2 = nocb+j2
            call tlf_exp(tmp, 21, ia1, ia2, s_mo, precomp, noca, nbf)
            s_ab(j1,j2) = tmp
         end do
      end do

    case (2)
      precomp = 1.0_dp
      do i = 1, noca
         precomp = precomp*s_mo(i,i)
      enddo
!     alpha determinant
      do i1 = 1, noca
         do i2 = 1, noca
            call tlf_exp(tmp1, 11, i1, i2, s_mo, precomp, noca, nbf)
            call tlf_exp(tmp2, 12, i1, i2, s_mo, precomp, noca, nbf)
            tmp = tmp1+tmp2
            if (i1/=i2) then
               tmp = -1.0_dp*tmp
            end if
            s_ij(i1,i2) = tmp
         end do
      end do

      precomp = 1.d+00
      do i = 1, noca-2
         precomp = precomp*s_mo(i,i)
      end do
!     1-2 det
      do j1 = 1, nvirb
         ia1 = nocb+j1
         do j2 = 1, nvirb
            ia2 = nocb+j2
            call tlf_exp(tmp1, 21, ia1, ia2, s_mo, precomp, noca, nbf)
            call tlf_exp(tmp2, 22, ia1, ia2, s_mo, precomp, noca, nbf)
            s_ab(j1,j2) = tmp1+tmp2
         end do
      end do
    case default
       error stop "Unknown TLF value (0,1,2)"
    end select

!   alpha determinant
    do i1 = 1, noca
       do j1 = 1, nvirb
          ia1 = nocb+j1
          call ov_exact(temp1 ,i1, i2, ia1, ia2, s_mo, 1, noc, 3)
          s_ia(i1,j1) = temp1
       end do
    end do

  end subroutine mrsf_tlf

  subroutine ov_exact(temp1, i1, i2, ia1, ia2, s_mo, ilow, noc, itype)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(out) :: temp1
    integer, intent(in) :: i1, i2, ia1, ia2
    real(kind=dp), intent(in), dimension(:,:) :: s_mo
    integer, intent(in) :: ilow, noc, itype

    real(kind=dp), dimension(noc*noc) :: ddet
    integer :: i, iipp, imax, imin, ipp

    select case (itype)
    case (1)
       imin = min(i1,i2)
       imax = max(i1,i2)
    !  (1,1) block
       do i = 1, imin-1
          do ipp = 1, imin-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+ilow-1,ipp+ilow-1)
          end do
       end do
    !  (1,2) block
       do i = 1, imin-1
          do ipp = imin, imax-2
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+ilow-1,ipp+1+ilow-1)
          end do
       end do
    !  (1,3) block
       do i = 1, imin-1
          do ipp = imax-1, noc-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+ilow-1,ipp+2+ilow-1)
          end do
       end do
    !  (2,1) block
       do i = imin, imax-2
          do ipp = 1, imin-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+1+ilow-1,ipp+ilow-1)
          end do
       end do
    !  (2,2) block
       do i = imin, imax-2
          do ipp = imin, imax-2
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+1+ilow-1,ipp+1+ilow-1)
          end do
       end do
    !  (2,3) block
       do i = imin, imax-2
          do ipp = imax-1, noc-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+1+ilow-1,ipp+2+ilow-1)
          end do
       end do
    !  (3,1) block
       do i = imax-1, noc-1
          do ipp = 1, imin-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+2+ilow-1,ipp+ilow-1)
          end do
       end do
    !  (3,2) block
       do i = imax-1, noc-1
          do ipp = imin, imax-2
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+2+ilow-1,ipp+1+ilow-1)
          end do
       end do
    !  (3,3) block
       do i = imax-1, noc-1
          do ipp = imax-1, noc-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+2+ilow-1,ipp+2+ilow-1)
          end do
       end do
    !  (1,4) block
       do i = 1, imin-1
          do ipp = noc, noc
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+ilow-1,i1+ilow-1)
          end do
       end do
    !  (2,4) block
       do i = imin, imax-2
          do ipp = noc, noc
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+1+ilow-1,i1+ilow-1)
          end do
       end do
    !  (3,4) block
       do i = imax-1, noc-1
          do ipp = noc, noc
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+2+ilow-1,i1+ilow-1)
          end do
       end do
    !  (4,1) block
       do i = noc, noc
          do ipp = 1, imin-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i2+ilow-1,ipp+ilow-1)
          end do
       end do
    !  (4,2) block
       do i = noc, noc
          do ipp = imin, imax-2
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i2+ilow-1,ipp+1+ilow-1)
          end do
       end do
    !  (4,3) block
       do i = noc, noc
          do ipp = imax-1, noc-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i2+ilow-1,ipp+2+ilow-1)
          end do
       end do
    !  (4,4) block
       do i = noc, noc
          do ipp = noc, noc
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i2+ilow-1,i1+ilow-1)
          end do
       end do
    !  Calculate alpha determinant
       temp1 = comp_det(ddet, noc)
       if (i1==i2) then
          return
       else if (i1/=i2) then
          temp1 = -1.0_dp*temp1
          return
       endif

    case (2)
    !  (1,1) block
       do i = 1, noc-1
          do ipp = 1, noc-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+ilow-1,ipp+ilow-1)
          end do
       end do
    !  (1,2) block
       ipp = noc
       do i = 1, noc-1
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i+ilow-1,ia2)
       end do
    !  (2,1) block
       i = noc
       do ipp = 1, noc-1
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(ia1,ipp+ilow-1)
       end do
    !  (2,2) block
       i = noc
       ipp = noc
       iipp = (ipp-1)*noc+i
       ddet(iipp) = s_mo(ia1,ia2)

    !  Calculate 2 det
       temp1 = comp_det(ddet, noc)
       return

    case (3)
    !  (1,1) block
       do i = 1, i1-1
          do ipp = 1, i1-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i,ipp)
          end do
       end do
    !  (1,2) block
       do i = 1, i1-1
          do ipp = i1, noc-2
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i,ipp+1)
          end do
       end do
    !  (1,3) block
       do i = 1, i1-1
          ipp = noc-1
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i,i1)
       end do
    !  (1,4) block
       do i = 1, i1-1
          ipp = noc
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i,ia1)
       end do
    !  (2,1) block
       do i = i1, noc-2
          do ipp = 1, i1-1
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+1,ipp)
          end do
       end do
    !  (2,2) block
       do i = i1, noc-2
          do ipp = i1, noc-2
             iipp = (ipp-1)*noc+i
             ddet(iipp) = s_mo(i+1,ipp+1)
          end do
       end do
    !  (2,3) block
       do i = i1, noc-2
          ipp = noc-1
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i+1,i1)
       end do
    !  (2,4) block
       do i = i1, noc-2
          ipp = noc
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i+1,ia1)
       end do
    !  (3,1) block
       i = noc-1
       do ipp = 1, i1-1
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i+1,ipp)
       end do
    !  (3,2) block
       i = noc-1
       do ipp = i1, noc-2
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i+1,ipp+1)
       end do
    !  (3,3) block
       i = noc-1
       ipp = noc-1
       iipp = (ipp-1)*noc+i
       ddet(iipp) = s_mo(i+1,i1)
    !  (3,4) block
       i = noc-1
       ipp = noc
       iipp = (ipp-1)*noc+i
       ddet(iipp) = s_mo(i+1,ia1)
    !  (4,1) block
       i = noc
       do ipp = 1, i1-1
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i+1,ipp)
       end do
    !  (4,2) block
       i = noc
       do ipp = i1, noc-2
          iipp = (ipp-1)*noc+i
          ddet(iipp) = s_mo(i+1,ipp+1)
       end do
    !  (4,3) block
       i = noc
       ipp = noc-1
       iipp = (ipp-1)*noc+i
       ddet(iipp) = s_mo(i+1,i1)
    !  (4,4) block
       i = noc
       ipp = noc
       iipp = (ipp-1)*noc+i
       ddet(iipp) = s_mo(i+1,ia1)

    !  Calculate alpha determinant
       temp1 = comp_det(ddet, noc)
    end select

  end subroutine ov_exact

  subroutine tlf_exp(ov,itype,i1,i2,s_mo,precomp,noca,nbf)

    use precision, only: dp
    implicit none

    real(kind=dp) :: ov, precomp
    integer :: i1, i2, itype, nbf, noca
    real(kind=dp), dimension(nbf,nbf) :: s_mo

    real(kind=dp) :: ov1, ov2
    integer :: ia1, ia2, l, lp

!   itype=11 : alpha 1st order
!   itype=21 : beta  1st order
!   itype=12 : alpha 2nd order
!   itype=22 : beta  2nd order

    select case (itype)
    case (11)
       ov = precomp*s_mo(i2,i1)/(s_mo(i1,i1)*s_mo(i2,i2))
       return

    case (21)
       ia1 = i1
       ia2 = i2
       ov = precomp*s_mo(ia1,ia2)
       return

    case (12)
       if (i1/=i2) then
          ov = 0.0_dp
          do l = 1, noca
            if (l/=i1 .and. l/=i2) then
               ov = ov+s_mo(i2,l)*s_mo(l,i1)/s_mo(l,l)
            end if
          end do
          ov = -1.0_dp*precomp*ov/(s_mo(i1,i1)*s_mo(i2,i2))
          return
       else
          ov = 0.0_dp
          do l = 1, noca-1
            if (l/=i1) then
              do lp = l+1, noca
                if (lp/=i1) then
                  ov = ov+s_mo(l,lp)*s_mo(lp,l)/(s_mo(l,l)*s_mo(lp,lp))
                end if
              end do
            end if
          end do
          ov = -1.0_dp*precomp*ov/s_mo(i1,i1)
          return
       end if

    case (22)
       ia1 = i1
       ia2 = i2
       if (ia1/=ia2) then
          ov = 0.0_dp
          do l = 1, noca-2
             ov = ov+s_mo(ia1,l)*s_mo(l,ia2)/s_mo(l,l)
          end do
          ov = -1.0_dp*precomp*ov
          return
       else
          ov1 = 0.d+00
          do l = 1, noca-2
             ov1 = ov1+s_mo(ia1,l)*s_mo(l,ia2)/s_mo(l,l)
          end do
          ov2 = 0.0_dp
          do l = 1, noca-3
            do lp = l+1, noca-2
               ov2 = ov2+s_mo(l,lp)*s_mo(lp,l)/(s_mo(l,l)*s_mo(lp,lp))
            end do
          end do
          ov2 = ov2*s_mo(ia1,ia1)
          ov = -1.0_dp*precomp*(ov1+ov2)
          return
       end if

    case default
       error stop "Unknown itype for tlf_exp"

    end select

  end subroutine tlf_exp
!>
!> @brief Compute derivative coupling vectors (DCV) using the finite difference method,
!>        typically denoted as d_IJ between states I and J.
!>
!>        d_IJ = <I| d/dR |J> = h_IJ / (E_I - E_J),
!>        where h_IJ is the nonadiabatic coupling, defined as
!>        h_IJ = <I| dH/dR |J>.
!>
!>        Note that R is an entire vector, and so is d_IJ.
!>
!>        This routine computes d_IJ^a = <I| d/da |J>, which is a single component of d_IJ by TLF.
!>        The component a can be either a geometric or a time derivative.
!>        The geometric derivative is used to construct the DCV, while
!>        the time derivative is used as NACME for nonadiabatic MD.
!>
!>        The d_IJ^a is computed by numerical differentiation using
!>        a first or second-order formula. In the case of first-order,
!>
!>        O_IJ = <I(a)|J(a+da)> = F,
!>        O_IJ = <I(a+da)|J(a)> = B,
!>        d_IJ^a = (F - B) / a,
!>        where O_IJ_F and O_IJ_B are the state overlaps.
!>
!>        This routine returns d_IJ^a * a = F - B.
!>        The denominator MUST be provided in subsequent calculations
!>        because it can be 2 * a in the case of second-order
!>        numerical differentiation formula. Also, a can be either a time
!>        or geometric parameter.
!>
!>        https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b01049,
!>
!>     @author Seunghoon Lee, Konstantin Komarov
!>
  subroutine get_dcv(nact, s_st, nstates)

    use precision, only: dp
    use io_constants, only: iw

    implicit none

    integer :: nstates
    real(kind=dp), dimension(nstates,*) :: nact, s_st

    integer :: i, j
    logical :: debug = .true.

    if (debug) write(iw, &
      fmt='(/x/,a,/29x," F              B          F - B")') &
      '  F = <I(a)|J(a+da)>,  B = <I(a+da)|J(a)>, where a is variable.'

    do i = 1, nstates
      do j = 1, nstates
        nact(i,j) = (s_st(i,j)-s_st(j,i))
      if (debug) write(iw, &
        fmt='(x,a,i0,a,i0,a,i0,a,i0,a,f12.8,a,f12.8,f12.8)') &
        "<S",i,"|S",j,"> and <S",j,"|S",i,"> = ",s_st(I,J), &
        " and", s_st(J,I), nact(i,j)
      end do
    end do

  end subroutine get_dcv

!>  This routine calculates the determinate of a square matrix.
!>  Gauss Elimination Method
!
!>  array    the matrix of order norder which is to be evaluated.
!>           this subprogram destroys the matrix array
!>  norder   the order of the square matrix to be evaluated.

  function comp_det(array, n) result(det)
    use precision, only: dp

    implicit none

    real(kind=dp) :: det
    real(kind=dp), intent(inout), dimension(n,n) :: array
    integer, intent(in) :: n

    real(kind=dp), dimension(n,n) :: work
    integer, dimension(n) :: num
    real(kind=dp) :: tmp, max
    integer i, k, l, m

    det = 1.0_dp
    do k = 1, n
       max = array(k,k)
       num(k) = k
       do i = k+1, n
          if(abs(max)<abs(array(i,k))) then
             max = array(i,k)
             num(k) = i
          end if
       end do
       if (num(k)/=k) then
          do l = k, n
             tmp = array(k,l)
             array(k,l) = array(num(k),l)
             array(num(k),l) = tmp
          end do
          det = -1.0_dp*det
       end if
       do m = k+1, n
          work(m,k) = array(m,k)/array(k,k)
          do l = k, n
             array(m,l) = array(m,l)-work(m,k)*array(k,l)
          end do
       end do !There we made matrix triangular!
    end do

    do i = 1, n
    det = det*array(i,i)
    end do

  end function comp_det

  subroutine print_nac(infos, state_overlap, nac)
    use precision, only: dp
    use io_constants, only: iw
    use types, only: information

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), dimension(:,:) :: state_overlap
    real(kind=dp), intent(in), dimension(:,:) :: nac

    integer :: ndtlf, max, imax, imin, i, j, nstates

    ndtlf = infos%tddft%tlf
    nstates = infos%tddft%nstate

    write (iw, fmt='(/5x,40(1h-)/&
               &5x,"state overlap integral between different"/ &
               &5x,"   time steps by using TLF(",i0,") approx"/ &
               &5x,"     (<phi^{i}(t-dt)|phi^{j}(t)>)"/ &
               &5x,40(1h-))') ndtlf
    max = 10
    imax = 0
    do
       imin = imax+1
       imax = imax+max
       if (imax>=nstates) imax = nstates
       write (iw,fmt='(5x,10(4x,i4,3x))') (i, i = imin, imax)
       do j = 1, nstates
          write (iw, fmt='(i5,10f11.6)') j,(state_overlap(j,i),i = imin,imax)
       end do
       if (imax>nstates) then
          cycle
       else
          exit
       end if
    end do

    write(iw, fmt='(/,3(/5x,a),/9x,a/,5x,a/)') &
          "---------------------------------", &
          "Derivative Coupling Term (a.u.)", &
          "by using finite difference approx", &
          "(<phi^{i}|d/dt|phi^{j}> = F - B)", &
          "---------------------------------"
    do j = 1, nstates
       write(iw, fmt='(i5,10f11.6)') j, (nac(j,i)*0.02418884254, i = 1, nstates)
    end do
    write (iw,*) " "

  end subroutine

end module get_state_overlap_mod
