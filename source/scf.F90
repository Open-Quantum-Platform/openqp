module scf
  use precision, only: dp

  implicit none

  character(len=*), parameter :: module_name = "scf"

  public :: scf_driver
  public :: fock_jk

contains

  subroutine scf_driver(basis, infos, molGrid)

  ! Main drirver for HF and DFT with RHF, ROHF and UHF
  ! DIIS is the main algorithm for SCF

     USE precision, only: dp
     use oqp_tagarray_driver
     use types, only: information
     use int2_compute, only: int2_compute_t, int2_fock_data_t, int2_rhf_data_t, int2_urohf_data_t
     use dft, only: dftexcor
     use mod_dft_molgrid, only: dft_grid_t
     use messages,  only: show_message, WITH_ABORT
     use guess, only: get_ab_initio_density, get_ab_initio_orbital
     use util, only: measure_time, e_charge_repulsion
     use printing, only: print_mo_range
     use mathlib, only: traceprod_sym_packed, matrix_invsqrt
     use mathlib, only: unpack_matrix
     use io_constants, only: IW
     use basis_tools, only: basis_set
     use scf_converger, only: scf_conv_result, scf_conv, &
             conv_cdiis, conv_ediis

     implicit none

     character(len=*), parameter :: subroutine_name = "scf_driver"

     type(basis_set), intent(in) :: basis
     type(information), target, intent(inout) :: infos
     type(dft_grid_t), intent(in) :: molGrid
     integer :: i, ii, iter, nschwz, nbf, nbf_tri, nbf2, ok, maxit
     real(kind=dp) :: ehf, ehf1, nenergy, etot, diffe, e_old, psinrm, &
                 scalefactor,vne, vnn, vee, vtot, virial, tkin
     real(kind=dp), allocatable :: tempvec(:), lwrk(:), lwrk2(:)
     real(kind=dp), allocatable, target :: smat_full(:,:), pdmat(:,:), pfock(:,:), rohf_bak(:)
     real(kind=dp), allocatable, target :: dold(:,:), fold(:,:)
     integer :: nfocks, diis_reset
!    For DIIS
     real(kind=dp) :: diis_error
     real(kind=dp), parameter :: ethr_cdiis_big = 2.0_dp
     real(kind=dp), parameter :: ethr_ediis = 1.0_dp
     integer :: diis_nfocks, diis_stat, maxdiis
     character(len=6), dimension(5) :: diis_name
!    Vshift
     real(kind=dp) :: vshift, H_U_gap, H_U_gap_crit
     logical :: vshift_last_iter

     type(scf_conv) :: conv
     class(scf_conv_result), allocatable :: conv_res
     character(16) :: scf_name = ""
     real(kind=dp) :: eexc, totele, totkin
     logical :: is_dft, diis_reset_condition
     integer :: scf_type
!    MOM
     logical :: do_mom, step_0_mom, do_mom_flag
     real(kind=dp), allocatable, dimension(:,:) :: mo_a_for_mom, mo_b_for_mom, work

     integer :: nelec, nelec_a, nelec_b
     integer, parameter :: scf_rhf  = 1, &
                           scf_uhf  = 2, &
                           scf_rohf = 3
     real(kind=dp), allocatable :: pfxc(:,:), qmat(:,:)

     type(int2_compute_t) :: int2_driver
     class(int2_fock_data_t), allocatable :: int2_data
!    pFON 
     logical :: do_pfon
     real(kind=dp) :: beta_pfon, start_temp, end_temp, temp_pfon
     real(kind=dp) :: electron_sum_a, electron_sum_b, pfon_start_temp 
     real(kind=dp), allocatable :: occ_a(:), occ_b(:)
     real(kind=dp) :: sum_occ_alpha, sum_occ_beta, cooling_rate
     real(kind=dp) :: pfon_cooling_rate
     real(kind=dp), parameter :: kB_HaK = 3.166811563e-6_dp
  ! tagarray
     real(kind=dp), contiguous, pointer :: &
       dmat_a(:), dmat_b(:), fock_a(:), fock_b(:), hcore(:), mo_b(:,:), &
       smat(:), tmat(:), mo_a(:,:), &
       mo_energy_b(:), mo_energy_a(:), mo_energy_a_for_mom(:)
     character(len=*), parameter :: tags_general(3) = (/ character(len=80) :: &
       OQP_SM, OQP_TM, OQP_Hcore /)
     character(len=*), parameter :: tags_alpha(4) = (/ character(len=80) :: &
       OQP_FOCK_A, OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
     character(len=*), parameter :: tags_beta(4) = (/ character(len=80) :: &
       OQP_FOCK_B, OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)

  !  Default values

  !  MOM settings
  !  Current MOM option works for both RHF and ROHF
     do_mom = infos%control%mom
  !  Vshift settings
  !  Current VSHIFT only works for ROHF
     vshift = infos%control%vshift
     vshift_last_iter=.false.
     H_U_gap_crit=0.02_dp

  !  pFON settings
     do_pfon = .false. 
     do_pfon = infos%control%pfon 
     start_temp = infos%control%pfon_start_temp
     cooling_rate = infos%control%pfon_cooling_rate
     if (start_temp <= 0.0_dp) then 
         start_temp = 2000.0_dp 
     end if 
     temp_pfon  = start_temp
     if (temp_pfon < 1.0_dp) temp_pfon = 1.0_dp

     beta_pfon = 1.0_dp / (kB_HaK * temp_pfon)

  !  DIIS options
  !  none IS NOT recommended!
  !  c-DIIS: Default commutator DIIS
  !  e-DIIS: Energy-Weighted DIIS
  !  a-DIIS: Augmented DIIS
  !  v-DIIS: Vshift with DIIS, which will be eventually become c-DIIS
  !  MOM: Maximum Overlap Method for better convergency

     diis_error = 2.0_dp
     diis_name = [character(len=6) :: "none", "c-DIIS","e-DIIS","a-DIIS","v-DIIS"]
  !  Read calculation metadata from infos
     is_dft = infos%control%hamilton >= 20
     select case (infos%control%scftype)
     case (1)
       scf_type = scf_rhf
     case (2)
       scf_type = scf_uhf
     case (3)
       scf_type = scf_rohf
     end select
  !  Total number of electrons
     nelec = infos%mol_prop%nelec
     nelec_a = infos%mol_prop%nelec_a
     nelec_b = infos%mol_prop%nelec_b
  !  Matrix size
     nbf = basis%nbf
     nbf_tri = nbf*(nbf+1)/2
     nbf2 = nbf*nbf
     maxit = infos%control%maxit
     maxdiis = infos%control%maxdiis
     diis_reset = infos%control%diis_reset_mod
  !  DFT
     if (.not.is_dft) then
         scalefactor = 1.0_dp
     else
         scalefactor = infos%dft%HFscale
     end if
  !  SCF Type
     select case (scf_type)
     case (scf_rhf)
       scf_name = "RHF"
       nfocks = 1
       diis_nfocks = 1
     case (scf_uhf)
       scf_name = "UHF"
       nfocks = 2
       diis_nfocks = 2
     case (scf_rohf)
       scf_name = "ROHF"
       nfocks = 2
       diis_nfocks = 1
     end select

  !  Now we are allocating dynamic memories of tag arrays
  !  Tag arrays
     call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
     call tagarray_get_data(infos%dat, OQP_Hcore, hcore)
     call tagarray_get_data(infos%dat, OQP_SM, smat)
     call tagarray_get_data(infos%dat, OQP_TM, tmat)

     call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
     call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
     call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
     call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
     call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

     if (nfocks > 1) then
       call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
       call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
       call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
       call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
       call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
     endif

     allocate(smat_full(nbf,nbf), pdmat(nbf_tri,nfocks), pfock(nbf_tri,nfocks), &
              qmat(nbf,nbf), &
              stat=ok, &
              source=0.0_dp)
  !  Alloc work matrices for MOM
     if (do_mom) then
        step_0_mom = .true.
        do_mom_flag=.false.
        allocate(mo_a_for_mom(nbf,nbf), &
                 mo_b_for_mom(nbf,nbf), &
                 mo_energy_a_for_mom(nbf), &
                 work(nbf,nbf), &
                 source=0.0_dp)
     end if
     if(ok/=0) call show_message('Cannot allocate memory for SCF',WITH_ABORT)
  !  Incremental Fock building for Direct SCF
     if (infos%control%scf_incremental /= 0) then
       allocate( dold(nbf_tri,nfocks), fold(nbf_tri,nfocks), &
                stat=ok, &
                source=0.0_dp)
       if(ok/=0) call show_message('Cannot allocate memory for SCF',WITH_ABORT)
     end if
  !  DFT
     if (is_dft) then
       allocate(pfxc(nbf_tri,nfocks), &
                stat=ok, &
                source=0.0_dp)
       if(ok/=0) call show_message('Cannot allocate memory for temporary vectors',WITH_ABORT)
     end if
  !  ROHF temporary arrays
     if (scf_type == scf_rohf .or. is_dft) then
       allocate(tempvec(nbf2), &
                stat=ok, &
                source=0.0_dp)
       if(ok/=0) call show_message('Cannot allocate memory for temporary vectors',WITH_ABORT)
     end if
     if (scf_type == scf_rohf) then
       allocate(lwrk(nbf), lwrk2(nbf2), rohf_bak(nbf_tri), &
                stat=ok, &
                source=0.0_dp)
       if(ok/=0) call show_message('Cannot allocate memory for ROHF temporaries',WITH_ABORT)
     end if

  !  Allocating dynamic memories done

     call measure_time(print_total=1, log_unit=iw)

     call matrix_invsqrt(smat, qmat, nbf)

  !  First Compute Nuclear-Nuclear energy
     nenergy = e_charge_repulsion(infos%atoms%xyz, infos%atoms%zn)

  !  During guess, the Hcore, Q nd Overlap matrices were formed.
  !  Using these, the initial orbitals (VEC) and density (Dmat) were subsequently computed.
  !  Now we are going to calculate ERI(electron repulsion integrals) to form a new FOCK
  !  matrix.

  !  Vshift settings
  !   vshift = infos%control%vshift
  !   vshift_last_iter=.false.
  !   H_U_gap_crit=0.02_dp

  !  Initialize ERI calculations
     call int2_driver%init(basis, infos)
     call int2_driver%set_screening()
     call flush(iw)

  !  The main Do loop of SCF
     select case (scf_type)
     case (scf_rhf)
         pdmat(:,1) = dmat_a
         allocate(int2_rhf_data_t :: int2_data)
         int2_data = int2_rhf_data_t(nfocks=1, d=pdmat, scale_exchange=scalefactor)
     case (scf_uhf, scf_rohf)
         pdmat(:,1) = dmat_a
         pdmat(:,2) = dmat_b
         allocate(int2_urohf_data_t :: int2_data)
         int2_data = int2_urohf_data_t(nfocks=2, d=pdmat, scale_exchange=scalefactor)
     end select
     call unpack_matrix(smat,smat_full,nbf,'U')

!    Variable DIIS options:
!      a) If we use VSHIFT option, the combination of e-DIIS and c-DIIS is used,
!      since e-DIIS can better reduce total energy than c-DIIS.
!      Once sufficiently converged, c-DIIS is used to finalize the SCF.
!      b) If VSHIFT is not set, c-DIIS is default.
!      c) If vdiis (diis_type=5) is chosen, the VSHIFT is initally turned on.
!      d) if MOM=.true., MOM turns on if DIIS error < mom_switch
     if (infos%control%diis_type == 5) then
        call conv%init(ldim=nbf, maxvec=maxdiis, &
             subconvergers=[conv_cdiis, conv_ediis, conv_cdiis], &
             thresholds=[ethr_cdiis_big, ethr_ediis, infos%control%vdiis_cdiis_switch], &
             overlap=smat_full, overlap_sqrt=qmat, &
             num_focks=diis_nfocks, verbose=1)
        if (infos%control%vshift == 0.0_dp) then
           infos%control%vshift=0.1_dp
           vshift=0.1_dp
           call show_message("")
           call show_message('Setting VSHIFT=0.1, since VDIIS is chosen without VSHIFT value.')
        endif
     elseif (infos%control%vshift /= 0.0_dp) then
        call conv%init(ldim=nbf, maxvec=maxdiis, &
             subconvergers=[conv_cdiis, conv_ediis, conv_cdiis], &
             thresholds=[ethr_cdiis_big, ethr_ediis, infos%control%vshift_cdiis_switch], &
             overlap=smat_full, overlap_sqrt=qmat, &
             num_focks=diis_nfocks, verbose=1)
     else
!    Normally, c-DIIS works best. But one can choose others (e-DIIS and a-DIIS).
        call conv%init(ldim=nbf, maxvec=maxdiis, &
             subconvergers=[infos%control%diis_type], &
             thresholds=[infos%control%diis_method_threshold], &
             overlap=smat_full, overlap_sqrt=qmat, &
             num_focks=diis_nfocks, verbose=1)
     endif

     eexc = 0.0_dp
     e_old = 0.0_dp

  !  SCF Options
     write(iw,'(/5X,"SCF options"/ &
                &5X,18("-")/ &
                &5X,"SCF type = ",A,5x,"MaxIT = ",I5/, &
                &5X,"MaxDIIS = ",I5,17x,"Conv = ",F14.10/, &
                &5X,"DIIS Type = ",A/, &
                &5X,"vDIIS_cDIIS_Switch = ",F8.5,3x,"vDIIS_vshift_Switch = ",F8.5/, &
                &5X,"DIIS Reset Mod = ",I5,10x,"DIIS Reset Conv = ",F12.8/, &
                &5X,"VShift = ",F8.5,15X,"VShift_cDIIS_Switch = ",F8.5)') &
                & scf_name, infos%control%maxit, &
                & infos%control%maxdiis, infos%control%conv, &
                & diis_name(infos%control%diis_type), &
                & infos%control%vdiis_cdiis_switch, infos%control%vdiis_vshift_switch, &
                & infos%control%diis_reset_mod, infos%control%diis_reset_conv, &
                & infos%control%vshift, infos%control%vshift_cdiis_switch
     write(iw,'(5X,"MOM = ",L5,21X,"MOM_Switch = ",F8.5)') &
                & infos%control%mom, infos%control%mom_switch 
     write(iw,'(5X,"pFON = ",L5,20X,"pFON Start Temp. = ",F9.2,/, &
               5X, "pFON Cooling Rate = ", F9.2)') &
               infos%control%pfon, infos%control%pfon_start_temp, &
               infos%control%pfon_cooling_rate

  !  Initial message
     write(IW,fmt="&
          &(/3x,'Direct SCF iterations begin.'/, &
          &  3x,93('='),/ &
          &  4x,'Iter',9x,'Energy',12x,'Delta E',9x,'Int Skip',5x,'DIIS Error',5x,'Shift',5x,'Method'/ &
          &  3x,93('='))")
     call flush(iw)

     do iter = 1, maxit

  !     The main SCF iteration loop

  !     pFON Cooling
        if (cooling_rate <= 0.0_dp) then
            cooling_rate = 50_dp 
        end if 
        if (do_pfon) then
            if ( (iter == maxit ) .or. (abs(diis_error) < 10.0_dp * infos%control%conv) ) then 
                temp_pfon = 0.0_dp 
            else 
                temp_pfon = temp_pfon - cooling_rate 
                if (temp_pfon < 1.0_dp) temp_pfon = 1.0_dp 
            end if 
            if (temp_pfon > 1.0e-12_dp) then 
                beta_pfon = 1.0_dp / (kB_HaK * temp_pfon)
            else 
                beta_pfon = 1.0e20_dp 
            end if
        end if 
 
        pfock = 0.0_dp

  !     Compute difference density matrix for incremental Fock build,
  !     which is the main advantage of direct SCF.
  !     It will provide much better ERI screening.

        if (infos%control%scf_incremental /= 0) then
          pdmat = pdmat - dold
        end if

        call int2_driver%run(int2_data, &
                cam=is_dft.and.infos%dft%cam_flag, &
                alpha=infos%dft%cam_alpha, &
                beta=infos%dft%cam_beta,&
                mu=infos%dft%cam_mu)
        nschwz = int2_driver%skipped

  !     Recover full Fock and density from difference matrices
        if (infos%control%scf_incremental /= 0) then
          pdmat = pdmat + dold
          int2_data%f(:,:,1) = int2_data%f(:,:,1) + fold
          fold = int2_data%f(:,:,1)
          dold = pdmat
        end if
  !     Scaling
        pfock(:,:) = 0.5_dp * int2_data%f(:,:,1)
        ii=0
        do i = 1, nbf
           ii = ii + i
           pfock(ii,1:nfocks) = 2.0_dp*pfock(ii,1:nfocks)
        end do

  !     Adding the skeleton H core to Fock for getting the new orbitals
        do i = 1, nfocks
          pfock(:,i)  = pfock(:,i) + hcore
        end do

  !     After this, we compute Energy.
        ehf = 0.0_dp
        ehf1 = 0.0_dp
        do i = 1, nfocks
          ehf1 = ehf1 + traceprod_sym_packed(pdmat(:,i),hcore,nbf)
          ehf = ehf + traceprod_sym_packed(pdmat(:,i),pfock(:,i),nbf)
        end do
        ehf = 0.5_dp*(ehf+ehf1)
        etot = ehf + nenergy
  !     DFT contribution
        if (is_dft) then
          if (scf_type == scf_rhf) then
            call dftexcor(basis,molgrid,1,pfxc,pfxc,mo_a,mo_a,nbf,nbf_tri,eexc,totele,totkin,infos)
          else if (scf_type == scf_uhf) then
             call dftexcor(basis,molgrid,2,pfxc(:,1),pfxc(:,2),mo_a,mo_b,nbf,nbf_tri,eexc,totele,totkin,infos)
          else if (scf_type == scf_rohf) then
!            ROHF does not have MO_B. So we copy MO_A to MO_B.
             mo_b = mo_a
             call dftexcor(basis,molgrid,2,pfxc(:,1),pfxc(:,2),mo_a,mo_b,nbf,nbf_tri,eexc,totele,totkin,infos)
          end if
          pfock = pfock + pfxc
          etot = etot + eexc
        end if

  !     Forming ROHF Fock by combing Alpha and Beta Focks.
        if (scf_type == scf_rohf) then
           rohf_bak = pfock(:,1)
           if (vshift_last_iter.eqv..true.) vshift = 0.0_dp
           call form_rohf_fock(pfock(:,1),pfock(:,2),tempvec, &
             mo_a,smat,lwrk2,nelec_a,nelec_b,nbf,vshift)
           pdmat(:,1) = pdmat(:,1) + pdmat(:,2)
        end if

  !     SCF Converger to get refined Fock Matrix
        call conv%add_data(f=pfock(:,1:diis_nfocks), &
                dens=pdmat(:,1:diis_nfocks), e=Etot)

  !     Run DIIS calculation, get DIIS error
        call conv%run(conv_res)
        diis_error = conv_res%error
  !     Checking the convergency
        diffe = etot-e_old
        write(iw,'(4x,i4.1,2x,f17.10,1x,f17.10,1x,i13,1x,f14.8,5x,f5.3,5x,a)') &
                iter, etot, diffe, nschwz, &
                diis_error, vshift, trim(conv_res%active_converger_name)
        call flush(iw)
  !     VDIIS option
        if ((infos%control%diis_type.eq.5) &
           .and.(diis_error < infos%control%vdiis_vshift_switch)) then
           vshift=0.0_dp
        elseif ((infos%control%diis_type.eq.5) &
           .and.(diis_error.ge.infos%control%vdiis_vshift_switch)) then
           vshift=infos%control%vshift
        endif

        e_old = etot

  !     Exit if convergence criteria achieved
        if ((abs(diis_error)<infos%control%conv).and.(vshift==0.0_dp)) then
           exit
        elseif ((abs(diis_error)<infos%control%conv).and.(vshift/=0.0_dp)) then
           write(iw,"(3x,64('-')/10x,'Performing a last SCF with zero VSHIFT.')")
           vshift_last_iter=.true.
        elseif (vshift_last_iter.eqv..true.) then
  !        Only for ROHF case
           call get_ab_initio_orbital(pfock(:,1),mo_a,mo_energy_a,qmat)
           exit
        endif

  !     Reset DIIS for difficult case

  !     VSHIFT=0 and slow cases
        diis_reset_condition=(((iter/diis_reset).ge.1) &
           .and.(modulo(iter,diis_reset).eq.0) &
           .and.(diis_error.gt.infos%control%diis_reset_conv) &
           .and.(infos%control%vshift==0.0_dp))
        if (diis_reset_condition) then
!          Resetting DIIS
           write(iw,"(3x,64('-')/10x,'Resetting DIIS.')")
           call conv_res%get_fock(matrix=pfock(:,1:diis_nfocks), istat=diis_stat)
           call conv%init(ldim=nbf, maxvec=maxdiis, &
             subconvergers=[conv_cdiis], &
             thresholds=[ethr_cdiis_big], &
             overlap=smat_full, overlap_sqrt=qmat, &
             num_focks=diis_nfocks, verbose=1)
!          After resetting DIIS, we need to skip SD
           call conv%add_data(f=pfock(:,1:diis_nfocks), &
                dens=pdmat(:,1:diis_nfocks), e=Etot)
           call conv%run(conv_res)
        else
  !        Form the interpolated the Fock/density matrix
           call conv_res%get_fock(matrix=pfock(:,1:diis_nfocks), istat=diis_stat)
        endif

  !     Calculate new orbitals and density.
  !
        if (int2_driver%pe%rank == 0) then
           call get_ab_initio_orbital(pfock(:,1),mo_a,mo_energy_a,qmat)

           if (scf_type == scf_uhf .and. nelec_b /= 0) then
     !        Only UHF has beta orbitals.
               call get_ab_initio_orbital(pfock(:,2),mo_b,mo_energy_b,qmat)
           end if
        end if
        if (scf_type == scf_uhf .and. nelec_b /= 0) then
            call int2_driver%pe%bcast(mo_b, size(mo_b))
            call int2_driver%pe%bcast(mo_energy_b, size(mo_energy_b))
        end if
        call int2_driver%pe%bcast(mo_a, size(mo_a))
        call int2_driver%pe%bcast(mo_energy_a, size(mo_energy_a))

        ! pFON section
        do_pfon = infos%control%pfon
        if (do_pfon) then
            if (.not. allocated(occ_a)) allocate(occ_a(nbf))
            if (.not. allocated(occ_b)) allocate(occ_b(nbf))

            select case (scf_type)
            case (scf_rhf)
                call pfon_occupations(mo_energy_a, nbf, nelec, occ_a, beta_pfon, scf_type)

            case (scf_uhf)
                call pfon_occupations(mo_energy_a, nbf, nelec_a, occ_a, beta_pfon, scf_type)
                if (nelec_b > 0) then
                    call pfon_occupations(mo_energy_b, nbf, nelec_b, occ_b, beta_pfon, scf_type)
                end if

            case (scf_rohf)
                call pfon_occupations(mo_energy_a, nbf, nelec_a, occ_a, beta_pfon, scf_type)
                if (nelec_b > 0) then
                    call pfon_occupations(mo_energy_a, nbf, nelec_b, occ_b, beta_pfon, scf_type)
                end if
            end select

            ! (Alpha)
            sum_occ_alpha = 0.0_dp
            do i = 1, nbf
                sum_occ_alpha = sum_occ_alpha + occ_a(i)
            end do

            ! (Beta)
            sum_occ_beta = 0.0_dp
            if (scf_type == scf_uhf .or. scf_type == scf_rohf) then
                do i = 1, nbf
                    sum_occ_beta = sum_occ_beta + occ_b(i)
                end do
            end if

            write(iw,'(T7, " pFON: Temp=",F9.2,", Beta=",ES11.4)') &
                 temp_pfon, beta_pfon

!            write(iw,'(" Start: ",F9.2,", END: Temp=",F9.2,", Elect Sum(a)=",F8.3,", Elect Sum(b)=",F8.3)') &
!                  start_temp ,end_temp, electron_sum_a, electron_sum_b
        end if


  !     MOM option works for RHF and ROHF
        if (do_mom .and. diis_error.lt.infos%control%mom_switch) do_mom_flag=.true.
        if (do_mom .and. do_mom_flag .and. .not. step_0_mom) then
           call mo_reorder(infos, mo_a_for_mom, mo_energy_a_for_mom, &
                           mo_a, mo_energy_a, smat_full)
        end if
        if (do_mom) then
           mo_a_for_mom = mo_a
           mo_energy_a_for_mom = mo_energy_a
        end if
        step_0_mom = .false.

  !     New density matrix in AO basis using MO.
        if (int2_driver%pe%rank == 0) then
            if (.not. do_pfon) then 
                call get_ab_initio_density(pdmat(:,1),mo_a,pdmat(:,2),mo_b,infos,basis)
            else 
                call build_pfon_density(pdmat, mo_a, mo_b, occ_a, occ_b, scf_type, nbf, nelec_a, nelec_b)
            end if
        end if 
        call int2_driver%pe%bcast(pdmat, size(pdmat))
 !     Checking the HOMO-LUMO gaps for predicting SCF convergency
        if ((iter > 10).and.(vshift==0.0_dp)) then
           select case (scf_type)
           case (scf_rhf)
              H_U_gap=mo_energy_a(nelec/2)-mo_energy_a(nelec/2-1)
           case (scf_uhf)
              H_U_gap=mo_energy_a(nelec/2)-mo_energy_a(nelec/2-1)
           case (scf_rohf)
              H_U_gap=mo_energy_a(nelec_a+1)-mo_energy_a(nelec_a)
           end select
        endif

  !  End of Iteration
     end do

     if (iter>maxit) then
       write(iw,"(3x,64('-')/10x,'SCF is not converged ....')")
       infos%mol_energy%SCF_converged=.false.
     else
       write(iw,"(3x,64('-')/10x,'SCF convergence achieved ....')")
       infos%mol_energy%SCF_converged=.true.
     end if

     write(iw,"(/' Final ',A,' energy is',F20.10,' after',I4,' iterations'/)") trim(scf_name), etot, iter

     if (is_dft) then
       write(iw,*)
       write(iw,"(' DFT: XC energy              = ',F20.10)") eexc
       write(iw,"(' DFT: total electron density = ',F20.10)") totele
       write(iw,"(' DFT: number of electrons    = ',I9,/)") nelec
     end if
  !
     if (scf_type == scf_uhf .and. nelec_b /= 0) then
         call int2_driver%pe%bcast(mo_b, size(mo_b))
         call int2_driver%pe%bcast(mo_energy_b, size(mo_energy_b))
     end if

     call int2_driver%pe%bcast(mo_a, size(mo_a))
     call int2_driver%pe%bcast(pdmat, size(pdmat))
     call int2_driver%pe%bcast(mo_energy_a, size(mo_energy_a))

     select case (scf_type)
     case (scf_rhf)
       fock_a = pfock(:,1)
       dmat_a = pdmat(:,1)
     case (scf_uhf)
       fock_a = pfock(:,1)
       fock_b = pfock(:,2)
       dmat_a = pdmat(:,1)
       dmat_b = pdmat(:,2)
     case (scf_rohf)
       fock_a = rohf_bak
       call mo_to_ao(fock_b, pfock(:,2), smat, mo_a, nbf, nbf)
       dmat_a = pdmat(:,1) - pdmat(:,2)
       dmat_b = pdmat(:,2)
       mo_b = mo_a
       mo_energy_b = mo_energy_a
     end select

  !  Print out the molecular orbitals
     call print_mo_range(basis, infos, mostart=1, moend=nbf)

  !  Print out the results
     psinrm = 0.0_dp
     tkin = 0.0_dp
     do i = 1, diis_nfocks
       psinrm = psinrm + traceprod_sym_packed(pdmat(:,i),smat,nbf)/nelec
       tkin = tkin + traceprod_sym_packed(pdmat(:,i),tmat,nbf)
     end do

  !  Writing out final results
     vne = ehf1 - tkin
     vee = etot - ehf1 - nenergy
     vnn = nenergy
     vtot = vne + vnn + vee
     virial = - vtot/tkin
     call print_scf_energy(psinrm, ehf1, nenergy, etot, vee, vne, vnn, vtot, tkin, virial)

  !  Save results to infos.
     infos%mol_energy%energy = etot
     infos%mol_energy%psinrm = psinrm
     infos%mol_energy%ehf1 = ehf1
     infos%mol_energy%vee = vee
     infos%mol_energy%nenergy = nenergy
     infos%mol_energy%vne = vne
     infos%mol_energy%vnn = vnn
     infos%mol_energy%vtot = vtot
     infos%mol_energy%tkin = tkin
     infos%mol_energy%virial = virial
     infos%mol_energy%energy = etot

  !  Clean up
     call int2_driver%clean()

     call measure_time(print_total=1, log_unit=iw)

  end subroutine scf_driver

  subroutine print_scf_energy(psinrm, ehf1, enuclear, etot, vee, vne, vnn, vtot, tkin, virial)
     use precision, only: dp
     use io_constants, only: iw
     implicit none
     real(kind=dp) :: psinrm, ehf1, enuclear, etot, vee, vne, vnn, vtot, tkin, virial

     write(iw,"(/10X,17('=')/10X,'Energy components'/10X,17('=')/)")
     write(iw,"('         Wavefunction normalization =',F19.10)") psinrm
     write(iw,*)
     write(iw,"('                One electron energy =',F19.10)") ehf1
     write(iw,"('                Two electron energy =',F19.10)") vee
     write(iw,"('           Nuclear repulsion energy =',F19.10)") enuclear
     write(iw,"(38X,18('-'))")
     write(iw,"('                       TOTAL energy =',F19.10)") etot
     write(iw,*)
     write(iw,"(' Electron-electron potential energy =',F19.10)") vee
     write(iw,"('  Nucleus-electron potential energy =',F19.10)") vne
     write(iw,"('   Nucleus-nucleus potential energy =',F19.10)") vnn
     write(iw,"(38X,18('-'))")
     write(iw,"('             TOTAL potential energy =',F19.10)") vtot
     write(iw,"('               TOTAL kinetic energy =',F19.10)") tkin
     write(iw,"('                 Virial ratio (V/T) =',F19.10)") virial
     write(iw,*)
  end subroutine print_scf_energy

!> @brief Form the ROHF Fock matrix in MO basis.
!>
!> This implementation is based on the work of M. F. Guest and V. Saunders,
!> as described in Mol. Phys. 28, 819 (1974).
!>
!> @brief Form the ROHF Fock matrix in the MO basis using the Guest-Saunders method.
!>
!> @detail
!> This subroutine transforms the alpha and beta Fock matrices from the AO to the MO basis,
!> constructs the ROHF Fock matrix using the Guest-Saunders method, and optionally applies
!> a level shift to the diagonal elements of the virtual-virtual block.
!>
!> @author Konstantin Komarov, 2023
!>
  subroutine form_rohf_fock(fock_a_ao, fock_b_ao, fock_mo, &
                            MOs, overlap_tri, work, &
                            nocca, noccb, nbf, vshift)

    use mathlib, only: orthogonal_transform_sym, &
                       orthogonal_transform2, &
                       triangular_to_full, &
                       unpack_matrix, &
                       pack_matrix

    implicit none

    real(kind=dp), intent(inout), dimension(:) :: fock_a_ao
    real(kind=dp), intent(inout), dimension(:) :: fock_b_ao
    real(kind=dp), intent(out), dimension(:) :: fock_mo
    real(kind=dp), intent(in), dimension(:,:) :: MOs
    real(kind=dp), intent(in), dimension(:) :: overlap_tri
    real(kind=dp), intent(out), dimension(:) :: work
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: vshift

    real(kind=dp), allocatable, dimension(:,:) :: &
          overlap, work_matrix, fock, fock_a, fock_b
    real(kind=dp) :: acc, aoo, avv, bcc, boo, bvv
    integer :: i, nbf_tri

    acc = 0.5_dp; aoo = 0.5_dp; avv = 0.5_dp
    bcc = 0.5_dp; boo = 0.5_dp; bvv = 0.5_dp
    nbf_tri = nbf * (nbf + 1) / 2

    ! Allocate full matrices
    allocate(overlap(nbf, nbf), &
             work_matrix(nbf, nbf), &
             fock(nbf, nbf), &
             fock_a(nbf, nbf), &
             fock_b(nbf, nbf), &
             source=0.0_dp)

    ! Transform alpha and beta Fock matrices to MO basis
    call orthogonal_transform_sym(nbf, nbf, fock_a_ao, MOs, nbf, fock_mo)
    fock_a_ao(:nbf_tri) = fock_mo(:nbf_tri)

    call orthogonal_transform_sym(nbf, nbf, fock_b_ao, MOs, nbf, fock_mo)
    fock_b_ao(:nbf_tri) = fock_mo(:nbf_tri)

    ! Unpack triangular matrices to full matrices
    call unpack_matrix(fock_a_ao, fock_a)
    call unpack_matrix(fock_b_ao, fock_b)

    ! Construct ROHF Fock matrix in MO basis using Guest-Saunders method
    associate ( na => nocca &
              , nb => noccb &
      )
      fock(1:nb, 1:nb) = acc * fock_a(1:nb, 1:nb) &
                       + bcc * fock_b(1:nb, 1:nb)
      fock(nb+1:na, nb+1:na) = aoo * fock_a(nb+1:na, nb+1:na) &
                             + boo * fock_b(nb+1:na, nb+1:na)
      fock(na+1:nbf, na+1:nbf) = avv * fock_a(na+1:nbf, na+1:nbf) &
                               + bvv * fock_b(na+1:nbf, na+1:nbf)
      fock(1:nb, nb+1:na) = fock_b(1:nb, nb+1:na)
      fock(nb+1:na, 1:nb) = fock_b(nb+1:na, 1:nb)
      fock(1:nb, na+1:nbf) = 0.5_dp * (fock_a(1:nb, na+1:nbf) &
                                     + fock_b(1:nb, na+1:nbf))
      fock(na+1:nbf, 1:nb) = 0.5_dp * (fock_a(na+1:nbf, 1:nb) &
                                     + fock_b(na+1:nbf, 1:nb))
      fock(nb+1:na, na+1:nbf) = fock_a(nb+1:na, na+1:nbf)
      fock(na+1:nbf, nb+1:na) = fock_a(na+1:nbf, nb+1:na)

      ! Apply Vshift to the diagonal
      do i = nb+1, na
        fock(i,i) = fock(i,i) + vshift * 0.5_dp
      end do
      do i = na+1, nbf
        fock(i,i) = fock(i,i) + vshift
      end do
    end associate

    ! Pack fock to fock_mo, make it triangular
    call pack_matrix(fock, fock_mo)

    call triangular_to_full(fock_mo, nbf, 'u')

    ! Back-transform ROHF Fock matrix to AO basis
    call unpack_matrix(overlap_tri, overlap)
    call dsymm('l', 'u', nbf, nbf, &
               1.0_dp, overlap, nbf, &
                       MOs, nbf, &
               0.0_dp, work, nbf)
    call orthogonal_transform2('t', nbf, nbf, work, nbf, fock, nbf, &
                               work_matrix, nbf, overlap)
    call pack_matrix(work_matrix, fock_a_ao)

  end subroutine form_rohf_fock

!> @brief Back-transform a symmetric operator `Fmo` expressed in
!>   the MO basis `V` to the AO basis
!> @detail compute the transformation:
!>   Fao = S*V * Fmo * (SV)^T
  subroutine mo_to_ao(fao, fmo, s, v, nmo, nbf)

    use mathlib, only: pack_matrix, unpack_matrix
    use oqp_linalg

    implicit none

    real(kind=dp), intent(in) :: fmo(*), s(*), v(*)
    real(kind=dp), intent(out) :: fao(*)
    integer, intent(in) :: nmo, nbf

    integer :: nbf2
    real(kind=dp), allocatable :: sv(:,:), ftmp(:,:), wrk(:,:)

    allocate(ftmp(nbf,nbf), sv(nbf,nmo), wrk(nbf,nbf))

    call unpack_matrix(fmo, ftmp)

    call unpack_matrix(s, wrk)
    ! compute S*V
    call dsymm('l', 'u', nbf, nmo, &
               1.0_dp, wrk, nbf, &
                       v,  nbf, &
               0.0_dp, sv, nbf)

    ! compute (S * V) * Fmo
    call dsymm('r', 'u', nbf, nmo, &
               1.0d0, ftmp, nbf, &
                      sv,  nbf, &
               0.0d0, wrk, nbf)

    ! compute ((S * V) * Fmo) * (S * V)^T
    call dgemm('n', 't', nbf, nbf, nmo, &
               1.0d0, wrk, nbf, &
                      sv, nbf, &
               0.0d0, ftmp, nbf)

    nbf2 = nbf*(nbf+1)/2
    call pack_matrix(ftmp, fao(:nbf2))

  end subroutine mo_to_ao


  subroutine fock_jk(basis, d, f, scalefactor, infos)
   use precision, only: dp
   use io_constants, only: iw
   use util, only: measure_time
   use basis_tools, only: basis_set
   use types, only: information
   use int2_compute, only: int2_compute_t, int2_fock_data_t, int2_rhf_data_t, int2_urohf_data_t

   implicit none

   type(basis_set), intent(in) :: basis
   type(information), intent(inout) :: infos
   real(kind=dp), optional, intent(in) :: scalefactor

   integer :: i, ii, nf, nschwz
   real(kind=dp) :: scalef
   real(kind=dp), target, intent(in) :: d(:,:)
   real(kind=dp), intent(inout) :: f(:,:)

     type(int2_compute_t) :: int2_driver
     class(int2_fock_data_t), allocatable :: int2_data

!  Initial Settings
   scalef = 1.0d0
   if (present(scalefactor)) scalef = scalefactor

   call measure_time(print_total=1, log_unit=iw)

   write(iw,"(/3x,'Form Two-Electron J and K Fock')")

!  Initialize ERI calculations
   call int2_driver%init(basis, infos)
   call int2_driver%set_screening()
   int2_data = int2_rhf_data_t(nfocks=1, d=d, scale_exchange=scalefactor)

   call flush(iw)

!  Constructing two electron Fock matrix
   call int2_driver%run(int2_data)
   nschwz = int2_driver%skipped

!  Scaling (everything except diagonal is halved)
   f =  0.5 * int2_data%f(:,:,1)

   do nf = 1, ubound(f,2)
     ii = 0
     do i = 1, basis%nbf
        ii = ii + i
        f(ii,nf) = 2*f(ii,nf)
     end do
   end do

   call int2_driver%clean()

 end subroutine fock_jk

 subroutine mo_reorder(infos, Va, Ea, Vb, Eb, Sq)
!  Va, Ea: MO and MO energies of Previous Point
!  Vb, Eb: MO and MO energies of Current Point
!  The Vb and Eb are reordered based on the orders of Va and Ea.
!  Sq: Square Form of Overlap Matrix in AO
   use precision, only: dp
   use io_constants, only: iw
   use types, only: information

   implicit none

   type(information), intent(inout) :: infos
   real(kind=dp), intent(inout), dimension(:,:) :: Va, Vb ! (nbf, nbf)
   real(kind=dp), intent(inout), dimension(:) :: Ea, Eb    ! (nbf)
   real(kind=dp), intent(in), dimension(:,:) :: Sq    ! (nbf, nbf)

   integer :: i, loc
   real(kind=dp) :: tmp, tmp2, tmp_abs
   real(kind=dp), allocatable, dimension(:,:) :: scr, smo
   logical, allocatable, dimension(:) :: pr
   integer, allocatable, dimension(:) :: locs
   integer :: na, nbf

   write(iw, fmt='(/," MOM is on: Reordering the New MO and MO energies....")')

   nbf = infos%basis%nbf
   na = infos%mol_prop%nelec_a

   allocate(scr(nbf,nbf), &
            smo(nbf,nbf), &
            source=0.0_dp)
   allocate(pr(nbf), source=.false.)
   allocate(locs(nbf), source=0)

   call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, Va, nbf, &
                       Sq, nbf, &
               0.0_dp, scr, nbf)
   call dgemm('n', 'n', nbf, nbf, nbf, &
               1.0_dp, scr, nbf, &
                       Vb, nbf, &
               0.0_dp, smo, nbf)
!  Now Smo is overlap matrix in MO, row and column corresponds to Old and New MOs.
   do i = 1, nbf
     smo(:,i) = smo(:,i) / norm2(smo(:,i))
   end do
!  Finding out the location of maximum value in column i
   do i = 1, nbf
     loc = maxloc(abs(smo(:nbf,i)), dim=1)
     locs(i) = loc
     pr(loc) = .true.
   enddo
!  Checking printouts
   if (.not.all(pr)) then
     write(iw, fmt='(/,"   Warning")')
     write(iw, fmt='(" Some orbitals are missing in reordering")')
     write(iw, advance='no', fmt='(" Their indices are:")')
     do i = 1, nbf
       if (.not.pr(i)) write(iw, advance='no', fmt='(I4,",")') i
     end do
     write(iw,*)
     write(iw, fmt='(/,"   Error. Stop")')
     write(iw, fmt='(" Some orbitals are missing in reordering")')
     write(iw,*)
   end if

   write(iw,fmt='(x,a)') 'Old MO Index ---> New MO Index'
   do i = 1, nbf
     tmp_abs = maxval(abs(smo(:nbf,i)), dim=1)
     loc = maxloc(abs(smo(:nbf,i)), dim=1)
     tmp = smo(loc,i)
     tmp2 = smo(i,i)

     if ((loc.ne.i).and.(i.le.na+1)) then
       write(iw, advance='no', &
        fmt='(5x,i4,14x,i4)') i, loc
       if (i == na-1) write(iw, advance='no',fmt='(2x,a)') ' HOMO'
       if (i == na) write(iw, advance='no',fmt='(2x,a)') ' LUMO'
       if (i /= loc .and. tmp_abs < 0.9d+00) then
         write(iw, fmt='(2x,a)') &
         ' rearranged, WARNING'
       elseif (i == loc .and. tmp_abs < 0.9d+00) then
         write(iw, fmt='(2x,a)')  &
         ' WARNING'
       elseif (i /= loc .and. tmp_abs > 0.9d+00) then
         write(iw, fmt='(2x,a)') &
         ' rearranged'
       else
         write(iw,*)
       end if
     endif
   enddo

!  Reordering MO and MO energies of Current point.
   call reorderMOs(Vb, Eb, Smo, nbf, nbf, 1, na+1)

 end subroutine mo_reorder

!> @brief      pFON Implementation in SCF Module
!> Author: Alireza Lashkaripour
!> Date: January 2025
!> Reference paper: https://doi.org/10.1063/1.478177
!> This subroutine incorporates the Partial Fractional Occupation Number (pFON) 
!> method into SCF calculations, ensuring smooth occupation numbers using 
!> Fermi-Dirac distribution. It dynamically adjusts temperature and beta 
!> factors to enhance SCF convergence, particularly for near-degenerate states.
 subroutine pfon_occupations(mo_energy, nbf, nelec, occ, beta_pfon, scf_type)
     use precision, only: dp 
     implicit none 

     integer, intent(in) :: nbf 
     integer, intent(in) :: nelec
     real(kind=dp), intent(in) :: beta_pfon
     real(kind=dp), intent(in) :: mo_energy(nbf)
     real(kind=dp), intent(inout) :: occ(nbf)
     integer, intent(in) :: scf_type ! 1,2,3 RHF,UHF,ROHF 
     real(kind=dp) :: eF, sum_occ
     integer :: i, i_homo, i_lumo
     real(kind=dp) :: tmp

     select case (scf_type)
     case(1) ! RHF 
         i_homo = max(1, nelec/2) 
     case(2) ! UHF 
         i_homo = max(1, nelec)
     case(3) ! ROHF 
         i_homo = max(1, nelec) 
     end select 

     i_lumo = i_homo + 1 
     if (i_lumo > nbf) i_lumo = nbf 

     ! Fermi level!
     eF = 0.5_dp * (mo_energy(i_homo) + mo_energy(i_lumo))

     ! pre-normalizrion occupation 
     do i = 1, nbf
        tmp = beta_pfon * (mo_energy(i) - eF)
        occ(i) = 1.0_dp / (1.0_dp + exp(tmp))
     end do 

     ! Re-normalization  
     sum_occ = 0.0_dp 
     do i = 1, nbf 
        sum_occ = sum_occ + occ(i) 
     end do 
     if (sum_occ < 1.0e-14_dp) then
        sum_occ = 1.0_dp 
     end if 
     do i = 1, nbf 
        occ(i) = occ(i) * (real(nelec,dp) / sum_occ)
     end do 
 end subroutine pfon_occupations
 
 ! pFON Density
 subroutine build_pfon_density(pdmat, mo_a, mo_b, occ_a, occ_b, scf_type, nbf, nelec_a, nelec_b)
    use precision, only: dp 
    use mathlib, only: pack_matrix
    implicit none 

    real(kind=dp), intent(inout) :: pdmat(:,:)
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), intent(in) :: occ_a(:), occ_b(:)
    integer, intent(in) :: nbf, scf_type, nelec_a, nelec_b

    real(kind=dp), allocatable :: dtmp(:,:), cmo(:,:)
    integer :: i, mu, nu 

    allocate(dtmp(nbf, nbf), source=0.0_dp)
    pdmat(:,1) = 0.0_dp 
    if (size(pdmat,2) > 1) pdmat(:,2) = 0.0_dp 

    ! build alpha density 
    do i = 1, nbf 
        if (occ_a(i) > 1.0e-14_dp) then 
            do mu = 1, nbf
                do nu = 1, nbf 
                    dtmp(mu,nu) = dtmp(mu,nu) + occ_a(i) * mo_a(mu,i)*mo_a(nu,i)
                end do 
            end do 
        end if 
    end do 

    call pack_matrix(dtmp, pdmat(:,1))

    ! bulid beta density
    if (scf_type == 2 .or. scf_type == 3) then 
        dtmp(:,:) = 0.0_dp 
        do i = 1, nbf 
            if (occ_b(i) > 1.0e-14_dp) then 
                do mu = 1, nbf 
                    do nu = 1, nbf 
                        dtmp(mu,nu) = dtmp(mu,nu) + occ_b(i)*mo_b(mu,i)*mo_b(nu,i)
                    end do 
                end do 
            end if 
        end do
        call pack_matrix(dtmp, pdmat(:,2))
    end if 
 end subroutine build_pfon_density 

!> @brief      This routine reorders orbitals to maximum overlap.
 subroutine reordermos(v,e,smo,l0,nbf,lr1,lr2)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(nbf,*) :: V
    real(kind=dp), intent(in), dimension(*) :: E
    real(kind=dp), intent(in), dimension(l0,*) :: Smo
    integer :: l0, nbf, lr1, lr2

    integer, allocatable, dimension(:) :: iwrk, iwrk2
    real(kind=dp), allocatable, dimension(:) :: wrk
    integer :: i, j, k, ip1
    real(kind=dp) :: ss, smax

    allocate(iwrk(l0), iwrk2(l0), source=0)
    allocate(wrk(l0), source=0.0_dp)

    do i = 1, l0
      smax = 0.0_dp
      iwrk(i) = 0
      do j = 1, l0
        do k = 1, i
          if(iwrk(k)==j) cycle
        end do
        ss = abs(smo(i,j))
        if(ss > smax) then
          smax = ss
          iwrk(i) = j
        end if
      end do
      if(smo(i, iwrk(i)) < 0.0_dp) then
        do j = 1, nbf
          v(j,iwrk(i)) = -v(j,iwrk(i))
        end do
      end if
    end do

    iwrk2 = iwrk

    do i = lr1, lr2
       j = iwrk(i)
       call dswap(nbf,v(1,i),1,v(1,j),1)
       ip1 = i+1
       do k = ip1, lr2
          if(iwrk(k)==i) iwrk(k) = j
       end do
    end do

    iwrk = iwrk2
    do i = lr1, lr2
       j = iwrk(i)
       call dswap(1,e(i),1,e(j),1)
       ip1 = i+1
       do k = ip1, lr2
          if(iwrk(k)==i) iwrk(k) = j
       end do
    end do
    return
 end subroutine

end module scf

