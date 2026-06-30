module tdhf_mrsf_energy_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_energy_mod"

contains

  subroutine tdhf_mrsf_energy_C(c_handle) bind(C, name="tdhf_mrsf_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    inf%tddft%umrsf = .false.
    call tdhf_mrsf_energy_with_restart(inf)
  end subroutine tdhf_mrsf_energy_C

  subroutine tdhf_umrsf_energy_C(c_handle) bind(C, name="tdhf_umrsf_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    logical :: previous_umrsf
    inf => oqp_handle_get_info(c_handle)
    previous_umrsf = inf%tddft%umrsf
    inf%tddft%umrsf = .true.
    call tdhf_mrsf_energy_with_restart(inf)
    inf%tddft%umrsf = previous_umrsf
  end subroutine tdhf_umrsf_energy_C

  ! Run the MRSF Davidson and, if it fails to converge, auto-restart with a
  ! larger subspace (maxvec) and more iterations (maxit_dav).  Re-invoking the
  ! driver reallocates a fresh, larger Krylov subspace, so no inner-loop state
  ! is reused.  The user's maxvec/maxit_dav are restored afterwards.
  subroutine tdhf_mrsf_energy_with_restart(infos)
    use types, only: information
    use io_constants, only: iw
    type(information), intent(inout) :: infos
    integer, parameter :: max_restarts = 2
    integer :: attempt, maxvec0, maxit0
    maxvec0 = infos%tddft%maxvec
    maxit0  = infos%control%maxit_dav
    do attempt = 0, max_restarts
      call tdhf_mrsf_energy(infos)
      if (infos%mol_energy%Davidson_converged) exit
      if (attempt < max_restarts) then
        infos%tddft%maxvec      = 2 * infos%tddft%maxvec
        infos%control%maxit_dav = 2 * infos%control%maxit_dav
        ! The energy routine closes the log on exit; reopen to record the restart.
        open(unit=iw, file=infos%log_filename, position="append")
        write(iw,'(/,2X,"MRSF Davidson not converged; auto-restart #",I0, &
                 &" with larger subspace (maxvec=",I0,", maxit_dav=",I0,")"/)') &
          attempt + 1, infos%tddft%maxvec, infos%control%maxit_dav
        close(iw)
      end if
    end do
    infos%tddft%maxvec      = maxvec0
    infos%control%maxit_dav = maxit0
  end subroutine tdhf_mrsf_energy_with_restart

  subroutine tdhf_mrsf_energy(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time

    use precision, only: dp
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_mrsf_data_t, int2_umrsf_data_t
    use tdhf_lib, only: sym_response_project, &
      int2_td_data_t
    use tdhf_lib, only: &
      iatogen, mntoia, rparedms, rpaeig, rpavnorm, &
      rpaechk, rpanewb, &
      rpaprint, inivec
    use tdhf_sf_lib, only: sfresvec, sfqvec, sfdmat, trfrmb, &
      get_transition_density, get_transitions, &
      get_transition_dipole, print_results, get_spin_square
    use tdhf_mrsf_lib, only: &
      mrinivec, mrsfcbc, umrsfcbc, mrsfmntoia, umrsfmntoia, mrsfesum, &
      mrsfqroesum, get_mrsf_transitions, &
      get_mrsf_transition_density, get_jacobi, umrsfssqu, mrsf_set_fp32
    use mathlib, only: orthogonal_transform, orthogonal_transform_sym, &
      unpack_matrix
    use oqp_linalg
    use int1, only: multipole_integrals
    use printing, only: print_module_info
    use iso_c_binding, only: c_f_pointer, c_int

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_energy"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: s_size, ok

    real(kind=dp), allocatable :: scr2(:),scr3(:)
    real(kind=dp), allocatable :: wrk1(:,:), qvec(:,:)
    real(kind=dp), allocatable :: sym_ritz(:,:)
    real(kind=dp), allocatable :: amo(:,:), wrk2(:,:)
    real(kind=dp), allocatable :: squared_S(:)
    real(kind=dp), allocatable :: amb(:,:), apb(:,:), smat_full(:,:)
    real(kind=dp), allocatable, target :: vl(:), vr(:)
    real(kind=dp), pointer :: vl_p(:,:), vr_p(:,:)
    real(kind=dp), allocatable :: xm(:), scr(:)
    real(kind=dp), allocatable :: bvec_mo(:,:), for_trnsf_b_vec(:,:)
    real(kind=dp), allocatable, dimension(:,:) :: fa, fb
    real(kind=dp), allocatable, dimension(:) :: rnorm
    real(kind=dp), allocatable, dimension(:) :: mo_energy_work_a, mo_energy_work_b
    real(kind=dp), allocatable, dimension(:,:,:,:) :: trden
    integer, allocatable, dimension(:,:) :: trans
    real(kind=dp), allocatable, target :: mrsf_density(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    real(kind=dp), allocatable, target :: fmrq1(:,:,:)
    real(kind=dp), allocatable :: dip(:,:,:), bvec_mo_tmp(:), eex(:)
    integer(c_int) , pointer :: ixcore_ptr(:)
    ! misc-excited-analysis: tagarray exposure of the MRSF densities / dipoles
    real(kind=dp), pointer :: trden_store(:,:,:), dip_store(:,:,:), dipao_store(:,:)
    real(kind=dp), allocatable :: mints_exp(:,:)
    real(kind=dp) :: com_exp(3)

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf2, xvec_dim
    integer :: mxvec, ist, jst, iend, nvec, novec
    integer :: iter, nv, iv, ivec
    integer :: diag_index, i
    integer :: mxiter
    logical :: tamm_dancoff
    integer :: imax
    integer :: ierr
    logical :: converged
    real(kind=dp) :: rc_save, rc_new
    real(kind=dp) :: mxerr, cnvtol, scale_exch
    real(kind=dp) :: spc_scale_coco, spc_scale_ovov, spc_scale_coov
    integer :: maxvec, mrst, nstates, target_state
    logical :: roref = .false.
    logical :: uhfref = .false.
    logical :: debug_mode

    type(int2_compute_t) :: int2_driver
    type(int2_mrsf_data_t), target :: int2_data_st
    type(int2_umrsf_data_t), target :: int2_udata_st

    type(int2_td_data_t), target :: int2_data_q

    logical :: dft = .false.
    integer :: scf_type, mol_mult

    logical :: umrsf

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), dmat_a(:), mo_A(:,:), mo_energy_a(:), &
      fock_b(:), dmat_b(:), mo_b(:,:), mo_energy_b(:), &
      smat(:), ta(:), tb(:), td_t(:,:), bvec_mo_out(:,:), &
      mrsf_energies(:)
    character(len=*), parameter :: tags_alloc(3) = (/ character(len=80) :: &
      OQP_td_bvec_mo, OQP_td_t, OQP_td_energies /)
    character(len=*), parameter :: tags_required(9) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B, OQP_SM /)

  ! Readings
  ! Files open
  ! 3. LOG: Write: Main output file
    open (unit=iw, file=infos%log_filename, position="append")

    umrsf = infos%tddft%umrsf
  !
    if (umrsf) then
      call print_module_info('UMRSF_TDHF_Energy','Computing Energy of UMRSF-TDDFT')
    else
      call print_module_info('MRSF_TDHF_Energy','Computing Energy of MRSF-TDDFT')
    end if

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

  ! Get Fortran pointer ixcore_ptr from C pointer
    if (.not. (infos%tddft%ixcore_len == 0)) &
    call c_f_pointer(infos%tddft%ixcore, ixcore_ptr, [infos%tddft%ixcore_len])

   ! Input parameters
    dft = infos%control%hamilton == 20 ! dft or hf
    mrst = infos%tddft%mult
    nstates = infos%tddft%nstate
    target_state = infos%tddft%target_state
    maxvec = infos%tddft%maxvec
    cnvtol = infos%tddft%cnvtol
    debug_mode = infos%tddft%debug_mode

    mol_mult = infos%mol_prop%mult
    if (umrsf) then
      if (mol_mult/=3) call show_message('UMRSF-TDDFT are available for UHF ref.&
          &with ONLY triplet multiplicity(mult=3)',with_abort)
    else
      if (mol_mult/=3) call show_message('MRSF-TDDFT are available for ROHF ref.&
          &with ONLY triplet multiplicity(mult=3)',with_abort)
    end if
    scf_type = infos%control%scftype
    if (.not. umrsf .and. scf_type==3) roref = .true.

    if (umrsf .and. scf_type/=2) then
      call show_message('UMRSF-TDDFT requires UHF reference (SCFTYPE=2).',with_abort)
    else if (umrsf) then
      uhfref = .true.
    end if

    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    s_size = (basis%nshell**2+basis%nshell)/2

    tamm_dancoff = .true.  ! tamm_dancoff: 0/1 means not doing/doing Tamm/Dancoff run

    nocca = infos%mol_prop%nelec_a
    nvira = nbf-nocca
    noccb = infos%mol_prop%nelec_b
    nvirb = nbf-noccb


    if (mrst==1 .or. mrst==3 ) then
      xvec_dim = nocca*nvirb
    else if (mrst==5) then
      xvec_dim = noccb*nvira
    end if


    if (mrst==1 ) then
      nstates = min(nstates, xvec_dim-1)
      mxvec = min(maxvec*nstates, xvec_dim-1, infos%control%maxit_dav*nstates)
    else if (mrst==3) then
      nstates = min(nstates, xvec_dim-3)
      mxvec = min(maxvec*nstates, xvec_dim-3, infos%control%maxit_dav*nstates)
    else if (mrst==5) then
      nstates = min(nstates, xvec_dim)
      mxvec = min(maxvec*nstates, xvec_dim, infos%control%maxit_dav*nstates)
    end if

    infos%tddft%nstate = nstates

    nvec = min(max(nstates,6), mxvec)

    call infos%dat%alloc_or_die(OQP_td_bvec_mo, (/xvec_dim, nstates/), bvec_mo_out, description=OQP_td_bvec_mo_comment)
    call infos%dat%alloc_or_die(OQP_td_t, (/ nbf2, 2 /), td_t, description=OQP_td_t_comment)
    call infos%dat%alloc_or_die(OQP_td_energies, (/ nstates /), mrsf_energies, description=OQP_td_energies_comment)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)


  ! Allocate temporary matrices for diagonalization
    allocate (fa(nbf,nbf), &
              fb(nbf,nbf), &
              stat=ok)
    if( ok/=0 ) call show_message('Cannot allocate memory',with_abort)

  ! Allocate TDDFT variables
    allocate(xm(xvec_dim),  &
             bvec_mo(xvec_dim,mxvec), &
             trden(nbf,nbf,nstates,nstates), &
             wrk1(nbf,nbf), &
             wrk2(nbf,nbf), &
             smat_full(nbf,nbf), &
             amo(xvec_dim,mxvec), &
             EEX(mxvec), &
             squared_S(nstates), &
             APB(mxvec,mxvec), &
             AMB(mxvec,mxvec), &
             VR(mxvec*mxvec), &
             VL(mxvec*mxvec), &
             for_trnsf_b_vec(mxvec,mxvec), & !
             dip(3,nstates,nstates), &
             scr(nbf2), &
             bvec_mo_tmp(xvec_dim), &
             scr2(mxvec*mxvec), &
             scr3(nocca), &
             qvec(xvec_dim,nstates), &
             RNORM(nstates), &
             source=0.0_dp,stat=ok)
    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)
    allocate(trans(xvec_dim,2), &
             source=0, stat=ok)
    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    mo_energy_work_a = mo_energy_a
    mo_energy_work_b = mo_energy_b
  ! MO rotations (Jacobi)
    if (umrsf) then
      call unpack_matrix(smat, smat_full, nbf, 'U')
      call get_jacobi(infos, mo_a, mo_energy_a, mo_b, mo_energy_b, smat_full, nocca, wrk1, wrk2, 0)
      call get_jacobi(infos, mo_a, mo_energy_a, mo_b, mo_energy_b, smat_full, nocca, wrk1, wrk2, 1)
    end if

    ta => td_t(:,1)
    tb => td_t(:,2)

    if (mrst==1 .or. mrst==3 ) then
      if (umrsf) then
        allocate(mrsf_density(nvec,11,nbf,nbf), &
                 source=0.0_dp, &
                 stat=ok)
      else
        allocate(mrsf_density(nvec,7,nbf,nbf), &
                 source=0.0_dp, &
                 stat=ok)
      end if
    else if( mrst==5  )then

      allocate(fmrq1(nbf,nbf,nvec), &
               source=0.0_dp, &
               stat=ok)

    end if

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    scale_exch = 1.0_dp
    if (infos%tddft%HFscale == -1.0_dp) &
          infos%tddft%HFscale = infos%dft%HFscale

    if (infos%dft%cam_flag) then
      if (infos%tddft%cam_alpha == -1.0_dp) &
            infos%tddft%cam_alpha = infos%dft%cam_alpha
      infos%tddft%HFscale = infos%tddft%cam_alpha
      if (infos%tddft%cam_beta == -1.0_dp) &
            infos%tddft%cam_beta = infos%dft%cam_beta
      if (infos%tddft%cam_mu == -1.0_dp) &
            infos%tddft%cam_mu = infos%dft%cam_mu
    end if
    if (dft) scale_exch = infos%tddft%HFscale
    ! Pure HF reference (no DFT functional): the effective exact-exchange scale
    ! is 1.0. Without a DFT functional infos%dft%HFscale is left at the -1.0
    ! sentinel, so the response HFscale (and hence the spin-pair coupling below)
    ! would inherit -1.0. The energy tolerates this (the fmrst2 rescale is
    ! skipped because spc == HFscale either way), but the MRSF gradient uses the
    ! spin-pair coupling values directly and needs the correct +1.0.
    if (.not. dft) infos%tddft%HFscale = 1.0_dp
    ! set spin-pair coupling
    if (infos%tddft%spc_coco==-1.0_dp) &
          infos%tddft%spc_coco = infos%tddft%HFscale
    if (infos%tddft%spc_ovov==-1.0_dp) &
          infos%tddft%spc_ovov = infos%tddft%HFscale
    if (infos%tddft%spc_coov==-1.0_dp) &
          infos%tddft%spc_coov = infos%tddft%HFscale

    if(debug_mode)then
      write(*,'(/,5x,"Input parameters:")')
      write(*,'(5x,"Number of states:                 ",1x,I0)') nstates
      write(*,'(5x,"Number of single excitations:     ",1x,I0)') xvec_dim
      write(*,'(5x,"Number of atomic orbitals:        ",1x,I0)') nbf
      write(*,'(5x,"Number of electrons:              ",1x,I0)') nocca+noccb
      write(*,'(5x,"Number of occupied alpha orbitals:",1x,I0)') nocca
      write(*,'(5x,"Number of occupied beta orbitals: ",1x,I0)') noccb
      write(*,'(5x,"Number of virtual alpha orbitals: ",1x,I0)') nvira
      write(*,'(5x,"Number of virtual beta orbitals:  ",1x,I0)') nvirb
      write(*,'(5x,"Maximum vectors:                  ",1x,I0)') mxvec
      write(*,'(5x,"Initial vectors:                  ",1x,I0)') nvec
      if (.not. (infos%tddft%ixcore_len == 0)) &
        write(*,'(5x,"Ixcore (MO index):                ",1x,I0)') ixcore_ptr
      write(*, '(/7x,"Fitting parameters for MRSF-TDDFT")')
      if (.not.infos%dft%cam_flag) then
        write(*, '(10x,"Exact HF exchange:")')
        write(*, '(5x,"Reference: |", t20, f6.3, t29, "|")') infos%dft%HFscale
        write(*, '(5x,"Response:  |", t20, f6.3, t29, "|")') infos%tddft%HFscale
      else
        write(*, '(10x,"CAM parametres:")')
        write(*, '(16x,"|   alpha   |    beta   |     mu    |")')
        write(*, '(5x,"Reference: |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
           infos%dft%cam_alpha, infos%dft%cam_beta, infos%dft%cam_mu
        write(*, '(5x,"Response:  |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
           infos%tddft%cam_alpha, infos%tddft%cam_beta, infos%tddft%cam_mu
      end if
      write(*, '(10x,"Spin-pair coupling parametres:")')
      write(*, '(16x,"|   CO-CO   |   OV-OV   |   CO-OV   |")')
      write(*, '(16x,"|", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
         infos%tddft%spc_coco, infos%tddft%spc_ovov, infos%tddft%spc_coov
    end if

    write(*,'(/,5x,46("="))')
    if (mrst==1) write(*,'(  5X,"Davidson algorithm for Singlet response states")')
    if (mrst==3) write(*,'(  5x,"Davidson algorithm for Triplet response states")')
    if (mrst==5) write(*,'(  5x,"Davidson algorithm for Quintet response states")')
    write(*,'(5x,46("="))')

    ! Loosen the 2e integral cutoff for the MRSF RESPONSE build only. The
    ! response is built on the converged orbitals, so it tolerates a far looser
    ! cutoff than the SCF default (5e-11). DEFAULT 1e-8 -- measured exact to the
    ! printed precision (<<1 ueV, far below the ~5e-5 regression tolerance and
    ! the iterative conv tolerance) -- removes integrals the response cannot
    ! resolve, cutting the integral COUNT (eval + digestion) for a modest free
    ! speedup that grows with system size. Override via env OQP_MRSF_RESP_CUTOFF
    ! (a.u.): set looser (e.g. 1e-7/1e-6) for more speed at ueV cost, or set to
    ! the SCF cutoff (5e-11) to recover the previous exact-tight behavior.
    ! max(SCF cutoff, requested) never goes tighter than the SCF integrals.
    ! Restored after the response so SCF / later steps are unaffected.
    ! Response 2e cutoff from [tdhf] resp_cutoff (infos%control%mrsf_resp_cutoff,
    ! default 1e-8). max(SCF cutoff, requested) never goes tighter than SCF.
    rc_save = infos%control%int2e_cutoff
    rc_new = infos%control%mrsf_resp_cutoff
    if (rc_new <= 0.0_dp) rc_new = 1.0e-8_dp
    infos%control%int2e_cutoff = max(rc_save, rc_new)

    ! FP32 response digestion from [tdhf] fp32 (infos%control%mrsf_fp32). Also
    ! reaches the z-vector gradient, which reuses the same process-global flag.
    call mrsf_set_fp32(int(infos%control%mrsf_fp32))

    ! Initialize ERI (Electron Repulsion Integrals) calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    call flush(iw)

  ! Prepare for ROHF
    if ((roref .and. .not. umrsf) .or. (uhfref .and. umrsf)) then
  !   Alpha
      call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, scr)

      ! shift Fock in MO basis here except MOs listed in ixcores
      if (.not. (infos%tddft%ixcore_len == 0)) then
        Do iter = 1, noccb
            if (.not. any(ixcore_ptr(1:infos%tddft%ixcore_len) == iter)) then
                diag_index = (iter + 1) * iter / 2
                scr(diag_index) = -1.0d6
            end if
        End Do
      end if

      call unpack_matrix(scr,fa)

  !   Beta
      call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, scr)
      call unpack_matrix(scr,fb)
    end if

    if (umrsf) then
      do i = 1, nbf
        mo_energy_work_a(i) = fa(i,i)
        mo_energy_work_b(i) = fb(i,i)
      end do
    end if


  ! Construct TD trial vector
    if (mrst==1 .or. mrst==3) then
      if (.not. umrsf) then
        call mrinivec(infos, mo_energy_work_a, mo_energy_work_a, bvec_mo, xm, nvec)
      else
        call mrinivec(infos, mo_energy_work_a, mo_energy_work_b, bvec_mo, xm, nvec)
      end if

    else if (mrst==5) then
      call inivec(mo_energy_a,mo_energy_a,bvec_mo,xm,noccb,nocca,nvec)
    end if

    ist = 1
    iend = nvec
    iter = 0
    mxiter = infos%control%maxit_dav
    ierr = 0

    do iter = 1, mxiter
      nv = iend-ist+1

      if( mrst==1 .or. mrst==3 ) then

        mrsf_density = 0.0_dp   ! bo2v, bo1v, bco1, bco2, o21v, co12, ball

      else if( mrst==5 ) then

        fmrq1 = 0.0_dp

      end if

      do ivec = ist, iend

        iv = ivec-ist+1

        if (mrst==1 .or. mrst==3) then

          call iatogen(bvec_mo(:,ivec), wrk1, nocca, noccb)
          if (umrsf) then
            call umrsfcbc(infos, mo_a, mo_b, wrk1, mrsf_density(iv,:,:,:))
          else
            call mrsfcbc(infos, mo_a, mo_b, wrk1, mrsf_density(iv,:,:,:))
          end if

        else if (mrst==5) then

          call iatogen(bvec_mo(:,ivec), wrk1, noccb, nocca)
          call orthogonal_transform('t', nbf, mo_a, wrk1, fmrq1(:,:,iv), wrk2)

        end if

      end do

      if (mrst==1 .or. mrst==3) then

        if (umrsf) then
          int2_udata_st = int2_umrsf_data_t( &
            d3 = mrsf_density(:iv,:,:,:), &
            tamm_dancoff = tamm_dancoff, &
            scale_exchange = scale_exch, &
            scale_coulomb = scale_exch)

          call int2_driver%run( &
            int2_udata_st, &
            cam = dft.and.infos%dft%cam_flag, &
            alpha = infos%tddft%cam_alpha, &
            alpha_coulomb = infos%tddft%cam_alpha, &
            beta = infos%tddft%cam_beta, &
            beta_coulomb = infos%tddft%cam_beta, &
            mu = infos%tddft%cam_mu)

          fmrst2 => int2_udata_st%f3(:,:,:,:,1) ! ado2v, ado1v, adco1, adco2, ao21v, aco12, agdlr

        else
          int2_data_st = int2_mrsf_data_t( &
            d3 = mrsf_density(:iv,:,:,:), &
            tamm_dancoff = tamm_dancoff, &
            scale_exchange = scale_exch, &
            scale_coulomb = scale_exch)

        call int2_driver%run( &
          int2_data_st, &
          cam = dft.and.infos%dft%cam_flag, &
          alpha = infos%tddft%cam_alpha, &
          alpha_coulomb = infos%tddft%cam_alpha, &
          beta = infos%tddft%cam_beta, &
          beta_coulomb = infos%tddft%cam_beta, &
          mu = infos%tddft%cam_mu)

        fmrst2 => int2_data_st%f3(:,:,:,:,1) ! ado2v, ado1v, adco1, adco2, ao21v, aco12, agdlr

        endif

        ! Scaling factor if triplet
        if (umrsf .and. mrst==3) then
          fmrst2(:,1:10,:,:) = -fmrst2(:,1:10,:,:)
        else if (mrst==3) then
          fmrst2(:,1:6,:,:) = -fmrst2(:,1:6,:,:)
        endif

        ! Spin pair coupling
        if (umrsf) then
          if (abs(infos%tddft%hfscale) > epsilon(1.0_dp)) then
            if (infos%tddft%spc_coco /= infos%tddft%hfscale) then
              spc_scale_coco = infos%tddft%spc_coco / infos%tddft%hfscale
              fmrst2(:,10,:,:) = fmrst2(:,10,:,:) * spc_scale_coco
            end if
            if (infos%tddft%spc_ovov /= infos%tddft%hfscale) then
              spc_scale_ovov = infos%tddft%spc_ovov / infos%tddft%hfscale
              fmrst2(:,9,:,:) = fmrst2(:,9,:,:) * spc_scale_ovov
            end if
            if (infos%tddft%spc_coov /= infos%tddft%hfscale) then
              spc_scale_coov = infos%tddft%spc_coov / infos%tddft%hfscale
              fmrst2(:,1:8,:,:) = fmrst2(:,1:8,:,:) * spc_scale_coov
            end if
          else if (infos%tddft%spc_coco /= 0.0_dp .or. &
                   infos%tddft%spc_ovov /= 0.0_dp .or. &
                   infos%tddft%spc_coov /= 0.0_dp) then
            call show_message('UMRSF-TDDFT spin-pair coupling overrides require nonzero HFscale.', with_abort)
          end if
        else
          if (infos%tddft%spc_coco /= infos%tddft%hfscale) &
             fmrst2(:,6,:,:) = fmrst2(:,6,:,:) * infos%tddft%spc_coco / infos%tddft%hfscale
          if (infos%tddft%spc_ovov /= infos%tddft%hfscale) &
             fmrst2(:,5,:,:) = fmrst2(:,5,:,:) * infos%tddft%spc_ovov / infos%tddft%hfscale
          if (infos%tddft%spc_coov /= infos%tddft%hfscale) &
             fmrst2(:,1:4,:,:) = fmrst2(:,1:4,:,:) * infos%tddft%spc_coov / infos%tddft%hfscale
        endif

      else if (mrst==5) then

        int2_data_q = int2_td_data_t( &
          d2=fmrq1(:,:,:iv), &
          int_apb = .false., &
          int_amb = .false., &
          tamm_dancoff = tamm_dancoff, &
          scale_exchange = scale_exch)
        call int2_driver%run( &
          int2_data_q, &
          cam = dft.and.infos%dft%cam_flag, &
          alpha = infos%tddft%cam_alpha, &
          beta = infos%tddft%cam_beta,&
          mu = infos%tddft%cam_mu)

      end if

      do ivec = ist, iend

        iv = ivec-ist+1

        if (mrst==1 .or. mrst==3) then

          ! Product (A-B)*X
          if (umrsf) then
            call umrsfmntoia(infos, fmrst2(iv,:,:,:), amo, mo_a, mo_b, ivec)
          else
            call mrsfmntoia(infos, fmrst2(iv,:,:,:), amo, mo_a, mo_b, ivec)
          end if

          call iatogen(bvec_mo(:,ivec), wrk1, nocca, noccb)

          call mrsfesum(infos, wrk1, fa, fb, amo, ivec)

        else if( mrst==5 )then

          call mntoia(int2_data_q%amb(:,:,iv,1), amo(:,ivec), mo_a, mo_b, noccb, nocca)

          ! Z(I+,A-)
          call iatogen(bvec_mo(:,ivec),wrk1,noccb,nocca)

          ! FB(I+,J+)*Z(J+,A-)
          call dgemm('n','n',noccb,nbf,noccb, &
                     1.0_dp,fb,nbf, &
                            wrk1,nbf, &
                     0.0_dp,wrk2,noccb)

          ! Z(I+,B-)*FA(B-,A-)
          call dgemm('n','n',noccb,nbf,nbf, &
                     1.0_dp,wrk1,nbf, &
                            fa,nbf, &
                    -1.0_dp,wrk2,noccb)

          call mrsfqroesum(wrk2,amo, &
                           nocca,noccb,nbf,ivec)

        end if
      end do

      vl_p(1:nvec, 1:nvec) => vl(1:nvec*nvec)
      vr_p(1:nvec, 1:nvec) => vr(1:nvec*nvec)
      call rparedms(bvec_mo,amo,amo,apb,amb,nvec,tamm_dancoff=.true.)
      call rpaeig(eex,vl_p,vr_p,apb,amb,scr2,tamm_dancoff=.true.)
      call rpavnorm(vr_p,vl_p,tamm_dancoff=.true.)
      call rpaechk(eex,nvec,nstates,imax,tamm_dancoff=.true.)

      for_trnsf_b_vec = vr_p
      call sfresvec(qvec,bvec_mo,amo,vr_p,eex,nvec,rnorm,nstates)
      call sfqvec(qvec,xm,eex,nstates)

!     Response-space symmetry blocking (no-op unless staged by pyoqp):
!     confine each root's update to the dominant irrep of its Ritz vector.
      sym_ritz = matmul(bvec_mo(:,1:nvec), vr_p(1:nvec,1:nstates))
      call sym_response_project(infos, sym_ritz, qvec, nstates)
      call rpaprint(eex, rnorm, cnvtol, iter, imax, nstates, do_neg=.true.)

      mxerr = maxval(rnorm)

!     Check convergence
      converged = mxerr<=cnvtol
      if (converged) exit

!     No space left for new vectors, exit
      if (nvec==mxvec) ierr = 1
      if (ierr/=0) exit

      call rpanewb(nstates,bvec_mo,qvec,novec,nvec,ierr,tamm_dancoff=.true.)

  !   ierr=1 nvec over mxvec: not converged case
      if (ierr/=0) exit

      ist = novec+1
      iend = nvec

    end do

    if (iter >= mxiter .and. .not. converged) ierr = -1

    select case (ierr)
    case (-1)
      write(*,'(/,2X,"MRSF-TD-DFT energies NOT CONVERGED after ",I4," iterations"/)') mxiter
      infos%mol_energy%Davidson_converged=.false.
    case (0)
      write(*,'(/,2X,"MRSF-TD-DFT energies converged in ",I4," iterations"/)') iter
      infos%mol_energy%Davidson_converged=.true.
    case (1)
      write(*,'(/,2X,"..something is wrong.. nvec = mxvec")')
      infos%mol_energy%Davidson_converged=.false.
    case (2)
      write(*,'(/,2x,"..something is wrong..  nvec > mxvec")')
      write(*,'(3x,"nvec/mxvec =",I4,"/",I4)') nvec, mxvec
      infos%mol_energy%Davidson_converged=.false.
    case (3)
      write(*,'(/,2x,"..something is wrong.. No vectors were added")')
      infos%mol_energy%Davidson_converged=.false.
    end select
    call flush(iw)

    call trfrmb(bvec_mo, for_trnsf_b_vec, nvec, nstates)

    select case (mrst)
      case(1)
        if (umrsf) then
          trden = 0.0_dp
        else
          do ist = 1, nstates
            do jst = ist, nstates
              call get_mrsf_transition_density(infos,trden(:,:,ist,jst), bvec_mo, ist, jst)
            end do
          end do
        end if

        if (umrsf) then
          do ist = 1, nstates
            call umrsfssqu(squared_S(ist), mo_a, mo_b, smat, wrk1, scr3, nbf, nbf2, &
                           xvec_dim, ist, nbf, bvec_mo, nocca, noccb, &
                           .true., .false.)
          end do
        else
          squared_S(:) = 0.0_dp
        end if
        call get_mrsf_transitions(trans, nocca, noccb, nbf)
        write(*,'(/,2x,35("="),/,2x,&
            &"Spin-adapted spin-flip excitations",/,2x,35("="))')
      case(3)
        if (umrsf) then
          trden = 0.0_dp
        else
          do ist = 1, nstates
            do jst = ist, nstates
              call get_mrsf_transition_density(infos, trden(:,:,ist,jst), bvec_mo, ist, jst)
            end do
          end do
        end if
        if (umrsf) then
          do ist = 1, nstates
            call umrsfssqu(squared_S(ist), mo_a, mo_b, smat, wrk1, scr3, nbf, nbf2, &
                           xvec_dim, ist, nbf, bvec_mo, nocca, noccb, &
                           .false., .true.)
          end do
        else
          squared_S(:) = 2.0_dp
        end if
        call get_mrsf_transitions(trans, nocca, noccb, nbf)
        write(*,'(/,2x,35("="),/,2x,&
            &"Spin-adapted spin-flip excitations",/,2x,35("="))')
      case(5)
        call get_transition_density(trden, bvec_mo, nbf, noccb, nocca, nstates)
        squared_S(:) = 6.0_dp
        call get_transitions(trans, noccb, nocca, nbf)
        write(*,'(/,2x,35("="),/,2x,&
            &"Beta -> Alpha spin-flip excitations",/,2x,35("="))')
      case default
        error stop "Unknown mrst value"
    end select

    call get_transition_dipole(basis, dip, mo_a, trden, nstates)

    ! --- misc-excited-analysis: expose the MRSF state-interaction transition /
    !     state-difference densities (alpha-MO basis), the transition dipoles,
    !     and the AO electric-dipole integrals for downstream Python analysis.
    !     Pure write-out; no physics above is altered (ported to the alloc_or_die
    !     tagarray API of current main).
    !     Skipped for UMRSF: there trden is set identically to zero above, so no
    !     genuine state-interaction densities exist and exposing them would
    !     publish misleading all-zero tags.
    if (.not. umrsf) then
      ! get_mrsf_transition_density / get_transition_dipole only populate the
      ! upper triangle (ist<=jst). Mirror it into the stored copies so reverse
      ! state pairs are correct: gamma^{j->i} = (gamma^{i->j})^T and the (real)
      ! transition dipole mu^{j->i} = mu^{i->j}. The live `dip`/`trden` arrays
      ! handed to print_results are left untouched.
      do jst = 1, nstates
        do ist = jst+1, nstates
          trden(:,:,ist,jst) = transpose(trden(:,:,jst,ist))
        end do
      end do

      call infos%dat%alloc_or_die(OQP_td_trans_density_mo, &
        (/ nbf, nbf, nstates*nstates /), trden_store, &
        description=OQP_td_trans_density_mo_comment)
      trden_store = reshape(trden(:,:,1:nstates,1:nstates), (/ nbf, nbf, nstates*nstates /))

      call infos%dat%alloc_or_die(OQP_td_trans_dipole, (/ 3, nstates, nstates /), &
        dip_store, description=OQP_td_trans_dipole_comment)
      dip_store = dip(:,1:nstates,1:nstates)
      do jst = 1, nstates
        do ist = jst+1, nstates
          dip_store(:,ist,jst) = dip_store(:,jst,ist)
        end do
      end do

      allocate(mints_exp(nbf2,3), source=0.0_dp)
      com_exp = basis%atoms%center(weight='mass')
      call multipole_integrals(basis, mints_exp, com_exp, 1)
      call infos%dat%alloc_or_die(OQP_td_dip_ao, (/ nbf2, 3 /), dipao_store, &
        description=OQP_td_dip_ao_comment)
      dipao_store = mints_exp
      deallocate(mints_exp)
    end if

    mrsf_energies = eex(1:nstates)
    bvec_mo_out = bvec_mo(:,1:nstates)
    infos%mol_energy%excited_energy = mrsf_energies(infos%tddft%target_state)
    call print_results(infos, bvec_mo, eex, trans, dip, squared_S, nstates)
    call flush(iw)

    call int2_driver%clean()
    infos%control%int2e_cutoff = rc_save

    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine tdhf_mrsf_energy


end module tdhf_mrsf_energy_mod
