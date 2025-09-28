module tdhf_mrsf_z_vector_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_z_vector_mod"

contains

  subroutine tdhf_mrsf_z_vector_C(c_handle) bind(C, name="tdhf_mrsf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_z_vector(inf)
  end subroutine tdhf_mrsf_z_vector_C

  subroutine tdhf_mrsf_z_vector(infos)
    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time

    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use tdhf_mrsf_lib, only: int2_mrsf_data_t
    use tdhf_lib, only: iatogen, mntoia
    use tdhf_sf_lib, only: sfrorhs, &
      sfromcal, sfrogen, sfrolhs, pcgrbpini, &
      pcgb, sfropcal, sfdmat
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform, &
            orthogonal_transform_sym
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix

    use tdhf_mrsf_lib, only: &
      mrinivec, mrsfcbc, mrsfxvec, mrsfsp, mrsfrowcal, &
      mrsfqrorhs, mrsfqropcal, mrsfqrowcal
    use oqp_linalg
    use printing, only: print_module_info


    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_z_vector"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: ok

    real(kind=dp), allocatable :: ab1_mo_a(:,:)
    real(kind=dp), allocatable :: ab1_mo_b(:,:)
    real(kind=dp), allocatable :: xm(:)
    real(kind=dp), pointer :: ab2(:,:,:)
    real(kind=dp), pointer :: ab1(:,:,:)
    real(kind=dp), allocatable :: fa(:,:), fb(:,:)
    real(kind=dp), pointer :: bvec(:,:,:)
    real(kind=dp), pointer :: wmo(:,:)
    real(kind=dp), allocatable :: bvec_mo_d(:,:)
    real(kind=dp), allocatable, target :: &
      fmrst1(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf_tri
    integer :: iter
    real(kind=dp) :: cnvtol, scale_exch, scale_exch2
    logical :: roref = .false.
    integer :: mrst

    ! Solver mode parameter - must be declared early
    integer, parameter :: solver_mode = 3  ! 1=Original CG, 2=CG+DIIS, 3=DIIS-only

    type(int2_compute_t) :: int2_driver
    type(int2_mrsf_data_t), allocatable, target :: int2_data_st
    type(int2_td_data_t), allocatable, target :: int2_data_q
    class(int2_td_data_t), allocatable, target :: int2_data
    type(dft_grid_t) :: molGrid

  ! scr data
    real(kind=dp), allocatable, target :: wrk1(:,:), wrk2(:,:), wrk3(:,:)
    real(kind=dp), pointer :: wrk1t(:)

  ! SF-TD Gradient data
    real(kind=dp), allocatable :: &
      rhs(:), lhs(:), xminv(:), xk(:), pk(:), errv(:), &
      hxa(:,:), hxb(:,:), tij(:,:), ppija(:,:), ppijb(:,:), tab(:,:)
    real(kind=dp), allocatable, target :: pa(:,:,:)
    integer :: nsocc, lzdim, xvec_dim

  ! General data
    real(kind=dp) :: alpha, error

    ! DIIS variables (solver_mode defined above)
    integer, parameter :: max_diis = 6
    integer :: diis_dim, diis_start
    real(kind=dp), allocatable :: diis_x(:,:), diis_err(:,:)
    real(kind=dp), allocatable :: diis_b(:,:), diis_c(:)
    real(kind=dp), allocatable :: diis_xv(:), diis_rv(:)
    logical :: use_diis, diis_active
    real(kind=dp) :: error_old
    integer :: diis_info

    ! CG+DIIS control
    real(kind=dp) :: step_denom
    integer :: no_progress_count, cg_restart_count

    logical :: dft
    integer :: scf_type, mol_mult, target_state

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), mo_a(:,:), mo_energy_a(:), &
      fock_b(:), mo_b(:,:), &
      td_p(:,:), td_t(:,:), ta(:), tb(:), td_abxc(:,:), &
      td_mrsf_den(:,:,:), bvec_mo(:,:), wao(:), mrsf_energies(:)
    character(len=*), parameter :: tags_alloc(4) = (/ character(len=80) :: &
      OQP_WAO, OQP_td_mrsf_density, OQP_td_p, OQP_td_abxc /)
    character(len=*), parameter :: tags_required(8) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_E_MO_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_VEC_MO_B, OQP_td_bvec_mo, OQP_td_t, &
      OQP_td_energies /)

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message(&
            'MRSF-TDDFT are available for ROHF/UHF ref.&
            &with ONLY triplet multiplicity(mult=3)', with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    dft = infos%control%hamilton == 20

    open (unit=iw, file=infos%log_filename, position="append")
    call print_module_info('MRSF_TDHF_Z_Vector','CG+DIIS Z-Vector Solver')

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    if (dft) call dft_initialize(infos, basis, molGrid)

    mrst = infos%tddft%mult
    cnvtol = infos%tddft%zvconv

    nocca = infos%mol_prop%nelec_A
    nvira = nbf-noccA
    noccb = infos%mol_prop%nelec_B
    nvirb = nbf-noccb
    nsocc = nocca-noccb
    lzdim = noccb*(nsocc+nvira)+nsocc*nvira

    ! DIIS settings 
    use_diis = (solver_mode >= 2)  ! DIIS for modes 2 and 3
    diis_start = 3  
    diis_dim = 0
    diis_active = .false.

    if(mrst==1 .or. mrst==3) then
      xvec_dim = nocca*nvirb
      allocate(&
        fmrst1(1,7,nbf,nbf), &
        bvec_mo_d(xvec_dim,1), &
        hxa(nbf,nocca), &
        hxb(nbf,nbf), &
        tij(nocca,nocca), &
        tab(nvirb,nvirb), &
        stat=ok, &
        source=0.0_dp)
    else if(mrst==5) then
      xvec_dim = noccb*nvira
      allocate(&
        hxa(nbf,nbf), &
        hxb(nbf,noccb), &
        tij(noccb,noccb), &
        tab(nvira,nvira), &
        stat=ok, &
        source=0.0_dp)
    endif
    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    allocate(&
      xminv(lzdim), &
      rhs(lzdim), &
      lhs(lzdim), &
      xm(lzdim), &
      xk(lzdim), &
      pk(lzdim), &
      errv(lzdim), &
      pa(nbf,nbf,2), &
      ppija(nocca,nocca), &
      ppijb(noccb,noccb), &
      fa(nbf,nbf), &
      fb(nbf,nbf), &
      ab1_mo_a(nocca,nvira), &
      ab1_mo_b(noccb,nvirb), &
      wrk1(nbf,nbf), &
      wrk2(nbf,nbf), &
      wrk3(nbf,nbf), &
      stat=ok, &
      source=0.0_dp)
    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    ! DIIS arrays
    if (use_diis) then
      allocate(&
        diis_x(lzdim, max_diis), &
        diis_err(lzdim, max_diis), &
        diis_b(max_diis+1, max_diis+1), &
        diis_c(max_diis+1), &
        diis_xv(lzdim), &
        diis_rv(lzdim), &
        stat=ok, source=0.0_dp)
      if( ok/=0 ) call show_message('Cannot allocate DIIS memory', with_abort)
    endif

    call infos%dat%remove_records(tags_alloc)
    call infos%dat%reserve_data(OQP_WAO, TA_TYPE_REAL64, nbf_tri, comment=OQP_WAO_comment)
    call infos%dat%reserve_data(OQP_td_mrsf_density, TA_TYPE_REAL64, nbf*nbf*7, (/7, nbf, nbf /), comment=OQP_td_mrsf_density)
    call infos%dat%reserve_data(OQP_td_p, TA_TYPE_REAL64, nbf_tri*2, (/ nbf_tri, 2 /), comment=OQP_td_p)
    call infos%dat%reserve_data(OQP_td_abxc, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_den)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, mrsf_energies)

    ta => td_t(:,1)
    tb => td_t(:,2)

    target_state = min(infos%tddft%target_state, infos%tddft%nstate)
    if (target_state /=infos%tddft%target_state) then
      write(*,'(/1x,66("-")&
               &/1x,"WARNING: Target state has been changed to the max available nstates"/&
               &/1x,66("-")/)')
    end if

    ! [Keep all original setup for density matrices and RHS calculation]
    if (mrst==1 .or. mrst==3 ) then
      call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))
      call sfdmat(bvec_mo_d(:,1), td_abxc, mo_a, ta, tb, nocca, noccb)
    else if (mrst==5 ) then
      call sfdmat(bvec_mo(:,target_state), td_abxc, mo_a, tb, ta, noccb, nocca)
    end if

    bvec(1:nbf,1:nbf,1:1) => td_abxc

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    write(*,'(/1x,71("-")&
             &/19x,"MRSF-DFT ENERGY GRADIENT CALCULATION"&
             &/1x,71("-")/)')

    write(iw,fmt='(5x,a/&
                  &5x,16("-")/&
                  &5x,a,x,i0,x,f17.10,x,"Hartree"/&
                  &5x,a,x,i0/&
                  &5x,a,x,e10.4/&
                  &5x,a,x,i0/&
                  &5x,a,x,i0/&
                  &5x,a,x,l1)') &
        'Z-vector options' &
      , 'Target state       is', target_state, infos%mol_energy%energy+mrsf_energies(target_state) &
      , 'Multiplicity       is', infos%tddft%mult &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv &
      , 'Solver mode        is', solver_mode &
      , 'Using DIIS          is', use_diis
    call flush(iw)

    ! [Keep all original ROHF and setup code]
    if( roref )then
        wrk1t(1:nbf*nbf) => wrk1
        call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, wrk1)
        call unpack_matrix(wrk1t, fa)
        call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, wrk1)
        call unpack_matrix(wrk1t, fb)
    end if

    call unpack_matrix(ta, pa(:,:,1))
    call unpack_matrix(tb, pa(:,:,2))

    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
       scale_exch = infos%dft%HFscale
       scale_exch2 = infos%tddft%HFscale
    end if

    if (mrst==1 .or. mrst==3 ) then
      int2_data_st = int2_mrsf_data_t( &
          d3 = fmrst1, &
          tamm_dancoff = .true., &
          scale_exchange = scale_exch2, &
          scale_coulomb = scale_exch2)
    else if( mrst==5  )then
      int2_data_q = int2_td_data_t( &
          d2=bvec, &
          int_apb = .false., &
          int_amb = .false., &
          tamm_dancoff = .true., &
          scale_exchange = scale_exch2)
    end if

    int2_data = int2_tdgrd_data_t( &
        d2 = pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    pa = pa*2
    call utddft_fxc( &
        basis = basis, &
        molGrid = molGrid, &
        isVecs = .true., &
        wfa = mo_a, &
        wfb = mo_b, &
        fxa = ab1(:,:,1:1), &
        fxb = ab1(:,:,2:2), &
        dxa = pa(:,:,1:1), &
        dxb = pa(:,:,2:2), &
        nmtx = 1, &
        threshold = 1.0d-15, &
        infos = infos)

    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

    ! [Keep all original RHS setup]
    if (mrst==1 .or. mrst==3) then
      call iatogen(bvec_mo(:,target_state), wrk1, nocca, noccb)
      call mrsfcbc(infos, mo_a, mo_a, wrk1, fmrst1(1,:,:,:))
      fmrst1(1,7,:,:) = td_abxc
      td_mrsf_den(1:7,:,:) = fmrst1(1,1:7,:,:)

      call int2_driver%run(int2_data_st, &
            cam = dft.and.infos%dft%cam_flag, &
            alpha = infos%tddft%cam_alpha, &
            alpha_coulomb = infos%tddft%cam_alpha, &
            beta = infos%tddft%cam_beta,&
            beta_coulomb = infos%tddft%cam_beta, &
            mu = infos%tddft%cam_mu)
      fmrst2 => int2_data_st%f3(:,:,:,:,1)

      if (mrst==3) fmrst2(:,1:6,:,:) = -1.0_dp*fmrst2(:,1:6,:,:)

      if (infos%tddft%spc_coco /= infos%tddft%hfscale) &
         fmrst2(:,6,:,:) = fmrst2(:,6,:,:) * infos%tddft%spc_coco / infos%tddft%hfscale
      if (infos%tddft%spc_ovov /= infos%tddft%hfscale) &
         fmrst2(:,5,:,:) = fmrst2(:,5,:,:) * infos%tddft%spc_ovov / infos%tddft%hfscale
      if (infos%tddft%spc_coov /= infos%tddft%hfscale) &
         fmrst2(:,1:4,:,:) = fmrst2(:,1:4,:,:) * infos%tddft%spc_coov / infos%tddft%hfscale

      call orthogonal_transform('n', nbf, mo_a, fmrst2(1,7,:,:), wrk2, wrk1)
      call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))
      call iatogen(bvec_mo_d(:,1), wrk3, nocca, noccb)

      call dgemm('n', 't', nbf, nocca, nbf, &
                 2.0_dp, wrk2, nbf, wrk3, nbf, 0.0_dp, hxa, nbf)
      call dgemm('t', 'n', nbf, nbf, nocca, &
                 2.0_dp, wrk2, nbf, wrk3, nbf, 0.0_dp, hxb, nbf)

      call mrsfsp(hxa, hxb, mo_a, mo_a, wrk3, fmrst2(1,:,:,:), nocca, noccb)

      call dgemm('n', 't', nocca, nocca, nvirb, &
                -1.0_dp, bvec_mo_d, nocca, bvec_mo_d, nocca, 0.0_dp, tij, nocca)
      call dgemm('t', 'n', nvirb, nvirb, nocca, &
                 1.0_dp, bvec_mo_d, nocca, bvec_mo_d, nocca, 0.0_dp, tab, nvirb)

      call sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, Tij, Tab, Fa, Fb, nocca, noccb)

    else if(mrst==5) then
      call int2_driver%run(int2_data_q, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%tddft%cam_alpha, &
            beta=infos%tddft%cam_beta,&
            mu=infos%tddft%cam_mu)

      call orthogonal_transform('n', nbf, mo_a, int2_data_q%amb(:,:,1,1), wrk2, wrk1)
      call iatogen(bvec_mo(:,target_state),wrk3,noccb,nocca)

      call dgemm('t', 'n', nbf, nbf, noccb, &
                 2.0_dp, wrk2, nbf, wrk3, nbf, 0.0_dp, hxa, nbf)
      call dgemm('n', 't', nbf, noccb, nbf, &
                 2.0_dp, wrk2, nbf, wrk3, nbf, 0.0_dp, hxb, nbf)

      call dgemm('n', 't', noccb, noccb, nvira, &
                -1.0_dp, bvec_mo(:,target_state), noccb, &
                         bvec_mo(:,target_state), noccb, &
                 0.0_dp, tij, noccb)
      call dgemm('t', 'n', nvira, nvira, noccb, &
                 1.0_dp, bvec_mo(:,target_state), noccb, &
                         bvec_mo(:,target_state), noccb, &
                 0.0_dp, tab, nvira)

      call mrsfqrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, tab, tij, fa, fb, nocca, noccb)
    end if

    select case (solver_mode)
    case (1)
        write(*,'(/3x,25("-")&
                 &/6x,"ORIGINAL CG Z-VECTOR"&
                 &/3x,25("-")/)')
    case (2)
        write(*,'(/3x,25("-")&
                 &/6x,"CG+DIIS Z-VECTOR"&
                 &/3x,25("-")/)')
    case (3)
        write(*,'(/3x,25("-")&
                 &/6x,"DIIS-ONLY Z-VECTOR"&
                 &/3x,25("-")/)')
    case default
        write(*,'(/3x,25("-")&
                 &/6x,"DEFAULT CG Z-VECTOR"&
                 &/3x,25("-")/)')
    end select

    call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)

    call sfrogen(wrk1, wrk2, xk, nocca, noccb)
    call orthogonal_transform('n', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
    call orthogonal_transform('n', nbf, mo_a, wrk2, pa(:,:,2), wrk3)

    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_td_data_t( &
        d2 = pa, &
        int_apb = .true., &
        int_amb = .true., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)
    ab2 => int2_data%amb(:,:,:,1)

    pa = pa*2
    call utddft_fxc( &
        basis = basis, &
        molGrid = molGrid, &
        isVecs = .true., &
        wfa = mo_a, &
        wfb = mo_b, &
        fxa = ab1(:,:,1:1), &
        fxb = ab1(:,:,2:2), &
        dxa = pa(:,:,1:1), &
        dxb = pa(:,:,2:2), &
        nmtx = 1, &
        threshold = 1.0d-15, &
        infos = infos)

    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

    wrk1 = 2*ab1(:,:,2) + ab2(:,:,2)
    call mntoia(wrk1, ab1_mo_b, mo_a, mo_a, noccb, noccb)

    call sfrolhs(lhs, xk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, nocca, noccb)

    call pcgrbpini(errv, pk, error, rhs, xminv, lhs)

    no_progress_count = 0
    cg_restart_count = 0

    write(iw,'(" Initial error =",3x,1p,e10.3,1x,"/",1p,e10.3)') error, cnvtol
    call flush(iw)

! ===============================
! MAIN ITERATION LOOP
! ===============================

    do iter = 1, infos%control%maxit_zv

      select case (solver_mode)
      case (1)
          ! Original CG method
          call sfrogen(wrk1, wrk2, pk, nocca, noccb)
          call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
          call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

          call int2_data%clean()
          deallocate(int2_data)
          int2_data = int2_tdgrd_data_t( &
              d2 = pa, &
              int_apb = .true., &
              int_amb = .false., &
              tamm_dancoff = .false., &
              scale_exchange = scale_exch)

          call int2_driver%run(int2_data, &
                cam=dft.and.infos%dft%cam_flag, &
                alpha=infos%dft%cam_alpha, &
                beta=infos%dft%cam_beta,&
                mu=infos%dft%cam_mu)
          ab1 => int2_data%apb(:,:,:,1)

          call symmetrize_matrix(pa(:,:,1), nbf)
          call symmetrize_matrix(pa(:,:,2), nbf)
          call utddft_fxc( &
              basis = basis, &
              molGrid = molGrid, &
              isVecs = .true., &
              wfa = mo_a, &
              wfb = mo_b, &
              fxa = ab1(:,:,1:1), &
              fxb = ab1(:,:,2:2), &
              dxa = pa(:,:,1:1), &
              dxb = pa(:,:,2:2), &
              nmtx = 1, &
              threshold = 1.0d-15, &
              infos = infos)

          call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)
          call mntoia(ab1(:,:,2), ab1_mo_b, mo_a, mo_a, noccb, noccb)

          call sfrolhs(lhs, pk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, nocca, noccb)

          ! Safe CG step calculation
          step_denom = dot_product(pk, lhs)
          if (abs(step_denom) > 1.0e-12_dp) then
              alpha = dot_product(errv, errv) / step_denom
              ! Limit step size for stability
              if (abs(alpha) > 1.0_dp) alpha = sign(1.0_dp, alpha)
          else
              alpha = 0.0_dp  ! Skip step if denominator too small
          endif

          if (abs(alpha) > 1.0e-15_dp) then
              xk = xk + pk * alpha
              errv = errv - alpha * lhs
          endif
          
          error = sqrt(dot_product(errv, errv))
          write(iw,'(" Iter",i4,": error =",1p,e10.3," alpha =",1p,e10.3)') iter, error, alpha
          
          if (error < cnvtol) exit
          
          call pcgb(pk, errv, xminv)

      case (2)
          ! CG+DIIS hybrid
          if (iter >= diis_start) then
              ! Try DIIS first
              call diis_extrapolate(xk, errv, diis_x, diis_err, diis_b, diis_c, &
                                    diis_xv, diis_rv, diis_dim, max_diis, lzdim, &
                                    error_old, diis_active, diis_info)
              
              if (diis_active) then
                  xk = diis_xv
                  errv = diis_rv
                  error = error_old
                  write(iw,'("   --> DIIS extrapolation: ",1p,e10.3)') error
                  if (error < cnvtol) exit
              else
                  ! Fall back to CG step
                  call execute_cg_step()
              endif
          else
              ! Initial CG steps
              call execute_cg_step()
          endif

      case (3)
          ! DIIS-only mode
          if (iter >= diis_start) then
              call diis_extrapolate(xk, errv, diis_x, diis_err, diis_b, diis_c, &
                                    diis_xv, diis_rv, diis_dim, max_diis, lzdim, &
                                    error_old, diis_active, diis_info)
              
              if (diis_active) then
                  xk = diis_xv
                  errv = diis_rv
                  error = error_old
                  write(iw,'("   --> DIIS extrapolation: ",1p,e10.3)') error
                  if (error < cnvtol) exit
              else
                  ! DIIS failed, take simple mixing step
                  xk = xk - 0.3_dp * errv * xminv
                  error = sqrt(dot_product(errv, errv))
                  write(iw,'("   --> DIIS failed, simple mixing: ",1p,e10.3)') error
                  if (error < cnvtol) exit
              endif
          else
              ! Build DIIS history with conservative CG steps
              call execute_cg_step()
          endif

      case default
          ! Default to original CG
          call execute_cg_step()
      end select

      call flush(iw)

    end do

    if (error > cnvtol) then
       infos%mol_energy%Z_Vector_converged=.false.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector not converged"&
             &/3x,24("-")/)')
    else
       infos%mol_energy%Z_Vector_converged=.true.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector converged"&
             &/3x,24("-")/)')
    endif

    call flush(iw)

    ! [Rest of subroutine remains exactly the same as original]
    if (mrst==1 .or. mrst==3) then
      call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)
    else if (mrst==5) then
      call mrsfqropcal(wrk1, wrk2, tab, tij, xk, nocca, noccb)
    end if

    call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
    call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_tdgrd_data_t( &
        d2 = pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    call symmetrize_matrix(pa(:,:,1), nbf)
    call symmetrize_matrix(pa(:,:,2), nbf)
    call pack_matrix(pa(:,:,1), td_p(:,1))
    call pack_matrix(pa(:,:,2), td_p(:,2))

    td_p = 0.5_dp*td_p

    call utddft_fxc( &
        basis = basis, &
        molGrid = molGrid, &
        isVecs = .true., &
        wfa = mo_a, &
        wfb = mo_b, &
        fxa = ab1(:,:,1:1), &
        fxb = ab1(:,:,2:2), &
        dxa = pa(:,:,1:1), &
        dxb = pa(:,:,2:2), &
        nmtx = 1, &
        threshold = 1.0d-15, &
        infos = infos)

    call dgemm('n', 'n', nbf, nocca, nbf,  &
               1.0_dp, ab1(:,:,1), nbf, mo_a, nbf, 0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', nocca, nocca, nbf,  &
               1.0_dp, mo_a,  nbf, wrk2,  nbf, 0.0_dp, ppija, nocca)

    call dgemm('n', 'n', nbf, noccb, nbf,  &
               1.0_dp, ab1(:,:,2), nbf, mo_a, nbf, 0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', noccb, noccb, nbf,  &
               1.0_dp, mo_a,  nbf, wrk2,  nbf, 0.0_dp, ppijb, noccb)

    wmo => wrk3
    wmo = 0
    if (mrst==1 .or. mrst==3) then
      call mrsfrowcal(wmo, mo_energy_a, fa, fb, xk, hxa, hxb, ppija, ppijb, nocca, noccb)
    else if (mrst==5) then
      call mrsfqrowcal(wmo, mo_energy_a, fa, fb, xk, hxa, hxb, ppija, ppijb, nocca, noccb)
    end if

    call orthogonal_transform('t', nbf, mo_a, wmo, wrk2, wrk1)
    call symmetrize_matrix(wrk2, nbf)
    call pack_matrix(wrk2, wao)
    wao = wao*0.5_dp
    wao = wao*0.5_dp

    call int2_driver%clean()
    if (dft) call dftclean(infos)
    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  contains

    subroutine execute_cg_step()
      implicit none
      
      call sfrogen(wrk1, wrk2, pk, nocca, noccb)
      call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
      call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

      call int2_data%clean()
      deallocate(int2_data)
      int2_data = int2_tdgrd_data_t( &
          d2 = pa, &
          int_apb = .true., &
          int_amb = .false., &
          tamm_dancoff = .false., &
          scale_exchange = scale_exch)

      call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)

      call symmetrize_matrix(pa(:,:,1), nbf)
      call symmetrize_matrix(pa(:,:,2), nbf)
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = pa(:,:,1:1), &
          dxb = pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)

      call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)
      call mntoia(ab1(:,:,2), ab1_mo_b, mo_a, mo_a, noccb, noccb)

      call sfrolhs(lhs, pk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, nocca, noccb)

      ! Conservative CG step for initialization/fallback
      step_denom = dot_product(pk, lhs)
      if (abs(step_denom) > 1.0e-12_dp) then
          if (solver_mode == 1) then
              ! Original CG formula
              alpha = dot_product(errv, errv) / step_denom
              if (abs(alpha) > 1.0_dp) alpha = sign(1.0_dp, alpha)
          else
              ! Conservative step for building DIIS history
              alpha = dot_product(errv, errv) / step_denom
              alpha = sign(min(abs(alpha), 0.3_dp), alpha)
          endif
      else
          alpha = 0.0_dp
      endif

      if (abs(alpha) > 1.0e-15_dp) then
          xk = xk + pk * alpha
          errv = errv - alpha * lhs
      endif
      
      error = sqrt(dot_product(errv, errv))
      write(iw,'(" Iter",i4,": error =",1p,e10.3," alpha =",1p,e10.3)') iter, error, alpha
      
      if (error < cnvtol) return
      
      call pcgb(pk, errv, xminv)
    end subroutine execute_cg_step

    subroutine diis_extrapolate(x, err, diis_x, diis_err, B_matrix, coeff, &
                               x_new, err_new, diis_dim, max_diis, vec_dim, &
                               error_new, success, info)
      implicit none
      integer, intent(in) :: max_diis, vec_dim
      integer, intent(inout) :: diis_dim
      real(kind=dp), intent(in) :: x(:), err(:)
      real(kind=dp), intent(inout) :: diis_x(:,:), diis_err(:,:)
      real(kind=dp), intent(inout) :: B_matrix(:,:), coeff(:)
      real(kind=dp), intent(out) :: x_new(:), err_new(:), error_new
      logical, intent(out) :: success
      integer, intent(out) :: info
      
      integer :: i, j, pos
      integer, allocatable :: ipiv(:)
      
      success = .false.
      
      ! Update DIIS dimension
      if (diis_dim < max_diis) then
          diis_dim = diis_dim + 1
          pos = diis_dim
      else
          ! Shift arrays
          do i = 1, max_diis-1
              diis_x(:,i) = diis_x(:,i+1)
              diis_err(:,i) = diis_err(:,i+1)
          enddo
          pos = max_diis
      endif
      
      ! Store current vectors
      diis_x(:,pos) = x
      diis_err(:,pos) = err
      
      if (diis_dim < 2) return  ! Need at least 2 vectors
      
      ! Build B matrix
      do i = 1, diis_dim
          do j = 1, diis_dim
              B_matrix(i,j) = dot_product(diis_err(:,i), diis_err(:,j))
          enddo
          B_matrix(i,diis_dim+1) = -1.0_dp
          B_matrix(diis_dim+1,i) = -1.0_dp
      enddo
      B_matrix(diis_dim+1,diis_dim+1) = 0.0_dp
      
      ! Add damping
      do i = 1, diis_dim
          B_matrix(i,i) = B_matrix(i,i) + 1.0e-12_dp
      enddo
      
      coeff = 0.0_dp
      coeff(diis_dim+1) = -1.0_dp
      
      allocate(ipiv(diis_dim+1))
      call dgesv(diis_dim+1, 1, B_matrix, diis_dim+1, ipiv, coeff, diis_dim+1, info)
      deallocate(ipiv)
      
      if (info == 0) then
          x_new = 0.0_dp
          err_new = 0.0_dp
          do i = 1, diis_dim
              x_new = x_new + coeff(i) * diis_x(:,i)
              err_new = err_new + coeff(i) * diis_err(:,i)
          enddo
          
          error_new = sqrt(dot_product(err_new, err_new))
          success = .true.
      endif
      
    end subroutine diis_extrapolate

  end subroutine tdhf_mrsf_z_vector

end module tdhf_mrsf_z_vector_mod