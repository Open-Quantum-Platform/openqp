module tdhf_sf_z_vector_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_sf_z_vector_mod"

contains

  subroutine tdhf_sf_z_vector_C(c_handle) bind(C, name="tdhf_sf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_sf_z_vector(inf)
  end subroutine tdhf_sf_z_vector_C


  subroutine tdhf_sf_z_vector(infos)

    use precision, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
    use tdhf_lib, only: iatogen, mntoia
    use tdhf_sf_lib, only: sfrorhs, &
      sfromcal, sfrogen, sfrolhs, pcgrbpini, &
      pcgb, sfropcal, sfrowcal, sfdmat
    use mod_dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform_sym, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix
    use oqp_linalg
    use printing, only: print_module_info
    use zvector_common, only: sanitize_zvector_preconditioner, &
      zv_opts_t, zv_read_opts, zv_prog_tau

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_sf_z_vector"
    real(kind=dp), parameter :: SF_ZVEC_DENOMINATOR_FLOOR = 1.0d-12
    type(zv_opts_t) :: zvo
    real(kind=dp) :: zv_rc_tight

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

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf_tri
    integer :: iter
    real(kind=dp) :: cnvtol, scale_exch, scale_exch2
    logical :: roref

    type(int2_compute_t) :: int2_driver
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
    integer :: nsocc, lzdim

  ! General data
    real(kind=dp) :: alpha, error, pap

    logical :: dft, zvector_breakdown
    integer :: scf_type, mol_mult

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), mo_a(:,:), mo_energy_a(:), td_abxc(:,:), &
      fock_b(:), mo_b(:,:), &
      wao(:), td_p(:,:), td_t(:,:), &
      ta(:), tb(:), bvec_mo(:,:), sf_energies(:)
    character(len=*), parameter :: tags_alloc(3) = (/ character(len=80) :: &
      OQP_WAO, OQP_td_p, OQP_td_abxc /)
    character(len=*), parameter :: tags_required(8) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_E_MO_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_VEC_MO_B, OQP_td_bvec_mo, OQP_td_t, &
      OQP_td_energies /)

    mol_mult = infos%mol_prop%mult
 !   if (.not. (mol_mult == 3 .or. mol_mult == 4)) then
 !     call show_message( &
 !       'SF-TDDFT only supports mult=3 (triplet) or mult=4 (quartet) references', &
 !       with_abort)
 !   end if 

    scf_type = infos%control%scftype
    roref = scf_type == 3

    ! A UHF reference has independent alpha and beta orbital spaces: its orbital response is
    ! the open-shell occ-vir problem of each spin, not the ROHF doc/socc/virt partition.
    if (scf_type == 2) then
      call tdhf_sf_z_vector_uhf(infos)
      return
    end if

    dft = infos%control%hamilton == 20

  ! Files open
  ! 3. LOG: Write: Main output file
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('SF_TDHF_Z_Vector','Solving Z-Vector for SF-TDDFT')
  ! Readings

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    if (dft) call dft_initialize(infos, basis, molGrid)

  ! Parameter it should be inputed later
  ! convergence tolerance in the iterative TD-DFT step.
    cnvtol = infos%tddft%zvconv
    ! Shared z-vector perf opt-ins (env OQP_SF_ZV_*); progressive screening
    ! default ON, zvconv override default off (see zvector_common).
    call zv_read_opts(zvo, "SF")
    if (zvo%conv_user > 0.0_dp) cnvtol = zvo%conv_user

    nocca = infos%mol_prop%nelec_A
    nvira = nbf-noccA
    noccb = infos%mol_prop%nelec_B
    nvirb = nbf-noccB
    nsocc = nocca-noccb
    lzdim = noccb*(nsocc+nvira)+nsocc*nvira

    allocate(&
  ! for Z-vector
      xminv(lzdim), &
      rhs(lzdim), &
      lhs(lzdim), &
      xm(lzdim), &
      xk(lzdim), &
      pk(lzdim), &
      errv(lzdim), &
  ! for gradient
      hxa(nbf,nocca), &
      hxb(nbf,nbf), &
      tij(nocca,nocca), &
      tab(nvirb,nvirb), &
      ppija(nocca,nocca), &
      ppijb(noccb,noccb), &
      pa(nbf,nbf,2), &
   ! Allocate TDDFT variables
      fa(nbf,nbf), &           ! Temporary matrix for diagonalization
      fb(nbf,nbf), &           ! Temporary matrix for diagonalization
      ab1_MO_a(nocca,nvirb), &
      ab1_MO_b(noccb,nvirb), &
!   For scratch
      wrk1(nbf,nbf), &
      wrk2(nbf,nbf), &
      wrk3(nbf,nbf), &
      stat=ok, &
      source=0.0_dp)

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    call infos%dat%alloc_or_die(OQP_WAO, (/ nbf_tri /), wao, description=OQP_WAO_comment)
    call infos%dat%alloc_or_die(OQP_td_p, (/ nbf_tri, 2 /), td_p, description=OQP_td_p_comment)
    call infos%dat%alloc_or_die(OQP_td_abxc, (/ nbf, nbf /), td_abxc, description=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, sf_energies)

    ! The energy stage clips nstate to the response-space size, so the target
    ! must lie within the vectors actually solved.
    if (infos%tddft%target_state > size(bvec_mo,2)) then
      write(*,'(2x,a,i0,a,i0,a)') 'Requested gradient state ', infos%tddft%target_state, &
        ' exceeds the ', size(bvec_mo,2), ' states solved in this response space.'
      call show_message('SF-TDDFT gradient target state lies outside the solved response space; '// &
                        'request a lower state or a larger basis.', with_abort)
    end if

    ta          => td_t(:,1)
    tb          => td_t(:,2)

    ! Save unrelaxed density matrices and the `b=A*x` vector for target state
    call sfdmat(bvec_mo(:,infos%tddft%target_state), td_abxc, mo_a, ta, tb, nocca, noccb)

  ! Initialize ERI calculations
    ! Progressive screening keeps init at the tight cutoff (full pair list) and
    ! ramps the run-time threshold per CG iteration; restore tight for the tail.
    zv_rc_tight = infos%control%int2e_cutoff
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    write(*,'(/1x,71("-")&
             &/19x,"SF-DFT ENERGY GRADIENT CALCULATION"&
             &/1x,71("-")/)')
    write(iw,fmt='(5x,a/&
                  &5x,16("-")/&
                  &5x,a,x,i0,x,f17.10,x,"Hartree"/&
                  &5x,a,x,e10.4/&
                  &5x,a,x,i0)') &
        'Z-vector options' &
      , 'Target state       is', infos%tddft%target_state, infos%mol_energy%energy+sf_energies(infos%tddft%target_state) &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv
    call flush(iw)

    bvec(1:nbf,1:nbf,1:1) => td_abxc

  ! Prepare for ROHF
    ! Fock matrices A and B
    if( roref )then
        wrk1t(1:nbf*nbf) => wrk1
  !   Alapha
      call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, wrk1)
      call unpack_matrix(wrk1t, fa)

  !   Beta
      call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, wrk1)
      call unpack_matrix(wrk1t, fb)
    end if

  ! Make density like part
    call unpack_matrix(ta, pa(:,:,1))
    call unpack_matrix(tb, pa(:,:,2))

  ! Initialize ERI calculations
    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
       scale_exch = infos%dft%HFscale    !> Reference HF exchange
       scale_exch2 = infos%tddft%HFscale !> Response HF exchange
    end if

    if (allocated(int2_data)) deallocate(int2_data)
    allocate(int2_data, source=int2_tdgrd_data_t(d2=pa, &
            int_apb=.true., &
            int_amb=.false., &
            tamm_dancoff=.false., &
            scale_exchange=scale_exch))

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    pa = pa*2
    call utddft_fxc(basis=basis, &
           molGrid=molGrid, &
           isVecs=.true., &
           wfa=MO_A, &
           wfb=MO_B, &
           fxa=ab1(:,:,1:1), &
           fxb=ab1(:,:,2:2), &
           dxa=pa(:,:,1:1), &
           dxb=pa(:,:,2:2), &
           nmtx=1, &
           !threshold=1.0d-15, &
           threshold=0.0d0, &
           infos=infos)

!   ALPHA: AO(M,N) -> MO(IA+)
    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

    call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

  ! Initialize ERI calculations
    call int2_data%clean()
    deallocate(int2_data)
    allocate(int2_data, source=int2_td_data_t(d2=bvec, &
            int_apb=.false., &
            int_amb=.false., &
            tamm_dancoff=.true., &
            scale_exchange=scale_exch2))

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%tddft%cam_alpha, &
            beta=infos%tddft%cam_beta,&
            mu=infos%tddft%cam_mu)
    ab2 => int2_data%amb(:,:,:,1)

    call orthogonal_transform('n', nbf, mo_a, ab2(:,:,1), wrk2, wrk1)

    call iatogen(bvec_mo(:,infos%tddft%target_state), wrk3, nocca, noccb)

    call dgemm('n', 't', nbf, nocca, nbf,  &
               2.0_dp, wrk2, nbf,  &
                       wrk3, nbf,  &
               0.0_dp, hxa,  nbf)
    call dgemm('t', 'n', nbf, nbf, nocca,  &
               2.0_dp, wrk2, nbf,  &
                       wrk3, nbf,  &
               0.0_dp, hxb,  nbf)

!   Unrelaxed difference density matries T_ij and T_ab
!     Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
    call dgemm('n', 't', nocca, nocca, nvirb,  &
              -1.0_dp, bvec_mo(:,infos%tddft%target_state), nocca,  &
                       bvec_mo(:,infos%tddft%target_state), nocca,  &
               0.0_dp, tij,     nocca)

    ! Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
    call dgemm('t', 'n', nvirb, nvirb, nocca,  &
               1.0_dp, bvec_mo(:,infos%tddft%target_state), nocca,  &
                       bvec_mo(:,infos%tddft%target_state), nocca,  &
               0.0_dp, tab,     nvirb)

    call sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                 Tij, Tab, Fa, Fb, nocca, noccb)

    write(*,'(/3x,25("-")&
             &/6x,"START Z-VECTOR LOOP"&
             &/3x,25("-")/)')
    call flush(iw)

    call run_sf_cg_zvector()
    if (zvo%prog_on) call int2_driver%set_cutoff(zv_rc_tight)


! -----------------------------------------------
    if (zvector_breakdown) then
       infos%mol_energy%Z_Vector_converged=.false.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector breakdown"&
             &/3x,24("-")/)')
    else if (error>cnvtol) then
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

    if (zvector_breakdown) then
      call int2_driver%clean()
      if (dft) call dftclean(infos)
      call measure_time(print_total=1, log_unit=iw)
      close(iw)
      return
    end if

    call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)

 !  Update density for alpha
    call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)

 !  Update density for beta
    call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

    call int2_data%clean()
    deallocate(int2_data)
    allocate(int2_data, source=int2_tdgrd_data_t(d2=pa, &
            int_apb=.true., int_amb=.false., tamm_dancoff=.false., &
            scale_exchange=scale_exch))

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

    call utddft_fxc(basis=basis, &
           molGrid=molGrid, &
           isVecs=.true., &
           wfa=MO_A, &
           wfb=MO_B, &
           fxa=ab1(:,:,1:1), &
           fxb=ab1(:,:,2:2), &
           dxa=pa(:,:,1:1), &
           dxb=pa(:,:,2:2), &
           nmtx=1, &
           !threshold=1.0d-15, &
           threshold=0.0d0, &
           infos=infos)

!   ALPHA AO(M,N) -> MO(I-,J-) ... LPPIJA
    call dgemm('n', 'n', nbf, nocca, nbf,  &
               1.0_dp, ab1(:,:,1), nbf,  &
                       mo_a, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', nocca, nocca, nbf,  &
               1.0_dp, mo_a,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppija, nocca)
!   BETA: AO(M,N) -> MO(I-,J-) ... LPPIJB
    call dgemm('n', 'n', nbf, noccb, nbf,  &
               1.0_dp, ab1(:,:,2), nbf,  &
                       mo_a, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', noccb, noccb, nbf,  &
               1.0_dp, mo_a,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppijb, noccb)

!   Calculate W (in MO basis)
    wmo => wrk3
    wmo = 0
    call sfrowcal(wmo,sf_energies(infos%tddft%target_state), &
                  mo_energy_a, fa, fb, bvec_mo(:,infos%tddft%target_state), xk, &
                  hxa, hxb, ppija, ppijb, &
                  nocca, noccb)

    call orthogonal_transform('t', nbf, mo_a, wmo, wrk2, wrk1)
    call symmetrize_matrix(wrk2, nbf)
    call pack_matrix(wrk2, wao)
    wao = wao*0.5_dp
!   ROHF, half one more time:
    wao = wao*0.5_dp

    call int2_driver%clean()

    if (dft) call dftclean(infos)

    call measure_time(print_total=1, log_unit=iw)
    close(iw)


  contains

    ! Preconditioned CG z-vector solve.  All state is reached by host
    ! association, so this is behaviorally identical to the inline version.
    subroutine run_sf_cg_zvector()
      call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)
      call sanitize_zvector_preconditioner(xm, xminv, iw, SF_ZVEC_DENOMINATOR_FLOOR, "SF")

      call pcgrbpini(errv, pk, error, rhs, xminv, lhs)
      zvector_breakdown = .false.
      if (.not. ieee_is_finite(error) .or. any(.not. ieee_is_finite(errv)) .or. &
          any(.not. ieee_is_finite(pk)) .or. any(.not. ieee_is_finite(lhs))) then
        zvector_breakdown = .true.
        write(*,'(/3x,24("-")&
              &/6x,"Z-Vector breakdown: non-finite initial PCG state"&
              &/3x,24("-")/)')
      end if

      if (infos%control%verbose >= 1) &
        write(*,'(" INITIAL ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') error, cnvtol

  ! -----------------------------------------------

      do iter = 1, infos%control%maxit_zv

        if (zvector_breakdown) exit

        call sfrogen(wrk1, wrk2, pk, nocca, noccb)
  !     Alpha
        call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
  !     Beta
        call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

  !     Progressive screening: loosen cutoff while the residual is large, pinned
  !     tight near convergence (zv_prog_tau); restored after the loop.
        if (zvo%prog_on) call int2_driver%set_cutoff(zv_prog_tau(zvo, error, zv_rc_tight))

  !     (A+B)*PK
        call int2_data%clean()
        deallocate(int2_data)
        allocate(int2_data, source=int2_tdgrd_data_t(d2=pa, &
                int_apb=.true., &
                int_amb=.false., &
                tamm_dancoff=.false., &
                scale_exchange=scale_exch))

        call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
        ab1 => int2_data%apb(:,:,:,1)

        !ab1 = ab1/2
        call symmetrize_matrix(pa(:,:,1), nbf)
        call symmetrize_matrix(pa(:,:,2), nbf)
        call utddft_fxc(basis=basis, &
               molGrid=molGrid, &
               isVecs=.true., &
               wfa=MO_A, &
               wfb=MO_B, &
               fxa=ab1(:,:,1:1), &
               fxb=ab1(:,:,2:2), &
               dxa=pa(:,:,1:1), &
               dxb=pa(:,:,2:2), &
               nmtx=1, &
               !threshold=1.0d-15, &
               threshold=0.0d0, &
               infos=infos)

  !     ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
        call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

        call mntoia(ab1(:,:,2), ab1_mo_b, mo_a, mo_a, noccb, noccb)

        call sfrolhs(lhs, pk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                     nocca, noccb)

        if (any(.not. ieee_is_finite(lhs)) .or. any(.not. ieee_is_finite(pk))) then
          zvector_breakdown = .true.
          write(*,'(" Z-Vector breakdown: non-finite SF PCG operator state at iter", I4)') iter
          exit
        end if

        pap = dot_product(pk, lhs)
        if (.not. ieee_is_finite(pap) .or. abs(pap) < SF_ZVEC_DENOMINATOR_FLOOR) then
          zvector_breakdown = .true.
          write(*,'(" Z-Vector breakdown: unsafe SF PCG denominator at iter", I4, 1x, 1p,e12.4)') iter, pap
          exit
        end if

        alpha = 1.0_dp / pap
        if (.not. ieee_is_finite(alpha)) then
          zvector_breakdown = .true.
          write(*,'(" Z-Vector breakdown: non-finite SF PCG alpha at iter", I4)') iter
          exit
        end if

        xk = xk + pk * alpha
        errv = errv - alpha*lhs
        if (any(.not. ieee_is_finite(xk)) .or. any(.not. ieee_is_finite(errv))) then
          zvector_breakdown = .true.
          write(*,'(" Z-Vector breakdown: non-finite SF PCG update at iter", I4)') iter
          exit
        end if

        error = dot_product(errv, errv)
        if (.not. ieee_is_finite(error)) then
          zvector_breakdown = .true.
          write(*,'(" Z-Vector breakdown: non-finite SF PCG residual at iter", I4)') iter
          exit
        end if
        if (infos%control%verbose >= 1) &
          write(*,'(" ITER#",I2," ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') &
          iter, error, cnvtol
        call flush(iw)

        if (error<cnvtol) exit

        call pcgb(pk, errv, xminv)
        if (any(.not. ieee_is_finite(pk))) then
          zvector_breakdown = .true.
          write(*,'(" Z-Vector breakdown: non-finite SF PCG search direction at iter", I4)') iter
          exit
        end if

      end do
    end subroutine run_sf_cg_zvector
  end subroutine tdhf_sf_z_vector

!###############################################################################
!> @brief Z-vector, relaxed difference density and energy-weighted density of SF-TDDFT
!>        (TDA) with a UHF/UKS triplet reference.
!>
!> The response energy of target state I with amplitudes X (alpha-occupied i, beta-virtual a) is
!>   omega = Tr(T^a F^a) + Tr(T^b F^b) - c_x sum_{ia,jb} X_ia X_jb (ij|ab),
!>   T^a_ij = -sum_a X_ia X_ja (alpha occ-occ),   T^b_ab = sum_i X_ia X_ib (beta virt-virt),
!> with the alpha holes in the alpha orbitals and the beta particles in the beta orbitals.
!> With Q^s_pq = sum_mu C^s_mu,p d(omega)/dC^s_mu,q, the orbital gradient of an occ-vir rotation
!> of spin s is R^s_ai = Q^s_ai - Q^s_ia + 2 G^s[T]_ai, where G^s[D] is the reference Fock
!> response (Coulomb of both spins, scaled same-spin exchange, f_xc).  The Z-vector solves the
!> open-shell orbital-Hessian equation M z = -R (cphf_solve_uhf).  The gradient then contracts
!> the relaxed density D^s = T^s + P^z_s/2, P^z_s = sum_ai z_ai (C_a C_i^T + C_i C_a^T), and the
!> energy-weighted density W^s = (Qt^s + Qt^s^T)/4 of the full Lagrangian (Qt includes the
!> z-explicit Fock terms and 2 G^s[D] on the occupied columns).
  subroutine tdhf_sf_z_vector_uhf(infos)

    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t, int2_tdgrd_data_t, iatogen, mntoia
    use tdhf_sf_lib, only: sfdmat
    use mod_dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: orthogonal_transform_sym, pack_matrix, unpack_matrix
    use oqp_linalg
    use printing, only: print_module_info
    use cphf_mod, only: cphf_solve_uhf

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_sf_z_vector_uhf"
    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    class(int2_td_data_t), allocatable, target :: int2_data
    type(dft_grid_t) :: molGrid

    integer :: nbf, nbf_tri, nocca, noccb, nvira, nvirb, la, lb, tgt, i, j, a, b, ij
    logical :: dft, zconv
    real(kind=dp) :: scale_exch, scale_exch2

    real(kind=dp), allocatable, target :: wrk1(:,:)
    real(kind=dp), pointer :: wrk1t(:)
    real(kind=dp), allocatable :: xmat(:,:), fa(:,:), fb(:,:), wrk2(:,:), wrk3(:,:), &
      tij(:,:), tab(:,:), hxa(:,:), hxb(:,:), qa(:,:), qb(:,:), hpta(:,:), hptb(:,:), &
      rhs(:,:), zvec(:,:), za(:,:), zb(:,:), wtot(:,:)
    real(kind=dp), allocatable, target :: pa(:,:,:)
    real(kind=dp), pointer :: ab1(:,:,:), ab2(:,:,:), bvec(:,:,:)

    real(kind=dp), contiguous, pointer :: &
      fock_a(:), fock_b(:), mo_a(:,:), mo_b(:,:), td_abxc(:,:), td_p(:,:), td_t(:,:), &
      wao(:), bvec_mo(:,:), sf_energies(:)
    character(len=*), parameter :: tags_required(8) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_VEC_MO_B, OQP_E_MO_A, OQP_td_bvec_mo, &
      OQP_td_t, OQP_td_energies /)

    dft = infos%control%hamilton == 20

    open (unit=iw, file=infos%log_filename, position="append")
    call print_module_info('SF_TDHF_Z_Vector','Solving Z-Vector for SF-TDDFT (UHF reference)')

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    la = nocca*nvira
    lb = noccb*nvirb

    call infos%dat%alloc_or_die(OQP_WAO, (/ nbf_tri /), wao, description=OQP_WAO_comment)
    call infos%dat%alloc_or_die(OQP_td_p, (/ nbf_tri, 2 /), td_p, description=OQP_td_p_comment)
    call infos%dat%alloc_or_die(OQP_td_abxc, (/ nbf, nbf /), td_abxc, description=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, sf_energies)

    ! cphf_solve_uhf builds the full-range open-shell orbital Hessian only
    if (dft) then
      if (infos%dft%cam_flag) call show_message('The UHF SF-TDDFT gradient does not support '// &
        'range-separated (CAM/LRC) exchange; use an ROHF reference.', with_abort)
    end if

    tgt = infos%tddft%target_state
    if (tgt < 1 .or. tgt > size(bvec_mo,2)) then
      write(*,'(2x,a,i0,a,i0,a)') 'Requested gradient state ', tgt, ' outside the ', &
        size(bvec_mo,2), ' states solved in this response space.'
      call show_message('SF-TDDFT gradient target state lies outside the solved response space; '// &
                        'request a lower state or a larger basis.', with_abort)
    end if

    write(iw,fmt='(5x,a/5x,16("-")/5x,a,x,i0,x,f17.10,x,"Hartree"/5x,a,x,e10.4/5x,a,x,i0)') &
        'Z-vector options' &
      , 'Target state       is', tgt, infos%mol_energy%energy+sf_energies(tgt) &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv
    call flush(iw)

    allocate(xmat(nbf,nbf), fa(nbf,nbf), fb(nbf,nbf), wrk1(nbf,nbf), wrk2(nbf,nbf), wrk3(nbf,nbf), &
             tij(nocca,nocca), tab(nvirb,nvirb), hxa(nbf,nocca), hxb(nbf,nbf), &
             qa(nbf,nbf), qb(nbf,nbf), hpta(nocca,nvira), hptb(noccb,nvirb), &
             rhs(la+lb,1), zvec(la+lb,1), za(nbf,nbf), zb(nbf,nbf), wtot(nbf,nbf), &
             pa(nbf,nbf,2), source=0.0_dp)

    ! Transition density D = C_a X C_b^T and unrelaxed difference densities (AO, packed in td_t)
    call sfdmat(bvec_mo(:,tgt), td_abxc, mo_a, td_t(:,1), td_t(:,2), nocca, noccb, mo_b=mo_b)
    call iatogen(bvec_mo(:,tgt), xmat, nocca, noccb)
    call dgemm('n','t', nocca,nocca,nvirb, -1.0_dp,bvec_mo(:,tgt),nocca, bvec_mo(:,tgt),nocca, &
               0.0_dp,tij,nocca)
    call dgemm('t','n', nvirb,nvirb,nocca,  1.0_dp,bvec_mo(:,tgt),nocca, bvec_mo(:,tgt),nocca, &
               0.0_dp,tab,nvirb)

    ! Reference Fock matrices in each spin's own MO basis
    wrk1t(1:nbf*nbf) => wrk1
    call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, wrk1)
    call unpack_matrix(wrk1t, fa)
    call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, wrk1)
    call unpack_matrix(wrk1t, fb)

    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
      scale_exch = infos%dft%HFscale     ! reference exact exchange
      scale_exch2 = infos%tddft%HFscale  ! response (spin-flip) exact exchange
    end if

    if (dft) call dft_initialize(infos, basis, molGrid)
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    ! 2 G[T]: reference Fock response to the unrelaxed difference densities
    call unpack_matrix(td_t(:,1), pa(:,:,1))
    call unpack_matrix(td_t(:,2), pa(:,:,2))
    call sf_uhf_fock_response(pa)
    call mntoia(ab1(:,:,1), hpta, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), hptb, mo_b, mo_b, noccb, noccb)

    ! Explicit orbital dependence of the spin-flip exchange term: M = -c_x K[D] (AO)
    call int2_data%clean()
    deallocate(int2_data)
    bvec(1:nbf,1:nbf,1:1) => td_abxc
    allocate(int2_data, source=int2_td_data_t(d2=bvec, int_apb=.false., int_amb=.false., &
             tamm_dancoff=.true., scale_exchange=scale_exch2))
    call int2_driver%run(int2_data, cam=dft.and.infos%dft%cam_flag, alpha=infos%tddft%cam_alpha, &
                         beta=infos%tddft%cam_beta, mu=infos%tddft%cam_mu)
    ab2 => int2_data%amb(:,:,:,1)
    ! wrk2 = C_a^T M C_b : rows alpha orbitals (holes), columns beta orbitals (particles)
    call dgemm('t','n', nbf,nbf,nbf, 1.0_dp,mo_a,nbf, ab2(:,:,1),nbf, 0.0_dp,wrk1,nbf)
    call dgemm('n','n', nbf,nbf,nbf, 1.0_dp,wrk1,nbf, mo_b,nbf, 0.0_dp,wrk2,nbf)
    call dgemm('n','t', nbf,nocca,nbf, 2.0_dp,wrk2,nbf, xmat,nbf, 0.0_dp,hxa,nbf)
    call dgemm('t','n', nbf,nbf,nocca, 2.0_dp,wrk2,nbf, xmat,nbf, 0.0_dp,hxb,nbf)

    ! Q of the explicit orbital dependence of omega (alpha: occupied columns; beta: virtual)
    qa(:,1:nocca) = hxa
    call dgemm('n','n', nbf,nocca,nocca, 2.0_dp,fa,nbf, tij,nocca, 1.0_dp,qa,nbf)
    qb(:,noccb+1:nbf) = hxb(:,noccb+1:nbf)
    call dgemm('n','n', nbf,nvirb,nvirb, 2.0_dp,fb(:,noccb+1:),nbf, tab,nvirb, 1.0_dp, &
               qb(:,noccb+1:),nbf)

    ! Right-hand side -R in the cphf_solve_uhf layout: alpha occ-vir then beta occ-vir,
    ! each occupied-index fastest
    ij = 0
    do a = nocca+1, nbf
      do i = 1, nocca
        ij = ij + 1
        rhs(ij,1) = -(hpta(i,a-nocca) + qa(a,i) - qa(i,a))
      end do
    end do
    do b = noccb+1, nbf
      do j = 1, noccb
        ij = ij + 1
        rhs(ij,1) = -(hptb(j,b-noccb) + qb(b,j) - qb(j,b))
      end do
    end do

    write(*,'(/3x,25("-")/6x,"START Z-VECTOR LOOP (UHF)"/3x,25("-")/)')
    call flush(iw)
    ! cphf_solve_uhf builds and releases its own DFT grid state
    if (dft) call dftclean(infos)
    call cphf_solve_uhf(infos, 1, rhs, zvec, tol=infos%tddft%zvconv, &
                        maxit=int(infos%control%maxit_zv), converged=zconv)
    if (dft) call dft_initialize(infos, basis, molGrid)
    infos%mol_energy%Z_Vector_converged = zconv
    if (zconv) then
      write(*,'(/3x,24("-")/6x,"Z-Vector converged"/3x,24("-")/)')
    else
      write(*,'(/3x,24("-")/6x,"Z-Vector not converged"/3x,24("-")/)')
    end if
    call flush(iw)

    ! z in nbf x nbf occ-vir blocks
    ij = 0
    do a = nocca+1, nbf
      do i = 1, nocca
        ij = ij + 1
        za(i,a) = zvec(ij,1)
      end do
    end do
    do b = noccb+1, nbf
      do j = 1, noccb
        ij = ij + 1
        zb(j,b) = zvec(ij,1)
      end do
    end do

    ! Relaxed difference densities D^s = T^s + P^z_s/2 (AO), with P^z_s = C (z + z^T) C^T
    wrk1 = za + transpose(za)
    call dgemm('n','n', nbf,nbf,nbf, 1.0_dp,mo_a,nbf, wrk1,nbf, 0.0_dp,wrk2,nbf)
    call dgemm('n','t', nbf,nbf,nbf, 0.5_dp,wrk2,nbf, mo_a,nbf, 0.0_dp,wrk3,nbf)
    call unpack_matrix(td_t(:,1), pa(:,:,1))
    pa(:,:,1) = pa(:,:,1) + wrk3
    wrk1 = zb + transpose(zb)
    call dgemm('n','n', nbf,nbf,nbf, 1.0_dp,mo_b,nbf, wrk1,nbf, 0.0_dp,wrk2,nbf)
    call dgemm('n','t', nbf,nbf,nbf, 0.5_dp,wrk2,nbf, mo_b,nbf, 0.0_dp,wrk3,nbf)
    call unpack_matrix(td_t(:,2), pa(:,:,2))
    pa(:,:,2) = pa(:,:,2) + wrk3
    call pack_matrix(pa(:,:,1), td_p(:,1))
    call pack_matrix(pa(:,:,2), td_p(:,2))

    ! Lagrangian Q: z-explicit Fock terms and 2 G[D] on the occupied columns
    call dgemm('n','n', nbf,nvira,nocca, 1.0_dp,fa,nbf, za(1:nocca,nocca+1:),nocca, 1.0_dp, &
               qa(:,nocca+1:),nbf)
    call dgemm('n','t', nbf,nocca,nvira, 1.0_dp,fa(:,nocca+1:),nbf, za(1:nocca,nocca+1:),nocca, &
               1.0_dp,qa,nbf)
    call dgemm('n','n', nbf,nvirb,noccb, 1.0_dp,fb,nbf, zb(1:noccb,noccb+1:),noccb, 1.0_dp, &
               qb(:,noccb+1:),nbf)
    call dgemm('n','t', nbf,noccb,nvirb, 1.0_dp,fb(:,noccb+1:),nbf, zb(1:noccb,noccb+1:),noccb, &
               1.0_dp,qb,nbf)

    call sf_uhf_fock_response(pa)
    call dgemm('t','n', nbf,nbf,nbf, 1.0_dp,mo_a,nbf, ab1(:,:,1),nbf, 0.0_dp,wrk1,nbf)
    call dgemm('n','n', nbf,nocca,nbf, 1.0_dp,wrk1,nbf, mo_a,nbf, 1.0_dp,qa,nbf)
    call dgemm('t','n', nbf,nbf,nbf, 1.0_dp,mo_b,nbf, ab1(:,:,2),nbf, 0.0_dp,wrk1,nbf)
    call dgemm('n','n', nbf,noccb,nbf, 1.0_dp,wrk1,nbf, mo_b,nbf, 1.0_dp,qb,nbf)

    ! W^s = (Q + Q^T)/4 in each spin's MO basis, summed in the AO basis
    wrk1 = 0.25_dp*(qa + transpose(qa))
    call dgemm('n','n', nbf,nbf,nbf, 1.0_dp,mo_a,nbf, wrk1,nbf, 0.0_dp,wrk2,nbf)
    call dgemm('n','t', nbf,nbf,nbf, 1.0_dp,wrk2,nbf, mo_a,nbf, 0.0_dp,wtot,nbf)
    wrk1 = 0.25_dp*(qb + transpose(qb))
    call dgemm('n','n', nbf,nbf,nbf, 1.0_dp,mo_b,nbf, wrk1,nbf, 0.0_dp,wrk2,nbf)
    call dgemm('n','t', nbf,nbf,nbf, 1.0_dp,wrk2,nbf, mo_b,nbf, 1.0_dp,wtot,nbf)
    ! The overlap gradient contracts eijden + 2*wao, where eijden is -W_ref packed with a halved
    ! diagonal; the response energy-weighted density follows the same convention.
    call pack_matrix(wtot, wao)
    ij = 0
    do i = 1, nbf
      ij = ij + i
      wao(ij) = 0.5_dp*wao(ij)
    end do
    wao = -0.5_dp*wao

    call int2_data%clean()
    deallocate(int2_data)
    call int2_driver%clean()
    if (dft) call dftclean(infos)
    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  contains

    !> ab1 = 2 G[d] (AO, per spin): reference Fock response to the spin densities d,
    !> Coulomb + scaled same-spin exchange from int2, f_xc on the grid.  d is doubled on exit.
    subroutine sf_uhf_fock_response(d)
      real(kind=dp), intent(inout), target :: d(:,:,:)
      if (allocated(int2_data)) then
        call int2_data%clean()
        deallocate(int2_data)
      end if
      allocate(int2_data, source=int2_tdgrd_data_t(d2=d, int_apb=.true., int_amb=.false., &
               tamm_dancoff=.false., scale_exchange=scale_exch))
      call int2_driver%run(int2_data, cam=dft.and.infos%dft%cam_flag, alpha=infos%dft%cam_alpha, &
                           beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)
      d = 2*d
      if (dft) call utddft_fxc(basis=basis, molGrid=molGrid, isVecs=.true., wfa=mo_a, wfb=mo_b, &
                               fxa=ab1(:,:,1:1), fxb=ab1(:,:,2:2), dxa=d(:,:,1:1), &
                               dxb=d(:,:,2:2), nmtx=1, threshold=0.0d0, infos=infos)
    end subroutine sf_uhf_fock_response

  end subroutine tdhf_sf_z_vector_uhf

end module tdhf_sf_z_vector_mod
