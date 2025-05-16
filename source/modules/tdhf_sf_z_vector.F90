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
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform_sym, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix
    use oqp_linalg
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_sf_z_vector"

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
    logical :: roref = .false.

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
    real(kind=dp) :: alpha, error

    logical :: dft
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
    if (.not. (mol_mult == 3 .or. mol_mult == 4)) then
      call show_message( &
        'SF-TDDFT only supports mult=3 (triplet) or mult=4 (quartet) references', &
        with_abort)
    end if 

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

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

    call infos%dat%remove_records(tags_alloc)

    call infos%dat%reserve_data(OQP_WAO, TA_TYPE_REAL64, nbf_tri, comment=OQP_WAO_comment)
    call infos%dat%reserve_data(OQP_td_p, TA_TYPE_REAL64, nbf_tri*2, (/ nbf_tri, 2 /), comment=OQP_td_p_comment)
    call infos%dat%reserve_data(OQP_td_abxc, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
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
    call tagarray_get_data(infos%dat, OQP_td_energies, sf_energies)

    ta          => td_t(:,1)
    tb          => td_t(:,2)

    ! Save unrelaxed density matrices and the `b=A*x` vector for target state
    call sfdmat(bvec_mo(:,infos%tddft%target_state), td_abxc, mo_a, ta, tb, nocca, noccb)

  ! Initialize ERI calculations
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

    int2_data = int2_tdgrd_data_t(d2=pa, &
            int_apb=.true., &
            int_amb=.false., &
            tamm_dancoff=.false., &
            scale_exchange=scale_exch)

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
    int2_data = int2_td_data_t(d2=bvec, &
            int_apb=.false., &
            int_amb=.false., &
            tamm_dancoff=.true., &
            scale_exchange=scale_exch2)

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

    call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)

    call pcgrbpini(errv, pk, error, rhs, xminv, lhs)

    write(*,'(" INITIAL ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') error, cnvtol

! -----------------------------------------------

    do iter = 1, infos%control%maxit_zv

      call sfrogen(wrk1, wrk2, pk, nocca, noccb)
!     Alpha
      call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
!     Beta
      call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

!     (A+B)*PK
      call int2_data%clean()
      deallocate(int2_data)
      int2_data = int2_tdgrd_data_t(d2=pa, &
              int_apb=.true., &
              int_amb=.false., &
              tamm_dancoff=.false., &
              scale_exchange=scale_exch)

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

      alpha = 1.0_dp/dot_product(pk, lhs)

      xk = xk + pk * alpha
      errv = errv - alpha*lhs

      error = dot_product(errv, errv)
      write(*,'(" ITER#",I2," ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') &
        iter, error, cnvtol
      call flush(iw)

      if (error<cnvtol) exit

      call pcgb(pk, errv, xminv)

    end do


! -----------------------------------------------
    if (error>cnvtol) then
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

    call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)

 !  Update density for alpha
    call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)

 !  Update density for beta
    call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_tdgrd_data_t(d2=pa, &
            int_apb=.true., int_amb=.false., tamm_dancoff=.false., &
            scale_exchange=scale_exch)

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

  end subroutine tdhf_sf_z_vector

end module tdhf_sf_z_vector_mod
