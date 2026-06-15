module tdhf_mrsf_gradient_mod

  use precision, only: dp
  use grd2, only: grd2_driver, grd2_compute_data_t
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use printing, only: print_module_info

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_gradient_mod"

  public tdhf_mrsf_gradient

  type, extends(grd2_compute_data_t) :: grd2_mrsf_compute_data_t
    real(kind=dp), pointer :: d2(:,:,:) => null()
    real(kind=dp), pointer :: p2(:,:,:) => null()
    real(kind=dp), pointer :: spc2(:,:,:) => null()
    ! Cartesian-effective (bfnrm-folded) copies + offsets for HARMONIC_ACTIVE:
    ! d/p (alpha,beta) and the seven spin-pair-coupling densities.
    real(kind=dp), allocatable :: d2a_c(:,:), d2b_c(:,:), p2a_c(:,:), p2b_c(:,:)
    real(kind=dp), allocatable :: ball_c(:,:), bo2v_c(:,:), bo1v_c(:,:), bco1_c(:,:), &
                                  bco2_c(:,:), o21v_c(:,:), co12_c(:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
    integer :: mrst = 1
    real(kind=dp), dimension(3) :: spcscale = [0.0_dp, 0.0_dp, 0.0_dp]
  contains
    procedure :: init => grd2_mrsf_compute_data_t_init
    procedure :: clean => grd2_mrsf_compute_data_t_clean
    procedure :: get_density => grd2_mrsf_compute_data_t_get_density
    procedure :: build_cart => grd2_mrsf_build_cart
  end type

  ! NAC amplitude-term (Phase 11): the bilinear generalisation of
  ! grd2_mrsf_compute_data_t for G_IJ = X_I^T (d_x A) X_J. The reference
  ! density d2 (state-independent) and the interstate relaxed difference
  ! density p2 (=p2^{IJ}) are shared scalars; the seven transition/amplitude
  ! channels are carried separately for states I (spcI) and J (spcJ). Every
  ! pure transition-density product f(a)*g(b) in get_density becomes the
  ! symmetrised bilinear 1/2[ f_I(a) g_J(b) + f_J(a) g_I(b) ], while the
  ! reference/relaxed (df1,dq1) terms keep the production form with p2->p2^{IJ}.
  ! At I=J this collapses bit-for-bit to grd2_mrsf_compute_data_t.
  type, extends(grd2_compute_data_t) :: grd2_mrsf_nac_compute_data_t
    real(kind=dp), pointer :: d2(:,:,:) => null()
    real(kind=dp), pointer :: p2(:,:,:) => null()
    real(kind=dp), pointer :: spcI(:,:,:) => null()
    real(kind=dp), pointer :: spcJ(:,:,:) => null()
    integer :: nbf = 0
    integer :: mrst = 1
    real(kind=dp), dimension(3) :: spcscale = [0.0_dp, 0.0_dp, 0.0_dp]
  contains
    procedure :: init => grd2_mrsf_nac_compute_data_t_init
    procedure :: clean => grd2_mrsf_nac_compute_data_t_clean
    procedure :: get_density => grd2_mrsf_nac_compute_data_t_get_density
  end type

contains

  subroutine tdhf_mrsf_gradient_C(c_handle) bind(C, name="tdhf_mrsf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_gradient(inf)
  end subroutine tdhf_mrsf_gradient_c

  subroutine tdhf_mrsf_gradient(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort

    use grd1, only: eijden, print_gradient
    use util, only: measure_time
    use mod_dft, only: dft_initialize, dftclean
    use mathlib, only: symmetrize_matrix
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
    use mathlib, only: unpack_matrix
    use tdhf_sf_gradient_mod, only: sf_1e_grad, sf_2e_grad
    use tdhf_sf_lib, only: mrsf_state_label

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_gradient"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: s_size

    integer :: nbf, nbf2
    integer :: mrst
    character(len=12) :: target_label
    character(len=16) :: method_name
    logical :: roref = .false.

    type(dft_grid_t) :: molGrid

  ! General data
    logical :: dft
    integer :: scf_type, mol_mult

    real(kind=dp), allocatable :: p(:,:,:), v(:,:,:), d(:,:,:), spc(:,:,:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:), td_mrsf_density(:,:,:), td_abxc(:,:), td_p(:,:)
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_td_abxc, OQP_td_p /)
    character(len=*), parameter :: tags_mrsf(1) = (/ character(len=80) :: &
      OQP_td_mrsf_density /)

    dft = infos%control%hamilton == 20
    if (dft) then
      method_name = 'MRSF-TDDFT'
    else
      method_name = 'MRSF-TDHF'
    end if

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message(&
            'MRSF requires a triplet ROHF/UHF internal reference (mult=3).', with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    mrst = infos%tddft%mult
    target_label = mrsf_state_label(mrst, infos%tddft%target_state)

  ! Files open
    open (unit=iw, file=infos%log_filename, position="append")
  !
    call print_module_info('MRSF_Grad','Computing Gradient of '//trim(method_name))
!
    write(iw,'(/5X,"Gradient options"/&
                &5X,18("-")/&
                &5X,"Physical target state: ",A/&
                &5X,"Internal response root: ",I8/)')&
                & trim(target_label), infos%tddft%target_state

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

   ! Input parameters
  ! Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    s_size = (basis%nshell**2+basis%nshell)/2

!   Compute 1e gradient
    call flush(iw)

    call sf_1e_grad(infos, basis)

    write(iw,"(' ..... End Of 1-Eelectron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    allocate(v(nbf,nbf,2), source=0.0d0)
    allocate(d(nbf,nbf,2), source=0.0d0)
    allocate(p(nbf,nbf,2), source=0.0d0)
    allocate(spc(7,nbf,nbf), source=0.0d0)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    if (mrst==1 .or. mrst==3) then
      call data_has_tags(infos%dat, tags_mrsf, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_density)
    endif

    call unpack_matrix(td_p(:,1), p(:,:,1))
    call unpack_matrix(td_p(:,2), p(:,:,2))

    call unpack_matrix(dmat_a, d(:,:,1))
    call unpack_matrix(dmat_b, d(:,:,2))

    v(:,:,1) = td_abxc
    if (mrst==1 .or. mrst==3) then
      spc(1:7,:,:) = td_mrsf_density
    end if

!   Compute xc gradient
    if (dft) then
      call dft_initialize(infos, basis, molGrid, verbose=.true.)

      call utddft_xc_gradient(basis=basis, &
           molGrid=molGrid, &
           dedft=infos%atoms%grad, &
           da=d(:,:,1), &
           db=d(:,:,2), &
           pa=p(:,:,1:1), &
           pb=p(:,:,2:2), &
           nmtx=1, &
           !threshold=1.0d-15, &
           threshold=0.0d0, &
           infos=infos)

      ! Lee Eq. (3.16) identifies P=T+Z as the relaxed linear XC probe.
      ! The call above contains the established ground-state and fixed-grid
      ! AO/basis terms, but it does not differentiate the normalized
      ! atom-centred quadrature weights or the motion of each grid slice with
      ! its owner atom.  Add only those two moving-grid contributions in a
      ! separate sweep.  Combining the ground-state and probe responses here
      ! would differentiate a different quantity and double count terms.
      block
        real(kind=dp), allocatable :: grid_correction(:,:), grid_d(:,:,:), &
                                      grid_p(:,:,:)
        allocate(grid_correction(3,infos%mol_prop%natom), &
                 grid_d(nbf,nbf,2), grid_p(nbf,nbf,2), source=0.0_dp)
        ! utddft_xc_gradient temporarily applies AO normalization factors to
        ! its density arguments.  Copies avoid a second normalization round
        ! trip changing the arrays used later by mrsf_2e_grad.
        grid_d = d
        grid_p = p
        call utddft_xc_gradient(basis=basis, &
             molGrid=molGrid, &
             dedft=grid_correction, &
             da=grid_d(:,:,1), &
             db=grid_d(:,:,2), &
             pa=grid_p(:,:,1:1), &
             pb=grid_p(:,:,2:2), &
             nmtx=1, &
             threshold=0.0_dp, &
             infos=infos, &
             include_ground_state=.false., &
             include_weight_derivative=.true., &
             weight_derivative_only=.true.)
        infos%atoms%grad = infos%atoms%grad + grid_correction

        ! The total excited-state energy also contains the ROKS reference
        ! energy.  Differentiate its finite XC quadrature separately: using a
        ! zero probe prevents the relaxed P response above from being mixed
        ! into the ground-state owner-motion term.
        grid_correction = 0.0_dp
        grid_p = 0.0_dp
        call utddft_xc_gradient(basis=basis, &
             molGrid=molGrid, &
             dedft=grid_correction, &
             da=grid_d(:,:,1), &
             db=grid_d(:,:,2), &
             pa=grid_p(:,:,1:1), &
             pb=grid_p(:,:,2:2), &
             nmtx=1, &
             threshold=0.0_dp, &
             infos=infos, &
             include_ground_state=.true., &
             include_weight_derivative=.true., &
             weight_derivative_only=.true.)
        infos%atoms%grad = infos%atoms%grad + grid_correction
        deallocate(grid_correction, grid_d, grid_p)
      end block

      call dftclean(infos)
      call measure_time(print_total=1, log_unit=iw)
      call flush(iw)
    end if

!   Compute 2e gradient
    if (mrst==1 .or. mrst==3) then
      call mrsf_2e_grad(basis, infos, d, p, spc, v(:,:,1))
    else if (mrst==5) then
      call sf_2e_grad(basis, infos, d, p, v(:,:,1))
    end if

    call print_gradient(infos)

!   Phase 11 self-test (opt-in): verify the bilinear NAC density type
!   reproduces this production 2e gradient bit-for-bit at I=J.
    block
      character(len=8) :: env_selftest
      integer :: env_len, env_stat
      call get_environment_variable("OQP_NAC_SELFTEST", env_selftest, &
                                    env_len, env_stat)
      if (env_stat == 0 .and. env_len > 0) call mrsf_nac_amp_selftest(infos)
    end block

!   Print timings
    call measure_time(print_total=1, log_unit=iw)

    close(iw)

  end subroutine tdhf_mrsf_gradient

!###############################################################################

  subroutine mrsf_nac_overlap_C(c_handle) bind(C, name="mrsf_nac_overlap")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_overlap(inf)
  end subroutine mrsf_nac_overlap_C

!> @brief Orbital-overlap (Pulay) contribution to the first-order NAC for all
!>   MRSF state pairs: d^ov_{IJ,Ac} = sum_uv (C gamma^IJ C^T)_uv <chi_u|d_{Ac} chi_v>,
!>   the interstate transition density contracted with the ket-half overlap
!>   derivative. This is the dominant part of the derivative coupling and is
!>   absent from any energy gradient; it complements the Hellmann-Feynman/CI
!>   part computed by gradient polarization in the Python driver. Result stored
!>   in OQP::nac_overlap, shape (3, natom, nstate, nstate).
  subroutine mrsf_nac_overlap(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use, intrinsic :: iso_c_binding, only: c_int32_t
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use grd1, only: der_overlap_matrix_ket, der_overlap_matrix
    use tdhf_mrsf_lib, only: get_mrsf_transition_density
    use mathlib, only: orthogonal_transform

    implicit none

    character(len=*), parameter :: subroutine_name = "mrsf_nac_overlap"
    character(len=*), parameter :: OQP_nac_overlap = "OQP::nac_overlap"
    character(len=*), parameter :: OQP_nac_trden = "OQP::nac_trden_mo"
    character(len=*), parameter :: OQP_nac_gamma = "OQP::nac_gamma_tlf"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    integer :: nbf, natom, nstate, i, j, ist, jst, mu, nu, c, a, ok
    integer :: noca, nocb, p, q, isp, jsq
    integer(c_int32_t) :: gstat, gtag_id
    character(len=80) :: tags_gamma(1)
    real(kind=dp), contiguous, pointer :: gam_tlf(:,:,:)
    logical :: have_custom
    real(kind=dp), allocatable :: dSket(:,:,:,:), dSfull(:,:,:,:), &
                                  trden(:,:), trden_ss(:,:), trden_ao(:,:), &
                                  tmp(:,:), gnorm(:,:), gnorm2(:,:)
    real(kind=dp), pointer :: nac_ov(:,:,:), trden_st(:,:,:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), bvec_mo(:,:)
    real(kind=dp) :: acc
    character(len=*), parameter :: tags_required(2) = (/ character(len=80) :: &
      OQP_VEC_MO_A, OQP_td_bvec_mo /)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    natom = ubound(infos%atoms%grad, 2)
    nstate = infos%tddft%nstate
    noca = infos%mol_prop%nelec_a
    nocb = infos%mol_prop%nelec_b

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)

    ! optional override: if OQP::nac_gamma_tlf is present (the TLF-overlap-
    ! consistent transition density, e.g. built in Python from the closed
    ! form), contract it instead of get_mrsf_transition_density's gamma.
    tags_gamma(1) = OQP_nac_gamma
    gstat = infos%dat%has_records(tags_gamma, gtag_id)
    have_custom = (gstat == ta_ok)
    if (have_custom) call tagarray_get_data(infos%dat, OQP_nac_gamma, gam_tlf)

    call infos%dat%remove_records((/ character(len=80) :: OQP_nac_overlap, &
                                                          OQP_nac_trden /))
    call infos%dat%reserve_data(OQP_nac_overlap, ta_type_real64, &
          3*natom*nstate*nstate, (/ 3*natom, nstate, nstate /))
    call tagarray_get_data(infos%dat, OQP_nac_overlap, nac_ov)
    nac_ov = 0.0_dp
    ! also export the MO interstate transition densities for diagnostics
    call infos%dat%reserve_data(OQP_nac_trden, ta_type_real64, &
          nbf*nbf*nstate*nstate, (/ nbf*nbf, nstate, nstate /))
    call tagarray_get_data(infos%dat, OQP_nac_trden, trden_st)
    trden_st = 0.0_dp

    allocate(dSket(nbf,nbf,3,natom), dSfull(nbf,nbf,3,natom), &
             trden(nbf,nbf), trden_ss(nbf,nbf), trden_ao(nbf,nbf), &
             tmp(nbf,nbf), gnorm(nbf,nbf), gnorm2(nbf,nbf), &
             source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)

    ! ket-half AO overlap derivative <chi_u | d_{A,c} chi_v> (frozen term)
    call der_overlap_matrix_ket(basis, dSket)
    ! full AO overlap derivative dS_uv/dR (constraint term)
    call der_overlap_matrix(basis, dSfull)

    do ist = 1, nstate
      do jst = 1, nstate
        if (ist == jst) cycle
        ! interstate transition density (MO), then -> AO: C * trden * C^T
        if (have_custom) then
          trden = reshape(gam_tlf(:, ist, jst), (/ nbf, nbf /))
        else
          call get_mrsf_transition_density(infos, trden, bvec_mo, ist, jst)
        end if
        trden_st(:, ist, jst) = reshape(trden, (/ nbf*nbf /))
        call orthogonal_transform('t', nbf, mo_a, trden, trden_ao, tmp)
        ! apply basis normalization (dSket is in unnormalized convention)
        do nu = 1, nbf
          do mu = 1, nbf
            gnorm(mu,nu) = trden_ao(mu,nu)*basis%bfnrm(mu)*basis%bfnrm(nu)
          end do
        end do
        ! Term (3) frozen S^x-half: sum_uv gnorm_uv dSket(u,v,c,A) on all blocks
        ! Skeleton S^[x] terms from eliminating the dependent U^x blocks
        ! (orthonormality U^x_pq + U^x_qp = -S^[x]_pq), contracted with the
        ! FULL overlap derivative * (-1):
        !   same-space blocks (doc-doc, socc-socc, virt-virt): U^x = -1/2 S^[x]
        !     -> weight 1/2 (constraint term),
        !   cross-space (lo,hi) blocks (doc-socc, doc-virt, socc-virt with the
        !     row in the LOWER space): U^x_(lo,hi) = -S^[x] - U^x_(hi,lo)
        !     -> weight 1; the -U^x_(hi,lo) part goes to the CPHF RHS
        !        (antisymmetrized gamma in build_mrsf_zvector_rhs).
        ! Spaces: doc=[1,nocb], socc=[nocb+1,noca], virt=[noca+1,nbf].
        trden_ss = 0.0_dp
        do q = 1, nbf
          jsq = merge(1, merge(2, 3, q<=noca), q<=nocb)
          do p = 1, nbf
            isp = merge(1, merge(2, 3, p<=noca), p<=nocb)
            if (isp == jsq) then
              trden_ss(p,q) = 0.5_dp*trden(p,q)
            else if (isp < jsq) then
              trden_ss(p,q) = trden(p,q)
            end if
          end do
        end do
        call orthogonal_transform('t', nbf, mo_a, trden_ss, trden_ao, tmp)
        do nu = 1, nbf
          do mu = 1, nbf
            gnorm2(mu,nu) = trden_ao(mu,nu)*basis%bfnrm(mu)*basis%bfnrm(nu)
          end do
        end do
        do a = 1, natom
          do c = 1, 3
            acc = 0.0_dp
            do nu = 1, nbf
              do mu = 1, nbf
                acc = acc + gnorm(mu,nu)*dSket(mu,nu,c,a) &                ! frozen
                          - gnorm2(mu,nu)*dSfull(mu,nu,c,a)               ! skeleton U^x elim.
              end do
            end do
            nac_ov((a-1)*3+c, ist, jst) = acc
          end do
        end do
      end do
    end do

    deallocate(dSket, trden, trden_ao, tmp, gnorm)

  end subroutine mrsf_nac_overlap

!###############################################################################

!> @brief The driver for the two electron gradient
  subroutine mrsf_2e_grad(basis, infos, d, p, spc, v)

    use basis_tools, only: basis_set
    use precision, only: dp
    use messages, only: show_message, WITH_ABORT
    use types, only: information
    use parallel, only: par_env_t
    use routec_bridge, only: routec_try_grad2_mrsf

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set) :: basis
    real(kind=dp), contiguous, target :: p(:,:,:), d(:,:,:), spc(:,:,:), v(:,:)

    logical :: urohf, dft
    character(len=16) :: method_name
    real(kind=dp) :: scale_exch  !> HF scale in Reference
    real(kind=dp) :: scale_exch2 !> HF scale in Response

    integer :: ok
    real(kind=dp), allocatable :: de(:,:)
    class(grd2_compute_data_t), allocatable :: gcomp

    ! Route-C external MRSF 2e-gradient seam
    type(par_env_t) :: pe
    integer :: ok_ext

    dft = infos%control%hamilton == 20 ! dft or hf
    if (dft) then
      method_name = 'MRSF-TDDFT'
    else
      method_name = 'MRSF-TDHF'
    end if
    urohf = infos%control%scftype >= 2

    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
      scale_exch = infos%dft%HFscale
      scale_exch2 = infos%tddft%HFscale
    end if

    allocate(de(3,ubound(infos%atoms%zn,1)), &
            source=0.0d0, &
            stat=ok)

    if(ok/=0) call show_message('cannot allocate memory', WITH_ABORT)

    write(*, '(/7x,"Fitting parameters for ",A)') trim(method_name)
    if (.not.infos%dft%cam_flag) then
      write(*, '(10x,"Exact HF exchange:")')
      write(*, '(5x,"Reference: |", t20, f6.3, t29, "|")') scale_exch
      write(*, '(5x,"Response:  |", t20, f6.3, t29, "|")') scale_exch2
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

    ! Route-C external MRSF DF 2e-gradient (inert unless $OQP_ROUTEC_GRAD_LIB
    ! exposes routec_grad2_mrsf). No CAM (the attenuated short-range pass stays
    ! native); d/p/spc passed RAW (before the init() sum/diff transform).
    ! Rank 0 computes, result broadcast; declined -> native.
    if (.not. (dft .and. infos%dft%cam_flag)) then
      call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)
      ok_ext = 0
      if (pe%rank == 0) then
        if (routec_try_grad2_mrsf(d, p, spc, infos%atoms%xyz, de, &
                                  basis%nbf, ubound(infos%atoms%zn,1), &
                                  scale_exch, scale_exch2, 1.0d0, &
                                  [infos%tddft%spc_coco, &
                                   infos%tddft%spc_ovov, &
                                   infos%tddft%spc_coov], &
                                  infos%tddft%mult)) ok_ext = 1
      end if
      call pe%bcast(ok_ext, 1)
      if (ok_ext == 1) then
        call pe%bcast(de, size(de))
        infos%atoms%grad = infos%atoms%grad + de
        return
      end if
      de = 0.0d0
    end if

    if (allocated(gcomp)) deallocate(gcomp)
    allocate(gcomp, source=grd2_mrsf_compute_data_t( d2 = d &
                                    , p2 = p &
                                    , spc2 = spc &
                                    , nbf = basis%nbf &
                                    , hfscale = scale_exch &
                                    , hfscale2 = scale_exch2 &
                                    , spcscale = [infos%tddft%spc_coco, &
                                                  infos%tddft%spc_ovov, &
                                                  infos%tddft%spc_coov] &
                                    , mrst = infos%tddft%mult ))

    call gcomp%init()

    select type (gcomp)
    class is (grd2_mrsf_compute_data_t)
      call gcomp%build_cart(basis)
    end select

    ! Opt in to the petite reduction: the density contracted here is the
    ! (totally symmetric) converged one, and the skeleton gradient this
    ! produces is projected afterwards by Molecule.symmetrize_gradient. The
    ! opt-in set and the projection set are the same four sites by
    ! construction. Callers that contract a NON-symmetric density -- the CPHF
    ! probes in fock_deriv, and hf_hessian's displaced-geometry resp_grad --
    ! deliberately do not opt in.
    call grd2_driver(infos, basis, de, gcomp, &
                     cam = dft.and.infos%dft%cam_flag, &
                     alpha = infos%tddft%cam_alpha, &
                     beta = infos%tddft%cam_beta, &
                     mu = infos%tddft%cam_mu, &
                     petite = .true.)

    infos%atoms%grad = infos%atoms%grad + de

    call gcomp%clean()

  end subroutine

!###############################################################################

  subroutine grd2_mrsf_compute_data_t_init(this)
    implicit none
    class(grd2_mrsf_compute_data_t), target, intent(inout) :: this

    call this%clean()

    this%d2(:,:,1) = this%d2(:,:,1) +   this%d2(:,:,2)
    this%d2(:,:,2) = this%d2(:,:,1) - 2*this%d2(:,:,2)

    this%p2(:,:,1) = this%p2(:,:,1) +   this%p2(:,:,2)
    this%p2(:,:,2) = this%p2(:,:,1) - 2*this%p2(:,:,2)

  end subroutine

!###############################################################################

!> @brief Cartesian-effective copies of the MRSF gradient densities (d/p alpha
!>   and beta + the seven spin-pair-coupling densities) for HARMONIC_ACTIVE.
!>   Call AFTER init (which combines the spin densities). The spc densities may
!>   be non-symmetric; the per-block expansion handles that.
  subroutine grd2_mrsf_build_cart(this, basis)
    class(grd2_mrsf_compute_data_t), intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, allocatable :: od(:)
    integer :: nc
    if (.not. HARMONIC_ACTIVE) return
    call mrsf_cart_one(basis, this%d2(:,:,1), this%d2a_c, this%cart_off, nc)
    call mrsf_cart_one(basis, this%d2(:,:,2), this%d2b_c, od, nc)
    call mrsf_cart_one(basis, this%p2(:,:,1), this%p2a_c, od, nc)
    call mrsf_cart_one(basis, this%p2(:,:,2), this%p2b_c, od, nc)
    call mrsf_cart_one(basis, this%spc2(7,:,:), this%ball_c, od, nc)
    call mrsf_cart_one(basis, this%spc2(1,:,:), this%bo2v_c, od, nc)
    call mrsf_cart_one(basis, this%spc2(2,:,:), this%bo1v_c, od, nc)
    call mrsf_cart_one(basis, this%spc2(3,:,:), this%bco1_c, od, nc)
    call mrsf_cart_one(basis, this%spc2(4,:,:), this%bco2_c, od, nc)
    call mrsf_cart_one(basis, this%spc2(5,:,:), this%o21v_c, od, nc)
    call mrsf_cart_one(basis, this%spc2(6,:,:), this%co12_c, od, nc)
  end subroutine grd2_mrsf_build_cart

  subroutine mrsf_cart_one(basis, m, m_cart, off, nc)
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: m(:,:)
    real(kind=dp), allocatable, intent(out) :: m_cart(:,:)
    integer, allocatable, intent(out) :: off(:)
    integer, intent(out) :: nc
    real(kind=dp), allocatable :: tmp(:,:)
    tmp = m
    call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
    call build_cart_density(basis, tmp, m_cart, off, nc)
  end subroutine mrsf_cart_one

!###############################################################################

  subroutine grd2_mrsf_compute_data_t_clean(this)
    implicit none
    class(grd2_mrsf_compute_data_t), target, intent(inout) :: this
  end subroutine

!###############################################################################

!> @brief This routine forms the product of density
!>        matrices for use in forming the two electron
!>        gradient. Valid for closed and open shell SCF.
  subroutine grd2_mrsf_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_mrsf_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: xcfact, xcfact2, coulfact, df1, dq1, dt2, bfn
    real(kind=dp) :: qfspcp1, qfspcp2, qfspcp3, sgnk
    real(kind=dp) :: db1, db2, dc1, dc2, dc3, dc4, dd1, dd2, dd3, dd4
    real(kind=dp), pointer, dimension(:,:) :: &
      ball, bo2v, bo1v, bco1, bco2, co12, o21v, &
      d2a, d2b, p2a, p2b
    logical :: usecart
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    integer :: i1, j1, k1, l1

    coulfact = 4*this%coulscale
    xcfact = this%hfscale
    xcfact2 = this%hfscale2
    qfspcp1 = this%spcscale(1)
    qfspcp2 = this%spcscale(2)
    qfspcp3 = this%spcscale(3)

!   The spin-pair-coupling channels belong to the FULL-ERI pass only.
!   With a range-separated (CAM) functional grd2_driver runs this digest
!   twice: pass 1 over the full operator, pass 2 over the erf-attenuated
!   short-range one.  It re-points coulscale/hfscale/hfscale2/mu per pass,
!   but spcscale is ours and it does not touch it -- so without this guard
!   the six SPC contractions (co12, o21v, bco1*bo2v, bco2*bo1v) are added in
!   BOTH passes.  The energy adds them in pass 1 ONLY: int2_mrsf_data_t_update's
!   cur_pass==2 branch digests slot 7 (ball) alone, so slots 1-6 never receive
!   an attenuated-pass contribution.  Adding them again against the attenuated
!   derivative integrals therefore puts a term in the gradient that has no
!   counterpart in the energy, and MRSF analytic gradients disagreed with
!   finite differences by up to ~1e-3 Hartree/Bohr for any range-separated
!   functional -- state-dependently, in proportion to the state's
!   closed->open (spin-pair) character.  Global hybrids never run pass 2,
!   which is why only CAM-type functionals were affected.
    if (this%attenuated) then
      qfspcp1 = 0.0_dp
      qfspcp2 = 0.0_dp
      qfspcp3 = 0.0_dp
    end if

    sgnk = 1.0_dp
    if (this%mrst==3) sgnk = -1.0_dp
    dabmax = 0

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      d2a => this%d2a_c;  d2b => this%d2b_c;  p2a => this%p2a_c;  p2b => this%p2b_c
      ball => this%ball_c;  bo2v => this%bo2v_c;  bo1v => this%bo1v_c
      bco1 => this%bco1_c;  bco2 => this%bco2_c;  o21v => this%o21v_c;  co12 => this%co12_c
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
    else
      d2a => this%d2(:,:,1);  d2b => this%d2(:,:,2)
      p2a => this%p2(:,:,1);  p2b => this%p2(:,:,2)
      ball => this%spc2(7,:,:)
      bo2v => this%spc2(1,:,:);  bo1v => this%spc2(2,:,:)
      bco1 => this%spc2(3,:,:);  bco2 => this%spc2(4,:,:)
      o21v => this%spc2(5,:,:);  co12 => this%spc2(6,:,:)
      loc = basis%ao_offset(id) - 1
      nbf = basis%naos(id)
    end if

    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i

      do j = 1, nbf(2)
        j1 = loc(2) + j

        do k = 1, nbf(3)
          k1 = loc(3) + k

          do l = 1, nbf(4)
            l1 = loc(4) + l
            df1 = (d2a(i1,j1)+p2a(i1,j1))*d2a(k1,l1) &
                +  d2a(i1,j1)                  *p2a(k1,l1)
            df1 = df1 * coulfact

            if (xcfact /= 0.0_dp .or. xcfact2 /= 0.0_dp) then
              dq1 = (d2a(i1,k1)+p2a(i1,k1))*d2a(j1,l1) &
                  +  d2a(i1,k1)                  *p2a(j1,l1) &
                  + (d2a(i1,l1)+p2a(i1,l1))*d2a(j1,k1) &
                  +  d2a(i1,l1)                  *p2a(j1,k1) &
                  + (d2b(i1,k1)+p2b(i1,k1))*d2b(j1,l1) &
                  +  d2b(i1,k1)                  *p2b(j1,l1) &
                  + (d2b(i1,l1)+p2b(i1,l1))*d2b(j1,k1) &
                  +  d2b(i1,l1)                  *p2b(j1,k1)
              dt2 = ball(i1,k1)*ball(j1,l1) &
                  + ball(k1,i1)*ball(l1,j1) &
                  + ball(i1,l1)*ball(j1,k1) &
                  + ball(l1,i1)*ball(k1,j1)

              df1 = df1-xcfact*dq1-xcfact2*2.0_dp*dt2
            end if

            if (qfspcp1 /= 0.0_dp) then
              db1 =  co12(i1,k1)*co12(l1,j1) &
                   + co12(i1,l1)*co12(k1,j1) &
                   + co12(j1,k1)*co12(l1,i1) &
                   + co12(j1,l1)*co12(k1,i1) &
                   + co12(l1,j1)*co12(i1,k1) &
                   + co12(k1,j1)*co12(i1,l1) &
                   + co12(l1,i1)*co12(j1,k1) &
                   + co12(k1,i1)*co12(j1,l1)

              df1 = df1 + sgnk*qfspcp1*db1
            end if

            if (qfspcp2 /= 0.0_dp) then
              db2 =  o21v(i1,k1)*o21v(l1,j1) &
                   + o21v(i1,l1)*o21v(k1,j1) &
                   + o21v(j1,k1)*o21v(l1,i1) &
                   + o21v(j1,l1)*o21v(k1,i1) &
                   + o21v(l1,j1)*o21v(i1,k1) &
                   + o21v(k1,j1)*o21v(i1,l1) &
                   + o21v(l1,i1)*o21v(j1,k1) &
                   + o21v(k1,i1)*o21v(j1,l1)

              df1 = df1 + sgnk*qfspcp2*db2
            end if

            if (qfspcp3 /= 0.0_dp) then
              dc1 =  bco1(i1,k1)*bo2v(j1,l1) &
                   + bco1(i1,l1)*bo2v(j1,k1) &
                   + bco1(j1,k1)*bo2v(i1,l1) &
                   + bco1(j1,l1)*bo2v(i1,k1) &
                   + bco1(l1,j1)*bo2v(k1,i1) &
                   + bco1(k1,j1)*bo2v(l1,i1) &
                   + bco1(l1,i1)*bo2v(k1,j1) &
                   + bco1(k1,i1)*bo2v(l1,j1)

              dc2 =  bco2(i1,k1)*bo1v(j1,l1) &
                   + bco2(i1,l1)*bo1v(j1,k1) &
                   + bco2(j1,k1)*bo1v(i1,l1) &
                   + bco2(j1,l1)*bo1v(i1,k1) &
                   + bco2(l1,j1)*bo1v(k1,i1) &
                   + bco2(k1,j1)*bo1v(l1,i1) &
                   + bco2(l1,i1)*bo1v(k1,j1) &
                   + bco2(k1,i1)*bo1v(l1,j1)

              dc3 =  bo2v(i1,k1)*bco1(j1,l1) &
                   + bo2v(i1,l1)*bco1(j1,k1) &
                   + bo2v(j1,k1)*bco1(i1,l1) &
                   + bo2v(j1,l1)*bco1(i1,k1) &
                   + bo2v(l1,j1)*bco1(k1,i1) &
                   + bo2v(k1,j1)*bco1(l1,i1) &
                   + bo2v(l1,i1)*bco1(k1,j1) &
                   + bo2v(k1,i1)*bco1(l1,j1)

              dc4 =  bo1v(i1,k1)*bco2(j1,l1) &
                   + bo1v(i1,l1)*bco2(j1,k1) &
                   + bo1v(j1,k1)*bco2(i1,l1) &
                   + bo1v(j1,l1)*bco2(i1,k1) &
                   + bo1v(l1,j1)*bco2(k1,i1) &
                   + bo1v(k1,j1)*bco2(l1,i1) &
                   + bo1v(l1,i1)*bco2(k1,j1) &
                   + bo1v(k1,i1)*bco2(l1,j1)

              dd1 =  bco1(i1,j1)*bo2v(l1,k1) &
                   + bco1(i1,j1)*bo2v(k1,l1) &
                   + bco1(j1,i1)*bo2v(l1,k1) &
                   + bco1(j1,i1)*bo2v(k1,l1) &
                   + bco1(l1,k1)*bo2v(i1,j1) &
                   + bco1(k1,l1)*bo2v(i1,j1) &
                   + bco1(l1,k1)*bo2v(j1,i1) &
                   + bco1(k1,l1)*bo2v(j1,i1)

              dd2 =  bco2(i1,j1)*bo1v(l1,k1) &
                   + bco2(i1,j1)*bo1v(k1,l1) &
                   + bco2(j1,i1)*bo1v(l1,k1) &
                   + bco2(j1,i1)*bo1v(k1,l1) &
                   + bco2(l1,k1)*bo1v(i1,j1) &
                   + bco2(k1,l1)*bo1v(i1,j1) &
                   + bco2(l1,k1)*bo1v(j1,i1) &
                   + bco2(k1,l1)*bo1v(j1,i1)

              dd3 =  bo2v(i1,j1)*bco1(l1,k1) &
                   + bo2v(i1,j1)*bco1(k1,l1) &
                   + bo2v(j1,i1)*bco1(l1,k1) &
                   + bo2v(j1,i1)*bco1(k1,l1) &
                   + bo2v(l1,k1)*bco1(i1,j1) &
                   + bo2v(k1,l1)*bco1(i1,j1) &
                   + bo2v(l1,k1)*bco1(j1,i1) &
                   + bo2v(k1,l1)*bco1(j1,i1)

              dd4 =  bo1v(i1,j1)*bco2(l1,k1) &
                   + bo1v(i1,j1)*bco2(k1,l1) &
                   + bo1v(j1,i1)*bco2(l1,k1) &
                   + bo1v(j1,i1)*bco2(k1,l1) &
                   + bo1v(l1,k1)*bco2(i1,j1) &
                   + bo1v(k1,l1)*bco2(i1,j1) &
                   + bo1v(l1,k1)*bco2(j1,i1) &
                   + bo1v(k1,l1)*bco2(j1,i1)

              df1  = df1 + sgnk*qfspcp3*(-dc1-dc2-dc3-dc4 &
                                         +dd1+dd2+dd3+dd4)
            end if

            dabmax = max(dabmax, abs(df1))
            bfn = 1.0_dp
            if (.not. usecart) bfn = product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i) = df1*bfn
          end do
        end do
      end do
    end do
  end subroutine grd2_mrsf_compute_data_t_get_density

!###############################################################################
!  Phase 11: NAC amplitude-term bilinear compute-data type
!###############################################################################

  subroutine grd2_mrsf_nac_compute_data_t_init(this)
    implicit none
    class(grd2_mrsf_nac_compute_data_t), target, intent(inout) :: this

    call this%clean()

    ! same alpha/beta -> total/spin recombination as the production type;
    ! applied to the shared reference (d2) and interstate relaxed (p2) densities
    this%d2(:,:,1) = this%d2(:,:,1) +   this%d2(:,:,2)
    this%d2(:,:,2) = this%d2(:,:,1) - 2*this%d2(:,:,2)

    this%p2(:,:,1) = this%p2(:,:,1) +   this%p2(:,:,2)
    this%p2(:,:,2) = this%p2(:,:,1) - 2*this%p2(:,:,2)

  end subroutine

!###############################################################################

  subroutine grd2_mrsf_nac_compute_data_t_clean(this)
    implicit none
    class(grd2_mrsf_nac_compute_data_t), target, intent(inout) :: this
  end subroutine

!###############################################################################

!> @brief Bilinear four-index density product for the NAC amplitude term
!>        G_IJ = X_I^T (d_x A) X_J. Identical in structure to
!>        grd2_mrsf_compute_data_t_get_density but every pure transition-density
!>        product f(a)*g(b) is replaced by 1/2[ f_I(a) g_J(b) + f_J(a) g_I(b) ].
!>        Collapses bit-for-bit to the production form at I=J.
  subroutine grd2_mrsf_nac_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_mrsf_nac_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: xcfact, xcfact2, coulfact, df1, dq1, dt2
    real(kind=dp) :: qfspcp1, qfspcp2, qfspcp3, sgnk
    real(kind=dp) :: db1, db2, dc1, dc2, dc3, dc4, dd1, dd2, dd3, dd4
    real(kind=dp), pointer, dimension(:,:) :: &
      ballI, bo2vI, bo1vI, bco1I, bco2I, co12I, o21vI, &
      ballJ, bo2vJ, bo1vJ, bco1J, bco2J, co12J, o21vJ
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    integer :: i1, j1, k1, l1

    ! state-I transition/amplitude channels
    ballI => this%spcI(7,:,:)
    bo2vI => this%spcI(1,:,:)
    bo1vI => this%spcI(2,:,:)
    bco1I => this%spcI(3,:,:)
    bco2I => this%spcI(4,:,:)
    o21vI => this%spcI(5,:,:)
    co12I => this%spcI(6,:,:)
    ! state-J transition/amplitude channels
    ballJ => this%spcJ(7,:,:)
    bo2vJ => this%spcJ(1,:,:)
    bo1vJ => this%spcJ(2,:,:)
    bco1J => this%spcJ(3,:,:)
    bco2J => this%spcJ(4,:,:)
    o21vJ => this%spcJ(5,:,:)
    co12J => this%spcJ(6,:,:)

    coulfact = 4*this%coulscale
    xcfact = this%hfscale
    xcfact2 = this%hfscale2
    qfspcp1 = this%spcscale(1)
    qfspcp2 = this%spcscale(2)
    qfspcp3 = this%spcscale(3)

    sgnk = 1.0_dp
    if (this%mrst==3) sgnk = -1.0_dp
    dabmax = 0
    loc = basis%ao_offset(id)-1

    nbf = basis%naos(id)

    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i

      do j = 1, nbf(2)
        j1 = loc(2) + j

        do k = 1, nbf(3)
          k1 = loc(3) + k

          do l = 1, nbf(4)
            l1 = loc(4) + l
            ! Coulomb + relaxed/reference: reference (d2) state-independent,
            ! relaxed (p2) already the interstate object -> production form.
            df1 = (this%d2(i1,j1,1)+this%p2(i1,j1,1))*this%d2(k1,l1,1) &
                +  this%d2(i1,j1,1)                  *this%p2(k1,l1,1)
            df1 = df1 * coulfact

            if (xcfact /= 0.0_dp .or. xcfact2 /= 0.0_dp) then
              dq1 = (this%d2(i1,k1,1)+this%p2(i1,k1,1))*this%d2(j1,l1,1) &
                  +  this%d2(i1,k1,1)                  *this%p2(j1,l1,1) &
                  + (this%d2(i1,l1,1)+this%p2(i1,l1,1))*this%d2(j1,k1,1) &
                  +  this%d2(i1,l1,1)                  *this%p2(j1,k1,1) &
                  + (this%d2(i1,k1,2)+this%p2(i1,k1,2))*this%d2(j1,l1,2) &
                  +  this%d2(i1,k1,2)                  *this%p2(j1,l1,2) &
                  + (this%d2(i1,l1,2)+this%p2(i1,l1,2))*this%d2(j1,k1,2) &
                  +  this%d2(i1,l1,2)                  *this%p2(j1,k1,2)
              ! channel-7 exchange (ball), symmetrised I<->J
              dt2 = 0.5_dp*(ballI(i1,k1)*ballJ(j1,l1) + ballJ(i1,k1)*ballI(j1,l1)) &
                  + 0.5_dp*(ballI(k1,i1)*ballJ(l1,j1) + ballJ(k1,i1)*ballI(l1,j1)) &
                  + 0.5_dp*(ballI(i1,l1)*ballJ(j1,k1) + ballJ(i1,l1)*ballI(j1,k1)) &
                  + 0.5_dp*(ballI(l1,i1)*ballJ(k1,j1) + ballJ(l1,i1)*ballI(k1,j1))

              df1 = df1-xcfact*dq1-xcfact2*2.0_dp*dt2
            end if

            if (qfspcp1 /= 0.0_dp) then
              db1 = 0.5_dp*(co12I(i1,k1)*co12J(l1,j1) + co12J(i1,k1)*co12I(l1,j1)) &
                  + 0.5_dp*(co12I(i1,l1)*co12J(k1,j1) + co12J(i1,l1)*co12I(k1,j1)) &
                  + 0.5_dp*(co12I(j1,k1)*co12J(l1,i1) + co12J(j1,k1)*co12I(l1,i1)) &
                  + 0.5_dp*(co12I(j1,l1)*co12J(k1,i1) + co12J(j1,l1)*co12I(k1,i1)) &
                  + 0.5_dp*(co12I(l1,j1)*co12J(i1,k1) + co12J(l1,j1)*co12I(i1,k1)) &
                  + 0.5_dp*(co12I(k1,j1)*co12J(i1,l1) + co12J(k1,j1)*co12I(i1,l1)) &
                  + 0.5_dp*(co12I(l1,i1)*co12J(j1,k1) + co12J(l1,i1)*co12I(j1,k1)) &
                  + 0.5_dp*(co12I(k1,i1)*co12J(j1,l1) + co12J(k1,i1)*co12I(j1,l1))

              df1 = df1 + sgnk*qfspcp1*db1
            end if

            if (qfspcp2 /= 0.0_dp) then
              db2 = 0.5_dp*(o21vI(i1,k1)*o21vJ(l1,j1) + o21vJ(i1,k1)*o21vI(l1,j1)) &
                  + 0.5_dp*(o21vI(i1,l1)*o21vJ(k1,j1) + o21vJ(i1,l1)*o21vI(k1,j1)) &
                  + 0.5_dp*(o21vI(j1,k1)*o21vJ(l1,i1) + o21vJ(j1,k1)*o21vI(l1,i1)) &
                  + 0.5_dp*(o21vI(j1,l1)*o21vJ(k1,i1) + o21vJ(j1,l1)*o21vI(k1,i1)) &
                  + 0.5_dp*(o21vI(l1,j1)*o21vJ(i1,k1) + o21vJ(l1,j1)*o21vI(i1,k1)) &
                  + 0.5_dp*(o21vI(k1,j1)*o21vJ(i1,l1) + o21vJ(k1,j1)*o21vI(i1,l1)) &
                  + 0.5_dp*(o21vI(l1,i1)*o21vJ(j1,k1) + o21vJ(l1,i1)*o21vI(j1,k1)) &
                  + 0.5_dp*(o21vI(k1,i1)*o21vJ(j1,l1) + o21vJ(k1,i1)*o21vI(j1,l1))

              df1 = df1 + sgnk*qfspcp2*db2
            end if

            if (qfspcp3 /= 0.0_dp) then
              dc1 = 0.5_dp*(bco1I(i1,k1)*bo2vJ(j1,l1) + bco1J(i1,k1)*bo2vI(j1,l1)) &
                  + 0.5_dp*(bco1I(i1,l1)*bo2vJ(j1,k1) + bco1J(i1,l1)*bo2vI(j1,k1)) &
                  + 0.5_dp*(bco1I(j1,k1)*bo2vJ(i1,l1) + bco1J(j1,k1)*bo2vI(i1,l1)) &
                  + 0.5_dp*(bco1I(j1,l1)*bo2vJ(i1,k1) + bco1J(j1,l1)*bo2vI(i1,k1)) &
                  + 0.5_dp*(bco1I(l1,j1)*bo2vJ(k1,i1) + bco1J(l1,j1)*bo2vI(k1,i1)) &
                  + 0.5_dp*(bco1I(k1,j1)*bo2vJ(l1,i1) + bco1J(k1,j1)*bo2vI(l1,i1)) &
                  + 0.5_dp*(bco1I(l1,i1)*bo2vJ(k1,j1) + bco1J(l1,i1)*bo2vI(k1,j1)) &
                  + 0.5_dp*(bco1I(k1,i1)*bo2vJ(l1,j1) + bco1J(k1,i1)*bo2vI(l1,j1))

              dc2 = 0.5_dp*(bco2I(i1,k1)*bo1vJ(j1,l1) + bco2J(i1,k1)*bo1vI(j1,l1)) &
                  + 0.5_dp*(bco2I(i1,l1)*bo1vJ(j1,k1) + bco2J(i1,l1)*bo1vI(j1,k1)) &
                  + 0.5_dp*(bco2I(j1,k1)*bo1vJ(i1,l1) + bco2J(j1,k1)*bo1vI(i1,l1)) &
                  + 0.5_dp*(bco2I(j1,l1)*bo1vJ(i1,k1) + bco2J(j1,l1)*bo1vI(i1,k1)) &
                  + 0.5_dp*(bco2I(l1,j1)*bo1vJ(k1,i1) + bco2J(l1,j1)*bo1vI(k1,i1)) &
                  + 0.5_dp*(bco2I(k1,j1)*bo1vJ(l1,i1) + bco2J(k1,j1)*bo1vI(l1,i1)) &
                  + 0.5_dp*(bco2I(l1,i1)*bo1vJ(k1,j1) + bco2J(l1,i1)*bo1vI(k1,j1)) &
                  + 0.5_dp*(bco2I(k1,i1)*bo1vJ(l1,j1) + bco2J(k1,i1)*bo1vI(l1,j1))

              dc3 = 0.5_dp*(bo2vI(i1,k1)*bco1J(j1,l1) + bo2vJ(i1,k1)*bco1I(j1,l1)) &
                  + 0.5_dp*(bo2vI(i1,l1)*bco1J(j1,k1) + bo2vJ(i1,l1)*bco1I(j1,k1)) &
                  + 0.5_dp*(bo2vI(j1,k1)*bco1J(i1,l1) + bo2vJ(j1,k1)*bco1I(i1,l1)) &
                  + 0.5_dp*(bo2vI(j1,l1)*bco1J(i1,k1) + bo2vJ(j1,l1)*bco1I(i1,k1)) &
                  + 0.5_dp*(bo2vI(l1,j1)*bco1J(k1,i1) + bo2vJ(l1,j1)*bco1I(k1,i1)) &
                  + 0.5_dp*(bo2vI(k1,j1)*bco1J(l1,i1) + bo2vJ(k1,j1)*bco1I(l1,i1)) &
                  + 0.5_dp*(bo2vI(l1,i1)*bco1J(k1,j1) + bo2vJ(l1,i1)*bco1I(k1,j1)) &
                  + 0.5_dp*(bo2vI(k1,i1)*bco1J(l1,j1) + bo2vJ(k1,i1)*bco1I(l1,j1))

              dc4 = 0.5_dp*(bo1vI(i1,k1)*bco2J(j1,l1) + bo1vJ(i1,k1)*bco2I(j1,l1)) &
                  + 0.5_dp*(bo1vI(i1,l1)*bco2J(j1,k1) + bo1vJ(i1,l1)*bco2I(j1,k1)) &
                  + 0.5_dp*(bo1vI(j1,k1)*bco2J(i1,l1) + bo1vJ(j1,k1)*bco2I(i1,l1)) &
                  + 0.5_dp*(bo1vI(j1,l1)*bco2J(i1,k1) + bo1vJ(j1,l1)*bco2I(i1,k1)) &
                  + 0.5_dp*(bo1vI(l1,j1)*bco2J(k1,i1) + bo1vJ(l1,j1)*bco2I(k1,i1)) &
                  + 0.5_dp*(bo1vI(k1,j1)*bco2J(l1,i1) + bo1vJ(k1,j1)*bco2I(l1,i1)) &
                  + 0.5_dp*(bo1vI(l1,i1)*bco2J(k1,j1) + bo1vJ(l1,i1)*bco2I(k1,j1)) &
                  + 0.5_dp*(bo1vI(k1,i1)*bco2J(l1,j1) + bo1vJ(k1,i1)*bco2I(l1,j1))

              dd1 = 0.5_dp*(bco1I(i1,j1)*bo2vJ(l1,k1) + bco1J(i1,j1)*bo2vI(l1,k1)) &
                  + 0.5_dp*(bco1I(i1,j1)*bo2vJ(k1,l1) + bco1J(i1,j1)*bo2vI(k1,l1)) &
                  + 0.5_dp*(bco1I(j1,i1)*bo2vJ(l1,k1) + bco1J(j1,i1)*bo2vI(l1,k1)) &
                  + 0.5_dp*(bco1I(j1,i1)*bo2vJ(k1,l1) + bco1J(j1,i1)*bo2vI(k1,l1)) &
                  + 0.5_dp*(bco1I(l1,k1)*bo2vJ(i1,j1) + bco1J(l1,k1)*bo2vI(i1,j1)) &
                  + 0.5_dp*(bco1I(k1,l1)*bo2vJ(i1,j1) + bco1J(k1,l1)*bo2vI(i1,j1)) &
                  + 0.5_dp*(bco1I(l1,k1)*bo2vJ(j1,i1) + bco1J(l1,k1)*bo2vI(j1,i1)) &
                  + 0.5_dp*(bco1I(k1,l1)*bo2vJ(j1,i1) + bco1J(k1,l1)*bo2vI(j1,i1))

              dd2 = 0.5_dp*(bco2I(i1,j1)*bo1vJ(l1,k1) + bco2J(i1,j1)*bo1vI(l1,k1)) &
                  + 0.5_dp*(bco2I(i1,j1)*bo1vJ(k1,l1) + bco2J(i1,j1)*bo1vI(k1,l1)) &
                  + 0.5_dp*(bco2I(j1,i1)*bo1vJ(l1,k1) + bco2J(j1,i1)*bo1vI(l1,k1)) &
                  + 0.5_dp*(bco2I(j1,i1)*bo1vJ(k1,l1) + bco2J(j1,i1)*bo1vI(k1,l1)) &
                  + 0.5_dp*(bco2I(l1,k1)*bo1vJ(i1,j1) + bco2J(l1,k1)*bo1vI(i1,j1)) &
                  + 0.5_dp*(bco2I(k1,l1)*bo1vJ(i1,j1) + bco2J(k1,l1)*bo1vI(i1,j1)) &
                  + 0.5_dp*(bco2I(l1,k1)*bo1vJ(j1,i1) + bco2J(l1,k1)*bo1vI(j1,i1)) &
                  + 0.5_dp*(bco2I(k1,l1)*bo1vJ(j1,i1) + bco2J(k1,l1)*bo1vI(j1,i1))

              dd3 = 0.5_dp*(bo2vI(i1,j1)*bco1J(l1,k1) + bo2vJ(i1,j1)*bco1I(l1,k1)) &
                  + 0.5_dp*(bo2vI(i1,j1)*bco1J(k1,l1) + bo2vJ(i1,j1)*bco1I(k1,l1)) &
                  + 0.5_dp*(bo2vI(j1,i1)*bco1J(l1,k1) + bo2vJ(j1,i1)*bco1I(l1,k1)) &
                  + 0.5_dp*(bo2vI(j1,i1)*bco1J(k1,l1) + bo2vJ(j1,i1)*bco1I(k1,l1)) &
                  + 0.5_dp*(bo2vI(l1,k1)*bco1J(i1,j1) + bo2vJ(l1,k1)*bco1I(i1,j1)) &
                  + 0.5_dp*(bo2vI(k1,l1)*bco1J(i1,j1) + bo2vJ(k1,l1)*bco1I(i1,j1)) &
                  + 0.5_dp*(bo2vI(l1,k1)*bco1J(j1,i1) + bo2vJ(l1,k1)*bco1I(j1,i1)) &
                  + 0.5_dp*(bo2vI(k1,l1)*bco1J(j1,i1) + bo2vJ(k1,l1)*bco1I(j1,i1))

              dd4 = 0.5_dp*(bo1vI(i1,j1)*bco2J(l1,k1) + bo1vJ(i1,j1)*bco2I(l1,k1)) &
                  + 0.5_dp*(bo1vI(i1,j1)*bco2J(k1,l1) + bo1vJ(i1,j1)*bco2I(k1,l1)) &
                  + 0.5_dp*(bo1vI(j1,i1)*bco2J(l1,k1) + bo1vJ(j1,i1)*bco2I(l1,k1)) &
                  + 0.5_dp*(bo1vI(j1,i1)*bco2J(k1,l1) + bo1vJ(j1,i1)*bco2I(k1,l1)) &
                  + 0.5_dp*(bo1vI(l1,k1)*bco2J(i1,j1) + bo1vJ(l1,k1)*bco2I(i1,j1)) &
                  + 0.5_dp*(bo1vI(k1,l1)*bco2J(i1,j1) + bo1vJ(k1,l1)*bco2I(i1,j1)) &
                  + 0.5_dp*(bo1vI(l1,k1)*bco2J(j1,i1) + bo1vJ(l1,k1)*bco2I(j1,i1)) &
                  + 0.5_dp*(bo1vI(k1,l1)*bco2J(j1,i1) + bo1vJ(k1,l1)*bco2I(j1,i1))

              df1  = df1 + sgnk*qfspcp3*(-dc1-dc2-dc3-dc4 &
                                         +dd1+dd2+dd3+dd4)
            end if

            dabmax = max(dabmax, abs(df1))
            ab(l,k,j,i) = df1*product(basis%bfnrm([i1,j1,k1,l1]))
          end do
        end do
      end do
    end do
  end subroutine grd2_mrsf_nac_compute_data_t_get_density

!###############################################################################

!> @brief Phase 11 self-test: at I=J the bilinear NAC compute-data type must
!>        reproduce the production grd2_mrsf_compute_data_t two-electron
!>        gradient bit-for-bit. Reads the same tagarrays the production gradient
!>        consumed (DM_A/B reference, td_p relaxed, td_mrsf_density channels)
!>        for the current target state, runs both grd2 paths in one process, and
!>        prints max|de_nac - de_prod|. Triggered by env OQP_NAC_SELFTEST.
  subroutine mrsf_nac_amp_selftest(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use mathlib, only: unpack_matrix

    implicit none

    character(len=*), parameter :: subroutine_name = "mrsf_nac_amp_selftest"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:), &
                                          td_mrsf_density(:,:,:), td_p(:,:)
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_td_p /)
    character(len=*), parameter :: tags_mrsf(1) = (/ character(len=80) :: &
      OQP_td_mrsf_density /)

    real(kind=dp), allocatable, target :: dA(:,:,:), pA(:,:,:), &
                                          dB(:,:,:), pB(:,:,:), spc(:,:,:)
    real(kind=dp), allocatable :: deP(:,:), deN(:,:)
    class(grd2_compute_data_t), allocatable :: gP
    class(grd2_compute_data_t), allocatable :: gN
    real(kind=dp) :: scale_exch, scale_exch2, dmax, gmax
    logical :: dft, do_cam
    integer :: nbf, mrst, natom, ok

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf   = basis%nbf
    natom = ubound(infos%atoms%zn,1)
    mrst  = infos%tddft%mult

    dft = infos%control%hamilton == 20
    scale_exch  = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
      scale_exch  = infos%dft%HFscale
      scale_exch2 = infos%tddft%HFscale
    end if
    do_cam = dft .and. infos%dft%cam_flag

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call data_has_tags(infos%dat, tags_mrsf, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_density)

    allocate(dA(nbf,nbf,2), pA(nbf,nbf,2), dB(nbf,nbf,2), pB(nbf,nbf,2), &
             spc(7,nbf,nbf), source=0.0_dp, stat=ok)
    if (ok/=0) call show_message('cannot allocate memory', WITH_ABORT)

    call unpack_matrix(td_p(:,1), pA(:,:,1))
    call unpack_matrix(td_p(:,2), pA(:,:,2))
    call unpack_matrix(dmat_a, dA(:,:,1))
    call unpack_matrix(dmat_b, dA(:,:,2))
    spc(1:7,:,:) = td_mrsf_density
    ! independent copies (init mutates d2/p2 in place)
    dB = dA
    pB = pA

    allocate(deP(3,natom), deN(3,natom), source=0.0_dp)

    ! Path A: production quadratic type
    gP = grd2_mrsf_compute_data_t( d2 = dA, p2 = pA, spc2 = spc, &
                                   nbf = nbf, hfscale = scale_exch, &
                                   hfscale2 = scale_exch2, &
                                   spcscale = [infos%tddft%spc_coco, &
                                               infos%tddft%spc_ovov, &
                                               infos%tddft%spc_coov], &
                                   mrst = mrst )
    call gP%init()
    call grd2_driver(infos, basis, deP, gP, &
                     cam = do_cam, alpha = infos%tddft%cam_alpha, &
                     beta = infos%tddft%cam_beta, mu = infos%tddft%cam_mu)
    call gP%clean()

    ! Path B: bilinear NAC type at I=J (spcI == spcJ == spc)
    gN = grd2_mrsf_nac_compute_data_t( d2 = dB, p2 = pB, spcI = spc, spcJ = spc, &
                                       nbf = nbf, hfscale = scale_exch, &
                                       hfscale2 = scale_exch2, &
                                       spcscale = [infos%tddft%spc_coco, &
                                                   infos%tddft%spc_ovov, &
                                                   infos%tddft%spc_coov], &
                                       mrst = mrst )
    call gN%init()
    call grd2_driver(infos, basis, deN, gN, &
                     cam = do_cam, alpha = infos%tddft%cam_alpha, &
                     beta = infos%tddft%cam_beta, mu = infos%tddft%cam_mu)
    call gN%clean()

    dmax = maxval(abs(deN - deP))
    gmax = maxval(abs(deP))
    write(iw,'(/5X,"=== Phase 11 NAC amplitude self-test (I=J) ===")')
    write(iw,'(5X,"production 2e-grad max |de|      = ",ES20.12)') gmax
    write(iw,'(5X,"max |de_nac(I=J) - de_prod|      = ",ES20.12)') dmax
    write(*, '(/5X,"=== Phase 11 NAC amplitude self-test (I=J) ===")')
    write(*, '(5X,"production 2e-grad max |de|      = ",ES20.12)') gmax
    write(*, '(5X,"max |de_nac(I=J) - de_prod|      = ",ES20.12)') dmax

    deallocate(dA, pA, dB, pB, spc, deP, deN)

  end subroutine mrsf_nac_amp_selftest

!###############################################################################

end module tdhf_mrsf_gradient_mod
