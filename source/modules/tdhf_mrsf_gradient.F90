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
    ! Omit the state-independent d2*d2 reference gradient.  The NAC amplitude
    ! kernel needs G(I,J)-G(SCF), and grd2 is linear in its four-index density,
    ! so forming that difference before the integral sweep avoids a second,
    ! identical derivative-integral pass.
    logical :: subtract_reference = .false.
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
    integer(c_int32_t) :: gtag_id
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
    have_custom = infos%dat%contains(tags_gamma, gtag_id)
    if (have_custom) call tagarray_get_data(infos%dat, OQP_nac_gamma, gam_tlf)

    call infos%dat%erase((/ character(len=80) :: OQP_nac_overlap, &
                                                          OQP_nac_trden /))
    call tagarray_reserve_data(infos%dat, OQP_nac_overlap, ta_type_real64, &
          3*natom*nstate*nstate, (/ 3*natom, nstate, nstate /))
    call tagarray_get_data(infos%dat, OQP_nac_overlap, nac_ov)
    nac_ov = 0.0_dp
    ! also export the MO interstate transition densities for diagnostics
    call tagarray_reserve_data(infos%dat, OQP_nac_trden, ta_type_real64, &
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

    ! optional export of the raw AO derivative-overlap matrices for the
    ! Python-side Lagrangian assembly gates (env NAC_DUMP_DS).
    ! NOTE: dSket/dSfull are in the UNNORMALIZED basis convention; the
    ! consumer must apply bfnrm(mu)*bfnrm(nu) exactly as done below.
    block
      character(len=16) :: ev_ds
      real(kind=dp), pointer :: dsk_out(:,:), dsf_out(:,:)
      call get_environment_variable('NAC_DUMP_DS', ev_ds)
      if (len_trim(ev_ds) > 0) then
        call infos%dat%erase((/ character(len=80) :: &
          'OQP::dbg_dsket', 'OQP::dbg_dsfull' /))
        call tagarray_reserve_data(infos%dat, 'OQP::dbg_dsket', ta_type_real64, &
              nbf*nbf*3*natom, (/ nbf*nbf, 3*natom /))
        call tagarray_reserve_data(infos%dat, 'OQP::dbg_dsfull', ta_type_real64, &
              nbf*nbf*3*natom, (/ nbf*nbf, 3*natom /))
        call tagarray_get_data(infos%dat, 'OQP::dbg_dsket', dsk_out)
        call tagarray_get_data(infos%dat, 'OQP::dbg_dsfull', dsf_out)
        do a = 1, natom
          do c = 1, 3
            do nu = 1, nbf
              do mu = 1, nbf
                dsk_out(mu+(nu-1)*nbf, (a-1)*3+c) = &
                  dSket(mu,nu,c,a)*basis%bfnrm(mu)*basis%bfnrm(nu)
                dsf_out(mu+(nu-1)*nbf, (a-1)*3+c) = &
                  dSfull(mu,nu,c,a)*basis%bfnrm(mu)*basis%bfnrm(nu)
              end do
            end do
          end do
        end do
      end if
    end block

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
            if (this%subtract_reference) then
              df1 = this%p2(i1,j1,1)*this%d2(k1,l1,1) &
                  + this%d2(i1,j1,1)*this%p2(k1,l1,1)
            else
              df1 = (this%d2(i1,j1,1)+this%p2(i1,j1,1))*this%d2(k1,l1,1) &
                  +  this%d2(i1,j1,1)                  *this%p2(k1,l1,1)
            end if
            df1 = df1 * coulfact

            if (xcfact /= 0.0_dp .or. xcfact2 /= 0.0_dp) then
              if (this%subtract_reference) then
                dq1 = this%p2(i1,k1,1)*this%d2(j1,l1,1) &
                    + this%d2(i1,k1,1)*this%p2(j1,l1,1) &
                    + this%p2(i1,l1,1)*this%d2(j1,k1,1) &
                    + this%d2(i1,l1,1)*this%p2(j1,k1,1) &
                    + this%p2(i1,k1,2)*this%d2(j1,l1,2) &
                    + this%d2(i1,k1,2)*this%p2(j1,l1,2) &
                    + this%p2(i1,l1,2)*this%d2(j1,k1,2) &
                    + this%d2(i1,l1,2)*this%p2(j1,k1,2)
              else
                dq1 = (this%d2(i1,k1,1)+this%p2(i1,k1,1))*this%d2(j1,l1,1) &
                    +  this%d2(i1,k1,1)                  *this%p2(j1,l1,1) &
                    + (this%d2(i1,l1,1)+this%p2(i1,l1,1))*this%d2(j1,k1,1) &
                    +  this%d2(i1,l1,1)                  *this%p2(j1,k1,1) &
                    + (this%d2(i1,k1,2)+this%p2(i1,k1,2))*this%d2(j1,l1,2) &
                    +  this%d2(i1,k1,2)                  *this%p2(j1,l1,2) &
                    + (this%d2(i1,l1,2)+this%p2(i1,l1,2))*this%d2(j1,k1,2) &
                    +  this%d2(i1,l1,2)                  *this%p2(j1,k1,2)
              end if
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
!  Phase 11: NAC amplitude-term driver (BP#1, explicit-ERI part)
!###############################################################################

!> @brief Build the seven per-state MRSF transition/amplitude channel densities
!>        spc(1:7,:,:) for state `ist` from the response amplitude bvec_mo(:,ist),
!>        EXACTLY mirroring the production build of OQP::td_mrsf_density
!>        (tdhf_mrsf_z_vector.F90:1675-1684):
!>          channels 1-6:  iatogen(bvec_mo(:,ist)) -> mrsfcbc(mo_a,mo_a,...)
!>          channel  7  :  ball := sfdmat( mrsfxvec(bvec_mo(:,ist)) )  [td_abxc]
!>        The triplet flip / spc rescale are NOT applied here (production applies
!>        them to the int2 OUTPUT Fock fmrst2, not to these INPUT densities); the
!>        grd2 get_density consumes the raw channels and applies sgnk/spcscale.
  subroutine mrsf_nac_amp_channels(infos, mo_a, bvec_mo, ist, spc)
    use precision, only: dp
    use types, only: information
    use tdhf_mrsf_lib, only: mrsfcbc, mrsfxvec
    use tdhf_sf_lib, only: sfdmat
    use tdhf_lib, only: iatogen
    use messages, only: show_message, with_abort

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), dimension(:,:) :: mo_a
    real(kind=dp), intent(in), dimension(:,:) :: bvec_mo
    integer, intent(in) :: ist
    real(kind=dp), intent(out), dimension(:,:,:) :: spc   ! (7,nbf,nbf)

    integer :: nbf, nocca, noccb, nbf_tri, ok
    real(kind=dp), allocatable :: wrk1(:,:), bvec_d(:), td_abxc(:,:), &
                                  ta(:), tb(:)

    nbf   = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    nbf_tri = nbf*(nbf+1)/2

    allocate(wrk1(nbf,nbf), td_abxc(nbf,nbf), &
             bvec_d(nocca*(nbf-noccb)), ta(nbf_tri), tb(nbf_tri), &
             source=0.0_dp, stat=ok)
    if (ok/=0) call show_message('Cannot allocate memory', with_abort)

    spc = 0.0_dp

    ! channels 1-6 (+ a mrsfcbc ball that channel 7 overwrites)
    call iatogen(bvec_mo(:,ist), wrk1, nocca, noccb)
    call mrsfcbc(infos, mo_a, mo_a, wrk1, spc(1:7,:,:))

    ! channel 7 (ball) = td_abxc from the mrsfxvec-folded amplitude (sfdmat)
    call mrsfxvec(infos, bvec_mo(:,ist), bvec_d)
    call sfdmat(bvec_d, td_abxc, mo_a, ta, tb, nocca, noccb)
    spc(7,:,:) = td_abxc

    deallocate(wrk1, td_abxc, bvec_d, ta, tb)
  end subroutine mrsf_nac_amp_channels

!###############################################################################

  subroutine mrsf_nac_amp_C(c_handle) bind(C, name="mrsf_nac_amp")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_amp(inf)
  end subroutine mrsf_nac_amp_C

!> Compatibility C entry for the explicit-ERI amplitude term of one ordered
!> state pair.  The resident production driver calls mrsf_nac_amp directly.
  subroutine mrsf_nac_amp_pair_C(c_handle, istate, jstate) &
      bind(C, name="mrsf_nac_amp_pair")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use, intrinsic :: iso_c_binding, only: c_int32_t
    type(oqp_handle_t) :: c_handle
    integer(c_int32_t), intent(in), value :: istate, jstate
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_amp(inf, int(istate), int(jstate))
  end subroutine mrsf_nac_amp_pair_C

!> @brief Phase 11 BP#1: the explicit-ERI part of the analytic amplitude term
!>        G_IJ = X_I^T (d_x A_2e) X_J for all MRSF state pairs (I/=J), stored as
!>        OQP::nac_amp (3, natom, nstate, nstate).  A caller may select one
!>        ordered pair; the resident driver supplies its gap-scaled left
!>        response y_IJ=X_I/(Omega_J-Omega_I) before calling this routine.
!>
!>        Mechanism: the validated bilinear engine grd2_mrsf_nac_compute_data_t
!>        is passed to grd2 once per requested pair with subtract_reference
!>        enabled.  Its density callback omits the state-independent d2*d2
!>        contribution before the derivative-integral sweep.  By linearity this
!>        is G_full-G_SCF, but it avoids both a second grd2 traversal and the
!>        cancellation that would result from subtracting two gradients after
!>        integration.
!>
!>        OQP::td_p must not be consumed here.  It is unlabelled scratch owned
!>        by the diagonal gradient Z-vector and can therefore contain the last
!>        requested root when NAC follows a gradient.  The interstate
!>        difference-density/Fock skeleton is evaluated by mrsf_nac_esum;
!>        reusing td_p here both makes the result call-history dependent and
!>        contaminates every pair with one diagonal state's relaxed density.
  subroutine mrsf_nac_amp(infos, only_istate, only_jstate)
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use mathlib, only: unpack_matrix

    implicit none

    character(len=*), parameter :: subroutine_name = "mrsf_nac_amp"
    character(len=*), parameter :: OQP_nac_amp = "OQP::nac_amp"

    type(information), target, intent(inout) :: infos
    integer, intent(in), optional :: only_istate, only_jstate
    type(basis_set), pointer :: basis

    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), bvec_mo(:,:)
    real(kind=dp), pointer :: nac_amp(:,:,:)
    character(len=*), parameter :: tags_required(4) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_VEC_MO_A, OQP_td_bvec_mo /)

    real(kind=dp), allocatable, target :: d0(:,:,:), pIJ(:,:,:), &
         spcI(:,:,:), spcJ(:,:,:), dcopy(:,:,:)
    real(kind=dp), allocatable, target :: dmat_a_owned(:), &
         dmat_b_owned(:), mo_a_owned(:,:), bvec_mo_owned(:,:)
    real(kind=dp), allocatable :: deFull(:,:)
    class(grd2_compute_data_t), allocatable :: gFull
    real(kind=dp) :: scale_exch, scale_exch2
    logical :: dft, do_cam
    integer :: nbf, mrst, natom, nstate, ist, jst, c, a, ok
    integer :: ist_first, ist_last, jst_first, jst_last

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf    = basis%nbf
    natom  = ubound(infos%atoms%grad, 2)
    nstate = infos%tddft%nstate
    mrst   = infos%tddft%mult

    if (present(only_istate) .neqv. present(only_jstate)) then
      call show_message('mrsf_nac_amp pair selection requires both states', &
                        with_abort)
    end if
    ist_first = 1
    ist_last = nstate
    jst_first = 1
    jst_last = nstate
    if (present(only_istate)) then
      if (only_istate < 1 .or. only_istate > nstate .or. &
          only_jstate < 1 .or. only_jstate > nstate .or. &
          only_istate == only_jstate) then
        call show_message('mrsf_nac_amp ordered pair is invalid', with_abort)
      end if
      ist_first = only_istate
      ist_last = only_istate
      jst_first = only_jstate
      jst_last = only_jstate
    end if

    if (mrst/=1 .and. mrst/=3) &
      call show_message('mrsf_nac_amp supports mult=1,3 only', with_abort)

    dft = infos%control%hamilton == 20
    scale_exch  = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
      scale_exch  = infos%dft%HFscale
      scale_exch2 = infos%tddft%HFscale
    end if
    do_cam = dft .and. infos%dft%cam_flag

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)

    ! Any output reserve may move TagArray storage and invalidate every cached
    ! input pointer.  Keep owned production inputs across the whole amplitude
    ! kernel, including nested gradient routines that may manage records.
    allocate(dmat_a_owned(size(dmat_a)), dmat_b_owned(size(dmat_b)), &
             mo_a_owned(size(mo_a,1),size(mo_a,2)), &
             bvec_mo_owned(size(bvec_mo,1),size(bvec_mo,2)), stat=ok)
    if (ok/=0) call show_message('cannot allocate memory', WITH_ABORT)
    dmat_a_owned = dmat_a
    dmat_b_owned = dmat_b
    mo_a_owned = mo_a
    bvec_mo_owned = bvec_mo
    dmat_a => dmat_a_owned
    dmat_b => dmat_b_owned
    mo_a => mo_a_owned
    bvec_mo => bvec_mo_owned

    call infos%dat%erase((/ character(len=80) :: OQP_nac_amp /))
    call tagarray_reserve_data(infos%dat, OQP_nac_amp, ta_type_real64, &
          3*natom*nstate*nstate, (/ 3*natom, nstate, nstate /))
    call tagarray_get_data(infos%dat, OQP_nac_amp, nac_amp)
    nac_amp = 0.0_dp

    allocate(d0(nbf,nbf,2), dcopy(nbf,nbf,2), pIJ(nbf,nbf,2), &
             spcI(7,nbf,nbf), spcJ(7,nbf,nbf), deFull(3,natom), &
             source=0.0_dp, stat=ok)
    if (ok/=0) call show_message('cannot allocate memory', WITH_ABORT)

    call unpack_matrix(dmat_a, d0(:,:,1))
    call unpack_matrix(dmat_b, d0(:,:,2))

    ! This engine owns its p2 input and keeps it zero.  In particular, never
    ! infer pair identity or freshness from the generic OQP::td_p record.

    do ist = ist_first, ist_last
      do jst = jst_first, jst_last
        if (ist == jst) cycle

        ! per-state channel densities (mirror production td_mrsf_density build)
        call mrsf_nac_amp_channels(infos, mo_a, bvec_mo, ist, spcI)
        call mrsf_nac_amp_channels(infos, mo_a, bvec_mo, jst, spcJ)

        ! Direct pair channel density with an owned zero p2.  init mutates
        ! d2/p2 in place, so retain a fresh copy for every requested pair.
        dcopy = d0
        deFull = 0.0_dp
        gFull = grd2_mrsf_nac_compute_data_t( d2 = dcopy, p2 = pIJ, &
                     spcI = spcI, spcJ = spcJ, nbf = nbf, &
                     subtract_reference = .true., &
                     hfscale = scale_exch, hfscale2 = scale_exch2, &
                     spcscale = [infos%tddft%spc_coco, &
                                 infos%tddft%spc_ovov, &
                                 infos%tddft%spc_coov], mrst = mrst )
        call gFull%init()
        call grd2_driver(infos, basis, deFull, gFull, &
                         cam = do_cam, alpha = infos%tddft%cam_alpha, &
                         beta = infos%tddft%cam_beta, mu = infos%tddft%cam_mu)
        call gFull%clean()

        ! gFull omits the d2*d2 reference density before the integral sweep,
        ! hence it is directly G_full-G_SCF without a cancellation-prone
        ! second derivative-integral pass.
        do a = 1, natom
          do c = 1, 3
            nac_amp((a-1)*3+c, ist, jst) = deFull(c,a)
          end do
        end do
      end do
    end do

    deallocate(d0, dcopy, pIJ, spcI, spcJ, deFull, &
               dmat_a_owned, dmat_b_owned, mo_a_owned, bvec_mo_owned)
  end subroutine mrsf_nac_amp

!###############################################################################

!> @brief CLOSED-FORM interstate `esum` term: the explicit 1e/Fock part of the
!>        matvec derivative X_I^T (dA/dR) X_J.
!>
!>   The MRSF matvec (tdhf_mrsf_energy.F90) builds A.x as
!>     amo = mrsfmntoia(2e generalized Fock)      -> differentiates to `ana2e`
!>         + mrsfesum(iatogen(x), fa, fb, amo)    -> the Fock part, THIS term
!>   with fa = mo_a^T FOCK_A mo_a, fb = mo_b^T FOCK_B mo_b. mrsfesum performs the
!>   TDA Fock contraction (F_vv X - X F_oo), so the coupling's Fock part is
!>     X_I^T mrsfesum(iatogen(X_J), fa, fb).
!>   Differentiating at FIXED MOs and FIXED (transported) amplitudes gives
!>     esum^x = Tr(fa^x . Gam_A) + Tr(fb^x . Gam_B)
!>            = Tr(P^IJ_a . dFOCK_A^skel/dR) + Tr(P^IJ_b . dFOCK_B^skel/dR)
!>   where Gam_A = tij (alpha occ-occ), Gam_B = tab (beta virt-virt) are exactly
!>   the interstate difference-density blocks from mrsf_interstate_tden, and
!>     P^IJ_a = mo_a(:,1:nocca) tij mo_a(:,1:nocca)^T
!>     P^IJ_b = mo_b(:,noccb+1:) tab mo_b(:,noccb+1:)^T.
!>
!>   CRITICAL (this is the p20 bug): dFOCK^skel is the PURE SKELETON Fock
!>   derivative at FROZEN reference density (DM_A/DM_B) --
!>     FOCK_A = h + J[Pa+Pb] - c_x K[Pa],  FOCK_B = h + J[Pa+Pb] - c_x K[Pb].
!>   p20_esum_grad.py instead injected P^IJ as td_p and ran the FULL gradient
!>   seam, which returns Tr(P^IJ dF/dR) PLUS P^IJ's own 2e response dG[P^IJ]
!>   PLUS the W.dS^x overlap term -> ~20x too large. Here the 1e piece uses only
!>   kinetic + nuclear-attraction (Hellmann-Feynman + Pulay) and DELIBERATELY
!>   omits grad_ee_overlap (that W.dS^x term belongs to d_ov, not to esum), and
!>   the 2e piece uses fock_deriv_contract_os at frozen densities.
!>
!>   NOTE: for a DFT reference the Fock also carries V_xc, whose derivative is
!>   NOT included here; validate against an HF reference first.
!>   Result -> OQP::nac_esum (3,natom).
  subroutine mrsf_nac_esum_C(c_handle, istate, jstate) &
      bind(C, name="mrsf_nac_esum")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use io_constants, only: iw
    use types, only: information
    use, intrinsic :: iso_c_binding, only: c_int64_t
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t), value :: istate, jstate
    type(information), pointer :: inf
    logical :: log_was_open

    inf => oqp_handle_get_info(c_handle)
    inquire(unit=iw, opened=log_was_open)
    if (.not. log_was_open) &
      open(unit=iw, file=inf%log_filename, position='append')
    call mrsf_nac_esum(inf, int(istate), int(jstate))
    if (.not. log_was_open) close(iw)
  end subroutine mrsf_nac_esum_C

  subroutine mrsf_nac_esum(infos, istate, jstate)
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use mathlib, only: unpack_matrix, pack_matrix, orthogonal_transform_sym
    use tdhf_mrsf_lib, only: mrsf_interstate_tden
    use fock_deriv_mod, only: fock_deriv_contract_os
    use grd1, only: grad_ee_kinetic, grad_en_hellman_feynman, grad_en_pulay, &
                    grad_ee_overlap
    use mod_dft, only: dft_initialize, dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient

    implicit none
    character(len=*), parameter :: subroutine_name = "mrsf_nac_esum"
    character(len=*), parameter :: OQP_nac_esum = "OQP::nac_esum"

    type(information), target, intent(inout) :: infos
    integer, intent(in) :: istate, jstate
    type(basis_set), pointer :: basis

    real(kind=dp), contiguous, pointer :: bvec_mo(:,:), mo_a(:,:), mo_b(:,:), &
                                          dmat_a(:), dmat_b(:)
    real(kind=dp), pointer :: out(:,:)
    real(kind=dp), allocatable :: tij(:,:), tab(:,:), pij_a(:,:), pij_b(:,:), &
                                  pa(:,:), pb(:,:), ptot(:,:), scr(:,:), &
                                  p1e(:), gx(:,:), g1(:,:)
    real(kind=dp), allocatable, target :: bvec_mo_owned(:,:), &
      mo_a_owned(:,:), mo_b_owned(:,:), dmat_a_owned(:), dmat_b_owned(:)
    real(kind=dp) :: hfscale, tol
    integer :: nbf, nbf2, nocca, noccb, nvirb, natom, ok
    character(len=80) :: tags_req(5)
    type(dft_grid_t) :: molGrid
    logical :: dft
    real(kind=dp), allocatable :: xcd(:,:,:), xcp(:,:,:), gxc(:,:)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb
    natom = ubound(infos%atoms%grad, 2)
    tol = log(10.0_dp)*infos%control%int2e_cutoff

    hfscale = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale = infos%dft%hfscale
    dft = infos%control%hamilton == 20

    tags_req(1) = OQP_td_bvec_mo
    tags_req(2) = OQP_VEC_MO_A
    tags_req(3) = OQP_VEC_MO_B
    tags_req(4) = OQP_DM_A
    tags_req(5) = OQP_DM_B
    call data_has_tags(infos%dat, tags_req, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)

    ! Pair/diagnostic output records are reserved repeatedly below.  Own every
    ! resident input that remains live after the first reserve so allocator
    ! movement cannot turn a valid calculation into a dangling-pointer read.
    allocate(bvec_mo_owned(size(bvec_mo,1),size(bvec_mo,2)), &
             mo_a_owned(size(mo_a,1),size(mo_a,2)), &
             mo_b_owned(size(mo_b,1),size(mo_b,2)), &
             dmat_a_owned(size(dmat_a)), dmat_b_owned(size(dmat_b)), stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)
    bvec_mo_owned = bvec_mo
    mo_a_owned = mo_a
    mo_b_owned = mo_b
    dmat_a_owned = dmat_a
    dmat_b_owned = dmat_b
    bvec_mo => bvec_mo_owned
    mo_a => mo_a_owned
    mo_b => mo_b_owned
    dmat_a => dmat_a_owned
    dmat_b => dmat_b_owned

    allocate(tij(nocca,nocca), tab(nvirb,nvirb), &
             pij_a(nbf,nbf), pij_b(nbf,nbf), &
             pa(nbf,nbf), pb(nbf,nbf), ptot(nbf,nbf), scr(nbf,nbf), &
             p1e(nbf2), gx(3,natom), g1(3,natom), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)

  ! (1) interstate difference-density blocks (SOMO-unfolded, symmetric)
    call mrsf_interstate_tden(infos, bvec_mo, istate, jstate, tij, tab)

  ! (2) to AO: P^IJ_a = mo_a(:,1:nocca) tij mo_a(:,1:nocca)^T
    call dgemm('n','n', nbf, nocca, nocca, &
               1.0_dp, mo_a, nbf, tij, nocca, 0.0_dp, scr, nbf)
    call dgemm('n','t', nbf, nbf, nocca, &
               1.0_dp, scr, nbf, mo_a, nbf, 0.0_dp, pij_a, nbf)
  !     P^IJ_b = mo_b(:,noccb+1:nbf) tab mo_b(:,noccb+1:nbf)^T
    scr = 0.0_dp
    call dgemm('n','n', nbf, nvirb, nvirb, &
               1.0_dp, mo_b(1,noccb+1), nbf, tab, nvirb, 0.0_dp, scr, nbf)
    call dgemm('n','t', nbf, nbf, nvirb, &
               1.0_dp, scr, nbf, mo_b(1,noccb+1), nbf, 0.0_dp, pij_b, nbf)

  ! Formal pair outputs used by the native response assembly.  These are part
  ! of the production ABI; production must not depend on env-gated dbg_* tags.
    block
      real(kind=dp), pointer :: pa_pair(:), pb_pair(:), pa_packed(:), pb_packed(:)
      call infos%dat%erase((/ character(len=80) :: &
        'OQP::nac_pij_a', 'OQP::nac_pij_b', &
        'OQP::nac_dm1_a', 'OQP::nac_dm1_b' /))
      call tagarray_reserve_data(infos%dat, 'OQP::nac_pij_a', ta_type_real64, &
            nbf*nbf, (/ nbf*nbf /), comment='ordered MRSF alpha AO pair density')
      call tagarray_reserve_data(infos%dat, 'OQP::nac_pij_b', ta_type_real64, &
            nbf*nbf, (/ nbf*nbf /), comment='ordered MRSF beta AO pair density')
      call tagarray_reserve_data(infos%dat, 'OQP::nac_dm1_a', ta_type_real64, &
            nbf2, (/ nbf2 /), comment='packed ordered MRSF alpha AO pair density')
      call tagarray_reserve_data(infos%dat, 'OQP::nac_dm1_b', ta_type_real64, &
            nbf2, (/ nbf2 /), comment='packed ordered MRSF beta AO pair density')
      call tagarray_get_data(infos%dat, 'OQP::nac_pij_a', pa_pair)
      call tagarray_get_data(infos%dat, 'OQP::nac_pij_b', pb_pair)
      call tagarray_get_data(infos%dat, 'OQP::nac_dm1_a', pa_packed)
      call tagarray_get_data(infos%dat, 'OQP::nac_dm1_b', pb_packed)
      pa_pair = reshape(pij_a, (/ nbf*nbf /))
      pb_pair = reshape(pij_b, (/ nbf*nbf /))
      call pack_matrix(pij_a, pa_packed)
      call pack_matrix(pij_b, pb_packed)
    end block

  ! Optional backward-compatible diagnostic aliases for forensic scripts.
    block
      character(len=16) :: ev_p
      real(kind=dp), pointer :: pa_out(:), pb_out(:)
      call get_environment_variable('NAC_DUMP_PIJ', ev_p)
      if (len_trim(ev_p) > 0) then
        call infos%dat%erase((/ character(len=80) :: &
          'OQP::dbg_pij_a', 'OQP::dbg_pij_b' /))
        call tagarray_reserve_data(infos%dat, 'OQP::dbg_pij_a', ta_type_real64, &
              nbf*nbf, (/ nbf*nbf /))
        call tagarray_reserve_data(infos%dat, 'OQP::dbg_pij_b', ta_type_real64, &
              nbf*nbf, (/ nbf*nbf /))
        call tagarray_get_data(infos%dat, 'OQP::dbg_pij_a', pa_out)
        call tagarray_get_data(infos%dat, 'OQP::dbg_pij_b', pb_out)
        pa_out = reshape(pij_a, (/ nbf*nbf /))
        pb_out = reshape(pij_b, (/ nbf*nbf /))
      end if
    end block

  ! (3) explicit 1e: Tr[(P^IJ_a + P^IJ_b) . dh/dR]
  !     kinetic + nuclear attraction ONLY -- grad_ee_overlap is EXCLUDED on
  !     purpose (the W.dS^x term is part of d_ov, and including it was one of
  !     the contaminations that made the p20 esum ~20x too big).
    scr = pij_a + pij_b
    call pack_matrix(scr, p1e)
    call grad_ee_kinetic(basis, p1e, gx)
    call grad_en_hellman_feynman(basis, infos%atoms%xyz, infos%atoms%zn, p1e, gx)
    call grad_en_pulay(basis, infos%atoms%xyz, infos%atoms%zn, p1e, gx)

  ! diagnostic split: export the 1e part alone so the 1e/2e balance can be checked
    block
      real(kind=dp), pointer :: o1(:,:)
      call infos%dat%erase((/ character(len=80) :: 'OQP::nac_esum_1e' /))
      call tagarray_reserve_data(infos%dat, 'OQP::nac_esum_1e', ta_type_real64, 3*natom, (/ 3, natom /))
      call tagarray_get_data(infos%dat, 'OQP::nac_esum_1e', o1)
      o1 = gx
    end block
  ! (4) explicit 2e at FROZEN reference density:
  !     Tr[P^IJ_a . dG^a[Pa,Pb]] + Tr[P^IJ_b . dG^b[Pa,Pb]]
    call unpack_matrix(dmat_a, pa)
    call unpack_matrix(dmat_b, pb)
    ptot = pa + pb
    g1 = 0.0_dp
    call fock_deriv_contract_os(infos, basis, ptot, pa, pij_a, hfscale, g1)
    gx = gx + g1
    g1 = 0.0_dp
    call fock_deriv_contract_os(infos, basis, ptot, pb, pij_b, hfscale, g1)
    gx = gx + g1

  ! (4b) explicit XC: Tr[P^IJ_a . dV_xc^a/dR] + Tr[P^IJ_b . dV_xc^b/dR] at the
  !      FROZEN reference density. This is the same production routine the
  !      validated MRSF gradient uses for its relaxed-density XC term, with the
  !      relaxed density replaced by the interstate probe P^IJ. No xa/xb is
  !      passed, so only the grad_v_xc branch runs (no f_xc kernel term, which
  !      belongs to the transition-density/2e side, not to the Fock part).
  !      NOTE two landmines: utddft_xc_gradient ACCUMULATES into dedft (start
  !      the scratch at zero) and its da/db/pa/pb are
  !      `intent(inout)` and get SCALED BY BASIS NORMS IN PLACE -> pass COPIES,
  !      otherwise pij_a/pij_b are silently corrupted for anything downstream.
    if (dft) then
      allocate(xcd(nbf,nbf,2), xcp(nbf,nbf,2), gxc(3,natom), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('Cannot allocate memory', with_abort)
      ! This is an interstate *linear probe*: exclude the probe-independent
      ! ground-state gradient and include the derivative of the normalized
      ! atom-centred fuzzy-cell weights.  The latter is essential here because
      ! the small unrelaxed probe is amplified by the near-degenerate pair-Z
      ! resolvent.  A zero-probe subtraction removes only the ground-state
      ! constant; it cannot recover the missing (dw/dR)*delta_P e_xc term.
      call dftclean(infos)
      call dft_initialize(infos, basis, molGrid, verbose=.false.)
      xcd(:,:,1) = pa
      xcd(:,:,2) = pb
      xcp(:,:,1) = pij_a
      xcp(:,:,2) = pij_b
      call utddft_xc_gradient(basis=basis, molGrid=molGrid, dedft=gxc, &
           da=xcd(:,:,1), db=xcd(:,:,2), &
           pa=xcp(:,:,1:1), pb=xcp(:,:,2:2), &
           nmtx=1, threshold=0.0d0, infos=infos, &
           include_ground_state=.false., &
           include_weight_derivative=.true.)
      call dftclean(infos)
      gx = gx + gxc
      deallocate(xcd, xcp, gxc)
    end if

  ! (4c) Fock-weighted interstate density term  -Tr[W^IJ . S^x]  -> OQP::nac_wsx
  !      (exported SEPARATELY; NOT added to gx, esum stays as validated).
  !      Differentiating the Fock part of the coupling, Tr[Gam . C^T F_AO C]:
  !        d/dR = Tr[Gam C^T dF_AO C]            <- esum (skeleton), validated
  !             + Tr[M . U^x],  M = Gam F_MO + F_MO Gam
  !      and orthonormality fixes the symmetric part U^x_sym = -1/2 S^x_MO, so
  !        resp  ⊃  -Tr[W . S^x],   W = 1/2 C M C^T.
  !      W is FOCK-WEIGHTED so it carries the orbital energies (core ~ -20 Ha) and
  !      is therefore LARGE — this is the piece that cancels the large skeleton
  !      esum. It is exactly the grad_ee_overlap contraction that was excluded
  !      from esum; it belongs to resp, not to d_ov.
    block
      real(kind=dp), allocatable :: fa(:,:), fb(:,:), mA(:,:), mB(:,:), &
                                    wao(:,:), wpk(:), gw(:,:), sc(:,:)
      real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:)
      real(kind=dp), pointer :: owsx(:,:)
      real(kind=dp), allocatable :: fpk(:)
      allocate(fa(nbf,nbf), fb(nbf,nbf), mA(nocca,nocca), mB(nvirb,nvirb), &
               wao(nbf,nbf), wpk(nbf2), gw(3,natom), sc(nbf,nbf), fpk(nbf2), &
               source=0.0_dp)
      call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
      ! F in MO basis (same construction the matvec uses)
      call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, fpk)
      call unpack_matrix(fpk, fa)
      call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, fpk)
      call unpack_matrix(fpk, fb)
      ! M_A = Gam_A F_oo + F_oo Gam_A   (alpha occ block)
      mA = matmul(tij, fa(1:nocca,1:nocca)) + matmul(fa(1:nocca,1:nocca), tij)
      ! M_B = Gam_B F_vv + F_vv Gam_B   (beta virt block)
      mB = matmul(tab, fb(noccb+1:nbf,noccb+1:nbf)) &
         + matmul(fb(noccb+1:nbf,noccb+1:nbf), tab)
      ! W = 1/2 [ C_o M_A C_o^T + C_v M_B C_v^T ]
      call dgemm('n','n', nbf, nocca, nocca, 0.5_dp, mo_a, nbf, mA, nocca, &
                 0.0_dp, sc, nbf)
      call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, sc, nbf, mo_a, nbf, &
                 0.0_dp, wao, nbf)
      sc = 0.0_dp
      call dgemm('n','n', nbf, nvirb, nvirb, 0.5_dp, mo_b(1,noccb+1), nbf, mB, &
                 nvirb, 0.0_dp, sc, nbf)
      call dgemm('n','t', nbf, nbf, nvirb, 1.0_dp, sc, nbf, mo_b(1,noccb+1), nbf, &
                 1.0_dp, wao, nbf)
      call pack_matrix(wao, wpk)
      gw = 0.0_dp
      call grad_ee_overlap(basis, wpk, gw)
      gw = -gw
      call infos%dat%erase((/ character(len=80) :: 'OQP::nac_wsx' /))
      call tagarray_reserve_data(infos%dat, 'OQP::nac_wsx', ta_type_real64, 3*natom, (/ 3, natom /))
      call tagarray_get_data(infos%dat, 'OQP::nac_wsx', owsx)
      owsx = gw
      deallocate(fa, fb, mA, mB, wao, wpk, gw, sc, fpk)
    end block

    call infos%dat%erase((/ character(len=80) :: OQP_nac_esum /))
    call tagarray_reserve_data(infos%dat, OQP_nac_esum, ta_type_real64, 3*natom, (/ 3, natom /))
    call tagarray_get_data(infos%dat, OQP_nac_esum, out)
    out = gx

    deallocate(tij, tab, pij_a, pij_b, pa, pb, ptot, scr, p1e, gx, g1, &
               bvec_mo_owned, mo_a_owned, mo_b_owned, dmat_a_owned, &
               dmat_b_owned)
  end subroutine mrsf_nac_esum

!###############################################################################

!> @brief Phase 12 (closed-form NAC): the interstate matvec derivative
!>        X_i^T (dA/dR) X_j by FORTRAN-LEVEL polarization of the analytic MRSF
!>        gradient. The excited-state gradient G(X) = E0^xi + Omega^xi(X) with
!>        Omega^xi(X) = X^T(dA/dR)X a pure quadratic form in the amplitude (Lee
!>        JCP 150,184111 Eq 3.21: no state eigenvalue in the assembly). Hence
!>          X_i^T(dA)X_j = 1/2[ G(X_i+X_j) - G(X_i) - G(X_j) + G(0) ]
!>        (coeffs sum to 0 so the ground gradient E0^xi cancels). Each G is the
!>        FULL existing z-vector + gradient with the amplitude injected in the
!>        target column IN FORTRAN (no Python set_tdhf_target reload). Result to
!>        OQP::nac_amp_polar (3,natom). Divide by (Om_j-Om_i) in Python for d_amp.
!>        NOTE: requires mrsf_nac_cphf_mode = .false. (standard RHS, not the
!>        overlap-gamma NAC-CPHF override) — caller must NOT set_mrsf_nac_cphf.
  subroutine mrsf_nac_polarize_C(c_handle, istate, jstate) &
      bind(C, name="mrsf_nac_polarize")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use, intrinsic :: iso_c_binding, only: c_int64_t
    type(oqp_handle_t) :: c_handle
    integer(c_int64_t), value :: istate, jstate
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_polarize(inf, int(istate), int(jstate))
  end subroutine mrsf_nac_polarize_C

  subroutine mrsf_nac_polarize(infos, istate, jstate)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use tdhf_mrsf_z_vector_mod, only: tdhf_mrsf_z_vector

    implicit none
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: istate, jstate

    character(len=*), parameter :: OQP_nac_amp_polar = "OQP::nac_amp_polar"
    real(kind=dp), contiguous, pointer :: bvec_mo(:,:)
    real(kind=dp), pointer :: out(:,:)
    real(kind=dp), allocatable :: save_bvec(:,:), Xi(:), Xj(:), accum(:,:)
    real(kind=dp) :: coef(4)
    integer :: natom, save_target, c, nx, nst
    character(len=80) :: tags_req(1)

    tags_req(1) = OQP_td_bvec_mo
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    nx  = size(bvec_mo, 1)
    nst = size(bvec_mo, 2)
    natom = ubound(infos%atoms%grad, 2)

    allocate(save_bvec(nx, nst), Xi(nx), Xj(nx), accum(3, natom))
    save_bvec = bvec_mo
    save_target = infos%tddft%target_state
    Xi = bvec_mo(:, istate)
    Xj = bvec_mo(:, jstate)
    accum = 0.0_dp
    coef = (/ 0.5_dp, -0.5_dp, -0.5_dp, 0.5_dp /)   ! X_i+X_j, X_i, X_j, 0
    infos%tddft%target_state = istate

    ! ---- Homogeneity diagnostic (env NAC_HOMOG) --------------------------
    ! Decisive test of the degree-2 premise behind polarization. Evaluate the
    ! analytic gradient G along scales s = 0,1,2,3 of the istate amplitude and
    ! dump each to OQP::nac_homog(3,natom,4). If G(X)=c+L(X)+X^T(dA)X (degree<=2)
    ! then with D_s = G(sX)-G(0): (D3-3 D1) = 6Q and (D2-2 D1) = 2Q, so the ratio
    ! (D3-3 D1)/(D2-2 D1) == 3 EXACTLY (componentwise). A ratio != 3 proves a
    ! degree>2 / non-polynomial term in the ground-config channel, which is the
    ! ONLY way polarization can be deficient while every operator is bilinear.
    block
      character(len=16) :: ev_homog
      real(kind=dp), pointer :: homog_out(:,:)
      real(kind=dp), allocatable :: homog(:,:,:)
      real(kind=dp) :: scl(4)
      call get_environment_variable('NAC_HOMOG', ev_homog)
      if (len_trim(ev_homog) > 0) then
        allocate(homog(3, natom, 4))
        scl = (/ 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp /)
        do c = 1, 4
          bvec_mo(:, istate) = scl(c) * Xi
          infos%atoms%grad = 0.0_dp
          call tdhf_mrsf_z_vector(infos)
          call tdhf_mrsf_gradient(infos)
          homog(:, :, c) = infos%atoms%grad
        end do
        bvec_mo = save_bvec
        infos%tddft%target_state = save_target
        infos%atoms%grad = 0.0_dp
        call infos%dat%erase((/ character(len=80) :: 'OQP::nac_homog' /))
        call tagarray_reserve_data(infos%dat, 'OQP::nac_homog', ta_type_real64, &
             3*natom*4, (/ 3*natom, 4 /))
        call tagarray_get_data(infos%dat, 'OQP::nac_homog', homog_out)
        homog_out = reshape(homog, (/ 3*natom, 4 /))
        write(iw,'(/5X,"=== NAC homogeneity dump (scales 0,1,2,3) done ===")')
        deallocate(homog, save_bvec, Xi, Xj, accum)
        return
      end if
    end block

    do c = 1, 4
      select case (c)
        case (1); bvec_mo(:, istate) = Xi + Xj
        case (2); bvec_mo(:, istate) = Xi
        case (3); bvec_mo(:, istate) = Xj
        case (4); bvec_mo(:, istate) = 0.0_dp
      end select
      infos%atoms%grad = 0.0_dp
      call tdhf_mrsf_z_vector(infos)
      call tdhf_mrsf_gradient(infos)
      accum = accum + coef(c) * infos%atoms%grad
    end do

    ! restore
    bvec_mo = save_bvec
    infos%tddft%target_state = save_target
    infos%atoms%grad = 0.0_dp

    call infos%dat%erase((/ character(len=80) :: OQP_nac_amp_polar /))
    call tagarray_reserve_data(infos%dat, OQP_nac_amp_polar, ta_type_real64, &
         3*natom, (/ 3, natom /))
    call tagarray_get_data(infos%dat, OQP_nac_amp_polar, out)
    out = accum

    write(iw,'(/5X,"=== NAC polarization X_",I0,"^T dA X_",I0," computed ===")') &
         istate, jstate
    deallocate(save_bvec, Xi, Xj, accum)
  end subroutine mrsf_nac_polarize

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

!###############################################################################
!> ROUTE A resident closed-form engine.  It evaluates
!>   M_pq = d[ytil^T A(C) X_J]/d theta_pq at frozen AO Fock
!> using the exact mrsfcbc/mrsfmntoia adjoint.  The two response vectors share
!> one batched ERI build; the frozen-Fock contribution is evaluated directly
!> in MO space, including the matvec's ixcore diagonal-overwrite semantics.
!> OQP::nac_ytil and OQP::nac_xstate are input folded vectors;
!> OQP::nac_mt_frozen is the column-major output matrix.
  subroutine mrsf_nac_wpair_C(c_handle, istate, jstate) &
      bind(C, name="mrsf_nac_wpair")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use, intrinsic :: iso_c_binding, only: c_int32_t
    type(oqp_handle_t) :: c_handle
    integer(c_int32_t), intent(in), value :: istate, jstate
    type(information), pointer :: inf
    interface
      subroutine mrsf_nac_wpair_impl(infos, istate, jstate)
        use types, only: information
        type(information), target, intent(inout) :: infos
        integer, intent(in) :: istate, jstate
      end subroutine mrsf_nac_wpair_impl
    end interface
    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_wpair_impl(inf, int(istate), int(jstate))
  end subroutine mrsf_nac_wpair_C

  subroutine mrsf_nac_wpair_impl(infos, istate, jstate)
    use oqp_tagarray_driver
    use types, only: information
    use precision, only: dp
    use messages, only: show_message, with_abort
    implicit none
    character(len=*), parameter :: module_name = "tdhf_mrsf_gradient_mod"
    character(len=*), parameter :: subroutine_name = "mrsf_nac_wpair_impl"
    character(len=*), parameter :: OQP_nac_ytil = "OQP::nac_ytil"
    character(len=*), parameter :: OQP_nac_xstate = "OQP::nac_xstate"
    character(len=*), parameter :: OQP_nac_mt = "OQP::nac_mt_frozen"
    character(len=*), parameter :: tags_required(2) = (/ character(len=80) :: &
      OQP_nac_ytil, OQP_nac_xstate /)
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: istate, jstate
    real(kind=dp), contiguous, pointer :: ytil(:), xstate(:)
    real(kind=dp), pointer :: mt_out(:)
    real(kind=dp), allocatable :: ytil_batch(:,:), xstate_batch(:,:), &
                                  mt_batch(:,:,:)
    integer :: nbf, xdim
    interface
      subroutine mrsf_nac_wpair_batch_impl(infos, ytil_batch, &
                                           xstate_batch, mt_batch)
        use types, only: information
        use precision, only: dp
        type(information), target, intent(inout) :: infos
        real(kind=dp), contiguous, intent(in) :: ytil_batch(:,:), &
                                                 xstate_batch(:,:)
        real(kind=dp), contiguous, intent(out) :: mt_batch(:,:,:)
      end subroutine mrsf_nac_wpair_batch_impl
    end interface

    if (istate == jstate) return
    nbf = infos%basis%nbf
    xdim = infos%mol_prop%nelec_a*(nbf-infos%mol_prop%nelec_b)
    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, &
                       with_abort)
    call tagarray_get_data(infos%dat, OQP_nac_ytil, ytil)
    call tagarray_get_data(infos%dat, OQP_nac_xstate, xstate)
    if (size(ytil) /= xdim .or. size(xstate) /= xdim) &
      call show_message('Invalid scalar MRSF NAC wpair vectors.', with_abort)
    allocate(ytil_batch(xdim,1), xstate_batch(xdim,1), &
             mt_batch(nbf,nbf,1))
    ytil_batch(:,1) = ytil
    xstate_batch(:,1) = xstate
    call mrsf_nac_wpair_batch_impl(infos, ytil_batch, xstate_batch, mt_batch)

    call infos%dat%erase((/ character(len=80) :: OQP_nac_mt /))
    call tagarray_reserve_data(infos%dat, OQP_nac_mt, ta_type_real64, nbf*nbf, &
                                (/ nbf*nbf /))
    call tagarray_get_data(infos%dat, OQP_nac_mt, mt_out)
    mt_out = reshape(mt_batch(:,:,1), (/ nbf*nbf /))
  end subroutine mrsf_nac_wpair_impl

!###############################################################################

!> Batch the frozen pair source for at most three canonical state pairs.  Each
!> pair contributes its folded X and Y density, so one integral traversal
!> evaluates 2*nrhs Fock-like sources.  mt_batch is kept outside TagArray;
!> the resident driver publishes one slice at the original pair-consumption
!> point, preserving the legacy record lifetime and metric-column streaming.
  subroutine mrsf_nac_wpair_batch_impl(infos, ytil_batch, xstate_batch, mt_batch)
    use oqp_tagarray_driver
    use types, only: information
    use precision, only: dp
    use messages, only: show_message, with_abort
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_mrsf_data_t, mrsfcbc, mrsfxvec, mrsfsp
    use tdhf_lib, only: iatogen
    use mathlib, only: orthogonal_transform_sym, unpack_matrix
    use oqp_linalg
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_int
    implicit none
    character(len=*), parameter :: module_name = "tdhf_mrsf_gradient_mod"
    character(len=*), parameter :: subroutine_name = &
      "mrsf_nac_wpair_batch_impl"
    integer, parameter :: max_batch_width = 3
    character(len=*), parameter :: tags_required(4) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_FOCK_B, OQP_VEC_MO_A, OQP_VEC_MO_B /)
    type(information), target, intent(inout) :: infos
    real(kind=dp), contiguous, intent(in) :: ytil_batch(:,:), xstate_batch(:,:)
    real(kind=dp), contiguous, intent(out) :: mt_batch(:,:,:)
    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(int2_mrsf_data_t), target :: int2_data_st
    real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:), &
                                          mo_a(:,:), mo_b(:,:)
    real(kind=dp), allocatable, target :: mrsf_density(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    real(kind=dp), allocatable :: wrk(:,:), xu(:,:), yu(:,:), &
                                  fa(:,:), fb(:,:), fpk(:), &
                                  gamma_a(:,:), gamma_b(:,:), mt(:,:), &
                                  ha_yx(:,:), hb_yx(:,:), &
                                  ha_xy(:,:), hb_xy(:,:), &
                                  xuvec(:), yuvec(:), &
                                  hx_tmp(:,:), hx_g(:,:), hx_f7(:,:)
    integer(c_int), pointer :: ixcore_ptr(:)
    real(kind=dp) :: scale_exch, hfs
    integer :: nbf, nbf2, noca, nocb, nvirb, xdim, mrst, k, ok
    integer :: ipair, nrhs, source_x, source_y
    logical :: dft

    if (infos%tddft%umrsf) &
      call show_message('Analytic NAC wpair is not available for UMRSF-TDDFT.', with_abort)
    mrst = infos%tddft%mult
    if (mrst /= 1 .and. mrst /= 3) &
      call show_message('Analytic NAC wpair requires singlet or triplet MRSF-TDDFT.', with_abort)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    noca = infos%mol_prop%nelec_a
    nocb = infos%mol_prop%nelec_b
    nvirb = nbf - nocb
    xdim = noca*nvirb
    nrhs = size(ytil_batch,2)
    if (nrhs < 1 .or. nrhs > max_batch_width .or. &
        size(ytil_batch,1) /= xdim .or. &
        size(xstate_batch,1) /= xdim .or. &
        size(xstate_batch,2) /= nrhs .or. &
        size(mt_batch,1) /= nbf .or. size(mt_batch,2) /= nbf .or. &
        size(mt_batch,3) /= nrhs) &
      call show_message('Invalid MRSF NAC wpair batch dimensions.', with_abort)
    dft = infos%control%hamilton == 20
    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%tddft%hfscale
    hfs = infos%tddft%hfscale

    if (abs(hfs) <= epsilon(1.0_dp) .and. &
        (infos%tddft%spc_coco /= 0.0_dp .or. &
         infos%tddft%spc_ovov /= 0.0_dp .or. &
         infos%tddft%spc_coov /= 0.0_dp)) &
      call show_message('MRSF-TDDFT spin-pair coupling overrides require nonzero HFscale.', with_abort)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)

    allocate(mrsf_density(2*nrhs,7,nbf,nbf), wrk(nbf,nbf), &
             xu(nbf,nbf), yu(nbf,nbf), xuvec(xdim), yuvec(xdim), &
             fa(nbf,nbf), fb(nbf,nbf), fpk(nbf2), &
             gamma_a(nbf,nbf), gamma_b(nbf,nbf), mt(nbf,nbf), &
             ha_yx(nbf,nbf), hb_yx(nbf,nbf), &
             ha_xy(nbf,nbf), hb_xy(nbf,nbf), &
             hx_tmp(nbf,nbf), hx_g(nbf,nbf), hx_f7(nbf,nbf), &
             source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)

    ! Raw frozen Fock in the unshifted MO basis.  The ixcore overwrite is a
    ! constant diagonal replacement, whose derivative is removed below.
    call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, fpk)
    call unpack_matrix(fpk, fa)
    call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, fpk)
    call unpack_matrix(fpk, fb)

    ! The ERI source uses the folded eigenvectors, exactly as the matvec does.
    do ipair = 1, nrhs
      source_x = 2*ipair - 1
      source_y = source_x + 1
      call iatogen(xstate_batch(:,ipair), wrk, noca, nocb)
      call mrsfcbc(infos, mo_a, mo_b, wrk, &
                   mrsf_density(source_x,:,:,:))
      call iatogen(ytil_batch(:,ipair), wrk, noca, nocb)
      call mrsfcbc(infos, mo_a, mo_b, wrk, &
                   mrsf_density(source_y,:,:,:))
    end do

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_data_st = int2_mrsf_data_t( &
      d3 = mrsf_density(:2*nrhs,:,:,:), &
      tamm_dancoff = .true., &
      scale_exchange = scale_exch, &
      scale_coulomb = scale_exch)
    call int2_driver%run(int2_data_st, &
      cam = dft .and. infos%dft%cam_flag, &
      alpha = infos%tddft%cam_alpha, &
      alpha_coulomb = infos%tddft%cam_alpha, &
      beta = infos%tddft%cam_beta, &
      beta_coulomb = infos%tddft%cam_beta, &
      mu = infos%tddft%cam_mu)
    fmrst2 => int2_data_st%f3(:,:,:,:,1)

    if (mrst == 3) fmrst2(:,1:6,:,:) = -fmrst2(:,1:6,:,:)
    if (abs(hfs) > epsilon(1.0_dp)) then
      if (infos%tddft%spc_coco /= hfs) &
        fmrst2(:,6,:,:) = fmrst2(:,6,:,:)*infos%tddft%spc_coco/hfs
      if (infos%tddft%spc_ovov /= hfs) &
        fmrst2(:,5,:,:) = fmrst2(:,5,:,:)*infos%tddft%spc_ovov/hfs
      if (infos%tddft%spc_coov /= hfs) &
        fmrst2(:,1:4,:,:) = fmrst2(:,1:4,:,:)*infos%tddft%spc_coov/hfs
    end if

    if (infos%tddft%ixcore_len > 0) &
      call c_f_pointer(infos%tddft%ixcore, ixcore_ptr, &
                       [infos%tddft%ixcore_len])

    do ipair = 1, nrhs
      source_x = 2*ipair - 1
      source_y = source_x + 1
      ! The derivative-side amplitudes are unfolded once.  Pre-unfolding the
      ! density source above would apply the SOMO fold twice.
      call mrsfxvec(infos, xstate_batch(:,ipair), xuvec)
      call iatogen(xuvec, xu, noca, nocb)
      call mrsfxvec(infos, ytil_batch(:,ipair), yuvec)
      call iatogen(yuvec, yu, noca, nocb)

      ! H(y,Kx) + H(x,Ky), including all six spin-pair channels.  The half
      ! factor is the adjoint symmetrisation of B(C)^T K B(C).
      call mrsf_nac_hx_side(mo_a, mo_b, yu, fmrst2(source_x,:,:,:), &
                            noca, nocb, ha_yx, hb_yx, &
                            hx_tmp, hx_g, hx_f7)
      call mrsf_nac_hx_side(mo_a, mo_b, xu, fmrst2(source_y,:,:,:), &
                            noca, nocb, ha_xy, hb_xy, &
                            hx_tmp, hx_g, hx_f7)
      mt = 0.5_dp*(ha_yx + hb_yx + ha_xy + hb_xy)

      ! Frozen-Fock part: d Tr[Gamma C^T F_AO C]/d theta.
      gamma_a = 0.0_dp
      gamma_b = 0.0_dp
      gamma_a(1:noca,1:noca) = -0.5_dp*( &
        matmul(yu(1:noca,nocb+1:nbf), transpose(xu(1:noca,nocb+1:nbf))) + &
        matmul(xu(1:noca,nocb+1:nbf), transpose(yu(1:noca,nocb+1:nbf))))
      gamma_b(nocb+1:nbf,nocb+1:nbf) = 0.5_dp*( &
        matmul(transpose(yu(1:noca,nocb+1:nbf)), &
               xu(1:noca,nocb+1:nbf)) + &
        matmul(transpose(xu(1:noca,nocb+1:nbf)), &
               yu(1:noca,nocb+1:nbf)))
      ! Restrict the Fock contractions to the nonzero occupied/virtual blocks.
      call dgemm('n', 'n', nbf, noca, noca, 2.0_dp, &
                 fa, nbf, gamma_a, nbf, 1.0_dp, mt, nbf)
      call dgemm('n', 'n', nbf, nvirb, nvirb, 2.0_dp, &
                 fb(1,nocb+1), nbf, gamma_b(nocb+1,nocb+1), nbf, &
                 1.0_dp, mt(1,nocb+1), nbf)

      if (infos%tddft%ixcore_len > 0) then
        do k = 1, nocb
          if (.not. any(ixcore_ptr(1:infos%tddft%ixcore_len) == k)) &
            mt(:,k) = mt(:,k) - 2.0_dp*fa(:,k)*gamma_a(k,k)
        end do
      end if
      mt_batch(:,:,ipair) = mt
    end do

    call int2_driver%clean()

  contains

    subroutine mrsf_nac_hx_side(ca, cb, v, f, nocc_a, nocc_b, ha, hb, &
                                tmp, g, f7)
      real(kind=dp), intent(in) :: ca(:,:), cb(:,:), v(:,:)
      real(kind=dp), intent(in), target :: f(:,:,:)
      integer, intent(in) :: nocc_a, nocc_b
      real(kind=dp), intent(out) :: ha(:,:), hb(:,:)
      real(kind=dp), intent(inout) :: tmp(:,:), g(:,:), f7(:,:)
      integer :: n, nvir

      n = size(ca,1)
      nvir = n - nocc_b
      ha = 0.0_dp
      hb = 0.0_dp

      ! General channel: g = C_alpha^T F_7 C_beta.
      f7 = f(7,:,:)
      call dgemm('t', 'n', n, n, n, 1.0_dp, ca, n, f7, n, &
                 0.0_dp, tmp, n)
      call dgemm('n', 'n', n, n, n, 1.0_dp, tmp, n, cb, n, &
                 0.0_dp, g, n)
      ! iatogen leaves v nonzero only in (1:nocc_a,nocc_b+1:n).
      ! Restrict both general-channel products to that rectangular block.
      call dgemm('n', 't', n, nocc_a, nvir, 2.0_dp, &
                 g(:,nocc_b+1:n), n, v(:,nocc_b+1:n), n, &
                 0.0_dp, ha, n)
      call dgemm('t', 'n', n, nvir, nocc_a, 2.0_dp, &
                 g, n, v(:,nocc_b+1:n), n, &
                 0.0_dp, hb(:,nocc_b+1:n), n)

      ! mrsfsp accumulates the channel 1:6 adjoint into the general part.
      call mrsfsp(ha, hb, ca, cb, v, f, nocc_a, nocc_b)
    end subroutine mrsf_nac_hx_side

  end subroutine mrsf_nac_wpair_batch_impl
