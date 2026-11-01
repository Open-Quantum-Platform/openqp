module soc_mrsf_mod

    
  use precision, only: dp  

  implicit none

  character(len=*), parameter :: module_name = "soc_mrsf_mod"


  !> Determinant expansion of an MRSF state in the sector representation
  !> described above compute_tdm's configuration algebra.
  !> The one-hole/one-particle sector is kept factorized: an MRSF core->virtual
  !> excitation populates two open/spin patterns per (c,w) with the same
  !> amplitude, and the spin ladder acts on the pattern only, so
  !>   hp(2c+sr, o, 2w+sv) = hc(sr, o, sv) * xcv(c, w)
  !> exactly (no dense 64*nc*nv tensor is ever formed).
  type sector_state_t
    integer :: nc = 0, nv = 0
    real(kind=dp), allocatable :: nn(:), hn(:,:), np(:,:)
    real(kind=dp) :: hc(0:1, 0:15, 0:1) = 0.0_dp
    real(kind=dp), allocatable :: xcv(:,:)
  end type sector_state_t

  logical, save :: open_tables_ready = .false.
  integer, save :: pop_tab(0:15), before_tab(0:3, 0:15), after_tab(0:3, 0:15)
  logical, save :: has_tab(0:3, 0:15)
  real(kind=dp), save :: lop_tab(0:3, 0:3, 0:15, 0:15), splus_open(0:15, 0:15)

  !> MO-basis SOC integrals, exported only in debug mode ([input] verbose=3,
  !> infos%tddft%debug_mode) for independent checks of H_SOC
  !> (tests/test_soc_mrsf_tdm.py): shape (nbf, nbf, 3) = (t, u, x/y/z), the real
  !> antisymmetric l_b(t,u) with L_b = -i l_b; the 2e record exists only when
  !> the mean-field 2e part was computed (soc_2e=1).
  character(len=*), parameter :: OQP_soc_lmo_1e = "OQP::soc_lmo_1e"
  character(len=*), parameter :: OQP_soc_lmo_2e = "OQP::soc_lmo_2e"

  private
  ! bind(C) entry points must stay public: GCC 16.2 hides PRIVATE ones (GCC PR126872)
  public :: soc_mrsf_C
  public soc_mrsf

contains

  subroutine soc_mrsf_C(c_handle) bind(C, name="soc_mrsf")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call soc_mrsf(inf)
  end subroutine soc_mrsf_C

!> @brief Compute spin-orbit coupling corrections for MRSF-TDDFT states
!> @details
!>  Driver for the MRSF SOC calculation. Performs the following steps:
!>    1. Compute 1e SOC AO integrals <mu|Z*L/r^3|nu> via Breit-Pauli operator
!>    2. Transform AO integrals to MO basis
!>    3. Optionally add 2e mean-field SOC correction (controlled by infos%control%soc_2e)
!>    4. Build spin-dependent transition density matrices from MRSF Davidson vectors
!>    5. Assemble the SOC Hamiltonian in the (singlet + 3*triplet) basis
!>    6. Diagonalize to obtain SOC-corrected adiabatic energies and eigenvectors
!>
!>  State ordering in the OpenQP SOC basis is:
!>    indices 1..ns        -> singlet states S0..S(ns-1)
!>    indices ns+1..ns+3nt -> triplet Ms sublevels T0(Ms=-1,0,+1), T1(...), ...
!>
!> @param[inout] infos  OQP information struct (basis, atoms, control, tagarray, log)
  subroutine soc_mrsf(infos)
    use io_constants, only: iw
    use types, only: information
    use oqp_tagarray_driver
    use precision, only: dp
    use printing, only: print_module_info
    use messages, only: show_message, with_abort
    use grd2_rys, only: soc2e_driver
    use mathlib, only: orthogonal_transform
    use parallel, only: par_env_t
   use physical_constants, only: alpha  => FINE_STRUCTURE, &
                                 ha2wn  => HA_TO_WAVENUM,  &
                                 ha2ev  => EV2HTREE
    implicit none

    character(len=*), parameter :: subroutine_name = "soc_mrsf"

    type(information), target, intent(inout) :: infos

    real(kind=dp), contiguous, pointer :: singlet_energies(:), triplet_energies(:)
    real(kind=dp), contiguous, pointer :: bvec_mo_s(:,:), bvec_mo_t(:,:), mo_a(:,:)
    real(kind=dp) :: e_ref
    integer :: ok, nbf, nbf2
    real(kind=dp), allocatable :: lx_ao(:), ly_ao(:), lz_ao(:)
    real(kind=dp), allocatable :: lx_2e_ao(:), ly_2e_ao(:), lz_2e_ao(:)
    integer :: ns, nt, ist, jst, ims, ims_i, ims_j, itemp, i, idx, j
    real(kind=dp), allocatable :: lx_mo(:,:), ly_mo(:,:), lz_mo(:,:)
    real(kind=dp), allocatable :: t00aa(:,:,:,:), t110aa(:,:,:,:), t11ab(:,:,:,:)
    integer :: nocca, noccb

    complex(kind=dp), allocatable :: hsoc(:,:), h1soc(:,:), h2soc(:,:)
    real(kind=dp), allocatable :: eval(:)
    complex(kind=dp), allocatable :: evec(:,:)
    real(kind=dp), parameter :: dfac  = alpha**2 / 2.0_dp * ha2wn  ! 5.8438 cm-1/a.u.
    real(kind=dp) :: re1e, im1e, re2e, im2e, abs12e
    character(len=7), dimension(3), parameter :: trip = ['(Ms=-1)', '(Ms= 0)', '(Ms=+1)']

    real(kind=dp), allocatable :: den_rohf(:,:)
    real(kind=dp), allocatable :: wao(:,:,:)
    real(kind=dp), allocatable :: lx_2e_mo(:,:), ly_2e_mo(:,:), lz_2e_mo(:,:)
    real(kind=dp), allocatable :: lx_12e_mo(:,:), ly_12e_mo(:,:), lz_12e_mo(:,:)

    logical :: do_2e_soc
    logical :: debug_soc_prints
    type(par_env_t) :: pe

    integer :: nstate_soc
    real(kind=dp), pointer :: eval_out(:)
    real(kind=dp), pointer :: evec_re_out(:,:), evec_im_out(:,:)
    real(kind=dp), pointer :: hsoc_re_out(:,:), hsoc_im_out(:,:)

    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    debug_soc_prints = (infos%control%verbose >= 3)

    do_2e_soc = (infos%control%soc_2e /= 0)

    ! The 2e SOC kernel (grd2_rys soc2e path) still scatters Cartesian
    ! component counts against spherical AO offsets: under ispher with pure
    ! shells it would silently corrupt memory. Abort until it is ported.
    block
      use constants, only: HARMONIC_ACTIVE
      if (do_2e_soc .and. HARMONIC_ACTIVE) then
        if (any(infos%basis%harmonic == 1)) &
          call show_message('soc_mrsf: 2e SOC is not yet available with '// &
                            'spherical-harmonic AOs; set soc_2e=0 or ispher=false', WITH_ABORT)
      end if
    end block

    if (pe%rank == 0) then
      open(unit=iw, file=infos%log_filename, position="append")
      if (do_2e_soc) then
        call print_module_info('SOC_MRSF (1e+2e)', 'Spin-Orbit Coupling: MRSF Energies')
      else
        call print_module_info('SOC_MRSF (1e)', 'Spin-Orbit Coupling: MRSF Energies')
      end if
    end if

!    write(iw, *) 'Do we 2e?', infos%control%soc_2e
    e_ref = infos%mol_energy%energy

    call data_has_tags(infos%dat, &
        (/ character(len=80) :: OQP_td_singlet_energies, OQP_td_triplet_energies /), &
        module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_td_singlet_energies, singlet_energies)
    call tagarray_get_data(infos%dat, OQP_td_triplet_energies, triplet_energies)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo_s, bvec_mo_s)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo_t, bvec_mo_t)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    nbf  = infos%basis%nbf
    nbf2 = nbf*(nbf+1)/2

    ns = size(singlet_energies)
    nt = size(triplet_energies)

    nstate_soc = ns + 3*nt

    ! Number of alpha/beta occupied MOs, needed for TDM flat index mapping
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b

    ! compute_tdm addresses the stock packed grid nocca*(nbf-noccb) per state.
    ! A vector of any other length (e.g. an extended response space) would be
    ! read with the states misaligned; refuse it instead of returning garbage.
    ! Only the total size is checked: the Python side may store the tag with
    ! its two dimensions in either order, the memory layout being the same.
    if (size(bvec_mo_s) /= nocca*(nbf-noccb)*ns .or. size(bvec_mo_t) /= nocca*(nbf-noccb)*nt) then
      if (pe%rank == 0) then
        write(iw,'(/,a,i0,a,i0,a,i0,a)') 'soc_mrsf: MRSF response vectors do not have the stock '// &
          'layout: expected ', nocca*(nbf-noccb), ' per state, found ', &
          size(bvec_mo_s)/max(ns,1), ' (singlet) / ', size(bvec_mo_t)/max(nt,1), ' (triplet)'
        call flush(iw)
        close(iw)
      end if
      call show_message('soc_mrsf: response-vector layout is not the stock MRSF grid', WITH_ABORT)
    end if

    ! --- Step 1: Compute SOC 1e AO integrals <mu|Z*L/r^3|nu> ---
    allocate(lx_ao(nbf2), ly_ao(nbf2), lz_ao(nbf2), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate AO L matrices', WITH_ABORT)
    call compute_soc_ao(infos, lx_ao, ly_ao, lz_ao)

    ! --- Step 2: Transform AO integrals to MO basis: L_MO = C^T * L_AO * C ---
    allocate(lx_mo(nbf,nbf), ly_mo(nbf,nbf), lz_mo(nbf,nbf), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate MO L matrices', WITH_ABORT)
    call ao2mo_soc(lx_ao, lx_mo, mo_a, nbf)
    call ao2mo_soc(ly_ao, ly_mo, mo_a, nbf)
    call ao2mo_soc(lz_ao, lz_mo, mo_a, nbf)

    allocate(lx_12e_mo(nbf, nbf), ly_12e_mo(nbf, nbf), lz_12e_mo(nbf, nbf), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate total MO L matrices', WITH_ABORT)

    ! --- Step 2b: compute SOC 2e AO integrals, AO2MO transformation ---
    if (do_2e_soc) then
      allocate(den_rohf(nbf, nbf), stat=ok)
      if (ok /= 0) call show_message('soc_mrsf: cannot allocate ROHF density', WITH_ABORT)
      allocate(lx_2e_mo(nbf, nbf), ly_2e_mo(nbf, nbf), lz_2e_mo(nbf, nbf), stat=ok)
      if (ok /= 0) call show_message('soc_mrsf: cannot allocate 2e MO L matrices', WITH_ABORT)

      den_rohf = 0.0_dp
      call dgemm('N','T', nbf, nbf, noccb, 1.0_dp, &
                 mo_a, nbf, mo_a, nbf, 0.0_dp, den_rohf, nbf)

      allocate(wao(3, nbf, nbf), stat=ok)
      if (ok /= 0) call show_message('soc_mrsf: cannot allocate 2e AO L matrices', WITH_ABORT)
      wao = 0.0_dp

      if (pe%rank == 0 .and. debug_soc_prints) then
        do i = 1, infos%basis%nshell
          write(iw,'(a,3i4)') 'shell ao_offset am:', i, &
            infos%basis%ao_offset(i), infos%basis%am(i)
        end do
        write(iw,'(a,2i4)') ' nocca, noccb = ', nocca, noccb
      end if

      call soc2e_driver(infos, infos%basis, den_rohf, wao)

      if (pe%rank == 0 .and. debug_soc_prints) then
        write(iw,'(/,a)') ' LX in AO (our wao)'
        do i = 1, nbf
          write(iw,'(*(f12.6))') (wao(1,i,j), j=1,i)
        end do
        write(iw,'(/,a)') ' LY in AO (our wao)'
        do i = 1, nbf
          write(iw,'(*(f12.6))') (wao(2,i,j), j=1,i)
        end do
        write(iw,'(/,a)') ' LZ in AO (our wao)'
        do i = 1, nbf
          write(iw,'(*(f12.6))') (wao(3,i,j), j=1,i)
        end do
        write(iw,'(/,a,3es14.6)') '  ||wao|| (Lx,Ly,Lz) = ', &
          sqrt(sum(wao(1,:,:)**2)), &
          sqrt(sum(wao(2,:,:)**2)), &
          sqrt(sum(wao(3,:,:)**2))
      end if

      deallocate(den_rohf)

      call orthogonal_transform('n', nbf, mo_a, wao(1,:,:), lx_2e_mo)
      call orthogonal_transform('n', nbf, mo_a, wao(2,:,:), ly_2e_mo)
      call orthogonal_transform('n', nbf, mo_a, wao(3,:,:), lz_2e_mo)
      deallocate(wao)

      lx_12e_mo = lx_mo + lx_2e_mo
      ly_12e_mo = ly_mo + ly_2e_mo
      lz_12e_mo = lz_mo + lz_2e_mo
    else
      lx_12e_mo = lx_mo
      ly_12e_mo = ly_mo
      lz_12e_mo = lz_mo
    end if

    deallocate(lx_ao, ly_ao, lz_ao)

    ! Debug-mode export of the MO integrals (independent regression of H_SOC
    ! in tests/test_soc_mrsf_tdm.py).  Gated like the other MRSF developer
    ! dumps (infos%tddft%debug_mode, [input] verbose=3): a production run
    ! does not retain these 3*nbf*nbf records.  A record whose export
    ! condition is false now is erased, so a molecule re-run with lower
    ! verbosity or soc_2e=0 does not keep a stale copy from an earlier run.
    block
      real(kind=dp), contiguous, pointer :: lout(:,:,:)
      integer :: ta_status
      if (infos%tddft%debug_mode) then
        ta_status = infos%dat%alloc(OQP_soc_lmo_1e, [nbf, nbf, 3], lout, &
          description="MO SOC 1e integrals l_b(t,u) (x,y,z), L_b = -i l_b, a.u. (no alpha^2/2)", &
          override=.true.)
        if (ta_status /= TA_OK) call show_message('soc_mrsf: cannot allocate OQP::soc_lmo_1e', WITH_ABORT)
        lout(:,:,1) = lx_mo; lout(:,:,2) = ly_mo; lout(:,:,3) = lz_mo
      else
        call infos%dat%erase((/ character(len=80) :: OQP_soc_lmo_1e /))
      end if
      if (infos%tddft%debug_mode .and. do_2e_soc) then
        ta_status = infos%dat%alloc(OQP_soc_lmo_2e, [nbf, nbf, 3], lout, &
          description="MO SOC mean-field 2e integrals (x,y,z), same convention", &
          override=.true.)
        if (ta_status /= TA_OK) call show_message('soc_mrsf: cannot allocate OQP::soc_lmo_2e', WITH_ABORT)
        lout(:,:,1) = lx_2e_mo; lout(:,:,2) = ly_2e_mo; lout(:,:,3) = lz_2e_mo
      else
        call infos%dat%erase((/ character(len=80) :: OQP_soc_lmo_2e /))
      end if
    end block

    ! --- Step 3: Build spin-dependent transition density matrices ---
    allocate(t00aa (ns, nt, nbf, nbf), &
             t110aa(nt, nt, nbf, nbf), &
             t11ab (nt, nt, nbf, nbf), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate TDM arrays', WITH_ABORT)
    call compute_tdm(bvec_mo_s, bvec_mo_t, nocca*(nbf-noccb), nocca, noccb, nbf, ns, nt, &
                     t00aa, t110aa, t11ab)

    ! --- Step 4: Assemble the 1e SOC Hamiltonian H_SOC ---
    allocate(hsoc(ns + 3*nt, ns + 3*nt), stat=ok)
    allocate(h1soc(ns + 3*nt, ns + 3*nt), stat=ok)
    allocate(h2soc(ns + 3*nt, ns + 3*nt), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate H_SOC matrix', WITH_ABORT)
    h2soc = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call compute_soc_matrix(t00aa, t110aa, t11ab, lx_mo, ly_mo, lz_mo, ns, nt, nbf, h1soc)
    if (do_2e_soc) call compute_soc_matrix(t00aa, t110aa, t11ab, lx_2e_mo, ly_2e_mo, lz_2e_mo, ns, nt, nbf, h2soc)
    call compute_soc_matrix(t00aa, t110aa, t11ab, lx_12e_mo, ly_12e_mo, lz_12e_mo, ns, nt, nbf, hsoc)
    deallocate(t00aa, t110aa, t11ab, lx_mo, ly_mo, lz_mo)

    ! --- Step 5: Print SOC coupling constants, separated into 1e and 2e parts ---
    if (pe%rank == 0) then
      write(iw,'(/,11x,89("-"))')
      write(iw,'(41x,a)') 'Absolute = sqrt(Re(1e+2e)**2+Im(1e+2e)**2)'
      write(iw,'(11x,89("-"))')
      write(iw,'(2x,a,4x,a,9x,a,6x,a,6x,a,6x,a,6x,a)') &
        'State_i', 'State_j', 'Re(1e)', 'Im(1e)', 'Re(2e)', 'Im(2e)', 'Absolute'

      ! S-T block
      do ist = 1, ns
        do jst = 1, nt
          do ims = 1, 3  ! Ms = -1, 0, +1
            idx = ns + (jst-1)*3 + ims

            re1e   = real (h1soc(ist, idx)) * dfac
            im1e   = aimag(h1soc(ist, idx)) * dfac
            re2e   = real (h2soc(ist, idx)) * dfac
            im2e   = aimag(h2soc(ist, idx)) * dfac
            abs12e = sqrt((re1e + re2e)**2 + (im1e + im2e)**2)

            write(iw,'(5x,a,i0,4x,"/",x,a,i0,a,x,4f12.4,f18.12)') &
              'S', ist-1, 'T', jst-1, trim(trip(ims)), re1e, im1e, re2e, im2e, abs12e
          end do
        end do
      end do

      ! T-T block
      do ist = 1, nt
        do jst = 1, nt
          do ims_i = 1, 3
            do ims_j = 1, 3
              i = ns + (ist-1)*3 + ims_i
              j = ns + (jst-1)*3 + ims_j

              re1e   = real (h1soc(i, j)) * dfac
              im1e   = aimag(h1soc(i, j)) * dfac
              re2e   = real (h2soc(i, j)) * dfac
              im2e   = aimag(h2soc(i, j)) * dfac
              abs12e = sqrt((re1e + re2e)**2 + (im1e + im2e)**2)

              write(iw,'(5x,a,i0,a,4x,"/",x,a,i0,a,x,4f12.4,f18.12)') &
                'T', ist-1, trim(trip(ims_i)), 'T', jst-1, trim(trip(ims_j)), &
                re1e, im1e, re2e, im2e, abs12e
            end do
          end do
        end do
      end do
    end if

    ! --- Step 6: Diagonalize H_SOC + excitation energies, print eigenvalues ---
    allocate(eval(ns + 3*nt), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate eigenvalue array', WITH_ABORT)
    allocate(evec(ns + 3*nt, ns + 3*nt), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate eigenvector array', WITH_ABORT)
    call diag_soc(hsoc, singlet_energies, triplet_energies, e_ref, ns, nt, eval, evec)

    if (pe%rank == 0) then
      call print_soc_eigenvalues(iw, eval, evec, singlet_energies, triplet_energies, e_ref, ns, nt)
      call print_soc_decomposition(iw, eval, evec, ns, nt)

      call infos%dat%alloc_or_die(OQP_soc_eval, (/nstate_soc/), eval_out, description=OQP_soc_eval_comment)
      call infos%dat%alloc_or_die(OQP_soc_evec_re, (/nstate_soc, nstate_soc/), evec_re_out, description=OQP_soc_evec_re_comment)
      call infos%dat%alloc_or_die(OQP_soc_evec_im, (/nstate_soc, nstate_soc/), evec_im_out, description=OQP_soc_evec_im_comment)
      call infos%dat%alloc_or_die(OQP_soc_hsoc_re, (/nstate_soc, nstate_soc/), hsoc_re_out, description=OQP_soc_hsoc_re_comment)
      call infos%dat%alloc_or_die(OQP_soc_hsoc_im, (/nstate_soc, nstate_soc/), hsoc_im_out, description=OQP_soc_hsoc_im_comment)

      eval_out    = eval
      evec_re_out = real(evec, kind=dp)
      evec_im_out = aimag(evec)
      hsoc_re_out = real(hsoc, kind=dp)
      hsoc_im_out = aimag(hsoc)

      write(iw,'(/,a)') 'SOC_MRSF done'
      call flush(iw)
      close(iw)
    end if

    deallocate(eval, evec)
    deallocate(lx_12e_mo, ly_12e_mo, lz_12e_mo)
    deallocate(hsoc)
    if (do_2e_soc) then
      deallocate(lx_2e_mo, ly_2e_mo, lz_2e_mo)
    end if

  end subroutine soc_mrsf


!> @brief Compute 1-electron SOC AO integrals using the Breit-Pauli operator
!> @details
!>  Evaluates <mu|Z_A * L_A / r_A^3|nu> for each atom A, where L_A is the
!>  angular momentum operator relative to nucleus A and Z_A is the (effective)
!>  nuclear charge. Loops over shell pairs (ii >= jj) and accumulates into
!>  packed lower-triangular arrays. Results are normalised with basis function norms.
!>
!>  Note: uses bare nuclear charges (ze = Z).
!> @param[inout] infos   OQP information struct (basis, atoms)
!> @param[out]   lx_ao   Lx AO integrals, packed lower-triangular (nbf*(nbf+1)/2)
!> @param[out]   ly_ao   Ly AO integrals, packed lower-triangular
!> @param[out]   lz_ao   Lz AO integrals, packed lower-triangular
subroutine compute_soc_ao(infos, lx_ao, ly_ao, lz_ao)
  use basis_tools,       only: basis_set, bas_norm_matrix
  use cart2sph,          only: cart2sph_mat
  use mod_1e_primitives, only: comp_soc_int1_prim, update_triang_matrix
  use mod_shell_tools,   only: shell_t, shpair_t
  use constants,         only: HARMONIC_ACTIVE, tol_int
  use precision,         only: dp
  use types,             only: information
  use parallel,          only: par_env_t

  implicit none

  type(information), target, intent(inout) :: infos
  real(kind=dp), intent(out) :: lx_ao(:), ly_ao(:), lz_ao(:)

  type(basis_set), pointer :: basis
  type(shell_t)   :: shi, shj
  type(shpair_t)  :: cntp

  integer, parameter :: blocksize = 28*28
  real(kind=dp) :: socblk(blocksize, 3)

  integer  :: ii, jj, ig, iat, iz, nat, nbf, mpi_ii
  real(kind=dp) :: ze, tol
  type(par_env_t) :: pe

  basis => infos%basis
  basis%atoms => infos%atoms
  nat  = size(infos%atoms%zn)
  nbf  = basis%nbf
  tol  = log(10.0_dp) * tol_int

  lx_ao = 0.0_dp
  ly_ao = 0.0_dp
  lz_ao = 0.0_dp

  call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

!$omp parallel &
!$omp   private(shi, shj, cntp, socblk, ii, jj, ig, iat, iz, ze, mpi_ii) &
!$omp   reduction(+:lx_ao, ly_ao, lz_ao)

  call cntp%alloc(basis)

!$omp barrier
  if (infos%mpiinfo%usempi) mpi_ii = 0

  do ii = basis%nshell, 1, -1
    if (infos%mpiinfo%usempi) then
      mpi_ii = mpi_ii + 1
      if (mod(mpi_ii, pe%size) /= pe%rank) cycle
    end if
    call shi%fetch_by_id(basis, ii)

!$omp do schedule(dynamic)
    do jj = 1, ii
      call shj%fetch_by_id(basis, jj)
      call cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
      if (cntp%numpairs == 0) cycle

      socblk = 0.0_dp

      do ig = 1, cntp%numpairs
        do iat = 1, nat
          iz = nint(infos%atoms%zn(iat))
          ze = real(iz, dp)
          call comp_soc_int1_prim(cntp, ig, infos%atoms%xyz(:,iat), ze, socblk)
        end do
      end do

      if (HARMONIC_ACTIVE .and. (shi%harmonic == 1 .or. shj%harmonic == 1)) then
        call cart2sph_mat(socblk(:,1), shj%ang, shj%harmonic, shi%ang, shi%harmonic, &
                          iandj=(shi%shid==shj%shid), antisym=.true.)
        call cart2sph_mat(socblk(:,2), shj%ang, shj%harmonic, shi%ang, shi%harmonic, &
                          iandj=(shi%shid==shj%shid), antisym=.true.)
        call cart2sph_mat(socblk(:,3), shj%ang, shj%harmonic, shi%ang, shi%harmonic, &
                          iandj=(shi%shid==shj%shid), antisym=.true.)
      end if
      call update_triang_matrix(shi, shj, socblk(:,1), lx_ao)
      call update_triang_matrix(shi, shj, socblk(:,2), ly_ao)
      call update_triang_matrix(shi, shj, socblk(:,3), lz_ao)

    end do
!$omp end do

  end do
!$omp end parallel

  call pe%allreduce(lx_ao, size(lx_ao))
  call pe%allreduce(ly_ao, size(ly_ao))
  call pe%allreduce(lz_ao, size(lz_ao))

  call bas_norm_matrix(lx_ao, basis%bfnrm, nbf)
  call bas_norm_matrix(ly_ao, basis%bfnrm, nbf)
  call bas_norm_matrix(lz_ao, basis%bfnrm, nbf)

end subroutine compute_soc_ao

!> @brief Print a packed SOC AO integral matrix in a compact block format
!> @details
!>  Writes the lower-triangular AO integral matrix to unit iw in blocks of
!>  NCOLS=5 columns, with basis function labels and row indices.
!>
!> @param[in]  iw     Log file unit
!> @param[in]  comp   Component label ('LX', 'LY', or 'LZ')
!> @param[in]  mat    Packed lower-triangular AO matrix (nbf*(nbf+1)/2)
!> @param[in]  nbf    Number of basis functions
!> @param[in]  basis  Basis set descriptor (used for bf_label)
subroutine print_soc_ao(iw, comp, mat, nbf, basis)
  use basis_tools, only: basis_set
  use precision,   only: dp
  implicit none

  integer,          intent(in) :: iw, nbf
  character(len=2), intent(in) :: comp       ! 'LX', 'LY', or 'LZ'
  real(kind=dp),    intent(in) :: mat(nbf*(nbf+1)/2)
  type(basis_set),  intent(in) :: basis

  integer, parameter :: NCOLS = 5
  integer :: i, j, jstart, jend, jend_row, idx

  write(iw, '(/,2x,a)') comp//'  AO INTEGRALS'

  jstart = 1
  do while (jstart <= nbf)
    jend = min(jstart + NCOLS - 1, nbf)

    ! column index header
    write(iw, '(/,17x)', advance='no')
    do j = jstart, jend
      write(iw, '(i11)', advance='no') j
    end do
    write(iw, '(/)')

    ! data rows (lower triangle only)
    do i = jstart, nbf
      jend_row = min(jend, i)
      if (jend_row < jstart) cycle
      write(iw, '(i5,2x,a8,2x)', advance='no') i, basis%bf_label(i)
      do j = jstart, jend_row
        idx = i*(i-1)/2 + j
        write(iw, '(f11.6)', advance='no') mat(idx)
      end do
      write(iw, *)
    end do

    jstart = jend + 1
  end do
  write(iw, *)

end subroutine print_soc_ao

!> @brief Transform a packed antisymmetric SOC AO matrix to the MO basis
!> @details
!>  Unpacks the lower-triangular AO integral into a full antisymmetric matrix
!>  (L(nu,mu) = -L(mu,nu), diagonal = 0), then applies the two-step MO
!>  transformation: tmp = L_AO * C, L_MO = C^T * tmp via DGEMM.
!>
!> @param[in]  l_tri  Packed lower-triangular AO integrals (nbf*(nbf+1)/2)
!> @param[out] l_mo   Full MO integral matrix (nbf x nbf)
!> @param[in]  cmo    MO coefficient matrix C(mu,p) (nbf x nbf)
!> @param[in]  nbf    Number of basis functions
subroutine ao2mo_soc(l_tri, l_mo, cmo, nbf)
  use precision, only: dp
  implicit none

  real(kind=dp), intent(in)  :: l_tri(nbf*(nbf+1)/2)  ! AO integrals, packed lower triangle
  real(kind=dp), intent(out) :: l_mo(nbf, nbf)         ! MO integrals, full matrix
  real(kind=dp), intent(in)  :: cmo(nbf, nbf)          ! MO coefficient matrix C(mu,p)
  integer,       intent(in)  :: nbf

  real(kind=dp), allocatable :: l_full(:,:), tmp(:,:)
  integer :: mu, nu, idx

  allocate(l_full(nbf,nbf), tmp(nbf,nbf))

  ! Unpack lower triangle into full antisymmetric matrix: L(nu,mu) = -L(mu,nu)
  l_full = 0.0_dp
  do mu = 1, nbf
    do nu = 1, mu-1
      idx = mu*(mu-1)/2 + nu
      l_full(mu,nu) =  l_tri(idx)
      l_full(nu,mu) = -l_tri(idx)
    end do
    ! diagonal is zero by antisymmetry
  end do

  ! tmp(mu,q) = sum_nu L^AO(mu,nu) * C(nu,q)
  call dgemm('N','N', nbf, nbf, nbf, &
             1.0_dp, l_full, nbf, &
                     cmo,   nbf, &
             0.0_dp, tmp,   nbf)

  ! L^MO(p,q) = sum_mu C(mu,p) * tmp(mu,q)  =  C^T * tmp
  call dgemm('T','N', nbf, nbf, nbf, &
             1.0_dp, cmo, nbf, &
                     tmp, nbf, &
             0.0_dp, l_mo, nbf)

  deallocate(l_full, tmp)

end subroutine ao2mo_soc

!> @brief Spin-component transition density matrices of the MRSF states
!> @details
!>  Every MRSF M_S = 0 state is expanded in determinants (configuration
!>  algebra below) and the one-particle transition densities
!>  D[P,Q] = <bra| a+_P a_Q |ket> over spin orbitals P = 2m+s are evaluated
!>  exactly; nothing is hand-coded per configuration class.  The arrays
!>  returned keep the convention of compute_soc_matrix (second index of the
!>  TDM pairs with the first index of the MO integral, l(t,u) <-> D(u,t)):
!>    t00aa (I,J,t,u) = <S_I  | a+_{u a} a_{t a} |T_J,0 >
!>    t110aa(I,J,t,u) = <T_I,0| a+_{u a} a_{t a} |T_J,0 >
!>    t11ab (I,J,t,u) = <T_I,0| a+_{u b} a_{t a} |T_J,+1>,  |T,+1> = S+|T,0>/sqrt2
!>  The former implementation built the core-open elements for the top core
!>  orbital only, transposed the core-core block and omitted every term of
!>  first order in the core->virtual amplitudes.
!>
!>  State map (two-SOMO MRSF, packed slot (i,a) = alpha-occupied i -> beta-virtual a):
!>    Psi = 2^-1/2 sum_{(i,a) not OO} X(i,a) [E+(i,a) + lam E-(i,a)]
!>        + 2^-1/2 X(O1,O1) [L + lam R] + [singlet] X(O2,O1) G + X(O1,O2) D
!>    E+(i,a) = a+(a b) a(i a) R+,  E-(i,a) = a+(a a) a(i b) R-,
!>    R+ = |closed O1a O2a>, R- = |closed O1b O2b>, L = E+(O1,O1), R = E+(O2,O2),
!>    G = E+(O2,O1), D = E+(O1,O2), lam = -1 (singlet) / +1 (triplet).
!>
!> @param[in]  bvec_s   Singlet Davidson vectors (xvec_dim x ns)
!> @param[in]  bvec_t   Triplet Davidson vectors (xvec_dim x nt)
!> @param[in]  xvec_dim nocca*(nbf-noccb)
!> @param[in]  nocca    Number of alpha occupied MOs
!> @param[in]  noccb    Number of beta  occupied MOs
!> @param[in]  nbf      Number of basis functions
!> @param[in]  ns, nt   Number of singlet/triplet states
!> @param[out] t00aa    Singlet-triplet TDM (ns x nt x nbf x nbf)
!> @param[out] t110aa   Triplet-triplet TDM, Ms=0 sector (nt x nt x nbf x nbf)
!> @param[out] t11ab    Triplet-triplet TDM, Ms=0/+1 sector (nt x nt x nbf x nbf)
subroutine compute_tdm(bvec_s, bvec_t, xvec_dim, nocca, noccb, nbf, ns, nt, &
                       t00aa, t110aa, t11ab)
  use precision, only: dp
  use messages, only: show_message, with_abort
  implicit none

  integer,       intent(in)  :: xvec_dim, nocca, noccb, nbf, ns, nt
  real(kind=dp), intent(in)  :: bvec_s(xvec_dim, ns)
  real(kind=dp), intent(in)  :: bvec_t(xvec_dim, nt)
  real(kind=dp), intent(out) :: t00aa (ns, nt, nbf, nbf)
  real(kind=dp), intent(out) :: t110aa(nt, nt, nbf, nbf)
  real(kind=dp), intent(out) :: t11ab (nt, nt, nbf, nbf)

  type(sector_state_t), allocatable :: sst(:), tt0(:), tt1(:)
  real(kind=dp), allocatable :: d(:,:)
  integer :: i, j, t, u

  if (nocca - noccb /= 2) &
    call show_message('soc_mrsf: two-SOMO MRSF requires nocca - noccb = 2', WITH_ABORT)
  if (xvec_dim /= nocca*(nbf - noccb)) &
    call show_message('soc_mrsf: response vectors do not have the stock MRSF layout', WITH_ABORT)

  call init_open_tables()

  allocate(sst(ns), tt0(nt), tt1(nt))
  do i = 1, ns
    call mrsf_sector_state(bvec_s(:, i), 1, nocca, noccb, nbf, sst(i))
  end do
  do j = 1, nt
    call mrsf_sector_state(bvec_t(:, j), 3, nocca, noccb, nbf, tt0(j))
    call spin_raise(tt0(j), tt1(j))
    tt1(j)%nn = tt1(j)%nn / sqrt(2.0_dp)
    tt1(j)%hn = tt1(j)%hn / sqrt(2.0_dp)
    tt1(j)%np = tt1(j)%np / sqrt(2.0_dp)
    tt1(j)%hc = tt1(j)%hc / sqrt(2.0_dp)
  end do

  allocate(d(0:2*nbf-1, 0:2*nbf-1))
  t00aa  = 0.0_dp
  t110aa = 0.0_dp
  t11ab  = 0.0_dp

  do j = 1, nt
    do i = 1, ns
      call transition_density(sst(i), tt0(j), d)
      do u = 1, nbf
        do t = 1, nbf
          t00aa(i, j, t, u) = d(2*(u-1), 2*(t-1))
        end do
      end do
    end do
  end do

  do j = 1, nt
    do i = 1, nt
      call transition_density(tt0(i), tt0(j), d)
      do u = 1, nbf
        do t = 1, nbf
          t110aa(i, j, t, u) = d(2*(u-1), 2*(t-1))
        end do
      end do
      call transition_density(tt0(i), tt1(j), d)
      do u = 1, nbf
        do t = 1, nbf
          t11ab(i, j, t, u) = d(2*(u-1)+1, 2*(t-1))
        end do
      end do
    end do
  end do

  deallocate(d, sst, tt0, tt1)

end subroutine compute_tdm

! ---------------------------------------------------------------------------
! Configuration algebra for two-SOMO MRSF states
!
! A determinant of an MRSF state (or of its S+- partners) has a doubly
! occupied closed set with at most one hole, any occupation of the four open
! spin orbitals (O1a, O1b, O2a, O2b) and at most one electron in the virtual
! set.  A state is stored as four sector tensors
!     nn(o)        no core hole, no virtual electron
!     hn(r, o)     core hole r,  no virtual electron
!     np(o, v)     no core hole, virtual electron v
!     hp(r, o, v)  core hole r,  virtual electron v  (= hc(sr,o,sv) xcv(c,w))
! with o = 0..15 the open-shell occupation bit pattern (bit 0 O1a, 1 O1b,
! 2 O2a, 3 O2b), r = 2c+s the core spin orbital (s = 0 alpha, 1 beta) and
! v = 2w+s the virtual spin orbital.  The basis vectors are occupation-number
! vectors in the fermion-mode order (virtual, open, core); a+_p a_q (p /= q)
! carries (-1)**N with N the number of occupied modes strictly between p and q.
! The spin-orbital index of the transition density is 2m+s, m = core 0..nc-1,
! O1 = nc, O2 = nc+1, virtuals nc+2.. ; the mode order coincides with it.
! ---------------------------------------------------------------------------

!> Bit tables of the four open modes and the one-body operators
!> lop(k,l,o1,o2) = <o1| a+_k a_l |o2>.
subroutine init_open_tables()
  implicit none
  integer :: k, l, o1, o2, lo, hi, between

  if (open_tables_ready) return
  do o2 = 0, 15
    pop_tab(o2) = popcnt(o2)
    do k = 0, 3
      has_tab(k, o2) = btest(o2, k)
      before_tab(k, o2) = popcnt(iand(o2, ishft(1, k) - 1))
      after_tab(k, o2) = popcnt(ishft(o2, -(k+1)))
    end do
  end do
  lop_tab = 0.0_dp
  do k = 0, 3
    do l = 0, 3
      do o2 = 0, 15
        if (.not. btest(o2, l)) cycle
        if (k == l) then
          lop_tab(k, l, o2, o2) = 1.0_dp
          cycle
        end if
        if (btest(o2, k)) cycle
        lo = min(k, l)
        hi = max(k, l)
        between = popcnt(iand(ishft(o2, -(lo+1)), ishft(1, hi-lo-1) - 1))
        o1 = ior(iand(o2, not(ishft(1, l))), ishft(1, k))
        if (mod(between, 2) == 1) then
          lop_tab(k, l, o1, o2) = -1.0_dp
        else
          lop_tab(k, l, o1, o2) = 1.0_dp
        end if
      end do
    end do
  end do
  splus_open = lop_tab(0, 1, :, :) + lop_tab(2, 3, :, :)
  open_tables_ready = .true.
end subroutine init_open_tables

pure function parity_sign(n) result(s)
  implicit none
  integer, intent(in) :: n
  real(kind=dp) :: s
  if (mod(n, 2) == 0) then
    s = 1.0_dp
  else
    s = -1.0_dp
  end if
end function parity_sign

subroutine alloc_sector_state(st, nc, nv)
  implicit none
  type(sector_state_t), intent(out) :: st
  integer, intent(in) :: nc, nv
  st%nc = nc
  st%nv = nv
  allocate(st%nn(0:15), st%hn(0:2*nc-1, 0:15), st%np(0:15, 0:2*nv-1), st%xcv(0:nc-1, 0:nv-1))
  st%nn = 0.0_dp
  st%hn = 0.0_dp
  st%np = 0.0_dp
  st%hc = 0.0_dp
  st%xcv = 0.0_dp
end subroutine alloc_sector_state

!> MRSF M_S = 0 state of one packed response vector x (slot (i,a) at
!> (a-noccb-1)*nocca + i), mult = 1 or 3.
subroutine mrsf_sector_state(x, mult, nocca, noccb, nbf, st)
  implicit none
  real(kind=dp), intent(in) :: x(:)
  integer, intent(in) :: mult, nocca, noccb, nbf
  type(sector_state_t), intent(out) :: st

  integer, parameter :: ref_plus = 5     ! O1a O2a  (bits 0 and 2)
  integer, parameter :: ref_minus = 10   ! O1b O2b  (bits 1 and 3)
  real(kind=dp), parameter :: isq2 = 1.0_dp / sqrt(2.0_dp)
  integer :: nc, nv, c, w, m, k, l, o, ra, rb
  real(kind=dp) :: lam, val

  nc = noccb
  nv = nbf - nocca
  lam = 1.0_dp
  if (mult == 1) lam = -1.0_dp
  call alloc_sector_state(st, nc, nv)

  ! core -> virtual: E+ puts (hole c alpha, pattern O1a O2a, particle w beta),
  ! E- puts (hole c beta, pattern O1b O2b, particle w alpha); the phases
  ! (-1)**(n_open + r) are +1 and -1 for every c, so the sector factorizes.
  do c = 0, nc - 1
    do w = 0, nv - 1
      st%xcv(c, w) = xslot(c, w + 2)
    end do
  end do
  st%hc(0, ref_plus, 1)  = isq2 * parity_sign(pop_tab(ref_plus))
  st%hc(1, ref_minus, 0) = lam * isq2 * parity_sign(pop_tab(ref_minus) + 1)

  ! core -> O_m
  do m = 0, 1
    do c = 0, nc - 1
      ra = 2*c
      rb = 2*c + 1
      val = xslot(c, m) * isq2
      k = 2*m + 1
      st%hn(ra, ior(ref_plus, ishft(1, k))) = st%hn(ra, ior(ref_plus, ishft(1, k))) &
        + parity_sign(after_tab(k, ref_plus) + ra) * val
      k = 2*m
      st%hn(rb, ior(ref_minus, ishft(1, k))) = st%hn(rb, ior(ref_minus, ishft(1, k))) &
        + lam * parity_sign(after_tab(k, ref_minus) + rb) * val
    end do
  end do

  ! O_m -> virtual
  do m = 0, 1
    do w = 0, nv - 1
      val = xslot(nc + m, w + 2) * isq2
      l = 2*m
      st%np(ieor(ref_plus, ishft(1, l)), 2*w + 1) = st%np(ieor(ref_plus, ishft(1, l)), 2*w + 1) &
        + parity_sign(before_tab(l, ref_plus)) * val
      l = 2*m + 1
      st%np(ieor(ref_minus, ishft(1, l)), 2*w) = st%np(ieor(ref_minus, ishft(1, l)), 2*w) &
        + lam * parity_sign(before_tab(l, ref_minus)) * val
    end do
  end do

  ! open -> open: L + lam R, and G, D for singlets
  do o = 0, 15
    st%nn(o) = st%nn(o) + xslot(nc, 0) * isq2 &
      * (lop_tab(1, 0, o, ref_plus) + lam * lop_tab(3, 2, o, ref_plus))
    if (mult == 1) then
      st%nn(o) = st%nn(o) + xslot(nc + 1, 0) * lop_tab(1, 2, o, ref_plus) &
                          + xslot(nc, 1)     * lop_tab(3, 0, o, ref_plus)
    end if
  end do

contains

  !> packed amplitude X(i, a) for 0-based alpha-occupied i and 0-based a-noccb
  pure function xslot(i0, a0) result(v)
    integer, intent(in) :: i0, a0
    real(kind=dp) :: v
    v = x(a0*nocca + i0 + 1)
  end function xslot

end subroutine mrsf_sector_state

!> S+ = sum_p a+(p alpha) a(p beta) applied to a sector state.
subroutine spin_raise(st, out)
  implicit none
  type(sector_state_t), intent(in) :: st
  type(sector_state_t), intent(out) :: out
  integer :: o1, o2, c, w

  call alloc_sector_state(out, st%nc, st%nv)
  out%xcv = st%xcv
  do o1 = 0, 15
    do o2 = 0, 15
      if (splus_open(o1, o2) == 0.0_dp) cycle
      out%nn(o1) = out%nn(o1) + splus_open(o1, o2) * st%nn(o2)
      out%hn(:, o1) = out%hn(:, o1) + splus_open(o1, o2) * st%hn(:, o2)
      out%np(o1, :) = out%np(o1, :) + splus_open(o1, o2) * st%np(o2, :)
      out%hc(:, o1, :) = out%hc(:, o1, :) + splus_open(o1, o2) * st%hc(:, o2, :)
    end do
  end do
  ! virtual electron (w beta) -> (w alpha)
  do w = 0, st%nv - 1
    out%np(:, 2*w) = out%np(:, 2*w) + st%np(:, 2*w + 1)
  end do
  out%hc(:, :, 0) = out%hc(:, :, 0) + st%hc(:, :, 1)
  ! core hole: a+(c alpha) a(c beta) moves the hole from (c alpha) to (c beta)
  do c = 0, st%nc - 1
    out%hn(2*c + 1, :) = out%hn(2*c + 1, :) + st%hn(2*c, :)
  end do
  out%hc(1, :, :) = out%hc(1, :, :) + st%hc(0, :, :)
end subroutine spin_raise

!> D(P,Q) = <bra| a+_P a_Q |ket>, P = 2m+s over all spin orbitals.
!> The hp sector enters only through the factorized form, so every contraction
!> over it is an O(nc nv) (or nc^2 / nv^2) product of the compact amplitudes.
subroutine transition_density(bra, ket, d)
  implicit none
  type(sector_state_t), intent(in) :: bra, ket
  real(kind=dp), intent(out) :: d(0:, 0:)

  integer :: nc, nv, nh, npv, o0, v0, p, q, a, b, o, k, pk, of, ox, r, m, c, w, sr, sp, sq, sv
  real(kind=dp) :: overlap, ph, s, sx
  real(kind=dp) :: gvv(0:1, 0:1), gcc(0:1, 0:1), goo(0:15, 0:15), acv(0:1, 0:1), avc(0:1, 0:1)
  real(kind=dp), allocatable :: mc(:,:), gam(:,:), mvv(:,:), mcc(:,:), qb(:,:), qk(:,:), ub(:,:), uk(:,:)

  nc = ket%nc
  nv = ket%nv
  nh = 2*nc
  npv = 2*nv
  o0 = 2*nc
  v0 = 2*nc + 4
  d = 0.0_dp

  ! contractions of the factorized hp sectors
  sx = sum(bra%xcv * ket%xcv)                                   ! sum_cw xb xk
  allocate(mvv(0:nv-1, 0:nv-1), mcc(0:nc-1, 0:nc-1))
  mvv = matmul(transpose(bra%xcv), ket%xcv)                     ! (w, w')
  mcc = matmul(bra%xcv, transpose(ket%xcv))                     ! (c, c')
  do sq = 0, 1
    do sp = 0, 1
      gvv(sp, sq) = sum(bra%hc(:, :, sp) * ket%hc(:, :, sq))
      gcc(sp, sq) = sum(bra%hc(sp, :, :) * ket%hc(sq, :, :))
    end do
  end do
  do m = 0, 15
    do o = 0, 15
      goo(o, m) = sum(bra%hc(:, o, :) * ket%hc(:, m, :))
    end do
  end do

  overlap = sum(bra%nn * ket%nn) + sum(bra%hn * ket%hn) + sum(bra%np * ket%np) &
          + sx * sum(bra%hc * ket%hc)

  ! virtual <- virtual (no intervening occupied mode)
  do q = 0, npv - 1
    do p = 0, npv - 1
      d(v0 + p, v0 + q) = sum(bra%np(:, p) * ket%np(:, q)) &
                        + mvv(p/2, q/2) * gvv(mod(p, 2), mod(q, 2))
    end do
  end do

  ! core <- core: ket hole b, bra hole a, phase (-1)**(|a-b|-1)
  allocate(mc(0:nh-1, 0:nh-1))
  do b = 0, nh - 1
    do a = 0, nh - 1
      mc(a, b) = sum(bra%hn(a, :) * ket%hn(b, :)) + mcc(a/2, b/2) * gcc(mod(a, 2), mod(b, 2))
    end do
  end do
  do b = 0, nh - 1
    do a = 0, nh - 1
      if (a == b) then
        d(a, a) = overlap - mc(a, a)
      else
        d(a, b) = parity_sign(abs(a - b) - 1) * mc(b, a)
      end if
    end do
  end do
  deallocate(mc)

  ! open <- open
  allocate(gam(0:15, 0:15))
  do m = 0, 15
    do o = 0, 15
      gam(o, m) = bra%nn(o) * ket%nn(m) + sum(bra%hn(:, o) * ket%hn(:, m)) &
                + sum(bra%np(o, :) * ket%np(m, :)) + sx * goo(o, m)
    end do
  end do
  do k = 0, 3
    do m = 0, 3
      d(o0 + k, o0 + m) = sum(gam * lop_tab(k, m, :, :))
    end do
  end do
  deallocate(gam)

  ! open <-> core, open <-> virtual.  qb(c,sp) = sum_w xb(c,w) np_ket(o,2w+sp)
  ! (and qk with bra/ket exchanged); ub(w,sr) = sum_c xb(c,w) hn_ket(2c+sr,o).
  allocate(qb(0:nc-1, 0:1), qk(0:nc-1, 0:1), ub(0:nv-1, 0:1), uk(0:nv-1, 0:1))
  do k = 0, 3
    pk = ishft(1, k)
    do o = 0, 15
      if (.not. has_tab(k, o)) then
        of = ior(o, pk)
        do sp = 0, 1
          qb(:, sp) = matmul(bra%xcv, ket%np(o, sp:npv-1:2))
        end do
        do sr = 0, 1
          uk(:, sr) = matmul(bra%hn(sr:nh-1:2, of), ket%xcv)
        end do
        ! O <- C
        do r = 0, nh - 1
          c = r/2
          sr = mod(r, 2)
          ph = parity_sign(after_tab(k, o) + r)
          s = bra%hn(r, of) * ket%nn(o) + sum(bra%hc(sr, of, :) * qb(c, :))
          d(o0 + k, r) = d(o0 + k, r) + ph * s
        end do
        ! O <- V
        ph = parity_sign(before_tab(k, o))
        do p = 0, npv - 1
          w = p/2
          sp = mod(p, 2)
          s = bra%nn(of) * ket%np(o, p) + sum(ket%hc(:, o, sp) * uk(w, :))
          d(o0 + k, v0 + p) = d(o0 + k, v0 + p) + ph * s
        end do
      else
        ox = ieor(o, pk)
        do sp = 0, 1
          qk(:, sp) = matmul(ket%xcv, bra%np(ox, sp:npv-1:2))
        end do
        do sr = 0, 1
          ub(:, sr) = matmul(ket%hn(sr:nh-1:2, o), bra%xcv)
        end do
        ! C <- O
        do r = 0, nh - 1
          c = r/2
          sr = mod(r, 2)
          ph = parity_sign(after_tab(k, o) + r)
          s = bra%nn(ox) * ket%hn(r, o) + sum(ket%hc(sr, o, :) * qk(c, :))
          d(r, o0 + k) = d(r, o0 + k) + ph * s
        end do
        ! V <- O
        ph = parity_sign(before_tab(k, o))
        do p = 0, npv - 1
          w = p/2
          sp = mod(p, 2)
          s = bra%np(ox, p) * ket%nn(o) + sum(bra%hc(:, ox, sp) * ub(w, :))
          d(v0 + p, o0 + k) = d(v0 + p, o0 + k) + ph * s
        end do
      end if
    end do
  end do
  deallocate(qb, qk, ub, uk)

  ! virtual <- core and core <- virtual: phase (-1)**(n_open + r) = (-1)**(n_open + sr)
  do sp = 0, 1
    do sr = 0, 1
      avc(sr, sp) = 0.0_dp
      acv(sr, sp) = 0.0_dp
      do o = 0, 15
        ph = parity_sign(pop_tab(o) + sr)
        avc(sr, sp) = avc(sr, sp) + ph * bra%hc(sr, o, sp) * ket%nn(o)
        acv(sr, sp) = acv(sr, sp) + ph * bra%nn(o) * ket%hc(sr, o, sp)
      end do
    end do
  end do
  do r = 0, nh - 1
    c = r/2
    sr = mod(r, 2)
    do p = 0, npv - 1
      w = p/2
      sp = mod(p, 2)
      d(v0 + p, r) = bra%xcv(c, w) * avc(sr, sp)
      d(r, v0 + p) = ket%xcv(c, w) * acv(sr, sp)
    end do
  end do
  deallocate(mvv, mcc)

end subroutine transition_density


subroutine compute_soc_matrix(t00aa, t110aa, t11ab, lx_mo, ly_mo, lz_mo, &
                               ns, nt, nbf, hsoc)
  !
  ! Assemble the full SOC Hamiltonian in the basis of MRSF spin-states.
  !
  ! OpenQP basis ordering:
  !   rows/cols 1..ns              : singlets S_I
  !   rows/cols ns+1..ns+3*nt      : triplets, grouped as
  !                                  (T_0,Ms=-1),(T_0,Ms=0),(T_0,Ms=+1),
  !                                  (T_1,Ms=-1),(T_1,Ms=0),(T_1,Ms=+1), ...
  !   index helper: itrp(J,Ms) = ns + (J-1)*3 + (Ms+2)
  !
  ! S-T block (only T00aa needed; T00bb = -T00aa by time reversal):
  !   <S_I| h_soc |T_J, Ms=0 > = sum_tu 2*celm_aa * T00aa(I,J,t,u)
  !   <S_I| h_soc |T_J, Ms=+1> = sum_tu celm_ba*(-sqrt2) * T00aa(I,J,t,u)
  !   <S_I| h_soc |T_J, Ms=-1> = sum_tu celm_ab*(+sqrt2) * T00aa(I,J,t,u)
  !
  ! T-T block (needs T110aa and T11ab; derived TDMs computed on the fly):
  !   T111aa   =  sqrt2*T11ab + T110aa
  !   T111bb   = -sqrt2*T11ab + T110aa  (= Tm11m1aa)
  !   T110bb   =  T110aa
  !   T1m1ba   =  T11ab                 (= Tm11m1bb)
  !
  !   <T_I,Ms=0 | h_soc |T_J,Ms=+1> = sum_tu celm_ba * T11ab(I,J,t,u)
  !   <T_I,Ms=0 | h_soc |T_J,Ms=-1> = sum_tu celm_ab * T11ab(I,J,t,u)
  !   <T_I,Ms=+1| h_soc |T_J,Ms=+1> = sum_tu (celm_aa*T111aa + celm_bb*T111bb)
  !   <T_I,Ms=0 | h_soc |T_J,Ms=0 > = sum_tu (celm_aa*T110aa + celm_bb*T110bb)
  !   <T_I,Ms=-1| h_soc |T_J,Ms=-1> = sum_tu (celm_aa*Tm11m1aa + celm_bb*Tm11m1bb)
  !
  ! Spin matrix elements (spnfac absorbed):
  !   celm_aa = (0, -Lz(t,u)/2)
  !   celm_bb = -celm_aa = (0, +Lz(t,u)/2)
  !   celm_ba = (Ly(t,u)/4, -Lx(t,u)/4)   [S- component; spnfac=0.5 already absorbed]
  !   celm_ab = (Ly(t,u)/4,  Lx(t,u)/4)   [S+ component]
  !
  ! Note: celm_ba here = sqrt(0.5)*0.5*(Ly-iLx).
  !       The extra sqrt(0.5) from the spin ladder operator S- is the spnfac.
  !
  ! Output: hsoc(nstate, nstate), nstate = ns + 3*nt, in Hartree.
  !
  use precision, only: dp
  implicit none
 
  real(kind=dp),    intent(in)  :: t00aa (ns, nt, nbf, nbf)
  real(kind=dp),    intent(in)  :: t110aa(nt, nt, nbf, nbf)
  real(kind=dp),    intent(in)  :: t11ab (nt, nt, nbf, nbf)
  real(kind=dp),    intent(in)  :: lx_mo(nbf, nbf)
  real(kind=dp),    intent(in)  :: ly_mo(nbf, nbf)
  real(kind=dp),    intent(in)  :: lz_mo(nbf, nbf)
  integer,          intent(in)  :: ns, nt, nbf
 
  complex(kind=dp), intent(out) :: hsoc(ns + 3*nt, ns + 3*nt)
 
  integer       :: ist, jst, it, iu, itrp_i, itrp_j
  complex(kind=dp) :: celm_aa, celm_bb, celm_ba, celm_ab
  real(kind=dp)    :: t111aa, t111bb, tm11m1aa, tm11m1bb, t110bb
  real(kind=dp), parameter :: sq2   = sqrt(2.0_dp)
  real(kind=dp), parameter :: sq05  = sqrt(0.5_dp)  
  real(kind=dp), parameter :: half  = 0.5_dp
  real(kind=dp), parameter :: quart = 0.25_dp
 
  hsoc = cmplx(0.0_dp, 0.0_dp, kind=dp)
 
  do iu = 1, nbf
    do it = 1, nbf
 
      celm_aa = cmplx( 0.0_dp,          -lz_mo(it,iu)*half,  kind=dp)
      celm_bb = cmplx( 0.0_dp,          +lz_mo(it,iu)*half,  kind=dp)  ! = -celm_aa
      celm_ba = cmplx(+ly_mo(it,iu)*half, -lx_mo(it,iu)*half, kind=dp)
      celm_ab = cmplx(-ly_mo(it,iu)*half, -lx_mo(it,iu)*half, kind=dp)
 
      do ist = 1, ns
        do jst = 1, nt
          ! --- S-T block: row=ist (singlet), col=ns+(jst-1)*3+Ms+2 ---
          ! Ms=0:
          hsoc(ist, ns+(jst-1)*3+2) = hsoc(ist, ns+(jst-1)*3+2) &
            + 2.0_dp * celm_aa * t00aa(ist,jst,it,iu)
          ! Ms=+1:
          hsoc(ist, ns+(jst-1)*3+3) = hsoc(ist, ns+(jst-1)*3+3) &
            + celm_ba * (-sq2) * t00aa(ist,jst,it,iu)
          ! Ms=-1:
          hsoc(ist, ns+(jst-1)*3+1) = hsoc(ist, ns+(jst-1)*3+1) &
            + celm_ab * (+sq2) * t00aa(ist,jst,it,iu)
        end do
      end do
 
      do ist = 1, nt
        do jst = 1, nt
          ! Derived TDMs (computed on the fly):
          t111aa    =  sq05 * t11ab(ist,jst,it,iu) + t110aa(ist,jst,it,iu)
          t111bb    = -sq05 * t11ab(ist,jst,it,iu) + t110aa(ist,jst,it,iu)
          t110bb    =  t110aa(ist,jst,it,iu)
          tm11m1aa  =  t111bb
          tm11m1bb  =  t111aa
 
          itrp_i = ns + (ist-1)*3
          itrp_j = ns + (jst-1)*3
 
          ! <T_I,Ms=0| h_soc |T_J,Ms=+1>: celm_ba * T11ab
          hsoc(itrp_i+2, itrp_j+3) = hsoc(itrp_i+2, itrp_j+3) &
            + celm_ba * t11ab(ist,jst,it,iu)
 
          ! <T_I,Ms=0| h_soc |T_J,Ms=-1>: celm_ab * T11ab (T1m1ba = T11ab)
          hsoc(itrp_i+2, itrp_j+1) = hsoc(itrp_i+2, itrp_j+1) &
            + celm_ab * t11ab(ist,jst,it,iu)
 
          ! <T_I,Ms=+1| h_soc |T_J,Ms=+1>: celm_aa*T111aa + celm_bb*T111bb
          hsoc(itrp_i+3, itrp_j+3) = hsoc(itrp_i+3, itrp_j+3) &
            + celm_aa * t111aa + celm_bb * t111bb

          ! <Ti,Ms=+1|Tj,Ms=0> = conjg(celm_ba) * T11ab(jst,ist,it,iu)
          hsoc(itrp_i+3, itrp_j+2) = hsoc(itrp_i+3, itrp_j+2) + conjg(celm_ba) * t11ab(jst,ist,it,iu) 

          ! <T_I,Ms=0| h_soc |T_J,Ms=0>: celm_aa*T110aa + celm_bb*T110bb
          hsoc(itrp_i+2, itrp_j+2) = hsoc(itrp_i+2, itrp_j+2) &
            + celm_aa * t110aa(ist,jst,it,iu) + celm_bb * t110bb

          ! <Ti,Ms=-1|Tj,Ms=0> = conjg(celm_ab) * T11ab(jst,ist,it,iu)
          hsoc(itrp_i+1, itrp_j+2) = hsoc(itrp_i+1, itrp_j+2) + conjg(celm_ab) * t11ab(jst,ist,it,iu)
 
          ! <T_I,Ms=-1| h_soc |T_J,Ms=-1>: celm_aa*Tm11m1aa + celm_bb*Tm11m1bb
          hsoc(itrp_i+1, itrp_j+1) = hsoc(itrp_i+1, itrp_j+1) &
            + celm_aa * tm11m1aa + celm_bb * tm11m1bb
 
        end do
      end do
 
    end do
  end do
 
  ! Fill lower triangle by Hermitian conjugation (hsoc should be Hermitian):
  ! H(j,i) = conjg(H(i,j)) for all i>j
  do ist = 1, ns + 3*nt
    do jst = ist+1, ns + 3*nt
      hsoc(jst, ist) = conjg(hsoc(ist, jst))
    end do
  end do
 
end subroutine compute_soc_matrix

!> @brief Diagonalize the SOC Hamiltonian and return adiabatic eigenvalues/eigenvectors
!> @details
!>  Builds the full (ns+3*nt) x (ns+3*nt) complex Hermitian matrix:
!>    diagonal  = (E_I - E_0) * ha2wn  [excitation energies in cm-1]
!>    off-diag  = hsoc(I,J)  * dfac    [SOC couplings in cm-1]
!>  where E_0 = min(singlet_energies(1), triplet_energies(1)).
!>  Diagonalizes via LAPACK zheev. The eigenvectors evec are returned in
!>  column-major order and used for the state decomposition print.
!>
!>  Note: uses explicit integer(4) arguments for LP64 LAPACK compatibility.
!>
!> @param[in]  hsoc                SOC Hamiltonian matrix in Hartree (ns+3*nt x ns+3*nt)
!> @param[in]  singlet_energies    MRSF singlet excitation energies (Hartree, rel. ROHF)
!> @param[in]  triplet_energies    MRSF triplet excitation energies (Hartree, rel. ROHF)
!> @param[in]  e_ref               ROHF reference energy (Hartree)
!> @param[in]  ns, nt              Number of singlet/triplet states
!> @param[out] eval                Adiabatic SOC eigenvalues (cm-1, rel. lowest state)
!> @param[out] evec                SOC eigenvectors (complex, column = adiabat)
subroutine diag_soc(hsoc, singlet_energies, triplet_energies, e_ref, ns, nt, eval, evec)
  use precision, only: dp
  use mathlib_types, only: blas_int
  use messages, only: show_message, WITH_ABORT
  use physical_constants, only: ha2wn => HA_TO_WAVENUM, &
                                 FINE_STRUCTURE
  implicit none

  complex(kind=dp), intent(in)  :: hsoc(ns+3*nt, ns+3*nt)
  real(kind=dp),    intent(in)  :: singlet_energies(ns)
  real(kind=dp),    intent(in)  :: triplet_energies(nt)
  real(kind=dp),    intent(in)  :: e_ref
  integer,          intent(in)  :: ns, nt
  real(kind=dp),    intent(out) :: eval(ns+3*nt)
  complex(kind=dp), intent(out) :: evec(ns+3*nt, ns+3*nt)

  integer :: nstate, ist, i, j, ioff
  integer(blas_int) :: nstate_, lwork_, info
  real(kind=dp) :: e0
  complex(kind=dp), allocatable :: work(:)
  complex(kind=dp) :: work_query(1)
  real(kind=dp),    allocatable :: rwork(:)

  real(kind=dp), parameter :: dfac = FINE_STRUCTURE**2 / 2.0_dp * ha2wn

  nstate  = ns + 3*nt
  nstate_ = int(nstate, blas_int)

  allocate(rwork(3*nstate))

  ! --- 1. Scale off-diagonal SOC elements to cm-1, fill diagonal with excitation energies ---
  do j = 1, nstate
    do i = 1, nstate
      evec(i,j) = hsoc(i,j) * dfac
    end do
  end do

  e0 = min(singlet_energies(1), triplet_energies(1))

  do ist = 1, ns
    evec(ist,ist) = cmplx((singlet_energies(ist) - e0)*ha2wn, 0.0_dp, kind=dp)
  end do

  do ist = 1, nt
    do j = 1, 3   ! Ms = -1, 0, +1 components share the same energy
      ioff = ns + (ist-1)*3 + j
      evec(ioff,ioff) = cmplx((triplet_energies(ist) - e0)*ha2wn, 0.0_dp, kind=dp)
    end do
  end do

  ! --- 2. Diagonalize via LAPACK zheev ---
  ! Workspace query
  ! nstate_ is the number of spin-orbit states -- small, and the same
  ! small-eigenproblem regime eigen_blas_threads describes.
  block
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use eigen, only: eigen_blas_scope_enter, eigen_blas_scope_exit
    integer(c_int64_t) :: nb
    nb = eigen_blas_scope_enter(int(nstate_))
    call zheev('V', 'U', nstate_, evec, nstate_, eval, work_query, -1_blas_int, rwork, info)
    lwork_ = int(real(work_query(1)), blas_int)
    allocate(work(lwork_))

    call zheev('V', 'U', nstate_, evec, nstate_, eval, work, lwork_, rwork, info)
    call eigen_blas_scope_exit(nb)
  end block

  if (info /= 0) then
    call show_message('(A,I0)', 'ZHEEV failed in diag_soc, info=', int(info), WITH_ABORT)
  end if

  deallocate(rwork, work)

end subroutine diag_soc



!> @brief Print SOC adiabatic eigenvalues and eigenvector decomposition table
!> @details
!>  Writes two blocks to iw:
!>    1. Eigenvalue table: state index, energy in cm-1, Hartree, and eV
!>    2. Eigenvector table: mixing coefficients printed in blocks of 5 columns;
!>       components with |c|^2 < 0.01 are suppressed for readability.
!>
!> @param[in]  iw                 Log file unit
!> @param[in]  eval               Adiabatic SOC eigenvalues (cm-1)
!> @param[in]  evec               SOC eigenvectors (complex, column = adiabat)
!> @param[in]  singlet_energies   MRSF singlet excitation energies (Hartree)
!> @param[in]  triplet_energies   MRSF triplet excitation energies (Hartree)
!> @param[in]  e_ref              ROHF reference energy (Hartree)
!> @param[in]  ns, nt             Number of singlet/triplet states
subroutine print_soc_eigenvalues(iw, eval, evec, singlet_energies, triplet_energies, e_ref, ns, nt)
  use precision, only: dp
  implicit none

  integer,          intent(in) :: iw, ns, nt
  real(kind=dp),    intent(in) :: eval(ns+3*nt)
  complex(kind=dp), intent(in) :: evec(ns+3*nt, ns+3*nt)
  real(kind=dp),    intent(in) :: singlet_energies(ns), triplet_energies(nt), e_ref

  integer :: nstate, ist, i, j, ncols, ioff, ms_idx
  real(kind=dp) :: e0, a, b, tmpmod
  real(kind=dp), parameter :: ha2wn = 219474.6_dp
  real(kind=dp), parameter :: ha2ev = 27.211386245988_dp

  nstate = ns + 3*nt
  e0     = min(singlet_energies(1), triplet_energies(1))

  ! Eigenvalues
  write(iw,'(/,11x,65("-"))')
  write(iw,'(11x,a)') 'SOC eigenvalues (adiabatic, SOC-corrected)'
  write(iw,'(11x,65("-"))')
  write(iw,'(a)') '  Non-SOC ground state is at 0 cm-1'
  write(iw,'()')
  write(iw,'(5x,a,12x,a,14x,a,14x,a)') 'State', 'cm-1', 'Hartree', 'eV'
  do ist = 1, nstate
    write(iw,'(5x,i4,2x,f14.4,2x,f18.10,2x,f14.6)') &
      ist, eval(ist), eval(ist)/ha2wn + e0, (eval(ist)/ha2wn + e0)*ha2ev
  end do

  ! Eigenvectors (5 columns at a time)
  write(iw,'(/,11x,a)') 'SOC eigenvectors (rows = diabatic states, cols = adiabats)'
  ncols = 5
  ioff  = 0
  do while (ioff < nstate)
    write(iw,'(/,10x)', advance='no')
    do j = ioff+1, min(ioff+ncols, nstate)
      write(iw,'(i6,8x)', advance='no') j
    end do
    write(iw,*)
    write(iw,'(8x,a)', advance='no') 'E(cm-1)'
    do j = ioff+1, min(ioff+ncols, nstate)
      write(iw,'(f10.2,4x)', advance='no') eval(j)
    end do
    write(iw,*)
    do i = 1, nstate
      if (i <= ns) then
        write(iw,'(4x,a1,i3,2x)', advance='no') 'S', i-1
      else
        ist    = (i - ns - 1)/3 + 1
        ms_idx = mod(i - ns - 1, 3)
        select case(ms_idx)
          case(0); write(iw,'(3x,a1,i3,a)', advance='no') 'T', ist-1, '(-1)'
          case(1); write(iw,'(3x,a1,i3,a)', advance='no') 'T', ist-1, '( 0)'
          case(2); write(iw,'(3x,a1,i3,a)', advance='no') 'T', ist-1, '(+1)'
        end select
      end if
      do j = ioff+1, min(ioff+ncols, nstate)
        a = real(evec(i,j)); b = aimag(evec(i,j))
        tmpmod = a**2 + b**2
        if (tmpmod >= 0.01_dp) then
          write(iw,'(f8.4,sp,f7.4,"i"," ")', advance='no') a, b
        else
          write(iw,'(16x)', advance='no')
        end if
      end do
      write(iw,*)
    end do
    ioff = ioff + ncols
  end do

end subroutine print_soc_eigenvalues

! 2e part

!> @brief Print per-state SOC decomposition: energy and top-3 diabatic contributions
!> @details
!>  For each adiabatic SOC state, computes the diabatic weights |c_I|^2 from
!>  the eigenvector matrix and identifies the three largest contributions.
!>  Output format (one line per state):
!>    index   cm-1   eV   Label1 (wt%)   Label2 (wt%)   Label3 (wt%)
!>  where labels are S<n> for singlets and T<n>(Ms) for triplet sublevels.
!>  Contributions below 0.1% are suppressed.
!>
!> @param[in]  iw    Log file unit
!> @param[in]  eval  Adiabatic SOC eigenvalues (cm-1)
!> @param[in]  evec  SOC eigenvectors (complex, column = adiabat)
!> @param[in]  ns    Number of singlet states
!> @param[in]  nt    Number of triplet states
subroutine print_soc_decomposition(iw, eval, evec, ns, nt)
  use precision, only: dp
  implicit none

  integer,          intent(in) :: iw, ns, nt
  real(kind=dp),    intent(in) :: eval(ns+3*nt)
  complex(kind=dp), intent(in) :: evec(ns+3*nt, ns+3*nt)

  integer, parameter :: ntop = 3
  integer  :: nstate, ist, i, j, ms_idx
  real(kind=dp) :: weight(ns+3*nt), tmpmod
  real(kind=dp), parameter :: ha2wn = 219474.6_dp
  real(kind=dp), parameter :: ha2ev = 27.211386245988_dp
  integer  :: idx(ntop)
  real(kind=dp) :: best(ntop)
  character(len=10) :: label

  nstate = ns + 3*nt

  write(iw,'(/,11x,65("-"))')
  write(iw,'(11x,a)') 'SOC state decomposition (top 3 diabatic contributions)'
  write(iw,'(11x,65("-"))')
  write(iw,'(/,2x,a,4x,a,8x,a,8x,a)') 'State', 'cm-1', 'eV', 'Composition'
  write(iw,*)

  do ist = 1, nstate
    ! compute weights
    do i = 1, nstate
      tmpmod = real(evec(i,ist))**2 + aimag(evec(i,ist))**2
      weight(i) = tmpmod
    end do

    ! find top-3 by simple selection
    idx  = 0
    best = -1.0_dp
    do j = 1, ntop
      do i = 1, nstate
        if (weight(i) > best(j)) then
          if (j == 1 .or. all(idx(1:j-1) /= i)) then
            best(j) = weight(i)
            idx(j)  = i
          end if
        end if
      end do
      ! mask already chosen
      if (idx(j) > 0) weight(idx(j)) = -1.0_dp
    end do

    ! restore weights for next iteration
    do i = 1, nstate
      tmpmod = real(evec(i,ist))**2 + aimag(evec(i,ist))**2
      weight(i) = tmpmod
    end do

    write(iw,'(2x,i4,2x,f10.2,2x,f10.6,2x)', advance='no') &
      ist, eval(ist), eval(ist)/ha2wn * ha2ev

    do j = 1, ntop
      i = idx(j)
      if (i == 0) exit
      if (weight(i) < 0.001_dp) exit
      if (i <= ns) then
        write(label,'(a1,i0)') 'S', i-1
      else
        ms_idx = mod(i - ns - 1, 3)
        select case(ms_idx)
          case(0); write(label,'(a1,i0,a)') 'T', (i-ns-1)/3, '(-1)'
          case(1); write(label,'(a1,i0,a)') 'T', (i-ns-1)/3, '( 0)'
          case(2); write(label,'(a1,i0,a)') 'T', (i-ns-1)/3, '(+1)'
        end select
      end if
      write(iw,'(a8,a,f5.1,a,a)', advance='no') &
        trim(label), ' (', weight(i)*100.0_dp, '%)', '   '
    end do
    write(iw,*)
  end do

  write(iw,'(11x,65("-"))')

end subroutine print_soc_decomposition





end module soc_mrsf_mod
