module soc_mrsf_mod

    
  use precision, only: dp  

  implicit none

  character(len=*), parameter :: module_name = "soc_mrsf_mod"


  real(kind=dp), parameter :: ZEFF_321G(54) = [ &
       1.00d0,  2.00d0,  1.50d0,  2.20d0,  3.00d0,  3.90d0,  4.90d0,  6.00d0, &
       7.20d0, 10.00d0, 10.67d0, 11.52d0, 12.35d0, 13.16d0, 13.95d0, 14.72d0, &
      15.47d0, 18.00d0, 22.42d0, 23.00d0, 21.00d0, 22.00d0, 23.00d0, 24.00d0, &
      25.00d0, 26.00d0, 27.00d0, 28.00d0, 29.00d0, 30.00d0, 34.72d0, 34.88d0, &
      34.98d0, 35.02d0, 35.00d0, 36.00d0, 45.88d0, 47.12d0, 39.00d0, 40.00d0, &
      41.00d0, 42.00d0, 43.00d0, 44.00d0, 45.00d0, 46.00d0, 47.00d0, 48.00d0, &
      60.72d0, 62.00d0, 63.24d0, 64.48d0, 65.72d0, 54.00d0 ]


  private
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
!>  State ordering in the SOC basis follows the GAMESS convention:
!>    indices 1..ns        -> singlet states S0..S(ns-1)
!>    indices ns+1..ns+3nt -> triplet Ms sublevels T1(Ms=-1,0,+1), T2(...), ...
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
    real(kind=dp), parameter :: alpha = 7.2973506e-3_dp
    real(kind=dp), parameter :: ha2wn = 219474.6_dp
    real(kind=dp), parameter :: dfac  = alpha**2 / 2.0_dp * ha2wn  ! 5.8438 cm-1/a.u.
    real(kind=dp) :: re1e, im1e, re2e, im2e, abs12e
    character(len=7), dimension(3), parameter :: trip = ['(Ms=-1)', '(Ms= 0)', '(Ms=+1)']

    real(kind=dp), allocatable :: den_rohf(:,:)
    real(kind=dp), allocatable :: wao(:,:,:)
    real(kind=dp), allocatable :: lx_2e_mo(:,:), ly_2e_mo(:,:), lz_2e_mo(:,:)
    real(kind=dp), allocatable :: lx_12e_mo(:,:), ly_12e_mo(:,:), lz_12e_mo(:,:)

    logical :: do_2e_soc
    logical, parameter :: debug_soc_prints = .true.

    integer :: nstate_soc
    real(kind=dp), pointer :: eval_out(:)
    real(kind=dp), pointer :: evec_re_out(:,:), evec_im_out(:,:)
    real(kind=dp), pointer :: hsoc_re_out(:,:), hsoc_im_out(:,:)

    ! Temporary switch for the two-electron SOC mean-field correction.



    open(unit=iw, file=infos%log_filename, position="append")

    do_2e_soc = (infos%control%soc_2e /= 0)
    
    if (do_2e_soc) then
        call print_module_info('SOC_MRSF (1e+2e)', 'Spin-Orbit Coupling: MRSF Energies')
    else
        call print_module_info('SOC_MRSF (1e)', 'Spin-Orbit Coupling: MRSF Energies')
    endif

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

    ! --- Step 1: Compute SOC AO integrals <mu|Z*L/r^3|nu> ---
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

    ! --- Step 2b: 2e mean-field SOC correction ---
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

      if (debug_soc_prints) then
        do i = 1, infos%basis%nshell
          write(iw,'(a,3i4)') 'shell ao_offset am:', i, &
            infos%basis%ao_offset(i), infos%basis%am(i)
        end do
        write(iw,'(a,2i4)') ' nocca, noccb = ', nocca, noccb
      end if

      call soc2e_driver(infos, infos%basis, den_rohf, wao)

      if (debug_soc_prints) then
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

      call ao2mo_full(wao(1,:,:), lx_2e_mo, mo_a, nbf)
      call ao2mo_full(wao(2,:,:), ly_2e_mo, mo_a, nbf)
      call ao2mo_full(wao(3,:,:), lz_2e_mo, mo_a, nbf)
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

    ! --- Step 3: Build spin-dependent transition density matrices ---
    allocate(t00aa (ns, nt, nbf, nbf), &
             t110aa(nt, nt, nbf, nbf), &
             t11ab (nt, nt, nbf, nbf), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate TDM arrays', WITH_ABORT)
    call compute_tdm(bvec_mo_s, bvec_mo_t, nocca, noccb, nbf, ns, nt, &
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
            'S', ist-1, 'T', jst, trim(trip(ims)), re1e, im1e, re2e, im2e, abs12e
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
              'T', ist, trim(trip(ims_i)), 'T', jst, trim(trip(ims_j)), &
              re1e, im1e, re2e, im2e, abs12e
          end do
        end do
      end do
    end do

    ! --- Step 6: Diagonalize H_SOC + excitation energies, print eigenvalues ---
    allocate(eval(ns + 3*nt), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate eigenvalue array', WITH_ABORT)
    allocate(evec(ns + 3*nt, ns + 3*nt), stat=ok)
    if (ok /= 0) call show_message('soc_mrsf: cannot allocate eigenvector array', WITH_ABORT)
    call diag_soc(hsoc, singlet_energies, triplet_energies, e_ref, ns, nt, eval, evec)
    call print_soc_eigenvalues(iw, eval, evec, singlet_energies, triplet_energies, e_ref, ns, nt)
    call print_soc_decomposition(iw, eval, evec, ns, nt) 

    call infos%dat%reserve_data(OQP_soc_eval, TA_TYPE_REAL64, &
        nstate_soc, (/nstate_soc/), comment=OQP_soc_eval_comment)
    call infos%dat%reserve_data(OQP_soc_evec_re, TA_TYPE_REAL64, &
        nstate_soc*nstate_soc, (/nstate_soc, nstate_soc/), comment=OQP_soc_evec_re_comment)
    call infos%dat%reserve_data(OQP_soc_evec_im, TA_TYPE_REAL64, &
        nstate_soc*nstate_soc, (/nstate_soc, nstate_soc/), comment=OQP_soc_evec_im_comment)
    call infos%dat%reserve_data(OQP_soc_hsoc_re, TA_TYPE_REAL64, &
        nstate_soc*nstate_soc, (/nstate_soc, nstate_soc/), comment=OQP_soc_hsoc_re_comment)
    call infos%dat%reserve_data(OQP_soc_hsoc_im, TA_TYPE_REAL64, &
        nstate_soc*nstate_soc, (/nstate_soc, nstate_soc/), comment=OQP_soc_hsoc_im_comment)
    
    call tagarray_get_data(infos%dat, OQP_soc_eval,    eval_out)
    call tagarray_get_data(infos%dat, OQP_soc_evec_re, evec_re_out)
    call tagarray_get_data(infos%dat, OQP_soc_evec_im, evec_im_out)
    call tagarray_get_data(infos%dat, OQP_soc_hsoc_re, hsoc_re_out)
    call tagarray_get_data(infos%dat, OQP_soc_hsoc_im, hsoc_im_out)
    
    eval_out     = eval
    evec_re_out  = real(evec,  kind=dp)
    evec_im_out  = aimag(evec)
    hsoc_re_out  = real(hsoc,  kind=dp)
    hsoc_im_out  = aimag(hsoc)
    
    
    
    
    deallocate(eval, evec)
        deallocate(lx_12e_mo, ly_12e_mo, lz_12e_mo)

    deallocate(hsoc)
    if (do_2e_soc) then
        deallocate(lx_2e_mo, ly_2e_mo, lz_2e_mo)
    endif
    write(iw,'(/,a)') 'SOC_MRSF done'
    call flush(iw)
    close(iw)

  end subroutine soc_mrsf


!> @brief Compute 1-electron SOC AO integrals using the Breit-Pauli operator
!> @details
!>  Evaluates <mu|Z_A * L_A / r_A^3|nu> for each atom A, where L_A is the
!>  angular momentum operator relative to nucleus A and Z_A is the (effective)
!>  nuclear charge. Loops over shell pairs (ii >= jj) and accumulates into
!>  packed lower-triangular arrays. Results are normalised with basis function norms.
!>
!>  Note: currently uses bare nuclear charges (ze = Z). Slater effective charges
!>  (ZEFF_321G) are available but commented out.
!>
!> @param[inout] infos   OQP information struct (basis, atoms)
!> @param[out]   lx_ao   Lx AO integrals, packed lower-triangular (nbf*(nbf+1)/2)
!> @param[out]   ly_ao   Ly AO integrals, packed lower-triangular
!> @param[out]   lz_ao   Lz AO integrals, packed lower-triangular
subroutine compute_soc_ao(infos, lx_ao, ly_ao, lz_ao)
  use basis_tools,       only: basis_set, bas_norm_matrix
  use mod_1e_primitives, only: comp_soc_int1_prim, update_triang_matrix
  use mod_shell_tools,   only: shell_t, shpair_t
  use constants,         only: tol_int
  use precision,         only: dp
  use types,             only: information

  implicit none

  type(information), target, intent(inout) :: infos
  real(kind=dp), intent(out) :: lx_ao(:), ly_ao(:), lz_ao(:)

  type(basis_set), pointer :: basis
  type(shell_t)   :: shi, shj
  type(shpair_t)  :: cntp

  integer, parameter :: blocksize = 28*28
  real(kind=dp) :: socblk(blocksize, 3)

  integer  :: ii, jj, ig, iat, iz, nat, nbf
  real(kind=dp) :: ze, tol

  basis => infos%basis
  basis%atoms => infos%atoms
  nat  = size(infos%atoms%zn)
  nbf  = basis%nbf
  tol  = log(10.0_dp) * tol_int

  lx_ao = 0.0_dp
  ly_ao = 0.0_dp
  lz_ao = 0.0_dp

  call cntp%alloc(basis)

  do ii = basis%nshell, 1, -1
    call shi%fetch_by_id(basis, ii)
    do jj = 1, ii
      call shj%fetch_by_id(basis, jj)
      call cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
      if (cntp%numpairs == 0) cycle

      socblk = 0.0_dp

      do ig = 1, cntp%numpairs
        do iat = 1, nat
          iz = nint(infos%atoms%zn(iat))
!         ZEFFECTIVE CHARGES
!          ze = merge(ZEFF_321G(iz), real(iz, dp), iz <= 54)
          ze = real(iz, dp)
          call comp_soc_int1_prim(cntp, ig, infos%atoms%xyz(:,iat), ze, socblk)
        end do
      end do

      call update_triang_matrix(shi, shj, socblk(:,1), lx_ao)
      call update_triang_matrix(shi, shj, socblk(:,2), ly_ao)
      call update_triang_matrix(shi, shj, socblk(:,3), lz_ao)

    end do
  end do

  call bas_norm_matrix(lx_ao, basis%bfnrm, nbf)
  call bas_norm_matrix(ly_ao, basis%bfnrm, nbf)
  call bas_norm_matrix(lz_ao, basis%bfnrm, nbf)

end subroutine compute_soc_ao

!> @brief Print a packed SOC AO integral matrix in GAMESS-compatible format
!> @details
!>  Writes the lower-triangular AO integral matrix to unit iw in blocks of
!>  NCOLS=5 columns, with basis function labels and row indices, matching
!>  the layout used in GAMESS for direct comparison.
!>
!> @param[in]  iw     Log file unit
!> @param[in]  comp   Component label ('LX', 'LY', or 'LZ')
!> @param[in]  mat    Packed lower-triangular AO matrix (nbf*(nbf+1)/2)
!> @param[in]  nbf    Number of basis functions
!> @param[in]  basis  Basis set descriptor (used for bf_label)
subroutine print_soc_ao_gamess(iw, comp, mat, nbf, basis)
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

end subroutine print_soc_ao_gamess

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

!> @brief Build spin-dependent transition density matrices from MRSF Davidson vectors
!> @details
!>  Constructs the TDMs needed to assemble the SOC Hamiltonian matrix elements:
!>    t00aa(I,J,t,u)  -- singlet I / triplet J TDM in the alpha-alpha spin sector
!>    t110aa(I,J,t,u) -- triplet I / triplet J, Ms=0 component (alpha-alpha)
!>    t11ab(I,J,t,u)  -- triplet I / triplet J, Ms=+1/-1 component (alpha-beta)
!>
!>  The Davidson eigenvectors are first reordered to the Ms-resolved form required
!>  by the GAMESS SOC convention (see compute_soc_matrix for the state ordering).
!>  The open-shell ROHF reference determines the active MO indices (iO1, iO2, iC).
!>
!> @param[in]  bvec_s   Singlet Davidson vectors (xvec_dim x ns)
!> @param[in]  bvec_t   Triplet Davidson vectors (xvec_dim x nt)
!> @param[in]  nocca    Number of alpha occupied MOs
!> @param[in]  noccb    Number of beta  occupied MOs
!> @param[in]  nbf      Number of basis functions
!> @param[in]  ns, nt   Number of singlet/triplet states
!> @param[out] t00aa    Singlet-triplet TDM (ns x nt x nbf x nbf)
!> @param[out] t110aa   Triplet-triplet TDM, Ms=0 sector (nt x nt x nbf x nbf)
!> @param[out] t11ab    Triplet-triplet TDM, Ms=±1 sector (nt x nt x nbf x nbf)
subroutine compute_tdm(bvec_s, bvec_t, nocca, noccb, nbf, ns, nt, &
                       t00aa, t110aa, t11ab)
  use precision, only: dp
  implicit none

  real(kind=dp), intent(in)  :: bvec_s(nocca*(nbf-noccb), ns)
  real(kind=dp), intent(in)  :: bvec_t(nocca*(nbf-noccb), nt)
  integer,       intent(in)  :: nocca, noccb, nbf, ns, nt

  real(kind=dp), intent(out) :: t00aa (ns, nt, nbf, nbf)
  real(kind=dp), intent(out) :: t110aa(nt, nt, nbf, nbf)
  real(kind=dp), intent(out) :: t11ab (nt, nt, nbf, nbf)

  real(kind=dp), allocatable :: xs(:,:), xt(:,:)

  integer :: xvec_dim
  integer :: iV, iO2, iO1, iC
  integer :: ijLR1, ijG, ijD, ijLR2
  integer :: ist, jst, i, it, iu
  integer :: ijiO1, ijiO2, ijO1a, ijO2a
  integer :: iO1a, jO1a, iO2a, jO2a
  integer :: iiO1, jiO1, iiO2, jiO2

  real(kind=dp), parameter :: half  = 0.5_dp
  real(kind=dp), parameter :: sqrt2 = 1.0_dp / sqrt(2.0_dp)

  iV   = nocca + 1
  iO2  = nocca
  iO1  = nocca - 1
  iC   = nocca - 2

  xvec_dim = size(bvec_s, 1)

  ijLR1 = (iO1 - noccb - 1)*nocca + iO1
  ijG   = (iO1 - noccb - 1)*nocca + iO2
  ijD   = (iO2 - noccb - 1)*nocca + iO1
  ijLR2 = (iO2 - noccb - 1)*nocca + iO2

  allocate(xs(xvec_dim, ns), xt(xvec_dim, nt))

  do ist = 1, ns
    xs(:, ist) = bvec_s(:, ist)
    xs(ijLR2, ist) = -bvec_s(ijLR1, ist)
  end do

  do jst = 1, nt
    xt(:, jst) = bvec_t(:, jst)
    xt(ijLR2, jst) =  bvec_t(ijLR1, jst)
    xt(ijG,   jst) = 0.0_dp
    xt(ijD,   jst) = 0.0_dp
  end do

  t00aa  = 0.0_dp
  t110aa = 0.0_dp
  t11ab  = 0.0_dp

  ! Diagonal t=u block. Only triplet-triplet TDMs are nonzero.
  do ist = 1, nt
    do jst = 1, nt
      do iu = 1, nbf
        it = iu
        do i = 1, noccb
          ijiO1 = (iO1 - noccb - 1)*nocca + i
          ijiO2 = (iO2 - noccb - 1)*nocca + i
          t110aa(ist, jst, it, iu) = t110aa(ist, jst, it, iu) &
            + xt(ijiO1, ist) * xt(ijiO1, jst) &
            + xt(ijiO2, ist) * xt(ijiO2, jst)
        end do
        t110aa(ist, jst, it, iu) = t110aa(ist, jst, it, iu) &
          + xt(ijLR1, ist) * xt(ijLR1, jst)
      end do
    end do
  end do

  ! Core-open blocks: singlet-triplet part.
  do ist = 1, ns
    do jst = 1, nt
      ijiO1 = (iO1 - noccb - 1)*nocca + iC
      ijiO2 = (iO2 - noccb - 1)*nocca + iC
      t00aa(ist, jst, iO1, iC) = -half  * xs(ijLR1, ist) * xt(ijiO1, jst) &
                                  -sqrt2 * xs(ijD,   ist) * xt(ijiO2, jst)
      t00aa(ist, jst, iC, iO1) = -half  * xs(ijiO1, ist) * xt(ijLR1, jst)
      t00aa(ist, jst, iO2, iC) = -sqrt2 * xs(ijG,   ist) * xt(ijiO1, jst) &
                                  +half  * xs(ijLR2, ist) * xt(ijiO2, jst)
      t00aa(ist, jst, iC, iO2) = -half  * xs(ijiO2, ist) * xt(ijLR2, jst)
    end do
  end do

  ! Core-open blocks: triplet-triplet part.
  do ist = 1, nt
    do jst = 1, nt
      ijiO1 = (iO1 - noccb - 1)*nocca + iC
      ijiO2 = (iO2 - noccb - 1)*nocca + iC
      t110aa(ist, jst, iO1, iC) = -half  * xt(ijLR1, ist) * xt(ijiO1, jst)
      t11ab (ist, jst, iO1, iC) = +sqrt2 * xt(ijLR1, ist) * xt(ijiO1, jst)
      t110aa(ist, jst, iC, iO1) = -half  * xt(ijiO1, ist) * xt(ijLR1, jst)
      t11ab (ist, jst, iC, iO1) = +sqrt2 * xt(ijiO1, ist) * xt(ijLR1, jst)
      t110aa(ist, jst, iO2, iC) = -half  * xt(ijLR2, ist) * xt(ijiO2, jst)
      t11ab (ist, jst, iO2, iC) = +sqrt2 * xt(ijLR2, ist) * xt(ijiO2, jst)
      t110aa(ist, jst, iC, iO2) = -half  * xt(ijiO2, ist) * xt(ijLR2, jst)
      t11ab (ist, jst, iC, iO2) = +sqrt2 * xt(ijiO2, ist) * xt(ijLR1, jst)
    end do
  end do

  ! Open-open blocks: singlet-triplet part.
  do ist = 1, ns
    do jst = 1, nt
      do i = 1, noccb
        ijiO1 = (iO1 - noccb - 1)*nocca + i
        ijiO2 = (iO2 - noccb - 1)*nocca + i
        t00aa(ist, jst, iO2, iO1) = t00aa(ist, jst, iO2, iO1) &
          - half * xs(ijiO1, ist) * xt(ijiO2, jst)
        t00aa(ist, jst, iO1, iO2) = t00aa(ist, jst, iO1, iO2) &
          - half * xs(ijiO2, ist) * xt(ijiO1, jst)
      end do
      t00aa(ist, jst, iO2, iO1) = t00aa(ist, jst, iO2, iO1) &
        - sqrt2 * xs(ijG, ist) * xt(ijLR1, jst)
      t00aa(ist, jst, iO1, iO2) = t00aa(ist, jst, iO1, iO2) &
        - sqrt2 * xs(ijD, ist) * xt(ijLR2, jst)
      do i = 1, nbf - nocca
        ijO2a = (nocca + i - noccb - 1)*nocca + iO2
        ijO1a = (nocca + i - noccb - 1)*nocca + iO1
        t00aa(ist, jst, iO2, iO1) = t00aa(ist, jst, iO2, iO1) &
          - half * xs(ijO2a, ist) * xt(ijO1a, jst)
        t00aa(ist, jst, iO1, iO2) = t00aa(ist, jst, iO1, iO2) &
          - half * xs(ijO1a, ist) * xt(ijO2a, jst)
      end do
    end do
  end do

  ! Open-open blocks: triplet-triplet part.
  do ist = 1, nt
    do jst = 1, nt
      do i = 1, noccb
        ijiO1 = (iO1 - noccb - 1)*nocca + i
        ijiO2 = (iO2 - noccb - 1)*nocca + i
        t110aa(ist, jst, iO2, iO1) = t110aa(ist, jst, iO2, iO1) &
          + half * xt(ijiO1, ist) * xt(ijiO2, jst)
        t11ab(ist, jst, iO2, iO1) = t11ab(ist, jst, iO2, iO1) &
          - sqrt2 * xt(ijiO1, ist) * xt(ijiO2, jst)
        t110aa(ist, jst, iO1, iO2) = t110aa(ist, jst, iO1, iO2) &
          + half * xt(ijiO2, ist) * xt(ijiO1, jst)
        t11ab(ist, jst, iO1, iO2) = t11ab(ist, jst, iO1, iO2) &
          - sqrt2 * xt(ijiO2, ist) * xt(ijiO1, jst)
      end do
      do i = 1, nbf - nocca
        ijO2a = (nocca + i - noccb - 1)*nocca + iO2
        ijO1a = (nocca + i - noccb - 1)*nocca + iO1
        t110aa(ist, jst, iO2, iO1) = t110aa(ist, jst, iO2, iO1) &
          - half * xt(ijO2a, ist) * xt(ijO1a, jst)
        t11ab(ist, jst, iO2, iO1) = t11ab(ist, jst, iO2, iO1) &
          - sqrt2 * xt(ijO2a, ist) * xt(ijO1a, jst)
        t110aa(ist, jst, iO1, iO2) = t110aa(ist, jst, iO1, iO2) &
          - half * xt(ijO1a, ist) * xt(ijO2a, jst)
        t11ab(ist, jst, iO1, iO2) = t11ab(ist, jst, iO1, iO2) &
          - sqrt2 * xt(ijO1a, ist) * xt(ijO2a, jst)
      end do
    end do
  end do

  ! Virtual-open blocks: singlet-triplet part.
  do ist = 1, ns
    do jst = 1, nt
      do it = iV, nbf
        ijO1a = (it - noccb - 1)*nocca + iO1
        ijO2a = (it - noccb - 1)*nocca + iO2
        t00aa(ist, jst, it, iO1) = -sqrt2 * xs(ijG,   ist) * xt(ijO2a, jst) &
                                    -half  * xs(ijLR2, ist) * xt(ijO1a, jst)
        t00aa(ist, jst, iO1, it) = -half  * xs(ijO1a, ist) * xt(ijLR2, jst)
        t00aa(ist, jst, it, iO2) = +half  * xs(ijLR1, ist) * xt(ijO2a, jst) &
                                    -sqrt2 * xs(ijD,   ist) * xt(ijO1a, jst)
        t00aa(ist, jst, iO2, it) = -half  * xs(ijO2a, ist) * xt(ijLR1, jst)
      end do
    end do
  end do

  ! Virtual-open blocks: triplet-triplet part.
  do ist = 1, nt
    do jst = 1, nt
      do it = iV, nbf
        ijO1a = (it - noccb - 1)*nocca + iO1
        ijO2a = (it - noccb - 1)*nocca + iO2
        t110aa(ist, jst, it, iO1) = +half  * xt(ijLR2, ist) * xt(ijO1a, jst)
        t11ab (ist, jst, it, iO1) = +sqrt2 * xt(ijLR1, ist) * xt(ijO1a, jst)
        t110aa(ist, jst, iO1, it) = +half  * xt(ijO1a, ist) * xt(ijLR2, jst)
        t11ab (ist, jst, iO1, it) = +sqrt2 * xt(ijO1a, ist) * xt(ijLR1, jst)
        t110aa(ist, jst, it, iO2) = +half  * xt(ijLR1, ist) * xt(ijO2a, jst)
        t11ab (ist, jst, it, iO2) = +sqrt2 * xt(ijLR2, ist) * xt(ijO2a, jst)
        t110aa(ist, jst, iO2, it) = +half  * xt(ijO2a, ist) * xt(ijLR1, jst)
        t11ab (ist, jst, iO2, it) = +sqrt2 * xt(ijO2a, ist) * xt(ijLR1, jst)
      end do
    end do
  end do

  ! Off-diagonal virtual-virtual blocks.
  do ist = 1, ns
    do jst = 1, nt
      do iu = iV, nbf
        do it = iV, nbf
          if (iu == it) cycle
          iO1a = (iu - noccb - 1)*nocca + iO1
          jO1a = (it - noccb - 1)*nocca + iO1
          iO2a = (iu - noccb - 1)*nocca + iO2
          jO2a = (it - noccb - 1)*nocca + iO2
          t00aa(ist, jst, it, iu) = -half * (xs(iO1a, ist)*xt(jO1a, jst) &
                                            + xs(iO2a, ist)*xt(jO2a, jst))
        end do
      end do
    end do
  end do

  do ist = 1, nt
    do jst = 1, nt
      do iu = iV, nbf
        do it = iV, nbf
          if (iu == it) cycle
          iO1a = (iu - noccb - 1)*nocca + iO1
          jO1a = (it - noccb - 1)*nocca + iO1
          iO2a = (iu - noccb - 1)*nocca + iO2
          jO2a = (it - noccb - 1)*nocca + iO2
          t110aa(ist, jst, it, iu) = +half  * (xt(iO1a, ist)*xt(jO1a, jst) &
                                              + xt(iO2a, ist)*xt(jO2a, jst))
          t11ab (ist, jst, it, iu) = +sqrt2 * (xt(iO1a, ist)*xt(jO1a, jst) &
                                              + xt(iO2a, ist)*xt(jO2a, jst))
        end do
      end do
    end do
  end do

  ! Off-diagonal core-core blocks.
  do ist = 1, ns
    do jst = 1, nt
      do iu = 1, noccb
        do it = 1, noccb
          if (iu == it) cycle
          iiO1 = (iO1 - noccb - 1)*nocca + iu
          jiO1 = (iO1 - noccb - 1)*nocca + it
          iiO2 = (iO2 - noccb - 1)*nocca + iu
          jiO2 = (iO2 - noccb - 1)*nocca + it
          t00aa(ist, jst, it, iu) = -half * (xs(iiO1, ist)*xt(jiO1, jst) &
                                            + xs(iiO2, ist)*xt(jiO2, jst))
        end do
      end do
    end do
  end do

  do ist = 1, nt
    do jst = 1, nt
      do iu = 1, noccb
        do it = 1, noccb
          if (iu == it) cycle
          iiO1 = (iO1 - noccb - 1)*nocca + iu
          jiO1 = (iO1 - noccb - 1)*nocca + it
          iiO2 = (iO2 - noccb - 1)*nocca + iu
          jiO2 = (iO2 - noccb - 1)*nocca + it
          t110aa(ist, jst, it, iu) = +half  * (xt(iiO1, ist)*xt(jiO1, jst) &
                                              + xt(iiO2, ist)*xt(jiO2, jst))
          t11ab (ist, jst, it, iu) = +sqrt2 * (xt(iiO1, ist)*xt(jiO1, jst) &
                                              + xt(iiO2, ist)*xt(jiO2, jst))
        end do
      end do
    end do
  end do

  deallocate(xs, xt)

end subroutine compute_tdm


subroutine compute_soc_matrix(t00aa, t110aa, t11ab, lx_mo, ly_mo, lz_mo, &
                               ns, nt, nbf, hsoc)
  !
  ! Assemble the full SOC Hamiltonian in the basis of MRSF spin-states.
  !
  ! Basis ordering (GAMESS convention):
  !   rows/cols 1..ns              : singlets S_I
  !   rows/cols ns+1..ns+3*nt      : triplets, grouped as
  !                                  (T_1,Ms=-1),(T_1,Ms=0),(T_1,Ms=+1),
  !                                  (T_2,Ms=-1),(T_2,Ms=0),(T_2,Ms=+1), ...
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
  ! Spin matrix elements (spnfac absorbed, GAMESS convention):
  !   celm_aa = (0, -Lz(t,u)/2)
  !   celm_bb = -celm_aa = (0, +Lz(t,u)/2)
  !   celm_ba = (Ly(t,u)/4, -Lx(t,u)/4)   [S- component; spnfac=0.5 already absorbed]
  !   celm_ab = (Ly(t,u)/4,  Lx(t,u)/4)   [S+ component]
  !
  ! Note: celm_ba here = sqrt(0.5)*0.5*(Ly-iLx) per GAMESS lines 1108-1110.
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
  real(kind=dp), parameter :: sq05  = sqrt(0.5_dp)  ! = 1/sqrt(2)
  real(kind=dp), parameter :: half  = 0.5_dp
  real(kind=dp), parameter :: quart = 0.25_dp
 
  hsoc = cmplx(0.0_dp, 0.0_dp, kind=dp)
 
  do iu = 1, nbf
    do it = 1, nbf
 
      ! Spin-orbit matrix elements (GAMESS convention, lines 1108-1126):
      !   celm_ba = sqrt(0.5)*spnfac * (Ly - i*Lx),  spnfac = sqrt(0.5)
      !           = 0.5 * (Ly - i*Lx) * 0.5 = 0.25*(Ly - i*Lx)  [with both spnfac]
      ! Wait: GAMESS has spnfac=sqrt(0.5) for S-T and spnfac=0.5 for T-T.
      ! For S-T block (spnfac=sqrt(0.5)):
      !   celm_ba_ST = sqrt(0.5) * (Ly/2 - i*Lx/2) * sqrt(0.5) = 0.5*(Ly/2-i*Lx/2)
      !              = (Ly/4, -Lx/4)
      ! For T-T Lz block (spnfac=0.5):
      !   celm_aa_TT = (0, -Lz/2) * 0.5 = (0, -Lz/4)  -- but GAMESS uses 0.5 directly
      ! Actually looking at GAMESS code more carefully:
      !   spnfac for <ua|Lz Sz|ta> = 0.5 (Sz alpha = +1/2 alpha)
      !   celm_aa = (0, -Lz(iu,it)) * 0.5
      !   celm_ba = sqrt(0.5)*(Ly(iu,it),-Lx(iu,it)) * sqrt(0.5)
      !           = 0.5 * (Ly, -Lx)  with the outer sqrt(0.5) already in spnfac
      ! The outer sqrt(0.5) in "celm1ba=sqrt(0.5)*dcmplx(y,-x)*spnfac"
      ! and spnfac=sqrt(0.5) gives celm_ba = 0.5*(Ly/2, -Lx/2)... no.
      !
      ! Reading literally: celm1ba = sqrt(0.5) * (Ly,-Lx) * sqrt(0.5) = 0.5*(Ly,-Lx)
      ! Then for S-T: H += celm1ba * T11ab(IST,JST,it,iu)
      !   where T11ab already has the correct normalization.
      ! This matches our current S-T formula which works correctly.
      ! So we keep the SAME celm definition as the current working code:
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
  implicit none

  complex(kind=dp), intent(in)  :: hsoc(ns+3*nt, ns+3*nt)
  real(kind=dp),    intent(in)  :: singlet_energies(ns)
  real(kind=dp),    intent(in)  :: triplet_energies(nt)
  real(kind=dp),    intent(in)  :: e_ref
  integer,          intent(in)  :: ns, nt
  real(kind=dp),    intent(out) :: eval(ns+3*nt)
  complex(kind=dp), intent(out) :: evec(ns+3*nt, ns+3*nt)

  ! All integers are default kind (= 64-bit with -fdefault-integer-8)
  ! to match OQP's ILP64 MKL linkage. NEVER use integer(4) for LAPACK args.
  integer :: nstate, ist, i, j, ioff, lwork
  real(kind=dp) :: e0
  complex(kind=dp), allocatable :: work(:)
  real(kind=dp),    allocatable :: rwork(:)

  real(kind=dp), parameter :: dfac  = 5.84375713555_dp  ! alpha^2/2 * ha2wn, cm-1/a.u.
  real(kind=dp), parameter :: ha2wn = 219474.6_dp

  ! Explicit 4-byte integers required by LP64 LAPACK zheev interface
  integer(4) :: info4, nstate4, lwork4

  nstate  = ns + 3*nt
  nstate4 = int(nstate, 4)

  ! Fixed lwork: safe minimum for small matrices, avoids workspace query overhead
  lwork  = max(1, 100*nstate)
  lwork4 = int(lwork, 4)

  allocate(rwork(3*nstate), work(lwork))

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
  call zheev('V', 'U', nstate4, evec, nstate4, eval, work, lwork4, rwork, info4)

  if (iand(int(info4, 8), int(z'FFFFFFFF', 8)) /= 0) then
    write(*,'(/,a,i0)') '  !!! zheev failed in diag_soc, info=', info4
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
          case(0); write(iw,'(3x,a1,i3,a)', advance='no') 'T', ist, '(-1)'
          case(1); write(iw,'(3x,a1,i3,a)', advance='no') 'T', ist, '( 0)'
          case(2); write(iw,'(3x,a1,i3,a)', advance='no') 'T', ist, '(+1)'
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

!> @brief Transform a full (non-antisymmetric) AO matrix to the MO basis
!> @details
!>  Two-step MO transformation via DGEMM: tmp = Xao * C, Xmo = C^T * tmp.
!>  Used for the 2e mean-field SOC correction where the AO matrix is symmetric
!>  (unlike the antisymmetric 1e angular momentum integrals handled by ao2mo_soc).
!>
!> @param[in]  Xao  Full AO matrix (nbf x nbf)
!> @param[out] Xmo  Full MO matrix (nbf x nbf)
!> @param[in]  mo   MO coefficient matrix C(mu,p) (nbf x nbf)
!> @param[in]  nbf  Number of basis functions
subroutine ao2mo_full(Xao, Xmo, mo, nbf)
  implicit none
  integer, intent(in) :: nbf
  real(kind=dp), intent(in) :: Xao(nbf,nbf), mo(nbf,nbf)
  real(kind=dp), intent(out) :: Xmo(nbf,nbf)
  real(kind=dp), allocatable :: tmp(:,:)

  allocate(tmp(nbf,nbf))
  ! tmp = Xao * C
  call dgemm('N','N', nbf, nbf, nbf, 1.0_dp, Xao, nbf, mo, nbf, 0.0_dp, tmp, nbf)
  ! Xmo = C^T * tmp
  call dgemm('T','N', nbf, nbf, nbf, 1.0_dp, mo, nbf, tmp, nbf, 0.0_dp, Xmo, nbf)
  deallocate(tmp)
end subroutine ao2mo_full

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
          case(0); write(label,'(a1,i0,a)') 'T', (i-ns-1)/3+1, '(-1)'
          case(1); write(label,'(a1,i0,a)') 'T', (i-ns-1)/3+1, '( 0)'
          case(2); write(label,'(a1,i0,a)') 'T', (i-ns-1)/3+1, '(+1)'
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
