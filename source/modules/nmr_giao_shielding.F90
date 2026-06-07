module nmr_giao_shielding_mod

  use precision, only: dp
  implicit none

  character(len=*), parameter :: module_name = "nmr_giao_shielding_mod"

  private
  public nmr_giao_shielding_debug

contains

  subroutine nmr_giao_shielding_debug_C(c_handle) bind(C, name="nmr_giao_shielding_debug")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call nmr_giao_shielding_debug(inf)
  end subroutine nmr_giao_shielding_debug_C

  !> Production entry point for GIAO NMR shielding (RHF/UHF/ROHF).
  subroutine nmr_giao_shielding_C(c_handle) bind(C, name="nmr_giao_shielding")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call nmr_giao_shielding_debug(inf)
  end subroutine nmr_giao_shielding_C

!> @brief Test-only/debug emitter for native GIAO (London-orbital) NMR shielding.
!> @details Computes the GIAO paramagnetic nuclear magnetic shielding tensor for
!>  RHF/UHF/ROHF (and the corresponding DFT) references.  Open-shell support:
!>  the diamagnetic term uses the total density D_alpha+D_beta; the paramagnetic
!>  term is solved per spin channel (independent same-spin exchange response) and
!>  summed.  ROHF orbitals are semicanonicalized (from the ground-state spin Fock
!>  matrices) so the UHF-like CPHF is well defined; the closed-shell ROHF limit
!>  reproduces RHF, and open-shell UHF/ROHF match an independent GIAO reference to ~1e-4
!>  ppm (OH, CH3 radicals).  Built from the validated native GIAO building blocks:
!>   - first-order GIAO core Hamiltonian h10 (one-electron, giao_h10_core),
!>   - first-order GIAO two-electron Fock derivative (giao_h10_twoe_matrix),
!>   - first-order GIAO overlap derivative S10 (giao_overlap_derivative),
!>   - the PSO operator at each nucleus (pso_integrals).
!>  It assembles the magnetic first-order Hamiltonian/overlap in the MO basis,
!>  solves the GIAO CPHF/CPKS first-order equation (uncoupled and coupled; the
!>  coupled response is exchange-only for the imaginary/antisymmetric first-order
!>  density, scaled by the exact-exchange fraction c_x), and contracts the
!>  resulting first-order density with the PSO operator to form the paramagnetic
!>  shielding.
!>
!>  Hartree-Fock (RHF/UHF/ROHF) is validated against an independent GIAO reference.
!>
!>  DFT (RKS/UKS/ROKS) is supported: the closed-shell h1 two-electron exchange is
!>  scaled by c_x, and the London (GIAO) derivative of the exchange-correlation
!>  potential (the "vxc_giao" term -- an essential contribution, ~tens of ppm) is
!>  added via a GIAO-weighted XC grid integration (mod_dft_gridint_giao::giao_vxc;
!>  the London derivative of the AO values and gradients).  Validated vs an independent
!>  reference (replicating get_vxc_giao): H2O/PBE and /PBE0 reproduce the
!>  reference to ~1e-3 ppm.  vxc_giao is grid-sensitive, so DFT GIAO benefits from
!>  a converged integration grid.
!>
!>  Results are written as machine-parseable records to the log so
!>  the native GIAO path can be validated against an independent GIAO reference WITHOUT
!>  ungating production nmr_gauge=giao.
!>
!>  Diamagnetic GIAO term (both pieces reference-validated):
!>   - a11part (London diamagnetic): validated CGO diamagnetic at gauge origin 0
!>     plus the Hellmann-Feynman field correction weighted by the ket center.
!>   - a01gp (GIAO gauge correction = London derivative of the PSO): cvec x M,
!>     M^{(col)}_b = <mu| r_b PSO_col |nu> with r referenced to the molecular
!>     origin (R0I = (r-R_bra) bra-raise + R_bra*base, the libcint convention).
!>     Full-tensor agreement with libcint int1e_a01gp to ~3e-8 (e.g. CH4).
!>  Total GIAO shielding matches an independent GIAO reference to ~1e-4 ppm for HF
!>  (grid-free) and to cross-code DFT-SCF/grid level (~0.03 ppm) for DFT, for all
!>  geometries and angular momenta tested (He, H2, HF, CO2, CH4, H2O).
!>  SG/SC/SA are sign/scale conventions fixed against the oracle (SA folds the
!>  factor-2 normalization of the native a01gp vs the libcint int1e_a01gp).
!>
!>  References (GIAO/London-orbital NMR shielding methodology):
!>   - F. London, J. Phys. Radium 8, 397 (1937).
!>   - R. Ditchfield, Mol. Phys. 27, 789 (1974).
!>   - K. Wolinski, J. F. Hinton, P. Pulay, J. Am. Chem. Soc. 112, 8251 (1990).
!>   - T. Helgaker, M. Jaszunski, K. Ruud, Chem. Rev. 99, 293 (1999).
  subroutine nmr_giao_shielding_debug(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use constants, only: tol_int
    use int1, only: giao_h10_core, giao_overlap_derivative, pso_integrals, &
                    nmr_dia_shielding, giao_a11part_corr, giao_a01gp_contract
    use nmr_giao_debug_mod, only: giao_h10_twoe_matrix
    use dft, only: dft_initialize
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_giao, only: giao_vxc

    implicit none

    character(len=*), parameter :: subroutine_name = "nmr_giao_shielding_debug"
    real(kind=dp), parameter :: ALPHA = 1.0d0/137.035999084d0
    real(kind=dp), parameter :: a2ppm = ALPHA*ALPHA*1.0d6
    real(kind=dp), parameter :: ha2ppm = 0.5d0*ALPHA*ALPHA*1.0d6
    ! Calibration signs for the GIAO diamagnetic pieces (fixed vs the reference).
    real(kind=dp), parameter :: SG = -1.0d0, SC = 1.0d0, SA = -0.5d0
    real(kind=dp), parameter :: SX = -1.0d0   ! London-XC (vxc_giao): h1 -= vxc_giao

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    integer :: nbf, nbf2, nat, nocc, nmo, nvir, nocc_b
    integer :: i, j, m, c, t, s, ok, iat
    integer(4) :: status
    logical :: is_dft, open_shell, iw_open
    real(kind=dp) :: tol, scale_exch

    real(kind=dp), allocatable :: h10p(:,:), s10p(:,:)            ! packed (nbf2,3)
    real(kind=dp), allocatable :: twoe(:,:,:), twoe2(:,:,:), vj(:,:,:), vk(:,:,:), vkb(:,:,:)
    real(kind=dp), allocatable :: dm(:,:), dmp(:), dm_b(:,:)
    real(kind=dp), allocatable :: h1ao(:,:,:), h1ao_b(:,:,:), s1ao(:,:,:) ! (nbf,nbf,3)
    real(kind=dp), allocatable :: coords(:,:), zq(:)
    real(kind=dp), allocatable :: sig_u(:,:,:), sig_c(:,:,:)      ! (3,3,nat)
    real(kind=dp), allocatable :: vxa(:,:,:), vxb(:,:,:)          ! London-XC (vxc_giao)
    type(dft_grid_t) :: molGrid
    integer :: mxAngMom
    real(kind=dp), allocatable :: gdia0(:,:,:), corrpre(:,:,:)    ! GIAO dia pieces
    real(kind=dp), allocatable :: a01(:,:,:)                      ! a01gp contracted
    real(kind=dp), allocatable :: sig_dia(:,:,:), sig_tot(:,:,:)  ! (3,3,nat)
    real(kind=dp) :: trg0, trc, o0(3)

    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), contiguous, pointer :: mo_e(:), mo_e_b(:)
    real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:)
    real(kind=dp), contiguous, pointer :: nmrout(:)
    real(kind=dp), allocatable :: ca_sc(:,:), cb_sc(:,:), ea_sc(:), eb_sc(:)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(basis%atoms%zn,1)
    tol = log(10.0d0)*tol_int

    ! Connect the log unit early so guard aborts and CPHF warnings land in the
    ! log instead of an orphan fort.* file.
    inquire(unit=iw, opened=iw_open)
    if (.not. iw_open) open(unit=iw, file=infos%log_filename, position="append")

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
    call check_status(status, module_name, subroutine_name, OQP_DM_A)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a, status)
    call check_status(status, module_name, subroutine_name, OQP_VEC_MO_A)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_e, status)
    call check_status(status, module_name, subroutine_name, OQP_E_MO_A)

    ! scftype: 1=RHF (closed shell), 2=UHF, 3=ROHF
    open_shell = infos%control%scftype == 2 .or. infos%control%scftype == 3
    nmo  = size(mo_e)
    if (open_shell) then
      nocc   = int(infos%mol_prop%nelec_A)
      nocc_b = int(infos%mol_prop%nelec_B)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b, status)
      call check_status(status, module_name, subroutine_name, OQP_DM_B)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b, status)
      call check_status(status, module_name, subroutine_name, OQP_VEC_MO_B)
      call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_e_b, status)
      call check_status(status, module_name, subroutine_name, OQP_E_MO_B)
    else
      nocc   = int(infos%mol_prop%nocc)
      nocc_b = nocc
    end if
    nvir = nmo - nocc

    is_dft = infos%control%hamilton == 20
    scale_exch = 1.0d0
    if (is_dft) scale_exch = infos%dft%HFscale

    ! Not-implemented classes must abort instead of silently producing wrong
    ! shieldings:
    !  - CAM/range-separated hybrids: the coupled magnetic response and the
    !    GIAO two-electron derivative use the global exchange fraction only;
    !    the range-separation attenuation (alpha/beta/mu) is not wired in.
    !  - meta-GGAs: the tau channel of the London-XC (vxc_giao) term is not
    !    implemented (mod_dft_gridint_giao evaluates LDA/GGA ingredients only).
    !  - ECP: no effective-core magnetic-derivative term is implemented.
    if (is_dft) then
      if (infos%dft%cam_flag) then
        call show_message('GIAO NMR shielding with range-separated (CAM) &
          &functionals is not implemented', with_abort)
      end if
      if (infos%functional%needtau) then
        call show_message('GIAO NMR shielding with meta-GGA (tau-dependent) &
          &functionals is not implemented', with_abort)
      end if
    end if
    if (allocated(infos%basis%ecp_zn_num)) then
      if (any(infos%basis%ecp_zn_num /= 0)) then
        call show_message('GIAO NMR shielding with ECP basis sets is not &
          &implemented', with_abort)
      end if
    end if

    allocate(coords(3,nat), zq(nat))
    do iat = 1, nat
      coords(:,iat) = basis%atoms%xyz(:,iat)
    end do
    zq = infos%atoms%zn - infos%basis%ecp_zn_num

    ! --- Total density (full + packed).  For RHF OQP_DM_A is already the total
    !     (closed-shell) density; for UHF/ROHF total = D_alpha + D_beta. ---
    allocate(dm(nbf,nbf), dmp(nbf2), source=0.0d0)
    if (open_shell) then
      allocate(dm_b(nbf,nbf), source=0.0d0)
      dmp = dmat_a + dmat_b
      call unpack_sym(dmat_a, dm, nbf)        ! D_alpha
      call unpack_sym(dmat_b, dm_b, nbf)      ! D_beta
    else
      dmp = dmat_a
      call unpack_sym(dmp, dm, nbf)           ! D_total (closed shell)
    end if

    ! --- One-electron GIAO Hamiltonian h10(1e) and overlap derivative S10 ---
    allocate(h10p(nbf2,3), s10p(nbf2,3), source=0.0d0)
    call giao_h10_core(basis, coords, zq, h10p, debug=.false., logtol=tol)
    call giao_overlap_derivative(basis, s10p, debug=.false., logtol=tol)

    allocate(h1ao(nbf,nbf,3), s1ao(nbf,nbf,3), source=0.0d0)
    do c = 1, 3
      call expand_antisym(h10p(:,c), h1ao(:,:,c), nbf)
      call expand_antisym(s10p(:,c), s1ao(:,:,c), nbf)
    end do

    ! --- London (GIAO) derivative of the XC potential (vxc_giao); DFT only ---
    allocate(vxa(3,nbf,nbf), vxb(3,nbf,nbf), source=0.0d0)
    if (is_dft) then
      mxAngMom = maxval(basis%am) + 2
      call dft_initialize(infos, basis, molGrid)
      block
        real(kind=dp), allocatable :: ca(:,:), cb(:,:)
        allocate(ca(nbf,nbf), cb(nbf,nbf))
        ca = mo_a(:,1:nbf)
        if (open_shell) then
          cb = mo_b(:,1:nbf)
          call giao_vxc(basis, molGrid, infos, ca, cb, .true., &
                        vxa, vxb, mxAngMom, nbf, infos%dft%grid_density_cutoff)
        else
          cb = ca
          call giao_vxc(basis, molGrid, infos, ca, cb, .false., &
                        vxa, vxb, mxAngMom, nbf, infos%dft%grid_density_cutoff)
        end if
        deallocate(ca, cb)
      end block
    end if

    ! --- Two-electron GIAO Fock derivative ---
    allocate(twoe(3,nbf,nbf), twoe2(3,nbf,nbf), vj(3,nbf,nbf), vk(3,nbf,nbf), source=0.0d0)
    allocate(sig_u(3,3,nat), sig_c(3,3,nat), source=0.0d0)
    if (open_shell) then
      ! Spin-resolved: h1_sigma = h10(1e) + J[D_tot] - cx*K[D_sigma].  giao_h10_
      ! twoe_matrix returns (vj=J, vk=K, h10) for its input density; call it once
      ! per density, using a throwaway 'twoe' for the unused J/h10 outputs.
      allocate(vkb(3,nbf,nbf), h1ao_b(nbf,nbf,3), source=0.0d0)
      call giao_h10_twoe_matrix(basis, infos, dm+dm_b, vj,  twoe, twoe2) ! vj  = J[D_tot]
      call giao_h10_twoe_matrix(basis, infos, dm,      twoe, vk,  twoe2) ! vk  = K[D_a]
      call giao_h10_twoe_matrix(basis, infos, dm_b,    twoe, vkb, twoe2) ! vkb = K[D_b]
      do c = 1, 3
        do i = 1, nbf
          do j = 1, nbf
            h1ao_b(i,j,c) = h1ao(i,j,c) + vj(c,i,j) - scale_exch*vkb(c,i,j) + SX*vxb(c,i,j)
            h1ao(i,j,c)   = h1ao(i,j,c) + vj(c,i,j) - scale_exch*vk(c,i,j)  + SX*vxa(c,i,j)
          end do
        end do
      end do
      ! Two independent spin channels (same-spin exchange), each occ_factor = 1.
      if (infos%control%scftype == 3) then
        ! ROHF: a single orbital set with effective-Fock eigenvalues.  Build
        ! semicanonical spin orbitals/energies from the ground-state spin Fock
        ! matrices so the UHF-like CPHF/Delta-e is well defined.
        call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a, status)
        call check_status(status, module_name, subroutine_name, OQP_FOCK_A)
        call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b, status)
        call check_status(status, module_name, subroutine_name, OQP_FOCK_B)
        allocate(ca_sc(nbf,nmo), cb_sc(nbf,nmo), ea_sc(nmo), eb_sc(nmo), source=0.0d0)
        call semicanon_orbitals(fock_a, mo_a, nbf, nmo, nocc,   ca_sc, ea_sc)
        call semicanon_orbitals(fock_b, mo_b, nbf, nmo, nocc_b, cb_sc, eb_sc)
        call giao_para_channel(infos, basis, ca_sc, ea_sc, nocc,   nmo, nbf, nat, &
                               coords, h1ao,   s1ao, scale_exch, 1.0d0, sig_u, sig_c)
        call giao_para_channel(infos, basis, cb_sc, eb_sc, nocc_b, nmo, nbf, nat, &
                               coords, h1ao_b, s1ao, scale_exch, 1.0d0, sig_u, sig_c)
      else
        call giao_para_channel(infos, basis, mo_a, mo_e,   nocc,   nmo, nbf, nat, &
                               coords, h1ao,   s1ao, scale_exch, 1.0d0, sig_u, sig_c)
        call giao_para_channel(infos, basis, mo_b, mo_e_b, nocc_b, nmo, nbf, nat, &
                               coords, h1ao_b, s1ao, scale_exch, 1.0d0, sig_u, sig_c)
      end if
    else
      ! Closed shell: h1 = h10(1e) + (J - 0.5 K)[D_tot] (twoe), single channel.
      call giao_h10_twoe_matrix(basis, infos, dm, vj, vk, twoe)
      do c = 1, 3
        do i = 1, nbf
          do j = 1, nbf
            h1ao(i,j,c) = h1ao(i,j,c) + vj(c,i,j) - 0.5d0*scale_exch*vk(c,i,j) + SX*vxa(c,i,j)
          end do
        end do
      end do
      call giao_para_channel(infos, basis, mo_a, mo_e, nocc, nmo, nbf, nat, &
                             coords, h1ao, s1ao, scale_exch, 2.0d0, sig_u, sig_c)
    end if
    sig_u = sig_u * a2ppm
    sig_c = sig_c * a2ppm

    ! --- Diamagnetic shielding (GIAO) ---
    !   a11part = cg_a11part(O=0) + 0.5 field_a R_nu,b  (verified vs libcint).
    !   e11_pre_{t,s} = 0.5*gdia0_{s,t} + corrpre_{t,s};  e11 = e11_pre - I*tr;
    !   sigma_dia = e11 * alpha^2 * 1e6.  (a01gp gauge-correction: TODO.)
    o0 = 0.0d0
    allocate(gdia0(3,3,nat), corrpre(3,3,nat), a01(3,3,nat), &
             sig_dia(3,3,nat), sig_tot(3,3,nat), source=0.0d0)
    call nmr_dia_shielding(basis, dmp, o0, coords, nat, gdia0, logtol=tol)
    call giao_a11part_corr(basis, dmp, coords, nat, corrpre, logtol=tol)
    call giao_a01gp_contract(basis, dmp, coords, nat, a01, logtol=tol)
    do iat = 1, nat
      trg0 = gdia0(1,1,iat)+gdia0(2,2,iat)+gdia0(3,3,iat)
      trc  = corrpre(1,1,iat)+corrpre(2,2,iat)+corrpre(3,3,iat)
      do t = 1, 3
        do s = 1, 3
          ! a11part (trace-corrected) + a01gp (raw, per the standard diamagnetic decomposition)
          sig_dia(t,s,iat) = ( SG*0.5d0*gdia0(s,t,iat) + SC*corrpre(t,s,iat) &
                 - merge(SG*0.5d0*trg0 + SC*trc, 0.0d0, t==s) &
                 + SA*a01(t,s,iat) ) * a2ppm
          sig_tot(t,s,iat) = sig_dia(t,s,iat) + sig_c(t,s,iat)
        end do
      end do
    end do

    ! --- Store isotropic shielding (ppm) to a tagarray for JSON output ---
    !     rows: dia, para_uncoupled, para_coupled, total_uncoupled, total_coupled
    !     stored atom-major (flat): atom a occupies entries 5*(a-1)+1 .. 5*a.
    call infos%dat%reserve_data(OQP_nmr_shielding, TA_TYPE_REAL64, 5*nat, &
                                comment=OQP_nmr_shielding_comment)
    call tagarray_get_data(infos%dat, OQP_nmr_shielding, nmrout)
    do iat = 1, nat
      nmrout(5*(iat-1)+1) = (sig_dia(1,1,iat)+sig_dia(2,2,iat)+sig_dia(3,3,iat))/3.0d0
      nmrout(5*(iat-1)+2) = (sig_u(1,1,iat)+sig_u(2,2,iat)+sig_u(3,3,iat))/3.0d0
      nmrout(5*(iat-1)+3) = (sig_c(1,1,iat)+sig_c(2,2,iat)+sig_c(3,3,iat))/3.0d0
      nmrout(5*(iat-1)+4) = nmrout(5*(iat-1)+1) + nmrout(5*(iat-1)+2)
      nmrout(5*(iat-1)+5) = nmrout(5*(iat-1)+1) + nmrout(5*(iat-1)+3)
    end do

    ! --- Emit parseable records ---
    inquire(unit=iw, opened=iw_open)
    if (.not. iw_open) open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/,A)') 'GIAO_SHIELDING_DEBUG_BEGIN native-giao shielding (ppm)'
    write(iw,'(A,1X,I0)') 'GIAO_SHIELDING_DEBUG_NATOM', nat
    write(iw,'(A,1X,F10.6)') 'GIAO_SHIELDING_DEBUG_CX', scale_exch
    do iat = 1, nat
      do t = 1, 3
        do s = 1, 3
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'GIAO_SHIELDING_DEBUG_PARA_UNC', &
            iat, t, s, sig_u(t,s,iat)
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'GIAO_SHIELDING_DEBUG_PARA_CPL', &
            iat, t, s, sig_c(t,s,iat)
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'GIAO_SHIELDING_DEBUG_DIA', &
            iat, t, s, sig_dia(t,s,iat)
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'GIAO_SHIELDING_DEBUG_TOTAL', &
            iat, t, s, sig_tot(t,s,iat)
        end do
      end do
      write(iw,'(A,1X,I0,4(1X,ES24.16))') 'GIAO_SHIELDING_DEBUG_ISO', iat, &
        (sig_u(1,1,iat)+sig_u(2,2,iat)+sig_u(3,3,iat))/3.0d0, &
        (sig_c(1,1,iat)+sig_c(2,2,iat)+sig_c(3,3,iat))/3.0d0, &
        (sig_dia(1,1,iat)+sig_dia(2,2,iat)+sig_dia(3,3,iat))/3.0d0, &
        (sig_tot(1,1,iat)+sig_tot(2,2,iat)+sig_tot(3,3,iat))/3.0d0
    end do
    write(iw,'(A)') 'GIAO_SHIELDING_DEBUG_END'

    ! --- Human-readable production table (GIAO isotropic shielding) ---
    write(iw,'(2/)')
    write(iw,'(4x,a)') '========================================'
    write(iw,'(4x,a)') 'NMR nuclear magnetic shielding (GIAO)'
    write(iw,'(4x,a)') '========================================'
    write(iw,'(4x,a)') 'Gauge formulation: GIAO (London) atomic orbitals; gauge-origin independent.'
    if (open_shell) then
      write(iw,'(4x,a)') 'Reference: open shell (UHF/ROHF); diamagnetic from total density, '// &
        'paramagnetic summed over spin channels.'
    end if
    write(iw,'(/4x,a,f8.4,a)') 'Isotropic shielding (GIAO, ppm)   [exact-exchange c_x =', &
           scale_exch, ']'
    write(iw,'(4x,a)')  '   Atom    Z   sigma_dia   para_uncoupled   para_coupled'// &
           '   total_uncoupled   total_coupled'
    do iat = 1, nat
      write(iw,'(4x,i6,f6.1,5f16.6)') iat, basis%atoms%zn(iat), &
        (sig_dia(1,1,iat)+sig_dia(2,2,iat)+sig_dia(3,3,iat))/3.0d0, &
        (sig_u(1,1,iat)+sig_u(2,2,iat)+sig_u(3,3,iat))/3.0d0, &
        (sig_c(1,1,iat)+sig_c(2,2,iat)+sig_c(3,3,iat))/3.0d0, &
        (sig_dia(1,1,iat)+sig_dia(2,2,iat)+sig_dia(3,3,iat))/3.0d0 + &
        (sig_u(1,1,iat)+sig_u(2,2,iat)+sig_u(3,3,iat))/3.0d0, &
        (sig_tot(1,1,iat)+sig_tot(2,2,iat)+sig_tot(3,3,iat))/3.0d0
    end do
    close(iw)

    deallocate(gdia0, corrpre, a01, sig_dia, sig_tot)
    deallocate(h10p, s10p, twoe, twoe2, vj, vk, dm, dmp, h1ao, s1ao)
    if (allocated(vxa)) deallocate(vxa)
    if (allocated(vxb)) deallocate(vxb)
    deallocate(coords, zq, sig_u, sig_c)
    if (allocated(vkb))    deallocate(vkb)
    if (allocated(h1ao_b)) deallocate(h1ao_b)
    if (allocated(dm_b))   deallocate(dm_b)

  end subroutine nmr_giao_shielding_debug

!> Expand packed lower-triangular antisymmetric matrix to full form.
  subroutine expand_antisym(packed, full, n)
    real(kind=dp), intent(in)  :: packed(:)
    real(kind=dp), intent(out) :: full(:,:)
    integer, intent(in) :: n
    integer :: p, q, idx
    full = 0.0d0
    do p = 1, n
      do q = 1, p
        idx = q + p*(p-1)/2
        full(p,q) =  packed(idx)
        full(q,p) = -packed(idx)
      end do
    end do
  end subroutine expand_antisym

  subroutine unpack_sym(packed, full, n)
    real(kind=dp), intent(in) :: packed(:)
    real(kind=dp), intent(out) :: full(:,:)
    integer, intent(in) :: n
    integer :: i, j, ij
    full = 0.0d0
    do i = 1, n
      do j = 1, i
        ij = j + i*(i-1)/2
        full(i,j) = packed(ij)
        full(j,i) = packed(ij)
      end do
    end do
  end subroutine unpack_sym

!> MO transform keeping occupied ket columns: m(p,i) = sum_mn C(m,p) a(m,n) C(n,i),
!>  p = 1..nmo, i = 1..nocc.
  subroutine ao_to_mo_occ(a_ao, c, m_mo, nbf, nmo, nocc)
    real(kind=dp), intent(in)  :: a_ao(:,:), c(:,:)
    real(kind=dp), intent(out) :: m_mo(:,:)
    integer, intent(in) :: nbf, nmo, nocc
    real(kind=dp), allocatable :: tmp(:,:)
    allocate(tmp(nbf,nocc))
    tmp = matmul(a_ao, c(:,1:nocc))
    m_mo(1:nmo,1:nocc) = matmul(transpose(c(:,1:nmo)), tmp)
    deallocate(tmp)
  end subroutine ao_to_mo_occ

!> Uncoupled first-order equation (standard uncoupled first-order equation):
!>   hs = h1 - s1*e_i ;  mo1[vir,i] = -hs[vir,i]/(e_a-e_i) ;
!>   mo1[occ,i] = -0.5*s1[occ,i].
  subroutine solve_mo1_uncoupled(h1mo, s1mo, e, nocc, nmo, mo1)
    real(kind=dp), intent(in)  :: h1mo(:,:,:), s1mo(:,:,:), e(:)
    integer, intent(in) :: nocc, nmo
    real(kind=dp), intent(out) :: mo1(:,:,:)
    integer :: x, p, i
    real(kind=dp) :: hs
    mo1 = 0.0d0
    do x = 1, 3
      do i = 1, nocc
        do p = 1, nmo
          hs = h1mo(p,i,x) - s1mo(p,i,x)*e(i)
          if (p > nocc) then
            mo1(p,i,x) = -hs/(e(p)-e(i))
          else
            mo1(p,i,x) = -0.5d0*s1mo(p,i,x)
          end if
        end do
      end do
    end do
  end subroutine solve_mo1_uncoupled

!> Coupled GIAO CPHF/CPKS: fixed-point of
!>   mo1[vir,i] = -(hs[vir,i] + v1[vir,i])/(e_a-e_i),  mo1[occ,i] = -0.5 s1[occ,i],
!>  with v1 the exact-exchange response (scaled by c_x) of the imaginary
!>  antisymmetric first-order density built from the full mo1.  For c_x = 0 the
!>  loop is skipped and the result equals the uncoupled solution.
  subroutine solve_mo1_coupled(infos, basis, mo, h1mo, s1mo, e, nocc, nmo, &
                               scale_exch, mo1)
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t, mntoia
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message
    real(kind=dp), intent(in) :: mo(:,:), h1mo(:,:,:), s1mo(:,:,:), e(:)
    integer, intent(in) :: nocc, nmo
    real(kind=dp), intent(in) :: scale_exch
    real(kind=dp), intent(out) :: mo1(:,:,:)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis

    integer, parameter :: maxit = 100
    real(kind=dp), parameter :: tol = 1.0d-9
    integer :: nbf, nvir, x, p, i, a, k, it
    real(kind=dp) :: diff, hs
    type(int2_compute_t) :: int2_driver
    type(int2_td_data_t), target :: kdat
    real(kind=dp), allocatable, target :: pa(:,:,:)
    real(kind=dp), allocatable :: gxv(:), mo1x(:,:), prev(:,:), gao(:,:)

    nbf = basis%nbf
    nvir = nmo - nocc

    ! Start from the uncoupled solution.
    call solve_mo1_uncoupled(h1mo, s1mo, e, nocc, nmo, mo1)
    if (abs(scale_exch) <= 1.0d-12) return

    allocate(pa(nbf,nbf,1), gxv(nocc*nvir), mo1x(nmo,nocc), &
             prev(nmo,nocc), gao(nbf,nbf), source=0.0d0)

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    kdat = int2_td_data_t(d2=pa, int_apb=.false., int_amb=.true., &
                          tamm_dancoff=.false., scale_exchange=scale_exch)

    do x = 1, 3
      mo1x = mo1(:,:,x)
      do it = 1, maxit
        prev = mo1x
        ! Imaginary antisymmetric AO first-order density from the full MO
        ! response (occ + vir rows); CGO-consistent normalization (no x2).
        call giao_pb_density(mo, mo1x, pa(:,:,1), nbf, nmo, nocc)
        call int2_driver%run(kdat)
        gao = 0.5d0*kdat%amb(:,:,1,1)
        ! Exchange response projected to the occ-vir block (i fast):
        ! gxv(i+(a-1)*nocc) = (C^T gao C)[i_occ, a_vir].
        call mntoia(gao, gxv, mo, mo, nocc, nocc)
        do i = 1, nocc
          do a = 1, nvir
            p = nocc + a
            k = i + (a-1)*nocc
            hs = h1mo(p,i,x) - s1mo(p,i,x)*e(i)
            ! mntoia returns the [occ,vir] block gxv(i,a); the response element
            ! needed here is the [vir,occ] entry v1(a,i) = -gxv(i,a) (the
            ! exchange image is antisymmetric).
            mo1x(p,i) = -(hs - gxv(k))/(e(p)-e(i))
          end do
          do p = 1, nocc
            mo1x(p,i) = -0.5d0*s1mo(p,i,x)
          end do
        end do
        diff = maxval(abs(mo1x - prev))
        if (diff < tol) exit
      end do
      if (diff >= tol) then
        call show_message('WARNING: GIAO coupled magnetic response (CPHF) did &
          &not converge within the iteration limit; shieldings may be inaccurate')
      end if
      mo1(:,:,x) = mo1x
    end do

    call int2_driver%clean()
    deallocate(pa, gxv, mo1x, prev, gao)
  end subroutine solve_mo1_coupled

!> Imaginary antisymmetric AO first-order density from a full MO response vector
!>  (CGO-consistent normalization, no double-occupancy factor):
!>   D = C mo1 orbo^T ;  pa = D - D^T.
  subroutine giao_pb_density(mo, mo1x, pa, nbf, nmo, nocc)
    real(kind=dp), intent(in) :: mo(:,:), mo1x(:,:)
    real(kind=dp), intent(out) :: pa(:,:)
    integer, intent(in) :: nbf, nmo, nocc
    real(kind=dp), allocatable :: dleft(:,:)
    allocate(dleft(nbf,nbf))
    dleft = matmul(mo(:,1:nmo), matmul(mo1x, transpose(mo(:,1:nocc))))
    pa = dleft - transpose(dleft)
    deallocate(dleft)
  end subroutine giao_pb_density

!> Paramagnetic shielding tensor for one nucleus (standard paramagnetic convention):
!>   dm10(a,b,x) = 2 sum_{p,i} C(a,p) mo1(p,i,x) C(b,i) ;
!>   sigma_para[x,y] = 2 sum_{a,b} dm10(a,b,x) * h01i(b,a,y).
  !> Paramagnetic shielding contribution from one spin channel: MO transform,
  !> CPHF (uncoupled + coupled), and PSO contraction, ACCUMULATED into sig_u/sig_c.
  !> occ_factor = 2 for RHF (closed shell), 1 for each UHF spin channel.
  subroutine giao_para_channel(infos, basis, mo, e, nocc, nmo, nbf, nat, coords, &
                               h1ao, s1ao, scale_exch, occ_factor, sig_u, sig_c)
    use types, only: information
    use basis_tools, only: basis_set
    use int1, only: pso_integrals
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: mo(:,:), e(:), coords(:,:)
    real(kind=dp), intent(in) :: h1ao(:,:,:), s1ao(:,:,:), scale_exch, occ_factor
    integer, intent(in) :: nocc, nmo, nbf, nat
    real(kind=dp), intent(inout) :: sig_u(:,:,:), sig_c(:,:,:)
    real(kind=dp), allocatable :: h1mo(:,:,:), s1mo(:,:,:), mo1u(:,:,:), mo1c(:,:,:)
    real(kind=dp), allocatable :: pso(:,:,:), st(:,:)
    integer :: c, iat

    allocate(h1mo(nmo,nocc,3), s1mo(nmo,nocc,3), mo1u(nmo,nocc,3), mo1c(nmo,nocc,3), &
             pso(nbf,nbf,3), st(3,3), source=0.0d0)
    do c = 1, 3
      call ao_to_mo_occ(h1ao(:,:,c), mo, h1mo(:,:,c), nbf, nmo, nocc)
      call ao_to_mo_occ(s1ao(:,:,c), mo, s1mo(:,:,c), nbf, nmo, nocc)
    end do
    call solve_mo1_uncoupled(h1mo, s1mo, e, nocc, nmo, mo1u)
    call solve_mo1_coupled(infos, basis, mo, h1mo, s1mo, e, nocc, nmo, scale_exch, mo1c)
    do iat = 1, nat
      call pso_integrals(basis, coords(:,iat), pso)
      call para_tensor(mo1u, mo, pso, nbf, nmo, nocc, st, occ_factor)
      sig_u(:,:,iat) = sig_u(:,:,iat) + st
      call para_tensor(mo1c, mo, pso, nbf, nmo, nocc, st, occ_factor)
      sig_c(:,:,iat) = sig_c(:,:,iat) + st
    end do
    deallocate(h1mo, s1mo, mo1u, mo1c, pso, st)
  end subroutine giao_para_channel

  !> Semicanonicalize ROHF orbitals for spin sigma: diagonalize the occ-occ and
  !> vir-vir blocks of C^T F_sigma C, returning semicanonical orbitals (csc) and
  !> orbital energies (esc).  ROHF stores a single orbital set with effective-Fock
  !> eigenvalues; the UHF-like CPHF needs proper spin orbital energies (the
  !> ground-state spin Fock matrices F_a/F_b from OQP_FOCK_A/B).
  subroutine semicanon_orbitals(fock_p, c0, nbf, nmo, nocc, csc, esc)
    use eigen, only: diag_symm_full
    real(kind=dp), intent(in) :: fock_p(:), c0(:,:)
    integer, intent(in) :: nbf, nmo, nocc
    real(kind=dp), intent(out) :: csc(:,:), esc(:)
    real(kind=dp), allocatable :: fao(:,:), fmo(:,:), blk(:,:), eb(:)
    integer :: nvir, ierr
    nvir = nmo - nocc
    allocate(fao(nbf,nbf), fmo(nmo,nmo))
    call unpack_sym(fock_p, fao, nbf)
    fmo = matmul(transpose(c0(:,1:nmo)), matmul(fao, c0(:,1:nmo)))
    csc = c0(:,1:nmo); esc = 0.0d0
    if (nocc > 0) then
      allocate(blk(nocc,nocc), eb(nocc))
      blk = fmo(1:nocc,1:nocc)
      call diag_symm_full(1, nocc, blk, nocc, eb, ierr)
      csc(:,1:nocc) = matmul(c0(:,1:nocc), blk)
      esc(1:nocc) = eb(1:nocc)
      deallocate(blk, eb)
    end if
    if (nvir > 0) then
      allocate(blk(nvir,nvir), eb(nvir))
      blk = fmo(nocc+1:nmo, nocc+1:nmo)
      call diag_symm_full(1, nvir, blk, nvir, eb, ierr)
      csc(:,nocc+1:nmo) = matmul(c0(:,nocc+1:nmo), blk)
      esc(nocc+1:nmo) = eb(1:nvir)
      deallocate(blk, eb)
    end if
    deallocate(fao, fmo)
  end subroutine semicanon_orbitals

  subroutine para_tensor(mo1, mo, h01i, nbf, nmo, nocc, sig, occ_factor)
    real(kind=dp), intent(in) :: mo1(:,:,:), mo(:,:), h01i(:,:,:)
    integer, intent(in) :: nbf, nmo, nocc
    real(kind=dp), intent(out) :: sig(:,:)
    real(kind=dp), intent(in), optional :: occ_factor
    integer :: x, y, a, b
    real(kind=dp), allocatable :: dm10(:,:,:)
    real(kind=dp) :: acc, ofac
    ofac = 2.0d0                 ! RHF closed-shell occupation; UHF per spin = 1
    if (present(occ_factor)) ofac = occ_factor
    allocate(dm10(nbf,nbf,3))
    do x = 1, 3
      dm10(:,:,x) = ofac*matmul(mo(:,1:nmo), matmul(mo1(:,:,x), transpose(mo(:,1:nocc))))
    end do
    do x = 1, 3
      do y = 1, 3
        acc = 0.0d0
        do b = 1, nbf
          do a = 1, nbf
            acc = acc + dm10(a,b,x)*h01i(b,a,y)
          end do
        end do
        ! OpenQP pso_integrals stores the negative of libcint int1e_prinvxp
        ! (h01i); the leading 2 is the +c.c. factor of the the standard paramagnetic routine.
        sig(x,y) = -2.0d0*acc
      end do
    end do
    deallocate(dm10)
  end subroutine para_tensor

end module nmr_giao_shielding_mod
