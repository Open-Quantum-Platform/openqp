module tdhf_mrsf_ekt_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_ekt_mod"

contains

  subroutine tdhf_mrsf_ekt_ip_C(c_handle) bind(C, name="tdhf_mrsf_ekt_ip")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_ekt(inf, .false.)
  end subroutine tdhf_mrsf_ekt_ip_C

  subroutine tdhf_mrsf_ekt_ea_C(c_handle) bind(C, name="tdhf_mrsf_ekt_ea")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_ekt(inf, .true.)
  end subroutine tdhf_mrsf_ekt_ea_C

  subroutine tdhf_mrsf_ekt(infos, electron_affinity)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information, energy_results
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time
    use precision, only: dp
    use mathlib, only: unpack_matrix, pack_matrix, orthogonal_transform_sym
    use eigen, only: diag_symm_full
    use printing, only: print_module_info
    use tdhf_mrsf_energy_mod, only: tdhf_mrsf_energy
    use tdhf_mrsf_z_vector_mod, only: tdhf_mrsf_z_vector
    use dft, only: dft_initialize
    use mod_dft_molgrid, only: dft_grid_t
    use scf_addons, only: calc_fock, scf_energy_t

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_ekt"
    type(information), target, intent(inout) :: infos
    logical, intent(in) :: electron_affinity

    type(basis_set), pointer :: basis
    integer :: nbf, nbf2, nroot, ok, i, nfocks, nschwz
    integer :: nactive
    integer, allocatable :: dom_mo_idx(:), dom_no_idx(:)
    ! EKT natural-orbital deflation threshold, matching GAMESS EKT-MRSF
    ! (thresh = 1.0d-2, sfgrad.src:3298).
    real(kind=dp), parameter :: ekt_occ_tol = 1.0e-2_dp
    real(kind=dp), parameter :: hartree_to_ev = 27.211386245988_dp
    real(kind=dp), allocatable :: fock_pack(:)
    type(dft_grid_t), target :: ea_molgrid
    type(scf_energy_t) :: ea_energy
    type(energy_results) :: saved_mol_energy
    real(kind=dp), allocatable :: density_alpha_ao(:), density_beta_ao(:)
    real(kind=dp), allocatable :: density_alpha_mo(:,:), density_beta_mo(:,:)
    real(kind=dp), allocatable :: density_ip_mo(:,:), density_ea_mo(:,:)
    real(kind=dp), allocatable :: smat(:)
    real(kind=dp), allocatable :: density_mo(:,:), fock_mo(:,:), lagrangian_mo(:,:)
    real(kind=dp), allocatable :: ea_density_ao(:,:), ea_fock_ao(:,:), saved_fock_ao(:,:)
    real(kind=dp), allocatable :: ekt_metric(:,:), ekt_operator(:,:)
    real(kind=dp), allocatable :: eig(:)
    real(kind=dp), allocatable :: orbitals(:,:), strengths(:), metric_norms(:)
    real(kind=dp), allocatable :: dom_mo_coeff(:), dom_no_coeff(:), dom_no_occ(:)
    real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:), dmat_a(:), dmat_b(:), mo_a(:,:), td_p(:,:), wao(:)
    real(kind=dp), contiguous, pointer :: density_store(:,:), lagrangian_store(:,:)
    real(kind=dp), contiguous, pointer :: fock_store(:,:), orbital_store(:,:), eig_store(:), strength_store(:)
    character(len=*), parameter :: tags_required(7) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_DM_A, OQP_DM_B, OQP_VEC_MO_A, OQP_td_p, OQP_WAO, OQP_td_bvec_mo /)
    character(len=*), parameter :: tags_alloc(6) = (/ character(len=80) :: &
      OQP_mrsf_ekt_density_mo, OQP_mrsf_ekt_lagrangian_mo, OQP_mrsf_ekt_fock_mo, &
      OQP_mrsf_ekt_orbitals_mo, OQP_mrsf_ekt_eigenvalues, OQP_mrsf_ekt_strengths /)

    ! Run the parent MRSF calculation and Z-vector first so the selected
    ! neutral state's relaxed density P and energy-weighted density /
    ! Lagrangian W are available.  The EKT step then builds the IP/EA
    ! generalized eigenproblem in the MRSF molecular-orbital basis.
    call tdhf_mrsf_energy(infos)
    call tdhf_mrsf_z_vector(infos)

    open(unit=iw, file=infos%log_filename, position="append")
    if (electron_affinity) then
      call print_module_info('MRSF_EKT_EA','Computing MRSF-EKT electron affinities')
    else
      call print_module_info('MRSF_EKT_IP','Computing MRSF-EKT ionization potentials')
    end if

    basis => infos%basis
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    ! Print/store all retained EKT Dyson roots available after the natural-
    ! orbital deflation, matching the GAMESS behavior.  The TDDFT NSTATE input
    ! still controls the parent MRSF calculation; it is not an EKT print cap.
    nroot = nbf
    select case (infos%control%scftype)
    case (1)
      nfocks = 1
    case default
      nfocks = 2
    end select

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)

    allocate(fock_pack(nbf2), &
             ea_density_ao(nbf2,nfocks), ea_fock_ao(nbf2,nfocks), saved_fock_ao(nbf2,nfocks), &
             density_alpha_ao(nbf2), density_beta_ao(nbf2), &
             density_alpha_mo(nbf,nbf), density_beta_mo(nbf,nbf), &
             density_ip_mo(nbf,nbf), density_ea_mo(nbf,nbf), &
             density_mo(nbf,nbf), fock_mo(nbf,nbf), &
             lagrangian_mo(nbf,nbf), ekt_metric(nbf,nbf), ekt_operator(nbf,nbf), &
             orbitals(nbf,nroot), strengths(nroot), smat(nbf2), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    ! MRSF-EKT uses the relaxed target-state one-particle density P and the
    ! relaxed energy-weighted density / Lagrangian W generated by the MRSF
    ! Z-vector.  OQP_td_p is the MRSF Z-vector relaxed *difference* density in
    ! the AO basis, spin-blocked (alpha,beta), in the same contravariant
    ! C*M*C^T convention as dmat and already 0.5-scaled in tdhf_mrsf_z_vector.
    ! Hence dmat_sigma + td_p(:,sigma) is the total relaxed one-particle
    ! density of the target (S0) response state for each spin.
    density_alpha_ao = dmat_a + td_p(:,1)
    density_beta_ao = dmat_b + td_p(:,2)

    ! AO -> MO transform of the one-particle DENSITY.  A density is
    ! contravariant: its natural-orbital (MO) representation is obtained with
    ! the metric-covariant form  D_MO = C^T S P S C, NOT C^T P C.  The SCF
    ! MOs satisfy C^T S C = I, so plain C^T P C mis-scales by the overlap and
    ! gives Tr(D_MO) != N_elec.  Operators (Fock, Lagrangian) keep the plain
    ! C^T A C transform below.
    call get_overlap_matrix(infos, nbf, smat)
    call density_ao_to_mo(nbf, nbf2, density_alpha_ao, smat, mo_a, density_alpha_mo)
    call density_ao_to_mo(nbf, nbf2, density_beta_ao,  smat, mo_a, density_beta_mo)
    density_ip_mo = density_alpha_mo
    density_ea_mo = density_beta_mo

    ! Acceptance checks (printed): AO electron counts Tr(P_sigma * S) and the
    ! MO-density traces must both equal nelec_sigma, and MO occupations must
    ! be physically bounded.
    call ekt_density_trace_checks(infos, nbf, nbf2, density_alpha_ao, &
         density_beta_ao, smat, density_alpha_mo, density_beta_mo, iw)

    call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, fock_pack)
    call unpack_matrix(fock_pack, fock_mo, nbf, 'U')

    ! ---------------------------------------------------------------------
    ! EKT-MRSF state-specific Lagrangian W~, reproduced exactly from the
    ! GAMESS-US reference implementation (sfgrad.src, the IF(MREKT) block;
    ! Pomogaev et al., JPCL 2021, eq 5):
    !
    !   W~_MO = -1/2 ( WMO + Lambda_MO )
    !
    ! where, in the MO basis and using the density metric C^T S M S C
    ! (a Lagrangian density is contravariant, like D):
    !   * Lambda = energy-weighted density Sum_k eps_k C_k C_k^T  (= eijden /
    !     GAMESS DAF record 36).  OpenQP's eijden halves the packed diagonal
    !     (a gradient-consumer convention, grd1.src:624); the GAMESS EKT path
    !     reads the raw record, so the half-diagonal is undone here.
    !   * WMO = the MRSF response Lagrangian.  OpenQP stores
    !     wao = 1/4 C WMO C^T (z-vector applies two halves), while GAMESS uses
    !     the bare WMO (records 419/429 = 1/2 C WMO C^T).  Hence
    !     WMO = 2 C^T S wao S C.
    ! ---------------------------------------------------------------------
    block
      use grd1, only: eijden
      real(kind=dp), allocatable :: lam_ao(:), lam_mo(:,:), wmo_mo(:,:)
      integer :: gok, kk, ijd
      allocate(lam_ao(nbf2), lam_mo(nbf,nbf), wmo_mo(nbf,nbf), &
               source=0.0_dp, stat=gok)
      if (gok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
      call eijden(lam_ao, nbf, infos)
      ijd = 0
      do kk = 1, nbf
        ijd = ijd + kk                       ! packed-upper diagonal index
        lam_ao(ijd) = 2.0_dp*lam_ao(ijd)     ! undo eijden's half-diagonal
      end do
      call density_ao_to_mo(nbf, nbf2, lam_ao, smat, mo_a, lam_mo)  ! C^T S Lambda S C
      call density_ao_to_mo(nbf, nbf2, wao,    smat, mo_a, wmo_mo)  ! C^T S wao S C
      wmo_mo = 2.0_dp*wmo_mo                  ! bare GAMESS WMO
      lagrangian_mo = -0.5_dp*(wmo_mo + lam_mo)
      lagrangian_mo = 0.5_dp*(lagrangian_mo + transpose(lagrangian_mo))
      deallocate(lam_ao, lam_mo, wmo_mo)
    end block

    if (electron_affinity) then
      ! EKT-EA: (F - W) * x = (I - P) * x * lambda.
      ! GAMESS rebuilds the Fock matrix for EA from the relaxed MRSF density
      ! (sfgrad.src:3269-3280, BUILDFOCK(..., MRDEA)) before forming F-W.
      ! Reproduce that path here instead of using the reference SCF Fock.
      ! GAMESS MRDEA uses the spin-averaged relaxed MRSF density
      ! P = 1/2*(Da+Db) for the EA complement metric, not beta alone.
      density_mo = 0.5_dp*(density_ip_mo + density_ea_mo)
      ekt_metric = -density_mo
      do i = 1, nbf
        ekt_metric(i,i) = ekt_metric(i,i) + 1.0_dp
      end do

      ! GAMESS MRDEA builds a closed-shell Fock from the spin-summed
      ! relaxed density: LDVAL = AO(Pavg), DSCAL(2), BUILDFOCK; internally
      ! its DFT path halves that total density for alpha=beta.
      ea_density_ao(:,1) = 0.5_dp*(density_alpha_ao + density_beta_ao)
      if (nfocks > 1) ea_density_ao(:,2) = ea_density_ao(:,1)
      saved_fock_ao(:,1) = fock_a
      if (nfocks > 1) then
        call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
        saved_fock_ao(:,2) = fock_b
      end if
      saved_mol_energy = infos%mol_energy

      nschwz = 0
      if (infos%control%hamilton >= 20) call dft_initialize(infos, basis, ea_molgrid)
      call calc_fock(basis, infos, ea_molgrid, ea_fock_ao, ea_energy, &
                     mo_a_in=mo_a, dens_in=ea_density_ao, nschwz=nschwz)

      ! calc_fock updates the stored SCF Fock and infos%mol_energy%energy as
      ! side effects; restore them because this EA-specific relaxed-density
      ! Fock is only an EKT operator ingredient, not a new SCF reference.
      fock_a = saved_fock_ao(:,1)
      if (nfocks > 1) fock_b = saved_fock_ao(:,2)
      infos%mol_energy = saved_mol_energy

      call orthogonal_transform_sym(nbf, nbf, ea_fock_ao(:,1), mo_a, nbf, fock_pack)
      call unpack_matrix(fock_pack, fock_mo, nbf, 'U')
      write(iw,'(/,2x,"EKT-EA: rebuilt relaxed-density Fock for EA operator (GAMESS MRDEA path)")')
      ekt_operator = fock_mo - lagrangian_mo
    else
      ! EKT: W * x = P * x * lambda.  Metric is the spin-summed relaxed
      ! density P = 1/2 (Da_MO + Db_MO) (GAMESS sfgrad.src:3115); for H2O S0
      ! this gives 5 occupied natural orbitals at occupation ~1.
      density_mo = 0.5_dp*(density_ip_mo + density_ea_mo)
      ekt_metric = density_mo
      ekt_operator = lagrangian_mo
    end if

    ! Natural-orbital occupation deflation + EKT generalized solve.
    ! The metric D_MO is non-diagonal, so a diag(D) > tol active-space
    ! selection is mathematically invalid.  Diagonalize the metric to obtain
    ! natural orbitals/occupations, keep occ > occ_tol, project W and D into
    ! the retained NO subspace, solve the generalized eigenproblem there,
    ! metric-normalize, and back-transform the Dyson orbitals.  Pole strengths
    ! come out <= 1 by construction.
    allocate(eig(nroot), metric_norms(nroot), dom_mo_coeff(nroot), dom_no_coeff(nroot), &
             dom_no_occ(nroot), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    allocate(dom_mo_idx(nroot), dom_no_idx(nroot), source=0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    ! Deflation threshold matches GAMESS EKT (thresh = 1.0d-2, sfgrad.src:3298):
    ! keep only natural orbitals with occupation > 1e-2, removing low-occupation
    ! null/relaxation channels.
    call solve_ekt_no_deflation(nbf, nroot, ekt_operator, ekt_metric, &
         ekt_occ_tol, eig, orbitals, strengths, metric_norms, dom_mo_idx, &
         dom_mo_coeff, dom_no_idx, dom_no_coeff, dom_no_occ, nactive, iw)

    call infos%dat%remove_records(tags_alloc)
    call infos%dat%reserve_data(OQP_mrsf_ekt_density_mo, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), &
        comment='MRSF-EKT metric density P in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_lagrangian_mo, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), &
        comment='MRSF-EKT Lagrangian W in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_fock_mo, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), &
        comment='MRSF-EKT Fock matrix in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_orbitals_mo, TA_TYPE_REAL64, nbf*nroot, (/ nbf, nroot /), &
        comment='MRSF-EKT Dyson-like orbital coefficients in MO basis')
    call infos%dat%reserve_data(OQP_mrsf_ekt_eigenvalues, TA_TYPE_REAL64, nroot, &
        comment='MRSF-EKT generalized-eigenproblem eigenvalues in Hartree')
    call infos%dat%reserve_data(OQP_mrsf_ekt_strengths, TA_TYPE_REAL64, nroot, &
        comment='MRSF-EKT pole strengths')

    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_density_mo, density_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_lagrangian_mo, lagrangian_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_fock_mo, fock_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_orbitals_mo, orbital_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_eigenvalues, eig_store)
    call tagarray_get_data(infos%dat, OQP_mrsf_ekt_strengths, strength_store)
    density_store = density_mo
    lagrangian_store = lagrangian_mo
    fock_store = fock_mo
    orbital_store = orbitals
    eig_store = eig
    strength_store = strengths

    if (electron_affinity) then
      write(iw,'(/,2x,"MRSF-EKT electron affinities (Dyson roots)")')
    else
      write(iw,'(/,2x,"MRSF-EKT ionization potentials (Dyson roots)")')
    end if
    write(iw,'(2x,"Rows index EKT Dyson roots/orbitals, not TDDFT excited states.")')
    ! eig holds the EKT eigenvalues epsilon of  W C = D C epsilon.  The
    ! electron binding energy (IP) of a detachment is -epsilon; print both
    ! the eigenvalue and the eBE in hartree and eV with the Dyson pole
    ! strength (norm of the Dyson vector, <= 1 for a physical root).
    write(iw,'(2x,"dyson",6x,"eig(ha)",8x,"eBE(ha)",8x,"eBE(eV)",7x,"metric",7x,"strength")')
    do i = 1, min(nroot, nactive)
      write(iw,'(2x,I5,2x,F14.6,2x,F14.6,2x,F12.4,2x,F12.6,2x,F12.6)') &
            i, eig(i), -eig(i), -eig(i)*hartree_to_ev, metric_norms(i), strengths(i)
    end do
    write(iw,'(/,2x,"--- EKT Dyson-root character diagnostics ---")')
    write(iw,'(2x,"dyson",2x,"dom_mo",2x,"mo_coeff",2x,"dom_no",2x,"no_coeff",2x,"no_occ",2x,"character",2x,"spurious")')
    do i = 1, min(nroot, nactive)
      write(iw,'(2x,I5,2x,I6,2x,F10.6,2x,I6,2x,F10.6,2x,F10.6,2x,A12,2x,L1)') &
            i, dom_mo_idx(i), dom_mo_coeff(i), dom_no_idx(i), dom_no_coeff(i), &
            dom_no_occ(i), ekt_root_character(dom_no_occ(i), electron_affinity), &
            ekt_is_spurious(-eig(i)*hartree_to_ev, metric_norms(i), strengths(i))
    end do
    write(iw,'(2x,"--------------------------------------")')

    if (electron_affinity) then
      write(iw,'(/,2x,"MRSF-EKT Dyson orbitals (AO-basis coefficients, GAMESS-style layout; columns = EA roots)")')
    else
      write(iw,'(/,2x,"MRSF-EKT Dyson orbitals (AO-basis coefficients, GAMESS-style layout; columns = IP roots)")')
    end if
    call print_dyson_orbitals(basis, mo_a, orbitals, eig, strengths, nbf, min(nroot, nactive))
    call flush(iw)

    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine tdhf_mrsf_ekt

  !> @brief Print Dyson orbitals in a GAMESS-like block.
  !> @detail The EKT solver returns Dyson-orbital coefficients in the OpenQP
  !>         MO basis.  For direct log-level comparison with GAMESS, transform
  !>         them back to the AO basis and print five roots per block with
  !>         ENERGY and STRENGTH header rows followed by AO basis-function
  !>         labels and coefficients.  This routine is output-only: it does not
  !>         change eigenvalues, strengths, root selection, or stored MO-basis
  !>         Dyson orbitals.
  subroutine print_dyson_orbitals(basis, mo_a, dyson_mo, eig, strengths, nbf, nroot)
    use precision, only: dp
    use io_constants, only: iw
    use basis_tools, only: basis_set

    implicit none

    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nbf, nroot
    real(kind=dp), intent(in) :: mo_a(nbf,nbf), dyson_mo(nbf,*), eig(*), strengths(*)

    integer, parameter :: mxlen = 5
    integer :: i, j, k, i0, i1
    real(kind=dp), allocatable :: dyson_ao(:,:)

    allocate(dyson_ao(nbf,nroot), source=0.0_dp)
    do i = 1, nroot
      do k = 1, nbf
        do j = 1, nbf
          dyson_ao(j,i) = dyson_ao(j,i) + mo_a(j,k)*dyson_mo(k,i)
        end do
      end do
    end do

    write(iw,'(/,10x,46("-"))')
    write(iw,'(12x,"EKT DYSON ORBITALS, ENERGIES AND NORMS")')
    write(iw,'(12x,"OpenQP AO-basis coefficients; GAMESS-style layout")')
    write(iw,'(10x,46("-"))')

    do i0 = 1, nroot, mxlen
      i1 = min(nroot, i0+mxlen-1)
      write(iw,'(/,14x,5i11)') (i, i=i0, i1)
      write(iw,'(5x,"ENERGY",3x,5f11.6)') (eig(i), i=i0, i1)
      write(iw,'(5x,"STRENGTH",1x,5f11.6)') (strengths(i), i=i0, i1)
      write(iw,'(14x,5(10x,a1))') ("A", i=i0, i1)
      do j = 1, nbf
        write(iw,'(i5,2x,a8,5f11.6)') j, basis%bf_label(j), (dyson_ao(j,i), i=i0, i1)
      end do
    end do

    deallocate(dyson_ao)

  end subroutine print_dyson_orbitals

  !> @brief Fetch the packed AO overlap matrix S (OQP_SM) into smat(nbf2).
  subroutine get_overlap_matrix(infos, nbf, smat)
    use types, only: information
    use oqp_tagarray_driver
    use precision, only: dp
    use messages, only: show_message, with_abort
    implicit none
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nbf
    real(kind=dp), intent(out) :: smat(:)
    real(kind=dp), contiguous, pointer :: s_p(:)
    character(len=*), parameter :: sn = "get_overlap_matrix"
    call data_has_tags(infos%dat, (/ character(len=80) :: OQP_SM /), &
         module_name, sn, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, s_p)
    smat(1:nbf*(nbf+1)/2) = s_p(1:nbf*(nbf+1)/2)
  end subroutine get_overlap_matrix

  !> @brief Metric-covariant AO->MO transform of a one-particle density:
  !>          D_MO = C^T (S P S) C
  !> @details A density is contravariant, so its MO (natural-orbital)
  !>          representation needs the overlap metric on both sides.  With
  !>          C^T S C = I this yields Tr(D_MO) = Tr(P S) = N_elec and MO
  !>          occupations bounded in [0,1] (per spin).  Plain C^T P C (used
  !>          for operators) is WRONG for a density.
  subroutine density_ao_to_mo(nbf, nbf2, p_ao_pack, s_pack, mo, d_mo)
    use precision, only: dp
    use mathlib, only: unpack_matrix, pack_matrix, orthogonal_transform_sym
    use messages, only: show_message, with_abort
    implicit none
    integer, intent(in) :: nbf, nbf2
    real(kind=dp), intent(in) :: p_ao_pack(:)   ! packed AO density (upper)
    real(kind=dp), intent(in) :: s_pack(:)      ! packed AO overlap (upper)
    real(kind=dp), intent(in) :: mo(nbf,nbf)    ! MO coefficients C (AO x MO)
    real(kind=dp), intent(out) :: d_mo(nbf,nbf) ! MO-basis density
    real(kind=dp), allocatable :: p_ao(:,:), s_ao(:,:), sps(:,:), tmp(:,:), sps_pack(:)
    integer :: ok
    allocate(p_ao(nbf,nbf), s_ao(nbf,nbf), sps(nbf,nbf), tmp(nbf,nbf), &
             sps_pack(nbf2), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    call unpack_matrix(p_ao_pack, p_ao, nbf, 'U')
    call unpack_matrix(s_pack, s_ao, nbf, 'U')
    tmp = matmul(s_ao, p_ao)        ! S P
    sps = matmul(tmp, s_ao)         ! S P S  (symmetric)
    sps = 0.5_dp*(sps + transpose(sps))
    call pack_matrix(sps, sps_pack, 'U')
    call orthogonal_transform_sym(nbf, nbf, sps_pack, mo, nbf, sps_pack)
    call unpack_matrix(sps_pack, d_mo, nbf, 'U')
    deallocate(p_ao, s_ao, sps, tmp, sps_pack)
  end subroutine density_ao_to_mo

  !> @brief Acceptance/trace checks for the EKT relaxed density (printed).
  !> @details Verifies Tr(P_sigma S) = nelec_sigma in the AO basis and that
  !>          the MO-density trace matches, with bounded MO occupations.
  subroutine ekt_density_trace_checks(infos, nbf, nbf2, pa_ao_pack, pb_ao_pack, &
       s_pack, da_mo, db_mo, iw)
    use types, only: information
    use precision, only: dp
    use mathlib, only: unpack_matrix
    use messages, only: show_message, with_abort
    implicit none
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nbf, nbf2, iw
    real(kind=dp), intent(in) :: pa_ao_pack(:), pb_ao_pack(:), s_pack(:)
    real(kind=dp), intent(in) :: da_mo(nbf,nbf), db_mo(nbf,nbf)
    real(kind=dp), allocatable :: pa(:,:), pb(:,:), s_ao(:,:)
    real(kind=dp) :: trPaS, trPbS, trDa, trDb, occ_min, occ_max
    integer :: i, j, nea, neb, ok
    allocate(pa(nbf,nbf), pb(nbf,nbf), s_ao(nbf,nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    call unpack_matrix(pa_ao_pack, pa, nbf, 'U')
    call unpack_matrix(pb_ao_pack, pb, nbf, 'U')
    call unpack_matrix(s_pack, s_ao, nbf, 'U')
    nea = int(infos%mol_prop%nelec_A)
    neb = int(infos%mol_prop%nelec_B)
    trPaS = 0.0_dp; trPbS = 0.0_dp
    do i = 1, nbf
      do j = 1, nbf
        trPaS = trPaS + pa(i,j)*s_ao(j,i)
        trPbS = trPbS + pb(i,j)*s_ao(j,i)
      end do
    end do
    trDa = 0.0_dp; trDb = 0.0_dp; occ_min = 1.0d30; occ_max = -1.0d30
    do i = 1, nbf
      trDa = trDa + da_mo(i,i)
      trDb = trDb + db_mo(i,i)
      occ_min = min(occ_min, da_mo(i,i), db_mo(i,i))
      occ_max = max(occ_max, da_mo(i,i), db_mo(i,i))
    end do
    write(iw,'(/,2x,"--- EKT relaxed-density acceptance checks ---")')
    write(iw,'(2x,"Tr(Pa S) = ",f12.6,"   Tr(Da_MO) = ",f12.6)') trPaS, trDa
    write(iw,'(2x,"Tr(Pb S) = ",f12.6,"   Tr(Db_MO) = ",f12.6)') trPbS, trDb
    write(iw,'(2x,"(AO and MO traces must agree; target-state singlet S0)")')
    write(iw,'(2x,"MO occupations range: [",f10.6,",",f10.6,"]")') occ_min, occ_max
    write(iw,'(2x,"--------------------------------------------",/)')
    call flush(iw)
    deallocate(pa, pb, s_ao)
  end subroutine ekt_density_trace_checks

  !> @brief EKT solve with natural-orbital occupation deflation.
  !> @details Solves  W X = D X eps  where the metric D (= ekt_metric) is the
  !>          relaxed one-particle density (IP) or its particle complement (EA),
  !>          and W (= ekt_operator) is the EKT Lagrangian.  Because D is
  !>          non-diagonal, deflation must be done on its NATURAL ORBITALS, not
  !>          its diagonal:
  !>            1. symmetrize D
  !>            2. diagonalize  D U = U n   (natural orbitals U, occupations n)
  !>            3. keep NOs with n_i > occ_tol
  !>            4. project  D_NO = U_keep^T D U_keep,  W_NO = U_keep^T W U_keep
  !>            5. solve the generalized problem in the retained NO space
  !>            6. metric-normalize eigenvectors:  X^T D_NO X = 1
  !>            7. back-transform Dyson orbitals  C = U_keep X
  !>            8. pole strength = X^T D_NO X  (<= 1 physical)
  subroutine solve_ekt_no_deflation(nbf, nroot, wmat, dmat, occ_tol, &
       eig, dyson_mo, strengths, metric_norms, dom_mo_idx, dom_mo_coeff, &
       dom_no_idx, dom_no_coeff, dom_no_occ, nkeep, iw)
    use precision, only: dp
    use eigen, only: diag_symm_full
    use messages, only: show_message, with_abort
    implicit none
    integer, intent(in) :: nbf, nroot, iw
    real(kind=dp), intent(in) :: wmat(nbf,nbf), dmat(nbf,nbf), occ_tol
    real(kind=dp), intent(out) :: eig(nroot)
    real(kind=dp), intent(out) :: dyson_mo(nbf,nroot)
    real(kind=dp), intent(out) :: strengths(nroot)
    real(kind=dp), intent(out) :: metric_norms(nroot)
    integer, intent(out) :: dom_mo_idx(nroot), dom_no_idx(nroot)
    real(kind=dp), intent(out) :: dom_mo_coeff(nroot), dom_no_coeff(nroot), dom_no_occ(nroot)
    integer, intent(out) :: nkeep
    real(kind=dp), allocatable :: dsym(:,:), nocc(:), uno(:,:), ukeep(:,:)
    integer, allocatable :: kept_idx(:)
    real(kind=dp), allocatable :: d_no(:,:), w_no(:,:), tmp(:,:)
    real(kind=dp), allocatable :: xvec(:,:), eps(:), cdys(:,:), d_times_x(:)
    real(kind=dp) :: tr_disc, occ_min, occ_max, xnorm, str, coeff_abs, best_abs
    integer :: i, j, m, ierr, ok, nout, best_idx, out_idx, n_skipped
    ! EKT eigenvalue print/retention threshold.  GAMESS skips near-zero EKT
    ! eigenvalues (|eps| <= thresh) when transferring solved EKT roots to its
    ! output (sfgrad.src:3411, thresh = 1.0d-2 Ha); mirror that here so the
    ! printed/stored Dyson ladder aligns with GAMESS without a one-root shift.
    ! Kept separate from the NO occupation tolerance (ekt_occ_tol) which only
    ! governs natural-orbital deflation.
    real(kind=dp), parameter :: ekt_eig_tol = 1.0e-2_dp

    ! (1) symmetrize the metric
    allocate(dsym(nbf,nbf), nocc(nbf), uno(nbf,nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    do j = 1, nbf
      do i = 1, nbf
        dsym(i,j) = 0.5_dp*(dmat(i,j) + dmat(j,i))
      end do
    end do

    ! (2) diagonalize the metric -> natural orbitals U, occupations n
    uno = dsym
    call diag_symm_full(1, nbf, uno, nbf, nocc, ierr)
    if (ierr /= 0) call show_message('EKT: NO diagonalization failed', WITH_ABORT)

    ! (3) keep NOs with physical occupation n_i > occ_tol
    m = 0
    tr_disc = 0.0_dp
    do i = 1, nbf
      if (nocc(i) > occ_tol) then
        m = m + 1
      else
        tr_disc = tr_disc + nocc(i)
      end if
    end do
    if (m == 0) call show_message('EKT: no NOs above occupation tolerance', WITH_ABORT)

    allocate(ukeep(nbf,m), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    allocate(kept_idx(m), source=0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    j = 0
    occ_min = 1.0d30; occ_max = -1.0d30
    do i = 1, nbf
      if (nocc(i) > occ_tol) then
        j = j + 1
        ukeep(:,j) = uno(:,i)
        kept_idx(j) = i
        occ_min = min(occ_min, nocc(i))
        occ_max = max(occ_max, nocc(i))
      end if
    end do

    ! (4) deflation diagnostics
    write(iw,'(/,2x,"--- EKT natural-orbital deflation ---")')
    write(iw,'(2x,"occ_tol = ",es10.2,"   retained NOs = ",i0," / ",i0)') occ_tol, m, nbf
    write(iw,'(2x,"discarded trace = ",f12.6,"   retained occ range = [",f10.6,",",f10.6,"]")') &
          tr_disc, occ_min, occ_max
    write(iw,'(2x,"eig_tol = ",es10.2," Ha   (skip EKT roots with |eps| <= eig_tol; GAMESS sfgrad.src)")') &
          ekt_eig_tol
    write(iw,'(2x,"-------------------------------------",/)')

    ! (5) project D and W into the retained NO space
    allocate(d_no(m,m), w_no(m,m), tmp(nbf,m), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    tmp  = matmul(dsym, ukeep)
    d_no = matmul(transpose(ukeep), tmp)
    tmp  = matmul(wmat, ukeep)
    w_no = matmul(transpose(ukeep), tmp)
    d_no = 0.5_dp*(d_no + transpose(d_no))
    w_no = 0.5_dp*(w_no + transpose(w_no))

    ! (6) solve the generalized problem in the retained NO space
    allocate(xvec(m,m), eps(m), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    call solve_symmetric_generalized(w_no, d_no, eps, xvec, m, ierr)
    if (ierr /= 0) call show_message('(A,I0)', &
         'EKT: NO-space generalized solve failed; info=', ierr, WITH_ABORT)

    ! (7) metric-normalize eigenvectors:  X^T D_NO X = 1
    do j = 1, m
      xnorm = 0.0_dp
      do i = 1, m
        xnorm = xnorm + xvec(i,j)*dot_product(d_no(i,:), xvec(:,j))
      end do
      if (xnorm > 1.0e-14_dp) xvec(:,j) = xvec(:,j)/sqrt(xnorm)
    end do

    ! (8) back-transform Dyson orbitals to the MO basis  C = U_keep X.
    !     Metric norms are X^T D_NO X after normalization (=1 for retained
    !     physical roots).  EKT pole strengths follow the GAMESS/Pomogaev
    !     definition X^T D_NO^2 X; weak Dyson roots therefore retain small
    !     strengths instead of being forced to 1 by metric normalization.
    allocate(cdys(nbf,m), d_times_x(m), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    cdys = matmul(ukeep, xvec)

    eig = 0.0_dp
    dyson_mo = 0.0_dp
    strengths = 0.0_dp
    metric_norms = 0.0_dp
    dom_mo_idx = 0
    dom_no_idx = 0
    dom_mo_coeff = 0.0_dp
    dom_no_coeff = 0.0_dp
    dom_no_occ = 0.0_dp
    ! Transfer solved EKT roots to the output arrays, skipping near-zero EKT
    ! eigenvalues the GAMESS way (|eps| <= ekt_eig_tol).  The source eigenpair
    ! index is j (xvec/cdys/eps); out_idx is a separate compressed counter into
    ! the output arrays so the printed/stored ladder matches GAMESS.
    out_idx = 0
    n_skipped = 0
    do j = 1, m
      if (abs(eps(j)) <= ekt_eig_tol) then
        n_skipped = n_skipped + 1
        cycle
      end if
      out_idx = out_idx + 1
      if (out_idx > nroot) exit
      eig(out_idx) = eps(j)
      dyson_mo(:,out_idx) = cdys(:,j)
      str = 0.0_dp
      do i = 1, m
        str = str + xvec(i,j)*dot_product(d_no(i,:), xvec(:,j))
      end do
      metric_norms(out_idx) = str
      d_times_x = matmul(d_no, xvec(:,j))
      strengths(out_idx) = dot_product(d_times_x, d_times_x)

      best_abs = -1.0_dp
      best_idx = 0
      do i = 1, nbf
        coeff_abs = abs(cdys(i,j))
        if (coeff_abs > best_abs) then
          best_abs = coeff_abs
          best_idx = i
        end if
      end do
      dom_mo_idx(out_idx) = best_idx
      if (best_idx > 0) dom_mo_coeff(out_idx) = cdys(best_idx,j)

      best_abs = -1.0_dp
      best_idx = 0
      do i = 1, m
        coeff_abs = abs(xvec(i,j))
        if (coeff_abs > best_abs) then
          best_abs = coeff_abs
          best_idx = i
        end if
      end do
      if (best_idx > 0) then
        dom_no_idx(out_idx) = kept_idx(best_idx)
        dom_no_coeff(out_idx) = xvec(best_idx,j)
        dom_no_occ(out_idx) = nocc(kept_idx(best_idx))
      end if
    end do
    nout = min(out_idx, nroot)
    ! The caller uses nkeep as nactive when printing rows/Dyson orbitals, so it
    ! must reflect the number of stored output roots (nout), not the retained
    ! NO-space dimension (m); otherwise a zero row is printed for skipped roots.
    nkeep = nout
    write(iw,'(2x,"EKT eigenvalue filter: skipped ",i0," near-zero root(s) ", &
         &"(|eps| <= ",es10.2," Ha); stored ",i0," / ",i0," roots")') &
          n_skipped, ekt_eig_tol, nout, nroot

    deallocate(dsym, nocc, uno, ukeep, kept_idx, d_no, w_no, tmp, xvec, eps, cdys, d_times_x)
  end subroutine solve_ekt_no_deflation

  pure logical function ekt_is_spurious(ebe_ev, metric_norm, pole_strength) result(flag)
    use precision, only: dp
    implicit none
    real(kind=dp), intent(in) :: ebe_ev, metric_norm, pole_strength
    flag = (abs(ebe_ev) > 1.0e4_dp) .or. (metric_norm < 0.0_dp) .or. &
           (metric_norm > 1.05_dp) .or. (pole_strength < 0.0_dp) .or. &
           (pole_strength > 1.05_dp)
  end function ekt_is_spurious

  pure function ekt_root_character(no_occ, electron_affinity) result(label)
    use precision, only: dp
    implicit none
    real(kind=dp), intent(in) :: no_occ
    logical, intent(in) :: electron_affinity
    character(len=12) :: label
    if (electron_affinity) then
      if (no_occ < 0.20_dp) then
        label = 'virtual'
      else if (no_occ > 0.80_dp) then
        label = 'occupied'
      else
        label = 'mixed'
      end if
    else
      if (no_occ > 0.80_dp) then
        label = 'occupied'
      else if (no_occ < 0.20_dp) then
        label = 'virtual'
      else
        label = 'mixed'
      end if
    end if
  end function ekt_root_character

  !> @brief Solve the EKT generalized eigenproblem  op * C = metric * C * lambda
  !> @details The EKT working equation (eq 1 of Park et al., JCTC 2024) is a
  !>          symmetric generalized eigenproblem in which the metric is the
  !>          relaxed one-particle density (IP) or its particle complement (EA).
  !>          That metric is symmetric positive semidefinite but in general NOT
  !>          diagonal: the MRSF Z-vector relaxation introduces off-diagonal
  !>          occupation couplings in the ground-state MO basis.  We therefore
  !>          reduce the problem by Loewdin (symmetric) orthogonalization built
  !>          from the full eigendecomposition of the metric, discarding the
  !>          null-space directions whose occupation eigenvalues fall below a
  !>          tolerance.  This is the rigorous replacement for a diagonal-only
  !>          metric scaling, which is exact only when the metric is diagonal.
  subroutine solve_symmetric_generalized(operator, metric, eigenvalues, eigenvectors, n, ierr)
    use precision, only: dp
    use eigen, only: diag_symm_full
    use messages, only: show_message, with_abort
    implicit none
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: operator(n,n)
    real(kind=dp), intent(in) :: metric(n,n)
    real(kind=dp), intent(out) :: eigenvalues(n), eigenvectors(n,n)
    integer, intent(out) :: ierr
    integer :: i, j, m, ok
    real(kind=dp), parameter :: metric_eval_tol = 1.0e-10_dp
    real(kind=dp), allocatable :: smat(:,:), seval(:), xorth(:,:)
    real(kind=dp), allocatable :: opsym(:,:), tmp(:,:), reduced(:,:)
    real(kind=dp), allocatable :: vec(:,:), eval_red(:)

    ierr = 0
    eigenvalues = 0.0_dp
    eigenvectors = 0.0_dp
    if (n <= 0) return

    ! Eigendecomposition of the (symmetrized) metric:  S = U diag(seval) U^T
    allocate(smat(n,n), seval(n), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    do j = 1, n
      do i = 1, n
        smat(i,j) = 0.5_dp*(metric(i,j) + metric(j,i))
      end do
    end do
    call diag_symm_full(1, n, smat, n, seval, ierr)
    if (ierr /= 0) return

    ! Keep only directions with positive occupation (rank of the metric).
    m = 0
    do i = 1, n
      if (seval(i) > metric_eval_tol) m = m + 1
    end do
    if (m == 0) then
      ierr = -1
      return
    end if

    ! Symmetric orthogonalizer  X = U_kept * diag(seval_kept)^(-1/2)   (n x m)
    allocate(xorth(n,m), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    j = 0
    do i = 1, n
      if (seval(i) > metric_eval_tol) then
        j = j + 1
        xorth(:,j) = smat(:,i) / sqrt(seval(i))
      end if
    end do

    ! Reduced operator in the orthonormal metric basis:  reduced = X^T op X
    allocate(opsym(n,n), tmp(n,m), reduced(m,m), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    do j = 1, n
      do i = 1, n
        opsym(i,j) = 0.5_dp*(operator(i,j) + operator(j,i))
      end do
    end do
    tmp = matmul(opsym, xorth)
    reduced = matmul(transpose(xorth), tmp)
    reduced = 0.5_dp*(reduced + transpose(reduced))

    ! Standard symmetric eigenproblem on the reduced operator.
    allocate(vec(m,m), eval_red(m), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    vec = reduced
    call diag_symm_full(1, m, vec, m, eval_red, ierr)
    if (ierr /= 0) return

    ! Map results back:  eigenvalues are lambda, eigenvectors C = X * vec.
    eigenvalues(1:m) = eval_red(1:m)
    eigenvectors(:,1:m) = matmul(xorth, vec)
  end subroutine solve_symmetric_generalized

end module tdhf_mrsf_ekt_mod
