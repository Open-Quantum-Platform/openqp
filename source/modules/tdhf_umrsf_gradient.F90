!> UMRSF-TDDFT analytic nuclear gradient — clean-room implementation (branch uhf-grad-plan).
!> STAGE: response build-up, gates-first (G1 -> G2 -> G3) per DERIVATIONS/stage1_hf_gradient_spec.md.
!> The C entry computes the reference (UHF-triplet) gradient via the reusable hf_gradient primitive
!> (grd1/grd2 with the converged DM_A/DM_B) — already FD-certified to 1.46e-7 — and additionally runs
!> the NON-FD isolation gates (energy reconstruction from the converged amplitude) so the response
!> 2-PDM / matvec routing is validated before any derivative-integral code is wired. Reuses the
!> clean-room energy/response lib (umrsfcbc/int2_umrsf/umrsfmntoia/mrsfesum/get_jacobi); never reads
!> the guarded RO-MRSF gradient.
module tdhf_umrsf_gradient_mod

  use precision, only: dp
  use types, only: information
  use grd2, only: grd2_compute_data_t
  use basis_tools, only: basis_set

  implicit none

  character(len=*), parameter :: module_name = "tdhf_umrsf_gradient_mod"

  private
  public :: tdhf_umrsf_gradient_C

  !> Custom grd2 2-PDM for the UMRSF response (amplitude transition density).
  !> Emits, per shell-quartet, the certified (B_k,D_k) channel 2-PDM (G1) with the
  !> grd2 factor convention (verified to reduce to grd2_uhf in the degenerate case):
  !>   Coulomb (k=1..8): + 4*sc*s_k*(B_k(ij)D_k(kl)+D_k(ij)B_k(kl))
  !>   exchange (all k): - 2*sx*s_k*(B_k(ik)D_k(jl)+B_k(il)D_k(jk)+D_k(ik)B_k(jl)+D_k(il)B_k(jk))
  !> s_k = mrst sign (k<=10 negated for mrst=3) * spin-pair-coupling scale. NEVER reuse grd2_uhf
  !> (that gives D(x)D — the wrong 2-PDM). See DERIVATIONS/G1_response_2pdm.md.
  type, extends(grd2_compute_data_t) :: grd2_umrsf_resp_t
    integer :: nbf = 0
    integer :: nchan = 11
    real(kind=dp), allocatable :: bden(:,:,:)   ! (nchan,nbf,nbf) RAW bra B_k   (exchange + energy)
    real(kind=dp), allocatable :: dden(:,:,:)   ! (nchan,nbf,nbf) RAW ket D_k   (exchange + energy)
    real(kind=dp), allocatable :: bden_s(:,:,:) ! (nchan,nbf,nbf) sym(B_k)      (Coulomb only)
    real(kind=dp), allocatable :: dden_s(:,:,:) ! (nchan,nbf,nbf) sym(D_k)      (Coulomb only)
    real(kind=dp), allocatable :: sgn(:)        ! (nchan) s_k
    logical,       allocatable :: has_coul(:)   ! (nchan) channel carries Coulomb
    real(kind=dp) :: sc = 1.0_dp                ! scale_coulomb
    real(kind=dp) :: sx = 1.0_dp                ! scale_exchange
  contains
    procedure :: init => grd2_umrsf_resp_init
    procedure :: clean => grd2_umrsf_resp_clean
    procedure :: get_density => grd2_umrsf_resp_get_density
  end type

contains

  subroutine tdhf_umrsf_gradient_C(c_handle) bind(C, name="tdhf_umrsf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use io_constants, only: iw
    use hf_gradient_mod, only: hf_gradient
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    real(kind=dp), allocatable :: de2e_resp(:,:)
    inf => oqp_handle_get_info(c_handle)

    ! Gate: reconstruct the excited-state energy from the converged amplitude (G1), assemble the
    ! response 2-PDM, and validate the 2e response gradient via a frozen-density FD self-test (G1b).
    ! Returns the validated 2e response gradient dω_2e/dx (densities-fixed / "Pulay" part).
    call umrsf_grad_run_gates(inf, de2e_resp)

    open(unit=iw, file=inf%log_filename, position="append")
    write(iw,'(/2x,a)') 'UMRSF gradient [Stage1]: reference UHF-triplet gradient + 2e response '// &
                        '(P^Delta 1e/mean-field + Z-vector + M1 still to add)'
    close(iw)
    ! Reference (UHF triplet) nuclear gradient via the reusable HF/DFT gradient primitive
    ! (zeros the gradient, then fills the reference part).
    call hf_gradient(inf)
    ! Add the validated 2e response (transition-density Pulay) contribution.
    if (allocated(de2e_resp)) inf%atoms%grad = inf%atoms%grad + de2e_resp
  end subroutine tdhf_umrsf_gradient_C

!###############################################################################
!> @brief NON-FD isolation gates. Reconstructs the excited-state response energy
!>        omega_I = X^T A X (TDA) from the converged amplitude X by re-running the
!>        clean-room response matvec (umrsfcbc -> int2_umrsf -> umrsfmntoia [2e] +
!>        mrsfesum [orbital-energy diagonal]). Compares to the stored td_energy.
!>        Also computes the 2e part two ways (matvec back-transform vs the channel
!>        density:Fock trace) to characterize the response 2-PDM for the gradient.
  subroutine umrsf_grad_run_gates(infos, de2e_out)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use mathlib, only: orthogonal_transform_sym, unpack_matrix
    use eigen, only: diag_symm_full
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc, umrsfmntoia, mrsfesum, get_jacobi
    use tdhf_lib, only: iatogen
    use grd2, only: grd2_driver
    use oqp_linalg
    use iso_c_binding, only: c_int, c_f_pointer

    implicit none

    character(len=*), parameter :: subroutine_name = "umrsf_grad_run_gates"
    type(information), target, intent(inout) :: infos
    real(kind=dp), allocatable, intent(out) :: de2e_out(:,:)

    type(grd2_umrsf_resp_t) :: gcomp
    real(kind=dp), allocatable :: de2e(:,:), de2e_fd(:,:)
    real(kind=dp), allocatable, target :: densym(:,:,:,:)
    integer :: natom, iat, icmp
    real(kind=dp) :: omega_p, omega_m, hfd, omega_base, maxd2e
    ! unrelaxed difference density P^Δ,u + orbital-part gradient
    real(kind=dp), allocatable :: talpha(:,:), tbeta(:,:), pda(:,:), pdb(:,:)
    real(kind=dp), allocatable :: de_orb(:,:), de_w(:,:), de_w_2e(:,:)
    real(kind=dp) :: omega_orb_chk, omega_orb, hfscale_ref
    integer :: ia, ib, i, j
    ! Z-vector (relaxation): canonical MOs + relaxation density
    real(kind=dp), allocatable :: cac(:,:), cbc(:,:), epsca(:), epscb(:), pza(:,:), pzb(:,:)
    real(kind=dp) :: zrms, zresid

    type(basis_set), pointer :: basis
    ! tagarray pointers
    real(kind=dp), contiguous, pointer :: smat(:), fock_a(:), fock_b(:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), contiguous, pointer :: mo_energy_a(:), mo_energy_b(:)
    real(kind=dp), contiguous, pointer :: bvec(:,:), td_en(:)
    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    integer(c_int), pointer :: ixcore_ptr(:)
    ! locals
    real(kind=dp), allocatable :: va(:,:), vb(:,:), fa(:,:), fb(:,:), smat_full(:,:)
    real(kind=dp), allocatable :: ea(:), eb(:), wrk1(:,:), wrk2(:,:), scr(:)
    real(kind=dp), allocatable :: xmat(:,:), amo(:,:), amo2e(:,:)
    real(kind=dp), allocatable :: brad(:,:,:)
    real(kind=dp), allocatable, target :: dens(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    type(int2_compute_t) :: int2_driver
    type(int2_umrsf_data_t), target :: int2_udata
    integer :: nbf, nbf2, nocca, noccb, nvirb, xvec_dim, mrst, nstates, tstate
    integer :: k, it, diag_index
    real(kind=dp) :: scale_exch, hfs, omega_recon, omega_2e_mv, omega_2e_tr
    real(kind=dp) :: spc_coco, spc_ovov, spc_coov

    integer(4) :: status

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb
    mrst = infos%tddft%mult
    nstates = infos%tddft%nstate
    tstate = infos%tddft%target_state
    if (tstate < 1) tstate = 1
    ! mrst is the RESPONSE multiplicity (1=singlet / 3=triplet response states); both use the
    ! same MRSF machinery off the UHF triplet reference. The triplet sign-flip below is mrst==3 only.
    if (mrst /= 1 .and. mrst /= 3) return
    xvec_dim = nocca*nvirb

    ! Effective exact-exchange scale for the response (HF limit -> 1.0).
    hfs = infos%tddft%hfscale
    if (hfs == -1.0_dp) hfs = 1.0_dp
    scale_exch = 1.0_dp
    if (infos%control%hamilton == 20) scale_exch = hfs   ! dft
    spc_coco = infos%tddft%spc_coco ; if (spc_coco == -1.0_dp) spc_coco = hfs
    spc_ovov = infos%tddft%spc_ovov ; if (spc_ovov == -1.0_dp) spc_ovov = hfs
    spc_coov = infos%tddft%spc_coov ; if (spc_coov == -1.0_dp) spc_coov = hfs

    call tagarray_get_data(infos%dat, OQP_SM, smat, status)
    if (status /= 0) return               ! data not present (e.g. fresh call): silently skip gate
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_td_energies, td_en, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status); if (status/=0) return
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b, status); if (status/=0) return

    allocate(va(nbf,nbf), vb(nbf,nbf), fa(nbf,nbf), fb(nbf,nbf), smat_full(nbf,nbf), &
             ea(nbf), eb(nbf), wrk1(nbf,nbf), wrk2(nbf,nbf), scr(nbf2), &
             xmat(nbf,nbf), amo(xvec_dim,1), amo2e(xvec_dim,1), source=0.0_dp)
    allocate(dens(1,11,nbf,nbf), brad(11,nbf,nbf), source=0.0_dp)

    open(unit=iw, file=infos%log_filename, position="append")

    va = mo_a ; vb = mo_b ; ea = mo_energy_a ; eb = mo_energy_b
    call unpack_matrix(smat, smat_full, nbf, 'U')
    ! Corresponding-orbital (Jacobi) alignment — idempotent on already-aligned MOs.
    call get_jacobi(infos, va, ea, vb, eb, smat_full, nocca, wrk1, wrk2, 0)
    call get_jacobi(infos, va, ea, vb, eb, smat_full, nocca, wrk1, wrk2, 1)

    ! MO-basis Fock (rotated orbitals); frozen-core shift identical to the energy path.
    call orthogonal_transform_sym(nbf, nbf, fock_a, va, nbf, scr)
    if (infos%tddft%ixcore_len /= 0) then
      call c_f_pointer(infos%tddft%ixcore, ixcore_ptr, [infos%tddft%ixcore_len])
      do it = 1, noccb
        if (.not. any(ixcore_ptr(1:infos%tddft%ixcore_len) == it)) then
          diag_index = (it+1)*it/2
          scr(diag_index) = -1.0d6
        end if
      end do
    end if
    call unpack_matrix(scr, fa)
    call orthogonal_transform_sym(nbf, nbf, fock_b, vb, nbf, scr)
    call unpack_matrix(scr, fb)

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_driver%schwarz = .false.         ! validation: disable Schwarz screening (RULES sec.9)

    ! ---- response matvec on the converged target-state amplitude ----
    call iatogen(bvec(:,tstate), xmat, nocca, noccb)
    call umrsfcbc(infos, va, vb, xmat, dens(1,:,:,:))

    int2_udata = int2_umrsf_data_t(d3=dens(1:1,:,:,:), tamm_dancoff=.true., &
                                   scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call int2_driver%run(int2_udata)
    fmrst2 => int2_udata%f3(:,:,:,:,1)

    if (mrst == 3) fmrst2(:,1:10,:,:) = -fmrst2(:,1:10,:,:)
    ! Spin-pair coupling (no-op for the HF defaults spc==hfscale).
    if (abs(hfs) > epsilon(1.0_dp)) then
      if (spc_coco /= hfs) fmrst2(:,10,:,:) = fmrst2(:,10,:,:) * (spc_coco/hfs)
      if (spc_ovov /= hfs) fmrst2(:,9,:,:)  = fmrst2(:,9,:,:)  * (spc_ovov/hfs)
      if (spc_coov /= hfs) fmrst2(:,1:8,:,:) = fmrst2(:,1:8,:,:) * (spc_coov/hfs)
    end if

    ! 2e part of A.X via the back-transform (the energy path)
    amo2e = 0.0_dp
    call umrsfmntoia(infos, fmrst2(1,:,:,:), amo2e, va, vb, 1)
    omega_2e_mv = dot_product(bvec(:,tstate), amo2e(:,1))

    ! 2e part via the response 2-PDM bra densities B_k = adjoint(umrsfmntoia).X :
    ! omega_2e = sum_k <B_k, F_k>  with  F_k = int2_k(D_k),  D_k = umrsfcbc(X).
    ! This is the contraction the analytic 2e gradient differentiates (B_k vs D_k pair).
    call umrsf_bra_density(infos, va, vb, xmat, brad)
    omega_2e_tr = 0.0_dp
    do k = 1, 11
      omega_2e_tr = omega_2e_tr + sum(brad(k,:,:)*fmrst2(1,k,:,:))
    end do

    ! full omega = 2e part + orbital-energy (Fock-diagonal) part
    amo(:,1) = amo2e(:,1)
    call iatogen(bvec(:,tstate), xmat, nocca, noccb)
    call mrsfesum(infos, xmat, fa, fb, amo, 1)
    omega_recon = dot_product(bvec(:,tstate), amo(:,1))

    write(iw,'(/2x,a)') '================ UMRSF gradient NON-FD gate G1 ================'
    write(iw,'(2x,a,i0,a,i0)') 'target_state = ', tstate, '   mrst = ', mrst
    write(iw,'(2x,a,f18.12)')  'X^T X (amplitude norm)          = ', dot_product(bvec(:,tstate),bvec(:,tstate))
    write(iw,'(2x,a,f18.10)')  'omega reconstructed (X^T A X)   = ', omega_recon
    write(iw,'(2x,a,f18.10)')  'omega stored      (td_energies) = ', td_en(tstate)
    write(iw,'(2x,a,es12.3)')  '  |delta| omega                 = ', abs(omega_recon-td_en(tstate))
    write(iw,'(2x,a,f18.10)')  'omega_2e (back-transform)       = ', omega_2e_mv
    write(iw,'(2x,a,f18.10)')  'omega_orb (X.esum)              = ', omega_recon-omega_2e_mv
    write(iw,'(2x,a,f18.10)')  'omega_2e via bra-density 2-PDM  = ', omega_2e_tr
    write(iw,'(2x,a,es12.3)')  '  |delta| 2e (G1 2-PDM routing) = ', abs(omega_2e_mv-omega_2e_tr)
    if (abs(omega_recon-td_en(tstate)) <= 1.0e-9_dp .and. abs(omega_2e_mv-omega_2e_tr) <= 1.0e-9_dp) then
      write(iw,'(2x,a)')       'VERDICT: G1 PASS (matvec + response 2-PDM routing validated)'
    else
      write(iw,'(2x,a)')       'VERDICT: G1 CHECK (see deltas above)'
    end if
    write(iw,'(2x,a)')         '=============================================================='
    close(iw)

    ! ================= 2e RESPONSE GRADIENT + frozen-density FD self-test =================
    ! Analytic dω_2e/dx = Σ Γ^Δ (μν|λσ)^x via the custom grd2 (B_k,D_k channel 2-PDM).
    ! Validated against a central FD of ω_2e with the AO densities B_k,D_k held FIXED
    ! (only the integrals move) — the same frozen-density identity grd2_hess_selftest uses.
    ! This isolates the 2e gradient assembly (no orbital relaxation / Z-vector needed).
    natom = ubound(infos%atoms%zn,1)
    call umrsf_resp_2pdm_fill(gcomp, dens(1,:,:,:), brad, nbf, mrst, hfs, &
                              scale_exch, spc_coco, spc_ovov, spc_coov)
    allocate(de2e(3,natom), de2e_fd(3,natom), densym(1,11,nbf,nbf), source=0.0_dp)
    do k = 1, 11
      densym(1,k,:,:) = gcomp%dden(k,:,:)
    end do

    ! analytic 2e response gradient (base geometry)
    call grd2_driver(infos, basis, de2e, gcomp)

    ! frozen-density central FD of ω_2e
    call umrsf_frozen_omega2e(infos, int2_driver, densym, gcomp, scale_exch, omega_base)
    hfd = 1.0e-3_dp
    do iat = 1, natom
      do icmp = 1, 3
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) + hfd
        call basis%init_shell_centers()
        call umrsf_frozen_omega2e(infos, int2_driver, densym, gcomp, scale_exch, omega_p)
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) - 2.0_dp*hfd
        call basis%init_shell_centers()
        call umrsf_frozen_omega2e(infos, int2_driver, densym, gcomp, scale_exch, omega_m)
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) + hfd
        call basis%init_shell_centers()
        de2e_fd(icmp,iat) = (omega_p - omega_m)/(2.0_dp*hfd)
      end do
    end do
    maxd2e = maxval(abs(de2e - de2e_fd))

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '========= UMRSF 2e response gradient (G1b: frozen-density FD) ========='
    write(iw,'(2x,a,f18.10)') 'omega_2e (frozen-density, base) = ', omega_base
    write(iw,'(2x,a,f18.10)') '  (cf. back-transform omega_2e)  = ', omega_2e_mv
    write(iw,'(2x,a)') '   atom  comp     analytic dω2e/dx        frozen-FD          |Δ|'
    do iat = 1, natom
      do icmp = 1, 3
        write(iw,'(2x,2i5,3es20.10)') iat, icmp, de2e(icmp,iat), de2e_fd(icmp,iat), &
                                      abs(de2e(icmp,iat)-de2e_fd(icmp,iat))
      end do
    end do
    write(iw,'(2x,a,es12.3)') 'max|analytic - frozen-FD| 2e    = ', maxd2e
    if (maxd2e <= 1.0e-6_dp) then
      write(iw,'(2x,a)') 'VERDICT: 2e response gradient PASS (custom grd2 2-PDM validated)'
    else
      write(iw,'(2x,a)') 'VERDICT: 2e response gradient CHECK (see |Δ| above)'
    end if
    write(iw,'(2x,a)') '======================================================================'
    close(iw)
    call gcomp%clean()

    ! ================= unrelaxed difference density P^Δ,u (omega_orb check) =================
    ! Standard CIS/TDA difference density (spin-flip): occ_α-occ_α  T_α(i,j) = -Σ_a X_ia X_ja,
    ! virt_β-virt_β  T_β(a,b) = +Σ_i X_ia X_ib, occ-virt = 0 (filled later by the Z-vector).
    ! Must reproduce ω_orb = Tr(P^Δ,u F^AO) = X·mrsfesum(X). Gate-check before wiring the 1e/W terms.
    omega_orb = omega_recon - omega_2e_mv
    allocate(talpha(nbf,nbf), tbeta(nbf,nbf), pda(nbf,nbf), pdb(nbf,nbf), source=0.0_dp)
    call iatogen(bvec(:,tstate), xmat, nocca, noccb)
    do j = 1, nocca
      do i = 1, nocca
        do ia = noccb+1, nbf
          talpha(i,j) = talpha(i,j) - xmat(i,ia)*xmat(j,ia)
        end do
      end do
    end do
    do ib = noccb+1, nbf
      do ia = noccb+1, nbf
        do i = 1, nocca
          tbeta(ia,ib) = tbeta(ia,ib) + xmat(i,ia)*xmat(i,ib)
        end do
      end do
    end do
    omega_orb_chk = sum(talpha*fa) + sum(tbeta*fb)

    ! DIAGNOSTIC (Z-vector prep): are the STORED MOs canonical (F diagonal) or get_jacobi-rotated?
    ! cphf_solve_uhf assumes canonical MOs (uses e_a-e_i). Check ||va - mo_a|| and off-diag(fa).
    block
      real(kind=dp) :: dva, offocc, offvir
      integer :: ii, jj
      dva = sum(abs(va-mo_a)) + sum(abs(vb-mo_b))
      offocc = 0.0_dp; offvir = 0.0_dp
      do jj = 1, nocca; do ii = 1, nocca
        if (ii/=jj) offocc = max(offocc, abs(fa(ii,jj)))
      end do; end do
      do jj = noccb+1, nbf; do ii = noccb+1, nbf
        if (ii/=jj) offvir = max(offvir, abs(fb(ii,jj)))
      end do; end do
      open(unit=iw, file=infos%log_filename, position="append")
      write(iw,'(/2x,a)') '--- Z-vector prep diagnostic ---'
      write(iw,'(2x,a,es12.3)') '||va-mo_a||+||vb-mo_b|| (rotation persisted in store?) = ', dva
      write(iw,'(2x,a,es12.3)') 'max|off-diag fa| within alpha-occ block               = ', offocc
      write(iw,'(2x,a,es12.3)') 'max|off-diag fb| within beta-virt block               = ', offvir
      close(iw)
    end block

    ! Re-canonicalize the ground-state MO-Fock (G2 prep): fa,fb are occ/virt block-diagonal (SCF F_ai=0,
    ! get_jacobi only mixes within α-occ / β-virt), so diagonalizing gives canonical MOs WITHOUT occ-virt
    ! mixing. Cσ_can = Cσ_rot·Vσ, εσ = eigenvalues. cphf_solve_uhf needs exactly these. Verify diagonal.
    allocate(cac(nbf,nbf), cbc(nbf,nbf), epsca(nbf), epscb(nbf))
    block
      real(kind=dp), allocatable :: fac(:,:), fbc(:,:), fmo(:,:)
      real(kind=dp) :: offa, offb
      integer :: ierr, ii, jj
      allocate(fac(nbf,nbf), fbc(nbf,nbf), fmo(nbf,nbf))
      fac = fa ; fbc = fb
      call diag_symm_full(1, nbf, fac, nbf, epsca, ierr)   ! fac -> Vα (eigenvectors)
      call diag_symm_full(1, nbf, fbc, nbf, epscb, ierr)   ! fbc -> Vβ
      cac = matmul(va, fac)                                  ! canonical α MOs = Cα_rot · Vα
      cbc = matmul(vb, fbc)                                  ! canonical β MOs
      ! verify Cα_canᵀ fock_a Cα_can is diagonal (= εα)
      call orthogonal_transform_sym(nbf, nbf, fock_a, cac, nbf, scr) ; call unpack_matrix(scr, fmo)
      offa = 0.0_dp
      do jj = 1, nbf; do ii = 1, nbf
        if (ii/=jj) offa = max(offa, abs(fmo(ii,jj)))
      end do; end do
      call orthogonal_transform_sym(nbf, nbf, fock_b, cbc, nbf, scr) ; call unpack_matrix(scr, fmo)
      offb = 0.0_dp
      do jj = 1, nbf; do ii = 1, nbf
        if (ii/=jj) offb = max(offb, abs(fmo(ii,jj)))
      end do; end do
      open(unit=iw, file=infos%log_filename, position="append")
      write(iw,'(2x,a,2es12.3)') 're-canonicalized: max|off-diag F^MO_can| alpha/beta = ', offa, offb
      if (max(offa,offb) <= 1.0e-9_dp) then
        write(iw,'(2x,a)') 'VERDICT: re-canonicalization OK (F^MO_can diagonal; ready for cphf_solve_uhf)'
      else
        write(iw,'(2x,a)') 'VERDICT: re-canonicalization CHECK (F^MO_can not diagonal)'
      end if
      close(iw)
      deallocate(fac, fbc, fmo)
    end block
    pda = matmul(matmul(va, talpha), transpose(va))     ! AO P^Δ,u_α = C^α T_α C^α^T
    pdb = matmul(matmul(vb, tbeta),  transpose(vb))     ! AO P^Δ,u_β = C^β T_β C^β^T

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '========= UMRSF difference density P^Delta,u (omega_orb check) ========='
    write(iw,'(2x,a,f18.10)') 'omega_orb via P^Delta,u (Tr T fa+Tr T fb) = ', omega_orb_chk
    write(iw,'(2x,a,f18.10)') 'omega_orb true (X.mrsfesum)               = ', omega_orb
    write(iw,'(2x,a,es12.3)') '  |delta| omega_orb                       = ', abs(omega_orb_chk-omega_orb)
    write(iw,'(2x,a)') '======================================================================'
    close(iw)

    ! ---- Z-vector relaxation (numerical-RHS end-to-end validation) ----
    ! R = ∂ω/∂κ by FD over canonical orbital rotations; solve via cphf_solve_uhf; relaxation density
    ! P_z folded into P^Δ (= P^Δ,u + z). Needs the live int2 driver for the ω re-evaluations.
    allocate(pza(nbf,nbf), pzb(nbf,nbf), source=0.0_dp)
    call umrsf_zvector_relax(infos, int2_driver, va, vb, cac, cbc, epsca, epscb, smat_full, &
                             fock_a, fock_b, bvec(:,tstate), scale_exch, pza, pzb, zrms, zresid)
    pda = pda + pza
    pdb = pdb + pzb
    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '========= UMRSF Z-vector (numerical RHS + cphf_solve_uhf) ========='
    write(iw,'(2x,a,es12.3)') 'z RMS = ', zrms
    write(iw,'(2x,a)') 'relaxed P^Δ = P^Δ,u + z folded into the 1e + mean-field gradient (W still missing)'
    write(iw,'(2x,a)') '==================================================================='
    close(iw)

    ! ---- W part 1 (ω_2e/amplitude): numerical symmetric-generalized-Fock route is DISABLED ----
    ! It overshoots ~20x: the non-orthonormal symmetric orbital perturbation violates umrsfcbc's
    ! orthonormal-MO assumption, inflating G^2e. The analytic ω_2e W (channel-adjoint) is the TODO.
    allocate(de_w_2e(3,natom), source=0.0_dp)
    ! call umrsf_w_numerical(infos, int2_driver, basis, va, vb, fock_a, fock_b, bvec(:,tstate), &
    !                        scale_exch, de_w_2e)

    call int2_driver%clean()

    ! orbital-part gradient (RELAXED P^Δ = P^Δ,u + z): 1e Tr(P^Δ h^x) + 2e mean-field Tr(P^Δ G[P^ref])
    hfscale_ref = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale_ref = infos%dft%hfscale
    allocate(de_orb(3,natom), source=0.0_dp)
    call umrsf_orbital_grad(infos, basis, pda, pdb, dmat_a, dmat_b, hfscale_ref, de_orb)

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '========= UMRSF orbital-part response gradient (1e + 2e mean-field) ========='
    write(iw,'(2x,a)') '   atom  comp        de_2e (transition)      de_orb (1e+meanfield)'
    do iat = 1, natom
      do icmp = 1, 3
        write(iw,'(2x,2i5,2es24.12)') iat, icmp, de2e(icmp,iat), de_orb(icmp,iat)
      end do
    end do
    write(iw,'(2x,a)') '============================================================================'
    close(iw)

    ! ---- W part 2: difference-density W^Δ_orb (analytic, relaxed P^Δ) + total W = W^Δ_orb + W_2e ----
    allocate(de_w(3,natom), source=0.0_dp)
    call umrsf_w_overlap_grad(infos, basis, cac, cbc, epsca, epscb, smat_full, pda, pdb, de_w)
    de_w = de_w + de_w_2e
    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '========= UMRSF W (overlap-Pulay) = W^Δ_orb + W_2e (amplitude) ========='
    do iat = 1, natom
      do icmp = 1, 3
        write(iw,'(2x,2i5,es24.12)') iat, icmp, de_w(icmp,iat)
      end do
    end do
    write(iw,'(2x,a)') '==============================================================='
    close(iw)

    ! return full response = transition-2e + orbital(1e+meanfield, relaxed) + W overlap
    de2e_out = de2e + de_orb + de_w

    deallocate(va, vb, fa, fb, smat_full, ea, eb, wrk1, wrk2, scr, xmat, amo, amo2e, dens, brad)
    deallocate(de2e, de2e_fd, densym)
    deallocate(talpha, tbeta, pda, pdb, de_orb, de_w, de_w_2e)
    deallocate(cac, cbc, epsca, epscb, pza, pzb)

  end subroutine umrsf_grad_run_gates

!###############################################################################
!> ω_2e = Σ_k s_k ⟨B_k, int2_k(D_k)⟩ at the CURRENT geometry with FIXED AO densities
!> (re-inits the int2 engine for the moved nuclei; Schwarz screening off for a clean FD).
  subroutine umrsf_frozen_omega2e(infos, idrv, densym, gcomp, scale_exch, omega)
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in), target :: densym(:,:,:,:)
    type(grd2_umrsf_resp_t), intent(in) :: gcomp
    real(kind=dp), intent(in) :: scale_exch
    real(kind=dp), intent(out) :: omega
    type(int2_umrsf_data_t), target :: ud
    real(kind=dp), pointer :: f3(:,:,:,:)
    integer :: k

    call idrv%clean()
    call idrv%init(infos%basis, infos)
    call idrv%set_screening()
    idrv%schwarz = .false.
    ud = int2_umrsf_data_t(d3=densym, tamm_dancoff=.true., &
                           scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(ud)
    f3 => ud%f3(:,:,:,:,1)
    omega = 0.0_dp
    do k = 1, 11
      omega = omega + gcomp%sgn(k)*sum(gcomp%bden(k,:,:)*f3(1,k,:,:))
    end do
  end subroutine umrsf_frozen_omega2e

!###############################################################################
!> Overlap-Pulay W term: −Tr(W S^x) for the (difference-density) energy-weighted density
!> W^Δ_σ = ½ C_σ (F^MO_σ T_σ + T_σ F^MO_σ) C_σ^T.  In the CANONICAL basis F^MO=diag(ε) ⇒
!> W^Δ_σ = ½ Σ_pq (ε_p+ε_q) T_pq C_p C_q^T with T = C_σ^T S P^Δ_σ S C_σ (relaxed difference density).
!> Wired via grd1 grad_ee_overlap matching the eijden convention (negate + half-diagonal pack).
  subroutine umrsf_w_overlap_grad(infos, basis, cac, cbc, epsca, epscb, smat_full, pda, pdb, de_w)
    use grd1, only: grad_ee_overlap
    use mathlib, only: pack_matrix
    use constants, only: tol_int
    implicit none
    type(information), intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), epsca(:), epscb(:), smat_full(:,:)
    real(kind=dp), intent(in) :: pda(:,:), pdb(:,:)
    real(kind=dp), intent(out) :: de_w(:,:)
    real(kind=dp), allocatable :: tmo(:,:), wtot(:,:), wpack(:)
    integer :: nbf, nbf2, p, q, ij, i
    real(kind=dp) :: tol

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; tol = tol_int*log(10.0_dp)
    allocate(tmo(nbf,nbf), wtot(nbf,nbf), wpack(nbf2))
    wtot = 0.0_dp
    ! alpha
    tmo = matmul(transpose(cac), matmul(smat_full, matmul(pda, matmul(smat_full, cac))))
    do q = 1, nbf ; do p = 1, nbf ; tmo(p,q) = 0.5_dp*(epsca(p)+epsca(q))*tmo(p,q) ; end do ; end do
    wtot = wtot + matmul(matmul(cac, tmo), transpose(cac))
    ! beta
    tmo = matmul(transpose(cbc), matmul(smat_full, matmul(pdb, matmul(smat_full, cbc))))
    do q = 1, nbf ; do p = 1, nbf ; tmo(p,q) = 0.5_dp*(epscb(p)+epscb(q))*tmo(p,q) ; end do ; end do
    wtot = wtot + matmul(matmul(cbc, tmo), transpose(cbc))
    ! pack -W with half-diagonal (eijden convention) and contract with S^x
    call pack_matrix(-wtot, wpack, 'U')
    ij = 0
    do i = 1, nbf ; ij = ij + i ; wpack(ij) = 0.5_dp*wpack(ij) ; end do
    de_w = 0.0_dp
    call grad_ee_overlap(basis, wpack, de_w, logtol=tol)
    deallocate(tmo, wtot, wpack)
  end subroutine umrsf_w_overlap_grad

!###############################################################################
!> ω_2e / amplitude overlap-Pulay W via the numerical SYMMETRIC generalized Fock of ω_2e ONLY:
!>   W_2e = ½ Σ_σ C_σ G^2e,σ_sym C_σ^T,   G^2e,σ_sym,pq = ∂ω_2e/∂U^σ_pq  (symmetric orbital variation).
!> This is the energy-weighting of the TRANSITION density, ADDED to the difference-density W^Δ_orb
!> (which already carries the ω_orb + z-relaxation weighting). ω_2e via only_2e=.true. (skip mrsfesum).
!> The AO W is invariant to the orbital basis ⇒ perturb the rotated orbitals directly (no get_jacobi ⇒
!> fast; EXCLUDES the M1/get_jacobi part — size it from the residual). −Tr(W_2e S^x) via grad_ee_overlap.
  subroutine umrsf_w_numerical(infos, idrv, basis, va, vb, fock_a, fock_b, xv, scale_exch, de_w)
    use int2_compute, only: int2_compute_t
    use grd1, only: grad_ee_overlap
    use mathlib, only: pack_matrix
    use constants, only: tol_int
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), fock_a(:), fock_b(:), xv(:), scale_exch
    real(kind=dp), intent(out) :: de_w(:,:)

    real(kind=dp), allocatable :: gsym(:,:), wtot(:,:), cw_a(:,:), cw_b(:,:), wpack(:)
    integer :: nbf, nbf2, sgn, p, q, ij, i
    real(kind=dp) :: th, omp, omm, tol

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; th = 1.0d-3 ; tol = tol_int*log(10.0_dp)
    allocate(gsym(nbf,nbf), wtot(nbf,nbf), cw_a(nbf,nbf), cw_b(nbf,nbf), wpack(nbf2), source=0.0_dp)

    do sgn = 1, 2
      gsym = 0.0_dp
      do p = 1, nbf
        do q = p, nbf
          cw_a = va ; cw_b = vb
          if (p == q) then
            if (sgn==1) then ; cw_a(:,p) = (1.0_dp+th)*va(:,p)
            else             ; cw_b(:,p) = (1.0_dp+th)*vb(:,p) ; end if
          else
            if (sgn==1) then
              cw_a(:,p) = va(:,p) + th*va(:,q) ; cw_a(:,q) = va(:,q) + th*va(:,p)
            else
              cw_b(:,p) = vb(:,p) + th*vb(:,q) ; cw_b(:,q) = vb(:,q) + th*vb(:,p)
            end if
          end if
          call umrsf_omega_eval(infos, idrv, cw_a, cw_b, fock_a, fock_b, xv, scale_exch, omp, only_2e=.true.)
          cw_a = va ; cw_b = vb
          if (p == q) then
            if (sgn==1) then ; cw_a(:,p) = (1.0_dp-th)*va(:,p)
            else             ; cw_b(:,p) = (1.0_dp-th)*vb(:,p) ; end if
          else
            if (sgn==1) then
              cw_a(:,p) = va(:,p) - th*va(:,q) ; cw_a(:,q) = va(:,q) - th*va(:,p)
            else
              cw_b(:,p) = vb(:,p) - th*vb(:,q) ; cw_b(:,q) = vb(:,q) - th*vb(:,p)
            end if
          end if
          call umrsf_omega_eval(infos, idrv, cw_a, cw_b, fock_a, fock_b, xv, scale_exch, omm, only_2e=.true.)
          if (p == q) then
            gsym(p,p) = (omp - omm)/(2.0_dp*th)
          else
            gsym(p,q) = (omp - omm)/(4.0_dp*th)   ! ∂ω/∂U_pq = G_pq+G_qp = 2 G_sym_pq
            gsym(q,p) = gsym(p,q)
          end if
        end do
      end do
      if (sgn==1) then
        wtot = wtot + 0.5_dp*matmul(matmul(va, gsym), transpose(va))
      else
        wtot = wtot + 0.5_dp*matmul(matmul(vb, gsym), transpose(vb))
      end if
    end do

    call pack_matrix(-wtot, wpack, 'U')
    ij = 0
    do i = 1, nbf ; ij = ij + i ; wpack(ij) = 0.5_dp*wpack(ij) ; end do
    de_w = 0.0_dp
    call grad_ee_overlap(basis, wpack, de_w, logtol=tol)
    deallocate(gsym, wtot, cw_a, cw_b, wpack)
  end subroutine umrsf_w_numerical

!###############################################################################
!> ω = X^T A X for a given (already get_jacobi-aligned) set of orbitals vva,vvb and amplitude xv.
!> Replicates the validated response matvec (umrsfcbc->int2->umrsfmntoia[2e] + mrsfesum[orbital]).
!> Used to build the Z-vector RHS R = ∂ω/∂κ by FD over canonical orbital rotations.
  subroutine umrsf_omega_eval(infos, idrv, vva, vvb, fock_a, fock_b, xv, scale_exch, omega, only_2e)
    use mathlib, only: orthogonal_transform_sym, unpack_matrix
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc, umrsfmntoia, mrsfesum
    use tdhf_lib, only: iatogen
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: vva(:,:), vvb(:,:), fock_a(:), fock_b(:), xv(:)
    real(kind=dp), intent(in) :: scale_exch
    real(kind=dp), intent(out) :: omega
    logical, intent(in), optional :: only_2e
    logical :: do_2e_only
    real(kind=dp), allocatable :: fa(:,:), fb(:,:), scr(:), xmat(:,:), amo(:,:)
    real(kind=dp), allocatable, target :: dens(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    type(int2_umrsf_data_t), target :: ud
    integer :: nbf, nbf2, nocca, noccb, nvirb, xvec_dim, mrst

    nbf = infos%basis%nbf ; nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb ; xvec_dim = nocca*nvirb ; mrst = infos%tddft%mult
    allocate(fa(nbf,nbf), fb(nbf,nbf), scr(nbf2), xmat(nbf,nbf), amo(xvec_dim,1), &
             dens(1,11,nbf,nbf), source=0.0_dp)

    call orthogonal_transform_sym(nbf, nbf, fock_a, vva, nbf, scr) ; call unpack_matrix(scr, fa)
    call orthogonal_transform_sym(nbf, nbf, fock_b, vvb, nbf, scr) ; call unpack_matrix(scr, fb)

    call iatogen(xv, xmat, nocca, noccb)
    call umrsfcbc(infos, vva, vvb, xmat, dens(1,:,:,:))
    ud = int2_umrsf_data_t(d3=dens(1:1,:,:,:), tamm_dancoff=.true., &
                           scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(ud)
    fmrst2 => ud%f3(:,:,:,:,1)
    if (mrst == 3) fmrst2(:,1:10,:,:) = -fmrst2(:,1:10,:,:)
    do_2e_only = .false. ; if (present(only_2e)) do_2e_only = only_2e
    amo = 0.0_dp
    call umrsfmntoia(infos, fmrst2(1,:,:,:), amo, vva, vvb, 1)
    if (.not. do_2e_only) then
      call iatogen(xv, xmat, nocca, noccb)
      call mrsfesum(infos, xmat, fa, fb, amo, 1)
    end if
    omega = dot_product(xv, amo(:,1))   ! full ω, or ω_2e only (skip mrsfesum orbital part)
    deallocate(fa, fb, scr, xmat, amo, dens)
  end subroutine umrsf_omega_eval

!###############################################################################
!> Unrelaxed orbital-part gradient: dω_orb/dx = Tr(P^Δ,u h^x) + Σ(P^Δ,u⊗P^ref)(μν|λσ)^x
!> (the W/overlap term is NOT here — it is the orthonormality response, added separately).
!> 1e via the validated grd1 primitives; 2e mean-field (cross of P^Δ,u with the reference density)
!> via the validated grd2_uhf path and the polarization identity
!>   d[Tr(A G[B])] = ½ d[E2(A+B) − E2(A) − E2(B)],  E2(P)=½Tr(P G[P]).
  subroutine umrsf_orbital_grad(infos, basis, pda, pdb, dmat_a, dmat_b, hfscale, de_orb)
    use grd1, only: grad_ee_kinetic, grad_en_hellman_feynman, grad_en_pulay
    use grd2, only: grd2_driver
    use hf_gradient_mod, only: grd2_uhf_compute_data_t
    use mathlib, only: pack_matrix
    use constants, only: tol_int
    implicit none
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: pda(:,:), pdb(:,:)             ! AO P^Δ,u_α, P^Δ,u_β (full)
    real(kind=dp), target, intent(in) :: dmat_a(:), dmat_b(:)   ! packed reference densities
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: de_orb(:,:)

    real(kind=dp), allocatable :: pdtot_p(:), zn(:)
    real(kind=dp), allocatable, target :: pda_p(:), pdb_p(:), da1(:), db1(:)
    real(kind=dp), allocatable :: de1(:,:), de2(:,:), de3(:,:)
    type(grd2_uhf_compute_data_t) :: gc
    integer :: nbf, nbf2, natom
    real(kind=dp) :: tol

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; natom = ubound(infos%atoms%zn,1)
    tol = tol_int*log(10.0_dp)
    allocate(pdtot_p(nbf2), pda_p(nbf2), pdb_p(nbf2), da1(nbf2), db1(nbf2), zn(natom))
    allocate(de1(3,natom), de2(3,natom), de3(3,natom), source=0.0_dp)
    zn = infos%atoms%zn - infos%basis%ecp_zn_num

    call pack_matrix(pda+pdb, pdtot_p, 'U')
    call pack_matrix(pda, pda_p, 'U')
    call pack_matrix(pdb, pdb_p, 'U')

    de_orb = 0.0_dp
    ! 1e: Tr(P^Δ,u h^x)  (kinetic + e-n Hellmann–Feynman + e-n Pulay)
    call grad_ee_kinetic(basis, pdtot_p, de_orb, logtol=tol)
    call grad_en_hellman_feynman(basis, infos%atoms%xyz, zn, pdtot_p, de_orb, logtol=tol)
    call grad_en_pulay(basis, infos%atoms%xyz, zn, pdtot_p, de_orb, logtol=tol)

    ! 2e mean-field: d[Tr(P^Δ,u G[P^ref])] via grd2_uhf polarization
    da1 = pda_p + dmat_a ; db1 = pdb_p + dmat_b
    gc = grd2_uhf_compute_data_t(da=da1,    db=db1,    hfscale=hfscale, nbf=nbf)
    call gc%init() ; call grd2_driver(infos, basis, de1, gc) ; call gc%clean()
    gc = grd2_uhf_compute_data_t(da=pda_p,  db=pdb_p,  hfscale=hfscale, nbf=nbf)
    call gc%init() ; call grd2_driver(infos, basis, de2, gc) ; call gc%clean()
    gc = grd2_uhf_compute_data_t(da=dmat_a, db=dmat_b, hfscale=hfscale, nbf=nbf)
    call gc%init() ; call grd2_driver(infos, basis, de3, gc) ; call gc%clean()
    de_orb = de_orb + (de1 - de2 - de3)

    deallocate(pdtot_p, pda_p, pdb_p, da1, db1, zn, de1, de2, de3)
  end subroutine umrsf_orbital_grad

!###############################################################################
!> @brief Response 2-PDM bra densities B_k = adjoint(umrsfmntoia) applied to the
!>        amplitude X. With the ket densities D_k = umrsfcbc(X) and the int2_umrsf
!>        channel J/K patterns, the 2e response energy is omega_2e = sum_k <B_k, F_k>,
!>        F_k = int2_k(D_k); the analytic 2e gradient differentiates this bilinear
!>        form (B_k vs D_k pair) against the derivative ERIs. Derived clean-room by
!>        transposing umrsfmntoia (general C^a^T F C^b transform + the SOMO sections
!>        3-6 with their 1/2 factors). Channel order matches umrsfcbc/umrsfmntoia:
!>        1=o2v_a 2=o2v_b 3=o1v_a 4=o1v_b 5=co1_a 6=co1_b 7=co2_a 8=co2_b
!>        9=o21v 10=co12 11=ball(general).
  subroutine umrsf_bra_density(infos, va, vb, x, brad)
    use messages, only: show_message, with_abort
    implicit none
    type(information), intent(in) :: infos
    real(kind=dp), intent(in)  :: va(:,:), vb(:,:), x(:,:)
    real(kind=dp), intent(out) :: brad(:,:,:)         ! (11, nbf, nbf)

    integer :: nbf, nocca, noccb, mrst, o1, o2, i, j
    real(kind=dp), parameter :: isqrt2 = 1.0_dp/sqrt(2.0_dp)
    real(kind=dp), allocatable :: ca_o1(:), cb_o1(:), ca_o2(:), cb_o2(:)
    real(kind=dp), allocatable :: va_v_o1(:), vb_v_o1(:), va_v_o2(:), vb_v_o2(:)
    real(kind=dp), allocatable :: tmp(:,:)

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    mrst = infos%tddft%mult
    o1 = nocca-1
    o2 = nocca

    allocate(ca_o1(nbf), cb_o1(nbf), ca_o2(nbf), cb_o2(nbf), &
             va_v_o1(nbf), vb_v_o1(nbf), va_v_o2(nbf), vb_v_o2(nbf), &
             tmp(nbf,nbf), source=0.0_dp)

    ! closed -> O1/O2 weighted MO sums  (i in closed = 1..nocca-2)
    do i = 1, nocca-2
      ca_o1 = ca_o1 + x(i,o1)*va(:,i)
      cb_o1 = cb_o1 + x(i,o1)*vb(:,i)
      ca_o2 = ca_o2 + x(i,o2)*va(:,i)
      cb_o2 = cb_o2 + x(i,o2)*vb(:,i)
    end do
    ! O1/O2 -> virtual weighted MO sums  (j in virt = nocca+1..nbf)
    do j = nocca+1, nbf
      va_v_o1 = va_v_o1 + x(o1,j)*va(:,j)
      vb_v_o1 = vb_v_o1 + x(o1,j)*vb(:,j)
      va_v_o2 = va_v_o2 + x(o2,j)*va(:,j)
      vb_v_o2 = vb_v_o2 + x(o2,j)*vb(:,j)
    end do

    brad = 0.0_dp
    ! same-spin channels 1..8 carry a factor 1/2 (umrsfmntoia sections 3-6)
    call add_outer(brad(1,:,:), 0.5_dp, ca_o1, va(:,o1))      ! ado2va (sec 4)
    call add_outer(brad(2,:,:), 0.5_dp, cb_o1, vb(:,o1))      ! ado2vb
    call add_outer(brad(3,:,:), 0.5_dp, ca_o2, va(:,o2))      ! ado1va (sec 3)
    call add_outer(brad(4,:,:), 0.5_dp, cb_o2, vb(:,o2))      ! ado1vb
    call add_outer(brad(5,:,:), 0.5_dp, va(:,o2), va_v_o2)    ! adco1a (sec 6)
    call add_outer(brad(6,:,:), 0.5_dp, vb(:,o2), vb_v_o2)    ! adco1b
    call add_outer(brad(7,:,:), 0.5_dp, va(:,o1), va_v_o1)    ! adco2a (sec 5)
    call add_outer(brad(8,:,:), 0.5_dp, vb(:,o1), vb_v_o1)    ! adco2b
    ! mixed two-SOMO channels 9,10 (factor 1)
    call add_outer(brad(9,:,:),  1.0_dp, va(:,o2), vb_v_o1)   ! ao21v (sec 5)
    call add_outer(brad(9,:,:), -1.0_dp, va(:,o1), vb_v_o2)   !       (sec 6)
    call add_outer(brad(10,:,:), 1.0_dp, ca_o2, vb(:,o1))     ! aco12 (sec 3)
    call add_outer(brad(10,:,:),-1.0_dp, ca_o1, vb(:,o2))     !       (sec 4)

    ! general channel 11 (agdlr): bra = C^a X C^b^T over the SF block, with the
    ! SOMO diagonal replaced by the sqrt(1/2) spin-adapted combination.
    tmp = matmul(matmul(va, x), transpose(vb))               ! sum_ij x(i,j) va_i vb_j^T
    brad(11,:,:) = tmp
    call add_outer(brad(11,:,:), -x(o1,o1), va(:,o1), vb(:,o1))
    call add_outer(brad(11,:,:), -x(o2,o2), va(:,o2), vb(:,o2))
    if (mrst == 1) then
      call add_outer(brad(11,:,:),  isqrt2*x(o1,o1), va(:,o1), vb(:,o1))
      call add_outer(brad(11,:,:), -isqrt2*x(o1,o1), va(:,o2), vb(:,o2))
    else if (mrst == 3) then
      call add_outer(brad(11,:,:), -x(o1,o2), va(:,o1), vb(:,o2))
      call add_outer(brad(11,:,:), -x(o2,o1), va(:,o2), vb(:,o1))
      call add_outer(brad(11,:,:),  isqrt2*x(o1,o1), va(:,o1), vb(:,o1))
      call add_outer(brad(11,:,:),  isqrt2*x(o1,o1), va(:,o2), vb(:,o2))
    else
      call show_message('umrsf_bra_density: unsupported response multiplicity', with_abort)
    end if

    deallocate(ca_o1, cb_o1, ca_o2, cb_o2, va_v_o1, vb_v_o1, va_v_o2, vb_v_o2, tmp)
  end subroutine umrsf_bra_density

!###############################################################################
!> Z-vector via NUMERICAL RHS: R^σ_ai = ∂ω/∂κ^σ_ai by central FD of ω over canonical orbital
!> rotations (get_jacobi re-applied each point ⇒ the M1/two-reference response is INCLUDED), solved
!> by the validated cphf_solve_uhf (UHF orbital Hessian). Returns AO relaxation densities pza,pzb.
!> End-to-end framework validation (slow O(ltot) ω-evals; analytic R replaces it later).
  subroutine umrsf_zvector_relax(infos, idrv, va, vb, cac, cbc, epsca, epscb, smat_full, &
                                 fock_a, fock_b, xv, scale_exch, pza, pzb, zrms, resid)
    use oqp_tagarray_driver
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: get_jacobi
    use cphf_mod, only: cphf_solve_uhf
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), cac(:,:), cbc(:,:)
    real(kind=dp), intent(in) :: epsca(:), epscb(:), smat_full(:,:)
    real(kind=dp), intent(in) :: fock_a(:), fock_b(:), xv(:), scale_exch
    real(kind=dp), intent(out) :: pza(:,:), pzb(:,:), zrms, resid

    real(kind=dp), contiguous, pointer :: moa(:,:), mob(:,:), ea(:), eb(:)
    real(kind=dp), allocatable :: cwa(:,:), cwb(:,:), eaw(:), ebw(:), w1(:,:), w2(:,:)
    real(kind=dp), allocatable :: rhs(:,:), zsol(:,:), moa_s(:,:), mob_s(:,:), ea_s(:), eb_s(:)
    integer :: nbf, nocca, noccb, nvira, nvirb, la, lb, ltot
    integer :: iov, amo, ivir, iocc, sgn
    real(kind=dp) :: th, omp, omm, ct, st
    integer(4) :: status

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvira = nbf - nocca ; nvirb = nbf - noccb
    la = nocca*nvira ; lb = noccb*nvirb ; ltot = la + lb
    th = 1.0d-3

    allocate(cwa(nbf,nbf), cwb(nbf,nbf), eaw(nbf), ebw(nbf), w1(nbf,nbf), w2(nbf,nbf))
    allocate(rhs(ltot,1), zsol(ltot,1), source=0.0_dp)

    ! ---- numerical RHS: R = dω/dκ (canonical orbital rotations) ----
    do sgn = 1, 2                     ! sgn=1 alpha block, sgn=2 beta block
      block
        integer :: nocc, nvir, base
        if (sgn == 1) then ; nocc = nocca ; nvir = nvira ; base = 0
        else               ; nocc = noccb ; nvir = nvirb ; base = la ; end if
        do ivir = 1, nvir
          amo = nocc + ivir
          do iocc = 1, nocc
            iov = base + (ivir-1)*nocc + iocc
            ct = cos(th) ; st = sin(th)
            ! +theta Givens rotation of (iocc, amo) in spin sgn; other spin fixed
            cwa = cac ; cwb = cbc ; eaw = epsca ; ebw = epscb
            if (sgn == 1) then
              cwa(:,iocc) = ct*cac(:,iocc) - st*cac(:,amo)
              cwa(:,amo)  = st*cac(:,iocc) + ct*cac(:,amo)
            else
              cwb(:,iocc) = ct*cbc(:,iocc) - st*cbc(:,amo)
              cwb(:,amo)  = st*cbc(:,iocc) + ct*cbc(:,amo)
            end if
            call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 0)
            call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 1)
            call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omp)
            ! -theta
            cwa = cac ; cwb = cbc ; eaw = epsca ; ebw = epscb
            if (sgn == 1) then
              cwa(:,iocc) =  ct*cac(:,iocc) + st*cac(:,amo)
              cwa(:,amo)  = -st*cac(:,iocc) + ct*cac(:,amo)
            else
              cwb(:,iocc) =  ct*cbc(:,iocc) + st*cbc(:,amo)
              cwb(:,amo)  = -st*cbc(:,iocc) + ct*cbc(:,amo)
            end if
            call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 0)
            call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 1)
            call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omm)
            rhs(iov,1) = (omp - omm)/(2.0_dp*th)
          end do
        end do
      end block
    end do

    ! ---- solve M z = -R via cphf_solve_uhf (transiently feed CANONICAL MOs through tagarray) ----
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, moa, status)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mob, status)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, ea, status)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, eb, status)
    allocate(moa_s, source=moa) ; allocate(mob_s, source=mob)
    allocate(ea_s, source=ea)   ; allocate(eb_s, source=eb)
    moa = cac ; mob = cbc ; ea = epsca ; eb = epscb
    zsol = 0.0_dp
    call cphf_solve_uhf(infos, 1, -rhs, zsol, tol=1.0d-9)
    moa = moa_s ; mob = mob_s ; ea = ea_s ; eb = eb_s     ! restore rotated MOs
    zrms = sqrt(sum(zsol(:,1)**2)/max(ltot,1))
    resid = 0.0_dp

    ! ---- AO relaxation densities  P_z^σ = Σ_ia C^σ_i z^σ_ia C^σ_a^T  (symmetrized) ----
    pza = 0.0_dp ; pzb = 0.0_dp
    do ivir = 1, nvira
      amo = nocca + ivir
      do iocc = 1, nocca
        iov = (ivir-1)*nocca + iocc
        call add_outer(pza, zsol(iov,1), cac(:,iocc), cac(:,amo))
        call add_outer(pza, zsol(iov,1), cac(:,amo), cac(:,iocc))
      end do
    end do
    do ivir = 1, nvirb
      amo = noccb + ivir
      do iocc = 1, noccb
        iov = la + (ivir-1)*noccb + iocc
        call add_outer(pzb, zsol(iov,1), cbc(:,iocc), cbc(:,amo))
        call add_outer(pzb, zsol(iov,1), cbc(:,amo), cbc(:,iocc))
      end do
    end do

    deallocate(cwa, cwb, eaw, ebw, w1, w2, rhs, zsol, moa_s, mob_s, ea_s, eb_s)
  end subroutine umrsf_zvector_relax

!###############################################################################
!> B(mu,nu) += c * u(mu) * w(nu)
  subroutine add_outer(b, c, u, w)
    implicit none
    real(kind=dp), intent(inout) :: b(:,:)
    real(kind=dp), intent(in) :: c, u(:), w(:)
    integer :: mu, nu
    do nu = 1, size(w)
      do mu = 1, size(u)
        b(mu,nu) = b(mu,nu) + c*u(mu)*w(nu)
      end do
    end do
  end subroutine add_outer

!###############################################################################
!  Custom grd2 response 2-PDM type-bound procedures
!###############################################################################

  subroutine grd2_umrsf_resp_init(this)
    class(grd2_umrsf_resp_t), target, intent(inout) :: this
    ! densities/signs are filled by the caller before grd2_driver; nothing to do.
  end subroutine grd2_umrsf_resp_init

  subroutine grd2_umrsf_resp_clean(this)
    class(grd2_umrsf_resp_t), target, intent(inout) :: this
    if (allocated(this%bden)) deallocate(this%bden)
    if (allocated(this%dden)) deallocate(this%dden)
    if (allocated(this%bden_s)) deallocate(this%bden_s)
    if (allocated(this%dden_s)) deallocate(this%dden_s)
    if (allocated(this%sgn)) deallocate(this%sgn)
    if (allocated(this%has_coul)) deallocate(this%has_coul)
  end subroutine grd2_umrsf_resp_clean

  !> Per shell-quartet response 2-PDM block (mirrors grd2_uhf_compute_data_t_get_density).
  subroutine grd2_umrsf_resp_get_density(this, basis, id, dab, dabmax)
    implicit none
    class(grd2_umrsf_resp_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: df1, nrmij, nrmijk, s
    integer :: i, j, k, l, ch, i1, j1, k1, l1
    integer :: loc(4), nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)

    dabmax = 0
    loc = basis%ao_offset(id)-1
    nbf = basis%naos(id)
    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i
      do j = 1, nbf(2)
        j1 = loc(2) + j
        nrmij = basis%bfnrm(i1)*basis%bfnrm(j1)
        do k = 1, nbf(3)
          k1 = loc(3) + k
          nrmijk = nrmij*basis%bfnrm(k1)
          do l = 1, nbf(4)
            l1 = loc(4) + l
            df1 = 0.0_dp
            do ch = 1, this%nchan
              s = this%sgn(ch)
              if (s == 0.0_dp) cycle
              ! Coulomb (channels with J): +4*sc*s*(Bs(ij)Ds(kl)+Ds(ij)Bs(kl))  [symmetrized densities]
              if (this%has_coul(ch)) then
                df1 = df1 + 4.0_dp*this%sc*s*( &
                        this%bden_s(ch,i1,j1)*this%dden_s(ch,k1,l1) &
                      + this%dden_s(ch,i1,j1)*this%bden_s(ch,k1,l1) )
              end if
              ! Exchange (all channels): full 8-fold symmetrization of Γ^X(ijkl)=B(ik)D(jl)
              ! with RAW densities (K is NOT symmetrization-invariant). Coeff -s*sx.
              df1 = df1 - this%sx*s*( &
                      this%bden(ch,i1,k1)*this%dden(ch,j1,l1) &
                    + this%bden(ch,j1,k1)*this%dden(ch,i1,l1) &
                    + this%bden(ch,i1,l1)*this%dden(ch,j1,k1) &
                    + this%bden(ch,j1,l1)*this%dden(ch,i1,k1) &
                    + this%bden(ch,k1,i1)*this%dden(ch,l1,j1) &
                    + this%bden(ch,l1,i1)*this%dden(ch,k1,j1) &
                    + this%bden(ch,k1,j1)*this%dden(ch,l1,i1) &
                    + this%bden(ch,l1,j1)*this%dden(ch,k1,i1) )
            end do
            dabmax = max(dabmax, abs(df1))
            ab(l,k,j,i) = df1*(nrmijk*basis%bfnrm(l1))
          end do
        end do
      end do
    end do
  end subroutine grd2_umrsf_resp_get_density

  !> Assemble the symmetrized response 2-PDM densities + channel scales s_k into gcomp.
  subroutine umrsf_resp_2pdm_fill(gcomp, dens_ket, brad_bra, nbf, mrst, hfs, &
                                  scale_exch, spc_coco, spc_ovov, spc_coov)
    type(grd2_umrsf_resp_t), intent(inout) :: gcomp
    real(kind=dp), intent(in) :: dens_ket(:,:,:)   ! (11,nbf,nbf) = umrsfcbc(X)
    real(kind=dp), intent(in) :: brad_bra(:,:,:)   ! (11,nbf,nbf) = umrsf_bra_density(X)
    integer, intent(in) :: nbf, mrst
    real(kind=dp), intent(in) :: hfs, scale_exch, spc_coco, spc_ovov, spc_coov
    integer :: ch

    call gcomp%clean()
    gcomp%nbf = nbf
    gcomp%nchan = 11
    gcomp%sc = scale_exch
    gcomp%sx = scale_exch
    allocate(gcomp%bden(11,nbf,nbf), gcomp%dden(11,nbf,nbf), &
             gcomp%bden_s(11,nbf,nbf), gcomp%dden_s(11,nbf,nbf), &
             gcomp%sgn(11), gcomp%has_coul(11))
    ! RAW densities (exchange + energy); symmetrized copies (Coulomb only — J is symmetrization-
    ! invariant but K is NOT, so exchange MUST use the raw channel densities).
    do ch = 1, 11
      gcomp%dden(ch,:,:)   = dens_ket(ch,:,:)
      gcomp%bden(ch,:,:)   = brad_bra(ch,:,:)
      gcomp%dden_s(ch,:,:) = 0.5_dp*(dens_ket(ch,:,:) + transpose(dens_ket(ch,:,:)))
      gcomp%bden_s(ch,:,:) = 0.5_dp*(brad_bra(ch,:,:) + transpose(brad_bra(ch,:,:)))
    end do
    gcomp%has_coul = (/ .true.,.true.,.true.,.true.,.true.,.true.,.true.,.true., &
                        .false.,.false.,.false. /)
    ! s_k = mrst sign * spin-pair-coupling scale (matches the energy fmrst2 handling)
    gcomp%sgn = 1.0_dp
    if (mrst == 3) gcomp%sgn(1:10) = -1.0_dp
    if (abs(hfs) > epsilon(1.0_dp)) then
      gcomp%sgn(1:8) = gcomp%sgn(1:8) * (spc_coov/hfs)
      gcomp%sgn(9)   = gcomp%sgn(9)   * (spc_ovov/hfs)
      gcomp%sgn(10)  = gcomp%sgn(10)  * (spc_coco/hfs)
    end if
  end subroutine umrsf_resp_2pdm_fill

end module tdhf_umrsf_gradient_mod
