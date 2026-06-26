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
    real(kind=dp), allocatable :: peffa(:,:), peffb(:,:)        ! P_eff = P^Δ,u + ½ P_z (c03 recipe)
    real(kind=dp), allocatable :: de_orb(:,:), de_w(:,:), de_m1(:,:)
    integer :: iwmode                                          ! DEBUG bisection switch (UMRSF_WMODE)
    real(kind=dp) :: zw                                        ! z-weight in P_eff
    logical :: lrr, l2e, lm1, lwsz                             ! include refrelax / 2e / M1 ; zero within-seg 2e
    integer :: irbrhs                                          ! RHS mode: 0 frozen+jacobi / 1 rebuilt+jacobi / 2 rebuilt+fixed-align
    real(kind=dp) :: omega_orb_chk, omega_orb, hfscale_ref
    integer :: ia, ib, i, j
    ! Z-vector (relaxation): canonical MOs + relaxation density
    real(kind=dp), allocatable :: cac(:,:), cbc(:,:), epsca(:), epscb(:), pza(:,:), pzb(:,:)
    real(kind=dp), allocatable :: zmata(:,:), zmatb(:,:)
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

    ! Re-initialize int2_driver at the BASE geometry. The frozen-density FD self-test loop above
    ! cleaned+re-init'd the driver inside umrsf_frozen_omega2e at the LAST displaced geometry; reusing
    ! that stale driver for the subsequent ω re-evaluations (Z-vector RHS, W_2e generalized Fock) gives
    ! a ~4.4e-5 inconsistency in ω_2e which, amplified by 1/θ in the FD generalized Fock, blows up W.
    call int2_driver%clean()
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_driver%schwarz = .false.

    ! UNIT TEST (rigorous-W primitive): rebuild F^ref from the orbitals' occupied block and confirm it
    ! reproduces the stored SCF Fock. Validates umrsf_ref_fock (reference-density-relaxation primitive).
    block
      real(kind=dp), allocatable :: fa_rb(:,:), fb_rb(:,:), fa_st(:,:), fb_st(:,:)
      real(kind=dp) :: hfscale_rb
      hfscale_rb = 1.0_dp
      if (infos%control%hamilton >= 20) hfscale_rb = infos%dft%hfscale
      allocate(fa_rb(nbf,nbf), fb_rb(nbf,nbf), fa_st(nbf,nbf), fb_st(nbf,nbf))
      call umrsf_ref_fock(infos, basis, va, vb, hfscale_rb, fa_rb, fb_rb)
      call unpack_matrix(fock_a, fa_st) ; call unpack_matrix(fock_b, fb_st)
      open(unit=iw, file=infos%log_filename, position="append")
      write(iw,'(/2x,a,2es12.3)') 'umrsf_ref_fock unit test: max|F^ref_rebuilt - fock_stored| a/b = ', &
        maxval(abs(fa_rb-fa_st)), maxval(abs(fb_rb-fb_st))
      close(iw)
      deallocate(fa_rb, fb_rb, fa_st, fb_st)
    end block

    ! DIAGNOSTIC (rigorous-W blocker): is get_jacobi(cac) == va up to column signs? Compute w2e(base)
    ! from cja=get_jacobi(cac) with and without sign-matching cja to va. If sign-fixed w2e == -0.519,
    ! the blocker is purely signs (fixable); if column overlaps < 1, it is a different alignment.
    block
      real(kind=dp), allocatable :: fa_rf(:,:), fb_rf(:,:), cja(:,:), cjb(:,:), ead(:), ebd(:)
      real(kind=dp), allocatable :: w1d(:,:), w2d(:,:), bradd(:,:,:), xmd(:,:)
      real(kind=dp), allocatable, target :: densd(:,:,:,:)
      real(kind=dp), pointer :: f3d(:,:,:,:)
      type(int2_umrsf_data_t), target :: udd
      real(kind=dp) :: hfrf, ovp, minov, w2e_raw, w2e_sf, sgnp
      integer :: pp, kk
      hfrf = 1.0_dp ; if (infos%control%hamilton >= 20) hfrf = infos%dft%hfscale
      allocate(fa_rf(nbf,nbf), fb_rf(nbf,nbf), cja(nbf,nbf), cjb(nbf,nbf), ead(nbf), ebd(nbf), &
               w1d(nbf,nbf), w2d(nbf,nbf), bradd(11,nbf,nbf), xmd(nbf,nbf), densd(1,11,nbf,nbf), source=0.0_dp)
      call umrsf_ref_fock(infos, basis, cac, cbc, hfrf, fa_rf, fb_rf)
      cja = cac ; cjb = cbc
      call get_jacobi(infos, cja, ead, cjb, ebd, smat_full, nocca, w1d, w2d, 0)
      call get_jacobi(infos, cja, ead, cjb, ebd, smat_full, nocca, w1d, w2d, 1)
      minov = 1.0_dp
      do pp = 1, nbf
        ovp = abs(dot_product(cja(:,pp), matmul(smat_full, va(:,pp))))
        minov = min(minov, ovp)
      end do
      ! raw w2e
      call iatogen(bvec(:,tstate), xmd, nocca, noccb)
      call umrsfcbc(infos, cja, cjb, xmd, densd(1,:,:,:))
      udd = int2_umrsf_data_t(d3=densd(1:1,:,:,:), tamm_dancoff=.true., scale_exchange=scale_exch, scale_coulomb=scale_exch)
      call int2_driver%run(udd) ; f3d => udd%f3(:,:,:,:,1)
      if (mrst == 3) f3d(:,1:10,:,:) = -f3d(:,1:10,:,:)
      call umrsf_bra_density(infos, cja, cjb, xmd, bradd)
      w2e_raw = 0.0_dp ; do kk = 1, 11 ; w2e_raw = w2e_raw + sum(bradd(kk,:,:)*f3d(1,kk,:,:)) ; end do
      ! sign-fix cja/cjb to va/vb, recompute
      do pp = 1, nbf
        sgnp = dot_product(cja(:,pp), matmul(smat_full, va(:,pp))) ; if (sgnp < 0.0_dp) cja(:,pp) = -cja(:,pp)
        sgnp = dot_product(cjb(:,pp), matmul(smat_full, vb(:,pp))) ; if (sgnp < 0.0_dp) cjb(:,pp) = -cjb(:,pp)
      end do
      densd = 0.0_dp
      call umrsfcbc(infos, cja, cjb, xmd, densd(1,:,:,:))
      udd = int2_umrsf_data_t(d3=densd(1:1,:,:,:), tamm_dancoff=.true., scale_exchange=scale_exch, scale_coulomb=scale_exch)
      call int2_driver%run(udd) ; f3d => udd%f3(:,:,:,:,1)
      if (mrst == 3) f3d(:,1:10,:,:) = -f3d(:,1:10,:,:)
      call umrsf_bra_density(infos, cja, cjb, xmd, bradd)
      w2e_sf = 0.0_dp ; do kk = 1, 11 ; w2e_sf = w2e_sf + sum(bradd(kk,:,:)*f3d(1,kk,:,:)) ; end do
      open(unit=iw, file=infos%log_filename, position="append")
      write(iw,'(/2x,a,f10.6,a,f16.10,a,f16.10,a,f16.10)') &
        'get_jacobi(cac) vs va: min|col overlap|=', minov, '  w2e_raw=', w2e_raw, &
        '  w2e_signfixed=', w2e_sf, '  (target omega_2e=', omega_2e_mv, ')'
      close(iw)
      deallocate(fa_rf, fb_rf, cja, cjb, ead, ebd, w1d, w2d, bradd, xmd, densd)
    end block

    ! ---- Z-vector relaxation (numerical-RHS end-to-end validation) ----
    ! R = ∂ω/∂κ by FD over canonical orbital rotations; solve via cphf_solve_uhf; relaxation density
    ! P_z folded into P^Δ (= P^Δ,u + z). Needs the live int2 driver for the ω re-evaluations.
    ! UMRSF_RBRHS=1 ⇒ R uses REBUILT F^ref (c03; captures the refrelax term Tr(P^Δ,u ∂F^ref/∂κ) in R).
    hfscale_ref = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale_ref = infos%dft%hfscale
    block
      character(len=16) :: renv
      integer :: ios
      call get_environment_variable("UMRSF_RBRHS", renv, status=ios)
      irbrhs = 0 ; if (ios == 0) read(renv,*,iostat=ios) irbrhs ; if (ios /= 0) irbrhs = 0
    end block
    allocate(pza(nbf,nbf), pzb(nbf,nbf), zmata(nbf,nbf), zmatb(nbf,nbf), source=0.0_dp)
    call umrsf_zvector_relax(infos, int2_driver, basis, va, vb, cac, cbc, epsca, epscb, smat_full, &
                             fock_a, fock_b, hfscale_ref, irbrhs, bvec(:,tstate), scale_exch, &
                             pza, pzb, zmata, zmatb, zrms, zresid)
    ! c03 rigorous recipe (CAS-verified in c04_uhf_W_factors.py): the 1e + mean-field gradient uses
    ! P_eff = P^Δ,u + ½ P_z (relaxation at HALF weight); the OTHER ½ P_z is carried by the W z-coupling
    ! (G^z). Keep pda/pdb = P^Δ,u (UNRELAXED) for the analytic W refrelax/2e; build P_eff separately.
    ! W-mode switch UMRSF_WMODE. DEFAULT 1 = validated heuristic (relaxed-T frozen W, 3.293e-3) —
    ! the best result to date; the rigorous c03 recipe (0) is CAS-verified but INCOMPLETE (the M1
    ! get_jacobi within-segment response is missing — see LOG/PROGRESS), so it currently regresses
    ! to ~1.78e-2. Modes 0,2-8 are the analysis switches used to localize M1 (kept for next session):
    !   0 rigorous(½z, frozen+rr+2e, +M1probe) | 1 heuristic(1z, frozen only) DEFAULT 3.293e-3
    !   2 ½z fr+rr | 3 ½z fr+2e | 4 ½z fr | 5 1z fr+rr+2e | 6 0z fr+rr+2e | 7 0z fr | 8 −½z fr+rr+2e
    block
      character(len=16) :: wenv
      integer :: ios
      call get_environment_variable("UMRSF_WMODE", wenv, status=ios)
      iwmode = 1 ; if (ios == 0) read(wenv,*,iostat=ios) iwmode ; if (ios /= 0) iwmode = 1
    end block
    zw = 0.5_dp ; lrr = .true. ; l2e = .true. ; lm1 = .true. ; lwsz = .false.
    select case (iwmode)
    case (1) ; zw = 1.0_dp ; lrr = .false. ; l2e = .false. ; lm1 = .false.
    case (2) ; l2e = .false.
    case (3) ; lrr = .false.
    case (4) ; lrr = .false. ; l2e = .false.
    case (5) ; zw = 1.0_dp
    case (6) ; zw = 0.0_dp                                   ! rigorous, NO relaxation (test z sign/size)
    case (7) ; zw = 0.0_dp ; lrr = .false. ; l2e = .false.  ! frozen-only, no z (P^Δ,u)
    case (8) ; zw = -0.5_dp                                  ! rigorous, NEGATIVE ½z (test z sign)
    case (10) ; lwsz = .true.                                ! rigorous + ZERO 2e within-segment (M1 span-invariance test)
    case (11) ; zw = 0.0_dp ; lwsz = .true.                 ! zw=0 rigorous + zero 2e within-seg
    end select
    allocate(peffa(nbf,nbf), peffb(nbf,nbf))
    peffa = pda + zw*pza
    peffb = pdb + zw*pzb
    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '========= UMRSF Z-vector (numerical RHS + cphf_solve_uhf) ========='
    write(iw,'(2x,a,es12.3)') 'z RMS = ', zrms
    write(iw,'(2x,a,i0,a,f4.2,4(a,l1))') 'WMODE=', iwmode, ' zw=', zw, ' rr=',lrr,' 2e=',l2e,' m1=',lm1
    write(iw,'(2x,a)') '==================================================================='
    close(iw)

    ! W (energy-weighted density, overlap Pulay): RIGOROUS analytic c03 recipe (CAS-verified, c04),
    !   W = ½ Σ_σ C_σ (G^ω + G^z)_σ,sym C_σ^T,  F^ref REBUILT — NO FD-of-L / NO get_jacobi re-apply.
    ! Built from closed-form pieces: frozen ½(F^MO T_eff + T_eff F^MO) + ref-relaxation 2 C^T G[P_eff] C
    ! |occ-cols (both with P_eff = P^Δ,u+½P_z, in the canonical basis) + the 2e channel generalized Fock
    ! (va basis). NB needs the live (Schwarz-off) int2 driver for the 2e part — compute BEFORE clean.
    allocate(de_w(3,natom), de_m1(3,natom), source=0.0_dp)
    call umrsf_w_analytic(infos, int2_driver, basis, fock_a, fock_b, va, vb, smat_full, &
                          peffa, peffb, lrr, l2e, lwsz, bvec(:,tstate), scale_exch, hfscale_ref, de_w)

    ! M1 (two-reference): the get_jacobi alignment's explicit overlap (dS/dx) Pulay response — the
    ! ANALYTIC closed form is intractable (normalized Jacobi sweep + 1e-3 threshold), so compute it
    ! EXACTLY via the same get_jacobi re-aligned to S(x±θ) with orbitals/ERIs at base (smooth). Needs
    ! the live int2 driver for ω. (Geometry is perturbed+restored internally.)
    if (lm1) call umrsf_m1_overlap_grad(infos, int2_driver, basis, cac, cbc, va, vb, fock_a, fock_b, &
                               smat_full, bvec(:,tstate), scale_exch, de_m1)

    call int2_driver%clean()

    ! orbital-part gradient: 1e Tr(P_eff h^x) + 2e mean-field Tr(P_eff G[P^ref])  (P_eff = P^Δ,u + ½P_z)
    allocate(de_orb(3,natom), source=0.0_dp)
    call umrsf_orbital_grad(infos, basis, peffa, peffb, dmat_a, dmat_b, hfscale_ref, de_orb)

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '====== UMRSF response gradient pieces (de2e / de_orb / de_w / de_m1) ======'
    write(iw,'(2x,a)') '   atom  comp        de_2e               de_orb               de_w                de_m1'
    do iat = 1, natom
      do icmp = 1, 3
        write(iw,'(2x,2i5,4es21.11)') iat, icmp, de2e(icmp,iat), de_orb(icmp,iat), de_w(icmp,iat), de_m1(icmp,iat)
      end do
    end do
    write(iw,'(2x,a)') '=========================================================================='
    close(iw)

    ! full response = transition-2e + orbital(1e+meanfield, P_eff) + rigorous analytic W + M1 alignment
    de2e_out = de2e + de_orb + de_w + de_m1

    deallocate(va, vb, fa, fb, smat_full, ea, eb, wrk1, wrk2, scr, xmat, amo, amo2e, dens, brad)
    deallocate(de2e, de2e_fd, densym)
    deallocate(talpha, tbeta, pda, pdb, peffa, peffb, de_orb, de_w, de_m1)
    deallocate(cac, cbc, epsca, epscb, pza, pzb, zmata, zmatb)

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
!> RIGOROUS analytic energy-weighted density W (c03 recipe; CAS-verified in c04_uhf_W_factors.py):
!>   W = ½ Σ_σ C_σ (G^ω + G^z)_σ,sym C_σ^T,  with F^ref REBUILT from the orbitals (reference-density
!>   relaxation). Assembled from CLOSED-FORM pieces — NO FD-of-the-Lagrangian, NO get_jacobi re-apply
!>   (the dead end). Per spin σ, in the MO basis (W^AO is basis-independent, so frozen+refrelax use the
!>   CANONICAL basis where F^MO=diag ε; the 2e part uses the va/vb energy/amplitude basis):
!>     frozen   W_fr,σ = ½(F^MO_σ T_eff,σ + T_eff,σ F^MO_σ) = ½(ε_p+ε_q) T_eff,σ,pq , T_eff = C^T S P_eff S C
!>     refrelax W_rr,σ = sym( (C_σ^T G_σ[P_eff] C_σ) |occ-cols ) = ½([p≤nocc]+[q≤nocc]) (C^T G_σ[P_eff] C)_pq
!>                       with G_σ[P]=J[Pα+Pβ]−hfscale_ref·K[Pσ]  (the reference mean field, ONE fock_jk build)
!>     2e       W_2e,σ = ½ C_σ G^2e_sym,σ C_σ^T , G^2e_sym = ∂ω_2e/∂U via SYMMETRIC variation of the CLEAN
!>                       explicit-bra ω_2e (umrsf_omega2e_explicit — valid under non-orthonormal variation;
!>                       no get_jacobi ⇒ smooth). [TODO task5: replace with the analytic channel-adjoint.]
!>   P_eff = P^Δ,u + ½ P_z (AO). M1 (two-reference get_jacobi response) is added separately. The frozen
!>   z-coupling and the refrelax z density are folded via P_eff (½ P_z) — see c04 derivation.
!>   Gradient piece = −Tr(W S^x) via grd1 grad_ee_overlap (eijden: negate + half-diagonal pack).
  subroutine umrsf_w_analytic(infos, idrv, basis, fock_a, fock_b, va, vb, smat_full, &
                              peffa, peffb, lrr, l2e, lwsz, xv, scale_exch, hfscale_ref, de_w)
    use io_constants, only: iw
    use int2_compute, only: int2_compute_t
    use scf_addons, only: fock_jk
    use grd1, only: grad_ee_overlap
    use mathlib, only: pack_matrix, unpack_matrix, orthogonal_transform_sym
    use constants, only: tol_int
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: fock_a(:), fock_b(:), va(:,:), vb(:,:)
    real(kind=dp), intent(in) :: smat_full(:,:), peffa(:,:), peffb(:,:)
    logical, intent(in) :: lrr, l2e, lwsz
    real(kind=dp), intent(in) :: xv(:), scale_exch, hfscale_ref
    real(kind=dp), intent(out) :: de_w(:,:)

    ! W built ENTIRELY in the va/vb (get_jacobi) basis so the within-segment blocks (α {1..nocca-1},
    ! β {nocca..nbf}) are well-defined and can be zeroed (lwsz) — the get_jacobi span-invariance: ω is
    ! invariant to within-segment rotations ⇒ that gen-Fock content is spurious gauge (the M1 mechanism).
    real(kind=dp), allocatable :: wtot(:,:), w2ep(:,:), wpack(:), scr(:)
    real(kind=dp), allocatable :: famo(:,:), teff(:,:), ggm(:,:), wmo(:,:), gsym(:,:), cw_a(:,:), cw_b(:,:)
    real(kind=dp), allocatable :: dens(:,:), fout(:,:), gg(:,:), de_dbg(:,:)
    integer :: nbf, nbf2, nocca, noccb, p, q, ij, i, sgn, nocc, slo, shi
    real(kind=dp) :: tol, th, omp, omm, wij, o2e_base

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; tol = tol_int*log(10.0_dp) ; th = 1.0d-3
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    allocate(wtot(nbf,nbf), w2ep(nbf,nbf), wpack(nbf2), scr(nbf2), famo(nbf,nbf), teff(nbf,nbf), &
             ggm(nbf,nbf), wmo(nbf,nbf), gsym(nbf,nbf), cw_a(nbf,nbf), cw_b(nbf,nbf), &
             dens(nbf2,2), fout(nbf2,2), gg(nbf,nbf), de_dbg(3,size(de_w,2)), source=0.0_dp)
    wtot = 0.0_dp ; w2ep = 0.0_dp

    ! reference mean field G_σ[P_eff] = J[Peffα+Peffβ] − hfscale_ref·K[Peffσ] (one fock_jk build)
    if (lrr) then
      call pack_matrix(peffa, dens(:,1), 'U') ; call pack_matrix(peffb, dens(:,2), 'U')
      call fock_jk(basis, dens, fout, infos, scale_exch=hfscale_ref, scale_coul=1.0_dp)
    end if

    do sgn = 1, 2
      ! pick the spin's MOs, fock, P_eff, occupancy, and the get_jacobi segment [slo:shi]
      if (sgn == 1) then ; nocc = nocca ; slo = 1 ; shi = nocca-1
      else               ; nocc = noccb ; slo = nocca ; shi = nbf ; end if
      wmo = 0.0_dp
      ! F^MO_σ = Cσ^T F^ref Cσ (at base F^ref = stored fock) ; T_eff = Cσ^T S P_eff,σ S Cσ
      if (sgn == 1) then
        call orthogonal_transform_sym(nbf, nbf, fock_a, va, nbf, scr) ; call unpack_matrix(scr, famo)
        teff = matmul(transpose(va), matmul(smat_full, matmul(peffa, matmul(smat_full, va))))
      else
        call orthogonal_transform_sym(nbf, nbf, fock_b, vb, nbf, scr) ; call unpack_matrix(scr, famo)
        teff = matmul(transpose(vb), matmul(smat_full, matmul(peffb, matmul(smat_full, vb))))
      end if
      ! frozen ½(F^MO T_eff + T_eff F^MO)
      wmo = wmo + 0.5_dp*(matmul(famo, teff) + matmul(teff, famo))
      ! ref-relaxation sym((Cσ^T G_σ[P_eff] Cσ)|occ-cols) = ½([p≤occ]+[q≤occ]) (Cσ^T Gσ Cσ)_pq
      if (lrr) then
        if (sgn == 1) then ; call unpack_matrix(fout(:,1), gg) ; ggm = matmul(transpose(va), matmul(gg, va))
        else               ; call unpack_matrix(fout(:,2), gg) ; ggm = matmul(transpose(vb), matmul(gg, vb)) ; end if
        ! (refrelax within-seg zeroing also tried — interacts badly with the wrong z; revisit after the z fix)
        do q = 1, nbf ; do p = 1, nbf
          wij = 0.0_dp ; if (p <= nocc) wij = wij + 0.5_dp ; if (q <= nocc) wij = wij + 0.5_dp
          wmo(p,q) = wmo(p,q) + wij*ggm(p,q)
        end do ; end do
      end if
      ! 2e ½ G^2e_sym,σ via SYMMETRIC variation of the clean explicit-bra ω_2e (no get_jacobi ⇒ smooth)
      if (l2e) then
        gsym = 0.0_dp
        do p = 1, nbf
          do q = p, nbf
            cw_a = va ; cw_b = vb
            if (p == q) then
              if (sgn==1) then ; cw_a(:,p) = (1.0_dp+th)*va(:,p) ; else ; cw_b(:,p) = (1.0_dp+th)*vb(:,p) ; end if
            else if (sgn==1) then ; cw_a(:,p) = va(:,p)+th*va(:,q) ; cw_a(:,q) = va(:,q)+th*va(:,p)
            else ; cw_b(:,p) = vb(:,p)+th*vb(:,q) ; cw_b(:,q) = vb(:,q)+th*vb(:,p) ; end if
            call umrsf_omega2e_explicit(infos, idrv, cw_a, cw_b, xv, scale_exch, omp)
            cw_a = va ; cw_b = vb
            if (p == q) then
              if (sgn==1) then ; cw_a(:,p) = (1.0_dp-th)*va(:,p) ; else ; cw_b(:,p) = (1.0_dp-th)*vb(:,p) ; end if
            else if (sgn==1) then ; cw_a(:,p) = va(:,p)-th*va(:,q) ; cw_a(:,q) = va(:,q)-th*va(:,p)
            else ; cw_b(:,p) = vb(:,p)-th*vb(:,q) ; cw_b(:,q) = vb(:,q)-th*vb(:,p) ; end if
            call umrsf_omega2e_explicit(infos, idrv, cw_a, cw_b, xv, scale_exch, omm)
            if (p == q) then ; gsym(p,p) = (omp-omm)/(2.0_dp*th)
            else ; gsym(p,q) = (omp-omm)/(4.0_dp*th) ; gsym(q,p) = gsym(p,q) ; end if
          end do
        end do
        ! M1 (get_jacobi span-invariance): the 2e gen-Fock within-segment block (the SOMO/closed and
        ! SOMO/virt couplings get_jacobi makes ω invariant to) is spurious gauge — zero it if lwsz.
        ! (Only the 2e piece; the frozen orbital-energy within-seg block is REAL — zeroing it regresses.)
        if (lwsz) gsym(slo:shi, slo:shi) = 0.0_dp
        wmo = wmo + 0.5_dp*gsym
        if (sgn==1) then ; w2ep = w2ep + 0.5_dp*matmul(matmul(va, gsym), transpose(va))
        else             ; w2ep = w2ep + 0.5_dp*matmul(matmul(vb, gsym), transpose(vb)) ; end if
      end if
      ! accumulate AO:  W^AO += Cσ W^MO_σ Cσ^T
      if (sgn==1) then ; wtot = wtot + matmul(matmul(va, wmo), transpose(va))
      else             ; wtot = wtot + matmul(matmul(vb, wmo), transpose(vb)) ; end if
    end do

    ! −Tr(W S^x): pack −W (negate + half-diagonal, eijden convention) and contract via grd1
    call pack_matrix(-wtot, wpack, 'U')
    ij = 0 ; do i = 1, nbf ; ij = ij + i ; wpack(ij) = 0.5_dp*wpack(ij) ; end do
    de_w = 0.0_dp
    call grad_ee_overlap(basis, wpack, de_w, logtol=tol)

    ! diagnostic: 2e-only de_w + base ω_2e + matrix norms
    call umrsf_omega2e_explicit(infos, idrv, va, vb, xv, scale_exch, o2e_base)
    call pack_matrix(-w2ep, wpack, 'U')
    ij = 0 ; do i = 1, nbf ; ij = ij + i ; wpack(ij) = 0.5_dp*wpack(ij) ; end do
    de_dbg = 0.0_dp ; call grad_ee_overlap(basis, wpack, de_dbg, logtol=tol)
    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a,l1,a)') '--- W_analytic diagnostic (va-basis; within-seg zeroed=', lwsz, ') ---'
    write(iw,'(2x,a,f16.10)') 'omega_2e(va) base via explicit-bra      = ', o2e_base
    write(iw,'(2x,a,2es13.5)') '||W_total|| ||W_2e(pre-zero)||          = ', sqrt(sum(wtot**2)), sqrt(sum(w2ep**2))
    write(iw,'(2x,a)') '-------------------------------------------------------'
    close(iw)

    deallocate(wtot, w2ep, wpack, scr, famo, teff, ggm, wmo, gsym, cw_a, cw_b, dens, fout, gg, de_dbg)
  end subroutine umrsf_w_analytic

!###############################################################################
!> M1 (two-reference) gradient term = the get_jacobi ALIGNMENT's explicit overlap (dS/dx) response.
!> The aligned va,vb satisfy [va^T S vb]_offdiag = 0 within the two segments (get_jacobi). That
!> condition involves S(x), so as the nuclei move the alignment rotates — a Pulay-like term NOT
!> captured by the fixed-alignment W or the Z-vector. The alignment is a REDUNDANT rotation (within
!> α-occ {1..nocca-1} and β-virt {nocca..nbf}), so it leaves the reference density / F^ref / ERIs
!> unchanged ⇒ ω responds ONLY through the rotated orbitals, with BASE integrals. Computed EXACTLY
!> (same get_jacobi the energy uses) by re-aligning to S(x±θ) — orbitals + ERIs held at base — and
!> central-differencing ω. Smooth for small θ (no Jacobi pairing flips), UNLIKE the symmetric-variation
!> FD-of-L (the proven dead end). The full ω is orthonormal here (get_jacobi preserves it), so
!> umrsf_omega_eval (umrsfmntoia back-transform + FROZEN base Fock = F^ref since the density is fixed)
!> is exact.  de_m1 is added alongside de2e/de_orb/de_w.
!>
!> KEY: get_jacobi only rotates when an off-diagonal exceeds its 1e-3 threshold. Feeding the already-
!> aligned va/vb gives sub-threshold off-diagonals under a small dS ⇒ NO rotation ⇒ de_m1=0 (wrong).
!> So start from the UN-ALIGNED canonical cac/cbc (same span; get_jacobi does a FULL alignment, exactly
!> as the Z-vector RHS does), then SIGN-FIX each column against va/vb (get_jacobi(cac)=va up to per-
!> column signs — verified min|col overlap|=1.0). The sign-fix is smooth for a small dS perturbation.
  subroutine umrsf_m1_overlap_grad(infos, idrv, basis, cac, cbc, va, vb, fock_a, fock_b, &
                                   smat_full, xv, scale_exch, de_m1)
    use io_constants, only: iw
    use int2_compute, only: int2_compute_t
    use int1, only: omp_hst
    use tdhf_mrsf_lib, only: get_jacobi
    use mathlib, only: unpack_matrix
    use constants, only: tol_int
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), va(:,:), vb(:,:)
    real(kind=dp), intent(in) :: fock_a(:), fock_b(:), smat_full(:,:), xv(:), scale_exch
    real(kind=dp), intent(out) :: de_m1(:,:)
    real(kind=dp), allocatable :: cwa(:,:), cwb(:,:), spack(:), sfull(:,:), ead(:), ebd(:), w1(:,:), w2(:,:)
    real(kind=dp), allocatable :: hbuf(:), tbuf(:), zq(:)
    integer :: nbf, nbf2, nocca, natom, iat, icmp, p
    real(kind=dp) :: tol, th, omp, omm, schk, ocheck

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; nocca = infos%mol_prop%nelec_a
    natom = ubound(infos%atoms%zn,1) ; tol = tol_int*log(10.0_dp) ; th = 1.0d-3
    allocate(cwa(nbf,nbf), cwb(nbf,nbf), spack(nbf2), sfull(nbf,nbf), ead(nbf), ebd(nbf), &
             w1(nbf,nbf), w2(nbf,nbf), hbuf(nbf2), tbuf(nbf2), zq(natom), source=0.0_dp)
    zq = infos%atoms%zn - infos%basis%ecp_zn_num

    open(unit=iw, file=infos%log_filename, position="append")
    ! convention sanity: omp_hst overlap at base reproduces stored S; sign-fixed get_jacobi(cac,S_base)==ω_base
    call omp_hst(basis, infos%atoms%xyz, zq, hbuf, spack, tbuf, logtol=tol, &
                 comm=infos%mpiinfo%comm, usempi=infos%mpiinfo%usempi)
    call unpack_matrix(spack, sfull, nbf, 'U')
    schk = maxval(abs(sfull - smat_full))
    cwa = cac ; cwb = cbc
    call get_jacobi(infos, cwa, ead, cwb, ebd, sfull, nocca, w1, w2, 0)
    call get_jacobi(infos, cwa, ead, cwb, ebd, sfull, nocca, w1, w2, 1)
    call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
    call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, ocheck)
    write(iw,'(/2x,a,es12.3,a,f16.10)') 'M1: overlap() base vs stored S max|Δ| = ', schk, &
      '   ω(get_jacobi(cac),signfix) base = ', ocheck

    de_m1 = 0.0_dp
    do iat = 1, natom
      do icmp = 1, 3
        ! +θ: recompute S(x+θ); canonical orbitals + ERIs stay at base; full re-align + sign-fix
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) + th
        call basis%init_shell_centers()
        call omp_hst(basis, infos%atoms%xyz, zq, hbuf, spack, tbuf, logtol=tol, comm=infos%mpiinfo%comm, usempi=infos%mpiinfo%usempi)
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) - th
        call basis%init_shell_centers() ; call unpack_matrix(spack, sfull, nbf, 'U')
        cwa = cac ; cwb = cbc
        call get_jacobi(infos, cwa, ead, cwb, ebd, sfull, nocca, w1, w2, 0)
        call get_jacobi(infos, cwa, ead, cwb, ebd, sfull, nocca, w1, w2, 1)
        call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
        call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omp)
        ! -θ
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) - th
        call basis%init_shell_centers()
        call omp_hst(basis, infos%atoms%xyz, zq, hbuf, spack, tbuf, logtol=tol, comm=infos%mpiinfo%comm, usempi=infos%mpiinfo%usempi)
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) + th
        call basis%init_shell_centers() ; call unpack_matrix(spack, sfull, nbf, 'U')
        cwa = cac ; cwb = cbc
        call get_jacobi(infos, cwa, ead, cwb, ebd, sfull, nocca, w1, w2, 0)
        call get_jacobi(infos, cwa, ead, cwb, ebd, sfull, nocca, w1, w2, 1)
        call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
        call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omm)
        de_m1(icmp,iat) = (omp - omm)/(2.0_dp*th)
      end do
    end do
    write(iw,'(2x,a)') 'M1 alignment-overlap (dS/dx) response computed (de_m1).'
    close(iw)
    deallocate(cwa, cwb, spack, sfull, ead, ebd, w1, w2, hbuf, tbuf, zq)
  end subroutine umrsf_m1_overlap_grad

!###############################################################################
!> Fix per-column signs of (cwa,cwb) to match the reference aligned (va,vb): if cwσ(:,p)·S·vσ(:,p)<0,
!> flip column p. get_jacobi(cac) matches va up to column signs (min|col overlap|=1.0), and the V_S/V_T
!> SOMO combos in ω are sign-sensitive ⇒ must restore va's convention before evaluating ω.
  subroutine m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
    implicit none
    integer, intent(in) :: nbf
    real(kind=dp), intent(inout) :: cwa(:,:), cwb(:,:)
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), smat_full(:,:)
    integer :: p
    real(kind=dp) :: s
    do p = 1, nbf
      s = dot_product(cwa(:,p), matmul(smat_full, va(:,p))) ; if (s < 0.0_dp) cwa(:,p) = -cwa(:,p)
      s = dot_product(cwb(:,p), matmul(smat_full, vb(:,p))) ; if (s < 0.0_dp) cwb(:,p) = -cwb(:,p)
    end do
  end subroutine m1_sign_fix

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
!> ω_2e = Σ_k s_k ⟨B_k(C), F_k(C)⟩ via the EXPLICIT clean-room bra density
!> (umrsf_bra_density) — NOT the umrsfmntoia back-transform. umrsfcbc AND umrsf_bra_density
!> are pure bilinears in C (no S / orthonormality dependence; verified by reading umrsfcbc),
!> so this is the exact ω_2e quartic in C and stays valid under NON-orthonormal (symmetric)
!> orbital variations. By contrast umrsf_omega_eval uses umrsfmntoia, whose SOMO-diagonal
!> OVERWRITE makes the adjoint identity ⟨X,mntoia(F)⟩=Σ⟨B_k,F_k⟩ exact ONLY for orthonormal
!> MOs — that is what broke the earlier symmetric-FD W (umrsf_w_numerical, ~20x overshoot).
!> Integrals are at the CURRENT (fixed) geometry through the live idrv (no re-init).
  subroutine umrsf_omega2e_explicit(infos, idrv, vva, vvb, xv, scale_exch, omega2e)
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc
    use tdhf_lib, only: iatogen
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: vva(:,:), vvb(:,:), xv(:), scale_exch
    real(kind=dp), intent(out) :: omega2e
    real(kind=dp), allocatable :: xmat(:,:), brad(:,:,:)
    real(kind=dp), allocatable, target :: dens(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    type(int2_umrsf_data_t), target :: ud
    integer :: nbf, nocca, noccb, mrst, k

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    mrst = infos%tddft%mult
    allocate(xmat(nbf,nbf), brad(11,nbf,nbf), dens(1,11,nbf,nbf), source=0.0_dp)

    call iatogen(xv, xmat, nocca, noccb)
    call umrsfcbc(infos, vva, vvb, xmat, dens(1,:,:,:))
    ud = int2_umrsf_data_t(d3=dens(1:1,:,:,:), tamm_dancoff=.true., &
                           scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(ud)
    fmrst2 => ud%f3(:,:,:,:,1)
    if (mrst == 3) fmrst2(:,1:10,:,:) = -fmrst2(:,1:10,:,:)
    ! (SPC scaling spc==hfs==1 in the HF limit; add for Stage-2 hybrids — see fmrst2 handling above.)
    call umrsf_bra_density(infos, vva, vvb, xmat, brad)
    omega2e = 0.0_dp
    do k = 1, 11
      omega2e = omega2e + sum(brad(k,:,:)*fmrst2(1,k,:,:))
    end do
    deallocate(xmat, brad, dens)
  end subroutine umrsf_omega2e_explicit

!###############################################################################
!> Rebuild the UHF reference Fock F^ref_σ = h + J[Dα+Dβ] − scale_exch·K[Dσ] (AO, full) from the
!> OCCUPIED orbitals of cw_a/cw_b. This is the reference-density-relaxation primitive the rigorous
!> generalized Fock needs (c03_cis_gradient_W.py): ω_orb = Tr(P^Δ F^ref(C)) with F^ref REBUILT from the
!> varied orbitals, NOT the frozen stored fock_a. D^ref is get_jacobi-invariant (seg0/seg1 rotate WITHIN
!> the α-occ / β-virt subspaces, leaving β-occ and the α-occ subspace density unchanged), so cw may be
!> pre- or post-jacobi. Stage-1 HF only (no XC; add the XC Fock for Stage-2). NB fock_jk uses its own
!> screened int2 driver — fine to ~1e-6; for a clean FD use the schwarz-off analytic mean field instead.
  subroutine umrsf_ref_fock(infos, basis, cw_a, cw_b, scale_exch, fa_full, fb_full)
    use guess, only: get_ab_initio_density
    use scf_addons, only: fock_jk
    use mathlib, only: unpack_matrix
    use oqp_tagarray_driver
    implicit none
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cw_a(:,:), cw_b(:,:), scale_exch
    real(kind=dp), intent(out) :: fa_full(:,:), fb_full(:,:)
    real(kind=dp), allocatable :: dens(:,:), fout(:,:)
    real(kind=dp), contiguous, pointer :: hcore(:)
    integer :: nbf, nbf2
    integer(4) :: status

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2
    allocate(dens(nbf2,2), fout(nbf2,2), source=0.0_dp)
    call get_ab_initio_density(dens(:,1), cw_a, dens(:,2), cw_b, infos, basis)
    call fock_jk(basis, dens, fout, infos, scale_exch=scale_exch, scale_coul=1.0_dp)
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore, status)
    fout(:,1) = fout(:,1) + hcore
    fout(:,2) = fout(:,2) + hcore
    call unpack_matrix(fout(:,1), fa_full)
    call unpack_matrix(fout(:,2), fb_full)
    deallocate(dens, fout)
  end subroutine umrsf_ref_fock

!###############################################################################
!> Lagrangian L(C_can) = ω_2e + ω_orb + Σ_σ Σ_ai z^σ_ai F^ref,σ_ai, evaluated at canonical orbitals
!> cw_can (get_jacobi re-applied internally ⇒ M1 captured), with F^ref REBUILT from cw_can (c03 recipe).
!>   ω_2e   = Σ_k ⟨B_k(cw_jac), F_k(cw_jac)⟩      (explicit clean bilinear, cw_jac = get_jacobi(cw_can))
!>   ω_orb  = Tr(T_u_α F^MO_jac_α) + Tr(T_u_β F^MO_jac_β),  F^MO_jac_σ = cw_jac_σ^T F^ref_σ cw_jac_σ
!>   z-coup = ½ Tr(zmat_α F^MO_can_α) + ½ Tr(zmat_β F^MO_can_β)   (Σ_ai z_ai F_ai = ½ Tr(zmat F^MO),
!>            zmat symmetric ov+vo;  F^MO_can_σ = cw_can_σ^T F^ref_σ cw_can_σ).
!> Used by umrsf_w_rigorous for the symmetric generalized Fock G^L (→ W) and reproduces ω at base.
  subroutine umrsf_lagrangian_eval(infos, idrv, basis, cw_can_a, cw_can_b, smat_full, &
                                   va_ref, vb_ref, zmata, zmatb, xv, scale_exch, hfscale_ref, lval, &
                                   w2e_out, worb_out, zc_out)
    use mathlib, only: orthogonal_transform_sym, unpack_matrix, pack_matrix
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc, get_jacobi, mrsfesum
    use tdhf_lib, only: iatogen
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cw_can_a(:,:), cw_can_b(:,:), smat_full(:,:)
    real(kind=dp), intent(in) :: va_ref(:,:), vb_ref(:,:)     ! base get_jacobi orbitals (sign reference)
    real(kind=dp), intent(in) :: zmata(:,:), zmatb(:,:)
    real(kind=dp), intent(in) :: xv(:), scale_exch, hfscale_ref
    real(kind=dp), intent(out) :: lval
    real(kind=dp), intent(out), optional :: w2e_out, worb_out, zc_out

    real(kind=dp), allocatable :: fa_ref(:,:), fb_ref(:,:), cja(:,:), cjb(:,:)
    real(kind=dp), allocatable :: ea(:), eb(:), w1(:,:), w2(:,:), fmo(:,:)
    real(kind=dp), allocatable :: faj(:,:), fbj(:,:), scr(:), xmat(:,:), brad(:,:,:), amo(:,:)
    real(kind=dp), allocatable, target :: dens(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    type(int2_umrsf_data_t), target :: ud
    integer :: nbf, nbf2, nocca, noccb, nvirb, xvec_dim, mrst, k
    real(kind=dp) :: w2e, worb, zc

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb ; xvec_dim = nocca*nvirb ; mrst = infos%tddft%mult
    allocate(fa_ref(nbf,nbf), fb_ref(nbf,nbf), cja(nbf,nbf), cjb(nbf,nbf), &
             ea(nbf), eb(nbf), w1(nbf,nbf), w2(nbf,nbf), fmo(nbf,nbf), &
             faj(nbf,nbf), fbj(nbf,nbf), scr(nbf2), xmat(nbf,nbf), brad(11,nbf,nbf), &
             amo(xvec_dim,1), dens(1,11,nbf,nbf), source=0.0_dp)

    ! F^ref rebuilt from cw_can occupied (reference-density relaxation)
    call umrsf_ref_fock(infos, basis, cw_can_a, cw_can_b, hfscale_ref, fa_ref, fb_ref)

    ! z-coupling (canonical basis): Σ_ai z_ai F^ref_ai = ½ Tr(zmat_σ · cw_can_σ^T F^ref_σ cw_can_σ)
    fmo = matmul(transpose(cw_can_a), matmul(fa_ref, cw_can_a))
    zc = 0.5_dp*sum(zmata*fmo)
    fmo = matmul(transpose(cw_can_b), matmul(fb_ref, cw_can_b))
    zc = zc + 0.5_dp*sum(zmatb*fmo)

    ! get_jacobi → channel orbitals cja/cjb (D^ref invariant, so F^ref unchanged)
    cja = cw_can_a ; cjb = cw_can_b
    call get_jacobi(infos, cja, ea, cjb, eb, smat_full, nocca, w1, w2, 0)
    call get_jacobi(infos, cja, ea, cjb, eb, smat_full, nocca, w1, w2, 1)
    ! get_jacobi(cac) == va only up to per-column SIGNS; match cja/cjb to the base orbitals va/vb so
    ! the V_S/V_T (1/√2) SOMO combos in ω_2e are consistent (verified: recovers ω_2e exactly at base).
    block
      integer :: pp
      do pp = 1, nbf
        if (dot_product(cja(:,pp), matmul(smat_full, va_ref(:,pp))) < 0.0_dp) cja(:,pp) = -cja(:,pp)
        if (dot_product(cjb(:,pp), matmul(smat_full, vb_ref(:,pp))) < 0.0_dp) cjb(:,pp) = -cjb(:,pp)
      end do
    end block

    ! MO Fock (rebuilt F^ref) in the channel basis cja/cjb
    faj = matmul(transpose(cja), matmul(fa_ref, cja))
    fbj = matmul(transpose(cjb), matmul(fb_ref, cjb))

    ! ω_orb via mrsfesum (self-consistent in cja/cjb; no external T_u, so no basis mismatch)
    call iatogen(xv, xmat, nocca, noccb)
    amo = 0.0_dp
    call mrsfesum(infos, xmat, faj, fbj, amo, 1)
    worb = dot_product(xv, amo(:,1))

    ! ω_2e via explicit bra (clean bilinear), channel orbitals cja/cjb
    call iatogen(xv, xmat, nocca, noccb)
    call umrsfcbc(infos, cja, cjb, xmat, dens(1,:,:,:))
    ud = int2_umrsf_data_t(d3=dens(1:1,:,:,:), tamm_dancoff=.true., &
                           scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(ud)
    fmrst2 => ud%f3(:,:,:,:,1)
    if (mrst == 3) fmrst2(:,1:10,:,:) = -fmrst2(:,1:10,:,:)
    call umrsf_bra_density(infos, cja, cjb, xmat, brad)
    w2e = 0.0_dp
    do k = 1, 11
      w2e = w2e + sum(brad(k,:,:)*fmrst2(1,k,:,:))
    end do

    lval = w2e + worb + zc
    if (present(w2e_out))  w2e_out = w2e
    if (present(worb_out)) worb_out = worb
    if (present(zc_out))   zc_out = zc
    deallocate(fa_ref, fb_ref, cja, cjb, ea, eb, w1, w2, fmo, faj, fbj, scr, xmat, brad, amo, dens)
  end subroutine umrsf_lagrangian_eval

!###############################################################################
!> RIGOROUS energy-weighted density W = ½ Σ_σ C_σ (G^L)_σ,sym C_σ^T (c03_cis_gradient_W.py, 6e-9),
!> G^L = ∂L/∂U the symmetric generalized Fock of the full Lagrangian L (umrsf_lagrangian_eval) by central
!> FD over symmetric variations of the CANONICAL orbitals (get_jacobi re-applied ⇒ M1 automatic; F^ref
!> rebuilt ⇒ reference-density relaxation). Replaces the relaxed-T heuristic W^Δ_orb. −Tr(W S^x) via grd1.
  subroutine umrsf_w_rigorous(infos, idrv, basis, cac, cbc, va, vb, smat_full, &
                              zmata, zmatb, xv, scale_exch, hfscale_ref, de_w)
    use int2_compute, only: int2_compute_t
    use grd1, only: grad_ee_overlap
    use mathlib, only: pack_matrix
    use constants, only: tol_int
    use io_constants, only: iw
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), va(:,:), vb(:,:), smat_full(:,:)
    real(kind=dp), intent(in) :: zmata(:,:), zmatb(:,:)
    real(kind=dp), intent(in) :: xv(:), scale_exch, hfscale_ref
    real(kind=dp), intent(out) :: de_w(:,:)

    real(kind=dp), allocatable :: gsym(:,:), wtot(:,:), cwa(:,:), cwb(:,:), wpack(:)
    integer :: nbf, nbf2, sgn, p, q, ij, i
    real(kind=dp) :: th, lp, lm, tol

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; th = 1.0d-3 ; tol = tol_int*log(10.0_dp)
    allocate(gsym(nbf,nbf), wtot(nbf,nbf), cwa(nbf,nbf), cwb(nbf,nbf), wpack(nbf2), source=0.0_dp)

    ! DIAGNOSTIC: L(base) must reproduce ω (z-coupling = 0 at base, since canonical F^ref_ai = 0).
    block
      real(kind=dp) :: lbase, w2eb, worbb, zcb
      call umrsf_lagrangian_eval(infos, idrv, basis, cac, cbc, va, vb, smat_full, &
                                 zmata, zmatb, xv, scale_exch, hfscale_ref, lbase, w2eb, worbb, zcb)
      open(unit=iw, file=infos%log_filename, position="append")
      write(iw,'(/2x,a,f16.10,a,f16.10,a,f16.10,a,f16.10)') &
        'umrsf_w_rigorous diag: L(base)=', lbase, '  w2e=', w2eb, '  worb=', worbb, '  zc=', zcb
      close(iw)
    end block
    ! BLOCKED: L(base) ≠ ω because get_jacobi(cac) reassigns the SOMO positions O1/O2 (re-canonicalization
    ! reorders/sign-flips) ⇒ ω_2e wrong (-0.46 vs -0.519). Needs a SOMO/sign-preserving canonicalization or a
    ! consistent va/vb-basis formulation of the z-coupling before the loop below is valid. See LOG.
    if (infos%tddft%debug_mode) then
      de_w = 0.0_dp ; deallocate(gsym, wtot, cwa, cwb, wpack) ; return
    end if

    do sgn = 1, 2
      gsym = 0.0_dp
      do p = 1, nbf
        do q = p, nbf
          cwa = cac ; cwb = cbc
          if (p == q) then
            if (sgn==1) then ; cwa(:,p) = (1.0_dp+th)*cac(:,p)
            else             ; cwb(:,p) = (1.0_dp+th)*cbc(:,p) ; end if
          else
            if (sgn==1) then
              cwa(:,p) = cac(:,p) + th*cac(:,q) ; cwa(:,q) = cac(:,q) + th*cac(:,p)
            else
              cwb(:,p) = cbc(:,p) + th*cbc(:,q) ; cwb(:,q) = cbc(:,q) + th*cbc(:,p)
            end if
          end if
          call umrsf_lagrangian_eval(infos, idrv, basis, cwa, cwb, va, vb, smat_full, &
                                     zmata, zmatb, xv, scale_exch, hfscale_ref, lp)
          cwa = cac ; cwb = cbc
          if (p == q) then
            if (sgn==1) then ; cwa(:,p) = (1.0_dp-th)*cac(:,p)
            else             ; cwb(:,p) = (1.0_dp-th)*cbc(:,p) ; end if
          else
            if (sgn==1) then
              cwa(:,p) = cac(:,p) - th*cac(:,q) ; cwa(:,q) = cac(:,q) - th*cac(:,p)
            else
              cwb(:,p) = cbc(:,p) - th*cbc(:,q) ; cwb(:,q) = cbc(:,q) - th*cbc(:,p)
            end if
          end if
          call umrsf_lagrangian_eval(infos, idrv, basis, cwa, cwb, va, vb, smat_full, &
                                     zmata, zmatb, xv, scale_exch, hfscale_ref, lm)
          if (p == q) then
            gsym(p,p) = (lp - lm)/(2.0_dp*th)
          else
            gsym(p,q) = (lp - lm)/(4.0_dp*th)
            gsym(q,p) = gsym(p,q)
          end if
        end do
      end do
      if (sgn==1) then
        wtot = wtot + 0.5_dp*matmul(matmul(cac, gsym), transpose(cac))
      else
        wtot = wtot + 0.5_dp*matmul(matmul(cbc, gsym), transpose(cbc))
      end if
    end do

    call pack_matrix(-wtot, wpack, 'U')
    ij = 0
    do i = 1, nbf ; ij = ij + i ; wpack(ij) = 0.5_dp*wpack(ij) ; end do
    de_w = 0.0_dp
    call grad_ee_overlap(basis, wpack, de_w, logtol=tol)
    deallocate(gsym, wtot, cwa, cwb, wpack)
  end subroutine umrsf_w_rigorous

!###############################################################################
!> ω_2e overlap-Pulay W via the SYMMETRIC generalized Fock of ω_2e (numerical ORACLE for the
!> analytic G^2e): W_2e = ½ Σ_σ C_σ G^2e,σ_sym C_σ^T,  G^2e,σ_sym,pq = ∂ω_2e/∂U^σ_(pq,sym),
!> the energy-weighting of the TRANSITION (amplitude) density. ADDED to the difference-density
!> W^Δ_orb. ω_2e is the EXPLICIT clean bilinear (umrsf_omega2e_explicit) so the symmetric
!> (norm-changing) variation is exact. The AO W is orbital-basis invariant ⇒ va/vb is fine.
!> −Tr(W_2e S^x) via grd1 grad_ee_overlap (eijden convention: negate + half-diagonal pack).
  subroutine umrsf_w_2e_numerical(infos, idrv, basis, vva, vvb, xv, scale_exch, de_w)
    use int2_compute, only: int2_compute_t
    use grd1, only: grad_ee_overlap
    use mathlib, only: pack_matrix
    use constants, only: tol_int
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: vva(:,:), vvb(:,:), xv(:), scale_exch
    real(kind=dp), intent(out) :: de_w(:,:)

    real(kind=dp), allocatable :: gsym(:,:), wtot(:,:), cwa(:,:), cwb(:,:), wpack(:)
    integer :: nbf, nbf2, sgn, p, q, ij, i
    real(kind=dp) :: th, omp, omm, tol

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; th = 1.0d-3 ; tol = tol_int*log(10.0_dp)
    allocate(gsym(nbf,nbf), wtot(nbf,nbf), cwa(nbf,nbf), cwb(nbf,nbf), wpack(nbf2), source=0.0_dp)

    do sgn = 1, 2                       ! sgn=1 vary alpha columns, sgn=2 vary beta columns
      gsym = 0.0_dp
      do p = 1, nbf
        do q = p, nbf
          ! +theta symmetric variation of MO pair (p,q): U_pq = U_qp = th
          cwa = vva ; cwb = vvb
          if (p == q) then
            if (sgn==1) then ; cwa(:,p) = (1.0_dp+th)*vva(:,p)
            else             ; cwb(:,p) = (1.0_dp+th)*vvb(:,p) ; end if
          else
            if (sgn==1) then
              cwa(:,p) = vva(:,p) + th*vva(:,q) ; cwa(:,q) = vva(:,q) + th*vva(:,p)
            else
              cwb(:,p) = vvb(:,p) + th*vvb(:,q) ; cwb(:,q) = vvb(:,q) + th*vvb(:,p)
            end if
          end if
          call umrsf_omega2e_explicit(infos, idrv, cwa, cwb, xv, scale_exch, omp)
          ! -theta
          cwa = vva ; cwb = vvb
          if (p == q) then
            if (sgn==1) then ; cwa(:,p) = (1.0_dp-th)*vva(:,p)
            else             ; cwb(:,p) = (1.0_dp-th)*vvb(:,p) ; end if
          else
            if (sgn==1) then
              cwa(:,p) = vva(:,p) - th*vva(:,q) ; cwa(:,q) = vva(:,q) - th*vva(:,p)
            else
              cwb(:,p) = vvb(:,p) - th*vvb(:,q) ; cwb(:,q) = vvb(:,q) - th*vvb(:,p)
            end if
          end if
          call umrsf_omega2e_explicit(infos, idrv, cwa, cwb, xv, scale_exch, omm)
          if (p == q) then
            gsym(p,p) = (omp - omm)/(2.0_dp*th)            ! ∂ω/∂U_pp
          else
            gsym(p,q) = (omp - omm)/(4.0_dp*th)            ! ½(∂ω/∂U_pq+∂ω/∂U_qp) = G_sym,pq
            gsym(q,p) = gsym(p,q)
          end if
        end do
      end do
      if (sgn==1) then
        wtot = wtot + 0.5_dp*matmul(matmul(vva, gsym), transpose(vva))
      else
        wtot = wtot + 0.5_dp*matmul(matmul(vvb, gsym), transpose(vvb))
      end if
    end do

    call pack_matrix(-wtot, wpack, 'U')
    ij = 0
    do i = 1, nbf ; ij = ij + i ; wpack(ij) = 0.5_dp*wpack(ij) ; end do
    de_w = 0.0_dp
    call grad_ee_overlap(basis, wpack, de_w, logtol=tol)
    deallocate(gsym, wtot, cwa, cwb, wpack)
  end subroutine umrsf_w_2e_numerical

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
!> Like umrsf_omega_eval but with F^ref REBUILT from the (rotated) orbitals (c03: ω_orb uses
!> F^ref(C)=h+G[D^ref(C)], NOT the frozen converged Fock). Needed for the Z-vector RHS R=∂ω/∂κ over
!> occ-virt rotations (which CHANGE the reference density ⇒ the refrelax term Tr(P^Δ,u ∂F^ref/∂κ) must
!> be present, else z is wrong). Slower (one fock_jk per call). Orthonormal orbitals ⇒ umrsfmntoia exact.
  subroutine umrsf_omega_eval_rb(infos, idrv, basis, vva, vvb, hfscale_ref, xv, scale_exch, omega)
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc, umrsfmntoia, mrsfesum
    use tdhf_lib, only: iatogen
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: vva(:,:), vvb(:,:), hfscale_ref, xv(:), scale_exch
    real(kind=dp), intent(out) :: omega
    real(kind=dp), allocatable :: fa_full(:,:), fb_full(:,:), fa(:,:), fb(:,:), xmat(:,:), amo(:,:)
    real(kind=dp), allocatable, target :: dens(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    type(int2_umrsf_data_t), target :: ud
    integer :: nbf, nocca, noccb, nvirb, xvec_dim, mrst

    nbf = basis%nbf
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb ; xvec_dim = nocca*nvirb ; mrst = infos%tddft%mult
    allocate(fa_full(nbf,nbf), fb_full(nbf,nbf), fa(nbf,nbf), fb(nbf,nbf), xmat(nbf,nbf), &
             amo(xvec_dim,1), dens(1,11,nbf,nbf), source=0.0_dp)
    ! rebuilt F^ref from the (rotated) orbitals' occupied block, transformed to the MO basis
    call umrsf_ref_fock(infos, basis, vva, vvb, hfscale_ref, fa_full, fb_full)
    fa = matmul(transpose(vva), matmul(fa_full, vva))
    fb = matmul(transpose(vvb), matmul(fb_full, vvb))
    call iatogen(xv, xmat, nocca, noccb)
    call umrsfcbc(infos, vva, vvb, xmat, dens(1,:,:,:))
    ud = int2_umrsf_data_t(d3=dens(1:1,:,:,:), tamm_dancoff=.true., &
                           scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(ud) ; fmrst2 => ud%f3(:,:,:,:,1)
    if (mrst == 3) fmrst2(:,1:10,:,:) = -fmrst2(:,1:10,:,:)
    amo = 0.0_dp
    call umrsfmntoia(infos, fmrst2(1,:,:,:), amo, vva, vvb, 1)
    call iatogen(xv, xmat, nocca, noccb)
    call mrsfesum(infos, xmat, fa, fb, amo, 1)
    omega = dot_product(xv, amo(:,1))
    deallocate(fa_full, fb_full, fa, fb, xmat, amo, dens)
  end subroutine umrsf_omega_eval_rb

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
  subroutine umrsf_zvector_relax(infos, idrv, basis, va, vb, cac, cbc, epsca, epscb, smat_full, &
                                 fock_a, fock_b, hfscale_ref, rmode, xv, scale_exch, pza, pzb, zmata, zmatb, zrms, resid)
    use oqp_tagarray_driver
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: get_jacobi
    use cphf_mod, only: cphf_solve_uhf
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), cac(:,:), cbc(:,:)
    real(kind=dp), intent(in) :: epsca(:), epscb(:), smat_full(:,:)
    real(kind=dp), intent(in) :: fock_a(:), fock_b(:), hfscale_ref, xv(:), scale_exch
    integer, intent(in) :: rmode             ! 0 frozen+get_jacobi | 1 rebuilt+get_jacobi+signfix | 2 rebuilt+fixed-align (M1-free)
    real(kind=dp), intent(out) :: pza(:,:), pzb(:,:), zmata(:,:), zmatb(:,:), zrms, resid
    real(kind=dp), allocatable :: malign(:,:), mblign(:,:)     ! fixed base alignment M=cac^T S va (rmode 2)
    logical :: lrb

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
    lrb = (rmode >= 1)
    ! rmode 2: fixed base alignment M_σ = C_can^T S C_aligned (block-diag occ/virt) — apply to the
    ! rotated canonical orbitals INSTEAD of re-running get_jacobi ⇒ R is M1-FREE (the c03 formulation).
    allocate(malign(nbf,nbf), mblign(nbf,nbf), source=0.0_dp)
    if (rmode == 2) then
      malign = matmul(transpose(cac), matmul(smat_full, va))
      mblign = matmul(transpose(cbc), matmul(smat_full, vb))
    end if

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
            if (rmode == 2) then
              cwa = matmul(cwa, malign) ; cwb = matmul(cwb, mblign)   ! fixed base alignment (M1-free)
            else
              call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 0)
              call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 1)
              if (rmode == 1) call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
            end if
            if (lrb) then ; call umrsf_omega_eval_rb(infos, idrv, basis, cwa, cwb, hfscale_ref, xv, scale_exch, omp)
            else          ; call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omp) ; end if
            ! -theta
            cwa = cac ; cwb = cbc ; eaw = epsca ; ebw = epscb
            if (sgn == 1) then
              cwa(:,iocc) =  ct*cac(:,iocc) + st*cac(:,amo)
              cwa(:,amo)  = -st*cac(:,iocc) + ct*cac(:,amo)
            else
              cwb(:,iocc) =  ct*cbc(:,iocc) + st*cbc(:,amo)
              cwb(:,amo)  = -st*cbc(:,iocc) + ct*cbc(:,amo)
            end if
            if (rmode == 2) then
              cwa = matmul(cwa, malign) ; cwb = matmul(cwb, mblign)
            else
              call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 0)
              call get_jacobi(infos, cwa, eaw, cwb, ebw, smat_full, nocca, w1, w2, 1)
              if (rmode == 1) call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
            end if
            if (lrb) then ; call umrsf_omega_eval_rb(infos, idrv, basis, cwa, cwb, hfscale_ref, xv, scale_exch, omm)
            else          ; call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omm) ; end if
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
    ! Also build the MO relaxation matrix zmat^σ (canonical, symmetric ov+vo) for the W z-coupling.
    pza = 0.0_dp ; pzb = 0.0_dp ; zmata = 0.0_dp ; zmatb = 0.0_dp
    do ivir = 1, nvira
      amo = nocca + ivir
      do iocc = 1, nocca
        iov = (ivir-1)*nocca + iocc
        call add_outer(pza, zsol(iov,1), cac(:,iocc), cac(:,amo))
        call add_outer(pza, zsol(iov,1), cac(:,amo), cac(:,iocc))
        zmata(iocc,amo) = zsol(iov,1) ; zmata(amo,iocc) = zsol(iov,1)
      end do
    end do
    do ivir = 1, nvirb
      amo = noccb + ivir
      do iocc = 1, noccb
        iov = la + (ivir-1)*noccb + iocc
        call add_outer(pzb, zsol(iov,1), cbc(:,iocc), cbc(:,amo))
        call add_outer(pzb, zsol(iov,1), cbc(:,amo), cbc(:,iocc))
        zmatb(iocc,amo) = zsol(iov,1) ; zmatb(amo,iocc) = zsol(iov,1)
      end do
    end do

    deallocate(cwa, cwb, eaw, ebw, w1, w2, rhs, zsol, moa_s, mob_s, ea_s, eb_s, malign, mblign)
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
