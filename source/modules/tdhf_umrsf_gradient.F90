!> UMRSF-TDDFT analytic nuclear gradient — clean-room implementation (branch uhf-grad-plan).
!>
!> VALIDATED SCOPE (2026-06-28, post-MILESTONE-A): PURE-HF response gradient (no XC). S1 gate PASS
!> (≤1e-5, RULES §11) EVERYWHERE — CH2(eq) 1.4e-7, SiH2 2.5e-8, H2CO 5.6e-7, C2H4-MEP(diradical) 4.6e-7,
!> stretched-CH2 5.7e-8. S3 also passes (CH2 1.5e-7). MILESTONE A (RULES §17) DONE: the energy's
!> get_jacobi (tdhf_mrsf_lib.F90) was unified onto the SAME cyclic+converged algorithm as umrsf_jacobi_smooth
!> (max|btt|<1e-12), so energy and analytic share ONE basin by construction (|δ vs td|→1e-15, smooth polish
!> a no-op). This CLOSED the distorted/diradical geometries that were gauge-limited (C2H4 1.19e-5→4.6e-7,
!> stretched-CH2 1.89e-4→5.7e-8) with NO regression (GAMESS-cross-checked; SCF ref unchanged). S2 remains
!> 1.3e-2 — NOT a gauge issue (analytic internally EXACT ≤1e-11): a confirmed STATE-TRACKING root-flip
!> (FD energy-rank vs analytic amplitude; sign-flip on the symmetry-breaking component) ⇒ §9-EXCLUDED,
!> pending the §11 character-following FD harness (independent of A/B). KNOWN-OPEN: (B) functionals NOT
!> implemented — the response is missing the XC kernel + grid-weight derivatives (MILESTONE B, after A).
!>
!> The C entry computes the reference (UHF-triplet) gradient via the reusable hf_gradient primitive
!> (grd1/grd2 with the converged DM_A/DM_B) — already FD-certified to 1.46e-7 — plus the response:
!> full-block Z-vector (oo+ov+vv) + full G^f (re-align, carries the dV/dC alignment Jacobian) +
!> W=½sym(G^f+G^z) + de_m1 (alignment overlap-Pulay). Reuses the clean-room energy/response lib
!> (umrsfcbc/int2_umrsf/umrsfmntoia/mrsfesum/get_jacobi); never reads the guarded RO-MRSF gradient.
module tdhf_umrsf_gradient_mod

  use precision, only: dp
  use types, only: information
  use grd2, only: grd2_compute_data_t
  use basis_tools, only: basis_set
  use mod_dft_molgrid, only: dft_grid_t

  implicit none

  character(len=*), parameter :: module_name = "tdhf_umrsf_gradient_mod"

  private
  public :: tdhf_umrsf_gradient_C

  !> ------- Stage-2 XC context (MILESTONE B / RULES §18) -------
  !> Set ONCE per gradient in umrsf_grad_run_gates (DFT runs only). The reference UKS XC kernel
  !> f_xc[ρ_ref] enters the response gradient ONLY through the reference-orbital relaxation (the
  !> MRSF response A-matrix has NO grid f_xc — energy is int2-only, so the f_xc·(X+Y)(X+Y) term is
  !> ABSENT). xc_meanfield_on toggles umrsf_meanfield's f_xc·P add-on (T3, propagates to the
  !> Z-vector Hessian / refrelax G^f / G^z / W). xc_refa/refb = reference density (defines the
  !> kernel). See DERIVATIONS/M2_xc_response.md.
  logical,                  save :: xc_meanfield_on = .false.
  type(dft_grid_t),         save :: xc_molgrid
  real(kind=dp), allocatable, save :: xc_refa(:,:), xc_refb(:,:)
  real(kind=dp),            save :: xc_thresh = 0.0_dp

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

  !> ------- ov-block Z-vector context (for the reusable pcg_optimize / minres_optimize) -------
  !> The matrix-free SPD solve of the ov–ov orbital Hessian A_{V,V} reuses OQP's stock solvers
  !> (source/pcg.F90 / source/minres.F90), whose matvec/precond callbacks take a single c_ptr `dat`
  !> context. This type carries everything umrsf_zov_matvec/umrsf_zov_precond need to apply A_{V,V}
  !> (= one umrsf_genfock_z over the ov DOFs) and the Jacobi preconditioner (ε_a−ε_i)⁻¹. The work
  !> arrays are members so the callbacks allocate nothing per matvec. See umrsf_zvector_iter.
  type :: umrsf_zov_ctx_t
    type(information), pointer :: infos => null()
    type(basis_set),   pointer :: basis => null()
    real(kind=dp), pointer :: cac(:,:) => null(), cbc(:,:) => null()
    real(kind=dp), pointer :: epsca(:) => null(), epscb(:) => null()
    real(kind=dp) :: hfscale_ref = 0.0_dp
    integer :: nbf = 0, ndof_ov = 0
    integer, allocatable :: dsp(:), dpr(:), dqr(:)        ! ov DOF maps (spin, p>q)
    real(kind=dp), allocatable :: pcinv(:)                ! ov Jacobi preconditioner 1/(ε_a−ε_i)
    real(kind=dp), allocatable :: za1(:,:), zb1(:,:), gza(:,:), gzb(:,:)  ! matvec scratch
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
    use mathlib, only: orthogonal_transform_sym, unpack_matrix, pack_matrix
    use eigen, only: diag_symm_full
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc, umrsfmntoia, mrsfesum, get_jacobi
    use tdhf_lib, only: iatogen
    use grd2, only: grd2_driver
    use grd1, only: grad_ee_overlap
    use constants, only: tol_int
    use oqp_linalg
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
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
    real(kind=dp), allocatable :: tua(:,:), tub(:,:)            ! standard-CIS T_u (SOMO-gate diagnostic only)
    real(kind=dp), allocatable :: peffa(:,:), peffb(:,:)        ! P_eff = P^Δ,u + ½ P_z (c03/c04 split)
    real(kind=dp), allocatable :: de_orb(:,:), de_w(:,:), de_m1(:,:), de_xc(:,:)
    real(kind=dp) :: omega_orb_chk, omega_orb, hfscale_ref
    real(kind=dp) :: omega_orb_tu, omega_orb_mine              ! SOMO gates: Tr(T_u F̃) / clean-room matvec
    real(kind=dp) :: dbg_zw                                     ! z weight in P_eff (c03/c04 split = 0.5)
    logical :: dbg_w2e, dbg_wrr, l_zov, l_gvt, l_m1             ! G̃ 2e/refrelax ; ablations: ov-only Z / V-transform G^f / M1
    logical :: l_zdense, l_zcmp                                 ! Z-vector solver: dense dgelss (UMRSF_ZDENSE) ; dense-vs-iter compare (UMRSF_ZCMP)
    logical :: l_g2efd, l_g2ecmp                                ! G̃ 2e: FD oracle (UMRSF_G2EFD) ; analytic-vs-FD compare (UMRSF_G2ECMP)
    logical :: dft_run, l_xck, l_xcg                           ! Stage-2 XC (§18): DFT run? f_xc kernel (T3)? diff-density XC grad (T2)?
    integer :: ia, ib, i, j
    ! Z-vector (relaxation): canonical MOs + relaxation density
    real(kind=dp), allocatable :: cac(:,:), cbc(:,:), epsca(:), epscb(:), pza(:,:), pzb(:,:)
    real(kind=dp), allocatable :: zmata(:,:), zmatb(:,:)

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
    real(kind=dp), allocatable :: xmat(:,:), amo(:,:), amo2e(:,:), xamp(:)
    real(kind=dp), allocatable :: brad(:,:,:)
    real(kind=dp), allocatable, target :: dens(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)
    type(int2_compute_t) :: int2_driver
    type(int2_umrsf_data_t), target :: int2_udata
    integer :: nbf, nbf2, nocca, noccb, nvirb, xvec_dim, mrst, nstates, tstate
    integer :: k, it, diag_index
    real(kind=dp) :: scale_exch, hfs, omega_recon, omega_2e_mv, omega_2e_tr, omega_eig
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
             xmat(nbf,nbf), amo(xvec_dim,1), amo2e(xvec_dim,1), xamp(xvec_dim), source=0.0_dp)
    allocate(dens(1,11,nbf,nbf), brad(11,nbf,nbf), source=0.0_dp)

    open(unit=iw, file=infos%log_filename, position="append")

    va = mo_a ; vb = mo_b ; ea = mo_energy_a ; eb = mo_energy_b
    call unpack_matrix(smat, smat_full, nbf, 'U')
    ! Corresponding-orbital (Jacobi) alignment — idempotent on already-aligned MOs.
    call get_jacobi(infos, va, ea, vb, eb, smat_full, nocca, wrk1, wrk2, 0)
    call get_jacobi(infos, va, ea, vb, eb, smat_full, nocca, wrk1, wrk2, 1)

    ! ---- SMOOTH/converged get_jacobi alignment (RULES §15 step 2.1; §17 aligner unification) ----
    ! POST-MILESTONE-A: get_jacobi (above) is now itself cyclic+converged (max|btt|<1e-12), so va,vb
    ! arrive ALREADY at the converged fixed point and this umrsf_jacobi_smooth call is a CONFIRMING
    ! NO-OP (polish size → ~1e-13). Kept as a defensive re-convergence + the authoritative residual/
    ! S-orthonormality gate. The analytic gradient needs the converged fixed point (within-seg btt → 0)
    ! so G^f = V G̃ Vᵀ is exact; by §17 the energy uses the SAME aligner ⇒ same basin by construction.
    block
      real(kind=dp), allocatable :: va_thr(:,:), vb_thr(:,:)
      real(kind=dp) :: off_thr, off_smooth, dva, orthoa, orthob
      allocate(va_thr, source=va) ; allocate(vb_thr, source=vb)
      ! off_thr/off_smooth = max within-seg |btt| (the get_jacobi STATIONARITY residual;
      ! c05: the off-diagonals themselves need NOT vanish, only btt → 0 at the fixed point).
      off_thr = within_seg_btt(va, vb, smat_full, nocca)
      call umrsf_jacobi_smooth(va, vb, smat_full, nocca, off_smooth)
      off_smooth = within_seg_btt(va, vb, smat_full, nocca)
      dva = sum(abs(va-va_thr)) + sum(abs(vb-vb_thr))
      orthoa = maxval(abs(matmul(transpose(va), matmul(smat_full, va)) - id_nbf(nbf)))
      orthob = maxval(abs(matmul(transpose(vb), matmul(smat_full, vb)) - id_nbf(nbf)))
      ! iw is already open (line above) — write directly; do NOT close (G1 gate below shares the bracket).
      write(iw,'(/2x,a)') '========= UMRSF gradient: SMOOTH get_jacobi alignment (§15 step 2.1) ========='
      write(iw,'(2x,a,es12.3)') 'within-seg max|btt| (stationarity) get_jacobi (§17 cyclic) = ', off_thr
      write(iw,'(2x,a,es12.3)') 'within-seg max|btt| (stationarity) SMOOTH (re-converged)   = ', off_smooth
      write(iw,'(2x,a,es12.3)') '||va_smooth - va_thr||_1 (polish size; →0 post-§17)        = ', dva
      write(iw,'(2x,a,2es12.3)') 'S-orthonormality max|CᵀSC−I| alpha/beta (smooth)        = ', orthoa, orthob
      if (off_smooth <= 1.0e-10_dp .and. max(orthoa,orthob) <= 1.0e-9_dp) then
        write(iw,'(2x,a)') 'VERDICT: smooth alignment CONVERGED (max|btt| → 0; S-orthonormal). '// &
                           'ω-reproduction = G1 gate below (omega_recon vs td_energies).'
      else
        write(iw,'(2x,a)') 'VERDICT: smooth alignment CHECK (see residuals above)'
      end if
      write(iw,'(2x,a)') '============================================================================='
      deallocate(va_thr, vb_thr)
    end block

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

    ! ---- RE-DIAGONALIZE A in the SMOOTH basis → genuine eigenvector xamp (§15 / c05) ----
    ! The stored bvec is the eigenvector in the energy's THRESHOLD basis (non-stationary here ~1.5e-6).
    ! The ov-only Z-vector + W machinery needs X to be a TRUE eigenvector of A in the SMOOTH basis.
    ! Dispatch: default = matrix-free Davidson (umrsf_track_amplitude_dav, ~tens of matvecs); UMRSF_TRKDENSE=1
    ! = the dense column-by-column oracle (umrsf_track_amplitude, nia matvecs + full diag); UMRSF_TRKCMP=1 =
    ! run BOTH and print the GATE max|x_dav - x_dense| + |dOmega| (must reproduce the dense xamp/omega <=1e-9).
    block
      character(len=24) :: e ; integer :: ios
      logical :: trk_dense, trk_cmp
      real(kind=dp), allocatable :: xamp_d(:) ; real(kind=dp) :: om_d, sgn_al
      trk_dense = .false. ; trk_cmp = .false.
      call get_environment_variable("UMRSF_TRKDENSE", e, status=ios) ; if (ios==0) trk_dense = (trim(e)=="1")
      call get_environment_variable("UMRSF_TRKCMP",   e, status=ios) ; if (ios==0) trk_cmp   = (trim(e)=="1")
      if (trk_cmp) then
        allocate(xamp_d(size(xamp)))
        call umrsf_track_amplitude(infos, int2_driver, va, vb, fa, fb, bvec(:,tstate), scale_exch, &
                                   hfs, spc_coco, spc_ovov, spc_coov, xamp_d, om_d)
        call umrsf_track_amplitude_dav(infos, int2_driver, va, vb, fa, fb, bvec(:,tstate), scale_exch, &
                                       hfs, spc_coco, spc_ovov, spc_coov, xamp, omega_eig)
        sgn_al = sign(1.0_dp, dot_product(xamp, xamp_d))     ! align global sign (both fixed to bvec_ref)
        write(iw,'(2x,a,2es12.3)') 'TRK GATE max|x_dav - x_dense| / |dOmega| (must <=1e-9) = ', &
          maxval(abs(sgn_al*xamp - xamp_d)), abs(omega_eig - om_d)
        deallocate(xamp_d)
      else if (trk_dense) then
        call umrsf_track_amplitude(infos, int2_driver, va, vb, fa, fb, bvec(:,tstate), scale_exch, &
                                   hfs, spc_coco, spc_ovov, spc_coov, xamp, omega_eig)
      else
        call umrsf_track_amplitude_dav(infos, int2_driver, va, vb, fa, fb, bvec(:,tstate), scale_exch, &
                                       hfs, spc_coco, spc_ovov, spc_coov, xamp, omega_eig)
      end if
    end block
    write(iw,'(/2x,a)') '----- smooth-basis amplitude re-diagonalization (genuine eigenvector) -----'
    write(iw,'(2x,a,f18.10)') 'omega (smooth-basis eigenvalue) = ', omega_eig
    write(iw,'(2x,a,f18.10)') 'omega stored (td_energies)      = ', td_en(tstate)
    write(iw,'(2x,a,es12.3)') '  |delta| vs td_energy          = ', abs(omega_eig-td_en(tstate))
    write(iw,'(2x,a,es12.3)') '  |overlap deficit| 1-|x·bvec|  = ', 1.0_dp-abs(dot_product(xamp,bvec(:,tstate)))
    ! The |delta| vs td_energy is the THRESHOLD-vs-SMOOTH alignment gap (OQP's energy uses the 1e-3
    ! threshold get_jacobi; the gradient uses the converged one). It is the irreducible §15-step-3
    ! "smoothness-permitting cross-check" residual; must be << 1e-5 (the S1 gate).
    if (abs(omega_eig-td_en(tstate)) <= 1.0e-5_dp) then
      write(iw,'(2x,a)') 'VERDICT: genuine smooth-basis eigenvector (threshold-vs-smooth gap << 1e-5 gate).'
    else
      write(iw,'(2x,a)') 'VERDICT: CHECK omega_eig vs td_energy delta (threshold/smooth gap too large).'
    end if

    ! ---- response matvec on the converged target-state amplitude ----
    call iatogen(xamp, xmat, nocca, noccb)
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
    omega_2e_mv = dot_product(xamp, amo2e(:,1))

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
    call iatogen(xamp, xmat, nocca, noccb)
    call mrsfesum(infos, xmat, fa, fb, amo, 1)
    omega_recon = dot_product(xamp, amo(:,1))

    write(iw,'(/2x,a)') '================ UMRSF gradient NON-FD gate G1 ================'
    write(iw,'(2x,a,i0,a,i0)') 'target_state = ', tstate, '   mrst = ', mrst
    write(iw,'(2x,a,f18.12)')  'X^T X (amplitude norm)          = ', dot_product(xamp,xamp)
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

    ! ===== SOMO-corrected unrelaxed difference density P_eff = sym(∂omega_orb/∂F̃) (MILESTONE C / M3) =====
    ! The MRSF orbital energy (mrsfesum, mrst=1/3) carries SOMO terms ∝ xlr=X(O1,O1) that the standard-CIS
    ! T_u OMITS ⇒ Tr(T_u F̃) ≠ omega_orb for SOMO-mixed (S2) states (gate gap ∝ xlr²; the S2 ~8e-3 error).
    ! omega_orb is LINEAR in F̃, so the correct unrelaxed difference density is P_eff = ∂omega_orb/∂F̃ (built
    ! by probing umrsf_orb_matvec with unit Focks; α occ-occ, β virt-virt like T_u; reduces to T_u when the
    ! SOMO-SOMO amplitudes vanish ⇒ S1/S3/non-SOMO untouched). P_eff REPLACES T_u as talpha/tbeta and
    ! propagates to pda/pdb (de_orb, refrelax), the frozen G̃ (gta=2 F̃ P_eff), G^f, the full-block Z, W.
    ! Model: DERIVATIONS/M3_s2_somo_diffdens.md, c09_peff_closure.py (≤1e-9), CAS c09_cas_peff.py.
    omega_orb = omega_recon - omega_2e_mv
    allocate(talpha(nbf,nbf), tbeta(nbf,nbf), pda(nbf,nbf), pdb(nbf,nbf), &
             tua(nbf,nbf), tub(nbf,nbf), source=0.0_dp)
    call iatogen(xamp, xmat, nocca, noccb)
    ! standard-CIS T_u (DIAGNOSTIC ONLY — its gate Tr(T_u F̃)−omega_orb is the SOMO tell: ~0 S1, ~2.4e-4 S2)
    do j = 1, nocca ; do i = 1, nocca ; do ia = noccb+1, nbf
      tua(i,j) = tua(i,j) - xmat(i,ia)*xmat(j,ia)
    end do ; end do ; end do
    do ib = noccb+1, nbf ; do ia = noccb+1, nbf ; do i = 1, nocca
      tub(ia,ib) = tub(ia,ib) + xmat(i,ia)*xmat(i,ib)
    end do ; end do ; end do
    omega_orb_tu = sum(tua*fa) + sum(tub*fb)
    ! SOMO-corrected difference density P_eff (the FIX) → talpha/tbeta
    call umrsf_build_peff(nbf, nocca, noccb, mrst, xmat, talpha, tbeta)
    omega_orb_chk = sum(talpha*fa) + sum(tbeta*fb)
    ! cross-check: my clean-room orbital matvec reproduces omega_orb (== the energy-path mrsfesum)
    call umrsf_orb_matvec(nbf, nocca, noccb, mrst, fa, fb, xmat, wrk1)
    omega_orb_mine = sum(xmat*wrk1)
    pda = matmul(matmul(va, talpha), transpose(va))     ! AO P_eff,α = C̃_α P_eff_α C̃_αᵀ
    pdb = matmul(matmul(vb, tbeta),  transpose(vb))     ! AO P_eff,β = C̃_β P_eff_β C̃_βᵀ

    ! Canonical SCF orbitals cac,cbc (= c05 Ca,Cb; the Z-vector frame) by diagonalizing the aligned
    ! MO-Fock fa,fb. These are occ/virt block-diagonal (SCF F_ai=0; the smooth get_jacobi only mixes
    ! within α-occ / β-virt), so the eigenvectors carry NO occ-virt mixing. cac is alignment-invariant.
    allocate(cac(nbf,nbf), cbc(nbf,nbf), epsca(nbf), epscb(nbf))
    block
      real(kind=dp), allocatable :: fac(:,:), fbc(:,:)
      integer :: ierr
      allocate(fac, source=fa) ; allocate(fbc, source=fb)
      call diag_symm_full(1, nbf, fac, nbf, epsca, ierr)   ! fac -> Vα (eigenvectors)
      call diag_symm_full(1, nbf, fbc, nbf, epscb, ierr)   ! fbc -> Vβ
      cac = matmul(va, fac) ; cbc = matmul(vb, fbc)
      deallocate(fac, fbc)
    end block

    hfscale_ref = 1.0_dp
    if (infos%control%hamilton >= 20) hfscale_ref = infos%dft%hfscale

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/2x,a)') '====== UMRSF response: c06 §16 CLOSED FORM (full G^f + FULL-BLOCK Z + W) ======'
    write(iw,'(2x,a,es12.3)') 'omega_orb gate |Tr(T_u  F̃) − omega_orb| (SOMO tell; S1~0, S2~2.4e-4) = ', abs(omega_orb_tu-omega_orb)
    write(iw,'(2x,a,es12.3)') 'omega_orb gate |Tr(P_eff F̃) − omega_orb| (the FIX; must → ~0)        = ', abs(omega_orb_chk-omega_orb)
    write(iw,'(2x,a,es12.3)') 'orbital matvec |X·esum_mine − omega_orb|  (clean-room == mrsfesum)    = ', abs(omega_orb_mine-omega_orb)

    ! Re-init int2 at the BASE geometry (the frozen-density 2e FD self-test left it displaced).
    call int2_driver%clean()
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_driver%schwarz = .false.

    allocate(pza(nbf,nbf), pzb(nbf,nbf), zmata(nbf,nbf), zmatb(nbf,nbf), source=0.0_dp)
    allocate(peffa(nbf,nbf), peffb(nbf,nbf), de_orb(3,natom), de_w(3,natom), de_m1(3,natom), &
             de_xc(3,natom), source=0.0_dp)

    ! Toggles: UMRSF_ZW (z weight in P_eff, default 0.5 = the c03/c04 split prefactor); UMRSF_W2E /
    ! UMRSF_WRR (include the 2e / refrelax pieces of the raw aligned G̃, default on). Ablations for the
    ! §16 fix: UMRSF_ZOV=1 → ov-only Z-vector (drops oo/vv → the wall); UMRSF_GVT=1 → G^f = V G̃ Vᵀ
    ! only (drops the dV/dC alignment Jacobian ΔG^f → the wall). Both off ⇒ the full c06 closed form.
    block
      character(len=16) :: e ; integer :: ios
      dbg_zw = 0.5_dp ; dbg_w2e = .true. ; dbg_wrr = .true. ; l_zov = .false. ; l_gvt = .false. ; l_m1 = .true.
      l_zdense = .false. ; l_zcmp = .false. ; l_g2efd = .false. ; l_g2ecmp = .false.
      call get_environment_variable("UMRSF_ZW", e, status=ios)
      if (ios==0) then ; read(e,*,iostat=ios) dbg_zw ; if (ios/=0) dbg_zw = 0.5_dp ; end if
      call get_environment_variable("UMRSF_W2E", e, status=ios) ; if (ios==0) dbg_w2e = (trim(e)/="0")
      call get_environment_variable("UMRSF_WRR", e, status=ios) ; if (ios==0) dbg_wrr = (trim(e)/="0")
      call get_environment_variable("UMRSF_ZOV", e, status=ios) ; if (ios==0) l_zov = (trim(e)=="1")
      call get_environment_variable("UMRSF_GVT", e, status=ios) ; if (ios==0) l_gvt = (trim(e)=="1")
      call get_environment_variable("UMRSF_M1",  e, status=ios) ; if (ios==0) l_m1  = (trim(e)/="0")
      ! Z-vector solver: default = matrix-free reduce+PCG (umrsf_zvector_iter). UMRSF_ZDENSE=1 →
      ! the dense dgelss oracle (umrsf_zvector_fullblock, rank-deficient-safe). UMRSF_ZCMP=1 → run BOTH
      ! and print max|z_iter − z_dense| (the perf-port GATE: reproduce the dense z to ≤1e-9).
      call get_environment_variable("UMRSF_ZDENSE", e, status=ios) ; if (ios==0) l_zdense = (trim(e)=="1")
      call get_environment_variable("UMRSF_ZCMP",   e, status=ios) ; if (ios==0) l_zcmp   = (trim(e)=="1")
      ! G̃ 2e gen-Fock: default = analytic (umrsf_g2e_analytic, 2 int2 builds). UMRSF_G2EFD=1 →
      ! the FD oracle umrsf_g2e_onesided (4·nbf² builds). UMRSF_G2ECMP=1 → run BOTH and print
      ! max|analytic − FD| per spin (the perf-port GATE: reproduce the FD g2e to ≤1e-9).
      call get_environment_variable("UMRSF_G2EFD",  e, status=ios) ; if (ios==0) l_g2efd  = (trim(e)=="1")
      call get_environment_variable("UMRSF_G2ECMP", e, status=ios) ; if (ios==0) l_g2ecmp = (trim(e)=="1")
    end block

    ! ---- Stage-2 XC context (RULES §18 / DERIVATIONS/M2_xc_response.md) ----
    ! DFT runs only. UMRSF_XCK (T3): add f_xc[ρ_ref]·P to the reference mean field ⇒ Z-vector Hessian /
    ! refrelax G^f / G^z / W. UMRSF_XCG (T2): add the difference-density XC gradient d/dR Tr(V_xc[ρ_ref]·P_eff).
    ! Both default ON for DFT (the response is otherwise missing all XC ⇒ the ~3.5e-2 BHHLYP S1 gap). The grid
    ! f_xc·(X+Y)(X+Y) term is ABSENT for MRSF (the energy A-matrix has no grid f_xc). Set the module XC context
    ! ONCE here (reference density + molGrid at the BASE geometry) for umrsf_meanfield (T3) + de_xc (T2).
    dft_run = (infos%control%hamilton == 20)
    l_xck = dft_run ; l_xcg = dft_run
    block
      character(len=16) :: e ; integer :: ios
      call get_environment_variable("UMRSF_XCK", e, status=ios) ; if (ios==0) l_xck = (trim(e)/="0")
      call get_environment_variable("UMRSF_XCG", e, status=ios) ; if (ios==0) l_xcg = (trim(e)/="0")
    end block
    xc_meanfield_on = .false.
    if (dft_run .and. (l_xck .or. l_xcg)) then
      call dft_initialize(infos, basis, xc_molgrid, verbose=.false.)
      if (allocated(xc_refa)) deallocate(xc_refa)
      if (allocated(xc_refb)) deallocate(xc_refb)
      allocate(xc_refa(nbf,nbf), xc_refb(nbf,nbf))
      call unpack_matrix(dmat_a, xc_refa, nbf, 'U')
      call unpack_matrix(dmat_b, xc_refb, nbf, 'U')
      xc_thresh = 0.0_dp
      xc_meanfield_on = l_xck       ! T3 inside umrsf_meanfield (refrelax / Z-Hessian / G^z / W)
    end if

    ! ============ c06 §16 closed form: G̃ → FULL G^f (re-align) → FULL-BLOCK Z → W ============
    ! THE FIX (both parts required; either alone leaves the ~1e-3 wall):
    !  (1) FULL-BLOCK Z-vector (oo+ov+vv) — the aligned 11-channel ω is NOT stationary to oo/vv rotations;
    !  (2) G^f carrying the dV/dC alignment Jacobian ΔG^f (via the numerical re-align, not V G̃ Vᵀ).
    ! W = Σ_σ C_σ ½sym(G^f+G^z)_σ C_σᵀ in the CANONICAL basis (full G^f). de_zexplicit folds into the
    ! P_eff = P^Δu + ½P_z orbital gradient (with the full-block P_z). Ablations: l_zov / l_gvt.
    block
      real(kind=dp), allocatable :: famoa(:,:), famob(:,:), ya(:,:), yb(:,:), tmp(:,:)
      real(kind=dp), allocatable :: gta(:,:), gtb(:,:), g2e(:,:), g2ea(:,:), g2eb(:,:)
      real(kind=dp), allocatable :: gfa(:,:), gfb(:,:), gza(:,:), gzb(:,:), wao(:,:), wpack(:)
      real(kind=dp) :: zrms, statio, tolw, wsa, wsb
      integer :: ij, ii, si, sj
      tolw = tol_int*log(10.0_dp)
      allocate(famoa(nbf,nbf), famob(nbf,nbf), ya(nbf,nbf), yb(nbf,nbf), tmp(nbf,nbf), &
               gta(nbf,nbf), gtb(nbf,nbf), g2e(nbf,nbf), g2ea(nbf,nbf), g2eb(nbf,nbf), &
               gfa(nbf,nbf), gfb(nbf,nbf), &
               gza(nbf,nbf), gzb(nbf,nbf), wao(nbf,nbf), wpack(nbf2), source=0.0_dp)

      ! ---- G̃_σ (raw aligned-basis generalized Fock = c04 factors, va/vb basis) ----
      ! frozen 2 F̃ T_u (F̃ = C̃ᵀ F^ref C̃) + 2e one-sided + refrelax 2(C̃ᵀ Y C̃)|occ, Y_σ = J[P^Δu]−hfscale·K[P^Δu_σ]
      call orthogonal_transform_sym(nbf, nbf, fock_a, va, nbf, scr) ; call unpack_matrix(scr, famoa)
      call orthogonal_transform_sym(nbf, nbf, fock_b, vb, nbf, scr) ; call unpack_matrix(scr, famob)
      call umrsf_meanfield(basis, infos, pda, pdb, hfscale_ref, ya, yb)
      gta = 2.0_dp*matmul(famoa, talpha)
      gtb = 2.0_dp*matmul(famob, tbeta)
      ! ---- G̃ 2e channel-adjoint: ANALYTIC (default, 2 int2 builds) / FD oracle / GATE ----
      if (dbg_w2e) then
        if (l_g2ecmp) then       ! GATE: analytic vs FD oracle (reproduce ≤1e-9), use analytic
          call umrsf_g2e_analytic(infos, int2_driver, va, vb, xamp, scale_exch, g2ea, g2eb)
          call umrsf_g2e_onesided(infos, int2_driver, va, vb, xamp, scale_exch, 1, g2e)
          write(iw,'(/2x,a,es12.3)') 'G2e GATE max|analytic − FD| α (must ≤1e-9) = ', maxval(abs(g2ea-g2e))
          call umrsf_g2e_onesided(infos, int2_driver, va, vb, xamp, scale_exch, 2, g2e)
          write(iw,'(2x,a,es12.3)')  'G2e GATE max|analytic − FD| β (must ≤1e-9) = ', maxval(abs(g2eb-g2e))
        else if (l_g2efd) then   ! FD oracle fallback (4·nbf² builds)
          call umrsf_g2e_onesided(infos, int2_driver, va, vb, xamp, scale_exch, 1, g2ea)
          call umrsf_g2e_onesided(infos, int2_driver, va, vb, xamp, scale_exch, 2, g2eb)
        else                     ! production: analytic, both spins one call
          call umrsf_g2e_analytic(infos, int2_driver, va, vb, xamp, scale_exch, g2ea, g2eb)
        end if
        gta = gta + g2ea ; gtb = gtb + g2eb
      end if
      if (dbg_wrr) then ; tmp = matmul(transpose(va), matmul(ya, va))
        gta(:,1:nocca) = gta(:,1:nocca) + 2.0_dp*tmp(:,1:nocca) ; end if
      if (dbg_wrr) then ; tmp = matmul(transpose(vb), matmul(yb, vb))
        gtb(:,1:noccb) = gtb(:,1:noccb) + 2.0_dp*tmp(:,1:noccb) ; end if

      ! ---- FULL G^f (canonical basis), carrying the dV/dC alignment Jacobian (M1) ----
      if (l_gvt) then
        ! ABLATION: G^f = V G̃ Vᵀ only (drop ΔG^f) — reproduces the wall.  V_σ = C_can,σᵀ S C̃_σ.
        tmp = matmul(transpose(cac), matmul(smat_full, va))     ! V_α
        gfa = matmul(tmp, matmul(gta, transpose(tmp)))
        tmp = matmul(transpose(cbc), matmul(smat_full, vb))     ! V_β
        gfb = matmul(tmp, matmul(gtb, transpose(tmp)))
      else
        call umrsf_genfock_full(infos, cac, cbc, va, vb, smat_full, gta, gtb, nocca, gfa, gfb)
      end if

      ! ---- FULL-BLOCK Z-vector: M z = −R, R = antisym(G^f) over all p>q (l_zov: ov-only ablation) ----
      ! Default = matrix-free reduce+PCG (umrsf_zvector_iter, ~tens of Fock builds). UMRSF_ZDENSE=1
      ! = the dense dgelss oracle (~ndof Fock builds, rank-deficient-safe). UMRSF_ZCMP=1 = run BOTH and
      ! print max|z_iter − z_dense| (the perf-port GATE; the gradient uses the iterative z).
      if (l_zcmp) then
        block
          real(kind=dp), allocatable :: pzad(:,:), pzbd(:,:), zmad(:,:), zmbd(:,:)
          real(kind=dp) :: zrd, std, dza, dzb
          allocate(pzad(nbf,nbf), pzbd(nbf,nbf), zmad(nbf,nbf), zmbd(nbf,nbf))
          call umrsf_zvector_fullblock(infos, basis, cac, cbc, epsca, epscb, gfa, gfb, &
                                       hfscale_ref, l_zov, pzad, pzbd, zmad, zmbd, zrd, std)
          call umrsf_zvector_iter(infos, basis, cac, cbc, epsca, epscb, gfa, gfb, &
                                  hfscale_ref, l_zov, pza, pzb, zmata, zmatb, zrms, statio)
          dza = maxval(abs(zmata - zmad)) ; dzb = maxval(abs(zmatb - zmbd))
          write(iw,'(2x,a,2es12.3)') 'Z-solver GATE max|z_iter − z_dense| α/β (must ≤1e-9)   = ', dza, dzb
          write(iw,'(2x,a,2es12.3)') 'Z-solver      dense statio / iter statio              = ', std, statio
          deallocate(pzad, pzbd, zmad, zmbd)
        end block
      else if (l_zdense) then
        call umrsf_zvector_fullblock(infos, basis, cac, cbc, epsca, epscb, gfa, gfb, &
                                     hfscale_ref, l_zov, pza, pzb, zmata, zmatb, zrms, statio)
      else
        call umrsf_zvector_iter(infos, basis, cac, cbc, epsca, epscb, gfa, gfb, &
                                hfscale_ref, l_zov, pza, pzb, zmata, zmatb, zrms, statio)
      end if

      ! ---- G^z (full-block) and W_ao = Σ_σ C_σ ½sym(G^f+G^z)_σ C_σᵀ ; de_w = −Tr(W S^x) ----
      call umrsf_genfock_z(infos, basis, cac, cbc, epsca, epscb, zmata, zmatb, hfscale_ref, gza, gzb)
      wao = matmul(cac, matmul(0.25_dp*(gfa+transpose(gfa)), transpose(cac))) &
          + matmul(cbc, matmul(0.25_dp*(gfb+transpose(gfb)), transpose(cbc))) &
          + matmul(cac, matmul(0.25_dp*(gza+transpose(gza)), transpose(cac))) &
          + matmul(cbc, matmul(0.25_dp*(gzb+transpose(gzb)), transpose(cbc)))
      de_w = 0.0_dp
      call pack_matrix(-wao, wpack, 'U')
      ij = 0 ; do ii = 1, nbf ; ij = ij + ii ; wpack(ij) = 0.5_dp*wpack(ij) ; end do
      call grad_ee_overlap(basis, wpack, de_w, logtol=tolw)

      ! ---- diagnostics (c06 cross-checks) ----
      ! within-segment antisym of G̃ — the c06 phenomenon (c05 toy ≈0; real MRSF ≈ 1e-3).
      wsa = 0.0_dp ; wsb = 0.0_dp
      do sj = 1, nocca-1 ; do si = 1, nocca-1
        wsa = max(wsa, abs(0.5_dp*(gta(si,sj)-gta(sj,si)))) ; end do ; end do
      do sj = nocca, nbf ; do si = nocca, nbf
        wsb = max(wsb, abs(0.5_dp*(gtb(si,sj)-gtb(sj,si)))) ; end do ; end do
      write(iw,'(/2x,a)') '----- §16 full-block Z + full G^f diagnostics -----'
      write(iw,'(2x,a,l1,a,l1)') 'ablations: ov-only Z = ', l_zov, '   V-transform G^f = ', l_gvt
      write(iw,'(2x,a,2es12.3)') 'within-seg antisym G̃ α/β (c06: ~1e-3 NONZERO)  = ', wsa, wsb
      block
        real(kind=dp) :: aoo, avv, boo, bvv
        integer :: rr, ss
        aoo = 0.0_dp ; avv = 0.0_dp ; boo = 0.0_dp ; bvv = 0.0_dp
        do ss = 1, nbf ; do rr = 1, nbf
          if (rr<=nocca .and. ss<=nocca) aoo = max(aoo, abs(gfa(rr,ss)-gfa(ss,rr)))
          if (rr>nocca  .and. ss>nocca ) avv = max(avv, abs(gfa(rr,ss)-gfa(ss,rr)))
          if (rr<=noccb .and. ss<=noccb) boo = max(boo, abs(gfb(rr,ss)-gfb(ss,rr)))
          if (rr>noccb  .and. ss>noccb ) bvv = max(bvv, abs(gfb(rr,ss)-gfb(ss,rr)))
        end do ; end do
        write(iw,'(2x,a,2es12.3)') 'canonical G^f antisym  α oo / β vv (≠0 ⇒ need full-block) = ', aoo, bvv
        write(iw,'(2x,a,2es12.3)') '          (α vv / β oo, ≈0 for canonical CIS)             = ', avv, boo
      end block
      write(iw,'(2x,a,es12.3)') 'full-block z RMS                               = ', zrms
      write(iw,'(2x,a,es12.3)') 'Z-vector stationarity ||antisym(G^f+G^z)||      = ', statio
      write(iw,'(2x,a,es12.3)') '||W_ao||                                       = ', sqrt(sum(wao**2))
      deallocate(famoa, famob, ya, yb, tmp, gta, gtb, g2e, gfa, gfb, gza, gzb, wao, wpack)
    end block

    ! de_m1 = the alignment's EXPLICIT-S response −∂(ω∘align)/∂S·∂S/∂x (the overlap-Pulay of the
    ! get_jacobi alignment): re-align the canonical orbitals at S(x±θ) (orbitals/ERIs/F^ref frozen at
    ! base; D^ref get_jacobi-invariant), central-difference ω. SMOOTH (converged) re-alignment, faithful
    ! to the model's de_explicit (which re-aligns). c06: NONZERO (~the residual after the fixed-alignment
    ! de2e+de_orb), unlike the c05 toy (within-seg invariant ⇒ de_m1≈0). Done while int2_driver is live.
    if (l_m1) call umrsf_m1_overlap_grad(infos, int2_driver, basis, cac, cbc, va, vb, &
                                         fock_a, fock_b, smat_full, xamp, scale_exch, de_m1)

    call int2_driver%clean()

    ! de_orb = umrsf_orbital_grad(P_eff = P^Δ,u + ½ P_z) = de_explicit_orbital + de_zexplicit
    ! (the ½ on P_z is the c03/c04 split; the OTHER ½ lives in the W z-coupling G^z).
    peffa = pda + dbg_zw*pza
    peffb = pdb + dbg_zw*pzb
    call umrsf_orbital_grad(infos, basis, peffa, peffb, dmat_a, dmat_b, hfscale_ref, de_orb)

    ! de_xc (T2, RULES §18): difference-density XC skeleton gradient  d/dR Tr(V_xc[ρ_ref]·P_eff)  — the XC
    ! analogue of de_orb's 2e mean-field (J−cK) part, with the SAME P_eff. utddft_xc_gradient with
    ! do_ground_state=.false. (reference XC grad already in hf_gradient/dftder ⇒ no double counting) and
    ! do_fxc=.false. (no transition-density term: MRSF A has no grid f_xc). dedft sign convention matches
    ! the energy (validated by per-term OQP-FD). Reuses the validated dftlib TD-XC gradient consumer.
    if (dft_run .and. l_xcg) then
      block
        real(kind=dp), allocatable :: pxa(:,:,:), pxb(:,:,:)
        allocate(pxa(nbf,nbf,1), pxb(nbf,nbf,1))
        pxa(:,:,1) = peffa ; pxb(:,:,1) = peffb
        call utddft_xc_gradient(basis, xc_molgrid, de_xc, xc_refa, xc_refb, pxa, pxb, &
                                nMtx=1, threshold=xc_thresh, infos=infos, do_ground_state=.false.)
        deallocate(pxa, pxb)
      end block
    end if

    write(iw,'(2x,a)') '   atom  comp        de_2e               de_orb               de_w                de_m1                de_xc'
    do iat = 1, natom ; do icmp = 1, 3
      write(iw,'(2x,2i5,5es21.11)') iat, icmp, de2e(icmp,iat), de_orb(icmp,iat), de_w(icmp,iat), &
                                    de_m1(icmp,iat), de_xc(icmp,iat)
    end do ; end do
    write(iw,'(2x,a)') '================================================================================='
    close(iw)

    ! full response = transition-2e (Pulay) + orbital(1e+mean-field, P_eff) + W + M1 (alignment-S Pulay)
    !                 + de_xc (Stage-2: XC kernel response; T3 folded into de_orb/de_w via umrsf_meanfield)
    de2e_out = de2e + de_orb + de_w + de_m1 + de_xc

    if (xc_meanfield_on .or. (dft_run .and. l_xcg)) then
      call dftclean(infos)
      if (allocated(xc_refa)) deallocate(xc_refa)
      if (allocated(xc_refb)) deallocate(xc_refb)
      xc_meanfield_on = .false.
    end if

    deallocate(va, vb, fa, fb, smat_full, ea, eb, wrk1, wrk2, scr, xmat, amo, amo2e, xamp, dens, brad)
    deallocate(de2e, de2e_fd, densym)
    deallocate(talpha, tbeta, tua, tub, pda, pdb, peffa, peffb, de_orb, de_w, de_m1, de_xc)
    deallocate(cac, cbc, epsca, epscb, pza, pzb, zmata, zmatb)

  end subroutine umrsf_grad_run_gates

!###############################################################################
!> @brief Orbital part of the UMRSF response matvec (A_orb·X)[i,j] on the occα×virβ block — clean-room
!>        transcription of the energy's mrsfesum (mrst=1/3; tdhf_mrsf_lib.F90 READ-OK, re-derived as
!>        DERIVATIONS/c09_somo_diffdens.py:mrsf_orb_matvec). LINEAR in the aligned-MO Fock (fij=F̃α occ-occ,
!>        fab=F̃β virt-virt) and in the amplitude wrk=X. umrsf_build_peff probes this with UNIT Focks to
!>        get P_eff = ∂omega_orb/∂F̃; calling it with the real F̃ reproduces omega_orb (== mrsfesum) — a gate.
!>        The SOMO terms (∝ xlr=X(O1,O1), O1=nocca-1 O2=nocca) are exactly what the standard-CIS T_u omits.
  subroutine umrsf_orb_matvec(nbf, nocca, noccb, mrst, fij, fab, wrk, w)
    implicit none
    integer, intent(in) :: nbf, nocca, noccb, mrst
    real(kind=dp), intent(in) :: fij(nbf,nbf), fab(nbf,nbf), wrk(nbf,nbf)
    real(kind=dp), intent(out) :: w(nbf,nbf)
    real(kind=dp), allocatable :: scr(:,:), tmp1(:,:)
    real(kind=dp) :: dumn, xlr
    real(kind=dp), parameter :: sqrt2 = 1.0_dp/sqrt(2.0_dp)
    integer :: lr1, lr2, i, j

    allocate(scr(nbf,nbf), tmp1(nbf,nbf), source=0.0_dp)
    lr1 = nocca-1 ; lr2 = nocca                       ! SOMOs (1-based, ⊂ both occα and virβ)
    scr = wrk ; scr(lr1,lr1) = 0.0_dp ; scr(lr2,lr2) = 0.0_dp
    xlr = wrk(lr1,lr1)
    w = 0.0_dp
    ! standard CIS on scr:  Σ_a scr(i,a)F̃β(a,j) − Σ_k F̃α(i,k)scr(k,j)   (i∈occα, j∈virβ)
    tmp1(1:nocca,noccb+1:nbf) = &
        matmul(scr(1:nocca,noccb+1:nbf), fab(noccb+1:nbf,noccb+1:nbf)) &
      - matmul(fij(1:nocca,1:nocca),     scr(1:nocca,noccb+1:nbf))
    if (mrst == 1) then
      do j = noccb+1, nbf ; do i = 1, nocca
        w(i,j) = tmp1(i,j)
        if (i==lr1) w(i,j) = w(i,j) + fab(j,lr1)*xlr*sqrt2
        if (i==lr2) w(i,j) = w(i,j) - fab(j,lr2)*xlr*sqrt2
        if (j==lr1) w(i,j) = w(i,j) - fij(i,lr1)*xlr*sqrt2
        if (j==lr2) w(i,j) = w(i,j) + fij(i,lr2)*xlr*sqrt2
      end do ; end do
      dumn = - dot_product(fij(lr1,1:nocca), scr(1:nocca,lr1)) &
             + dot_product(fij(lr2,1:nocca), scr(1:nocca,lr2)) &
             + dot_product(fab(lr1,noccb+1:nbf), scr(lr1,noccb+1:nbf)) &
             - dot_product(fab(lr2,noccb+1:nbf), scr(lr2,noccb+1:nbf))
      w(lr1,lr1) = dumn*sqrt2 + xlr*(fab(lr1,lr1)+fab(lr2,lr2)-fij(lr1,lr1)-fij(lr2,lr2))*0.5_dp
      w(lr2,lr2) = 0.0_dp
    else                                              ! mrst == 3 (triplet response)
      do j = noccb+1, nbf ; do i = 1, nocca
        w(i,j) = tmp1(i,j)
        if (i==lr1) w(i,j) = w(i,j) + fab(j,lr1)*xlr*sqrt2
        if (i==lr2) w(i,j) = w(i,j) + fab(j,lr2)*xlr*sqrt2
        if (j==lr1) w(i,j) = w(i,j) - fij(i,lr1)*xlr*sqrt2
        if (j==lr2) w(i,j) = w(i,j) - fij(i,lr2)*xlr*sqrt2
      end do ; end do
      dumn = - dot_product(fij(lr1,1:nocca), scr(1:nocca,lr1)) &
             - dot_product(fij(lr2,1:nocca), scr(1:nocca,lr2)) &
             + dot_product(fab(lr1,noccb+1:nbf), scr(lr1,noccb+1:nbf)) &
             + dot_product(fab(lr2,noccb+1:nbf), scr(lr2,noccb+1:nbf))
      w(lr1,lr1) = dumn*sqrt2 + xlr*(fab(lr1,lr1)+fab(lr2,lr2)-fij(lr1,lr1)-fij(lr2,lr2))*0.5_dp
      w(lr2,lr1) = 0.0_dp ; w(lr1,lr2) = 0.0_dp ; w(lr2,lr2) = 0.0_dp
    end if
    deallocate(scr, tmp1)
  end subroutine umrsf_orb_matvec

!###############################################################################
!> @brief P_eff = sym(∂omega_orb/∂F̃) — the SOMO-corrected unrelaxed difference density that REPLACES the
!>        standard-CIS T_u (talpha/tbeta) for SOMO-mixed states (closes S2). omega_orb = Σ X·orb_matvec(F̃,X)
!>        is LINEAR in F̃, so P_raw_σ[p,q] = omega_orb evaluated with F̃_σ = unit E_pq (the only blocks F̃
!>        enters: α occ-occ p,q∈1:nocca; β virt-virt p,q∈noccb+1:nbf). P_eff = sym(P_raw) (physical density;
!>        Tr(sym(P)·F̃)=Tr(P·F̃) for symmetric F̃ ⇒ gate Tr(P_eff F̃)=omega_orb preserved). Reduces to T_u
!>        when X(O1,O1)=X(O2,O2)=0. Model: DERIVATIONS/c09_peff_closure.py:Praw (≤1e-9), CAS c09_cas_peff.py.
  subroutine umrsf_build_peff(nbf, nocca, noccb, mrst, xmat, peffa, peffb)
    implicit none
    integer, intent(in) :: nbf, nocca, noccb, mrst
    real(kind=dp), intent(in) :: xmat(nbf,nbf)
    real(kind=dp), intent(out) :: peffa(nbf,nbf), peffb(nbf,nbf)
    real(kind=dp), allocatable :: epq(:,:), zero(:,:), w(:,:)
    integer :: p, q

    allocate(epq(nbf,nbf), zero(nbf,nbf), w(nbf,nbf), source=0.0_dp)
    peffa = 0.0_dp ; peffb = 0.0_dp
    ! α: ∂omega_orb/∂F̃α[p,q], nonzero only on the occ-occ block
    do q = 1, nocca ; do p = 1, nocca
      epq = 0.0_dp ; epq(p,q) = 1.0_dp
      call umrsf_orb_matvec(nbf, nocca, noccb, mrst, epq, zero, xmat, w)
      peffa(p,q) = sum(xmat*w)
    end do ; end do
    ! β: ∂omega_orb/∂F̃β[p,q], nonzero only on the virt-virt block
    do q = noccb+1, nbf ; do p = noccb+1, nbf
      epq = 0.0_dp ; epq(p,q) = 1.0_dp
      call umrsf_orb_matvec(nbf, nocca, noccb, mrst, zero, epq, xmat, w)
      peffb(p,q) = sum(xmat*w)
    end do ; end do
    ! symmetrize → physical difference density (block structure preserved: α occ-occ, β virt-virt)
    peffa = 0.5_dp*(peffa + transpose(peffa))
    peffb = 0.5_dp*(peffb + transpose(peffb))
    deallocate(epq, zero, w)
  end subroutine umrsf_build_peff

!###############################################################################
!> SMOOTH (converged) corresponding-orbital alignment — clean-room port of the c05 MODEL
!> (DERIVATIONS/c05_uhf_full_gradient.py:get_jacobi_align), faithful to the energy's get_jacobi
!> 2×2 angle law + segments (tdhf_mrsf_lib.F90, READ-OK) but run to CONVERGENCE: NO 1e-3 threshold,
!> NO min-|θ| early exit — cyclic sweeps until max|btt| < tol. RULES §15: the threshold get_jacobi is
!> non-smooth and was THE contamination of the earlier numerical-RHS Z-vector; the analytic gradient
!> path needs the converged alignment so the within-segment generalized-Fock blocks vanish (then the
!> M1 V-transform G^f = V G̃ Vᵀ is exact). seg0 rotates ALPHA columns {1..nocca-1} (closed+O1); seg1
!> rotates BETA columns {nocca..nbf} (O2+virt). s_mo = vaᵀ S vb (columns normalized once). 2×2 law per
!> pair (i,j): aa=s(i,i) bb=s(j,j) cc=s(i,j) dd=s(j,i); att=½(aa²+bb²−cc²−dd²); seg0 btt=aa·dd−bb·cc,
!> seg1 btt=aa·cc−bb·dd; θ=½ atan2(btt,att). seg0 rotates va cols(i,j) + s_mo ROWS(i,j); seg1 rotates
!> vb cols(i,j) + s_mo COLS(i,j). Sign-fix (flip β column so s_mo(i,i)≥0). Start from the threshold-
!> aligned va,vb (already in the right basin + sign convention) ⇒ polishing preserves ω (span-invariance)
!> and reaches the EXACT fixed point. Returns the post-convergence max within-seg |btt| in offmax.
  subroutine umrsf_jacobi_smooth(va, vb, smat_full, nocca, offmax)
    implicit none
    real(kind=dp), intent(inout) :: va(:,:), vb(:,:)
    real(kind=dp), intent(in) :: smat_full(:,:)
    integer, intent(in) :: nocca
    real(kind=dp), intent(out) :: offmax
    real(kind=dp), allocatable :: s_mo(:,:), ri(:), rj(:)
    integer :: nbf, i, j, p, sweep, slo, shi, seg
    integer, parameter :: maxsweep = 500
    real(kind=dp), parameter :: tol = 1.0e-12_dp
    real(kind=dp) :: aa, bb, cc, dd, att, btt, th, ct, st, off, nrm

    nbf = size(va,1)
    allocate(s_mo(nbf,nbf), ri(nbf), rj(nbf))
    ! s_mo = vaᵀ S vb ; normalize columns once (faithful to get_jacobi)
    s_mo = matmul(transpose(va), matmul(smat_full, vb))
    do i = 1, nbf
      nrm = max(norm2(s_mo(:,i)), 1.0e-10_dp)
      s_mo(:,i) = s_mo(:,i)/nrm
    end do

    do seg = 0, 1
      if (seg == 0) then ; slo = 1 ; shi = nocca-1
      else               ; slo = nocca ; shi = nbf ; end if
      do sweep = 1, maxsweep
        off = 0.0_dp
        do i = slo, shi-1
          do j = i+1, shi
            aa = s_mo(i,i) ; bb = s_mo(j,j) ; cc = s_mo(i,j) ; dd = s_mo(j,i)
            att = 0.5_dp*(aa*aa + bb*bb - cc*cc - dd*dd)
            if (seg == 0) then ; btt = aa*dd - bb*cc
            else               ; btt = aa*cc - bb*dd ; end if
            off = max(off, abs(btt))
            th = 0.5_dp*atan2(btt, att)
            ct = cos(th) ; st = sin(th)
            if (seg == 0) then
              ri = va(:,i) ; rj = va(:,j)            ! rotate va columns (i,j)
              va(:,i) = ct*ri + st*rj ; va(:,j) = -st*ri + ct*rj
              ri = s_mo(i,:) ; rj = s_mo(j,:)        ! rotate s_mo ROWS (i,j)
              s_mo(i,:) = ct*ri + st*rj ; s_mo(j,:) = -st*ri + ct*rj
            else
              ri = vb(:,i) ; rj = vb(:,j)            ! rotate vb columns (i,j)
              vb(:,i) = ct*ri + st*rj ; vb(:,j) = -st*ri + ct*rj
              ri = s_mo(:,i) ; rj = s_mo(:,j)        ! rotate s_mo COLUMNS (i,j)
              s_mo(:,i) = ct*ri + st*rj ; s_mo(:,j) = -st*ri + ct*rj
            end if
          end do
        end do
        if (off < tol) exit
      end do
    end do

    ! sign fix: flip β column p so s_mo(p,p) ≥ 0 (faithful to check_sign / c05)
    do p = 1, nbf
      if (s_mo(p,p) < 0.0_dp) then
        vb(:,p) = -vb(:,p)
        s_mo(:,p) = -s_mo(:,p)
      end if
    end do

    ! authoritative post-convergence residual: max within-segment |s_mo off-diagonal|
    offmax = 0.0_dp
    do seg = 0, 1
      if (seg == 0) then ; slo = 1 ; shi = nocca-1
      else               ; slo = nocca ; shi = nbf ; end if
      do i = slo, shi
        do j = slo, shi
          if (i /= j) offmax = max(offmax, abs(s_mo(i,j)))
        end do
      end do
    end do
    deallocate(s_mo, ri, rj)
  end subroutine umrsf_jacobi_smooth

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
    use mathlib, only: unpack_matrix
    use constants, only: tol_int
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), va(:,:), vb(:,:)
    real(kind=dp), intent(in) :: fock_a(:), fock_b(:), smat_full(:,:), xv(:), scale_exch
    real(kind=dp), intent(out) :: de_m1(:,:)
    real(kind=dp), allocatable :: cwa(:,:), cwb(:,:), spack(:), sfull(:,:)
    real(kind=dp), allocatable :: hbuf(:), tbuf(:), zq(:)
    integer :: nbf, nbf2, nocca, natom, iat, icmp, p
    real(kind=dp) :: tol, th, omp, omm, schk, ocheck, dum

    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2 ; nocca = infos%mol_prop%nelec_a
    natom = ubound(infos%atoms%zn,1) ; tol = tol_int*log(10.0_dp) ; th = 1.0d-3
    allocate(cwa(nbf,nbf), cwb(nbf,nbf), spack(nbf2), sfull(nbf,nbf), &
             hbuf(nbf2), tbuf(nbf2), zq(natom), source=0.0_dp)
    zq = infos%atoms%zn - infos%basis%ecp_zn_num

    open(unit=iw, file=infos%log_filename, position="append")
    ! convention sanity: omp_hst overlap at base reproduces stored S; sign-fixed get_jacobi(cac,S_base)==ω_base
    call omp_hst(basis, infos%atoms%xyz, zq, hbuf, spack, tbuf, logtol=tol, &
                 comm=infos%mpiinfo%comm, usempi=infos%mpiinfo%usempi)
    call unpack_matrix(spack, sfull, nbf, 'U')
    schk = maxval(abs(sfull - smat_full))
    cwa = cac ; cwb = cbc
    call umrsf_jacobi_smooth(cwa, cwb, sfull, nocca, dum)
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
        call umrsf_jacobi_smooth(cwa, cwb, sfull, nocca, dum)
        call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
        call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omp)
        ! -θ
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) - th
        call basis%init_shell_centers()
        call omp_hst(basis, infos%atoms%xyz, zq, hbuf, spack, tbuf, logtol=tol, comm=infos%mpiinfo%comm, usempi=infos%mpiinfo%usempi)
        infos%atoms%xyz(icmp,iat) = infos%atoms%xyz(icmp,iat) + th
        call basis%init_shell_centers() ; call unpack_matrix(spack, sfull, nbf, 'U')
        cwa = cac ; cwb = cbc
        call umrsf_jacobi_smooth(cwa, cwb, sfull, nocca, dum)
        call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
        call umrsf_omega_eval(infos, idrv, cwa, cwb, fock_a, fock_b, xv, scale_exch, omm)
        de_m1(icmp,iat) = (omp - omm)/(2.0_dp*th)
      end do
    end do
    write(iw,'(2x,a)') 'M1 alignment-overlap (dS/dx) response computed (de_m1).'
    close(iw)
    deallocate(cwa, cwb, spack, sfull, hbuf, tbuf, zq)
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
!> Reference UHF mean field of an arbitrary AO density (P_a,P_b full): Y_σ = J[P_a+P_b] − hfscale·K[P_σ]
!> (one fock_jk build, screening tol from infos). Used for the refrelaxation generalized-Fock blocks
!> (G̃ uses P^Δ,u; G^z uses ½P_z). fock_jk's own screened driver is fine for the density-cross terms.
  subroutine umrsf_meanfield(basis, infos, pa_full, pb_full, hfscale, ya_full, yb_full)
    use scf_addons, only: fock_jk
    use mathlib, only: pack_matrix, unpack_matrix
    use mod_dft_gridint_fxc, only: utddft_fxc
    implicit none
    type(basis_set), intent(inout) :: basis
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: pa_full(:,:), pb_full(:,:), hfscale
    real(kind=dp), intent(out) :: ya_full(:,:), yb_full(:,:)
    real(kind=dp), allocatable :: dens(:,:), fout(:,:)
    integer :: nbf, nbf2
    nbf = basis%nbf ; nbf2 = nbf*(nbf+1)/2
    allocate(dens(nbf2,2), fout(nbf2,2), source=0.0_dp)
    call pack_matrix(pa_full, dens(:,1), 'U') ; call pack_matrix(pb_full, dens(:,2), 'U')
    call fock_jk(basis, dens, fout, infos, scale_exch=hfscale, scale_coul=1.0_dp)
    call unpack_matrix(fout(:,1), ya_full) ; call unpack_matrix(fout(:,2), yb_full)
    deallocate(dens, fout)
    ! Stage-2 (RULES §18): add the reference UKS XC kernel response f_xc[ρ_ref]·P (T3). The reference
    ! mean field becomes G_σ[P] = J[P] − hfscale·K[P_σ] + (f_xc·P)_σ (collinear UKS f_xc, spin-conserving
    ! P). One grid pass per call; propagates to the Z-vector Hessian, refrelax G^f, G^z and W.
    if (xc_meanfield_on) then
      block
        real(kind=dp), allocatable :: fxa(:,:,:), fxb(:,:,:), dxa(:,:,:), dxb(:,:,:)
        allocate(fxa(nbf,nbf,1), fxb(nbf,nbf,1), dxa(nbf,nbf,1), dxb(nbf,nbf,1), source=0.0_dp)
        dxa(:,:,1) = pa_full ; dxb(:,:,1) = pb_full
        call utddft_fxc(basis, xc_molgrid, .false., xc_refa, xc_refb, &
                        fxa, fxb, dxa, dxb, 1, xc_thresh, infos)
        ya_full = ya_full + fxa(:,:,1) ; yb_full = yb_full + fxb(:,:,1)
        deallocate(fxa, fxb, dxa, dxb)
      end block
    end if
  end subroutine umrsf_meanfield

!###############################################################################
!> Full ONE-SIDED 2e generalized Fock  G2e_pq = ∂ω_2e/∂U^ispin_pq  in the aligned basis va,vb,
!> via central FD of the EXACT explicit-bra ω_2e (umrsf_omega2e_explicit; the clean bilinear in C —
!> valid under non-orthonormal column variation, unlike the umrsfmntoia back-transform). Column
!> variation dC^ispin_q += th·C^ispin_p (one-sided, NOT symmetric ⇒ the full non-sym gen-Fock that R
!> needs; W later takes ½sym). This realizes the c04/c05 "2e channel-adjoint" factor to FD precision
!> (c05 verified analytic == this FD oracle to ~5e-11). th = 1e-5 (FD truncation ~1e-10; int2 clean,
!> Schwarz off). Cost = 2·nbf² int2 builds/spin — paid once per gradient.
  subroutine umrsf_g2e_onesided(infos, idrv, va, vb, xamp, scale_exch, ispin, g2e)
    use int2_compute, only: int2_compute_t
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), xamp(:), scale_exch
    integer, intent(in) :: ispin
    real(kind=dp), intent(out) :: g2e(:,:)
    real(kind=dp), allocatable :: cwa(:,:), cwb(:,:)
    integer :: nbf, p, q
    real(kind=dp) :: th, omp, omm
    nbf = size(va,1) ; th = 1.0d-5
    allocate(cwa(nbf,nbf), cwb(nbf,nbf))
    g2e = 0.0_dp
    do q = 1, nbf
      do p = 1, nbf
        cwa = va ; cwb = vb
        if (ispin == 1) then ; cwa(:,q) = va(:,q) + th*va(:,p)
        else                 ; cwb(:,q) = vb(:,q) + th*vb(:,p) ; end if
        call umrsf_omega2e_explicit(infos, idrv, cwa, cwb, xamp, scale_exch, omp)
        cwa = va ; cwb = vb
        if (ispin == 1) then ; cwa(:,q) = va(:,q) - th*va(:,p)
        else                 ; cwb(:,q) = vb(:,q) - th*vb(:,p) ; end if
        call umrsf_omega2e_explicit(infos, idrv, cwa, cwb, xamp, scale_exch, omm)
        g2e(p,q) = (omp - omm)/(2.0_dp*th)
      end do
    end do
    deallocate(cwa, cwb)
  end subroutine umrsf_g2e_onesided

!###############################################################################
!> ANALYTIC 11-channel 2e generalized Fock  g2e^σ_pq = ∂ω_2e/∂U^σ_pq  (BOTH spins, one call),
!> replacing the FD oracle umrsf_g2e_onesided (4·nbf² int2 builds → 2 builds). Derivation +
!> CAS/FD gate: DERIVATIONS/{c10_g2e_analytic.py,M4_g2e_analytic.md} (worst 5.25e-12 ≤1e-9).
!>   ω_2e = Σ_k ⟨B_k, F_k⟩,  F_k = s_k int2[D]_k,  D=umrsfcbc, B=umrsf_bra_density,  int2_umrsf
!>   is CHANNEL-DIAGONAL & SELF-ADJOINT (J/K of 8-fold-sym ERIs) ⇒ for the ket-derivative term
!>   ⟨B_k, s_k int2[∂D_k]⟩ = ⟨GB_k, ∂D_k⟩ with GB_k = s_k int2[B]_k (precomputed once).
!>   g2e = V_σᵀ Γ_σ ,  Γ_σ = Σ_k [⟨F_k,∂B_k/∂V_σ⟩ + ⟨GB_k,∂D_k/∂V_σ⟩].
!> Each channel density is V_l A_k V_rᵀ (sparse amplitude core A_k; ket core A, bra core AB;
!> spin pairing aa=1,3,5,7 / bb=2,4,6,8 / ab=9,10,11). Gradient of ⟨M,V_l A V_rᵀ⟩:
!>   ab: Γa += M Vb Aᵀ, Γb += Mᵀ Va A ;  aa: Γa += M Va Aᵀ + Mᵀ Va A ;  bb: analogous.
!> s_k = mrst sign (fmrst2(:,1:10) flip for mrst==3); scale_exchange=scale_coulomb=scale_exch
!> (mirrors umrsf_omega2e_explicit exactly — the FD oracle's integrand).
  subroutine umrsf_g2e_analytic(infos, idrv, va, vb, xamp, scale_exch, g2ea, g2eb)
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc
    use tdhf_lib, only: iatogen
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), xamp(:), scale_exch
    real(kind=dp), intent(out) :: g2ea(:,:), g2eb(:,:)
    real(kind=dp), allocatable :: xmat(:,:), fk(:,:,:), gbk(:,:,:)
    real(kind=dp), allocatable :: acore(:,:,:), abcore(:,:,:), ga(:,:), gb(:,:)
    real(kind=dp), allocatable :: mc(:,:), vacc(:,:)
    real(kind=dp), allocatable, target :: densd(:,:,:,:), densb(:,:,:,:)
    type(int2_umrsf_data_t), target :: udd, udb
    real(kind=dp), parameter :: isqrt2 = 1.0_dp/sqrt(2.0_dp)
    integer :: nbf, nocca, noccb, mrst, o1, o2, i, j, k, pass, ls, rs
    integer :: lspin(11), rspin(11)
    !                   ch:  1  2  3  4  5  6  7  8  9 10 11
    lspin = (/ 1,2,1,2,1,2,1,2,1,1,1 /)   ! 1=α 2=β  (left  spin)
    rspin = (/ 1,2,1,2,1,2,1,2,2,2,2 /)   !          (right spin)

    nbf = size(va,1)
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    mrst = infos%tddft%mult ; o1 = nocca-1 ; o2 = nocca
    allocate(xmat(nbf,nbf), fk(nbf,nbf,11), gbk(nbf,nbf,11), &
             acore(nbf,nbf,11), abcore(nbf,nbf,11), ga(nbf,nbf), gb(nbf,nbf), &
             mc(nbf,nbf), vacc(nbf,nbf), &
             densd(1,11,nbf,nbf), densb(1,11,nbf,nbf), source=0.0_dp)

    ! ---- channel densities D=umrsfcbc(ket), B=umrsf_bra_density(bra) ----
    call iatogen(xamp, xmat, nocca, noccb)
    call umrsfcbc(infos, va, vb, xmat, densd(1,:,:,:))
    call umrsf_bra_density(infos, va, vb, xmat, densb(1,:,:,:))

    ! ---- F_k = s_k int2[D]_k  and  GB_k = s_k int2[B]_k  (2 int2 builds) ----
    udd = int2_umrsf_data_t(d3=densd(1:1,:,:,:), tamm_dancoff=.true., &
                            scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(udd)
    do k = 1, 11 ; fk(:,:,k) = udd%f3(1,k,:,:,1) ; end do
    if (mrst == 3) fk(:,:,1:10) = -fk(:,:,1:10)
    udb = int2_umrsf_data_t(d3=densb(1:1,:,:,:), tamm_dancoff=.true., &
                            scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(udb)
    do k = 1, 11 ; gbk(:,:,k) = udb%f3(1,k,:,:,1) ; end do
    if (mrst == 3) gbk(:,:,1:10) = -gbk(:,:,1:10)

    ! ---- amplitude cores A_k (ket), AB_k (bra) : density_k = V_l A_k V_rᵀ ----
    do j = nocca+1, nbf                                   ! virt-row channels
      acore(o2,j,1)=xmat(o2,j) ; acore(o1,j,3)=xmat(o1,j)               ! aa  bo2va/bo1va
      acore(o2,j,2)=xmat(o2,j) ; acore(o1,j,4)=xmat(o1,j)               ! bb
      acore(o1,j,9)=xmat(o2,j) ; acore(o2,j,9)=-xmat(o1,j)              ! ab  o21v
      abcore(o2,j,5)=0.5_dp*xmat(o2,j) ; abcore(o1,j,7)=0.5_dp*xmat(o1,j)   ! aa  adco1a/adco2a
      abcore(o2,j,6)=0.5_dp*xmat(o2,j) ; abcore(o1,j,8)=0.5_dp*xmat(o1,j)   ! bb
      abcore(o2,j,9)=xmat(o1,j) ; abcore(o1,j,9)=-xmat(o2,j)            ! ab  ao21v
    end do
    do i = 1, nocca-2                                     ! closed-col channels
      acore(i,o1,5)=xmat(i,o1) ; acore(i,o2,7)=xmat(i,o2)              ! aa  bco1a/bco2a
      acore(i,o1,6)=xmat(i,o1) ; acore(i,o2,8)=xmat(i,o2)              ! bb
      acore(i,o2,10)=xmat(i,o1) ; acore(i,o1,10)=-xmat(i,o2)           ! ab  co12
      abcore(i,o1,1)=0.5_dp*xmat(i,o1) ; abcore(i,o2,3)=0.5_dp*xmat(i,o2)   ! aa  ado2va/ado1va
      abcore(i,o1,2)=0.5_dp*xmat(i,o1) ; abcore(i,o2,4)=0.5_dp*xmat(i,o2)   ! bb
      abcore(i,o1,10)=xmat(i,o2) ; abcore(i,o2,10)=-xmat(i,o1)         ! ab  aco12
    end do
    ! channel 11 (ab): SF block; ket A11 excludes the SOMO 2×2, bra AB11 is full then corrected
    do j = noccb+1, nbf ; do i = 1, nocca
      abcore(i,j,11) = xmat(i,j)
      if (.not. ((i==o1 .or. i==o2) .and. (j==o1 .or. j==o2))) acore(i,j,11) = xmat(i,j)
    end do ; end do
    abcore(o1,o1,11) = abcore(o1,o1,11) - xmat(o1,o1)
    abcore(o2,o2,11) = abcore(o2,o2,11) - xmat(o2,o2)
    if (mrst == 1) then
      acore(o2,o1,11)=xmat(o2,o1) ; acore(o1,o2,11)=xmat(o1,o2)
      acore(o1,o1,11)=isqrt2*xmat(o1,o1) ; acore(o2,o2,11)=-isqrt2*xmat(o1,o1)
      abcore(o1,o1,11)=abcore(o1,o1,11)+isqrt2*xmat(o1,o1)
      abcore(o2,o2,11)=abcore(o2,o2,11)-isqrt2*xmat(o1,o1)
    else if (mrst == 3) then
      acore(o1,o1,11)=isqrt2*xmat(o1,o1) ; acore(o2,o2,11)=isqrt2*xmat(o1,o1)
      abcore(o1,o2,11)=abcore(o1,o2,11)-xmat(o1,o2)
      abcore(o2,o1,11)=abcore(o2,o1,11)-xmat(o2,o1)
      abcore(o1,o1,11)=abcore(o1,o1,11)+isqrt2*xmat(o1,o1)
      abcore(o2,o2,11)=abcore(o2,o2,11)+isqrt2*xmat(o1,o1)
    end if

    ! ---- Γa, Γb = Σ_k [bra pass (M=F_k, core=AB_k) + ket pass (M=GB_k, core=A_k)] ----
    ga = 0.0_dp ; gb = 0.0_dp
    do k = 1, 11
      ls = lspin(k) ; rs = rspin(k)
      do pass = 1, 2
        if (pass == 1) then ; mc = fk(:,:,k) ; vacc = abcore(:,:,k)   ! bra
        else                ; mc = gbk(:,:,k) ; vacc = acore(:,:,k)   ! ket
        end if
        if (ls==1 .and. rs==1) then          ! aa
          ga = ga + matmul(mc, matmul(va, transpose(vacc))) &
                  + matmul(transpose(mc), matmul(va, vacc))
        else if (ls==2 .and. rs==2) then     ! bb
          gb = gb + matmul(mc, matmul(vb, transpose(vacc))) &
                  + matmul(transpose(mc), matmul(vb, vacc))
        else                                 ! ab
          ga = ga + matmul(mc, matmul(vb, transpose(vacc)))
          gb = gb + matmul(transpose(mc), matmul(va, vacc))
        end if
      end do
    end do

    g2ea = matmul(transpose(va), ga)
    g2eb = matmul(transpose(vb), gb)
    deallocate(xmat, fk, gbk, acore, abcore, ga, gb, mc, vacc, densd, densb)
  end subroutine umrsf_g2e_analytic

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
!> Spin-flip TDA response matvec  ax = A·xin  in the (aligned) basis vva,vvb, with the MO-basis
!> reference Fock famo/fbmo ALREADY built (= vᵀ fock v with the frozen-core shift, identical to the
!> energy path).  ax = 2e [umrsfcbc → int2_umrsf (J/K, mrst sign, SPC) → umrsfmntoia] + orbital
!> [mrsfesum]. Reused to (i) build the full A column-by-column for the smooth-basis re-diagonalization
!> and (ii) reconstruct ω = xᵀ A x. Mirrors the inline energy-path matvec exactly (clean-room).
  subroutine umrsf_response_Ax(infos, idrv, vva, vvb, famo, fbmo, xin, scale_exch, &
                               hfs, spc_coco, spc_ovov, spc_coov, ax)
    use int2_compute, only: int2_compute_t
    use tdhf_mrsf_lib, only: int2_umrsf_data_t, umrsfcbc, umrsfmntoia, mrsfesum
    use tdhf_lib, only: iatogen
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: vva(:,:), vvb(:,:), famo(:,:), fbmo(:,:), xin(:)
    real(kind=dp), intent(in) :: scale_exch, hfs, spc_coco, spc_ovov, spc_coov
    real(kind=dp), intent(out) :: ax(:)
    real(kind=dp), allocatable :: xmat(:,:), amo(:,:)
    real(kind=dp), allocatable, target :: dens(:,:,:,:)
    real(kind=dp), pointer :: f3(:,:,:,:)
    type(int2_umrsf_data_t), target :: ud
    integer :: nbf, nocca, noccb, nvirb, xvec_dim, mrst

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb ; xvec_dim = nocca*nvirb ; mrst = infos%tddft%mult
    allocate(xmat(nbf,nbf), amo(xvec_dim,1), dens(1,11,nbf,nbf), source=0.0_dp)

    call iatogen(xin, xmat, nocca, noccb)
    call umrsfcbc(infos, vva, vvb, xmat, dens(1,:,:,:))
    ud = int2_umrsf_data_t(d3=dens(1:1,:,:,:), tamm_dancoff=.true., &
                           scale_exchange=scale_exch, scale_coulomb=scale_exch)
    call idrv%run(ud)
    f3 => ud%f3(:,:,:,:,1)
    if (mrst == 3) f3(:,1:10,:,:) = -f3(:,1:10,:,:)
    if (abs(hfs) > epsilon(1.0_dp)) then
      if (spc_coco /= hfs) f3(:,10,:,:) = f3(:,10,:,:) * (spc_coco/hfs)
      if (spc_ovov /= hfs) f3(:,9,:,:)  = f3(:,9,:,:)  * (spc_ovov/hfs)
      if (spc_coov /= hfs) f3(:,1:8,:,:) = f3(:,1:8,:,:) * (spc_coov/hfs)
    end if
    amo = 0.0_dp
    call umrsfmntoia(infos, f3(1,:,:,:), amo, vva, vvb, 1)
    call iatogen(xin, xmat, nocca, noccb)
    call mrsfesum(infos, xmat, famo, fbmo, amo, 1)
    ax = amo(:,1)
    deallocate(xmat, amo, dens)
  end subroutine umrsf_response_Ax

!###############################################################################
!> Re-diagonalize the spin-flip TDA response matrix A in the SMOOTH-aligned basis and overlap-track
!> to the stored amplitude bvec_ref → the genuine eigenvector xamp (+ eigenvalue omega_eig).
!> RULES §15 / c05: the ov-only Z-vector + W machinery is exact ONLY when X is a TRUE eigenvector of
!> A in the alignment basis. The stored bvec is the eigenvector in the ENERGY's THRESHOLD basis and is
!> non-stationary in the converged-smooth basis (Rayleigh quotient ~1.5e-6 above the eigenvalue),
!> which would corrupt R. Building A column-by-column via the energy matvec (nia int2 builds, cheap)
!> and diagonalizing recovers the stationary X. State-following by max |eigenvector·bvec_ref| overlap
!> (RULES §11), sign-fixed to bvec_ref.
  subroutine umrsf_track_amplitude(infos, idrv, va, vb, famo, fbmo, bvec_ref, scale_exch, &
                                   hfs, spc_coco, spc_ovov, spc_coov, xamp, omega_eig)
    use eigen, only: diag_symm_full
    use int2_compute, only: int2_compute_t
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), famo(:,:), fbmo(:,:), bvec_ref(:)
    real(kind=dp), intent(in) :: scale_exch, hfs, spc_coco, spc_ovov, spc_coov
    real(kind=dp), intent(out) :: xamp(:), omega_eig
    real(kind=dp), allocatable :: amat(:,:), ek(:), axk(:), ev(:)
    integer :: nbf, nocca, noccb, nvirb, nia, k, ktrack, ierr
    real(kind=dp) :: ov, ovmax, sgn

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb ; nia = nocca*nvirb
    allocate(amat(nia,nia), ek(nia), axk(nia), ev(nia), source=0.0_dp)

    do k = 1, nia
      ek = 0.0_dp ; ek(k) = 1.0_dp
      call umrsf_response_Ax(infos, idrv, va, vb, famo, fbmo, ek, scale_exch, &
                             hfs, spc_coco, spc_ovov, spc_coov, axk)
      amat(:,k) = axk
    end do
    block ; use io_constants, only: iw
      write(iw,'(2x,a,es12.3)') 'response A asymmetry ||A−Aᵀ|| (SOMO-overwrite probe) = ', &
        maxval(abs(amat - transpose(amat)))
    end block
    amat = 0.5_dp*(amat + transpose(amat))     ! symmetrize (TDA A is symmetric; kills matvec noise)
    call diag_symm_full(1, nia, amat, nia, ev, ierr)   ! columns of amat → eigenvectors

    ! state-following: pick the eigenvector with max |overlap| to the stored bvec_ref; sign-fix to it
    ovmax = -1.0_dp ; ktrack = 1
    do k = 1, nia
      ov = abs(dot_product(amat(:,k), bvec_ref))
      if (ov > ovmax) then ; ovmax = ov ; ktrack = k ; end if
    end do
    sgn = sign(1.0_dp, dot_product(amat(:,ktrack), bvec_ref))
    xamp = sgn*amat(:,ktrack)
    omega_eig = ev(ktrack)
    deallocate(amat, ek, axk, ev)
  end subroutine umrsf_track_amplitude

!###############################################################################
!> MATRIX-FREE Davidson replacement for umrsf_track_amplitude — the perf fix for the thymine wall.
!> The dense path builds A by nia=nocca*nvirb umrsf_response_Ax calls + a full diag_symm_full (nia~3910
!> at thymine, the wall). Here the SAME response A is diagonalized matrix-free: seed the subspace with
!> bvec_ref (the energy eigenvector, ~1.5e-6 from the smooth-basis one) and converge the SINGLE nearest
!> root by Davidson, applying A via ONE umrsf_response_Ax per new subspace vector (~tens, not nia).
!> STATE-FOLLOW: each iteration picks the Ritz pair whose Ritz vector has max |overlap with bvec_ref|
!> (= max |g.Y(:,j)|, g = V^T bvec_ref) — the SAME criterion as the dense path (lines ~1788-1794),
!> restricted to the converged subspace. Returns the genuine smooth-basis eigenvector xamp (sign-fixed
!> to bvec_ref) + eigenvalue omega_eig, reproducing the dense result to the Davidson tolerance. Jacobi
!> preconditioner = (theta - (eps_b - eps_i))^-1 (orbital gap; DOF k packed column-major over
!> i=1..nocca, a=noccb+1..nbf, per iatogen). Tunables: UMRSF_TRKTOL (residual norm, 1e-10),
!> UMRSF_TRKMAXSUB (max subspace before thick-restart, default min(nia,60)).
  subroutine umrsf_track_amplitude_dav(infos, idrv, va, vb, famo, fbmo, bvec_ref, scale_exch, &
                                       hfs, spc_coco, spc_ovov, spc_coov, xamp, omega_eig)
    use io_constants, only: iw
    use eigen, only: diag_symm_full
    use int2_compute, only: int2_compute_t
    implicit none
    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: idrv
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), famo(:,:), fbmo(:,:), bvec_ref(:)
    real(kind=dp), intent(in) :: scale_exch, hfs, spc_coco, spc_ovov, spc_coov
    real(kind=dp), intent(out) :: xamp(:), omega_eig
    real(kind=dp), allocatable :: vsub(:,:), avsub(:,:), hsub(:,:), hcopy(:,:), theta(:), yy(:,:)
    real(kind=dp), allocatable :: adiag(:), xr(:), axr(:), rr(:), tt(:), gg(:)
    integer :: nbf, nocca, noccb, nvirb, nia, i, a, k, j, m, mmax, jsel, ierr, iter, maxit
    real(kind=dp) :: rtol, nrm, rnorm, theta_sel, ov, ovmax, denom, sgn
    real(kind=dp), parameter :: pcfloor = 1.0e-3_dp

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvirb = nbf - noccb ; nia = nocca*nvirb

    ! orbital-gap diagonal of A (preconditioner); DOF k packed column-major over (i=1..nocca, a=noccb+1..nbf)
    allocate(adiag(nia))
    k = 0
    do a = noccb+1, nbf ; do i = 1, nocca
      k = k + 1 ; adiag(k) = fbmo(a,a) - famo(i,i)
    end do ; end do

    rtol = 1.0e-10_dp ; mmax = min(nia, 60) ; maxit = 300
    block
      character(len=24) :: e ; integer :: ios ; real(kind=dp) :: rv ; integer :: iv
      call get_environment_variable("UMRSF_TRKTOL", e, status=ios)
      if (ios==0) then ; read(e,*,iostat=ios) rv ; if (ios==0 .and. rv>0.0_dp) rtol = rv ; end if
      call get_environment_variable("UMRSF_TRKMAXSUB", e, status=ios)
      if (ios==0) then ; read(e,*,iostat=ios) iv ; if (ios==0 .and. iv>1) mmax = min(nia, iv) ; end if
    end block

    allocate(vsub(nia,mmax), avsub(nia,mmax), hsub(mmax,mmax), hcopy(mmax,mmax), theta(mmax), &
             yy(mmax,mmax), xr(nia), axr(nia), rr(nia), tt(nia), gg(mmax), source=0.0_dp)

    ! seed the subspace with the (normalized) stored amplitude
    nrm = sqrt(sum(bvec_ref**2)) ; vsub(:,1) = bvec_ref / nrm
    call umrsf_response_Ax(infos, idrv, va, vb, famo, fbmo, vsub(:,1), scale_exch, &
                           hfs, spc_coco, spc_ovov, spc_coov, avsub(:,1))
    m = 1 ; theta_sel = 0.0_dp ; rnorm = huge(1.0_dp) ; iter = 0
    davidson: do iter = 1, maxit
      ! subspace matrix H = V^T A V (symmetric)
      do j = 1, m ; do i = 1, m ; hsub(i,j) = dot_product(vsub(:,i), avsub(:,j)) ; end do ; end do
      hcopy = 0.0_dp ; hcopy(1:m,1:m) = 0.5_dp*(hsub(1:m,1:m) + transpose(hsub(1:m,1:m)))
      call diag_symm_full(1, m, hcopy, mmax, theta, ierr)     ! cols of hcopy(1:m,1:m) -> eigvecs
      yy(1:m,1:m) = hcopy(1:m,1:m)
      ! state-follow: ritz pair with max |overlap to bvec_ref| = max |g.Y(:,j)|, g = V^T bvec_ref
      do i = 1, m ; gg(i) = dot_product(vsub(:,i), bvec_ref) ; end do
      ovmax = -1.0_dp ; jsel = 1
      do j = 1, m
        ov = abs(dot_product(gg(1:m), yy(1:m,j)))
        if (ov > ovmax) then ; ovmax = ov ; jsel = j ; end if
      end do
      theta_sel = theta(jsel)
      ! Ritz vector x = V Y(:,jsel) ; A x = AV Y(:,jsel) ; residual r = A x - theta x
      xr = 0.0_dp ; axr = 0.0_dp
      do i = 1, m ; xr = xr + yy(i,jsel)*vsub(:,i) ; axr = axr + yy(i,jsel)*avsub(:,i) ; end do
      rr = axr - theta_sel*xr ; rnorm = sqrt(sum(rr**2))
      if (rnorm <= rtol) exit davidson
      ! thick-restart: subspace full -> collapse to the tracked Ritz vector
      if (m == mmax) then ; vsub(:,1) = xr ; avsub(:,1) = axr ; m = 1 ; end if
      ! preconditioned correction t = r / (theta - adiag) (floored)
      do k = 1, nia
        denom = theta_sel - adiag(k) ; if (abs(denom) < pcfloor) denom = sign(pcfloor, denom)
        tt(k) = rr(k) / denom
      end do
      ! orthonormalize t against V(:,1:m) (modified Gram-Schmidt, twice)
      do j = 1, 2 ; do i = 1, m ; tt = tt - dot_product(tt, vsub(:,i))*vsub(:,i) ; end do ; end do
      nrm = sqrt(sum(tt**2))
      if (nrm < 1.0e-12_dp) exit davidson     ! correction exhausted -> converged in subspace
      tt = tt / nrm ; m = m + 1 ; vsub(:,m) = tt
      call umrsf_response_Ax(infos, idrv, va, vb, famo, fbmo, vsub(:,m), scale_exch, &
                             hfs, spc_coco, spc_ovov, spc_coov, avsub(:,m))
    end do davidson

    ! sign-fix to bvec_ref (matches the dense path)
    sgn = sign(1.0_dp, dot_product(xr, bvec_ref)) ; xamp = sgn*xr ; omega_eig = theta_sel
    write(iw,'(2x,a,i0,a,es10.2,a,es10.2)') 'track Davidson: iters = ', iter, &
      ', resid = ', rnorm, ', overlap deficit 1-|x.bvec| = ', 1.0_dp - abs(dot_product(xamp, bvec_ref))
    deallocate(vsub, avsub, hsub, hcopy, theta, yy, adiag, xr, axr, rr, tt, gg)
  end subroutine umrsf_track_amplitude_dav

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
!> Z-vector via the ANALYTIC RHS (§15 / c05): R^σ = (G^f_σ)_ov antisym in the CANONICAL SCF basis,
!> with G^f_σ = V_σ G̃_σ V_σᵀ (the M1 V-transform of the aligned-basis raw generalized Fock G̃).
!> Solve M z = −R via the validated UHF spin-coupled CPHF Hessian (cphf_solve_uhf), feeding the
!> CANONICAL MOs/energies through tagarray transiently (cphf reads them from there). Ra,Rb are the
!> ov blocks (Ra(iocc,ivir) = G^f_{ia} − G^f_{ai}; the ∂ω/∂κ convention = the old numerical RHS, so
!> the cphf call passes −R exactly as the validated numerical path did). Returns AO relaxation
!> densities pza,pzb (symmetric) + the MO relaxation matrices zmata,zmatb (canonical, ov+vo).
  subroutine umrsf_zvector_analytic(infos, cac, cbc, epsca, epscb, ra, rb, pza, pzb, zmata, zmatb, zrms)
    use oqp_tagarray_driver
    use cphf_mod, only: cphf_solve_uhf
    implicit none
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), epsca(:), epscb(:)
    real(kind=dp), intent(in) :: ra(:,:), rb(:,:)      ! ov blocks (nocca,nvira) / (noccb,nvirb)
    real(kind=dp), intent(out) :: pza(:,:), pzb(:,:), zmata(:,:), zmatb(:,:), zrms
    real(kind=dp), contiguous, pointer :: moa(:,:), mob(:,:), ea(:), eb(:)
    real(kind=dp), allocatable :: rhs(:,:), zsol(:,:), moa_s(:,:), mob_s(:,:), ea_s(:), eb_s(:)
    integer :: nbf, nocca, noccb, nvira, nvirb, la, lb, ltot, ivir, iocc, iov, amo
    integer(4) :: status

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    nvira = nbf - nocca ; nvirb = nbf - noccb
    la = nocca*nvira ; lb = noccb*nvirb ; ltot = la + lb
    allocate(rhs(ltot,1), zsol(ltot,1), source=0.0_dp)

    ! pack R into the cphf occ-major layout  iov = base + (ivir-1)*nocc + iocc
    do ivir = 1, nvira ; do iocc = 1, nocca
      rhs((ivir-1)*nocca + iocc, 1) = ra(iocc, ivir)
    end do ; end do
    do ivir = 1, nvirb ; do iocc = 1, noccb
      rhs(la + (ivir-1)*noccb + iocc, 1) = rb(iocc, ivir)
    end do ; end do

    ! solve M z = −R with the canonical MOs transiently in tagarray (cphf reads them there)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, moa, status)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mob, status)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, ea, status)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, eb, status)
    allocate(moa_s, source=moa) ; allocate(mob_s, source=mob)
    allocate(ea_s, source=ea)   ; allocate(eb_s, source=eb)
    moa = cac ; mob = cbc ; ea = epsca ; eb = epscb
    call cphf_solve_uhf(infos, 1, -rhs, zsol, tol=1.0d-10)
    moa = moa_s ; mob = mob_s ; ea = ea_s ; eb = eb_s
    zrms = sqrt(sum(zsol(:,1)**2)/max(ltot,1))

    ! AO relaxation densities P_z^σ = Σ_ia z^σ_ia (C_i C_a^T + C_a C_i^T) + MO zmat (ov+vo)
    pza = 0.0_dp ; pzb = 0.0_dp ; zmata = 0.0_dp ; zmatb = 0.0_dp
    do ivir = 1, nvira ; amo = nocca + ivir ; do iocc = 1, nocca
      iov = (ivir-1)*nocca + iocc
      call add_outer(pza, zsol(iov,1), cac(:,iocc), cac(:,amo))
      call add_outer(pza, zsol(iov,1), cac(:,amo), cac(:,iocc))
      zmata(iocc,amo) = zsol(iov,1) ; zmata(amo,iocc) = zsol(iov,1)
    end do ; end do
    do ivir = 1, nvirb ; amo = noccb + ivir ; do iocc = 1, noccb
      iov = la + (ivir-1)*noccb + iocc
      call add_outer(pzb, zsol(iov,1), cbc(:,iocc), cbc(:,amo))
      call add_outer(pzb, zsol(iov,1), cbc(:,amo), cbc(:,iocc))
      zmatb(iocc,amo) = zsol(iov,1) ; zmatb(amo,iocc) = zsol(iov,1)
    end do ; end do
    deallocate(rhs, zsol, moa_s, mob_s, ea_s, eb_s)
  end subroutine umrsf_zvector_analytic

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
!> Analytic z-coupling generalized Fock G^z (the orbital-Hessian ACTION, also the final G^z for W):
!>   G^z_σ,pq = ε^σ_p z^σ_pq + [q ≤ nocc_σ] (C_σᵀ G_σ[P_z] C_σ)_pq ,  P_z,σ = C_σ z^σ C_σᵀ (AO),
!>   G_σ[P_z] = J[P_z,α+P_z,β] − hfscale·K[P_z,σ]   (ONE umrsf_meanfield / fock_jk build of P_z).
!> This is EXACTLY the gen-Fock of the z-coupling Σ_pq z_pq F^ref_pq(C) for ANY symmetric z (ov OR
!> full-block oo+ov+vv): the frozen part ε_p z_pq + the reference-density-relaxation mean field, the
!> latter nonzero only for occupied columns q (rotating a virtual column does not change D^ref). The
!> 2(Cᵀ G[½P_z] C) form used inside umrsf_meanfield(½P_z) ×2 = (Cᵀ G[P_z] C) by linearity. zmata/zmatb
!> are symmetric MO matrices in the canonical basis cac/cbc. Used per trial z to build the full-block
!> Hessian M, and once at the solved z to assemble G^z for the energy-weighted density W.
  subroutine umrsf_genfock_z(infos, basis, cac, cbc, epsca, epscb, zmata, zmatb, hfscale_ref, gza, gzb)
    implicit none
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), epsca(:), epscb(:)
    real(kind=dp), intent(in) :: zmata(:,:), zmatb(:,:), hfscale_ref
    real(kind=dp), intent(out) :: gza(:,:), gzb(:,:)
    real(kind=dp), allocatable :: pza(:,:), pzb(:,:), ya(:,:), yb(:,:), tmp(:,:)
    integer :: nbf, nocca, noccb, pp, qq

    nbf = basis%nbf ; nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b
    allocate(pza(nbf,nbf), pzb(nbf,nbf), ya(nbf,nbf), yb(nbf,nbf), tmp(nbf,nbf), source=0.0_dp)
    ! AO relaxation density P_z,σ = C_σ z^σ C_σᵀ  (z symmetric off-diagonal in the MO basis)
    pza = matmul(cac, matmul(zmata, transpose(cac)))
    pzb = matmul(cbc, matmul(zmatb, transpose(cbc)))
    ! reference mean field G_σ[P_z] (one fock_jk build); 2·G[½P_z] = G[P_z]
    call umrsf_meanfield(basis, infos, 0.5_dp*pza, 0.5_dp*pzb, hfscale_ref, ya, yb)
    ! frozen orbital-energy part ε_p z_pq (canonical F^MO = diag ε)
    do qq = 1, nbf ; do pp = 1, nbf
      gza(pp,qq) = epsca(pp)*zmata(pp,qq)
      gzb(pp,qq) = epscb(pp)*zmatb(pp,qq)
    end do ; end do
    ! reference-density-relaxation mean field, occupied columns only
    tmp = matmul(transpose(cac), matmul(ya, cac)) ; gza(:,1:nocca) = gza(:,1:nocca) + 2.0_dp*tmp(:,1:nocca)
    tmp = matmul(transpose(cbc), matmul(yb, cbc)) ; gzb(:,1:noccb) = gzb(:,1:noccb) + 2.0_dp*tmp(:,1:noccb)
    deallocate(pza, pzb, ya, yb, tmp)
  end subroutine umrsf_genfock_z

!###############################################################################
!> FULL generalized Fock G^f of f = ω∘align at the CANONICAL orbitals cac/cbc — it CARRIES the
!> alignment Jacobian dV/dC (the M1 term the V-transform G^f = V G̃ Vᵀ drops). c06 §16: the real
!> 11-channel ω is NOT within-segment invariant, so ΔG^f = G^f − V G̃ Vᵀ ≠ 0 (concentrated in the
!> oo-α / vv-β antisym blocks) and is REQUIRED in both the Z-vector RHS and W (dropping it leaves the
!> ~1e-3 wall). Computed by the re-alignment response (c06 Gf_fdalign / the FD-oracle, CLEAN — no
!> int2/fock in the loop): for each canonical rotation (p,q,spin), one-sided perturb cac/cbc, re-align
!> with the SMOOTH (converged) get_jacobi + sign-fix to va/vb, form the induced ALIGNED-basis rotation
!> dU_σ = C̃_σᵀ S dC̃_σ, and contract with the raw aligned gen-Fock G̃ (gta/gtb):
!>   G^f_pq = Σ_rs G̃a_rs dUa_rs + Σ_rs G̃b_rs dUb_rs       (both spins respond to one spin's rotation).
!> Central difference (th=1e-4). Returns G^f in the CANONICAL (cac/cbc) basis. (In-model: unseeded
!> smooth re-align reproduces the FD-oracle G^f and closes the gradient to ~5e-9.)
  subroutine umrsf_genfock_full(infos, cac, cbc, va, vb, smat_full, gta, gtb, nocca, gfa, gfb)
    implicit none
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), va(:,:), vb(:,:), smat_full(:,:)
    real(kind=dp), intent(in) :: gta(:,:), gtb(:,:)
    integer, intent(in) :: nocca
    real(kind=dp), intent(out) :: gfa(:,:), gfb(:,:)
    real(kind=dp), allocatable :: cwa(:,:), cwb(:,:), ctap(:,:), ctbp(:,:), dua(:,:), dub(:,:)
    integer :: nbf, p, q, isp
    real(kind=dp) :: th, dum

    nbf = size(cac,1) ; th = 1.0d-4
    block   ! UMRSF_TH: re-align FD step for G^f (diagnostic — S2's large within-seg antisym stresses it)
      character(len=24) :: e ; integer :: ios ; real(kind=dp) :: thv
      call get_environment_variable("UMRSF_TH", e, status=ios)
      if (ios==0) then ; read(e,*,iostat=ios) thv ; if (ios==0 .and. thv>0.0_dp) th = thv ; end if
    end block
    allocate(cwa(nbf,nbf), cwb(nbf,nbf), ctap(nbf,nbf), ctbp(nbf,nbf), dua(nbf,nbf), dub(nbf,nbf), &
             source=0.0_dp)
    gfa = 0.0_dp ; gfb = 0.0_dp
    do isp = 1, 2
      do q = 1, nbf
        do p = 1, nbf
          ! +th : one-sided perturbation of canonical column q in the direction of column p
          cwa = cac ; cwb = cbc
          if (isp == 1) then ; cwa(:,q) = cac(:,q) + th*cac(:,p)
          else               ; cwb(:,q) = cbc(:,q) + th*cbc(:,p) ; end if
          call umrsf_jacobi_smooth(cwa, cwb, smat_full, nocca, dum)
          call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
          ctap = cwa ; ctbp = cwb
          ! -th
          cwa = cac ; cwb = cbc
          if (isp == 1) then ; cwa(:,q) = cac(:,q) - th*cac(:,p)
          else               ; cwb(:,q) = cbc(:,q) - th*cbc(:,p) ; end if
          call umrsf_jacobi_smooth(cwa, cwb, smat_full, nocca, dum)
          call m1_sign_fix(cwa, cwb, va, vb, smat_full, nbf)
          ! induced aligned-basis rotation dU_σ = C̃_σᵀ S dC̃_σ (base aligned C̃ = va/vb)
          dua = matmul(transpose(va), matmul(smat_full, ctap - cwa)) / (2.0_dp*th)
          dub = matmul(transpose(vb), matmul(smat_full, ctbp - cwb)) / (2.0_dp*th)
          if (isp == 1) then ; gfa(p,q) = sum(gta*dua) + sum(gtb*dub)
          else               ; gfb(p,q) = sum(gta*dua) + sum(gtb*dub) ; end if
        end do
      end do
    end do
    deallocate(cwa, cwb, ctap, ctbp, dua, dub)
  end subroutine umrsf_genfock_full

!###############################################################################
!> FULL-BLOCK Z-vector (c06 §16, the decisive fix): solve M z = −R over ALL off-diagonal orbital
!> rotations p>q (oo + ov + vv, both spins), not just ov. The aligned 11-channel ω is NOT stationary
!> to oo/vv rotations (canonical-C G^f antisym: α-oo, β-vv ≠ 0), so the standard ov-only Z-vector
!> leaves the ~1e-3 wall. R = antisym(G^f) over the p>q pairs. M = the spin-coupled orbital Hessian =
!> the antisym of the z-coupling gen-Fock (umrsf_genfock_z), built DENSE column-by-column and solved
!> by dgelss (SVD least-squares: M is indefinite over oo/vv, and RANK-DEFICIENT for SOMO-degenerate
!> systems like linear molecules' π_x/π_y ⇒ min-norm; reduces to the exact dgesv solve when full-rank).
!> ovonly=.true. restricts the DOFs to ov pairs
!> (the ablation that reproduces the wall). Returns AO relaxation densities pza/pzb and symmetric MO
!> z-matrices zmata/zmatb (canonical basis), zrms, and statio = ||antisym(G^f+G^z(z))|| (the solve
!> residual; → 0 confirms M z = −R, i.e. the full-block stationarity the ov-only Z cannot reach).
  subroutine umrsf_zvector_fullblock(infos, basis, cac, cbc, epsca, epscb, gfa, gfb, &
                                     hfscale_ref, ovonly, pza, pzb, zmata, zmatb, zrms, statio)
    implicit none
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(in) :: cac(:,:), cbc(:,:), epsca(:), epscb(:), gfa(:,:), gfb(:,:)
    real(kind=dp), intent(in) :: hfscale_ref
    logical, intent(in) :: ovonly
    real(kind=dp), intent(out) :: pza(:,:), pzb(:,:), zmata(:,:), zmatb(:,:), zrms, statio
    integer, allocatable :: dsp(:), dpr(:), dqr(:)
    real(kind=dp), allocatable :: mmat(:,:), rhs(:,:), za1(:,:), zb1(:,:), gza(:,:), gzb(:,:)
    integer :: nbf, nocca, noccb, ndof, k, j, p, q, sp, nocc, info
    logical :: keep
    real(kind=dp) :: rr
    external :: dgelss     ! ILP64 LAPACK: OQP_BLAS_INT=8 ⇒ default integer matches; call directly (cf. resp.F90)

    nbf = basis%nbf ; nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b

    ! ---- enumerate the rotation DOFs (sp, p>q): full-block, or ov-only (p virt, q occ) if ovonly ----
    ndof = 0
    do sp = 1, 2
      nocc = nocca ; if (sp == 2) nocc = noccb
      do p = 1, nbf ; do q = 1, p-1
        keep = .true. ; if (ovonly) keep = (p > nocc .and. q <= nocc)
        if (keep) ndof = ndof + 1
      end do ; end do
    end do
    allocate(dsp(ndof), dpr(ndof), dqr(ndof), source=0)
    allocate(mmat(ndof,ndof), rhs(ndof,1), za1(nbf,nbf), zb1(nbf,nbf), &
             gza(nbf,nbf), gzb(nbf,nbf), source=0.0_dp)
    k = 0
    do sp = 1, 2
      nocc = nocca ; if (sp == 2) nocc = noccb
      do p = 1, nbf ; do q = 1, p-1
        keep = .true. ; if (ovonly) keep = (p > nocc .and. q <= nocc)
        if (keep) then ; k = k + 1 ; dsp(k) = sp ; dpr(k) = p ; dqr(k) = q ; end if
      end do ; end do
    end do

    ! ---- RHS  R_k = antisym(G^f)  ;  rhs = −R ----
    do k = 1, ndof
      if (dsp(k) == 1) then ; rhs(k,1) = -(gfa(dpr(k),dqr(k)) - gfa(dqr(k),dpr(k)))
      else                  ; rhs(k,1) = -(gfb(dpr(k),dqr(k)) - gfb(dqr(k),dpr(k))) ; end if
    end do

    ! ---- dense Hessian M: column k = antisym of the z-coupling gen-Fock for unit z at DOF k ----
    do k = 1, ndof
      za1 = 0.0_dp ; zb1 = 0.0_dp
      if (dsp(k) == 1) then ; za1(dpr(k),dqr(k)) = 1.0_dp ; za1(dqr(k),dpr(k)) = 1.0_dp
      else                  ; zb1(dpr(k),dqr(k)) = 1.0_dp ; zb1(dqr(k),dpr(k)) = 1.0_dp ; end if
      call umrsf_genfock_z(infos, basis, cac, cbc, epsca, epscb, za1, zb1, hfscale_ref, gza, gzb)
      do j = 1, ndof
        if (dsp(j) == 1) then ; mmat(j,k) = gza(dpr(j),dqr(j)) - gza(dqr(j),dpr(j))
        else                  ; mmat(j,k) = gzb(dpr(j),dqr(j)) - gzb(dqr(j),dpr(j)) ; end if
      end do
    end do

    ! ---- solve M z = −R via SVD least-squares (dgelss): rank-revealing + minimum-norm ----
    ! M is indefinite over oo/vv (so not CG) and can be RANK-DEFICIENT for SOMO-degenerate systems
    ! (e.g. linear molecules whose two SOMOs are the degenerate π_x/π_y pair: rotating one into the
    ! other is a zero-energy mode ⇒ a null direction; a plain LU/dgesv blows up). The null space is
    ! pure symmetry gauge (R has no component along it), so the min-norm SVD solution is the physical
    ! relaxation. For full-rank M (CH2, butadiene, …) dgelss returns the exact dgesv solution.
    block
      use io_constants, only: iw
      real(kind=dp), allocatable :: svals(:), work(:)
      real(kind=dp) :: wq(1), rcond
      integer :: rank, lwork
      allocate(svals(ndof))
      rcond = 1.0e-9_dp                                   ! drop σ ≤ rcond·σ_max (the degenerate null space)
      call dgelss(ndof, ndof, 1, mmat, ndof, rhs, ndof, svals, rcond, rank, wq, -1, info)
      lwork = max(int(wq(1)), 1) ; allocate(work(lwork))
      call dgelss(ndof, ndof, 1, mmat, ndof, rhs, ndof, svals, rcond, rank, work, lwork, info)
      if (info /= 0) write(iw,'(2x,a,i0)') 'umrsf_zvector_fullblock: dgelss info = ', info
      if (rank < ndof) write(iw,'(2x,a,i0,a,i0,a,es10.2)') &
        'umrsf_zvector_fullblock: M rank-deficient (SOMO degeneracy) rank ', rank, ' / ', ndof, &
        '; σ_min(kept)/σ_max = ', svals(rank)/max(svals(1), tiny(1.0_dp))
      deallocate(svals, work)
    end block

    ! ---- unpack z → symmetric MO matrices + AO relaxation densities ----
    zmata = 0.0_dp ; zmatb = 0.0_dp
    do k = 1, ndof
      if (dsp(k) == 1) then
        zmata(dpr(k),dqr(k)) = rhs(k,1) ; zmata(dqr(k),dpr(k)) = rhs(k,1)
      else
        zmatb(dpr(k),dqr(k)) = rhs(k,1) ; zmatb(dqr(k),dpr(k)) = rhs(k,1)
      end if
    end do
    pza = matmul(cac, matmul(zmata, transpose(cac)))
    pzb = matmul(cbc, matmul(zmatb, transpose(cbc)))
    zrms = sqrt(sum(rhs(:,1)**2)/max(ndof,1))

    ! ---- stationarity residual ||antisym(G^f + G^z(z))|| over the DOFs (→ 0 ⇒ M z = −R solved) ----
    call umrsf_genfock_z(infos, basis, cac, cbc, epsca, epscb, zmata, zmatb, hfscale_ref, gza, gzb)
    statio = 0.0_dp
    do k = 1, ndof
      if (dsp(k) == 1) then
        rr = (gfa(dpr(k),dqr(k)) - gfa(dqr(k),dpr(k))) + (gza(dpr(k),dqr(k)) - gza(dqr(k),dpr(k)))
      else
        rr = (gfb(dpr(k),dqr(k)) - gfb(dqr(k),dpr(k))) + (gzb(dpr(k),dqr(k)) - gzb(dqr(k),dpr(k)))
      end if
      statio = max(statio, abs(rr))
    end do

    deallocate(dsp, dpr, dqr, mmat, rhs, za1, zb1, gza, gzb)
  end subroutine umrsf_zvector_fullblock

!###############################################################################
!> MATRIX-FREE iterative full-block Z-vector — the perf replacement for umrsf_zvector_fullblock,
!> REUSING OQP's stock matrix-free linear solvers (source/pcg.F90, source/minres.F90) instead of a
!> bespoke GMRES. Solves the IDENTICAL system M z = -R (M = the spin-coupled orbital Hessian = antisym
!> of the z-coupling gen-Fock umrsf_genfock_z over the p>q DOFs; R = antisym(G^f)) but WITHOUT forming
!> M densely — the dense path needs ndof (~5112 at nbf=72, ~21462 at thymine) Fock builds, ONE per
!> column. M z is applied MATRIX-FREE via a single umrsf_genfock_z (one fock_jk + one f_xc grid pass
!> for DFT).
!>   STRUCTURE — M is BLOCK LOWER-TRIANGULAR. tmp = C^T G[1/2 P_z] C is symmetric (ya symmetric, same
!>   C both sides), so the refrelax (occ-column) part contributes ZERO antisym on the oo/vv readouts
!>   => the oo/vv DOF rows of M are PURE DIAGONAL eps_p-eps_q (M_{D,D}=diag, M_{D,V}=0, D=oo+vv); the
!>   ov rows DO couple to D via the density (M_{V,D}!=0) and the ov-ov block A_{V,V} is SYMMETRIC (but
!>   INDEFINITE: the excited-state/SOMO orbital Hessian has negative eigenvalues). So the solve splits
!>   EXACTLY:
!>     (1) z_D = (eps_p-eps_q)^-1 b_D                     [diagonal divide, no matvec]
!>     (2) A_{V,V} z_V = b_V - M_{V,D} z_D                [sym. indefinite => minres_optimize]
!>   where M_{V,D} z_D = the ov-antisym readout of ONE umrsf_genfock_z applied to z_D. (2) is solved by
!>   minres_optimize (matvec umrsf_zov_matvec = the ov-restricted umrsf_genfock_z, preconditioner
!>   umrsf_zov_precond = the floored (eps_a-eps_i)^-1 Jacobi diagonal) — MINRES (Paige-Saunders) is
!>   residual-minimizing and stable on indefinite A, converging like the old GMRES (~16 matvecs);
!>   pcg_optimize (UMRSF_ZPCG=1) STALLS here (CG needs SPD) and is kept only for a verified-SPD
!>   reference. Both callbacks read a umrsf_zov_ctx_t passed through the solvers' c_ptr `dat`.
!> Returns IDENTICAL outputs to umrsf_zvector_fullblock (pza/pzb, zmata/zmatb, zrms, statio =
!> ||antisym(G^f+G^z(z))|| = the M z = -R residual). For the SOMO-degenerate rank-deficient case
!> (linear diradicals: eps_p-eps_q -> 0 on a symmetry pair) the dense SVD min-norm path (UMRSF_ZDENSE=1)
!> is still preferred. Tunables: UMRSF_ZTOL (relative residual, 1e-11), UMRSF_ZMAXIT (max iters),
!> UMRSF_ZPCG (1 => pcg instead of minres; only valid if A_{V,V} is SPD).
  subroutine umrsf_zvector_iter(infos, basis, cac, cbc, epsca, epscb, gfa, gfb, &
                                hfscale_ref, ovonly, pza, pzb, zmata, zmatb, zrms, statio)
    use io_constants, only: iw
    use zvector_common, only: sanitize_zvector_preconditioner
    use pcg_mod, only: pcg_optimize
    use minres_mod, only: minres_optimize
    implicit none
    type(information), target, intent(inout) :: infos
    type(basis_set), target, intent(inout) :: basis
    real(kind=dp), target, intent(in) :: cac(:,:), cbc(:,:), epsca(:), epscb(:)
    real(kind=dp), intent(in) :: gfa(:,:), gfb(:,:)
    real(kind=dp), intent(in) :: hfscale_ref
    logical, intent(in) :: ovonly
    real(kind=dp), intent(out) :: pza(:,:), pzb(:,:), zmata(:,:), zmatb(:,:), zrms, statio
    type(umrsf_zov_ctx_t), target :: ctx
    integer, allocatable :: dsp(:), dpr(:), dqr(:)
    logical, allocatable :: isov(:)
    real(kind=dp), allocatable :: bvec(:), zvec(:), diagm(:), pcinv(:), rhsov(:)
    real(kind=dp), allocatable :: za1(:,:), zb1(:,:), gza(:,:), gzb(:,:)
    integer :: nbf, nocca, noccb, ndof, ndofov, k, p, q, sp, nocc, kk, mxit, iters
    logical :: keep, use_minres
    real(kind=dp) :: rtol, bnorm, errout, cgit, relres, rr
    character(len=8) :: sname

    nbf = basis%nbf ; nocca = infos%mol_prop%nelec_a ; noccb = infos%mol_prop%nelec_b

    ! ---- enumerate ALL p>q DOFs (full-block, or ov-only if ovonly) + mark the ov subset ----
    ndof = 0
    do sp = 1, 2
      nocc = nocca ; if (sp == 2) nocc = noccb
      do p = 1, nbf ; do q = 1, p-1
        keep = .true. ; if (ovonly) keep = (p > nocc .and. q <= nocc)
        if (keep) ndof = ndof + 1
      end do ; end do
    end do
    allocate(dsp(ndof), dpr(ndof), dqr(ndof), source=0)
    allocate(isov(ndof)) ; isov = .false.
    allocate(bvec(ndof), zvec(ndof), diagm(ndof), pcinv(ndof), source=0.0_dp)
    allocate(za1(nbf,nbf), zb1(nbf,nbf), gza(nbf,nbf), gzb(nbf,nbf), source=0.0_dp)
    k = 0
    do sp = 1, 2
      nocc = nocca ; if (sp == 2) nocc = noccb
      do p = 1, nbf ; do q = 1, p-1
        keep = .true. ; if (ovonly) keep = (p > nocc .and. q <= nocc)
        if (keep) then
          k = k + 1 ; dsp(k) = sp ; dpr(k) = p ; dqr(k) = q
          isov(k) = (p > nocc .and. q <= nocc)
        end if
      end do ; end do
    end do
    ndofov = count(isov)

    ! ---- RHS b = -R = -antisym(G^f) ; diagonal = eps_p - eps_q (the frozen part of M) ----
    do k = 1, ndof
      if (dsp(k) == 1) then
        bvec(k)  = -(gfa(dpr(k),dqr(k)) - gfa(dqr(k),dpr(k)))
        diagm(k) =   epsca(dpr(k)) - epsca(dqr(k))
      else
        bvec(k)  = -(gfb(dpr(k),dqr(k)) - gfb(dqr(k),dpr(k)))
        diagm(k) =   epscb(dpr(k)) - epscb(dqr(k))
      end if
    end do
    call sanitize_zvector_preconditioner(diagm, pcinv, iw, 1.0e-12_dp, 'UMRSF-Z')

    ! ---- (1) oo/vv block: pure-diagonal solve z_D = b_D / (eps_p - eps_q) ----
    zvec = 0.0_dp
    do k = 1, ndof
      if (.not. isov(k)) zvec(k) = pcinv(k) * bvec(k)
    end do

    ! ---- build the ov context (data for the pcg/minres matvec + Jacobi preconditioner) ----
    ctx%infos => infos ; ctx%basis => basis
    ctx%cac => cac ; ctx%cbc => cbc ; ctx%epsca => epsca ; ctx%epscb => epscb
    ctx%hfscale_ref = hfscale_ref ; ctx%nbf = nbf ; ctx%ndof_ov = ndofov
    allocate(ctx%dsp(ndofov), ctx%dpr(ndofov), ctx%dqr(ndofov), source=0)
    allocate(ctx%pcinv(ndofov), source=0.0_dp)
    allocate(ctx%za1(nbf,nbf), ctx%zb1(nbf,nbf), ctx%gza(nbf,nbf), ctx%gzb(nbf,nbf), source=0.0_dp)
    allocate(rhsov(max(ndofov,1)), source=0.0_dp)
    kk = 0
    do k = 1, ndof
      if (isov(k)) then
        kk = kk + 1
        ctx%dsp(kk) = dsp(k) ; ctx%dpr(kk) = dpr(k) ; ctx%dqr(kk) = dqr(k)
        ctx%pcinv(kk) = pcinv(k) ; rhsov(kk) = bvec(k)
      end if
    end do

    ! ---- (2) ov solve  A_{V,V} z_V = b_V - M_{V,D} z_D  via minres_optimize (pcg if SPD) ----
    ! A_{V,V} is SYMMETRIC but INDEFINITE (the excited-state/SOMO orbital Hessian has negative
    ! eigenvalues), so minres_optimize (Paige-Saunders, residual-minimizing, robust to indefiniteness)
    ! is the default — it converges like the old GMRES (~16 matvecs). pcg_optimize (UMRSF_ZPCG=1) STALLS
    ! here (CG requires SPD) and is kept only for a verified-SPD reference.
    rtol = 1.0e-11_dp ; mxit = min(ndofov, 5000) ; use_minres = .true.
    iters = 0 ; relres = 0.0_dp
    if (ndofov > 0) then
      ! correction M_{V,D} z_D = ov-antisym readout of umrsf_genfock_z applied to the oo/vv solution
      za1 = 0.0_dp ; zb1 = 0.0_dp
      do k = 1, ndof
        if (.not. isov(k)) then
          if (dsp(k) == 1) then ; za1(dpr(k),dqr(k)) = zvec(k) ; za1(dqr(k),dpr(k)) = zvec(k)
          else                  ; zb1(dpr(k),dqr(k)) = zvec(k) ; zb1(dqr(k),dpr(k)) = zvec(k) ; end if
        end if
      end do
      call umrsf_genfock_z(infos, basis, cac, cbc, epsca, epscb, za1, zb1, hfscale_ref, gza, gzb)
      kk = 0
      do k = 1, ndof
        if (isov(k)) then
          kk = kk + 1
          if (dsp(k) == 1) then ; rhsov(kk) = rhsov(kk) - (gza(dpr(k),dqr(k)) - gza(dqr(k),dpr(k)))
          else                  ; rhsov(kk) = rhsov(kk) - (gzb(dpr(k),dqr(k)) - gzb(dqr(k),dpr(k))) ; end if
        end if
      end do

      block
        character(len=24) :: e ; integer :: ios ; real(kind=dp) :: rv ; integer :: iv
        call get_environment_variable("UMRSF_ZTOL", e, status=ios)
        if (ios==0) then ; read(e,*,iostat=ios) rv ; if (ios==0 .and. rv>0.0_dp) rtol = rv ; end if
        call get_environment_variable("UMRSF_ZMAXIT", e, status=ios)
        if (ios==0) then ; read(e,*,iostat=ios) iv ; if (ios==0 .and. iv>0) mxit = iv ; end if
        call get_environment_variable("UMRSF_ZPCG", e, status=ios)
        if (ios==0 .and. trim(e)=="1") use_minres = .false.
      end block

      bnorm = sqrt(sum(rhsov**2))
      if (bnorm > tiny(1.0_dp)) then
        ! pcg/minres tol is on the residual NORM (absolute); target the relative tolerance rtol
        if (use_minres) then
          call minres_optimize(rhsov, umrsf_zov_matvec, umrsf_zov_precond, ctx, mxit, &
                               tol=rtol*bnorm, err=errout, iters=iters)
        else
          call pcg_optimize(rhsov, umrsf_zov_matvec, umrsf_zov_precond, ctx, mxit, &
                            tol=rtol*bnorm, err=errout, cgiters=cgit)
          iters = int(cgit)
        end if
        relres = errout / bnorm
      end if
      ! scatter z_V back into the full z
      kk = 0
      do k = 1, ndof
        if (isov(k)) then ; kk = kk + 1 ; zvec(k) = rhsov(kk) ; end if
      end do
    end if

    ! ---- unpack z -> symmetric MO matrices + AO relaxation densities (identical to dense) ----
    zmata = 0.0_dp ; zmatb = 0.0_dp
    do k = 1, ndof
      if (dsp(k) == 1) then ; zmata(dpr(k),dqr(k)) = zvec(k) ; zmata(dqr(k),dpr(k)) = zvec(k)
      else                  ; zmatb(dpr(k),dqr(k)) = zvec(k) ; zmatb(dqr(k),dpr(k)) = zvec(k) ; end if
    end do
    pza = matmul(cac, matmul(zmata, transpose(cac)))
    pzb = matmul(cbc, matmul(zmatb, transpose(cbc)))
    zrms = sqrt(sum(zvec**2)/max(ndof,1))

    ! ---- stationarity residual ||antisym(G^f + G^z(z))|| = max|M z - b| (-> 0 => M z = -R solved) ----
    call umrsf_genfock_z(infos, basis, cac, cbc, epsca, epscb, zmata, zmatb, hfscale_ref, gza, gzb)
    statio = 0.0_dp
    do k = 1, ndof
      if (dsp(k) == 1) then
        rr = (gfa(dpr(k),dqr(k)) - gfa(dqr(k),dpr(k))) + (gza(dpr(k),dqr(k)) - gza(dqr(k),dpr(k)))
      else
        rr = (gfb(dpr(k),dqr(k)) - gfb(dqr(k),dpr(k))) + (gzb(dpr(k),dqr(k)) - gzb(dqr(k),dpr(k)))
      end if
      statio = max(statio, abs(rr))
    end do
    sname = 'PCG' ; if (use_minres) sname = 'MINRES'
    write(iw,'(2x,3a,i0,a,i0,a,i0,a,es10.2)') 'iterative Z (', trim(sname), &
      '): ndof = ', ndof, ', ov-DOFs = ', ndofov, ', matvecs = ', iters, ', final rel-resid = ', relres

    deallocate(dsp, dpr, dqr, isov, bvec, zvec, diagm, pcinv, rhsov, za1, zb1, gza, gzb)
  end subroutine umrsf_zvector_iter

!###############################################################################
!> ov-block matvec for the reduced SPD Z-vector solve (pcg_optimize / minres_optimize `update`
!> callback). y = A_{V,V} x : embed x into the ov entries of a symmetric MO z, apply ONE
!> umrsf_genfock_z, read the ov antisym back = (eps_a-eps_i) x + 2 (C^T G[C z_V C^T] C)_{ai}. The
!> umrsf_zov_ctx_t is recovered from the solver's c_ptr `dat`. NB: no intent on the dummies — the
!> signature must match the abstract pcg_matvec / minres_matvec interfaces exactly.
  subroutine umrsf_zov_matvec(y, x, dat)
    use iso_c_binding, only: c_ptr, c_f_pointer
    real(kind=dp) :: y(:)
    real(kind=dp) :: x(:)
    type(c_ptr) :: dat
    type(umrsf_zov_ctx_t), pointer :: c
    integer :: kk
    call c_f_pointer(dat, c)
    c%za1 = 0.0_dp ; c%zb1 = 0.0_dp
    do kk = 1, c%ndof_ov
      if (c%dsp(kk) == 1) then ; c%za1(c%dpr(kk),c%dqr(kk)) = x(kk) ; c%za1(c%dqr(kk),c%dpr(kk)) = x(kk)
      else                     ; c%zb1(c%dpr(kk),c%dqr(kk)) = x(kk) ; c%zb1(c%dqr(kk),c%dpr(kk)) = x(kk) ; end if
    end do
    call umrsf_genfock_z(c%infos, c%basis, c%cac, c%cbc, c%epsca, c%epscb, &
                         c%za1, c%zb1, c%hfscale_ref, c%gza, c%gzb)
    do kk = 1, c%ndof_ov
      if (c%dsp(kk) == 1) then ; y(kk) = c%gza(c%dpr(kk),c%dqr(kk)) - c%gza(c%dqr(kk),c%dpr(kk))
      else                     ; y(kk) = c%gzb(c%dpr(kk),c%dqr(kk)) - c%gzb(c%dqr(kk),c%dpr(kk)) ; end if
    end do
  end subroutine umrsf_zov_matvec

!###############################################################################
!> Jacobi preconditioner for the ov SPD solve: y = (eps_a - eps_i)^-1 x (floored). pcg/minres
!> `precond` callback (no intent on the dummies — must match the abstract interface).
  subroutine umrsf_zov_precond(y, x, dat)
    use iso_c_binding, only: c_ptr, c_f_pointer
    real(kind=dp) :: y(:)
    real(kind=dp) :: x(:)
    type(c_ptr) :: dat
    type(umrsf_zov_ctx_t), pointer :: c
    call c_f_pointer(dat, c)
    y = c%pcinv * x
  end subroutine umrsf_zov_precond

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
!> n×n identity matrix (small helper for S-orthonormality probes).
  pure function id_nbf(n) result(idm)
    implicit none
    integer, intent(in) :: n
    real(kind=dp) :: idm(n,n)
    integer :: i
    idm = 0.0_dp
    do i = 1, n ; idm(i,i) = 1.0_dp ; end do
  end function id_nbf

!###############################################################################
!> max within-segment |btt| = the get_jacobi STATIONARITY residual (c05 align_offdiag_btt).
!> s_mo = vaᵀ S vb (columns normalized); seg0 {1..nocca-1} btt=aa·dd−bb·cc, seg1 {nocca..nbf}
!> btt=aa·cc−bb·dd with aa=s(i,i) bb=s(j,j) cc=s(i,j) dd=s(j,i). 0 at the converged alignment.
  function within_seg_btt(va, vb, smat_full, nocca) result(m)
    implicit none
    real(kind=dp), intent(in) :: va(:,:), vb(:,:), smat_full(:,:)
    integer, intent(in) :: nocca
    real(kind=dp) :: m
    real(kind=dp), allocatable :: s_mo(:,:)
    integer :: nbf, i, j, slo, shi, seg
    real(kind=dp) :: aa, bb, cc, dd, btt
    nbf = size(va,1)
    allocate(s_mo(nbf,nbf))
    s_mo = matmul(transpose(va), matmul(smat_full, vb))
    do i = 1, nbf ; s_mo(:,i) = s_mo(:,i)/max(norm2(s_mo(:,i)),1.0e-10_dp) ; end do
    m = 0.0_dp
    do seg = 0, 1
      if (seg==0) then ; slo=1 ; shi=nocca-1 ; else ; slo=nocca ; shi=nbf ; end if
      do i = slo, shi-1
        do j = i+1, shi
          aa = s_mo(i,i) ; bb = s_mo(j,j) ; cc = s_mo(i,j) ; dd = s_mo(j,i)
          if (seg==0) then ; btt = aa*dd - bb*cc ; else ; btt = aa*cc - bb*dd ; end if
          m = max(m, abs(btt))
        end do
      end do
    end do
    deallocate(s_mo)
  end function within_seg_btt

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
