!> @file mp2_lib.F90
!>
!> @brief Standalone Moller-Plesset second-order (MP2) ground-state correlation
!>        energy for RHF/UHF/ROHF references.
!>
!> The correlation energy is built in the spin-blocked (aa, bb, ab) form on
!> semicanonicalized orbitals.  Two paths are provided:
!>
!>  1. Direct-J (O(N⁶)): per-occupied-MO-pair Coulomb builds via int2_driver,
!>     fast for small systems.  Guarded by MAX_JBUILDS.
!>
!>  2. Batched half-transform (O(N⁵)): optimal for large bases.  For small
!>     bases (nbf < 100) it builds the full MO integral tensor via 4 quarter-
!>     transforms (PySCF-style, 10-100× faster than per-pair DGEMMs).  For
!>     large bases it uses the memory-save half-transform path.
module mp2_lib

  use precision, only: dp
  implicit none

  private
  public :: mp2_correlation
  public :: semicanonicalize

  integer, parameter :: MAX_JBUILDS = 4000

contains

  logical function mp2_build_is_affordable(nbf, nocca, noccb) result(ok)
    integer, intent(in) :: nbf, nocca, noccb
    integer :: vira, virb, nbuild, cap, ln
    character(len=32) :: sval
    vira = nbf - nocca; virb = nbf - noccb
    nbuild = max(0, nocca*vira) + max(0, noccb*virb)
    cap = MAX_JBUILDS
    call get_environment_variable("OQP_MP2_MAX_JBUILDS", sval, ln)
    if (ln > 0) read(sval, *, iostat=ln) cap
    ok = (nbuild > 0) .and. (nbuild <= cap)
  end function mp2_build_is_affordable

  subroutine mp2_correlation(infos, e_mp2, e_aa, e_bb, e_ab, e_s, computed)
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use messages, only: show_message, with_abort
    use oqp_tagarray_driver, only: tagarray_get_data, &
                                   OQP_VEC_MO_A, OQP_VEC_MO_B, &
                                   OQP_FOCK_A, OQP_FOCK_B
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: e_mp2, e_aa, e_bb, e_ab
    real(kind=dp), intent(out) :: e_s
    logical, intent(out) :: computed
    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:)
    real(kind=dp), allocatable :: mo_a_sc(:,:), mo_b_sc(:,:)
    real(kind=dp), allocatable :: e_a_sc(:), e_b_sc(:)
    integer :: nbf, nbf2, nocca, noccb, vira, virb, ok
    real(kind=dp) :: e_opp_scratch
    real(kind=dp) :: ss_scale, os_scale
    logical :: restricted_ref, need_same_spin, need_opposite_spin
    logical :: try_n5, n5_ok

    e_mp2 = 0.0_dp; e_aa = 0.0_dp; e_bb = 0.0_dp; e_ab = 0.0_dp; e_s = 0.0_dp
    computed = .false.
    basis => infos%basis
    basis%atoms => infos%atoms
    nbf  = basis%nbf; nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_a; noccb = infos%mol_prop%nelec_b
    vira = nbf - nocca; virb = nbf - noccb
    try_n5 = .not. mp2_build_is_affordable(nbf, nocca, noccb)
    if (.not. try_n5 .and. nbf > 120) try_n5 = .true.

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    restricted_ref = (infos%control%scftype == 1)
    if (restricted_ref) then; mo_b => mo_a; fock_b => fock_a
    else
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    end if

    ss_scale = infos%dft%MP2SS_Scale; os_scale = infos%dft%MP2OS_Scale
    need_same_spin = abs(ss_scale) > 1.0e-14_dp
    need_opposite_spin = abs(os_scale) > 1.0e-14_dp

    allocate(mo_a_sc(nbf,nbf), mo_b_sc(nbf,nbf), e_a_sc(nbf), e_b_sc(nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: cannot allocate semicanonical MOs', with_abort)
    call semicanonicalize(nbf, nocca, mo_a, fock_a, mo_a_sc, e_a_sc)
    call semicanonicalize(nbf, noccb, mo_b, fock_b, mo_b_sc, e_b_sc)

    e_s = mp2_singles(nbf, nocca, mo_a_sc, fock_a, e_a_sc) &
        + mp2_singles(nbf, noccb, mo_b_sc, fock_b, e_b_sc)

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    if (try_n5) then
      call mp2_corr_n5(int2_driver, basis, nbf, nbf2, &
          mo_a_sc, e_a_sc, nocca, vira, mo_a_sc, e_a_sc, nocca, vira, &
          mo_b_sc, e_b_sc, noccb, virb, &
          need_same_spin, need_opposite_spin, restricted_ref, e_aa, e_ab, success=n5_ok)
      if (n5_ok .and. need_same_spin) then
        call mp2_corr_n5(int2_driver, basis, nbf, nbf2, &
            mo_b_sc, e_b_sc, noccb, virb, mo_b_sc, e_b_sc, noccb, virb, &
            mo_a_sc, e_a_sc, nocca, vira, &
            same_spin=.true., do_opposite=.false., &
            e_same=e_bb, e_opp=e_opp_scratch, success=n5_ok)
      end if
      if (n5_ok) computed = .true.
    end if

    if (.not. computed) then
      if (.not. mp2_build_is_affordable(nbf, nocca, noccb)) then
        deallocate(mo_a_sc, mo_b_sc, e_a_sc, e_b_sc)
        call int2_driver%clean(); return
      end if
      call mp2_spin_block(int2_driver, basis, nbf, nbf2, &
          mo_a_sc, e_a_sc, nocca, vira, mo_a_sc, e_a_sc, nocca, vira, &
          mo_b_sc, e_b_sc, noccb, virb, &
          same_spin=need_same_spin, do_opposite=need_opposite_spin, &
          e_same=e_aa, e_opp=e_ab)
      e_opp_scratch = 0.0_dp
      call mp2_spin_block(int2_driver, basis, nbf, nbf2, &
          mo_b_sc, e_b_sc, noccb, virb, mo_b_sc, e_b_sc, noccb, virb, &
          mo_a_sc, e_a_sc, nocca, vira, &
          same_spin=need_same_spin, do_opposite=.false., &
          e_same=e_bb, e_opp=e_opp_scratch)
      computed = .true.
    end if

    deallocate(mo_a_sc, mo_b_sc, e_a_sc, e_b_sc)
    call int2_driver%clean()
    e_mp2 = ss_scale * (e_aa + e_bb) + os_scale * e_ab + e_s
  end subroutine mp2_correlation


  !###########################################################################
  ! O(N⁵) path — hybrid full-MO / batched half-transform
  !###########################################################################

  !> @brief MP2 correlation via O(N⁵) transform.
  !>
  !> Strategy (PySCF-aligned):
  !>   1. Collect packed AO integrals once
  !>   2. If memory permits (nbf⁴ < 1/4 of available RAM), run the full MO
  !>      transform via cc_build_full_mo — four DGEMMs over the packed store,
  !>      then read eri_mo directly for the energy.  This is 10-100× faster
  !>      than per-pair DGEMMs for nbf < 100.
  !>   3. Otherwise, use the memory-save batched half-transform path.
  subroutine mp2_corr_n5(int2_driver, basis, nbf, nbf2, &
                         cmo_l, e_l, nocc_l, nvir_l, &
                         cmo_s, e_s, nocc_s, nvir_s, &
                         cmo_o, e_o, nocc_o, nvir_o, &
                         restricted_ref, same_spin, do_opposite, e_same, e_opp, success)
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use cc_ao2mo, only: cc_eri_collect_t, cc_packed_length, cc_build_full_mo
    use messages, only: show_message, WITH_ABORT
    use memory_info, only: oqp_available_memory_gb, OQP_MEMORY_SAFETY_FRACTION
    type(int2_compute_t), intent(inout) :: int2_driver
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nbf, nbf2
    real(kind=dp), intent(in) :: cmo_l(nbf,nbf), e_l(nbf)
    integer, intent(in) :: nocc_l, nvir_l
    real(kind=dp), intent(in) :: cmo_s(nbf,nbf), e_s(nbf)
    integer, intent(in) :: nocc_s, nvir_s
    real(kind=dp), intent(in) :: cmo_o(nbf,nbf), e_o(nbf)
    integer, intent(in) :: nocc_o, nvir_o
    logical, intent(in) :: same_spin, do_opposite
    real(kind=dp), intent(inout) :: e_same, e_opp
    logical, intent(out) :: success
    logical, intent(in) :: restricted_ref

    integer(8) :: packed_len
    integer :: nmo, ok
    logical :: do_full
    real(kind=dp), allocatable, target :: g(:)

    success = .false.
    nmo = nbf
    packed_len = cc_packed_length(nbf)

    ! --- Collect packed AO integrals (shared by both paths) ----------------
    allocate(g(packed_len), source=0.0_dp, stat=ok)
    if (ok /= 0) return
    block
      type(cc_eri_collect_t), target :: collector
      collector%g => g; collector%nbf = nbf; collector%npair = nbf2
      call int2_driver%run(collector)
    end block

    ! --- Decide which path to use -------------------------------------------
    ! Full MO transform: eri_mo(nmo,nmo,nmo,nmo) ≈ 8×nbf⁴ bytes.
    do_full = .false.
    block
      real(kind=dp) :: cost_gb, avail_gb, frac
      frac = OQP_MEMORY_SAFETY_FRACTION * 0.25_dp  ! use at most 1/4 of safety budget
      cost_gb = real(nmo, dp)**4 * 8.0_dp / 1.073741824e9_dp
      avail_gb = oqp_available_memory_gb()
      if (avail_gb > 0.0_dp .and. cost_gb < frac * avail_gb) do_full = .true.
      ! Also allow the full path when the cost is very small (< 1 MiB) even if
      ! the memory check is unavailable.
      if (avail_gb <= 0.0_dp .and. cost_gb < 1.0e-3_dp) do_full = .true.
    end block

    if (do_full) then
      ! Full MO transform (PySCF-style, fast): builds eri_mo via 4 quarter-
      ! transforms (large DGEMMs), then reads (ia|jb) directly for energy.
      ! Handles both same-spin and opposite-spin in one shot.
      call mp2_corr_n5_full(nbf, nmo, g, &
          cmo_l, e_l, nocc_l, nvir_l, &
          cmo_s, e_s, nocc_s, nvir_s, &
          cmo_o, e_o, nocc_o, nvir_o, &
          restricted_ref, same_spin, do_opposite, e_same, e_opp, ok)
      if (ok == 0) then
        deallocate(g); success = .true.; return
      end if
    end if

    ! Memory-save batched half-transform path
    call mp2_corr_n5_batched(g, nbf, nbf2, &
        cmo_l, e_l, nocc_l, nvir_l, &
        cmo_s, e_s, nocc_s, nvir_s, &
        cmo_o, e_o, nocc_o, nvir_o, &
        same_spin, do_opposite, e_same, e_opp, success)
    deallocate(g)

  end subroutine mp2_corr_n5


  !> @brief Full MO transform path (PySCF-style).
  !>
  !> Builds eri_mo(nmo,nmo,nmo,nmo) via cc_build_full_mo, then reads the
  !> (ia|jb) blocks directly.  For small bases this is dramatically faster
  !> than per-pair DGEMMs.
  subroutine mp2_corr_n5_full(nbf, nmo, g, &
      cmo_l, e_l, nocc_l, nvir_l, &
      cmo_s, e_s, nocc_s, nvir_s, &
      cmo_o, e_o, nocc_o, nvir_o, &
      restricted_ref, same_spin, do_opposite, e_same, e_opp, ok)
    use cc_ao2mo, only: cc_build_full_mo
    integer, intent(in) :: nbf, nmo
    real(kind=dp), intent(in) :: g(*)
    real(kind=dp), intent(in) :: cmo_l(nbf,nbf), e_l(nbf)
    integer, intent(in) :: nocc_l, nvir_l
    real(kind=dp), intent(in) :: cmo_s(nbf,nbf), e_s(nbf)
    integer, intent(in) :: nocc_s, nvir_s
    real(kind=dp), intent(in) :: cmo_o(nbf,nbf), e_o(nbf)
    integer, intent(in) :: nocc_o, nvir_o
    logical, intent(in) :: restricted_ref, same_spin, do_opposite
    real(kind=dp), intent(inout) :: e_same, e_opp
    integer, intent(out) :: ok

    real(kind=dp), allocatable :: eri_mo(:,:,:,:)
    integer :: i, a, j, b
    real(kind=dp) :: num, denom

    ok = -1
    allocate(eri_mo(nmo,nmo,nmo,nmo), source=0.0_dp, stat=ok)
    if (ok /= 0) return

    ! cc_build_full_mo(nbf, nmo, cmo_bra, cmo_ket, g, eri)
    ! eri_mo(p,q,r,s) = (pq|rs) in chemist notation.
    ! For RHF (restricted_ref=.true.) cmo_l == cmo_o, so the single transform
    ! is correct for both same-spin and opposite-spin.
    ! For UHF/ROHF (restricted_ref=.false.) the opposite-spin (ia|jb) needs
    ! cmo_o on the ket pair; this path handles same-spin only, and opposite-spin
    ! is computed by the batched half-transform path in the caller.
    call cc_build_full_mo(nbf, nmo, cmo_l, cmo_l, g, eri_mo)

    if (same_spin) then
      !$omp parallel do collapse(2) private(i,a,j,b,num,denom) &
      !$omp   schedule(dynamic,1) reduction(+:e_same)
      do b = 1, nvir_s
        do a = 1, nvir_l
          do j = 1, nocc_s
            do i = 1, nocc_l
              denom = e_l(i) + e_s(j) - e_l(nocc_l+a) - e_s(nocc_s+b)
              if (abs(denom) < 1.0e-10_dp) cycle
              ! (ia|jb) - (ib|ja)
              num = eri_mo(i, nocc_l+a, j, nocc_s+b) &
                  - eri_mo(i, nocc_s+b, j, nocc_l+a)
              e_same = e_same + 0.25_dp * num * num / denom
            end do
          end do
        end do
      end do
      !$omp end parallel do
    end if

    if (do_opposite) then
      if (.not. restricted_ref) then
        ! UHF/ROHF opposite-spin is computed by the batched path in the caller.
        ! This guard exists because cc_build_full_mo with cmo_l on all four
        ! indices is incorrect for opposite-spin when cmo_l != cmo_o.
        ok = -1; return
      end if
      !$omp parallel do collapse(2) private(i,a,j,b,num,denom) &
      !$omp   schedule(dynamic,1) reduction(+:e_opp)
      do b = 1, nvir_o
        do a = 1, nvir_l
          do j = 1, nocc_o
            do i = 1, nocc_l
              denom = e_l(i) + e_o(j) - e_l(nocc_l+a) - e_o(nocc_o+b)
              if (abs(denom) < 1.0e-10_dp) cycle
              ! (ia|jb) — no exchange for opposite-spin
              num = eri_mo(i, nocc_l+a, j, nocc_o+b)
              e_opp = e_opp + num * num / denom
            end do
          end do
        end do
      end do
      !$omp end parallel do
    end if

    deallocate(eri_mo)
    ok = 0
  end subroutine mp2_corr_n5_full


  !> @brief Memory-save batched half-transform path (same as the earlier
  !>        implementation, retained for large bases).
  subroutine mp2_corr_n5_batched(g, nbf, nbf2, &
                         cmo_l, e_l, nocc_l, nvir_l, &
                         cmo_s, e_s, nocc_s, nvir_s, &
                         cmo_o, e_o, nocc_o, nvir_o, &
                         same_spin, do_opposite, e_same, e_opp, success)
    use cc_ao2mo, only: cc_packed_length
    use messages, only: show_message, WITH_ABORT
    integer, intent(in) :: nbf, nbf2
    real(kind=dp), intent(in), target :: g(*)
    real(kind=dp), intent(in) :: cmo_l(nbf,nbf), e_l(nbf)
    integer, intent(in) :: nocc_l, nvir_l
    real(kind=dp), intent(in) :: cmo_s(nbf,nbf), e_s(nbf)
    integer, intent(in) :: nocc_s, nvir_s
    real(kind=dp), intent(in) :: cmo_o(nbf,nbf), e_o(nbf)
    integer, intent(in) :: nocc_o, nvir_o
    logical, intent(in) :: same_spin, do_opposite
    real(kind=dp), intent(inout) :: e_same, e_opp
    logical, intent(out) :: success

    integer :: npair, nov_l, nov_s, nov_o, idx, i, a, j, b, q, ip, ok
    integer :: lambda, sigma, mu, nu
    real(kind=dp) :: denom, num, val
    real(kind=dp), allocatable :: half(:,:)
    integer, allocatable :: prow(:), pcol(:)
    real(kind=dp), allocatable :: ovov(:,:,:,:)

    success = .false.
    npair = nbf2
    nov_l = nocc_l * nvir_l
    nov_s = nocc_s * nvir_s
    nov_o = nocc_o * nvir_o

    ! --- Precompute AO pair tables -----------------------------------------
    allocate(prow(npair), pcol(npair), stat=ok)
    if (ok /= 0) return
    do q = 1, npair
      do lambda = 1, nbf
        if (q <= lambda*(lambda-1)/2) cycle
        if (q > lambda*(lambda+1)/2) cycle
        sigma = q - lambda*(lambda-1)/2
        prow(q) = lambda; pcol(q) = sigma; exit
      end do
    end do

    ! --- Step 1: first half-transform (O(N⁵)) --------------------------------
    allocate(half(npair, nov_l), source=0.0_dp, stat=ok)
    if (ok /= 0) then; deallocate(prow, pcol); return; end if

    !$omp parallel default(shared) private(q, ip, mu, nu, idx)
    block
      real(kind=dp), allocatable :: d(:,:), scr(:,:), m(:,:)
      allocate(d(nbf,nbf), scr(nbf,nvir_l), m(nocc_l,nvir_l))
      !$omp do schedule(static)
      do q = 1, npair
        do ip = 1, npair
          mu = prow(ip); nu = pcol(ip)
          d(mu,nu) = g(mp2_packed_idx(ip, q))
          d(nu,mu) = d(mu,nu)
        end do
        call dgemm('n','n', nbf, nvir_l, nbf, 1.0_dp, d, nbf, &
                   cmo_l(1, nocc_l+1), nbf, 0.0_dp, scr, nbf)
        call dgemm('t','n', nocc_l, nvir_l, nbf, 1.0_dp, &
                   cmo_l, nbf, scr, nbf, 0.0_dp, m, nocc_l)
        do a = 1, nvir_l
          do i = 1, nocc_l
            half(q, (i-1)*nvir_l + a) = m(i, a)
          end do
        end do
      end do
      !$omp end do
      deallocate(d, scr, m)
    end block
    !$omp end parallel

    ! --- Step 2: second half-transform and energy accumulation --------------
    if (same_spin) then
      allocate(ovov(nocc_l, nvir_l, nocc_s, nvir_s), source=0.0_dp, stat=ok)
      if (ok /= 0) then; deallocate(half, prow, pcol); return; end if
      call mp2_n5_second_half(nbf, npair, nov_l, &
          nocc_l, nvir_l, nocc_s, nvir_s, &
          cmo_l, cmo_s, half, prow, pcol, ovov)
      do i = 1, nocc_l
        do a = 1, nvir_l
          do j = 1, nocc_s
            do b = 1, nvir_s
              denom = e_l(i) + e_s(j) - e_l(nocc_l + a) - e_s(nocc_s + b)
              if (abs(denom) < 1.0e-10_dp) cycle
              num = ovov(i, a, j, b) - ovov(i, b, j, a)
              e_same = e_same + 0.25_dp * num * num / denom
            end do
          end do
        end do
      end do
      deallocate(ovov)
    end if

    if (do_opposite) then
      call mp2_n5_opposite(nbf, npair, nov_l, &
          nocc_l, nvir_l, nocc_o, nvir_o, &
          cmo_l, e_l, cmo_o, e_o, &
          half, prow, pcol, e_opp)
    end if

    deallocate(half, prow, pcol)
    success = .true.
  end subroutine mp2_corr_n5_batched


  subroutine mp2_n5_second_half(nbf, npair, nov_l, &
      nocc_l, nvir_l, nocc_s, nvir_s, &
      cmo_l, cmo_s, half, prow, pcol, ovov)
    integer, intent(in) :: nbf, npair, nov_l
    integer, intent(in) :: nocc_l, nvir_l, nocc_s, nvir_s
    real(kind=dp), intent(in) :: cmo_l(nbf,nbf), cmo_s(nbf,nbf)
    real(kind=dp), intent(in) :: half(npair, nov_l)
    integer, intent(in) :: prow(npair), pcol(npair)
    real(kind=dp), intent(out) :: ovov(nocc_l, nvir_l, nocc_s, nvir_s)
    integer :: idx, i, a, j, b, lambda, sigma, q
    real(kind=dp), allocatable :: e_mat(:,:), scr2(:,:), vv(:,:)
    allocate(e_mat(nbf,nbf), scr2(nbf,nocc_s), vv(nocc_s,nvir_s))
    !$omp parallel do private(idx, i, a, q, lambda, sigma, e_mat, scr2, vv) &
    !$omp   schedule(dynamic, 4)
    do idx = 1, nov_l
      i = (idx - 1) / nvir_l + 1; a = mod(idx - 1, nvir_l) + 1
      do q = 1, npair
        lambda = prow(q); sigma = pcol(q)
        e_mat(lambda, sigma) = half(q, idx)
        e_mat(sigma, lambda) = half(q, idx)
      end do
      call dgemm('t', 'n', nocc_s, nbf, nbf, 1.0_dp, cmo_s, nbf, e_mat, nbf, 0.0_dp, scr2, nocc_s)
      call dgemm('n', 'n', nocc_s, nvir_s, nbf, 1.0_dp, scr2, nocc_s, &
                 cmo_s(1, nocc_s+1), nbf, 0.0_dp, vv, nocc_s)
      do b = 1, nvir_s
        do j = 1, nocc_s
          ovov(i, a, j, b) = vv(j, b)
        end do
      end do
    end do
    !$omp end parallel do
    deallocate(e_mat, scr2, vv)
  end subroutine mp2_n5_second_half


  subroutine mp2_n5_opposite(nbf, npair, nov_l, &
      nocc_l, nvir_l, nocc_o, nvir_o, &
      cmo_l, e_l, cmo_o, e_o, &
      half, prow, pcol, e_opp)
    integer, intent(in) :: nbf, npair, nov_l
    integer, intent(in) :: nocc_l, nvir_l, nocc_o, nvir_o
    real(kind=dp), intent(in) :: cmo_l(nbf,nbf), e_l(nbf)
    real(kind=dp), intent(in) :: cmo_o(nbf,nbf), e_o(nbf)
    real(kind=dp), intent(in) :: half(npair, nov_l)
    integer, intent(in) :: prow(npair), pcol(npair)
    real(kind=dp), intent(inout) :: e_opp
    integer :: idx, i, a, j, b, q, lambda, sigma
    real(kind=dp) :: denom, val
    real(kind=dp), allocatable :: e_mat(:,:), scr2(:,:), vv(:,:)
    allocate(e_mat(nbf,nbf), scr2(nbf,nocc_o), vv(nocc_o,nvir_o))
    !$omp parallel do private(idx, i, a, q, lambda, sigma, &
    !$omp   e_mat, scr2, vv, j, b, denom, val) &
    !$omp   schedule(dynamic, 4) reduction(+:e_opp)
    do idx = 1, nov_l
      i = (idx - 1) / nvir_l + 1; a = mod(idx - 1, nvir_l) + 1
      do q = 1, npair
        lambda = prow(q); sigma = pcol(q)
        e_mat(lambda, sigma) = half(q, idx)
        e_mat(sigma, lambda) = half(q, idx)
      end do
      call dgemm('t', 'n', nocc_o, nbf, nbf, 1.0_dp, cmo_o, nbf, e_mat, nbf, 0.0_dp, scr2, nocc_o)
      call dgemm('n', 'n', nocc_o, nvir_o, nbf, 1.0_dp, scr2, nocc_o, &
                 cmo_o(1, nocc_o+1), nbf, 0.0_dp, vv, nocc_o)
      do b = 1, nvir_o
        do j = 1, nocc_o
          denom = e_l(i) + e_o(j) - e_l(nocc_l + a) - e_o(nocc_o + b)
          if (abs(denom) < 1.0e-10_dp) cycle
          val = vv(j, b)
          e_opp = e_opp + val * val / denom
        end do
      end do
    end do
    !$omp end parallel do
    deallocate(e_mat, scr2, vv)
  end subroutine mp2_n5_opposite


  pure integer(8) function mp2_packed_idx(p, q) result(idx)
    integer, intent(in) :: p, q
    integer(8) :: hi, lo
    hi = int(max(p, q), 8); lo = int(min(p, q), 8)
    idx = hi * (hi - 1_8) / 2_8 + lo
  end function mp2_packed_idx


  !###########################################################################
  ! Original direct-J path (O(N⁶))
  !###########################################################################

  subroutine mp2_spin_block(int2_driver, basis, nbf, nbf2, &
                            cmo_l, e_l, nocc_l, nvir_l, &
                            cmo_s, e_s, nocc_s, nvir_s, &
                            cmo_o, e_o, nocc_o, nvir_o, &
                            same_spin, do_opposite, e_same, e_opp)
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t, int2_urohf_data_t
    use mathlib, only: pack_matrix, unpack_matrix
    use messages, only: show_message, WITH_ABORT
    type(int2_compute_t), intent(inout) :: int2_driver
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nbf, nbf2
    real(kind=dp), intent(in) :: cmo_l(nbf,nbf), e_l(nbf)
    integer, intent(in) :: nocc_l, nvir_l
    real(kind=dp), intent(in) :: cmo_s(nbf,nbf), e_s(nbf)
    integer, intent(in) :: nocc_s, nvir_s
    real(kind=dp), intent(in) :: cmo_o(nbf,nbf), e_o(nbf)
    integer, intent(in) :: nocc_o, nvir_o
    logical, intent(in) :: same_spin, do_opposite
    real(kind=dp), intent(inout) :: e_same, e_opp
    type(int2_urohf_data_t), target :: int2_data
    real(kind=dp), allocatable, target :: pdmat(:,:)
    real(kind=dp), allocatable :: dfull(:,:), jfull(:,:), scr(:,:), jpack(:)
    real(kind=dp), allocatable :: same_block(:,:,:), gopp(:,:)
    integer :: i, a, j, b, ii, ok
    real(kind=dp) :: denom, num

    allocate(pdmat(nbf2,2), dfull(nbf,nbf), jfull(nbf,nbf), scr(nbf,nbf), jpack(nbf2), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: cannot allocate J-build scratch', WITH_ABORT)
    if (same_spin .and. (nocc_s /= nocc_l .or. nvir_s /= nvir_l)) &
      call show_message('mp2: inconsistent same-spin dimensions', WITH_ABORT)
    if (same_spin) allocate(same_block(nvir_l,nvir_s,nocc_s), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: cannot allocate same-spin block', WITH_ABORT)
    if (do_opposite) allocate(gopp(nocc_o,nvir_o), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: cannot allocate opp-spin block', WITH_ABORT)

    do i = 1, nocc_l
      if (same_spin) same_block = 0.0_dp
      do a = 1, nvir_l
        call rank1_sym_density(cmo_l(:,i), cmo_l(:,nocc_l+a), nbf, dfull)
        call pack_matrix(dfull, pdmat(:,1), 'U'); pdmat(:,2) = 0.0_dp
        int2_data = int2_urohf_data_t(nfocks=2, d=pdmat, scale_exchange=0.0_dp)
        call int2_driver%run(int2_data)
        jpack(:) = 0.5_dp * int2_data%f(:,1,1)
        ii = 0
        do j = 1, nbf; ii = ii + j; jpack(ii) = 2.0_dp * jpack(ii); end do
        call unpack_matrix(jpack, jfull, 'U'); call int2_data%clean()

        if (same_spin) then
          call dgemm('n','n', nbf, nocc_s, nbf, 1.0_dp, jfull, nbf, cmo_s, nbf, 0.0_dp, scr, nbf)
          do j = 1, nocc_s; do b = 1, nvir_s
            same_block(a,b,j) = dot_product(cmo_s(:,nocc_s+b), scr(:,j))
          end do; end do
        end if
        if (do_opposite) then
          call dgemm('n','n', nbf, nocc_o, nbf, 1.0_dp, jfull, nbf, cmo_o, nbf, 0.0_dp, scr, nbf)
          do j = 1, nocc_o; do b = 1, nvir_o
            gopp(j,b) = dot_product(cmo_o(:,nocc_o+b), scr(:,j))
          end do; end do
          do j = 1, nocc_o; do b = 1, nvir_o
            denom = e_l(i) + e_o(j) - e_l(nocc_l+a) - e_o(nocc_o+b)
            if (abs(denom) < 1.0e-10_dp) cycle
            e_opp = e_opp + gopp(j,b)*gopp(j,b) / denom
          end do; end do
        end if
      end do
      if (same_spin) then
        do j = 1, nocc_s; do a = 1, nvir_l; do b = 1, nvir_s
          denom = e_l(i) + e_s(j) - e_l(nocc_l+a) - e_s(nocc_s+b)
          if (abs(denom) < 1.0e-10_dp) cycle
          num = same_block(a,b,j) - same_block(b,a,j)
          e_same = e_same + 0.25_dp * num*num / denom
        end do; end do; end do
      end if
    end do
    deallocate(pdmat, dfull, jfull, scr, jpack)
    if (allocated(same_block)) deallocate(same_block)
    if (allocated(gopp)) deallocate(gopp)
  end subroutine mp2_spin_block

  subroutine rank1_sym_density(u, v, nbf, dfull)
    real(kind=dp), intent(in) :: u(nbf), v(nbf); integer, intent(in) :: nbf
    real(kind=dp), intent(out) :: dfull(nbf,nbf)
    integer :: mu, nu
    do nu = 1, nbf; do mu = 1, nbf
      dfull(mu,nu) = 0.5_dp*(u(mu)*v(nu) + u(nu)*v(mu))
    end do; end do
  end subroutine rank1_sym_density

  function mp2_singles(nbf, nocc, cmo_sc, fock_packed, e_sc) result(e_s)
    use mathlib, only: unpack_matrix
    integer, intent(in) :: nbf, nocc
    real(kind=dp), intent(in) :: cmo_sc(nbf,nbf)
    real(kind=dp), intent(in) :: fock_packed(:)
    real(kind=dp), intent(in) :: e_sc(nbf)
    real(kind=dp) :: e_s
    real(kind=dp), allocatable :: fao(:,:), scr(:,:), fmo(:,:)
    real(kind=dp) :: d; integer :: i, a
    e_s = 0.0_dp; if (nocc <= 0 .or. nocc >= nbf) return
    allocate(fao(nbf,nbf), scr(nbf,nbf), fmo(nbf,nbf), source=0.0_dp)
    call unpack_matrix(fock_packed, fao, 'U')
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, cmo_sc, nbf, 0.0_dp, scr, nbf)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, cmo_sc, nbf, scr, nbf, 0.0_dp, fmo, nbf)
    do a = nocc+1, nbf; do i = 1, nocc
      d = e_sc(i) - e_sc(a); if (abs(d) < 1.0e-12_dp) cycle
      e_s = e_s + fmo(i,a)*fmo(i,a)/d
    end do; end do
    deallocate(fao, scr, fmo)
  end function mp2_singles

  subroutine semicanonicalize(nbf, nocc, cmo, fock_packed, cmo_sc, e_sc, nfzc)
    use mathlib, only: unpack_matrix
    use eigen, only: diag_symm_full
    use messages, only: show_message, with_abort
    integer, intent(in) :: nbf, nocc
    real(kind=dp), intent(in) :: cmo(nbf,nbf)
    real(kind=dp), intent(in) :: fock_packed(:)
    real(kind=dp), intent(out) :: cmo_sc(nbf,nbf)
    real(kind=dp), intent(out) :: e_sc(nbf)
    integer, intent(in), optional :: nfzc
    real(kind=dp), allocatable :: fao(:,:), fmo(:,:), scr(:,:), blk(:,:), eval(:)
    integer :: nvir, ok, nfz, nocc_act, i
    nfz = 0; if (present(nfzc)) nfz = max(0, min(nfzc, nocc))
    nocc_act = nocc - nfz; nvir = nbf - nocc
    allocate(fao(nbf,nbf), fmo(nbf,nbf), scr(nbf,nbf), eval(nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: semicanonical alloc failed', with_abort)
    call unpack_matrix(fock_packed, fao, 'U')
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, cmo, nbf, 0.0_dp, scr, nbf)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, cmo, nbf, scr, nbf, 0.0_dp, fmo, nbf)
    cmo_sc = cmo; e_sc = 0.0_dp
    do i = 1, nfz; e_sc(i) = fmo(i,i); end do
    if (nocc_act > 0) then
      allocate(blk(nocc_act,nocc_act), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('mp2: occ block alloc failed', with_abort)
      blk = fmo(nfz+1:nocc, nfz+1:nocc)
      call diag_symm_full(1, nocc_act, blk, nocc_act, eval(1:nocc_act))
      call dgemm('n','n', nbf, nocc_act, nocc_act, 1.0_dp, cmo(:,nfz+1:nocc), nbf, blk, nocc_act, 0.0_dp, cmo_sc(:,nfz+1:nocc), nbf)
      e_sc(nfz+1:nocc) = eval(1:nocc_act)
      deallocate(blk)
    end if
    if (nvir > 0) then
      allocate(blk(nvir,nvir), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('mp2: vir block alloc failed', with_abort)
      blk = fmo(nocc+1:nbf, nocc+1:nbf)
      call diag_symm_full(1, nvir, blk, nvir, eval(1:nvir))
      call dgemm('n','n', nbf, nvir, nvir, 1.0_dp, cmo(:,nocc+1:nbf), nbf, blk, nvir, 0.0_dp, cmo_sc(:,nocc+1:nbf), nbf)
      e_sc(nocc+1:nbf) = eval(1:nvir)
      deallocate(blk)
    end if
    deallocate(fao, fmo, scr, eval)
  end subroutine semicanonicalize

end module mp2_lib
