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
!>  2. Batched half-transform (O(N⁵)): collect packed AO integrals once, then
!>     two half-transforms produce the (ia|jb) block.  No per-pair J-build
!>     limit.  Used automatically when nocc×nvir exceeds MAX_JBUILDS or nbf
!>     exceeds 120 (O(N⁵) is always cheaper for large bases).
module mp2_lib

  use precision, only: dp
  implicit none

  private
  public :: mp2_correlation
  !> Reused by the open-shell coupled-cluster path, which needs the same
  !> semicanonical basis before its denominators are defined.
  public :: semicanonicalize

  !> Default guard on the number of per-MO-pair Coulomb builds the correlation
  !> assembly performs (nocc*nvir over both spins); overridable at run time via
  !> OQP_MP2_MAX_JBUILDS.  Prevents an accidental large number of J-builds.
  integer, parameter :: MAX_JBUILDS = 4000

contains

  !> .true. when the per-MO-pair Coulomb-build count for this reference is within
  !> MAX_JBUILDS (overridable via OQP_MP2_MAX_JBUILDS).
  logical function mp2_build_is_affordable(nbf, nocca, noccb) result(ok)
    integer, intent(in) :: nbf, nocca, noccb
    integer :: vira, virb, nbuild, cap, ln
    character(len=32) :: sval

    vira = nbf - nocca
    virb = nbf - noccb
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
    !> Second-order singles energy.  Zero for RHF and for canonical UHF;
    !> non-zero for ROHF, where it is a required part of E(2).
    real(kind=dp), intent(out) :: e_s
    logical, intent(out) :: computed

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver

    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:)
    ! Semicanonical orbitals/energies (occ-occ and vir-vir Fock blocks
    ! diagonalized) -- required for a standard ROHF/UHF MP2.
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

    nbf  = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    vira = nbf - nocca
    virb = nbf - noccb

    ! Decide whether the O(N⁵) path should be attempted.
    ! try_n5 = true when direct-J is unaffordable, or when the basis is large
    ! enough that O(N⁵) is always preferable.
    try_n5 = .not. mp2_build_is_affordable(nbf, nocca, noccb)
    if (.not. try_n5 .and. nbf > 120) try_n5 = .true.

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    restricted_ref = (infos%control%scftype == 1)
    if (restricted_ref) then
      mo_b => mo_a
      fock_b => fock_a
    else
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    end if

    ss_scale = infos%dft%MP2SS_Scale
    os_scale = infos%dft%MP2OS_Scale
    need_same_spin = abs(ss_scale) > 1.0e-14_dp
    need_opposite_spin = abs(os_scale) > 1.0e-14_dp

    ! Semicanonicalize each spin so the MP2 denominators use canonical orbital
    ! energies (Fock occ-occ / vir-vir sub-blocks diagonalized).
    allocate(mo_a_sc(nbf,nbf), mo_b_sc(nbf,nbf), e_a_sc(nbf), e_b_sc(nbf), &
             source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: cannot allocate semicanonical MOs', with_abort)
    call semicanonicalize(nbf, nocca, mo_a, fock_a, mo_a_sc, e_a_sc)
    call semicanonicalize(nbf, noccb, mo_b, fock_b, mo_b_sc, e_b_sc)

    ! Second-order singles.
    e_s = mp2_singles(nbf, nocca, mo_a_sc, fock_a, e_a_sc) &
        + mp2_singles(nbf, noccb, mo_b_sc, fock_b, e_b_sc)

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    if (try_n5) then
      ! O(N⁵) batched half-transform path.  Falls back to direct-J when the
      ! packed AO integral allocation fails.
      call mp2_corr_n5(int2_driver, basis, nbf, nbf2, &
          mo_a_sc, e_a_sc, nocca, vira, &
          mo_a_sc, e_a_sc, nocca, vira, &
          mo_b_sc, e_b_sc, noccb, virb, &
          need_same_spin, need_opposite_spin, e_aa, e_ab, &
          success=n5_ok)
      if (n5_ok .and. need_same_spin) then
        ! Beta same-spin block (opposite-spin already counted from alpha side)
        call mp2_corr_n5(int2_driver, basis, nbf, nbf2, &
            mo_b_sc, e_b_sc, noccb, virb, &
            mo_b_sc, e_b_sc, noccb, virb, &
            mo_a_sc, e_a_sc, nocca, vira, &
            same_spin=.true., do_opposite=.false., &
            e_same=e_bb, e_opp=e_opp_scratch, &
            success=n5_ok)
      end if
      if (n5_ok) computed = .true.
    end if

    if (.not. computed) then
      ! Direct-J (O(N⁶)) path: original per-(i,a) J-builds.
      if (.not. mp2_build_is_affordable(nbf, nocca, noccb)) then
        deallocate(mo_a_sc, mo_b_sc, e_a_sc, e_b_sc)
        call int2_driver%clean()
        return
      end if

      ! Same-spin alpha block + opposite-spin block.
      call mp2_spin_block(int2_driver, basis, nbf, nbf2, &
                          mo_a_sc, e_a_sc, nocca, vira, &
                          mo_a_sc, e_a_sc, nocca, vira, &
                          mo_b_sc, e_b_sc, noccb, virb, &
                          same_spin=need_same_spin, do_opposite=need_opposite_spin, &
                          e_same=e_aa, e_opp=e_ab)

      ! Same-spin beta block (opposite-spin already counted above).
      e_opp_scratch = 0.0_dp
      call mp2_spin_block(int2_driver, basis, nbf, nbf2, &
                          mo_b_sc, e_b_sc, noccb, virb, &
                          mo_b_sc, e_b_sc, noccb, virb, &
                          mo_a_sc, e_a_sc, nocca, vira, &
                          same_spin=need_same_spin, do_opposite=.false., &
                          e_same=e_bb, e_opp=e_opp_scratch)

      computed = .true.
    end if

    deallocate(mo_a_sc, mo_b_sc, e_a_sc, e_b_sc)

    call int2_driver%clean()

    ! The spin-component scales are defined for the doubles components; the
    ! singles term is neither same- nor opposite-spin, so it enters unscaled.
    e_mp2 = ss_scale * (e_aa + e_bb) + os_scale * e_ab + e_s

  end subroutine mp2_correlation


  !###########################################################################
  ! O(N⁵) batched half-transform path
  !###########################################################################

  !> @brief MP2 correlation via batched half-transform (O(N⁵)).
  !>
  !> Collects packed AO integrals once, then two half-transforms over the
  !> (occupied, virtual) MO block produce (ia|jb) without per-pair J-builds.
  !> The old O(N⁶) direct-J path has no per-pair limit and raises the
  !> effective MAX_JBUILDS to infinity.
  !>
  !> @param[in]  int2_driver  initialised two-electron driver
  !> @param[in]  basis        basis set
  !> @param[in]  nbf, nbf2    basis size and packed length
  !> @param[in]  cmo_l, e_l   left-spin semicanonical MOs and energies
  !> @param[in]  nocc_l, nvir_l  left-spin occupied/virtual counts
  !> @param[in]  cmo_s, e_s   same-spin partner MOs and energies
  !> @param[in]  nocc_s, nvir_s  same-spin occ/vir counts
  !> @param[in]  cmo_o, e_o   opposite-spin partner MOs and energies
  !> @param[in]  nocc_o, nvir_o  opposite-spin occ/vir counts
  !> @param[in]  same_spin    .true. to include same-spin aa/bb contribution
  !> @param[in]  do_opposite  .true. to include opposite-spin ab contribution
  !> @param[out] e_same       same-spin energy accumulated
  !> @param[out] e_opp        opposite-spin energy accumulated
  !> @param[out] success      .true. if the n5 path completed
  subroutine mp2_corr_n5(int2_driver, basis, nbf, nbf2, &
                         cmo_l, e_l, nocc_l, nvir_l, &
                         cmo_s, e_s, nocc_s, nvir_s, &
                         cmo_o, e_o, nocc_o, nvir_o, &
                         same_spin, do_opposite, e_same, e_opp, success)

    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use cc_ao2mo, only: cc_eri_collect_t, cc_packed_length
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
    logical, intent(out) :: success

    integer(8) :: packed_len
    integer :: npair, nov_l, nov_s, nov_o, idx, i, a, j, b, q, ip, ok
    integer :: lambda, sigma, mu, nu
    real(kind=dp) :: denom, num, val

    ! Packed AO integrals (TARGET because cc_eri_collect_t%g points here)
    real(kind=dp), allocatable, target :: g(:)

    ! Half-transform intermediate: (i a | lambda sigma) for every AO ket pair
    real(kind=dp), allocatable :: half(:,:)

    ! Precomputed pair tables
    integer, allocatable :: prow(:), pcol(:)

    ! ovov(i,a,j,b) for same-spin antisymmetrisation
    real(kind=dp), allocatable :: ovov(:,:,:,:)

    success = .false.

    npair = nbf2  ! nbf*(nbf+1)/2
    packed_len = cc_packed_length(nbf)
    nov_l = nocc_l * nvir_l
    nov_s = nocc_s * nvir_s
    nov_o = nocc_o * nvir_o

    ! --- Step 0: allocate packed AO integrals and collect -------------------
    allocate(g(packed_len), source=0.0_dp, stat=ok)
    if (ok /= 0) return  ! fall back to direct-J

    block
      type(cc_eri_collect_t), target :: collector
      collector%g => g
      collector%nbf = nbf
      collector%npair = npair
      call int2_driver%run(collector)
    end block

    ! --- Precompute AO pair tables -----------------------------------------
    allocate(prow(npair), pcol(npair), stat=ok)
    if (ok /= 0) then; deallocate(g); return; end if

    ! Compute the pair tables (trivial O(npair) work)
    do q = 1, npair
      do lambda = 1, nbf
        if (q <= lambda*(lambda-1)/2) cycle
        if (q > lambda*(lambda+1)/2) cycle
        sigma = q - lambda*(lambda-1)/2
        prow(q) = lambda
        pcol(q) = sigma
        exit
      end do
    end do

    ! --- Step 1: first half-transform (O(N⁵)) --------------------------------
    ! For each AO ket pair q = pair(lambda,sigma):
    !   m(i,a) = (i a | lambda sigma)
    !          = sum_{mu,nu} C_mu_i * C_nu_a * g(pair(mu,nu), pair(lambda,sigma))
    ! Store in half(q, idx(i,a)).
    allocate(half(npair, nov_l), source=0.0_dp, stat=ok)
    if (ok /= 0) then; deallocate(g, prow, pcol); return; end if

    !$omp parallel default(shared) private(q, ip, mu, nu, idx)
    block
      real(kind=dp), allocatable :: d(:,:), scr(:,:), m(:,:)
      allocate(d(nbf,nbf), scr(nbf,nvir_l), m(nocc_l,nvir_l))
      !$omp do schedule(static)
      do q = 1, npair
        ! Unpack g column q into d(mu,nu)
        !$omp simd
        do ip = 1, npair
          mu = prow(ip); nu = pcol(ip)
          d(mu,nu) = g(mp2_packed_idx(ip, q))
          d(nu,mu) = d(mu,nu)
        end do

        ! scr_na = sum_mu d(mu,nu) * cmo_l(mu, nocc_l + a)  (all a)
        call dgemm('n','n', nbf, nvir_l, nbf, 1.0_dp, d, nbf, &
                   cmo_l(1, nocc_l+1), nbf, 0.0_dp, scr, nbf)
        ! m(i,a) = sum_nu cmo_l(nu,i) * scr_na(nu,a)
        call dgemm('t','n', nocc_l, nvir_l, nbf, 1.0_dp, &
                   cmo_l, nbf, scr, nbf, 0.0_dp, m, nocc_l)

        ! Store into half
        do a = 1, nvir_l
          do i = 1, nocc_l
            idx = (i-1)*nvir_l + a
            half(q, idx) = m(i, a)
          end do
        end do
      end do
      !$omp end do
      deallocate(d, scr, m)
    end block
    !$omp end parallel

    ! Free AO integrals — no longer needed after first half
    deallocate(g)

    ! --- Step 2: second half-transform and energy accumulation --------------
    ! Same-spin needs the full ovov(i,a,j,b) for antisymmetrisation
    if (same_spin) then
      allocate(ovov(nocc_l, nvir_l, nocc_s, nvir_s), source=0.0_dp, stat=ok)
      if (ok /= 0) then
        deallocate(half, prow, pcol); return
      end if

      call mp2_n5_second_half(nbf, npair, nov_l, &
          nocc_l, nvir_l, nocc_s, nvir_s, &
          cmo_l, cmo_s, half, prow, pcol, ovov)

      ! Antisymmetrised accumulation
      ! E_aa = 0.25 * sum_{ijab} ((ia|jb) - (ib|ja))^2 / D
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

    ! Opposite-spin: accumulate on the fly (no antisymmetrisation needed)
    if (do_opposite) then
      call mp2_n5_opposite(nbf, npair, nov_l, &
          nocc_l, nvir_l, nocc_o, nvir_o, &
          cmo_l, e_l, cmo_o, e_o, &
          half, prow, pcol, e_opp)
    end if

    deallocate(half, prow, pcol)
    success = .true.

  end subroutine mp2_corr_n5


  !> @brief Second half-transform: (ia|jb) from (ia|lambda sigma).
  !>
  !> For each (i,a) of the left spin, reconstruct the (ia|λσ) matrix from
  !> the packed half-transform intermediate, then transform the ket indices
  !> to MO basis.  Fills ovov(i,a,j,b) = (ia|jb).
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
      i = (idx - 1) / nvir_l + 1
      a = mod(idx - 1, nvir_l) + 1

      ! Reconstruct full (ia|lambda,sigma) from packed half
      !$omp simd
      do q = 1, npair
        lambda = prow(q); sigma = pcol(q)
        e_mat(lambda, sigma) = half(q, idx)
        e_mat(sigma, lambda) = half(q, idx)
      end do

      ! Transform ket: (ia|jb)
      ! scr2(j, sigma) = sum_lambda C_s(j, lambda) * e_mat(lambda, sigma)
      call dgemm('t', 'n', nocc_s, nbf, nbf, 1.0_dp, &
                 cmo_s, nbf, e_mat, nbf, 0.0_dp, scr2, nocc_s)
      ! vv(j, b) = sum_sigma scr2(j, sigma) * C_s(sigma, nocc_s+b)
      call dgemm('n', 'n', nocc_s, nvir_s, nbf, 1.0_dp, &
                 scr2, nocc_s, cmo_s(1, nocc_s+1), nbf, 0.0_dp, vv, nocc_s)

      do b = 1, nvir_s
        do j = 1, nocc_s
          ovov(i, a, j, b) = vv(j, b)
        end do
      end do
    end do
    !$omp end parallel do

    deallocate(e_mat, scr2, vv)

  end subroutine mp2_n5_second_half


  !> @brief Opposite-spin contribution (ia|jb)^2 / denom, no antisym.
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
      i = (idx - 1) / nvir_l + 1
      a = mod(idx - 1, nvir_l) + 1

      ! Reconstruct full (ia|lambda,sigma) from packed half
      !$omp simd
      do q = 1, npair
        lambda = prow(q); sigma = pcol(q)
        e_mat(lambda, sigma) = half(q, idx)
        e_mat(sigma, lambda) = half(q, idx)
      end do

      ! Transform ket to opposite-spin MO basis
      call dgemm('t', 'n', nocc_o, nbf, nbf, 1.0_dp, &
                 cmo_o, nbf, e_mat, nbf, 0.0_dp, scr2, nocc_o)
      call dgemm('n', 'n', nocc_o, nvir_o, nbf, 1.0_dp, &
                 scr2, nocc_o, cmo_o(1, nocc_o+1), nbf, 0.0_dp, vv, nocc_o)

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


  !> @brief Triangular packed index for the symmetrised pair-pair store.
  !>
  !> For AO pair p = pair(mu,nu) and q = pair(lambda,sigma), returns
  !> g_offset such that the canonical 1/8-of-tensor value (mu nu | lambda sigma)
  !> with p >= q is stored at g(g_offset).
  pure integer(8) function mp2_packed_idx(p, q) result(idx)
    integer, intent(in) :: p, q
    integer(8) :: hi, lo
    hi = int(max(p, q), 8)
    lo = int(min(p, q), 8)
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
    ! same_block(a,b,j) = (i a | j b) for one occupied i.
    real(kind=dp), allocatable :: same_block(:,:,:)
    real(kind=dp), allocatable :: gopp(:,:)
    integer :: i, a, j, b, ii, ok
    real(kind=dp) :: denom, num

    allocate(pdmat(nbf2,2), dfull(nbf,nbf), jfull(nbf,nbf), scr(nbf,nbf), &
             jpack(nbf2), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: cannot allocate J-build scratch', WITH_ABORT)
    if (same_spin .and. (nocc_s /= nocc_l .or. nvir_s /= nvir_l)) then
      call show_message('mp2: inconsistent same-spin dimensions', WITH_ABORT)
    end if
    if (same_spin) then
      allocate(same_block(nvir_l,nvir_s,nocc_s), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('mp2: cannot allocate same-spin block', WITH_ABORT)
    end if
    if (do_opposite) then
      allocate(gopp(nocc_o,nvir_o), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('mp2: cannot allocate opp-spin block', WITH_ABORT)
    end if

    do i = 1, nocc_l
      if (same_spin) same_block = 0.0_dp
      do a = 1, nvir_l
        ! Rank-1 symmetric AO density of MO_i (x) MO_(occ+a) (left spin).
        call rank1_sym_density(cmo_l(:,i), cmo_l(:,nocc_l+a), nbf, dfull)
        call pack_matrix(dfull, pdmat(:,1), 'U')
        pdmat(:,2) = 0.0_dp

        int2_data = int2_urohf_data_t(nfocks=2, d=pdmat, scale_exchange=0.0_dp)
        call int2_driver%run(int2_data)
        ! The packed Fock accumulator stores OFF-DIAGONAL elements doubled and
        ! the DIAGONAL untouched (same convention as scf_addons::fock_jk): the
        ! true matrix is 0.5*f off-diagonal, f on the diagonal.  Recover it
        ! before unpacking so J(D) = (mu nu | i a) has the correct magnitude.
        jpack(:) = 0.5_dp * int2_data%f(:,1,1)
        ii = 0
        do j = 1, nbf
          ii = ii + j
          jpack(ii) = 2.0_dp * jpack(ii)
        end do
        call unpack_matrix(jpack, jfull, 'U')
        call int2_data%clean()

        if (same_spin) then
          ! Same-spin: (i a | j b) = MO_j^T J MO_(occ+b), all j,b of left spin.
          ! scr = J * C_occ(left)  -> (nbf, nocc_s)
          call dgemm('n','n', nbf, nocc_s, nbf, 1.0_dp, jfull, nbf, &
                     cmo_s, nbf, 0.0_dp, scr, nbf)
          do j = 1, nocc_s
            do b = 1, nvir_s
              ! (i a | j b) = sum_mu C_(occ+b)(mu) * scr(mu,j)
              same_block(a,b,j) = dot_product(cmo_s(:,nocc_s+b), scr(:,j))
            end do
          end do
        end if

        if (do_opposite) then
          ! Opposite spin: (i a | j b) with j,b on the OTHER spin.
          call dgemm('n','n', nbf, nocc_o, nbf, 1.0_dp, jfull, nbf, &
                     cmo_o, nbf, 0.0_dp, scr, nbf)
          do j = 1, nocc_o
            do b = 1, nvir_o
              gopp(j,b) = dot_product(cmo_o(:,nocc_o+b), scr(:,j))
            end do
          end do
          do j = 1, nocc_o
            do b = 1, nvir_o
              denom = e_l(i) + e_o(j) - e_l(nocc_l+a) - e_o(nocc_o+b)
              if (abs(denom) < 1.0e-10_dp) cycle
              e_opp = e_opp + gopp(j,b)*gopp(j,b) / denom
            end do
          end do
        end if
      end do
      ! Same-spin contraction with antisymmetrized integrals.
      if (same_spin) then
        do j = 1, nocc_s
          do a = 1, nvir_l
            do b = 1, nvir_s
              denom = e_l(i) + e_s(j) - e_l(nocc_l+a) - e_s(nocc_s+b)
              if (abs(denom) < 1.0e-10_dp) cycle
              ! <ij||ab> = (ia|jb) - (ib|ja)
              num = same_block(a,b,j) - same_block(b,a,j)
              e_same = e_same + 0.25_dp * num*num / denom
            end do
          end do
        end do
      end if
    end do

    deallocate(pdmat, dfull, jfull, scr, jpack)
    if (allocated(same_block)) deallocate(same_block)
    if (allocated(gopp)) deallocate(gopp)

  end subroutine mp2_spin_block

  subroutine rank1_sym_density(u, v, nbf, dfull)
    real(kind=dp), intent(in) :: u(nbf), v(nbf)
    integer, intent(in) :: nbf
    real(kind=dp), intent(out) :: dfull(nbf,nbf)
    integer :: mu, nu
    do nu = 1, nbf
      do mu = 1, nbf
        dfull(mu,nu) = 0.5_dp*(u(mu)*v(nu) + u(nu)*v(mu))
      end do
    end do
  end subroutine rank1_sym_density

!> @brief Second-order singles energy in the semicanonical basis.
!>
!>   E_S = sum_ia |f_ia|^2 / (e_i - e_a)
!>
!> @param[in] cmo_sc       semicanonical MO coefficients for this spin
!> @param[in] fock_packed  the same spin's Fock matrix, packed, in the AO basis
!> @param[in] e_sc         semicanonical orbital energies for this spin
!>
!> Returns zero whenever the occupied-virtual Fock block vanishes, which
!> covers RHF and canonical UHF; only an ROHF reference makes it non-zero.
  function mp2_singles(nbf, nocc, cmo_sc, fock_packed, e_sc) result(e_s)

    use mathlib, only: unpack_matrix

    integer, intent(in) :: nbf, nocc
    real(kind=dp), intent(in) :: cmo_sc(nbf,nbf)
    real(kind=dp), intent(in) :: fock_packed(:)
    real(kind=dp), intent(in) :: e_sc(nbf)
    real(kind=dp) :: e_s

    real(kind=dp), allocatable :: fao(:,:), scr(:,:), fmo(:,:)
    real(kind=dp) :: d
    integer :: i, a

    e_s = 0.0_dp
    if (nocc <= 0 .or. nocc >= nbf) return

    allocate(fao(nbf,nbf), scr(nbf,nbf), fmo(nbf,nbf), source=0.0_dp)

    call unpack_matrix(fock_packed, fao, 'U')
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, cmo_sc, nbf, 0.0_dp, scr, nbf)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, cmo_sc, nbf, scr, nbf, 0.0_dp, fmo, nbf)

    do a = nocc+1, nbf
      do i = 1, nocc
        d = e_sc(i) - e_sc(a)
        if (abs(d) < 1.0e-12_dp) cycle
        e_s = e_s + fmo(i,a)*fmo(i,a)/d
      end do
    end do

    deallocate(fao, scr, fmo)

  end function mp2_singles

!###############################################################################

!> Rotate the occupied and virtual blocks to the semicanonical basis, in which
!> the Fock matrix is diagonal within each of them.
!>
!> @param[in] nfzc  frozen core orbitals, excluded from the occupied rotation.
!>                  Optional, default 0 (rotate the whole occupied space).
!>
!> A frozen core has to be taken out BEFORE the rotation, not after it.  The
!> occupied-occupied rotation mixes core with valence, so dropping the first
!> nfzc columns of an already-rotated set discards a different subspace than
!> the intended core -- and, since alpha and beta are rotated separately, an
!> ROHF reference ends up freezing a different spatial orbital in each spin.
!> Restricting the rotation to the active window keeps the correlated space
!> equal to the span of the reference orbitals nfzc+1..nbf, which is what
!> every other code means by frozen-core CC.  The correlation energy is
!> invariant to the choice of basis *within* that space, so the active-window
!> rotation is free to make the denominators well defined.
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

    nfz = 0
    if (present(nfzc)) nfz = max(0, min(nfzc, nocc))
    nocc_act = nocc - nfz

    nvir = nbf - nocc
    allocate(fao(nbf,nbf), fmo(nbf,nbf), scr(nbf,nbf), eval(nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: semicanonical alloc failed', with_abort)

    ! F_mo = C^T F_ao C
    call unpack_matrix(fock_packed, fao, 'U')
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, cmo, nbf, 0.0_dp, scr, nbf)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, cmo, nbf, scr, nbf, 0.0_dp, fmo, nbf)

    cmo_sc = cmo
    e_sc = 0.0_dp

    ! Frozen core keeps its reference orbitals; report their diagonal Fock
    ! elements so e_sc is meaningful over the whole range even though the
    ! correlated methods never read this part.
    do i = 1, nfz
      e_sc(i) = fmo(i,i)
    end do

    ! --- occupied-occupied block, over the correlated window only ---
    if (nocc_act > 0) then
      allocate(blk(nocc_act,nocc_act), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('mp2: occ block alloc failed', with_abort)
      blk = fmo(nfz+1:nocc, nfz+1:nocc)
      call diag_symm_full(1, nocc_act, blk, nocc_act, eval(1:nocc_act))
      ! blk now holds eigenvectors (columns); rotate occupied MOs.
      call dgemm('n','n', nbf, nocc_act, nocc_act, 1.0_dp, cmo(:,nfz+1:nocc), nbf, &
                 blk, nocc_act, 0.0_dp, cmo_sc(:,nfz+1:nocc), nbf)
      e_sc(nfz+1:nocc) = eval(1:nocc_act)
      deallocate(blk)
    end if

    ! --- virtual-virtual block ---
    if (nvir > 0) then
      allocate(blk(nvir,nvir), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('mp2: vir block alloc failed', with_abort)
      blk = fmo(nocc+1:nbf, nocc+1:nbf)
      call diag_symm_full(1, nvir, blk, nvir, eval(1:nvir))
      call dgemm('n','n', nbf, nvir, nvir, 1.0_dp, cmo(:,nocc+1:nbf), nbf, &
                 blk, nvir, 0.0_dp, cmo_sc(:,nocc+1:nbf), nbf)
      e_sc(nocc+1:nbf) = eval(1:nvir)
      deallocate(blk)
    end if

    deallocate(fao, fmo, scr, eval)

  end subroutine semicanonicalize

end module mp2_lib
