!> @file mp2_lib.F90
!>
!> @brief Standalone Moller-Plesset second-order (MP2) ground-state correlation
!>        energy for RHF/UHF/ROHF references.
!>
!> The correlation energy is built in the spin-blocked (aa, bb, ab) form on
!> semicanonicalized orbitals, reusing the validated two-electron driver
!> (`int2_compute`) via per-occupied-MO-pair Coulomb builds -- so no full O(N^4)
!> MO integral tensor is ever stored.  A ROHF reference is semicanonicalized first
!> (occ-occ and vir-vir Fock blocks diagonalized) so the canonical MP2 amplitude
!> denominators are well defined.  Validated to 1e-8 Ha against PySCF UMP2.
module mp2_lib

  use precision, only: dp
  implicit none

  private
  public :: mp2_correlation

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

  subroutine mp2_correlation(infos, e_mp2, e_aa, e_bb, e_ab, computed)

    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use messages, only: show_message, with_abort
    use oqp_tagarray_driver, only: tagarray_get_data, &
                                   OQP_VEC_MO_A, OQP_VEC_MO_B, &
                                   OQP_FOCK_A, OQP_FOCK_B

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: e_mp2, e_aa, e_bb, e_ab
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

    e_mp2 = 0.0_dp; e_aa = 0.0_dp; e_bb = 0.0_dp; e_ab = 0.0_dp
    computed = .false.

    basis => infos%basis
    basis%atoms => infos%atoms

    nbf  = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    vira = nbf - nocca
    virb = nbf - noccb

    if (.not. mp2_build_is_affordable(nbf, nocca, noccb)) return

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
    ! energies (Fock occ-occ / vir-vir sub-blocks diagonalized).  For a UHF
    ! reference this is a no-op (orbitals already canonical); for ROHF it yields
    ! the standard ROHF-MP2 amplitudes.
    allocate(mo_a_sc(nbf,nbf), mo_b_sc(nbf,nbf), e_a_sc(nbf), e_b_sc(nbf), &
             source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: cannot allocate semicanonical MOs', with_abort)
    call semicanonicalize(nbf, nocca, mo_a, fock_a, mo_a_sc, e_a_sc)
    call semicanonicalize(nbf, noccb, mo_b, fock_b, mo_b_sc, e_b_sc)

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    ! Same-spin alpha block + opposite-spin block share the alpha (i,a) Coulomb
    ! builds, so they are accumulated together.
    call mp2_spin_block(int2_driver, basis, nbf, nbf2, &
                        mo_a_sc, e_a_sc, nocca, vira, &    ! "left"  = alpha occ/vir
                        mo_a_sc, e_a_sc, nocca, vira, &    ! same-spin partner = alpha
                        mo_b_sc, e_b_sc, noccb, virb, &    ! opposite-spin partner = beta
                        same_spin=need_same_spin, do_opposite=need_opposite_spin, &
                        e_same=e_aa, e_opp=e_ab)

    ! Same-spin beta block (opposite-spin already counted once above).
    e_opp_scratch = 0.0_dp
    call mp2_spin_block(int2_driver, basis, nbf, nbf2, &
                        mo_b_sc, e_b_sc, noccb, virb, &
                        mo_b_sc, e_b_sc, noccb, virb, &
                        mo_a_sc, e_a_sc, nocca, vira, &
                        same_spin=need_same_spin, do_opposite=.false., &
                        e_same=e_bb, e_opp=e_opp_scratch)

    deallocate(mo_a_sc, mo_b_sc, e_a_sc, e_b_sc)

    call int2_driver%clean()

    e_mp2 = ss_scale * (e_aa + e_bb) + os_scale * e_ab
    computed = .true.

  end subroutine mp2_correlation

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
      ! Same-spin contraction with antisymmetrized integrals.  Only the
      ! virtual-virtual blocks for the current occupied i are retained, avoiding
      ! the previous (nocc*nvir)^2 tensor.
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

  subroutine semicanonicalize(nbf, nocc, cmo, fock_packed, cmo_sc, e_sc)

    use mathlib, only: unpack_matrix
    use eigen, only: diag_symm_full
    use messages, only: show_message, with_abort

    integer, intent(in) :: nbf, nocc
    real(kind=dp), intent(in) :: cmo(nbf,nbf)
    real(kind=dp), intent(in) :: fock_packed(:)
    real(kind=dp), intent(out) :: cmo_sc(nbf,nbf)
    real(kind=dp), intent(out) :: e_sc(nbf)

    real(kind=dp), allocatable :: fao(:,:), fmo(:,:), scr(:,:), blk(:,:), eval(:)
    integer :: nvir, ok

    nvir = nbf - nocc
    allocate(fao(nbf,nbf), fmo(nbf,nbf), scr(nbf,nbf), eval(nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('mp2: semicanonical alloc failed', with_abort)

    ! F_mo = C^T F_ao C
    call unpack_matrix(fock_packed, fao, 'U')
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, cmo, nbf, 0.0_dp, scr, nbf)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, cmo, nbf, scr, nbf, 0.0_dp, fmo, nbf)

    cmo_sc = cmo
    e_sc = 0.0_dp

    ! --- occupied-occupied block ---
    if (nocc > 0) then
      allocate(blk(nocc,nocc), source=0.0_dp, stat=ok)
      if (ok /= 0) call show_message('mp2: occ block alloc failed', with_abort)
      blk = fmo(1:nocc,1:nocc)
      call diag_symm_full(1, nocc, blk, nocc, eval(1:nocc))
      ! blk now holds eigenvectors (columns); rotate occupied MOs.
      call dgemm('n','n', nbf, nocc, nocc, 1.0_dp, cmo(:,1:nocc), nbf, &
                 blk, nocc, 0.0_dp, cmo_sc(:,1:nocc), nbf)
      e_sc(1:nocc) = eval(1:nocc)
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
