!> @brief Fortran engine for the determinant-space bookkeeping and the mean-field
!> Fock build of pyoqp caspt2_dyall.py (CASPT2 / NEVPT2 / MRMP2 / MCQDPT2).
!>
!> pyoqp is a driver; liboqp computes.  These are the pieces of the PT2 path that
!> were still evaluated in Python or NumPy after the QDPT2 direct engine
!> (qdpt2_kernel.F90), the determinant Hamiltonian (fci_hamiltonian.F90) and the
!> RDMs (rdm_kernel.F90) had been moved down.  Every one of them keeps its Python
!> implementation in caspt2_dyall.py as the numerical pin and as the fallback for
!> when this symbol is absent, so the dense reference path is unchanged.
!>
!>  * pt2_effective_fock    -- F = h + J(D) - 1/2 K(D), the closed+active
!>                             generalized Fock whose diagonal is the
!>                             semicanonical orbital energy `eps`.
!>  * pt2_h0_dyall_active   -- the core-dressed active one-electron block of
!>                             Dyall's H0, h_tu + sum_i [2(tu|ii) - (ti|iu)].
!>  * pt2_external_indices  -- the external (first-order) determinant positions:
!>                             every determinant that is not a pure CAS
!>                             configuration (core doubly occupied, no virtual).
!>  * pt2_diagonal_h0       -- E0(D) = sum_p eps_p n_p(D), the determinant-basis
!>                             diagonal of a one-electron H0 that is already
!>                             diagonal in the MO basis (the `h0=fock` zeroth
!>                             order).  The caller has already established that
!>                             H0 really is diagonal; this only evaluates it.
!>  * pt2_occupation_blocks -- the H0-invariant core+virtual occupation blocks
!>                             over which Dyall's H0 is exactly block-diagonal.
!>                             Returns a permutation of the external space in
!>                             ascending occupation-signature order together with
!>                             the block boundaries; the caller still verifies
!>                             numerically that the partition closes before using
!>                             it, exactly as before.
!>
!> Determinant keys are the fci.py integers: alpha bits 0..norb-1, beta bits
!> norb..2*norb-1, occupied spatial orbitals as set bits.  This needs
!> 2*norb <= 62, which the callers check.
!>
!> All arrays cross the boundary in the caller's C order (last index fastest).
!> h1e, D, F and the ERI tensor are symmetric under the index exchanges that
!> distinguish the two layouts -- (pq|rs) = (qp|rs) = (pq|sr) = (rs|pq) -- so the
!> C-order buffer and the Fortran column-major view are the same tensor and no
!> repacking is needed anywhere below.
module pt2_kernel_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: pt2_effective_fock, pt2_h0_dyall_active, pt2_external_indices
  public :: pt2_diagonal_h0, pt2_occupation_blocks

contains

  !> Closed+active mean-field (generalized) Fock  F = h + J - 1/2 K, with
  !>     J[p,q] = sum_rs D[r,s] (pq|rs),   K[p,q] = sum_rs D[r,s] (pr|sq).
  !>
  !> J is a single DGEMV against the ERI tensor viewed as an (nbf^2 x nbf^2)
  !> matrix.  K contracts the two *middle* Fortran axes, so it is written as an
  !> explicit loop nest ordered to walk `eri` contiguously on the innermost
  !> index; that is the same O(nbf^4) work the NumPy einsum does, without the
  !> nbf^4 temporary einsum materializes for the transposed operand.
  !>
  !> @param[in]  nbf   number of orbitals
  !> @param[in]  h1e   MO core Hamiltonian, C-order [nbf,nbf]
  !> @param[in]  eri   MO ERIs (chemist (pq|rs)), C-order [nbf,nbf,nbf,nbf]
  !> @param[in]  dmat  spin-summed 1-RDM, C-order [nbf,nbf]
  !> @param[out] fock  generalized Fock, C-order [nbf,nbf]
  subroutine pt2_effective_fock(nbf, h1e, eri, dmat, fock) &
      bind(C, name="pt2_effective_fock")
    integer(c_int32_t), value :: nbf
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: h1e(0:*), eri(0:*), dmat(0:*)
    real(dp), intent(inout) :: fock(0:*)

    ! Local copies in the default (ILP64) integer kind the linked BLAS uses.
    integer :: n, n2, p, q, r, s
    integer(i8) :: base_pr, off
    real(dp) :: drs
    real(dp), allocatable :: jmat(:), kmat(:,:)

    n = int(nbf)
    if (n <= 0) return
    n2 = n * n

    allocate(jmat(0:n2 - 1))
    allocate(kmat(0:n - 1, 0:n - 1))

    ! J(q,p) = sum_(s,r) D(s,r) eri(s,r,q,p): the leading n2 axes are contracted,
    ! so this is y = A^T x with A the (n2 x n2) view of eri.
    call dgemv('T', n2, n2, 1.0_dp, eri, n2, dmat, 1, 0.0_dp, jmat, 1)

    ! K(q,p) = sum_(r,s) D(s,r) eri(q,s,r,p).  Innermost loop over q runs along
    ! the fastest Fortran axis of both eri(:,s,r,p) and kmat(:,p).
    kmat = 0.0_dp
    do p = 0, n - 1
      do r = 0, n - 1
        do s = 0, n - 1
          ! D[r,s], not D[s,r]: the two coincide for every density this is
          ! called with, but the einsum being reproduced reads D[r,s].
          drs = dmat(int(r, i8) * int(n, i8) + int(s, i8))
          if (drs == 0.0_dp) cycle
          ! eri(0,s,r,p) is at the flat offset ((p*n + r)*n + s)*n
          base_pr = ((int(p, i8) * int(n, i8) + int(r, i8)) * int(n, i8) &
                     + int(s, i8)) * int(n, i8)
          do q = 0, n - 1
            kmat(q, p) = kmat(q, p) + drs * eri(base_pr + int(q, i8))
          end do
        end do
      end do
    end do

    do p = 0, n - 1
      do q = 0, n - 1
        off = int(p, i8) * int(n, i8) + int(q, i8)
        fock(off) = h1e(off) + jmat(int(off)) - 0.5_dp * kmat(q, p)
      end do
    end do

    deallocate(jmat, kmat)
  end subroutine pt2_effective_fock

  !> Core-dressed active one-electron block of Dyall's zeroth-order Hamiltonian:
  !>     hact[t,u] = h1e[T,U] + sum_{i<ncore} [ 2 (TU|ii) - (Ti|iU) ],
  !> with T = ncore + t.  Replaces the Python triple loop over (t,u,i).
  !>
  !> @param[in]  nbf    number of orbitals
  !> @param[in]  ncore  inactive (doubly occupied) orbitals
  !> @param[in]  nact   active orbitals
  !> @param[in]  h1e    MO core Hamiltonian, C-order [nbf,nbf]
  !> @param[in]  eri    MO ERIs (chemist (pq|rs)), C-order [nbf,nbf,nbf,nbf]
  !> @param[out] hact   dressed active block, C-order [nact,nact]
  subroutine pt2_h0_dyall_active(nbf, ncore, nact, h1e, eri, hact) &
      bind(C, name="pt2_h0_dyall_active")
    integer(c_int32_t), value :: nbf, ncore, nact
    real(dp), intent(in) :: h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: hact(0:*)

    integer :: n, nc, na, t, u, i, tt, uu
    integer(i8) :: n_i, acc_off
    real(dp) :: acc

    n = int(nbf); nc = int(ncore); na = int(nact)
    if (na <= 0) return
    n_i = int(n, i8)

    do t = 0, na - 1
      tt = nc + t
      do u = 0, na - 1
        uu = nc + u
        acc = h1e(int(tt, i8) * n_i + int(uu, i8))
        do i = 0, nc - 1
          ! (TU|ii) = eri_C[tt,uu,i,i]; (Ti|iU) = eri_C[tt,i,i,uu]
          acc = acc + 2.0_dp * eri(((int(tt, i8) * n_i + int(uu, i8)) * n_i &
                                    + int(i, i8)) * n_i + int(i, i8)) &
                    - eri(((int(tt, i8) * n_i + int(i, i8)) * n_i &
                           + int(i, i8)) * n_i + int(uu, i8))
        end do
        acc_off = int(t, i8) * int(na, i8) + int(u, i8)
        hact(acc_off) = acc
      end do
    end do
  end subroutine pt2_h0_dyall_active

  !> Positions of the external (first-order) determinants: every determinant that
  !> is NOT a pure CAS configuration, i.e. that does not have every core orbital
  !> doubly occupied and every virtual orbital empty, in both spins.
  !>
  !> @param[in]  norb   spatial orbitals (needs 2*norb <= 62)
  !> @param[in]  ncore  inactive orbitals
  !> @param[in]  nact   active orbitals
  !> @param[in]  ndet   number of determinants
  !> @param[in]  dets   determinant keys, [ndet]
  !> @param[out] ext    external positions in ascending order, [ndet] capacity
  !> @return            the number of external determinants written to `ext`
  function pt2_external_indices(norb, ncore, nact, ndet, dets, ext) result(next) &
      bind(C, name="pt2_external_indices")
    integer(c_int32_t), value :: norb, ncore, nact
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(0:*)
    integer(i8), intent(inout) :: ext(0:*)
    integer(i8) :: next

    integer :: no, nc, na, i
    integer(i8) :: k, det, a, b, core_mask, virt_mask, full_mask

    no = int(norb); nc = int(ncore); na = int(nact)
    next = 0_i8

    core_mask = 0_i8
    do i = 0, nc - 1
      core_mask = ior(core_mask, ishft(1_i8, i))
    end do
    virt_mask = 0_i8
    do i = nc + na, no - 1
      virt_mask = ior(virt_mask, ishft(1_i8, i))
    end do
    full_mask = ishft(1_i8, no) - 1_i8

    do k = 0_i8, ndet - 1_i8
      det = dets(k)
      a = iand(det, full_mask)
      b = ishft(det, -no)
      if (iand(a, core_mask) == core_mask .and. &
          iand(b, core_mask) == core_mask .and. &
          iand(a, virt_mask) == 0_i8 .and. &
          iand(b, virt_mask) == 0_i8) cycle
      ext(next) = k
      next = next + 1_i8
    end do
  end function pt2_external_indices

  !> Determinant-basis diagonal of a one-electron H0 that is diagonal in the MO
  !> basis:  E0(D) = sum_p eps_p n_p(D), n_p = (alpha occ) + (beta occ).
  !>
  !> The caller is responsible for having established that H0 is genuinely
  !> diagonal (no two-electron part, no off-diagonal one-electron part); this
  !> routine only evaluates the resulting diagonal.
  !>
  !> @param[in]  norb  spatial orbitals (needs 2*norb <= 62)
  !> @param[in]  ndet  number of determinants
  !> @param[in]  dets  determinant keys, [ndet]
  !> @param[in]  eps   MO-diagonal of H0, [norb]
  !> @param[out] diag  determinant diagonal, [ndet]
  subroutine pt2_diagonal_h0(norb, ndet, dets, eps, diag) &
      bind(C, name="pt2_diagonal_h0")
    integer(c_int32_t), value :: norb
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(0:*)
    real(dp), intent(in) :: eps(0:*)
    real(dp), intent(inout) :: diag(0:*)

    integer :: no, p
    integer(i8) :: k, det, bit, abit
    real(dp) :: e

    no = int(norb)

    do k = 0_i8, ndet - 1_i8
      diag(k) = 0.0_dp
    end do

    ! Orbital-major, determinant-minor: the same loop order NumPy uses, so each
    ! diag(k) still accumulates over p in ascending order and the result is
    ! bit-identical -- but the inner loop is now a long unit-stride pass the
    ! compiler can vectorise, instead of a short branchy one per determinant.
    do p = 0, no - 1
      bit = ishft(1_i8, p)
      abit = ishft(bit, no)
      e = eps(p)
      do k = 0_i8, ndet - 1_i8
        det = dets(k)
        if (iand(det, bit) /= 0_i8) diag(k) = diag(k) + e
        if (iand(det, abit) /= 0_i8) diag(k) = diag(k) + e
      end do
    end do
  end subroutine pt2_diagonal_h0

  !> Partition the external determinants into the core+virtual occupation blocks
  !> that Dyall's H0 conserves.
  !>
  !> Dyall H0 is sum_p eps_p n_p on the inactive/virtual orbitals plus an operator
  !> that lives entirely inside the active window, so it cannot move an electron
  !> in or out of a core or virtual orbital: the full core+virtual occupation
  !> pattern (per spin) is invariant, and determinants sharing it form a closed
  !> block.  Determinants are sorted by that signature; `order` comes back as the
  !> permutation of external-space positions in ascending signature order and
  !> `starts` holds the block boundaries, so block b is
  !> order(starts(b) : starts(b+1)-1).
  !>
  !> The sort is a shell sort carrying the permutation, matching the stable-order
  !> idiom of rdm_kernel.F90: no scratch allocation beyond the signature copy, and
  !> the block order is deterministic run to run.  The caller still measures that
  !> the off-block Hamiltonian elements vanish before trusting the partition.
  !>
  !> @param[in]  norb    spatial orbitals (needs 2*norb <= 62)
  !> @param[in]  ncore   inactive orbitals
  !> @param[in]  nact    active orbitals
  !> @param[in]  next    number of external determinants
  !> @param[in]  dets    determinant keys of the FULL space, [ndet]
  !> @param[in]  ext     external positions into `dets`, [next]
  !> @param[out] order   permutation of 0..next-1 in signature order, [next]
  !> @param[out] starts  block start offsets, [next+1] capacity
  !> @return             the number of blocks
  function pt2_occupation_blocks(norb, ncore, nact, next, dets, ext, order, starts) &
      result(nblock) bind(C, name="pt2_occupation_blocks")
    integer(c_int32_t), value :: norb, ncore, nact
    integer(i8), value :: next
    integer(i8), intent(in) :: dets(0:*), ext(0:*)
    integer(i8), intent(inout) :: order(0:*), starts(0:*)
    integer(i8) :: nblock

    integer :: no, nc, na, i
    integer(i8) :: k, act_mask, frozen, keep, gap, j, skey, spos
    integer(i8), allocatable :: sig(:)

    no = int(norb); nc = int(ncore); na = int(nact)
    nblock = 0_i8
    if (next <= 0_i8) then
      starts(0) = 0_i8
      return
    end if

    act_mask = 0_i8
    do i = nc, nc + na - 1
      act_mask = ior(act_mask, ishft(1_i8, i))
    end do
    frozen = iand(ishft(1_i8, no) - 1_i8, not(act_mask))
    keep = ior(frozen, ishft(frozen, no))          ! both spins, one mask

    allocate(sig(0:next - 1))
    do k = 0_i8, next - 1_i8
      sig(k) = iand(dets(ext(k)), keep)
      order(k) = k
    end do

    ! Shell sort of `sig` carrying `order` alongside.  The comparison is on the
    ! COMPOSITE key (signature, original position), which is a total order
    ! because the positions are distinct.  That makes the result identical to
    ! numpy's stable argsort even though shell sort is not itself stable, so the
    ! intra-block member order -- and with it the per-block eigenproblem the
    ! caller then solves -- is bit-for-bit what the Python fallback produces.
    gap = next / 2_i8
    do while (gap > 0_i8)
      do k = gap, next - 1_i8
        skey = sig(k)
        spos = order(k)
        j = k
        do while (j >= gap)
          if (sig(j - gap) < skey) exit
          if (sig(j - gap) == skey .and. order(j - gap) < spos) exit
          sig(j) = sig(j - gap)
          order(j) = order(j - gap)
          j = j - gap
        end do
        sig(j) = skey
        order(j) = spos
      end do
      gap = gap / 2_i8
    end do

    starts(0) = 0_i8
    nblock = 1_i8
    do k = 1_i8, next - 1_i8
      if (sig(k) /= sig(k - 1_i8)) then
        starts(nblock) = k
        nblock = nblock + 1_i8
      end if
    end do
    starts(nblock) = next

    deallocate(sig)
  end function pt2_occupation_blocks

end module pt2_kernel_mod
