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
!$ use omp_lib
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  ! Radix width of the occupation-signature sort in pt2_occupation_blocks.
  integer, parameter :: rbits = 8
  integer, parameter :: rbuck = 2**rbits
  integer(i8), parameter :: rmask = int(rbuck - 1, i8)

  public :: pt2_effective_fock, pt2_h0_dyall_active, pt2_external_indices
  public :: pt2_diagonal_h0, pt2_occupation_blocks, pt2_minor_transform

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
    !
    ! This nest is the only O(nbf^4) work in the routine and each `p` owns its
    ! own column of `kmat`, so the columns are independent and nothing has to be
    ! reduced across threads.  Splitting on `p` also keeps each thread on its own
    ! contiguous nbf^3 slab of `eri`.  Every kmat(:,p) still accumulates over
    ! (r,s) in the same ascending order as the serial nest, so the result does
    ! not depend on the thread count.
    kmat = 0.0_dp
    !$omp parallel do default(shared) schedule(static) if(int(n, i8)**3 >= 32768_i8) &
    !$omp   private(p, r, s, q, drs, base_pr)
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
    !$omp end parallel do

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
    integer(i8) :: k, kb, khi, det, bit, abit
    real(dp) :: e

    ! Determinants per thread chunk.  Large enough that the inner orbital sweep
    ! still vectorises over a long run, small enough that a chunk's slice of
    ! `diag` stays in L1 across all `no` passes over it.
    integer(i8), parameter :: dblock = 4096_i8

    no = int(norb)

    ! Orbital-major, determinant-minor: the same loop order NumPy uses, so each
    ! diag(k) still accumulates over p in ascending order and the result is
    ! bit-identical -- but the inner loop is a long unit-stride pass the
    ! compiler can vectorise, instead of a short branchy one per determinant.
    !
    ! Blocking the determinant axis and handing the blocks to OpenMP keeps that
    ! order exactly (within a block the orbital loop is still the outer one and
    ! still ascends), while making the sweep parallel: no diag(k) is ever
    ! touched by two threads, so there is no reduction and no dependence of the
    ! answer on the thread count.
    !$omp parallel do default(shared) schedule(static) if(ndet >= 8192_i8) &
    !$omp   private(kb, khi, p, k, det, bit, abit, e)
    do kb = 0_i8, ndet - 1_i8, dblock
      khi = min(ndet - 1_i8, kb + dblock - 1_i8)
      do k = kb, khi
        diag(k) = 0.0_dp
      end do
      do p = 0, no - 1
        bit = ishft(1_i8, p)
        abit = ishft(bit, no)
        e = eps(p)
        do k = kb, khi
          det = dets(k)
          if (iand(det, bit) /= 0_i8) diag(k) = diag(k) + e
          if (iand(det, abit) /= 0_i8) diag(k) = diag(k) + e
        end do
      end do
    end do
    !$omp end parallel do
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
  !> The sort is an LSD radix sort on the signature, in 8-bit digits.  Radix is
  !> the right shape here twice over: the signature is a small non-negative
  !> integer (it is `det` masked down to the inactive+virtual bits, so it never
  !> exceeds 2^62 and is never negative), and a counting pass is inherently
  !> STABLE -- equal signatures come out in ascending original position, which is
  !> exactly what numpy's `argsort(kind="stable")` produces.  The earlier
  !> revision used a shell sort and had to compare the composite key
  !> (signature, original position) to recover that order; stability is now a
  !> property of the algorithm rather than of the comparator, but the requirement
  !> behind it is unchanged and is the reason radix was chosen over an
  !> introsort: the intra-block member order feeds a per-block eigendecomposition,
  !> so any reordering inside a block would perturb the printed energies.
  !>
  !> All the digit histograms are built in one sweep, and a pass whose digit is
  !> constant over the whole array is skipped -- which most of them are, since
  !> `keep` masks out the whole active window and the signature is sparse.  The
  !> caller still measures that the off-block Hamiltonian elements vanish before
  !> trusting the partition.
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

    integer :: no, nc, na, i, d, nbit, npass, ip, nz
    integer(i8) :: k, act_mask, frozen, keep
    logical :: staged
    integer(i8), allocatable :: sig(:), sig2(:), ord2(:), cnt(:,:)

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

    allocate(sig(0:next - 1), sig2(0:next - 1), ord2(0:next - 1))
    do k = 0_i8, next - 1_i8
      sig(k) = iand(dets(ext(k)), keep)
      order(k) = k
    end do

    ! Only the digits a signature can actually reach are worth a pass, and
    ! `keep` bounds those exactly.
    nbit = 0
    do i = 0, 62
      if (btest(keep, i)) nbit = i + 1
    end do
    npass = max(1, (nbit + rbits - 1) / rbits)

    allocate(cnt(0:rbuck - 1, 0:npass - 1))
    cnt = 0_i8
    do k = 0_i8, next - 1_i8
      do ip = 0, npass - 1
        d = int(iand(ishft(sig(k), -ip * rbits), rmask))
        cnt(d, ip) = cnt(d, ip) + 1_i8
      end do
    end do

    ! `staged` tracks which pair of buffers currently holds the data, so the
    ! passes ping-pong without copying between them.
    staged = .false.
    do ip = 0, npass - 1
      nz = 0
      do d = 0, rbuck - 1
        if (cnt(d, ip) /= 0_i8) nz = nz + 1
      end do
      if (nz <= 1) cycle              ! digit constant: the pass is a no-op
      if (staged) then
        call radix_pass(next, ip * rbits, cnt(0, ip), sig2, ord2, sig, order)
      else
        call radix_pass(next, ip * rbits, cnt(0, ip), sig, order, sig2, ord2)
      end if
      staged = .not. staged
    end do
    if (staged) then
      do k = 0_i8, next - 1_i8
        sig(k) = sig2(k)
        order(k) = ord2(k)
      end do
    end if
    deallocate(sig2, ord2, cnt)

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

  !> The single-spin orbital-rotation transform on the string space:
  !> T[t,s] = det( R[t_occ, s_occ] ), the Loewdin minors of `R`.
  !>
  !> This is the MS/XMS rotation's inner loop and it was the most expensive
  !> Python left anywhere on the PT2 path -- not because the arithmetic is
  !> large (nstr^2 determinants of an nelec x nelec matrix is ~4 MFlop at
  !> norb=10) but because it was `np.linalg.det` called once per string PAIR:
  !> 63504 LAPACK round trips through the interpreter for 4 MFlop of work.
  !>
  !> The determinant is therefore evaluated inline rather than by calling
  !> LAPACK per minor, and the elimination reproduces `dgetf2` step for step --
  !> IDAMAX pivoting taking the FIRST row of maximal magnitude, the row swap,
  !> the reciprocal-multiply scaling (with dgetf2's own fallback to division
  !> below the safe minimum), the right-looking rank-1 update, and finally the
  !> diagonal product accumulated in ascending k.  That is the same
  !> factorization `np.linalg.det` gets from `dgetrf`, which at these sizes
  !> dispatches straight to the unblocked `dgetf2`.  Measured against numpy on
  !> random unitary R the two agree to one ulp (max 1.7e-16 absolute on O(1)
  !> minors at norb=10) rather than bit for bit -- close enough that nothing
  !> printed moves, but not identical, so do not rely on exact equality.  The
  !> zero-pivot exit IS exact on both sides: a structurally singular minor
  !> returns a bit-zero determinant.
  !>
  !> The `nelec` rows of `R` a given `t` selects are gathered once and reused
  !> across all `nstr` values of `s`, so the strided read of `R` happens
  !> nstr times rather than nstr^2.
  !>
  !> @param[in]  norb   spatial orbitals
  !> @param[in]  nelec  electrons of this spin (the minor is nelec x nelec)
  !> @param[in]  nstr   number of single-spin strings
  !> @param[in]  rmat   the rotation, C-order [norb,norb]
  !> @param[in]  occ    occupied orbital indices per string, C-order [nstr,nelec]
  !> @param[out] tmat   C-order [nstr,nstr]
  subroutine pt2_minor_transform(norb, nelec, nstr, rmat, occ, tmat) &
      bind(C, name="pt2_minor_transform")
    integer(c_int32_t), value :: norb, nelec
    integer(i8), value :: nstr
    real(dp), intent(in) :: rmat(0:*)
    integer(i8), intent(in) :: occ(0:*)
    real(dp), intent(inout) :: tmat(0:*)

    integer :: no, ne, x, y, k, ipiv, ib
    integer(i8) :: t, s, tb, sb, ob
    real(dp) :: amax, piv, det, sgn, v, sfmin
    real(dp), allocatable :: a(:,:), rows(:,:)

    no = int(norb); ne = int(nelec)
    if (nstr <= 0_i8) return

    ! The 0 x 0 minor has determinant 1, which is what numpy returns too.
    if (ne <= 0) then
      do t = 0_i8, nstr * nstr - 1_i8
        tmat(t) = 1.0_dp
      end do
      return
    end if

    sfmin = tiny(1.0_dp)                 ! dlamch('S')

    ! nstr^2 independent minors, and `t` selects a disjoint row band of `tmat`
    ! (tmat(t*nstr : t*nstr + nstr - 1)), so the string pairs parallelise with
    ! no sharing at all -- the factorisation of one minor never reads another's.
    ! The elimination scratch has to be per-thread, so it is allocated INSIDE
    ! the region: `private` on an allocatable hands each thread its own
    ! unallocated copy, which is exactly what is wanted here, and allocating
    ! once per thread rather than once per minor keeps it off the inner path.
    !$omp parallel default(shared) if(nstr * nstr >= 4096_i8) &
    !$omp   private(a, rows, t, s, tb, sb, ob, x, y, k, ipiv, ib, amax, piv, &
    !$omp           det, sgn, v)
    allocate(a(0:ne - 1, 0:ne - 1), rows(0:ne - 1, 0:no - 1))
    !$omp do schedule(static)
    do t = 0_i8, nstr - 1_i8
      tb = t * int(ne, i8)
      do x = 0, ne - 1
        ib = int(occ(tb + int(x, i8)))
        do y = 0, no - 1
          rows(x, y) = rmat(int(ib, i8) * int(no, i8) + int(y, i8))
        end do
      end do

      ob = t * nstr
      do s = 0_i8, nstr - 1_i8
        sb = s * int(ne, i8)
        do y = 0, ne - 1
          ib = int(occ(sb + int(y, i8)))
          do x = 0, ne - 1
            a(x, y) = rows(x, ib)
          end do
        end do

        sgn = 1.0_dp
        det = 1.0_dp
        do k = 0, ne - 1
          ipiv = k                        ! idamax: first row of maximal |.|
          amax = abs(a(k, k))
          do x = k + 1, ne - 1
            if (abs(a(x, k)) > amax) then
              amax = abs(a(x, k))
              ipiv = x
            end if
          end do
          piv = a(ipiv, k)
          if (piv == 0.0_dp) then
            det = 0.0_dp
            exit
          end if
          if (ipiv /= k) then
            do y = 0, ne - 1
              v = a(k, y); a(k, y) = a(ipiv, y); a(ipiv, y) = v
            end do
            sgn = -sgn
          end if
          if (abs(piv) >= sfmin) then
            v = 1.0_dp / piv
            do x = k + 1, ne - 1
              a(x, k) = a(x, k) * v
            end do
          else
            do x = k + 1, ne - 1
              a(x, k) = a(x, k) / piv
            end do
          end if
          do y = k + 1, ne - 1
            v = a(k, y)
            do x = k + 1, ne - 1
              a(x, y) = a(x, y) - a(x, k) * v
            end do
          end do
          det = det * piv
        end do
        tmat(ob + s) = sgn * det
      end do
    end do
    !$omp end do
    deallocate(a, rows)
    !$omp end parallel
  end subroutine pt2_minor_transform

  !> One stable LSD counting pass: distribute (key, val) into `dst` by the
  !> `rbits`-wide digit of the key at `shift`, using the pre-built histogram.
  !>
  !> Stability is the whole point (see pt2_occupation_blocks): the source is
  !> walked in ascending order and each bucket cursor only advances, so records
  !> sharing a digit keep their relative order.
  subroutine radix_pass(n, shift, cnt, src_key, src_val, dst_key, dst_val)
    integer(i8), value :: n
    integer, value :: shift
    integer(i8), intent(in) :: cnt(0:*), src_key(0:*), src_val(0:*)
    integer(i8), intent(inout) :: dst_key(0:*), dst_val(0:*)

    integer(i8) :: off(0:rbuck - 1), k, pos
    integer :: d

    off(0) = 0_i8
    do d = 1, rbuck - 1
      off(d) = off(d - 1) + cnt(d - 1)
    end do
    do k = 0_i8, n - 1_i8
      d = int(iand(ishft(src_key(k), -shift), rmask))
      pos = off(d)
      off(d) = pos + 1_i8
      dst_key(pos) = src_key(k)
      dst_val(pos) = src_val(k)
    end do
  end subroutine radix_pass

end module pt2_kernel_mod
