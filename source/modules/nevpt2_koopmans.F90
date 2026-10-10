!> @brief Fortran engine for the dominant Koopmans intermediates of pyoqp
!> nevpt2_sc.py (strongly contracted NEVPT2, Dyall zeroth order).
!>
!> pyoqp is a driver; liboqp computes.  Every routine here keeps its Python
!> implementation in nevpt2_sc.py as the numerical pin and as the fallback for
!> when the symbol is absent.
!>
!> Measured on synthetic operands of the shape the module actually builds
!> (best-of-5, engine and NumPy interleaved per kernel so a drifting machine
!> load hits both sides of each ratio equally -- absolute numbers on a shared
!> machine are not comparable across runs, the ratios are):
!>
!>                      nact=4        nact=6        nact=8      NumPy at 8
!>     f3ca/f3ac         1.6x          1.0x          2.1x         35.3 ms
!>     a16               7.8x          4.0x          2.1x          6.4 ms
!>     a22               6.8x          3.8x          1.8x          7.4 ms
!>     hdm3             15.5x         33.3x         13.3x          3.7 ms
!>     a7                7.9x          5.3x          3.1x          1.6 ms
!>     a13               2.4x          3.1x          1.5x          0.5 ms
!>     a9                5.8x          2.1x          1.9x          0.5 ms
!>     a12               1.9x          1.8x          1.4x          0.2 ms
!>
!> The ratios fall with nact because the NumPy side is itself BLAS-bound by then
!> (see `_ein`); what the engine removes is the temporaries and the dispatch, not
!> the flops.  hdm3 is the outlier because it is not a contraction at all -- see
!> its own header.
!>
!> Still in NumPy, and deliberately: a17, a19, a23, a25, a3, k27 and hdm2, which
!> together cost ~1.2 ms at nact=8.
!>
!> `h2e` is the PHYSICIST-ordered active two-electron tensor (chemist (pq|rs)
!> transposed (0,2,1,3)); `dm3`/`dm4` are the spin-free RDMs in the PySCF
!> make_dm1234 convention.  All arrays cross in the caller's C order (last index
!> fastest) and every offset below is computed explicitly in that order, so the
!> code never relies on an index symmetry to reconcile the two layouts.
!>
!> The two eri-folded 4-pdm intermediates reduce to plain GEMMs once the einsum
!> transposes are folded into the index expression.  Working the two
!> `.transpose(argsort(...))` calls of _f3ca_f3ac through by hand gives
!>
!>     f3ca[o0,o1,o2,o3,o4,o5] = sum_ijk h2e[k,o5,i,j] dm4[o0,o1,o2,o3,o4,j,k,i]
!>     f3ac[o0,o1,o2,o3,o4,o5] = sum_ijk h2e[i,j,k,o4] dm4[o0,o1,o2,o3,j,o5,i,k]
!>
!> In f3ca the contracted indices are exactly the three fastest axes of dm4 and
!> the free indices exactly the five slowest, so the whole thing is ONE DGEMM of
!> (nact x nact^3) against (nact^3 x nact^5) with no repacking at all.  f3ac
!> interleaves o5 among the contracted axes, so it is done as nact^4 small GEMMs
!> over a repacked nact^3 x nact panel.
!>
!> a16 and a22 are both driven by the observation that their dm3 terms all carry
!> the same spectator indices in dm3's SLOWEST axes, so dm3 is a stack of
!> contiguous nact^3 (a16) or nact^4 (a22) blocks; a12 and a13 do the same for Sir, and hdm3/a9/a7
!> close out Sij and Srs.  See each routine's own header for how its terms are laid
!> onto that stack.
!>
module nevpt2_koopmans_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private
  ! bind(C) entry points must stay public: GCC 16.2 hides PRIVATE ones (GCC PR126872)
  public :: nevpt2_a12, nevpt2_a13, nevpt2_a7, nevpt2_a9, &
            nevpt2_hdm3

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: nevpt2_f3ca_f3ac, nevpt2_a16, nevpt2_a22, nevpt2_hdm1
  public :: nevpt2_a3, nevpt2_a17, nevpt2_a19, nevpt2_a23, nevpt2_a25, nevpt2_k27, nevpt2_hdm2
  public :: nevpt2_sijrs, nevpt2_srs, nevpt2_sij

contains

  !> The two eri-folded 4-pdm intermediates f3ca and f3ac.
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  h2e   physicist-ordered active ERIs, C-order [n,n,n,n]
  !> @param[in]  dm4   spin-free active 4-RDM, C-order [n,n,n,n,n,n,n,n]
  !> @param[out] f3ca  C-order [n,n,n,n,n,n]
  !> @param[out] f3ac  C-order [n,n,n,n,n,n]
  subroutine nevpt2_f3ca_f3ac(nact, h2e, dm4, f3ca, f3ac) &
      bind(C, name="nevpt2_f3ca_f3ac")
    integer(c_int32_t), value :: nact
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: h2e(0:*), dm4(0:*)
    real(dp), intent(inout) :: f3ca(0:*), f3ac(0:*)

    ! Local copies in the default (ILP64) integer kind the linked BLAS uses.
    integer :: n, n2, n3, i, j, k, o4, o5
    integer(i8) :: n_i, n3_i, n4_i, n5_i, q, qbase, m
    real(dp), allocatable :: bmat(:,:), amat(:,:), spanel(:,:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n_i = int(n, i8)
    n3_i = n_i**3
    n4_i = n_i**4
    n5_i = n_i**5

    ! ---- f3ca: one DGEMM, F(o5, P) = sum_m B(o5,m) DM(m, P) ----------------
    ! B(o5, m) = h2e[k, o5, i, j] with m = i + n*k + n^2*j, matching the flat
    ! order of dm4's three fastest axes (j, k, i) -> offset j*n^2 + k*n + i.
    allocate(bmat(0:n - 1, 0:n3 - 1))
    do j = 0, n - 1
      do k = 0, n - 1
        do i = 0, n - 1
          m = int(j, i8) * int(n2, i8) + int(k, i8) * n_i + int(i, i8)
          do o5 = 0, n - 1
            bmat(o5, m) = h2e(((int(k, i8) * n_i + int(o5, i8)) * n_i &
                               + int(i, i8)) * n_i + int(j, i8))
          end do
        end do
      end do
    end do
    ! dm4 viewed column-major as (n^3 x n^5): within one P block the three
    ! fastest C axes are contiguous, which is exactly the leading dimension.
    call dgemm('N', 'N', n, int(n5_i), n3, 1.0_dp, bmat, n, dm4, n3, &
               0.0_dp, f3ca, n)
    deallocate(bmat)

    ! ---- f3ac: n^4 small GEMMs --------------------------------------------
    ! A(o4, m) = h2e[i, j, k, o4] with m = k + n*i + n^2*j.
    allocate(amat(0:n - 1, 0:n3 - 1))
    do j = 0, n - 1
      do i = 0, n - 1
        do k = 0, n - 1
          m = int(j, i8) * int(n2, i8) + int(i, i8) * n_i + int(k, i8)
          do o4 = 0, n - 1
            amat(o4, m) = h2e(((int(i, i8) * n_i + int(j, i8)) * n_i &
                               + int(k, i8)) * n_i + int(o4, i8))
          end do
        end do
      end do
    end do

    allocate(spanel(0:n3 - 1, 0:n - 1))
    do q = 0_i8, n4_i - 1_i8
      qbase = q * n4_i
      ! Sp(m, o5) = dm4[Q, j, o5, i, k], m = k + n*i + n^2*j
      do j = 0, n - 1
        do o5 = 0, n - 1
          do i = 0, n - 1
            do k = 0, n - 1
              m = int(j, i8) * int(n2, i8) + int(i, i8) * n_i + int(k, i8)
              spanel(m, o5) = dm4(qbase + int(j, i8) * n3_i &
                                  + int(o5, i8) * int(n2, i8) &
                                  + int(i, i8) * n_i + int(k, i8))
            end do
          end do
        end do
      end do
      ! C(o5, o4) = sum_m Sp(m,o5) A(o4,m); the output block f3ac[Q, o4, o5]
      ! is contiguous with o5 fastest, i.e. exactly column-major (o5, o4).
      call dgemm('T', 'T', n, n, n3, 1.0_dp, spanel, n3, amat, n, &
                 0.0_dp, f3ac(q * int(n2, i8)), n)
    end do

    deallocate(amat, spanel)
  end subroutine nevpt2_f3ca_f3ac

  !> The a16 Koopmans intermediate (the Sr subspace).
  !>
  !> Every dm3 term carries the same spectator triple (r,p,q) in dm3's three
  !> SLOWEST axes, so dm3 is a single (nact^3 x nact^3) matrix whose columns are
  !> the spectators and whose rows are the contiguous block that a16[p,q,r,:,:,:]
  !> is built from.  Naming the three block axes h1..h3 fastest-first, the terms
  !> contract
  !>
  !>     h3        terms 1+11         h2,h1     terms 5 and 12
  !>     h2        term 2             h3,h2     term 6
  !>     h1        terms 3+13         h1,h3     term 7
  !>     h1,h2,h3  the fdm2 fold
  !>
  !> so three permuted copies of dm3 -- P6=(h3,h2,h1), P2=(h2,h3,h1) and
  !> P7=(h1,h3,h2) -- together with dm3 itself put the contracted axes leading in
  !> every case, and the whole intermediate is eight GEMMs over the full stack.
  !> The four two-index terms are the n^8 ones; the first revision of this
  !> routine evaluated them as scalar loop nests over one nact^3 block at a time,
  !> which is what made it lose to NumPy's einsum at nact=8.
  !>
  !> Two pairs of terms also collapse before any arithmetic happens --
  !>
  !>     -sum_i h1e[i,b] G[i,a,c] + sum_ij h2e[j,b,i,j] G[i,a,c]
  !>   = sum_i ( sum_j h2e[j,b,i,j] - h1e[i,b] ) G[i,a,c]
  !>
  !> and likewise for the pair carrying `c` -- so the thirteen einsum calls of
  !> the Python become eight GEMMs plus the f3 fold-ins.
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  h1e   active one-electron matrix, C-order [n,n]
  !> @param[in]  h2e   physicist-ordered active ERIs, C-order [n,n,n,n]
  !> @param[in]  dm3   spin-free active 3-RDM, C-order [n,n,n,n,n,n]
  !> @param[in]  f3ca  C-order [n,n,n,n,n,n]
  !> @param[in]  f3ac  C-order [n,n,n,n,n,n]
  !> @param[out] a16   C-order [n,n,n,n,n,n]
  subroutine nevpt2_a16(nact, h1e, h2e, dm3, f3ca, f3ac, a16) &
      bind(C, name="nevpt2_a16")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm3(0:*)
    real(dp), intent(in) :: f3ca(0:*), f3ac(0:*)
    real(dp), intent(inout) :: a16(0:*)

    ! OpenQP is ILP64, so the default integer is 8 bytes and every offset here
    ! (at most n^6) is far inside its range; the BLAS dimension arguments below
    ! therefore take the default kind, as the linked ILP64 BLAS expects.
    integer :: n, n2, n3, n4, n5, n6
    integer :: p, q, r, a, b, c, i, j, u, v, y, z, s, sb, ob, rp
    integer :: h1x, h2x, h3x
    real(dp) :: acc, val
    real(dp), allocatable :: u1(:,:), u2(:,:), h1t(:,:)
    real(dp), allocatable :: b5(:,:), b6(:,:), b7(:,:), b12(:,:), bf(:,:)
    real(dp), allocatable :: p6(:), p2(:), p7(:)
    real(dp), allocatable :: ra(:), rb(:), rc(:), rd(:), fd(:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n4 = n3 * n
    n5 = n4 * n
    n6 = n5 * n

    ! u1(i,b) = sum_j h2e[j,b,i,j] - h1e[i,b]   (terms 1 + 11 folded)
    ! u2(i,c) = sum_j h2e[j,c,i,j] - h1e[c,i]   (terms 3 + 13 folded)
    allocate(u1(0:n - 1, 0:n - 1), u2(0:n - 1, 0:n - 1), h1t(0:n - 1, 0:n - 1))
    do b = 0, n - 1
      do i = 0, n - 1
        acc = 0.0_dp
        do j = 0, n - 1
          acc = acc + h2e(((j * n + b) * n + i) * n + j)
        end do
        u1(i, b) = acc - h1e(i * n + b)
        u2(i, b) = acc - h1e(b * n + i)
        h1t(i, b) = h1e(i * n + b)      ! h1e read the other way round
      end do
    end do

    ! The four two-index ERI panels.  Row index = the contracted pair in the
    ! order the matching dm3 copy stores it; column index = the pair of new
    ! output indices in the order the result buffer wants them.
    allocate(b5(0:n2 - 1, 0:n2 - 1), b6(0:n2 - 1, 0:n2 - 1), &
             b7(0:n2 - 1, 0:n2 - 1), b12(0:n2 - 1, 0:n2 - 1))
    do v = 0, n - 1
      do u = 0, n - 1
        do y = 0, n - 1
          do z = 0, n - 1
            ! b5 (i+nk, b+na) = h2e[k,b,i,a]
            b5(z + n * y, u + n * v) = h2e(((y * n + u) * n + z) * n + v)
            ! b6 (j+nk, b+na) = h2e[k,b,a,j]
            b6(z + n * y, u + n * v) = h2e(((y * n + u) * n + v) * n + z)
            ! b7 (i+nj, c+nb) = h2e[c,b,i,j]
            b7(z + n * y, u + n * v) = h2e(((u * n + v) * n + z) * n + y)
            ! b12(k+nj, c+na) = h2e[c,j,k,a]
            b12(z + n * y, u + n * v) = h2e(((u * n + y) * n + z) * n + v)
          end do
        end do
      end do
    end do

    ! bf(i+nk+n2 j, b) = h2e[k,b,i,j] -- the fdm2 fold contracts all three
    ! block axes, so its spectator triple is (r,p,a), not (r,p,q).
    allocate(bf(0:n3 - 1, 0:n - 1))
    do i = 0, n - 1
      do y = 0, n - 1
        do j = 0, n - 1
          do b = 0, n - 1
            bf(i + n * y + n2 * j, b) = h2e(((y * n + b) * n + i) * n + j)
          end do
        end do
      end do
    end do

    ! ---- the three permuted copies of dm3 ----------------------------------
    allocate(p6(0:n6 - 1), p2(0:n6 - 1), p7(0:n6 - 1))
    do s = 0, n3 - 1
      sb = s * n3
      do h3x = 0, n - 1
        do h2x = 0, n - 1
          do h1x = 0, n - 1
            val = dm3(sb + h1x + n * h2x + n2 * h3x)
            p6(sb + h3x + n * h2x + n2 * h1x) = val
            p2(sb + h2x + n * h3x + n2 * h1x) = val
            p7(sb + h1x + n * h3x + n2 * h2x) = val
          end do
        end do
      end do
    end do

    ! ---- one GEMM per term, over the whole spectator stack -----------------
    ! Terms that land on the same flat offset accumulate into one buffer.
    allocate(ra(0:n6 - 1), rb(0:n6 - 1), rc(0:n6 - 1), rd(0:n6 - 1), &
             fd(0:n4 - 1))

    ! ra: +u1[b,i] G[i,a,c] (1+11), -h2e[k,b,a,j] G[j,k,c] (6),
    !     -h2e[k,b,i,a] G[c,k,i] (5)
    call dgemm('T', 'N', n, n5, n, 1.0_dp, u1, n, p6, n, 0.0_dp, ra, n)
    call dgemm('T', 'N', n2, n4, n2, -1.0_dp, b6, n2, p6, n2, 1.0_dp, ra, n2)
    call dgemm('T', 'N', n2, n4, n2, -1.0_dp, b5, n2, dm3, n2, 1.0_dp, ra, n2)
    ! rb: +h1e[i,a] G[b,i,c] (2)
    call dgemm('T', 'N', n, n5, n, 1.0_dp, h1t, n, p2, n, 0.0_dp, rb, n)
    ! rc: +u2[c,i] G[b,a,i] (3+13), -h2e[c,j,k,a] G[b,j,k] (12)
    call dgemm('T', 'N', n, n5, n, 1.0_dp, u2, n, dm3, n, 0.0_dp, rc, n)
    call dgemm('T', 'N', n2, n4, n2, -1.0_dp, b12, n2, dm3, n2, 1.0_dp, rc, n2)
    ! rd: +h2e[c,b,i,j] G[j,a,i] (7)
    call dgemm('T', 'N', n2, n4, n2, 1.0_dp, b7, n2, p7, n2, 0.0_dp, rd, n2)
    ! fd(b, (r,p,a)) = sum_kij h2e[k,b,i,j] dm3[r,p,a,j,k,i]
    call dgemm('T', 'N', n, n3, n3, 1.0_dp, bf, n3, dm3, n3, 0.0_dp, fd, n)

    deallocate(p6, p2, p7)

    ! ---- gather, one spectator triple at a time ----------------------------
    do s = 0, n3 - 1
      r = s / n2
      rp = s - r * n2
      p = rp / n
      q = rp - p * n
      sb = s * n3
      ob = ((p * n + q) * n + r) * n3
      do a = 0, n - 1
        do b = 0, n - 1
          do c = 0, n - 1
            ! terms 9 and 10: + f3ac[r,p,q,b,a,c] - f3ca[r,p,q,b,a,c]
            ! term 4:         - f3ca[r,p,a,c,q,b]
            val = ra(sb + b + n * a + n2 * c) &
                + rb(sb + a + n * b + n2 * c) &
                + rc(sb + c + n * a + n2 * b) &
                + rd(sb + c + n * b + n2 * a) &
                + f3ac(sb + b * n2 + a * n + c) &
                - f3ca(sb + b * n2 + a * n + c) &
                - f3ca(((((r * n + p) * n + a) * n + c) * n + q) * n + b)
            a16(ob + a * n2 + b * n + c) = val
          end do
        end do
      end do
      ! term 8: a16[p,q,r,a,b,c] += fdm2[p,r,a,b] on the q == c plane
      do a = 0, n - 1
        do b = 0, n - 1
          a16(ob + a * n2 + b * n + q) = a16(ob + a * n2 + b * n + q) &
                                       + fd(b + n * (((r * n + p) * n + a)))
        end do
      end do
    end do

    deallocate(u1, u2, h1t, b5, b6, b7, b12, bf, ra, rb, rc, rd, fd)
  end subroutine nevpt2_a16

  !> The a22 Koopmans intermediate (the Si subspace).
  !>
  !> Every one of the ~20 einsum terms of the Python contracts a dm3 whose two
  !> SLOWEST axes are the same spectator pair (k,i) -- `kipjac`, `kibjpc`,
  !> `kiqcpr`, ... -- so dm3 is one stack of n^2 contiguous n^4 blocks and each
  !> term is a single large GEMM over the whole stack, provided the contracted
  !> axes sit fastest inside a block.  Naming dm3's four block axes g1..g4 from
  !> fastest to slowest (g1 = dm3's last C axis), the terms contract
  !>
  !>     g1        term 3+5          g1,g2     terms 4, 9
  !>     g2        term 2            g2,g4     term 8
  !>     g4        term 1            g1,g4     term 10
  !>     g1,g2,g3  X2                g3,g4     term 14
  !>     g1,g4,g2  X1 (= both fdm2 folds)
  !>
  !> which four permuted copies of dm3 cover: TB=(g2,g4,g1,g3) serves terms 2
  !> and 8, TC=(g1,g4,g2,g3) serves term 10 and X1, TD=(g3,g4,g1,g2) serves
  !> term 14, TE=(g4,g1,g2,g3) serves term 1, and dm3 itself serves terms 3+5,
  !> 4, 9 and X2.  The repacks are n^6 each; the GEMMs they enable are n^8.
  !>
  !> Two pairs of terms fold before any arithmetic.  With
  !>
  !>     w[c,p] = h1e[c,p] - sum_q h2e[q,c,p,q]
  !>
  !>       + sum_p h1e[c,p] G[b,j,a,p] - sum_pq h2e[q,c,p,q] G[b,j,a,p]
  !>     = + sum_p w[c,p] G[b,j,a,p]
  !>
  !> and identically inside the second fdm2 fold, so twenty einsum calls become
  !> ten GEMMs.  The two `for i in range(...)` diagonal inserts of the Python are
  !> the j == a and j == b planes of the output and are applied after the gather.
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  h1e   active one-electron matrix, C-order [n,n]
  !> @param[in]  h2e   physicist-ordered active ERIs, C-order [n,n,n,n]
  !> @param[in]  dm2   spin-free active 2-RDM, C-order [n,n,n,n]
  !> @param[in]  dm3   spin-free active 3-RDM, C-order [n,n,n,n,n,n]
  !> @param[in]  f3ca  C-order [n,n,n,n,n,n]
  !> @param[in]  f3ac  C-order [n,n,n,n,n,n]
  !> @param[out] a22   C-order [n,n,n,n,n,n]
  subroutine nevpt2_a22(nact, h1e, h2e, dm2, dm3, f3ca, f3ac, a22) &
      bind(C, name="nevpt2_a22")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm2(0:*), dm3(0:*)
    real(dp), intent(in) :: f3ca(0:*), f3ac(0:*)
    real(dp), intent(inout) :: a22(0:*)

    ! OpenQP is ILP64, so the default integer is 8 bytes and every offset here
    ! (at most n^6) is far inside its range; the BLAS dimension arguments below
    ! therefore take the default kind, as the linked ILP64 BLAS expects.
    integer :: n, n2, n3, n4, n5, n6
    integer :: p, q, r, u, v, y, z, i, j, k, a, b, c, aa, bb, cc, ki, m
    integer :: g1, g2, g3, g4
    integer :: kb, ob, d2b, xb
    real(dp) :: acc, val, hjb
    real(dp), allocatable :: h1t(:,:), w(:,:), yblk(:,:)
    real(dp), allocatable :: b4(:,:), b8(:,:), b9(:,:), b10(:,:), b14(:,:)
    real(dp), allocatable :: bx1(:,:), bx2(:,:)
    real(dp), allocatable :: tb(:), tc(:), td(:), te(:)
    real(dp), allocatable :: r1(:), r2(:), r34(:), r8(:), r9(:), r10(:), r14(:)
    real(dp), allocatable :: x1(:), x2(:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n4 = n3 * n
    n5 = n4 * n
    n6 = n5 * n

    ! ---- operands that do not depend on the RDMs ---------------------------
    allocate(h1t(0:n - 1, 0:n - 1), w(0:n - 1, 0:n - 1))
    do b = 0, n - 1
      do p = 0, n - 1
        h1t(p, b) = h1e(p * n + b)            ! h1e read the other way round
      end do
    end do
    do c = 0, n - 1
      do p = 0, n - 1
        acc = 0.0_dp
        do q = 0, n - 1
          acc = acc + h2e(((q * n + c) * n + p) * n + q)
        end do
        w(p, c) = h1e(c * n + p) - acc
      end do
    end do

    ! The five two-index ERI panels.  Their row index is the contracted pair in
    ! the order the matching dm3 copy stores it, their column index is the pair
    ! of new output indices in the order the output block wants them.
    allocate(b4(0:n2 - 1, 0:n2 - 1), b8(0:n2 - 1, 0:n2 - 1), &
             b9(0:n2 - 1, 0:n2 - 1), b10(0:n2 - 1, 0:n2 - 1), &
             b14(0:n2 - 1, 0:n2 - 1))
    do v = 0, n - 1
      do u = 0, n - 1
        do y = 0, n - 1
          do z = 0, n - 1
            ! b4 (r+nq, c+na) = h2e[c,q,r,a]
            b4(z + n * y, u + n * v) = h2e(((u * n + y) * n + z) * n + v)
            ! b8 (p+nq, b+na) = h2e[p,q,a,b]
            b8(y + n * z, u + n * v) = h2e(((y * n + z) * n + v) * n + u)
            ! b9 (r+np, c+nb) = h2e[p,c,r,b]
            b9(z + n * y, u + n * v) = h2e(((y * n + u) * n + z) * n + v)
            ! b10(r+nq, c+nb) = h2e[c,q,r,b]
            b10(z + n * y, u + n * v) = h2e(((u * n + y) * n + z) * n + v)
            ! b14(r+np, b+nj) = h2e[p,j,r,b]
            b14(z + n * y, u + n * v) = h2e(((y * n + v) * n + z) * n + u)
          end do
        end do
      end do
    end do

    allocate(bx1(0:n3 - 1, 0:n - 1), bx2(0:n3 - 1, 0:n - 1))
    do p = 0, n - 1
      do q = 0, n - 1
        do r = 0, n - 1
          do u = 0, n - 1
            ! bx1(r+nq+n2 p, x) = h2e[p,q,r,x]
            bx1(r + n * q + n2 * p, u) = h2e(((p * n + q) * n + r) * n + u)
            ! bx2(p+nr+n2 q, c) = h2e[r,c,p,q]
            bx2(p + n * r + n2 * q, u) = h2e(((r * n + u) * n + p) * n + q)
          end do
        end do
      end do
    end do

    ! ---- the four permuted copies of dm3 -----------------------------------
    allocate(tb(0:n6 - 1), tc(0:n6 - 1), td(0:n6 - 1), te(0:n6 - 1))
    do ki = 0, n2 - 1
      kb = ki * n4
      do g4 = 0, n - 1
        do g3 = 0, n - 1
          do g2 = 0, n - 1
            do g1 = 0, n - 1
              val = dm3(kb + g1 + n * g2 + n2 * g3 + n3 * g4)
              tb(kb + g2 + n * g4 + n2 * g1 + n3 * g3) = val
              tc(kb + g1 + n * g4 + n2 * g2 + n3 * g3) = val
              td(kb + g3 + n * g4 + n2 * g1 + n3 * g2) = val
              te(kb + g4 + n * g1 + n2 * g2 + n3 * g3) = val
            end do
          end do
        end do
      end do
    end do

    ! ---- one GEMM per term, over the whole (k,i) stack ----------------------
    ! Every call is C = alpha * A^T B: A holds the ERI/one-electron panel with
    ! the contracted index leading, B is the dm3 copy that stores that same
    ! index leading, and the result buffer is read back below with the flat
    ! offset the gather needs.
    allocate(r1(0:n6 - 1), r2(0:n6 - 1), r34(0:n6 - 1), r8(0:n6 - 1), &
             r9(0:n6 - 1), r10(0:n6 - 1), r14(0:n6 - 1), &
             x1(0:n4 - 1), x2(0:n4 - 1))

    ! term 1   -h1e[p,b] dm3[k,i,p,j,a,c]
    call dgemm('T', 'N', n, n5, n, -1.0_dp, h1t, n, te, n, 0.0_dp, r1, n)
    ! term 2   -h1e[p,a] dm3[k,i,b,j,p,c]
    call dgemm('T', 'N', n, n5, n, -1.0_dp, h1t, n, tb, n, 0.0_dp, r2, n)
    ! terms 3+5 +w[c,p] dm3[k,i,b,j,a,p]
    call dgemm('T', 'N', n, n5, n, 1.0_dp, w, n, dm3, n, 0.0_dp, r34, n)
    ! term 4   +h2e[c,q,r,a] dm3[k,i,b,j,q,r]  -- same flat offset as 3+5, so
    ! it accumulates straight into r34 instead of costing another n^6 buffer.
    call dgemm('T', 'N', n2, n4, n2, 1.0_dp, b4, n2, dm3, n2, 1.0_dp, r34, n2)
    ! term 8   -h2e[p,q,a,b] dm3[k,i,q,j,p,c]
    call dgemm('T', 'N', n2, n4, n2, -1.0_dp, b8, n2, tb, n2, 0.0_dp, r8, n2)
    ! term 9   +h2e[p,c,r,b] dm3[k,i,a,j,p,r]
    call dgemm('T', 'N', n2, n4, n2, 1.0_dp, b9, n2, dm3, n2, 0.0_dp, r9, n2)
    ! term 10  +h2e[c,q,r,b] dm3[k,i,q,j,a,r]
    call dgemm('T', 'N', n2, n4, n2, 1.0_dp, b10, n2, tc, n2, 0.0_dp, r10, n2)
    ! term 14  +2 h2e[p,j,r,b] dm3[k,i,p,r,a,c]
    call dgemm('T', 'N', n2, n4, n2, 2.0_dp, b14, n2, td, n2, 0.0_dp, r14, n2)
    ! X1(x,c) = h2e[p,q,r,x] dm3[k,i,q,c,p,r] -- the Python builds this same
    ! object twice, once as the j==a insert and once inside the j==b insert.
    call dgemm('T', 'N', n, n3, n3, 1.0_dp, bx1, n3, tc, n3, 0.0_dp, x1, n)
    ! X2(c,a) = h2e[r,c,p,q] dm3[k,i,a,q,r,p]
    call dgemm('T', 'N', n, n3, n3, 1.0_dp, bx2, n3, dm3, n3, 0.0_dp, x2, n)

    deallocate(tb, tc, td, te)

    ! ---- gather, one (k,i) block at a time ---------------------------------
    allocate(yblk(0:n - 1, 0:n - 1))
    do ki = 0, n2 - 1
      k = ki / n
      i = ki - k * n
      kb = ki * n4
      d2b = ki * n2
      xb = ki * n2

      ! Y(a,c): the whole second fdm2 fold, inserted on the j == b plane.
      do c = 0, n - 1
        do a = 0, n - 1
          acc = 0.0_dp
          do p = 0, n - 1
            acc = acc + h1t(p, a) * dm2(d2b + p * n + c) &
                      - w(p, c) * dm2(d2b + a * n + p)
          end do
          ! -h2e[c,q,r,a] dm2[k,i,q,r]: the 2-RDM block is flat in exactly the
          ! (r + n q) order b4 already stores.
          do m = 0, n2 - 1
            acc = acc - b4(m, c + n * a) * dm2(d2b + m)
          end do
          yblk(a, c) = acc + x1(xb + a + n * c) - x2(xb + c + n * a)
        end do
      end do

      do j = 0, n - 1
        ob = ((i * n + j) * n + k) * n3
        do a = 0, n - 1
          do b = 0, n - 1
            hjb = 2.0_dp * h1e(j * n + b)    ! term 13, +2 h1e[j,b] dm2[k,i,a,c]
            do c = 0, n - 1
              val = r1(kb + b + n * c + n2 * a + n3 * j) &
                  + r2(kb + a + n * b + n2 * c + n3 * j) &
                  + r34(kb + c + n * a + n2 * j + n3 * b) &
                  + r8(kb + b + n * a + n2 * c + n3 * j) &
                  + r9(kb + c + n * b + n2 * j + n3 * a) &
                  + r10(kb + c + n * b + n2 * a + n3 * j) &
                  + r14(kb + b + n * j + n2 * c + n3 * a) &
                  - f3ac(kb + j + n * b + n2 * c + n3 * a) &
                  - f3ac(kb + c + n * a + n2 * j + n3 * b) &
                  + f3ca(kb + c + n * a + n2 * j + n3 * b) &
                  + hjb * dm2(d2b + a * n + c)
              a22(ob + a * n2 + b * n + c) = val
            end do
          end do
        end do
        ! the two diagonal inserts: a22[:,x,:,x,:,:] -= X1 and
        ! a22[:,x,:,:,x,:] += 2 Y, i.e. the j == a and j == b planes
        do bb = 0, n - 1
          do cc = 0, n - 1
            a22(ob + j * n2 + bb * n + cc) = a22(ob + j * n2 + bb * n + cc) &
                                           - x1(xb + bb + n * cc)
          end do
        end do
        do aa = 0, n - 1
          do cc = 0, n - 1
            a22(ob + aa * n2 + j * n + cc) = a22(ob + aa * n2 + j * n + cc) &
                                           + 2.0_dp * yblk(aa, cc)
          end do
        end do
      end do
    end do

    deallocate(h1t, w, yblk, b4, b8, b9, b10, b14, bx1, bx2)
    deallocate(r1, r2, r34, r8, r9, r10, r14, x1, x2)
  end subroutine nevpt2_a22

  !> The hdm3 "hole" 3-RDM of the Sij subspace.
  !>
  !> Despite the nine einsum calls this is not a contraction at all: every term
  !> carries at least one Kronecker delta, so each is a plane of the output
  !> rather than a sum over anything.  NumPy has to materialise the full n^6
  !> outer product for each of them; here the deltas are loop bounds, and the
  !> whole intermediate is one n^6 transpose of dm3 plus six n^5 planes and two
  !> n^4 lines.  Written out for the output element [p,q,r,a,b,c] the nine terms
  !> are
  !>
  !>     -dm3[b,q,a,p,c,r]                                       (always)
  !>     -hdm2[q,r,a,c]   on p == b        -hdm2[p,q,a,c]  on b == r
  !>     +2 hdm2[p,r,a,c] on b == q        +2 dm2[b,q,c,r] on a == p
  !>     +2 dm2[b,q,a,p]  on c == r        -dm2[b,q,c,p]   on a == r
  !>     -4 dm1[b,q]      on a == p and c == r
  !>     +2 dm1[b,q]      on a == r and p == c
  !>
  !> `hdm1` is an argument of the Python but never appears in its body, so it is
  !> not passed here.
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  dm1   spin-free active 1-RDM, C-order [n,n]
  !> @param[in]  dm2   spin-free active 2-RDM, C-order [n,n,n,n]
  !> @param[in]  dm3   spin-free active 3-RDM, C-order [n,n,n,n,n,n]
  !> @param[in]  hdm2  the hole 2-RDM, C-order [n,n,n,n]
  !> @param[out] hdm3  C-order [n,n,n,n,n,n]
  subroutine nevpt2_hdm3(nact, dm1, dm2, dm3, hdm2, hdm3) &
      bind(C, name="nevpt2_hdm3")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: dm1(0:*), dm2(0:*), dm3(0:*), hdm2(0:*)
    real(dp), intent(inout) :: hdm3(0:*)

    integer :: n, n2, n3, p, q, r, a, b, c, ob, qr

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n

    ! The one term with no delta.  `c` runs innermost so the STORE is
    ! sequential; the dm3 read is then stride n, which is the better trade.
    do p = 0, n - 1
      do q = 0, n - 1
        do r = 0, n - 1
          ob = ((p * n + q) * n + r) * n3
          do a = 0, n - 1
            do b = 0, n - 1
              qr = (((b * n + q) * n + a) * n + p) * n2 + r
              do c = 0, n - 1
                hdm3(ob + a * n2 + b * n + c) = -dm3(qr + c * n)
              end do
            end do
          end do
        end do
      end do
    end do

    do p = 0, n - 1
      do q = 0, n - 1
        do r = 0, n - 1
          ob = ((p * n + q) * n + r) * n3
          do a = 0, n - 1
            do c = 0, n - 1
              hdm3(ob + a * n2 + p * n + c) = hdm3(ob + a * n2 + p * n + c) &
                - hdm2(((q * n + r) * n + a) * n + c)
              hdm3(ob + a * n2 + r * n + c) = hdm3(ob + a * n2 + r * n + c) &
                - hdm2(((p * n + q) * n + a) * n + c)
              hdm3(ob + a * n2 + q * n + c) = hdm3(ob + a * n2 + q * n + c) &
                + 2.0_dp * hdm2(((p * n + r) * n + a) * n + c)
            end do
          end do
          do b = 0, n - 1
            do c = 0, n - 1
              hdm3(ob + p * n2 + b * n + c) = hdm3(ob + p * n2 + b * n + c) &
                + 2.0_dp * dm2(((b * n + q) * n + c) * n + r)
              hdm3(ob + r * n2 + b * n + c) = hdm3(ob + r * n2 + b * n + c) &
                - dm2(((b * n + q) * n + c) * n + p)
            end do
          end do
          do a = 0, n - 1
            do b = 0, n - 1
              hdm3(ob + a * n2 + b * n + r) = hdm3(ob + a * n2 + b * n + r) &
                + 2.0_dp * dm2(((b * n + q) * n + a) * n + p)
            end do
          end do
          do b = 0, n - 1
            hdm3(ob + p * n2 + b * n + r) = hdm3(ob + p * n2 + b * n + r) &
              - 4.0_dp * dm1(b * n + q)
            hdm3(ob + r * n2 + b * n + p) = hdm3(ob + r * n2 + b * n + p) &
              + 2.0_dp * dm1(b * n + q)
          end do
        end do
      end do
    end do
  end subroutine nevpt2_hdm3

  !> The a9 Koopmans intermediate (the Sij subspace).
  !>
  !> Six of the nine terms are one-index contractions against hdm2, and they
  !> collapse into two applications of a single folded panel: with
  !>
  !>     v[t,x] = h1e[t,x] + 2 sum_i h2e[i,t,i,x] - sum_j h2e[t,j,j,x]
  !>
  !> terms 1+2+3 are sum_t v[t,b] hdm2[p,q,a,t] and terms 5+6+8 are the same
  !> panel again as sum_t v[t,a] hdm2[p,q,t,b] -- the two folds are literally
  !> the same function of the integrals, which is why v is built once.  What is
  !> left is one n^6 term against hdm2 and the two n^7 terms against hdm3.
  !>
  !> The hdm3 terms contract axes that are neither leading nor contiguous
  !> (1,4,5 and 2,3,5), so each gets a permuted copy that brings its three
  !> contracted axes to the front; those two n^6 repacks are what turn the pair
  !> into single GEMMs over the whole array.  `hdm1` is an argument of the
  !> Python but never appears in its body, so it is not passed here.
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  h1e   active one-electron matrix, C-order [n,n]
  !> @param[in]  h2e   physicist-ordered active ERIs, C-order [n,n,n,n]
  !> @param[in]  hdm2  the hole 2-RDM, C-order [n,n,n,n]
  !> @param[in]  hdm3  the hole 3-RDM, C-order [n,n,n,n,n,n]
  !> @param[out] a9    C-order [n,n,n,n]
  subroutine nevpt2_a9(nact, h1e, h2e, hdm2, hdm3, a9) bind(C, name="nevpt2_a9")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), hdm2(0:*), hdm3(0:*)
    real(dp), intent(inout) :: a9(0:*)

    integer :: n, n2, n3, n4, n6
    integer :: p, q, a, b, i, j, k, t, x
    real(dp) :: s2, s3
    real(dp), allocatable :: v1(:,:), bc(:,:), bd(:,:), be(:,:)
    real(dp), allocatable :: pb(:), pd(:), pe(:), ra(:), rb(:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n4 = n3 * n
    n6 = n4 * n2

    allocate(v1(0:n - 1, 0:n - 1))
    do t = 0, n - 1
      do x = 0, n - 1
        s2 = 0.0_dp
        s3 = 0.0_dp
        do i = 0, n - 1
          s2 = s2 + h2e(((i * n + t) * n + i) * n + x)
          s3 = s3 + h2e(((t * n + i) * n + i) * n + x)
        end do
        v1(t, x) = h1e(t * n + x) + 2.0_dp * s2 - s3
      end do
    end do

    ! bc(i+nj, b+na) = h2e[i,j,b,a]
    allocate(bc(0:n2 - 1, 0:n2 - 1))
    do i = 0, n - 1
      do j = 0, n - 1
        do b = 0, n - 1
          do a = 0, n - 1
            bc(i + n * j, b + n * a) = h2e(((i * n + j) * n + b) * n + a)
          end do
        end do
      end do
    end do

    ! bd(j+ni+n2 k, x) and be(i+nj+n2 k, x) are the same ERI panel read with the
    ! two different contracted-index orders the two hdm3 copies store.
    allocate(bd(0:n3 - 1, 0:n - 1), be(0:n3 - 1, 0:n - 1))
    do i = 0, n - 1
      do j = 0, n - 1
        do k = 0, n - 1
          do x = 0, n - 1
            bd(j + n * i + n2 * k, x) = h2e(((i * n + j) * n + k) * n + x)
            be(i + n * j + n2 * k, x) = h2e(((i * n + j) * n + k) * n + x)
          end do
        end do
      end do
    end do

    ! pb(t,b,q,p) from hdm2[p,q,t,b]: term B contracts hdm2's THIRD axis, so
    ! the two fastest have to be swapped before it can lead a GEMM.
    allocate(pb(0:n4 - 1))
    do p = 0, n - 1
      do q = 0, n - 1
        do t = 0, n - 1
          do b = 0, n - 1
            pb(t + n * b + n2 * q + n3 * p) = hdm2(((p * n + q) * n + t) * n + b)
          end do
        end do
      end do
    end do

    allocate(pd(0:n6 - 1), pe(0:n6 - 1))
    ! pd(j,i,k, a,q,p) from hdm3[p,k,q,a,i,j]
    do p = 0, n - 1
      do k = 0, n - 1
        do q = 0, n - 1
          do a = 0, n - 1
            do i = 0, n - 1
              do j = 0, n - 1
                pd(j + n * i + n2 * k + n3 * (a + n * q + n2 * p)) = &
                  hdm3(((((p * n + k) * n + q) * n + a) * n + i) * n + j)
              end do
            end do
          end do
        end do
      end do
    end do
    ! pe(i,j,k, b,q,p) from hdm3[p,q,k,j,b,i]
    do p = 0, n - 1
      do q = 0, n - 1
        do k = 0, n - 1
          do j = 0, n - 1
            do b = 0, n - 1
              do i = 0, n - 1
                pe(i + n * j + n2 * k + n3 * (b + n * q + n2 * p)) = &
                  hdm3(((((p * n + q) * n + k) * n + j) * n + b) * n + i)
              end do
            end do
          end do
        end do
      end do
    end do

    allocate(ra(0:n4 - 1), rb(0:n4 - 1))
    ! ra: terms 1+2+3, then term 7, then term 4 -- all on flat b+na+n2 q+n3 p
    call dgemm('T', 'N', n, n3, n, 1.0_dp, v1, n, hdm2, n, 0.0_dp, ra, n)
    call dgemm('T', 'N', n2, n2, n2, -1.0_dp, bc, n2, hdm2, n2, 1.0_dp, ra, n2)
    call dgemm('T', 'N', n, n3, n3, -1.0_dp, bd, n3, pd, n3, 1.0_dp, ra, n)
    ! rb: terms 5+6+8, then term 9 -- both on flat a+nb+n2 q+n3 p
    call dgemm('T', 'N', n, n3, n, 1.0_dp, v1, n, pb, n, 0.0_dp, rb, n)
    call dgemm('T', 'N', n, n3, n3, -1.0_dp, be, n3, pe, n3, 1.0_dp, rb, n)

    do p = 0, n - 1
      do q = 0, n - 1
        do a = 0, n - 1
          do b = 0, n - 1
            a9(((p * n + q) * n + a) * n + b) = ra(b + n * a + n2 * q + n3 * p) &
                                              + rb(a + n * b + n2 * q + n3 * p)
          end do
        end do
      end do
    end do

    deallocate(v1, bc, bd, be, pb, pd, pe, ra, rb)
  end subroutine nevpt2_a9

  !> The rm2/rm3 "reduced" RDMs and the a7 intermediate (the Srs subspace).
  !>
  !> rm2 is returned as well as consumed: the caller uses it for the Srs norm,
  !> exactly as the Python `_a7` returns the pair.
  !>
  !> rm2 and rm3 are delta-dressed transposes of dm2/dm3 -- no contraction --
  !> so they are built by direct scatter.  a7 itself is five terms, of which the
  !> two n^7 ones contract rm3's axes (2,3,4) and (2,3,5); as in a9 those sit in
  !> the middle of the array, so each takes a permuted copy that brings its
  !> three contracted axes to the front.  Both then share the SAME ERI panel,
  !> h2e[k,x,i,j], because the two terms differ only in which free index plays
  !> the role of `x`.
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  h1e   active one-electron matrix, C-order [n,n]
  !> @param[in]  h2e   physicist-ordered active ERIs, C-order [n,n,n,n]
  !> @param[in]  dm1   spin-free active 1-RDM, C-order [n,n]
  !> @param[in]  dm2   spin-free active 2-RDM, C-order [n,n,n,n]
  !> @param[in]  dm3   spin-free active 3-RDM, C-order [n,n,n,n,n,n]
  !> @param[out] rm2   C-order [n,n,n,n]
  !> @param[out] a7    C-order [n,n,n,n]
  subroutine nevpt2_a7(nact, h1e, h2e, dm1, dm2, dm3, rm2, a7) &
      bind(C, name="nevpt2_a7")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm1(0:*), dm2(0:*), dm3(0:*)
    real(dp), intent(inout) :: rm2(0:*), a7(0:*)

    integer :: n, n2, n3, n4, n6
    integer :: p, q, a, b, i, j, k, l, m, u, x
    real(dp) :: val
    real(dp), allocatable :: ha(:,:), b5(:,:), b3(:,:)
    real(dp), allocatable :: rm3(:), q1(:), p3(:), p4(:), ra(:), rb(:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n4 = n3 * n
    n6 = n4 * n2

    ! rm2[i,j,k,l] = dm2[i,l,j,k] - delta[j,l] dm1[i,k]
    do i = 0, n - 1
      do j = 0, n - 1
        do k = 0, n - 1
          do l = 0, n - 1
            val = dm2(((i * n + l) * n + j) * n + k)
            if (j == l) val = val - dm1(i * n + k)
            rm2(((i * n + j) * n + k) * n + l) = val
          end do
        end do
      end do
    end do

    ! rm3[i,j,k,l,m,u] = dm3[i,u,j,m,k,l] - delta[j,u] dm2[i,m,k,l]
    !                    - delta[k,m] rm2[i,j,l,u] - delta[k,u] rm2[i,j,m,l]
    allocate(rm3(0:n6 - 1))
    do i = 0, n - 1
      do j = 0, n - 1
        do k = 0, n - 1
          do l = 0, n - 1
            do m = 0, n - 1
              do u = 0, n - 1
                val = dm3(((((i * n + u) * n + j) * n + m) * n + k) * n + l)
                if (j == u) val = val - dm2(((i * n + m) * n + k) * n + l)
                if (k == m) val = val - rm2(((i * n + j) * n + l) * n + u)
                if (k == u) val = val - rm2(((i * n + j) * n + m) * n + l)
                rm3(((((i * n + j) * n + k) * n + l) * n + m) * n + u) = val
              end do
            end do
          end do
        end do
      end do
    end do

    allocate(ha(0:n - 1, 0:n - 1))
    do i = 0, n - 1
      do x = 0, n - 1
        ha(i, x) = h1e(x * n + i)          ! ha(i,x) = h1e[x,i]
      end do
    end do

    ! b5(j+ni, b+na) = h2e[b,a,i,j]
    allocate(b5(0:n2 - 1, 0:n2 - 1))
    do i = 0, n - 1
      do j = 0, n - 1
        do b = 0, n - 1
          do a = 0, n - 1
            b5(j + n * i, b + n * a) = h2e(((b * n + a) * n + i) * n + j)
          end do
        end do
      end do
    end do

    ! b3(j+ni+n2 k, x) = h2e[k,x,i,j] -- shared by both n^7 terms
    allocate(b3(0:n3 - 1, 0:n - 1))
    do i = 0, n - 1
      do j = 0, n - 1
        do k = 0, n - 1
          do x = 0, n - 1
            b3(j + n * i + n2 * k, x) = h2e(((k * n + x) * n + i) * n + j)
          end do
        end do
      end do
    end do

    ! q1(i,a,q,p) from rm2[p,q,i,a]: term 1 contracts rm2's THIRD axis.
    allocate(q1(0:n4 - 1))
    do p = 0, n - 1
      do q = 0, n - 1
        do i = 0, n - 1
          do a = 0, n - 1
            q1(i + n * a + n2 * q + n3 * p) = rm2(((p * n + q) * n + i) * n + a)
          end do
        end do
      end do
    end do

    ! p3(j,i,k, x,q,p) from rm3[p,q,k,i,j,x]
    ! p4(j,i,k, x,q,p) from rm3[p,q,k,i,x,j]
    allocate(p3(0:n6 - 1), p4(0:n6 - 1))
    do p = 0, n - 1
      do q = 0, n - 1
        do k = 0, n - 1
          do i = 0, n - 1
            do j = 0, n - 1
              do x = 0, n - 1
                p3(j + n * i + n2 * k + n3 * (x + n * q + n2 * p)) = &
                  rm3(((((p * n + q) * n + k) * n + i) * n + j) * n + x)
                p4(j + n * i + n2 * k + n3 * (x + n * q + n2 * p)) = &
                  rm3(((((p * n + q) * n + k) * n + i) * n + x) * n + j)
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(rm3)

    allocate(ra(0:n4 - 1), rb(0:n4 - 1))
    ! ra: terms 1, 5 and 3 -- all on flat b+na+n2 q+n3 p
    call dgemm('T', 'N', n, n3, n, -1.0_dp, ha, n, q1, n, 0.0_dp, ra, n)
    call dgemm('T', 'N', n2, n2, n2, -1.0_dp, b5, n2, rm2, n2, 1.0_dp, ra, n2)
    call dgemm('T', 'N', n, n3, n3, -1.0_dp, b3, n3, p3, n3, 1.0_dp, ra, n)
    ! rb: terms 2 and 4 -- both on flat a+nb+n2 q+n3 p
    call dgemm('T', 'N', n, n3, n, -1.0_dp, ha, n, rm2, n, 0.0_dp, rb, n)
    call dgemm('T', 'N', n, n3, n3, -1.0_dp, b3, n3, p4, n3, 1.0_dp, rb, n)

    do p = 0, n - 1
      do q = 0, n - 1
        do a = 0, n - 1
          do b = 0, n - 1
            a7(((p * n + q) * n + a) * n + b) = ra(b + n * a + n2 * q + n3 * p) &
                                              + rb(a + n * b + n2 * q + n3 * p)
          end do
        end do
      end do
    end do

    deallocate(ha, b5, b3, q1, p3, p4, ra, rb)
  end subroutine nevpt2_a7

  !> The a12 Koopmans intermediate (the Sir subspace).
  !>
  !> Same shape as a22, one dimension down: all six terms take dm2/dm3 whose
  !> two SLOWEST axes are the spectator pair, here `qp` -- note the order, the
  !> output is indexed [p,q,a,b] but the RDMs are indexed [q,p,...], so the
  !> gather below reads the spectator as S = q*n + p.  Naming the dm3 block
  !> axes g1..g4 fastest-first, term 4 contracts g1,g2,g3 (dm3 as it stands)
  !> and term 3 contracts g1,g2,g4 (one permuted copy).  On the dm2 side term
  !> 5 contracts both block axes, terms 2+6 the fastest, and term 1 the other
  !> one (one permuted copy).  Terms 2 and 6 fold first:
  !>
  !>     -sum_i h1e[b,i] G[a,i] + sum_ij h2e[j,b,i,j] G[a,i]
  !>   = sum_i ( sum_j h2e[j,b,i,j] - h1e[b,i] ) G[a,i]
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  h1e   active one-electron matrix, C-order [n,n]
  !> @param[in]  h2e   physicist-ordered active ERIs, C-order [n,n,n,n]
  !> @param[in]  dm2   spin-free active 2-RDM, C-order [n,n,n,n]
  !> @param[in]  dm3   spin-free active 3-RDM, C-order [n,n,n,n,n,n]
  !> @param[out] a12   C-order [n,n,n,n]
  subroutine nevpt2_a12(nact, h1e, h2e, dm2, dm3, a12) &
      bind(C, name="nevpt2_a12")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm2(0:*), dm3(0:*)
    real(dp), intent(inout) :: a12(0:*)

    integer :: n, n2, n3, n4, n6
    integer :: p, q, a, b, i, j, k, x, s, e1, e2
    real(dp) :: acc
    real(dp), allocatable :: u(:,:), ha(:,:), b5(:,:), b3(:,:), b4(:,:)
    real(dp), allocatable :: q2(:), p3(:), ra(:), rb(:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n4 = n3 * n
    n6 = n4 * n2

    allocate(u(0:n - 1, 0:n - 1), ha(0:n - 1, 0:n - 1))
    do b = 0, n - 1
      do i = 0, n - 1
        acc = 0.0_dp
        do j = 0, n - 1
          acc = acc + h2e(((j * n + b) * n + i) * n + j)
        end do
        u(i, b) = acc - h1e(b * n + i)
        ha(i, b) = h1e(i * n + b)
      end do
    end do

    ! b5(k+nj, b+na) = h2e[b,j,k,a]
    allocate(b5(0:n2 - 1, 0:n2 - 1))
    do j = 0, n - 1
      do k = 0, n - 1
        do b = 0, n - 1
          do a = 0, n - 1
            b5(k + n * j, b + n * a) = h2e(((b * n + j) * n + k) * n + a)
          end do
        end do
      end do
    end do

    ! b3(k+ni+n2 j, a) = h2e[i,j,k,a];  b4(i+nk+n2 j, b) = h2e[k,b,i,j]
    allocate(b3(0:n3 - 1, 0:n - 1), b4(0:n3 - 1, 0:n - 1))
    do i = 0, n - 1
      do j = 0, n - 1
        do k = 0, n - 1
          do x = 0, n - 1
            b3(k + n * i + n2 * j, x) = h2e(((i * n + j) * n + k) * n + x)
            b4(i + n * k + n2 * j, x) = h2e(((k * n + x) * n + i) * n + j)
          end do
        end do
      end do
    end do

    ! q2(e2,e1,S) from dm2[S,e2,e1];  p3(k,i,j, b,S) from dm3[S,j,b,i,k]
    allocate(q2(0:n4 - 1), p3(0:n6 - 1))
    do s = 0, n2 - 1
      do e2 = 0, n - 1
        do e1 = 0, n - 1
          q2(e2 + n * e1 + n2 * s) = dm2((s * n + e2) * n + e1)
        end do
      end do
      do j = 0, n - 1
        do b = 0, n - 1
          do i = 0, n - 1
            do k = 0, n - 1
              p3(k + n * i + n2 * j + n3 * (b + n * s)) = &
                dm3((((s * n + j) * n + b) * n + i) * n + k)
            end do
          end do
        end do
      end do
    end do

    allocate(ra(0:n4 - 1), rb(0:n4 - 1))
    ! ra: terms 2+6, then 5, then 4 -- all on flat b+na+n2 S
    call dgemm('T', 'N', n, n3, n, 1.0_dp, u, n, dm2, n, 0.0_dp, ra, n)
    call dgemm('T', 'N', n2, n2, n2, -1.0_dp, b5, n2, dm2, n2, 1.0_dp, ra, n2)
    call dgemm('T', 'N', n, n3, n3, -1.0_dp, b4, n3, dm3, n3, 1.0_dp, ra, n)
    ! rb: terms 1 then 3 -- both on flat a+nb+n2 S
    call dgemm('T', 'N', n, n3, n, 1.0_dp, ha, n, q2, n, 0.0_dp, rb, n)
    call dgemm('T', 'N', n, n3, n3, 1.0_dp, b3, n3, p3, n3, 1.0_dp, rb, n)

    do p = 0, n - 1
      do q = 0, n - 1
        s = q * n + p                      ! the RDMs are indexed [q,p,...]
        do a = 0, n - 1
          do b = 0, n - 1
            a12(((p * n + q) * n + a) * n + b) = ra(b + n * a + n2 * s) &
                                               + rb(a + n * b + n2 * s)
          end do
        end do
      end do
    end do

    deallocate(u, ha, b5, b3, b4, q2, p3, ra, rb)
  end subroutine nevpt2_a12

  !> The a13 Koopmans intermediate (the Sir subspace).
  !>
  !> Unlike a12 this one has only ONE spectator: the RDMs are indexed [q,...]
  !> and `p` sits inside them, so the natural unit of work is a q block rather
  !> than a (q,p) block.  Only the two n^7 terms (5 and 6, the dm3 ones) are
  !> worth a GEMM; each contracts three axes that are not adjacent in dm3, so
  !> each takes a permuted copy.  Everything else is n^6 or smaller and is
  !> evaluated in the gather.
  !>
  !> Three of the twelve terms fold away before any arithmetic.  With
  !>
  !>     v[b,i] = h1e[b,i] - sum_l h2e[l,b,i,l]
  !>
  !> terms 3 and 10 are sum_i v[b,i] dm2[q,i,a,p], and terms 4 and 12 -- which
  !> both live on the a == p plane together with term 11 -- collapse into
  !>
  !>     W[b] = -2 sum_i v[b,i] dm1[q,i]
  !>            -2 sum_klm h2e[m,b,k,l] dm2[q,l,m,k]
  !>
  !> added on that plane alone.
  !>
  !> @param[in]  nact  active orbitals
  !> @param[in]  h1e   active one-electron matrix, C-order [n,n]
  !> @param[in]  h2e   physicist-ordered active ERIs, C-order [n,n,n,n]
  !> @param[in]  dm1   spin-free active 1-RDM, C-order [n,n]
  !> @param[in]  dm2   spin-free active 2-RDM, C-order [n,n,n,n]
  !> @param[in]  dm3   spin-free active 3-RDM, C-order [n,n,n,n,n,n]
  !> @param[out] a13   C-order [n,n,n,n]
  subroutine nevpt2_a13(nact, h1e, h2e, dm1, dm2, dm3, a13) &
      bind(C, name="nevpt2_a13")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm1(0:*), dm2(0:*), dm3(0:*)
    real(dp), intent(inout) :: a13(0:*)

    integer :: n, n2, n3, n4, n5, n6
    integer :: p, q, a, b, i, j, k, l, m, x, y
    integer :: d1b, d2b
    real(dp) :: acc, val
    real(dp), allocatable :: v(:,:), b5(:,:), b6(:,:), b7(:,:), b8(:,:)
    real(dp), allocatable :: p5(:), p6(:), q7(:), r5(:), r6(:), r7(:), w(:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n4 = n3 * n
    n5 = n4 * n
    n6 = n5 * n

    allocate(v(0:n - 1, 0:n - 1))
    do b = 0, n - 1
      do i = 0, n - 1
        acc = 0.0_dp
        do l = 0, n - 1
          acc = acc + h2e(((l * n + b) * n + i) * n + l)
        end do
        v(i, b) = h1e(b * n + i) - acc
      end do
    end do

    ! b5(k+ni+n2 j, a) = h2e[i,j,k,a];  b6(i+nk+n2 j, b) = h2e[k,b,i,j]
    allocate(b5(0:n3 - 1, 0:n - 1), b6(0:n3 - 1, 0:n - 1))
    do i = 0, n - 1
      do j = 0, n - 1
        do k = 0, n - 1
          do x = 0, n - 1
            b5(k + n * i + n2 * j, x) = h2e(((i * n + j) * n + k) * n + x)
            b6(i + n * k + n2 * j, x) = h2e(((k * n + x) * n + i) * n + j)
          end do
        end do
      end do
    end do

    ! p5(k,i,j, p,b, q) from dm3[q,b,j,p,i,k]   (term 5, contract i,j,k)
    ! p6(i,k,j, p,a, q) from dm3[q,j,a,p,k,i]   (term 6, contract i,j,k)
    ! One sweep builds both: x plays b for p5 and j for p6, y the other way.
    allocate(p5(0:n6 - 1), p6(0:n6 - 1))
    do q = 0, n - 1
      do x = 0, n - 1
        do y = 0, n - 1
          do p = 0, n - 1
            do i = 0, n - 1
              do k = 0, n - 1
                p5(k + n * i + n2 * y + n3 * (p + n * x + n2 * q)) = &
                  dm3(((((q * n + x) * n + y) * n + p) * n + i) * n + k)
                p6(i + n * k + n2 * x + n3 * (p + n * y + n2 * q)) = &
                  dm3(((((q * n + x) * n + y) * n + p) * n + k) * n + i)
              end do
            end do
          end do
        end do
      end do
    end do

    ! Terms 7 and 8 are only n^6, but n^6 of SCALAR work is what a NumPy
    ! einsum does with a BLAS call, so they get GEMMs too -- leaving them in
    ! the gather is what made the first revision of this routine lose to
    ! NumPy at nact=8.  Term 8 contracts dm2's two fastest block axes (dm2 as
    ! it stands) and lands on term 5's flat offset; term 7 contracts the two
    ! slowest and needs the one dm2 repack below.
    ! b7(l+nm, b+na) = h2e[b,l,m,a];  b8(m+nk, a+np) = 2 h2e[k,p,m,a]
    allocate(b7(0:n2 - 1, 0:n2 - 1), b8(0:n2 - 1, 0:n2 - 1))
    do l = 0, n - 1
      do m = 0, n - 1
        do x = 0, n - 1
          do y = 0, n - 1
            b7(l + n * m, x + n * y) = h2e(((x * n + l) * n + m) * n + y)
            b8(m + n * l, y + n * x) = 2.0_dp * h2e(((l * n + x) * n + m) * n + y)
          end do
        end do
      end do
    end do

    ! q7(e2,e3,e1,q) from dm2[q,e3,e2,e1]
    allocate(q7(0:n4 - 1))
    do q = 0, n - 1
      do k = 0, n - 1
        do l = 0, n - 1
          do m = 0, n - 1
            q7(l + n * k + n2 * m + n3 * q) = dm2(((q * n + k) * n + l) * n + m)
          end do
        end do
      end do
    end do

    allocate(r5(0:n4 - 1), r6(0:n4 - 1), r7(0:n4 - 1))
    call dgemm('T', 'N', n, n3, n3, -1.0_dp, b5, n3, p5, n3, 0.0_dp, r5, n)
    call dgemm('T', 'N', n2, n2, n2, 1.0_dp, b8, n2, dm2, n2, 1.0_dp, r5, n2)
    call dgemm('T', 'N', n, n3, n3, 1.0_dp, b6, n3, p6, n3, 0.0_dp, r6, n)
    call dgemm('T', 'N', n2, n2, n2, 1.0_dp, b7, n2, q7, n2, 0.0_dp, r7, n2)
    deallocate(p5, p6, q7)

    allocate(w(0:n - 1))
    do q = 0, n - 1
      d1b = q * n
      d2b = q * n3

      ! W[b]: terms 4, 11 and 12, which only reach the a == p plane
      do b = 0, n - 1
        acc = 0.0_dp
        do i = 0, n - 1
          acc = acc - 2.0_dp * v(i, b) * dm1(d1b + i)
        end do
        do k = 0, n - 1
          do l = 0, n - 1
            do m = 0, n - 1
              acc = acc - 2.0_dp * h2e(((m * n + b) * n + k) * n + l) &
                        * dm2(d2b + (l * n + m) * n + k)
            end do
          end do
        end do
        w(b) = acc
      end do

      do p = 0, n - 1
        do a = 0, n - 1
          do b = 0, n - 1
            val = r5(a + n * p + n2 * b + n3 * q) &
                + r6(b + n * p + n2 * a + n3 * q) &
                + r7(b + n * a + n2 * p + n3 * q)
            do i = 0, n - 1
              ! term 1     -h1e[i,a] dm2[q,b,i,p]
              ! terms 3+10 +v[b,i]   dm2[q,i,a,p]
              val = val - h1e(i * n + a) * dm2(d2b + (b * n + i) * n + p) &
                        + v(i, b) * dm2(d2b + (i * n + a) * n + p)
            end do
            ! term 2 +2 h1e[p,a] dm1[q,b]
            val = val + 2.0_dp * h1e(p * n + a) * dm1(d1b + b)
            do m = 0, n - 1
              ! term 9 -2 h2e[b,p,m,a] dm1[q,m]
              val = val - 2.0_dp * h2e(((b * n + p) * n + m) * n + a) &
                        * dm1(d1b + m)
            end do
            if (a == p) val = val + w(b)
            a13(((p * n + q) * n + a) * n + b) = val
          end do
        end do
      end do
    end do

    deallocate(v, b5, b6, b7, b8, r5, r6, r7, w)
  end subroutine nevpt2_a13

  subroutine nevpt2_hdm1(nact, dm1, hdm1) bind(C, name="nevpt2_hdm1")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: dm1(0:*)
    real(dp), intent(inout) :: hdm1(0:*)
    integer :: n, p, a
    n = int(nact)
    if (n <= 0) return
    do a = 0, n - 1
      do p = 0, n - 1
        hdm1(p + a*n) = -dm1(a + p*n)
      end do
      hdm1(a + a*n) = hdm1(a + a*n) + 2.0_dp
    end do
  end subroutine nevpt2_hdm1

  subroutine nevpt2_hdm2(nact, dm1, dm2, hdm2) bind(C, name="nevpt2_hdm2")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: dm1(0:*), dm2(0:*)
    real(dp), intent(inout) :: hdm2(0:*)
    integer :: n, i, j, k, l
    n = int(nact)
    if (n <= 0) return
    !$omp parallel do collapse(2) default(none) shared(n, dm1, dm2, hdm2) private(i,j,k,l)
    do j = 0, n - 1
      do i = 0, n - 1
        do l = 0, n - 1
          do k = 0, n - 1
            hdm2(((i*n + j)*n + k)*n + l) = dm2(((k*n + i)*n + l)*n + j)
            if (j == k) then
              hdm2(((i*n + j)*n + k)*n + l) = hdm2(((i*n + j)*n + k)*n + l) + dm1(l*n + i)
            end if
            if (i == k) then
              hdm2(((i*n + j)*n + k)*n + l) = hdm2(((i*n + j)*n + k)*n + l) - 2.0_dp * dm1(l*n + j)
            end if
            if (j == l) then
              hdm2(((i*n + j)*n + k)*n + l) = hdm2(((i*n + j)*n + k)*n + l) - 2.0_dp * dm1(k*n + i)
            end if
            if (i == l .and. j == k) then
              hdm2(((i*n + j)*n + k)*n + l) = hdm2(((i*n + j)*n + k)*n + l) - 2.0_dp
            end if
            if (i == k .and. j == l) then
              hdm2(((i*n + j)*n + k)*n + l) = hdm2(((i*n + j)*n + k)*n + l) + 4.0_dp
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine nevpt2_hdm2

  subroutine nevpt2_a17(nact, h1e, h2e, dm2, dm3, a17) &
      bind(C, name="nevpt2_a17")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm2(0:*), dm3(0:*)
    real(dp), intent(inout) :: a17(0:*)
    integer :: n, a, b, c, p, i, j, k
    real(dp) :: s
    n = int(nact)
    if (n <= 0) return
    do a = 0, n - 1
      do b = 0, n - 1
        do c = 0, n - 1
          do p = 0, n - 1
            s = 0.0_dp
            do i = 0, n - 1
              s = s - h1e(p + i*n) * dm2((((c*n + a)*n + b)*n + i))
            end do
            do i = 0, n - 1
              do j = 0, n - 1
                do k = 0, n - 1
                  s = s - h2e(((k*n + p)*n + i)*n + j) &
                      * dm3((((c*n + a)*n + b)*n*n + j*n + k)*n + i)
                end do
              end do
            end do
            a17((((a*n + b)*n + c)*n + p)) = s
          end do
        end do
      end do
    end do
  end subroutine nevpt2_a17

  subroutine nevpt2_a19(nact, h1e, h2e, dm1, dm2, a19) &
      bind(C, name="nevpt2_a19")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm1(0:*), dm2(0:*)
    real(dp), intent(inout) :: a19(0:*)
    integer :: n, a, p, i, j, k
    real(dp) :: s
    n = int(nact)
    if (n <= 0) return
    do a = 0, n - 1
      do p = 0, n - 1
        s = 0.0_dp
        do i = 0, n - 1
          s = s - h1e(p + i*n) * dm1(a + i*n)
        end do
        do i = 0, n - 1
          do j = 0, n - 1
            do k = 0, n - 1
              s = s - h2e(((k*n + p)*n + i)*n + j) &
                  * dm2((((a*n + j)*n + k)*n + i))
            end do
          end do
        end do
        a19(p + a*n) = s
      end do
    end do
  end subroutine nevpt2_a19

  subroutine nevpt2_a23(nact, h1e, h2e, dm1, dm2, dm3, a23) &
      bind(C, name="nevpt2_a23")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm1(0:*), dm2(0:*), dm3(0:*)
    real(dp), intent(inout) :: a23(0:*)
    integer :: n, a, b, c, p, i, j, k
    real(dp) :: s
    n = int(nact)
    if (n <= 0) return
    do a = 0, n - 1
      do b = 0, n - 1
        do c = 0, n - 1
          do p = 0, n - 1
            s = 0.0_dp
            do i = 0, n - 1
              s = s - h1e(i + p*n) * dm2((((c*n + a)*n + i)*n + b))
            end do
            do i = 0, n - 1
              do j = 0, n - 1
                do k = 0, n - 1
                  s = s - h2e(((p*n + i)*n + j)*n + k) &
                      * dm3((((c*n + a)*n + j)*n*n + b*n + i)*n + k)
                end do
              end do
            end do
            s = s + 2.0_dp * h1e(b + p*n) * dm1(c + a*n)
            do i = 0, n - 1
              do k = 0, n - 1
                s = s + 2.0_dp * h2e(((p*n + i)*n + b)*n + k) &
                    * dm2((((c*n + a)*n + i)*n + k))
              end do
            end do
            a23((((a*n + b)*n + c)*n + p)) = s
          end do
        end do
      end do
    end do
  end subroutine nevpt2_a23

  subroutine nevpt2_a25(nact, h1e, h2e, dm1, dm2, a25) &
      bind(C, name="nevpt2_a25")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm1(0:*), dm2(0:*)
    real(dp), intent(inout) :: a25(0:*)
    integer :: n, a, p, i, j, k
    real(dp) :: s
    n = int(nact)
    if (n <= 0) return
    do a = 0, n - 1
      do p = 0, n - 1
        s = 0.0_dp
        do i = 0, n - 1
          s = s - h1e(p + i*n) * dm1(a + i*n)
        end do
        do i = 0, n - 1
          do j = 0, n - 1
            do k = 0, n - 1
              s = s - h2e(((p*n + i)*n + j)*n + k) &
                  * dm2((((j*n + a)*n + i)*n + k))
            end do
          end do
        end do
        s = s + 2.0_dp * h1e(p + a*n)
        do i = 0, n - 1
          do j = 0, n - 1
            s = s + 2.0_dp * h2e(((p*n + i)*n + a)*n + j) &
                * dm1(i + j*n)
          end do
        end do
        a25(p + a*n) = s
      end do
    end do
  end subroutine nevpt2_a25

  subroutine nevpt2_a3(nact, h1e, h2e, dm1, dm2, hdm1, a3) &
      bind(C, name="nevpt2_a3")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm1(0:*), dm2(0:*), hdm1(0:*)
    real(dp), intent(inout) :: a3(0:*)
    integer :: n, p, a, i, j, k
    real(dp) :: s
    n = int(nact)
    if (n <= 0) return
    do p = 0, n - 1
      do a = 0, n - 1
        s = 0.0_dp
        do i = 0, n - 1
          s = s + h1e(i + a*n) * hdm1(i + p*n)
        end do
        do i = 0, n - 1
          do k = 0, n - 1
            s = s + 2.0_dp * h2e(((i*n + p)*n + k)*n + a) * dm1(i + k*n)
          end do
        end do
        do i = 0, n - 1
          do j = 0, n - 1
            do k = 0, n - 1
              s = s - h2e(((i*n + j)*n + k)*n + a) &
                  * dm2((((j*n + p)*n + i)*n + k))
            end do
          end do
        end do
        a3(p + a*n) = s
      end do
    end do
  end subroutine nevpt2_a3

  subroutine nevpt2_k27(nact, h1e, h2e, dm1, dm2, k27) &
      bind(C, name="nevpt2_k27")
    integer(c_int32_t), value :: nact
    real(dp), intent(in) :: h1e(0:*), h2e(0:*), dm1(0:*), dm2(0:*)
    real(dp), intent(inout) :: k27(0:*)
    integer :: n, p, a, i, j, k
    real(dp) :: s
    n = int(nact)
    if (n <= 0) return
    do p = 0, n - 1
      do a = 0, n - 1
        s = 0.0_dp
        do i = 0, n - 1
          s = s - h1e(a + i*n) * dm1(p + i*n)
        end do
        do i = 0, n - 1
          do j = 0, n - 1
            do k = 0, n - 1
              s = s - h2e(((i*n + a)*n + j)*n + k) &
                  * dm2((((p*n + k)*n + i)*n + j))
            end do
          end do
        end do
        do j = 0, n - 1
          do i = 0, n - 1
            s = s + h2e(((i*n + a)*n + j)*n + i) * dm1(p + j*n)
          end do
        end do
        k27(p + a*n) = s
      end do
    end do
  end subroutine nevpt2_k27

! nevpt2_sijrs — Fortran implementation of the Sijrs subspace
! (core,core;virt,virt) doubly-occupied→virtual SC-NEVPT2 correction.
!
! norm and energy for the Sijrs subspace (i/core, j/core, a/virt, b/virt):
!
!   norm_i = sum_jab  g_iajb * (2*g_iajb - g_ibja)
!   e_i    = sum_jab  g_iajb / (e_i - e_a + e_j - e_b) * (2*g_iajb - g_ibja)
!
! g(i,a,j,b) = eri_mo(i,a,j,b) in chemist notation (core,virt,core,virt).
! e_core(i/j) and e_virt(a/b) are semicanonical orbital energies.
!
subroutine nevpt2_sijrs(ncore, nvirt, g, e_core, e_virt, norm, energy) &
    bind(C, name="nevpt2_sijrs")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int32_t), value :: ncore, nvirt
  ! C-order [ncore, nvirt, ncore, nvirt] passed as flat 1-D
  real(c_double), intent(in)  :: g(0:*)
  real(c_double), intent(in)  :: e_core(0:*), e_virt(0:*)
  real(c_double), intent(out) :: norm, energy

  real(c_double), parameter :: eps = 1.0e-12_c_double
  integer :: i, j, a, b
  real(c_double) :: gi, gj_ab_swapped, theta, denom, t2i

  if (ncore <= 0 .or. nvirt <= 0) then
    norm = 0.0_c_double
    energy = 0.0_c_double
    return
  end if

  norm = 0.0_c_double
  energy = 0.0_c_double

  !$omp parallel do collapse(4) default(none) &
  !$omp shared(ncore, nvirt, g, e_core, e_virt) &
  !$omp private(i, j, a, b, gi, gj_ab_swapped, theta, denom, t2i) &
  !$omp reduction(+:norm, energy)
  do i = 0, ncore - 1
    do j = 0, ncore - 1
      do a = 0, nvirt - 1
        do b = 0, nvirt - 1
          ! C-order flat offset: ((i*nvirt + a)*ncore + j)*nvirt + b
          gi = g(((i*nvirt + a)*ncore + j)*nvirt + b)
          gj_ab_swapped = g(((i*nvirt + b)*ncore + j)*nvirt + a)
          theta = 2.0_c_double * gi - gj_ab_swapped
          norm = norm + gi * theta
          denom = e_core(i) - e_virt(a) + e_core(j) - e_virt(b)
          if (abs(denom) > eps) then
            t2i = gi / denom
            energy = energy + t2i * theta
          end if
        end do
      end do
    end do
  end do
  !$omp end parallel do
end subroutine nevpt2_sijrs


! nevpt2_srs — Fortran implementation of the Srs subspace final contraction.
!
!   norm = 0.5 * sum_{pqrsab}  h2e_v(r,s,q,p) * h2e_v(r,s,b,a) * rm2(p,q,b,a)
!   h    = 0.5 * sum_{pqrsab}  h2e_v(r,s,q,p) * h2e_v(r,s,b,a) * a7(p,q,a,b)
!   diff(r,s) = e_virt(r) + e_virt(s)  (computed in Python, used in _norm_to_energy)
!
! h2e_v is (nvirt, ncas, nvirt, ncas) = phys-ordered ERI block (V,A,V,A).
! rm2 and a7 are (ncas, ncas, ncas, ncas) from _a7.
!
subroutine nevpt2_srs(nact, nvirt, h2e_v, rm2, a7, norm, energy) &
    bind(C, name="nevpt2_srs")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int32_t), value :: nact, nvirt
  ! All arrays are C-order, passed as flat 1-D:
  !   h2e_v: [nvirt, nact, nvirt, nact]
  !   rm2, a7: [nact, nact, nact, nact]
  !   norm, energy: [nvirt, nvirt]
  real(c_double), intent(in)  :: h2e_v(0:*)
  real(c_double), intent(in)  :: rm2(0:*), a7(0:*)
  real(c_double), intent(out) :: norm(0:*), energy(0:*)

  real(c_double), parameter :: half = 0.5_c_double
  integer :: p, q, r, s, a, b
  real(c_double) :: h2e_rsqp, h2e_rsba, rm2_pqba, a7_pqab
  ! stride helpers for h2e_v
  integer :: snva, sna2
  snva = nvirt * nact   ! stride for dim 2 of h2e_v (nvirt * nact)
  sna2 = snva * nvirt   ! stride for dim 3 of h2e_v (nvirt * nact * nvirt)

  norm(0:nvirt*nvirt - 1) = 0.0_c_double
  energy(0:nvirt*nvirt - 1) = 0.0_c_double
  if (nvirt <= 0 .or. nact <= 0) return

  !$omp parallel do collapse(2) default(none) &
  !$omp shared(nact, nvirt, h2e_v, rm2, a7, norm, energy, snva, sna2) &
  !$omp private(r, s, p, q, a, b, h2e_rsqp, h2e_rsba, rm2_pqba, a7_pqab)
  do r = 0, nvirt - 1
    do s = 0, nvirt - 1
      do p = 0, nact - 1
        do q = 0, nact - 1
          ! h2e_v[r, q, s, p] = h2e_phys[r, s, q, p]  (after transpose(0,2,1,3))
          ! C-order offset: ((r*nact + q)*nvirt + s)*nact + p
          h2e_rsqp = h2e_v(((r*nact + q)*nvirt + s)*nact + p)
          do a = 0, nact - 1
            do b = 0, nact - 1
              ! h2e_v[r, b, s, a] = h2e_phys[r, s, b, a]  (after transpose(0,2,1,3))
              ! C-order offset: ((r*nact + b)*nvirt + s)*nact + a
              h2e_rsba = h2e_v(((r*nact + b)*nvirt + s)*nact + a)
              ! rm2[p, q, b, a]  C-order: (((p*nact + q)*nact + b)*nact + a
              rm2_pqba = rm2(((p*nact + q)*nact + b)*nact + a)
              ! a7[p, q, a, b]  C-order
              a7_pqab  = a7(((p*nact + q)*nact + a)*nact + b)
              norm(r*nvirt + s)   = norm(r*nvirt + s)   + half * h2e_rsqp * h2e_rsba * rm2_pqba
              energy(r*nvirt + s) = energy(r*nvirt + s) + half * h2e_rsqp * h2e_rsba * a7_pqab
            end do
          end do
        end do
      end do
    end do
  end do
  !$omp end parallel do
end subroutine nevpt2_srs


! nevpt2_sij — Fortran implementation of the Sij subspace final contraction.
!
!   norm = 0.5 * sum_{pqijab}  h2e_v(q,p,i,j) * h2e_v(b,a,i,j) * hdm2(p,q,a,b)
!   h    = 0.5 * sum_{pqijab}  h2e_v(q,p,i,j) * h2e_v(b,a,i,j) * a9(p,q,a,b)
!   diff(i,j) = e_core(i) + e_core(j)  (computed in Python, used in _norm_to_energy)
!
! h2e_v is (ncas, ncore, ncas, ncore) = phys-ordered ERI block (A,C,A,C).
! hdm2 and a9 are (ncas, ncas, ncas, ncas).
!
subroutine nevpt2_sij(nact, ncore, h2e_v, hdm2, a9, norm, energy) &
    bind(C, name="nevpt2_sij")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int32_t), value :: nact, ncore
  ! All arrays are C-order, passed as flat 1-D:
  !   h2e_v: [nact, ncore, nact, ncore]
  !   hdm2, a9: [nact, nact, nact, nact]
  !   norm, energy: [ncore, ncore]
  real(c_double), intent(in)  :: h2e_v(0:*)
  real(c_double), intent(in)  :: hdm2(0:*), a9(0:*)
  real(c_double), intent(out) :: norm(0:*), energy(0:*)

  real(c_double), parameter :: half = 0.5_c_double
  integer :: p, q, a, b, i, j
  real(c_double) :: h2e_qpij, h2e_baij, hdm2_pqab, a9_pqab

  norm(0:ncore*ncore - 1) = 0.0_c_double
  energy(0:ncore*ncore - 1) = 0.0_c_double
  if (ncore <= 0 .or. nact <= 0) return

  !$omp parallel do collapse(2) default(none) &
  !$omp shared(nact, ncore, h2e_v, hdm2, a9, norm, energy) &
  !$omp private(i, j, p, q, a, b, h2e_qpij, h2e_baij, hdm2_pqab, a9_pqab)
  do i = 0, ncore - 1
    do j = 0, ncore - 1
      do p = 0, nact - 1
        do q = 0, nact - 1
          ! h2e_v[q, i, p, j] = h2e_phys[q, p, i, j]  (after transpose(0,2,1,3))
          ! C-order offset: (((q*ncore + i)*nact + p)*ncore + j)
          h2e_qpij = h2e_v((((q*ncore + i)*nact + p)*ncore + j))
          do a = 0, nact - 1
            do b = 0, nact - 1
              ! h2e_v[b, i, a, j] = h2e_phys[b, a, i, j]  (after transpose(0,2,1,3))
              ! C-order offset: (((b*ncore + i)*nact + a)*ncore + j)
              h2e_baij = h2e_v((((b*ncore + i)*nact + a)*ncore + j))
              ! hdm2[p, q, a, b]  C-order: (((p*nact + q)*nact + a)*nact + b)
              hdm2_pqab = hdm2((((p*nact + q)*nact + a)*nact + b))
              ! a9[p, q, a, b]  C-order
              a9_pqab   = a9((((p*nact + q)*nact + a)*nact + b))
              norm(i*ncore + j)   = norm(i*ncore + j)   + half * h2e_qpij * h2e_baij * hdm2_pqab
              energy(i*ncore + j) = energy(i*ncore + j) + half * h2e_qpij * h2e_baij * a9_pqab
            end do
          end do
        end do
      end do
    end do
  end do
  !$omp end parallel do
end subroutine nevpt2_sij


! nevpt2_sijr — Fortran implementation of the Sijr subspace (virt,core,core).
!
! h2e_v: C-order [nvirt, nact, ncore, ncore] from blocks['Sijr']
!   = phys(eri_mo[V, C, A, C]) → eri_mo[V, A, C, C]
!
! norm[r,j,i] = sum_{p,a} h2e_v[r,p,j,i] * (2*h2e_v[r,a,j,i] - h2e_v[r,a,i,j]) * hdm1[p,a]
!
! Same expression with a3 replacing hdm1 for energy.
!
subroutine nevpt2_sijr(nvirt, ncore, nact, h2e_v, hdm1, a3, norm, energy) &
    bind(C, name="nevpt2_sijr")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int32_t), value :: nvirt, ncore, nact
  ! h2e_v: [nvirt, nact, ncore, ncore]; hdm1, a3: [nact, nact]
  real(c_double), intent(in)  :: h2e_v(0:*), hdm1(0:*), a3(0:*)
  real(c_double), intent(out) :: norm(0:*), energy(0:*)

  integer :: r, i, j, p, a, nc2, nr2
  real(c_double) :: h2e_rpji, t

  nc2 = ncore * ncore
  nr2 = nvirt * nc2
  norm(0:nr2 - 1) = 0.0_c_double
  energy(0:nr2 - 1) = 0.0_c_double
  if (ncore <= 0 .or. nvirt <= 0 .or. nact <= 0) return

  !$omp parallel do collapse(3) default(none) &
  !$omp shared(nvirt, ncore, nact, h2e_v, hdm1, a3, norm, energy, nc2) &
  !$omp private(r, i, j, p, a, h2e_rpji, t)
  do r = 0, nvirt - 1
    do j = 0, ncore - 1
      do i = 0, ncore - 1
        do p = 0, nact - 1
          h2e_rpji = h2e_v(((r*nact + p)*ncore + j)*ncore + i)
          t = 0.0_c_double
          do a = 0, nact - 1
            ! 2*h2e_v[r,a,j,i] - h2e_v[r,a,i,j]
            t = t + (2.0_c_double * h2e_v(((r*nact + a)*ncore + j)*ncore + i) &
                    - h2e_v(((r*nact + a)*ncore + i)*ncore + j)) * hdm1(p*nact + a)
          end do
          norm(r*nc2 + j*ncore + i) = norm(r*nc2 + j*ncore + i) + h2e_rpji * t
        end do
        do p = 0, nact - 1
          h2e_rpji = h2e_v(((r*nact + p)*ncore + j)*ncore + i)
          t = 0.0_c_double
          do a = 0, nact - 1
            t = t + (2.0_c_double * h2e_v(((r*nact + a)*ncore + j)*ncore + i) &
                    - h2e_v(((r*nact + a)*ncore + i)*ncore + j)) * a3(p*nact + a)
          end do
          energy(r*nc2 + j*ncore + i) = energy(r*nc2 + j*ncore + i) + h2e_rpji * t
        end do
      end do
    end do
  end do
  !$omp end parallel do

  ! Symmetrise (r,i,j) += (r,j,i) and halve diagonal
  do r = 0, nvirt - 1
    do i = 0, ncore - 1
      do j = i + 1, ncore - 1
        norm(r*nc2 + i*ncore + j) = norm(r*nc2 + i*ncore + j) + norm(r*nc2 + j*ncore + i)
        norm(r*nc2 + j*ncore + i) = norm(r*nc2 + i*ncore + j)
        energy(r*nc2 + i*ncore + j) = energy(r*nc2 + i*ncore + j) + energy(r*nc2 + j*ncore + i)
        energy(r*nc2 + j*ncore + i) = energy(r*nc2 + i*ncore + j)
      end do
      norm(r*nc2 + i*ncore + i) = norm(r*nc2 + i*ncore + i) * 0.5_c_double
      energy(r*nc2 + i*ncore + i) = energy(r*nc2 + i*ncore + i) * 0.5_c_double
    end do
  end do
end subroutine nevpt2_sijr


! nevpt2_srsi — Fortran implementation of the Srsi subspace (virt,virt,core).
!
! h2e_v: C-order [nvirt, nvirt, ncore, nact] from blocks['Srsi']
!   = phys(eri_mo[V, C, V, A]) → eri_mo[V, V, C, A]
!
! norm[r,s,i] = sum_{p,a} h2e_v[r,s,i,p] * (2*h2e_v[r,s,i,a] - h2e_v[s,r,i,a]) * dm1[p,a]
!
! Same expression with k27 replacing dm1 for energy.
!
subroutine nevpt2_srsi(nvirt, ncore, nact, h2e_v, dm1, k27, norm, energy) &
    bind(C, name="nevpt2_srsi")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int32_t), value :: nvirt, ncore, nact
  ! h2e_v: [nvirt, nvirt, ncore, nact]; dm1, k27: [nact, nact]
  real(c_double), intent(in)  :: h2e_v(0:*), dm1(0:*), k27(0:*)
  real(c_double), intent(out) :: norm(0:*), energy(0:*)

  integer :: r, s, i, p, a, nvnc, nv2, nv2nc
  real(c_double) :: h2e_rsip, t

  nvnc = nvirt * ncore
  nv2 = nvirt * nvirt
  nv2nc = nv2 * ncore
  norm(0:nv2nc - 1) = 0.0_c_double
  energy(0:nv2nc - 1) = 0.0_c_double
  if (ncore <= 0 .or. nvirt <= 0 .or. nact <= 0) return

  !$omp parallel do collapse(3) default(none) &
  !$omp shared(nvirt, ncore, nact, h2e_v, dm1, k27, norm, energy, nvnc) &
  !$omp private(r, s, i, p, a, h2e_rsip, t)
  do r = 0, nvirt - 1
    do s = 0, nvirt - 1
      do i = 0, ncore - 1
        do p = 0, nact - 1
          h2e_rsip = h2e_v(((r*nvirt + s)*ncore + i)*nact + p)
          t = 0.0_c_double
          do a = 0, nact - 1
            ! 2*h2e_v[r,s,i,a] - h2e_v[s,r,i,a]
            t = t + (2.0_c_double * h2e_v(((r*nvirt + s)*ncore + i)*nact + a) &
                    - h2e_v(((s*nvirt + r)*ncore + i)*nact + a)) * dm1(p*nact + a)
          end do
          norm((r*nvirt + s)*ncore + i) = norm((r*nvirt + s)*ncore + i) + h2e_rsip * t
        end do
        do p = 0, nact - 1
          h2e_rsip = h2e_v(((r*nvirt + s)*ncore + i)*nact + p)
          t = 0.0_c_double
          do a = 0, nact - 1
            t = t + (2.0_c_double * h2e_v(((r*nvirt + s)*ncore + i)*nact + a) &
                    - h2e_v(((s*nvirt + r)*ncore + i)*nact + a)) * k27(p*nact + a)
          end do
          energy((r*nvirt + s)*ncore + i) = energy((r*nvirt + s)*ncore + i) + h2e_rsip * t
        end do
      end do
    end do
  end do
  !$omp end parallel do

  ! Symmetrise (s,r,i) += (r,s,i) and halve diagonal virt pairs
  do r = 0, nvirt - 1
    do s = r + 1, nvirt - 1
      do i = 0, ncore - 1
        norm((s*nvirt + r)*ncore + i) = norm((s*nvirt + r)*ncore + i) + norm((r*nvirt + s)*ncore + i)
        norm((r*nvirt + s)*ncore + i) = norm((s*nvirt + r)*ncore + i)
        energy((s*nvirt + r)*ncore + i) = energy((s*nvirt + r)*ncore + i) + energy((r*nvirt + s)*ncore + i)
        energy((r*nvirt + s)*ncore + i) = energy((s*nvirt + r)*ncore + i)
      end do
    end do
    do i = 0, ncore - 1
      norm((r*nvirt + r)*ncore + i) = norm((r*nvirt + r)*ncore + i) * 0.5_c_double
      energy((r*nvirt + r)*ncore + i) = energy((r*nvirt + r)*ncore + i) * 0.5_c_double
    end do
  end do
end subroutine nevpt2_srsi


! nevpt2_sir — Fortran implementation of the Sir subspace (virt,core).
!
! h2e_v1: C-order [nvirt, nact, ncore, nact] from blocks['Sir'][0]
!   = phys(eri_mo[V, C, A, A]) → eri_mo[V, A, C, A]
! h2e_v2: C-order [nvirt, nact, nact, ncore] from blocks['Sir'][1]
!   = phys(eri_mo[V, A, A, C]) → eri_mo[V, A, A, C]
! h1e_v: C-order [nvirt, ncore] from blocks['Sir'][2]
!
! norm[i,r] = sum_{pqab}
!   2.0 * v1[r,p,i,q] * v1[r,a,i,b] * dm2[q,p,a,b]
! - 1.0 * v1[r,p,i,q] * v2[r,a,b,i] * dm2[q,p,a,b]
! - 1.0 * v2[r,p,q,i] * v1[r,a,i,b] * dm2[q,p,a,b]
! + 2.0 * v2[r,a,q,i] * v2[r,a,b,i] * dm1[q,b]
! - 1.0 * v2[r,p,q,i] * v2[r,a,b,i] * dm2[q,b,a,p]
! + 1.0 * v2[r,p,q,i] * v2[r,a,a,i] * dm1[q,p]
! + 4.0 * v1[r,p,i,q] * h1e_v[r,i] * dm1[q,p]
! - 2.0 * v2[r,p,q,i] * h1e_v[r,i] * dm1[q,p]
! + 2.0 * h1e_v[r,i] * h1e_v[r,i]
!
! Same expression with a12/a13 replacing dm2/dm1 for energy.
!
subroutine nevpt2_sir(nvirt, ncore, nact, &
    h2e_v1, h2e_v2, h1e_v, dm1, dm2, a12, a13, norm, energy) &
    bind(C, name="nevpt2_sir")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int32_t), value :: nvirt, ncore, nact
  ! h2e_v1: [nvirt, nact, ncore, nact]; h2e_v2: [nvirt, nact, nact, ncore]
  ! h1e_v: [nvirt, ncore]
  ! dm1, dm2, a12, a13: active-space intermediates
  real(c_double), intent(in)  :: h2e_v1(0:*), h2e_v2(0:*), h1e_v(0:*)
  real(c_double), intent(in)  :: dm1(0:*), dm2(0:*), a12(0:*), a13(0:*)
  real(c_double), intent(out) :: norm(0:*), energy(0:*)

  integer(c_int32_t) :: r, i, p, q, a, b, nr
  real(c_double) :: h1e_ri, nsum, esum, v2_diag_ri
  real(c_double) :: v1_rpiq, v2_rpqi, v1_raib, v2_rabi, v2_raqi
  real(c_double) :: dm1_qp, dm2_qpab, dm2_qbap
  real(c_double) :: a12_qpab, a13_pqab

  nr = nvirt * ncore
  norm(0:nr - 1) = 0.0_c_double
  energy(0:nr - 1) = 0.0_c_double
  if (ncore <= 0 .or. nvirt <= 0 .or. nact <= 0) return

  !$omp parallel do collapse(2) default(none) &
  !$omp shared(nvirt, ncore, nact, h2e_v1, h2e_v2, h1e_v, &
  !$omp        dm1, dm2, a12, a13, norm, energy) &
  !$omp private(r, i, p, q, a, b, h1e_ri, nsum, esum, v2_diag_ri, &
  !$omp         v1_rpiq, v2_rpqi, v1_raib, v2_rabi, v2_raqi, &
  !$omp         dm1_qp, dm2_qpab, dm2_qbap, a12_qpab, a13_pqab)
  do r = 0, nvirt - 1
    do i = 0, ncore - 1
      h1e_ri = h1e_v(r * ncore + i)
      nsum = 2.0_c_double * h1e_ri * h1e_ri  ! T9
      esum = 0.0_c_double

      ! Precompute sum_a v2[r,a,a,i] for term 6
      v2_diag_ri = 0.0_c_double
      do a = 0, nact - 1
        v2_diag_ri = v2_diag_ri + h2e_v2(((r * nact + a) * nact + a) * ncore + i)
      end do

      ! Terms 1,2,3,5,6,7,8 (p,q,a,b) + Energy E1-E4
      ! T6,T7,T8 only need p,q (not a,b), computed inside the 4-index loop
      ! T4 needs a,q,b (no p), handled in a separate loop below
      do p = 0, nact - 1
        do q = 0, nact - 1
          dm1_qp = dm1(q * nact + p)
          v1_rpiq = h2e_v1(((r * nact + p) * ncore + i) * nact + q)
          v2_rpqi = h2e_v2(((r * nact + p) * nact + q) * ncore + i)

          ! T7 + T8 + T6 (p,q-dependent, no a,b needed)
          nsum = nsum + (4.0_c_double * v1_rpiq - 2.0_c_double * v2_rpqi) * h1e_ri * dm1_qp
          nsum = nsum + v2_rpqi * v2_diag_ri * dm1_qp  ! T6

          do a = 0, nact - 1
            do b = 0, nact - 1
              v1_raib = h2e_v1(((r * nact + a) * ncore + i) * nact + b)
              v2_rabi = h2e_v2(((r * nact + a) * nact + b) * ncore + i)

              dm2_qpab = dm2(((q * nact + p) * nact + a) * nact + b)
              dm2_qbap = dm2(((q * nact + b) * nact + a) * nact + p)

              ! T1: 2*v1_rpiq*v1_raib*dm2_qpab
              nsum = nsum + 2.0_c_double * v1_rpiq * v1_raib * dm2_qpab
              ! T2: -v1_rpiq*v2_rabi*dm2_qpab
              nsum = nsum - v1_rpiq * v2_rabi * dm2_qpab
              ! T3: -v2_rpqi*v1_raib*dm2_qpab
              nsum = nsum - v2_rpqi * v1_raib * dm2_qpab
              ! T5: -v2_rpqi*v2_rabi*dm2_qbap
              nsum = nsum - v2_rpqi * v2_rabi * dm2_qbap

              ! Energy terms (a12 for dm2-based, a13 for dm1-based)
              a12_qpab = a12(((q * nact + p) * nact + a) * nact + b)
              a13_pqab = a13(((p * nact + q) * nact + a) * nact + b)

              ! E1: 2*v1_rpiq*v1_raib*a12_qpab
              esum = esum + 2.0_c_double * v1_rpiq * v1_raib * a12_qpab
              ! E2: -v1_rpiq*v2_rabi*a12_qpab
              esum = esum - v1_rpiq * v2_rabi * a12_qpab
              ! E3: -v2_rpqi*v1_raib*a12_qpab
              esum = esum - v2_rpqi * v1_raib * a12_qpab
              ! E4: +v2_rpqi*v2_rabi*a13_pqab
              esum = esum + v2_rpqi * v2_rabi * a13_pqab
            end do
          end do
        end do
      end do

      ! Term 4: 2.0 * sum_{a,q,b} v2[r,a,q,i] * v2[r,a,b,i] * dm1[q,b]
      ! (no p dependence, separate loop avoids O(nact^5) work)
      do a = 0, nact - 1
        do q = 0, nact - 1
          v2_raqi = h2e_v2(((r * nact + a) * nact + q) * ncore + i)
          do b = 0, nact - 1
            v2_rabi = h2e_v2(((r * nact + a) * nact + b) * ncore + i)
            nsum = nsum + 2.0_c_double * v2_raqi * v2_rabi * dm1(q * nact + b)
          end do
        end do
      end do

      norm(r * ncore + i) = nsum
      energy(r * ncore + i) = esum
    end do
  end do
  !$omp end parallel do
end subroutine nevpt2_sir

end module nevpt2_koopmans_mod
