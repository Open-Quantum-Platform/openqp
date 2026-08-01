!> @brief Fortran engine for the dominant Koopmans intermediates of pyoqp
!> nevpt2_sc.py (strongly contracted NEVPT2, Dyall zeroth order).
!>
!> pyoqp is a driver; liboqp computes.  Measured on synthetic operands of the
!> shape the module actually builds (best-of-5, isolated microbenchmark, so no
!> shared-machine wall-clock), the NumPy `optimize=True` einsum path spends
!>
!>     nact:              4        6        8
!>     f3ca/f3ac      0.23 ms   5.4 ms   41.1 ms     (18% of the module at 8)
!>     a16            0.76 ms   6.7 ms   42.7 ms     (19%)
!>     a22            0.93 ms   6.1 ms   35.3 ms     (16%)
!>     everything else                   ~104 ms
!>
!> so these three carry ~53% of the intermediate cost and grow fastest.  They are
!> ported here; the remaining intermediates stay in NumPy for now and every one of
!> these keeps its Python implementation as the numerical pin and as the fallback
!> when the symbol is absent.
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
module nevpt2_koopmans_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: nevpt2_f3ca_f3ac, nevpt2_a16

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
  !> Nine of the thirteen terms contract dm3 only over its three FASTEST axes
  !> while (r,p,q) ride along as spectators, so the whole intermediate is
  !> evaluated block-wise: for each (p,q,r) the nact^3 sub-block dm3[r,p,q,:,:,:]
  !> is contiguous and so is the output block a16[p,q,r,:,:,:].  Two pairs of
  !> terms also collapse before any arithmetic happens --
  !>
  !>     -sum_i h1e[i,b] G[i,a,c] + sum_ij h2e[j,b,i,j] G[i,a,c]
  !>   = sum_i ( sum_j h2e[j,b,i,j] - h1e[i,b] ) G[i,a,c]
  !>
  !> and likewise for the pair carrying `c` -- so the thirteen einsum calls of
  !> the Python become nine block contractions plus the f3 fold-ins.
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

    integer :: n, n2, n3, p, q, r, a, b, c, i, j, k, mm
    integer(i8) :: n_i, n2_i, n3_i, gbase, obase, fbase
    real(dp) :: acc, hv
    real(dp), allocatable :: u1(:,:), u2(:,:), gb(:), ob(:)
    real(dp), allocatable :: fdm2(:,:,:,:), db(:)

    n = int(nact)
    if (n <= 0) return
    n2 = n * n
    n3 = n2 * n
    n_i = int(n, i8)
    n2_i = n_i * n_i
    n3_i = n2_i * n_i

    ! u1[b,i] = sum_j h2e[j,b,i,j] - h1e[i,b]   (terms 1 + 11 folded)
    ! u2[c,i] = sum_j h2e[j,c,i,j] - h1e[c,i]   (terms 3 + 13 folded)
    allocate(u1(0:n - 1, 0:n - 1), u2(0:n - 1, 0:n - 1))
    do b = 0, n - 1
      do i = 0, n - 1
        acc = 0.0_dp
        do j = 0, n - 1
          acc = acc + h2e(((int(j, i8) * n_i + int(b, i8)) * n_i &
                           + int(i, i8)) * n_i + int(j, i8))
        end do
        u1(b, i) = acc - h1e(int(i, i8) * n_i + int(b, i8))
        u2(b, i) = acc - h1e(int(b, i8) * n_i + int(i, i8))
      end do
    end do

    ! fdm2[p,r,a,b] = sum_kij h2e[k,b,i,j] dm3[r,p,a,j,k,i]
    allocate(fdm2(0:n - 1, 0:n - 1, 0:n - 1, 0:n - 1))
    allocate(db(0:n3 - 1))
    do p = 0, n - 1
      do r = 0, n - 1
        do a = 0, n - 1
          gbase = ((int(r, i8) * n_i + int(p, i8)) * n_i + int(a, i8)) * n3_i
          do mm = 0, n3 - 1
            db(mm) = dm3(gbase + int(mm, i8))
          end do
          do b = 0, n - 1
            acc = 0.0_dp
            do j = 0, n - 1
              do k = 0, n - 1
                do i = 0, n - 1
                  acc = acc + h2e(((int(k, i8) * n_i + int(b, i8)) * n_i &
                                   + int(i, i8)) * n_i + int(j, i8)) &
                            * db(j * n2 + k * n + i)
                end do
              end do
            end do
            fdm2(p, r, a, b) = acc
          end do
        end do
      end do
    end do
    deallocate(db)

    ! ---- block loop over the spectator triple (p,q,r) ----------------------
    allocate(gb(0:n3 - 1), ob(0:n3 - 1))
    do p = 0, n - 1
      do q = 0, n - 1
        do r = 0, n - 1
          gbase = ((int(r, i8) * n_i + int(p, i8)) * n_i + int(q, i8)) * n3_i
          obase = ((int(p, i8) * n_i + int(q, i8)) * n_i + int(r, i8)) * n3_i
          do mm = 0, n3 - 1
            gb(mm) = dm3(gbase + int(mm, i8))
          end do
          ! terms 9 and 10: + f3ac[r,p,q,b,a,c] - f3ca[r,p,q,b,a,c]
          do a = 0, n - 1
            do b = 0, n - 1
              do c = 0, n - 1
                fbase = int(b, i8) * n2_i + int(a, i8) * n_i + int(c, i8)
                ob(a * n2 + b * n + c) = f3ac(gbase + fbase) - f3ca(gbase + fbase)
              end do
            end do
          end do

          do a = 0, n - 1
            do b = 0, n - 1
              do c = 0, n - 1
                acc = ob(a * n2 + b * n + c)
                do i = 0, n - 1
                  ! (1+11) sum_i u1[b,i] G[i,a,c]
                  acc = acc + u1(b, i) * gb(i * n2 + a * n + c)
                  ! (2) + sum_i h1e[i,a] G[b,i,c]
                  acc = acc + h1e(int(i, i8) * n_i + int(a, i8)) &
                            * gb(b * n2 + i * n + c)
                  ! (3+13) sum_i u2[c,i] G[b,a,i]
                  acc = acc + u2(c, i) * gb(b * n2 + a * n + i)
                end do
                do k = 0, n - 1
                  do i = 0, n - 1
                    ! (5) - sum_ki h2e[k,b,i,a] G[c,k,i]
                    acc = acc - h2e(((int(k, i8) * n_i + int(b, i8)) * n_i &
                                     + int(i, i8)) * n_i + int(a, i8)) &
                              * gb(c * n2 + k * n + i)
                    ! (6) - sum_kj h2e[k,b,a,j] G[j,k,c]   (i plays j)
                    acc = acc - h2e(((int(k, i8) * n_i + int(b, i8)) * n_i &
                                     + int(a, i8)) * n_i + int(i, i8)) &
                              * gb(i * n2 + k * n + c)
                    ! (7) + sum_ij h2e[c,b,i,j] G[j,a,i]   (k plays j)
                    acc = acc + h2e(((int(c, i8) * n_i + int(b, i8)) * n_i &
                                     + int(i, i8)) * n_i + int(k, i8)) &
                              * gb(k * n2 + a * n + i)
                    ! (12) - sum_jk h2e[c,j,k,a] G[b,j,k]  (k plays j, i plays k)
                    acc = acc - h2e(((int(c, i8) * n_i + int(k, i8)) * n_i &
                                     + int(i, i8)) * n_i + int(a, i8)) &
                              * gb(b * n2 + k * n + i)
                  end do
                end do
                ob(a * n2 + b * n + c) = acc
              end do
            end do
          end do

          ! term 8: a16[p,q,r,a,b,c] += fdm2[p,r,a,b] when q == c
          do a = 0, n - 1
            do b = 0, n - 1
              ob(a * n2 + b * n + q) = ob(a * n2 + b * n + q) + fdm2(p, r, a, b)
            end do
          end do

          ! term 4: a16[p,q,r,a,b,c] -= f3ca[r,p,a,c,q,b]
          do a = 0, n - 1
            do b = 0, n - 1
              do c = 0, n - 1
                fbase = ((((int(r, i8) * n_i + int(p, i8)) * n_i + int(a, i8)) &
                          * n_i + int(c, i8)) * n_i + int(q, i8)) * n_i + int(b, i8)
                ob(a * n2 + b * n + c) = ob(a * n2 + b * n + c) - f3ca(fbase)
              end do
            end do
          end do

          do mm = 0, n3 - 1
            a16(obase + int(mm, i8)) = ob(mm)
          end do
        end do
      end do
    end do

    deallocate(u1, u2, gb, ob, fdm2)
  end subroutine nevpt2_a16

end module nevpt2_koopmans_mod
