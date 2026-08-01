!> @brief Fortran engine for the fixed-CI half of the analytic CASSCF orbital
!> Hessian of pyoqp casscf_hessian.py.
!>
!> The Python assembly walks the non-redundant rotation pairs and, for each,
!> materialises the one-index derivative integrals as a full nbf^4 tensor and
!> contracts them against the 2-RDM:
!>
!>     T^(l) = one-index derivative of g along K^(l)          O(nbf^4) storage
!>     Z^(l) = (D t - t D)^T + 2 sum_bcd G[n,b,c,d] T^(l)[m,b,c,d]   O(nbf^5)
!>     B[:,l] = Z^(l)[p_k,q_k] - Z^(l)[q_k,p_k]
!>
!> so the loop costs npar * O(nbf^5).  With a large basis and a small active
!> space -- the regime the dense-spectrum analytic Hessian is actually good
!> for -- that dominates everything else in the module.
!>
!> Neither the nbf^5 nor the nbf^4 is necessary, because T^(l) is sparse: the
!> generator K^(l) = e_PQ - e_QP applied to each of the four slots leaves only
!> eight slabs of g,
!>
!>     T[a,b,c,d] = d(a,Q)g[P,b,c,d] - d(a,P)g[Q,b,c,d]
!>                + d(b,Q)g[a,P,c,d] - d(b,P)g[a,Q,c,d]
!>                + d(c,Q)g[a,b,P,d] - d(c,P)g[a,b,Q,d]
!>                + d(d,Q)g[a,b,c,P] - d(d,P)g[a,b,c,Q].
!>
!> Contracting slab by slab, the slot-1 term collapses onto a single
!> pair-independent intermediate and the other three become nbf^4 GEMMs:
!>
!>     Y[n,a]      = sum_bcd G[n,b,c,d] g[a,b,c,d]      ONCE, not per pair
!>     A[X,Y][m,n] = sum_cd  g[m,X,c,d] G[n,Y,c,d]
!>     B[X,Y][m,n] = sum_bd  g[m,b,X,d] G[n,b,Y,d]
!>     C[X,Y][m,n] = sum_bc  g[m,b,c,X] G[n,b,c,Y]
!>
!>     C1[m,n] = A[P,Q]-A[Q,P] + B[P,Q]-B[Q,P] + C[P,Q]-C[Q,P]
!>     C1[Q,n] += Y[n,P],   C1[P,n] -= Y[n,Q]
!>
!> so the per-pair cost drops from O(nbf^5) to O(nbf^4) and no nbf^4 temporary
!> is ever formed.  The A contraction needs no gather at all: g[:,X,:,:] is
!> already a valid column-major (nbf^2 x nbf) matrix with leading dimension
!> nbf^3.  B and C fix an interior index, so those two need their operands
!> gathered into contiguous buffers -- nbf^3 per gather, negligible beside the
!> nbf^4 GEMM it feeds.
!>
!> The single contraction above (rather than the four the Python originally
!> evaluated) is justified separately: T inherits the full eight-fold symmetry
!> of the real MO ERIs, so all four coincide.  See casscf_hessian._z_matrix and
!> tests/test_casscf_zmatrix.py.
!>
!> The folded active derivative integrals f^(l), g^(l) that the CI-relaxation
!> half consumes are produced here too, evaluated directly from the eight slabs
!> rather than by slicing a materialised T.
!>
!> All arrays cross the boundary in the caller's C order (last index fastest).
module casscf_hess_bmat_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: casscf_hess_bmat

contains

  !> Fixed-CI Hessian columns B and the folded active derivative integrals.
  !>
  !> @param[in]  nbf    orbitals
  !> @param[in]  ncore  inactive orbitals
  !> @param[in]  nact   active orbitals
  !> @param[in]  npar   non-redundant rotation pairs
  !> @param[in]  pairs  pair list, C-order [npar,2] as (p,q)
  !> @param[in]  dmat   full spin-summed 1-RDM, C-order [nbf,nbf]
  !> @param[in]  rdm2   full spin-summed chemist 2-RDM, C-order [nbf,nbf,nbf,nbf]
  !> @param[in]  h1e    MO core Hamiltonian, C-order [nbf,nbf]
  !> @param[in]  eri    MO ERIs (chemist), C-order [nbf,nbf,nbf,nbf]
  !> @param[out] bmat   fixed-CI Hessian columns, C-order [npar,npar]
  !> @param[out] fder   folded 1e derivative integrals, C-order [npar,nact,nact]
  !> @param[out] gder   folded 2e derivative integrals, C-order
  !>                    [npar,nact,nact,nact,nact]
  subroutine casscf_hess_bmat(nbf, ncore, nact, npar, pairs, dmat, rdm2, &
                              h1e, eri, bmat, fder, gder) &
      bind(C, name="casscf_hess_bmat")
    integer(c_int32_t), value :: nbf, ncore, nact, npar
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    integer(c_int32_t), intent(in) :: pairs(0:*)
    real(dp), intent(in) :: dmat(0:*), rdm2(0:*), h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: bmat(0:*), fder(0:*), gder(0:*)

    ! Default (ILP64) integer kind, matching the linked BLAS.
    integer :: n, nc, na, np, n2, l, k, p, q, i, j, m, t, u, v, w
    integer(i8) :: n3
    real(dp) :: acc, zpq, zqp
    real(dp), allocatable :: ymat(:,:), cmat(:,:), tmat(:,:), w1(:,:), w2(:,:)
    real(dp), allocatable :: dm(:,:), gb(:,:), rb(:,:), gc(:,:), rc(:,:)
    real(dp), allocatable :: gbp(:,:), rbp(:,:), gcp(:,:), rcp(:,:)
    integer :: last_p
    real(dp), allocatable :: zmat(:,:)

    n = int(nbf)
    nc = int(ncore)
    na = int(nact)
    np = int(npar)
    n2 = n * n
    n3 = int(n, i8) * int(n, i8) * int(n, i8)
    if (n <= 0 .or. np <= 0) return

    allocate(ymat(0:n-1, 0:n-1), cmat(0:n-1, 0:n-1), tmat(0:n-1, 0:n-1))
    allocate(w1(0:n-1, 0:n-1), w2(0:n-1, 0:n-1), zmat(0:n-1, 0:n-1))
    allocate(dm(0:n-1, 0:n-1))
    allocate(gb(0:n2-1, 0:n-1), rb(0:n2-1, 0:n-1))
    allocate(gc(0:n2-1, 0:n-1), rc(0:n2-1, 0:n-1))
    ! Separate buffers for the p-indexed gathers so they can be reused across
    ! the run of pairs that share a p (see the loop below).
    allocate(gbp(0:n2-1, 0:n-1), rbp(0:n2-1, 0:n-1))
    allocate(gcp(0:n2-1, 0:n-1), rcp(0:n2-1, 0:n-1))
    last_p = -1

    do i = 0, n - 1
      do j = 0, n - 1
        dm(i, j) = dmat(int(i, i8)*int(n, i8) + int(j, i8))
      end do
    end do

    ! ---- Y[n,a] = sum_bcd G[n,b,c,d] g[a,b,c,d], once for all pairs.
    ! Both buffers are already column-major ((b,c,d), leading index) matrices
    ! with leading dimension nbf^3, so this is one GEMM with no repacking.
    ! ymat(a,n) holds Y[n,a].
    call dgemm('T', 'N', n, n, int(n3), 1.0_dp, eri, int(n3), &
               rdm2, int(n3), 0.0_dp, ymat, n)

    do l = 0, np - 1
      p = int(pairs(2*l))
      q = int(pairs(2*l + 1))

      ! The four p-indexed gathers are invariant across a run of pairs sharing
      ! the same p, and _nonredundant_pairs emits p outermost -- (active,
      ! inactive) then (virtual, inactive) then (virtual, active) -- so those
      ! runs are consecutive.  Re-gathering them per pair was half of this
      ! routine's memory traffic, which is what kept it at 11-16 GFlop/s while
      ! doing 8x fewer flops than the NumPy it replaced.
      if (p /= last_p) then
        call gather_third(eri, p, n, gbp)
        call gather_third(rdm2, p, n, rbp)
        call gather_fourth(eri, p, n, gcp)
        call gather_fourth(rdm2, p, n, rcp)
        last_p = p
      end if

      ! ---- C1 = A[P,Q] - A[Q,P]
      ! g[:,X,:,:] is a column-major (nbf^2 x nbf) matrix, base X*nbf^2 and
      ! leading dimension nbf^3 -- passed by element reference, no gather.
      call dgemm('T', 'N', n, n, n2, 1.0_dp, eri(int(p, i8)*int(n2, i8)), int(n3), &
                 rdm2(int(q, i8)*int(n2, i8)), int(n3), 0.0_dp, cmat, n)
      call dgemm('T', 'N', n, n, n2, -1.0_dp, eri(int(q, i8)*int(n2, i8)), int(n3), &
                 rdm2(int(p, i8)*int(n2, i8)), int(n3), 1.0_dp, cmat, n)

      ! ---- C1 += B[P,Q] - B[Q,P],  B[X,Y][m,n] = sum_bd g[m,b,X,d] G[n,b,Y,d]
      call gather_third(rdm2, q, n, rb)
      call dgemm('T', 'N', n, n, n2, 1.0_dp, gbp, n2, rb, n2, 1.0_dp, cmat, n)
      call gather_third(eri, q, n, gb)
      call dgemm('T', 'N', n, n, n2, -1.0_dp, gb, n2, rbp, n2, 1.0_dp, cmat, n)

      ! ---- C1 += C[P,Q] - C[Q,P],  C[X,Y][m,n] = sum_bc g[m,b,c,X] G[n,b,c,Y]
      call gather_fourth(rdm2, q, n, rc)
      call dgemm('T', 'N', n, n, n2, 1.0_dp, gcp, n2, rc, n2, 1.0_dp, cmat, n)
      call gather_fourth(eri, q, n, gc)
      call dgemm('T', 'N', n, n, n2, -1.0_dp, gc, n2, rcp, n2, 1.0_dp, cmat, n)

      ! ---- slot-1 slabs: the only pair dependence is which row they land on
      do j = 0, n - 1
        cmat(q, j) = cmat(q, j) + ymat(p, j)
        cmat(p, j) = cmat(p, j) - ymat(q, j)
      end do

      ! ---- t = [h, K] for this pair (casscf_hessian._one_index_derivative_h)
      tmat = 0.0_dp
      do i = 0, n - 1
        tmat(i, q) = tmat(i, q) + h1e(int(i, i8)*int(n, i8) + int(p, i8))
        tmat(i, p) = tmat(i, p) - h1e(int(i, i8)*int(n, i8) + int(q, i8))
      end do
      do j = 0, n - 1
        tmat(p, j) = tmat(p, j) - h1e(int(q, i8)*int(n, i8) + int(j, i8))
        tmat(q, j) = tmat(q, j) + h1e(int(p, i8)*int(n, i8) + int(j, i8))
      end do

      ! ---- Z[m,n] = (D t - t D)[n,m] + 2 C1[m,n]
      call dgemm('N', 'N', n, n, n, 1.0_dp, dm, n, tmat, n, 0.0_dp, w1, n)
      call dgemm('N', 'N', n, n, n, 1.0_dp, tmat, n, dm, n, 0.0_dp, w2, n)
      do m = 0, n - 1
        do j = 0, n - 1
          zmat(m, j) = (w1(j, m) - w2(j, m)) + 2.0_dp * cmat(m, j)
        end do
      end do

      ! ---- B[k,l] = Z[p_k,q_k] - Z[q_k,p_k]
      do k = 0, np - 1
        i = int(pairs(2*k))
        j = int(pairs(2*k + 1))
        zpq = zmat(i, j)
        zqp = zmat(j, i)
        bmat(int(k, i8)*int(np, i8) + int(l, i8)) = zpq - zqp
      end do

      ! ---- folded active derivative integrals, straight from the slabs
      do t = 0, na - 1
        do u = 0, na - 1
          do v = 0, na - 1
            do w = 0, na - 1
              gder((((int(l, i8)*int(na, i8) + int(t, i8))*int(na, i8) &
                      + int(u, i8))*int(na, i8) + int(v, i8))*int(na, i8) &
                      + int(w, i8)) = &
                  tval(eri, n, p, q, nc + t, nc + u, nc + v, nc + w)
            end do
          end do
        end do
      end do
      do t = 0, na - 1
        do u = 0, na - 1
          acc = tmat(nc + t, nc + u)
          do i = 0, nc - 1
            acc = acc + 2.0_dp * tval(eri, n, p, q, nc + t, nc + u, i, i) &
                      - tval(eri, n, p, q, nc + t, i, i, nc + u)
          end do
          fder((int(l, i8)*int(na, i8) + int(t, i8))*int(na, i8) + int(u, i8)) = acc
        end do
      end do
    end do

    deallocate(ymat, cmat, tmat, w1, w2, zmat, dm, gb, rb, gc, rc)
    deallocate(gbp, rbp, gcp, rcp)

  contains

    !> out((b,d), m) = src[m,b,X,d]: fixes the third index of a C-order
    !> [n,n,n,n] buffer and makes the contracted pair contiguous.
    subroutine gather_third(src, x, nn, out)
      real(dp), intent(in) :: src(0:*)
      integer, intent(in) :: x, nn
      real(dp), intent(inout) :: out(0:,0:)
      integer :: mm, bb, dd
      integer(i8) :: base
      do mm = 0, nn - 1
        do bb = 0, nn - 1
          base = int(mm, i8)*int(nn, i8)**3 + int(bb, i8)*int(nn, i8)**2 &
               + int(x, i8)*int(nn, i8)
          do dd = 0, nn - 1
            out(bb*nn + dd, mm) = src(base + int(dd, i8))
          end do
        end do
      end do
    end subroutine gather_third

    !> out((b,c), m) = src[m,b,c,X]: fixes the fourth index.
    subroutine gather_fourth(src, x, nn, out)
      real(dp), intent(in) :: src(0:*)
      integer, intent(in) :: x, nn
      real(dp), intent(inout) :: out(0:,0:)
      integer :: mm, bb, cc
      integer(i8) :: base
      do mm = 0, nn - 1
        do bb = 0, nn - 1
          base = int(mm, i8)*int(nn, i8)**3 + int(bb, i8)*int(nn, i8)**2 &
               + int(x, i8)
          do cc = 0, nn - 1
            out(bb*nn + cc, mm) = src(base + int(cc, i8)*int(nn, i8))
          end do
        end do
      end do
    end subroutine gather_fourth

    !> One element of the one-index derivative tensor T[a,b,c,d] for the pair
    !> (P,Q), evaluated from the eight slabs instead of a stored tensor.
    pure function tval(g, nn, pp, qq, a, b, c, d) result(val)
      real(dp), intent(in) :: g(0:*)
      integer, intent(in) :: nn, pp, qq, a, b, c, d
      real(dp) :: val
      val = 0.0_dp
      if (a == qq) val = val + gel(g, nn, pp, b, c, d)
      if (a == pp) val = val - gel(g, nn, qq, b, c, d)
      if (b == qq) val = val + gel(g, nn, a, pp, c, d)
      if (b == pp) val = val - gel(g, nn, a, qq, c, d)
      if (c == qq) val = val + gel(g, nn, a, b, pp, d)
      if (c == pp) val = val - gel(g, nn, a, b, qq, d)
      if (d == qq) val = val + gel(g, nn, a, b, c, pp)
      if (d == pp) val = val - gel(g, nn, a, b, c, qq)
    end function tval

    !> C-order element g[a,b,c,d] of an [n,n,n,n] buffer.
    pure function gel(g, nn, a, b, c, d) result(val)
      real(dp), intent(in) :: g(0:*)
      integer, intent(in) :: nn, a, b, c, d
      real(dp) :: val
      val = g(((int(a, i8)*int(nn, i8) + int(b, i8))*int(nn, i8) &
               + int(c, i8))*int(nn, i8) + int(d, i8))
    end function gel

  end subroutine casscf_hess_bmat

end module casscf_hess_bmat_mod
