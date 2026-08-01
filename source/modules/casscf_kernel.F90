!> @brief Fortran engine for the CASSCF generalized Fock matrix and the
!> orbital-rotation gradient of pyoqp casscf.py.
!>
!> Every macroiteration of the orbital optimizer evaluates the gradient once
!> per finite-difference displacement (2*npar of them for the FD Hessian), so
!> this build sits on the hottest path of CASSCF.  One bind(C) entry point
!> replaces the Python assembly (the Python implementation remains as the
!> numerical pin):
!>
!>  * casscf_gfock_grad -- weighted generalized Fock F and the gradient g
!>
!> The Python reference forms the full spin-summed 2-RDM
!>
!>     G[p,q,r,s] = D[p,q] D[r,s] - 1/2 D[p,s] D[r,q]     (active block
!>                                                         replaced by Gamma)
!>
!> as an nbf^4 array and contracts it against the ERIs, which costs O(nbf^5)
!> flops and O(nbf^4) memory.  Neither is necessary.  Writing G as its
!> separable part plus a correction that is nonzero only when all four indices
!> are active,
!>
!>     G = G_sep(D) + dG,   dG[t,u,v,w] = Gamma[t,u,v,w] - G_sep(D)[t,u,v,w],
!>
!> the separable part contracts through the ordinary Coulomb/exchange matrices
!>
!>     J[n,q] = sum_rs (nq|rs) D[r,s],   K[n,s] = sum_qr (nq|rs) D[r,q],
!>     F_sep  = D (h + J - 1/2 K),
!>
!> which is O(nbf^4) and needs only nbf^2 storage, while dG contributes a
!> single small DGEMM of size (nbf x nact^3)(nact^3 x nact).  The result is
!> the same matrix to round-off; the nbf^5 term and the nbf^4 allocation are
!> simply gone.
!>
!> For SA-CASSCF the caller hands over one (gamma, Gamma) pair per root and the
!> weighted sum sum_I w_I F_I is accumulated here.  G_sep is quadratic in D, so
!> the per-root loop is not obviously redundant -- but the term quadratic in
!> the root-dependent (active) part of D has all four indices active, which is
!> exactly the block dG overwrites, so it never reaches F.  On the surviving
!> blocks the build is affine in (gamma, Gamma) and averaging the RDMs first
!> gives the same F *provided the weights sum to 1*: the core-core term picks
!> up sum_I w_I one way and (sum_I w_I)^2 the other.  CASSCF normalizes its
!> weights, but the per-root loop is kept so the engine does not silently
!> depend on that; nroot is small and each pass is only O(nbf^4).
!>
!> All arrays cross the boundary in the caller's C order (last index fastest).
!> Because D, J, K and h are symmetric, their C-order buffers coincide with
!> the Fortran column-major views, and the column-major product written into
!> `fock` is already the C-order F the caller expects.
module casscf_kernel_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: casscf_gfock_grad, casscf_effective_fock

contains

  !> Closed+active mean-field matrix  W = h + J - K/2  for a spin-summed 1-RDM.
  !>
  !>     J[j,q] = sum_rs (j q|rs) D[r,s],   K[j,s] = sum_qr (j q|rs) D[r,q]
  !>
  !> Shared by the generalized-Fock build below (where W is the inner factor of
  !> the separable part) and by the canonicalization Fock `casscf_effective_fock`
  !> -- they are the same matrix, so this is written once.
  !>
  !> `dmat` is the Fortran (n,n) density, `wmat` is returned as wmat(i,j) = W[i,j];
  !> `jmat`/`kmat` are caller-supplied scratch of the same shape.  `h1e` and
  !> `eri` are the caller's C-order buffers (D, J, K and h are symmetric, so
  !> their C-order buffers coincide with the column-major views).
  subroutine build_jkw(n, dmat, h1e, eri, jmat, kmat, wmat)
    integer, intent(in) :: n
    real(dp), intent(in) :: dmat(0:n-1, 0:n-1)
    real(dp), intent(in) :: h1e(0:*), eri(0:*)
    real(dp), intent(out) :: jmat(0:n-1, 0:n-1), kmat(0:n-1, 0:n-1)
    real(dp), intent(out) :: wmat(0:n-1, 0:n-1)

    integer :: i, j, n2
    integer(i8) :: off, n3

    n2 = n * n
    n3 = int(n, i8) * int(n, i8) * int(n, i8)

    ! The ERI buffer viewed column-major is M((r,s),(j,q)) with the pair (r,s)
    ! contiguous, so J is one DGEMV against the flattened D.  Vector arguments
    ! are passed by element reference so no rank-mismatch temporary is made.
    call dgemv('T', n2, n2, 1.0_dp, eri, n2, dmat(0, 0), 1, 0.0_dp, jmat(0, 0), 1)

    ! For fixed j the block is contiguous as (s,(q,r)) with s fastest, so each
    ! column of K is a DGEMV against the same flattened D.  The ERI block is
    ! passed by element reference so the following columns stay contiguous
    ! through sequence association.
    do j = 0, n - 1
      off = int(j, i8) * n3
      call dgemv('N', n, n2, 1.0_dp, eri(off), n, dmat(0, 0), 1, &
                 0.0_dp, kmat(0, j), 1)
    end do

    ! The DGEMV outputs land transposed relative to W (jmat(q,j) = J[j,q],
    ! kmat(s,j) = K[j,s]), so they are read back transposed here.
    do i = 0, n - 1
      do j = 0, n - 1
        wmat(i, j) = h1e(int(i, i8)*int(n, i8) + int(j, i8)) &
                   + jmat(j, i) - 0.5_dp * kmat(j, i)
      end do
    end do
  end subroutine build_jkw

  !> Closed+active mean-field Fock used to canonicalize the CASSCF orbitals.
  !>
  !> Replaces casscf.py `_effective_fock`, whose two O(nbf^5)-shaped einsums
  !>
  !>     J = einsum("rs,pqrs->pq", D, eri),  K = einsum("rs,prsq->pq", D, eri)
  !>
  !> build the very same J and K the generalized-Fock engine already forms (the
  !> K contractions differ by a transpose of D, which is symmetric here), so
  !> this entry is a thin bind(C) wrapper over the shared `build_jkw`.
  !>
  !> @param[in]  nbf   number of orbitals
  !> @param[in]  dmat  spin-summed 1-RDM, C-order [nbf,nbf]
  !> @param[in]  h1e   MO core Hamiltonian, C-order [nbf,nbf]
  !> @param[in]  eri   MO ERIs (chemist (pq|rs)), C-order [nbf,nbf,nbf,nbf]
  !> @param[out] fout  h + J - K/2, C-order [nbf,nbf]
  subroutine casscf_effective_fock(nbf, dmat, h1e, eri, fout) &
      bind(C, name="casscf_effective_fock")
    integer(c_int32_t), value :: nbf
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: dmat(0:*), h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: fout(0:*)

    real(dp), allocatable :: dwork(:,:), jmat(:,:), kmat(:,:), wmat(:,:)
    integer :: n, i, j

    n = int(nbf)
    if (n <= 0) return

    allocate(dwork(0:n-1, 0:n-1), jmat(0:n-1, 0:n-1), kmat(0:n-1, 0:n-1))
    allocate(wmat(0:n-1, 0:n-1))

    ! D is symmetric, so its C-order buffer is already the column-major view;
    ! the gather is written out anyway so the engine does not depend on it.
    do i = 0, n - 1
      do j = 0, n - 1
        dwork(i, j) = dmat(int(i, i8)*int(n, i8) + int(j, i8))
      end do
    end do

    call build_jkw(n, dwork, h1e, eri, jmat, kmat, wmat)

    do i = 0, n - 1
      do j = 0, n - 1
        fout(int(i, i8)*int(n, i8) + int(j, i8)) = wmat(i, j)
      end do
    end do

    deallocate(dwork, jmat, kmat, wmat)
  end subroutine casscf_effective_fock

  !> Weighted generalized Fock F[p,q] and orbital gradient g.
  !>
  !> @param[in]  nbf      number of orbitals
  !> @param[in]  ncore    inactive (doubly occupied) orbitals
  !> @param[in]  nact     active orbitals
  !> @param[in]  nroot    number of state-average roots (1 for state-specific)
  !> @param[in]  weights  state-average weights, [nroot]
  !> @param[in]  gamma    active 1-RDMs, C-order [nroot,nact,nact]
  !> @param[in]  gamma2   active 2-RDMs, C-order [nroot,nact,nact,nact,nact]
  !> @param[in]  h1e      MO core Hamiltonian, C-order [nbf,nbf]
  !> @param[in]  eri      MO ERIs (chemist (pq|rs)), C-order [nbf,nbf,nbf,nbf]
  !> @param[out] fock     generalized Fock, C-order [nbf,nbf]
  !> @param[out] grad     orbital gradient over the non-redundant pairs, in
  !>                      the ordering of casscf.py `_nonredundant_pairs`:
  !>                      (active,inactive), (virtual,inactive), (virtual,active)
  subroutine casscf_gfock_grad(nbf, ncore, nact, nroot, weights, gamma, gamma2, &
                               h1e, eri, fock, grad) &
      bind(C, name="casscf_gfock_grad")
    integer(c_int32_t), value :: nbf, ncore, nact, nroot
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: weights(0:*), gamma(0:*), gamma2(0:*)
    real(dp), intent(in) :: h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: fock(0:*), grad(0:*)

    ! Local copies in the default (ILP64) integer kind the linked BLAS uses.
    integer :: n, nc, na, nr, na3, n2
    integer :: i, j, p, t, u, v, m, root, k
    integer(i8) :: off, n3
    real(dp) :: w, dsep
    real(dp), allocatable :: dmat(:,:), jmat(:,:), kmat(:,:), wmat(:,:)
    real(dp), allocatable :: fwork(:,:), eact(:,:), dgmat(:,:)

    n = int(nbf)
    nc = int(ncore)
    na = int(nact)
    nr = int(nroot)
    na3 = na * na * na
    n2 = n * n
    n3 = int(n, i8) * int(n, i8) * int(n, i8)

    do i = 0, n2 - 1
      fock(i) = 0.0_dp
    end do

    allocate(dmat(0:n-1, 0:n-1), jmat(0:n-1, 0:n-1), kmat(0:n-1, 0:n-1))
    allocate(wmat(0:n-1, 0:n-1), fwork(0:n-1, 0:n-1))
    if (na > 0) allocate(eact(0:n-1, 0:na3-1), dgmat(0:na3-1, 0:na-1))

    ! ERI slice over the active triple (t,u,v): root independent, gathered
    ! once.  eact(j, (t,u,v)) = (j t | u v).  The j loop is outermost so each
    ! pass reads one contiguous nbf^3 block rather than striding through all
    ! of them per element.
    if (na > 0) then
      do j = 0, n - 1
        off = int(j, i8) * n3
        do t = 0, na - 1
          do u = 0, na - 1
            do v = 0, na - 1
              eact(j, (t*na + u)*na + v) = &
                  eri(off + int(nc + t, i8)*int(n2, i8) &
                          + int(nc + u, i8)*int(n, i8) + int(nc + v, i8))
            end do
          end do
        end do
      end do
    end if

    do root = 0, nr - 1
      w = weights(root)

      ! ---- D for this root: 2 on the inactive diagonal, gamma on the active block
      dmat = 0.0_dp
      do i = 0, nc - 1
        dmat(i, i) = 2.0_dp
      end do
      do t = 0, na - 1
        do u = 0, na - 1
          dmat(nc + t, nc + u) = gamma(int(root, i8)*int(na, i8)*int(na, i8) &
                                       + int(t, i8)*int(na, i8) + int(u, i8))
        end do
      end do

      ! ---- W = h + J - K/2 (the same mean-field matrix casscf_effective_fock
      !      returns; the shared builder keeps the two from drifting apart).
      !      W is in fact symmetric for a physically valid ERI, but nothing
      !      below relies on that: wmat(i,j) is W[i,j] as written.
      call build_jkw(n, dmat, h1e, eri, jmat, kmat, wmat)

      ! ---- F_sep[m,j] = sum_q D[m,q] W[j,q]
      ! Computed as W D so that the column-major result (j,m) lands exactly on
      ! the caller's C-order element F[m,j] with no transpose.
      call dgemm('N', 'N', n, n, n, 1.0_dp, wmat, n, dmat, n, 0.0_dp, fwork, n)

      ! ---- all-active correction  F[m,j] += sum_tuv dG[m,t,u,v] (j t|u v)
      if (na > 0) then
        do m = 0, na - 1
          do t = 0, na - 1
            do u = 0, na - 1
              do v = 0, na - 1
                dsep = dmat(nc + m, nc + t) * dmat(nc + u, nc + v) &
                     - 0.5_dp * dmat(nc + m, nc + v) * dmat(nc + u, nc + t)
                dgmat((t*na + u)*na + v, m) = &
                    gamma2((((int(root, i8)*int(na, i8) + int(m, i8))*int(na, i8) &
                            + int(t, i8))*int(na, i8) + int(u, i8))*int(na, i8) &
                            + int(v, i8)) - dsep
              end do
            end do
          end do
        end do
        call dgemm('N', 'N', n, na, na3, 1.0_dp, eact, n, dgmat, na3, &
                   1.0_dp, fwork(0, nc), n)
      end if

      ! ---- accumulate w_root * F_root.  fwork's column-major (j,m) buffer is
      !      already the C-order F[m,j], so this is a flat weighted add.
      do m = 0, n - 1
        do j = 0, n - 1
          fock(int(m, i8)*int(n, i8) + int(j, i8)) = &
              fock(int(m, i8)*int(n, i8) + int(j, i8)) + w * fwork(j, m)
        end do
      end do
    end do

    ! ---- gradient over the non-redundant pairs, in casscf.py's ordering:
    !      g = 2 (F[q,p] - F[p,q]) for each pair (p,q).
    k = 0
    do t = nc, nc + na - 1
      do i = 0, nc - 1
        grad(k) = 2.0_dp * (fock(int(i, i8)*int(n, i8) + int(t, i8)) &
                          - fock(int(t, i8)*int(n, i8) + int(i, i8)))
        k = k + 1
      end do
    end do
    do p = nc + na, n - 1
      do i = 0, nc - 1
        grad(k) = 2.0_dp * (fock(int(i, i8)*int(n, i8) + int(p, i8)) &
                          - fock(int(p, i8)*int(n, i8) + int(i, i8)))
        k = k + 1
      end do
    end do
    do p = nc + na, n - 1
      do t = nc, nc + na - 1
        grad(k) = 2.0_dp * (fock(int(t, i8)*int(n, i8) + int(p, i8)) &
                          - fock(int(p, i8)*int(n, i8) + int(t, i8)))
        k = k + 1
      end do
    end do

    deallocate(dmat, jmat, kmat, wmat, fwork)
    if (allocated(eact)) deallocate(eact, dgmat)
  end subroutine casscf_gfock_grad

end module casscf_kernel_mod
