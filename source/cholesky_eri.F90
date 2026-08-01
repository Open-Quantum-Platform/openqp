!> @file cholesky_eri.F90
!>
!> @brief Pivoted Cholesky decomposition of the two-electron integrals.
!>
!> The two-electron integral matrix, indexed by compound pairs,
!>
!>   M(P,Q) = (mu nu | lambda sigma),   P = pair(mu,nu), Q = pair(lambda,sigma)
!>
!> is symmetric positive semi-definite, so it factorises as
!>
!>   M(P,Q) = sum_J L(P,J) L(Q,J)
!>
!> and its numerical rank is far below its dimension: the vectors needed to
!> reach a given accuracy grow roughly linearly in the basis size rather than
!> quadratically.  That turns an N^4 object into an N^3 one, with the
!> truncation controlled by a tolerance rather than by an auxiliary basis --
!> the same factorisation density fitting produces, obtained from the
!> integrals themselves and needing no fitting set.
!>
!> The algorithm is the standard pivoted Cholesky, matching the formulation in
!> Open-Quantum-Platform/AFQMC (`cholesky.F90`): repeatedly take the largest
!> remaining diagonal as the pivot, form its column of the residual, scale by
!> the square root of the pivot, and deflate.
!>
!> It is written differently here for one reason.  That version materialises
!> the residual as a dense (nbf^2)^2 array and updates it in place, which
!> needs the whole N^4 matrix resident -- exactly the cost the factorisation
!> exists to avoid.  Storing only the accumulated vectors and the running
!> diagonal gives the same numbers with O(npair * nchol) memory, so the
!> decomposition itself never becomes the bottleneck it is meant to remove.
module cholesky_eri

  use precision, only: dp

  implicit none

  private
  public :: cholesky_eri_decompose
  public :: cholesky_eri_max_vectors
  public :: cholesky_packed_index
  public :: cholesky_transform_vv

contains

!###############################################################################

!> Offset of the pair-pair element (P,Q) in a packed triangular store.
!> int64 throughout -- the packed length passes 2^31 near nbf = 400.
  pure integer(8) function cholesky_packed_index(p, q) result(idx)
    integer, intent(in) :: p, q
    integer(8) :: hi, lo
    hi = int(max(p,q), 8)
    lo = int(min(p,q), 8)
    idx = hi*(hi-1)/2 + lo
  end function cholesky_packed_index

!###############################################################################

!> A safe upper bound on the number of vectors to allocate for.
!>
!> The rank is bounded by the matrix dimension, but in practice it lands near
!> a small multiple of nbf.  Callers size their buffer with this and the
!> decomposition stops early on tolerance; @p factor is how many multiples of
!> nbf to allow before the bound is simply npair.
  pure integer function cholesky_eri_max_vectors(nbf, factor) result(n)
    integer, intent(in) :: nbf
    integer, intent(in) :: factor
    integer :: npair
    npair = nbf*(nbf+1)/2
    n = min(npair, max(1, factor)*nbf)
  end function cholesky_eri_max_vectors

!###############################################################################

!> @brief Factorise the packed AO integrals into Cholesky vectors.
!>
!> @param[in]  nbf      number of basis functions
!> @param[in]  g        packed AO integrals, canonical eighth, as built by
!>                      cc_ao2mo's collecting consumer
!> @param[in]  tol      stop once the largest remaining diagonal falls below
!>                      this; the reconstruction error is bounded by it
!> @param[in]  maxchol  capacity of @p lvec
!> @param[out] lvec     the vectors, lvec(P,J), only 1:nchol filled
!> @param[out] nchol    how many were produced
!> @param[out] err      largest remaining diagonal on exit -- the achieved
!>                      accuracy, which the caller should report
!> @param[out] truncated .true. if it ran out of capacity before reaching tol
!>
!> Only the running diagonal and the accepted vectors are kept, so the working
!> set is O(npair * nchol) rather than the O(npair^2) a residual matrix would
!> need.
  subroutine cholesky_eri_decompose(nbf, g, tol, maxchol, lvec, nchol, err, &
                                    truncated)

    integer, intent(in) :: nbf
    real(dp), intent(in) :: g(*)
    real(dp), intent(in) :: tol
    integer, intent(in) :: maxchol
    real(dp), intent(out) :: lvec(:,:)
    integer, intent(out) :: nchol
    real(dp), intent(out) :: err
    logical, intent(out) :: truncated

    real(dp), allocatable :: diag(:), col(:)
    integer :: npair, p, j, pivot
    real(dp) :: vmax, scale_

    npair = nbf*(nbf+1)/2
    nchol = 0
    truncated = .false.

    allocate(diag(npair), col(npair))

    ! The diagonal is the Schwarz quantity (mu nu | mu nu), positive by
    ! construction; a negative entry would mean the integrals are not the
    ! matrix this factorisation assumes.
    do p = 1, npair
      diag(p) = g(cholesky_packed_index(p, p))
    end do

    err = maxval(diag)

    do
      pivot = maxloc(diag, dim=1)
      vmax = diag(pivot)
      err = vmax

      if (vmax < tol) exit
      if (nchol >= maxchol) then
        truncated = .true.
        exit
      end if

      nchol = nchol + 1

      ! Column of the residual at the pivot: the original column minus what
      ! the vectors accepted so far already account for.
      !$omp parallel do default(shared) private(p) schedule(static)
      do p = 1, npair
        col(p) = g(cholesky_packed_index(p, pivot))
      end do
      !$omp end parallel do

      do j = 1, nchol-1
        call daxpy(npair, -lvec(pivot,j), lvec(1,j), 1, col, 1)
      end do

      scale_ = 1.0_dp/sqrt(vmax)
      !$omp parallel do default(shared) private(p) schedule(static)
      do p = 1, npair
        lvec(p,nchol) = col(p)*scale_
        ! Deflate.  Clamping at zero keeps rounding on a nearly exhausted
        ! diagonal from turning into a spurious negative pivot.
        diag(p) = max(diag(p) - lvec(p,nchol)*lvec(p,nchol), 0.0_dp)
      end do
      !$omp end parallel do

      diag(pivot) = 0.0_dp
    end do

    deallocate(diag, col)

  end subroutine cholesky_eri_decompose

!###############################################################################

!> @brief Transform the AO Cholesky vectors into the virtual-virtual MO block.
!>
!>   B(a,c,J) = sum_(mu,nu) C(mu,a) C(nu,c) L(pair(mu,nu), J)
!>
!> so that (ac|bd) = sum_J B(a,c,J) B(b,d,J).  This is the whole point of the
!> factorisation: the v^4 ladder integrals are never stored, only this v^2 x
!> nchol object, and the blocks the ladder wants are rebuilt from it.
!>
!> @param[in]  nbf    basis functions
!> @param[in]  nv     virtual orbitals
!> @param[in]  cv     virtual MO coefficients, (nbf, nv)
!> @param[in]  lvec   AO Cholesky vectors, (npair, nchol)
!> @param[in]  nchol  number of vectors
!> @param[out] bvv    (nv*nv, nchol), indexed (a + (c-1)*nv, J)
  subroutine cholesky_transform_vv(nbf, nv, cv, lvec, nchol, bvv)

    integer, intent(in) :: nbf, nv, nchol
    real(dp), intent(in) :: cv(nbf, nv)
    real(dp), intent(in) :: lvec(:,:)
    real(dp), intent(out) :: bvv(nv*nv, nchol)

    real(dp), allocatable :: lao(:,:), scr(:,:)
    integer :: j, mu, nu, ip

    !$omp parallel default(shared) private(j, mu, nu, ip, lao, scr)
    allocate(lao(nbf,nbf), scr(nbf,nv))
    !$omp do schedule(static)
    do j = 1, nchol
      ! Unpack the vector into a full symmetric AO matrix.
      do mu = 1, nbf
        do nu = 1, mu
          ip = mu*(mu-1)/2 + nu
          lao(mu,nu) = lvec(ip, j)
          lao(nu,mu) = lao(mu,nu)
        end do
      end do
      call dgemm('n','n', nbf, nv, nbf, 1.0_dp, lao, nbf, cv, nbf, 0.0_dp, scr, nbf)
      call dgemm('t','n', nv, nv, nbf, 1.0_dp, cv, nbf, scr, nbf, 0.0_dp, &
                 bvv(1,j), nv)
    end do
    !$omp end do
    deallocate(lao, scr)
    !$omp end parallel

  end subroutine cholesky_transform_vv

end module cholesky_eri
