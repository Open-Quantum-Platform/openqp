!> @brief Fortran engine for the CASSCF orbital rotation C <- C exp(K) of
!> pyoqp casscf.py.
!>
!> Every trial step, every backtracking bisection and every one of the 2*npar
!> finite-difference displacements that build the FD orbital Hessian applies
!> one orbital rotation, so this is the single most frequently called piece of
!> the CASSCF optimizer -- on the default `hessian = fd` path it runs
!> 2*npar + O(1) times per macroiteration.  The Python it replaces is
!>
!>     K = _kappa_matrix(vec, pairs, nbf)      # antisymmetric, npar Python loop
!>     C_new = C @ scipy.linalg.expm(K)
!>
!> and both halves move here: `casscf_orbital_rotate` builds K from the
!> non-redundant pair list, exponentiates it, and applies it to the orbitals in
!> one call.  The Python assembly remains as the numerical pin
!> (tests/test_casscf_orbrot.py).
!>
!> The exponential is scaling-and-squaring with the degree-13 diagonal Pade
!> approximant (Higham, SIAM J. Matrix Anal. Appl. 26 (2005) 1179), the same
!> algorithm scipy.linalg.expm uses for its top-order branch:
!>
!>     s = max(0, ceil(log2(||K||_1 / theta13)))    theta13 = 5.371920351148152
!>     A = K / 2^s
!>     A2 = A A,  A4 = A2 A2,  A6 = A2 A4
!>     U = A [ A6 (b13 A6 + b11 A4 + b9 A2) + b7 A6 + b5 A4 + b3 A2 + b1 I ]
!>     V =     A6 (b12 A6 + b10 A4 + b8 A2) + b6 A6 + b4 A4 + b2 A2 + b0 I
!>     (V - U) X = (V + U),   exp(K) = X^(2^s)
!>
!> Only the degree-13 branch is implemented: scipy additionally switches down
!> to Pade 3/5/7/9 for small norms purely to save flops, and using the higher
!> order there is if anything more accurate, never less.  CASSCF caps the
!> rotation at `max_rotation_norm` (default 0.2--0.5), so in practice s = 0 and
!> the whole rotation is 6 GEMMs plus one DGESV.  Pinned against
!> scipy.linalg.expm on random antisymmetric matrices to ~1e-14 relative.
!>
!> Arrays cross in the caller's C order (last index fastest).  A C-order
!> [n,n] buffer read as a Fortran (n,n) column-major array is the transpose,
!> so the orbital coefficients are gathered index-explicitly rather than
!> relying on that coincidence -- unlike the symmetric matrices in
!> casscf_kernel.F90, neither C nor exp(K) is symmetric.
module casscf_orbrot_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  !> Higham's theta_13: the largest ||A||_1 for which the degree-13 Pade
  !> approximant is backward stable to double-precision round-off.
  real(dp), parameter :: theta13 = 5.371920351148152_dp

  public :: casscf_orbital_rotate, casscf_expm

contains

  !> exp(a) for a general real n x n matrix, in place of `a`.
  !>
  !> `a` is a Fortran column-major (n,n) array in the ordinary mathematical
  !> sense (a(i,j) = A[i,j]); the caller is responsible for any C-order
  !> transposition.  Returns 0, or the LAPACK `info` if the Pade solve failed.
  function expm_inplace(n, a) result(info)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(0:n-1, 0:n-1)
    integer :: info

    real(dp), parameter :: b(0:13) = [ &
        64764752532480000.0_dp, 32382376266240000.0_dp, 7771770303897600.0_dp, &
        1187353796428800.0_dp,    129060195264000.0_dp,   10559470521600.0_dp, &
        670442572800.0_dp,           33522128640.0_dp,       1323241920.0_dp, &
        40840800.0_dp,                   960960.0_dp,            16380.0_dp, &
        182.0_dp,                             1.0_dp ]

    real(dp), allocatable :: a2(:,:), a4(:,:), a6(:,:), u(:,:), v(:,:), w(:,:)
    integer, allocatable :: ipiv(:)
    real(dp) :: anorm, colsum, scal
    integer :: i, j, s, k

    info = 0
    if (n <= 0) return

    ! ---- ||A||_1 = max_j sum_i |a(i,j)|
    anorm = 0.0_dp
    do j = 0, n - 1
      colsum = 0.0_dp
      do i = 0, n - 1
        colsum = colsum + abs(a(i, j))
      end do
      if (colsum > anorm) anorm = colsum
    end do

    ! ---- scaling: bring the norm inside the degree-13 Pade radius
    s = 0
    if (anorm > theta13) then
      s = ceiling(log(anorm / theta13) / log(2.0_dp))
      if (s < 0) s = 0
      scal = 1.0_dp / (2.0_dp ** s)
      a = a * scal
    end if

    allocate(a2(0:n-1, 0:n-1), a4(0:n-1, 0:n-1), a6(0:n-1, 0:n-1))
    allocate(u(0:n-1, 0:n-1), v(0:n-1, 0:n-1), w(0:n-1, 0:n-1))

    call dgemm('N', 'N', n, n, n, 1.0_dp, a,  n, a,  n, 0.0_dp, a2, n)
    call dgemm('N', 'N', n, n, n, 1.0_dp, a2, n, a2, n, 0.0_dp, a4, n)
    call dgemm('N', 'N', n, n, n, 1.0_dp, a2, n, a4, n, 0.0_dp, a6, n)

    ! ---- U = A [ A6 (b13 A6 + b11 A4 + b9 A2) + b7 A6 + b5 A4 + b3 A2 + b1 I ]
    w = b(13) * a6 + b(11) * a4 + b(9) * a2
    call dgemm('N', 'N', n, n, n, 1.0_dp, a6, n, w, n, 0.0_dp, u, n)
    u = u + b(7) * a6 + b(5) * a4 + b(3) * a2
    do i = 0, n - 1
      u(i, i) = u(i, i) + b(1)
    end do
    call dgemm('N', 'N', n, n, n, 1.0_dp, a, n, u, n, 0.0_dp, w, n)
    u = w

    ! ---- V = A6 (b12 A6 + b10 A4 + b8 A2) + b6 A6 + b4 A4 + b2 A2 + b0 I
    w = b(12) * a6 + b(10) * a4 + b(8) * a2
    call dgemm('N', 'N', n, n, n, 1.0_dp, a6, n, w, n, 0.0_dp, v, n)
    v = v + b(6) * a6 + b(4) * a4 + b(2) * a2
    do i = 0, n - 1
      v(i, i) = v(i, i) + b(0)
    end do

    ! ---- solve (V - U) X = (V + U);  a <- X
    a = v + u
    v = v - u
    allocate(ipiv(n))
    ! DGESV takes the default (ILP64) integer kind the linked LAPACK was built
    ! with, which is what `n` and `ipiv` already are.
    call dgesv(n, n, v, n, ipiv, a, n, info)
    deallocate(ipiv)
    if (info /= 0) then
      deallocate(a2, a4, a6, u, v, w)
      return
    end if

    ! ---- undo the scaling: X <- X^(2^s)
    do k = 1, s
      call dgemm('N', 'N', n, n, n, 1.0_dp, a, n, a, n, 0.0_dp, w, n)
      a = w
    end do

    deallocate(a2, a4, a6, u, v, w)
  end function expm_inplace

  !> exp(K) for a general real matrix handed over in the caller's C order.
  !>
  !> @param[in]  n     matrix dimension
  !> @param[in]  kin   K, C-order [n,n]
  !> @param[out] eout  exp(K), C-order [n,n]
  !> @return           0 on success, LAPACK info on a failed Pade solve
  !>
  !> Exposed so the Python pin can compare the engine's exponential against
  !> scipy.linalg.expm directly, without going through an orbital rotation.
  function casscf_expm(n, kin, eout) result(info) &
      bind(C, name="casscf_expm")
    integer(c_int32_t), value :: n
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: kin(0:*)
    real(dp), intent(inout) :: eout(0:*)
    integer(c_int32_t) :: info

    real(dp), allocatable :: kf(:,:)
    integer :: nn, i, j, ierr

    nn = int(n)
    info = 0_c_int32_t
    if (nn <= 0) return

    allocate(kf(0:nn-1, 0:nn-1))
    do i = 0, nn - 1
      do j = 0, nn - 1
        kf(i, j) = kin(int(i, i8)*int(nn, i8) + int(j, i8))
      end do
    end do

    ierr = expm_inplace(nn, kf)
    info = int(ierr, c_int32_t)
    if (ierr == 0) then
      do i = 0, nn - 1
        do j = 0, nn - 1
          eout(int(i, i8)*int(nn, i8) + int(j, i8)) = kf(i, j)
        end do
      end do
    end if
    deallocate(kf)
  end function casscf_expm

  !> Orbital rotation C_new = C exp(K(vec)) over the non-redundant pairs.
  !>
  !> @param[in]  nbf   number of orbitals
  !> @param[in]  npar  non-redundant rotation pairs
  !> @param[in]  pairs pair list, C-order [npar,2] as (p,q), the ordering of
  !>                   casscf.py `_nonredundant_pairs`
  !> @param[in]  vec   rotation amplitudes, [npar]
  !> @param[in]  cin   orbital coefficients, C-order [nbf,nbf] (columns are MOs)
  !> @param[out] cout  rotated coefficients, C-order [nbf,nbf]
  !> @return           0 on success, LAPACK info on a failed Pade solve
  !>
  !> K is built exactly as `_kappa_matrix` does, with `+=` / `-=` accumulation
  !> so a repeated pair sums rather than overwrites (the Python pin has the
  !> same semantics and one test exercises it).
  function casscf_orbital_rotate(nbf, npar, pairs, vec, cin, cout) result(info) &
      bind(C, name="casscf_orbital_rotate")
    integer(c_int32_t), value :: nbf, npar
    integer(c_int32_t), intent(in) :: pairs(0:*)
    real(dp), intent(in) :: vec(0:*), cin(0:*)
    real(dp), intent(inout) :: cout(0:*)
    integer(c_int32_t) :: info

    real(dp), allocatable :: kf(:,:), cf(:,:), rf(:,:)
    integer :: n, np, i, j, k, p, q, ierr
    real(dp) :: val

    n = int(nbf)
    np = int(npar)
    info = 0_c_int32_t
    if (n <= 0) return

    allocate(kf(0:n-1, 0:n-1), cf(0:n-1, 0:n-1), rf(0:n-1, 0:n-1))

    ! ---- K[p,q] += vec[k],  K[q,p] -= vec[k]
    kf = 0.0_dp
    do k = 0, np - 1
      p = int(pairs(2*k))
      q = int(pairs(2*k + 1))
      val = vec(k)
      kf(p, q) = kf(p, q) + val
      kf(q, p) = kf(q, p) - val
    end do

    ierr = expm_inplace(n, kf)
    info = int(ierr, c_int32_t)
    if (ierr /= 0) then
      deallocate(kf, cf, rf)
      return
    end if

    ! ---- gather C from the caller's C-order buffer into cf(i,j) = C[i,j]
    do i = 0, n - 1
      do j = 0, n - 1
        cf(i, j) = cin(int(i, i8)*int(n, i8) + int(j, i8))
      end do
    end do

    ! ---- R = C exp(K)
    call dgemm('N', 'N', n, n, n, 1.0_dp, cf, n, kf, n, 0.0_dp, rf, n)

    do i = 0, n - 1
      do j = 0, n - 1
        cout(int(i, i8)*int(n, i8) + int(j, i8)) = rf(i, j)
      end do
    end do

    deallocate(kf, cf, rf)
  end function casscf_orbital_rotate

end module casscf_orbrot_mod
