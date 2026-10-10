!> @brief Fortran engine for the AO -> MO integral transformation of pyoqp
!> `fci._transform_integrals`.
!>
!> Every CASSCF gradient evaluation, every CASCI/FCI energy and every CASPT2 /
!> NEVPT2 assembly starts by transforming the one-electron core Hamiltonian and
!> the dense AO ERI tensor into the current MO basis.  The Python reference is
!> a single `np.einsum("up,vq,wr,xs,uvwx->pqrs", ...)` with `optimize=True`;
!> it stays in the tree as the numerical pin and as the fallback when this
!> symbol is unavailable.
!>
!> Two bind(C) entry points:
!>
!>  * mo_transform_h1e -- C^T h C, one pair of DGEMMs
!>  * mo_transform_eri -- the four quarter transformations
!>
!> The four-index transform is done as four DGEMMs and nothing else.  Writing
!> the AO tensor as a matrix over (first index) x (remaining three), each
!> quarter transformation is
!>
!>     T(J, s) = sum_x A(x, J) C(s, x)
!>
!> which is DGEMM('T','T', nbf^3, nbf, nbf).  Taking the result transposed is
!> the whole trick: it leaves the next index to be contracted sitting in the
!> leading position, so the index roll between steps is free and no explicit
!> nbf^4 transpose is ever materialized.  After four passes the indices have
!> cycled exactly once and the result is back in the caller's layout.
!>
!> That matters because the naive four-GEMM chain loses to numpy: measured at
!> nbf=24..58 an explicit chain with `ascontiguousarray` rolls between steps
!> runs at 0.48-0.78x of the einsum, since the copies cost more than the
!> einsum's own bookkeeping.  The transposed-output formulation removes them.
!>
!> Memory: one internal nbf^4 buffer, not two.  The passes ping-pong between it
!> and the caller's output array (in -> work -> out -> work -> out), which is
!> why the step count has to be even.
!>
!> All arrays cross the boundary in the caller's C order (last index fastest).
!> A C-order [nbf,nbf] `coeff` with coeff[u,p] (AO u, MO p) is the Fortran
!> column-major matrix C(p,u), i.e. already the transposed coefficient matrix
!> the contractions above want, so no transpose of C is needed either.
module mo_transform_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: mo_transform_h1e, mo_transform_eri

contains

  !> MO core Hamiltonian h1e = C^T hcore C.
  !>
  !> @param[in]  nbf    number of basis functions (== number of MOs)
  !> @param[in]  hcore  AO core Hamiltonian, C-order [nbf,nbf] (symmetric)
  !> @param[in]  coeff  MO coefficients, C-order [nbf,nbf], coeff[u,p]
  !> @param[out] h1e    MO core Hamiltonian, C-order [nbf,nbf]
  subroutine mo_transform_h1e(nbf, hcore, coeff, h1e) &
      bind(C, name="mo_transform_h1e")
    integer(c_int32_t), value :: nbf
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: hcore(0:*), coeff(0:*)
    real(dp), intent(inout) :: h1e(0:*)

    ! Default (ILP64) integer kind, matching the linked BLAS.
    integer :: n
    real(dp), allocatable :: work(:,:)

    n = int(nbf)
    if (n <= 0) return
    allocate(work(0:n-1, 0:n-1))

    ! Fortran views: Cf(p,u) = coeff[u,p], Hf(v,u) = hcore[u,v] = Hf(u,v).
    ! work(p,v) = sum_u Cf(p,u) Hf(u,v)
    call dgemm('N', 'N', n, n, n, 1.0_dp, coeff(0), n, hcore(0), n, &
               0.0_dp, work(0,0), n)
    ! h1e(p,q) = sum_v work(p,v) Cf(q,v).  Symmetric, so the C-order buffer
    ! and the column-major result coincide.
    call dgemm('N', 'T', n, n, n, 1.0_dp, work(0,0), n, coeff(0), n, &
               0.0_dp, h1e(0), n)

    deallocate(work)
  end subroutine mo_transform_h1e

  !> MO two-electron integrals (chemists' (pq|rs)) from the AO tensor.
  !>
  !> eri_mo[p,q,r,s] = sum_{uvwx} C[u,p] C[v,q] C[w,r] C[x,s] eri_ao[u,v,w,x]
  !>
  !> @param[in]  nbf     number of basis functions (== number of MOs)
  !> @param[in]  eri_ao  AO ERIs, C-order [nbf,nbf,nbf,nbf]
  !> @param[in]  coeff   MO coefficients, C-order [nbf,nbf], coeff[u,p]
  !> @param[out] eri_mo  MO ERIs, C-order [nbf,nbf,nbf,nbf]
  !> @return     0 on success, -1 when the internal buffer cannot be allocated
  !>             (the caller then falls back to the Python reference)
  function mo_transform_eri(nbf, eri_ao, coeff, eri_mo) result(status) &
      bind(C, name="mo_transform_eri")
    integer(c_int32_t), value :: nbf
    real(dp), intent(in) :: eri_ao(0:*), coeff(0:*)
    real(dp), intent(inout) :: eri_mo(0:*)
    integer(i8) :: status

    integer :: n, n3, step, ierr
    integer(i8) :: n4
    real(dp), allocatable :: work(:)

    status = 0_i8
    n = int(nbf)
    if (n <= 0) return
    n3 = n * n * n
    n4 = int(n, i8) * int(n3, i8)

    allocate(work(0:n4-1), stat=ierr)
    if (ierr /= 0) then
      status = -1_i8
      return
    end if

    ! Pass 1 reads the caller's input, so it is peeled out of the loop; after
    ! it the passes alternate work -> eri_mo -> work -> eri_mo.  Four passes
    ! is even, so the last one lands in eri_mo.
    !
    ! Each pass is T(J,s) = sum_x A(x,J) Cf(s,x), i.e. A^T C^T with
    ! A treated as (K=n) x (M=n^3) and Cf as (N=n) x (K=n).
    call dgemm('T', 'T', n3, n, n, 1.0_dp, eri_ao(0), n, coeff(0), n, &
               0.0_dp, work(0), n3)
    do step = 2, 4
      if (mod(step, 2) == 0) then
        call dgemm('T', 'T', n3, n, n, 1.0_dp, work(0), n, coeff(0), n, &
                   0.0_dp, eri_mo(0), n3)
      else
        call dgemm('T', 'T', n3, n, n, 1.0_dp, eri_mo(0), n, coeff(0), n, &
                   0.0_dp, work(0), n3)
      end if
    end do

    deallocate(work)
  end function mo_transform_eri

end module mo_transform_mod
