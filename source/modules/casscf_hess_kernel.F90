!> @brief Fortran engine for the CI-relaxation half of the analytic CASSCF
!> orbital Hessian of pyoqp casscf_hessian.py.
!>
!> The analytic Hessian replaces 2*npar CI solves with one dense active-space
!> diagonalization, and profiling puts ~78% of what remains in a single loop:
!> for every rotation pair k, `_apply_active_operator` applies the derivative
!> operator to the reference CI vector and projects it on the active eigenbasis,
!>
!>     sigma^(k)_a = sum_tu f^(k)_tu W_tua
!>                 + 1/2 sum_tub stack[t,u,a,b] x^(k)_tub
!>                 - 1/2 sum_tw (sum_u g^(k)_tuuw) W_twa
!>     x^(k)_tub   = sum_vw g^(k)_tuvw W_vwb
!>     amp[k,j]    = sum_a sigma^(k)_a vecs[a,j]
!>
!> with W_tua = (E_tu|c>)_a.  In Python that is npar separate einsums over the
!> (nact,nact,ndet,ndet) excitation stack, and the profile shows the time going
!> to array reshapes rather than to arithmetic.  Every one of these
!> contractions is a matrix product and the pair index k is just another free
!> column index, so the whole loop becomes a handful of DGEMMs with k batched:
!>
!>   1. sigma  = W . f^(k)                    one GEMM over all k
!>   2. x      = W . g^(k)                    one GEMM over all k
!>   3. sigma -= 1/2 W . gtr^(k)              one GEMM over all k
!>   4. sigma += 1/2 stack^T . x              nact^2 GEMMs, accumulating
!>   5. amp    = vecs . sigma                 one GEMM over all k
!>
!> Step 4 is written as one GEMM per (t,u) so the stack is consumed in place:
!> each (t,u) block is already a contiguous ndet x ndet matrix stored as (b,a),
!> and the matching slice of x is a strided submatrix, which sequence
!> association expresses as a leading dimension.  Transposing the stack into
!> ((t,u,b),a) order instead would cost a full nact^2*ndet^2 copy -- the very
!> memory traffic this kernel exists to remove.
!>
!> The x buffer is npar*nact^2*ndet doubles, which outgrows memory for the
!> larger active spaces the dense-spectrum path still admits, so the pair index
!> is processed in chunks sized to a fixed budget.  Results do not depend on
!> the chunking: the pairs are independent columns.
!>
!> All arrays cross the boundary in the caller's C order (last index fastest).
!> A C-order [m,n] buffer is the Fortran column-major (n,m) matrix, which is
!> what lets each contraction above be a GEMM with no repacking.
module casscf_hess_kernel_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  !> Budget for the per-chunk x buffer, in doubles (128 MiB).
  integer(i8), parameter :: x_budget = 16777216_i8

  public :: casscf_hess_amp

contains

  !> Projected derivative-operator amplitudes for the CI-relaxation Hessian.
  !>
  !> @param[in]  nact   active orbitals
  !> @param[in]  ndet   determinants in the active space
  !> @param[in]  npar   non-redundant rotation pairs
  !> @param[in]  stack  spin-summed excitation matrices, C-order
  !>                    [nact,nact,ndet,ndet]; (E_tu)_{a,b}
  !> @param[in]  fder   folded one-body derivative integrals, C-order
  !>                    [npar,nact,nact]
  !> @param[in]  gder   folded two-body derivative integrals, C-order
  !>                    [npar,nact,nact,nact,nact]
  !> @param[in]  wmat   E_tu applied to the reference CI vector, C-order
  !>                    [nact,nact,ndet]
  !> @param[in]  vecs   active-Hamiltonian eigenvectors, C-order [ndet,ndet];
  !>                    column j is eigenvector j
  !> @param[out] amp    projected amplitudes, C-order [npar,ndet]
  subroutine casscf_hess_amp(nact, ndet, npar, stack, fder, gder, wmat, vecs, amp) &
      bind(C, name="casscf_hess_amp")
    integer(c_int32_t), value :: nact, npar
    integer(i8), value :: ndet
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: stack(0:*), fder(0:*), gder(0:*), wmat(0:*), vecs(0:*)
    real(dp), intent(inout) :: amp(0:*)

    ! Default (ILP64) integer kind, matching the linked BLAS.
    integer :: na, na2, na4, np, nd, nb, ncols, kbase, nchunk
    integer :: chunk, kk, t, u, w, tu
    integer(i8) :: need
    real(dp) :: acc
    real(dp), allocatable :: sigma(:,:), xbuf(:), gtr(:)

    na = int(nact)
    np = int(npar)
    nd = int(ndet)
    na2 = na * na
    na4 = na2 * na2
    if (np <= 0 .or. nd <= 0 .or. na <= 0) return

    ! Chunk the pair index so the x buffer stays inside its budget, but never
    ! fewer than one pair per chunk.
    need = int(na2, i8) * int(nd, i8)
    nb = int(max(1_i8, min(int(np, i8), x_budget / max(1_i8, need))))
    nchunk = (np + nb - 1) / nb

    allocate(sigma(0:nd-1, 0:nb-1))
    allocate(xbuf(0:need*int(nb, i8) - 1_i8))
    allocate(gtr(0:int(na2, i8)*int(nb, i8) - 1_i8))

    do chunk = 0, nchunk - 1
      kbase = chunk * nb
      ncols = min(nb, np - kbase)

      ! ---- 1. sigma(a,k) = sum_(t,u) W(a,(t,u)) f^(k)((t,u))
      ! C-order wmat[(t,u),a] is the Fortran (a,(t,u)) matrix and C-order
      ! fder[k,(t,u)] is the Fortran ((t,u),k) matrix, so the pair index rides
      ! along as the free column index and the whole chunk is one call.
      call dgemm('N', 'N', nd, ncols, na2, 1.0_dp, &
                 wmat, nd, fder(int(kbase, i8)*int(na2, i8)), na2, &
                 0.0_dp, sigma, nd)

      ! ---- 3. sigma -= 1/2 sum_(t,w) W(a,(t,w)) gtr^(k)(t,w),
      !         gtr^(k)_tw = sum_u g^(k)_tuuw
      do kk = 0, ncols - 1
        do t = 0, na - 1
          do w = 0, na - 1
            acc = 0.0_dp
            do u = 0, na - 1
              acc = acc + gder((((int(kbase + kk, i8)*int(na, i8) + int(t, i8)) &
                                 *int(na, i8) + int(u, i8)) &
                                 *int(na, i8) + int(u, i8)) &
                                 *int(na, i8) + int(w, i8))
            end do
            gtr(int(t*na + w, i8) + int(kk, i8)*int(na2, i8)) = acc
          end do
        end do
      end do
      call dgemm('N', 'N', nd, ncols, na2, -0.5_dp, &
                 wmat, nd, gtr, na2, 1.0_dp, sigma, nd)

      ! ---- 2. x(b,(k,t,u)) = sum_(v,w) W(b,(v,w)) g^(k)((t,u),(v,w))
      call dgemm('N', 'N', nd, ncols*na2, na2, 1.0_dp, &
                 wmat, nd, gder(int(kbase, i8)*int(na4, i8)), na2, &
                 0.0_dp, xbuf, nd)

      ! ---- 4. sigma(a,k) += 1/2 sum_(t,u) sum_b stack[t,u,a,b] x(b,(k,t,u))
      ! The (t,u) stack block is contiguous and stored as (b,a), so 'T'
      ! contracts over b; the matching x columns are strided by na2*nd, which
      ! becomes the leading dimension.  Both operands are passed by element
      ! reference so no temporary is materialised.
      do tu = 0, na2 - 1
        call dgemm('T', 'N', nd, ncols, nd, 0.5_dp, &
                   stack(int(tu, i8)*int(nd, i8)*int(nd, i8)), nd, &
                   xbuf(int(tu, i8)*int(nd, i8)), na2*nd, &
                   1.0_dp, sigma, nd)
      end do

      ! ---- 5. amp(j,k) = sum_a vecs(j,a) sigma(a,k)
      call dgemm('N', 'N', nd, ncols, nd, 1.0_dp, &
                 vecs, nd, sigma, nd, &
                 0.0_dp, amp(int(kbase, i8)*int(nd, i8)), nd)
    end do

    deallocate(sigma, xbuf, gtr)
  end subroutine casscf_hess_amp

end module casscf_hess_kernel_mod
