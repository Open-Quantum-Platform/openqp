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
!> Step 4 has two forms and picks between them on size.  Written as one GEMM
!> per (t,u) the stack is consumed in place: each (t,u) block is already a
!> contiguous ndet x ndet matrix stored as (b,a), and the matching slice of x
!> is a strided submatrix, which sequence association expresses as a leading
!> dimension.  Written as ((t,u,b), a) the whole of step 4 collapses into ONE
!> GEMM of inner dimension nact^2*ndet -- and getting there is NOT a transpose:
!> the (t,u) block is stored as (b,a) already, so the re-blocking is a copy of
!> contiguous ndet-long runs.
!>
!> The merged form removes nact^2 - 1 GEMM calls and their separate
!> accumulations into sigma, and pays one extra pass over the stack.  Measured
!> against the per-block loop in one process, best of many, nact 3-8 and
!> ndet 9-400: 1.35x at ndet=9, 1.26x at ndet=36 for every nact, 1.18-1.20x at
!> ndet=60, 1.10-1.12x at ndet=100, and 1.00x or below once the re-blocked
!> stack passes ~2 MiB -- past that the extra pass costs what the merge saves.
!> The per-block loop is kept above that threshold, and below a handful of
!> pairs, which is the other way the copy fails to amortise.  Both forms are
!> one call per chunk of the pair index; the re-blocking is done once for the
!> whole call, and its buffer is bounded by the same threshold.
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

  !> Largest re-blocked excitation stack worth building, in doubles (2 MiB).
  !> This is the measured crossover of the two step-4 forms, and it doubles as
  !> the memory bound on the copy: the buffer never exceeds this.
  integer(i8), parameter :: sblk_max = 262144_i8

  !> Fewest rotation pairs that amortise the re-blocking copy.  The copy is
  !> independent of npar while the work it feeds is linear in it, so below a
  !> handful of pairs the merged form loses: measured 0.52x at npar=2, 1.0x
  !> around npar=16, 1.19x at 64, 1.27x at 131.
  integer, parameter :: merge_min_pairs = 16

  public :: casscf_hess_amp, casscf_hess_wmat, casscf_hess_relax

contains

  !> W_tua = (E_tu |c>)_a, the excitation operators applied to one CI vector.
  !>
  !> Every CI-relaxation amplitude is contracted against this, so it is built
  !> once per averaged root.  The C-order stack [nact,nact,ndet,ndet] viewed
  !> column-major is the matrix (b, (t,u,a)) with leading dimension ndet, so
  !> the whole thing is a single DGEMV against the CI vector -- no per-block
  !> loop and no repacking.
  !>
  !> @param[in]  nact   active orbitals
  !> @param[in]  ndet   determinants
  !> @param[in]  stack  excitation matrices, C-order [nact,nact,ndet,ndet]
  !> @param[in]  civec  CI vector, [ndet]
  !> @param[out] wmat   C-order [nact,nact,ndet]
  subroutine casscf_hess_wmat(nact, ndet, stack, civec, wmat) &
      bind(C, name="casscf_hess_wmat")
    integer(c_int32_t), value :: nact
    integer(i8), value :: ndet
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: stack(0:*), civec(0:*)
    real(dp), intent(inout) :: wmat(0:*)

    integer :: na, nd, ncols

    na = int(nact)
    nd = int(ndet)
    if (na <= 0 .or. nd <= 0) return
    ncols = na * na * nd

    call dgemv('T', nd, ncols, 1.0_dp, stack, nd, civec, 1, 0.0_dp, wmat, 1)
  end subroutine casscf_hess_wmat

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
    integer :: chunk, kk, t, u, w, tu, b
    integer(i8) :: need, nblk
    real(dp) :: acc
    logical :: merged
    real(dp), allocatable :: sigma(:,:), xbuf(:), gtr(:), sblk(:)

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

    ! Re-block the stack into ((t,u,b), a) once for the whole call, when that
    ! is the cheaper of the two step-4 forms (see the module header).  Both
    ! sides of the copy are contiguous ndet-long runs.
    nblk = int(na2, i8) * int(nd, i8) * int(nd, i8)
    merged = nblk <= sblk_max .and. np >= merge_min_pairs
    if (merged) then
      allocate(sblk(0:nblk - 1_i8))
      do tu = 0, na2 - 1
        do b = 0, nd - 1
          sblk(int(tu, i8)*int(nd, i8) + int(b, i8)*int(na2, i8)*int(nd, i8) : &
               int(tu, i8)*int(nd, i8) + int(b, i8)*int(na2, i8)*int(nd, i8) &
               + int(nd, i8) - 1_i8) = &
              stack(int(tu, i8)*int(nd, i8)*int(nd, i8) + int(b, i8)*int(nd, i8) : &
                    int(tu, i8)*int(nd, i8)*int(nd, i8) + int(b, i8)*int(nd, i8) &
                    + int(nd, i8) - 1_i8)
        end do
      end do
    end if

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
      ! x is already ((t,u,b), k) with leading dimension na2*nd, so against the
      ! re-blocked stack this is one GEMM contracting over (t,u,b).  Otherwise
      ! one GEMM per (t,u): the block is contiguous and stored as (b,a), so 'T'
      ! contracts over b, and the matching x columns are strided by na2*nd,
      ! which becomes the leading dimension.  Both operands are passed by
      ! element reference so no temporary is materialised.
      if (merged) then
        call dgemm('T', 'N', nd, ncols, na2*nd, 0.5_dp, &
                   sblk, na2*nd, xbuf, na2*nd, 1.0_dp, sigma, nd)
      else
        do tu = 0, na2 - 1
          call dgemm('T', 'N', nd, ncols, nd, 0.5_dp, &
                     stack(int(tu, i8)*int(nd, i8)*int(nd, i8)), nd, &
                     xbuf(int(tu, i8)*int(nd, i8)), na2*nd, &
                     1.0_dp, sigma, nd)
        end do
      end if

      ! ---- 5. amp(j,k) = sum_a vecs(j,a) sigma(a,k)
      call dgemm('N', 'N', nd, ncols, nd, 1.0_dp, &
                 vecs, nd, sigma, nd, &
                 0.0_dp, amp(int(kbase, i8)*int(nd, i8)), nd)
    end do

    deallocate(sigma, xbuf, gtr)
    if (allocated(sblk)) deallocate(sblk)
  end subroutine casscf_hess_amp

  !> CI-relaxation accumulation for one averaged root.
  !>
  !> Builds the response weights
  !>
  !>     factor_j = 2 coef_j / (E_I - eps_j)
  !>
  !> and accumulates `hess += sum_j factor_j amp[:,j] amp[:,j]^T`, which is one
  !> DGEMM against a column-scaled copy rather than the ndet-long Python loop
  !> plus a masked NumPy product.
  !>
  !> `coef_j` follows casscf_hessian exactly: eigenstates that ARE the averaged
  !> root are skipped; couplings *within* the averaged set carry 0.5*(w_I - w_J)
  !> so equal weights cancel over the two ordered visits; everything else
  !> carries w_I.
  !>
  !> Degeneracy is not silently absorbed.  A denominator below `degen_tol` with
  !> a genuinely non-zero coupling means the state-averaged objective is not
  !> smooth and no orbital Hessian exists there, so the routine refuses:
  !> it returns j+1 and the caller raises.  A vanishing coupling at the same
  !> denominator is the spin/symmetry-forbidden case and is simply skipped.
  !>
  !> @param[in]     npar, ndet, navg, iavg  sizes and which averaged root this is
  !> @param[in]     ovl        squared overlaps, C-order [navg,ndet]
  !> @param[in]     weights    state-average weights, [navg]
  !> @param[in]     eps        active-Hamiltonian eigenvalues, [ndet]
  !> @param[in]     e_i        <c|H|c> for this root
  !> @param[in]     amp        projected amplitudes, C-order [npar,ndet]
  !> @param[inout]  hess       accumulated into, C-order [npar,npar]
  !> @return 0 on success, or j+1 for the offending degenerate eigenstate.
  function casscf_hess_relax(npar, ndet, navg, iavg, ovl, weights, eps, &
                             e_i, degen_tol, noise_tol, amp, hess) result(info) &
      bind(C, name="casscf_hess_relax")
    integer(c_int32_t), value :: npar, navg, iavg
    integer(i8), value :: ndet
    real(dp), value :: e_i, degen_tol, noise_tol
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: ovl(0:*), weights(0:*), eps(0:*), amp(0:*)
    real(dp), intent(inout) :: hess(0:*)
    integer(i8) :: info

    integer :: np, na, ia, m, mm, k
    integer(i8) :: j, nd
    real(dp) :: coef, denom, best, amax, w_i
    real(dp), allocatable :: factors(:), bmat(:,:)
    logical :: any_nz

    info = 0_i8
    np = int(npar)
    na = int(navg)
    ia = int(iavg)
    nd = ndet
    if (np <= 0 .or. nd <= 0_i8) return
    w_i = weights(ia)

    allocate(factors(0:nd-1))
    factors = 0.0_dp
    any_nz = .false.

    do j = 0_i8, nd - 1_i8
      ! ovl is C-order [navg, ndet]: entry (a, j) sits at a*ndet + j.
      if (ovl(int(ia, i8)*nd + j) > 0.5_dp) cycle   ! this eigenstate IS root I
      m = 0
      best = ovl(j)
      do mm = 1, na - 1
        if (ovl(int(mm, i8)*nd + j) > best) then
          best = ovl(int(mm, i8)*nd + j)
          m = mm
        end if
      end do
      if (best > 0.5_dp) then
        coef = 0.5_dp * (w_i - weights(m))
      else
        coef = w_i
      end if
      if (coef == 0.0_dp) cycle
      denom = e_i - eps(j)
      if (abs(denom) < degen_tol) then
        amax = 0.0_dp
        do k = 0, np - 1
          amax = max(amax, abs(amp(int(k, i8)*nd + j)))
        end do
        if (amax * abs(coef) < noise_tol) cycle    ! forbidden partner: exact zero
        info = j + 1_i8                            ! genuine degeneracy: refuse
        deallocate(factors)
        return
      end if
      factors(j) = 2.0_dp * coef / denom
      any_nz = .true.
    end do

    if (any_nz) then
      ! hess += (amp . diag(factors)) . amp^T.  C-order amp[k,j] is the Fortran
      ! (j,k) matrix, so scaling columns of the Fortran view scales rows of the
      ! caller's -- exactly the diag(factors) product -- and one DGEMM closes it.
      allocate(bmat(0:nd-1, 0:np-1))
      do k = 0, np - 1
        do j = 0_i8, nd - 1_i8
          bmat(j, k) = amp(int(k, i8)*nd + j) * factors(j)
        end do
      end do
      call dgemm('T', 'N', np, np, int(nd), 1.0_dp, bmat, int(nd), &
                 amp, int(nd), 1.0_dp, hess, np)
      deallocate(bmat)
    end if

    deallocate(factors)
  end function casscf_hess_relax

end module casscf_hess_kernel_mod
