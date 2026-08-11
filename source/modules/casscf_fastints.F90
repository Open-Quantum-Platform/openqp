!> Partial AO->MO integral transformations for the CASSCF hot path.
!>
!> One CASSCF gradient evaluation needs three things from the two-electron
!> integrals:
!>
!>   * the active-block MO ERIs `(tu|vw)` for the CI solve            (nact^4)
!>   * the one-general-index slice `(j t|u v)` for the generalized
!>     Fock's all-active correction                                (nbf nact^3)
!>   * Coulomb/exchange contractions `J(D)`/`K(D)` against densities
!>     (the frozen-core fold and the per-root mean field W)
!>
!> The original driver produced all three from a full nbf^4 MO transformation
!> (`mo_transform_eri`, 8 nbf^5 flops) per evaluation -- and the TRAH
!> converger's finite-difference Hessian-vector product calls the evaluation
!> twice per micro-iteration, which made the transform the single dominant
!> cost of every CASSCF run beyond minimal basis sets.
!>
!> Here the first two come from a shared partial transformation that contracts
!> three indices down to the active window (2 nbf^4 nact + O(nbf^3 nact^2)
!> flops -- nbf/(4 nact) times cheaper than the full transform), and the J/K
!> contractions run directly over the AO tensor with AO-basis densities, which
!> is what `build_jkw` already does over MO quantities.  The full transform
!> survives at the two once-per-run call sites (the final CI and the
!> canonicalization Fock) and in the analytic-Hessian build, which genuinely
!> consumes the full tensor.
!>
!> All tensors cross these routines in the caller's conventions: the AO ERI is
!> C-order [nbf,nbf,nbf,nbf] chemist (ab|cd), `cbuf` is the driver's orbital
!> buffer (Fortran (MO p, AO a), i.e. the C-order [AO, MO] coefficient matrix),
!> and the outputs land exactly in the layouts `fci_solve` and
!> `casscf_gfock_grad_w` read.
module casscf_fastints_mod
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: cas_partial_eri, cas_ao_density, cas_wmo_from_ao, cas_dmo_to_ao
  public :: cas_act_exchange_tensor, cas_act_fock
  public :: cas_fused_core_sweep, cas_partial_eri_from_t1

contains

  !> Core J/K and the pass-1 quarter transform in ONE blocked sweep of the AO
  !> ERI tensor.
  !>
  !> One gradient evaluation reads the nbf^4 tensor three times -- the core
  !> Coulomb DGEMV, the core exchange DGEMVs, and the first quarter
  !> transformation -- and each of those runs at the memory bandwidth, not the
  !> flop ceiling.  All three are row-range decomposable over the leading AO
  !> index, so this routine walks the tensor once in slabs and applies the
  !> three contractions per slab while it is still cache-resident: the DRAM
  !> traffic of the evaluation's dominant term drops ~3x and each piece is
  !> still its optimal BLAS call.
  !>
  !> Outputs exactly match the separate calls (same BLAS calls on the same
  !> operands, only reordered across disjoint output ranges):
  !>   jmat(q,j) = J[j,q], kmat(s,j) = K[j,s]   (build_jkw's layouts)
  !>   fao(a,b)  = h + J - K/2                  (the core mean field)
  !>   t1        = pass-1 tensor, C-order [t1,a,b,c] (cas_partial_eri's)
  subroutine cas_fused_core_sweep(n, na, nc, eri_ao, hcore, cbuf, dcore, &
                                  jmat, kmat, fao, t1)
    integer, intent(in) :: n, na, nc
    real(dp), intent(in) :: eri_ao(0:*), hcore(0:*)
    real(dp), contiguous, intent(in) :: cbuf(0:, 0:)
    real(dp), intent(in) :: dcore(0:n-1, 0:n-1)
    real(dp), intent(inout) :: jmat(0:n-1, 0:n-1), kmat(0:n-1, 0:n-1)
    real(dp), intent(out) :: fao(0:n-1, 0:n-1)
    real(dp), intent(inout) :: t1(0:*)

    integer :: a0, blk, nblk, i, j, n2
    integer(i8) :: off, n3

    ! ~32 MiB slabs: big enough that the pass-1 DGEMM keeps a useful shape,
    ! small enough that the second and third pass over a slab can come from
    ! the last-level cache rather than DRAM on current server parts.
    n2 = n * n
    n3 = int(n, i8) * int(n, i8) * int(n, i8)
    blk = max(1, int(32_i8 * 1024_i8 * 1024_i8 / max(8_i8 * n3, 1_i8)))

    do a0 = 0, n - 1, blk
      nblk = min(blk, n - a0)
      off = int(a0, i8) * n3
      ! J: columns (j,q) with j in the slab are the contiguous column range
      call dgemv('T', n2, nblk * n, 1.0_dp, eri_ao(off), n2, dcore(0, 0), 1, &
                 0.0_dp, jmat(0, a0), 1)
      ! K: one column per slab row, same block
      do j = a0, a0 + nblk - 1
        call dgemv('N', n, n2, 1.0_dp, eri_ao(int(j, i8) * n3), n, &
                   dcore(0, 0), 1, 0.0_dp, kmat(0, j), 1)
      end do
      ! pass 1: the slab's rows of every t1 column (ldc = n^3 keeps the
      ! caller's full-tensor layout)
      call dgemm('T', 'T', nblk * n2, na, n, 1.0_dp, eri_ao(off), n, &
                 cbuf(nc, 0), n, 0.0_dp, t1(int(a0, i8) * int(n2, i8)), &
                 int(n3))
    end do

    ! the same mean-field assembly build_jkw performs
    do i = 0, n - 1
      do j = 0, n - 1
        fao(i, j) = hcore(int(i, i8) * int(n, i8) + int(j, i8)) &
                  + jmat(j, i) - 0.5_dp * kmat(j, i)
      end do
    end do
  end subroutine cas_fused_core_sweep

  !> Passes 2..4 of `cas_partial_eri`, starting from a pass-1 tensor the
  !> caller already owns (the fused sweep above produces it).
  subroutine cas_partial_eri_from_t1(n, na, nc, cbuf, t1, t2, t3, eri_act, eact)
    integer, intent(in) :: n, na, nc
    real(dp), contiguous, intent(in) :: cbuf(0:, 0:)
    real(dp), intent(in) :: t1(0:*)
    real(dp), intent(inout) :: t2(0:*), t3(0:*)
    real(dp), intent(inout) :: eri_act(0:*), eact(0:*)

    integer :: na3

    na3 = na * na * na
    call dgemm('T', 'T', n * n * na, na, n, 1.0_dp, t1(0), n, cbuf(nc, 0), n, &
               0.0_dp, t2(0), n * n * na)
    call dgemm('T', 'T', n * na * na, na, n, 1.0_dp, t2(0), n, cbuf(nc, 0), n, &
               0.0_dp, t3(0), n * na * na)
    call dgemm('T', 'T', na3, na, n, 1.0_dp, t3(0), n, cbuf(nc, 0), n, &
               0.0_dp, eri_act(0), na3)
    call dgemm('N', 'N', n, na3, n, 1.0_dp, cbuf(0, 0), n, t3(0), n, &
               0.0_dp, eact(0), n)
  end subroutine cas_partial_eri_from_t1

  !> The exchange-shaped active slice  Y[v,w,a,b] = (a v | b w)  (v,w active;
  !> a,b AO), C-order, built from the pass-1 intermediate of
  !> `cas_partial_eri` (`t1`, C-order [w,a,b,c] = (ab|c w)).
  !>
  !> Together with the pass-2 intermediate (`t2`, [v,w,a,b] = (ab|v w)) this
  !> closes the active-density mean field without touching the nbf^4 AO tensor
  !> again: J(D_act) and K(D_act) become DGEMVs over these O(nbf^2 nact^2)
  !> slices (see `cas_act_fock`), so one gradient evaluation streams the AO
  !> tensor only for the CORE J/K and the shared pass-1 transform -- the
  !> per-state work no longer scales with the AO tensor at all.
  !>
  !> `ytr` is n^3*nact scratch for the index shuffle.
  subroutine cas_act_exchange_tensor(n, na, nc, cbuf, t1, ytr, yten)
    integer, intent(in) :: n, na, nc
    real(dp), contiguous, intent(in) :: cbuf(0:, 0:)
    real(dp), intent(in) :: t1(0:*)
    real(dp), intent(inout) :: ytr(0:*)
    real(dp), intent(out) :: yten(0:*)

    integer(i8) :: w, a, b, c, nn, na8, base_w, base_wa, base_wab

    nn = int(n, i8)
    na8 = int(na, i8)

    ! ytr[w,a,c,b] <- t1[w,a,b,c]: move the to-be-contracted index (b) last
    do w = 0_i8, na8 - 1_i8
      base_w = w * nn * nn * nn
      do a = 0_i8, nn - 1_i8
        base_wa = base_w + a * nn * nn
        do b = 0_i8, nn - 1_i8
          base_wab = base_wa + b * nn
          do c = 0_i8, nn - 1_i8
            ytr(base_wa + c * nn + b) = t1(base_wab + c)
          end do
        end do
      end do
    end do

    ! contract the (now last) index b with the active rows:
    ! yten C-order [v, w, a, c] = sum_b (ab|c w) C[b, nc+v] = (a v | c w)
    call dgemm('T', 'T', na * n * n, na, n, 1.0_dp, ytr(0), n, cbuf(nc, 0), n, &
               0.0_dp, yten(0), na * n * n)
  end subroutine cas_act_exchange_tensor

  !> Per-state AO mean field  F = F_core + J(D_act) - K(D_act)/2  from the
  !> active-slice tensors -- no nbf^4 stream.
  !>
  !>   J(D_act)[a,b] = sum_vw (ab|vw) gamma[v,w]   (DGEMV over t2)
  !>   K(D_act)[a,b] = sum_vw (av|bw) gamma[v,w]   (DGEMV over yten)
  !>
  !> `gamma` is the active 1-RDM (C-order [na,na], symmetric), `fcore` the
  !> core mean field from `build_jkw`/`cas_wmo_from_ao` (Fortran (n,n)), and
  !> `fout` the state mean field (Fortran (n,n)).  `jv`/`kv` are n^2 scratch.
  subroutine cas_act_fock(n, na, t2, yten, gamma, fcore, jv, kv, fout)
    integer, intent(in) :: n, na
    real(dp), intent(in) :: t2(0:*), yten(0:*)
    real(dp), intent(in) :: gamma(0:*)
    real(dp), intent(in) :: fcore(0:n-1, 0:n-1)
    real(dp), intent(inout) :: jv(0:*), kv(0:*)
    real(dp), intent(out) :: fout(0:n-1, 0:n-1)

    integer :: a, b, n2, na2

    n2 = n * n
    na2 = na * na

    ! t2 C-order [v,w,a,b] viewed column-major is M((a,b), (v,w)) with the
    ! (a,b) pair contiguous -- J is one DGEMV against the flattened gamma.
    call dgemv('N', n2, na2, 1.0_dp, t2(0), n2, gamma(0), 1, 0.0_dp, jv(0), 1)
    ! yten C-order [v,w,a,c]: same shape, same contraction -> K
    call dgemv('N', n2, na2, 1.0_dp, yten(0), n2, gamma(0), 1, 0.0_dp, kv(0), 1)

    ! jv/kv flat index (a*n + b) is the C-order [a,b] element; F is symmetric
    ! in exact arithmetic but assembled exactly as written.
    do a = 0, n - 1
      do b = 0, n - 1
        fout(a, b) = fcore(a, b) + jv(int(a, i8) * n + b) &
                   - 0.5_dp * kv(int(a, i8) * n + b)
      end do
    end do
  end subroutine cas_act_fock

  !> AO-basis image  D_ao = C D_mo C^T  of a full MO-space density.
  !> `dmo` is C-order [n,n] and symmetric (so its column-major view coincides);
  !> `tmp` is (n,n) scratch; `dao` is Fortran (n,n) (symmetric on exit).
  subroutine cas_dmo_to_ao(n, cbuf, dmo, tmp, dao)
    integer, intent(in) :: n
    real(dp), contiguous, intent(in) :: cbuf(0:, 0:)
    real(dp), intent(in) :: dmo(0:*)
    real(dp), intent(inout) :: tmp(0:n-1, 0:n-1)
    real(dp), intent(out) :: dao(0:n-1, 0:n-1)

    ! tmp(a, q) = sum_p C[a,p] D[p,q] with cbuf(p,a) = C[a,p]
    call dgemm('T', 'N', n, n, n, 1.0_dp, cbuf(0, 0), n, dmo(0), n, &
               0.0_dp, tmp(0, 0), n)
    ! dao(a, b) = sum_q tmp(a,q) C[b,q] = sum_q tmp(a,q) cbuf(q,b)
    call dgemm('N', 'N', n, n, n, 1.0_dp, tmp(0, 0), n, cbuf(0, 0), n, &
               0.0_dp, dao(0, 0), n)
  end subroutine cas_dmo_to_ao

  !> Active-space ERIs and the (j t|u v) slice from one shared 3-index
  !> partial transformation.
  !>
  !> Chain (each pass contracts the current C-order LAST index and prepends
  !> the new MO index in front, exactly the `mo_transform_eri` trick with
  !> rectangular coefficient blocks):
  !>
  !>   AO [a,b,c,d]            (ab|cd)
  !>   -> [t1,a,b,c]           (ab|c t1)      cost 2 n^4 na
  !>   -> [t2,t1,a,b]          (ab|t2 t1)     cost 2 n^3 na^2
  !>   -> [t3,t2,t1,a]         (a t3|t2 t1)   cost 2 n^2 na^3
  !>   -> eri_act [t4,t3,t2,t1] = (t4 t3|t2 t1)          cost 2 n na^4
  !>   -> eact    (j,(t3,t2,t1)) = (j t3|t2 t1)          cost 2 n^2 na^3
  !>
  !> @param[in]  n        basis functions
  !> @param[in]  na       active orbitals
  !> @param[in]  nc       inactive orbitals (active window starts at row nc)
  !> @param[in]  eri_ao   AO ERIs, C-order [n,n,n,n]
  !> @param[in]  cbuf     orbitals, Fortran (n, n): cbuf(p, a) = C[a, p]
  !> @param[out] t1       scratch, >= n^3 * na
  !> @param[out] t2       scratch, >= n^2 * na^2  (may alias unused tails)
  !> @param[out] t3       scratch, >= n * na^3
  !> @param[out] eri_act  active MO ERIs, C-order [na,na,na,na]
  !> @param[out] eact     Fortran (n, na^3): eact(j,(t*na+u)*na+v) = (j t|u v)
  subroutine cas_partial_eri(n, na, nc, eri_ao, cbuf, t1, t2, t3, eri_act, eact)
    integer, intent(in) :: n, na, nc
    real(dp), intent(in) :: eri_ao(0:*)
    real(dp), contiguous, intent(in) :: cbuf(0:, 0:)
    real(dp), intent(inout) :: t1(0:*), t2(0:*), t3(0:*)
    real(dp), intent(inout) :: eri_act(0:*), eact(0:*)

    integer :: n3, na3

    n3 = n * n * n
    na3 = na * na * na

    ! pass 1: contract d -> t1 over the active rows of cbuf
    call dgemm('T', 'T', n3, na, n, 1.0_dp, eri_ao(0), n, cbuf(nc, 0), n, &
               0.0_dp, t1(0), n3)
    ! pass 2: contract c
    call dgemm('T', 'T', n * n * na, na, n, 1.0_dp, t1(0), n, cbuf(nc, 0), n, &
               0.0_dp, t2(0), n * n * na)
    ! pass 3: contract b
    call dgemm('T', 'T', n * na * na, na, n, 1.0_dp, t2(0), n, cbuf(nc, 0), n, &
               0.0_dp, t3(0), n * na * na)
    ! pass 4a: contract a over the active rows -> (t4 t3|t2 t1)
    call dgemm('T', 'T', na3, na, n, 1.0_dp, t3(0), n, cbuf(nc, 0), n, &
               0.0_dp, eri_act(0), na3)
    ! pass 4b: contract a over ALL MOs, MO index leading -> eact(j, tuv)
    call dgemm('N', 'N', n, na3, n, 1.0_dp, cbuf(0, 0), n, t3(0), n, &
               0.0_dp, eact(0), n)
  end subroutine cas_partial_eri

  !> AO-basis density  D = 2 C_core^T C_core + C_act^T gamma C_act.
  !>
  !> `gamma` is the active 1-RDM, C-order [na,na] (symmetric); pass na=0 to get
  !> the bare core density.  `dout` is Fortran (n,n) (symmetric, so its C-order
  !> view coincides).  `tmp` is (na, n) scratch.
  subroutine cas_ao_density(n, nc, na, cbuf, gamma, dout, tmp)
    integer, intent(in) :: n, nc, na
    real(dp), contiguous, intent(in) :: cbuf(0:, 0:)
    real(dp), intent(in) :: gamma(0:*)
    real(dp), intent(out) :: dout(0:n-1, 0:n-1)
    real(dp), intent(inout) :: tmp(0:*)

    if (nc > 0) then
      call dgemm('T', 'N', n, n, nc, 2.0_dp, cbuf(0, 0), n, cbuf(0, 0), n, &
                 0.0_dp, dout(0, 0), n)
    else
      dout = 0.0_dp
    end if
    if (na > 0) then
      ! tmp(t, a) = sum_u gamma[t,u] cbuf(nc+u, a); gamma symmetric, so the
      ! C-order buffer is its own column-major view.
      call dgemm('N', 'N', na, n, na, 1.0_dp, gamma(0), na, cbuf(nc, 0), n, &
                 0.0_dp, tmp(0), na)
      call dgemm('T', 'N', n, n, na, 1.0_dp, cbuf(nc, 0), n, tmp(0), na, &
                 1.0_dp, dout(0, 0), n)
    end if
  end subroutine cas_ao_density

  !> MO mean field  W = C^T F_ao C  with  F_ao = h + J(D) - K(D)/2 built in the
  !> AO basis.  `fao` is the (n,n) AO mean field on exit (the caller reuses it
  !> for the frozen-core energy), `wmo` is Fortran (n,n): wmo(p,q) = W[p,q].
  !> `js`, `ks`, `tmp` are (n,n) scratch.
  subroutine cas_wmo_from_ao(n, eri_ao, hcore, cbuf, dao, js, ks, fao, tmp, wmo)
    use casscf_kernel_mod, only: build_jkw
    integer, intent(in) :: n
    real(dp), intent(in) :: eri_ao(0:*), hcore(0:*)
    real(dp), contiguous, intent(in) :: cbuf(0:, 0:)
    real(dp), intent(in) :: dao(0:n-1, 0:n-1)
    real(dp), intent(inout) :: js(0:n-1, 0:n-1), ks(0:n-1, 0:n-1)
    real(dp), intent(out) :: fao(0:n-1, 0:n-1)
    real(dp), intent(inout) :: tmp(0:n-1, 0:n-1)
    real(dp), intent(out) :: wmo(0:n-1, 0:n-1)

    ! build_jkw is basis-agnostic: with AO h, AO ERIs and an AO density it
    ! returns the AO mean field h + J - K/2.
    call build_jkw(n, dao, hcore, eri_ao, js, ks, fao)

    ! W = C^T F_ao C: tmp(p, b) = sum_a cbuf(p, a) fao(a, b)
    call dgemm('N', 'N', n, n, n, 1.0_dp, cbuf(0, 0), n, fao(0, 0), n, &
               0.0_dp, tmp(0, 0), n)
    call dgemm('N', 'T', n, n, n, 1.0_dp, tmp(0, 0), n, cbuf(0, 0), n, &
               0.0_dp, wmo(0, 0), n)
  end subroutine cas_wmo_from_ao

end module casscf_fastints_mod
