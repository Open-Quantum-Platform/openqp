!> @brief Native assembly of the exact (SA-)CASSCF orbital-rotation Hessian.
!>
!> This is the Fortran counterpart of
!> `pyoqp/oqp/library/casscf_hessian.py::analytic_orbital_hessian`, and it is
!> orchestration, not new numerics: every expensive step already had a native
!> kernel, and the Python function existed only to call them in order.  What
!> moves here is the ordering plus the array plumbing between them --
!>
!>   * the state-averaged full-space 1-/2-RDMs (`_full_rdms`), built per root
!>     from `rdm1_spatial` / `rdm2_spatial`;
!>   * the fixed-CI Hessian columns and the folded derivative integrals
!>     (`casscf_hess_bmat`), and the `0.5 (B + B^T)` symmetrization;
!>   * the inactive-Fock fold of the *undisplaced* integrals (`_fold_active`);
!>   * the active-space Hamiltonian of those folded integrals -- which is just
!>     the CASCI Hamiltonian, so `fci_spin_orbital_integrals` +
!>     `fci_dense_build` build it, exactly as the Python routes to
!>     `fci._build_dense_hamiltonian` rather than owning a second
!>     implementation;
!>   * its dense spectrum (`oqp_dsyevd_f`), the averaged-root overlaps, and the
!>     per-root `casscf_hess_wmat` / `casscf_hess_amp` / `casscf_hess_relax`
!>     chain that accumulates the CI-relaxation part.
!>
!> The excitation-matrix stack is the one genuinely persistent object: it is
!> `nact^2 x ndet^2` doubles and depends only on the active space, so it is
!> built once per run and cached in the context (the Python `lru_cache`s it for
!> the same reason).
!>
!> Refusals, preserved exactly
!> ---------------------------
!> Two conditions must NOT return a number:
!>
!>   * an excitation stack beyond `CAS_HESS_MAX_STACK` (2 GiB, the Python's
!>     `_MAX_STACK_BYTES`) -- the dense-spectrum path does not scale there;
!>   * a genuine root degeneracy with non-zero orbital coupling, which
!>     `casscf_hess_relax` already refuses by returning `j+1`.  A finite Hessian
!>     there would be wrong: the state-averaged objective is not smooth.
!>
!> Both come back as negative statuses, and the caller (`casscf_driver`) turns
!> them into a declined run so the Python path re-raises its own message
!> verbatim.  Nothing here silently absorbs either case.
!>
!> Array conventions are the stack's usual ones: everything crosses in C order,
!> and "C-order [m,n] IS the Fortran (n,m) matrix" is used only where the
!> operand is symmetric.  The one asymmetric operand is the eigenvector matrix:
!> `casscf_hess_amp` wants the C-order buffer `vecs[a][j]` = component a of
!> eigenvector j, while LAPACK writes eigenvectors as Fortran COLUMNS, so the
!> transpose is taken explicitly (`vt = transpose(vecs)`).  Getting that
!> backwards is silent -- both shapes are square.
module casscf_anhess_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use casscf_hess_bmat_mod, only: casscf_hess_bmat
  use casscf_hess_kernel_mod, only: casscf_hess_amp, casscf_hess_wmat, &
                                    casscf_hess_relax
  use casscf_exc_stack_mod, only: casscf_excitation_stack
  use fci_setup_mod, only: fci_spin_orbital_integrals
  use fci_hamiltonian_mod, only: fci_sort_dets, fci_dense_build, oqp_dsyevd_f
  use rdm_kernel_mod, only: rdm1_spatial, rdm2_spatial
  ! DGEMM/DGEMV through the ILP64 wrapper layer (AGENTS.md rule 1).
  use oqp_linalg
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: cas_anhess_ctx_t, cas_anhess_init, cas_anhess_free, cas_anhess_build
  public :: CAS_HESS_OK, CAS_HESS_ERR_INPUT, CAS_HESS_ERR_ALLOC, &
            CAS_HESS_ERR_EIGEN, CAS_HESS_ERR_DEGEN, CAS_HESS_ERR_STACK

  !> Mirrors casscf_hessian.py: `_MAX_STACK_BYTES`, `_DEGENERACY_TOL`,
  !> `_COUPLING_NOISE`.  Fixed constants there, so fixed constants here.
  integer(i8), parameter :: CAS_HESS_MAX_STACK = 2_i8 * 1024_i8**3
  real(dp), parameter :: CAS_HESS_DEGEN_TOL = 1.0e-8_dp
  real(dp), parameter :: CAS_HESS_NOISE = 1.0e-8_dp

  integer(i8), parameter :: CAS_HESS_OK        =  0_i8
  integer(i8), parameter :: CAS_HESS_ERR_INPUT = -1_i8
  integer(i8), parameter :: CAS_HESS_ERR_ALLOC = -2_i8
  integer(i8), parameter :: CAS_HESS_ERR_EIGEN = -5_i8
  integer(i8), parameter :: CAS_HESS_ERR_DEGEN = -8_i8   !< non-smooth objective
  integer(i8), parameter :: CAS_HESS_ERR_STACK = -9_i8   !< active space too large

  !> Every buffer one analytic Hessian build needs, allocated once per run.
  !>
  !> `stack` is the expensive one (nact^2 * ndet^2 doubles) and depends only on
  !> the active space, so it survives every build.  The rest are scratch and are
  !> kept only to avoid an allocate/deallocate storm across macroiterations.
  type :: cas_anhess_ctx_t
    integer :: n = 0, nc = 0, na = 0, npar = 0, nstate = 0, nthreads = 1
    integer(i8) :: ndet = 0
    real(dp), allocatable :: stack(:)              !< [na,na,ndet,ndet] C-order
    real(dp), allocatable :: dmat(:), rdm2(:)      !< SA D [n,n], G [n,n,n,n]
    real(dp), allocatable :: gam1(:), gam2(:)      !< per-root active RDMs
    real(dp), allocatable :: gacc(:)               !< SA active 2-RDM block
    real(dp), allocatable :: bmat(:), fder(:), gder(:)
    real(dp), allocatable :: f0(:), g0(:), hspin(:), gspin(:)
    real(dp), allocatable :: hact(:), eps(:)
    real(dp), allocatable :: vecs(:,:), vt(:,:)
    real(dp), allocatable :: ovl(:,:), cir(:,:)
    real(dp), allocatable :: wmat(:), amp(:), cvec(:), hc(:)
    integer(i8), allocatable :: skeys(:), sperm(:)
  end type cas_anhess_ctx_t

contains

  !> Allocate the context and build the excitation stack for this active space.
  !>
  !> Returns CAS_HESS_ERR_STACK when the dense-spectrum stack would exceed the
  !> memory guard -- the same refusal the Python raises, and the caller turns it
  !> back into that message.
  function cas_anhess_init(ctx, n, nc, na, nalpha, nbeta, npar, nstate, &
                           nthreads, ndet, dets) result(status)
    type(cas_anhess_ctx_t), intent(inout) :: ctx
    integer, intent(in) :: n, nc, na, nalpha, nbeta, npar, nstate, nthreads
    integer(i8), intent(in) :: ndet
    ! assumed-SIZE throughout, as everywhere else in this stack: these arrays
    ! reach bind(C) kernels whose own dummies are assumed-size.
    integer(i8), intent(in) :: dets(0:*)
    integer(i8) :: status

    integer :: ierr, ns
    integer(i8) :: n2, n4, na2, na4, nd, nstk, ns2, ns4

    status = CAS_HESS_OK
    if (n <= 0 .or. na <= 0 .or. npar <= 0 .or. nstate <= 0 .or. ndet <= 0_i8) then
      status = CAS_HESS_ERR_INPUT
      return
    end if
    if (2 * na > 62) then          ! determinant encoding must fit int64
      status = CAS_HESS_ERR_INPUT
      return
    end if
    ! `nalpha`/`nbeta` are carried only so the caller states which active-space
    ! occupation the cached stack belongs to; the determinant list itself is
    ! what the kernel consumes.
    if (nalpha < 0 .or. nbeta < 0 .or. nalpha > na .or. nbeta > na) then
      status = CAS_HESS_ERR_INPUT
      return
    end if

    nd = ndet
    na2 = int(na, i8)**2
    na4 = int(na, i8)**4
    n2 = int(n, i8)**2
    n4 = int(n, i8)**4
    nstk = 8_i8 * na2 * nd * nd
    if (nstk > CAS_HESS_MAX_STACK .or. nstk < 0_i8) then
      status = CAS_HESS_ERR_STACK
      return
    end if

    ctx%n = n
    ctx%nc = nc
    ctx%na = na
    ctx%npar = npar
    ctx%nstate = nstate
    ctx%nthreads = max(1, nthreads)
    ctx%ndet = nd
    ns = 2 * na
    ns2 = int(ns, i8)**2
    ns4 = int(ns, i8)**4

    call cas_anhess_free(ctx)
    allocate(ctx%stack(0:na2*nd*nd - 1_i8), &
             ctx%dmat(0:n2 - 1_i8), ctx%rdm2(0:n4 - 1_i8), &
             ctx%gam1(0:na2 - 1_i8), ctx%gam2(0:na4 - 1_i8), &
             ctx%gacc(0:na4 - 1_i8), &
             ctx%bmat(0:int(npar, i8)**2 - 1_i8), &
             ctx%fder(0:int(npar, i8)*na2 - 1_i8), &
             ctx%gder(0:int(npar, i8)*na4 - 1_i8), &
             ctx%f0(0:na2 - 1_i8), ctx%g0(0:na4 - 1_i8), &
             ctx%hspin(0:ns2 - 1_i8), ctx%gspin(0:ns4 - 1_i8), &
             ctx%hact(0:nd*nd - 1_i8), ctx%eps(0:nd - 1_i8), &
             ctx%vecs(0:nd-1, 0:nd-1), ctx%vt(0:nd-1, 0:nd-1), &
             ctx%ovl(0:nd-1, 0:nstate-1), ctx%cir(0:nd-1, 0:nstate-1), &
             ctx%wmat(0:na2*nd - 1_i8), ctx%amp(0:int(npar, i8)*nd - 1_i8), &
             ctx%cvec(0:nd - 1_i8), ctx%hc(0:nd - 1_i8), &
             ctx%skeys(nd), ctx%sperm(nd), stat=ierr)
    if (ierr /= 0) then
      status = CAS_HESS_ERR_ALLOC
      call cas_anhess_free(ctx)
      return
    end if

    ! The stack is written in full (including its structural zeros) by the
    ! kernel, which clears it itself in parallel; nothing to pre-zero here.
    if (casscf_excitation_stack(int(na, c_int32_t), nd, dets, ctx%stack) &
        /= 0_i8) then
      status = CAS_HESS_ERR_ALLOC
      call cas_anhess_free(ctx)
      return
    end if
  end function cas_anhess_init

  subroutine cas_anhess_free(ctx)
    type(cas_anhess_ctx_t), intent(inout) :: ctx
    if (allocated(ctx%stack)) deallocate(ctx%stack)
    if (allocated(ctx%dmat)) deallocate(ctx%dmat)
    if (allocated(ctx%rdm2)) deallocate(ctx%rdm2)
    if (allocated(ctx%gam1)) deallocate(ctx%gam1)
    if (allocated(ctx%gam2)) deallocate(ctx%gam2)
    if (allocated(ctx%gacc)) deallocate(ctx%gacc)
    if (allocated(ctx%bmat)) deallocate(ctx%bmat)
    if (allocated(ctx%fder)) deallocate(ctx%fder)
    if (allocated(ctx%gder)) deallocate(ctx%gder)
    if (allocated(ctx%f0)) deallocate(ctx%f0)
    if (allocated(ctx%g0)) deallocate(ctx%g0)
    if (allocated(ctx%hspin)) deallocate(ctx%hspin)
    if (allocated(ctx%gspin)) deallocate(ctx%gspin)
    if (allocated(ctx%hact)) deallocate(ctx%hact)
    if (allocated(ctx%eps)) deallocate(ctx%eps)
    if (allocated(ctx%vecs)) deallocate(ctx%vecs)
    if (allocated(ctx%vt)) deallocate(ctx%vt)
    if (allocated(ctx%ovl)) deallocate(ctx%ovl)
    if (allocated(ctx%cir)) deallocate(ctx%cir)
    if (allocated(ctx%wmat)) deallocate(ctx%wmat)
    if (allocated(ctx%amp)) deallocate(ctx%amp)
    if (allocated(ctx%cvec)) deallocate(ctx%cvec)
    if (allocated(ctx%hc)) deallocate(ctx%hc)
    if (allocated(ctx%skeys)) deallocate(ctx%skeys)
    if (allocated(ctx%sperm)) deallocate(ctx%sperm)
  end subroutine cas_anhess_free

  !> One analytic Hessian at the orbitals whose MO integrals are `h1e`/`eri`.
  !>
  !> @param[in]  pairs   non-redundant rotation pairs, C-order [npar,2]
  !> @param[in]  weights state-average weights, [nstate]
  !> @param[in]  roots   0-based CI root indices, [nstate]
  !> @param[in]  nroot   columns of `ci`
  !> @param[in]  dets    determinant list in CI ordering
  !> @param[in]  h1e     MO one-electron integrals, C-order [n,n]
  !> @param[in]  eri     MO chemist ERIs, C-order [n,n,n,n]
  !> @param[in]  ci      CI vectors, C-order [ndet,nroot]
  !> @param[out] hess    C-order [npar,npar]
  function cas_anhess_build(ctx, pairs, weights, roots, nroot, dets, h1e, eri, &
                            ci, hess) result(status)
    type(cas_anhess_ctx_t), intent(inout) :: ctx
    ! bind(C) kernels downstream take assumed-SIZE dummies, so these are
    ! assumed-size too and reach them by sequence association.
    integer(c_int32_t), intent(in) :: pairs(0:*)
    real(dp), intent(in) :: weights(0:*)
    integer, intent(in) :: roots(0:*)
    integer, intent(in) :: nroot
    integer(i8), intent(in) :: dets(0:*)
    real(dp), intent(in) :: h1e(0:*), eri(0:*), ci(0:*)
    real(dp), intent(inout) :: hess(0:*)
    integer(i8) :: status

    integer :: n, nc, na, npar, nstate, k, r, i, j, p, q, ierr
    integer(i8) :: nd, cap, rc, idx, na2, na4, n2
    real(dp) :: w, val, s1, s2, nrm, e_i

    status = CAS_HESS_OK
    n = ctx%n
    nc = ctx%nc
    na = ctx%na
    npar = ctx%npar
    nstate = ctx%nstate
    nd = ctx%ndet
    na2 = int(na, i8)**2
    na4 = int(na, i8)**4
    n2 = int(n, i8)**2

    ! ---------------------------------------------------- state-averaged RDMs
    ! `casscf._full_rdms` per root, accumulated with the state-average weights
    ! in root order.  The separable part of G is nonlinear in D, so it must be
    ! formed per root and then averaged -- averaging D first is a different
    ! (wrong) object.
    ctx%dmat = 0.0_dp
    ctx%rdm2 = 0.0_dp
    ctx%gacc = 0.0_dp
    cap = nd * int(2 * na, i8)**2 + 1_i8
    do k = 0, nstate - 1
      r = roots(k)
      w = weights(k)
      do idx = 0_i8, nd - 1_i8
        ctx%cvec(idx) = ci(idx * int(nroot, i8) + int(r, i8))
      end do
      call rdm1_spatial(int(na, c_int32_t), nd, dets, ctx%cvec, ctx%gam1)
      rc = rdm2_spatial(int(na, c_int32_t), nd, dets, ctx%cvec, cap, ctx%gam2, &
                        0_c_int32_t)
      if (rc /= 0_i8) then
        status = CAS_HESS_ERR_ALLOC
        return
      end if
      call accumulate_rdms(ctx, n, nc, na, w)
    end do
    ! the active block of G is the averaged active 2-RDM, not the separable form
    do i = 0, na - 1
      do j = 0, na - 1
        do p = 0, na - 1
          do q = 0, na - 1
            ctx%rdm2(((int(nc + i, i8)*int(n, i8) + int(nc + j, i8))*int(n, i8) &
                      + int(nc + p, i8))*int(n, i8) + int(nc + q, i8)) = &
                ctx%gacc(((int(i, i8)*int(na, i8) + int(j, i8))*int(na, i8) &
                          + int(p, i8))*int(na, i8) + int(q, i8))
          end do
        end do
      end do
    end do

    ! ------------------------------ part 1: fixed-CI Hessian + derivative fold
    call casscf_hess_bmat(int(n, c_int32_t), int(nc, c_int32_t), &
                          int(na, c_int32_t), int(npar, c_int32_t), pairs, &
                          ctx%dmat, ctx%rdm2, h1e, eri, ctx%bmat, ctx%fder, &
                          ctx%gder)
    do i = 0, npar - 1
      do j = 0, npar - 1
        hess(int(i, i8)*int(npar, i8) + int(j, i8)) = &
            0.5_dp * (ctx%bmat(int(i, i8)*int(npar, i8) + int(j, i8)) &
                    + ctx%bmat(int(j, i8)*int(npar, i8) + int(i, i8)))
      end do
    end do

    ! ------------------------------------- part 2: CI relaxation, undisplaced
    ! `_fold_active(h1e, eri, ncore, nact)`: the two core reductions are kept
    ! as SEPARATE sums, matching the two NumPy einsums the reference forms.
    do i = 0, na - 1
      do j = 0, na - 1
        s1 = 0.0_dp
        do k = 0, nc - 1
          s1 = s1 + eri(((int(nc + i, i8)*int(n, i8) + int(nc + j, i8)) &
                          *int(n, i8) + int(k, i8))*int(n, i8) + int(k, i8))
        end do
        s2 = 0.0_dp
        do k = 0, nc - 1
          s2 = s2 + eri(((int(nc + i, i8)*int(n, i8) + int(k, i8)) &
                          *int(n, i8) + int(k, i8))*int(n, i8) + int(nc + j, i8))
        end do
        val = h1e(int(nc + i, i8)*int(n, i8) + int(nc + j, i8)) + 2.0_dp * s1
        ctx%f0(int(i, i8)*int(na, i8) + int(j, i8)) = val - s2
      end do
    end do
    do i = 0, na - 1
      do j = 0, na - 1
        do p = 0, na - 1
          do q = 0, na - 1
            ctx%g0(((int(i, i8)*int(na, i8) + int(j, i8))*int(na, i8) &
                    + int(p, i8))*int(na, i8) + int(q, i8)) = &
                eri(((int(nc + i, i8)*int(n, i8) + int(nc + j, i8))*int(n, i8) &
                     + int(nc + p, i8))*int(n, i8) + int(nc + q, i8))
          end do
        end do
      end do
    end do

    ! The determinant-basis matrix of `sum f E + 1/2 sum g (EE - dE)` IS the
    ! CASCI Hamiltonian of the folded integrals, so the validated CI builder
    ! produces it -- no second implementation of the same operator.
    call fci_spin_orbital_integrals(int(na, c_int32_t), ctx%f0, ctx%g0, &
                                    ctx%hspin, ctx%gspin, &
                                    int(ctx%nthreads, c_int32_t))
    call fci_sort_dets(nd, dets, ctx%skeys, ctx%sperm)
    call fci_dense_build(2 * na, nd, dets, ctx%skeys, ctx%sperm, ctx%hspin, &
                         ctx%gspin, 0.0_dp, ctx%hact, ctx%nthreads)

    ! hact is symmetric, so its C-order buffer and its column-major view are
    ! the same matrix and DSYEVD may read the copy in place.
    do idx = 0_i8, nd - 1_i8
      ctx%vecs(0:nd-1, idx) = ctx%hact(idx*nd : idx*nd + nd - 1_i8)
    end do
    call oqp_dsyevd_f(int(nd), ctx%vecs, ctx%eps, ierr)
    if (ierr /= 0) then
      status = CAS_HESS_ERR_EIGEN
      return
    end if
    ! `casscf_hess_amp` wants the C-order buffer vecs[a][j]; LAPACK left the
    ! eigenvectors in Fortran COLUMNS, so this transpose is mandatory and
    ! silent if omitted (every shape here is square).
    ctx%vt = transpose(ctx%vecs)

    ! ovl(j, i) is the C-order [nstate, ndet] buffer casscf_hess_relax reads.
    do k = 0, nstate - 1
      r = roots(k)
      do idx = 0_i8, nd - 1_i8
        ctx%cir(idx, k) = ci(idx * int(nroot, i8) + int(r, i8))
      end do
    end do
    call dgemm('T', 'N', int(nd), nstate, int(nd), 1.0_dp, ctx%vecs, int(nd), &
               ctx%cir, int(nd), 0.0_dp, ctx%ovl, int(nd))
    ctx%ovl = ctx%ovl * ctx%ovl

    do k = 0, nstate - 1
      r = roots(k)
      do idx = 0_i8, nd - 1_i8
        ctx%cvec(idx) = ci(idx * int(nroot, i8) + int(r, i8))
      end do
      nrm = sqrt(sum(ctx%cvec * ctx%cvec))
      if (nrm <= 0.0_dp) then
        status = CAS_HESS_ERR_INPUT
        return
      end if
      ctx%cvec = ctx%cvec / nrm
      ! hact is symmetric so 'N' and 'T' agree; e_i = <c|H|c>
      call dgemv('N', int(nd), int(nd), 1.0_dp, ctx%hact, int(nd), ctx%cvec, 1, &
                 0.0_dp, ctx%hc, 1)
      e_i = sum(ctx%cvec * ctx%hc)

      call casscf_hess_wmat(int(na, c_int32_t), nd, ctx%stack, ctx%cvec, ctx%wmat)
      call casscf_hess_amp(int(na, c_int32_t), nd, int(npar, c_int32_t), &
                           ctx%stack, ctx%fder, ctx%gder, ctx%wmat, ctx%vt, &
                           ctx%amp)
      rc = casscf_hess_relax(int(npar, c_int32_t), nd, int(nstate, c_int32_t), &
                             int(k, c_int32_t), ctx%ovl, weights, ctx%eps, &
                             e_i, CAS_HESS_DEGEN_TOL, CAS_HESS_NOISE, &
                             ctx%amp, hess)
      if (rc /= 0_i8) then
        ! A genuine degeneracy with live coupling.  Refuse -- a finite number
        ! here would be wrong -- and let the caller surface the Python message.
        status = CAS_HESS_ERR_DEGEN
        return
      end if
    end do
  end function cas_anhess_build

  !> `casscf._full_rdms` for one root, accumulated with weight `w`.
  !>
  !> D gets 2 on the inactive diagonal and the active 1-RDM on the active
  !> block; G is the separable form `D_pq D_rs - 1/2 D_ps D_rq` everywhere
  !> (its active block is overwritten by the averaged active 2-RDM afterwards,
  !> which is why `gacc` accumulates separately in the same root order).
  subroutine accumulate_rdms(ctx, n, nc, na, w)
    type(cas_anhess_ctx_t), intent(inout) :: ctx
    integer, intent(in) :: n, nc, na
    real(dp), intent(in) :: w

    real(dp), allocatable :: dr(:,:)
    integer :: i, j, p, q
    integer(i8) :: base_p, base_pq, base_pqr
    real(dp) :: val

    allocate(dr(0:n-1, 0:n-1))
    dr = 0.0_dp
    do i = 0, nc - 1
      dr(i, i) = 2.0_dp
    end do
    do i = 0, na - 1
      do j = 0, na - 1
        dr(nc + i, nc + j) = ctx%gam1(int(i, i8)*int(na, i8) + int(j, i8))
      end do
    end do
    do i = 0, n - 1
      do j = 0, n - 1
        ctx%dmat(int(i, i8)*int(n, i8) + int(j, i8)) = &
            ctx%dmat(int(i, i8)*int(n, i8) + int(j, i8)) + w * dr(i, j)
      end do
    end do

    !$omp parallel do collapse(2) default(shared) private(p, q, i, j, val, &
    !$omp   base_p, base_pq, base_pqr) if(int(n, i8)**4 >= 262144_i8)
    do p = 0, n - 1
      do q = 0, n - 1
        base_p = int(p, i8) * int(n, i8)**3
        base_pq = base_p + int(q, i8) * int(n, i8)**2
        do i = 0, n - 1
          base_pqr = base_pq + int(i, i8) * int(n, i8)
          do j = 0, n - 1
            val = dr(p, q) * dr(i, j) - 0.5_dp * (dr(p, j) * dr(i, q))
            ctx%rdm2(base_pqr + int(j, i8)) = &
                ctx%rdm2(base_pqr + int(j, i8)) + w * val
          end do
        end do
      end do
    end do
    !$omp end parallel do

    do i = 0, na - 1
      do j = 0, na - 1
        do p = 0, na - 1
          do q = 0, na - 1
            ctx%gacc(((int(i, i8)*int(na, i8) + int(j, i8))*int(na, i8) &
                      + int(p, i8))*int(na, i8) + int(q, i8)) = &
                ctx%gacc(((int(i, i8)*int(na, i8) + int(j, i8))*int(na, i8) &
                          + int(p, i8))*int(na, i8) + int(q, i8)) &
                + w * ctx%gam2(((int(i, i8)*int(na, i8) + int(j, i8)) &
                                *int(na, i8) + int(p, i8))*int(na, i8) &
                                + int(q, i8))
          end do
        end do
      end do
    end do
    deallocate(dr)
  end subroutine accumulate_rdms

end module casscf_anhess_mod
