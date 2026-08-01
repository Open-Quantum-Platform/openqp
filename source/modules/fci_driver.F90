!> @brief One bind(C) entry point for a complete determinant-CI solve.
!>
!> The wavefunction stack on this branch exposed ~15 fine-grained kernels that a
!> Python loop orchestrated: for one CASSCF macroiteration the boundary was
!> crossed once for the AO->MO transform, once for the core fold, once for the
!> spin-orbital expansion, once per Davidson application, once per subspace
!> diagonalization and once per spin diagnostic -- hundreds of crossings for a
!> single CI solve, with the iteration control, the solver dispatch and the root
!> selection all living in Python.  Every other method in OpenQP exposes one C
!> entry point per method (`hf_energy(inf)`, `mp2_energy(inf)`); this restores
!> that shape for the CI layer, which is the piece that was blocking a single
!> `casscf_energy(inf)`.
!>
!> `fci_solve` runs active-space setup, the core fold, the spin-orbital
!> expansion, determinant generation, the dense/Davidson dispatch, the Davidson
!> iteration, the spin diagnostics and root selection without returning to the
!> caller.  It is deliberately **handle-free**, taking plain arrays like
!> `qdpt2_stream_kernel` rather than `struct oqp_handle_t *`:
!>
!>   * its dominant caller is the CASSCF microiteration, whose MO integrals are
!>     a *trial point* produced from a candidate orbital set inside a line
!>     search.  Those integrals do not live in `infos%dat` and putting them
!>     there would add a full nbf^4 copy per evaluation to save nothing;
!>   * the end state `casscf_energy(inf)` is itself Fortran, and will hold
!>     `h1e`/`eri` in local arrays.  A handle-taking `fci_solve` would have to
!>     be unwrapped again to be callable from it -- i.e. it would block the very
!>     thing it is meant to enable.
!>
!> The handle is read once, above this layer, by whatever obtains `OQP::Hcore`,
!> `OQP::AO_ERI` and `OQP::VEC_MO_A`; from there down the CI path is pure
!> array-in/array-out.
!>
!> **Option schema.**  `[cas]`/`[ci]` are Python-side sections with no
!> representation in `source/types.F90`, and there are ~111 keys across the
!> wavefunction sections.  Bolting them onto the shared `control_parameters`
!> struct would put a method's private options into every method's namespace, so
!> instead only what the compute path actually reads is packed into two small
!> fixed-layout arrays.  The authoritative schema is the parameter block below;
!> `include/oqp.h` documents it and `pyoqp/oqp/library/fci.py` mirrors it, and
!> `tests/test_fci_solve.py::test_option_schema_matches_fortran` parses this
!> file to assert the mirror is exact rather than merely intended.
!>
!> Everything the schema leaves out stays in Python on purpose: option parsing,
!> the validation that produces user-facing error messages, and the run
!> metadata.  All of it is integer bookkeeping executed once per solve, none of
!> it is on the compute path, and keeping it out means a native solve raises
!> exactly the messages the Python one always did.  Any negative status is a
!> "could not do it" signal, and the caller re-runs the Python `solve_fci`,
!> which is still the numerical pin.
!>
!> Array conventions follow the rest of the CI stack: everything crosses in the
!> caller's C order (last index fastest).  `h1e`/`eri` are the full MO-basis
!> integrals over `norb` orbitals; `civecs` comes back C-order [ndet, nroot], so
!> root k is `civecs[:, k]` exactly as `solve_fci` returns it.
module fci_driver_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use fci_hamiltonian_mod, only: fci_sort_dets, fci_dense_build, &
                                 fci_diag_build, fci_matvec_apply, oqp_dsyevd_f
  use fci_setup_mod, only: fci_spin_orbital_integrals, fci_fold_core, &
                           fci_spin_square
  ! DGEMM/DGEMV through the ILP64 wrapper layer (AGENTS.md rule 1).  `only:`
  ! cannot reach the renamed re-exports, so the whole module is used; on an
  ! ILP64 build (OQP_BLAS_INT /= 4) it is empty and the calls bind the raw
  ! 8-byte-integer symbols, both of which are registered in
  ! OQP_ACCELERATE_ILP64_SYMS.  The scalar reductions use Fortran intrinsics
  ! rather than DDOT/DNRM2 so no BLAS *function* has to be declared external in
  ! a build where the wrapper module contributes no interface.
  use oqp_linalg
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: fci_solve
  ! Reused by the CASSCF driver (casscf_driver.F90), which holds `h1e`/`eri`
  ! in local arrays and calls `fci_solve` directly -- exactly the in-process
  ! Fortran caller this entry point was shaped handle-free for.  It packs the
  ! same `iopt`/`dopt`, so it takes the schema from here rather than
  ! re-declaring index constants that could drift.
  public :: build_determinants
  public :: FCI_I_NORB, FCI_I_NACT, FCI_I_NCORE, FCI_I_NALPHA, FCI_I_NBETA
  public :: FCI_I_NROOT, FCI_I_SOLVER, FCI_I_MAXITER, FCI_I_SUBSPACE
  public :: FCI_I_MULT, FCI_I_MAXMEMORY, FCI_I_NTHREADS, FCI_I_WANT_S2
  public :: FCI_NIOPT, FCI_D_ECORE, FCI_D_EIG_TOL, FCI_D_CUTOFF, FCI_NDOPT
  public :: FCI_MAX_NSPIN

  ! ------------------------------------------------------------------ schema
  ! AUTHORITATIVE OPTION SCHEMA -- mirrored in include/oqp.h and in
  ! pyoqp/oqp/library/fci.py (_FCI_IOPT / _FCI_DOPT), and checked by
  ! tests/test_fci_solve.py.  Indices are 0-based, as seen from C.
  !
  ! iopt (int32):
  integer, parameter :: FCI_I_NORB       = 0   ! MOs spanned by h1e/eri
  integer, parameter :: FCI_I_NACT       = 1   ! active orbitals
  integer, parameter :: FCI_I_NCORE      = 2   ! frozen doubly-occupied orbitals
  integer, parameter :: FCI_I_NALPHA     = 3   ! active alpha electrons
  integer, parameter :: FCI_I_NBETA      = 4   ! active beta electrons
  integer, parameter :: FCI_I_NROOT      = 5   ! roots requested and returned
  integer, parameter :: FCI_I_SOLVER     = 6   ! 0 auto, 1 dense, 2 davidson
  integer, parameter :: FCI_I_MAXITER    = 7   ! Davidson iteration cap
  integer, parameter :: FCI_I_SUBSPACE   = 8   ! Davidson subspace cap, 0 = auto
  integer, parameter :: FCI_I_MULT       = 9   ! target multiplicity, 0 = any
  integer, parameter :: FCI_I_MAXMEMORY  = 10  ! working-set budget, MiB
  integer, parameter :: FCI_I_NTHREADS   = 11  ! OpenMP threads for the kernels
  integer, parameter :: FCI_I_WANT_S2    = 12  ! 1 = also return <S^2> per root
  integer, parameter :: FCI_NIOPT        = 13
  !
  ! dopt (double):
  integer, parameter :: FCI_D_ECORE      = 0   ! scalar added to every root
  integer, parameter :: FCI_D_EIG_TOL    = 1   ! eigenpair residual tolerance
  integer, parameter :: FCI_D_CUTOFF     = 2   ! integral screening cutoff
  integer, parameter :: FCI_NDOPT        = 3

  ! Status codes.  Success is the number of roots written (>= 0); every
  ! negative value means "the caller should fall back to the Python solver",
  ! which then raises the established user-facing message.
  integer(i8), parameter :: FCI_ERR_INPUT    = -1_i8  ! unsupported/invalid sizes
  integer(i8), parameter :: FCI_ERR_ALLOC    = -2_i8  ! out of memory
  integer(i8), parameter :: FCI_ERR_BUDGET   = -3_i8  ! over the max_memory budget
  integer(i8), parameter :: FCI_ERR_DAVIDSON = -4_i8  ! did not converge
  integer(i8), parameter :: FCI_ERR_SPIN     = -5_i8  ! no/too few roots of the target spin
  integer(i8), parameter :: FCI_ERR_EIGEN    = -6_i8  ! LAPACK failure or bad residual

  ! Determinants are packed into one 64-bit key (alpha bits 0..nact-1, beta
  ! bits nact..2*nact-1), so the active space is capped at 31 orbitals.
  integer, parameter :: FCI_MAX_NSPIN = 62

contains

  !> Complete CI solve: active space -> determinants -> eigenpairs -> roots.
  !>
  !> @param[in]  iopt     integer options, FCI_NIOPT entries (schema above)
  !> @param[in]  dopt     real options, FCI_NDOPT entries (schema above)
  !> @param[in]  active   active MO indices, 0-based, [nact]
  !> @param[in]  core     core MO indices, 0-based, [ncore] (may be a dummy
  !>                      single element when ncore == 0)
  !> @param[in]  h1e      full MO one-electron integrals, C-order [norb,norb]
  !> @param[in]  eri      full MO (pq|rs), C-order [norb,norb,norb,norb]
  !> @param[out] energies nroot root energies, already offset by dopt(ECORE)
  !> @param[out] civecs   CI vectors, C-order [ndet, nroot]
  !> @param[out] s2       <S^2> per returned root, written only when WANT_S2
  !> @return     the number of roots written, or a negative status
  function fci_solve(iopt, dopt, active, core, h1e, eri, &
                     energies, civecs, s2) result(status) bind(C, name="fci_solve")
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    integer(c_int32_t), intent(in) :: iopt(0:*)
    real(dp), intent(in) :: dopt(0:*)
    integer(c_int32_t), intent(in) :: active(0:*), core(0:*)
    real(dp), intent(in) :: h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: energies(0:*), civecs(0:*), s2(0:*)
    integer(i8) :: status

    integer :: norb, nact, ncore, na, nb, nroot, solver, maxiter, subspace
    integer :: mult, maxmem, nthreads, want_s2, nspin, ierr
    real(dp) :: ecore, eig_tol, cutoff
    integer(i8) :: ndet, dense_bytes, budget_bytes

    integer(i8), allocatable :: dets(:), skeys(:), sperm(:)
    real(dp), allocatable :: h_act(:), eri_act(:), hspin(:), gspin(:)
    real(dp), allocatable :: evals(:), evecs(:)
    real(dp), allocatable :: s2_win(:)
    integer, allocatable :: keep(:)
    real(dp) :: eref(1)
    integer :: nsel, k, r

    status = 0_i8

    norb     = int(iopt(FCI_I_NORB))
    nact     = int(iopt(FCI_I_NACT))
    ncore    = int(iopt(FCI_I_NCORE))
    na       = int(iopt(FCI_I_NALPHA))
    nb       = int(iopt(FCI_I_NBETA))
    nroot    = int(iopt(FCI_I_NROOT))
    solver   = int(iopt(FCI_I_SOLVER))
    maxiter  = int(iopt(FCI_I_MAXITER))
    subspace = int(iopt(FCI_I_SUBSPACE))
    mult     = int(iopt(FCI_I_MULT))
    maxmem   = int(iopt(FCI_I_MAXMEMORY))
    nthreads = int(iopt(FCI_I_NTHREADS))
    want_s2  = int(iopt(FCI_I_WANT_S2))
    ecore    = dopt(FCI_D_ECORE)
    eig_tol  = dopt(FCI_D_EIG_TOL)
    cutoff   = dopt(FCI_D_CUTOFF)

    nspin = 2 * nact
    if (nact <= 0 .or. norb <= 0 .or. nspin > FCI_MAX_NSPIN) then
      status = FCI_ERR_INPUT
      return
    end if
    if (na < 0 .or. nb < 0 .or. na > nact .or. nb > nact) then
      status = FCI_ERR_INPUT
      return
    end if
    if (nroot < 1 .or. maxiter < 1 .or. eig_tol <= 0.0_dp) then
      status = FCI_ERR_INPUT
      return
    end if

    ndet = binom(nact, na) * binom(nact, nb)
    if (ndet < 1_i8 .or. int(nroot, i8) > ndet) then
      status = FCI_ERR_INPUT
      return
    end if

    ! ---- determinants, active integrals, spin-orbital expansion
    allocate(dets(ndet), skeys(ndet), sperm(ndet), stat=ierr)
    if (ierr /= 0) then
      status = FCI_ERR_ALLOC
      return
    end if
    call build_determinants(nact, na, nb, ndet, dets)
    call fci_sort_dets(ndet, dets, skeys, sperm)

    allocate(h_act(0:int(nact, i8)**2 - 1_i8), &
             eri_act(0:int(nact, i8)**4 - 1_i8), &
             hspin(0:int(nspin, i8)**2 - 1_i8), &
             gspin(0:int(nspin, i8)**4 - 1_i8), stat=ierr)
    if (ierr /= 0) then
      status = FCI_ERR_ALLOC
      return
    end if

    call gather_active_eri(norb, nact, active, eri, eri_act)
    eref(1) = ecore
    call fci_fold_core(int(norb, c_int32_t), int(nact, c_int32_t), &
                       int(ncore, c_int32_t), active, core, h1e, eri, &
                       h_act, eref)
    ecore = eref(1)
    call fci_spin_orbital_integrals(int(nact, c_int32_t), h_act, eri_act, &
                                    hspin, gspin, int(nthreads, c_int32_t))

    ! ---- solver dispatch, matching solve_fci's "auto" rule exactly
    dense_bytes = 8_i8 * ndet * ndet
    budget_bytes = int(max(1, maxmem), i8) * 1024_i8 * 1024_i8
    if (solver == 0) then
      solver = merge(1, 2, dense_bytes <= budget_bytes)
    end if
    if (solver == 1 .and. dense_bytes > budget_bytes) then
      status = FCI_ERR_BUDGET
      return
    end if

    if (solver == 1) then
      call solve_dense(ndet, nspin, dets, skeys, sperm, hspin, gspin, cutoff, &
                       nroot, mult, eig_tol, maxiter, nthreads, &
                       evals, evecs, s2_win, keep, nsel, status)
    else
      call solve_iterative(ndet, nspin, dets, skeys, sperm, hspin, gspin, &
                           cutoff, nroot, mult, eig_tol, maxiter, subspace, &
                           budget_bytes, nthreads, &
                           evals, evecs, s2_win, keep, nsel, status)
    end if
    if (status < 0_i8) return
    if (nsel < nroot) then
      status = FCI_ERR_SPIN
      return
    end if

    ! ---- outputs.  `evecs` is column-major (ndet, *); the caller's buffer is
    ! C-order [ndet, nroot], so the write is an explicit transpose -- the two
    ! layouts coincide only for symmetric matrices, which this is not.
    do k = 1, nroot
      r = keep(k)
      energies(k - 1) = evals(r) + ecore
      call scatter_column(ndet, evecs, r, nroot, k, civecs)
    end do

    if (want_s2 /= 0) then
      if (allocated(s2_win)) then
        ! the spin diagnostics already ran on the root window; the per-root
        ! accumulation is independent, so the subset values are bit-identical
        do k = 1, nroot
          s2(k - 1) = s2_win(keep(k))
        end do
      else
        status = fci_spin_square(int(nact, c_int32_t), ndet, &
                                 int(nroot, c_int32_t), dets, civecs, &
                                 int(na - nb, c_int32_t), s2, &
                                 int(nthreads, c_int32_t))
        if (status /= 0_i8) then
          status = FCI_ERR_ALLOC
          return
        end if
      end if
    end if

    status = int(nroot, i8)
  end function fci_solve

  ! ==================================================================== dense
  !> Explicit ndet x ndet build followed by root extraction, mirroring
  !> solve_fci's dense branch: for `target_spin=any` the lowest `nroot` roots
  !> come from a Davidson pass over the already-built matrix (a full dsyevd
  !> produces all ndet eigenpairs to return a handful), and the full solve is
  !> the fallback and the path taken whenever a spin filter needs a window.
  subroutine solve_dense(ndet, nspin, dets, skeys, sperm, hspin, gspin, cutoff, &
                         nroot, mult, eig_tol, maxiter, nthreads, &
                         evals, evecs, s2_win, keep, nsel, status)
    integer(i8), intent(in) :: ndet
    integer, intent(in) :: nspin, nroot, mult, maxiter, nthreads
    integer(i8), intent(in) :: dets(*), skeys(*), sperm(*)
    real(dp), intent(in) :: hspin(*), gspin(*), cutoff, eig_tol
    real(dp), allocatable, intent(out) :: evals(:), evecs(:), s2_win(:)
    integer, allocatable, intent(out) :: keep(:)
    integer, intent(out) :: nsel
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: hmat(:), work(:)
    integer :: ierr, k, limit, nact
    logical :: have_roots, hsym
    real(dp) :: worst

    nact = nspin / 2
    allocate(hmat(ndet * ndet), stat=ierr)
    if (ierr /= 0) then
      status = FCI_ERR_ALLOC
      return
    end if
    call fci_dense_build(nspin, ndet, dets, skeys, sperm, hspin, gspin, &
                         cutoff, hmat, nthreads)
    ! fci_dense_build writes one symmetrized value into both triangles, so this
    ! holds by construction -- but the matvec below applies H^T in its place and
    ! would be silently wrong, not crash, if it ever stopped holding.  One
    ! streaming pass, against the build that just ran.
    hsym = is_symmetric(ndet, hmat)

    have_roots = .false.
    if (mult == 0) then
      ! Measured crossover of the whole dense path (see _lowest_dense_roots):
      ! below ndet=256, or when the requested window is a large fraction of
      ! the spectrum, the full solve already wins.
      if (ndet >= 256_i8 .and. int(8 * nroot, i8) < ndet) then
        allocate(evals(nroot), evecs(ndet * int(nroot, i8)), stat=ierr)
        if (ierr /= 0) then
          status = FCI_ERR_ALLOC
          return
        end if
        call ci_davidson(1, hmat, hsym, nspin, ndet, dets, skeys, sperm, &
                         hspin, gspin, cutoff, nroot, 0.1_dp * eig_tol, &
                         maxiter, 0, nthreads, evals, evecs, ierr)
        if (ierr == 0) then
          ! _davidson can also bail out early, so re-check rather than let the
          ! residual assertion below turn a solvable case into a hard failure
          worst = worst_residual(ndet, hmat, hsym, nthreads, nroot, evals, &
                                 evecs)
          have_roots = worst <= eig_tol
        end if
        if (.not. have_roots) deallocate(evals, evecs)
      end if
    end if

    if (.not. have_roots) then
      allocate(evals(ndet), evecs(ndet * ndet), stat=ierr)
      if (ierr /= 0) then
        status = FCI_ERR_ALLOC
        return
      end if
      ! dsyevd overwrites its input; the residual check below needs `hmat`
      evecs(1:ndet * ndet) = hmat(1:ndet * ndet)
      call oqp_dsyevd_f(int(ndet), evecs, evals, ierr)
      if (ierr /= 0) then
        status = FCI_ERR_EIGEN
        return
      end if
    end if

    if (mult == 0) then
      allocate(keep(nroot))
      do k = 1, nroot
        keep(k) = k
      end do
      nsel = nroot
    else
      ! Spin-classify only a low-root prefix and grow it on demand: the
      ! diagnostics cost is per-vector, and classifying all ndet eigenvectors
      ! of a large dense solve dominates a filter that consumes a handful.
      limit = int(min(ndet, int(max(4 * nroot, nroot + 12), i8)))
      do
        call classify_prefix(nact, ndet, dets, evecs, limit, nspin, nthreads, &
                             s2_win, ierr)
        if (ierr /= 0) then
          status = FCI_ERR_ALLOC
          return
        end if
        call select_by_multiplicity(s2_win, limit, mult, nroot, keep, nsel)
        if (nsel >= nroot) exit
        if (int(limit, i8) >= ndet) then
          status = FCI_ERR_SPIN
          return
        end if
        limit = int(min(ndet, 2_i8 * int(limit, i8)))
      end do
    end if

    ! Every returned root is verified against the matrix it came from, exactly
    ! as the Python driver does before handing roots back.
    worst = 0.0_dp
    do k = 1, nsel
      worst = max(worst, worst_residual(ndet, hmat, hsym, nthreads, 1, &
                                        evals(keep(k):keep(k)), &
                                        evecs((keep(k) - 1) * ndet + 1:)))
    end do
    if (worst > eig_tol) status = FCI_ERR_EIGEN
  end subroutine solve_dense

  ! ================================================================ iterative
  !> Matrix-free Davidson with the same root-window growth the Python driver
  !> uses when a spin filter is active.
  subroutine solve_iterative(ndet, nspin, dets, skeys, sperm, hspin, gspin, &
                             cutoff, nroot, mult, eig_tol, maxiter, subspace, &
                             budget_bytes, nthreads, &
                             evals, evecs, s2_win, keep, nsel, status)
    integer(i8), intent(in) :: ndet, budget_bytes
    integer, intent(in) :: nspin, nroot, mult, maxiter, subspace, nthreads
    integer(i8), intent(in) :: dets(*), skeys(*), sperm(*)
    real(dp), intent(in) :: hspin(*), gspin(*), cutoff, eig_tol
    real(dp), allocatable, intent(out) :: evals(:), evecs(:), s2_win(:)
    integer, allocatable, intent(out) :: keep(:)
    integer, intent(out) :: nsel
    integer(i8), intent(inout) :: status

    real(dp) :: dummy(1)
    integer :: solve_nroot, eff_sub, max_sub, ierr, k
    integer(i8) :: work_bytes
    integer :: nact

    nact = nspin / 2
    solve_nroot = nroot
    do
      eff_sub = 0
      if (subspace > 0) eff_sub = max(subspace, solve_nroot)
      max_sub = subspace_limit(ndet, solve_nroot, eff_sub)
      work_bytes = 8_i8 * ndet * int(2 * max_sub + 4 * solve_nroot, i8)
      if (work_bytes > budget_bytes) then
        status = FCI_ERR_BUDGET
        return
      end if

      if (allocated(evals)) deallocate(evals)
      if (allocated(evecs)) deallocate(evecs)
      allocate(evals(solve_nroot), evecs(ndet * int(solve_nroot, i8)), stat=ierr)
      if (ierr /= 0) then
        status = FCI_ERR_ALLOC
        return
      end if
      call ci_davidson(2, dummy, .false., nspin, ndet, dets, skeys, sperm, &
                       hspin, gspin, &
                       cutoff, solve_nroot, eig_tol, maxiter, eff_sub, &
                       nthreads, evals, evecs, ierr)
      if (ierr /= 0) then
        if (ierr == 2) then
          status = FCI_ERR_ALLOC
        else if (ierr == 3) then
          status = FCI_ERR_EIGEN
        else
          status = FCI_ERR_DAVIDSON
        end if
        return
      end if

      if (mult == 0) then
        allocate(keep(nroot))
        do k = 1, nroot
          keep(k) = k
        end do
        nsel = nroot
        return
      end if

      call classify_prefix(nact, ndet, dets, evecs, solve_nroot, nspin, &
                           nthreads, s2_win, ierr)
      if (ierr /= 0) then
        status = FCI_ERR_ALLOC
        return
      end if
      call select_by_multiplicity(s2_win, solve_nroot, mult, nroot, keep, nsel)
      if (nsel >= nroot) return
      if (int(solve_nroot, i8) >= ndet) then
        status = FCI_ERR_SPIN
        return
      end if
      solve_nroot = int(min(ndet, int(max(solve_nroot + 1, 2 * solve_nroot), i8)))
    end do
  end subroutine solve_iterative

  ! ================================================================= Davidson
  !> Block Davidson for the lowest `nroot` eigenpairs, reproducing
  !> fci.py::_davidson step for step.
  !>
  !> `mode` 1 applies an explicit dense matrix through DGEMM, `mode` 2 the
  !> matrix-free determinant enumeration; the sort of the determinant keys is
  !> hoisted out of the iteration, which the per-call bind(C) kernels could not
  !> do (they re-sorted on every application).
  !>
  !> The correction acceptance test is deliberately scale-free: a preconditioned
  !> correction is r/(theta - D), so its norm shrinks with the residual, and an
  !> absolute cutoff on it silently discards every correction once the residual
  !> gets small -- with eig_tol=1e-10 that stalled at a residual of 3.7e-8 and
  !> returned unconverged without saying so.  What matters is whether the
  !> direction survives orthogonalization, so the correction is normalized
  !> first and the test applied to what is left.
  !>
  !> ierr: 0 ok, 1 not converged, 2 allocation failure, 3 eigensolver failure.
  subroutine ci_davidson(mode, hmat, hsym, nspin, ndet, dets, skeys, sperm, &
                         hspin, gspin, cutoff, nroot, tol, max_iter, &
                         max_subspace_in, nthreads, theta_out, ritz_out, ierr)
    integer, intent(in) :: mode, nspin, nroot, max_iter, max_subspace_in, nthreads
    logical, intent(in) :: hsym
    integer(i8), intent(in) :: ndet
    real(dp), intent(in) :: hmat(*)
    integer(i8), intent(in) :: dets(*), skeys(*), sperm(*)
    real(dp), intent(in) :: hspin(*), gspin(*), cutoff, tol
    real(dp), intent(out) :: theta_out(*), ritz_out(*)
    integer, intent(out) :: ierr

    real(dp), allocatable :: basis(:,:), sigma(:,:), diag(:), sub(:), subv(:)
    real(dp), allocatable :: theta(:), ritz(:,:), resid(:,:), corr(:,:), vec(:)
    real(dp), allocatable :: proj(:), scratch(:), dense(:)
    real(dp), allocatable :: rnorm(:)
    integer, allocatable :: order(:)
    integer :: max_sub, nbas, ncur, iter, i, j, nadd, nnew, tiny_cap, info
    real(dp) :: scale, nrm, denom

    ierr = 0
    max_sub = subspace_limit(ndet, nroot, max_subspace_in)

    ! Tiny spaces: a direct dense diagonalization is cheaper and more robust
    ! than iterating, and it also covers nroot close to ndet.
    tiny_cap = max(2 * nroot + 20, 40)
    if (ndet <= int(tiny_cap, i8)) then
      allocate(dense(ndet * ndet), theta(ndet), stat=info)
      if (info /= 0) then
        ierr = 2
        return
      end if
      if (mode == 1) then
        dense(1:ndet * ndet) = hmat(1:ndet * ndet)
      else
        call fci_dense_build(nspin, ndet, dets, skeys, sperm, hspin, gspin, &
                             cutoff, dense, nthreads)
      end if
      call oqp_dsyevd_f(int(ndet), dense, theta, info)
      if (info /= 0) then
        ierr = 3
        return
      end if
      theta_out(1:nroot) = theta(1:nroot)
      ritz_out(1:ndet * int(nroot, i8)) = dense(1:ndet * int(nroot, i8))
      return
    end if

    allocate(diag(ndet), basis(ndet, max_sub), sigma(ndet, max_sub), &
             sub(max_sub * max_sub), subv(max_sub * max_sub), theta(max_sub), &
             ritz(ndet, nroot), resid(ndet, nroot), corr(ndet, nroot), &
             vec(ndet), proj(max_sub), rnorm(nroot), order(nroot), &
             scratch(2_i8 * ndet * int(max_sub, i8)), stat=info)
    if (info /= 0) then
      ierr = 2
      return
    end if

    call fci_diag_build(nspin, ndet, dets, hspin, gspin, diag)

    ! Initial guess: unit vectors on the nroot lowest diagonal elements.  They
    ! are already orthonormal, so the reference's QR of them is the identity up
    ! to column signs, which the subspace does not see.  Ties go to the lower
    ! index, matching a stable argsort.
    call lowest_indices(ndet, diag, nroot, order)
    basis(:, 1:nroot) = 0.0_dp
    do i = 1, nroot
      basis(order(i), i) = 1.0_dp
    end do
    nbas = nroot
    call apply_h(mode, hmat, hsym, nspin, ndet, dets, skeys, sperm, hspin, &
                 gspin, cutoff, nthreads, nbas, basis(:, 1:nbas), &
                 sigma(:, 1:nbas), scratch)

    do iter = 1, max_iter
      ! subspace matrix, symmetrized
      call dgemm('T', 'N', nbas, nbas, int(ndet), 1.0_dp, basis, int(ndet), &
                 sigma, int(ndet), 0.0_dp, sub, nbas)
      do j = 1, nbas
        do i = j + 1, nbas
          sub((j - 1) * nbas + i) = 0.5_dp * (sub((j - 1) * nbas + i) &
                                            + sub((i - 1) * nbas + j))
          sub((i - 1) * nbas + j) = sub((j - 1) * nbas + i)
        end do
      end do
      subv(1:nbas * nbas) = sub(1:nbas * nbas)
      call oqp_dsyevd_f(nbas, subv, theta, info)
      if (info /= 0) then
        ierr = 3
        return
      end if

      ! Ritz vectors and residuals for the lowest nroot
      call dgemm('N', 'N', int(ndet), nroot, nbas, 1.0_dp, basis, int(ndet), &
                 subv, nbas, 0.0_dp, ritz, int(ndet))
      call dgemm('N', 'N', int(ndet), nroot, nbas, 1.0_dp, sigma, int(ndet), &
                 subv, nbas, 0.0_dp, resid, int(ndet))
      do i = 1, nroot
        resid(:, i) = resid(:, i) - theta(i) * ritz(:, i)
        rnorm(i) = sqrt(sum(resid(:, i) * resid(:, i)))
      end do
      if (maxval(rnorm) < tol) then
        theta_out(1:nroot) = theta(1:nroot)
        call pack_columns(ndet, nroot, ritz, ritz_out)
        return
      end if

      ! preconditioned corrections for the unconverged roots
      nadd = 0
      do i = 1, nroot
        if (rnorm(i) < tol) cycle
        nadd = nadd + 1
        do j = 1, int(ndet)
          denom = theta(i) - diag(j)
          if (abs(denom) < 1.0e-12_dp) denom = 1.0e-12_dp
          corr(j, nadd) = resid(j, i) / denom
        end do
      end do

      ! restart from the current Ritz vectors when the subspace would overflow
      if (nadd == 0 .or. nbas + nadd > max_sub) then
        call orthonormalize(ndet, nroot, ritz)
        basis(:, 1:nroot) = ritz(:, 1:nroot)
        nbas = nroot
        call apply_h(mode, hmat, hsym, nspin, ndet, dets, skeys, sperm, hspin, &
                     gspin, cutoff, nthreads, nbas, basis(:, 1:nbas), &
                     sigma(:, 1:nbas), scratch)
        cycle
      end if

      ncur = nbas
      nnew = 0
      do i = 1, nadd
        vec(:) = corr(:, i)
        scale = sqrt(sum(vec * vec))
        if (scale <= 0.0_dp) cycle
        vec(:) = vec(:) / scale
        ! Orthogonalize twice for numerical stability, against the basis as it
        ! grows -- a new direction must also be independent of the ones already
        ! accepted in this sweep.
        do j = 1, 2
          call dgemv('T', int(ndet), ncur, 1.0_dp, basis, int(ndet), vec, 1, &
                     0.0_dp, proj, 1)
          call dgemv('N', int(ndet), ncur, -1.0_dp, basis, int(ndet), proj, 1, &
                     1.0_dp, vec, 1)
        end do
        nrm = sqrt(sum(vec * vec))
        if (nrm > 1.0e-6_dp) then
          ncur = ncur + 1
          basis(:, ncur) = vec(:) / nrm
          nnew = nnew + 1
        end if
      end do
      if (nnew == 0) then
        theta_out(1:nroot) = theta(1:nroot)
        call pack_columns(ndet, nroot, ritz, ritz_out)
        return
      end if
      call apply_h(mode, hmat, hsym, nspin, ndet, dets, skeys, sperm, hspin, &
                   gspin, cutoff, nthreads, nnew, basis(:, nbas + 1:ncur), &
                   sigma(:, nbas + 1:ncur), scratch)
      nbas = ncur
    end do

    ierr = 1
  end subroutine ci_davidson

  !> Y = 0.5(H + H^T) X for a block of `m` column-major vectors.
  !>
  !> `hsym` says the caller has VERIFIED (not assumed) that `hmat` is exactly
  !> symmetric; see is_symmetric and sym_matvec.
  subroutine apply_h(mode, hmat, hsym, nspin, ndet, dets, skeys, sperm, hspin, &
                     gspin, cutoff, nthreads, m, x, y, scratch)
    integer, intent(in) :: mode, nspin, nthreads, m
    logical, intent(in) :: hsym
    integer(i8), intent(in) :: ndet
    real(dp), intent(in) :: hmat(*)
    integer(i8), intent(in) :: dets(*), skeys(*), sperm(*)
    real(dp), intent(in) :: hspin(*), gspin(*), cutoff
    real(dp), intent(in) :: x(ndet, m)
    real(dp), intent(out) :: y(ndet, m)
    real(dp), intent(inout) :: scratch(*)
    integer(i8) :: j
    integer :: k

    if (mode == 1) then
      ! Up to two columns the transposed DGEMV (see sym_matvec) is the cheapest
      ! form: 2.26 ms per 4900^2 column against 6.10 ms for the DGEMM of the
      ! same shape, which reads `hmat` once no matter how many columns it is
      ! given and so only pays off from three columns on.
      if (hsym .and. m <= 2) then
        do k = 1, m
          call sym_matvec(ndet, hmat, nthreads, x(1, k), y(1, k))
        end do
      else if (m == 1) then
        ! A single-column DGEMM is a matrix-vector product spelled badly:
        ! Accelerate takes a much worse path for N=1 than DGEMV does.  Pass the
        ! element reference, not a rank-1 slice.
        call dgemv('N', int(ndet), int(ndet), 1.0_dp, hmat, int(ndet), &
                   x(1, 1), 1, 0.0_dp, y(1, 1), 1)
      else
        call dgemm('N', 'N', int(ndet), m, int(ndet), 1.0_dp, hmat, int(ndet), &
                   x, int(ndet), 0.0_dp, y, int(ndet))
      end if
      return
    end if
    ! The enumeration kernel indexes vector blocks in the caller's C order
    ! (vector index fastest), while the Davidson algebra wants column-major
    ! (determinant index fastest); repack across the call.  The repack is
    ! ndet*m doubles against an O(ndet * nspin^4) enumeration.
    do k = 1, m
      do j = 1, ndet
        scratch((j - 1) * m + k) = x(j, k)
      end do
    end do
    call fci_matvec_apply(nspin, ndet, dets, skeys, sperm, hspin, gspin, &
                          cutoff, m, scratch, scratch(ndet * int(m, i8) + 1), &
                          nthreads)
    do k = 1, m
      do j = 1, ndet
        y(j, k) = scratch(ndet * int(m, i8) + (j - 1) * m + k)
      end do
    end do
  end subroutine apply_h

  !> y = H x for a dense SYMMETRIC H held column-major, `nblk`-way threaded.
  !>
  !> The transposed DGEMV, not the plain one.  H = H^T so the two compute the
  !> same product, but not at the same speed: DGEMV('N') is the axpy form, which
  !> sweeps the whole of `y` once per column of H, while DGEMV('T') is the dot
  !> form, which streams each column into a single accumulator.  On a 4900^2
  !> matrix (192 MB) that is 7.46 ms against 3.25 ms -- 26 GB/s against 59 GB/s,
  !> i.e. the transposed form runs at the machine's single-core read ceiling and
  !> the plain one at less than half of it.  This is the whole of the dense
  !> Davidson's cost, and it is exactly where the Python driver's advantage came
  !> from: numpy's `A @ x` on a C-order matrix IS the transposed call.
  !>
  !> Splitting the COLUMNS across threads then buys the rest: each thread owns a
  !> disjoint slice of `y`, every y(j) stays one uninterrupted dot product over
  !> the full row range, and the read is spread over the cores instead of
  !> saturating one -- 2.26 ms, 85 GB/s at eight blocks.  Because no y(j) is ever
  !> split, the result is bit-identical to the serial call for every block count
  !> (verified 2..8 on a 1997^2 case: zero elements differ), so the thread count
  !> cannot move a printed energy.
  !>
  !> Accelerate does not thread DGEMV itself here -- the serial 59 GB/s is the
  !> same rate a single-threaded `sum()` reaches -- so this does not oversubscribe
  !> a BLAS that is already parallel.
  subroutine sym_matvec(ndet, hmat, nthreads, x, y)
    integer(i8), intent(in) :: ndet
    real(dp), intent(in) :: hmat(*)
    integer, intent(in) :: nthreads
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)

    integer :: nblk, ib
    integer(i8) :: lo, hi, chunk

    nblk = 1
    ! Below ~64k columns per block the OpenMP region costs more than the sweep
    ! it splits, so small matrices stay on the single call.
!$  nblk = max(1, min(max(1, nthreads), int(ndet / 512_i8)))

    if (nblk <= 1) then
      call dgemv('T', int(ndet), int(ndet), 1.0_dp, hmat, int(ndet), &
                 x, 1, 0.0_dp, y, 1)
      return
    end if

    chunk = (ndet + int(nblk, i8) - 1_i8) / int(nblk, i8)
    !$omp parallel do num_threads(nblk) schedule(static, 1) default(shared) &
    !$omp   private(ib, lo, hi)
    do ib = 0, nblk - 1
      lo = int(ib, i8) * chunk + 1_i8
      hi = min(ndet, lo + chunk - 1_i8)
      if (lo > hi) cycle
      ! element reference, never a rank-1 slice: the columns after `lo` have to
      ! stay contiguous for DGEMV to walk them
      call dgemv('T', int(ndet), int(hi - lo + 1_i8), 1.0_dp, &
                 hmat((lo - 1_i8) * ndet + 1_i8), int(ndet), &
                 x, 1, 0.0_dp, y(lo), 1)
    end do
    !$omp end parallel do
  end subroutine sym_matvec

  !> .true. when the dense Hamiltonian is symmetric to the last bit.
  !>
  !> `sym_matvec` applies H^T where the algorithm asks for H, which is the same
  !> operator only if H = H^T exactly.  `fci_dense_build` guarantees that -- it
  !> assigns one symmetrized expression into BOTH triangles -- but the failure
  !> mode if that ever stops holding is a wrong energy with no crash and no
  !> bounds violation, because every shape involved is square (the same trap
  !> that a transposed eigenvector operand set elsewhere on this branch).  So it
  !> is checked rather than inferred, once per dense solve.
  !>
  !> The scan is tiled because the naive nest reads one operand down a column
  !> and the other across a row; at 64 x 64 both stay in cache and the whole
  !> check is one streaming pass -- ~3 ms at ndet=4900, against the seconds the
  !> build it follows already spent.
  logical function is_symmetric(n, a) result(sym)
    integer(i8), intent(in) :: n
    real(dp), intent(in) :: a(*)
    integer(i8), parameter :: tile = 64_i8
    integer(i8) :: jb, ib, j, i, jhi, ihi

    sym = .true.
    do jb = 1_i8, n, tile
      jhi = min(n, jb + tile - 1_i8)
      do ib = jb, n, tile
        ihi = min(n, ib + tile - 1_i8)
        do j = jb, jhi
          do i = max(ib, j + 1_i8), ihi
            if (a((j - 1_i8) * n + i) /= a((i - 1_i8) * n + j)) then
              sym = .false.
              return
            end if
          end do
        end do
      end do
    end do
  end function is_symmetric

  ! ================================================================== helpers

  !> Determinant keys in CI order: alpha strings outer, beta strings inner,
  !> each string in itertools.combinations (lexicographic) order.
  subroutine build_determinants(norb, na, nb, ndet, dets)
    integer, intent(in) :: norb, na, nb
    integer(i8), intent(in) :: ndet
    integer(i8), intent(out) :: dets(ndet)
    integer(i8), allocatable :: astr(:), bstr(:)
    integer(i8) :: nas, nbs, ia, ib, k

    nas = binom(norb, na)
    nbs = binom(norb, nb)
    allocate(astr(nas), bstr(nbs))
    call string_basis(norb, na, nas, astr)
    call string_basis(norb, nb, nbs, bstr)
    k = 0_i8
    do ia = 1_i8, nas
      do ib = 1_i8, nbs
        k = k + 1_i8
        dets(k) = ior(astr(ia), shiftl(bstr(ib), norb))
      end do
    end do
    deallocate(astr, bstr)
  end subroutine build_determinants

  !> Occupation bit patterns for all C(norb, nelec) orbital subsets, in
  !> lexicographic (itertools.combinations) order.
  subroutine string_basis(norb, nelec, nstr, out)
    integer, intent(in) :: norb, nelec
    integer(i8), intent(in) :: nstr
    integer(i8), intent(out) :: out(nstr)
    integer :: c(max(nelec, 1)), i, j
    integer(i8) :: idx, key

    if (nelec == 0) then
      out(1) = 0_i8
      return
    end if
    do i = 1, nelec
      c(i) = i - 1
    end do
    idx = 0_i8
    do
      idx = idx + 1_i8
      key = 0_i8
      do i = 1, nelec
        key = ibset(key, c(i))
      end do
      out(idx) = key
      if (idx >= nstr) exit
      i = nelec
      do while (i >= 1)
        if (c(i) /= norb - nelec + i - 1) exit
        i = i - 1
      end do
      if (i < 1) exit
      c(i) = c(i) + 1
      do j = i + 1, nelec
        c(j) = c(j - 1) + 1
      end do
    end do
  end subroutine string_basis

  pure function binom(n, k) result(c)
    integer, intent(in) :: n, k
    integer(i8) :: c
    integer :: i
    c = 0_i8
    if (k < 0 .or. k > n) return
    c = 1_i8
    do i = 1, min(k, n - k)
      c = c * int(n - min(k, n - k) + i, i8) / int(i, i8)
    end do
  end function binom

  !> eri_act[p,q,r,s] = eri[active(p),active(q),active(r),active(s)], both in
  !> C order.  The innermost index is unit-stride in the destination.
  subroutine gather_active_eri(norb, nact, active, eri, eri_act)
    integer, intent(in) :: norb, nact
    integer(c_int32_t), intent(in) :: active(0:*)
    real(dp), intent(in) :: eri(0:*)
    real(dp), intent(out) :: eri_act(0:*)
    integer :: p, q, r, s
    integer(i8) :: n, n2, n3, m, m2, m3, src, dst

    n = int(norb, i8)
    n2 = n * n
    n3 = n2 * n
    m = int(nact, i8)
    m2 = m * m
    m3 = m2 * m
    do p = 0, nact - 1
      do q = 0, nact - 1
        do r = 0, nact - 1
          src = int(active(p), i8) * n3 + int(active(q), i8) * n2 &
                + int(active(r), i8) * n
          dst = int(p, i8) * m3 + int(q, i8) * m2 + int(r, i8) * m
          do s = 0, nact - 1
            eri_act(dst + int(s, i8)) = eri(src + int(active(s), i8))
          end do
        end do
      end do
    end do
  end subroutine gather_active_eri

  !> Davidson subspace cap, matching _davidson_subspace_limit.
  pure function subspace_limit(ndet, nroot, requested) result(lim)
    integer(i8), intent(in) :: ndet
    integer, intent(in) :: nroot, requested
    integer :: lim, want
    want = requested
    if (want <= 0) want = max(4 * nroot, 2 * nroot + 20)
    lim = int(min(ndet, int(want, i8)))
  end function subspace_limit

  !> Indices of the `k` smallest entries of `v`, ascending, ties to the lower
  !> index.  A selection rather than a full argsort: k is the root count, so
  !> this is O(k*n) against O(n log n) plus an n-sized permutation.
  subroutine lowest_indices(n, v, k, idx)
    integer(i8), intent(in) :: n
    real(dp), intent(in) :: v(n)
    integer, intent(in) :: k
    integer, intent(out) :: idx(k)
    logical, allocatable :: taken(:)
    integer(i8) :: i, best
    integer :: j

    allocate(taken(n))
    taken = .false.
    do j = 1, k
      best = 0_i8
      do i = 1_i8, n
        if (taken(i)) cycle
        if (best == 0_i8) then
          best = i
        else if (v(i) < v(best)) then
          best = i
        end if
      end do
      idx(j) = int(best)
      taken(best) = .true.
    end do
    deallocate(taken)
  end subroutine lowest_indices

  !> Modified Gram-Schmidt orthonormalization of `m` column-major vectors.
  subroutine orthonormalize(n, m, a)
    integer(i8), intent(in) :: n
    integer, intent(in) :: m
    real(dp), intent(inout) :: a(n, m)
    integer :: i, j
    real(dp) :: r, nrm
    do i = 1, m
      do j = 1, i - 1
        r = dot_product(a(:, j), a(:, i))
        a(:, i) = a(:, i) - r * a(:, j)
      end do
      nrm = sqrt(sum(a(:, i) * a(:, i)))
      if (nrm > 0.0_dp) a(:, i) = a(:, i) / nrm
    end do
  end subroutine orthonormalize

  !> max_k || A v_k - w_k v_k || over `m` column-major eigenpairs.
  function worst_residual(n, a, asym, nthreads, m, w, v) result(worst)
    integer(i8), intent(in) :: n
    integer, intent(in) :: m, nthreads
    logical, intent(in) :: asym
    real(dp), intent(in) :: a(*), w(*), v(*)
    real(dp) :: worst
    real(dp), allocatable :: r(:)
    integer :: k
    integer(i8) :: j

    allocate(r(n))
    worst = 0.0_dp
    do k = 1, m
      ! rank-1 slices must never be handed to BLAS as matrices; pass the
      ! element reference so the following columns stay contiguous
      if (asym) then
        call sym_matvec(n, a, nthreads, v((int(k, i8) - 1_i8) * n + 1_i8), r)
      else
        call dgemv('N', int(n), int(n), 1.0_dp, a, int(n), &
                   v((int(k, i8) - 1_i8) * n + 1_i8), 1, 0.0_dp, r, 1)
      end if
      do j = 1_i8, n
        r(j) = r(j) - w(k) * v((int(k, i8) - 1_i8) * n + j)
      end do
      worst = max(worst, sqrt(sum(r * r)))
    end do
    deallocate(r)
  end function worst_residual

  !> <S^2> for the first `limit` column-major eigenvectors.
  subroutine classify_prefix(nact, ndet, dets, evecs, limit, nspin, nthreads, &
                             s2_win, ierr)
    integer, intent(in) :: nact, limit, nspin, nthreads
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: evecs(*)
    real(dp), allocatable, intent(inout) :: s2_win(:)
    integer, intent(out) :: ierr
    real(dp), allocatable :: cbuf(:)
    integer(i8) :: j, rc
    integer :: k, na, nb

    ierr = 0
    if (allocated(s2_win)) deallocate(s2_win)
    allocate(s2_win(limit), cbuf(ndet * int(limit, i8)), stat=ierr)
    if (ierr /= 0) return
    ! fci_spin_square reads C-order [ndet, nvec]; the eigenvectors are
    ! column-major, so this is a real transpose, not a reinterpretation
    do k = 1, limit
      do j = 1_i8, ndet
        cbuf((j - 1_i8) * int(limit, i8) + int(k, i8)) = &
            evecs((int(k, i8) - 1_i8) * ndet + j)
      end do
    end do
    ! ms2 is carried by the determinant list itself; recover it from any
    ! determinant (all share the same alpha/beta electron counts)
    na = popcnt(iand(dets(1), shiftl(1_i8, nact) - 1_i8))
    nb = popcnt(shiftr(dets(1), nact))
    rc = fci_spin_square(int(nact, c_int32_t), ndet, int(limit, c_int32_t), &
                         dets, cbuf, int(na - nb, c_int32_t), s2_win, &
                         int(nthreads, c_int32_t))
    if (rc /= 0_i8) ierr = 1
    deallocate(cbuf)
  end subroutine classify_prefix

  !> Indices of the first `nroot` roots whose multiplicity matches `mult`.
  !> multiplicity = round(sqrt(1 + 4 <S^2>)), floored at 1, as in Python.
  subroutine select_by_multiplicity(s2_win, limit, mult, nroot, keep, nsel)
    real(dp), intent(in) :: s2_win(*)
    integer, intent(in) :: limit, mult, nroot
    integer, allocatable, intent(inout) :: keep(:)
    integer, intent(out) :: nsel
    integer :: k, m

    if (allocated(keep)) deallocate(keep)
    allocate(keep(nroot))
    nsel = 0
    do k = 1, limit
      m = max(1, nint(sqrt(max(0.0_dp, 1.0_dp + 4.0_dp * s2_win(k)))))
      if (m /= mult) cycle
      nsel = nsel + 1
      keep(nsel) = k
      if (nsel == nroot) exit
    end do
  end subroutine select_by_multiplicity

  !> Copy `m` column-major vectors of length n into a flat column-major buffer.
  subroutine pack_columns(n, m, src, dst)
    integer(i8), intent(in) :: n
    integer, intent(in) :: m
    real(dp), intent(in) :: src(n, m)
    real(dp), intent(out) :: dst(*)
    integer :: k
    integer(i8) :: j
    do k = 1, m
      do j = 1_i8, n
        dst((int(k, i8) - 1_i8) * n + j) = src(j, k)
      end do
    end do
  end subroutine pack_columns

  !> Write column-major eigenvector `col` into slot `slot` of a C-order
  !> [ndet, nroot] output buffer.
  subroutine scatter_column(ndet, evecs, col, nroot, slot, out)
    integer(i8), intent(in) :: ndet
    real(dp), intent(in) :: evecs(*)
    integer, intent(in) :: col, nroot, slot
    real(dp), intent(out) :: out(0:*)
    integer(i8) :: j, base
    base = (int(col, i8) - 1_i8) * ndet
    do j = 1_i8, ndet
      out((j - 1_i8) * int(nroot, i8) + int(slot, i8) - 1_i8) = evecs(base + j)
    end do
  end subroutine scatter_column

end module fci_driver_mod
