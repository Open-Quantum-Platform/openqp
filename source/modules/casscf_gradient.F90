!> @brief Analytic state-specific CASSCF nuclear gradient.
!>
!> This is the derivative of the CASSCF energy with respect to the NUCLEAR
!> COORDINATES.  It is a different object from `casscf_gfock_grad`
!> (`casscf_kernel.F90`), which is the derivative with respect to the ORBITAL
!> ROTATION parameters and is what the orbital optimizer drives to zero.  The
!> orbital-rotation gradient vanishing is a PRECONDITION for the expression
!> implemented here, not an approximation of it.
!>
!> Energy expression and conventions
!> ---------------------------------
!> Chemist notation throughout, `(pq|rs)` with the charge distributions `pq`
!> and `rs`.  Spin-summed density matrices,
!>
!>     D_pq       = sum_sigma <a^+_{p sigma} a_{q sigma}>
!>     Gamma_pqrs = sum_{sigma tau} <a^+_{p sigma} a^+_{r tau} a_{s tau} a_{q sigma}>
!>
!> so that
!>
!>     E = sum_pq D_pq h_pq + 1/2 sum_pqrs Gamma_pqrs (pq|rs) + V_NN.
!>
!> For CASSCF with inactive orbitals `i,j`, active `t,u,v,w` and virtual `a,b`,
!> `D` is 2 on the inactive diagonal and the CI 1-RDM `gamma_tu` on the active
!> block, and `Gamma` is the separable form
!>
!>     Gamma^sep_pqrs = D_pq D_rs - 1/2 D_ps D_rq
!>
!> everywhere EXCEPT the all-active block, which is the CI 2-RDM `Gamma_tuvw`.
!> These are exactly the conventions of `casscf.py::_full_rdms` and of the
!> generalized-Fock engine `casscf_kernel.F90`, whose
!>
!>     F_pq = sum_r D_pr h_qr + sum_rst Gamma_prst (qr|st)
!>
!> is reused here rather than re-derived.
!>
!> The gradient
!> ------------
!> With `x` a nuclear coordinate and the superscript `x` a derivative of the
!> AO integrals at fixed MO coefficients,
!>
!>     dE/dx = sum_{mu nu} D^AO_{mu nu} h^x_{mu nu}
!>           + 1/2 sum_{mu nu la si} Gamma^AO_{mu nu la si} (mu nu|la si)^x
!>           - sum_{mu nu} X^AO_{mu nu} S^x_{mu nu}
!>           + dV_NN/dx,
!>
!>     D^AO = C D C^T,   X^AO = C F C^T,   Gamma^AO = (C C C C) Gamma.
!>
!> There is NO orbital-response (Z-vector) term and no CI-response term.  That
!> is a consequence of full variationality, and it holds only under the
!> stationarity conditions checked below:
!>
!>   * the CI vector is an eigenvector of the active Hamiltonian at the current
!>     orbitals, so the energy is stationary with respect to every CI parameter
!>     (this is true of an excited root as well -- a stationary point need not
!>     be a minimum);
!>   * the orbital-rotation gradient vanishes over the non-redundant
!>     inactive-active, inactive-virtual and active-virtual pairs, which is what
!>     the optimizer converged;
!>   * rotations WITHIN the inactive, active and virtual blocks leave the energy
!>     invariant (for the active block because the CI space spans the full
!>     active CAS), so the energy is stationary there too.
!>
!> Together these make `F` symmetric, and `X^AO = C F C^T` is then the ordinary
!> energy-weighted density: for a closed-shell RHF reference `F_ij = 2 eps_i
!> delta_ij` and the term reduces to the familiar `-2 sum_i eps_i S^x_ii` that
!> `grd1::eijden` builds.  `|g_orb|` and `max|F - F^T|` are returned to the
!> caller so a run that was NOT converged cannot quietly produce a gradient
!> that looks plausible.
!>
!> This is the STATE-SPECIFIC gradient: the root the orbitals were optimized
!> for, i.e. `[casscf] root`.  It is valid only when the optimized objective is
!> that single root's energy.  A state-averaged run optimizes `sum_I w_I E_I`,
!> for which
!>
!>   * the gradient of the AVERAGED objective needs the weight-averaged RDMs
!>     (a small extension of what is here), and
!>   * the gradient of an INDIVIDUAL averaged state needs a coupled-perturbed /
!>     Z-vector solution, because the individual state's energy is NOT
!>     stationary with respect to the orbital rotations -- the averaged sum is.
!>
!> Neither is implemented, and a state-averaged run is REFUSED here rather than
!> silently given the state-specific formula.
!>
!> Two-particle density in the AO basis
!> ------------------------------------
!> The derivative-ERI driver `grd2` asks for a two-particle density block per
!> shell quartet and assumes it carries the full eight-fold permutational
!> symmetry of `(mu nu|la si)`.  Writing `Gamma^AO` as its separable part plus
!> an all-active correction,
!>
!>     Gamma^AO = Gamma^sep(D^AO) + dGamma,
!>     dGamma_{tuvw} = Gamma_tuvw - Gamma^sep_tuvw   (active MO indices only),
!>
!> the separable part is ALREADY what the Hartree-Fock kernel builds from a
!> density matrix -- symmetrized over `mu <-> nu`,
!>
!>     Gamma^sep = D_{mu nu} D_{la si}
!>               - 1/4 (D_{mu la} D_{nu si} + D_{nu la} D_{mu si}),
!>
!> which is `4 x` the Hartree-Fock `get_density` expression at `hfscale = 1`.
!> So the CASSCF two-particle density is the Hartree-Fock one evaluated at the
!> CASSCF AO density, plus a correction confined to the active space.
!>
!> The correction is NOT separable, and an nbf^4 AO tensor is neither necessary
!> nor affordable.  Symmetrized, `dGamma` is a symmetric matrix `M` over the
!> composite index `(tu)`, of dimension nact^2, so its eigendecomposition
!>
!>     M_{(tu),(vw)} = sum_k lambda_k V^k_{tu} V^k_{vw}
!>
!> turns the correction into a short sum of SEPARABLE products of AO matrices,
!>
!>     dGamma_{mu nu la si} = sum_k lambda_k A^k_{mu nu} A^k_{la si},
!>     A^k = C_act V^k C_act^T,
!>
!> with at most nact(nact+1)/2 nonzero terms (the antisymmetric half of the
!> composite space is annihilated by the symmetrization).  Each `A^k` is an
!> ordinary symmetric AO matrix, so it goes through the same bfnrm folding and
!> spherical->Cartesian expansion (`build_cart_density`) the density does, and
!> each term is manifestly eight-fold symmetric.  Storage is
!> nact(nact+1)/2 x nbf_cart^2 instead of nbf_cart^4.
module casscf_gradient_mod

  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use precision, only: dp
  use types, only: information
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use grd2, only: grd2_driver, grd2_compute_data_t
  use casscf_driver_mod, only: cas_ctx_t, cas_context_init, cas_evaluate
  use fci_hamiltonian_mod, only: oqp_dsyevd_f
  ! DGEMM through the ILP64 wrapper layer (AGENTS.md rule 1); `only:` cannot be
  ! used here for the same reason as in fci_driver.F90.
  use oqp_linalg

  implicit none
  private

  integer, parameter :: i8 = c_int64_t

  character(len=*), parameter :: module_name = "casscf_gradient_mod"

  public :: casscf_gradient

  !> `info` slots handed back to the caller (all double):
  integer, parameter :: CAS_G_GNORM  = 0  !< |g_orb|, the orbital-rotation gradient
  integer, parameter :: CAS_G_ENERGY = 1  !< energy of the differentiated root
  integer, parameter :: CAS_G_FASYM  = 2  !< max |F_pq - F_qp|, stationarity probe
  integer, parameter :: CAS_G_NVEC   = 3  !< active-correction factorization rank
  integer, parameter :: CAS_G_NINFO  = 4

  ! Status codes.  0 is success.  Codes at or below -20 are this module's own;
  ! anything between -1 and -9 is propagated unchanged from `cas_context_init`
  ! (see the CAS_ERR_* block in casscf_driver.F90) so the caller can report the
  ! same cause it would report for an energy run.
  integer(i8), parameter :: CASG_ERR_STATE_AVERAGE = -20_i8  !< nstate /= 1
  integer(i8), parameter :: CASG_ERR_DFT           = -21_i8  !< non-HF Hamiltonian
  integer(i8), parameter :: CASG_ERR_ALLOC         = -22_i8  !< out of memory
  integer(i8), parameter :: CASG_ERR_EIGEN         = -23_i8  !< LAPACK failure
  integer(i8), parameter :: CASG_ERR_EVALUATE      = -24_i8  !< CI / Fock build failed

  !> Two-particle density for the derivative-ERI driver.
  !>
  !> `d2` is the CASSCF spin-summed AO one-particle density; `av(k,:,:)` and
  !> `lam(k)` are the factorized all-active correction described in the module
  !> header.  The `_cart` copies are the bfnrm-folded, Cartesian-expanded
  !> versions the derivative kernels contract under HARMONIC_ACTIVE, built by
  !> `init` exactly as the Hartree-Fock path builds its own.
  type, extends(grd2_compute_data_t) :: grd2_cas_compute_data_t
    real(kind=dp), allocatable :: d2(:,:)
    real(kind=dp), allocatable :: av(:,:,:)
    real(kind=dp), allocatable :: lam(:)
    real(kind=dp), allocatable :: d2_cart(:,:)
    real(kind=dp), allocatable :: av_cart(:,:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
    integer :: nbf_cart = 0
    integer :: nvec = 0
  contains
    procedure :: init => grd2_cas_init
    procedure :: clean => grd2_cas_clean
    procedure :: get_density => grd2_cas_get_density
    procedure :: build_cart => grd2_cas_build_cart
  end type

contains

!###############################################################################

  !> C-bound entry point.  One call produces the complete nuclear gradient in
  !> `infos%atoms%grad`.
  !>
  !> The option block is the SAME `iopt`/`dopt` schema `casscf_energy` takes,
  !> so a gradient run is configured by the arrays that configured the energy
  !> run it differentiates; only the optimizer-specific entries are ignored.
  !>
  !> @param[in]     c_handle  the OQP handle (`OQP::Hcore`, `OQP::AO_ERI` and
  !>                          `OQP::VEC_MO_A` are read; `infos%atoms%grad` is
  !>                          overwritten)
  !> @param[in]     iopt      integer options, CAS_NIOPT entries
  !> @param[in]     dopt      real options, CAS_NDOPT entries
  !> @param[in]     weights   state weights; must be the single entry 1.0
  !> @param[in]     roots     the 0-based CI root to differentiate, one entry
  !> @param[out]    info      CAS_G_NINFO diagnostics (see the slot constants)
  !> @return        0, or a negative status
  function casscf_gradient_C(c_handle, iopt, dopt, weights, roots, info) &
      result(status) bind(C, name="casscf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    integer(c_int32_t), intent(in) :: iopt(0:*)
    real(dp), intent(in) :: dopt(0:*)
    real(dp), intent(in) :: weights(0:*)
    integer(c_int32_t), intent(in) :: roots(0:*)
    real(dp), intent(inout) :: info(0:*)
    integer(i8) :: status
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    status = casscf_gradient(inf, iopt, dopt, weights, roots, info)
  end function casscf_gradient_C

!###############################################################################

  !> The run itself, reachable from an in-process Fortran caller.
  function casscf_gradient(infos, iopt, dopt, weights, roots, info) &
      result(status)
    use io_constants, only: iw
    use grd1, only: print_gradient
    use util, only: measure_time
    use printing, only: print_module_info

    type(information), target, intent(inout) :: infos
    integer(c_int32_t), intent(in) :: iopt(0:*)
    real(dp), intent(in) :: dopt(0:*)
    real(dp), intent(in) :: weights(0:*)
    integer(c_int32_t), intent(in) :: roots(0:*)
    real(dp), intent(inout) :: info(0:*)
    integer(i8) :: status

    type(cas_ctx_t), target :: ctx
    type(basis_set), pointer :: basis
    real(dp), allocatable :: cbuf(:,:), gorb(:)
    real(dp), allocatable :: dmo(:,:), fmo(:,:), dao(:,:), xao(:,:), scr(:,:)
    real(dp), allocatable :: dm_packed(:), wm_packed(:)
    type(grd2_cas_compute_data_t) :: gcomp
    integer :: n, nc, na, k, p, q, nbf_tri, ierr
    integer(i8) :: rc
    real(dp) :: objective, asym

    status = 0_i8
    do k = 0, CAS_G_NINFO - 1
      info(k) = 0.0_dp
    end do

    ! ---- context: the same one the optimizer converged in
    call cas_context_init(infos, iopt, dopt, weights, roots, ctx, cbuf, status)
    if (status < 0_i8) return

    ! A state-averaged objective is a different derivative; refuse rather than
    ! hand back the state-specific formula evaluated on averaged quantities.
    if (ctx%nstate /= 1) then
      status = CASG_ERR_STATE_AVERAGE
      return
    end if
    ! The energy expression above has no exchange-correlation term, and CASSCF
    ! requires Hartree-Fock integrals anyway.
    if (infos%control%hamilton >= 20) then
      status = CASG_ERR_DFT
      return
    end if

    n  = ctx%n
    nc = ctx%nc
    na = ctx%na
    nbf_tri = n * (n + 1) / 2

    open (unit=IW, file=infos%log_filename, position="append")
    call print_module_info('CASSCF_Gradient', &
                           'Analytic State-Specific CASSCF Nuclear Gradient')

    ! ---- one energy/gradient evaluation at the converged orbitals: CI, the
    ! per-root RDMs, the generalized Fock and the orbital-rotation gradient
    allocate(gorb(0:ctx%npar-1), stat=ierr)
    if (ierr /= 0) then
      status = CASG_ERR_ALLOC
      close(iw)
      return
    end if
    call cas_evaluate(ctx, cbuf, objective, gorb, rc)
    if (rc < 0_i8) then
      status = CASG_ERR_EVALUATE
      close(iw)
      return
    end if

    info(CAS_G_GNORM) = sqrt(sum(gorb * gorb))
    info(CAS_G_ENERGY) = objective

    ! ---- MO-basis one-particle density and the (symmetrized) generalized Fock
    allocate(dmo(0:n-1, 0:n-1), fmo(0:n-1, 0:n-1), dao(0:n-1, 0:n-1), &
             xao(0:n-1, 0:n-1), scr(0:n-1, 0:n-1), stat=ierr)
    if (ierr /= 0) then
      status = CASG_ERR_ALLOC
      close(iw)
      return
    end if

    dmo = 0.0_dp
    do k = 0, nc - 1
      dmo(k, k) = 2.0_dp
    end do
    do p = 0, na - 1
      do q = 0, na - 1
        dmo(nc + p, nc + q) = ctx%g1(int(p, i8)*int(na, i8) + int(q, i8))
      end do
    end do

    ! `ctx%fock` is the C-order buffer F[p][q].  At a stationary point F is
    ! symmetric; the residual asymmetry is reported and then averaged away, so
    ! the energy-weighted density is symmetric by construction (the packed
    ! representation the overlap-derivative kernel takes requires it).
    asym = 0.0_dp
    do p = 0, n - 1
      do q = 0, n - 1
        asym = max(asym, abs(ctx%fock(int(p, i8)*int(n, i8) + int(q, i8)) &
                           - ctx%fock(int(q, i8)*int(n, i8) + int(p, i8))))
        fmo(p, q) = 0.5_dp * (ctx%fock(int(p, i8)*int(n, i8) + int(q, i8)) &
                            + ctx%fock(int(q, i8)*int(n, i8) + int(p, i8)))
      end do
    end do
    info(CAS_G_FASYM) = asym

    ! ---- AO transforms.  `cbuf(p, u)` is C_{u p} (MO index first), so
    ! M^AO = C M C^T is  (C^T)(M C)  with the transpose on the left factor.
    call dgemm('N', 'N', n, n, n, 1.0_dp, dmo, n, cbuf, n, 0.0_dp, scr, n)
    call dgemm('T', 'N', n, n, n, 1.0_dp, cbuf, n, scr, n, 0.0_dp, dao, n)
    call dgemm('N', 'N', n, n, n, 1.0_dp, fmo, n, cbuf, n, 0.0_dp, scr, n)
    call dgemm('T', 'N', n, n, n, 1.0_dp, cbuf, n, scr, n, 0.0_dp, xao, n)

    allocate(dm_packed(nbf_tri), wm_packed(nbf_tri), stat=ierr)
    if (ierr /= 0) then
      status = CASG_ERR_ALLOC
      close(iw)
      return
    end if
    call pack_symmetric(n, dao, dm_packed)
    ! The overlap-derivative kernel ADDS  sum_{mu nu} W_{mu nu} S^x_{mu nu},
    ! and the gradient term is  -sum X^AO S^x, so the negated Lagrangian is
    ! what goes in -- the same sign convention `grd1::eijden` uses for RHF.
    call pack_symmetric(n, -xao, wm_packed)

    ! ---- AO one-particle density for the two-electron contraction (1-based,
    ! the indexing `grd2` uses through `basis%ao_offset`)
    allocate(gcomp%d2(n, n), stat=ierr)
    if (ierr /= 0) then
      status = CASG_ERR_ALLOC
      close(iw)
      return
    end if
    do p = 0, n - 1
      do q = 0, n - 1
        gcomp%d2(p + 1, q + 1) = dao(p, q)
      end do
    end do

    ! ---- factorized all-active two-particle correction
    call build_active_correction(ctx, cbuf, gcomp, status)
    if (status < 0_i8) then
      close(iw)
      return
    end if
    info(CAS_G_NVEC) = real(gcomp%nvec, dp)

    basis => infos%basis
    basis%atoms => infos%atoms

    write(iw,"(/' ..... Beginning CASSCF Gradient Calculation...'/)")
    write(iw,"(' Differentiated root                 ',I8)") ctx%roots(0)
    write(iw,"(' Orbital-rotation gradient norm      ',ES16.8)") info(CAS_G_GNORM)
    write(iw,"(' Generalized Fock asymmetry          ',ES16.8)") info(CAS_G_FASYM)
    write(iw,"(' Active 2-RDM correction vectors     ',I8)") gcomp%nvec
    call flush(iw)

    ! ---- one-electron, nuclear-repulsion and overlap (Pulay) terms
    call cas_1e_grad(infos, basis, dm_packed, wm_packed)
    write(iw,"(' ..... End Of 1-Electron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    ! ---- two-electron term
    call cas_2e_grad(infos, basis, gcomp, status)
    if (status < 0_i8) then
      close(iw)
      return
    end if
    write(iw, fmt="(' ...... End Of 2-Electron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    call print_gradient(infos)
    call measure_time(print_total=1, log_unit=iw)
    close(iw)

    call gcomp%clean()
  end function casscf_gradient

!###############################################################################

  !> Pack the upper triangle of a symmetric matrix in the layout the gradient
  !> kernels read (`grd1::prepare_grad_density` unpacks it with the same
  !> default convention `mathlib::pack_matrix` writes).
  subroutine pack_symmetric(n, a, ap)
    use mathlib, only: pack_matrix
    integer, intent(in) :: n
    real(dp), intent(in) :: a(0:n-1, 0:n-1)
    real(dp), intent(out) :: ap(:)
    real(dp), allocatable :: sym(:,:)
    integer :: i, j

    allocate(sym(n, n))
    do i = 0, n - 1
      do j = 0, n - 1
        sym(i + 1, j + 1) = 0.5_dp * (a(i, j) + a(j, i))
      end do
    end do
    call pack_matrix(sym, ap)
    deallocate(sym)
  end subroutine pack_symmetric

!###############################################################################

  !> Build the factorized all-active two-particle correction and load the
  !> compute object with everything `grd2` will contract.
  !>
  !> `dGamma_{tuvw} = Gamma_tuvw - (gamma_tu gamma_vw - 1/2 gamma_tw gamma_vu)`
  !> is symmetrized over the full eight-element permutation group of
  !> `(tu|vw)`, which makes it a symmetric nact^2 x nact^2 matrix; its
  !> eigendecomposition gives the separable terms
  !> `dGamma = sum_k lambda_k A^k (x) A^k` with `A^k = C_act V^k C_act^T`.
  !>
  !> The symmetrization annihilates the antisymmetric half of the composite
  !> index, so nact(nact-1)/2 eigenvalues are numerically zero and are dropped
  !> against a relative threshold; the remaining rank is at most
  !> nact(nact+1)/2.
  subroutine build_active_correction(ctx, cbuf, gcomp, status)
    type(cas_ctx_t), intent(in) :: ctx
    real(dp), contiguous, intent(in) :: cbuf(0:,0:)
    type(grd2_cas_compute_data_t), intent(inout) :: gcomp
    integer(i8), intent(inout) :: status

    real(dp), parameter :: DROP_REL = 1.0e-14_dp

    real(dp), allocatable :: dg(:,:,:,:), m2(:,:), eval(:), vmat(:,:)
    real(dp), allocatable :: cact(:,:), amat(:,:), tmpa(:,:)
    integer :: n, nc, na, na2, t, u, v, w, k, nkeep, ierr, idx
    integer(i8) :: na_i8
    real(dp) :: gsep, thr, big

    n = ctx%n
    nc = ctx%nc
    na = ctx%na
    na2 = na * na
    na_i8 = int(na, i8)

    allocate(dg(0:na-1, 0:na-1, 0:na-1, 0:na-1), m2(na2, na2), eval(na2), &
             stat=ierr)
    if (ierr /= 0) then
      status = CASG_ERR_ALLOC
      return
    end if

    ! dGamma in the CI's own index order (`ctx%g2` is C-order [na,na,na,na],
    ! chemist (tu|vw); `ctx%g1` is C-order [na,na]).
    do t = 0, na - 1
      do u = 0, na - 1
        do v = 0, na - 1
          do w = 0, na - 1
            gsep = ctx%g1(int(t, i8)*na_i8 + int(u, i8)) &
                 * ctx%g1(int(v, i8)*na_i8 + int(w, i8)) &
                 - 0.5_dp * ctx%g1(int(t, i8)*na_i8 + int(w, i8)) &
                          * ctx%g1(int(v, i8)*na_i8 + int(u, i8))
            dg(t, u, v, w) = ctx%g2((((int(t, i8)*na_i8 + int(u, i8))*na_i8 &
                                      + int(v, i8))*na_i8 + int(w, i8))) - gsep
          end do
        end do
      end do
    end do

    ! Symmetrize over the eight permutations the ERI has.  The two-particle
    ! density only carries four of them on its own (pair exchange and the
    ! simultaneous within-pair swap), and `grd2` contracts one representative
    ! per orbit with a multiplicity factor, so the missing ones must be
    ! averaged in explicitly or the correction would be double- or
    ! under-counted for mixed quartets.
    do t = 0, na - 1
      do u = 0, na - 1
        do v = 0, na - 1
          do w = 0, na - 1
            m2(t*na + u + 1, v*na + w + 1) = 0.125_dp * &
                (dg(t, u, v, w) + dg(u, t, v, w) + dg(t, u, w, v) &
               + dg(u, t, w, v) + dg(v, w, t, u) + dg(w, v, t, u) &
               + dg(v, w, u, t) + dg(w, v, u, t))
          end do
        end do
      end do
    end do
    deallocate(dg)

    ! Columns of `m2` come back as the eigenvectors, ascending in `eval`.
    call oqp_dsyevd_f(na2, m2, eval, ierr)
    if (ierr /= 0) then
      status = CASG_ERR_EIGEN
      return
    end if

    big = 0.0_dp
    do k = 1, na2
      big = max(big, abs(eval(k)))
    end do
    thr = DROP_REL * max(big, 1.0_dp)
    nkeep = 0
    do k = 1, na2
      if (abs(eval(k)) > thr) nkeep = nkeep + 1
    end do

    gcomp%nbf = n
    gcomp%nvec = nkeep
    allocate(gcomp%lam(max(nkeep, 1)), gcomp%av(max(nkeep, 1), n, n), &
             cact(n, na), vmat(na, na), tmpa(n, na), amat(n, n), stat=ierr)
    if (ierr /= 0) then
      status = CASG_ERR_ALLOC
      return
    end if

    ! Active MO coefficients in the AO basis: cact(u, t) = C_{u, nc+t}.
    do t = 0, na - 1
      do u = 0, n - 1
        cact(u + 1, t + 1) = cbuf(nc + t, u)
      end do
    end do

    idx = 0
    do k = 1, na2
      if (abs(eval(k)) <= thr) cycle
      idx = idx + 1
      gcomp%lam(idx) = eval(k)
      do t = 0, na - 1
        do u = 0, na - 1
          vmat(t + 1, u + 1) = m2(t*na + u + 1, k)
        end do
      end do
      ! A^k = C_act V^k C_act^T
      call dgemm('N', 'N', n, na, na, 1.0_dp, cact, n, vmat, na, 0.0_dp, tmpa, n)
      call dgemm('N', 'T', n, n, na, 1.0_dp, tmpa, n, cact, n, 0.0_dp, amat, n)
      gcomp%av(idx, :, :) = amat
    end do

    deallocate(m2, eval, cact, vmat, tmpa, amat)
  end subroutine build_active_correction

!###############################################################################

  !> Nuclear-repulsion, overlap (Pulay), kinetic, nuclear-attraction and ECP
  !> terms.  Structurally identical to `hf_gradient_mod::hf_1e_grad`; the only
  !> difference is where the two densities come from -- the CASSCF spin-summed
  !> AO density and the CASSCF energy-weighted density replace the SCF ones.
  subroutine cas_1e_grad(infos, basis, dm_packed, wm_packed)
    use constants, only: tol_int
    use grd1, only: grad_nn, grad_ee_overlap, grad_ee_kinetic, &
                    grad_en_hellman_feynman, grad_en_pulay, grad_1e_ecp

    type(information), intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(dp), intent(inout) :: dm_packed(:), wm_packed(:)

    real(dp) :: tol

    tol = tol_int * log(10.0_dp)

    associate(grad => infos%atoms%grad, &
              xyz  => infos%atoms%xyz, &
              zn   => infos%atoms%zn - infos%basis%ecp_zn_num)

      grad = 0.0_dp

      call grad_nn(infos%atoms, infos%basis%ecp_zn_num)

      ! Energy-weighted density (already negated by the caller)
      call grad_ee_overlap(basis, wm_packed, grad, logtol=tol)

      call grad_ee_kinetic(basis, dm_packed, grad, logtol=tol)
      call grad_en_hellman_feynman(basis, xyz, zn, dm_packed, grad, logtol=tol)
      call grad_en_pulay(basis, xyz, zn, dm_packed, grad, logtol=tol)
      call grad_1e_ecp(infos, basis, xyz, dm_packed, grad, logtol=tol)

    end associate
  end subroutine cas_1e_grad

!###############################################################################

  !> Two-electron term: contract the CASSCF two-particle density against the
  !> derivative ERIs.
  !>
  !> The symmetry petite list is deliberately NOT opted into.  It is valid only
  !> for a totally symmetric contracted density, and the state-specific CASSCF
  !> two-particle density of an arbitrary `[casscf] root` carries no such
  !> guarantee -- a spatially degenerate root in an abelian subgroup is the
  !> obvious counterexample.  The cost of declining is speed; the cost of
  !> opting in wrongly is a wrong gradient.
  subroutine cas_2e_grad(infos, basis, gcomp, status)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    type(grd2_cas_compute_data_t), intent(inout) :: gcomp
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: de(:,:)
    integer :: ierr

    allocate(de(3, ubound(infos%atoms%zn, 1)), source=0.0_dp, stat=ierr)
    if (ierr /= 0) then
      status = CASG_ERR_ALLOC
      return
    end if

    gcomp%hfscale = 1.0_dp
    gcomp%coulscale = 1.0_dp
    call gcomp%init()
    call gcomp%build_cart(basis)
    ! `cam=.false.` explicitly rather than by inheritance. The range-separated
    ! path runs `grd2_driver_gen` TWICE with different Coulomb/exchange scales,
    ! and the all-active correction is neither Coulomb nor exchange -- it would
    ! be added unscaled in both passes. CASSCF requires Hartree-Fock integrals
    ! (the caller already refused a DFT Hamiltonian), so this only pins that.
    call grd2_driver(infos, basis, de, gcomp, cam=.false.)
    infos%atoms%grad = infos%atoms%grad + de

    deallocate(de)
  end subroutine cas_2e_grad

!###############################################################################

  !> Fold the basis normalization into the density and the correction factors
  !> and, under HARMONIC_ACTIVE, expand both to the Cartesian-effective basis
  !> the derivative kernels iterate.  Each `A^k` is an ordinary symmetric AO
  !> matrix, so it takes exactly the same route the density does.
  subroutine grd2_cas_init(this)
    class(grd2_cas_compute_data_t), target, intent(inout) :: this
    this%nbf_cart = 0
  end subroutine grd2_cas_init

  !> The Cartesian build needs the basis, which `init` does not receive; it is
  !> driven from `cas_2e_grad` through this entry instead.
  subroutine grd2_cas_build_cart(this, basis)
    class(grd2_cas_compute_data_t), intent(inout) :: this
    type(basis_set), intent(in) :: basis

    real(dp), allocatable :: tmp(:,:), expanded(:,:)
    integer, allocatable :: dummy_off(:)
    integer :: k, nc_tmp

    if (.not. HARMONIC_ACTIVE) return

    tmp = this%d2
    call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
    call build_cart_density(basis, tmp, this%d2_cart, this%cart_off, &
                            this%nbf_cart)

    if (this%nvec > 0) then
      allocate(this%av_cart(this%nvec, this%nbf_cart, this%nbf_cart))
      do k = 1, this%nvec
        tmp = this%av(k, :, :)
        call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
        call build_cart_density(basis, tmp, expanded, dummy_off, nc_tmp)
        this%av_cart(k, :, :) = expanded
        deallocate(expanded, dummy_off)
      end do
    end if
  end subroutine grd2_cas_build_cart

  subroutine grd2_cas_clean(this)
    class(grd2_cas_compute_data_t), target, intent(inout) :: this
    if (allocated(this%d2)) deallocate(this%d2)
    if (allocated(this%av)) deallocate(this%av)
    if (allocated(this%lam)) deallocate(this%lam)
    if (allocated(this%d2_cart)) deallocate(this%d2_cart)
    if (allocated(this%av_cart)) deallocate(this%av_cart)
    if (allocated(this%cart_off)) deallocate(this%cart_off)
  end subroutine grd2_cas_clean

!###############################################################################

  !> Two-particle density for one shell quartet.
  !>
  !> The separable part is the Hartree-Fock expression evaluated at the CASSCF
  !> AO density -- it is `4 x` the symmetrized `Gamma^sep` of the module header
  !> and reduces to `grd2_rhf_compute_data_t_get_density` term for term.  The
  !> all-active correction is added as `4 sum_k lambda_k A^k_{ij} A^k_{kl}`
  !> before the normalization factor, because both pieces are four-index
  !> products of AO-basis objects and take the same `bfnrm` weight.
  subroutine grd2_cas_get_density(this, basis, id, dab, dabmax)

    class(grd2_cas_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: coulfact, xcfact, df1, dq1, corr, bfn
    logical :: do_exchange, usecart
    integer :: i, j, k, l, v
    integer :: loc(4), nbf(4)
    integer :: i1, j1, k1, l1
    real(kind=dp), pointer :: ab(:,:,:,:)
    real(kind=dp), pointer :: d2(:,:), av(:,:,:)

    coulfact = 4 * this%coulscale
    xcfact = this%hfscale
    do_exchange = xcfact /= 0.0_dp

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      d2 => this%d2_cart
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
    else
      d2 => this%d2
      loc = basis%ao_offset(id) - 1
      nbf = basis%naos(id)
    end if
    av => null()
    if (this%nvec > 0) then
      if (usecart) then
        av => this%av_cart
      else
        av => this%av
      end if
    end if

    dabmax = 0
    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i

      do j = 1, nbf(2)
        j1 = loc(2) + j

        do k = 1, nbf(3)
          k1 = loc(3) + k

          do l = 1, nbf(4)
            l1 = loc(4) + l
            df1 = coulfact*d2(i1,j1)*d2(k1,l1)
            if (do_exchange) then
              dq1 = d2(i1,k1)*d2(j1,l1) &
                  + d2(i1,l1)*d2(j1,k1)
              df1 = df1-xcfact*dq1
            end if
            if (this%nvec > 0) then
              corr = 0.0_dp
              ! `av` is stored with the factor index FIRST so this inner loop
              ! walks two contiguous runs; it is the hot loop of the build.
              do v = 1, this%nvec
                corr = corr + this%lam(v)*av(v,i1,j1)*av(v,k1,l1)
              end do
              df1 = df1 + 4*corr
            end if
            dabmax = max(dabmax, abs(df1))
            bfn = 1.0_dp
            if (.not. usecart) bfn = product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i) = df1*bfn
          end do
        end do
      end do
    end do
  end subroutine grd2_cas_get_density

end module casscf_gradient_mod
