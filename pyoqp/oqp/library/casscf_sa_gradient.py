"""SA-CASSCF nuclear gradients: the weighted objective, and individual roots.

Two derivatives, not one, and neither is the state-specific formula of
:mod:`oqp.library.casscf_gradient`:

* **Path A**, ``[casscf] gradient_state = averaged`` -- the gradient of the
  optimized objective ``L = sum_I w_I E_I``.  ``L`` *is* stationary with
  respect to every wavefunction parameter, so this is the ordinary
  Hellmann-Feynman-plus-Pulay expression evaluated on weight-averaged RDMs.  It
  is what a state-averaged geometry optimization should follow, and it is not
  the gradient of any physical state.

* **Path B**, ``[casscf] gradient_state = <root>`` -- the gradient of one
  averaged root.  ``E_J`` is stationary against CI variation (its vector is an
  eigenvector) and **not** stationary against orbital rotation: only the
  weighted sum is.  Repairing that needs a Lagrangian whose orbital multipliers
  drag the CI in, because the SA orbital gradient depends on the CI vectors.
  The resulting coupled orbital+CI Z-vector reduces, after eliminating the CI
  multipliers in closed form, to

      H~ z = -g^J

  with ``H~`` the SA orbital Hessian *including* CI relaxation -- which is
  exactly the matrix :func:`oqp.library.casscf_hessian.analytic_orbital_hessian`
  already assembles, because that Hessian is the Schur complement of the
  coupled (orbital, CI) Hessian.  No second Hessian implementation exists here,
  and the CI multipliers are rebuilt from the same ``amp`` amplitudes.

The derivation, the redundancy analysis behind the Pulay term, and the list of
refusals are in ``docs/sa_casscf_gradients.md``.  What follows implements it.

Where the work happens
----------------------
Everything up to the derivative integrals is numpy over objects the CASSCF
stack already builds (MO integrals, per-root RDMs, the analytic Hessian and its
amplitudes).  The derivative-integral contraction is the new liboqp entry point
``casscf_ao_gradient`` (``source/modules/casscf_ao_gradient.F90``), which takes
the two-particle density in the *factorized* form every piece of this
construction naturally has -- a short list of symmetric AO matrices -- rather
than as a dense ``nbf^4`` AO tensor, which ``grd2`` could not expand to the
Cartesian-effective basis anyway.

Scope, honestly
---------------
Path B needs the MO ERI tensor (for the Hessian) and a dense ``nbf^4`` MO
two-particle density (for the effective generalized Fock), so it carries the
same ``O(nbf^4)`` memory and validation-grade active-space limits as
``[casscf] hessian = analytic``.  Large active spaces are declined by the
excitation-stack guard rather than approximated.
"""
from __future__ import annotations

import numpy as np

from oqp.library.casscf import (
    _as_f64c,
    _casscf_options,
    _full_rdms,
    _generalized_fock,
    _kappa_matrix,
    _lib_backend,
    _lib_gfock_grad,
    _log,
    _nonredundant_pairs,
    _solve_active_rdms,
)
from oqp.library.casscf_hessian import analytic_orbital_hessian
from oqp.library.fci import (
    _symmetric_eigh,
    _transform_integrals,
    _unpack_lower_triangle,
    contiguous_active_space,
    settings_from_casci_config,
)
from oqp.library.rdm import make_rdm1_spatial, make_rdm2_spatial

__all__ = [
    "averaged_objective_densities",
    "check_densities",
    "relaxed_state_densities",
    "resolve_gradient_state",
    "sa_casscf_analytic_gradient",
    "sa_casscf_gradient",
    "state_average_plan",
    "state_average_requested",
]

#: ``info`` slots of ``casscf_ao_gradient`` (mirrors the CAS_AG_* block in
#: include/oqp.h and source/modules/casscf_ao_gradient.F90).
_CAS_AG_NINFO = 2

#: Relative cutoff below which an eigenvalue of the symmetrized all-active
#: two-particle correction is a structural zero.  Same constant, same reason,
#: as ``DROP_REL`` in ``casscf_gradient.F90``: the eight-fold symmetrization
#: annihilates the antisymmetric half of the composite index.
_FACTOR_DROP_REL = 1.0e-14

#: Floor on the accepted SA orbital-rotation gradient norm, and the multiple of
#: the run's own declared threshold above which the point is only *reported* as
#: loosely converged.  Same policy as the state-specific gradient: the whole
#: construction assumes the SA stationarity condition holds, and the error is
#: first order in |g_SA|.
_STATIONARITY_FLOOR = 1.0e-4
_STATIONARITY_WARN = 1.0e2

#: Default Z-vector conditioning tolerance ([casscf] zvector_tol).  Used
#: relative to the largest |eigenvalue| of H~ for the null-space cutoff, and to
#: the norm of the right-hand side for the null-space leakage and the residual.
_ZVECTOR_TOL = 1.0e-8
#: Default root-degeneracy refusal threshold ([casscf] zvector_degeneracy_tol).
_ZVECTOR_DEGEN_TOL = 1.0e-8
#: Absolute floor on the retained |eigenvalue|, so a numerically zero Hessian
#: cannot produce an enormous z through a relative test alone.
_ZVECTOR_ABS_FLOOR = 1.0e-12
#: Tolerance on max |F_eff - F_eff^T| relative to max |F_eff|.  A non-vanishing
#: antisymmetric part means the Lagrangian is not stationary against some
#: redundant rotation and the Pulay term is not defined (see the docs).
_FOCK_ASYM_TOL = 1.0e-7
#: Tolerance on the two-particle energy computed from the dense MO density
#: versus the factorized AO density handed to liboqp.  These are two
#: independent assemblies of the same object; a mismatch is an index or
#: convention error, not a physical condition.
_FACTORIZATION_TOL = 1.0e-8


# --------------------------------------------------------------------- selection
def state_average_requested(mol, settings) -> bool:
    """True when the run optimized a state-averaged objective."""
    method = str(mol.config["input"]["method"]).strip().lower().replace("_", "-")
    return (method in {"sa-casscf", "sacasscf"}
            or bool(getattr(settings, "state_average_enabled", False)))


def resolve_gradient_state(config, roots):
    """Resolve ``[casscf] gradient_state`` against the averaged root list.

    Returns ``None`` for the weighted objective (Path A) or the 0-based CI root
    index for an individual state (Path B).  Raises ``ValueError`` with the
    available roots named whenever the selector cannot be honoured -- naming a
    root that is not averaged is the mistake this catches, and answering it
    with some other root would be silently wrong.
    """
    raw = (config.get("casscf", {}) or {}).get("gradient_state", "averaged")
    text = str(raw).strip().lower()
    if text in {"", "averaged", "average", "sa", "weighted", "objective"}:
        return None
    try:
        state = int(text)
    except (TypeError, ValueError):
        raise ValueError(
            f"[casscf] gradient_state must be 'averaged' or one of the "
            f"state-averaged roots {sorted(int(r) for r in roots)}; got "
            f"{raw!r}.") from None
    if state not in [int(r) for r in roots]:
        raise ValueError(
            f"[casscf] gradient_state = {state} is not one of the "
            f"state-averaged roots {sorted(int(r) for r in roots)}. An "
            "individual-state gradient is only defined for a root the "
            "objective actually averaged over: the Z-vector response is built "
            "from that root's own orbital gradient and from the CI couplings "
            "the state average constrains.")
    return state


# ------------------------------------------------------------- density algebra
def _separable_rdm2(dmat):
    """``Gamma^sep(D)_pqrs = D_pq D_rs - 1/2 D_ps D_rq``."""
    return (np.einsum("pq,rs->pqrs", dmat, dmat)
            - 0.5 * np.einsum("ps,rq->pqrs", dmat, dmat))


def _cross_rdm2(amat, bmat):
    """Bilinear polarization ``Gamma^sep(A+B) - Gamma^sep(A) - Gamma^sep(B)``.

    This is the one-index transform of a separable density (the product rule
    applied to ``D (x) D``) and the inactive-Fock folding of a transition
    density, written once."""
    return (np.einsum("pq,rs->pqrs", amat, bmat)
            + np.einsum("pq,rs->pqrs", bmat, amat)
            - 0.5 * (np.einsum("ps,rq->pqrs", amat, bmat)
                     + np.einsum("ps,rq->pqrs", bmat, amat)))


def _one_index_transform_rdm2(kmat, rdm2):
    """One-index transform of a four-index tensor by the generator ``K``.

    Same derivation on every slot -- ``(oit X)_pq = sum_m (K_pm X_mq + K_qm
    X_pm)`` generalized -- so this is the tensor counterpart of the commutator
    ``[K, D]`` the one-particle density takes."""
    return (np.einsum("pm,mqrs->pqrs", kmat, rdm2)
            + np.einsum("qm,pmrs->pqrs", kmat, rdm2)
            + np.einsum("rm,pqms->pqrs", kmat, rdm2)
            + np.einsum("sm,pqrm->pqrs", kmat, rdm2))


def _symmetrize_pair_tensor(tensor):
    """Average a four-index tensor over the eight-fold ERI permutation group.

    A two-particle density carries only four of the eight permutations on its
    own (pair exchange and the simultaneous within-pair swap).  ``grd2``
    contracts one representative per orbit with a multiplicity factor, so the
    missing ones must be averaged in explicitly or a mixed quartet would be
    counted wrong.  The averaging is invisible to any contraction against a
    genuinely eight-fold-symmetric object, which is why the dense MO density
    used for the effective Fock keeps the unaveraged form."""
    t = np.asarray(tensor, dtype=float)
    return 0.125 * (t
                    + t.transpose(1, 0, 2, 3)
                    + t.transpose(0, 1, 3, 2)
                    + t.transpose(1, 0, 3, 2)
                    + t.transpose(2, 3, 0, 1)
                    + t.transpose(3, 2, 0, 1)
                    + t.transpose(2, 3, 1, 0)
                    + t.transpose(3, 2, 1, 0))


def _factorize_active_block(block, ncore, nact, nbf):
    """Factorize an all-active four-index block into separable matrix products.

    After the eight-fold symmetrization the block is a symmetric ``nact^2 x
    nact^2`` matrix over the composite index ``(t,u)``, so its
    eigendecomposition turns it into ``sum_k lam_k M^k (x) M^k`` with ``M^k``
    ordinary symmetric matrices -- the same construction
    ``casscf_gradient.F90::build_active_correction`` performs, and for the same
    reason: each ``M^k`` takes the bfnrm folding and the spherical-to-Cartesian
    expansion a density takes, and each term is manifestly eight-fold
    symmetric.  Returns ``(lams, mats)`` with ``mats`` embedded in the full
    orbital space.
    """
    m2 = _symmetrize_pair_tensor(block).reshape(nact * nact, nact * nact)
    evals, evecs = _symmetric_eigh(m2)
    big = float(np.max(np.abs(evals))) if evals.size else 0.0
    thr = _FACTOR_DROP_REL * max(big, 1.0)
    keep = np.abs(evals) > thr
    lams = evals[keep]
    mats = []
    act = slice(ncore, ncore + nact)
    for col in np.nonzero(keep)[0]:
        mat = np.zeros((nbf, nbf))
        mat[act, act] = evecs[:, col].reshape(nact, nact)
        mats.append(mat)
    return lams, mats


def _two_particle_energy_factorized(sep_terms, low_terms, eri):
    """``1/2 sum Gamma (pq|rs)`` from the factorized form, in any one basis.

    Used only as the cross-check against the dense assembly: the two are
    independent constructions of the same tensor, and the scalar is invariant
    under the eight-fold symmetrization that separates them."""
    total = 0.0
    for coef, pmat in sep_terms:
        jmat = np.einsum("ls,mnls->mn", pmat, eri, optimize=True)
        direct = float(np.einsum("mn,mn->", pmat, jmat))
        exch = float(np.einsum("ms,ln,mnls->", pmat, pmat, eri, optimize=True))
        total += coef * (direct - 0.5 * exch)
    for lam, amat in low_terms:
        jmat = np.einsum("ls,mnls->mn", amat, eri, optimize=True)
        total += lam * float(np.einsum("mn,mn->", amat, jmat))
    return 0.5 * total


# ------------------------------------------------------------------- transition
def _symmetrized_transition_rdms(bra, ket, dets, nact):
    """Symmetrized transition 1-/2-RDMs of two determinant-space vectors.

    ``sym T(a,b) = 1/2 [T(a+b) - T(a) - T(b)]`` is exact because every RDM is a
    quadratic form in the CI vector, and it gives
    ``1/2 (<a|E|b> + <b|E|a>)`` -- which is what a contraction against the
    symmetric folded integrals sees.  So the validated ``make_rdm*_spatial``
    engines cover the transition case and no transition-RDM kernel is needed.
    """
    total = bra + ket
    g1 = 0.5 * (make_rdm1_spatial(total, dets, nact)
                - make_rdm1_spatial(bra, dets, nact)
                - make_rdm1_spatial(ket, dets, nact))
    g2 = 0.5 * (make_rdm2_spatial(total, dets, nact)
                - make_rdm2_spatial(bra, dets, nact)
                - make_rdm2_spatial(ket, dets, nact))
    return g1, g2


# --------------------------------------------------------------------- Z-vector
def _solve_zvector(hess, rhs, tol):
    """Solve ``H~ z = rhs`` through the symmetric spectrum, fail-closed.

    ``H~`` is symmetric and generally indefinite -- a converged SA solution
    need not be a minimum of the SA objective in every direction -- so a
    Cholesky or a positive-definite iterative solver is the wrong tool.  The
    eigendecomposition also makes the two conditions that must be checked
    explicit rather than hidden in a residual:

    * a direction whose curvature is numerically zero is dropped, and
    * if the right-hand side has a non-negligible component along a dropped
      direction, the answer does not exist: the SA objective is flat along it
      while ``E_J`` is not, so ``E_J`` depends on which member of a continuum of
      equivalent SA solutions the optimizer stopped at.

    Exact redundancies of the SA objective are harmless (``dL_SA/dtheta``
    vanishes identically along them, so ``d2 L_SA/dtheta dx`` does too); the
    leakage test is what separates them from accidental near-null directions.
    """
    evals, evecs = _symmetric_eigh(hess)
    scale = float(np.max(np.abs(evals))) if evals.size else 0.0
    cutoff = max(tol * max(scale, 1.0e-30), _ZVECTOR_ABS_FLOOR)
    proj = evecs.T @ rhs
    keep = np.abs(evals) > cutoff
    rhs_norm = float(np.linalg.norm(rhs))
    leak = float(np.max(np.abs(proj[~keep]))) if np.any(~keep) else 0.0
    zvec = evecs[:, keep] @ (proj[keep] / evals[keep])
    residual = float(np.linalg.norm(hess @ zvec - rhs))
    info = {
        "n_dropped": int(np.count_nonzero(~keep)),
        "min_kept": float(np.min(np.abs(evals[keep]))) if np.any(keep) else 0.0,
        "max_eig": scale,
        "null_leak": leak,
        "residual": residual,
        "rhs_norm": rhs_norm,
        "z_norm": float(np.linalg.norm(zvec)),
    }
    ref = max(rhs_norm, 1.0)
    if leak > tol * ref:
        raise ValueError(
            "SA-CASSCF Z-vector: the state gradient has a component "
            f"({leak:.3e}) in the numerical null space of the state-averaged "
            f"orbital Hessian (cutoff {cutoff:.3e} on {info['n_dropped']} "
            "direction(s)). The individual-state gradient is not defined "
            "there: the averaged objective is flat along that rotation while "
            "the individual state is not, so the state energy depends on which "
            "of a continuum of equivalent SA solutions the optimization "
            "stopped at.")
    if residual > tol * ref:
        raise ValueError(
            f"SA-CASSCF Z-vector did not solve: |H z + g| = {residual:.3e} "
            f"against |g| = {rhs_norm:.3e}. The state-averaged orbital Hessian "
            "is too ill-conditioned for a reliable response at this geometry.")
    return zvec, info


def _ci_response_vectors(resp, weights, roots, zvec, ci, degen_tol):
    """CI Lagrange multipliers, reduced to one vector per averaged root.

    ``zeta_I = sum_{M != I, w_M != w_I} |M> (A^{MI} . z) / (E_I - eps_M)``
    with ``A^{MI}_k = <M| dHeff/dkappa_k |c_I>`` the ``amp`` array the Hessian
    build already produced.  The skip rule is literally the ``coef == 0`` skip
    of ``casscf_hess_relax``: an equal-weight partner contributes a pair of
    terms that cancel between ``zeta_I`` and ``zeta_K``, and evaluating them
    would be a catastrophic cancellation over ``+-(E_I - E_K)`` exactly when
    two averaged roots approach each other.

    A vanishing denominator with a live coupling is refused, not regularized --
    same condition, same message shape, as the Hessian builder.
    """
    eps = resp["eps"]
    vecs = resp["vecs"]
    ovl = resp["ovl"]
    noise = resp["coupling_noise"]
    ndet = eps.size
    zetas = []
    for i_avg, (w_i, r) in enumerate(zip(weights, roots)):
        amp = resp["amp"][i_avg]              # (npar, ndet)
        e_i = float(resp["e_avg"][i_avg])
        proj = amp.T @ zvec                   # (ndet,)
        factors = np.zeros(ndet)
        for j in range(ndet):
            if ovl[i_avg, j] > 0.5:
                continue                      # this eigenstate IS root I
            m = int(np.argmax(ovl[:, j]))
            if ovl[m, j] > 0.5:
                coef = float(w_i) - float(weights[m])
            else:
                coef = float(w_i)
            if coef == 0.0:
                continue
            denom = e_i - float(eps[j])
            if abs(denom) < degen_tol:
                if float(np.max(np.abs(amp[:, j]))) < noise:
                    continue                  # spin/symmetry-forbidden partner
                raise ValueError(
                    "SA-CASSCF Z-vector: root degeneracy with non-zero orbital "
                    f"coupling (|E_I - E_J| = {abs(denom):.2e}); the CI "
                    "response is singular and the state-averaged objective is "
                    "not smooth here.")
            factors[j] = proj[j] / denom
        zeta = vecs @ factors
        cvec = ci[:, int(r)]
        cvec = cvec / np.linalg.norm(cvec)
        zeta = zeta - cvec * float(cvec @ zeta)   # <zeta_I|c_I> = 0 exactly
        zetas.append(zeta)
    return zetas


# ------------------------------------------------------------------ liboqp call
def _ao_gradient_backend():
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_ao_gradient"):
        return None
    return lib, ffi


def _run_ao_gradient(mol, nbf, dm_ao, xm_ao, sep_terms, low_terms):
    """Contract the derivative integrals against the factorized AO densities."""
    backend = _ao_gradient_backend()
    if backend is None:
        raise ValueError(
            "The compiled OpenQP library has no casscf_ao_gradient entry "
            "point; rebuild liboqp to use SA-CASSCF analytic gradients.")
    lib, ffi = backend

    nsep = len(sep_terms)
    nvec = len(low_terms)
    sepc = _as_f64c(np.array([c for c, _ in sep_terms], dtype=np.float64))
    sepm = _as_f64c(np.array([m for _, m in sep_terms], dtype=np.float64)
                    .reshape(max(nsep, 1), nbf, nbf) if nsep else
                    np.zeros((1, nbf, nbf)))
    lam = _as_f64c(np.array([c for c, _ in low_terms], dtype=np.float64)
                   if nvec else np.zeros(1))
    avm = _as_f64c(np.array([m for _, m in low_terms], dtype=np.float64)
                   .reshape(max(nvec, 1), nbf, nbf) if nvec else
                   np.zeros((1, nbf, nbf)))
    dmc = _as_f64c(np.ascontiguousarray(dm_ao))
    xmc = _as_f64c(np.ascontiguousarray(xm_ao))
    info = np.zeros(_CAS_AG_NINFO, dtype=np.float64)

    status = int(lib.casscf_ao_gradient(
        mol.data._data, int(nbf),
        ffi.cast("double *", dmc.ctypes.data),
        ffi.cast("double *", xmc.ctypes.data),
        int(nsep),
        ffi.cast("double *", sepc.ctypes.data),
        ffi.cast("double *", sepm.ctypes.data),
        int(nvec),
        ffi.cast("double *", lam.ctypes.data),
        ffi.cast("double *", avm.ctypes.data),
        ffi.cast("double *", info.ctypes.data)))
    if status < 0:
        reason = {
            -30: "the basis size or factor count disagreed with the handle",
            -31: "the gradient driver could not allocate its working arrays",
        }.get(status, f"the native AO-density gradient declined (status {status})")
        raise ValueError(f"SA-CASSCF gradient failed: {reason}.")

    natom = int(mol.data["natom"])
    return np.asarray(mol.get_grad(), dtype=float).reshape((natom, 3)).copy()


# ------------------------------------------------------- basis-free density core
#
# Everything below the AO transform is ordinary linear algebra over MO-basis
# objects, and it is where the whole derivation lives.  Keeping it in functions
# that take arrays rather than a `mol` is what lets it be validated on its own,
# against an abstract SA-CASSCF model whose "geometry" is a scalar parameter in
# the integrals -- no basis set, no derivative integrals, no Pulay term, so a
# disagreement there is a defect in the response algebra and nothing else.
# `tests/test_casscf_sa_gradient.py` runs exactly that protocol before any
# molecule is touched.
def _gfock_and_grad(nbf, ncore, nact, weights, gammas, gammas2, h1e, eri, pairs):
    """Weighted generalized Fock and orbital gradient, engine or numpy.

    The numpy branch is the definition -- ``_generalized_fock`` of the full
    RDMs, then ``2 (F_qp - F_pq)`` over the pairs -- and is both the fallback
    and the pin for the O(nbf^4) engine."""
    built = _lib_gfock_grad(nbf, ncore, nact, np.asarray(weights, dtype=float),
                            gammas, gammas2, h1e, eri, len(pairs))
    if built is not None:
        return built
    dmat = np.zeros((nbf, nbf))
    rdm2 = np.zeros((nbf, nbf, nbf, nbf))
    for w, g1, g2 in zip(weights, gammas, gammas2):
        d_r, g_r = _full_rdms(g1, g2, ncore, nact, nbf)
        dmat += float(w) * d_r
        rdm2 += float(w) * g_r
    fock = _generalized_fock(dmat, rdm2, h1e, eri)
    grad = np.array([fock[q, p] - fock[p, q] for (p, q) in pairs]) * 2.0
    return fock, grad


def averaged_objective_densities(ncore, nact, nbf, weights, gammas, gammas2):
    """Path A densities: the weighted objective, which needs no response.

    The one structural point is that ``Gamma^sep`` is quadratic in ``D`` and so
    does not commute with weight averaging -- except that the only block where
    two root-dependent factors meet is the all-active one, which the averaged
    active 2-RDM overwrites anyway.  Every surviving block is affine in
    ``gamma_I``, so the averaged density has exactly the separable-plus-active
    structure the state-specific path exploits and the factorization is the
    same size.
    """
    gamma_sa = np.einsum("k,kpq->pq", weights, gammas)
    gamma2_sa = np.einsum("k,kpqrs->pqrs", weights, gammas2)
    d_sa, g_sa = _full_rdms(gamma_sa, gamma2_sa, ncore, nact, nbf)
    # The all-active block of Gamma^sep(D) involves only D's active block, so
    # the correction is built at nact^4 rather than sliced out of nbf^4.
    lams, mats = _factorize_active_block(
        gamma2_sa - _separable_rdm2(gamma_sa), ncore, nact, nbf)
    return {
        "d_eff": d_sa,
        "g_eff": g_sa,
        "sep_terms": [(1.0, d_sa)],
        "low_terms": list(zip(lams, mats)),
    }


def relaxed_state_densities(h1e, eri, ncore, nact, active_nelec, pairs, weights,
                            roots, target, ci, dets, gammas, gammas2, *,
                            tol=_ZVECTOR_TOL, degen_tol=_ZVECTOR_DEGEN_TOL,
                            response="full"):
    """Path B densities: one averaged root, through the coupled Z-vector.

    Returns the relaxed one- and two-particle densities, the effective
    generalized Fock they define, the factorized two-particle density for the
    contraction kernel, and the diagnostics that decide whether any of it is
    admissible.

    ``response`` is a test-only lever with no input keyword behind it, because
    the only way to prove a response term is being APPLIED rather than merely
    computed is to remove it and watch the answer break:

    * ``"full"``   -- the coupled orbital + CI Z-vector (production);
    * ``"orbital"``-- the orbital multipliers only, CI multipliers dropped.
      This is the plausible-looking mistake, and it is the one the effective
      Fock's active-active asymmetry detects: the docs claim that block is
      annihilated by CI stationarity and by nothing else, so dropping the CI
      multipliers must break it;
    * ``"none"``   -- no response at all, i.e. the state-specific formula
      applied to an averaged run. Off by a first-order term.
    """
    if response not in {"full", "orbital", "none"}:
        raise ValueError("response must be 'full', 'orbital' or 'none'")
    nbf = int(h1e.shape[0])
    act = slice(ncore, ncore + nact)
    npar = len(pairs)
    k_target = [int(r) for r in roots].index(int(target))

    # ---- state J's own densities and its own orbital gradient.  The engine
    # that builds the weighted Fock builds this one too: the weights are the
    # delta vector on J, so nothing about the construction differs from the
    # averaged case except which RDMs are summed.
    delta_w = np.zeros(len(roots))
    delta_w[k_target] = 1.0
    _fock_j, grad_j = _gfock_and_grad(nbf, ncore, nact, delta_w, gammas,
                                      gammas2, h1e, eri, pairs)
    d_j, g_j = _full_rdms(gammas[k_target], gammas2[k_target], ncore, nact, nbf)

    hess, resp = analytic_orbital_hessian(
        h1e, eri, ncore, nact, active_nelec, pairs, weights, roots, ci,
        degeneracy_tol=degen_tol, return_response=True)

    # ---- the target root must be resolved from every other CI root.  At an
    # exact crossing the individual adiabatic energy is not differentiable,
    # whatever the weights, and no response formulation repairs that; the
    # averaged objective stays well defined there and is the honest answer.
    e_target = float(resp["e_avg"][k_target])
    others = np.abs(resp["eps"] - e_target)
    others[int(np.argmin(others))] = np.inf          # the root itself
    root_gap = float(np.min(others))
    if root_gap < degen_tol:
        raise ValueError(
            f"SA-CASSCF individual-state gradient refused: root {target} is "
            f"degenerate with another CI root to {root_gap:.3e} Eh. The "
            "adiabatic energy is not differentiable at a crossing; the "
            "averaged-objective gradient (gradient_state=averaged) is still "
            "well defined there.")

    if response == "none":
        zvec = np.zeros(npar)
        zinfo = {"n_dropped": 0, "min_kept": 0.0, "max_eig": 0.0,
                 "null_leak": 0.0, "residual": 0.0,
                 "rhs_norm": float(np.linalg.norm(grad_j)), "z_norm": 0.0}
    else:
        zvec, zinfo = _solve_zvector(hess, -np.asarray(grad_j, dtype=float), tol)
    if response == "full":
        zetas = _ci_response_vectors(resp, weights, roots, zvec, ci, degen_tol)
    else:
        zetas = [np.zeros(ci.shape[0]) for _ in roots]

    d_sa = resp["D_sa"]
    g_sa = resp["G_sa"]
    gamma_sa = d_sa[act, act]
    delta_sa = g_sa[act, act, act, act] - _separable_rdm2(gamma_sa)

    # ---- orbital relaxation: the one-index transform of the AVERAGED
    # densities by K(z).  A directional derivative of the fixed-density energy
    # functional is that functional on transformed densities, and the transform
    # is a derivation, so it reaches the separable part by the product rule.
    kmat = _kappa_matrix(zvec, pairs, nbf)
    d_z = kmat @ d_sa - d_sa @ kmat
    g_z = _one_index_transform_rdm2(kmat, g_sa)

    # ---- CI relaxation: transition densities between zeta_I and c_I.  These
    # fold into the full orbital space WITHOUT the inactive occupancy, because
    # <zeta_I|c_I> = 0 kills the core-energy term; what survives of the
    # inactive-Fock folding is exactly a cross two-particle density against the
    # inactive one-particle density.
    gamma_ci = np.zeros((nact, nact))
    gamma2_ci = np.zeros((nact,) * 4)
    for i_avg, (w_i, r) in enumerate(zip(weights, roots)):
        cvec = ci[:, int(r)]
        cvec = cvec / np.linalg.norm(cvec)
        t1, t2 = _symmetrized_transition_rdms(zetas[i_avg], cvec, dets, nact)
        gamma_ci += 2.0 * float(w_i) * t1
        gamma2_ci += 2.0 * float(w_i) * t2

    d_ci = np.zeros((nbf, nbf))
    d_ci[act, act] = gamma_ci
    d_in = np.zeros((nbf, nbf))
    for i in range(ncore):
        d_in[i, i] = 2.0
    g_ci = _cross_rdm2(d_ci, d_in)
    g_ci[act, act, act, act] = gamma2_ci

    # ---- relaxed densities, and the effective generalized Fock they define
    d_eff = d_j + d_z + d_ci
    g_eff = g_j + g_z + g_ci
    fock = _generalized_fock(d_eff, g_eff, h1e, eri)

    # ---- factorized two-particle density.  Every bilinear piece is the
    # polarization of the separable or low-rank forms the contraction kernel
    # already understands, so the relaxed case adds no new AO kernel:
    #   Gamma^cross(A,B) = Gamma^sep(A+B) - Gamma^sep(A) - Gamma^sep(B)
    #   N (x) M + M (x) N = (M+N)(x)(M+N) - M(x)M - N(x)N
    delta_j = gammas2[k_target] - _separable_rdm2(gammas[k_target])
    lam_j, mat_j = _factorize_active_block(delta_j, ncore, nact, nbf)
    lam_sa, mat_sa = _factorize_active_block(delta_sa, ncore, nact, nbf)
    lam_ci, mat_ci = _factorize_active_block(gamma2_ci, ncore, nact, nbf)

    sep_terms = [(1.0, d_j),
                 (1.0, d_z + d_sa), (-1.0, d_z), (-1.0, d_sa),
                 (1.0, d_ci + d_in), (-1.0, d_ci), (-1.0, d_in)]
    low_terms = list(zip(lam_j, mat_j)) + list(zip(lam_ci, mat_ci))
    for lam, mat in zip(lam_sa, mat_sa):
        nmat = kmat @ mat - mat @ kmat
        low_terms += [(lam, mat + nmat), (-lam, mat), (-lam, nmat)]

    return {
        "d_eff": d_eff,
        "g_eff": g_eff,
        "fock": fock,
        "sep_terms": sep_terms,
        "low_terms": low_terms,
        "zvec": zvec,
        "zetas": zetas,
        "diagnostics": {
            "gnorm_state": float(np.linalg.norm(grad_j)),
            "root_gap": root_gap,
            "zvector": zinfo,
            "zeta_norm": float(max(np.linalg.norm(z) for z in zetas)),
            "n_par": int(npar),
        },
    }


def check_densities(d_mo, g_mo, fock_mo, sep_terms, low_terms, eri):
    """Two silent-failure guards, run before any derivative integral is touched.

    * ``max |F - F^T|`` -- the Pulay term is ``-sum (C sym(F) C^T) S^x`` only
      because the antisymmetric part of ``F`` is annihilated by the
      stationarity of the Lagrangian in every *redundant* rotation.  In the
      active-active block that holds only if the CI multipliers were solved
      (see the docs), so a residual asymmetry is a defect in the response
      solution, not a tolerance to widen.
    * the two-particle energy from the dense density against the same quantity
      from the factorized one.  These are independent assemblies of a single
      tensor and only the factorized one is contracted, so agreement is what
      says the factorization and its index conventions are right.  The scalar
      is invariant under the eight-fold symmetrization that separates them,
      which is why the dense form is kept unaveraged.
    """
    fock_mo = np.asarray(fock_mo, dtype=float)
    scale = max(float(np.max(np.abs(fock_mo))), 1.0)
    asym = float(np.max(np.abs(fock_mo - fock_mo.T)))
    if asym > _FOCK_ASYM_TOL * scale:
        raise ValueError(
            "SA-CASSCF gradient refused: the effective generalized Fock is not "
            f"symmetric (max |F - F^T| = {asym:.3e}, scale {scale:.3e}). The "
            "Pulay term is only the energy-weighted density when the "
            "Lagrangian is stationary against every redundant orbital "
            "rotation, so this is a defect in the response solution, not a "
            "tolerance to widen.")

    e2_dense = 0.5 * float(np.einsum("pqrs,pqrs->", g_mo, eri, optimize=True))
    e2_factor = _two_particle_energy_factorized(sep_terms, low_terms, eri)
    mismatch = abs(e2_dense - e2_factor)
    if mismatch > _FACTORIZATION_TOL * max(abs(e2_dense), 1.0):
        raise ValueError(
            "SA-CASSCF gradient refused: the factorized two-particle density "
            f"does not reproduce the dense one ({e2_factor:.12f} vs "
            f"{e2_dense:.12f}). This is an internal consistency failure, not a "
            "property of the input.")
    return {"fock_asymmetry": asym,
            "two_particle_energy": e2_dense,
            "factorization_mismatch": mismatch}


# ------------------------------------------------------------------------ driver
def sa_casscf_gradient(mol, *, response="full"):
    """Nuclear gradient of a state-averaged CASSCF run.

    ``[casscf] gradient_state`` selects between the weighted objective and one
    averaged root.  Returns ``(grad, diagnostics)`` with ``grad`` of shape
    ``(natom, 3)``; the same gradient is left in the native buffer
    ``mol.get_grad()`` reads.
    """
    settings = settings_from_casci_config(mol.config)
    options = _casscf_options(mol.config)
    cfg = mol.config.get("casscf", {}) or {}

    if not state_average_requested(mol, settings):
        raise ValueError(
            "sa_casscf_gradient is for a state-averaged run; a state-specific "
            "CASSCF gradient is oqp.library.casscf_gradient.")
    if mol.config["input"]["functional"]:
        raise ValueError("CASSCF requires HF integrals; unset input.functional")

    nbf = int(mol.data.get_basis()["nbf"])
    enuc = float(mol.mol_energy.nenergy)
    nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
    ncore, nact, active_nelec, _plan = contiguous_active_space(
        nbf, nelec, settings, "SA-CASSCF")

    weights, roots = state_average_plan(settings)
    target = resolve_gradient_state(mol.config, roots)
    tol = float(cfg.get("zvector_tol", _ZVECTOR_TOL))
    degen_tol = float(cfg.get("zvector_degeneracy_tol", _ZVECTOR_DEGEN_TOL))

    # ---- integrals and orbitals at the converged solution
    hcore_ao = _unpack_lower_triangle(
        np.asarray(mol.data["OQP::Hcore"], dtype=float), nbf)
    eri_ao = np.asarray(mol.data["OQP::AO_ERI"], dtype=float).reshape(
        (nbf, nbf, nbf, nbf), order="F")
    coeff = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float).reshape(
        (nbf, nbf)).T
    h1e, eri = _transform_integrals(hcore_ao, eri_ao, coeff)

    # ---- the CI at those orbitals, and the per-root RDMs
    energies, ci, dets, _d_sa, _g_sa, gammas, gammas2 = _solve_active_rdms(
        h1e, eri, ncore, nact, active_nelec, enuc, settings, weights, roots,
        with_g=False)

    pairs = _nonredundant_pairs(ncore, nact, nbf)
    fock_sa, grad_sa = _gfock_and_grad(nbf, ncore, nact, weights, gammas,
                                       gammas2, h1e, eri, pairs)
    gnorm = float(np.linalg.norm(grad_sa))
    declared = float(options.gradient_norm_tol)
    limit = max(_STATIONARITY_FLOOR, _STATIONARITY_WARN * declared)
    if not np.isfinite(gnorm) or gnorm > limit:
        raise ValueError(
            "SA-CASSCF gradient refused: the state-averaged orbital-rotation "
            f"gradient at the reported orbitals is {gnorm:.3e}, above the "
            f"acceptance limit {limit:.3e}. Both the averaged-objective "
            "gradient and the Z-vector response assume that gradient vanishes; "
            "the error is first order in it. Converge the orbital optimization "
            "(tighten [casscf] gradient_norm_tol or raise "
            "max_macro_iterations).")

    diagnostics = {
        "gnorm_sa": gnorm,
        "weights": np.asarray(weights, dtype=float),
        "roots": [int(r) for r in roots],
        "energies": np.asarray(energies, dtype=float),
        "objective": float(np.dot(weights, np.asarray(energies)[roots])),
        "target": target,
    }

    if target is None:
        built = averaged_objective_densities(ncore, nact, nbf, weights, gammas,
                                             gammas2)
        built["fock"] = np.asarray(fock_sa, dtype=float)
        label = "averaged"
    else:
        built = relaxed_state_densities(
            h1e, eri, ncore, nact, active_nelec, pairs, weights, roots, target,
            ci, dets, gammas, gammas2, tol=tol, degen_tol=degen_tol,
            response=response)
        diagnostics.update(built["diagnostics"])
        label = f"root {target}"

    diagnostics.update(check_densities(
        built["d_eff"], built["g_eff"], built["fock"], built["sep_terms"],
        built["low_terms"], eri))
    diagnostics["n_sep_factors"] = len(built["sep_terms"])
    diagnostics["n_low_factors"] = len(built["low_terms"])
    diagnostics["label"] = label

    # ---- to the AO basis, and into the derivative-integral contraction
    fock = np.asarray(built["fock"], dtype=float)
    d_ao = coeff @ built["d_eff"] @ coeff.T
    x_ao = coeff @ (0.5 * (fock + fock.T)) @ coeff.T
    sep_ao = [(c, coeff @ m @ coeff.T) for c, m in built["sep_terms"]]
    low_ao = [(c, coeff @ m @ coeff.T) for c, m in built["low_terms"]]
    grad = _run_ao_gradient(mol, nbf, d_ao, x_ao, sep_ao, low_ao)
    return grad, diagnostics


def state_average_plan(settings):
    """The weights/roots the energy run optimized, resolved the same way.

    Deliberately a copy of ``CASSCF._state_average_plan``'s logic rather than a
    call into it: the gradient must differentiate the objective that was
    actually optimized, so the resolution is pinned by
    ``tests/test_casscf_sa_gradient.py`` against the energy path instead of
    being coupled to a method on the energy class."""
    roots = list(getattr(settings, "state_average_target_roots", ()) or ())
    nstate = int(getattr(settings, "state_average_nstate", 0) or 0)
    if not roots:
        nstate = nstate or int(settings.nroot)
        roots = list(range(max(1, nstate)))
    if getattr(settings, "state_average_equal_weights", False):
        weights = None
    else:
        weights = getattr(settings, "state_average_weights", None)
    if weights is None or len(weights) != len(roots):
        weights = np.full(len(roots), 1.0 / len(roots))
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()
    return weights, [int(r) for r in roots]


# -------------------------------------------------------------------- entry point
def sa_casscf_analytic_gradient(mol):
    """Public entry: gradient array shaped ``(1, natom, 3)`` plus a log block."""
    grad, diag = sa_casscf_gradient(mol)
    target = diag["target"]
    _log(mol, "")
    _log(mol, "   ==============================================")
    _log(mol, "   PyOQP: analytic SA-CASSCF nuclear gradient")
    _log(mol, "   ==============================================")
    if target is None:
        _log(mol, f"   {'differentiated quantity':<34}"
                  f"{'weighted objective':>20}")
    else:
        _log(mol, f"   {'differentiated root':<34}{target:>20d}")
        _log(mol, f"   {'energy of that root':<34}"
                  f"{diag['energies'][target]:>20.10f}")
    _log(mol, f"   {'objective value':<34}{diag['objective']:>20.10f}")
    _log(mol, f"   {'SA orbital gradient |g_SA|':<34}"
              f"{diag['gnorm_sa']:>20.3e}")
    if target is not None:
        z = diag["zvector"]
        _log(mol, f"   {'state orbital gradient |g_J|':<34}"
                  f"{diag['gnorm_state']:>20.3e}")
        _log(mol, f"   {'nearest CI root gap':<34}{diag['root_gap']:>20.3e}")
        _log(mol, f"   {'Z-vector norm':<34}{z['z_norm']:>20.3e}")
        _log(mol, f"   {'Z-vector residual':<34}{z['residual']:>20.3e}")
        _log(mol, f"   {'Hessian null directions':<34}"
                  f"{z['n_dropped']:>20d}")
        _log(mol, f"   {'null-space leakage':<34}{z['null_leak']:>20.3e}")
        _log(mol, f"   {'smallest kept curvature':<34}"
                  f"{z['min_kept']:>20.3e}")
        _log(mol, f"   {'CI relaxation |zeta|':<34}"
                  f"{diag['zeta_norm']:>20.3e}")
    _log(mol, f"   {'effective Fock asymmetry':<34}"
              f"{diag['fock_asymmetry']:>20.3e}")
    _log(mol, f"   {'2-particle factorization check':<34}"
              f"{diag['factorization_mismatch']:>20.3e}")
    _log(mol, f"   {'two-particle density factors':<34}"
              f"{diag['n_sep_factors'] + diag['n_low_factors']:>20d}")
    _log(mol, "   ==============================================")

    natom = int(mol.data["natom"])
    return np.array([grad.copy()]).reshape((1, natom, 3))
