"""Native CASSCF / SA-CASSCF for OpenQP.

State-specific and state-averaged complete-active-space SCF, built on OpenQP's
compiled core (AO integrals, Hcore) and the validated determinant CI solver in
:mod:`oqp.library.fci`.  The CASSCF energy is the CASCI energy evaluated at the
orbitals that make the orbital gradient vanish.

Algorithm
---------
At fixed orbitals the active-space CI is solved with the verified determinant
solver; the spin-summed 1-/2-RDMs build the full generalized Fock matrix

    F[p,q] = sum_r D[p,r] h[q,r] + sum_rst Gamma[p,r,s,t] (qr|st),

and the orbital-rotation gradient is ``g[p,q] = 2 (F[q,p] - F[p,q])`` over the
non-redundant inactive-active / inactive-virtual / active-virtual pairs.  Both
are assembled by the Fortran engine in ``casscf_kernel.F90``, which contracts
the separable part of Gamma through ordinary J/K matrices and so never forms
the nbf^4 2-RDM the formula above suggests; the Python assembly below stays as
its numerical pin (see ``tests/test_casscf_gfock.py``).  The
orbitals are rotated ``C <- C exp(K)`` with a Newton step built from a numerical
orbital Hessian, made descent-safe by flooring the Hessian eigenvalues.  This is
robust for the small/medium active spaces this validation-grade path targets.

The rotation itself is also native (``casscf_orbrot.F90``): one call builds K
from the non-redundant pair list, exponentiates it by scaling-and-squaring with
the degree-13 Pade approximant, and applies it to the orbitals.  It is the most
frequently executed step of the optimizer -- the default finite-difference
Hessian alone applies 2*n_par of them per macroiteration -- and replacing the
``C @ scipy.linalg.expm(_kappa_matrix(...))`` expression measured x5.2-6.5 on
interleaved microbenchmarks.  ``_kappa_matrix`` below stays as its pin.
Canonicalization's mean-field Fock is native too and deliberately reuses the
generalized-Fock kernel's own J/K builder rather than duplicating it.

For SA-CASSCF the same machinery is used with the weight-averaged RDMs, so the
state-averaged energy ``sum_I w_I E_I`` is stationary.

This module intentionally keeps a small, readable surface and reuses the
first-class ``[cas]`` / ``[ci]`` / ``[state_average]`` options plus a ``[casscf]``
optimizer block.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import time

import numpy as np
from scipy.linalg import expm

from oqp.utils.file_utils import print_module_banner

import oqp
from oqp.library.fci import (
    _active_space,
    _determinants,
    _transform_integrals,
    _unpack_lower_triangle,
    fci_spin_diagnostics,
    settings_from_casci_config,
    solve_fci,
    _symmetric_eigh,
    _lib_backend,
    _as_f64c,
)
from oqp.library.rdm import (
    make_rdm1_spatial,
    make_rdm2_spatial,
    natural_orbital_occupations,
)

_BOOL_TRUE = {"true", "t", "yes", "y", "1", "on"}
_BOOL_FALSE = {"false", "f", "no", "n", "0", "off", ""}


# --------------------------------------------------------------------------- options
@dataclass
class CASSCFOptions:
    max_macro_iterations: int = 50
    optimizer: str = "newton"
    gradient_norm_tol: float = 1.0e-6
    energy_decrease_tol: float = 1.0e-9
    step_norm_tol: float = 1.0e-7
    max_rotation_norm: float = 0.5
    level_shift: float = 1.0e-3
    canonicalize: bool = True
    root: int = 0


def _as_bool(value, label):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in _BOOL_TRUE:
        return True
    if text in _BOOL_FALSE:
        return False
    raise ValueError(f"{label} must be a boolean")


def _casscf_options(config: dict) -> CASSCFOptions:
    raw = config.get("casscf", {}) or {}
    opt = str(raw.get("optimizer", "newton")).strip().lower()
    if opt in {"newton", "newton-ah", "augmented-hessian", "microiteration"}:
        opt = "newton"
    elif opt in {"powell"}:
        opt = "powell"
    else:
        raise ValueError("casscf.optimizer must be 'newton' or 'powell'")
    return CASSCFOptions(
        max_macro_iterations=int(raw.get("max_macro_iterations", 50)),
        optimizer=opt,
        gradient_norm_tol=float(raw.get("gradient_norm_tol", 1.0e-6)),
        energy_decrease_tol=float(raw.get("energy_decrease_tol", 1.0e-9)),
        step_norm_tol=float(raw.get("step_norm_tol", 1.0e-7)),
        max_rotation_norm=float(raw.get("max_rotation_norm", 0.5)),
        level_shift=float(raw.get("level_shift", 1.0e-3)),
        canonicalize=_as_bool(raw.get("canonicalize", True), "casscf.canonicalize"),
        root=int(raw.get("root", 0)),
    )


# --------------------------------------------------------------------------- log
def _log(mol, text: str = "") -> None:
    with open(mol.log, "a") as handle:
        handle.write(text + "\n")


# --------------------------------------------------------------------------- RDMs / Fock
def _full_rdm1(gamma, ncore, nact, nbf):
    """Full spin-summed 1-RDM: 2 on the inactive diagonal, gamma on the active block."""
    active = list(range(ncore, ncore + nact))
    D = np.zeros((nbf, nbf))
    for i in range(ncore):
        D[i, i] = 2.0
    D[np.ix_(active, active)] = gamma
    return D


def _full_rdms(gamma, Gamma_act, ncore, nact, nbf):
    """Full spin-summed 1-/2-RDM (chemist order; E2 = 0.5 sum (pq|rs) G[pqrs])."""
    active = list(range(ncore, ncore + nact))
    D = _full_rdm1(gamma, ncore, nact, nbf)
    G = np.einsum("pq,rs->pqrs", D, D) - 0.5 * np.einsum("ps,rq->pqrs", D, D)
    G[np.ix_(active, active, active, active)] = Gamma_act
    return D, G


def _generalized_fock(D, G, h1e, eri):
    return D @ h1e + np.einsum("mqrs,nqrs->mn", G, eri)


def _averaged_rdm2(gammas, Gammas, weights, ncore, nact, nbf):
    """Weight-averaged full 2-RDM from the per-root active RDMs."""
    G = np.zeros((nbf, nbf, nbf, nbf))
    for w, gamma, Gamma in zip(weights, gammas, Gammas):
        G += w * _full_rdms(gamma, Gamma, ncore, nact, nbf)[1]
    return G


# --------------------------------------------------------------------------- Fortran Fock/gradient
def _gfock_backend():
    """liboqp (lib, ffi) for the Fortran CASSCF Fock/gradient engine, or None."""
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_gfock_grad"):
        return None
    return lib, ffi


def _lib_gfock_grad(nbf, ncore, nact, weights, gammas, Gammas, h1e, eri, npair):
    """Weighted generalized Fock and orbital gradient through the Fortran engine.

    The engine reproduces ``_generalized_fock(*_full_rdms(...))`` followed by
    ``2 (F[q,p] - F[p,q])`` over ``_nonredundant_pairs``, but never forms the
    nbf^4 2-RDM: the separable part of G contracts through ordinary J/K
    matrices in O(nbf^4), leaving only an nact^3-sized correction.  Returns
    ``(F, grad)``, or ``None`` when the engine is unavailable."""
    backend = _gfock_backend()
    if backend is None:
        return None
    lib, ffi = backend
    w = _as_f64c(weights)
    g1 = _as_f64c(gammas)
    g2 = _as_f64c(Gammas)
    hh = _as_f64c(h1e)
    vv = _as_f64c(eri)
    fock = np.zeros((nbf, nbf), dtype=np.float64)
    grad = np.zeros(npair, dtype=np.float64)
    lib.casscf_gfock_grad(
        int(nbf), int(ncore), int(nact), int(w.size),
        ffi.cast("double *", w.ctypes.data),
        ffi.cast("double *", g1.ctypes.data),
        ffi.cast("double *", g2.ctypes.data),
        ffi.cast("double *", hh.ctypes.data),
        ffi.cast("double *", vv.ctypes.data),
        ffi.cast("double *", fock.ctypes.data),
        ffi.cast("double *", grad.ctypes.data))
    return fock, grad


def _nonredundant_pairs(ncore, nact, nbf):
    inactive = list(range(ncore))
    active = list(range(ncore, ncore + nact))
    virtual = list(range(ncore + nact, nbf))
    pairs = [(t, i) for t in active for i in inactive]
    pairs += [(a, i) for a in virtual for i in inactive]
    pairs += [(a, t) for a in virtual for t in active]
    return pairs


def _kappa_matrix(vec, pairs, nbf):
    """Antisymmetric rotation generator K from the non-redundant pair list.

    Kept as the numerical pin for the Fortran engine in ``casscf_orbrot.F90``,
    which builds the same matrix (with the same ``+=`` accumulation semantics
    for a repeated pair) and exponentiates it without returning it."""
    K = np.zeros((nbf, nbf))
    for (p, q), val in zip(pairs, vec):
        K[p, q] += val
        K[q, p] -= val
    return K


# --------------------------------------------------------------------- rotation
_PAIRS_CACHE: dict = {}


def _pairs_array(pairs):
    """int32 [npar,2] view of the pair list, cached on the list's identity.

    The optimizer builds one pair list per run and then rotates orbitals
    2*npar + O(1) times per macroiteration, so re-marshalling the list of
    tuples on every rotation would cost more Python than the rotation saves.
    A single-entry cache is enough (and the identity check keeps it correct if
    a caller does swap lists)."""
    key = id(pairs)
    hit = _PAIRS_CACHE.get(key)
    if hit is not None and hit[0] is pairs:
        return hit[1]
    arr = np.ascontiguousarray(np.asarray(pairs, dtype=np.int32).reshape(-1, 2))
    _PAIRS_CACHE.clear()
    _PAIRS_CACHE[key] = (pairs, arr)
    return arr


def _rotate_backend():
    """liboqp (lib, ffi) for the Fortran orbital-rotation engine, or None."""
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_orbital_rotate"):
        return None
    return lib, ffi


def _orbital_rotate(C, vec, pairs, nbf):
    """``C @ expm(_kappa_matrix(vec, pairs, nbf))`` through the Fortran engine.

    This is the most frequently executed step of the whole CASSCF optimizer:
    every trial step, every backtracking bisection and every one of the 2*npar
    finite-difference displacements behind the default FD Hessian applies one
    rotation.  The engine builds K, exponentiates it by scaling-and-squaring
    with the degree-13 Pade approximant, and applies it in a single call; the
    NumPy/scipy expression below stays as the numerical pin and as the
    fallback when the symbol is unavailable."""
    backend = _rotate_backend()
    if backend is not None:
        lib, ffi = backend
        pr = _pairs_array(pairs)
        v = _as_f64c(vec)
        cin = _as_f64c(C)
        out = np.zeros((nbf, nbf), dtype=np.float64)
        info = lib.casscf_orbital_rotate(
            int(nbf), int(pr.shape[0]),
            ffi.cast("int32_t *", pr.ctypes.data),
            ffi.cast("double *", v.ctypes.data),
            ffi.cast("double *", cin.ctypes.data),
            ffi.cast("double *", out.ctypes.data))
        if int(info) == 0:
            return out
    return C @ expm(_kappa_matrix(vec, pairs, nbf))


def _effective_fock_backend():
    """liboqp (lib, ffi) for the Fortran canonicalization Fock, or None."""
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_effective_fock"):
        return None
    return lib, ffi


def _lib_effective_fock(h1e, eri, D):
    """``h + J - K/2`` through the Fortran engine, or None when unavailable.

    Reuses the very J/K builder the generalized-Fock engine already contains
    (``casscf_kernel.F90`` ``build_jkw``) rather than duplicating it."""
    backend = _effective_fock_backend()
    if backend is None:
        return None
    lib, ffi = backend
    nbf = int(h1e.shape[0])
    d = _as_f64c(D)
    hh = _as_f64c(h1e)
    vv = _as_f64c(eri)
    out = np.zeros((nbf, nbf), dtype=np.float64)
    lib.casscf_effective_fock(
        nbf,
        ffi.cast("double *", d.ctypes.data),
        ffi.cast("double *", hh.ctypes.data),
        ffi.cast("double *", vv.ctypes.data),
        ffi.cast("double *", out.ctypes.data))
    return out


# --------------------------------------------------------------------------- CASCI inside CASSCF
def _solve_active_rdms(h1e, eri, ncore, nact, active_nelec, enuc, settings,
                       weights, roots, with_g=True):
    """Solve the active CI at the current orbitals and collect its RDMs.

    Returns ``(energies, coeffs, dets, D, G, gammas, Gammas)`` where ``D``/``G``
    are the weight-averaged full 1-/2-RDMs and ``gammas``/``Gammas`` are the
    per-root *active* RDMs stacked along axis 0.  ``with_g=False`` skips the
    nbf^4 ``G`` entirely (it returns ``None``), which is what the Fortran
    Fock/gradient engine wants -- it rebuilds the separable part of G on the
    fly from the active RDMs instead of being handed it."""
    h_act, eri_act, _nelec, ecore, _meta = _active_space(
        h1e, eri, (ncore + active_nelec[0], ncore + active_nelec[1]), enuc, settings
    )
    nroot = max(1, int(max(roots)) + 1)
    energies, coeffs = solve_fci(
        h_act, eri_act, active_nelec,
        ecore=ecore, nroot=nroot, max_det=settings.max_det,
        max_memory=settings.max_memory, eig_tol=settings.eig_tol,
        solver=settings.solver, davidson_maxiter=settings.davidson_maxiter,
        davidson_subspace=settings.davidson_subspace, target_spin=settings.target_spin,
        active_section="[cas]", ci_section="[ci]",
    )
    dets = _determinants(nact, active_nelec)
    nbf = h1e.shape[0]
    gammas = np.zeros((len(roots), nact, nact))
    Gammas = np.zeros((len(roots),) + (nact,) * 4)
    D = np.zeros((nbf, nbf))
    G = np.zeros((nbf, nbf, nbf, nbf)) if with_g else None
    for k, (w, r) in enumerate(zip(weights, roots)):
        gamma = make_rdm1_spatial(coeffs[:, r], dets, nact)
        Gamma = make_rdm2_spatial(coeffs[:, r], dets, nact)
        gammas[k] = gamma
        Gammas[k] = Gamma
        if with_g:
            Dr, Gr = _full_rdms(gamma, Gamma, ncore, nact, nbf)
            G += w * Gr
        else:
            Dr = _full_rdm1(gamma, ncore, nact, nbf)
        D += w * Dr
    return np.asarray(energies), coeffs, dets, D, G, gammas, Gammas


def _solve_active(h1e, eri, ncore, nact, active_nelec, enuc, settings, weights, roots):
    """Solve the active CI at the current orbitals; return (energies, averaged D, G)."""
    energies, coeffs, dets, D, G, _g1, _g2 = _solve_active_rdms(
        h1e, eri, ncore, nact, active_nelec, enuc, settings, weights, roots
    )
    return energies, coeffs, dets, D, G


# --------------------------------------------------------------------------- optimizer
def _optimize(mol, hcore_ao, eri_ao, enuc, coeff, ncore, nact, active_nelec,
              settings, weights, roots, options):
    nbf = coeff.shape[1]
    pairs = _nonredundant_pairs(ncore, nact, nbf)
    npar = len(pairs)
    obj_weights = np.asarray(weights, dtype=float)
    obj_roots = list(roots)

    # The Fortran Fock/gradient engine replaces the O(nbf^5) einsum over the
    # explicit nbf^4 2-RDM; probing once keeps the per-evaluation path free of
    # attribute lookups, and `None` falls back to the Python reference.
    native_gfock = _gfock_backend() is not None

    def evaluate(C):
        h1e, eri = _transform_integrals(hcore_ao, eri_ao, C)
        energies, coeffs, dets, D, G, gammas, Gammas = _solve_active_rdms(
            h1e, eri, ncore, nact, active_nelec, enuc, settings, obj_weights,
            obj_roots, with_g=not native_gfock,
        )
        objective = float(np.dot(obj_weights, energies[obj_roots]))
        built = None
        if native_gfock:
            built = _lib_gfock_grad(nbf, ncore, nact, obj_weights, gammas,
                                    Gammas, h1e, eri, npar)
        if built is None:
            if G is None:   # engine went away between the probe and the call
                G = _averaged_rdm2(gammas, Gammas, obj_weights, ncore, nact, nbf)
            F = _generalized_fock(D, G, h1e, eri)
            grad = np.array([2.0 * (F[q, p] - F[p, q]) for (p, q) in pairs])   # sign matches C->C exp(K)
        else:
            _F, grad = built
        return objective, grad, energies, coeffs

    C = np.array(coeff, dtype=float, copy=True)

    # Converger dispatch ([casscf] converger = twophase | ah | diis | auto).
    # The key is read with a dict-get default so inputs without it (and explicit
    # 'twophase') take the unchanged two-phase production path below.
    if hasattr(mol, "_casscf_converger_trace"):
        del mol._casscf_converger_trace   # stale trace from a previous run on this mol
    cfg = mol.config.get("casscf", {}) or {}
    converger = str(cfg.get("converger", "twophase")).strip().lower()
    # Orbital-Hessian backend ([casscf] hessian = fd | analytic, dict-get
    # default): None keeps the FD builder (byte-identical default behaviour).
    hess_fn = _hessian_provider(cfg, hcore_ao, eri_ao, ncore, nact, active_nelec,
                                pairs, obj_weights, obj_roots, settings)
    if converger not in ("", "twophase", "two-phase", "default"):
        from oqp.library.casscf_convergers import run_converger
        return run_converger(converger, mol, C, evaluate, pairs, nbf, options,
                             obj_weights, obj_roots, hess_fn=hess_fn)

    # Phase 1: damped (level-shifted) Newton to the nearest stationary point.
    # This is the original optimizer body and is left bit-for-bit unchanged, so
    # every system the production optimizer already converges is reproduced
    # exactly (same energy, same iteration count, same gradient).
    C, energies, coeffs, history, converged, it, last_curv = _floored_newton_loop(
        C, evaluate, pairs, nbf, options, hess_fn=hess_fn)
    objective = float(np.dot(obj_weights, energies[obj_roots]))
    total_it = it

    # Phase 2: curvature-gated saddle escape (safety net, no-op when Phase 1
    # already found a minimum).  A level-shifted Newton step *floors* negative
    # Hessian eigenvalues, so it can stop AT a spurious CASSCF saddle and report
    # "converged" there.  If the converged point has a *deep* negative orbital
    # curvature (below ``saddle_curv_tol``), step down that mode and re-converge;
    # accept the new point only if the energy strictly drops by more than
    # ``saddle_egain_tol``.  The threshold is set deep enough that the shallow
    # symmetry-breaking instabilities of genuine reference states (which are not
    # worth escaping) never trigger it.
    if converged and options.optimizer == "newton":
        C, energies, coeffs, history, converged, total_it = _escape_saddles(
            C, energies, coeffs, last_curv, evaluate, pairs, nbf, options,
            history, total_it, obj_weights, obj_roots, hess_fn=hess_fn)

    return C, energies, coeffs, history, converged, total_it


def _hessian_provider(cfg, hcore_ao, eri_ao, ncore, nact, active_nelec, pairs,
                      weights, roots, settings):
    """Resolve ``[casscf] hessian = fd (default) | analytic`` to a builder.

    Returns ``None`` for the default -- callers then run the unchanged
    finite-difference path (``_fd_orbital_hessian``) -- or a callable
    ``hess_fn(C, ci_coeffs)`` from :mod:`oqp.library.casscf_hessian` that
    assembles the exact orbital Hessian with zero extra CI solves."""
    choice = str(cfg.get("hessian", "fd")).strip().lower()
    if choice in ("", "fd", "finite-difference", "finite_difference",
                  "numerical", "default"):
        return None
    if choice not in ("analytic", "exact"):
        raise ValueError("casscf.hessian must be 'fd' (default) or 'analytic'")
    from oqp.library.casscf_hessian import make_hessian_provider
    return make_hessian_provider(hcore_ao, eri_ao, ncore, nact, active_nelec,
                                 pairs, weights, roots, settings=settings)


def _floored_newton_loop(C, evaluate, pairs, nbf, options, hess_fn=None):
    """Original level-shifted Newton optimizer body (unchanged behaviour).

    Additionally returns ``last_curv`` = the eigendecomposition ``(w, U)`` of the
    orbital Hessian built by the *final* Newton step, so the caller can test for a
    spurious saddle without recomputing a Hessian (zero extra FCI solves).

    ``hess_fn(C, coeffs)`` optionally supplies the orbital Hessian (the
    ``[casscf] hessian = analytic`` backend); ``None`` keeps the original
    finite-difference builder, bit-for-bit.
    """
    objective, grad, energies, coeffs = evaluate(C)
    history = [(0, objective, 0.0, float(np.linalg.norm(grad)), 0.0)]
    converged = False
    last_curv = None
    it = 0
    for it in range(1, options.max_macro_iterations + 1):
        gnorm = float(np.linalg.norm(grad))
        if gnorm < options.gradient_norm_tol:
            converged = True
            break
        if options.optimizer == "newton":
            hess = hess_fn(C, coeffs) if hess_fn is not None else None
            step, last_curv = _newton_step(C, grad, pairs, nbf, evaluate, options,
                                           return_curv=True, hess=hess)
        else:
            step = _powell_step(C, grad, pairs, nbf, evaluate, options)
        obj_old = objective
        accepted_step = step
        for _bt in range(12):
            Cn = _orbital_rotate(C, accepted_step, pairs, nbf)
            obj_new, grad_new, energies_new, coeffs_new = evaluate(Cn)
            if obj_new <= obj_old + 1.0e-12:
                break
            accepted_step = 0.5 * accepted_step
        C, objective, grad, energies, coeffs = Cn, obj_new, grad_new, energies_new, coeffs_new
        step_norm = float(np.linalg.norm(accepted_step))
        history.append((it, objective, objective - obj_old, float(np.linalg.norm(grad)), step_norm))
        if abs(objective - obj_old) < options.energy_decrease_tol and step_norm < options.step_norm_tol:
            converged = float(np.linalg.norm(grad)) < options.gradient_norm_tol
            break
    return C, energies, coeffs, history, converged, it, last_curv


def _escape_saddles(C, energies, coeffs, last_curv, evaluate, pairs, nbf, options,
                    history, total_it, obj_weights, obj_roots,
                    saddle_curv_tol=2.5e-2, saddle_egain_tol=1.0e-3,
                    max_escapes=8, hess_fn=None):
    """Detect and descend a deep spurious CASSCF saddle, then re-converge.

    The curvature test reuses ``last_curv`` = the eigendecomposition of the
    Hessian already built by the final Newton step, so the common case (a genuine
    minimum) costs ZERO extra FCI solves -- the phase just inspects ``w.min()``
    and returns.  Only a *deep* negative curvature (``< -saddle_curv_tol``)
    triggers a step, and the escaped point is accepted only when it strictly
    lowers the energy by more than ``saddle_egain_tol``; otherwise the original
    (converged) point is kept.  The threshold sits well below the shallow
    symmetry-breaking instabilities of genuine reference states, so the phase is
    a guaranteed no-op for every case the damped-Newton loop already handles and
    can never regress one."""
    if last_curv is None:
        return C, energies, coeffs, history, True, total_it
    objective = float(np.dot(obj_weights, energies[obj_roots]))
    w, U = last_curv
    escapes = 0
    while float(w.min()) < -saddle_curv_tol and escapes < max_escapes:
        vneg = U[:, int(np.argmin(w))]
        # line-search the negative-curvature direction (both signs / magnitudes)
        best_obj, best_C = None, None
        for amp in (0.3, 0.2, 0.1, -0.1, -0.2, -0.3):
            Cn = _orbital_rotate(C, amp * vneg, pairs, nbf)
            on = float(np.dot(obj_weights, evaluate(Cn)[2][obj_roots]))
            if best_obj is None or on < best_obj:
                best_obj, best_C = on, Cn
        # re-converge from the kicked geometry with the damped-Newton loop
        C2, en2, co2, hist2, conv2, it2, curv2 = _floored_newton_loop(
            best_C, evaluate, pairs, nbf, options, hess_fn=hess_fn)
        obj2 = float(np.dot(obj_weights, en2[obj_roots]))
        total_it += it2
        for h in hist2[1:]:
            history.append((total_it, h[1], h[2], h[3], h[4]))
        if not conv2 or obj2 >= objective - saddle_egain_tol or curv2 is None:
            break   # no genuine improvement -> keep the original converged point
        C, energies, coeffs, objective = C2, en2, co2, obj2
        w, U = curv2
        escapes += 1
    return C, energies, coeffs, history, True, total_it


def _fd_orbital_hessian(C, evaluate, pairs, nbf, hh=1.0e-4):
    """Symmetrized finite-difference orbital Hessian from the analytic gradient.

    Cost: 2*n_par gradient evaluations, each one active-space CI solve plus a
    generalized-Fock build.  Shared by the default Newton step and the 'ah'
    converger in :mod:`oqp.library.casscf_convergers`.  ``[casscf] hessian =
    analytic`` replaces this builder with the exact CI-solve-free Hessian of
    :mod:`oqp.library.casscf_hessian`; the default stays FD."""
    npar = len(pairs)
    hess = np.zeros((npar, npar))
    for k in range(npar):
        e = np.zeros(npar)
        e[k] = hh
        _, gp, _, _ = evaluate(_orbital_rotate(C, e, pairs, nbf))
        _, gm, _, _ = evaluate(_orbital_rotate(C, -e, pairs, nbf))
        hess[:, k] = (gp - gm) / (2.0 * hh)
    return 0.5 * (hess + hess.T)


def _newton_step(C, grad, pairs, nbf, evaluate, options, return_curv=False,
                 hess=None):
    if hess is None:
        hess = _fd_orbital_hessian(C, evaluate, pairs, nbf)
    w, U = _symmetric_eigh(hess)
    curv = (w.copy(), U)                       # raw (unfloored) curvature
    wf = np.where(w > options.level_shift, w, options.level_shift)
    step = -(U @ ((U.T @ grad) / wf))
    nrm = np.linalg.norm(step)
    if nrm > options.max_rotation_norm:
        step *= options.max_rotation_norm / nrm
    if return_curv:
        return step, curv
    return step


def _powell_step(C, grad, pairs, nbf, evaluate, options):
    # simple scaled-gradient fallback (no curvature); robust but slow
    step = -grad
    nrm = np.linalg.norm(step)
    if nrm > options.max_rotation_norm:
        step *= options.max_rotation_norm / nrm
    return step


# --------------------------------------------------------------------------- driver
class CASSCF:
    """State-specific (``method=casscf``) and state-averaged (``method=sa-casscf``)
    CASSCF energy."""

    def __init__(self, mol):
        self.mol = mol
        self.settings = settings_from_casci_config(mol.config)
        self.options = _casscf_options(mol.config)
        self.method = str(mol.config["input"]["method"]).strip().lower()

    def energy(self, ref_energy=None):
        mol = self.mol
        settings = self.settings
        options = self.options
        if settings.integral_backend != "native":
            raise ValueError("casscf integral_backend must be native")
        if mol.config["input"]["functional"]:
            raise ValueError("CASSCF requires HF integrals; unset input.functional")
        if mol.config["scf"]["type"] != "rhf":
            raise ValueError("CASSCF currently supports a closed-shell RHF reference")
        if int(mol.data["nelec_A"]) != int(mol.data["nelec_B"]):
            raise ValueError("CASSCF currently supports closed-shell singlets")

        nbf = int(mol.data.get_basis()["nbf"])
        oqp.fci_ao_integrals(mol)
        hcore_ao = _unpack_lower_triangle(np.asarray(mol.data["OQP::Hcore"], dtype=float), nbf)
        coeff = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T
        eri_ao = np.asarray(mol.data["OQP::AO_ERI"], dtype=float).reshape(
            (nbf, nbf, nbf, nbf), order="F"
        )
        enuc = float(mol.mol_energy.nenergy)

        ncore = int(settings.frozen_core)
        nact = int(settings.active_orbitals)
        nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
        active_nelec = (nelec[0] - ncore, nelec[1] - ncore)
        if nact <= 0 or ncore + nact > nbf:
            raise ValueError("CASSCF needs a valid [cas] active_orbitals / frozen_core")
        if min(active_nelec) < 0:
            raise ValueError("frozen_core exceeds the available electron count")

        weights, roots, sa_enabled = self._state_average_plan(settings, options)

        t0 = time.time()
        coeff_opt, energies, coeffs, history, converged, niter = _optimize(
            mol, hcore_ao, eri_ao, enuc, coeff, ncore, nact, active_nelec,
            settings, weights, roots, options,
        )
        if options.canonicalize:
            coeff_opt = self._canonicalize(coeff_opt, hcore_ao, eri_ao, enuc,
                                           ncore, nact, active_nelec, settings)
        # commit optimized orbitals in place (Fortran core reads this buffer)
        target = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float)
        mol.data["OQP::VEC_MO_A"][...] = np.ascontiguousarray(coeff_opt.T.reshape(target.shape))

        # final CI on the optimized orbitals (all requested roots)
        h1e, eri = _transform_integrals(hcore_ao, eri_ao, coeff_opt)
        energies, coeffs, dets, _D, _G = _solve_active(
            h1e, eri, ncore, nact, active_nelec, enuc, settings, weights, roots
        )
        s2, mult = fci_spin_diagnostics(coeffs, dets, nact, active_nelec)

        nroot_report = max(1, int(settings.nroot))
        report_energies = [float(energies[r]) for r in range(min(nroot_report, len(energies)))]
        if sa_enabled:
            sa_energy = float(np.dot(weights, energies[roots]))
            report_energies = [float(energies[r]) for r in roots]
        mol.energies = report_energies
        mol.mol_energy.energy = float(report_energies[0])
        mol.data["OQP::CASSCF_ENERGIES"] = _as_f64c(report_energies)

        self._write_log(ref_energy, ncore, nact, active_nelec, settings, options,
                        history, converged, niter, energies, s2, mult,
                        weights, roots, sa_enabled, time.time() - t0)
        return mol.energies

    # ---- helpers
    def _state_average_plan(self, settings, options):
        sa_enabled = (self.method in {"sa-casscf", "sacasscf"}) or bool(
            getattr(settings, "state_average_enabled", False)
        )
        if not sa_enabled:
            return np.array([1.0]), [int(options.root)], False
        roots = list(getattr(settings, "state_average_target_roots", ()) or ())
        nstate = int(getattr(settings, "state_average_nstate", 0) or 0)
        if not roots:
            nstate = nstate or max(2, int(settings.nroot))
            roots = list(range(nstate))
        weights = getattr(settings, "state_average_weights", None)
        if weights is None or len(weights) != len(roots):
            weights = np.full(len(roots), 1.0 / len(roots))
        weights = np.asarray(weights, dtype=float)
        weights = weights / weights.sum()
        return weights, [int(r) for r in roots], True

    def _canonicalize(self, coeff, hcore_ao, eri_ao, enuc, ncore, nact, active_nelec, settings):
        """Diagonalize the (inactive, active, virtual) blocks of the generalized
        Fock so inactive/virtual carry orbital energies; active stays natural."""
        h1e, eri = _transform_integrals(hcore_ao, eri_ao, coeff)
        energies, coeffs, dets, D, G = _solve_active(
            h1e, eri, ncore, nact, active_nelec, enuc, settings,
            np.array([1.0]), [0],
        )
        nbf = coeff.shape[1]
        # effective 1-particle Fock for canonicalization: F^I + F^A (closed+active mean field)
        Feff = self._effective_fock(h1e, eri, D, ncore, nact)
        Cnew = np.array(coeff, dtype=float, copy=True)
        for block in (range(ncore), range(ncore, ncore + nact), range(ncore + nact, nbf)):
            idx = list(block)
            if len(idx) < 2:
                continue
            sub = Feff[np.ix_(idx, idx)]
            _w, vec = _symmetric_eigh(0.5 * (sub + sub.T))
            Cnew[:, idx] = coeff[:, idx] @ vec
        return Cnew

    @staticmethod
    def _effective_fock(h1e, eri, D, ncore, nact):
        """Closed+active mean-field Fock h + J - K/2 for canonicalization.

        The Fortran engine reuses the J/K builder already inside the
        generalized-Fock kernel; the einsum reference below is the numerical
        pin and the fallback."""
        built = _lib_effective_fock(h1e, eri, D)
        if built is not None:
            return built
        return CASSCF._effective_fock_reference(h1e, eri, D)

    @staticmethod
    def _effective_fock_reference(h1e, eri, D):
        # closed+active mean-field Fock  h + sum_q D[q,q'] [ (pq|q's) - 0.5 (pq'|qs) ]
        J = np.einsum("rs,pqrs->pq", D, eri)
        K = np.einsum("rs,prsq->pq", D, eri)
        return h1e + J - 0.5 * K

    def _write_log(self, ref_energy, ncore, nact, active_nelec, settings, options,
                   history, converged, niter, energies, s2, mult,
                   weights, roots, sa_enabled, wall):
        mol = self.mol
        label = "SA-CASSCF" if sa_enabled else "CASSCF"
        molname = mol.config["input"].get("system", "") if isinstance(mol.config["input"], dict) else ""
        ne = active_nelec[0] + active_nelec[1]
        print_module_banner(
            mol, label,
            "State-Averaged Complete Active Space SCF" if sa_enabled
            else "Complete Active Space Self-Consistent Field")
        _log(mol)
        _log(mol, "   ==============================================")
        _log(mol, f"   PyOQP: {label} (state-averaged) " if sa_enabled else
                  "   PyOQP: Complete Active Space Self-Consistent Field")
        _log(mol, "   ==============================================")
        _log(mol)
        _log(mol, f"   PyOQP method:                       {self.method}")
        _log(mol, f"   PyOQP reference:                    closed-shell RHF")
        if ref_energy:
            _log(mol, f"   PyOQP RHF reference energy:         {ref_energy[0]:18.10f}")
        _log(mol, f"   PyOQP active space:                 CAS({ne},{nact})")
        _log(mol, f"   PyOQP inactive (core) orbitals:     {ncore}")
        _log(mol, f"   PyOQP active orbitals:              {nact}")
        _log(mol, f"   PyOQP optimizer:                    {options.optimizer}")
        # Converger trace: set only by the [casscf] converger framework, so the
        # default (two-phase) log stays byte-identical.
        for line in getattr(mol, "_casscf_converger_trace", ()) or ():
            _log(mol, f"   PyOQP {line}")
        _log(mol, f"   PyOQP CI solver:                    {settings.solver}")
        if sa_enabled:
            _log(mol, f"   PyOQP state-average roots:          {', '.join(str(r) for r in roots)}")
            _log(mol, f"   PyOQP state-average weights:        {', '.join(f'{w:.4f}' for w in weights)}")
        _log(mol)
        _log(mol, "   --- macro iterations ---")
        _log(mol, f"   {'it':>4} {'E(objective)':>20} {'dE':>14} {'|g_orb|':>12} {'|step|':>12}")
        for (it, e, de, g, st) in history:
            de_s = "" if it == 0 else f"{de:14.2e}"
            st_s = "" if it == 0 else f"{st:12.2e}"
            _log(mol, f"   {it:4d} {e:20.10f} {de_s:>14} {g:12.3e} {st_s:>12}")
        _log(mol)
        _log(mol, f"   PyOQP CASSCF converged:             {'yes' if converged else 'no'}")
        _log(mol, f"   PyOQP CASSCF macro iterations:      {niter} / {options.max_macro_iterations}")
        _log(mol, f"   PyOQP CASSCF final |g_orb|:         {history[-1][3]:.3e}")
        if sa_enabled:
            _log(mol, f"   PyOQP state-average energy:         {float(np.dot(weights, energies[roots])):18.10f}")
        _log(mol)
        _log(mol, f"   PyOQP {label} energies")
        report = roots if sa_enabled else list(range(min(max(1, int(settings.nroot)), len(energies))))
        for r in report:
            _log(mol, f"   PyOQP state {r:<3d} {float(energies[r]):18.10f}    "
                      f"<S^2> {float(s2[r]):8.6f}   multiplicity {int(mult[r])}")
        _log(mol)
        _log(mol, f"   PyOQP timing:                       {wall:.2f} s")
        _log(mol)
