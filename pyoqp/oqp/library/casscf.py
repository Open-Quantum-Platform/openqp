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

One call per run
----------------
The whole optimizer runs inside liboqp behind a single ``casscf_energy(inf)``
entry point (``source/modules/casscf_driver.F90``), the shape every other
OpenQP method already has (``hf_energy(inf)``, ``mp2_energy(inf)``).  What used
to be a Python loop calling ~20 fine-grained kernels -- six crossings per
energy/gradient evaluation, 2*n_par of those per macroiteration for the
finite-difference Hessian -- is one crossing: measured 43,317 -> 1 for
H2O/cc-pVDZ CAS(6,6) and 3,188 -> 1 for the H2O/STO-3G CAS(4,4) example, at an
identical converged energy.

Every ``[casscf] converger`` (``twophase``, ``ah``, ``trah``, ``diis``,
``auto``) and both orbital-Hessian builders (``hessian = fd | analytic``) are
covered.  ``trah`` is the odd one out: it runs the shared trust-region core of
``source/trah_core.F90`` -- the same optimizer the SCF ``converger_type=trah``
uses -- over a Hessian-vector product that never assembles the orbital
Hessian.  The
``ah`` + ``analytic`` combination is the one that used to cost most: the
converger's own control flow, the augmented-Hessian model step and the whole
analytic-Hessian assembly each crossed the boundary per macroiteration.

Everything below stays: the option parsing and validation that produces the
user-facing messages, the state-average plan, and the log -- the driver returns
the macroiteration table and the converger counters, and ``_write_log`` /
``_converger_trace`` format them unchanged.  The whole Python optimizer,
including :mod:`oqp.library.casscf_convergers` and
:mod:`oqp.library.casscf_hessian`, also remains as the numerical pin and as the
fallback whenever the driver declines -- the established ``_gfock_backend()``
pattern.

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
    _determinants,
    active_space_plan,
    solve_active_ci,
    resolve_ci_solve,
    _transform_integrals,
    _unpack_lower_triangle,
    fci_spin_diagnostics,
    settings_from_casci_config,
    solve_fci,
    _symmetric_eigh,
    _lib_backend,
    _as_f64c,
    _fci_lib_threads,
    _FCI_MAX_NSPIN,
    _FCI_SOLVER_CODE,
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


# ------------------------------------------------------- native CASSCF driver
# Mirror of the option schema in source/modules/casscf_driver.F90 (CAS_I_* /
# CAS_D_* parameter block).  The Fortran file is authoritative; the names below
# are matched against it by tests/test_casscf_energy.py, so the two cannot
# drift.  The CI half of the run is configured with fci_solve's own schema,
# which the Fortran imports rather than re-declares.
_CAS_IOPT = (
    "ncore", "nact", "nalpha", "nbeta", "nstate", "nroot", "solver", "maxiter",
    "subspace", "mult", "maxmemory", "nthreads", "maxmacro", "optimizer",
    "canonical", "maxescape", "maxhist", "converger", "hessian", "ah_micro",
    "ah_reject", "diis_space", "diis_start", "auto_stag",
)
_CAS_DOPT = (
    "enuc", "eig_tol", "cutoff", "grad_tol", "ener_tol", "step_tol", "maxrot",
    "shift", "fd_step", "saddle_c", "saddle_e", "ah_start", "ah_maxtr",
    "ah_mintr", "ah_sadc", "ah_sade",
)
_CAS_IOPT_INDEX = {name: i for i, name in enumerate(_CAS_IOPT)}
_CAS_DOPT_INDEX = {name: i for i, name in enumerate(_CAS_DOPT)}
_CAS_OPTIMIZER_CODE = {"newton": 0, "powell": 1}

# Converger / Hessian spellings the driver can run, mapped to its codes.  The
# accepted spellings are exactly ``casscf_convergers._CONV_ALIASES`` and
# ``_hessian_provider``'s -- anything else must come back as ``None`` so the
# Python path runs and raises the established message.  Note ``analytical`` is
# deliberately absent: ``_hessian_provider`` rejects it too.
# ``trah`` used to be a loose alias for ``ah``.  It is now its own converger --
# the shared trust-region core of ``source/trah_core.F90`` over a matrix-free
# CASSCF Hessian-vector product -- and ``ah`` (dense Hessian, augmented-Hessian
# eigenproblem) keeps every spelling it had except that one.
_CAS_CONVERGER_CODE = {
    "twophase": 0, "two-phase": 0, "2phase": 0,
    "ah": 1, "augmented-hessian": 1, "augmentedhessian": 1,
    "diis": 2, "auto": 3, "trah": 4, "": 4, "default": 4,
}
_CAS_CONVERGER_NAME = {0: "twophase", 1: "ah", 2: "diis", 3: "auto", 4: "trah"}
_CAS_HESSIAN_CODE = {
    "": 0, "fd": 0, "finite-difference": 0, "finite_difference": 0,
    "numerical": 0, "default": 0, "analytic": 1, "exact": 1,
}
# Raw spellings for which casscf._optimize does NOT enter the converger
# framework, hence for which no trace is printed.  The framework canonicalizes
# more spellings than this (``2phase`` runs the same two-phase code but through
# ``run_converger``, and therefore does print a trace), so the test is on the
# raw string, not on the resolved code.
# Spellings whose run does NOT go through the converger framework, and for
# which therefore no converger trace is printed.  `""`/`"default"` are absent
# on purpose: they now mean `trah`, which always goes through the framework.
_CAS_NO_FRAMEWORK = ("twophase", "two-phase")
# Length of the driver's `stats` block (CAS_NSTATS in casscf_driver.F90).
_CAS_NSTATS = 10

# Fixed algorithmic constants of the two-phase path.  They are hard-coded in
# the Python functions below (not config keys), so they are passed as the
# fixed values the Python uses rather than becoming new options.
_CAS_FD_STEP = 1.0e-4              # _fd_orbital_hessian hh
_CAS_MAX_ESCAPES = 8               # _escape_saddles max_escapes
_CAS_SADDLE_CURV_TOL = 2.5e-2      # _escape_saddles saddle_curv_tol
_CAS_SADDLE_EGAIN_TOL = 1.0e-3     # _escape_saddles saddle_egain_tol


def _casscf_energy_backend():
    """liboqp ``(lib, ffi)`` for the one-call CASSCF driver, or ``None``."""
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_energy"):
        return None
    return lib, ffi


def _native_converger_codes(cfg):
    """``(converger, hessian, framework)`` driver codes for ``[casscf]``.

    ``None`` means the driver must decline and the Python optimizer must run:
    that is how an unrecognized ``converger`` / ``hessian`` spelling keeps
    producing its established error message, which is raised by
    ``casscf_convergers.run_converger`` / ``_hessian_provider`` and not here.
    ``framework`` says whether ``casscf._optimize`` would have gone through the
    converger framework, i.e. whether a converger trace is printed.

    Read with ``dict.get`` defaults, exactly as the optimizer reads them, so an
    input without the keys takes the native default path."""
    raw_conv = str(cfg.get("converger", "trah")).strip().lower()
    raw_hess = str(cfg.get("hessian", "fd")).strip().lower()
    converger = _CAS_CONVERGER_CODE.get(raw_conv)
    hessian = _CAS_HESSIAN_CODE.get(raw_hess)
    if converger is None or hessian is None:
        return None
    return converger, hessian, raw_conv not in _CAS_NO_FRAMEWORK


def _converger_trace(converger, hessian, framework, nevals, stats):
    """The log trace ``casscf_convergers.run_converger`` would have written.

    Formatting is presentation, so it stays here; the driver returns the
    counters.  Line order and field widths match ``run_converger`` exactly,
    including that the ``ah stagnation`` note comes from ``_ah_converge`` and
    therefore also appears under ``auto``."""
    if not framework:
        return []
    lines = [f"{'converger:':<30}{_CAS_CONVERGER_NAME[converger]}"]
    if hessian == 1:
        lines.append(f"{'orbital hessian:':<30}analytic")
    if converger in (1, 3) and stats["stagnated"]:
        lines.append(f"{'ah stagnation:':<30}trust region collapsed / no descent")
    if converger == 2:
        lines.append(f"{'diis extrapolations used:':<30}"
                     f"{stats['diis_used']}/{stats['diis_tried']}")
    elif converger == 3:
        if stats["auto_code"] == 0:
            lines.append(f"{'auto:':<30}ah converged; no fallback")
        else:
            reason = ("stagnated" if stats["auto_code"] == 1
                      else "hit the macroiteration cap")
            lines.append(f"{'auto:':<30}ah {reason} after {stats['auto_it']} "
                         "macroiterations; falling back to twophase")
    lines.append(f"{'converger CI evaluations:':<30}{nevals}")
    if hessian == 1:
        lines.append(f"{'analytic hessian builds:':<30}{stats['nhess']}")
    return lines


def _lib_casscf_energy(mol, settings, options, nbf, ncore, nact, active_nelec,
                       enuc, weights, roots, cfg=None, converger=0, hessian=0):
    """A complete CASSCF run inside liboqp, or ``None`` to use the Python path.

    One boundary crossing replaces the entire orchestration loop: the
    macroiteration loop, its backtracking line search, the 2*n_par
    finite-difference orbital Hessian (the dominant cost -- 2*n_par CI solves
    per macroiteration, each of which was itself six crossings), the
    level-shifted Newton step, the curvature-gated saddle escape,
    canonicalization, the commit of the optimized orbitals and the final CI
    with its spin diagnostics all happen without returning here.

    Everything that is not on the compute path stays here: option parsing, the
    validation that produces the user-facing messages (``resolve_ci_solve`` is
    still called, and still raises them), the state-average plan and the log.
    The driver returns the macroiteration table so ``_write_log`` formats it
    unchanged.

    The ``ah`` / ``diis`` / ``auto`` convergers and ``hessian = analytic`` run
    inside the driver too; ``cfg`` is the ``[casscf]`` block and the options
    those convergers read are parsed here, with the same helpers and therefore
    the same error messages the Python framework would raise -- and only for
    the converger that actually reads them, so an unused bad value keeps being
    ignored exactly as it is today.

    Returns ``(energies, s2, multiplicity, history, converged, niter, nevals,
    stats)`` or ``None`` -- a missing symbol, an unsupported option
    combination, or any negative driver status -- in which case the caller runs
    the Python optimizer, which remains the numerical pin.  Two of those
    negative statuses are refusals rather than fallbacks (a root degeneracy
    with live coupling, and an excitation stack past the memory guard); the
    Python path re-raises their messages, which is why they come back as
    ``None`` like everything else."""
    backend = _casscf_energy_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if options.optimizer not in _CAS_OPTIMIZER_CODE:
        return None
    if nact <= 0 or 2 * nact > _FCI_MAX_NSPIN:
        return None
    if nbf - ncore - nact < 0:
        return None

    # The optimizer is written around contiguous inactive/active/virtual
    # blocks (`_nonredundant_pairs`, `_full_rdm1`), so an explicit
    # `[cas] active_orbital_indices` selection is not something the driver can
    # represent; decline and let the Python path deal with it.
    plan = active_space_plan(
        nbf, (ncore + active_nelec[0], ncore + active_nelec[1]), settings)
    if (tuple(plan.active) != tuple(range(ncore, ncore + nact))
            or tuple(plan.core) != tuple(range(ncore))
            or tuple(plan.nelec) != tuple(active_nelec)
            or int(plan.norb) != nbf):
        return None

    nroot_ci = max(1, int(max(roots)) + 1)
    # Same call `solve_active_ci` makes, so the CI options the driver runs with
    # are the ones the Python path would have used -- and so the established
    # error messages (max_det, max_memory, nroot vs ndet) still come from here.
    spec = resolve_ci_solve(
        plan.nact, plan.nelec,
        ecore=enuc, nroot=nroot_ci,
        max_det=settings.max_det, max_memory=settings.max_memory,
        eig_tol=settings.eig_tol, integral_cutoff=0.0,
        solver=settings.solver, davidson_maxiter=settings.davidson_maxiter,
        davidson_subspace=(max(settings.davidson_subspace, nroot_ci)
                           if settings.davidson_subspace else 0),
        target_spin=settings.target_spin,
        active_section="[cas]", ci_section="[ci]",
    )
    # Canonicalization re-solves the CI for one root only; the driver carries a
    # single Davidson subspace cap, so an explicit cap would apply the
    # nroot-widened value to that solve too.  Only reachable with an explicit
    # `[ci] davidson_subspace` on an iterative multi-root solve.
    if (spec.solver == "davidson" and settings.davidson_subspace
            and nroot_ci > 1):
        return None

    npar = nact * ncore + (nbf - ncore - nact) * (ncore + nact)
    if npar <= 0:
        return None

    if hessian == 1 and (
        tuple(getattr(settings, "active_orbital_indices", ()) or ())
        or tuple(getattr(settings, "core_orbital_indices", ()) or ())
    ):
        return None      # make_hessian_provider raises for this; let it

    maxmacro = max(0, int(options.max_macro_iterations))
    # Upper bound on the rows an optimization can append: one seed row plus one
    # per macroiteration, and one full re-convergence per accepted escape --
    # doubled because `auto` runs an AH optimization with its escape phase and
    # then a two-phase one with its own.
    maxhist = maxmacro * (2 + 2 * _CAS_MAX_ESCAPES) + 2 * _CAS_MAX_ESCAPES + 4

    iopt = np.zeros(len(_CAS_IOPT), dtype=np.int32)
    iopt[_CAS_IOPT_INDEX["ncore"]] = ncore
    iopt[_CAS_IOPT_INDEX["nact"]] = nact
    iopt[_CAS_IOPT_INDEX["nalpha"]] = spec.nalpha
    iopt[_CAS_IOPT_INDEX["nbeta"]] = spec.nbeta
    iopt[_CAS_IOPT_INDEX["nstate"]] = len(roots)
    iopt[_CAS_IOPT_INDEX["nroot"]] = spec.nroot
    iopt[_CAS_IOPT_INDEX["solver"]] = _FCI_SOLVER_CODE[spec.solver]
    iopt[_CAS_IOPT_INDEX["maxiter"]] = spec.davidson_maxiter
    iopt[_CAS_IOPT_INDEX["subspace"]] = spec.davidson_subspace
    iopt[_CAS_IOPT_INDEX["mult"]] = spec.target_multiplicity or 0
    iopt[_CAS_IOPT_INDEX["maxmemory"]] = spec.max_memory
    iopt[_CAS_IOPT_INDEX["nthreads"]] = _fci_lib_threads()
    iopt[_CAS_IOPT_INDEX["maxmacro"]] = maxmacro
    iopt[_CAS_IOPT_INDEX["optimizer"]] = _CAS_OPTIMIZER_CODE[options.optimizer]
    iopt[_CAS_IOPT_INDEX["canonical"]] = 1 if options.canonicalize else 0
    iopt[_CAS_IOPT_INDEX["maxescape"]] = _CAS_MAX_ESCAPES
    iopt[_CAS_IOPT_INDEX["maxhist"]] = maxhist
    iopt[_CAS_IOPT_INDEX["converger"]] = converger
    iopt[_CAS_IOPT_INDEX["hessian"]] = hessian

    dopt = np.zeros(len(_CAS_DOPT), dtype=np.float64)
    dopt[_CAS_DOPT_INDEX["enuc"]] = spec.ecore
    dopt[_CAS_DOPT_INDEX["eig_tol"]] = spec.eig_tol
    dopt[_CAS_DOPT_INDEX["cutoff"]] = spec.integral_cutoff
    dopt[_CAS_DOPT_INDEX["grad_tol"]] = options.gradient_norm_tol
    dopt[_CAS_DOPT_INDEX["ener_tol"]] = options.energy_decrease_tol
    dopt[_CAS_DOPT_INDEX["step_tol"]] = options.step_norm_tol
    dopt[_CAS_DOPT_INDEX["maxrot"]] = options.max_rotation_norm
    dopt[_CAS_DOPT_INDEX["shift"]] = options.level_shift
    dopt[_CAS_DOPT_INDEX["fd_step"]] = _CAS_FD_STEP
    dopt[_CAS_DOPT_INDEX["saddle_c"]] = _CAS_SADDLE_CURV_TOL
    dopt[_CAS_DOPT_INDEX["saddle_e"]] = _CAS_SADDLE_EGAIN_TOL

    # Converger options.  Parsed with the framework's own helpers -- and only
    # for the converger that reads them, so `[casscf] diis_space = nonsense`
    # keeps being ignored under `converger=ah` exactly as it is today.
    cfg = cfg or {}
    if converger in (1, 3, 4):
        # `trah` reuses the trust-region keys (`ah_start_trust_radius`,
        # `ah_max_micro`) rather than introducing a second spelling for the
        # same two numbers; the parser stays the single one in
        # `casscf_convergers._ah_params`.
        from oqp.library.casscf_convergers import _ah_params
        par = _ah_params(cfg, options)
        iopt[_CAS_IOPT_INDEX["ah_micro"]] = par.max_micro
        iopt[_CAS_IOPT_INDEX["ah_reject"]] = par.max_rejects
        dopt[_CAS_DOPT_INDEX["ah_start"]] = par.start_trust
        dopt[_CAS_DOPT_INDEX["ah_maxtr"]] = par.max_trust
        dopt[_CAS_DOPT_INDEX["ah_mintr"]] = par.min_trust
        dopt[_CAS_DOPT_INDEX["ah_sadc"]] = par.saddle_curv_tol
        dopt[_CAS_DOPT_INDEX["ah_sade"]] = par.saddle_egain_tol
    if converger == 2:
        from oqp.library.casscf_convergers import _cfg_int
        iopt[_CAS_IOPT_INDEX["diis_space"]] = max(2, _cfg_int(cfg, "diis_space", 8))
        iopt[_CAS_IOPT_INDEX["diis_start"]] = max(2, _cfg_int(cfg, "diis_start", 2))
    if converger == 3:
        from oqp.library.casscf_convergers import _cfg_int
        iopt[_CAS_IOPT_INDEX["auto_stag"]] = max(
            1, _cfg_int(cfg, "auto_stagnation", 3))

    w = _as_f64c(np.asarray(weights, dtype=np.float64))
    r = np.ascontiguousarray(np.asarray(roots, dtype=np.int32))
    energies = np.zeros(spec.nroot, dtype=np.float64)
    s2 = np.zeros(spec.nroot, dtype=np.float64)
    history = np.zeros((maxhist, 5), dtype=np.float64)
    stats = np.zeros(_CAS_NSTATS, dtype=np.int32)

    status = int(lib.casscf_energy(
        mol.data._data,
        ffi.cast("int32_t *", iopt.ctypes.data),
        ffi.cast("double *", dopt.ctypes.data),
        ffi.cast("double *", w.ctypes.data),
        ffi.cast("int32_t *", r.ctypes.data),
        ffi.cast("double *", energies.ctypes.data),
        ffi.cast("double *", s2.ctypes.data),
        ffi.cast("double *", history.ctypes.data),
        ffi.cast("int32_t *", stats.ctypes.data)))
    if status < 0:
        return None

    nrows = int(stats[0])
    rows = [(int(round(history[i, 0])), float(history[i, 1]),
             float(history[i, 2]), float(history[i, 3]), float(history[i, 4]))
            for i in range(nrows)]
    multiplicity = np.rint(np.sqrt(np.maximum(0.0, 1.0 + 4.0 * s2))).astype(np.int64)
    multiplicity = np.maximum(multiplicity, 1)
    counters = {
        "nhess": int(stats[4]),
        "diis_used": int(stats[5]),
        "diis_tried": int(stats[6]),
        "stagnated": bool(stats[7]),
        "auto_code": int(stats[8]),
        "auto_it": int(stats[9]),
    }
    return (energies, _as_f64c(s2), multiplicity, rows,
            bool(stats[2]), int(stats[1]), int(stats[3]), counters)


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
    plan = active_space_plan(
        h1e.shape[0], (ncore + active_nelec[0], ncore + active_nelec[1]), settings
    )
    nroot = max(1, int(max(roots)) + 1)
    # One crossing for the whole CI solve.  The Python driver used to run the
    # frozen-core fold, the active gather, the spin-orbital expansion, every
    # Davidson application and the root selection as separate calls -- hundreds
    # of crossings per macroiteration, and this is the innermost thing CASSCF
    # does.  solve_active_ci falls back to that driver when the engine declines.
    energies, coeffs, _s2 = solve_active_ci(
        h1e, eri, plan, enuc, settings,
        nroot=nroot, want_s2=False, use_target_spin=True,
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
            F, grad = built
        # Everything an orbital converger may need beyond (objective, grad) that
        # this closure already has in hand.  Only ``converger = trah`` reads it
        # (for the full-space gradient matrix behind its matrix-free
        # Hessian-vector product and for its preconditioner); attaching it to
        # the closure keeps every other converger's call signature unchanged.
        evaluate.state = {"F": F, "D": D, "h1e": h1e, "eri": eri}
        return objective, grad, energies, coeffs

    C = np.array(coeff, dtype=float, copy=True)

    # Converger dispatch ([casscf] converger = twophase | ah | diis | auto).
    # The key is read with a dict-get default so inputs without it (and explicit
    # 'twophase') take the unchanged two-phase production path below.
    if hasattr(mol, "_casscf_converger_trace"):
        del mol._casscf_converger_trace   # stale trace from a previous run on this mol
    cfg = mol.config.get("casscf", {}) or {}
    converger = str(cfg.get("converger", "trah")).strip().lower()
    # Orbital-Hessian backend ([casscf] hessian = fd | analytic, dict-get
    # default): None keeps the FD builder (byte-identical default behaviour).
    hess_fn = _hessian_provider(cfg, hcore_ao, eri_ao, ncore, nact, active_nelec,
                                pairs, obj_weights, obj_roots, settings)
    # An active space with no inactive and no virtual orbitals -- H4/STO-3G
    # CAS(4,4) with frozen_core=0 is the shipped example -- has NO non-redundant
    # rotations at all.  There is nothing to optimize: CASSCF is CASCI at these
    # orbitals, and the orbital gradient is the empty vector, hence converged.
    # Handled here, before any converger sees it, because it is a property of
    # the active space rather than of the optimizer; `trah` in particular
    # divides by |d|^2 in its trust-boundary root and would raise
    # ZeroDivisionError on the empty direction.  The native driver reaches the
    # same conclusion its own way, declining with CAS_ERR_INPUT
    # (casscf_driver.F90, `ctx%npar <= 0`), which is what routes these runs
    # here in the first place.
    if npar == 0:
        energies, coeffs, _dets, _D, _G, _gam, _Gam = _solve_active_rdms(
            *_transform_integrals(hcore_ao, eri_ao, C), ncore, nact,
            active_nelec, enuc, settings, obj_weights, obj_roots, with_g=False,
        )
        objective = float(np.dot(obj_weights, energies[obj_roots]))
        history = [(0, objective, 0.0, 0.0, 0.0)]
        return C, energies, coeffs, history, True, 0

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
        # One call for the whole run when the default two-phase / FD-Hessian
        # path is selected (casscf_driver.F90).  It reads Hcore / AO_ERI /
        # VEC_MO_A off the handle itself, runs the optimizer, canonicalizes,
        # commits the orbitals and solves the final CI, so nothing below the
        # option parsing crosses the boundary.  `None` means the driver
        # declined and the Python optimizer below runs, unchanged.
        native = None
        cas_cfg = mol.config.get("casscf", {}) or {}
        codes = _native_converger_codes(cas_cfg)
        if codes is not None:
            native = _lib_casscf_energy(
                mol, settings, options, nbf, ncore, nact, active_nelec, enuc,
                weights, roots, cfg=cas_cfg, converger=codes[0],
                hessian=codes[1])
        if native is not None:
            if hasattr(mol, "_casscf_converger_trace"):
                del mol._casscf_converger_trace   # stale trace from a previous run
            energies, s2, mult, history, converged, niter, _nev, _cnt = native
            trace = _converger_trace(codes[0], codes[1], codes[2], _nev, _cnt)
            if trace:
                mol._casscf_converger_trace = trace
        else:
            # Only the Python optimizer needs these NumPy views of the handle's
            # records; the driver reads the records itself, and unpacking the
            # Hcore triangle here is an O(nbf^2) Python loop.
            hcore_ao = _unpack_lower_triangle(
                np.asarray(mol.data["OQP::Hcore"], dtype=float), nbf)
            coeff = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T
            eri_ao = np.asarray(mol.data["OQP::AO_ERI"], dtype=float).reshape(
                (nbf, nbf, nbf, nbf), order="F"
            )
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

        # An exhausted macroiteration budget is a failure, not a result.  The
        # orbitals have already been committed and the CI re-solved on them,
        # so everything downstream -- CASPT2/NEVPT2/QDPT2, which take this as
        # their reference by default, and every gradient-driven runtype --
        # would otherwise proceed from a NONSTATIONARY point and report
        # energies that look ordinary.  The full log above is emitted first so
        # the macroiteration history and final |g_orb| are visible.
        #
        # `max_macro_iterations = 0` is the explicitly supported fixed-orbital
        # scaffold (CASCI on the reference orbitals): nothing was optimized,
        # so there is nothing to have failed.
        if not converged and options.max_macro_iterations > 0:
            raise RuntimeError(
                f"{label} did not converge: {niter} of "
                f"{options.max_macro_iterations} macroiterations used, final "
                f"|g_orb| = {history[-1][3]:.3e} against "
                f"gradient_norm_tol = {options.gradient_norm_tol:.3e}. "
                "Raise [casscf] max_macro_iterations, loosen "
                "gradient_norm_tol, or try another [casscf] converger."
            )
