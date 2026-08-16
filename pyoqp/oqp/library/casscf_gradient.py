"""Analytic state-specific CASSCF nuclear gradient.

Thin Python shell over the ``casscf_gradient`` entry point in
``source/modules/casscf_gradient.F90``: this module resolves the active space
and the CI options exactly as :mod:`oqp.library.casscf` does for the energy,
enforces the validity conditions of the state-specific formula, makes the one
call, and reports the stationarity diagnostics the driver hands back.

What "state-specific" buys
--------------------------
At a converged CASSCF solution the energy is stationary with respect to every
wavefunction parameter -- the CI coefficients (the vector is an eigenvector of
the active Hamiltonian) and the orbital rotations (non-redundant ones because
the optimizer drove them to zero, redundant ones by invariance).  The nuclear
gradient is then the plain Hellmann-Feynman-plus-Pulay expression, with no
coupled-perturbed step:

    dE/dx = sum D h^x + 1/2 sum Gamma (..|..)^x - sum X S^x + dV_NN/dx.

This is what makes the analytic gradient cheap here: it costs one CI solve plus
one pass of derivative integrals, independent of the number of nuclei, where
the central-difference alternative costs 6*natom full CASSCF runs.

What it does NOT cover
----------------------
A state-averaged run optimizes ``sum_I w_I E_I``, and the two things one might
want from it are separate extensions:

* the gradient of the AVERAGED objective, which needs the weight-averaged RDMs
  in the same expression -- a small change, but a different quantity;
* the gradient of an INDIVIDUAL averaged state, which is NOT variational: only
  the weighted sum is stationary against orbital rotation, so an individual
  state needs a Lagrangian / Z-vector response solve.

Neither is implemented.  A state-averaged run is refused rather than handed the
state-specific formula, which would be silently wrong by exactly the response
term that is missing.

Non-stationary starting points are refused for the same reason: the expression
above is only the derivative of the CASSCF energy when the orbital-rotation
gradient is zero and the CI vector is stationary.  Both ``|g_orb|`` and the
generalized-Fock asymmetry are checked before a gradient is returned.
"""
from __future__ import annotations

import numpy as np

from oqp.library.casscf import (
    _as_f64c,
    _casscf_options,
    _cas_gradient_wavefunction_arrays,
    _lib_backend,
    _log,
    _CAS_IOPT_INDEX,
    _CAS_DOPT_INDEX,
    _CAS_FD_STEP,
)
from oqp.library.fci import (
    contiguous_active_space,
    settings_from_casci_config,
)

#: ``info`` slots of the Fortran entry point (mirrors the CAS_G_* block in
#: include/oqp.h and source/modules/casscf_gradient.F90).
_CAS_G_GNORM = 0
_CAS_G_ENERGY = 1
_CAS_G_FASYM = 2
_CAS_G_NVEC = 3
_CAS_G_NINFO = 4

#: Driver status codes this module turns into user-facing messages.  -1..-9 are
#: the shared casscf_energy codes; -20 and below are gradient-specific.
_CASG_MESSAGES = {
    -20: "the analytic CASSCF gradient is state-specific and the run is "
         "state-averaged",
    -21: "the analytic CASSCF gradient requires Hartree-Fock integrals",
    -22: "the CASSCF gradient could not allocate its working arrays",
    -23: "the active 2-RDM factorization eigensolve failed",
    -24: "the CASSCF gradient could not re-evaluate the CI at the converged "
         "orbitals",
}

#: Floor on the accepted orbital-rotation gradient norm.  The nuclear gradient
#: error from a non-stationary point is first order in |g_orb|, so a run
#: converged to the default [casscf] gradient_norm_tol = 1e-6 has ample margin;
#: this only rejects the case where the optimizer plainly did not converge.
_STATIONARITY_FLOOR = 1.0e-4
#: Multiple of the run's own declared threshold above which the point is
#: reported as loosely converged (but still differentiated).
_STATIONARITY_WARN = 1.0e2
#: Maximum accepted ``max |F_pq - F_qp|`` in Hartree.  This independently
#: verifies CI and redundant active-active stationarity, including a full
#: active space where there are no non-redundant orbital rotations and
#: therefore ``|g_orb|`` is identically zero.
_FOCK_ASYMMETRY_LIMIT = 1.0e-6


def _casscf_gradient_backend():
    """liboqp ``(lib, ffi)`` for the analytic gradient, or ``None``."""
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_gradient"):
        return None
    return lib, ffi


def _state_average_requested(mol, settings):
    method = str(mol.config["input"]["method"]).strip().lower().replace("_", "-")
    return (method in {"sa-casscf", "sacasscf"}
            or bool(getattr(settings, "state_average_enabled", False)))


def _enforce_stationarity(gnorm, fasym, declared):
    """Require orbital and CI stationarity for the response-free derivative."""
    limit = max(_STATIONARITY_FLOOR, _STATIONARITY_WARN * declared)
    if not np.isfinite(gnorm) or gnorm > limit:
        raise ValueError(
            "Analytic CASSCF gradient refused: the orbital-rotation gradient "
            f"at the reported orbitals is {gnorm:.3e}, above the acceptance "
            f"limit {limit:.3e}. The state-specific gradient expression has no "
            "orbital-response term and is only the energy derivative at a "
            "stationary point. Converge the orbital optimization (tighten "
            "[casscf] gradient_norm_tol or raise max_macro_iterations)."
        )
    if not np.isfinite(fasym) or fasym > _FOCK_ASYMMETRY_LIMIT:
        raise ValueError(
            "Analytic CASSCF gradient refused: the generalized Fock asymmetry "
            f"max|F_pq-F_qp| is {fasym:.3e}, above the acceptance limit "
            f"{_FOCK_ASYMMETRY_LIMIT:.3e} Hartree. The response-free "
            "state-specific expression requires a stationary CI vector and "
            "active-active orbital block. Tighten [ci] eig_tol or improve the "
            "CI eigensolver convergence."
        )
    return limit


def casscf_analytic_gradient(mol):
    """Analytic nuclear gradient of the CASSCF root selected by ``[casscf] root``.

    Returns a ``(1, natom, 3)`` array and leaves the same gradient in the
    native buffer ``mol.get_grad()`` reads, so every gradient-driven runtype
    (``grad``, ``optimize``, ``ts``, ``mep``, ``irc``) consumes it
    through the ordinary path.

    Raises ``ValueError`` with the reason whenever the state-specific formula
    does not apply -- a state-averaged objective, a non-Hartree-Fock
    Hamiltonian, a non-contiguous active space, or orbitals that are not
    stationary.
    """
    settings = settings_from_casci_config(mol.config)
    options = _casscf_options(mol.config)

    if _state_average_requested(mol, settings):
        raise ValueError(
            "Analytic CASSCF gradients are state-specific. A state-averaged "
            "run optimizes sum_I w_I E_I, whose individual state gradients "
            "need a Lagrangian/Z-vector response that is not implemented; run "
            "method=casscf with [casscf] root, or use a numerical gradient."
        )

    backend = _casscf_gradient_backend()
    if backend is None:
        raise ValueError(
            "The compiled OpenQP library has no casscf_gradient entry point; "
            "rebuild liboqp to use analytic CASSCF gradients."
        )
    lib, ffi = backend

    nbf = int(mol.data.get_basis()["nbf"])
    enuc = float(mol.mol_energy.nenergy)
    nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
    ncore, nact, active_nelec, _plan = contiguous_active_space(
        nbf, nelec, settings, "CASSCF")

    root = int(options.root)
    packed = _cas_gradient_wavefunction_arrays(
        settings, options, nbf, ncore, nact, active_nelec, enuc, [root])
    if packed is None:
        raise ValueError(
            "The analytic CASSCF gradient needs the contiguous "
            "inactive/active/virtual partition the native driver runs on; "
            "this active space or CI option combination is not supported. "
            "Use a numerical gradient instead."
        )
    iopt, dopt, _spec = packed

    # The optimizer-only entries are unused by the gradient, but the shared
    # schema still validates a few of them, so give them their defaults rather
    # than leaving zeros that would read as invalid.
    iopt[_CAS_IOPT_INDEX["maxmacro"]] = 0
    iopt[_CAS_IOPT_INDEX["maxhist"]] = 2
    dopt[_CAS_DOPT_INDEX["fd_step"]] = _CAS_FD_STEP

    weights = _as_f64c(np.array([1.0], dtype=np.float64))
    roots = np.ascontiguousarray(np.array([root], dtype=np.int32))
    info = np.zeros(_CAS_G_NINFO, dtype=np.float64)

    status = int(lib.casscf_gradient(
        mol.data._data,
        ffi.cast("int32_t *", iopt.ctypes.data),
        ffi.cast("double *", dopt.ctypes.data),
        ffi.cast("double *", weights.ctypes.data),
        ffi.cast("int32_t *", roots.ctypes.data),
        ffi.cast("double *", info.ctypes.data)))
    if status < 0:
        reason = _CASG_MESSAGES.get(
            status, "the native CASSCF gradient driver declined (status "
                    f"{status})")
        raise ValueError(f"Analytic CASSCF gradient failed: {reason}.")

    gnorm = float(info[_CAS_G_GNORM])
    fasym = float(info[_CAS_G_FASYM])
    declared = float(options.gradient_norm_tol)
    _enforce_stationarity(gnorm, fasym, declared)

    _log(mol, "")
    _log(mol, "   ==============================================")
    _log(mol, "   PyOQP: analytic CASSCF nuclear gradient")
    _log(mol, "   ==============================================")
    _log(mol, f"   {'differentiated root':<34}{root}")
    _log(mol, f"   {'energy of that root':<34}{info[_CAS_G_ENERGY]:>20.10f}")
    _log(mol, f"   {'orbital gradient norm |g_orb|':<34}{gnorm:>20.3e}")
    _log(mol, f"   {'generalized Fock asymmetry':<34}"
              f"{fasym:>20.3e}")
    _log(mol, f"   {'active 2-RDM correction rank':<34}"
              f"{int(round(float(info[_CAS_G_NVEC]))):>20d}")
    if gnorm > _STATIONARITY_WARN * declared:
        _log(mol, "   NOTE: |g_orb| exceeds 100x the requested convergence "
                  "threshold;")
        _log(mol, "         the gradient error is first order in |g_orb|.")
    _log(mol, "   ==============================================")

    natom = int(mol.data["natom"])
    grad = np.asarray(mol.get_grad(), dtype=float).reshape((natom, 3))
    return np.array([grad.copy()]).reshape((1, natom, 3))
