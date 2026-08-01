"""Selectable orbital-optimization convergers for the native CASSCF.

The macroiteration engine in :mod:`oqp.library.casscf` minimizes the (state-
averaged) CASCI energy over orbital rotations ``C <- C exp(K)``; every energy /
gradient evaluation is one active-space CI solve plus a generalized-Fock build.
This module provides interchangeable *convergers* for that outer loop, selected
with the ``[casscf] converger`` key.  The key is read with ``dict.get`` defaults
(no schema entry is required), and an input without it runs the unchanged
two-phase production path in ``casscf.py`` -- this module is not even imported
on the default path.

Convergers
----------
``twophase`` (default)
    The existing production scheme, unchanged: eigenvalue-floored (level-
    shifted) Newton with backtracking, followed by the curvature-gated saddle
    escape.  Selecting it explicitly routes through the very same functions.

``ah``
    Trust-region augmented Hessian (TRAH-style).  Each macroiteration builds
    the orbital Hessian and solves the level-shifted augmented-Hessian
    eigenproblem in the Hessian eigenbasis,

        [[0, a g^T], [a g, H]] v = lambda v,    x = v[1:] / (a v[0]),

    so the step solves ``(H - lambda) x = -g`` with ``lambda <= min(0, w_min)``:
    a positive-definite shifted Hessian, hence a descent direction even through
    negative curvature.  Microiterations bisect the scale ``a`` until ``|x|``
    matches the trust radius; the radius adapts on the predicted-vs-actual
    energy ratio; uphill trial steps are rejected and re-solved at a smaller
    radius (model re-solves cost no CI evaluations).  A gradient-converged
    point with deep negative orbital curvature is escaped along the lowest
    Hessian mode and re-converged with the same AH loop.  The objective /
    gradient closure already carries the state-average weights, so the scheme
    works unchanged for SA-CASSCF.

``diis``
    Orbital-gradient DIIS (Pulay) acceleration over the existing quasi-Newton
    (or Powell, following ``[casscf] optimizer``) step.  Accumulated-rotation /
    gradient pairs are extrapolated; the extrapolated point is *evaluated* and
    accepted only when its true energy beats the plain step, so acceleration
    can only match or improve the unaccelerated iteration.  The accumulated
    rotation is the BCH-truncated sum of the accepted steps -- exact bookkeeping
    is not required because every candidate is re-evaluated variationally.

``auto``
    ``ah`` with a stagnation watchdog (no objective decrease beyond
    ``energy_decrease_tol`` for ``auto_stagnation`` consecutive
    macroiterations, or a collapsed trust region).  On stagnation or a missed
    macroiteration cap it falls back to the two-phase converger, restarted from
    the best (lowest-objective) point the AH loop reached; AH accepts only
    non-increasing steps, so the fallback never starts above the initial point
    and at worst reproduces today's default behaviour.

Honest cost statement
---------------------
By default ``ah`` reuses the same symmetrized finite-difference Hessian the
default Newton step already builds (``casscf._fd_orbital_hessian``): 2*n_par
gradient evaluations -- i.e. 2*n_par active-space CI solves -- per
macroiteration, the same per-iteration cost as the default converger.  The AH
microiterations and trial-step rejections operate on that fixed quadratic
model and are CI-free.

``[casscf] hessian = analytic`` swaps that builder for the exact MCSCF
orbital-rotation Hessian of :mod:`oqp.library.casscf_hessian` (fixed-CI part
from the SA RDMs and derivative integrals + CI-relaxation part from one dense
diagonalization of the active Hamiltonian).  Each analytic build costs zero
``evaluate`` calls -- the per-macroiteration CI cost drops from ``2*n_par + 1``
solves to ~1 solve plus one dense active-space diagonalization (roughly one
extra CI-solve-equivalent, reported separately in the trace as ``analytic
hessian builds``).  The dense-spectrum relaxation term limits this to
small/medium active spaces (a guard raises beyond); very large systems still
need an iterative Hessian-vector product before ``ah`` scales.  With the key
absent the FD path runs unchanged.

Why this module has no Fortran engine (measured)
------------------------------------------------
The two kernels the CASSCF optimizer leans on were moved to Fortran
(``casscf_kernel.F90``, ``casscf_hess_kernel.F90``); this module was profiled
afterwards and deliberately left in NumPy.  Everything here is control flow
around those kernels, and the only real arithmetic it owns is the dense
eigendecomposition of the augmented-Hessian matrix inside ``_ah_model_step``
(``solve``).  H2O CAS(6,6), ``converger=ah``, ``hessian=analytic`` -- the
configuration that gives this module its largest possible share, since the
analytic Hessian removes the 2*n_par CI solves that otherwise swamp
everything.  Wall-clock split, measured at the sparse-derivative-slab Hessian:

    basis         n_par   AH model step        hessian builds     model share
    cc-pVDZ        140    0.066 s / 16 macro   0.876 s / 16 macro     5.1 %
    aug-cc-pVDZ    276    0.094 s /  4 macro   2.239 s /  4 macro     3.2 %

The share still *falls* as the problem grows, which is the load-bearing
number: doubling n_par costs the model step 5.7x (its dense eigensolve is
O(n_par^3)) but costs the Hessian build 10.2x.  On the default FD path the
same 64 eigensolves are 0.073 s of a 74.8 s run (0.10 %), because 2*n_par CI
solves per macroiteration dwarf everything else.  The other convergers own
less still: ``_diis_coefficients`` is 0.004 s of a 3.8 s ``converger=diis``
run (0.1 %), and ``twophase`` only delegates to ``casscf.py``.

Note how that share moved.  Against the previous Hessian builder the model
step was 1.1 % / 0.3 % on the same two systems; the two commits that made the
analytic Hessian 8-22x faster are what lifted it to 5.1 % / 3.2 %.  So the
trigger for revisiting this is explicit: another ~5-10x off the Hessian build
makes ``_ah_model_step`` the leading term, and the paragraph below is then the
place to start.

The obvious replacement was prototyped before it was declined, and it does not
survive contact with this optimizer.  ``[[0, a g^T], [a g, diag(w)]]`` is an
*arrowhead* matrix, so its lowest eigenpair is the root of the scalar secular
equation ``l - a^2 sum_i ge_i^2/(l - w_i) = 0`` below ``min{w_i: ge_i != 0}``,
with ``x_i = ge_i/(l - w_i)``.  A safeguarded Newton iteration on that root
replaces the O(n_par^3) eigensolve with O(n_par) per step and benchmarks ~4x
faster at n_par=140, 50-80x at n_par=300.  Over 720 randomized cases it gets
``l`` right to 2.7e-15 relative -- and the *step* only to ~3e-4, because
``x_i`` divides by ``l - w_i``, which is itself ~1e-11: an eigenvalue accurate
to 1e-15 absolute is not an ``x`` accurate to anything much.  Worse, in the
regime this optimizer converges into -- small gradient, one deep negative mode
-- the root lies closer to a Hessian eigenvalue than double precision can
resolve, so the iteration lands *on* the pole: in 115 of 180 such cases it
returned a finite, wrong step where the dense path correctly reports no
reference component and the caller falls back to ``_lowest_mode_step``.  Doing
this properly needs the LAPACK ``dlaed4`` treatment -- solve for the offset
``d = l - w_k`` from the nearest pole so the small denominator is never formed
by cancellation.  That is delicate numerical work, and at this module's
measured share it buys at most the 5.1 % above -- on a step whose output feeds
an optimizer whose macroiteration count is already accumulation-order
sensitive.  So the dense solve stays: it is the correct, better-tested code,
and it is not what this optimizer is waiting on.

OpenTrustRegion note: the compiled core ships an OTR bridge
(``source/otr_interface.F90``), but its callbacks are hard-wired to the SCF
Fock/density machinery (``trah_converger``) and it is not linked in builds
configured with ``ENABLE_OPENTRAH=OFF``; there is no generic Python-facing
orbital-rotation entry point.  The ``ah`` converger here is therefore a native
NumPy implementation of the same trust-region AH idea; if a generic OTR
binding (objective + gradient + Hessian-vector callbacks) is exposed to Python
later, it can replace ``_ah_model_step``/``_ah_inner`` behind the same
``converger=ah`` key.

``[casscf]`` keys read here (all optional, ``dict.get`` defaults)
-----------------------------------------------------------------
converger              twophase | ah | diis | auto        (default twophase)
hessian                fd | analytic (read in casscf.py)  (default fd)
ah_start_trust_radius  initial trust radius               (default 0.2)
ah_max_trust_radius    trust-radius ceiling               (default max_rotation_norm)
ah_min_trust_radius    trust-radius floor / stagnation    (default 1.0e-6)
ah_max_micro           AH scale-bisection microiterations (default 32)
ah_max_rejects         uphill-step rejections per macro   (default 6)
ah_saddle_curv_tol     deep-negative-curvature threshold  (default 2.5e-2)
ah_saddle_egain_tol    strict gain to accept an escape    (default 1.0e-3)
diis_space             stored rotation/gradient pairs     (default 8)
diis_start             pairs required before extrapolating (default 2)
auto_stagnation        stalled macroiterations before falling back (default 3)

The shared stopping criteria (``gradient_norm_tol``, ``energy_decrease_tol``,
``step_norm_tol``, ``max_macro_iterations``, ``max_rotation_norm``) are the
existing ``[casscf]`` options and keep their meaning.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from oqp.library.fci import _symmetric_eigh
from scipy.linalg import expm

from oqp.library.casscf import (
    _escape_saddles,
    _fd_orbital_hessian,
    _floored_newton_loop,
    _kappa_matrix,
    _newton_step,
    _powell_step,
)

# Fixed algorithmic constants (deliberately not config keys: standard
# trust-region ratio thresholds and factors, plus the acceptance tolerance the
# two-phase backtracking already uses).
_ACCEPT_TOL = 1.0e-12          # non-increase tolerance for accepting a trial step
_RHO_SHRINK = 0.25             # actual/predicted ratio below which the radius shrinks
_RHO_GROW = 0.75               # ratio above which the radius may grow
_GROW_FACTOR = 2.0
_REJECT_FACTOR = 0.25          # radius factor on a rejected (uphill) step
_PRED_FLOOR = 1.0e-13          # ignore ratio updates when the model decrease is noise
_V0_TOL = 1.0e-10              # AH reference-component cutoff (pure saddle mode)
_ESCAPE_AMPS = (0.3, 0.2, 0.1, -0.1, -0.2, -0.3)   # same kick line search as casscf
_MAX_ESCAPES = 8

_CONV_ALIASES = {
    "twophase": "twophase", "two-phase": "twophase", "2phase": "twophase",
    "default": "twophase", "": "twophase",
    "ah": "ah", "trah": "ah", "augmented-hessian": "ah", "augmentedhessian": "ah",
    "diis": "diis",
    "auto": "auto",
}


def _cfg_float(cfg, key, default):
    try:
        return float(cfg.get(key, default))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"casscf.{key} must be a number") from exc


def _cfg_int(cfg, key, default):
    try:
        return int(str(cfg.get(key, default)).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"casscf.{key} must be an integer") from exc


# --------------------------------------------------------------------------- entry
def run_converger(name, mol, C, evaluate, pairs, nbf, options, obj_weights, obj_roots,
                  hess_fn=None):
    """Dispatch one CASSCF orbital optimization to the selected converger.

    Same contract as the tail of ``casscf._optimize``: returns
    ``(C, energies, coeffs, history, converged, total_iterations)`` with
    ``history`` rows ``(it, objective, dE, |g|, |step|)``.

    ``hess_fn(C, coeffs)`` optionally supplies the orbital Hessian (the
    ``[casscf] hessian = analytic`` backend resolved by ``casscf._optimize``);
    ``None`` keeps the finite-difference builder everywhere, unchanged.
    """
    canonical = _CONV_ALIASES.get(str(name).strip().lower())
    if canonical is None:
        raise ValueError(
            f"[casscf] converger '{name}' is not recognized; "
            "choose twophase (default), ah, diis, or auto")

    cfg = {}
    if mol is not None and isinstance(getattr(mol, "config", None), dict):
        cfg = mol.config.get("casscf", {}) or {}

    ncall = [0]

    def counted(Cm):
        ncall[0] += 1
        return evaluate(Cm)

    trace = []
    if mol is not None:
        mol._casscf_converger_trace = trace
    trace.append(f"{'converger:':<30}{canonical}")

    nhess = [0]
    counted_hess = None
    if hess_fn is not None:
        trace.append(f"{'orbital hessian:':<30}analytic")

        def counted_hess(Cm, coeffs):
            nhess[0] += 1
            return hess_fn(Cm, coeffs)

    if canonical == "twophase":
        result = _twophase_optimize(C, counted, pairs, nbf, options,
                                    obj_weights, obj_roots, hess_fn=counted_hess)
    elif canonical == "ah":
        result, _stagnated = _ah_converge(C, counted, pairs, nbf, options,
                                          _ah_params(cfg, options), trace,
                                          hess_fn=counted_hess)
    elif canonical == "diis":
        result = _diis_optimize(C, counted, pairs, nbf, options, cfg,
                                obj_weights, obj_roots, trace,
                                hess_fn=counted_hess)
    else:  # auto
        result = _auto_optimize(C, counted, pairs, nbf, options, cfg,
                                obj_weights, obj_roots, trace,
                                hess_fn=counted_hess)

    trace.append(f"{'converger CI evaluations:':<30}{ncall[0]}")
    if hess_fn is not None:
        trace.append(f"{'analytic hessian builds:':<30}{nhess[0]}")
    return result


# --------------------------------------------------------------------------- twophase
def _twophase_optimize(C, evaluate, pairs, nbf, options, obj_weights, obj_roots,
                       hess_fn=None):
    """The production two-phase scheme, called through the framework.

    Identical code path to ``casscf._optimize``'s default branch (the same
    functions in the same order); kept here so 'auto' and direct callers can
    reuse it."""
    C, energies, coeffs, history, converged, it, last_curv = _floored_newton_loop(
        C, evaluate, pairs, nbf, options, hess_fn=hess_fn)
    total_it = it
    if converged and options.optimizer == "newton":
        C, energies, coeffs, history, converged, total_it = _escape_saddles(
            C, energies, coeffs, last_curv, evaluate, pairs, nbf, options,
            history, total_it, obj_weights, obj_roots, hess_fn=hess_fn)
    return C, energies, coeffs, history, converged, total_it


# --------------------------------------------------------------------------- AH
@dataclass
class _AHParams:
    start_trust: float
    max_trust: float
    min_trust: float
    max_micro: int
    max_rejects: int
    saddle_curv_tol: float
    saddle_egain_tol: float


def _ah_params(cfg, options) -> _AHParams:
    # <= 0 means "auto" (the schema default), i.e. follow max_rotation_norm;
    # this keeps the dynamic default intact now that the config schema always
    # materializes the key.
    raw_max = _cfg_float(cfg, "ah_max_trust_radius", 0.0)
    max_trust = raw_max if raw_max > 0.0 else max(float(options.max_rotation_norm), 1.0e-3)
    start = min(_cfg_float(cfg, "ah_start_trust_radius", 0.2), max_trust)
    return _AHParams(
        start_trust=start,
        max_trust=max_trust,
        min_trust=_cfg_float(cfg, "ah_min_trust_radius", 1.0e-6),
        max_micro=max(1, _cfg_int(cfg, "ah_max_micro", 32)),
        max_rejects=max(0, _cfg_int(cfg, "ah_max_rejects", 6)),
        saddle_curv_tol=_cfg_float(cfg, "ah_saddle_curv_tol", 2.5e-2),
        saddle_egain_tol=_cfg_float(cfg, "ah_saddle_egain_tol", 1.0e-3),
    )


def _ah_model_step(grad, w, U, trust, max_micro):
    """Level-shifted augmented-Hessian step, constrained to the trust region.

    Works in the Hessian eigenbasis (H = U diag(w) U^T).  For a scale ``a`` the
    lowest eigenpair (lambda, v) of ``[[0, a g^T], [a g, diag(w)]]`` yields
    ``x = v[1:]/(a v[0])`` solving ``(H - lambda) x = -g`` with
    ``lambda <= min(0, w_min)``.  Larger ``a`` means a more negative shift and
    a shorter step, so the microiterations first grow ``a`` until the step is
    inside the region, then bisect in log-``a`` to land in [0.8, 1.0]*trust.
    All of this is dense linear algebra on the fixed model: no CI solves.

    Returns ``(step, shift, predicted_decrease, n_micro)``; ``step`` is None
    when the AH eigenvector has (numerically) no reference component -- the
    pure negative-curvature case the caller handles by stepping along the
    lowest Hessian mode.
    """
    ge = U.T @ grad
    n = ge.size

    def solve(alpha):
        A = np.zeros((n + 1, n + 1))
        A[0, 1:] = alpha * ge
        A[1:, 0] = alpha * ge
        A[1:, 1:][np.diag_indices(n)] = w
        vals, vecs = _symmetric_eigh(A)
        v = vecs[:, 0]
        if abs(v[0]) < _V0_TOL:
            return float(vals[0]), None
        return float(vals[0]), v[1:] / (alpha * v[0])

    def model(x):
        return float(ge @ x + 0.5 * np.dot(w * x, x))

    nmic = 1
    shift, x = solve(1.0)
    if x is None:
        return None, shift, 0.0, nmic
    if float(np.linalg.norm(x)) <= trust:
        return U @ x, shift, model(x), nmic

    # step longer than the radius: grow alpha until inside, then bisect
    a_lo, a_hi = 1.0, 1.0
    x_in, s_in = None, shift
    while nmic < max_micro:
        a_hi *= 4.0
        nmic += 1
        shift, x = solve(a_hi)
        if x is None:
            return None, shift, 0.0, nmic
        if float(np.linalg.norm(x)) <= trust:
            x_in, s_in = x, shift
            break
        a_lo = a_hi
    if x_in is None:
        # microiteration budget exhausted while still outside: hard rescale
        x = x * (trust / float(np.linalg.norm(x)))
        return U @ x, shift, model(x), nmic
    while nmic < max_micro and float(np.linalg.norm(x_in)) < 0.8 * trust:
        a_mid = math.sqrt(a_lo * a_hi)
        nmic += 1
        shift, x = solve(a_mid)
        if x is None:
            break
        if float(np.linalg.norm(x)) <= trust:
            x_in, s_in, a_hi = x, shift, a_mid
        else:
            a_lo = a_mid
    return U @ x_in, s_in, model(x_in), nmic


def _lowest_mode_step(grad, w, U, trust):
    """Saddle-escape trial step along the lowest Hessian mode (downhill sign)."""
    u = U[:, 0]
    overlap = float(grad @ u)
    sgn = -1.0 if overlap > 0.0 else 1.0
    step = (sgn * trust) * u
    pred = float(grad @ step) + 0.5 * float(w[0]) * trust * trust
    return step, pred


def _ah_inner(C, evaluate, pairs, nbf, options, params, stagnation_break=0,
              hess_fn=None):
    """One trust-region AH macroiteration loop (no saddle-escape phase).

    Returns ``(C, energies, coeffs, history, converged, it, curv, stagnated)``
    with ``curv`` the raw eigendecomposition of the Hessian built by the final
    step (same convention as the two-phase loop: zero extra CI solves for the
    caller's curvature test).  ``hess_fn(C, coeffs)`` optionally replaces the
    finite-difference Hessian builder (``[casscf] hessian = analytic``)."""
    tol_g = options.gradient_norm_tol
    trust = params.start_trust
    objective, grad, energies, coeffs = evaluate(C)
    history = [(0, objective, 0.0, float(np.linalg.norm(grad)), 0.0)]
    converged = False
    stagnated = False
    curv = None
    best_obj = objective
    stall = 0
    it = 0
    for it in range(1, options.max_macro_iterations + 1):
        gnorm = float(np.linalg.norm(grad))
        if gnorm < tol_g:
            converged = True
            break

        if hess_fn is not None:
            hess = hess_fn(C, coeffs)                         # zero CI solves
        else:
            hess = _fd_orbital_hessian(C, evaluate, pairs, nbf)   # 2*n_par CI solves
        w, U = _symmetric_eigh(hess)
        curv = (w.copy(), U)

        # trial step + rejection loop (model re-solves are CI-free; each trial
        # energy costs one CI solve)
        obj_old = objective
        accepted = False
        de = 0.0
        step = None
        pred = 0.0
        Cn, obj_new, grad_new, en_new, co_new = C, objective, grad, energies, coeffs
        for _rej in range(params.max_rejects + 1):
            step, _shift, pred, _nmic = _ah_model_step(grad, w, U, trust, params.max_micro)
            if step is None:
                step, pred = _lowest_mode_step(grad, w, U, trust)
            Cn = C @ expm(_kappa_matrix(step, pairs, nbf))
            obj_new, grad_new, en_new, co_new = evaluate(Cn)
            de = obj_new - obj_old
            if de <= _ACCEPT_TOL:
                accepted = True
                break
            trust = max(_REJECT_FACTOR * trust, params.min_trust)
            if trust <= params.min_trust:
                break
        if not accepted:
            stagnated = True
            break

        # adaptive radius from the predicted-vs-actual decrease ratio
        step_norm = float(np.linalg.norm(step))
        if pred < -_PRED_FLOOR:
            rho = de / pred
            if rho < _RHO_SHRINK:
                trust = max(0.5 * trust, params.min_trust)
            elif rho > _RHO_GROW and step_norm >= 0.8 * trust:
                trust = min(_GROW_FACTOR * trust, params.max_trust)

        C, objective, grad, energies, coeffs = Cn, obj_new, grad_new, en_new, co_new
        history.append((it, objective, objective - obj_old,
                        float(np.linalg.norm(grad)), step_norm))
        if abs(objective - obj_old) < options.energy_decrease_tol and \
                step_norm < options.step_norm_tol:
            converged = float(np.linalg.norm(grad)) < tol_g
            break
        if stagnation_break:
            if objective < best_obj - options.energy_decrease_tol:
                best_obj = objective
                stall = 0
            else:
                stall += 1
                if stall >= stagnation_break:
                    stagnated = True
                    break
    return C, energies, coeffs, history, converged, it, curv, stagnated


def _curvature_escape(C, energies, coeffs, curv, evaluate, pairs, nbf,
                      history, total_it, reconverge, curv_tol, egain_tol):
    """Curvature-gated saddle escape, generic over the re-convergence loop.

    Mirrors ``casscf._escape_saddles`` (kick along the lowest Hessian mode,
    line-search both signs, re-converge, accept only a strict energy gain) but
    re-converges with the caller's own loop so the 'ah' converger stays AH
    end-to-end.  ``curv`` is reused from the final step: the common case (a
    genuine minimum) costs zero extra CI solves."""
    if curv is None:
        return C, energies, coeffs, history, True, total_it
    objective = history[-1][1]
    w, U = curv
    escapes = 0
    while float(w.min()) < -curv_tol and escapes < _MAX_ESCAPES:
        vneg = U[:, int(np.argmin(w))]
        best_obj, best_C = None, None
        for amp in _ESCAPE_AMPS:
            Cn = C @ expm(_kappa_matrix(amp * vneg, pairs, nbf))
            on = evaluate(Cn)[0]
            if best_obj is None or on < best_obj:
                best_obj, best_C = on, Cn
        C2, en2, co2, hist2, conv2, it2, curv2 = reconverge(best_C)
        for h in hist2[1:]:
            history.append((total_it + h[0], h[1], h[2], h[3], h[4]))
        total_it += it2
        obj2 = hist2[-1][1]
        if not conv2 or obj2 >= objective - egain_tol or curv2 is None:
            break   # no genuine improvement -> keep the converged point
        C, energies, coeffs, objective = C2, en2, co2, obj2
        w, U = curv2
        escapes += 1
    return C, energies, coeffs, history, True, total_it


def _ah_converge(C, evaluate, pairs, nbf, options, params, trace, stagnation_break=0,
                 hess_fn=None):
    """Full 'ah' converger: trust-region AH loop + curvature-gated escape.

    Returns ``((C, energies, coeffs, history, converged, total_it), stagnated)``.
    """
    C, energies, coeffs, history, converged, it, curv, stagnated = _ah_inner(
        C, evaluate, pairs, nbf, options, params, stagnation_break=stagnation_break,
        hess_fn=hess_fn)
    total_it = it
    if converged:
        def reconverge(Ck):
            return _ah_inner(Ck, evaluate, pairs, nbf, options, params,
                             hess_fn=hess_fn)[:7]
        C, energies, coeffs, history, converged, total_it = _curvature_escape(
            C, energies, coeffs, curv, evaluate, pairs, nbf, history, total_it,
            reconverge, params.saddle_curv_tol, params.saddle_egain_tol)
    if stagnated:
        trace.append(f"{'ah stagnation:':<30}trust region collapsed / no descent")
    return (C, energies, coeffs, history, converged, total_it), stagnated


# --------------------------------------------------------------------------- DIIS
def _diis_coefficients(gradients):
    """Pulay coefficients minimizing |sum_i c_i g_i| with sum c_i = 1.

    Drops the oldest vectors while the bordered B matrix is ill-conditioned.
    Returns the coefficient vector for the *last* ``len(coef)`` stored entries,
    or None when no stable extrapolation exists."""
    gs = list(gradients)
    while len(gs) >= 2:
        n = len(gs)
        B = np.empty((n + 1, n + 1))
        for i in range(n):
            for j in range(i, n):
                B[i, j] = B[j, i] = float(np.dot(gs[i], gs[j]))
        B[:n, n] = 1.0
        B[n, :n] = 1.0
        B[n, n] = 0.0
        rhs = np.zeros(n + 1)
        rhs[n] = 1.0
        try:
            if not np.all(np.isfinite(B)) or np.linalg.cond(B) > 1.0e14:
                raise np.linalg.LinAlgError("ill-conditioned DIIS matrix")
            coef = np.linalg.solve(B, rhs)[:n]
        except np.linalg.LinAlgError:
            gs.pop(0)
            continue
        if np.all(np.isfinite(coef)):
            return coef
        gs.pop(0)
    return None


def _diis_optimize(C, evaluate, pairs, nbf, options, cfg, obj_weights, obj_roots, trace,
                   hess_fn=None):
    """Orbital-gradient DIIS over the existing quasi-Newton/two-phase step.

    Each macroiteration takes the production step (Newton or Powell, following
    ``[casscf] optimizer``) with the production backtracking, records the
    (accumulated-rotation, gradient) pair, and -- once ``diis_start`` pairs
    exist -- evaluates the Pulay-extrapolated point, keeping whichever of the
    two candidates has the lower true energy (one extra CI solve per
    extrapolated macroiteration).  Converged points get the same
    curvature-gated saddle escape as the default converger."""
    nspace = max(2, _cfg_int(cfg, "diis_space", 8))
    nstart = max(2, _cfg_int(cfg, "diis_start", 2))

    C0 = np.array(C, dtype=float, copy=True)
    T = np.zeros(len(pairs))
    objective, grad, energies, coeffs = evaluate(C)
    history = [(0, objective, 0.0, float(np.linalg.norm(grad)), 0.0)]
    store = []          # (accumulated rotation, gradient) pairs
    converged = False
    last_curv = None
    n_extrap = 0
    n_used = 0
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
            Cn = C @ expm(_kappa_matrix(accepted_step, pairs, nbf))
            obj_new, grad_new, en_new, co_new = evaluate(Cn)
            if obj_new <= obj_old + _ACCEPT_TOL:
                break
            accepted_step = 0.5 * accepted_step
        T_new = T + accepted_step

        store.append((np.array(T_new, copy=True), np.array(grad_new, copy=True)))
        if len(store) > nspace:
            store.pop(0)
        if len(store) >= nstart:
            coef = _diis_coefficients([g for _t, g in store])
            if coef is not None:
                sub = store[-len(coef):]
                T_x = np.zeros_like(T_new)
                for c, (t_i, _g_i) in zip(coef, sub):
                    T_x += c * t_i
                C_x = C0 @ expm(_kappa_matrix(T_x, pairs, nbf))
                obj_x, grad_x, en_x, co_x = evaluate(C_x)
                n_extrap += 1
                if obj_x < obj_new - _ACCEPT_TOL:
                    n_used += 1
                    Cn, obj_new, grad_new, en_new, co_new = C_x, obj_x, grad_x, en_x, co_x
                    T_new = T_x
                    store[-1] = (np.array(T_x, copy=True), np.array(grad_x, copy=True))

        step_norm = float(np.linalg.norm(T_new - T))
        C, objective, grad, energies, coeffs, T = Cn, obj_new, grad_new, en_new, co_new, T_new
        history.append((it, objective, objective - obj_old,
                        float(np.linalg.norm(grad)), step_norm))
        if abs(objective - obj_old) < options.energy_decrease_tol and \
                step_norm < options.step_norm_tol:
            converged = float(np.linalg.norm(grad)) < options.gradient_norm_tol
            break

    total_it = it
    if converged and options.optimizer == "newton":
        C, energies, coeffs, history, converged, total_it = _escape_saddles(
            C, energies, coeffs, last_curv, evaluate, pairs, nbf, options,
            history, total_it, obj_weights, obj_roots, hess_fn=hess_fn)
    trace.append(f"{'diis extrapolations used:':<30}{n_used}/{n_extrap}")
    return C, energies, coeffs, history, converged, total_it


# --------------------------------------------------------------------------- auto
def _auto_optimize(C, evaluate, pairs, nbf, options, cfg, obj_weights, obj_roots, trace,
                   hess_fn=None):
    """'ah' with automatic fallback to 'twophase' on stagnation.

    The AH loop accepts only non-increasing steps, so its final point is its
    best point; the fallback restarts the unchanged two-phase converger from
    there with the full macroiteration budget.  If AH made no progress the
    fallback trajectory is exactly today's default, so 'auto' never diverges
    and at worst matches the default converger's result (total iteration count
    then includes the abandoned AH attempts)."""
    stagnation = max(1, _cfg_int(cfg, "auto_stagnation", 3))
    (C1, en1, co1, hist, conv1, total_it), stagnated = _ah_converge(
        C, evaluate, pairs, nbf, options, _ah_params(cfg, options), trace,
        stagnation_break=stagnation, hess_fn=hess_fn)
    if conv1:
        trace.append(f"{'auto:':<30}ah converged; no fallback")
        return C1, en1, co1, hist, conv1, total_it

    reason = "stagnated" if stagnated else "hit the macroiteration cap"
    trace.append(f"{'auto:':<30}ah {reason} after {total_it} macroiterations; "
                 "falling back to twophase")
    C2, en2, co2, hist2, conv2, it2, curv2 = _floored_newton_loop(
        C1, evaluate, pairs, nbf, options, hess_fn=hess_fn)
    offset = total_it
    for h in hist2[1:]:
        hist.append((offset + h[0], h[1], h[2], h[3], h[4]))
    total_it = offset + it2
    if conv2 and options.optimizer == "newton":
        C2, en2, co2, hist, conv2, total_it = _escape_saddles(
            C2, en2, co2, curv2, evaluate, pairs, nbf, options,
            hist, total_it, obj_weights, obj_roots, hess_fn=hess_fn)
    return C2, en2, co2, hist, conv2, total_it
