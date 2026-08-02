"""Selectable orbital-optimization convergers for the native CASSCF.

The macroiteration engine in :mod:`oqp.library.casscf` minimizes the (state-
averaged) CASCI energy over orbital rotations ``C <- C exp(K)``; every energy /
gradient evaluation is one active-space CI solve plus a generalized-Fock build.
This module provides interchangeable *convergers* for that outer loop, selected
with the ``[casscf] converger`` key.  The key is read with ``dict.get`` defaults
(no schema entry is required), and an input without it runs the unchanged
two-phase production path in ``casscf.py`` -- this module is not even imported
on the default path.

Where these convergers now RUN
------------------------------
All four of them run inside liboqp: ``casscf_driver.F90`` implements
``_ah_inner`` / ``_curvature_escape`` / ``_ah_converge``, ``_diis_optimize`` and
``_auto_optimize`` natively, on top of the same three kernels this module
calls (``casscf_ah.F90``).  ``casscf_energy`` is therefore ONE boundary
crossing for ``converger = ah`` as well, where the loop below was ~650.

The Python below is unchanged and stays for two reasons, both load-bearing:
it is the numerical pin the native loops are validated against in the same
process (``tests/test_casscf_energy.py``), and it is the fallback for every
run the driver declines -- an unrecognized option spelling, a non-sequential
active space, a Hessian refusal.  The option parsing in ``_ah_params`` /
``_cfg_int`` is also still the only parser: ``casscf.py`` calls it to fill the
driver's option block, so the error messages have exactly one source.

Note that moving the loop is an architectural change, not a performance one --
the measurement section below was taken before the port and is still the
honest statement of what this module's arithmetic costs.

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

``trah``
    Trust-region augmented Hessian, matrix-free.  Runs the SAME optimizer core
    the SCF ``converger_type=trah`` uses (``source/trah_core.F90``): the macro
    trust-region loop, the Steihaug-Toint preconditioned CG and the stability
    check are written once, and each method supplies a gradient and a
    Hessian-vector product.  CASSCF's product is the central difference of the
    ANALYTIC orbital gradient along the requested direction -- two CI solves,
    against the ``2*n_par`` an assembled FD Hessian costs and against the
    ``casscf_hess_bmat`` contraction an assembled analytic one costs -- with the
    Baker-Campbell-Hausdorff (frame) term subtracted in closed form so the
    product is the true symmetric Hessian.  The CI is re-solved at both
    displaced points, so CI relaxation is exact.  ``hessian=analytic`` switches
    it to an assembled Hessian with free matrix-vector products instead, which
    is the cheaper choice only while ``n_par`` stays small.

    Note: ``trah`` used to be an accepted *alias* for ``ah``.  It is now its own
    converger; ``ah`` keeps every other spelling it had.

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

What this module's Fortran engine is worth (measured)
-----------------------------------------------------
This module now has one (``casscf_ah.F90``: the AH model step, the lowest-mode
step and the DIIS coefficients), but the measurement below is what sets
expectations for it, so it is kept verbatim.  It was taken when the module was
still pure NumPy, to decide whether to write an engine at all; the answer on
performance grounds was no, and that answer was correct and still is.  The
engine exists because pyoqp is to be a driver and liboqp is to compute, as it
already does for HF/DFT/TDDFT -- not because the profile asked for it.  What
the port actually bought, interleaved against this same module's NumPy
reference (which remains in place as the numerical pin and the fallback):
``_ah_model_step`` x1.04-1.40, ``_diis_coefficients`` x1.9-2.5, and
``_lowest_mode_step`` 4.2x SLOWER -- the last because it is ~3 us of NumPy on
one eigenvector column against ~13 us of fixed cffi marshalling, on a function
called at most once per run.  End to end the difference is below what a shared
machine can resolve.  The numbers below explain why all of that was the
predictable outcome.

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

That conclusion survived the port: ``casscf_ah.F90`` solves the bordered
eigenproblem densely with LAPACK ``DSYEVD``, which is what the NumPy path was
already doing through ``fci._symmetric_eigh`` -> ``oqp_dsyevd``.  Moving it to
Fortran therefore saves the per-microiteration matrix assembly and the Python
round-trip, not the eigensolve -- which is exactly why the measured gain is
x1.04-1.40 and shrinks as ``n_par`` grows.  Anyone tempted to reach for the
arrowhead shortcut inside the engine should read the paragraph above first.

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
converger              twophase | ah | trah | diis | auto (default twophase)
hessian                fd | analytic (read in casscf.py)  (default fd)
ah_start_trust_radius  initial trust radius               (default 0.2)
ah_max_trust_radius    trust-radius ceiling               (default max_rotation_norm)
ah_min_trust_radius    trust-radius floor / stagnation    (default 1.0e-6)
ah_max_micro           AH scale-bisection microiterations (default 32)
ah_max_rejects         uphill-step rejections per macro   (default 6)
ah_saddle_curv_tol     deep-negative-curvature threshold  (default 2.5e-2)
ah_saddle_egain_tol    strict gain to accept an escape    (default 1.0e-3)
                       (`trah` reuses ah_start_trust_radius, ah_max_trust_radius
                        and ah_max_micro rather than spelling them twice)
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

from oqp.library.casscf import (
    _escape_saddles,
    _fd_orbital_hessian,
    _floored_newton_loop,
    _kappa_matrix,
    _lib_effective_fock,
    _newton_step,
    _orbital_rotate,
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
    "ah": "ah", "augmented-hessian": "ah", "augmentedhessian": "ah",
    "diis": "diis",
    "auto": "auto",
    "trah": "trah",
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
            "choose twophase (default), ah, trah, diis, or auto")

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
    elif canonical == "trah":
        result = _trah_optimize(C, evaluate, counted, pairs, nbf, options,
                                _ah_params(cfg, options), trace)
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


def _ah_backend(symbol):
    """liboqp (lib, ffi) exporting ``symbol``, or None."""
    from oqp.library.fci import _lib_backend

    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, symbol):
        return None
    return lib, ffi


def _lib_ah_model_step(grad, w, U, trust, max_micro):
    """AH model step through the Fortran engine.

    Returns the same ``(step, shift, pred, nmic)`` tuple as the NumPy
    reference (``step`` None for the no-reference-component case), or ``None``
    when the engine is unavailable or LAPACK failed -- the caller then runs
    the reference."""
    backend = _ah_backend("casscf_ah_model_step")
    if backend is None:
        return None
    lib, ffi = backend
    npar = int(np.asarray(grad).size)
    g = np.ascontiguousarray(grad, dtype=np.float64)
    ww = np.ascontiguousarray(w, dtype=np.float64)
    uu = np.ascontiguousarray(U, dtype=np.float64)
    step = np.zeros(npar, dtype=np.float64)
    shift = np.zeros(1, dtype=np.float64)
    pred = np.zeros(1, dtype=np.float64)
    nmic = np.zeros(1, dtype=np.int32)
    status = int(lib.casscf_ah_model_step(
        npar,
        ffi.cast("double *", g.ctypes.data),
        ffi.cast("double *", ww.ctypes.data),
        ffi.cast("double *", uu.ctypes.data),
        float(trust), int(max_micro), float(_V0_TOL),
        ffi.cast("double *", step.ctypes.data),
        ffi.cast("double *", shift.ctypes.data),
        ffi.cast("double *", pred.ctypes.data),
        ffi.cast("int32_t *", nmic.ctypes.data)))
    if status == 1:
        return None, float(shift[0]), 0.0, int(nmic[0])
    if status != 0:
        return None                       # LAPACK failure -> Python reference
    return step, float(shift[0]), float(pred[0]), int(nmic[0])


def _ah_model_step(grad, w, U, trust, max_micro):
    """Level-shifted augmented-Hessian step, constrained to the trust region.

    Dispatches to the Fortran engine in ``casscf_ah.F90``; the NumPy
    implementation below is the numerical pin and the fallback.  The engine is
    a one-for-one transcription, including the microiteration budget and the
    fall-through cases, and is pinned to the reference over randomized cases
    (see ``tests/test_casscf_convergers.py``)."""
    built = _lib_ah_model_step(grad, w, U, trust, max_micro)
    if built is not None:
        return built
    return _ah_model_step_reference(grad, w, U, trust, max_micro)


def _ah_model_step_reference(grad, w, U, trust, max_micro):
    """NumPy reference for :func:`_ah_model_step` (the numerical pin).

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
    """Saddle-escape trial step along the lowest Hessian mode (downhill sign).

    Fortran engine first; the NumPy expression below is the numerical pin."""
    backend = _ah_backend("casscf_lowest_mode_step")
    if backend is not None:
        lib, ffi = backend
        npar = int(np.asarray(grad).size)
        g = np.ascontiguousarray(grad, dtype=np.float64)
        ww = np.ascontiguousarray(w, dtype=np.float64)
        uu = np.ascontiguousarray(U, dtype=np.float64)
        step = np.zeros(npar, dtype=np.float64)
        pred = np.zeros(1, dtype=np.float64)
        lib.casscf_lowest_mode_step(
            npar,
            ffi.cast("double *", g.ctypes.data),
            ffi.cast("double *", ww.ctypes.data),
            ffi.cast("double *", uu.ctypes.data), float(trust),
            ffi.cast("double *", step.ctypes.data),
            ffi.cast("double *", pred.ctypes.data))
        return step, float(pred[0])
    return _lowest_mode_step_reference(grad, w, U, trust)


def _lowest_mode_step_reference(grad, w, U, trust):
    """NumPy reference for :func:`_lowest_mode_step` (the numerical pin)."""
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
            Cn = _orbital_rotate(C, step, pairs, nbf)
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
            Cn = _orbital_rotate(C, amp * vneg, pairs, nbf)
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
_DIIS_CONDMAX = 1.0e14         # bordered-B conditioning ceiling


# --------------------------------------------------------------------------- trah
# Guards of the shared TRAH core (source/trah_core.F90); fixed constants there,
# so fixed constants here.
_TRAH_STAB_EIG_TOL = 1.0e-4    # Hessian eigenvalue below -this = unstable point
_TRAH_PRED_FLOOR = 1.0e-11     # model decrease below this is FP noise
_TRAH_DELTA_MIN = 1.0e-4       # trust radius below this = collapsed
_TRAH_GTOL_FP = 1.0e-4         # |g| below this at a collapse counts as converged
_TRAH_STAB_STEP = 1.0e-3       # step above this at small |g| = saddle escape
_TRAH_FD_STEP = 1.0e-4         # displacement norm of the Hessian-vector product


def _trah_hess_vec(C, ev, pairs, nbf, v, gmat, fd_step=_TRAH_FD_STEP):
    """``H.v`` without ever assembling H.

    The central difference of the ANALYTIC orbital gradient along ``v`` -- the
    same object ``_fd_orbital_hessian`` builds column by column, taken along one
    direction instead of all ``npar`` of them, so it costs 2 CI solves rather
    than ``2*npar``.  The CI is re-solved at both displaced points, so the CI
    relaxation is in there exactly.

    ``g`` is the gradient in the DISPLACED frame, so ``dg/dkappa`` is the second
    derivative only up to the Baker-Campbell-Hausdorff term of
    ``exp(K(kappa)) exp(K(delta))``.  That term is exactly antisymmetric in the
    two directions and of order ``|g|`` -- it is why ``_fd_orbital_hessian``
    symmetrizes at all.  Matrix-free there is no transpose to average with, so
    it is subtracted in closed form with ``R = K(v) G - G K(v)`` and the FULL
    antisymmetric gradient matrix ``G`` (the commutator leaks into the
    active-active block, whose gradient is non-zero for SA-CASSCF)."""
    nv = float(np.linalg.norm(v))
    if nv <= 0.0:
        return np.zeros_like(v)
    t = fd_step / nv
    _, gp, _, _ = ev(_orbital_rotate(C, t * v, pairs, nbf))
    _, gm, _, _ = ev(_orbital_rotate(C, -t * v, pairs, nbf))
    jv = (gp - gm) / (2.0 * t)
    K = _kappa_matrix(v, pairs, nbf)
    R = K @ gmat - gmat @ K
    corr = np.array([R[q, p] - R[p, q] for (p, q) in pairs])
    return jv - 0.25 * corr


def _trah_gradient_matrix(F):
    """Full antisymmetric orbital-gradient matrix ``G[p,q] = 2(F[q,p]-F[p,q])``."""
    return 2.0 * (np.asarray(F).T - np.asarray(F))


def _trah_hdiag(state, pairs):
    """Preconditioner diagonal ``2 (n_q - n_p) (eps_p - eps_q)``.

    ``eps`` is the diagonal of the closed+active mean-field Fock (the matrix
    canonicalization already uses) and ``n`` the state-averaged occupations; for
    a virtual/inactive rotation this is the familiar ``4(eps_a - eps_i)``.  It
    only sets the CG's metric, so an approximation costs iterations and never
    accuracy -- the exact diagonal would cost an assembled Hessian, which is the
    whole thing this converger avoids."""
    D = state["D"]
    feff = _lib_effective_fock(state["h1e"], state["eri"], D)
    if feff is None:                      # engine unavailable: unit metric
        return np.ones(len(pairs))
    eps = np.diag(feff)
    occ = np.diag(D)
    h = np.array([2.0 * (occ[q] - occ[p]) * (eps[p] - eps[q]) for (p, q) in pairs])
    return np.maximum(h, 1.0e-3)


def _trah_boundary(p, d, delta):
    """Positive root ``tau`` of ``|p + tau d| = delta``."""
    a = float(d @ d)
    b = 2.0 * float(p @ d)
    c = float(p @ p) - delta * delta
    disc = max(b * b - 4.0 * a * c, 0.0)
    return (-b + math.sqrt(disc)) / (2.0 * a)


def _trah_steihaug(hv_fn, g, hdiag, delta, nmic):
    """Steihaug-Toint preconditioned CG for the trust-region subproblem."""
    n = g.size
    p = np.zeros(n)
    hp = np.zeros(n)
    r = g.copy()
    m = np.maximum(hdiag, 1.0e-6)
    y = r / m
    d = -y
    ry = float(r @ y)
    rnorm0 = float(np.linalg.norm(r))
    used = 0
    for k in range(1, max(1, nmic) + 1):
        used = k
        hd = hv_fn(d)
        curv = float(d @ hd)
        if curv <= 0.0:                        # negative curvature -> boundary
            tau = _trah_boundary(p, d, delta)
            p = p + tau * d
            hp = hp + tau * hd
            break
        alpha = ry / curv
        if float(p @ p) + 2.0 * alpha * float(p @ d) \
                + alpha * alpha * float(d @ d) >= delta * delta:
            tau = _trah_boundary(p, d, delta)
            p = p + tau * d
            hp = hp + tau * hd
            break
        p = p + alpha * d
        hp = hp + alpha * hd
        r = r + alpha * hd
        rnorm = float(np.linalg.norm(r))
        if rnorm <= min(0.1, math.sqrt(rnorm0)) * rnorm0 or rnorm < 1.0e-10:
            break
        y = r / m
        ry_new = float(r @ y)
        beta = ry_new / ry
        d = -y + beta * d
        ry = ry_new
    pred = -(float(g @ p) + 0.5 * float(p @ hp))
    return p, pred, used


def _trah_lowest_eig(hv_fn, hdiag, n, nmic, seed=987654321):
    """Lowest eigenpair of the orbital Hessian by matrix-free Davidson.

    The stability check: at a point where |g| ~ 0 a negative eigenvalue means a
    saddle, and its mode is the escape direction.  Random starts (fixed seed, so
    the run is reproducible) so a symmetry-breaking mode is found even where the
    gradient cannot see it."""
    mmax = min(max(int(nmic), 8), 40)
    rng = np.random.RandomState(seed)
    V = np.zeros((n, mmax))
    W = np.zeros((n, mmax))
    m = 0
    for _ in range(min(4, mmax)):
        tc = 2.0 * rng.random_sample(n) - 1.0
        for kk in range(m):
            tc = tc - float(V[:, kk] @ tc) * V[:, kk]
        nv = float(np.linalg.norm(tc))
        if nv > 1.0e-8:
            V[:, m] = tc / nv
            m += 1
    theta = 0.0
    u = np.zeros(n)
    mw = 0
    for _ in range(mmax):
        while mw < m:
            W[:, mw] = hv_fn(V[:, mw])
            mw += 1
        Tm = V[:, :m].T @ W[:, :m]
        w, Z = _symmetric_eigh(Tm)
        theta = float(w[0])
        u = V[:, :m] @ Z[:, 0]
        r = W[:, :m] @ Z[:, 0] - theta * u
        if float(np.linalg.norm(r)) < 1.0e-5 or m == mmax:
            break
        den = hdiag - theta
        den = np.where(np.abs(den) < 1.0e-6, np.sign(den) * 1.0e-6 + 1.0e-12, den)
        tc = r / den
        for i in range(m):
            tc = tc - float(V[:, i] @ tc) * V[:, i]
        nv = float(np.linalg.norm(tc))
        if nv < 1.0e-8:
            break
        V[:, m] = tc / nv
        m += 1
    nu = float(np.linalg.norm(u))
    if nu > 0.0:
        u = u / nu
    return theta, u


def _trah_optimize(C, raw_evaluate, evaluate, pairs, nbf, options, params, trace):
    """``[casscf] converger = trah``: matrix-free trust-region augmented Hessian.

    The NumPy pin for, and fallback of, ``casscf_driver.F90``'s ``trah_optimize``
    over the shared core in ``source/trah_core.F90``: the same macro loop, the
    same Steihaug-Toint CG, the same stability check and the same
    Hessian-vector product.  The orbital Hessian is never assembled, so
    ``casscf_hess_bmat`` -- the dominant cost of an assembled-Hessian run -- is
    never called; each CG iteration costs two CI solves instead."""
    tol_g = options.gradient_norm_tol
    npar = len(pairs)
    nmic = max(1, params.max_micro)
    delta = params.start_trust
    # the core's ceiling, not `max_rotation_norm`; see the note in
    # casscf_driver.F90::trah_optimize
    dmax = max(4.0, 8.0 * delta)

    objective, grad, energies, coeffs = evaluate(C)
    state = raw_evaluate.state
    gmat = _trah_gradient_matrix(state["F"])
    hdiag = _trah_hdiag(state, pairs)
    history = [(0, objective, 0.0, float(np.linalg.norm(grad)), 0.0)]
    converged = False
    nmatvec = [0]
    it = 0

    def hv(v):
        nmatvec[0] += 1
        return _trah_hess_vec(C, evaluate, pairs, nbf, v, gmat)

    for it in range(1, options.max_macro_iterations + 1):
        gnorm = float(np.linalg.norm(grad))
        step, pred, _used = _trah_steihaug(hv, grad, hdiag, delta, nmic)
        snorm = float(np.linalg.norm(step))

        # converged only at a genuine minimum: small gradient AND no escape step
        if gnorm < tol_g and snorm < _TRAH_STAB_STEP:
            lam, vmin = _trah_lowest_eig(hv, hdiag, npar, nmic)
            if lam < -_TRAH_STAB_EIG_TOL:
                # A saddle escape is a step, not a converger event: the native
                # driver returns only the history rows and the counters, so the
                # trace must stay identical between the two arms.
                C = _orbital_rotate(C, 0.1 * vmin, pairs, nbf)
                objective, grad, energies, coeffs = evaluate(C)
                state = raw_evaluate.state
                gmat = _trah_gradient_matrix(state["F"])
                hdiag = _trah_hdiag(state, pairs)
                delta = params.start_trust
                continue
            converged = True
            break

        # The model can no longer predict a meaningful energy reduction: near a
        # minimum pred ~ 0.5 g.p ~ gnorm^2, so it underflows the floor while
        # gnorm is still ~1e-7.  The energy is converged and the orbitals are as
        # tight as floating point allows.  (The Fortran core keeps stepping here
        # to tighten the gradient quadratically before it stops; this reference
        # pins the ENERGY, not the trajectory, so it simply stops.)
        if pred <= _TRAH_PRED_FLOOR and gnorm < _TRAH_GTOL_FP \
                and snorm < _TRAH_STAB_STEP:
            converged = True
            break

        obj_old = objective
        Cn = _orbital_rotate(C, step, pairs, nbf)
        obj_new, grad_new, en_new, co_new = evaluate(Cn)
        rho = (obj_old - obj_new) / pred if pred > 0.0 else -1.0
        accepted = rho > 0.1

        # halving line search along the rejected direction
        if not accepted and pred > _TRAH_PRED_FLOOR:
            fac = 0.5
            for _ls in range(5):
                Ct = _orbital_rotate(C, fac * step, pairs, nbf)
                o2, g2, e2, c2 = evaluate(Ct)
                if o2 < obj_old - 1.0e-12:
                    step = fac * step
                    snorm = fac * snorm
                    Cn, obj_new, grad_new, en_new, co_new = Ct, o2, g2, e2, c2
                    rho = (obj_old - obj_new) / (pred * fac)
                    accepted = True
                    break
                fac *= 0.5

        if accepted:
            C, objective, grad, energies, coeffs = Cn, obj_new, grad_new, en_new, co_new
            state = raw_evaluate.state
            gmat = _trah_gradient_matrix(state["F"])
            hdiag = _trah_hdiag(state, pairs)
            history.append((it, objective, objective - obj_old,
                            float(np.linalg.norm(grad)), snorm))

        if rho < _RHO_SHRINK:
            delta = 0.25 * delta
        elif rho > _RHO_GROW and snorm > 0.8 * delta:
            delta = min(_GROW_FACTOR * delta, dmax)

        if delta < _TRAH_DELTA_MIN:
            converged = gnorm < _TRAH_GTOL_FP
            break

    return C, energies, coeffs, history, converged, it


def _diis_coefficients(gradients):
    """Pulay coefficients minimizing |sum_i c_i g_i| with sum c_i = 1.

    Drops the oldest vectors while the bordered B matrix is ill-conditioned.
    Returns the coefficient vector for the *last* ``len(coef)`` stored entries,
    or None when no stable extrapolation exists.

    Fortran engine first; the NumPy implementation below is the numerical pin
    and the fallback.  The engine forms the Gram matrix once and slices a
    trailing sub-block per retry, rather than rebuilding B from scratch each
    time the oldest vector is dropped."""
    backend = _ah_backend("casscf_diis_coeffs")
    if backend is not None:
        lib, ffi = backend
        gm = np.ascontiguousarray(np.asarray(gradients, dtype=np.float64))
        if gm.ndim == 2 and gm.shape[0] >= 2 and gm.shape[1] > 0:
            nvec, npar = gm.shape
            coef = np.zeros(nvec, dtype=np.float64)
            nused = np.zeros(1, dtype=np.int32)
            lib.casscf_diis_coeffs(
                int(nvec), int(npar),
                ffi.cast("double *", gm.ctypes.data), float(_DIIS_CONDMAX),
                ffi.cast("double *", coef.ctypes.data),
                ffi.cast("int32_t *", nused.ctypes.data))
            n = int(nused[0])
            return coef[:n] if n > 0 else None
    return _diis_coefficients_reference(gradients)


def _diis_coefficients_reference(gradients):
    """NumPy reference for :func:`_diis_coefficients` (the numerical pin)."""
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
            if not np.all(np.isfinite(B)) or np.linalg.cond(B) > _DIIS_CONDMAX:
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
            Cn = _orbital_rotate(C, accepted_step, pairs, nbf)
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
                C_x = _orbital_rotate(C0, T_x, pairs, nbf)
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
