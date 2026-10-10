"""Equivalence tests for the one-call CI driver (``fci_solve``, fci_driver.F90).

``solve_fci`` remains the numerical pin: every native solve must reproduce the
Python driver's energies, ``<S^2>`` values and CI vectors (up to the sign a
CI vector is only defined to) in the *same build*, against the *same*
integrals -- a hard-coded reference constant would only pin the build against
itself.

Also asserts that the fixed-layout option schema in ``fci_driver.F90``, the
copy in ``include/oqp.h`` and the mirror in ``fci.py`` are the same list, so
the three cannot drift apart silently.

All numerical tests skip when liboqp lacks ``fci_solve`` (pre-rebuild checkouts).
"""
import os
import re

import numpy as np
import pytest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_FORTRAN = os.path.join(_REPO, "source", "modules", "fci_driver.F90")
_HEADER = os.path.join(_REPO, "include", "oqp.h")


def _driver_available() -> bool:
    try:
        import oqp
    except Exception:
        return False
    if not bool(getattr(oqp, "BACKEND_AVAILABLE", False)):
        return False
    from oqp import lib
    return hasattr(lib, "fci_solve")


requires_driver = pytest.mark.skipif(
    not _driver_available(), reason="liboqp has no fci_solve (rebuild needed)"
)


# ------------------------------------------------------------------- schema
def _fortran_schema(prefix):
    """(name, index) pairs for the FCI_<prefix>_* parameters, in file order."""
    with open(_FORTRAN, encoding="utf-8") as handle:
        text = handle.read()
    pattern = re.compile(
        rf"^\s*integer, parameter :: FCI_{prefix}_(\w+)\s*=\s*(\d+)", re.M)
    return [(m.group(1).lower(), int(m.group(2))) for m in pattern.finditer(text)]


def _header_schema(prefix):
    with open(_HEADER, encoding="utf-8") as handle:
        text = handle.read()
    pattern = re.compile(rf"^\s*FCI_{prefix}_(\w+)\s*=\s*(\d+)", re.M)
    return [(m.group(1).lower(), int(m.group(2))) for m in pattern.finditer(text)]


def test_released_option_lengths_are_frozen():
    """v1.3.0/v1.3.1 published FCI_NIOPT=14, FCI_NDOPT=3.

    Those are the lengths a binary compiled against the released oqp.h
    allocates, and `fci_solve` is the symbol it calls.  Reading iopt[14],
    iopt[15] or dopt[3] there is an out-of-bounds read of the caller's memory,
    and adjacent nonzero memory reads as an irrep request.  The three copies of
    that promise -- Fortran, header, Python -- have to agree and must not move.
    """
    from oqp.library.fci import _FCI_NDOPT_V1, _FCI_NIOPT_V1

    assert (_FCI_NIOPT_V1, _FCI_NDOPT_V1) == (14, 3)
    for path in (_FORTRAN, _HEADER):
        text = open(path, encoding="utf-8").read()
        niopt = int(re.search(r"FCI_NIOPT_V1\s*=\s*(\d+)", text).group(1))
        ndopt = int(re.search(r"FCI_NDOPT_V1\s*=\s*(\d+)", text).group(1))
        assert (niopt, ndopt) == (14, 3), f"{path} moved the released lengths"


@requires_driver
def test_released_entry_point_reads_only_the_released_layout():
    """`fci_solve` still works when handed exactly 14 iopt and 3 dopt entries.

    This is the call an unrecompiled v1.3.x consumer makes.  It must return the
    same roots as the extended entry point given the same options -- if the old
    symbol had grown a dependency on the new slots, the two would disagree (or
    the short call would read past its arrays).
    """
    import numpy as np
    from oqp import ffi, lib
    from oqp.library import fci as F

    plan, settings = _plan(4, (2, 2), frozen_core=0, active_orbitals=4,
                           nroot=3, solver="dense")
    h, v = _system(4, seed=11)
    reference = F.solve_active_ci(
        h, v, plan, 0.0, settings, nroot=3, want_s2=False,
        use_target_spin=False)
    assert reference is not None
    ref_energies = reference[0]

    spec = F.resolve_ci_solve(
        plan.nact, plan.nelec, ecore=0.0, nroot=3,
        max_det=settings.max_det, max_memory=settings.max_memory,
        eig_tol=settings.eig_tol, solver="dense")

    # exactly the released lengths -- nothing allocated beyond them
    iopt = np.zeros(F._FCI_NIOPT_V1, dtype=np.int32)
    dopt = np.zeros(F._FCI_NDOPT_V1, dtype=np.float64)
    iopt[F._FCI_IOPT_INDEX["norb"]] = plan.norb
    iopt[F._FCI_IOPT_INDEX["nact"]] = plan.nact
    iopt[F._FCI_IOPT_INDEX["ncore"]] = plan.ncore
    iopt[F._FCI_IOPT_INDEX["nalpha"]] = spec.nalpha
    iopt[F._FCI_IOPT_INDEX["nbeta"]] = spec.nbeta
    iopt[F._FCI_IOPT_INDEX["nroot"]] = spec.nroot
    iopt[F._FCI_IOPT_INDEX["solver"]] = F._FCI_SOLVER_CODE["dense"]
    iopt[F._FCI_IOPT_INDEX["maxiter"]] = spec.davidson_maxiter
    iopt[F._FCI_IOPT_INDEX["maxmemory"]] = spec.max_memory
    iopt[F._FCI_IOPT_INDEX["nthreads"]] = 1
    dopt[F._FCI_DOPT_INDEX["eig_tol"]] = spec.eig_tol

    active = np.ascontiguousarray(plan.active, dtype=np.int32)
    core = (np.ascontiguousarray(plan.core, dtype=np.int32) if plan.core
            else np.zeros(1, dtype=np.int32))
    h_c = F._as_f64c(h)
    v_c = F._as_f64c(v)
    ndet = len(F._determinants(plan.nact, (spec.nalpha, spec.nbeta)))
    energies = np.zeros(spec.nroot, dtype=np.float64)
    civecs = np.zeros((ndet, spec.nroot), dtype=np.float64)
    s2 = np.zeros(spec.nroot, dtype=np.float64)

    status = int(lib.fci_solve(
        ffi.cast("int32_t *", iopt.ctypes.data),
        ffi.cast("double *", dopt.ctypes.data),
        ffi.cast("int32_t *", active.ctypes.data),
        ffi.cast("int32_t *", core.ctypes.data),
        ffi.cast("double *", h_c.ctypes.data),
        ffi.cast("double *", v_c.ctypes.data),
        ffi.cast("double *", energies.ctypes.data),
        ffi.cast("double *", civecs.ctypes.data),
        ffi.cast("double *", s2.ctypes.data)))

    assert status == spec.nroot
    assert np.allclose(energies, ref_energies, rtol=0.0, atol=1.0e-10)


def test_option_schema_matches_fortran():
    """The Python mirror of the iopt/dopt layout is the Fortran one."""
    from oqp.library.fci import _FCI_DOPT, _FCI_IOPT

    for prefix, mirror in (("I", _FCI_IOPT), ("D", _FCI_DOPT)):
        fortran = _fortran_schema(prefix)
        assert [n for n, _ in fortran] == list(mirror), (
            f"fci.py _FCI_{prefix}OPT does not match fci_driver.F90")
        assert [i for _, i in fortran] == list(range(len(mirror)))


def test_option_schema_matches_header():
    """include/oqp.h documents exactly the Fortran layout."""
    from oqp.library.fci import _FCI_DOPT, _FCI_IOPT

    header_i = [n for n, _ in _header_schema("I")]
    header_d = [n for n, _ in _header_schema("D")]
    assert header_i == list(_FCI_IOPT)
    assert header_d == list(_FCI_DOPT)


def test_schema_counts_are_consistent():
    from oqp.library.fci import _FCI_DOPT, _FCI_IOPT

    text = open(_FORTRAN, encoding="utf-8").read()
    niopt = int(re.search(r"FCI_NIOPT\s*=\s*(\d+)", text).group(1))
    ndopt = int(re.search(r"FCI_NDOPT\s*=\s*(\d+)", text).group(1))
    assert niopt == len(_FCI_IOPT)
    assert ndopt == len(_FCI_DOPT)


# -------------------------------------------------------------- numerics
def _system(norb, seed):
    rng = np.random.default_rng(seed)
    h = rng.normal(size=(norb, norb)) * 0.2
    h = 0.5 * (h + h.T) + np.diag(np.linspace(-4.0, 4.0, norb))
    v = rng.normal(size=(norb, norb, norb, norb)) * 0.02
    v = v + v.transpose(1, 0, 2, 3)
    v = v + v.transpose(0, 1, 3, 2)
    v = v + v.transpose(2, 3, 0, 1)
    return np.ascontiguousarray(h), np.ascontiguousarray(v)


def _both_arms(h, v, plan, settings, ecore, **kwargs):
    """(native, python) results of solve_active_ci in this one build."""
    from oqp.library import fci as F

    native = F.solve_active_ci(h, v, plan, ecore, settings, **kwargs)
    saved = F._fci_solve_backend
    F._fci_solve_backend = lambda: None
    try:
        reference = F.solve_active_ci(h, v, plan, ecore, settings, **kwargs)
    finally:
        F._fci_solve_backend = saved
    return native, reference


def _plan(norb, nelec_total, **kw):
    from oqp.library.fci import FCISettings, active_space_plan

    settings = FCISettings(max_det=100000, max_memory=4096, eig_tol=1.0e-10, **kw)
    return active_space_plan(norb, nelec_total, settings), settings


@requires_driver
@pytest.mark.parametrize(
    "norb,nelec,frozen,nact,nroot,solver",
    [
        (4, (2, 2), 0, 4, 1, "dense"),
        (4, (2, 2), 0, 4, 3, "dense"),
        (6, (3, 3), 0, 6, 1, "dense"),
        (6, (3, 3), 0, 6, 4, "dense"),
        (6, (3, 3), 0, 6, 1, "davidson"),
        (6, (3, 3), 0, 6, 3, "davidson"),
        (8, (5, 5), 2, 6, 2, "dense"),
        (8, (5, 5), 2, 6, 2, "davidson"),
    ],
)
def test_native_solve_matches_python_driver(norb, nelec, frozen, nact, nroot, solver):
    plan, settings = _plan(
        norb, nelec, frozen_core=frozen, active_orbitals=nact, nroot=nroot,
        solver=solver)
    h, v = _system(norb, seed=norb * 17 + nroot)
    (e_n, c_n, s2_n), (e_p, c_p, s2_p) = _both_arms(
        h, v, plan, settings, 1.25,
        nroot=nroot, want_s2=True, use_target_spin=False,
        active_section="[cas]", ci_section="[ci]")

    assert np.allclose(e_n, e_p, rtol=0.0, atol=1.0e-10)
    assert np.allclose(s2_n, s2_p, rtol=0.0, atol=1.0e-9)
    # A CI vector is defined only up to a sign, so both engines impose the SAME
    # canonical phase on the way out -- canonical_phase() in fci_driver.F90 and
    # canonicalize_ci_phase() in fci.py.  The overlap is therefore required to
    # be +1, not |1|: an unsigned check here is what let the two drift apart
    # and let a thread-count change flip the sign of a multistate CASPT2
    # off-diagonal.
    overlap = np.sum(c_n * c_p, axis=0)
    assert np.allclose(overlap, 1.0, atol=1.0e-8)


@requires_driver
@pytest.mark.parametrize("solver", ["dense", "davidson"])
@pytest.mark.parametrize("target", ["singlet", "triplet"])
def test_native_target_spin_matches_python_driver(solver, target):
    plan, settings = _plan(
        6, (3, 3), frozen_core=0, active_orbitals=6, nroot=1, solver=solver,
        target_spin=target)
    h, v = _system(6, seed=99)
    (e_n, c_n, _), (e_p, c_p, _) = _both_arms(
        h, v, plan, settings, 0.0,
        nroot=1, want_s2=False, use_target_spin=True,
        active_section="[cas]", ci_section="[ci]")
    assert np.allclose(e_n, e_p, rtol=0.0, atol=1.0e-10)
    # signed, for the reason above: spin filtering must not lose the phase
    assert np.allclose(np.sum(c_n * c_p, axis=0), 1.0, atol=1.0e-8)


def test_canonical_ci_phase_makes_the_largest_amplitude_positive():
    """The convention, stated directly."""
    from oqp.library.fci import canonicalize_ci_phase

    vecs = np.array([
        [0.1, -0.2],
        [-0.7, 0.3],
        [0.3, -0.9],
    ])
    out = canonicalize_ci_phase(vecs.copy())
    # column 0: |c| peaks at row 1 with -0.7 -> flipped
    assert np.allclose(out[:, 0], [-0.1, 0.7, -0.3])
    # column 1: |c| peaks at row 2 with -0.9 -> flipped
    assert np.allclose(out[:, 1], [0.2, -0.3, 0.9])
    # and applying it again changes nothing
    assert np.allclose(canonicalize_ci_phase(out.copy()), out)


def test_canonical_ci_phase_breaks_a_round_off_tie_by_determinant_index():
    """Symmetry-equivalent amplitudes must not let round-off pick the phase.

    Two determinants related by symmetry carry equal weight in exact
    arithmetic and differ by ~1e-16 in a real solve.  An exact argmax over
    |c| would follow that noise to whichever row happened to come out
    fractionally larger -- and if the two have opposite signs, the vector
    flips.  The tie window has to send both orderings to the same answer.
    """
    from oqp.library.fci import canonicalize_ci_phase

    lo = np.array([[-0.5], [0.5 + 1.0e-16], [0.1]])
    hi = np.array([[-0.5 - 1.0e-16], [0.5], [0.1]])
    out_lo = canonicalize_ci_phase(lo.copy())
    out_hi = canonicalize_ci_phase(hi.copy())
    # row 0 is the lowest index inside the window, and it is negative in both
    assert out_lo[0, 0] > 0.0 and out_hi[0, 0] > 0.0
    assert np.allclose(out_lo, out_hi, atol=1.0e-15)


def test_canonical_ci_phase_leaves_a_zero_vector_alone():
    from oqp.library.fci import canonicalize_ci_phase

    zeros = np.zeros((4, 2))
    assert np.array_equal(canonicalize_ci_phase(zeros.copy()), zeros)


@requires_driver
@pytest.mark.parametrize("solver", ["dense", "davidson"])
def test_native_roots_obey_the_canonical_phase(solver):
    """The engine applies the convention, not just the Python fallback."""
    from oqp.library import fci as F

    plan, settings = _plan(6, (3, 3), frozen_core=0, active_orbitals=6,
                           nroot=3, solver=solver)
    h, v = _system(6, seed=31)
    _energies, civecs, _ = F.solve_active_ci(
        h, v, plan, 0.0, settings, nroot=3, want_s2=False,
        use_target_spin=False, active_section="[cas]", ci_section="[ci]")

    magnitudes = np.abs(civecs)
    window = magnitudes.max(axis=0) * (1.0 - F._CI_PHASE_TIE_RTOL)
    lead = np.argmax(magnitudes >= window, axis=0)
    assert np.all(civecs[lead, np.arange(civecs.shape[1])] > 0.0)


@requires_driver
def test_native_civectors_are_returned_in_c_order():
    """civecs is C-order [ndet, nroot]: root k is column k, contiguous rows.

    A C-order [m, n] buffer is the Fortran (n, m) matrix, so the engine has to
    transpose on the way out; getting that backwards is silent when ndet and
    nroot happen to be equal, so pin it with nroot != ndet and an explicit
    per-root residual against the Hamiltonian the roots came from.
    """
    from oqp.library import fci as F

    plan, settings = _plan(4, (2, 2), frozen_core=0, active_orbitals=4,
                           nroot=3, solver="dense")
    h, v = _system(4, seed=5)
    energies, civecs, _ = F.solve_active_ci(
        h, v, plan, 0.0, settings, nroot=3, want_s2=False,
        use_target_spin=False)
    h_act, eri_act, nelec, ecore = F.apply_active_space(h, v, plan, 0.0)
    dets = F._determinants(plan.nact, nelec)
    hs, gs = F._spin_orbital_integrals(h_act, eri_act)
    ham = F._build_dense_hamiltonian(dets, None, hs, gs, 2 * plan.nact, 0.0)

    assert civecs.shape == (len(dets), 3)
    for root in range(3):
        vec = civecs[:, root]
        residual = ham @ vec - (energies[root] - ecore) * vec
        assert np.linalg.norm(residual) < 1.0e-8


@requires_driver
def test_native_declines_beyond_the_determinant_key_width():
    """An active space wider than a 64-bit determinant key falls back."""
    from oqp.library import fci as F

    plan, settings = _plan(40, (20, 20), frozen_core=8, active_orbitals=32,
                           nroot=1, solver="davidson")
    spec = F.resolve_ci_solve(plan.nact, plan.nelec, nroot=1, max_det=10**18,
                              max_memory=8, solver="davidson")
    assert F._lib_fci_solve(
        np.zeros((40, 40)), np.zeros((40,) * 4), plan, spec,
        nthreads=1, want_s2=False, use_target_spin=False) is None
