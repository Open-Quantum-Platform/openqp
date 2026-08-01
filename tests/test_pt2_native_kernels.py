"""Equivalence pin for the Fortran PT2 engines (pt2_kernel.F90, nevpt2_koopmans.F90).

These kernels replace arithmetic that used to run in Python/NumPy on the PT2
execution path:

* ``pt2_kernel.F90``       -- the mean-field Fock build, the Dyall active
                              dressing, and the determinant-space bookkeeping
                              (external space, diagonal H0, occupation blocks)
                              of ``caspt2_dyall.py``;
* ``nevpt2_koopmans.F90``  -- the two eri-folded 4-pdm intermediates and the
                              ``a16`` Koopmans intermediate of ``nevpt2_sc.py``.

Each Python implementation stays in the tree precisely so it can serve as this
pin, and as the fallback when liboqp predates the symbol.  Every test below
evaluates the same function twice -- once through liboqp, once with the backend
probe forced to ``None`` -- so the two paths are compared directly rather than
against a hard-coded constant.

The integer-valued results (external indices, occupation blocks) must agree
EXACTLY: the block partition feeds a per-block eigendecomposition, so a
different intra-block member order would perturb the printed energies even
though it is mathematically equivalent.  The floating-point results are held to
round-off rather than to a loose tolerance that would hide a real discrepancy.
"""
import numpy as np
import pytest

import oqp.library.caspt2_dyall as cp
import oqp.library.nevpt2_sc as sc


def _has(name, probe):
    backend = probe()
    return backend is not None and hasattr(backend[0], name)


def _without_backend(mod, probe_name, fn):
    """Evaluate ``fn`` with the module's native-backend probe disabled."""
    saved = getattr(mod, probe_name)
    setattr(mod, probe_name, lambda: None)
    try:
        return fn()
    finally:
        setattr(mod, probe_name, saved)


def _both(mod, probe_name, fn):
    return fn(), _without_backend(mod, probe_name, fn)


def _symmetric_eri(nbf, rng):
    """Random ERIs carrying the full 8-fold permutational symmetry."""
    a = rng.standard_normal((nbf, nbf, nbf, nbf)) * 0.1
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    return a + a.transpose(2, 3, 0, 1)


# --------------------------------------------------------------------- pt2_kernel
pt2_only = pytest.mark.skipif(
    not _has("pt2_effective_fock", cp._pt2_lib),
    reason="liboqp pt2_kernel engine unavailable",
)


@pt2_only
@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (12, 3, 6)])
def test_effective_fock_matches_python(nbf, ncore, nact):
    rng = np.random.default_rng(11 + nbf)
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _symmetric_eri(nbf, rng)
    D = rng.standard_normal((nbf, nbf))
    D = D + D.T

    native, python = _both(cp, "_pt2_lib", lambda: cp._effective_fock(h1e, eri, D))
    assert native.shape == python.shape == (nbf, nbf)
    assert np.max(np.abs(native - python)) < 1e-11


@pt2_only
@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (12, 3, 6)])
def test_h0_integrals_dyall_matches_python(nbf, ncore, nact):
    rng = np.random.default_rng(23 + nbf)
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _symmetric_eri(nbf, rng)
    eps = rng.standard_normal(nbf)

    native, python = _both(
        cp, "_pt2_lib",
        lambda: cp._h0_integrals(h1e, eri, eps, ncore, nact, "dyall")[0])
    assert np.max(np.abs(native - python)) < 1e-11
    # the inactive/virtual diagonal is untouched by the dressing
    assert np.allclose(np.diag(native)[ncore + nact:], eps[ncore + nact:])


@pt2_only
@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (10, 2, 4)])
def test_external_indices_match_python_exactly(nbf, ncore, nact):
    from oqp.library.fci import _determinants

    dets = _determinants(nbf, (ncore + 2, ncore + 2))
    native, python = _both(
        cp, "_pt2_lib", lambda: cp._external_indices(dets, ncore, nact, nbf))
    assert np.array_equal(np.asarray(native), np.asarray(python))
    # the CAS configurations are exactly the complement
    assert 0 < len(native) < len(dets)


@pt2_only
@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (10, 2, 4)])
def test_diagonal_zeroth_order_matches_python(nbf, ncore, nact):
    from oqp.library.fci import _determinants

    rng = np.random.default_rng(31 + nbf)
    dets = _determinants(nbf, (ncore + 2, ncore + 2))
    eps = rng.standard_normal(nbf)
    h0_1e = np.diag(eps)
    g0 = np.zeros((nbf,) * 4)

    native, python = _both(
        cp, "_pt2_lib",
        lambda: cp._diagonal_zeroth_order(dets, h0_1e, g0, nbf))
    assert native is not None and python is not None
    assert np.max(np.abs(native - python)) < 1e-11


@pt2_only
@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (10, 2, 4)])
def test_occupation_blocks_match_python_exactly(nbf, ncore, nact):
    from oqp.library.fci import _determinants

    dets = _determinants(nbf, (ncore + 2, ncore + 2))
    ext = cp._external_indices(dets, ncore, nact, nbf)

    native, python = _both(
        cp, "_pt2_lib",
        lambda: cp._occupation_blocks(dets, ext, ncore, nact, nbf))
    assert len(native) == len(python)
    for a, b in zip(native, python):
        # exact, including the member order inside each block
        assert np.array_equal(np.asarray(a), np.asarray(b))
    # the blocks partition the external space
    assert sum(len(g) for g in native) == len(ext)
    assert np.array_equal(np.sort(np.concatenate(native)), np.arange(len(ext)))


# ---------------------------------------------------------------- Koopmans
koopmans_only = pytest.mark.skipif(
    not _has("nevpt2_f3ca_f3ac", sc._koopmans_lib),
    reason="liboqp nevpt2_koopmans engine unavailable",
)


def _koopmans_operands(n, seed):
    rng = np.random.default_rng(seed)
    h1e = rng.standard_normal((n, n))
    h1e = 0.5 * (h1e + h1e.T)
    h2e = rng.standard_normal((n, n, n, n)) * 0.1
    h2e = h2e + h2e.transpose(3, 2, 1, 0)      # real physicist-ordered ERIs
    dm3 = rng.standard_normal((n,) * 6) * 0.1
    dm4 = rng.standard_normal((n,) * 8) * 0.1
    return h1e, h2e, dm3, dm4


@koopmans_only
@pytest.mark.parametrize("n", [3, 4, 5])
def test_f3ca_f3ac_match_python(n):
    _h1e, h2e, _dm3, dm4 = _koopmans_operands(n, 41 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._f3ca_f3ac(h2e, dm4))
    for a, b in zip(native, python):
        assert a.shape == (n,) * 6
        assert np.max(np.abs(a - b)) < 1e-11


@koopmans_only
@pytest.mark.parametrize("n", [3, 4, 5])
def test_a16_matches_python(n):
    h1e, h2e, dm3, dm4 = _koopmans_operands(n, 53 + n)
    f3ca, f3ac = sc._f3ca_f3ac(h2e, dm4)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._a16(h1e, h2e, dm3, f3ca, f3ac))
    assert native.shape == (n,) * 6
    assert np.max(np.abs(native - python)) < 1e-11
