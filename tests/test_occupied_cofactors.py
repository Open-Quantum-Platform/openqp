"""Installed mrsf_tlf cofactors against signed minors and its previous routine."""

import ctypes
from decimal import Decimal, localcontext
from itertools import permutations

import numpy as np
import pytest
import oqp.runtime as runtime


def _cofactors(s):
    """Independent Leibniz determinants, with no LU or production determinant."""
    n = len(s)
    result = np.empty((n, n))
    with localcontext() as ctx:
        ctx.prec = 80
        for i in range(n):
            for j in range(n):
                a = [[Decimal(float(s[r, c])) for c in range(n) if c != j]
                     for r in range(n) if r != i]
                value = Decimal(0)
                for p in permutations(range(n - 1)):
                    inversions = sum(p[r] > p[c] for r in range(n - 1)
                                     for c in range(r + 1, n - 1))
                    term = Decimal(-1 if inversions % 2 else 1)
                    for r, c in enumerate(p):
                        term *= a[r][c]
                    value += term
                result[i, j] = float(value * (-1 if (i + j) % 2 else 1))
    return result


def _call(s, *, force_batch=False):
    root, source = runtime.resolve_oqp_root()
    lib = ctypes.CDLL(runtime.library_path(root, source))
    # Both entry points are mandatory: missing symbols fail rather than skip.
    fn = (lib.occupied_cofactors_batch_selftest if force_batch
          else lib.occupied_cofactors_selftest)
    p = ctypes.POINTER(ctypes.c_double)
    fn.argtypes = [ctypes.c_int64, p, p, p, ctypes.POINTER(ctypes.c_int)]
    fn.restype = None
    s = np.array(s, dtype=np.float64, order="F")
    original = s.copy()
    actual, previous = np.empty_like(s), np.empty_like(s)
    used = ctypes.c_int()
    fn(len(s), s.ctypes.data_as(p), actual.ctypes.data_as(p),
       previous.ctypes.data_as(p), ctypes.byref(used))
    assert np.isfinite(actual).all()
    np.testing.assert_array_equal(s, original)
    return actual, previous, bool(used.value)


@pytest.mark.parametrize("s,batch", [
    ([[0.0]], False),
    ([[2.0]], False),
    ([[0.0, 0.7], [0.4, 0.2]], True),
    ([[1.0, 0.2, -0.1], [0.3, 0.9, 0.4], [-0.2, 0.1, 1.1]], True),
    ([[1., 2., 3.], [4., 5., 6.], [5., 7., 9.]], False),
    ([[1., 0.2], [0.3, 0.06000000000001]], False),
    ([[1e-200, 0.0], [0.0, 1e-200]], False),
    ([[1e-310, 0.2], [0.0, 1.0]], False),
    ([[1.0, 0.0], [1e-310, 1.0]], False),
    ([[0.0, 0.2], [0.0, 1.0]], False),
    ([[1.0, 0.0], [0.0, 2e-8]], False),
    ([[1.0, 0.0], [0.0, 1e-8]], False),
    ([[1e-153, 0.0], [0.0, 1e-153]], True),
    ([[1e-154, 0.0], [0.0, 1e-154]], False),
    ([[2.0**-511, 0.0], [0.0, 2.0**-511]], True),
    ([[2.0**-512, 0.0], [0.0, 2.0**-512]], False),
    ([[-1.0, 0.2], [0.3, 0.9]], True),
    ([[2e154, 0.0], [0.0, 2e154]], False),
    ([[2.0, 0.0], [0.0, 2.0]], True),
    ([[2.1, 0.0], [0.0, 2.1]], False),
    ([[1e100, 0., 0.], [1e-150, 1e100, 0.], [0., 0., 1e100]], False),
    ([[1e-100, 0., 0.], [0., 1e-150, 0.], [0., 0., 1e-160]], False),
])
@pytest.mark.parametrize("force_batch", [False, True], ids=["production", "guarded"])
def test_actual_mrsf_tlf(s, batch, force_batch):
    s = np.array(s)
    actual, previous, used = _call(s, force_batch=force_batch)
    if force_batch:
        assert used is batch
    else:
        assert not used
        np.testing.assert_array_equal(actual, previous)
    # The original comp_det classifies a 1x1 subnormal pivot as zero; retain
    # that convention, while the finite O(1) cofactors must still survive.
    expected = _cofactors(s)
    if len(s) == 2:
        expected[np.abs(expected) <= np.finfo(float).tiny] = 0.0
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=1e-315)
    np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-315)


def test_repeated_dimensions_and_nonsymmetric_occupied_overlap():
    for n in (6, 2, 0, 5, 1, 3, 6):
        s = np.eye(n) + np.arange(n*n).reshape(n, n) * 0.002
        actual, previous, used = _call(s, force_batch=True)
        assert used is (n > 1)
        np.testing.assert_allclose(actual, _cofactors(s), rtol=2e-13, atol=1e-15)
        np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-15)


def test_direct_wrappers_return_pivots_and_solve():
    """The mrsf_tlf row-exchange case must differ from a transposed answer."""
    actual, previous, used = _call([[0.0, 0.7], [0.4, 0.2]], force_batch=True)
    assert used
    np.testing.assert_allclose(actual, [[0.2, -0.4], [-0.7, 0.0]], atol=1e-15)
    np.testing.assert_allclose(actual, previous, atol=1e-15)


@pytest.mark.parametrize("rcond,batch", [(0.016, False), (0.020, False), (0.024, True)])
def test_accuracy_guard_boundary_on_nonsymmetric_overlap(rcond, batch):
    # For this triangular matrix, ||S||_1=1 and ||S^-1||_1=1.25/d,
    # hence rcond=0.8*d. At 0.020 the per-entry floor also rejects every
    # cofactor; at 0.024 at least one cofactor retains the LU result.
    s = np.array([[1.0, 0.25], [0.0, rcond / 0.8]])
    actual, previous, used = _call(s, force_batch=True)
    assert used is batch
    np.testing.assert_allclose(actual, _cofactors(s), rtol=2e-13, atol=1e-315)
    np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-315)
    if not batch:
        np.testing.assert_array_equal(actual, previous)


def test_near_dependent_nonsymmetric_overlap_uses_exact_minors():
    # The former sqrt(eps) guard admitted this rcond=2e-8 input, although
    # its eps/rcond estimate exceeds the cofactor tolerance by four orders.
    s = np.array([[1.0, 0.25], [0.0, 2.5e-8]])
    actual, previous, used = _call(s, force_batch=True)
    assert not used
    np.testing.assert_array_equal(actual, previous)
    np.testing.assert_allclose(actual, _cofactors(s), rtol=2e-13, atol=1e-315)


def test_accuracy_guard_accounts_for_matrix_dimension():
    s = np.eye(3)
    s[:2, :2] = [[1.0, 0.25], [0.0, 0.030]]
    # The same rcond=0.024 accepted at n=2 exceeds the n=3 error budget.
    _, _, used_small = _call(s[:2, :2], force_batch=True)
    assert used_small
    actual, previous, used = _call(s, force_batch=True)
    assert not used
    np.testing.assert_array_equal(actual, previous)
    np.testing.assert_allclose(actual, _cofactors(s), rtol=2e-13, atol=1e-315)


def test_small_pivot_product_does_not_accumulate_logarithm_error():
    n = 20
    pivot = 1e-15
    actual, previous, used = _call(pivot * np.eye(n), force_batch=True)
    assert used
    with localcontext() as ctx:
        ctx.prec = 80
        diagonal = float(Decimal(pivot) ** (n - 1))
    np.testing.assert_allclose(np.diag(actual), diagonal, rtol=2e-13, atol=1e-315)
    np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-315)
    np.testing.assert_array_equal(actual - np.diag(np.diag(actual)), np.zeros((n, n)))


def test_tiny_minor_pivot_preserves_previous_permutation():
    # The old minor ordering produces a tiny pivot even though the natural
    # Leibniz cofactor is subnormal. Keep the existing numerical contract.
    s = np.eye(3)
    s[1, 0] = s[2, 1] = 1e-155
    actual, previous, used = _call(s, force_batch=True)
    assert not used
    assert _cofactors(s)[0, 2] > 0.0
    assert actual[0, 2] == previous[0, 2] == 0.0
    np.testing.assert_array_equal(actual, previous)

    # A normal determinant may scale a subnormal inverse entry back into a
    # normal cofactor. Detect the inverse entry before that multiplication.
    s = 2.0 * np.eye(30)
    s[1, 0] = s[2, 1] = 1e-155
    actual, previous, used = _call(s, force_batch=True)
    assert not used
    assert actual[0, 2] == previous[0, 2] == 0.0
    np.testing.assert_array_equal(actual, previous)
    np.testing.assert_array_equal(np.diag(actual), np.full(30, 2.0**29))


@pytest.mark.parametrize("s", [
    np.array([-5, 0, -3, 1, 5, 1, 2, 5, 1], dtype=float).reshape(3, 3, order="F") / 3,
    np.array([[1.0, 0.2, -0.3], [0.1, 1.0, 1.0], [0.4, 0.2, 0.2]]),
])
@pytest.mark.parametrize("perturbation", [0.0, 2.0**-40, 2.0**-28])
def test_cancellation_sensitive_cofactor_preserves_previous_minor(s, perturbation):
    s = s.copy()
    s[1, 2] += perturbation
    # Both fixtures clear the matrix-wide admission rule. The (1,1) minor
    # cancels exactly, or is small after perturbation, within a sound matrix.
    rcond = 1.0 / np.linalg.cond(s, 1)
    assert rcond > 8 * len(s) * np.finfo(float).eps / 2e-13
    actual, previous, used = _call(s, force_batch=True)
    # Keep the old minor's floating-point convention, including any vendor
    # residue: mathematical zero alone is not its numerical contract.
    assert not used
    assert actual[0, 0] == previous[0, 0]
    np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-315)


@pytest.mark.parametrize("perturbation", [0.0, 2.0**-40, 2.0**-28])
def test_hybrid_lu_and_cancellation_minor(perturbation):
    s = np.array([[1.0, 1.0, 0.0], [0.0, 0.5, 0.7], [1.0, 0.5, 0.7]])
    s[1, 2] += perturbation
    actual, previous, used = _call(s, force_batch=True)
    assert used
    assert actual[0, 0] == previous[0, 0]
    np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-315)


@pytest.mark.parametrize("n", [6, 7])
def test_production_size_dispatch_boundary(n):
    s = np.eye(n) + np.arange(n*n).reshape(n, n) * 0.002
    assert not np.array_equal(s, s.T)
    actual, previous, used = _call(s)
    assert used is (n > 6)
    expected = _cofactors(s)
    assert not np.allclose(expected, expected.T, rtol=2e-13, atol=1e-15)
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=1e-15)
    np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-15)
    if n <= 6:
        np.testing.assert_array_equal(actual, previous)


def test_repeated_production_dispatch_dimensions():
    for n in (6, 7, 6, 0, 1, 7, 2):
        s = np.eye(n) + np.arange(n*n).reshape(n, n) * 0.002
        actual, previous, used = _call(s)
        assert used is (n > 6)
        np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-15)
        if n <= 6:
            np.testing.assert_array_equal(actual, previous)


def test_debug_export_owns_its_temporary_container(monkeypatch):
    monkeypatch.setenv("NAC_DUMP_MINORS", "1")
    retained = []
    for n in (6, 7, 6):
        s = np.eye(n) + np.arange(n*n).reshape(n, n) * 0.002
        actual, previous, used = _call(s)
        assert used is (n > 6)
        assert actual.flags.owndata and previous.flags.owndata
        assert not np.shares_memory(actual, previous)
        np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-15)
        if n <= 6:
            np.testing.assert_array_equal(actual, previous)
        retained.append((actual, actual.copy(), previous, previous.copy()))
    # Every native container is already deleted, and later calls used fresh owners.
    for actual, saved_actual, previous, saved_previous in retained:
        np.testing.assert_array_equal(actual, saved_actual)
        np.testing.assert_array_equal(previous, saved_previous)


def _persistent_call(sizes):
    root, source = runtime.resolve_oqp_root()
    lib = ctypes.CDLL(runtime.library_path(root, source))
    fn = lib.occupied_cofactors_persistent_selftest  # mandatory export
    p = ctypes.POINTER(ctypes.c_double)
    pi64 = ctypes.POINTER(ctypes.c_int64)
    pi = ctypes.POINTER(ctypes.c_int)
    fn.argtypes = [ctypes.c_int64, ctypes.c_int64, pi64, p, p, p, pi, pi]
    fn.restype = None
    sizes = np.array(sizes, dtype=np.int64)
    nmax = int(max(sizes))
    inputs = np.zeros((nmax, nmax, len(sizes)), dtype=np.float64, order="F")
    for step, n in enumerate(sizes):
        inputs[:n, :n, step] = np.eye(n) + np.arange(n*n).reshape(n, n) * 0.002
    saved = inputs.copy()
    actual = np.full_like(inputs, np.nan, order="F")
    previous = np.full_like(inputs, np.nan, order="F")
    checks = np.full(6, -1, dtype=np.int32)
    status = ctypes.c_int(-1)
    fn(nmax, len(sizes), sizes.ctypes.data_as(pi64), inputs.ctypes.data_as(p),
       actual.ctypes.data_as(p), previous.ctypes.data_as(p),
       checks.ctypes.data_as(pi), ctypes.byref(status))
    np.testing.assert_array_equal(inputs, saved)
    return inputs, actual, previous, checks, status.value


@pytest.mark.parametrize("sizes", [(6, 7, 3, 6), (7, 3, 6, 7)])
def test_debug_records_replace_in_one_owner_and_recover(monkeypatch, sizes):
    monkeypatch.setenv("NAC_DUMP_MINORS", "1")
    inputs, actual, previous, checks, status = _persistent_call(sizes)
    assert status == 0, f"persistent cofactor lifetime check failed: status={status}, checks={checks}"
    np.testing.assert_array_equal(checks, [1, len(sizes), 4*len(sizes), len(sizes), len(sizes), 2])
    assert actual.flags.owndata and previous.flags.owndata
    assert not np.shares_memory(actual, previous)
    for step, n in enumerate(sizes):
        # Checked after all replacements and both native owner deletions.
        expected = np.zeros_like(actual[:, :, step])
        expected[:n, :n] = _cofactors(inputs[:n, :n, step])
        np.testing.assert_allclose(actual[:, :, step], expected, rtol=2e-13, atol=1e-15)
        np.testing.assert_allclose(actual[:, :, step], previous[:, :, step], rtol=2e-13, atol=1e-15)
        if n <= 6:
            np.testing.assert_array_equal(actual[:, :, step], previous[:, :, step])


def test_persistent_debug_invalid_dimension_reports_failure_and_cleans(monkeypatch):
    monkeypatch.setenv("NAC_DUMP_MINORS", "1")
    _, actual, previous, checks, status = _persistent_call((6, 1, 7))
    assert status == 1
    np.testing.assert_array_equal(checks, [0, 0, 0, 0, 0, 2])
    np.testing.assert_array_equal(actual, np.zeros_like(actual))
    np.testing.assert_array_equal(previous, np.zeros_like(previous))


def test_persistent_debug_missing_records_cleans_live_owner(monkeypatch):
    monkeypatch.delenv("NAC_DUMP_MINORS", raising=False)
    inputs, actual, previous, checks, status = _persistent_call((6, 7))
    assert status == 4  # one computed call, then the absent-record contract fails
    np.testing.assert_array_equal(checks, [1, 1, 0, 0, 0, 2])
    np.testing.assert_array_equal(actual[:, :, 1], np.zeros((7, 7)))
    np.testing.assert_array_equal(previous[:, :, 1], np.zeros((7, 7)))
    np.testing.assert_array_equal(actual[:, :, 0], previous[:, :, 0])
    np.testing.assert_allclose(actual[:6, :6, 0], _cofactors(inputs[:6, :6, 0]),
                               rtol=2e-13, atol=1e-15)


@pytest.mark.parametrize("control", ["unrecognized", "missing-record," + "x"*40])
def test_persistent_debug_control_rejects_invalid_selector(monkeypatch, control):
    monkeypatch.setenv("NAC_DUMP_MINORS", "1")
    monkeypatch.setenv("OQP_COFACTOR_LIFETIME_CONTROL", control)
    _, actual, previous, checks, status = _persistent_call((6, 7, 3, 6))
    assert status == 8
    np.testing.assert_array_equal(checks, [0, 0, 0, 0, 0, 2])
    np.testing.assert_array_equal(actual, np.zeros_like(actual))
    np.testing.assert_array_equal(previous, np.zeros_like(previous))
