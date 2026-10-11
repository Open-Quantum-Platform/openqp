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


def _call(s):
    root, source = runtime.resolve_oqp_root()
    lib = ctypes.CDLL(runtime.library_path(root, source))
    fn = lib.occupied_cofactors_selftest  # Missing shim is a failure, never a skip.
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
    ([[1.0, 0.0], [0.0, 2e-8]], True),
    ([[1.0, 0.0], [0.0, 1e-8]], False),
    ([[1e-153, 0.0], [0.0, 1e-153]], True),
    ([[1e-154, 0.0], [0.0, 1e-154]], False),
    ([[2e154, 0.0], [0.0, 2e154]], False),
    ([[1e-100, 0., 0.], [0., 1e-150, 0.], [0., 0., 1e-160]], False),
])
def test_actual_mrsf_tlf(s, batch):
    s = np.array(s)
    actual, previous, used = _call(s)
    assert used is batch
    # The original comp_det classifies a 1x1 subnormal pivot as zero; retain
    # that convention, while the finite O(1) cofactors must still survive.
    expected = _cofactors(s)
    if len(s) == 2:
        expected[np.abs(expected) <= np.finfo(float).tiny] = 0.0
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=1e-315)
    np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-315)


def test_repeated_dimensions_and_nonsymmetric_occupied_overlap():
    for n in (6, 2, 5, 1, 3, 6):
        s = np.eye(n) + np.arange(n*n).reshape(n, n) * 0.002
        actual, previous, used = _call(s)
        assert used is (n != 1)
        np.testing.assert_allclose(actual, _cofactors(s), rtol=2e-13, atol=1e-15)
        np.testing.assert_allclose(actual, previous, rtol=2e-13, atol=1e-15)


def test_direct_wrappers_return_pivots_and_solve():
    """The mrsf_tlf row-exchange case must differ from a transposed answer."""
    actual, previous, used = _call([[0.0, 0.7], [0.4, 0.2]])
    assert used
    np.testing.assert_allclose(actual, [[0.2, -0.4], [-0.7, 0.0]], atol=1e-15)
    np.testing.assert_allclose(actual, previous, atol=1e-15)
