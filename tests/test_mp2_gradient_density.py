"""Pin the production contraction with non-symmetric, signed operands.

The reference constructs both MO blocks and applies the original four-index
transformation, independently of the occupied-index batching in production.
Missing native support is an error; this test must never pass by skipping.
"""
import numpy as np
import oqp


def _dense_reference(c, l2):
    n = c.shape[0]
    no = l2.shape[0]
    gmo = np.zeros((n, n, n, n))
    gmo[:no, no:, :no, no:] = 2 * l2.transpose(0, 2, 1, 3)
    gmo[no:, :no, no:, :no] = 2 * l2.transpose(2, 0, 3, 1)
    g = np.einsum("mp,nq,lr,ks,pqrs->mnlk", c, c, c, c, gmo,
                  optimize=True)
    permutations = ((0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2),
                    (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 0, 1),
                    (2, 3, 1, 0), (3, 2, 1, 0))
    return sum(g.transpose(p) for p in permutations) / 8


def test_compact_density_matches_dense_and_resizes_in_one_process():
    rng = np.random.default_rng(710470)
    first = None
    # Grow and shrink n, no, and nv. No symmetry is imposed on c or l2.
    for n, no in ((3, 1), (7, 3), (4, 2), (7, 1), (3, 0), (3, 3), (3, 1)):
        if (n, no) == (3, 1) and first is not None:
            c, l2, prior = first
        else:
            c = np.asfortranarray(rng.normal(size=(n, n)) / n)
            l2 = np.asfortranarray(rng.normal(size=(no, no, n-no, n-no)))
        g = np.empty((n, n, n, n), order="F")
        oqp.lib.mp2_gradient_density_test(
            n, no, oqp.ffi.cast("double *", c.ctypes.data),
            oqp.ffi.cast("double *", l2.ctypes.data),
            oqp.ffi.cast("double *", g.ctypes.data))
        np.testing.assert_allclose(g, _dense_reference(c, l2),
                                   atol=2e-14, rtol=2e-14)
        if (n, no) == (3, 1):
            if first is None:
                first = c.copy(order="F"), l2.copy(order="F"), g.copy()
            else:
                np.testing.assert_array_equal(g, prior)
