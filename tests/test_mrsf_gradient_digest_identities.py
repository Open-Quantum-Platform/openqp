"""Algebraic identities the MRSF two-particle-density digests rely on.

`grd2_mrsf_compute_data_t_get_density` in
`source/modules/tdhf_mrsf_gradient.F90` used to spell out eight `dc`/`dd`
blocks. Most of that was arithmetic the engine threw away:

  * `dc3`/`dc4`/`dd3`/`dd4` hold the same multiset of products as
    `dc1`/`dc2`/`dd1`/`dd2` -- only the order of the two factors differs -- so
    `-dc1-dc2-dc3-dc4+dd1+dd2+dd3+dd4` is `2*(-dc1-dc2+dd1+dd2)`.
  * `dd1`/`dd2` are products of index-symmetric sums, so sixteen multiplies
    collapse to four.

These are exact identities in real arithmetic. The digest now relies on them, so
if anyone re-derives those blocks and the identity stops holding, the engine is
silently wrong rather than merely slow. This test states the identities directly
on random matrices; it needs no build.
"""

import unittest

import numpy as np

SEED = 20260919
NDIM = 8
TOL = 1e-12


def _matrices(count):
    rng = np.random.default_rng(SEED)
    return [rng.standard_normal((NDIM, NDIM)) for _ in range(count)]


def _dc(a, b, i, j, k, l):
    """The eight-group `dc` block, single-state form."""
    return (
        a[i, k] * b[j, l] + a[i, l] * b[j, k] + a[j, k] * b[i, l]
        + a[j, l] * b[i, k] + a[l, j] * b[k, i] + a[k, j] * b[l, i]
        + a[l, i] * b[k, j] + a[k, i] * b[l, j]
    )


def _dd(a, b, i, j, k, l):
    """The eight-group `dd` block, single-state form."""
    return (
        a[i, j] * b[l, k] + a[i, j] * b[k, l] + a[j, i] * b[l, k]
        + a[j, i] * b[k, l] + a[l, k] * b[i, j] + a[k, l] * b[i, j]
        + a[l, k] * b[j, i] + a[k, l] * b[j, i]
    )


def _dd_factored(a, b, i, j, k, l):
    """What the engine computes instead of `_dd`."""
    return (a[i, j] + a[j, i]) * (b[l, k] + b[k, l]) + (
        a[l, k] + a[k, l]
    ) * (b[i, j] + b[j, i])


def _indices():
    for i in range(NDIM):
        for j in range(NDIM):
            for k in range(NDIM):
                for l in range(NDIM):
                    yield i, j, k, l


class DigestIdentities(unittest.TestCase):
    def _assert_close(self, lhs, rhs, scale, what):
        worst = max(abs(a - b) for a, b in zip(lhs, rhs))
        self.assertGreater(scale, 0.0, f"{what}: degenerate test data")
        self.assertLess(
            worst / scale, TOL,
            f"{what}: worst relative deviation {worst / scale:.3e}",
        )

    def test_dc3_repeats_dc1(self):
        """Swapping the two factor matrices reproduces the same block."""
        bco1, bo2v = _matrices(2)
        lhs, rhs, scale = [], [], 0.0
        for idx in _indices():
            a = _dc(bco1, bo2v, *idx)
            lhs.append(a)
            rhs.append(_dc(bo2v, bco1, *idx))
            scale = max(scale, abs(a))
        self._assert_close(lhs, rhs, scale, "dc3 vs dc1")

    def test_dd3_repeats_dd1(self):
        bco1, bo2v = _matrices(2)
        lhs, rhs, scale = [], [], 0.0
        for idx in _indices():
            a = _dd(bco1, bo2v, *idx)
            lhs.append(a)
            rhs.append(_dd(bo2v, bco1, *idx))
            scale = max(scale, abs(a))
        self._assert_close(lhs, rhs, scale, "dd3 vs dd1")

    def test_dd_factorisation(self):
        bco1, bo2v = _matrices(2)
        lhs, rhs, scale = [], [], 0.0
        for idx in _indices():
            a = _dd(bco1, bo2v, *idx)
            lhs.append(a)
            rhs.append(_dd_factored(bco1, bo2v, *idx))
            scale = max(scale, abs(a))
        self._assert_close(lhs, rhs, scale, "dd factorisation")

    def test_db_groups_5_to_8_repeat_1_to_4(self):
        """The single-state db blocks carry the same doubling."""
        co12, = _matrices(1)
        groups = [(0, 2, 3, 1), (0, 3, 2, 1), (1, 2, 3, 0), (1, 3, 2, 0)]
        lhs, rhs, scale = [], [], 0.0
        for idx in _indices():
            pick = lambda t: (idx[t[0]], idx[t[1]])
            first4 = sum(
                co12[pick((g[0], g[1]))] * co12[pick((g[2], g[3]))]
                for g in groups
            )
            full = first4 + sum(
                co12[pick((g[2], g[3]))] * co12[pick((g[0], g[1]))]
                for g in groups
            )
            lhs.append(full)
            rhs.append(2.0 * first4)
            scale = max(scale, abs(full))
        self._assert_close(lhs, rhs, scale, "db doubling")

if __name__ == "__main__":
    unittest.main()
