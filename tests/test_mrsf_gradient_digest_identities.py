"""Algebraic identities the MRSF two-particle-density digests rely on.

`grd2_mrsf_compute_data_t_get_density` and its NAC counterpart in
`source/modules/tdhf_mrsf_gradient.F90` used to spell out eight `dc`/`dd`
blocks and, in the NAC case, eight I<->J-symmetrised groups per channel. Most
of that was arithmetic the engine threw away:

  * `dc3`/`dc4`/`dd3`/`dd4` hold the same multiset of products as
    `dc1`/`dc2`/`dd1`/`dd2` -- only the order of the two factors differs -- so
    `-dc1-dc2-dc3-dc4+dd1+dd2+dd3+dd4` is `2*(-dc1-dc2+dd1+dd2)`.
  * `dd1`/`dd2` are products of index-symmetric sums, so sixteen multiplies (two
    in the non-NAC form) collapse to four.
  * In the NAC digest each channel's groups 5-8 repeat groups 1-4, because
    `0.5*(I(A)J(B)+J(A)I(B))` and `0.5*(I(B)J(A)+J(B)I(A))` are the same pair of
    products.

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
    """The eight-group `dc` block, single-state (non-NAC) form."""
    return (
        a[i, k] * b[j, l] + a[i, l] * b[j, k] + a[j, k] * b[i, l]
        + a[j, l] * b[i, k] + a[l, j] * b[k, i] + a[k, j] * b[l, i]
        + a[l, i] * b[k, j] + a[k, i] * b[l, j]
    )


def _dd(a, b, i, j, k, l):
    """The eight-group `dd` block, single-state (non-NAC) form."""
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


def _dc_nac(a_i, a_j, b_i, b_j, i, j, k, l):
    """The NAC `dc` block: `_dc` symmetrised over I<->J."""
    return 0.5 * (
        _dc_pairs(a_i, b_j, i, j, k, l) + _dc_pairs(a_j, b_i, i, j, k, l)
    )


def _dc_pairs(a, b, i, j, k, l):
    return _dc(a, b, i, j, k, l)


def _db_nac(m_i, m_j, i, j, k, l):
    """The NAC `db` block as originally written: eight symmetrised groups."""
    groups = [
        ((i, k), (l, j)), ((i, l), (k, j)), ((j, k), (l, i)), ((j, l), (k, i)),
        ((l, j), (i, k)), ((k, j), (i, l)), ((l, i), (j, k)), ((k, i), (j, l)),
    ]
    return sum(
        0.5 * (m_i[p] * m_j[q] + m_j[p] * m_i[q]) for p, q in groups
    )


def _db_nac_halved(m_i, m_j, i, j, k, l):
    """What the engine computes instead of `_db_nac`."""
    groups = [
        ((i, k), (l, j)), ((i, l), (k, j)), ((j, k), (l, i)), ((j, l), (k, i)),
    ]
    return sum(m_i[p] * m_j[q] + m_j[p] * m_i[q] for p, q in groups)


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

    def test_nac_dc3_repeats_dc1(self):
        """Same identity survives the I<->J symmetrisation of the NAC digest."""
        bco1_i, bco1_j, bo2v_i, bo2v_j = _matrices(4)
        lhs, rhs, scale = [], [], 0.0
        for idx in _indices():
            a = _dc_nac(bco1_i, bco1_j, bo2v_i, bo2v_j, *idx)
            lhs.append(a)
            rhs.append(_dc_nac(bo2v_i, bo2v_j, bco1_i, bco1_j, *idx))
            scale = max(scale, abs(a))
        self._assert_close(lhs, rhs, scale, "NAC dc3 vs dc1")

    def test_nac_db_groups_5_to_8_repeat_1_to_4(self):
        co12_i, co12_j = _matrices(2)
        lhs, rhs, scale = [], [], 0.0
        for idx in _indices():
            a = _db_nac(co12_i, co12_j, *idx)
            lhs.append(a)
            rhs.append(_db_nac_halved(co12_i, co12_j, *idx))
            scale = max(scale, abs(a))
        self._assert_close(lhs, rhs, scale, "NAC db halving")


if __name__ == "__main__":
    unittest.main()
