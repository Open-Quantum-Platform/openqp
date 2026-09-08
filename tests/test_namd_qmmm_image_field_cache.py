"""NAMD_QMMM periodic image field: the cached active-state gradient is only
reused for the same state at the same geometry, and after a surface hop the
force assembly re-iterates the field for the new state (seeded from the
previous state's charges) instead of evaluating it in the old state's field.

Pure bookkeeping test: the electronic-structure calls are stubbed out, so it
needs only the importable Python layer."""
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

try:
    from oqp.library import namd as namd_mod
    _HAVE = True
except Exception:  # pragma: no cover - uncompiled backend
    _HAVE = False


def _bare(active, r_all, mm_units=None):
    nd = object.__new__(namd_mod.NAMD_QMMM)
    nd.mol = SimpleNamespace(data={"OQP::partial_charges": np.zeros(3)},
                             energies={active: -1.0, 1: -1.0, 2: -0.9},
                             config={"input": {"method": "tdhf"}, "properties": {}})
    nd.driver = SimpleNamespace(espf_full=True)
    nd._u = None
    nd.active = active
    nd.r_all = np.array(r_all, dtype=float)
    nd._img_ctx = None
    nd._grad_cache = None
    nd._q_img = np.array([0.1, -0.2, 0.1])
    return nd


@unittest.skipUnless(_HAVE, "compiled OpenQP backend unavailable")
class TestImageFieldCache(unittest.TestCase):
    def test_cache_hit_needs_same_state_and_geometry(self):
        r = np.arange(9.0).reshape(3, 3)
        g = np.ones((3, 3)); q = np.array([0.3, -0.5, 0.2])
        nd = _bare(1, r)
        nd._grad_cache = (1, g, q, r.copy())
        out = nd._qm_gradient()
        np.testing.assert_array_equal(out, g)
        np.testing.assert_array_equal(nd.mol.data["OQP::partial_charges"], q)
        self.assertIsNone(nd._grad_cache)          # consumed once
        # a different state or a moved geometry must NOT return the cached gradient
        for active, geom in ((2, r), (1, r + 1e-9)):
            nd = _bare(active, geom)
            nd._grad_cache = (1, g, q, r.copy())
            with mock.patch.object(namd_mod, "Gradient") as grad_cls, \
                    mock.patch.object(namd_mod.oqp, "grad_esp_qmmm_excited", create=True):
                nd.mol.grads = {active: np.zeros(9)}
                nd.mol.data["OQP::ESPF_GRAD"] = np.zeros(9)
                nd._qm_gradient()
                grad_cls.assert_called_once()

    def test_total_force_reiterates_field_after_hop_only(self):
        r = np.arange(9.0).reshape(3, 3)
        nd = _bare(1, r)
        nd._img_ctx = {"geom": r.copy()}
        nd._grad_cache = (1, np.zeros((3, 3)), np.zeros(3), r.copy())
        calls = []
        def refine(q0):                             # stub: records the seed, caches like the real one
            calls.append(np.array(q0))
            nd._grad_cache = (nd.active, np.zeros((3, 3)), np.zeros(3), r.copy())
            return np.zeros(3)
        nd._refine_image_field = refine
        nd._total_force_espf = lambda potmm, gqm, pchg: (gqm, 0.0)
        nd._apply_odp_to_force_energy = lambda f, e: (f, e)
        with mock.patch.object(namd_mod, "dump_log"):
            nd._total_force(np.zeros(3))
            self.assertEqual(calls, [])             # same state: cached gradient reused
            nd.active = 2                           # surface hop 1 -> 2 at the same geometry
            nd._total_force(np.zeros(3))
            self.assertEqual(len(calls), 1)
            np.testing.assert_array_equal(calls[0], nd._q_img)   # seeded from the previous field
            nd.r_all = r + 1.0                      # geometry moved without an electronic step
            with self.assertRaises(RuntimeError):
                nd._total_force(np.zeros(3))
        nd_np = _bare(1, r)                         # non-periodic: no context, plain path
        nd_np._grad_cache = (1, np.zeros((3, 3)), np.zeros(3), r.copy())
        nd_np._total_force_espf = lambda potmm, gqm, pchg: (gqm, 0.0)
        nd_np._apply_odp_to_force_energy = lambda f, e: (f, e)
        nd_np._total_force(np.zeros(3))


if __name__ == "__main__":
    unittest.main()
