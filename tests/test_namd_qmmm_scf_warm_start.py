"""NAMD QM/MM embedded SCF: each MD step starts from the previous step's
orbitals, and solves that start from converged orbitals use the second-order
converger first.

Background (DNA thymine MRSF NAMD, ROHF triplet reference): a fresh Hueckel
guess every step started the SCF ~1.3 Hartree above the solution, and a DIIS
refill of already-converged orbitals could swap an occupied and an open
orbital (+0.2 Hartree), stall, and fall down the escalation ladder into TRAH.
"""
import unittest
from unittest import mock

import numpy as np


def _runtime_available():
    try:
        import oqp.library.namd  # noqa: F401
        return True
    except Exception:
        return False


def _pack_lower(m):
    n = m.shape[0]
    return np.array([m[i, j] for i in range(n) for j in range(i + 1)], dtype=float)


class _FakeData(dict):
    def __init__(self, nbf):
        super().__init__()
        self.nbf = nbf

    def get_basis(self):
        return {"nbf": self.nbf}


class _FakeMol:
    def __init__(self, data):
        self.data = data


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestWarmSolvesStartWithSoscf(unittest.TestCase):

    def _sp(self, result=True, exc=None):
        seen = []

        class SP:
            converger_type = "diis"

            def _run_scf(self):
                seen.append(self.converger_type)
                if exc is not None:
                    raise exc
                return result

        return SP(), seen

    def test_warm_solve_uses_soscf_and_restores_the_primary(self):
        from oqp.library.namd import NAMD_QMMM
        sp, seen = self._sp()
        NAMD_QMMM._embedded_scf(sp, warm=True)
        NAMD_QMMM._embedded_scf(sp)
        self.assertEqual(seen, ["soscf", "diis"])
        self.assertEqual(sp.converger_type, "diis")

    def test_primary_is_restored_when_the_solve_raises(self):
        from oqp.library.namd import NAMD_QMMM
        sp, _ = self._sp(exc=RuntimeError("solver failure"))
        with self.assertRaises(RuntimeError):
            NAMD_QMMM._embedded_scf(sp, warm=True)
        self.assertEqual(sp.converger_type, "diis")

    def test_unconverged_warm_solve_still_stops_the_run(self):
        from oqp.library.namd import NAMD_QMMM
        sp, _ = self._sp(result=False)
        with self.assertRaisesRegex(RuntimeError, "did not converge"):
            NAMD_QMMM._embedded_scf(sp, warm=True)
        self.assertEqual(sp.converger_type, "diis")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestStepStartsFromPreviousOrbitals(unittest.TestCase):
    """_start_scf_orbitals with the integral and guess back ends mocked."""

    NBF = 6

    def _owner(self, ready, basis_nbf=None):
        from oqp.library.namd import NAMD_QMMM
        n = self.NBF
        rng = np.random.default_rng(3)
        c0, _ = np.linalg.qr(rng.normal(size=(n, n)))          # orthonormal for S = 1
        a = rng.normal(size=(n, n)) * 0.03
        s1 = np.eye(n) + 0.5 * (a + a.T)                        # overlap at the new geometry
        data = _FakeData(basis_nbf or n)
        data["OQP::VEC_MO_A"] = np.ascontiguousarray(c0.T).copy()
        data["OQP::VEC_MO_B"] = np.zeros((n, n))
        data["OQP::E_MO_A"] = np.linspace(-1.0, 1.0, n)
        data["OQP::E_MO_B"] = np.zeros(n)
        data["OQP::SM"] = _pack_lower(s1)
        owner = NAMD_QMMM.__new__(NAMD_QMMM)
        owner.mol = _FakeMol(data)
        if ready:
            owner._scf_orbitals_ready = True
        return owner, c0, s1

    def _patches(self):
        return (mock.patch("oqp.library.set_basis"), mock.patch("oqp.library.namd.ints_1e"),
                mock.patch("oqp.guess_json"), mock.patch("oqp.library.guess"),
                mock.patch("oqp.library.namd.dump_log"))

    def test_first_step_uses_the_configured_guess(self):
        owner, _, _ = self._owner(ready=False)
        sp = mock.Mock()
        p = self._patches()
        with p[0], p[1], p[2] as gj, p[3], p[4]:
            self.assertFalse(owner._start_scf_orbitals(sp))
            gj.assert_not_called()
        sp._prep_guess.assert_called_once_with()

    def test_later_steps_reuse_orthonormalised_previous_orbitals(self):
        owner, c0, s1 = self._owner(ready=True)
        sp = mock.Mock()
        p = self._patches()
        with p[0] as sb, p[1] as i1, p[2] as gj, p[3] as gs, p[4]:
            self.assertTrue(owner._start_scf_orbitals(sp))
            sb.assert_called_once(); i1.assert_called_once(); gj.assert_called_once(); gs.assert_not_called()
        sp._prep_guess.assert_not_called()
        d = owner.mol.data
        c = np.asarray(d["OQP::VEC_MO_A"]).reshape(self.NBF, self.NBF).T
        np.testing.assert_allclose(c.T @ s1 @ c, np.eye(self.NBF), atol=1e-12)      # orthonormal in the new metric
        self.assertGreater(np.abs(np.diag(c0.T @ s1 @ c)).min(), 0.95)               # each orbital keeps its character
        np.testing.assert_array_equal(d["OQP::VEC_MO_B"], d["OQP::VEC_MO_A"])
        np.testing.assert_array_equal(d["OQP::E_MO_B"], d["OQP::E_MO_A"])
        # control: the stored orbitals were not orthonormal in the new metric
        self.assertGreater(np.abs(c0.T @ s1 @ c0 - np.eye(self.NBF)).max(), 1e-3)

    def test_orbital_count_mismatch_falls_back_to_a_fresh_guess(self):
        owner, _, _ = self._owner(ready=True, basis_nbf=self.NBF + 1)
        sp = mock.Mock()
        p = self._patches()
        with p[0], p[1], p[2] as gj, p[3] as gs, p[4]:
            self.assertFalse(owner._start_scf_orbitals(sp))
            gs.assert_called_once(); gj.assert_not_called()


if __name__ == "__main__":
    unittest.main()
