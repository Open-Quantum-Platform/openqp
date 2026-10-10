"""Periodic ESPF QM/MM NAMD: the reference-density QM-image charge loop accepts
a field whose ESPF charges only fluctuate at the SCF noise floor.

Background (DNA thymine MRSF NAMD, ROHF triplet reference): near a nearly
degenerate open/closed orbital pair, SOSCF/TRAH stop at |g| ~ 1e-8 and the
charges of successive image iterations differ by 1e-6 to 1e-5 e while the SCF
energy no longer moves.  Nine trajectories stopped with "QM-image charge
self-consistency did not converge in 50 iterations" although nothing changed
any more; the same checkpoints replayed on another machine passed those steps.
"""
import re
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NAMD = ROOT / "pyoqp" / "oqp" / "library" / "namd.py"
DRIVER = ROOT / "pyoqp" / "oqp" / "library" / "qmmm_driver.py"

DRV = types.SimpleNamespace(IMAGE_STAGNANT_MINITER=10, IMAGE_TOL_STAGNANT=1e-4, IMAGE_ETOL=1e-7)
# reference SCF energies of image iterations 4-6 at the e25 step-44 noise floor
NOISE = [-456.3375043006, -456.3375043002, -456.3375043003]


def _helper():
    try:
        from oqp.library.namd import _image_field_stagnant
        return _image_field_stagnant
    except Exception:
        return None


class TestImageFieldStagnationSource(unittest.TestCase):
    def test_driver_keeps_the_strict_tolerance_and_defines_the_fallback(self):
        src = DRIVER.read_text()
        self.assertIn("    IMAGE_TOL = 1e-7\n", src)
        self.assertIn("    IMAGE_STAGNANT_MINITER = 10\n", src)
        self.assertIn("    IMAGE_TOL_STAGNANT = 1e-4\n", src)
        self.assertIn("    IMAGE_ETOL = 1e-7\n", src)

    def test_both_reference_loops_fall_back_only_after_the_strict_test(self):
        src = NAMD.read_text()
        pattern = (r"e_hist\.append\(float\(mol\.get_scf_energy\(\)\)\)\s*\n"
                   r"\s*if delta < self\.driver\.IMAGE_TOL:\s*\n(?:.*\n){3}"
                   r"\s*if _image_field_stagnant\(it, delta, e_hist, self\.driver\):[\s\S]{0,900}?"
                   r"converged = True\s*\n\s*break\s*\n"
                   r"\s*q_prev = 0\.5 \* \(q_new \+ q_prev\) if it > 6 else q_new")
        self.assertEqual(len(re.findall(pattern, src)), 2)


@unittest.skipUnless(_helper() is not None, "compiled OpenQP runtime unavailable")
class TestImageFieldStagnant(unittest.TestCase):
    def setUp(self):
        self.stagnant = _helper()

    def test_noise_floor_after_ten_iterations_is_accepted(self):
        self.assertTrue(self.stagnant(9, 2.7e-6, NOISE, DRV))

    def test_not_before_the_minimum_iteration_count(self):
        # a step that converges normally keeps the strict tolerance
        self.assertFalse(self.stagnant(8, 2.7e-6, NOISE, DRV))

    def test_charges_still_moving_are_rejected(self):
        self.assertFalse(self.stagnant(20, 3.5e-4, NOISE, DRV))

    def test_moving_energy_is_rejected(self):
        # a switch between two SCF solutions moves the energy far beyond 1e-8
        self.assertFalse(self.stagnant(20, 2.7e-6, [-456.3375043006, -456.3375750506, -456.3375043003], DRV))

    def test_needs_three_energies(self):
        self.assertFalse(self.stagnant(20, 2.7e-6, NOISE[:2], DRV))

    def test_only_the_last_three_energies_count(self):
        self.assertTrue(self.stagnant(12, 2.7e-6, [-456.30] + NOISE, DRV))


class TestFallbackSource(unittest.TestCase):
    def test_fallback_keeps_the_scf_input_charges(self):
        src = NAMD.read_text()
        blocks = re.findall(r"if _image_field_stagnant\(it, delta, e_hist, self\.driver\):([\s\S]*?)converged = True", src)
        self.assertEqual(len(blocks), 2)
        for block in blocks:
            self.assertNotIn("q_prev = q_new", block)


class _FakeData(dict):
    def __init__(self, nbf):
        super().__init__()
        self.nbf = nbf

    def get_basis(self):
        return {"nbf": self.nbf}


@unittest.skipUnless(_helper() is not None, "compiled OpenQP runtime unavailable")
class TestReferenceImageLoop(unittest.TestCase):
    """NAMD_QMMM._electronic_qmmm's reference-density image loop with the
    electronic structure stubbed: scripted ESPF charges and SCF energies, and an
    identity image matrix, so the MM potential each SCF sees is its input field."""

    NAT, NBF = 2, 2
    BASE = np.array([0.30, -0.30])
    E0 = -456.3375043003

    def _noise(self, n, amp=3e-6):
        return [self.BASE + amp * np.array([(-1) ** k, (-1) ** (k + 1)]) for k in range(n)]

    def _run(self, charges, energies, maxiter=50, tol_stagnant=1e-4, etol=1e-7):
        import oqp.library.namd as namd
        nat, nbf = self.NAT, self.NBF
        data = _FakeData(nbf)
        data["natom"] = nat
        calls = {"scf": 0, "potmm": [], "logs": []}
        mol = types.SimpleNamespace(data=data, config={"input": {"method": "hf"}},
                                    get_hcore=lambda: np.zeros(nbf * (nbf + 1) // 2),
                                    set_hcore=lambda h: None,
                                    get_scf_energy=lambda: energies[calls["scf"] - 1])
        ewald = types.SimpleNamespace(qm_image_matrix=lambda pos: (np.eye(nat), np.zeros((nat, nat, 3))))
        driver = types.SimpleNamespace(espf_full=True, IMAGE_MAXITER=maxiter, IMAGE_TOL=1e-7,
                                       IMAGE_STAGNANT_MINITER=10, IMAGE_TOL_STAGNANT=tol_stagnant, IMAGE_ETOL=etol,
                                       _ewald=lambda: ewald,
                                       _qm_center_positions_bohr=lambda: np.zeros((nat, 3)))
        obj = namd.NAMD_QMMM.__new__(namd.NAMD_QMMM)
        obj.mol, obj.driver = mol, driver
        obj._embedding_field = lambda: (np.zeros(nat), np.zeros((nat, nat)))
        obj._start_scf_orbitals = lambda sp: False
        obj._geometry_key = lambda: "geometry"

        def scf(sp, warm=False):
            calls["potmm"].append(np.array(data["OQP::POTMM"], dtype=float))
            calls["scf"] += 1
        obj._embedded_scf = scf
        seed = {}

        def refine(q):
            seed["q"] = np.array(q, dtype=float)
            return "refined potential"
        obj._refine_image_field = refine

        def form_esp_charges(m):
            data["OQP::partial_charges"] = charges[calls["scf"] - 1]
        drv_mod = types.SimpleNamespace(unpack_lower_tri_single=lambda a, n: np.zeros((n, n)),
                                        unpack_lower_tri_multi=lambda a, n, k: np.zeros((k, n, n)),
                                        pack_lower_tri_single=lambda h: h)
        with mock.patch.dict(sys.modules, {"oqp.library.qmmm_driver": drv_mod}), \
                mock.patch.object(namd.oqp, "espf_op_corr", lambda m: data.__setitem__("OQP::ESPF_CORR", None), create=True), \
                mock.patch.object(namd.oqp, "form_esp_charges", form_esp_charges, create=True), \
                mock.patch.object(namd, "ints_1e", lambda m: None), \
                mock.patch.object(namd, "SinglePoint", mock.MagicMock()), \
                mock.patch.object(namd, "LastStep", mock.MagicMock()), \
                mock.patch.object(namd, "dump_log", lambda m, title="", **kw: calls["logs"].append(title)):
            result = obj._electronic_qmmm(with_overlap=False)
        return result, seed.get("q"), calls

    def test_noise_floor_is_accepted_with_the_field_the_scf_saw(self):
        n = 50
        charges = self._noise(n)
        _, seed, calls = self._run(charges, [self.E0 + 2e-10 * (-1) ** k for k in range(n)])
        self.assertEqual(calls["scf"], 10)
        self.assertTrue(any("accepted on stagnation after 10 iterations" in t for t in calls["logs"]))
        # the active-state loop starts from the input field of the last SCF, not its output charges
        np.testing.assert_array_equal(seed, calls["potmm"][-1])
        self.assertGreater(float(np.abs(seed - charges[9]).max()), 1e-6)

    def test_two_state_alternation_at_the_e25_step_45_size_is_accepted(self):
        """e25 step 45: every SCF converges to the same solution, the energy
        alternates within 1.3e-9 Hartree and max |dq| against the damped field
        settles at 1.56e-5 e.  Accepted at 1e-4 e; the earlier 1e-5 e raised."""
        n = 50
        charges = self._noise(n, amp=1.04e-5)
        energies = [self.E0 + 6.5e-10 * (-1) ** k for k in range(n)]
        _, seed, calls = self._run(charges, energies)
        self.assertEqual(calls["scf"], 10)
        self.assertTrue(any("accepted on stagnation after 10 iterations" in t for t in calls["logs"]))
        with self.assertRaisesRegex(RuntimeError, "did not stagnate"):
            self._run(charges, energies, tol_stagnant=1e-5)

    def test_energy_scatter_of_trah_exits_at_the_e52_step_741_size_is_accepted(self):
        """e52 step 741: every image iteration ends in TRAH, the charges meet the
        fallback and the SCF energies scatter by 3.7e-8 Hartree.  Accepted at
        1e-7 Hartree; the earlier 1e-8 Hartree raised."""
        n = 50
        charges = self._noise(n)
        energies = [self.E0 + 1.85e-8 * (-1) ** k for k in range(n)]
        _, _, calls = self._run(charges, energies)
        self.assertEqual(calls["scf"], 10)
        self.assertTrue(any("accepted on stagnation after 10 iterations" in t for t in calls["logs"]))
        with self.assertRaisesRegex(RuntimeError, "did not stagnate"):
            self._run(charges, energies, etol=1e-8)

    def test_moving_energy_still_raises(self):
        n = 50
        with self.assertRaisesRegex(RuntimeError, "did not stagnate"):
            self._run(self._noise(n), [self.E0 + 1e-6 * (-1) ** k for k in range(n)])

    def test_strict_convergence_is_unchanged(self):
        n = 50
        charges = [self.BASE + 0.01] + [self.BASE] * (n - 1)
        _, seed, calls = self._run(charges, [self.E0] * n)
        self.assertEqual(calls["scf"], 3)
        np.testing.assert_array_equal(seed, self.BASE)
        self.assertFalse(any("stagnation" in t for t in calls["logs"]))

    def test_no_fallback_before_ten_iterations(self):
        n = 9
        with self.assertRaisesRegex(RuntimeError, "did not converge in 9 iterations"):
            self._run(self._noise(n), [self.E0] * n, maxiter=9)


if __name__ == "__main__":
    unittest.main()
