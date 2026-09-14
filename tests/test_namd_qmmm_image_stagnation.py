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
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAMD = ROOT / "pyoqp" / "oqp" / "library" / "namd.py"
DRIVER = ROOT / "pyoqp" / "oqp" / "library" / "qmmm_driver.py"

DRV = types.SimpleNamespace(IMAGE_STAGNANT_MINITER=10, IMAGE_TOL_STAGNANT=1e-5, IMAGE_ETOL=1e-8)
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
        self.assertIn("    IMAGE_TOL_STAGNANT = 1e-5\n", src)
        self.assertIn("    IMAGE_ETOL = 1e-8\n", src)

    def test_both_reference_loops_fall_back_only_after_the_strict_test(self):
        src = NAMD.read_text()
        pattern = (r"e_hist\.append\(float\(mol\.get_scf_energy\(\)\)\)\s*\n"
                   r"\s*if delta < self\.driver\.IMAGE_TOL:\s*\n(?:.*\n){3}"
                   r"\s*if _image_field_stagnant\(it, delta, e_hist, self\.driver\):[\s\S]{0,500}?"
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
        self.assertFalse(self.stagnant(20, 3.5e-5, NOISE, DRV))

    def test_moving_energy_is_rejected(self):
        # a switch between two SCF solutions moves the energy far beyond 1e-8
        self.assertFalse(self.stagnant(20, 2.7e-6, [-456.3375043006, -456.3375750506, -456.3375043003], DRV))

    def test_needs_three_energies(self):
        self.assertFalse(self.stagnant(20, 2.7e-6, NOISE[:2], DRV))

    def test_only_the_last_three_energies_count(self):
        self.assertTrue(self.stagnant(12, 2.7e-6, [-456.30] + NOISE, DRV))


if __name__ == "__main__":
    unittest.main()
