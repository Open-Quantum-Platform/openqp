"""The Duschinsky orthogonality gate must not be disabled by a non-finite tolerance."""

import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load_vibronic_module():
    path = ROOT / "pyoqp" / "oqp" / "library" / "vibronic.py"
    spec = importlib.util.spec_from_file_location("openqp_vibronic_orthogonality", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class OrthogonalityToleranceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()

    def create(self, rotation, tolerance):
        return self.vibronic.HarmonicVibronicModel.create(
            [1000.0, 1200.0],
            [1000.0, 1200.0],
            rotation,
            [0.0, 0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic two-mode model",
            orthogonality_tolerance=tolerance,
        )

    def test_non_finite_tolerances_are_rejected(self):
        # A shear is far from orthogonal; NaN used to wave it through because
        # both "tolerance < 0" and "residual > tolerance" are false for NaN.
        shear = [[1.0, 0.5], [0.0, 1.0]]
        for tolerance in (float("nan"), float("inf")):
            with self.subTest(tolerance=tolerance):
                with self.assertRaisesRegex(ValueError, "orthogonality_tolerance"):
                    self.create(shear, tolerance)

    def test_finite_tolerances_keep_their_meaning(self):
        with self.assertRaisesRegex(ValueError, "orthogonality residual"):
            self.create([[1.0, 0.5], [0.0, 1.0]], 1.0e-8)
        model = self.create([[0.0, 1.0], [1.0, 0.0]], 1.0e-8)
        self.assertLessEqual(model.orthogonality_residual, 1.0e-8)


if __name__ == "__main__":
    unittest.main()
