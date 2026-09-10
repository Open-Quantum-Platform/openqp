"""Thermal and state-cap inputs of the vibronic spectrum must be finite."""

import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load_vibronic_module():
    path = ROOT / "pyoqp" / "oqp" / "library" / "vibronic.py"
    spec = importlib.util.spec_from_file_location("openqp_vibronic_nonfinite", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class NonFiniteVibronicInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()
        cls.model = cls.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [900.0],
            [[1.0]],
            [0.3],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic one-mode model",
        )

    def test_non_finite_temperature_is_rejected(self):
        # NaN used to pass "< 0" and produce NaN populations and intensities.
        for temperature in (float("nan"), float("inf")):
            with self.subTest(temperature=temperature):
                with self.assertRaisesRegex(
                    ValueError, "temperature_kelvin must be non-negative and finite"
                ):
                    self.vibronic.harmonic_vibronic_spectrum(
                        self.model,
                        electronic_origin_cm1=20000.0,
                        origin_kind="zero_zero",
                        max_final_quanta=2,
                        temperature_kelvin=temperature,
                        max_initial_quanta=2,
                        minimum_thermal_population=0.0,
                        minimum_franck_condon_completeness=0.0,
                    )

    def test_non_finite_state_cap_is_rejected(self):
        # A NaN or +inf cap was never exceeded, so the basis size was unbounded.
        for cap in (float("nan"), float("inf")):
            with self.subTest(max_states=cap):
                with self.assertRaisesRegex(
                    ValueError, "state-enumeration bounds must be non-negative and finite"
                ):
                    self.vibronic.enumerate_vibrational_states(2, 3, max_states=cap)


if __name__ == "__main__":
    unittest.main()
