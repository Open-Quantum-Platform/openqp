"""Final-state completeness of FC/HT vibronic spectra measured on HT strength."""

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_vibronic_module():
    path = ROOT / "pyoqp" / "oqp" / "library" / "vibronic.py"
    spec = importlib.util.spec_from_file_location("openqp_vibronic_ht_completeness", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class HerzbergTellerCompletenessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()
        # Identical surfaces: all Franck-Condon weight lies in n -> n.
        cls.model = cls.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [1000.0],
            [[1.0]],
            [0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="shared normal-mode phase",
        )
        cls.omega = 1000.0 * cls.vibronic.CM1_TO_HARTREE

    def spectrum(self, **options):
        transition = self.vibronic.ElectronicTransitionMoment.create(
            (1.0, 0.0, 0.0),
            electronic_state="S1",
            electronic_phase_convention="shared electronic-state phase",
        )
        derivative = self.vibronic.ExcitedStatePropertyDerivative.create(
            np.array([2.0, 0.0, 0.0]).reshape(3, 1),
            property_kind="transition_dipole",
            coordinate_basis="ground_normal",
            coordinate_unit="sqrt(me)*bohr",
            property_unit="e*bohr",
            electronic_state="S1",
            electronic_phase_convention="shared electronic-state phase",
            coordinate_phase_convention="shared normal-mode phase",
            provenance="synthetic analytic derivative",
            electronic_state_role="target_excited_state",
        )
        return self.vibronic.harmonic_vibronic_spectrum(
            self.model,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            transition=transition,
            transition_dipole_derivative=derivative,
            normalization="none",
            **options,
        )

    def test_fc_complete_truncation_that_drops_the_ht_fundamental_is_rejected(self):
        # max_final_quanta=0 keeps all FC weight, but dmu/dQ=2 puts
        # 4/(2 omega) ~ 439 of ~440 strength units into the omitted 0 -> 1 line.
        with self.assertRaisesRegex(ValueError, "transition-strength completeness"):
            self.spectrum(max_final_quanta=0)
        partial = self.spectrum(
            max_final_quanta=0, minimum_franck_condon_completeness=0.0
        )
        self.assertAlmostEqual(partial.franck_condon_completeness, 1.0, places=13)
        self.assertAlmostEqual(
            partial.transition_strength_completeness,
            1.0 / (1.0 + 4.0 / (2.0 * self.omega)),
            places=12,
        )

    def test_thermal_closure_carries_the_two_n_plus_one_factor(self):
        # A linear dipole couples n only to n and n +- 1, so finals up to
        # max_initial_quanta + 1 are complete and the ratio must be exactly one.
        # At 800 K the n >= 1 states hold ~16% of the population, so a closure
        # written for n = 0 alone would miss this by far more than rounding.
        result = self.spectrum(
            max_final_quanta=5, temperature_kelvin=800.0, max_initial_quanta=4
        )
        self.assertGreater(result.retained_thermal_population, 0.999)
        self.assertAlmostEqual(result.transition_strength_completeness, 1.0, places=11)
        self.assertAlmostEqual(result.franck_condon_completeness, 1.0, places=11)


if __name__ == "__main__":
    unittest.main()
