"""Broadening-grid size and tracking-overlap roundoff in the vibronic tools."""

import importlib.util
from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_vibronic_module():
    path = ROOT / "pyoqp" / "oqp" / "library" / "vibronic.py"
    spec = importlib.util.spec_from_file_location("openqp_vibronic_grid_overlap", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class GridAndOverlapBounds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()
        cls.model = cls.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [1000.0],
            [[1.0]],
            [0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic grid test",
        )

    def spectrum(self, **options):
        return self.vibronic.harmonic_vibronic_spectrum(
            self.model,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            max_final_quanta=8,
            fwhm_cm1=0.001,
            normalization="none",
            **options,
        )

    def test_automatic_grid_is_bounded_before_allocation(self):
        # An 8000 cm-1 line span at a 0.001 cm-1 FWHM would need ~2e8 points.
        with mock.patch.object(
            self.vibronic.np, "linspace", side_effect=AssertionError("grid was allocated")
        ):
            with self.assertRaisesRegex(ValueError, "max_grid_points"):
                self.spectrum()
        explicit = np.linspace(19990.0, 20010.0, 2001)
        self.assertEqual(self.spectrum(grid_cm1=explicit).wavenumbers_cm1.size, 2001)

    def test_tracking_overlap_roundoff_above_one_is_accepted(self):
        def derivative(overlap):
            return self.vibronic.finite_difference_excited_state_property(
                [[1.0, 0.0, 0.0]],
                [[0.8, 0.0, 0.0]],
                [0.01],
                property_kind="state_dipole",
                coordinate_basis="excited_normal",
                coordinate_unit="sqrt(amu)*bohr",
                property_unit="e*bohr",
                electronic_state="S1",
                electronic_phase_convention="central-state gauge",
                coordinate_phase_convention="synthetic modes",
                state_tracking_overlaps=[[overlap, 0.999]],
                minimum_overlap=0.99,
            )

        # MRSFTrackedPropertySnapshot accepts matched_overlap up to 1 + 1e-8.
        result = derivative(1.0000000000000002)
        np.testing.assert_allclose(np.asarray(result.values).reshape(-1), [10.0, 0.0, 0.0])
        with self.assertRaisesRegex(ValueError, "state_tracking_overlaps"):
            derivative(1.001)


if __name__ == "__main__":
    unittest.main()
