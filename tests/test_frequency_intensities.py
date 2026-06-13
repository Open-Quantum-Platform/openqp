import importlib.util
import sys
import types
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_frequency_module():
    stub_names = ("oqp", "oqp.utils", "oqp.utils.constants")
    saved_modules = {name: sys.modules.get(name) for name in stub_names}

    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp

        utils = types.ModuleType("oqp.utils")
        utils.__path__ = []
        sys.modules["oqp.utils"] = utils

        constants = types.ModuleType("oqp.utils.constants")
        constants.SPEED_OF_LIGHT = 2.99792458e10
        constants.ATMOS = 101.325
        constants.BOHR = 0.52917721090299996e-10
        constants.FREQ_TO_INV_CM = 5140.489195376594
        constants.AMU_to_KG = 1.66053886E-27
        constants.J_TO_AU = 1 / (4.184 * 627.509541 * 1000.0)
        constants.GAS_CONSTANT = 8.3144621
        constants.PLANCK_CONSTANT = 6.62606957e-34
        constants.BOLTZMANN_CONSTANT = 1.3806488e-23
        constants.AVOGADRO_CONSTANT = 6.0221415e23
        sys.modules["oqp.utils.constants"] = constants

        spec = importlib.util.spec_from_file_location(
            "frequency_under_test",
            ROOT / "pyoqp/oqp/library/frequency.py",
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules["frequency_under_test"] = module
        spec.loader.exec_module(module)
        return module
    finally:
        for name, module in saved_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


class TestVibrationalIntensities(unittest.TestCase):
    def test_infrared_intensities_project_cartesian_dipole_derivatives_onto_modes(self):
        frequency = load_frequency_module()
        modes = np.array([[1.0, 0.0], [0.0, 2.0]])
        dipole_derivatives = np.array(
            [
                [1.0, 3.0],
                [2.0, 0.0],
                [0.0, 4.0],
            ]
        )

        intensities, mode_dipoles = frequency.infrared_intensities(dipole_derivatives, modes)

        np.testing.assert_allclose(mode_dipoles, [[1.0, 2.0, 0.0], [6.0, 0.0, 8.0]])
        np.testing.assert_allclose(
            intensities,
            frequency.IR_INTENSITY_CONVERSION_KM_MOL * np.array([5.0, 100.0]),
        )

    def test_raman_activities_use_isotropic_and_anisotropic_polarizability_derivatives(self):
        frequency = load_frequency_module()
        modes = np.array([[1.0, 0.0], [0.0, 1.0]])
        polarizability_derivatives = np.zeros((3, 3, 2))
        polarizability_derivatives[:, :, 0] = np.eye(3) * 2.0
        polarizability_derivatives[:, :, 1] = np.diag([1.0, -1.0, 0.0])

        activities, mode_polarizabilities = frequency.raman_activities(polarizability_derivatives, modes)

        np.testing.assert_allclose(mode_polarizabilities[0], np.eye(3) * 2.0)
        np.testing.assert_allclose(mode_polarizabilities[1], np.diag([1.0, -1.0, 0.0]))
        np.testing.assert_allclose(activities, [180.0, 21.0])
    def test_linear_diatomic_retains_one_vibrational_mode(self):
        frequency = load_frequency_module()
        coord = np.array([0.0, 0.0, -0.7, 0.0, 0.0, 0.7])
        mass = np.array([1.0, 1.0])
        hessian = np.zeros((6, 6))
        hessian[2, 2] = 1.0
        hessian[5, 5] = 1.0
        hessian[2, 5] = -1.0
        hessian[5, 2] = -1.0

        freqs, modes, inertia = frequency.normal_mode(coord, mass, hessian)

        self.assertEqual(freqs.shape, (1,))
        self.assertEqual(modes.shape, (1, 6))
        self.assertGreater(freqs[0], 0.0)
        self.assertEqual(inertia.shape, (3,))


if __name__ == "__main__":
    unittest.main()
