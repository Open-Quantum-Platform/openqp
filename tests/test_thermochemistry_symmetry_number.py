"""Rigid-rotor thermochemistry: symmetry number, linear rotors, Gibbs sign.

Three defects are pinned here, all of which changed printed numbers:

1. The rotational symmetry number sigma was absent from the partition function,
   so S_rot -- and every S and G derived from it -- was too large by R*ln(sigma).
   Water: 11.8271 vs the correct 10.4497 cal/mol/K.

2. Linear molecules printed -inf. A linear rotor has one vanishing principal
   moment, and both the rotational entropy (via rt[0], which numpy sorts
   ascending and is therefore that infinity) and Grimme's free-rotor average
   (via mean(rc), likewise infinite) propagated it. Fixing only the first still
   leaves the TOTAL entropy at -inf, so both are covered.

3. G was assembled as H + TS against a printed formula reading G = H - TS.

The totals are checked against tabulated standard entropies, which is the only
external gate available for this path -- no shipped reference pins any
thermochemistry number.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

ANGSTROM_TO_BOHR = 1.8897259886
# Ha -> cal/mol/K at 298.15 K
HARTREE_TS_TO_CAL = 627.5094740631 / 298.15 * 1000.0
GAS_CONSTANT_CAL = 1.98720425


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
            "frequency_thermo_under_test",
            ROOT / "pyoqp/oqp/library/frequency.py",
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules["frequency_thermo_under_test"] = module
        spec.loader.exec_module(module)
        return module
    finally:
        for name, module in saved_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


def load_symmetry_detect_module():
    spec = importlib.util.spec_from_file_location(
        "symmetry_detect_thermo_under_test",
        ROOT / "pyoqp/oqp/library/symmetry_detect.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["symmetry_detect_thermo_under_test"] = module
    spec.loader.exec_module(module)
    return module


def thermo(frequency, atoms, mass, geometry_angstrom, freqs, sigma, linear):
    """Run thermal_analysis on a geometry given in Angstrom."""
    coords = np.asarray(geometry_angstrom, dtype=float) * ANGSTROM_TO_BOHR
    mass = np.asarray(mass, dtype=float)
    ncoord = 3 * len(atoms)
    _, _, inertia = frequency.normal_mode(
        coords.ravel(), mass, np.zeros((ncoord, ncoord)))
    return frequency.thermal_analysis(
        energy=0.0, atoms=atoms, mass=mass,
        freqs=np.asarray(freqs, dtype=float), inertia=inertia,
        temperature=298.15, linear=linear, sigma=sigma, mult=1)


def total_entropy_cal(data):
    return (data['st_el'] + data['st_trans']
            + data['st_rot'] + data['st_vib']) * HARTREE_TS_TO_CAL


WATER = ([8, 1, 1], [15.9949, 1.00783, 1.00783],
         [[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]],
         [1595.0, 3657.0, 3756.0])

CARBON_DIOXIDE = ([8, 6, 8], [15.9949, 12.0, 15.9949],
                  [[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]],
                  [667.0, 667.0, 1333.0, 2349.0])

CARBON_MONOXIDE = ([6, 8], [12.0, 15.9949],
                   [[0.0, 0.0, 0.0], [0.0, 0.0, 1.128]],
                   [2143.0])


class TestRotationalSymmetryNumber(unittest.TestCase):
    """sigma is the order of the rotational subgroup: count det = +1 operations."""

    def test_textbook_symmetry_numbers(self):
        detect = load_symmetry_detect_module()

        def ring(radius, count):
            return [[radius * np.cos(i * np.pi / 3),
                     radius * np.sin(i * np.pi / 3), 0.0] for i in range(count)]

        cases = [
            ('water C2v', [8, 1, 1], WATER[2], 2),
            ('methane Td', [6, 1, 1, 1, 1],
             [[0, 0, 0], [0.6276, 0.6276, 0.6276], [-0.6276, -0.6276, 0.6276],
              [-0.6276, 0.6276, -0.6276], [0.6276, -0.6276, -0.6276]], 12),
            ('benzene D6h', [6] * 6 + [1] * 6,
             ring(1.39, 6) + ring(2.47, 6), 12),
            # Linear molecules need no special case: the seed set collapses to
            # the two permutation classes a linear geometry admits.
            ('carbon dioxide D-inf-h', [8, 6, 8], CARBON_DIOXIDE[2], 2),
            ('carbon monoxide C-inf-v', [6, 8], CARBON_MONOXIDE[2], 1),
        ]
        for name, charges, geometry, expected in cases:
            with self.subTest(molecule=name):
                coords = np.asarray(geometry, dtype=float) * ANGSTROM_TO_BOHR
                self.assertEqual(
                    detect.rotational_symmetry_number(charges, coords), expected)

    def test_high_order_axes_are_recovered_by_group_closure(self):
        """The element survey stops at order 8, but closure rebuilds the rest.

        Two perpendicular C2 axes -- or two mirror planes -- separated by pi/n
        multiply to C_n, so a ring with any second symmetry element regenerates
        its full C_n even though the survey recorded only a divisor.
        """
        detect = load_symmetry_detect_module()

        def ring(n, radius=1.4):
            return [[radius * np.cos(2 * np.pi * i / n),
                     radius * np.sin(2 * np.pi * i / n), 0.0] for i in range(n)]

        for n, expected in ((9, 18), (10, 20), (11, 22), (12, 24)):
            with self.subTest(ring=f'D{n}h'):
                coords = np.asarray(ring(n), dtype=float) * ANGSTROM_TO_BOHR
                self.assertEqual(
                    detect.rotational_symmetry_number([6] * n, coords), expected)

    def test_a_single_atom_is_sigma_one(self):
        detect = load_symmetry_detect_module()
        self.assertEqual(
            detect.rotational_symmetry_number([10], [[0.0, 0.0, 0.0]]), 1)

    def test_failure_falls_back_to_one_rather_than_inventing_symmetry(self):
        detect = load_symmetry_detect_module()
        # Degenerate input the detector cannot handle: must not raise, and must
        # not claim symmetry. sigma=1 reproduces the old (over-counted) entropy,
        # which is why callers print the value they got.
        self.assertEqual(
            detect.rotational_symmetry_number([], np.zeros((0, 3))), 1)


class TestSymmetryNumberEntersTheEntropy(unittest.TestCase):
    def test_water_rotational_entropy_drops_by_r_ln_two(self):
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = WATER

        without = thermo(frequency, atoms, mass, geometry, freqs, 1, False)
        with_sigma = thermo(frequency, atoms, mass, geometry, freqs, 2, False)

        drop = (without['st_rot'] - with_sigma['st_rot']) * HARTREE_TS_TO_CAL
        self.assertAlmostEqual(drop, GAS_CONSTANT_CAL * np.log(2.0), places=3)
        self.assertAlmostEqual(
            with_sigma['st_rot'] * HARTREE_TS_TO_CAL, 10.4497, places=3)

    def test_only_the_rotational_term_moves(self):
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = WATER

        without = thermo(frequency, atoms, mass, geometry, freqs, 1, False)
        with_sigma = thermo(frequency, atoms, mass, geometry, freqs, 2, False)

        for key in ('zpe', 'u_trans', 'u_rot', 'u_vib', 'pv',
                    'st_el', 'st_trans', 'st_vib'):
            self.assertEqual(without[key], with_sigma[key], msg=key)

    def test_reported_sigma_is_what_was_applied(self):
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = WATER
        data = thermo(frequency, atoms, mass, geometry, freqs, 2, False)
        self.assertEqual(data['sigma'], 2)
        self.assertFalse(data['linear'])


class TestLinearRotorsAreFinite(unittest.TestCase):
    """Both -inf sources are covered; fixing only the rotational one is not enough."""

    def test_carbon_dioxide_entropy_terms_are_all_finite(self):
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = CARBON_DIOXIDE
        data = thermo(frequency, atoms, mass, geometry, freqs, 2, True)

        self.assertTrue(np.isfinite(data['st_rot']), 'rotational term is -inf')
        self.assertTrue(np.isfinite(data['st_vib']),
                        'free-rotor average still carries the infinite '
                        'rotational constant')
        self.assertTrue(np.isfinite(total_entropy_cal(data)))

    def test_linear_rotational_energy_is_rt_not_three_halves_rt(self):
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = CARBON_DIOXIDE
        linear = thermo(frequency, atoms, mass, geometry, freqs, 2, True)
        # Deliberately inconsistent: a linear molecule declared nonlinear, so
        # that u_rot can be compared across the two branches. That drives the
        # nonlinear entropy formula onto a vanishing moment and legitimately
        # divides by zero -- contained here rather than allowed to escape as a
        # RuntimeWarning from the entropy code, which is alarming to read in a
        # PR that exists to fix -inf entropies. Only u_rot is asserted on.
        with np.errstate(divide='ignore', invalid='ignore'):
            treated_as_nonlinear = thermo(
                frequency, atoms, mass, geometry, freqs, 2, False)
        self.assertAlmostEqual(
            linear['u_rot'] * 1.5, treated_as_nonlinear['u_rot'], places=12)

    def test_diatomic_is_finite(self):
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = CARBON_MONOXIDE
        data = thermo(frequency, atoms, mass, geometry, freqs, 1, True)
        self.assertTrue(np.isfinite(total_entropy_cal(data)))


class TestLinearRotorsAreOrientationIndependent(unittest.TestCase):
    """The vanishing moment must be selected against on the INERTIA.

    `eigh` returns it as exactly 0.0 only when the molecular axis happens to lie
    along a coordinate axis. Tilt the molecule and it comes back as a tiny value
    of either sign, and a test on the derived rotational constants lets it
    through:

        CO2 along (1,1,1)   moment -5.4e-31  ->  S_vib = NaN
        CO2 tilted in xy    moment +5.3e-15  ->  S_rot = -62.26 cal/mol/K

    Both used to pass an isfinite-and-positive filter. The threshold now matches
    the one normal_mode uses to decide `linear`.
    """

    ORIENTATIONS = {
        'along z': [0.0, 0.0, 1.0],
        'along (1,1,1)': [1.0, 1.0, 1.0],
        'tilted in xy': [0.9319, 0.3626, 0.0],
        'arbitrary axis': [0.31, -0.77, 0.55],
    }

    def test_entropy_does_not_depend_on_the_molecular_axis(self):
        frequency = load_frequency_module()
        _, mass, _, freqs = CARBON_DIOXIDE

        reference = None
        for name, axis in self.ORIENTATIONS.items():
            with self.subTest(orientation=name):
                unit = np.asarray(axis, dtype=float)
                unit = unit / np.linalg.norm(unit)
                geometry = np.array([-1.16 * unit, [0.0, 0.0, 0.0], 1.16 * unit])
                data = thermo(frequency, [8, 6, 8], mass, geometry, freqs, 2, True)

                self.assertTrue(np.isfinite(data['st_rot']), name)
                self.assertTrue(np.isfinite(data['st_vib']), name)
                total = total_entropy_cal(data)
                if reference is None:
                    reference = total
                self.assertAlmostEqual(total, reference, places=9,
                                       msg=f'{name} disagrees with "along z"')

    def test_a_nonlinear_molecule_keeps_all_three_moments(self):
        """The selection must not quietly drop a moment from a real top."""
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = WATER
        data = thermo(frequency, atoms, mass, geometry, freqs, 2, False)
        self.assertAlmostEqual(data['st_rot'] * HARTREE_TS_TO_CAL, 10.4497,
                               places=3)


class TestAgainstTabulatedStandardEntropies(unittest.TestCase):
    """The external gate. Internal identities cannot catch a missing sigma."""

    def test_standard_entropies(self):
        frequency = load_frequency_module()
        cases = [
            ('water', WATER, 2, False, 45.10),
            ('carbon dioxide', CARBON_DIOXIDE, 2, True, 51.07),
            ('carbon monoxide', CARBON_MONOXIDE, 1, True, 47.21),
        ]
        for name, (atoms, mass, geometry, freqs), sigma, linear, reference in cases:
            with self.subTest(molecule=name):
                data = thermo(frequency, atoms, mass, geometry,
                              freqs, sigma, linear)
                self.assertAlmostEqual(
                    total_entropy_cal(data), reference, delta=0.25,
                    msg=f'{name}: S deviates from the tabulated value')

    def test_water_without_sigma_misses_the_tabulated_value(self):
        """Guard against the fix silently regressing to sigma = 1."""
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = WATER
        without = thermo(frequency, atoms, mass, geometry, freqs, 1, False)
        self.assertGreater(total_entropy_cal(without) - 45.10, 1.0)


class TestGibbsSign(unittest.TestCase):
    def test_gibbs_correction_subtracts_the_entropy(self):
        import re
        source = (ROOT / 'pyoqp/oqp/utils/file_utils.py').read_text()
        self.assertRegex(source, r'g_el\s*=\s*h_el\s*-\s*st\b')
        self.assertNotRegex(source, r'g_el\s*=\s*h_el\s*\+\s*st\b')
        # The printed formula and the arithmetic must agree.
        self.assertIn('G = H - TS', source)


if __name__ == '__main__':
    unittest.main()
