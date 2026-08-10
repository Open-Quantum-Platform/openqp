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


class TestTheCheapScreenNeverChangesSigma(unittest.TestCase):
    """sigma is computed on every Hessian analysis, including `hess.read`.

    The detector's element survey builds a direction per atom pair and dedupes
    each against every direction kept so far, so a geometry with no symmetry --
    where nothing dedupes away -- costs O(N^4): 1.7 s for a random 50-atom
    geometry and 7.1 s for a 75-atom one, both to return sigma = 1. A
    same-element/same-radius screen answers those without the detector.

    It is a necessary condition, so it can only skip work: if no atom has a
    same-element partner at the same radius, every operation fixes every atom,
    and for a nonlinear geometry the only proper one that does so is identity.
    """

    @staticmethod
    def random_geometry(natom):
        rng = np.random.default_rng(natom)
        charges = list(rng.integers(1, 9, size=natom))
        coords = (rng.random((natom, 3)) * 8.0 - 4.0) * ANGSTROM_TO_BOHR
        return charges, coords

    @staticmethod
    def c2_symmetric_geometry(half):
        """`half` random atoms plus their images under C2 about z: sigma = 2."""
        rng = np.random.default_rng(7)
        charges = list(rng.integers(1, 9, size=half))
        xyz = (rng.random((half, 3)) * 8.0 - 4.0) * ANGSTROM_TO_BOHR
        return charges + charges, np.vstack(
            [xyz, xyz * np.array([-1.0, -1.0, 1.0])])

    @staticmethod
    def spy_on_the_detector(detect):
        """Count detector calls without suppressing it.

        Raising from the stub would prove nothing: rotational_symmetry_number
        wraps the call in `except Exception: return 1`, so a stub that raised
        would give the same answer whether or not the screen fired.
        """
        real = detect.enumerate_full_group
        calls = []

        def spy(*args, **kwargs):
            calls.append(1)
            return real(*args, **kwargs)

        detect.enumerate_full_group = spy
        return calls

    def test_a_geometry_with_no_partners_answers_without_the_detector(self):
        detect = load_symmetry_detect_module()
        calls = self.spy_on_the_detector(detect)
        charges, coords = self.random_geometry(50)

        self.assertEqual(detect.rotational_symmetry_number(charges, coords), 1)
        self.assertEqual(calls, [])

    def test_a_symmetric_geometry_of_the_same_size_still_reaches_it(self):
        """The control: without this, a screen that always fired would pass.

        Same atom count as the case above, but every atom now has a partner at
        its own radius, so the screen must decline and the detector must run.
        """
        detect = load_symmetry_detect_module()
        calls = self.spy_on_the_detector(detect)
        charges, coords = self.c2_symmetric_geometry(25)

        self.assertFalse(detect._every_atom_is_its_own_class(
            np.asarray(charges, dtype=float), coords, 1.0e-5))
        self.assertEqual(detect.rotational_symmetry_number(charges, coords), 2)
        self.assertEqual(len(calls), 1)

    def test_nearly_linear_geometries_are_left_to_the_detector(self):
        """Both reported by chatgpt-codex-connector, against successive versions.

        The match is approximate, so 'every atom is fixed' does not force the
        identity for a geometry that is *almost* collinear: a rotation by theta
        moves an atom by 2*sin(theta/2)*d(atom, axis), so an atom close enough
        to the axis is fixed to within tolerance.

        The second case is the one that defeats `_is_linear`: its axis comes
        from the first atom with a non-negligible position, and here that atom
        sits 0.1*tolerance from the center, so the axis points along x and every
        other atom reads as a huge residual -- even though all four are within
        0.3*tolerance of the z axis. Measuring the off-axis extent instead needs
        no axis estimate.
        """
        detect = load_symmetry_detect_module()
        tolerance = 1.0e-5
        charges = np.asarray([1.0, 1.0, 1.0, 1.0])
        cases = {
            # Largest residual 1.0276 * tolerance; a C8 about z moves every
            # atom by at most 8.5e-6 bohr.
            'off-axis by about one tolerance': np.array([
                [0.350000 * tolerance, 0.0, -3.25],
                [0.175000 * tolerance, 0.0, -1.25],
                [-1.108333 * tolerance, 0.0, 0.75],
                [0.583333 * tolerance, 0.0, 3.75],
            ]),
            # A C2 about z moves the farthest atom by 0.6 * tolerance.
            'first atom sits at the center': np.array([
                [0.1 * tolerance, 0.0, 0.0],
                [0.3 * tolerance, 0.0, 1.0],
                [-0.2 * tolerance, 0.0, 2.0],
                [-0.2 * tolerance, 0.0, -3.0],
            ]),
        }
        for name, coords in cases.items():
            with self.subTest(geometry=name):
                # The control: `_is_linear` at its own threshold does NOT catch
                # these, which is why the screen used to fire on them.
                centered = coords - np.einsum(
                    'i,ij->j', charges, coords) / float(np.sum(charges))
                self.assertIsNone(detect._is_linear(centered, tolerance))

                self.assertFalse(detect._every_atom_is_its_own_class(
                    charges, coords, tolerance))

    def test_the_screen_declines_on_every_symmetric_molecule_here(self):
        detect = load_symmetry_detect_module()
        cases = {
            'water': ([8, 1, 1], WATER[2]),
            'carbon dioxide': ([8, 6, 8], CARBON_DIOXIDE[2]),
            'benzene': ([6] * 6 + [1] * 6,
                        [[1.39 * np.cos(i * np.pi / 3),
                          1.39 * np.sin(i * np.pi / 3), 0.0] for i in range(6)]
                        + [[2.47 * np.cos(i * np.pi / 3),
                            2.47 * np.sin(i * np.pi / 3), 0.0]
                           for i in range(6)]),
        }
        for name, (charges, geometry) in cases.items():
            with self.subTest(molecule=name):
                coords = np.asarray(geometry, dtype=float) * ANGSTROM_TO_BOHR
                self.assertFalse(detect._every_atom_is_its_own_class(
                    np.asarray(charges, dtype=float), coords, 1.0e-5))

    def test_carbon_monoxide_is_left_to_the_detector(self):
        """Linear geometries are excluded from the screen by construction.

        A rotation about the molecular axis fixes every atom without being the
        identity, so 'nothing can be permuted' does not imply sigma = 1 there.
        """
        detect = load_symmetry_detect_module()
        charges, geometry, _ = (CARBON_MONOXIDE[0], CARBON_MONOXIDE[2],
                                CARBON_MONOXIDE[3])
        coords = np.asarray(geometry, dtype=float) * ANGSTROM_TO_BOHR
        self.assertFalse(detect._every_atom_is_its_own_class(
            np.asarray(charges, dtype=float), coords, 1.0e-5))


class TestAnUnusableToleranceCannotInventSymmetry(unittest.TestCase):
    """Reported by chatgpt-codex-connector on #320.

    `_match_permutation` rejects a partner when `dist[j] > tolerance`, and every
    comparison against NaN is false, so a NaN tolerance makes the detector
    accept geometrically invalid operations. Measured on 40 random asymmetric
    geometries: 12 returned sigma = 2 instead of 1, i.e. G wrong by
    R*T*ln(2) = 0.41 kcal/mol.
    """

    # One of the 12, kept verbatim so the scenario is reproducible.
    CHARGES = [6, 4, 6]
    COORDS = np.array([[-1.4854231, -1.5451678, 0.61338209],
                       [1.41382842, -1.19288346, -1.12792545],
                       [0.86633854, -0.11720132, -0.33911228]])

    def test_the_detector_really_was_fooled_by_it(self):
        """The control: without this, the guards below could be pinning nothing.

        Run against the matcher directly, since the detector entry points now
        refuse the tolerance before they reach it.
        """
        detect = load_symmetry_detect_module()
        charges = np.asarray(self.CHARGES, dtype=float)
        centered = self.COORDS - np.einsum(
            'i,ij->j', charges, self.COORDS) / float(np.sum(charges))
        # A quarter turn about z is not a symmetry operation of this geometry.
        quarter_turn = detect._rotation_matrix(
            np.array([0.0, 0.0, 1.0]), np.pi / 2.0)
        self.assertIsNone(detect._match_permutation(
            charges, centered, centered @ quarter_turn.T, 1.0e-5))
        self.assertIsNotNone(detect._match_permutation(
            charges, centered, centered @ quarter_turn.T, float('nan')))

    def test_a_non_finite_or_non_positive_tolerance_returns_one(self):
        detect = load_symmetry_detect_module()
        for value in (float('nan'), float('inf'), 0.0, -1.0e-5, 'loose'):
            with self.subTest(tolerance=value):
                self.assertEqual(
                    detect.rotational_symmetry_number(
                        self.CHARGES, self.COORDS, tolerance=value), 1)

    def test_the_detector_entry_points_reject_it_rather_than_detect(self):
        """A strict-symmetry run detects before the input checker runs, so the
        guard has to live in the detector too or the user gets a misleading
        point-group mismatch instead of a tolerance diagnostic."""
        detect = load_symmetry_detect_module()
        for value in (float('nan'), float('inf'), 0.0):
            with self.subTest(tolerance=value):
                with self.assertRaises(ValueError):
                    detect.detect_point_group(self.CHARGES, self.COORDS,
                                              tolerance=value)
                with self.assertRaises(ValueError):
                    detect.enumerate_full_group(self.CHARGES, self.COORDS,
                                                tolerance=value)

    def test_a_usable_tolerance_is_untouched(self):
        detect = load_symmetry_detect_module()
        coords = np.asarray(WATER[2], dtype=float) * ANGSTROM_TO_BOHR
        self.assertEqual(
            detect.rotational_symmetry_number([8, 1, 1], coords,
                                              tolerance=1.0e-5), 2)


class TestPositionalArgumentsKeepTheirMeaning(unittest.TestCase):
    """`sigma` was inserted between `linear` and `mult`, re-aiming every
    positional call: thermal_analysis(..., 298.15, False, 1) meant mult=1 and
    became sigma=1 with mult=0, so st_el = R*T*ln(0) = -inf."""

    def test_sigma_is_keyword_only(self):
        import inspect
        frequency = load_frequency_module()
        parameters = inspect.signature(frequency.thermal_analysis).parameters
        self.assertEqual(parameters['sigma'].kind,
                         inspect.Parameter.KEYWORD_ONLY)
        positional = [name for name, p in parameters.items()
                      if p.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD]
        self.assertEqual(positional,
                         ['energy', 'atoms', 'mass', 'freqs', 'inertia',
                          'temperature', 'linear', 'mult',
                          'freq_scale_factor', 'freq_cutoff'])

    def test_a_pre_existing_positional_call_still_sets_the_multiplicity(self):
        frequency = load_frequency_module()
        atoms, mass, geometry, freqs = WATER
        coords = np.asarray(geometry, dtype=float) * ANGSTROM_TO_BOHR
        _, _, inertia = frequency.normal_mode(
            coords.ravel(), np.asarray(mass, dtype=float),
            np.zeros((3 * len(atoms), 3 * len(atoms))))

        data = frequency.thermal_analysis(
            0.0, atoms, np.asarray(mass, dtype=float),
            np.asarray(freqs, dtype=float), inertia, 298.15, False, 1)

        # mult = 1 -> st_el = R*T*ln(1) = 0. Pre-fix this was -inf.
        self.assertTrue(np.isfinite(data['st_el']))
        self.assertEqual(data['st_el'], 0.0)
        self.assertEqual(data['sigma'], 1)


if __name__ == '__main__':
    unittest.main()
