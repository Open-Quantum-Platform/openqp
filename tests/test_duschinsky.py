import importlib.util
import subprocess
import sys
import textwrap
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_duschinsky_module():
    name = "openqp_duschinsky_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/library/duschinsky.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def pair_spring_hessian(geometry, spring_constants):
    """Independent Cartesian Hessian for pair-distance harmonic potentials."""

    geometry = np.asarray(geometry, dtype=float)
    hessian = np.zeros((3 * len(geometry), 3 * len(geometry)))
    for (first, second), force_constant in spring_constants.items():
        direction = geometry[first] - geometry[second]
        direction /= np.linalg.norm(direction)
        block = force_constant * np.outer(direction, direction)
        first_slice = slice(3 * first, 3 * first + 3)
        second_slice = slice(3 * second, 3 * second + 3)
        hessian[first_slice, first_slice] += block
        hessian[second_slice, second_slice] += block
        hessian[first_slice, second_slice] -= block
        hessian[second_slice, first_slice] -= block
    return hessian


def axis_angle_rotation(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    cross = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return (
        np.eye(3) * np.cos(angle)
        + (1.0 - np.cos(angle)) * np.outer(axis, axis)
        + np.sin(angle) * cross
    )


class FakeMolecule:
    """Minimal stand-in for the public OpenQP Molecule accessors."""

    def __init__(self, geometry, masses, atoms, hessian):
        self._geometry = np.asarray(geometry, dtype=float)
        self._masses = np.asarray(masses, dtype=float)
        self._atoms = np.asarray(atoms)
        self._hessian = np.asarray(hessian, dtype=float)

    def get_system(self):
        return self._geometry.reshape(-1).copy()

    def get_mass(self):
        return self._masses.copy()

    def get_atoms(self):
        return self._atoms.copy()

    def get_hess(self):
        return self._hessian.copy()


class TestProjectedNormalModes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.duschinsky = load_duschinsky_module()

    def test_diatomic_has_five_external_and_one_physical_mode(self):
        geometry = np.array([[-0.7, 0.0, 0.0], [0.7, 0.0, 0.0]])
        masses = np.array([2.0, 2.0])
        hessian = pair_spring_hessian(geometry, {(0, 1): 3.0})

        modes = self.duschinsky.projected_normal_modes(
            geometry, masses, hessian
        )

        self.assertTrue(modes.is_linear)
        self.assertEqual(modes.external_rank, 5)
        self.assertEqual(modes.force_constants.shape, (1,))
        self.assertAlmostEqual(modes.force_constants[0], 3.0, places=12)
        np.testing.assert_allclose(
            modes.mass_weighted_modes.T @ modes.mass_weighted_modes,
            np.eye(1),
            atol=1.0e-13,
        )
        self.assertLess(modes.external_contamination, 1.0e-13)

    def test_nonlinear_triatomic_retains_three_physical_modes(self):
        geometry = np.array(
            [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [-0.2, 0.95, 0.3]]
        )
        masses = np.array([16.0, 1.0, 2.0])
        hessian = pair_spring_hessian(
            geometry, {(0, 1): 1.1, (0, 2): 1.7, (1, 2): 2.3}
        )

        modes = self.duschinsky.projected_normal_modes(
            geometry, masses, hessian
        )

        self.assertFalse(modes.is_linear)
        self.assertEqual(modes.external_rank, 6)
        self.assertEqual(modes.force_constants.shape, (3,))
        self.assertTrue(np.all(modes.force_constants > 0.0))
        np.testing.assert_allclose(
            modes.mass_weighted_modes.T @ modes.mass_weighted_modes,
            np.eye(3),
            atol=2.0e-13,
        )
        coordinate_masses = np.repeat(masses, 3)
        np.testing.assert_allclose(
            modes.cartesian_modes.T
            @ (coordinate_masses[:, None] * modes.cartesian_modes),
            np.eye(3),
            atol=2.0e-13,
        )
        self.assertLess(modes.external_contamination, 2.0e-13)

    def test_negative_force_constant_is_retained(self):
        geometry = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
        masses = np.ones(2)
        hessian = pair_spring_hessian(geometry, {(0, 1): -0.4})

        modes = self.duschinsky.projected_normal_modes(
            geometry, masses, hessian
        )

        self.assertAlmostEqual(modes.force_constants[0], -0.8, places=12)

    def test_external_rank_is_independent_of_coordinate_length_unit(self):
        geometry = np.array(
            [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [-0.2, 0.95, 0.3]]
        )
        masses = np.array([16.0, 1.0, 2.0])
        hessian = pair_spring_hessian(
            geometry, {(0, 1): 1.1, (0, 2): 1.7, (1, 2): 2.3}
        )

        ordinary = self.duschinsky.projected_normal_modes(
            geometry, masses, hessian
        )
        tiny_length_unit = self.duschinsky.projected_normal_modes(
            geometry * 1.0e-12, masses, hessian
        )

        self.assertEqual(ordinary.external_rank, 6)
        self.assertEqual(tiny_length_unit.external_rank, 6)
        np.testing.assert_allclose(
            tiny_length_unit.force_constants,
            ordinary.force_constants,
            atol=2.0e-14,
        )

    def test_asymmetric_hessian_fails_closed(self):
        geometry = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
        hessian = np.zeros((6, 6))
        hessian[0, 1] = 1.0e-4

        with self.assertRaisesRegex(ValueError, "not symmetric"):
            self.duschinsky.projected_normal_modes(
                geometry, np.ones(2), hessian
            )

    def test_nonpositive_mass_is_rejected(self):
        geometry = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])

        with self.assertRaisesRegex(ValueError, "must be positive"):
            self.duschinsky.projected_normal_modes(
                geometry, np.array([1.0, 0.0]), np.zeros((6, 6))
            )


class TestModeGaugeAndDiagnostics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.duschinsky = load_duschinsky_module()

    def test_phase_and_degenerate_subspace_alignment(self):
        ground = np.eye(4)[:, :3]
        angle = 0.63
        excited = np.column_stack(
            (
                np.cos(angle) * ground[:, 0] + np.sin(angle) * ground[:, 1],
                -np.sin(angle) * ground[:, 0] + np.cos(angle) * ground[:, 1],
                -ground[:, 2],
            )
        )

        alignment = self.duschinsky.align_mode_gauge(
            ground,
            excited,
            np.array([1.0, 1.0, 3.0]),
            np.array([2.0, 2.0, 4.0]),
        )

        self.assertEqual(alignment.aligned_degenerate_blocks, ((0, 2),))
        np.testing.assert_allclose(alignment.modes, ground, atol=2.0e-15)
        np.testing.assert_allclose(
            alignment.transformation.T @ alignment.transformation,
            np.eye(3),
            atol=2.0e-15,
        )
        np.testing.assert_allclose(
            alignment.overlap_with_ground, np.eye(3), atol=2.0e-15
        )

    def test_nondegenerate_mixing_is_not_rotated_away(self):
        ground = np.eye(3)
        angle = 0.4
        excited = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )

        alignment = self.duschinsky.align_mode_gauge(
            ground,
            excited,
            np.array([1.0, 2.0, 3.0]),
            np.array([1.1, 2.1, 3.1]),
        )

        np.testing.assert_allclose(alignment.modes, excited, atol=1.0e-15)
        self.assertGreater(abs(alignment.overlap_with_ground[0, 1]), 0.3)

    def test_phase_uses_largest_overlap_when_diagonal_overlap_is_zero(self):
        ground = np.eye(3)
        excited = np.column_stack((-ground[:, 1], ground[:, 0], ground[:, 2]))

        alignment = self.duschinsky.align_mode_gauge(
            ground,
            excited,
            np.array([1.0, 2.0, 3.0]),
            np.array([1.1, 2.1, 3.1]),
        )

        self.assertGreater(alignment.overlap_with_ground[0, 1], 0.0)
        self.assertGreater(alignment.overlap_with_ground[1, 0], 0.0)

    def test_orthogonality_diagnostics_are_two_sided(self):
        proper = axis_angle_rotation([1.0, 2.0, -1.0], 0.7)
        proper_diagnostics = self.duschinsky.orthogonality_diagnostics(proper)
        self.assertLess(proper_diagnostics.max_abs, 5.0e-16)
        np.testing.assert_allclose(
            proper_diagnostics.singular_values, np.ones(3), atol=5.0e-16
        )
        self.assertAlmostEqual(proper_diagnostics.determinant, 1.0, places=14)

        nonorthogonal = np.diag([1.0, 0.8])
        diagnostics = self.duschinsky.orthogonality_diagnostics(nonorthogonal)
        self.assertAlmostEqual(diagnostics.left_max_abs, 0.36, places=14)
        self.assertAlmostEqual(diagnostics.right_max_abs, 0.36, places=14)
        np.testing.assert_allclose(diagnostics.singular_values, [1.0, 0.8])


class TestDuschinskyTransformation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.duschinsky = load_duschinsky_module()
        cls.ground_geometry = np.array(
            [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [-0.2, 0.95, 0.3]]
        )
        cls.masses = np.array([16.0, 1.0, 2.0])
        cls.ground_hessian = pair_spring_hessian(
            cls.ground_geometry,
            {(0, 1): 1.1, (0, 2): 1.7, (1, 2): 2.3},
        )

    def test_identical_surfaces_give_identity_and_zero_shift(self):
        result = self.duschinsky.compute_duschinsky(
            self.ground_geometry,
            self.ground_geometry,
            self.masses,
            self.ground_hessian,
            self.ground_hessian,
            orthogonality_tolerance=1.0e-12,
        )

        np.testing.assert_allclose(result.J, np.eye(3), atol=3.0e-14)
        np.testing.assert_allclose(result.K, np.zeros(3), atol=1.0e-14)
        self.assertLess(result.orthogonality.max_abs, 5.0e-14)
        self.assertLess(result.mass_weighted_rmsd, 1.0e-14)
        np.testing.assert_allclose(result.eckart_residual, 0.0, atol=1.0e-14)

    def test_rigid_rotation_and_translation_leave_J_and_K_invariant(self):
        rotation_ground_to_moving = axis_angle_rotation([1.0, -2.0, 0.5], 0.91)
        moving = self.ground_geometry @ rotation_ground_to_moving + np.array(
            [2.3, -1.7, 0.8]
        )
        cartesian_rotation = np.kron(
            np.eye(len(self.masses)), rotation_ground_to_moving.T
        )
        moving_hessian = (
            cartesian_rotation
            @ self.ground_hessian
            @ cartesian_rotation.T
        )

        result = self.duschinsky.compute_duschinsky(
            self.ground_geometry,
            moving,
            self.masses,
            self.ground_hessian,
            moving_hessian,
            orthogonality_tolerance=1.0e-11,
        )

        self.assertAlmostEqual(
            np.linalg.det(result.rotation_excited_to_ground), 1.0, places=13
        )
        np.testing.assert_allclose(
            result.aligned_excited_geometry,
            self.ground_geometry,
            atol=3.0e-14,
        )
        np.testing.assert_allclose(result.J, np.eye(3), atol=3.0e-13)
        np.testing.assert_allclose(result.K, np.zeros(3), atol=2.0e-13)
        np.testing.assert_allclose(result.eckart_residual, 0.0, atol=2.0e-13)

    def test_diatomic_shift_has_correct_mass_weighting_and_coordinate_identity(self):
        ground = np.array([[-0.7, 0.0, 0.0], [0.7, 0.0, 0.0]])
        excited = np.array([[-0.9, 0.0, 0.0], [0.9, 0.0, 0.0]])
        masses = np.array([2.0, 2.0])
        ground_hessian = pair_spring_hessian(ground, {(0, 1): 2.0})
        excited_hessian = pair_spring_hessian(excited, {(0, 1): 3.0})

        result = self.duschinsky.compute_duschinsky(
            ground, excited, masses, ground_hessian, excited_hessian
        )

        np.testing.assert_allclose(result.J, [[1.0]], atol=2.0e-14)
        self.assertAlmostEqual(abs(result.K[0]), 0.4, places=13)

        ground_coordinate = np.array([0.17])
        coordinate_masses = np.repeat(masses, 3)
        point = ground + (
            result.ground_modes.mass_weighted_modes @ ground_coordinate
            / np.sqrt(coordinate_masses)
        ).reshape(-1, 3)
        excited_coordinate = (
            result.excited_modes.mass_weighted_modes.T
            @ (
                np.sqrt(coordinate_masses)
                * (point - result.aligned_excited_geometry).reshape(-1)
            )
        )
        np.testing.assert_allclose(
            excited_coordinate,
            result.J @ ground_coordinate + result.K,
            atol=2.0e-14,
        )

    def test_deformed_nonlinear_geometry_obeys_the_coordinate_identity(self):
        excited = self.ground_geometry.copy()
        excited[1] += np.array([0.3, 0.4, -0.2])
        excited_hessian = pair_spring_hessian(
            excited, {(0, 1): 1.4, (0, 2): 1.2, (1, 2): 2.0}
        )
        result = self.duschinsky.compute_duschinsky(
            self.ground_geometry,
            excited,
            self.masses,
            self.ground_hessian,
            excited_hessian,
        )

        ground_coordinate = np.array([0.13, -0.27, 0.08])
        coordinate_masses = np.repeat(self.masses, 3)
        point = self.ground_geometry + (
            result.ground_modes.mass_weighted_modes @ ground_coordinate
            / np.sqrt(coordinate_masses)
        ).reshape(-1, 3)
        excited_coordinate = (
            result.excited_modes.mass_weighted_modes.T
            @ (
                np.sqrt(coordinate_masses)
                * (point - result.aligned_excited_geometry).reshape(-1)
            )
        )

        np.testing.assert_allclose(
            excited_coordinate,
            result.J @ ground_coordinate + result.K,
            atol=3.0e-14,
        )
        self.assertGreater(result.orthogonality.max_abs, 1.0e-3)
        self.assertLess(result.orthogonality.singular_values[-1], 0.999)

    def test_linear_to_nonlinear_change_is_rejected(self):
        ground = np.array([[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        excited = ground.copy()
        excited[1, 1] = 0.2
        masses = np.array([1.0, 16.0, 1.0])

        with self.assertRaisesRegex(ValueError, "different physical-mode counts"):
            self.duschinsky.compute_duschinsky(
                ground,
                excited,
                masses,
                np.zeros((9, 9)),
                np.zeros((9, 9)),
            )

    def test_optional_orthogonality_threshold_fails_closed(self):
        excited = self.ground_geometry.copy()
        excited[1] += np.array([0.3, 0.4, -0.2])
        excited_hessian = pair_spring_hessian(
            excited, {(0, 1): 1.4, (0, 2): 1.2, (1, 2): 2.0}
        )

        with self.assertRaisesRegex(ValueError, "orthogonality residual"):
            self.duschinsky.compute_duschinsky(
                self.ground_geometry,
                excited,
                self.masses,
                self.ground_hessian,
                excited_hessian,
                orthogonality_tolerance=1.0e-14,
            )


class TestPublicDuschinskyAPI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.duschinsky = load_duschinsky_module()
        cls.geometry = np.array(
            [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [-0.2, 0.95, 0.3]]
        )
        cls.masses = np.array([16.0, 1.0, 2.0])
        cls.atoms = np.array([8, 1, 1])
        cls.hessian = pair_spring_hessian(
            cls.geometry,
            {(0, 1): 1.1, (0, 2): 1.7, (1, 2): 2.3},
        )

    def test_explicit_state_arrays_support_a_linear_molecule(self):
        geometry = np.array([[-0.7, 0.0, 0.0], [0.7, 0.0, 0.0]])
        masses = np.array([2.0, 2.0])
        atoms = np.array([1, 1])
        hessian = pair_spring_hessian(geometry, {(0, 1): 3.0})

        result = self.duschinsky.compute_duschinsky_from_arrays(
            geometry,
            masses,
            hessian,
            geometry,
            masses.copy(),
            hessian.copy(),
            ground_atoms=atoms,
            excited_atoms=atoms.copy(),
        )

        self.assertTrue(result.ground_modes.is_linear)
        np.testing.assert_allclose(result.J, [[1.0]], atol=2.0e-14)
        np.testing.assert_allclose(result.K, [0.0], atol=1.0e-14)

    def test_two_molecules_use_geometry_mass_atom_and_hessian_accessors(self):
        ground = FakeMolecule(
            self.geometry, self.masses, self.atoms, self.hessian
        )
        excited = FakeMolecule(
            self.geometry, self.masses, self.atoms, self.hessian
        )

        result = self.duschinsky.compute_duschinsky_from_molecules(
            ground, excited, orthogonality_tolerance=1.0e-12
        )

        self.assertFalse(result.ground_modes.is_linear)
        np.testing.assert_allclose(result.J, np.eye(3), atol=3.0e-14)
        np.testing.assert_allclose(result.K, np.zeros(3), atol=1.0e-14)

    def test_explicit_hessians_replace_empty_molecule_results(self):
        empty_hessian = np.zeros(0)
        ground = FakeMolecule(
            self.geometry, self.masses, self.atoms, empty_hessian
        )
        excited = FakeMolecule(
            self.geometry, self.masses, self.atoms, empty_hessian
        )

        result = self.duschinsky.compute_duschinsky_from_molecules(
            ground,
            excited,
            ground_hessian=self.hessian,
            excited_hessian=self.hessian,
        )

        np.testing.assert_allclose(result.J, np.eye(3), atol=3.0e-14)

    def test_empty_molecule_hessian_fails_with_actionable_message(self):
        molecule = FakeMolecule(
            self.geometry, self.masses, self.atoms, np.zeros(0)
        )

        with self.assertRaisesRegex(ValueError, "pass ground_hessian explicitly"):
            self.duschinsky.compute_duschinsky_from_molecules(
                molecule, molecule
            )

    def test_mass_mismatch_is_rejected_before_mode_analysis(self):
        excited_masses = self.masses.copy()
        excited_masses[-1] += 0.01

        with self.assertRaisesRegex(ValueError, "atomic mass mismatch"):
            self.duschinsky.compute_duschinsky_from_arrays(
                self.geometry,
                self.masses,
                self.hessian,
                self.geometry,
                excited_masses,
                self.hessian,
                ground_atoms=self.atoms,
                excited_atoms=self.atoms,
            )

    def test_atom_identity_or_order_mismatch_is_rejected(self):
        excited_atoms = self.atoms[[1, 0, 2]].copy()

        with self.assertRaisesRegex(ValueError, "atom identity or order mismatch"):
            self.duschinsky.compute_duschinsky_from_arrays(
                self.geometry,
                self.masses,
                self.hessian,
                self.geometry,
                self.masses,
                self.hessian,
                ground_atoms=self.atoms,
                excited_atoms=excited_atoms,
            )

    def test_public_library_import_exports_both_helpers(self):
        package_root = ROOT / "pyoqp/oqp"
        script = textwrap.dedent(
            f"""
            import sys
            import types

            oqp = types.ModuleType("oqp")
            oqp.__path__ = [{str(package_root)!r}]
            sys.modules["oqp"] = oqp

            for leaf in ("guess", "ints_1e", "ints_2e", "set_basis", "project_basis"):
                module = types.ModuleType(f"oqp.library.{{leaf}}")
                module.__all__ = []
                sys.modules[module.__name__] = module

            odp = types.ModuleType("oqp.library.odp")
            odp.ODPUmbrella = object
            odp.odp_wham = lambda: None
            odp.write_odp_wham = lambda: None
            sys.modules[odp.__name__] = odp

            from oqp.library import (
                compute_duschinsky_from_arrays,
                compute_duschinsky_from_molecules,
            )

            assert callable(compute_duschinsky_from_arrays)
            assert callable(compute_duschinsky_from_molecules)
            """
        )
        subprocess.run([sys.executable, "-c", script], check=True)


if __name__ == "__main__":
    unittest.main()
