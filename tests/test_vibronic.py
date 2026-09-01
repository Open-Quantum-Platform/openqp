import importlib.util
import json
from math import factorial
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_vibronic_module():
    name = "openqp_vibronic_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/library/vibronic.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class TestHarmonicOverlap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()

    def model(self, *, frequency=1000.0, displacement=0.0):
        return self.vibronic.HarmonicVibronicModel.create(
            [frequency],
            [frequency],
            [[1.0]],
            [displacement],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic one-mode phase",
        )

    def test_identical_surfaces_have_kronecker_delta_overlaps(self):
        engine = self.vibronic.HarmonicOverlapEngine(self.model())
        for initial in range(5):
            for final in range(5):
                expected = 1.0 if initial == final else 0.0
                self.assertAlmostEqual(
                    engine.overlap((initial,), (final,)), expected, places=13
                )

    def test_displaced_equal_frequency_ground_progression_is_poisson(self):
        frequency = 1200.0
        displacement = 8.0
        model = self.model(frequency=frequency, displacement=displacement)
        engine = self.vibronic.HarmonicOverlapEngine(model)
        huang_rhys = (
            frequency
            * self.vibronic.CM1_TO_HARTREE
            * displacement**2
            / 2.0
        )
        for final in range(6):
            expected = (
                np.exp(-huang_rhys)
                * huang_rhys**final
                / factorial(final)
            )
            self.assertAlmostEqual(
                engine.overlap((0,), (final,)) ** 2, expected, places=13
            )

    def test_frequency_change_zero_zero_factor_matches_analytic_limit(self):
        model = self.vibronic.HarmonicVibronicModel.create(
            [800.0],
            [1800.0],
            [[1.0]],
            [0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic one-mode phase",
        )
        engine = self.vibronic.HarmonicOverlapEngine(model)
        expected = 2.0 * np.sqrt(800.0 * 1800.0) / (800.0 + 1800.0)
        self.assertAlmostEqual(engine.overlap((0,), (0,)) ** 2, expected, places=13)

    def test_amu_coordinate_displacement_is_converted_to_electron_mass(self):
        model = self.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [1000.0],
            [[1.0]],
            [0.2],
            coordinate_unit="sqrt(amu)*bohr",
            coordinate_phase_convention="synthetic amu coordinate",
        )
        engine = self.vibronic.HarmonicOverlapEngine(model)
        huang_rhys = (
            1000.0
            * self.vibronic.CM1_TO_HARTREE
            * self.vibronic.AMU_TO_ELECTRON_MASS
            * 0.2**2
            / 2.0
        )
        self.assertAlmostEqual(
            engine.overlap((0,), (0,)) ** 2,
            np.exp(-huang_rhys),
            places=13,
        )

    def test_excited_mode_phase_change_preserves_fc_factors(self):
        first = self.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [1200.0],
            [[1.0]],
            [4.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="first phase",
        )
        second = self.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [1200.0],
            [[-1.0]],
            [-4.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="opposite excited phase",
        )
        first_engine = self.vibronic.HarmonicOverlapEngine(first)
        second_engine = self.vibronic.HarmonicOverlapEngine(second)
        for final in range(5):
            self.assertAlmostEqual(
                first_engine.overlap((0,), (final,)) ** 2,
                second_engine.overlap((0,), (final,)) ** 2,
                places=13,
            )

    def test_nonorthogonal_duschinsky_matrix_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "orthogonality residual"):
            self.vibronic.HarmonicVibronicModel.create(
                [1000.0, 1200.0],
                [900.0, 1300.0],
                [[1.0, 0.0], [0.0, 0.98]],
                [0.0, 0.0],
                coordinate_unit="sqrt(amu)*bohr",
                coordinate_phase_convention="synthetic",
            )

    def test_two_mode_rotation_preserves_overlap_normalization(self):
        angle = 0.37
        rotation = np.array(
            [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        )
        model = self.vibronic.HarmonicVibronicModel.create(
            [900.0, 1400.0],
            [900.0, 1400.0],
            rotation,
            [0.0, 0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="rotated synthetic modes",
        )
        engine = self.vibronic.HarmonicOverlapEngine(model)
        states = self.vibronic.enumerate_vibrational_states(2, 12)
        total = sum(engine.overlap((0, 0), state) ** 2 for state in states)
        self.assertGreater(total, 0.999999)
        self.assertLessEqual(total, 1.0 + 2.0e-12)


class TestVibronicSpectrumAndHT(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()
        cls.model = cls.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [1000.0],
            [[1.0]],
            [0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="shared normal-mode phase",
        )

    def transition(self, value=(1.0, 0.0, 0.0)):
        return self.vibronic.ElectronicTransitionMoment.create(
            value,
            electronic_state="S1",
            electronic_phase_convention="shared electronic-state phase",
        )

    def derivative(self, value):
        return self.vibronic.ExcitedStatePropertyDerivative.create(
            np.asarray(value).reshape(3, 1),
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

    def test_condon_spectrum_and_grid_are_normalized(self):
        result = self.vibronic.harmonic_vibronic_spectrum(
            self.model,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            max_final_quanta=3,
            transition=self.transition(),
            fwhm_cm1=80.0,
            normalization="sum",
        )
        self.assertAlmostEqual(result.lines[0].raw_strength, 1.0, places=13)
        self.assertTrue(all(line.raw_strength < 1.0e-28 for line in result.lines[1:]))
        integration = getattr(np, "trapezoid", np.trapz if hasattr(np, "trapz") else None)
        self.assertIsNotNone(integration)
        self.assertAlmostEqual(
            integration(result.intensity, result.wavenumbers_cm1), 1.0, places=12
        )
        self.assertAlmostEqual(result.franck_condon_completeness, 1.0, places=13)

    def test_nonzero_temperature_requires_explicit_initial_truncation(self):
        with self.assertRaisesRegex(ValueError, "max_initial_quanta is required"):
            self.vibronic.harmonic_vibronic_spectrum(
                self.model,
                electronic_origin_cm1=20000.0,
                origin_kind="zero_zero",
                max_final_quanta=2,
                temperature_kelvin=300.0,
            )

    def test_thermal_population_is_checked_against_exact_partition_function(self):
        result = self.vibronic.harmonic_vibronic_spectrum(
            self.model,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            max_final_quanta=8,
            temperature_kelvin=1000.0,
            max_initial_quanta=8,
            minimum_thermal_population=0.999,
        )
        self.assertGreater(result.retained_thermal_population, 0.999)

    def test_incomplete_final_basis_is_rejected_before_normalization(self):
        displaced = self.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [900.0],
            [[1.0]],
            [50.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="strongly displaced synthetic mode",
        )
        with self.assertRaisesRegex(ValueError, "Franck-Condon completeness"):
            self.vibronic.harmonic_vibronic_spectrum(
                displaced,
                electronic_origin_cm1=20000.0,
                origin_kind="zero_zero",
                max_final_quanta=0,
                normalization="sum",
            )

        partial = self.vibronic.harmonic_vibronic_spectrum(
            displaced,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            max_final_quanta=0,
            normalization="sum",
            minimum_franck_condon_completeness=0.0,
        )
        self.assertLess(partial.franck_condon_completeness, 0.01)

    def test_ht_fundamental_matches_ladder_operator_limit(self):
        derivative = self.derivative((2.0, 0.0, 0.0))
        engine = self.vibronic.HarmonicOverlapEngine(self.model)
        moment = self.vibronic.vibronic_transition_moment(
            engine,
            (0,),
            (1,),
            self.transition((0.0, 0.0, 0.0)),
            derivative,
        )
        omega = 1000.0 * self.vibronic.CM1_TO_HARTREE
        self.assertAlmostEqual(moment[0].real, 2.0 / np.sqrt(2.0 * omega), places=12)
        self.assertAlmostEqual(moment[1].real, 0.0, places=13)

    def test_ht_derivative_requires_explicit_zero_transition_dipole(self):
        with self.assertRaisesRegex(ValueError, "explicit equilibrium transition dipole"):
            self.vibronic.harmonic_vibronic_spectrum(
                self.model,
                electronic_origin_cm1=20000.0,
                origin_kind="zero_zero",
                max_final_quanta=1,
                transition_dipole_derivative=self.derivative((1.0, 0.0, 0.0)),
            )

    def test_electronic_phase_mismatch_is_rejected(self):
        wrong_derivative = self.vibronic.ExcitedStatePropertyDerivative.create(
            [[1.0], [0.0], [0.0]],
            property_kind="transition_dipole",
            coordinate_basis="ground_normal",
            coordinate_unit="sqrt(me)*bohr",
            property_unit="e*bohr",
            electronic_state="S1",
            electronic_phase_convention="different electronic phase",
            coordinate_phase_convention="shared normal-mode phase",
            provenance="synthetic derivative",
            electronic_state_role="target_excited_state",
        )
        with self.assertRaisesRegex(ValueError, "electronic phases"):
            self.vibronic.harmonic_vibronic_spectrum(
                self.model,
                electronic_origin_cm1=20000.0,
                origin_kind="zero_zero",
                max_final_quanta=1,
                transition=self.transition(),
                transition_dipole_derivative=wrong_derivative,
            )

    def test_normal_coordinate_phase_mismatch_is_rejected(self):
        wrong_derivative = self.vibronic.ExcitedStatePropertyDerivative.create(
            [[1.0], [0.0], [0.0]],
            property_kind="transition_dipole",
            coordinate_basis="ground_normal",
            coordinate_unit="sqrt(me)*bohr",
            property_unit="e*bohr",
            electronic_state="S1",
            electronic_phase_convention="shared electronic-state phase",
            coordinate_phase_convention="opposite normal-mode phase",
            provenance="synthetic derivative",
            electronic_state_role="target_excited_state",
        )
        with self.assertRaisesRegex(ValueError, "normal-coordinate phases"):
            self.vibronic.harmonic_vibronic_spectrum(
                self.model,
                electronic_origin_cm1=20000.0,
                origin_kind="zero_zero",
                max_final_quanta=1,
                transition=self.transition(),
                transition_dipole_derivative=wrong_derivative,
            )

    def test_reference_state_property_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "target_excited_state"):
            self.vibronic.ExcitedStatePropertyDerivative.create(
                np.ones((3, 1)),
                property_kind="state_dipole",
                coordinate_basis="excited_normal",
                coordinate_unit="sqrt(amu)*bohr",
                property_unit="e*bohr",
                electronic_state="ROHF reference",
                electronic_phase_convention="reference electronic phase",
                coordinate_phase_convention="reference normal phase",
                provenance="ROHF property",
                electronic_state_role="reference_state",
            )

    def test_excited_state_role_cannot_be_omitted(self):
        with self.assertRaises(TypeError):
            self.vibronic.ExcitedStatePropertyDerivative.create(
                np.ones((3, 1)),
                property_kind="state_dipole",
                coordinate_basis="excited_normal",
                coordinate_unit="sqrt(amu)*bohr",
                property_unit="e*bohr",
                electronic_state="S1",
                electronic_phase_convention="phase-invariant expectation value",
                coordinate_phase_convention="tracked S1 normal modes",
                provenance="unspecified calculation",
            )


class TestPropertyDerivativeObservables(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()

    def property(self, values, kind, unit):
        return self.vibronic.ExcitedStatePropertyDerivative.create(
            values,
            property_kind=kind,
            coordinate_basis="excited_normal",
            coordinate_unit="sqrt(amu)*bohr",
            property_unit=unit,
            electronic_state="S1",
            electronic_phase_convention="phase-invariant expectation value",
            coordinate_phase_convention="tracked S1 normal modes",
            provenance="synthetic analytic derivative",
            electronic_state_role="target_excited_state",
        )

    def test_ir_intensity_is_squared_mode_dipole_norm(self):
        derivative = self.property(
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 3.0, 0.0]],
            "state_dipole",
            "e*bohr",
        )
        result = self.vibronic.excited_state_ir_intensities(derivative)
        np.testing.assert_allclose(
            result.intensities_km_mol,
            self.vibronic.IR_INTENSITY_CONVERSION_KM_MOL
            * np.array([9.0, 9.0, 0.0]),
        )
        self.assertAlmostEqual(
            self.vibronic.IR_INTENSITY_CONVERSION_KM_MOL,
            974.88011,
            places=5,
        )

    def test_raman_isotropic_and_anisotropic_analytic_limits(self):
        values = np.zeros((3, 3, 2))
        values[:, :, 0] = 2.0 * np.eye(3)
        values[:, :, 1] = np.diag([1.0, -1.0, 0.0])
        result = self.vibronic.excited_state_raman_activities(
            self.property(values, "polarizability", "bohr^3")
        )
        np.testing.assert_allclose(result.activities_au, [180.0, 21.0])
        np.testing.assert_allclose(result.depolarization_ratio_polarized, [0.0, 0.75])
        np.testing.assert_allclose(
            result.depolarization_ratio_unpolarized, [0.0, 6.0 / 7.0]
        )

    def test_finite_difference_aligns_transition_dipole_phases(self):
        plus = np.array([[[1.2, 0.0, 0.0]]]).reshape(1, 3)
        minus = np.array([[[-0.8, 0.0, 0.0]]]).reshape(1, 3)
        derivative = self.vibronic.finite_difference_excited_state_property(
            plus,
            minus,
            [0.1],
            property_kind="transition_dipole",
            coordinate_basis="ground_normal",
            coordinate_unit="sqrt(me)*bohr",
            property_unit="e*bohr",
            electronic_state="S1",
            electronic_phase_convention="overlap-aligned S1",
            coordinate_phase_convention="ground normal-mode phases",
            state_tracking_overlaps=[[0.99, 0.98]],
            transition_phase_factors=[[1.0, -1.0]],
        )
        np.testing.assert_allclose(derivative.values[:, 0], [2.0, 0.0, 0.0])
        self.assertAlmostEqual(derivative.minimum_state_overlap, 0.98)

    def test_finite_difference_rejects_untracked_state(self):
        with self.assertRaisesRegex(ValueError, "state-tracking overlap"):
            self.vibronic.finite_difference_excited_state_property(
                [[1.0, 0.0, 0.0]],
                [[0.8, 0.0, 0.0]],
                [0.1],
                property_kind="state_dipole",
                coordinate_basis="excited_normal",
                coordinate_unit="sqrt(amu)*bohr",
                property_unit="e*bohr",
                electronic_state="S1",
                electronic_phase_convention="phase-invariant expectation value",
                coordinate_phase_convention="tracked S1 normal modes",
                state_tracking_overlaps=[[0.95, 0.4]],
            )

    def test_finite_difference_rejects_nonphysical_overlap(self):
        with self.assertRaisesRegex(ValueError, "must lie in"):
            self.vibronic.finite_difference_excited_state_property(
                [[1.0, 0.0, 0.0]],
                [[0.8, 0.0, 0.0]],
                [0.1],
                property_kind="state_dipole",
                coordinate_basis="excited_normal",
                coordinate_unit="sqrt(amu)*bohr",
                property_unit="e*bohr",
                electronic_state="S1",
                electronic_phase_convention="phase-invariant expectation value",
                coordinate_phase_convention="tracked S1 normal modes",
                state_tracking_overlaps=[[1.01, 0.99]],
            )

    def test_mrsf_analytic_provider_hook_requires_complete_target_state_record(self):
        class Provider:
            def get_mrsf_analytic_property_derivative(
                self, *, property_kind, electronic_state
            ):
                return {
                    "values": [[1.0], [0.0], [0.0]],
                    "property_kind": property_kind,
                    "coordinate_basis": "excited_normal",
                    "coordinate_unit": "sqrt(amu)*bohr",
                    "property_unit": "e*bohr",
                    "electronic_state": electronic_state,
                    "electronic_state_role": "target_excited_state",
                    "electronic_phase_convention": (
                        "phase-invariant expectation value"
                    ),
                    "coordinate_phase_convention": "tracked S1 normal modes",
                    "electronic_method": "MRSF-TDDFT",
                    "derivative_order": 1,
                    "provenance": "analytic MRSF-TDDFT dipole derivative",
                }

        derivative = self.vibronic.mrsf_analytic_property_derivative_from_provider(
            Provider(), property_kind="state_dipole", electronic_state="S1"
        )
        self.assertEqual(derivative.electronic_state_role, "target_excited_state")
        self.assertIn("analytic", derivative.provenance)
        json.dumps(derivative.to_dict())

    def test_mrsf_provider_hook_never_falls_back_to_reference_property(self):
        class ReferenceOnlyProvider:
            def get_dipole_derivatives(self):
                raise AssertionError("reference accessor must not be called")

        with self.assertRaisesRegex(TypeError, "generic ROHF/SCF"):
            self.vibronic.mrsf_analytic_property_derivative_from_provider(
                ReferenceOnlyProvider(),
                property_kind="state_dipole",
                electronic_state="S1",
            )

    def test_ir_and_raman_records_are_json_compatible(self):
        ir = self.vibronic.excited_state_ir_intensities(
            self.property([[1.0], [0.0], [0.0]], "state_dipole", "e*bohr")
        )
        polarizability = np.zeros((3, 3, 1))
        polarizability[:, :, 0] = np.eye(3)
        raman = self.vibronic.excited_state_raman_activities(
            self.property(polarizability, "polarizability", "bohr^3")
        )
        ir_record = ir.to_dict()
        raman_record = raman.to_dict()
        self.assertEqual(ir_record["mode_dipole_derivative_unit"], "e/sqrt(amu)")
        self.assertEqual(raman_record["activity_unit"], "bohr^4/amu")
        json.dumps(ir_record)
        json.dumps(raman_record)


class TestResonanceRamanAndPublicAPI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()
        cls.model = cls.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [1000.0],
            [[1.0]],
            [0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="RR synthetic normal modes",
        )
        cls.transition = cls.vibronic.ElectronicTransitionMoment.create(
            [1.0, 0.0, 0.0],
            electronic_state="S1",
            electronic_phase_convention="RR synthetic electronic phase",
        )

    def test_constant_transition_dipole_gives_zero_fundamental_for_identical_surfaces(self):
        result = self.vibronic.resonance_raman_fc_ht(
            self.model,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            incident_wavenumber_cm1=19500.0,
            damping_cm1=100.0,
            transition=self.transition,
            transition_dipole_derivative=None,
            max_intermediate_quanta=3,
        )
        self.assertLess(result.tensor_norm_squared[0], 1.0e-24)

    def test_first_order_ht_gives_nonzero_resonance_raman_fundamental(self):
        derivative = self.vibronic.ExcitedStatePropertyDerivative.create(
            [[0.02], [0.0], [0.0]],
            property_kind="transition_dipole",
            coordinate_basis="ground_normal",
            coordinate_unit="sqrt(me)*bohr",
            property_unit="e*bohr",
            electronic_state="S1",
            electronic_phase_convention="RR synthetic electronic phase",
            coordinate_phase_convention="RR synthetic normal modes",
            provenance="synthetic derivative",
            electronic_state_role="target_excited_state",
        )
        result = self.vibronic.resonance_raman_fc_ht(
            self.model,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            incident_wavenumber_cm1=19500.0,
            damping_cm1=100.0,
            transition=self.transition,
            transition_dipole_derivative=derivative,
            max_intermediate_quanta=3,
        )
        self.assertGreater(result.tensor_norm_squared[0], 0.0)
        self.assertIn("resonant-term", result.approximation)
        self.assertGreater(result.intermediate_transition_completeness, 0.999)
        self.assertEqual(result.to_dict()["tensor_unit"], "(e*bohr)^2/cm^-1")

    def test_incomplete_resonance_raman_sum_is_rejected(self):
        displaced = self.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [900.0],
            [[1.0]],
            [40.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="RR synthetic normal modes",
        )
        with self.assertRaisesRegex(ValueError, "transition-strength completeness"):
            self.vibronic.resonance_raman_fc_ht(
                displaced,
                electronic_origin_cm1=20000.0,
                origin_kind="zero_zero",
                incident_wavenumber_cm1=19500.0,
                damping_cm1=100.0,
                transition=self.transition,
                transition_dipole_derivative=None,
                max_intermediate_quanta=0,
            )

        partial = self.vibronic.resonance_raman_fc_ht(
            displaced,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            incident_wavenumber_cm1=19500.0,
            damping_cm1=100.0,
            transition=self.transition,
            transition_dipole_derivative=None,
            max_intermediate_quanta=0,
            minimum_intermediate_completeness=0.0,
        )
        self.assertLess(partial.intermediate_transition_completeness, 0.1)

    def test_json_cli_writes_lines_and_normalized_spectrum(self):
        data = {
            "ground_frequencies_cm1": [1000.0],
            "excited_frequencies_cm1": [1000.0],
            "J": [[1.0]],
            "K": [3.0],
            "coordinate_unit": "sqrt(me)*bohr",
            "coordinate_phase_convention": "CLI synthetic phase",
            "electronic_origin_cm1": 20000.0,
            "origin_kind": "zero_zero",
            "max_final_quanta": 4,
            "fwhm_cm1": 50.0,
        }
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input.json"
            target = Path(directory) / "output.json"
            source.write_text(json.dumps(data))
            self.assertEqual(self.vibronic.main([str(source), str(target)]), 0)
            output = json.loads(target.read_text())
        self.assertEqual(output["model"], "multidimensional harmonic Franck-Condon")
        self.assertEqual(len(output["lines"]), 5)
        self.assertEqual(output["normalization"], "sum")

    def test_public_library_import_exports_vibronic_api(self):
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
                HarmonicVibronicModel,
                harmonic_vibronic_spectrum,
                excited_state_ir_intensities,
                excited_state_raman_activities,
            )

            assert HarmonicVibronicModel
            assert callable(harmonic_vibronic_spectrum)
            assert callable(excited_state_ir_intensities)
            assert callable(excited_state_raman_activities)
            """
        )
        subprocess.run([sys.executable, "-c", script], check=True)


if __name__ == "__main__":
    unittest.main()
