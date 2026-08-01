"""Pure-Python regressions for gauge-resolved frozen NAC comparisons."""

from __future__ import annotations

import configparser
from pathlib import Path
import tempfile
import unittest

import numpy as np

from tools.nac_lagrangian.nac_reference_gate import (
    NACGateError,
    compare_files,
    compare_payloads,
)


ROOT = Path(__file__).resolve().parents[1]
GATE_DIR = ROOT / "tools" / "nac_lagrangian"


def antisymmetric_dcv(natom: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    dcv = np.zeros((3, 3, natom, 3))
    for istate in range(3):
        for jstate in range(istate + 1, 3):
            block = rng.standard_normal((natom, 3))
            dcv[istate, jstate] = block
            dcv[jstate, istate] = -block
    return dcv


def apply_state_gauge(dcv: np.ndarray, signs: tuple[int, ...]) -> np.ndarray:
    gauged = np.array(dcv, copy=True)
    for istate in range(dcv.shape[0]):
        for jstate in range(dcv.shape[1]):
            gauged[istate, jstate] *= signs[istate] * signs[jstate]
    return gauged


def nacv_from(dcv: np.ndarray, excited_energies: np.ndarray) -> np.ndarray:
    nacv = np.zeros_like(dcv)
    for istate in range(dcv.shape[0]):
        for jstate in range(dcv.shape[1]):
            nacv[istate, jstate] = (
                excited_energies[jstate] - excited_energies[istate]
            ) * dcv[istate, jstate]
    return nacv


class NACReferenceGateTests(unittest.TestCase):
    def test_h2o_convergence_resolves_state_gauge(self):
        reference = antisymmetric_dcv(natom=3, seed=1)
        signs = (1, -1, 1)
        candidate = apply_state_gauge(reference, signs)
        noise = np.zeros_like(candidate)
        noise[0, 1, 0, 0] = 2.0e-7
        noise[1, 0, 0, 0] = -2.0e-7
        candidate += noise
        result = compare_payloads(
            {"dcv": reference},
            {"dcv": candidate},
            component_atol=3.0e-7,
            label="H2O",
        )
        self.assertEqual(result.state_signs, signs)
        self.assertAlmostEqual(result.max_component_error, 2.0e-7)

    def test_ethylene_rejects_independent_pair_flip(self):
        reference = antisymmetric_dcv(natom=6, seed=2)
        candidate = apply_state_gauge(reference, (1, -1, 1))
        candidate[1, 2] *= -1.0
        candidate[2, 1] *= -1.0
        with self.assertRaisesRegex(NACGateError, "not a state gauge"):
            compare_payloads(
                {"dcv": reference},
                {"dcv": candidate},
                component_atol=1.0e-12,
                label="ethylene",
            )

    def test_acrolein_enforces_component_tolerance(self):
        reference = antisymmetric_dcv(natom=8, seed=3)
        candidate = np.array(reference, copy=True)
        candidate[0, 2, 7, 2] += 2.0e-4
        candidate[2, 0, 7, 2] -= 2.0e-4
        with self.assertRaisesRegex(NACGateError, "max component error"):
            compare_payloads(
                {"dcv": reference},
                {"dcv": candidate},
                component_atol=1.0e-4,
                label="Acrolein",
            )

    def test_npz_flags_energies_and_gap_identity(self):
        dcv = antisymmetric_dcv(natom=3, seed=4)
        excited = np.array([-76.3, -76.0, -75.9])
        payload = {
            "dcv": dcv,
            "nacv": nacv_from(dcv, excited),
            "energies": np.concatenate(([-76.6], excited)),
            "flags": np.array(["computed"] * 18),
        }
        with tempfile.TemporaryDirectory() as tmp:
            reference = Path(tmp) / "reference.npz"
            candidate = Path(tmp) / "candidate.npz"
            np.savez(reference, **payload)
            np.savez(candidate, **payload)
            result = compare_files(
                reference,
                candidate,
                component_atol=0.0,
                energy_atol=0.0,
                require_flags=True,
            )
        self.assertEqual(result.max_component_error, 0.0)
        self.assertEqual(result.max_energy_error, 0.0)

    def test_noncomputed_worker_flag_is_rejected(self):
        dcv = antisymmetric_dcv(natom=3, seed=5)
        with self.assertRaisesRegex(NACGateError, "unsuccessful flags"):
            compare_payloads(
                {"dcv": dcv, "flags": np.array(["computed"])},
                {"dcv": dcv, "flags": np.array(["computed", "failed"])},
                component_atol=0.0,
            )

    def test_restart_and_analytic_success_flags_are_accepted(self):
        dcv = antisymmetric_dcv(natom=3, seed=6)
        compare_payloads(
            {"dcv": dcv, "flags": np.array([b"loaded"])},
            {"dcv": dcv, "flags": np.array(["analytic-v3-zvector"])},
            component_atol=0.0,
            require_flags=True,
        )

    def test_energy_tolerance_requires_energies_in_both_artifacts(self):
        dcv = antisymmetric_dcv(natom=3, seed=7)
        with self.assertRaisesRegex(NACGateError, "requires energies"):
            compare_payloads(
                {"dcv": dcv, "energies": np.arange(3.0)},
                {"dcv": dcv},
                component_atol=0.0,
                energy_atol=1.0e-8,
            )

    def test_acrolein_inputs_preserve_gamess_bohr_geometry(self):
        original_bohr = np.array(
            [
                [-6.9913814, -2.0220633, -0.0574336],
                [-4.8213949, -0.8772009, -0.0045365],
                [-2.2317045, -1.8893978, 0.0785222],
                [-8.7466350, -0.9479836, 0.0817211],
                [-7.3824334, -3.9017900, 0.3044129],
                [-4.6757554, 1.1527992, 0.1730773],
                [-0.3289985, -0.9804003, -0.0309508],
                [-2.2770514, -4.3034240, -0.2653938],
            ]
        )
        bohr_radius_angstrom = 0.529177210903
        for filename, expected_type, expected_dx in (
            ("Acrolein_S2_tlf2_dx5e4.inp", "numerical", 5.0e-4),
            ("Acrolein_S2_tlf2_dx25e5.inp", "numerical", 2.5e-4),
            ("Acrolein_S2_tlf2_analytic.inp", "analytical", None),
        ):
            with self.subTest(filename=filename):
                parser = configparser.ConfigParser(interpolation=None)
                loaded = parser.read(GATE_DIR / filename)
                self.assertEqual(len(loaded), 1)
                rows = parser["input"]["system"].strip().splitlines()
                parsed = np.array(
                    [[float(value) for value in row.split()[1:4]] for row in rows]
                )
                np.testing.assert_allclose(
                    parsed / bohr_radius_angstrom,
                    original_bohr,
                    atol=1.0e-9,
                )
                self.assertEqual(parser["input"]["runtype"], "energy")
                self.assertEqual(parser["input"]["functional"], "bhhlyp")
                self.assertEqual(parser["tdhf"].getint("target"), 3)
                self.assertEqual(parser["tdhf"].getint("tlf"), 2)
                self.assertLessEqual(parser["scf"].getfloat("conv"), 1.0e-10)
                self.assertLessEqual(parser["tdhf"].getfloat("conv"), 1.0e-10)
                self.assertEqual(parser["nac"]["type"], expected_type)
                if expected_dx is None:
                    self.assertNotIn("dx", parser["nac"])
                else:
                    self.assertAlmostEqual(
                        parser["nac"].getfloat("dx"), expected_dx
                    )


if __name__ == "__main__":
    unittest.main()
