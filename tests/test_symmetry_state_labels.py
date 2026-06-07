import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relpath):
    spec = importlib.util.spec_from_file_location(name, str(ROOT / relpath))
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


WATER_CHARGES = [8, 1, 1]
WATER_COORDS = np.array(
    [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ]
)
WATER_SHELLS = [(0, 0), (0, 0), (0, 1), (1, 0), (2, 0)]
WATER_NAO = 7


def water_salcs():
    symmetry = load_module(
        "openqp_symmetry_state_salc_helper",
        "pyoqp/oqp/library/symmetry.py",
    )
    detect = load_module(
        "openqp_symmetry_detect_state_helper",
        "pyoqp/oqp/library/symmetry_detect.py",
    )
    detection = detect.detect_point_group(WATER_CHARGES, WATER_COORDS, tolerance=1.0e-6)
    u, labels = symmetry.build_symmetry_adapted_transform(
        WATER_SHELLS, detection["operations"], detection["character_table"]
    )
    return symmetry, detection, u, labels


class TestAssignStateIrreps(unittest.TestCase):
    def setUp(self):
        self.symmetry, self.detection, self.u, self.labels = water_salcs()

    def single_amplitude(self, n_occ, n_vir, i, a):
        x = np.zeros((1, n_occ, n_vir))
        x[0, i, a] = 1.0
        return x

    def test_single_excitation_label_is_pair_product(self):
        n_occ, n_vir = 5, 2
        occ, vir = self.u[:, :n_occ], self.u[:, n_occ:]
        for i in range(n_occ):
            for a in range(n_vir):
                result = self.symmetry.assign_state_irreps(
                    self.single_amplitude(n_occ, n_vir, i, a),
                    occ, vir, np.eye(WATER_NAO), WATER_SHELLS,
                    self.detection["operations"],
                    self.detection["character_table"],
                    matrix_key="matrix_input_frame",
                )
                expected = self.symmetry.product_irrep(
                    [self.labels[i], self.labels[n_occ + a]],
                    self.detection["character_table"],
                )
                self.assertEqual(result["labels"][0], expected)
                self.assertLess(result["max_deviation"], 1.0e-10)

    def test_mixed_symmetry_amplitude_is_flagged(self):
        n_occ, n_vir = 5, 2
        occ, vir = self.u[:, :n_occ], self.u[:, n_occ:]
        # Combine two single excitations with different pair products.
        def product(i, a):
            return self.symmetry.product_irrep(
                [self.labels[i], self.labels[n_occ + a]],
                self.detection["character_table"],
            )

        pairs = [(i, a) for i in range(n_occ) for a in range(n_vir)]
        first = pairs[0]
        second = next(p for p in pairs if product(*p) != product(*first))
        x = np.zeros((1, n_occ, n_vir))
        x[0, first[0], first[1]] = 1.0 / np.sqrt(2.0)
        x[0, second[0], second[1]] = 1.0 / np.sqrt(2.0)
        result = self.symmetry.assign_state_irreps(
            x, occ, vir, np.eye(WATER_NAO), WATER_SHELLS,
            self.detection["operations"], self.detection["character_table"],
            matrix_key="matrix_input_frame",
        )
        self.assertEqual(result["labels"][0], "mixed")

    def test_reference_labels_multiply_in(self):
        n_occ, n_vir = 5, 2
        occ, vir = self.u[:, :n_occ], self.u[:, n_occ:]
        result = self.symmetry.assign_state_irreps(
            self.single_amplitude(n_occ, n_vir, 0, 0),
            occ, vir, np.eye(WATER_NAO), WATER_SHELLS,
            self.detection["operations"], self.detection["character_table"],
            reference_labels=["b1", "b2"],
            matrix_key="matrix_input_frame",
        )
        transition = result["transition_labels"][0]
        expected = self.symmetry.product_irrep(
            [transition, "b1", "b2"], self.detection["character_table"]
        )
        self.assertEqual(result["labels"][0], expected)


class TestMoleculeStateLabelWiring(unittest.TestCase):
    def _build_molecule(self, td_type, na, nb):
        symmetry, _detection, u, labels = water_salcs()
        detect = sys.modules["openqp_symmetry_detect_state_helper"]
        metadata_tests = load_module(
            "openqp_symmetry_metadata_tests_for_states",
            "tests/test_symmetry_metadata.py",
        )
        molecule_mod = metadata_tests.load_molecule_module()

        library_pkg = types.ModuleType("oqp.library")
        library_pkg.__path__ = [str(ROOT / "pyoqp/oqp/library")]
        library_pkg.symmetry_detect = detect
        library_pkg.symmetry = symmetry
        sys.modules["oqp.library"] = library_pkg
        sys.modules["oqp.library.symmetry_detect"] = detect
        sys.modules["oqp.library.symmetry"] = symmetry

        mol = object.__new__(molecule_mod.Molecule)
        mol.config = {
            "symmetry": {"enabled": "true"},
            "scf": {"type": "rhf" if td_type in ("tda", "rpa") else "rohf"},
            "tdhf": {"type": td_type, "mult": 1},
        }
        mol.symmetry_metadata = {}
        mol.get_atoms = lambda: np.array(WATER_CHARGES, dtype=float)
        mol.get_system = lambda: WATER_COORDS.ravel().copy()
        mol.initialize_symmetry_metadata()

        vec_mo = u.T.copy()
        packed = np.eye(WATER_NAO)[np.tril_indices(WATER_NAO)]

        class _Data(dict):
            def get_basis(self):
                return {
                    "centers": np.array([at for at, _ in WATER_SHELLS]),
                    "angs": np.array([l for _, l in WATER_SHELLS]),
                    "nbf": WATER_NAO,
                }

        mol.data = _Data()
        mol.data["OQP::SM"] = packed
        mol.data["OQP::VEC_MO_A"] = vec_mo.ravel()
        mol.data["OQP::VEC_MO_B"] = vec_mo.ravel()
        mol.data["nelec_A"] = na
        mol.data["nelec_B"] = nb
        mol.data["nocc"] = max(na, nb)
        return mol, labels

    @staticmethod
    def pack_states(amplitudes):
        """(n_states, n_occ, n_vir) -> Fortran-layout flat buffer."""
        return np.concatenate(
            [x.T.ravel() for x in amplitudes]
        )

    def test_tda_single_excitations_labeled(self):
        mol, labels = self._build_molecule("tda", 5, 5)
        symmetry = sys.modules["oqp.library.symmetry"]
        table = mol.symmetry_metadata["detection"]["character_table"]
        n_occ, n_vir = 5, 2

        x1 = np.zeros((n_occ, n_vir))
        x1[4, 0] = 1.0  # HOMO -> LUMO
        x2 = np.zeros((n_occ, n_vir))
        x2[3, 0] = 1.0
        mol.data["OQP::td_bvec_mo"] = self.pack_states([x1, x2])

        result = mol.label_excited_states()
        self.assertEqual(result["status"], "ok")
        expected1 = symmetry.product_irrep([labels[4], labels[5]], table)
        expected2 = symmetry.product_irrep([labels[3], labels[5]], table)
        self.assertEqual(result["labels"], [expected1, expected2])
        self.assertEqual(result["terms"][0], f"1{expected1.upper()}")
        self.assertEqual(
            mol.symmetry_metadata["state_labels"]["labels"], result["labels"]
        )

    def test_mrsf_includes_somo_reference_product(self):
        # Triplet reference: na=5, nb=3 -> SOMOs are alpha MOs 3 and 4.
        mol, labels = self._build_molecule("mrsf", 5, 3)
        symmetry = sys.modules["oqp.library.symmetry"]
        table = mol.symmetry_metadata["detection"]["character_table"]
        n_occ, n_vir = 5, 4  # occ alpha x vir beta

        x = np.zeros((n_occ, n_vir))
        x[4, 1] = 1.0  # alpha MO 4 -> beta MO 4 (index nb+1)
        mol.data["OQP::td_bvec_mo"] = self.pack_states([x])

        result = mol.label_excited_states()
        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["reference_labels"], labels[3:5])
        transition = symmetry.product_irrep([labels[4], labels[3 + 1]], table)
        expected = symmetry.product_irrep(
            [transition] + labels[3:5], table
        )
        self.assertEqual(result["labels"][0], expected)

    def test_unsupported_td_type_skips(self):
        mol, _ = self._build_molecule("tda", 5, 5)
        mol.config["tdhf"]["type"] = "exotic"
        self.assertIsNone(mol.label_excited_states())
        self.assertEqual(
            mol.symmetry_metadata["state_labels"]["status"],
            "skipped_unsupported_td_type_exotic",
        )

    def test_state_labels_written_to_log(self):
        import tempfile

        mol, labels = self._build_molecule("tda", 5, 5)
        n_occ, n_vir = 5, 2
        x = np.zeros((n_occ, n_vir))
        x[4, 0] = 1.0
        mol.data["OQP::td_bvec_mo"] = self.pack_states([x])
        mol.data["OQP::td_energies"] = np.array([0.31])

        with tempfile.TemporaryDirectory() as tmp:
            mol.log = str(Path(tmp) / "run.log")
            mol.label_excited_states()
            with open(mol.log, encoding="utf-8") as fin:
                content = fin.read()
        self.assertIn("excited-state symmetry labels", content)
        self.assertIn("state   1", content)
        self.assertIn("0.31000000", content)


if __name__ == "__main__":
    unittest.main()
