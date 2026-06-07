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


def random_rotation(seed):
    rng = np.random.RandomState(seed)
    q, _ = np.linalg.qr(rng.randn(3, 3))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


WATER_CHARGES = [8, 1, 1]
WATER_COORDS = np.array(
    [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ]
)
# O 1s, 2s, 2p; H 1s each -> 7 AOs.
WATER_SHELLS = [(0, 0), (0, 0), (0, 1), (1, 0), (2, 0)]
WATER_NAO = 7


class TestShellBlockRepresentation(unittest.TestCase):
    def setUp(self):
        self.symmetry = load_module(
            "openqp_symmetry_shell_block_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )

    def test_p_through_g_blocks_form_a_representation(self):
        m1 = random_rotation(1)
        m2 = random_rotation(2)
        for l in (1, 2, 3, 4):
            b12 = self.symmetry._shell_block(l, m1 @ m2)
            b1 = self.symmetry._shell_block(l, m1)
            b2 = self.symmetry._shell_block(l, m2)
            self.assertLess(float(np.max(np.abs(b12 - b1 @ b2))), 1.0e-11)

    def test_f_and_g_blocks_reduce_to_component_signs(self):
        for l in (3, 4):
            for signs in [(-1, -1, 1), (1, -1, -1), (-1, -1, -1)]:
                block = self.symmetry._shell_block(l, np.diag(signs).astype(float))
                expected = np.diag(self.symmetry._component_signs(l, signs)).astype(float)
                self.assertTrue(np.allclose(block, expected), f"l={l}, signs={signs}")

    def test_f_and_g_blocks_preserve_shell_overlap_metric(self):
        # <x^a y^b z^c | x^a' y^b' z^c'> over a Gaussian is the product of
        # 1D moments (p-1)!! for even p (zero for odd), normalized.
        def metric(l):
            monos = self.symmetry._CART_MONOMIALS[l]
            norms = self.symmetry._monomial_norms(l)
            df = self.symmetry._double_factorial
            size = len(monos)
            s = np.zeros((size, size))
            for i, (a, b, c) in enumerate(monos):
                for j, (d, e, f) in enumerate(monos):
                    if (a + d) % 2 or (b + e) % 2 or (c + f) % 2:
                        continue
                    s[i, j] = (
                        df(a + d - 1) * df(b + e - 1) * df(c + f - 1)
                    ) * norms[i] * norms[j]
            return s

        for l in (2, 3, 4):
            m = random_rotation(5 + l)
            block = self.symmetry._shell_block(l, m)
            s = metric(l)
            deviation = np.max(np.abs(block.T @ s @ block - s))
            self.assertLess(float(deviation), 1.0e-11, f"l={l}")

    def test_d_block_preserves_shell_overlap_metric(self):
        # Normalized Cartesian d functions are not mutually orthogonal
        # (<xx|yy> = 1/3), so rotations preserve this metric, not the identity.
        m = random_rotation(3)
        block = self.symmetry._shell_block(2, m)
        metric = np.eye(6)
        metric[:3, :3] = (np.eye(3) * 2.0 + np.ones((3, 3))) / 3.0
        deviation = np.max(np.abs(block.T @ metric @ block - metric))
        self.assertLess(float(deviation), 1.0e-12)

    def test_d_block_reduces_to_component_signs_for_sign_ops(self):
        signs = (-1, 1, -1)  # C2y
        block = self.symmetry._shell_block(2, np.diag(signs).astype(float))
        expected = np.diag(self.symmetry._component_signs(2, signs)).astype(float)
        self.assertTrue(np.allclose(block, expected))


class TestFShellLabeling(unittest.TestCase):
    """End-to-end SALC + MO labeling with an f shell present."""

    def test_water_with_f_shell_labels_cleanly(self):
        symmetry = load_module(
            "openqp_symmetry_fshell_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        detect = load_module(
            "openqp_symmetry_detect_fshell_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        detection = detect.detect_point_group(
            WATER_CHARGES, WATER_COORDS, tolerance=1.0e-6
        )
        shells = [(0, 0), (0, 3), (1, 0), (2, 0)]  # O: s+f, H: s each
        n_ao = 1 + 10 + 1 + 1

        u, labels = symmetry.build_symmetry_adapted_transform(
            shells, detection["operations"], detection["character_table"]
        )
        self.assertEqual(u.shape, (n_ao, n_ao))
        deviation = np.max(np.abs(u.T @ u - np.eye(n_ao)))
        self.assertLess(float(deviation), 1.0e-12)
        self.assertTrue(set(labels) <= {"a1", "a2", "b1", "b2"})

        result = symmetry.assign_mo_irreps(
            u, np.eye(n_ao), shells,
            detection["operations"], detection["character_table"],
            matrix_key="matrix_input_frame",
        )
        self.assertEqual(result["labels"], labels)
        self.assertLess(result["max_deviation"], 1.0e-10)


class TestRotatedFrameLabeling(unittest.TestCase):
    """MO labels must be frame-independent via matrix_input_frame."""

    def setUp(self):
        self.symmetry = load_module(
            "openqp_symmetry_rotated_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        self.detect = load_module(
            "openqp_symmetry_detect_rotated_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )

    def test_water_labels_invariant_under_input_rotation(self):
        rotation = random_rotation(42)
        rotated_coords = WATER_COORDS @ rotation.T

        reference = self.detect.detect_point_group(
            WATER_CHARGES, WATER_COORDS, tolerance=1.0e-6
        )
        rotated = self.detect.detect_point_group(
            WATER_CHARGES, rotated_coords, tolerance=1.0e-6
        )
        self.assertEqual(rotated["point_group"], "c2v")
        self.assertEqual(rotated["abelian_subgroup"], "c2v")

        # SALCs built in each standard frame, expressed in the local AO basis.
        u_ref, labels_ref = self.symmetry.build_symmetry_adapted_transform(
            WATER_SHELLS, reference["operations"], reference["character_table"]
        )

        # Express the reference standard-frame SALCs in the rotated input
        # frame: p coefficients transform with W = R_ref^T then R_rot maps
        # input axes; combined per-shell block W.
        r_ref = np.asarray(reference["orientation"])
        r_rot = np.asarray(rotated["orientation"])
        # AO function defined along reference-standard axes, re-expressed in
        # the rotated input frame: rows = rotated-input components.
        w3 = rotation @ r_ref.T
        w = np.eye(WATER_NAO)
        w[2:5, 2:5] = self.symmetry._shell_block(1, w3)
        mos_in_rotated_frame = w @ u_ref

        result = self.symmetry.assign_mo_irreps(
            mos_in_rotated_frame,
            np.eye(WATER_NAO),
            WATER_SHELLS,
            rotated["operations"],
            rotated["character_table"],
            matrix_key="matrix_input_frame",
        )
        self.assertLess(result["max_deviation"], 1.0e-8)
        # a1/a2 labels are orientation-independent; b1/b2 may swap with the
        # in-plane axis convention, so compare modulo that swap.
        def normalize(seq):
            return ["b" if lbl in ("b1", "b2") else lbl for lbl in seq]

        self.assertEqual(normalize(result["labels"]), normalize(labels_ref))
        from collections import Counter

        self.assertEqual(Counter(result["labels"]), Counter(labels_ref))


class TestMoleculeMoLabelWiring(unittest.TestCase):
    def _build_molecule(self, scf_type="rhf"):
        detect = load_module(
            "oqp_symmetry_detect_molabel_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        symmetry = load_module(
            "oqp_symmetry_molabel_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        metadata_tests = load_module(
            "openqp_symmetry_metadata_tests_for_molabel",
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
        mol.config = {"symmetry": {"enabled": "true"}, "scf": {"type": scf_type}}
        mol.symmetry_metadata = {}
        mol.get_atoms = lambda: np.array(WATER_CHARGES, dtype=float)
        mol.get_system = lambda: WATER_COORDS.ravel().copy()
        mol.initialize_symmetry_metadata()

        detection = mol.symmetry_metadata["detection"]
        u, labels = symmetry.build_symmetry_adapted_transform(
            WATER_SHELLS, detection["operations"], detection["character_table"]
        )
        # MOs in the input frame == standard frame here (no rotation applied),
        # Fortran stores MOs column-wise => flat buffer rows are MOs.
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
        return mol, labels

    def test_label_molecular_orbitals_rhf(self):
        mol, expected = self._build_molecule("rhf")
        result = mol.label_molecular_orbitals()
        self.assertIsNotNone(result)
        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["alpha"]["labels"], expected)
        self.assertNotIn("beta", result)
        self.assertEqual(
            mol.symmetry_metadata["mo_labels"]["alpha"]["labels"], expected
        )

    def test_label_molecular_orbitals_uhf_labels_both_spins(self):
        mol, expected = self._build_molecule("uhf")
        result = mol.label_molecular_orbitals()
        self.assertEqual(result["alpha"]["labels"], expected)
        self.assertEqual(result["beta"]["labels"], expected)

    def test_labeling_disabled_when_symmetry_off(self):
        mol, _ = self._build_molecule("rhf")
        mol.symmetry_metadata = {}
        self.assertIsNone(mol.label_molecular_orbitals())

    def test_label_mo_false_skips(self):
        mol, _ = self._build_molecule("rhf")
        mol.symmetry_metadata["label_mo"] = False
        self.assertIsNone(mol.label_molecular_orbitals())
        self.assertNotIn("mo_labels", mol.symmetry_metadata)

    def test_product_irrep_c2v_table(self):
        symmetry = load_module(
            "openqp_symmetry_product_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        table = {
            "a1": [1, 1, 1, 1],
            "a2": [1, 1, -1, -1],
            "b1": [1, -1, 1, -1],
            "b2": [1, -1, -1, 1],
        }
        self.assertEqual(symmetry.product_irrep([], table), "a1")
        self.assertEqual(symmetry.product_irrep(["b1", "b2"], table), "a2")
        self.assertEqual(symmetry.product_irrep(["b1", "b1"], table), "a1")
        self.assertEqual(symmetry.product_irrep(["a1", "b2"], table), "b2")
        self.assertEqual(symmetry.product_irrep(["a1", "mixed"], table), "mixed")

    def test_closed_shell_state_is_totally_symmetric(self):
        mol, labels = self._build_molecule("rhf")
        mol.data["nelec_A"] = 5
        mol.data["nelec_B"] = 5
        mol.mult = 1
        result = mol.label_molecular_orbitals()
        self.assertEqual(result["scf_state"]["irrep"], "a1")
        self.assertEqual(result["scf_state"]["term"], "1A1")

    def test_open_shell_state_carries_somo_irrep(self):
        mol, labels = self._build_molecule("rohf")
        mol.data["nelec_A"] = 5
        mol.data["nelec_B"] = 4
        mol.mult = 2
        result = mol.label_molecular_orbitals()
        somo = result["alpha"]["labels"][4]
        self.assertEqual(result["scf_state"]["irrep"], somo)
        self.assertEqual(result["scf_state"]["term"], f"2{somo.upper()}")

    def test_state_label_in_log(self):
        import tempfile

        mol, _ = self._build_molecule("rhf")
        mol.data["nelec_A"] = 5
        mol.data["nelec_B"] = 5
        mol.mult = 1
        with tempfile.TemporaryDirectory() as tmp:
            mol.log = str(Path(tmp) / "run.log")
            mol.label_molecular_orbitals()
            with open(mol.log, encoding="utf-8") as fin:
                content = fin.read()
        self.assertIn("SCF state:        1A1", content)

    def test_mo_labels_written_to_log(self):
        import tempfile

        mol, expected = self._build_molecule("rhf")
        with tempfile.TemporaryDirectory() as tmp:
            mol.log = str(Path(tmp) / "run.log")
            mol.data["OQP::E_MO_A"] = np.linspace(-20.0, 2.0, WATER_NAO)
            mol.label_molecular_orbitals()
            with open(mol.log, encoding="utf-8") as fin:
                content = fin.read()
        self.assertIn("MO symmetry labels", content)
        self.assertIn("point group:      c2v", content)
        self.assertIn("alpha MOs:", content)
        for label in set(expected):
            self.assertIn(label, content)

    def test_h_shells_are_skipped_gracefully(self):
        mol, _ = self._build_molecule("rhf")

        class _HData(dict):
            def get_basis(self):
                return {
                    "centers": np.array([0]),
                    "angs": np.array([5]),
                    "nbf": 21,
                }

        data = _HData()
        data["OQP::SM"] = np.zeros(21 * 22 // 2)
        mol.data = data
        mol.label_molecular_orbitals()
        self.assertEqual(
            mol.symmetry_metadata["mo_labels"]["status"],
            "skipped_unsupported_shells_beyond_g",
        )


if __name__ == "__main__":
    unittest.main()
