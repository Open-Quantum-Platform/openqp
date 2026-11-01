"""Custom basis paths must never become Molden filename components."""
import ast
import hashlib
import ntpath
from pathlib import Path
from types import SimpleNamespace
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "pyoqp/oqp/library/single_point.py"

def pack_name(basis, log="job.log", idx=1):
    tree = ast.parse(SOURCE.read_text())
    method = next(node for cls in tree.body if isinstance(cls, ast.ClassDef)
                  for node in cls.body if isinstance(node, ast.FunctionDef)
                  and node.name == "pack_molden_name")
    scope = {"hashlib": hashlib}
    exec(compile(ast.Module(body=[method], type_ignores=[]), str(SOURCE), "exec"), scope)
    obj = SimpleNamespace(basis=basis, mol=SimpleNamespace(log=log, idx=idx))
    return scope["pack_molden_name"](obj, "scf", "rhf", "hf")

class PortableMoldenNamesTests(unittest.TestCase):
    def test_standard_basis_filename_is_unchanged(self):
        self.assertEqual(pack_name("6-31g*"), "job_scf_rhf_hf_6-31gs.molden")

    def test_custom_path_does_not_escape_output_directory(self):
        for basis in ("file:../basis.json", r"file:C:\basis\custom.json", "file:/tmp/basis.json"):
            with self.subTest(basis=basis):
                name = pack_name(basis, log="/output/job.log")
                self.assertEqual(Path(name).parent, Path("/output"))
                self.assertTrue(Path(name).name.startswith("job_scf_rhf_hf_custom-"))
                self.assertFalse(any(c in Path(name).name for c in '<>:"/\\|?*'))

    def test_custom_labels_are_stable_and_distinct(self):
        self.assertEqual(pack_name("file:a.json"), pack_name("file:a.json"))
        self.assertNotEqual(pack_name("file:a.json"), pack_name("file:b.json"))

    def test_failed_windows_dkh_path_stays_below_max_path(self):
        directory = r"C:\GitLab-Runner\builds\Nc2j9q4px\1\open-quantum-platform\internal\openqp\openqp_test_tmp_2026-10-03_17-28-00\HBr_RHF_DKH2_UNCONTRACTED_LEGACY_ENERGY__86bf0d3cbb7e"
        log = directory + r"\HBr_RHF_DKH2_UNCONTRACTED_LEGACY_ENERGY.log"
        name = pack_name("file:x2c-tzvpall_uncontracted_hbr.json", log)
        self.assertEqual(ntpath.dirname(name), directory)
        self.assertLess(len(name), 260)
        self.assertNotIn(":", ntpath.basename(name))

    def test_indexed_frames_keep_their_suffix(self):
        self.assertIn("_2_scf_rhf_hf_custom-", pack_name("file:a.json", idx=2))

if __name__ == "__main__":
    unittest.main()
