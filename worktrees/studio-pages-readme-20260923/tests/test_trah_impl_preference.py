"""Guards for the default TRAH implementation preference."""

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
TYPES = ROOT / "source" / "types.F90"
SCF = ROOT / "source" / "scf.F90"


class TRAHImplementationPreferenceTests(unittest.TestCase):
    def test_python_auto_prefers_native(self):
        src = OQPDATA.read_text()

        self.assertIn('"auto": 1', src)
        self.assertIn('"native": 1', src)
        self.assertIn('"otr": 0', src)
        self.assertIn("molecule.control.trh_impl = 1", src)

        auto_block = re.search(
            r"if trh_choice == 'auto':(?P<body>.*?)(?=\n\s*natom =)",
            src,
            re.S,
        )
        self.assertIsNotNone(auto_block)
        self.assertNotIn("runtype", auto_block.group("body"))
        self.assertNotIn("method", auto_block.group("body"))

    def test_fortran_default_prefers_native(self):
        src = TYPES.read_text()
        self.assertRegex(src, r"trh_impl\s*=\s*1\s*!< TRAH solver: 1=native")

    def test_otr_remains_explicit_second_choice(self):
        src = SCF.read_text()
        self.assertIn("infos%control%trh_impl == 1", src)
        self.assertIn("call trah_native_run", src)
        self.assertIn("explicit trh_impl=otr", src)
        self.assertIn("call init_trah_solver", src)
        self.assertIn("call run_trah_solver", src)


if __name__ == "__main__":
    unittest.main()
