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
        # Only the default value is pinned, not the trailing comment: trh_impl no
        # longer selects an implementation, so its wording is free to change.
        self.assertRegex(src, r"trh_impl\s*=\s*1\b")

    def test_canonicalization_is_not_gated_on_trh_impl(self):
        """Post-TRAH canonicalization must not depend on trh_impl.

        run_trah always runs the native solver now, so gating the Fock
        diagonalization on control%trh_impl == 1 meant a non-Python caller
        using include/oqp.h with trh_impl != 1 still entered TRAH but skipped
        canonicalization, handing rotated orbitals and orbital energies to
        post-SCF properties and analytic gradients.
        """
        src = SCF.read_text()
        self.assertNotIn("control%trh_impl", src)
        self.assertIn("if (scf_type /= scf_rohf) then", src)



if __name__ == "__main__":
    unittest.main()
