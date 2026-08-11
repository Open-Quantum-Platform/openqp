"""The XOR encoding staged for correlated methods must be the real product.

``Molecule.stage_mo_irreps`` hands Fortran one integer per irrep -- the
bitmask of operations whose character is -1 -- so that a consumer can get the
irrep of a determinant, or of an amplitude block, with one ``ieor`` per
occupied orbital instead of a table lookup per pair.

That shortcut is only valid if the encoding is faithful (distinct irreps get
distinct codes) and if XOR really is the direct product. Both are properties
of the character tables, not of the staging code, so they are checked here
against every table the detector can produce -- a non-abelian table, or one
with a character outside +-1, would silently make the shortcut wrong rather
than fail loudly.

Pure stdlib: reads the tables and ``product_irrep`` directly, so it runs in a
source-only checkout.
"""

import importlib.util
import itertools
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


detect = _load("_openqp_symmetry_detect", "pyoqp/oqp/library/symmetry_detect.py")
symmetry = _load("_openqp_symmetry", "pyoqp/oqp/library/symmetry.py")


def _codes(table):
    """The encoding stage_mo_irreps writes, kept in step with it by hand."""
    out = {}
    for name, row in table.items():
        code = 0
        for op, chi in enumerate(row):
            if chi < 0:
                code |= (1 << op)
        out[name] = code
    return out


class MoIrrepEncodingTests(unittest.TestCase):
    def test_every_table_has_plus_minus_one_characters(self):
        for group, table in detect.CHARACTER_TABLES.items():
            for name, row in table.items():
                for chi in row:
                    self.assertIn(
                        chi, (1, -1),
                        "%s/%s has character %r; the bitmask encoding assumes "
                        "+-1 and would be wrong for this table" % (group, name, chi),
                    )

    def test_encoding_is_faithful(self):
        for group, table in detect.CHARACTER_TABLES.items():
            codes = _codes(table)
            self.assertEqual(
                len(set(codes.values())), len(table),
                "%s: two irreps share a code, so a consumer could not tell "
                "them apart: %r" % (group, codes),
            )

    def test_xor_reproduces_the_direct_product(self):
        for group, table in detect.CHARACTER_TABLES.items():
            codes = _codes(table)
            inverse = {code: name for name, code in codes.items()}
            for a, b in itertools.product(table, repeat=2):
                want = symmetry.product_irrep([a, b], table)
                got = inverse.get(codes[a] ^ codes[b])
                self.assertEqual(
                    got, want,
                    "%s: %s x %s is %s but XOR gives %s" % (group, a, b, want, got),
                )

    def test_totally_symmetric_irrep_encodes_to_zero(self):
        """A closed-shell determinant must come out totally symmetric, and it
        only does if the identity of the product group is code 0."""
        for group, table in detect.CHARACTER_TABLES.items():
            codes = _codes(table)
            zeros = [name for name, code in codes.items() if code == 0]
            self.assertEqual(
                len(zeros), 1,
                "%s: expected exactly one all-symmetric irrep, got %r"
                % (group, zeros),
            )
            for name in table:
                self.assertEqual(
                    codes[name] ^ codes[zeros[0]], codes[name],
                    "%s: %s is not the identity of the product" % (group, zeros[0]),
                )


if __name__ == "__main__":
    unittest.main()
