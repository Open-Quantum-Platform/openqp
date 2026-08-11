"""Integration test for determinant/CI-root point-group classification.

Drives the ``fci_symmetry_selftest`` bind(C) harness in
``tests/fortran/fci_symmetry_selftest.F90``. That harness needs no molecule --
it builds determinants and CI vectors by hand -- so the properties it pins are
exact rather than tolerance-based:

  * a closed-shell determinant is totally symmetric whatever the orbital
    irreps are (the XOR cancels because every orbital enters twice);
  * a b2 -> b1 single excitation carries a2, the C2v product;
  * a CI vector confined to one irrep classifies with purity 1, and a
    deliberate 50/50 mixture reports 0.5, so a caller filtering on irrep can
    tell a symmetry eigenstate from a mixture;
  * selection returns matching roots in energy order and rejects impure ones.

Skipped unless a built OpenQP is importable, like the other selftest drivers.
"""

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELFTEST_OUT = Path("/tmp/fci_symmetry_selftest.out")


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        import oqp  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "built OpenQP runtime not importable")
class FciSymmetrySelfTest(unittest.TestCase):
    def test_determinant_and_root_irreps(self):
        import oqp

        if SELFTEST_OUT.exists():
            SELFTEST_OUT.unlink()

        try:
            oqp.ffi.cdef("void fci_symmetry_selftest(void);")
        except Exception:
            # cdef is process-global; a second declaration is harmless.
            pass
        lib = oqp.ffi.dlopen(str(_library_path()))
        lib.fci_symmetry_selftest()

        self.assertTrue(SELFTEST_OUT.exists(), "selftest produced no output")
        results = {}
        for line in SELFTEST_OUT.read_text().splitlines():
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            results[key.strip()] = value.strip()

        for key in (
            "closed_shell_is_totally_symmetric",
            "b2_to_b1_excitation_is_a2",
            "pure_roots_classify_with_purity_1",
            "mixed_root_reports_purity_half",
            "selection_skips_impure_root",
        ):
            self.assertEqual(results.get(key), "T", "%s failed: %s" % (key, results))
        self.assertEqual(results.get("ALL_PASS"), "T", results)


def _library_path():
    root = Path(os.environ["OPENQP_ROOT"])
    for pattern in ("lib/liboqp.dylib", "lib/liboqp.so", "lib64/liboqp.so"):
        candidate = root / pattern
        if candidate.exists():
            return candidate
    matches = sorted(root.glob("lib*/liboqp.*"))
    if not matches:
        raise unittest.SkipTest("liboqp not found under OPENQP_ROOT")
    return matches[0]


if __name__ == "__main__":
    unittest.main()
