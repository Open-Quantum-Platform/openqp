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
import tempfile
import unittest
from pathlib import Path
from unittest import mock


def _runtime_available():
    try:
        import oqp  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "built OpenQP runtime not importable")
class FciSymmetrySelfTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="oqp_fci_symmetry_")
        self.addCleanup(tmp.cleanup)
        self.tmpdir = Path(tmp.name)
        # The native self-test writes its report to $OQP_SELFTEST_DIR
        # (source/io/selftest_output.F90); keep it private to this test.
        env = mock.patch.dict(os.environ, {"OQP_SELFTEST_DIR": str(self.tmpdir)})
        env.start()
        self.addCleanup(env.stop)

    def test_determinant_and_root_irreps(self):
        import oqp

        try:
            oqp.ffi.cdef("void fci_symmetry_selftest(void);")
        except Exception:
            # cdef is process-global; a second declaration is harmless.
            pass
        lib = oqp.ffi.dlopen(str(_library_path()))
        lib.fci_symmetry_selftest()
        report = self.tmpdir / "fci_symmetry_selftest.out"

        self.assertTrue(report.exists(), "selftest produced no output")
        results = {}
        for line in report.read_text().splitlines():
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
    # The liboqp that `import oqp` resolved, not $OPENQP_ROOT: the variable is
    # only a fallback in oqp.runtime, and reading it here tested whatever tree
    # it happened to name (usually the unbuilt source tree, so the test skipped).
    from oqp.runtime import library_path, resolve_oqp_root

    root, suffix = resolve_oqp_root()
    path = library_path(root, suffix)
    if path is None:
        raise AssertionError(f"oqp resolved runtime root {root} but it has no liboqp.{suffix}")
    return path


if __name__ == "__main__":
    unittest.main()
