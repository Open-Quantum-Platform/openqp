"""Dependency-light configuration tests for the outer QM/MM MD driver."""

import unittest

try:
    from oqp.library.qmmm_md import _qm_subcall_config
    _HAVE_QMMM_MD = True
except Exception:  # pragma: no cover - optional OpenMM/OpenQP runtime
    _HAVE_QMMM_MD = False


@unittest.skipUnless(_HAVE_QMMM_MD, "OpenMM or compiled OpenQP unavailable")
class TestQMMMMDConfig(unittest.TestCase):
    def test_md_runtype_is_normalized_for_qm_subcalls(self):
        outer = {
            "input.runtype": "md",
            "input.method": "dftb",
            "scf.type": "rhf",
        }

        inner = _qm_subcall_config(outer)

        self.assertEqual(inner["input.runtype"], "energy")
        self.assertEqual(inner["input.method"], "dftb")
        self.assertEqual(outer["input.runtype"], "md")

    def test_mol_mode_keeps_none_config(self):
        self.assertIsNone(_qm_subcall_config(None))


if __name__ == "__main__":
    unittest.main()
