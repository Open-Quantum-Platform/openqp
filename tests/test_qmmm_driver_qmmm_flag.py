"""Regression guard: config-mode QM/MM sets input.qmmm_flag for the native engine.

The Z-vector solvers form the relaxed ESPF charges of the target state, which
the MM electrostatic-embedding forces use, only when ``input.qmmm_flag`` is
set.  ``OpenQpQMMM(..., oqp_cfg=...)`` without it ran the embedded SCF and
returned excited-state forces built from the reference-state charges, silently
inconsistent with the energy (finite-difference check in
tests/test_restricted_td_relaxed_density_consumers.py, ``..._without_qmmm_flag``).

``_with_qmmm_flag`` sets the flag for electrostatic-type embeddings of the
native methods, rejects an explicit false, and leaves mechanical embedding and
the tight-binding backends alone.  Needs OpenMM (the driver imports it).
"""
import os
import unittest

OPENMM_REQUIRED = os.getenv("OQP_REQUIRE_OPENMM_TESTS") == "1"

try:
    from oqp.library.qmmm_driver import _with_qmmm_flag
    _HAVE = True
except Exception:  # optional OpenMM / uncompiled backend
    _HAVE = False


class QMMMFlagBuildGate(unittest.TestCase):
    @unittest.skipUnless(OPENMM_REQUIRED, "set OQP_REQUIRE_OPENMM_TESTS=1 where OpenMM is installed")
    def test_driver_importable(self):
        self.assertTrue(_HAVE, "oqp.library.qmmm_driver must import (OQP_REQUIRE_OPENMM_TESTS=1)")


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class QMMMFlagGuardTests(unittest.TestCase):
    def test_absent_flag_is_set_for_electrostatic_embeddings(self):
        for embedding in ("electrostatic", "espf", "espf_full", "split"):
            cfg = {"input.method": "tdhf"}
            out = _with_qmmm_flag(cfg, embedding)
            self.assertIs(out, cfg)
            self.assertEqual(cfg["input.qmmm_flag"], "True", embedding)

    def test_explicit_true_is_kept(self):
        for value in (True, "True", "true", "1", "on"):
            cfg = {"input.method": "hf", "input.qmmm_flag": value}
            _with_qmmm_flag(cfg, "electrostatic")
            self.assertEqual(cfg["input.qmmm_flag"], value)

    def test_explicit_false_is_rejected(self):
        for value in (False, "False", "0", "off"):
            with self.assertRaisesRegex(ValueError, "needs input.qmmm_flag=True"):
                _with_qmmm_flag({"input.method": "tdhf", "input.qmmm_flag": value},
                                "electrostatic")

    def test_mechanical_and_tight_binding_are_untouched(self):
        cfg = {"input.method": "tdhf"}
        _with_qmmm_flag(cfg, "mechanical")
        self.assertNotIn("input.qmmm_flag", cfg)
        cfg = {"input.method": "dftb"}
        _with_qmmm_flag(cfg, "electrostatic")
        self.assertNotIn("input.qmmm_flag", cfg)
        cfg = {"input.method": "dftb", "input.qmmm_flag": False}
        _with_qmmm_flag(cfg, "electrostatic")     # TB backends: not rejected
        self.assertIs(cfg["input.qmmm_flag"], False)


if __name__ == "__main__":
    unittest.main()
