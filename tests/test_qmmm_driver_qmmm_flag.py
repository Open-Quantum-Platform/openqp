"""Regression guard: config-mode QM/MM sets input.qmmm_flag for the native engine.

The Z-vector solvers form the relaxed ESPF charges of the target state, which
the MM electrostatic-embedding forces use, only when ``input.qmmm_flag`` is
set.  ``OpenQpQMMM(..., oqp_cfg=...)`` without it ran the embedded SCF and
returned excited-state forces built from the reference-state charges, silently
inconsistent with the energy (finite-difference check in
tests/test_restricted_td_relaxed_density_consumers.py, ``..._without_qmmm_flag``).

``_with_qmmm_flag`` sets the flag for electrostatic-type embeddings of the
native methods, rejects an explicit false, and leaves mechanical embedding and
the tight-binding backends alone.  In mol mode (``OpenQpQMMM(mol=...)``) the
caller owns the molecule's configuration, so ``_require_mol_qmmm_flag`` raises
the same error instead of changing it.  Needs OpenMM (the driver imports it).
"""
import os
import unittest

OPENMM_REQUIRED = os.getenv("OQP_REQUIRE_OPENMM_TESTS") == "1"

try:
    import oqp.library.qmmm_driver as _drv
    from oqp.library.qmmm_driver import OpenQpQMMM
    _HAVE = True
except Exception:  # optional OpenMM / uncompiled backend
    _HAVE = False
# Looked up separately so that the constructor test below still runs (and
# fails) on a build without the guard instead of skipping the whole module.
_with_qmmm_flag = getattr(_drv, "_with_qmmm_flag", None) if _HAVE else None
_require_mol_qmmm_flag = getattr(_drv, "_require_mol_qmmm_flag", None) if _HAVE else None


class QMMMFlagBuildGate(unittest.TestCase):
    @unittest.skipUnless(OPENMM_REQUIRED, "set OQP_REQUIRE_OPENMM_TESTS=1 where OpenMM is installed")
    def test_driver_importable(self):
        self.assertTrue(_HAVE, "oqp.library.qmmm_driver must import (OQP_REQUIRE_OPENMM_TESTS=1)")


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class QMMMFlagGuardTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(_with_qmmm_flag, "qmmm_driver._with_qmmm_flag is missing")

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


class _FakeMol:
    def __init__(self, method, flag=None):
        self.config = {"input": {"method": method}}
        if flag is not None:
            self.config["input"]["qmmm_flag"] = flag


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class MolModeQMMMFlagTests(unittest.TestCase):
    def test_missing_or_false_flag_is_rejected(self):
        self.assertIsNotNone(_require_mol_qmmm_flag, "qmmm_driver._require_mol_qmmm_flag is missing")
        for flag in (None, False, "False", "0"):
            for embedding in ("electrostatic", "espf", "espf_full", "split"):
                with self.assertRaisesRegex(ValueError, "needs input.qmmm_flag=True"):
                    _require_mol_qmmm_flag(_FakeMol("tdhf", flag), embedding)

    def test_true_flag_mechanical_and_tight_binding_pass(self):
        self.assertIsNotNone(_require_mol_qmmm_flag, "qmmm_driver._require_mol_qmmm_flag is missing")
        for flag in (True, "True", "1"):
            _require_mol_qmmm_flag(_FakeMol("tdhf", flag), "electrostatic")
        _require_mol_qmmm_flag(_FakeMol("tdhf"), "mechanical")
        _require_mol_qmmm_flag(_FakeMol("dftb"), "electrostatic")

    def test_driver_constructor_rejects_mol_without_flag(self):
        # The guard runs in OpenQpQMMM.__init__ before the molecule is used,
        # and the supplied molecule is not modified.
        import openmm.app as app
        import openmm.unit as unit
        from openmm import Vec3
        top = app.Topology()
        res = top.addResidue("HOH", top.addChain())
        o = top.addAtom("O", app.element.oxygen, res)
        h1 = top.addAtom("H1", app.element.hydrogen, res)
        h2 = top.addAtom("H2", app.element.hydrogen, res)
        top.addBond(o, h1)
        top.addBond(o, h2)
        pos = unit.Quantity([Vec3(0, 0, 0), Vec3(0.09572, 0, 0),
                             Vec3(-0.024, 0.0927, 0)], unit.nanometer)
        mol = _FakeMol("tdhf")
        with self.assertRaisesRegex(ValueError, "needs input.qmmm_flag=True"):
            OpenQpQMMM(pos, top, app.ForceField("tip3p.xml"), [0, 1, 2], mol=mol,
                       Cutoff=app.NoCutoff, Embedding="electrostatic")
        self.assertNotIn("qmmm_flag", mol.config["input"])


if __name__ == "__main__":
    unittest.main()
