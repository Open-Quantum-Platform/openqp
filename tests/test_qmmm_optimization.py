"""QM/MM geometry optimisation (runtype = optimize with qmmm_flag = true).

Bookkeeping that needs no compiled backend: the input checker admits plain
optimisation and still rejects the reaction-path drivers, validates the
movable-shell radius, and the driver maps [optimize] istate onto the gradient
root and selects whole MM residues within the radius.  The behavioural check
(the optimisation actually lowers the QM/MM energy and reduces the gradient)
lives in the shipped example examples/QMMM/ala-dipeptide_RHF-QMMM-OPT-linkatom.
"""
import importlib.util
import sys
import types
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def _load_checker():
    if "oqp.utils.mpi_utils" not in sys.modules:      # the checker imports it at module level
        mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
        class MPIManager:
            use_mpi = False
            size = 1
        mpi_utils.MPIManager = MPIManager
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    spec = importlib.util.spec_from_file_location(
        "input_checker_qmmm_opt_under_test", ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


class TestCheckerAdmitsPlainQMMMOptimisation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chk = _load_checker()

    def _diags(self, runtype, optimize=None, method="hf", qmmm_extra=None):
        cfg = {"input": {"runtype": runtype, "qmmm_flag": True, "method": method, "basis": "6-31g",
                         "system": "ala.pdb 9 10 17 18 19", "charge": 0},
               "scf": {"type": "rhf", "multiplicity": 1},
               "optimize": {"istate": 0, **(optimize or {})},
               "qmmm": {"pdb_file": "ala.pdb", "qm_atoms": "8,9,16,17,18",
                        "forcefield_files": "amber14-all.xml", "cutoff": "NoCutoff", **(qmmm_extra or {})}}
        if method == "tdhf":
            cfg["scf"] = {"type": "rohf", "multiplicity": 3}
            cfg["tdhf"] = {"type": "mrsf", "nstate": 2, "multiplicity": 1}
            cfg["optimize"]["istate"] = 1
        report = self.chk.CheckReport()
        self.chk.check_input_values(cfg, raise_error=False, emit=False) if False else None
        for fn in ("_check_qmmm_driver_options", "_check_optimize"):
            if hasattr(self.chk, fn):
                getattr(self.chk, fn)(cfg, report)
        return [(d.severity, d.path) for d in report.diagnostics
                if "qmmm" in d.path or d.path in ("optimize.lib", "input.method", "optimize.freeze")]

    def test_optimize_is_admitted_and_paths_are_not(self):
        self.assertNotIn(("ERROR", "input.qmmm_flag"), self._diags("optimize"))
        for rt in ("meci", "mecp", "ts", "irc", "neb", "mep"):
            self.assertIn(("ERROR", "input.qmmm_flag"), self._diags(rt), rt)

    def test_optimize_keeps_the_periodic_and_embedding_controls(self):
        """The optimiser builds the same driver as runtype=md and passes every
        periodic/embedding control through, so they must not be flagged as
        ignored the way they are for a single point."""
        cfg = {"qmmm_radius": 4.0}
        for extra in ({"cutoff": "PME", "ewald_tol": "1e-6"}, {"mm_charge_width": "0.7"},
                      {"lj_switch": True, "h_lj": True}):
            diags = self._diags("optimize", cfg, qmmm_extra=extra)
            self.assertFalse([d for d in diags if d[1] == "qmmm.cutoff"], extra)
        self.assertIn(("ERROR", "qmmm.embedding"),
                      self._diags("optimize", cfg, qmmm_extra={"cutoff": "PME", "embedding": "split"}))

    def test_method_and_constraints_are_gated(self):
        for method in ("hf", "tdhf", "dftb"):
            self.assertNotIn(("ERROR", "input.method"), self._diags("optimize", method=method))
        for method in ("casscf", "mp2"):
            self.assertIn(("ERROR", "input.method"), self._diags("optimize", method=method))
        self.assertIn(("ERROR", "optimize.freeze"), self._diags("optimize", {"freeze": "distance(1,2)"}))

    def test_radius_and_lib_are_validated(self):
        self.assertIn(("ERROR", "optimize.qmmm_radius"), self._diags("optimize", {"qmmm_radius": -1.0}))
        self.assertNotIn(("ERROR", "optimize.qmmm_radius"), self._diags("optimize", {"qmmm_radius": 4.0}))
        self.assertIn(("WARNING", "optimize.lib"), self._diags("optimize", {"lib": "geometric"}))


try:
    from oqp.library.qmmm_opt import QMMM_Opt
    _HAVE = True
except Exception:  # pragma: no cover - OpenMM / backend missing
    _HAVE = False


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class TestMovableSetAndState(unittest.TestCase):
    def _bare(self, radius, box_bohr=None):
        import openmm.app as app
        from openmm import Vec3, unit
        o = object.__new__(QMMM_Opt)
        top = app.Topology(); ch = top.addChain()
        pos = []
        # residue A: two QM atoms at the origin plus two MM atoms of the SAME
        # residue (a covalent cut inside a residue), one 1.5 A and one 4.5 A away
        ra = top.addResidue("QMA", ch)
        for nm, p in (("C1", (0, 0, 0)), ("C2", (1.2, 0, 0)), ("CA", (-1.5, 0, 0)), ("CB", (-4.5, 0, 0))):
            top.addAtom(nm, app.element.carbon, ra); pos.append(p)
        # residue B: a water 3 A away, residue C: a water 8 A away
        for name, x in (("HOH", 3.0), ("HOH", 8.0)):
            rb = top.addResidue(name, ch)
            for nm, el, dp in (("O", app.element.oxygen, 0.0), ("H1", app.element.hydrogen, 0.6), ("H2", app.element.hydrogen, -0.6)):
                top.addAtom(nm, el, rb); pos.append((x, dp, 0.0))
        o.pdb = SimpleNamespace(topology=top, positions=unit.Quantity([Vec3(*p) for p in pos], unit.angstrom))
        o.qm_atoms = np.array([0, 1]); o.driver = SimpleNamespace(_box_lengths_bohr=lambda: box_bohr)
        return o

    def test_whole_residues_within_radius_move(self):
        # atoms: 0,1 QM | 2 CA (1.5 A) 3 CB (4.5 A) same residue | 4-6 water at 3 A | 7-9 water at 8 A
        o = self._bare(0.0)
        np.testing.assert_array_equal(QMMM_Opt._movable_atoms(o, 0.0), [0, 1])
        # distances to the NEAREST QM atom: CA 1.5, near-water O 1.8 (from C2), CB 4.5, far-water O 6.8
        np.testing.assert_array_equal(QMMM_Opt._movable_atoms(o, 1.6), [0, 1, 2])          # the link host only
        np.testing.assert_array_equal(QMMM_Opt._movable_atoms(o, 2.0), [0, 1, 2, 4, 5, 6]) # + the near water, whole
        np.testing.assert_array_equal(QMMM_Opt._movable_atoms(o, 5.0), [0, 1, 2, 3, 4, 5, 6])   # + CB, not the far water
        np.testing.assert_array_equal(QMMM_Opt._movable_atoms(o, 7.0), np.arange(10))      # everything
        o = self._bare(0.0, box_bohr=np.array([10.0, 10.0, 10.0]) * 1.8897259886)         # 10 A cell
        np.testing.assert_array_equal(QMMM_Opt._movable_atoms(o, 2.5),                    # far water is 2 A away through the face
                                      [0, 1, 2, 4, 5, 6, 7, 8, 9])

    def test_list_parsers_accept_the_api_and_namd_forms(self):
        from oqp.library.qmmm_md import _parse_int_list, _parse_str_list
        self.assertEqual(_parse_int_list("0 1 2"), [0, 1, 2])           # job.qmmm(qm_atoms=[0,1,2]) serialises to this
        self.assertEqual(_parse_int_list("0-3, 8 9"), [0, 1, 2, 3, 8, 9])
        self.assertEqual(_parse_str_list("amber14-all.xml amber14/tip3p.xml"), ["amber14-all.xml", "amber14/tip3p.xml"])
        self.assertEqual(_parse_str_list("a.xml,b.xml"), ["a.xml", "b.xml"])
        import os, tempfile
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "my forcefield.xml"); open(p, "w").close()
            self.assertEqual(_parse_str_list(p), [p])                  # an existing path with a space stays whole

    def test_istate_is_the_gradient_root(self):
        src = (ROOT / "pyoqp" / "oqp" / "library" / "qmmm_opt.py").read_text()
        self.assertIn('mol.config.setdefault("properties", {})["grad"] = [self.istate]', src)


if __name__ == "__main__":
    unittest.main()
