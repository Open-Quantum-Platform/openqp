"""A dynamics run must not repeat the MO coefficient table for every step.

The SCF prints the whole table on every call and a trajectory calls the SCF at
least once per step, so a 100-step QM/MM NAMD run of an 18-atom QM region wrote
405 tables and 83 MB of log.  ``[scf] verbose = 0`` now suppresses the table in
``source/printing.F90``, and ``runtype = md`` / ``namd`` default to it.
"""
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]

try:
    from oqp.molecule import Molecule
    _HAVE = True
except Exception:  # pragma: no cover - uncompiled backend
    _HAVE = False


def _quiet(runtype, scf=None):
    mol = object.__new__(Molecule)
    mol.config = {"input": {"runtype": runtype}}
    if scf is not None:
        mol.config["scf"] = dict(scf)
    Molecule._quiet_orbitals_in_dynamics(mol)
    return mol.config.get("scf", {}).get("verbose", 1)


@unittest.skipUnless(_HAVE, "compiled OpenQP backend unavailable")
class TestDynamicsOrbitalPrinting(unittest.TestCase):
    def test_dynamics_defaults_to_silent_orbitals(self):
        for runtype in ("md", "namd", "MD", " NAMD "):
            self.assertEqual(_quiet(runtype), 0, runtype)
            self.assertEqual(_quiet(runtype, {"verbose": 1}), 0, runtype)

    def test_an_explicit_request_for_detail_is_kept(self):
        for v in (2, 3):
            self.assertEqual(_quiet("namd", {"verbose": v}), v)
        self.assertEqual(_quiet("namd", {"verbose": 0}), 0)

    def test_every_other_runtype_still_prints(self):
        for runtype in ("energy", "grad", "optimize", "hess", "soc", ""):
            self.assertEqual(_quiet(runtype), 1, runtype)
            self.assertEqual(_quiet(runtype, {"verbose": 1}), 1, runtype)

    def test_the_fortran_gate_exists(self):
        src = (ROOT / "source" / "printing.F90").read_text()
        head = src[src.index("subroutine print_mo_range"):src.index("end subroutine print_mo_range")]
        self.assertIn("infos%control%verbose < 1", head)
        self.assertLess(head.index("verbose < 1"), head.index("Molecular Orbitals and Energies"))


if __name__ == "__main__":
    unittest.main()
