"""A dynamics run must not repeat the MO coefficient table for every step.

The SCF prints the whole table on every call and a trajectory calls the SCF at
least once per step, so a 100-step QM/MM NAMD run of an 18-atom QM region wrote
405 tables and 83 MB of log.  ``[scf] verbose = 0`` now suppresses the table in
``source/printing.F90``, and ``runtype = md`` / ``namd`` default to it.
"""
import os
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DECK = """[input]
system=
   1   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.740000000
charge=0
method=hf
basis=6-31g
runtype=energy
[guess]
type=huckel
save_mol=false
[scf]
type=rhf
multiplicity=1
conv=1.0e-8
{verbose}"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False

try:
    from oqp.molecule import Molecule
    from oqp.utils.input_parser import OQPConfigParser
    from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
    _HAVE = True
except Exception:  # pragma: no cover - uncompiled backend
    _HAVE = False


def _verbose_after(runtype, verbose=None):
    """Run the parser-level default exactly as Molecule.get_config does."""
    parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA, allow_no_value=True)
    parser.set("input", "runtype", runtype)
    if verbose is not None:
        parser.set("scf", "verbose", str(verbose))
    Molecule._quiet_orbitals_in_dynamics(parser)
    return parser.get("scf", "verbose")


@unittest.skipUnless(_HAVE, "compiled OpenQP backend unavailable")
class TestDynamicsOrbitalPrinting(unittest.TestCase):
    def test_dynamics_defaults_to_silent_orbitals(self):
        for runtype in ("md", "namd", "MD", " NAMD "):
            self.assertEqual(_verbose_after(runtype), "0", runtype)
            self.assertEqual(_verbose_after(runtype, 1), "0", runtype)

    def test_an_explicit_request_for_detail_is_kept(self):
        for v in (2, 3):
            self.assertEqual(_verbose_after("namd", v), str(v))
        self.assertEqual(_verbose_after("namd", 0), "0")

    def test_every_other_runtype_still_prints(self):
        for runtype in ("energy", "grad", "optimize", "hess", "soc"):
            self.assertEqual(_verbose_after(runtype), "1", runtype)
            self.assertEqual(_verbose_after(runtype, 2), "2", runtype)

    def test_the_config_mode_md_path_applies_the_same_default(self):
        """QMMM_MD(oqp_cfg=...) rewrites runtype to 'energy' before the molecule
        is built, so the parser-level default cannot see the dynamics."""
        src = (ROOT / "pyoqp" / "oqp" / "library" / "qmmm_md.py").read_text()
        block = src[src.index("Force the internal QM runtype"):]
        block = block[:block.index("self.oqp_cfg = qm_cfg")]
        self.assertIn("scf.verbose", block)
        self.assertIn("'0'", block)
        self.assertLess(block.index("= 'energy'"), block.index("scf.verbose"))

    def test_the_fortran_gate_exists(self):
        src = (ROOT / "source" / "printing.F90").read_text()
        head = src[src.index("subroutine print_mo_range"):src.index("end subroutine print_mo_range")]
        self.assertIn("infos%control%verbose < 1", head)
        self.assertLess(head.index("verbose < 1"), head.index("Molecular Orbitals and Energies"))


@unittest.skipUnless(_HAVE and _runtime_available(),
                     "compiled OpenQP runtime unavailable")
class TestOrbitalTableIsActuallySuppressed(unittest.TestCase):
    """The behavioural half: run the SCF both ways and count the tables.

    This is the arm that fails if the Fortran gate is removed; the checks above
    only pin the plumbing that decides the verbosity."""

    def _run(self, verbose):
        from oqp.pyoqp import Runner
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "h2.inp"
            inp.write_text(DECK.format(verbose="" if verbose is None else f"verbose={verbose}\n"))
            log = Path(tmp) / "h2.log"
            runner = Runner(project="h2", input_file=str(inp), log=str(log),
                            silent=1, usempi=False)
            runner.run()
            text = log.read_text(errors="replace")
        energy = float(text.rsplit("TOTAL energy =", 1)[1].split()[0])
        return text.count("Molecular Orbitals and Energies"), energy

    def test_verbose_zero_removes_the_table_and_changes_nothing_else(self):
        n_default, e_default = self._run(None)
        n_quiet, e_quiet = self._run(0)
        self.assertEqual(n_default, 1)          # control: the table is there by default
        self.assertEqual(n_quiet, 0)            # the feature
        self.assertEqual(e_default, e_quiet)    # and it is numerically inert


if __name__ == "__main__":
    unittest.main()
