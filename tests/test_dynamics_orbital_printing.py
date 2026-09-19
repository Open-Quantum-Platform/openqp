"""A dynamics run must not repeat the MO coefficient table for every step.

The SCF prints the whole table on every call and a trajectory calls the SCF at
least once per step, so a 100-step QM/MM NAMD run of an 18-atom QM region wrote
405 tables and 83 MB of log.  ``verbose = 0`` now suppresses the table in
``source/printing.F90``, and ``runtype = md`` / ``namd`` default to it.  The
default level lists orbital energies only; ``verbose >= 2`` adds coefficients.
"""
import os
import re
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
runtype={runtype}
{functional}{input_verbose}
[guess]
type=huckel
save_mol=false
[scf]
type=rhf
multiplicity=1
conv=1.0e-8
{verbose}"""

# One row of the coefficient table: index, basis-function label, coefficients.
COEFFICIENT_ROW = re.compile(r"^\s*\d+\s+H\s+\d+\s+S\s+-?\d+\.\d{10}", re.M)


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


def _verbose_after(runtype, verbose=None, input_verbose=None):
    """Run the parser-level default exactly as Molecule.get_config does."""
    parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA, allow_no_value=True)
    parser.set("input", "runtype", runtype)
    if verbose is not None:
        parser.set("scf", "verbose", str(verbose))
    if input_verbose is not None:
        parser.set("input", "verbose", str(input_verbose))
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

    def test_an_explicit_global_level_is_kept(self):
        """``[input] verbose`` is the global switch: asking for more there keeps
        the legacy ``[scf]`` value at its default, and the resolver takes 2."""
        from oqp.utils.log_format import resolve_verbosity
        self.assertEqual(_verbose_after("namd", input_verbose=2), "1")
        self.assertEqual(resolve_verbosity({"input": {"verbose": 2},
                                            "scf": {"verbose": 1}}), 2)

    def test_every_other_runtype_still_prints(self):
        for runtype in ("energy", "grad", "optimize", "hess", "soc"):
            self.assertEqual(_verbose_after(runtype), "1", runtype)
            self.assertEqual(_verbose_after(runtype, 2), "2", runtype)

    def test_one_level_reaches_every_native_gate(self):
        """Both spellings resolve to one level, which also drives the MRSF
        developer dumps that no input could switch on before."""
        from types import SimpleNamespace
        from oqp.molecule.oqpdata import OQPData
        data = OQPData.__new__(OQPData)
        data._data = SimpleNamespace(control=SimpleNamespace(verbose=None),
                                     tddft=SimpleNamespace(debug_mode=None))
        data.set_input_verbose(1)
        data.set_scf_verbose(0)
        self.assertEqual(data._data.control.verbose, 0)
        self.assertFalse(data._data.tddft.debug_mode)
        data.set_input_verbose(3)
        self.assertEqual(data._data.control.verbose, 3)
        self.assertTrue(data._data.tddft.debug_mode)

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
        self.assertIn("energies_only=(infos%control%verbose < 2)", head)


@unittest.skipUnless(_HAVE and _runtime_available(),
                     "compiled OpenQP runtime unavailable")
class TestOrbitalTableIsActuallySuppressed(unittest.TestCase):
    """The behavioural half: run the SCF both ways and count the tables.

    This is the arm that fails if the Fortran gate is removed; the checks above
    only pin the plumbing that decides the verbosity."""

    def _log(self, verbose=None, section="scf", runtype="energy", functional=""):
        from oqp.pyoqp import Runner
        line = "" if verbose is None else f"verbose={verbose}\n"
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "h2.inp"
            inp.write_text(DECK.format(
                runtype=runtype,
                functional=f"functional={functional}\n" if functional else "",
                input_verbose=line if section == "input" else "",
                verbose=line if section == "scf" else ""))
            log = Path(tmp) / "h2.log"
            runner = Runner(project="h2", input_file=str(inp), log=str(log),
                            silent=1, usempi=False)
            runner.run()
            return log.read_text(errors="replace")

    def _run(self, verbose, section="scf"):
        text = self._log(verbose, section)
        # "Final RHF energy is" is written at every level; the energy
        # components block that carries "TOTAL energy =" is level 1 and up.
        energy = float(text.rsplit("Final RHF energy is", 1)[1].split()[0])
        return text.count("Molecular Orbitals and Energies"), energy

    def test_verbose_zero_removes_the_table_and_changes_nothing_else(self):
        n_default, e_default = self._run(None)
        n_quiet, e_quiet = self._run(0)
        self.assertEqual(n_default, 1)          # control: the table is there by default
        self.assertEqual(n_quiet, 0)            # the feature
        self.assertEqual(e_default, e_quiet)    # and it is numerically inert

    def test_input_verbose_is_the_same_switch(self):
        n_quiet, e_quiet = self._run(0, section="input")
        _, e_default = self._run(None)
        self.assertEqual(n_quiet, 0)
        self.assertEqual(e_quiet, e_default)

    def test_default_lists_energies_and_detail_adds_coefficients(self):
        normal = self._log(None)
        detailed = self._log(2, section="input")
        self.assertEqual(normal.count("Molecular Orbitals and Energies"), 1)
        self.assertEqual(len(COEFFICIENT_ROW.findall(normal)), 0)
        self.assertGreater(len(COEFFICIENT_ROW.findall(detailed)), 0)

    def test_functional_references_are_written_once(self):
        """The SCF and the gradient both set the functional and grid up; before
        the fix each set-up repeated the LibXC header, the grid description and
        every functional reference."""
        text = self._log(None, runtype="grad", functional="bhhlyp")
        announced = set(re.findall(r"^The .* functional will be used with a coefficient.*$",
                                   text, re.M))
        self.assertGreater(len(announced), 0)
        self.assertEqual(text.count("The functional has been described"), len(announced))
        self.assertEqual(text.count("The libXC interfaces are described"), 1)
        self.assertEqual(len(re.findall(r"Lebedev grid-based DFT options|Standard Grid", text)), 1)

    def test_appended_evaluations_share_one_description(self):
        """QM/MM optimisation and dynamics build a Runner per geometry and
        append to one log; the set-up is described in that log once."""
        from oqp.pyoqp import Runner
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "h2.inp"
            inp.write_text(DECK.format(runtype="energy", functional="functional=bhhlyp\n",
                                       input_verbose="", verbose=""))
            log = str(Path(tmp) / "h2.log")
            os.chdir(tmp)
            try:
                Runner(project="h2", input_file=str(inp), log=log, silent=1, usempi=False).run()
                Runner(project="h2", input_file=str(inp), log=log, silent=1, usempi=False,
                       append_log=True).run()
            finally:
                os.chdir(cwd)
            text = Path(log).read_text(errors="replace")
        self.assertEqual(text.count("OpenQP: Open Quantum Platform"), 1)   # control: one run log
        self.assertEqual(text.count("Final RHF energy is"), 2)             # both evaluations ran
        self.assertEqual(text.count("The libXC interfaces are described"), 1)
        self.assertEqual(len(re.findall(r"Lebedev grid-based DFT options|Standard Grid", text)), 1)

    def test_each_log_gets_one_description_however_runs_are_ordered(self):
        """Runners built before any of them runs, runs interleaved across two
        logs, and an evaluation appended to the first log: each log file carries
        its own description exactly once."""
        from oqp.pyoqp import Runner
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "h2.inp"
            inp.write_text(DECK.format(runtype="energy", functional="functional=bhhlyp\n",
                                       input_verbose="", verbose=""))
            log_a, log_b = str(Path(tmp) / "a.log"), str(Path(tmp) / "b.log")
            os.chdir(tmp)
            try:
                first = Runner(project="a", input_file=str(inp), log=log_a, silent=1, usempi=False)
                second = Runner(project="b", input_file=str(inp), log=log_b, silent=1, usempi=False)
                first.run()
                second.run()
                Runner(project="a", input_file=str(inp), log=log_a, silent=1, usempi=False,
                       append_log=True).run()
            finally:
                os.chdir(cwd)
            a = Path(log_a).read_text(errors="replace")
            b = Path(log_b).read_text(errors="replace")
        self.assertEqual(a.count("Final RHF energy is"), 2)    # control: both evaluations ran
        self.assertEqual(b.count("Final RHF energy is"), 1)
        for text in (a, b):
            self.assertEqual(text.count("The libXC interfaces are described"), 1)
            self.assertGreater(text.count("The functional has been described"), 0)

    def test_a_recreated_log_gets_its_description_again(self):
        """Two separate runs that reuse one log path in one process: the second
        run starts the file afresh, so it must describe the set-up again."""
        from oqp.pyoqp import Runner
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "h2.inp"
            inp.write_text(DECK.format(runtype="energy", functional="functional=bhhlyp\n",
                                       input_verbose="", verbose=""))
            log = str(Path(tmp) / "h2.log")
            os.chdir(tmp)
            try:
                Runner(project="h2", input_file=str(inp), log=log, silent=1, usempi=False).run()
                Runner(project="h2", input_file=str(inp), log=log, silent=1, usempi=False).run()
            finally:
                os.chdir(cwd)
            text = Path(log).read_text(errors="replace")
        self.assertEqual(text.count("Final RHF energy is"), 1)    # control: the file was started afresh
        self.assertEqual(text.count("The libXC interfaces are described"), 1)
        self.assertGreater(text.count("The functional has been described"), 0)
        self.assertEqual(len(re.findall(r"Lebedev grid-based DFT options|Standard Grid", text)), 1)

    def test_every_run_describes_its_functional(self):
        """The once-per-run records are reset at the start of a run, so a second
        run in the same Python process still documents its functional."""
        for _ in range(2):
            text = self._log(None, functional="bhhlyp")
            self.assertEqual(text.count("The libXC interfaces are described"), 1)
            self.assertGreater(text.count("The functional has been described"), 0)


if __name__ == "__main__":
    unittest.main()
