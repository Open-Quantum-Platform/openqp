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

    def _diags(self, runtype, optimize=None, method="hf", qmmm_extra=None, input_extra=None):
        cfg = {"input": {"runtype": runtype, "qmmm_flag": True, "method": method, "basis": "6-31g",
                         "system": "ala.pdb 9 10 17 18 19", "charge": 0, **(input_extra or {})},
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
                if "qmmm" in d.path or d.path in ("optimize.lib", "optimize.maxit", "input.method",
                                                  "input.d4", "optimize.freeze")]

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

    def test_blank_required_settings_are_missing(self):
        for key in ("pdb_file", "qm_atoms", "forcefield_files"):
            self.assertIn(("ERROR", f"qmmm.{key}"), self._diags("optimize", qmmm_extra={key: ""}))

    def test_unsupported_options_are_rejected_not_ignored(self):
        self.assertIn(("ERROR", "optimize.maxit"), self._diags("optimize", {"maxit": 0}))
        self.assertIn(("ERROR", "qmmm.qm_atoms_xyz"), self._diags("optimize", qmmm_extra={"qm_atoms_xyz": "qm.xyz"}))
        self.assertIn(("ERROR", "qmmm.qm_list"), self._diags("optimize", qmmm_extra={"qm_list": "0,1,2"}))
        self.assertIn(("ERROR", "input.d4"), self._diags("optimize", input_extra={"d4": True}))
        self.assertNotIn(("ERROR", "input.d4"), self._diags("optimize", input_extra={"d4": False}))

    def test_radius_and_lib_are_validated(self):
        self.assertIn(("ERROR", "optimize.qmmm_radius"), self._diags("optimize", {"qmmm_radius": -1.0}))
        for bad in (float("nan"), float("inf"), "nan", "inf"):
            self.assertIn(("ERROR", "optimize.qmmm_radius"), self._diags("optimize", {"qmmm_radius": bad}), bad)
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

    def test_driver_honours_init_scf_and_resolves_the_deck_path(self):
        src = (ROOT / "pyoqp" / "oqp" / "library" / "qmmm_opt.py").read_text()
        self.assertIn("self.driver._reuse_orbitals = not self.init_scf", src)
        self.assertIn("self.pdb = app.PDBFile(self._resolve_aux_file(pdb_file))", src)
        self.assertIn("max(1, int(opt.get(\"maxit\", 30)))", src)
        self.assertIn('self._forcefield_paths(qmmm_cfg.get("forcefield_files", ""))', src)
        self.assertIn('self.coordsys = self._resolve_coordsys(eng.get("coordsys", "auto"))', src)
        self.assertIn('return "cartesian"', src)                       # auto means Cartesian

    def test_istate_is_the_gradient_root(self):
        src = (ROOT / "pyoqp" / "oqp" / "library" / "qmmm_opt.py").read_text()
        self.assertIn('mol.config.setdefault("properties", {})["grad"] = [self.istate]', src)


def _runtime_available():
    import os
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        import openmm  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_HAVE and _runtime_available(), "OpenMM or compiled OpenQP runtime unavailable")
class TestOptimisationPublishesItsResult(unittest.TestCase):
    """Two evaluations of the shipped example through Runner, from a directory
    that is not the deck's: the deck-relative PDB and force-field resolution,
    and Runner.results()['energy'] carrying the QM/MM objective."""

    def test_runner_result_holds_the_qmmm_energy(self):
        import os, tempfile
        from oqp.pyoqp import Runner
        deck = ROOT / "examples" / "QMMM" / "ala-dipeptide_RHF-QMMM-OPT-linkatom.inp"
        text = deck.read_text().replace("maxit=12", "maxit=2").replace("save_mol=true", "save_mol=false")
        natoms = sum(1 for line in (deck.parent / "ala.pdb").read_text().splitlines()
                     if line.startswith(("ATOM", "HETATM")))
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "opt.inp"
            # the deck names ala.pdb relative to itself; run from elsewhere
            inp.write_text(text.replace("pdb_file=ala.pdb", f"pdb_file={deck.parent / 'ala.pdb'}")
                               .replace("system=ala.pdb", f"system={deck.parent / 'ala.pdb'}"))
            cwd = os.getcwd(); os.chdir(tmp)
            try:
                r = Runner(project="opt", input_file=str(inp), log=str(Path(tmp) / "opt.log"),
                           silent=1, usempi=False)
                r.run()
                res = r.results()
            finally:
                os.chdir(cwd)
        e = res["energy"]
        self.assertEqual(len(e), 1)
        self.assertTrue(np.isfinite(e[0]))
        self.assertLess(e[0], -168.0)                       # the QM/MM total, not a fragment or MM piece
        self.assertEqual(r.mol.qmmm_optimization["evaluations"], 2)
        self.assertFalse(r.mol.qmmm_optimization["recovery"])       # the deck switches it off
        self.assertAlmostEqual(r.mol.qmmm_optimization["energy_hartree"], e[0], places=12)
        # the gradient published with that energy, shaped like the published
        # atoms: the objective's derivative for each real QM atom and a zero
        # row for the link hydrogen; the full-system gradient stays on the driver
        self.assertEqual(len(res["grad"]), 1)
        g = np.asarray(res["grad"][0])
        nqm = len(r.qmmm_opt.qm_atoms)
        self.assertEqual(g.shape, (len(res["atoms"]), 3))
        self.assertTrue(np.all(g[nqm:] == 0.0))
        self.assertAlmostEqual(float(np.abs(g[:nqm]).max()), r.mol.qmmm_optimization["max_grad"], places=10)
        self.assertEqual(np.asarray(r.qmmm_opt.gradient_full).shape, (natoms, 3))
        np.testing.assert_allclose(g[:nqm], np.asarray(r.qmmm_opt.gradient_full)[r.qmmm_opt.qm_atoms], rtol=0, atol=0)

    def test_recovery_stage_restarts_an_unconverged_search(self):
        """maxit=1 cannot converge; with the default auto_recovery the search
        restarts from its best geometry for recovery_maxit more evaluations."""
        import os, tempfile
        from oqp.pyoqp import Runner
        deck = ROOT / "examples" / "QMMM" / "ala-dipeptide_RHF-QMMM-OPT-linkatom.inp"
        text = (deck.read_text().replace("maxit=12", "maxit=1").replace("save_mol=true", "save_mol=false")
                .replace("auto_recovery=false", "auto_recovery=true\nrecovery_maxit=1\nrecovery_trust=0.02"))
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "rec.inp"
            inp.write_text(text.replace("pdb_file=ala.pdb", f"pdb_file={deck.parent / 'ala.pdb'}")
                               .replace("system=ala.pdb", f"system={deck.parent / 'ala.pdb'}"))
            cwd = os.getcwd(); os.chdir(tmp)
            try:
                r = Runner(project="rec", input_file=str(inp), log=str(Path(tmp) / "rec.log"),
                           silent=1, usempi=False)
                r.run()
                log = (Path(tmp) / "rec.log").read_text()
            finally:
                os.chdir(cwd)
        info = r.mol.qmmm_optimization
        self.assertTrue(info["recovery"])
        self.assertEqual(info["evaluations"], 2)
        self.assertIn("recovery selected after the iteration limit", log)
        self.assertIn("trust=0.020", log)

    def test_electronic_failure_restarts_from_the_best_geometry(self):
        """An SCF that does not converge at the second trial geometry rejects
        it and enters the recovery stage from the first (best) geometry."""
        import os, tempfile
        from unittest import mock
        from oqp.pyoqp import Runner
        deck = ROOT / "examples" / "QMMM" / "ala-dipeptide_RHF-QMMM-OPT-linkatom.inp"
        text = (deck.read_text().replace("maxit=12", "maxit=2").replace("save_mol=true", "save_mol=false")
                .replace("auto_recovery=false", "auto_recovery=true\nrecovery_maxit=1"))
        real = QMMM_Opt._energy_force
        calls = [0]

        def flaky(self, X):
            calls[0] += 1
            if calls[0] == 2:
                raise RuntimeError("SCF did not converge in 200 iterations")
            return real(self, X)

        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "ef.inp"
            inp.write_text(text.replace("pdb_file=ala.pdb", f"pdb_file={deck.parent / 'ala.pdb'}")
                               .replace("system=ala.pdb", f"system={deck.parent / 'ala.pdb'}"))
            cwd = os.getcwd(); os.chdir(tmp)
            try:
                with mock.patch.object(QMMM_Opt, "_energy_force", flaky):
                    r = Runner(project="ef", input_file=str(inp), log=str(Path(tmp) / "ef.log"),
                               silent=1, usempi=False)
                    r.run()
                log = (Path(tmp) / "ef.log").read_text()
            finally:
                os.chdir(cwd)
        info = r.mol.qmmm_optimization
        self.assertTrue(info["electronic_failure"])
        self.assertTrue(info["recovery"])
        self.assertEqual(info["evaluations"], 2)          # the rejected trial geometry is not counted
        self.assertIn("electronic solver did not converge", log)
        self.assertIn("recovery selected after electronic non-convergence", log)


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class TestForceFieldPaths(unittest.TestCase):
    """[qmmm] forcefield_files: a deck-relative file whose name contains a
    space is one file, resolved before any splitting; lists split on commas
    or whitespace and each entry is resolved the same way."""

    def _opt(self, deck_dir):
        import types
        o = QMMM_Opt.__new__(QMMM_Opt)
        o.mol = types.SimpleNamespace(input_file=str(Path(deck_dir) / "x.inp"))
        return o

    def test_single_file_with_space_is_kept_whole(self):
        import os, tempfile
        with tempfile.TemporaryDirectory() as deck_dir, tempfile.TemporaryDirectory() as cwd:
            (Path(deck_dir) / "my forcefield.xml").write_text("<ForceField/>")
            here = os.getcwd(); os.chdir(cwd)
            try:
                o = self._opt(deck_dir)
                self.assertEqual(o._forcefield_paths("my forcefield.xml"),
                                 [str(Path(deck_dir) / "my forcefield.xml")])
                self.assertEqual(o._forcefield_paths("amber14-all.xml amber14/tip3p.xml"),
                                 ["amber14-all.xml", "amber14/tip3p.xml"])
                self.assertEqual(o._forcefield_paths("amber14-all.xml,amber14/tip3p.xml"),
                                 ["amber14-all.xml", "amber14/tip3p.xml"])
                self.assertEqual(o._forcefield_paths(""), [])
            finally:
                os.chdir(here)


class TestPublishedEnergyBelongsToThisRun(unittest.TestCase):
    """Molecule.get_results() publishes the QM/MM optimisation objective only
    for an optimisation; a summary left on a reused Molecule by an earlier
    optimisation must not replace a later run's energy, and Runner.run()
    clears it at the start of every calculation."""

    def _get_results(self):
        import ast
        from types import SimpleNamespace
        path = ROOT / "pyoqp" / "oqp" / "molecule" / "molecule.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "Molecule")
        fn = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "get_results")
        ns = {"np": np}
        exec(compile(ast.Module(body=[fn], type_ignores=[]), str(path), "exec"), ns)

        def fake(runtype, method="hf", energies=(-10.0,)):
            return SimpleNamespace(
                config={"input": {"runtype": runtype, "method": method}},
                mol_energy=SimpleNamespace(energy=-10.0), energies=list(energies),
                symmetry_metadata={}, data={"OQP::td_energies": [0.0]},
                get_atoms=lambda: np.array([1]), get_system=lambda: np.zeros(3),
                has_grad=lambda: False, get_grad=lambda: [], get_nac=lambda: [],
                get_soc=lambda: [], get_hess=lambda: [], get_mrsf_ekt_results=lambda: {},
                get_state_tracking=lambda: None, explicit_scf_props=lambda: [],
                qmmm_optimization={"converged": True, "energy_hartree": -12.5, "evaluations": 3,
                                   "recovery": False, "rms_grad": 1e-5, "max_grad": 2e-5,
                                   "output": "x_opt.pdb", "movable_atoms": [0, 1],
                                   "grad_fragment": [[0.1, 0.2, 0.3]]})
        return ns["get_results"], fake

    def test_optimisation_publishes_the_objective(self):
        get_results, fake = self._get_results()
        out = get_results(fake("optimize"))
        self.assertEqual(out["energy"], -12.5)
        self.assertEqual(out["qmmm_optimization"]["evaluations"], 3)
        self.assertNotIn("movable_atoms", out["qmmm_optimization"])
        self.assertEqual(out["grad"], [[0.1, 0.2, 0.3]])            # the objective's gradient
        self.assertNotIn("grad_fragment", out["qmmm_optimization"])

    def test_dftb_excited_state_optimisation_keeps_the_objective(self):
        # istate=1: mol.energies = [nan, objective]; the DFTB block would publish nan
        get_results, fake = self._get_results()
        out = get_results(fake("optimize", method="dftb", energies=(float("nan"), -12.5)))
        self.assertEqual(out["energy"], -12.5)

    def test_stale_summary_does_not_leak_into_another_runtype(self):
        get_results, fake = self._get_results()
        out = get_results(fake("energy"))
        self.assertEqual(out["energy"], -10.0)
        self.assertNotIn("qmmm_optimization", out)
        self.assertNotIn("grad", out)

    def test_runner_clears_the_summary_each_run(self):
        src = (ROOT / "pyoqp" / "oqp" / "pyoqp.py").read_text()
        run = src[src.index("    def run(self, test_mod=False):"):]
        run = run[:run.index('run_type = self.mol.config["input"]["runtype"]')]
        self.assertIn("self.mol.qmmm_optimization = None", run)


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class TestSelectionMatchesMolecule(unittest.TestCase):
    """The QM Molecule (from [input] system) must be the [qmmm] selection in
    topology order plus the link hydrogens."""

    def _opt(self, z_mol, qm, nlink):
        import types
        import openmm.app as app
        o = QMMM_Opt.__new__(QMMM_Opt)
        o.pdb = app.PDBFile(str(ROOT / "examples" / "QMMM" / "ala.pdb"))
        o.qm_atoms = np.array(qm, dtype=int)
        o.driver = types.SimpleNamespace(link_atoms=[object()] * nlink)
        o.mol = types.SimpleNamespace(get_atoms2=lambda what: np.array(z_mol, dtype=float))
        return o

    def _z(self, idx):
        import openmm.app as app
        atoms = list(app.PDBFile(str(ROOT / "examples" / "QMMM" / "ala.pdb")).topology.atoms())
        return [atoms[i].element.atomic_number for i in idx]

    def test_matching_layout_passes(self):
        qm = [8, 9, 16, 17, 18]
        self._opt(self._z(qm) + [1], qm, 1)._validate_qm_molecule_layout()

    def test_one_based_copy_of_the_selection_is_rejected(self):
        qm = [8, 9, 16, 17, 18]
        wrong = self._z([7, 8, 15, 16, 17]) + [1]           # ala.pdb has 19 atoms: shift down, not up
        if wrong == self._z(qm) + [1]:
            self.skipTest("shifted selection happens to have the same elements")
        with self.assertRaisesRegex(ValueError, "1-based"):
            self._opt(wrong, qm, 1)._validate_qm_molecule_layout()

    def test_wrong_atom_count_is_rejected(self):
        qm = [8, 9, 16, 17, 18]
        with self.assertRaisesRegex(ValueError, "does not match"):
            self._opt(self._z(qm), qm, 1)._validate_qm_molecule_layout()      # link H missing
        with self.assertRaisesRegex(ValueError, "beyond"):
            self._opt(self._z(qm) + [1], [8, 9, 16, 17, 100000], 1)._validate_qm_molecule_layout()


class TestRecoveryAndGradientAreWired(unittest.TestCase):
    def test_source(self):
        src = (ROOT / "pyoqp" / "oqp" / "library" / "qmmm_opt.py").read_text()
        for key in ("auto_recovery", "recovery_maxit", "recovery_trust"):
            self.assertIn(f'eng.get("{key}"', src)
        self.assertIn("mol.grads = ", src)
        self.assertIn("self._last_gradient = -f", src)
        self.assertIn("math.isfinite(self.radius)", src)


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class TestConstraintsOnMovableAtoms(unittest.TestCase):
    """[qmmm] rigidwater / constraints hold on the movable MM atoms."""

    QM = [8, 9, 16, 17, 18]

    def _opt(self, movable, cfg):
        import openmm.app as app
        o = QMMM_Opt.__new__(QMMM_Opt)
        o.pdb = app.PDBFile(str(ROOT / "examples" / "QMMM" / "ala.pdb"))
        o.forcefield = app.ForceField("amber14-all.xml")
        o.qm_atoms = np.array(self.QM, dtype=int)
        o.movable = np.array(sorted(movable), dtype=int)
        return o, o._constraint_pairs(cfg)

    def test_nothing_requested_nothing_held(self):
        o, pairs = self._opt(self.QM, {"constraints": "None", "rigidwater": False})
        self.assertEqual(pairs, [])

    def test_hbonds_exclude_qm_and_pull_in_partners(self):
        import openmm.app as app
        top = app.PDBFile(str(ROOT / "examples" / "QMMM" / "ala.pdb")).topology
        qm = set(self.QM)
        heavy, hyd = next((b[0].index, b[1].index) if b[1].element.symbol == "H" else (b[1].index, b[0].index)
                          for b in top.bonds()
                          if {b[0].element.symbol, b[1].element.symbol} >= {"H"} and len({b[0].element.symbol, b[1].element.symbol}) == 2
                          and b[0].index not in qm and b[1].index not in qm)
        o, pairs = self._opt(self.QM + [heavy], {"constraints": "HBonds", "rigidwater": False})
        self.assertIn(tuple(sorted((heavy, hyd))), [tuple(sorted(p)) for p in pairs])
        self.assertIn(hyd, o.movable)                                  # partner pulled in
        self.assertTrue(all(i not in qm and j not in qm for i, j in pairs))
        self.assertTrue(all(i in o.movable and j in o.movable for i, j in pairs))

    def test_constrained_group_split_across_the_cell_is_made_whole(self):
        import types
        L = 1.8                                                      # nm
        o = QMMM_Opt.__new__(QMMM_Opt)
        o.driver = types.SimpleNamespace(_box_lengths_bohr=lambda: [L / 0.052917721067] * 3)
        X = np.zeros((4, 3))
        X[0] = [1.75, 0.5, 0.5]                                     # O at the +x face
        X[1] = [1.75 + 0.09572 - L, 0.5, 0.5]                       # its H stored on the -x face
        X[2] = [1.75 - 0.024, 0.5 + 0.0927, 0.5]                    # the other H, whole
        X[3] = [0.3, 0.3, 0.3]                                      # unconstrained
        pairs = [(0, 1), (0, 2), (1, 2)]
        Y = o._unwrap_constrained(X, pairs)
        self.assertAlmostEqual(np.linalg.norm(Y[1] - Y[0]), 0.09572, places=9)
        np.testing.assert_allclose(Y[[0, 2, 3]], X[[0, 2, 3]], atol=0)
        o.driver = types.SimpleNamespace(_box_lengths_bohr=lambda: None)   # no cell: unchanged
        np.testing.assert_allclose(o._unwrap_constrained(X, pairs), X, atol=0)

    def test_unknown_constraint_name_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "HBonds"):
            self._opt(self.QM, {"constraints": "bonds-please", "rigidwater": False})


@unittest.skipUnless(_HAVE and _runtime_available(), "OpenMM or compiled OpenQP runtime unavailable")
class TestConstrainedAndReevaluatedResults(unittest.TestCase):
    def _run(self, text, patch_energy=None):
        import os, tempfile
        from unittest import mock
        from oqp.pyoqp import Runner
        deck = ROOT / "examples" / "QMMM" / "ala-dipeptide_RHF-QMMM-OPT-linkatom.inp"
        text = text.replace("pdb_file=ala.pdb", f"pdb_file={deck.parent / 'ala.pdb'}").replace(
            "system=ala.pdb", f"system={deck.parent / 'ala.pdb'}")
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "c.inp"
            inp.write_text(text)
            cwd = os.getcwd(); os.chdir(tmp)
            try:
                self._files = []
                ctx = mock.patch.object(QMMM_Opt, "_energy_force", patch_energy) if patch_energy else None
                if ctx:
                    ctx.start()
                try:
                    r = Runner(project="c", input_file=str(inp), log=str(Path(tmp) / "c.log"), silent=1, usempi=False)
                    r.run()
                    self._files = sorted(os.listdir(tmp))
                finally:
                    if ctx:
                        ctx.stop()
            finally:
                os.chdir(cwd)
        return r

    def _deck(self):
        deck = ROOT / "examples" / "QMMM" / "ala-dipeptide_RHF-QMMM-OPT-linkatom.inp"
        return deck.read_text().replace("save_mol=true", "save_mol=false")

    def test_hbond_distances_of_movable_mm_atoms_are_held(self):
        text = (self._deck().replace("maxit=12", "maxit=3").replace("qmmm_radius=0.0", "qmmm_radius=3.0\nqmmm_output=held.pdb")
                .replace("cutoff=NoCutoff", "cutoff=NoCutoff\nconstraints=HBonds"))
        r = self._run(text)
        o = r.qmmm_opt
        self.assertEqual(r.mol.qmmm_optimization["output"], "held.pdb")      # [optimize] qmmm_output honoured
        self.assertIn("held.pdb", self._files)
        self.assertGreater(len(o.frozen_pairs), 0)
        self.assertEqual(r.mol.qmmm_optimization["constraints"], len(o.frozen_pairs))
        import openmm.unit as unit
        X0 = np.array(o.pdb.positions.value_in_unit(unit.nanometer))
        X1 = np.asarray(o.positions_nm)
        moved = float(np.abs(X1[o.movable] - X0[o.movable]).max())
        self.assertGreater(moved, 1e-6)                                 # the search did move atoms
        for i, j in o.frozen_pairs:
            self.assertAlmostEqual(np.linalg.norm(X1[i] - X1[j]), np.linalg.norm(X0[i] - X0[j]), delta=1e-7)

    def test_rejected_final_trial_restores_the_reported_geometry(self):
        # auto_recovery off; the second trial moves mol, then the SCF 'fails':
        # the first geometry is reported, so mol must be brought back to it
        real = QMMM_Opt._energy_force
        calls = [0]

        def failing(self, X):
            calls[0] += 1
            if calls[0] == 2:
                real(self, X)                           # the driver/mol move to the trial
                raise RuntimeError("SCF did not converge in 200 iterations")
            return real(self, X)

        text = self._deck().replace("maxit=12", "maxit=3")
        r = self._run(text, patch_energy=failing)
        o, info = r.qmmm_opt, r.mol.qmmm_optimization
        self.assertTrue(info["electronic_failure"])
        self.assertEqual(calls[0], 3)                    # the re-evaluation at the reported geometry
        nqm = len(o.qm_atoms)
        mol_xyz = np.asarray(r.mol.get_system(), dtype=float).reshape(-1, 3)[:nqm]
        np.testing.assert_allclose(mol_xyz, o.history[-1]["x"].reshape(-1, 3), atol=1e-8)

    def test_same_coordinates_different_evaluation_is_reevaluated(self):
        # maxit=1 then a recovery restart from the best (= first) geometry: the
        # second evaluation is at identical coordinates but reports 1e-7 more,
        # so the first entry is reported while the driver holds the second;
        # the third (re-)evaluation returns E + 3e-7, which must be published
        real = QMMM_Opt._energy_force
        calls = [0]

        def shifted(self, X):
            calls[0] += 1
            e, f = real(self, X)
            return (e + 1e-7, f) if calls[0] == 2 else ((e + 3e-7, f) if calls[0] == 3 else (e, f))

        text = (self._deck().replace("maxit=12", "maxit=1")
                .replace("auto_recovery=false", "auto_recovery=true\nrecovery_maxit=1"))
        r = self._run(text, patch_energy=shifted)
        o, info = r.qmmm_opt, r.mol.qmmm_optimization
        self.assertTrue(info["recovery"])
        self.assertTrue(np.array_equal(o.history[0]["x"], o.history[1]["x"]))   # the same coordinates
        self.assertEqual(calls[0], 3)
        self.assertAlmostEqual(info["energy_hartree"] - o.history[0]["e"], 3e-7, delta=1e-9)

    def test_published_energy_is_the_reevaluated_one(self):
        # force 'best is not last': the second evaluation reports +1 Hartree, so
        # the first geometry is reported and re-evaluated; the re-evaluation
        # returns E + 5e-7, and that is what must be published
        real = QMMM_Opt._energy_force
        calls = [0]

        def shifted(self, X):
            calls[0] += 1
            e, f = real(self, X)
            if calls[0] == 2:
                return e + 1.0, f
            if calls[0] == 3:
                return e + 5e-7, f
            return e, f

        text = self._deck().replace("maxit=12", "maxit=2")
        r = self._run(text, patch_energy=shifted)
        info, hist = r.mol.qmmm_optimization, r.qmmm_opt.history
        self.assertEqual(calls[0], 3)
        self.assertAlmostEqual(info["energy_hartree"] - hist[0]["e"], 5e-7, delta=1e-9)
        self.assertEqual(r.mol.energies[0], info["energy_hartree"])


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class TestVirtualSites(unittest.TestCase):
    """TIP4P-Ew water: the M site (element None) is not an optimisation
    coordinate, and its position is rebuilt from the moved O and H atoms."""

    def _water4(self):
        import openmm as mm
        import openmm.app as app
        pdb = app.PDBFile(str(ROOT / "examples" / "QMMM" / "formaldehyde_water.pdb"))
        mod = app.Modeller(pdb.topology, pdb.positions)
        mod.delete([r for r in mod.topology.residues() if r.name != "HOH"])
        ff = app.ForceField("amber14/tip4pew.xml")
        mod.addExtraParticles(ff)
        system = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff, rigidWater=False)
        sim = app.Simulation(mod.topology, system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
        return mod, system, sim

    def _opt(self, mod, system, sim, qm):
        import types
        o = QMMM_Opt.__new__(QMMM_Opt)
        o.pdb = types.SimpleNamespace(topology=mod.topology, positions=mod.positions)
        o.qm_atoms = np.array(qm, dtype=int)
        o.driver = types.SimpleNamespace(mm_systems={"sys0": system, "sim0": sim},
                                         _box_lengths_bohr=lambda: None)
        return o

    def test_movable_set_has_no_virtual_site(self):
        mod, system, sim = self._water4()
        atoms = list(mod.topology.atoms())
        self.assertTrue(any(a.element is None for a in atoms))
        qm = [a.index for a in list(mod.topology.residues())[0].atoms() if a.element is not None]
        o = self._opt(mod, system, sim, qm)
        mv = o._movable_atoms(20.0)                       # every water within reach
        self.assertGreater(len(mv), len(qm))
        self.assertTrue(all(atoms[i].element is not None for i in mv))
        self.assertEqual([int(atoms[i].element.atomic_number) for i in mv][:1], [8])   # symbols build

    def test_virtual_site_follows_its_parents(self):
        import openmm.unit as unit
        mod, system, sim = self._water4()
        atoms = list(mod.topology.atoms())
        o = self._opt(mod, system, sim, [0, 1, 2])
        X = np.array(mod.positions.value_in_unit(unit.nanometer))
        ep = next(a.index for a in atoms if a.element is None)
        parents = [a.index for a in atoms[ep].residue.atoms() if a.element is not None]
        X2 = X.copy(); X2[parents] += np.array([0.05, -0.02, 0.01])   # translate the real atoms
        Y = o._with_virtual_sites(X2)
        np.testing.assert_allclose(Y[ep] - X[ep], [0.05, -0.02, 0.01], atol=1e-9)
        np.testing.assert_allclose(Y[parents], X2[parents], atol=0)


class TestStateAndCoordinatesForQmmmOptimisation(unittest.TestCase):
    def _report(self, method="hf", istate=0, coordsys="auto"):
        from oqp.utils import input_checker as chk
        cfg = {"input": {"runtype": "optimize", "qmmm_flag": True, "method": method, "basis": "6-31g",
                         "system": "ala.pdb 9 10 17 18 19", "charge": 0},
               "optimize": {"lib": "oqp", "istate": istate}, "oqp": {"coordsys": coordsys},
               "qmmm": {"pdb_file": "ala.pdb", "qm_atoms": "8,9,16,17,18", "forcefield_files": "amber14-all.xml"}}
        report = chk.CheckReport()
        chk._check_optimize(cfg, report)
        return [(d.severity, d.path) for d in report.diagnostics]

    def test_checker_rejects_negative_or_ground_tdhf_state(self):
        self.assertIn(("ERROR", "optimize.istate"), self._report(istate=-1))
        self.assertIn(("ERROR", "optimize.istate"), self._report(method="tdhf", istate=0))
        self.assertIn(("ERROR", "optimize.istate"), self._report(method="tdhf", istate=-1))
        self.assertNotIn(("ERROR", "optimize.istate"), self._report(istate=0))
        self.assertNotIn(("ERROR", "optimize.istate"), self._report(method="tdhf", istate=2))

    def test_checker_and_driver_reject_swapmo(self):
        from oqp.utils import input_checker as chk
        cfg = {"input": {"runtype": "optimize", "qmmm_flag": True, "method": "hf", "basis": "6-31g",
                         "system": "ala.pdb 9 10 17 18 19", "charge": 0},
               "optimize": {"lib": "oqp", "istate": 0}, "guess": {"swapmo": "5 6"},
               "qmmm": {"pdb_file": "ala.pdb", "qm_atoms": "8,9,16,17,18", "forcefield_files": "amber14-all.xml"}}
        report = chk.CheckReport()
        chk._check_optimize(cfg, report)
        self.assertIn(("ERROR", "guess.swapmo"), [(d.severity, d.path) for d in report.diagnostics])
        cfg["guess"]["swapmo"] = ""
        report = chk.CheckReport()
        chk._check_optimize(cfg, report)
        self.assertNotIn(("ERROR", "guess.swapmo"), [(d.severity, d.path) for d in report.diagnostics])
        if _HAVE:
            for bad in ("5 6", [5, 6]):
                with self.assertRaisesRegex(ValueError, "swapmo"):
                    QMMM_Opt._reject_swapmo(bad)
            for ok in ("", [], None):
                QMMM_Opt._reject_swapmo(ok)

    def test_checker_rejects_dlc_and_ric(self):
        for cs in ("dlc", "ric", "internal"):
            self.assertIn(("ERROR", "oqp.coordsys"), self._report(coordsys=cs), cs)
        for cs in ("auto", "cartesian", "cart", "tric"):
            self.assertNotIn(("ERROR", "oqp.coordsys"), self._report(coordsys=cs), cs)

    @unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
    def test_driver_guards(self):
        self.assertEqual(QMMM_Opt._resolve_coordsys("auto"), "cartesian")
        self.assertEqual(QMMM_Opt._resolve_coordsys("TRIC"), "tric")
        for cs in ("dlc", "ric"):
            with self.assertRaisesRegex(ValueError, "translations and rotations"):
                QMMM_Opt._resolve_coordsys(cs)
        QMMM_Opt._validate_istate(0, "hf")
        QMMM_Opt._validate_istate(1, "tdhf")
        for istate, method in ((-1, "hf"), (0, "tdhf"), (-2, "tdhf")):
            with self.assertRaises(ValueError):
                QMMM_Opt._validate_istate(istate, method)


@unittest.skipUnless(_HAVE and _runtime_available(), "OpenMM or compiled OpenQP runtime unavailable")
class TestTip4pBoxOptimisation(unittest.TestCase):
    """Formaldehyde in five TIP4P-Ew waters, whose M sites have no element:
    the QM/MM driver builds (its link-atom scan used to dereference every
    atom's element), the optimisation with movable waters runs, and each M
    site it reports sits where OpenMM places it from the moved O and H."""

    def test_optimisation_with_movable_tip4p_waters(self):
        import os, shutil, tempfile
        import openmm as mm
        import openmm.app as app
        import openmm.unit as unit
        from oqp.pyoqp import Runner
        ex = ROOT / "examples" / "QMMM"
        pdb = app.PDBFile(str(ex / "formaldehyde_water.pdb"))
        ff = app.ForceField(str(ex / "formaldehyde.xml"), "amber14/tip4pew.xml")
        mod = app.Modeller(pdb.topology, pdb.positions)
        mod.addExtraParticles(ff)
        deck = """[input]
system=box4.pdb 1 2 3 4
charge=0
runtype=optimize
basis=sto-3g
method=hf
qmmm_flag=True
[scf]
type=rhf
multiplicity=1
[optimize]
istate=0
maxit=2
qmmm_radius=3.0
[oqp]
auto_recovery=false
[qmmm]
pdb_file=box4.pdb
forcefield_files=formaldehyde.xml amber14/tip4pew.xml
qm_atoms=0-3
cutoff=NoCutoff
embedding=electrostatic
"""
        with tempfile.TemporaryDirectory() as tmp:
            with open(Path(tmp) / "box4.pdb", "w") as fh:
                app.PDBFile.writeFile(mod.topology, mod.positions, fh, keepIds=True)
            shutil.copy(ex / "formaldehyde.xml", Path(tmp) / "formaldehyde.xml")
            (Path(tmp) / "opt4.inp").write_text(deck)
            cwd = os.getcwd(); os.chdir(tmp)
            try:
                r = Runner(project="opt4", input_file="opt4.inp", log="opt4.log", silent=1, usempi=False)
                r.run()
            finally:
                os.chdir(cwd)
        o = r.qmmm_opt
        atoms = list(o.pdb.topology.atoms())
        ep = [a.index for a in atoms if a.element is None]
        self.assertEqual(len(ep), 5)
        self.assertTrue(np.isfinite(r.mol.energies[0]))
        mv = set(r.mol.qmmm_optimization["movable_atoms"])
        self.assertFalse(mv & set(ep))                                   # M sites are not coordinates
        moved_waters = [a.residue for a in atoms if a.index in mv and a.residue.name == "HOH"]
        self.assertTrue(moved_waters)                                    # the 3 A shell reaches some water
        # every reported M site is OpenMM's placement from the reported O/H
        system = ff.createSystem(o.pdb.topology, nonbondedMethod=app.NoCutoff, rigidWater=False)
        ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
        X = np.asarray(o.positions_nm)
        ctx.setPositions(unit.Quantity(X, unit.nanometer))
        ctx.computeVirtualSites()
        P = np.asarray(ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        np.testing.assert_allclose(X[ep], P[ep], atol=1e-9)


class TestNativeControlsCheckedForQmmmOptimisation(unittest.TestCase):
    def test_recovery_controls_validated_whatever_lib_says(self):
        from oqp.utils import input_checker as chk
        for lib in ("oqp", "geometric", "scipy"):
            cfg = {"input": {"runtype": "optimize", "qmmm_flag": True, "method": "hf", "basis": "6-31g",
                             "system": "ala.pdb 9 10 17 18 19", "charge": 0},
                   "optimize": {"lib": lib, "istate": 0}, "oqp": {"recovery_trust": -1.0},
                   "qmmm": {"pdb_file": "ala.pdb", "qm_atoms": "8,9,16,17,18", "forcefield_files": "amber14-all.xml"}}
            report = chk.CheckReport()
            chk._check_optimize(cfg, report)
            self.assertIn(("ERROR", "oqp.recovery_trust"), [(d.severity, d.path) for d in report.diagnostics], lib)


if __name__ == "__main__":
    unittest.main()
