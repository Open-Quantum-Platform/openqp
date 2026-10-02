"""Checkpoint/restart of ground-state QM/MM MD and the snapshot hand-over to
later runs (QM/MM MD and QM/MM surface hopping).

The protocol these serve: equilibrate (classical NPT, then QM/MM NVT), take
phase-space snapshots along the equilibrated trajectory, start one
surface-hopping trajectory from each.  Before this, QM/MM MD wrote positions
only and could not be continued, and QM/MM NAMD always drew fresh
Maxwell-Boltzmann velocities, so no trajectory could inherit an equilibrated
state.
"""
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "QMMM"

_SPEC = importlib.util.spec_from_file_location(
    "_md_snapshot_direct", ROOT / "pyoqp" / "oqp" / "utils" / "md_snapshot.py")
md_snapshot = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = md_snapshot
_SPEC.loader.exec_module(md_snapshot)


# ---------------------------------------------------------------- format ---

def _state(natom=4):
    rng = np.random.default_rng(5)
    return dict(positions_nm=rng.normal(size=(natom, 3)),
                velocities_nm_ps=rng.normal(size=(natom, 3)),
                masses_dalton=np.array([12.011, 15.999, 1.008, 1.008])[:natom])


def test_snapshot_round_trip_and_optional_fields(tmp_path):
    path = tmp_path / "s.npz"
    state = _state()
    md_snapshot.write_snapshot(path, **state)
    back = md_snapshot.read_snapshot(path)
    np.testing.assert_array_equal(back["positions_nm"], state["positions_nm"])
    np.testing.assert_array_equal(back["velocities_nm_ps"], state["velocities_nm_ps"])
    assert back["box_nm"] is None and back["integrator_velocities_nm_ps"] is None
    assert back["step"] == 0 and back["natom"] == 4

    leapfrog = state["velocities_nm_ps"] - 0.1
    md_snapshot.write_snapshot(path, box_nm=[3.0, 3.1, 3.2], step=40, time_ps=0.02,
                               integrator_velocities_nm_ps=leapfrog,
                               qm_atoms=[3, 0], **state)
    back = md_snapshot.read_snapshot(path)
    np.testing.assert_allclose(back["box_nm"], [3.0, 3.1, 3.2])
    np.testing.assert_array_equal(back["integrator_velocities_nm_ps"], leapfrog)
    assert back["step"] == 40 and list(back["qm_atoms"]) == [0, 3]
    assert not list(tmp_path.glob(".snapshot-*"))          # no temporary left behind


def test_snapshot_rejects_bad_input_and_foreign_files(tmp_path):
    state = _state()
    with pytest.raises(ValueError, match="velocities must have shape"):
        md_snapshot.write_snapshot(tmp_path / "a.npz", **{
            **state, "velocities_nm_ps": np.zeros((3, 3))})
    with pytest.raises(ValueError, match="NaN or infinity"):
        bad = dict(state)
        bad["positions_nm"] = state["positions_nm"].copy()
        bad["positions_nm"][0, 0] = np.nan
        md_snapshot.write_snapshot(tmp_path / "a.npz", **bad)
    assert not (tmp_path / "a.npz").exists()                # nothing half-written

    with pytest.raises(FileNotFoundError):
        md_snapshot.read_snapshot(tmp_path / "missing.npz")
    np.savez(tmp_path / "other.npz", coordinates=np.zeros(3))  # e.g. a NAMD checkpoint
    with pytest.raises(ValueError, match="not an OpenQP MD snapshot"):
        md_snapshot.read_snapshot(tmp_path / "other.npz")


def test_snapshot_of_another_system_is_refused(tmp_path):
    path = tmp_path / "s.npz"
    state = _state()
    md_snapshot.write_snapshot(path, **state)
    snap = md_snapshot.read_snapshot(path)
    md_snapshot.check_snapshot_matches(snap, state["masses_dalton"])
    # held atoms / virtual sites are massless on one side: compared by count only
    md_snapshot.check_snapshot_matches(snap, [12.011, 0.0, 1.008, 1.008])
    with pytest.raises(ValueError, match="has 4 atoms, the system has 3"):
        md_snapshot.check_snapshot_matches(snap, state["masses_dalton"][:3])
    with pytest.raises(ValueError, match="different topology or atom order"):
        md_snapshot.check_snapshot_matches(snap, state["masses_dalton"][::-1])


def test_numbered_snapshot_names():
    f = md_snapshot.numbered_snapshot_path
    assert f("run.restart.npz", 100) == "run.snapshot.00000100.npz"
    assert f("qmmm_md.restart.npz", 7) == "qmmm_md.snapshot.00000007.npz"
    assert f("/a/b/state.npz", 12) == "/a/b/state.snapshot.00000012.npz"


# ------------------------------------------------------------- the driver ---

DECK = """[input]
system = {pdb} 1 2 3 4
functional = bhhlyp
basis = 6-31g
method = hf
runtype = md
qmmm_flag = true

[scf]
type = rhf
conv = 1e-9

[qmmm]
pdb_file = {pdb}
forcefield_files = {ff} {tip}
qm_atoms = 0-3
cutoff = NoCutoff
n_steps = {nsteps}
timestep = 0.5
ensemble = {ensemble}
trajectory_format = {fmt}
trajectory_file = {name}.{fmt}
log_file = {name}.dat
energy_file = {name}.npz

[md]
restart_file = {name}.restart.npz
{md}
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        import openmm  # noqa: F401
        from oqp.library.qmmm_md import QMMM_MD  # noqa: F401
        return (EXAMPLES / "formaldehyde_water.pdb").exists()
    except Exception:
        return False


def _run(name, nsteps, md="velocity = zero\ncommon_control_keys = velocity",
         ensemble="nve", fmt="pdb"):
    from oqp.library.qmmm_md import QMMM_MD
    deck = Path(f"{name}.{nsteps}.inp")
    deck.write_text(DECK.format(pdb=EXAMPLES / "formaldehyde_water.pdb",
                                ff=EXAMPLES / "formaldehyde.xml",
                                tip=EXAMPLES / "tip3p.xml", nsteps=nsteps,
                                ensemble=ensemble, fmt=fmt, name=name, md=md))
    driver = QMMM_MD(oqp_cfg=str(deck))
    return driver, driver.run()


@unittest.skipUnless(_runtime_available(), "OpenQP runtime with OpenMM required")
class TestRestartAndSnapshots(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._cwd = os.getcwd()
        os.chdir(self._tmp.name)

    def tearDown(self):
        os.chdir(self._cwd)
        self._tmp.cleanup()

    def test_a_restarted_nve_run_is_the_uninterrupted_run(self):
        """Stop at step 2, continue to step 4: energies, log rows and frames
        equal those of one 4-step run, with nothing repeated."""
        import openmm.app as app
        import openmm.unit as u
        full = _run("full", 4)[1]
        _run("part", 2)
        zero = "velocity = zero\ncommon_control_keys = velocity"
        cont = _run("part", 4, md=zero + "\nrestart = true")[1]

        self.assertEqual(list(cont["step"]), [0, 1, 2, 3, 4])
        np.testing.assert_allclose(cont["E_tot"], full["E_tot"], rtol=0, atol=1e-4)
        np.testing.assert_allclose(cont["E_kin"], full["E_kin"], rtol=0, atol=1e-4)
        rows = np.loadtxt("part.dat", delimiter=",", skiprows=1)
        np.testing.assert_allclose(rows[:, 0], [0.0, 0.0005, 0.001, 0.0015, 0.002])
        with open("part.dat") as stream:
            self.assertEqual(sum(line.startswith("#") for line in stream), 1)
        a, b = app.PDBFile("part.pdb"), app.PDBFile("full.pdb")
        self.assertEqual(a.getNumFrames(), 5)
        for frame in range(5):
            np.testing.assert_allclose(
                np.array(a.getPositions(frame=frame).value_in_unit(u.nanometer)),
                np.array(b.getPositions(frame=frame).value_in_unit(u.nanometer)),
                atol=2e-4)                                   # PDB precision
        final = md_snapshot.read_snapshot("part.restart.npz")
        self.assertEqual(final["step"], 4)

    def test_restart_drops_output_written_after_the_checkpoint(self):
        """A run that died between checkpoints left rows and frames past the
        checkpoint; continuing must not keep them."""
        import openmm.app as app
        zero = "velocity = zero\ncommon_control_keys = velocity"
        full = _run("crash", 3, md=zero + "\nrestart_interval = 2")[1]
        # put back the checkpoint of step 2, as if the run had died at step 3
        # (the end-of-run checkpoint replaced it)
        _run("ref", 2)
        os.replace("ref.restart.npz", "crash.restart.npz")
        cont = _run("crash", 4, md=zero + "\nrestart = true")[1]
        self.assertEqual(list(cont["step"]), [0, 1, 2, 3, 4])
        rows = np.loadtxt("crash.dat", delimiter=",", skiprows=1)
        np.testing.assert_allclose(rows[:, 0], [0.0, 0.0005, 0.001, 0.0015, 0.002])
        self.assertEqual(app.PDBFile("crash.pdb").getNumFrames(), 5)
        np.testing.assert_allclose(cont["E_tot"][:3], full["E_tot"][:3], atol=1e-4)

    @staticmethod
    def _frame_steps(path):
        with open(path) as stream:
            return [int(line.split()[3]) for line in stream
                    if line.startswith("REMARK OPENQP STEP")]

    def test_a_trajectory_started_at_a_checkpoint_is_resumed_by_step(self):
        """A restart whose trajectory file is missing starts the file at the
        checkpoint, not at step 0.  A later restart counted frames as if the
        file began at step 0, kept a frame written after its checkpoint and
        then wrote that step a second time."""
        import openmm.app as app
        zero = "velocity = zero\ncommon_control_keys = velocity"
        _run("seg", 2)
        self.assertEqual(self._frame_steps("seg.pdb"), [0, 1, 2])
        os.remove("seg.pdb")
        _run("seg", 5, md=zero + "\nrestart = true")
        self.assertEqual(self._frame_steps("seg.pdb"), [2, 3, 4, 5])
        # as if that run had died after step 5 with its last checkpoint at 4
        _run("ref", 4)
        os.replace("ref.restart.npz", "seg.restart.npz")
        _run("seg", 6, md=zero + "\nrestart = true")
        self.assertEqual(self._frame_steps("seg.pdb"), [2, 3, 4, 5, 6])
        self.assertEqual(app.PDBFile("seg.pdb").getNumFrames(), 5)

        # a checkpoint older than every frame on file: nothing to append to
        _run("old", 1)
        os.replace("old.restart.npz", "seg.restart.npz")
        _run("seg", 3, md=zero + "\nrestart = true")
        self.assertEqual(self._frame_steps("seg.pdb"), [1, 2, 3])

    def test_a_dcd_started_at_a_checkpoint_is_resumed_by_step(self):
        zero = "velocity = zero\ncommon_control_keys = velocity"

        def header():
            with open("seg.dcd", "rb") as stream:
                return [int(v) for v in np.frombuffer(stream.read(20)[8:20], dtype="<i4")]

        _run("seg", 2, fmt="dcd")
        os.remove("seg.dcd")
        _run("seg", 5, md=zero + "\nrestart = true", fmt="dcd")
        self.assertEqual(header(), [4, 2, 1])            # frames 2..5, first step 2
        _run("seg", 6, md=zero + "\nrestart = true", fmt="dcd")
        self.assertEqual(header()[0], 5)                 # 2..6, nothing repeated
        # checkpoint at step 4 while the file runs to step 6: refuse, a DCD
        # cannot be shortened (this used to append a second frame for step 5)
        _run("ref", 4)
        os.replace("ref.restart.npz", "seg.restart.npz")
        before = Path("seg.dcd").read_bytes()
        with self.assertRaisesRegex(ValueError, "starts at step 2 and holds 5 frames"):
            _run("seg", 7, md=zero + "\nrestart = true", fmt="dcd")
        self.assertEqual(Path("seg.dcd").read_bytes(), before)

    def test_energy_table_lands_in_the_file_a_restart_reloads(self):
        """np.savez renames by a case-sensitive suffix rule; "e.NPZ" went to
        "e.NPZ.npz", the restart looked for "e.NPZ", found nothing and wrote a
        table holding only the continued segment."""
        from oqp.library.qmmm_md import QMMM_MD
        zero = "velocity = zero\ncommon_control_keys = velocity"
        for given, written in (("e.NPZ", "e.NPZ"), ("f.Npz", "f.Npz"),
                               ("g", "g.npz"), ("h.npy", "h.npz"), ("i.dat2", "i.dat2.npz")):
            name = given.split(".")[0]
            for nsteps, md in ((2, zero), (4, zero + "\nrestart = true")):
                deck = Path(f"{name}.{nsteps}.inp")
                deck.write_text(DECK.format(
                    pdb=EXAMPLES / "formaldehyde_water.pdb",
                    ff=EXAMPLES / "formaldehyde.xml", tip=EXAMPLES / "tip3p.xml",
                    nsteps=nsteps, ensemble="nve", fmt="pdb", name=name,
                    md=md).replace(f"energy_file = {name}.npz", f"energy_file = {given}"))
                driver = QMMM_MD(oqp_cfg=str(deck))
                self.assertEqual(driver._energy_npz_path(), written)
                driver.run()
            self.assertFalse(Path(written + ".npz").exists(), given)
            with np.load(written) as table:
                self.assertEqual(list(table["step"]), [0, 1, 2, 3, 4], given)

    def test_restart_with_another_time_step_keeps_the_physical_velocities(self):
        """The integrator's velocities sit half a step behind the positions,
        by half of the time step they were written with.  Installed unchanged
        under another time step they are not the velocities of the checkpoint."""
        from oqp.library.qmmm_md import QMMM_MD
        import openmm.unit as u
        zero = "velocity = zero\ncommon_control_keys = velocity"
        _run("dt", 4)
        saved = md_snapshot.read_snapshot("dt.restart.npz")
        self.assertAlmostEqual(saved["timestep_ps"], 0.0005)
        self.assertGreater(np.abs(saved["velocities_nm_ps"]).max(), 1e-3)

        def restarted(dt_fs):
            deck = Path(f"dt.{dt_fs}.inp")
            deck.write_text(DECK.format(
                pdb=EXAMPLES / "formaldehyde_water.pdb",
                ff=EXAMPLES / "formaldehyde.xml", tip=EXAMPLES / "tip3p.xml",
                nsteps=6, ensemble="nve", fmt="pdb", name="dt",
                md=zero + "\nrestart = true").replace(
                    "timestep = 0.5", f"timestep = {dt_fs}"))
            driver = QMMM_MD(oqp_cfg=str(deck))
            driver.setup()
            context = driver.simulation_md.context
            leapfrog = np.asarray(context.getState(getVelocities=True).getVelocities(
                asNumpy=True).value_in_unit(u.nanometer / u.picosecond))
            masses = np.asarray(saved["masses_dalton"])
            moving = masses > 0.0
            onstep = leapfrog.copy()
            onstep[moving] += (0.5 * dt_fs * 1e-3 * driver._force_kj_nm[moving]
                               / masses[moving, None])
            return leapfrog, onstep

        same_leapfrog, _ = restarted(0.5)                     # exact continuation
        np.testing.assert_array_equal(same_leapfrog, saved["integrator_velocities_nm_ps"])
        for dt_fs in (0.25, 1.0):
            leapfrog, onstep = restarted(dt_fs)
            # the velocities AT the checkpoint positions are the saved ones
            # (QM atoms: the rigid waters also get a constraint projection)
            np.testing.assert_allclose(
                onstep[:4], saved["velocities_nm_ps"][:4], atol=1e-6)
            # ... which the raw integrator velocities of the old step are not
            self.assertGreater(
                np.abs(leapfrog - saved["integrator_velocities_nm_ps"]).max(), 1e-4)

    def test_restart_validation(self):
        zero = "velocity = zero\ncommon_control_keys = velocity"
        with self.assertRaises(FileNotFoundError):
            _run("none", 2, md=zero + "\nrestart = true")
        _run("done", 2)
        with self.assertRaisesRegex(ValueError, "nstep is the total length"):
            _run("done", 2, md=zero + "\nrestart = true")
        before = Path("done.pdb").read_text()
        with self.assertRaisesRegex(ValueError, "give one of them, not both"):
            _run("done", 4, md="restart = true\nsnapshot = done.restart.npz")
        self.assertEqual(Path("done.pdb").read_text(), before)   # untouched

    def test_a_restarted_driver_can_be_run_again(self):
        """run() completes the continued PDB (END record) and closes it; a
        second run() on the same driver must reopen it, not write to a closed
        file."""
        import openmm.app as app
        zero = "velocity = zero\ncommon_control_keys = velocity"
        _run("again", 2)
        driver, _ = _run("again", 3, md=zero + "\nrestart = true")
        self.assertEqual(app.PDBFile("again.pdb").getNumFrames(), 4)
        driver.n_steps = 1
        data = driver.run()                      # one more step, same driver
        self.assertEqual(list(data["step"]), [0, 1, 2, 3, 4])
        self.assertEqual(app.PDBFile("again.pdb").getNumFrames(), 5)
        with open("again.pdb") as stream:
            text = stream.read()
        self.assertEqual(text.count("\nEND\n"), 1)
        del driver

    def test_a_run_started_from_a_snapshot_follows_the_run_that_wrote_it(self):
        """NVE: a snapshot taken at step 2 and used as the start of a new run
        must reproduce steps 3 and 4 of the run it came from.

        The snapshot stores velocities AT its positions, while the leapfrog
        integrator wants them half a step behind; handing them over unshifted
        adds half a kick to the first step and the two trajectories separate
        at once."""
        zero = "velocity = zero\ncommon_control_keys = velocity"
        source = _run("src", 4, md=zero + "\nsnapshot_interval = 2")[1]
        child = _run("child", 2, md="snapshot = src.snapshot.00000002.npz")[1]
        self.assertEqual(list(child["step"]), [0, 1, 2])            # a new run
        np.testing.assert_allclose(child["E_pot"], source["E_pot"][2:], rtol=0, atol=2e-3)
        np.testing.assert_allclose(child["E_kin"], source["E_kin"][2:], rtol=0, atol=2e-3)

    def test_restart_finds_the_energy_table_under_the_name_numpy_gave_it(self):
        """energy_file without .npz is written as <name>.npz by np.savez; the
        restart must reload that file, or the rows before the checkpoint are
        lost when the continued run rewrites it."""
        from oqp.library.qmmm_md import QMMM_MD
        zero = "velocity = zero\ncommon_control_keys = velocity"

        def run(nsteps, extra=""):
            deck = Path(f"bare.{nsteps}.inp")
            deck.write_text(DECK.format(
                pdb=EXAMPLES / "formaldehyde_water.pdb", ff=EXAMPLES / "formaldehyde.xml",
                tip=EXAMPLES / "tip3p.xml", nsteps=nsteps, ensemble="nve", fmt="pdb",
                name="bare", md=zero + extra).replace("energy_file = bare.npz",
                                                      "energy_file = energies"))
            return QMMM_MD(oqp_cfg=str(deck)).run()

        first = run(2)
        self.assertTrue(os.path.isfile("energies.npz"))
        self.assertFalse(os.path.isfile("energies"))
        data = run(4, "\nrestart = true")
        self.assertEqual(list(data["step"]), [0, 1, 2, 3, 4])
        np.testing.assert_allclose(data["E_tot"][:3], first["E_tot"], rtol=0, atol=1e-9)
        with np.load("energies.npz") as saved:
            self.assertEqual(list(saved["step"]), [0, 1, 2, 3, 4])

    def test_a_driver_kept_alive_cannot_damage_the_continued_trajectory(self):
        """A script that restarts a run usually still holds the first driver.

        OpenMM's PDBReporter writes its END record when it is garbage
        collected, at its own old file offset -- into the middle of the file
        the restarted driver has appended to.  The driver's own reporter
        completes and closes the file at the end of run() instead."""
        import gc
        import openmm.app as app
        zero = "velocity = zero\ncommon_control_keys = velocity"
        first_driver, _ = _run("keep", 2)                  # kept alive on purpose
        with open("keep.pdb") as stream:
            self.assertEqual(stream.read().count("\nEND\n"), 1)   # complete after run()
        self.assertFalse(any(type(r).__name__ == "PDBReporter"
                             for r in first_driver.simulation_md.reporters))
        second_driver, data = _run("keep", 4, md=zero + "\nrestart = true")
        self.assertEqual(list(data["step"]), [0, 1, 2, 3, 4])
        continued = Path("keep.pdb").read_text()
        del first_driver                                   # now it is collected
        gc.collect()
        self.assertEqual(Path("keep.pdb").read_text(), continued)   # untouched
        self.assertEqual(app.PDBFile("keep.pdb").getNumFrames(), 5)
        self.assertEqual(continued.count("\nEND\n"), 1)
        del second_driver

    def test_velocities_given_at_the_start_are_on_position_velocities(self):
        """velocity=zero means zero kinetic energy at step 0, and a velocity
        file means exactly its kinetic energy.

        The integrator keeps velocities half a step behind the positions, so
        installing the given ones unshifted added half a kick: a run started
        "at rest" reported a nonzero kinetic energy at step 0."""
        zero = "velocity = zero\ncommon_control_keys = velocity"
        data = _run("rest", 1, md=zero)[1]
        self.assertAlmostEqual(float(data["E_kin"][0]), 0.0, places=9)

        from oqp.library.qmmm_md import AU_VELOCITY_TO_NM_PER_PS
        rng = np.random.default_rng(3)
        v_au = np.zeros((19, 3))
        v_au[:4] = rng.normal(scale=2.0e-4, size=(4, 3))      # QM atoms only: no constraints
        np.savetxt("v.dat", v_au)
        driver, data = _run("file", 1, md="velocity = v.dat\ncommon_control_keys = velocity")
        masses = driver._sys0_masses
        expected = 0.5 * np.sum(masses[:, None] * (v_au * AU_VELOCITY_TO_NM_PER_PS) ** 2)
        self.assertAlmostEqual(float(data["E_kin"][0]), expected, delta=1e-6 * expected)
        del driver

    def test_trajectory_interval_sets_the_output_cadence(self):
        import openmm.app as app
        md = ("velocity = zero\ntrajectory_interval = 2\n"
              "common_control_keys = velocity,trajectory_interval")
        data = _run("cad", 4, md=md)[1]
        self.assertEqual(list(data["step"]), [0, 1, 2, 3, 4])       # table: every step
        self.assertEqual(app.PDBFile("cad.pdb").getNumFrames(), 3)    # frames 0, 2, 4
        rows = np.loadtxt("cad.dat", delimiter=",", skiprows=1)
        np.testing.assert_allclose(rows[:, 0], [0.0, 0.001, 0.002])

    def test_a_snapshot_is_never_overwritten_by_the_run_it_starts(self):
        from oqp.library.qmmm_md import QMMM_MD
        _run("orig", 2, md="velocity = zero\ncommon_control_keys = velocity\n"
                           "snapshot_interval = 2")
        before = Path("orig.snapshot.00000002.npz").read_bytes()

        def deck(name, md):
            path = Path(f"{name}.inp")
            path.write_text(DECK.format(
                pdb=EXAMPLES / "formaldehyde_water.pdb", ff=EXAMPLES / "formaldehyde.xml",
                tip=EXAMPLES / "tip3p.xml", nsteps=2, ensemble="nve", fmt="pdb",
                name=name, md=md))
            return str(path)

        # the checkpoint destination is the snapshot itself
        with self.assertRaisesRegex(ValueError, "restart_file=.* would overwrite the input"):
            QMMM_MD(oqp_cfg=deck("orig", "snapshot = orig.restart.npz"))
        # a numbered snapshot of this run's own family
        with self.assertRaisesRegex(ValueError, "would overwrite the starting point"):
            QMMM_MD(oqp_cfg=deck("orig", "snapshot = orig.snapshot.00000002.npz\n"
                                         "snapshot_interval = 2"))
        self.assertEqual(Path("orig.snapshot.00000002.npz").read_bytes(), before)
        # a different run name is fine, and the input survives the whole run
        QMMM_MD(oqp_cfg=deck("next", "snapshot = orig.snapshot.00000002.npz\n"
                                     "snapshot_interval = 2")).run()
        self.assertEqual(Path("orig.snapshot.00000002.npz").read_bytes(), before)
        self.assertTrue(os.path.isfile("next.snapshot.00000002.npz"))

    def test_no_output_may_land_on_an_input_or_on_another_output(self):
        """Every destination is checked, not only the checkpoint: an energy
        table, a trajectory or a log named like the snapshot would replace or
        truncate the starting point of every later trajectory."""
        from oqp.library.qmmm_md import QMMM_MD
        _run("src", 2, md="velocity = zero\ncommon_control_keys = velocity\n"
                          "snapshot_interval = 2")
        snapshot = Path("src.snapshot.00000002.npz")
        before = snapshot.read_bytes()
        pdb = EXAMPLES / "formaldehyde_water.pdb"
        pdb_before = pdb.read_bytes()

        def build(**names):
            text = DECK.format(
                pdb=pdb, ff=EXAMPLES / "formaldehyde.xml", tip=EXAMPLES / "tip3p.xml",
                nsteps=1, ensemble="nve", fmt="pdb", name="out",
                md=f"snapshot = {snapshot}")
            for key, value in names.items():
                text = re.sub(rf"(?m)^{key} = .*$", f"{key} = {value}", text)
            deck = Path("clash.inp")
            deck.write_text(text)
            return QMMM_MD(oqp_cfg=str(deck))

        for key in ("energy_file", "trajectory_file", "log_file"):
            with self.subTest(output=key):
                with self.assertRaisesRegex(ValueError, "would overwrite the input"):
                    build(**{key: snapshot})
        # np.savez would append .npz: the name actually written is what counts
        with self.assertRaisesRegex(ValueError, "would overwrite the input"):
            build(energy_file="src.snapshot.00000002")
        with self.assertRaisesRegex(ValueError, "pdb_file"):
            build(trajectory_file=pdb)
        with self.assertRaisesRegex(ValueError, "outputs must be distinct"):
            build(trajectory_file="same.out", log_file="same.out")
        # every force-field file is protected, not only the last one named,
        # and so is the deck itself
        import shutil
        first, last = Path("first.xml"), Path("last.xml")
        shutil.copy(EXAMPLES / "formaldehyde.xml", first)
        shutil.copy(EXAMPLES / "tip3p.xml", last)
        for victim in (first, last, Path("clash.inp")):
            # (energy_file is written as <name>.npz, so it cannot land on these)
            for key in ("restart_file", "trajectory_file", "log_file"):
                with self.subTest(victim=victim.name, output=key):
                    kept = None if victim.name == "clash.inp" else victim.read_bytes()
                    with self.assertRaisesRegex(ValueError, "would overwrite the input"):
                        build(forcefield_files=f"{first} {last}", **{key: victim})
                    if kept is not None:
                        self.assertEqual(victim.read_bytes(), kept)
        # the velocity file and the QM geometry override are inputs too
        natom = sum(line.startswith(("ATOM", "HETATM"))
                    for line in pdb.read_text().splitlines())
        velocity = Path("start.vel")
        np.savetxt(velocity, np.zeros((natom, 3)))
        geometry = Path("qm.xyz")
        geometry.write_text("4\n\nC 0 0 0\nO 0 0 1.2\nH 0.9 0 -0.5\nH -0.9 0 -0.5\n")

        def build_reading(victim, key):
            text = DECK.format(
                pdb=pdb, ff=EXAMPLES / "formaldehyde.xml", tip=EXAMPLES / "tip3p.xml",
                nsteps=1, ensemble="nve", fmt="pdb", name="out",
                md=f"velocity = {velocity}\ncommon_control_keys = velocity")
            text = text.replace("cutoff = NoCutoff",
                                f"cutoff = NoCutoff\nqm_atoms_xyz = {geometry}")
            text = re.sub(rf"(?m)^{key} = .*$", f"{key} = {victim}", text)
            Path("reads.inp").write_text(text)
            return QMMM_MD(oqp_cfg="reads.inp")

        for victim, label in ((velocity, "velocity"), (geometry, "qm_atoms_xyz")):
            for key in ("restart_file", "trajectory_file", "log_file"):
                with self.subTest(victim=victim.name, output=key):
                    kept = victim.read_bytes()
                    with self.assertRaisesRegex(
                            ValueError, f"would overwrite the input .*{label}"):
                        build_reading(victim, key)
                    self.assertEqual(victim.read_bytes(), kept)
        build_reading("elsewhere.out", "log_file")           # control: deck is valid
        self.assertEqual(snapshot.read_bytes(), before)
        self.assertEqual(pdb.read_bytes(), pdb_before)
        build().run()                                    # distinct names still run
        self.assertEqual(snapshot.read_bytes(), before)

    def test_restart_keeps_the_energy_table_when_the_text_log_is_gone(self):
        zero = "velocity = zero\ncommon_control_keys = velocity"
        first = _run("nolog", 2)[1]
        os.remove("nolog.dat")
        data = _run("nolog", 4, md=zero + "\nrestart = true")[1]
        self.assertEqual(list(data["step"]), [0, 1, 2, 3, 4])
        np.testing.assert_allclose(data["E_tot"][:3], first["E_tot"], rtol=0, atol=1e-9)
        with np.load("nolog.npz") as saved:
            self.assertEqual(list(saved["step"]), [0, 1, 2, 3, 4])
        with open("nolog.dat") as stream:                 # a new log, with its header
            self.assertTrue(stream.readline().startswith("#"))

    def test_sectioned_md_controls_are_honoured_without_a_marker(self):
        """[md] nstep / velocity / ensemble in a plain .inp deck used to be
        ignored: the driver took [qmmm] n_steps (1000) and Maxwell velocities."""
        from oqp.library.qmmm_md import QMMM_MD
        deck = Path("plain.inp")
        deck.write_text(DECK.format(
            pdb=EXAMPLES / "formaldehyde_water.pdb", ff=EXAMPLES / "formaldehyde.xml",
            tip=EXAMPLES / "tip3p.xml", nsteps=1000, ensemble="nve", fmt="pdb",
            name="plain", md="nstep = 3\ndt = 0.25\nvelocity = zero\nensemble = nvt"))
        driver = QMMM_MD(oqp_cfg=str(deck))
        import openmm.unit as u
        self.assertEqual(driver.n_steps, 3)
        self.assertAlmostEqual(driver.timestep.value_in_unit(u.femtoseconds), 0.25)
        self.assertEqual(driver.velocity_source, "zero")
        self.assertEqual(driver.ensemble, "nvt")
        # [qmmm] spellings still work when [md] does not set the control
        deck.write_text(DECK.format(
            pdb=EXAMPLES / "formaldehyde_water.pdb", ff=EXAMPLES / "formaldehyde.xml",
            tip=EXAMPLES / "tip3p.xml", nsteps=7, ensemble="nve", fmt="pdb",
            name="plain", md="snapshot_interval = 0"))
        self.assertEqual(QMMM_MD(oqp_cfg=str(deck)).n_steps, 7)

    def test_a_restarted_thermostat_does_not_replay_its_noise(self):
        nvt = "seed = 7\ncommon_control_keys = seed"
        first = _run("noise", 2, md=nvt, ensemble="nvt")[0]
        self.assertEqual(first.integrator_seed, first.random_seed)
        again = _run("noise", 4, md=nvt + "\nrestart = true", ensemble="nvt")[0]
        self.assertEqual(again.random_seed, first.random_seed)       # same user seed
        self.assertNotEqual(again.integrator_seed, first.integrator_seed)
        from oqp.utils.md_controls import continuation_seed
        self.assertEqual(again.integrator_seed, continuation_seed(again.random_seed, 2))
        self.assertEqual(again.simulation_md.integrator.getRandomNumberSeed(),
                         again.integrator_seed)
        del first, again

    def test_dcd_trajectories_are_appended(self):
        zero = "velocity = zero\ncommon_control_keys = velocity"
        _run("d", 2, fmt="dcd")
        _run("d", 4, md=zero + "\nrestart = true", fmt="dcd")
        with open("d.dcd", "rb") as stream:
            frames = int(np.frombuffer(stream.read(12)[8:12], dtype="<i4")[0])
        self.assertEqual(frames, 5)

    def test_numbered_snapshots_start_a_new_run(self):
        """snapshot_interval writes phase-space points; snapshot= starts a new
        trajectory (step 0) from one, with its positions and velocities."""
        import openmm.unit as u
        nvt = "seed = 7\ncommon_control_keys = seed\nsnapshot_interval = 2"
        _run("eq", 4, md=nvt, ensemble="nvt")
        self.assertTrue(os.path.isfile("eq.snapshot.00000002.npz"))
        self.assertTrue(os.path.isfile("eq.snapshot.00000004.npz"))
        self.assertFalse(os.path.isfile("eq.snapshot.00000000.npz"))
        snap = md_snapshot.read_snapshot("eq.snapshot.00000004.npz")
        self.assertEqual(snap["step"], 4)
        self.assertGreater(np.abs(snap["velocities_nm_ps"]).max(), 0.0)
        # the stored on-step velocities differ from the integrator's half-step ones
        self.assertFalse(np.allclose(snap["velocities_nm_ps"],
                                     snap["integrator_velocities_nm_ps"]))

        from oqp.library.qmmm_md import QMMM_MD
        deck = Path("prod.inp")
        deck.write_text(DECK.format(
            pdb=EXAMPLES / "formaldehyde_water.pdb", ff=EXAMPLES / "formaldehyde.xml",
            tip=EXAMPLES / "tip3p.xml", nsteps=1, ensemble="nve", fmt="pdb",
            name="prod", md="snapshot = eq.snapshot.00000004.npz"))
        d = QMMM_MD(oqp_cfg=str(deck))
        d.setup()
        state = d.simulation_md.context.getState(getPositions=True, getVelocities=True)
        self.assertEqual(d.simulation_md.currentStep, 0)          # a NEW run
        np.testing.assert_allclose(
            np.asarray(state.getPositions(asNumpy=True).value_in_unit(u.nanometer)),
            snap["positions_nm"], atol=1e-6)
        got = np.asarray(state.getVelocities(asNumpy=True).value_in_unit(
            u.nanometer / u.picosecond))
        # the integrator receives them half a kick EARLIER than stored: its
        # velocities live half a step behind the positions (QM atoms carry no
        # constraint, so the relation is exact for them)
        dt_ps = 0.0005
        expected = (snap["velocities_nm_ps"][:4]
                    - 0.5 * dt_ps * d._force_kj_nm[:4] / d._sys0_masses[:4, None])
        np.testing.assert_allclose(got[:4], expected, atol=1e-5)
        self.assertFalse(np.allclose(got[:4], snap["velocities_nm_ps"][:4], atol=1e-5))
        d._log_handle.close()

        with self.assertRaisesRegex(ValueError, "cannot be combined with snapshot"):
            deck.write_text(deck.read_text().replace(
                "snapshot = eq", "velocity = zero\ncommon_control_keys = velocity\nsnapshot = eq"))
            QMMM_MD(oqp_cfg=str(deck))

    def test_snapshot_carries_the_cell(self):
        """An equilibrated cell (classical NPT) must reach the QM/MM contexts."""
        from oqp.library.qmmm_md import QMMM_MD
        import openmm.app as app
        import openmm.unit as u
        pdb = app.PDBFile(str(EXAMPLES / "water_dimer.pdb"))
        xyz = np.array(pdb.positions.value_in_unit(u.nanometer))
        masses = [a.element.mass.value_in_unit(u.dalton) for a in pdb.topology.atoms()]
        md_snapshot.write_snapshot("classical.npz", positions_nm=xyz + 0.01,
                                   velocities_nm_ps=np.zeros_like(xyz),
                                   masses_dalton=masses, box_nm=[2.9, 3.05, 3.1])
        deck = Path("w.inp")
        deck.write_text(DECK.format(
            pdb=EXAMPLES / "water_dimer.pdb", ff="", tip=EXAMPLES / "tip3p.xml",
            nsteps=1, ensemble="nve", fmt="pdb", name="w",
            md="snapshot = classical.npz").replace("1 2 3 4", "1 2 3")
            .replace("qm_atoms = 0-3", "qm_atoms = 0-2")
            .replace("cutoff = NoCutoff", "cutoff = PME"))
        d = QMMM_MD(oqp_cfg=str(deck))
        d.setup()
        np.testing.assert_allclose(d._box_lengths_nm(), [2.9, 3.05, 3.1])
        for key in ("sim0", "simew"):
            box = d.mm_systems[key].context.getState().getPeriodicBoxVectors()
            np.testing.assert_allclose(
                [box[i][i].value_in_unit(u.nanometer) for i in range(3)],
                [2.9, 3.05, 3.1])
        d._log_handle.close()
        # a snapshot of another system is refused
        md_snapshot.write_snapshot("wrong.npz", positions_nm=xyz[:5],
                                   velocities_nm_ps=np.zeros((5, 3)),
                                   masses_dalton=masses[:5])
        deck.write_text(deck.read_text().replace("classical.npz", "wrong.npz"))
        with self.assertRaisesRegex(ValueError, "has 5 atoms, the PDB has 6"):
            QMMM_MD(oqp_cfg=str(deck))


@unittest.skipUnless(_runtime_available(), "OpenQP runtime with OpenMM required")
class TestSurfaceHoppingStartsFromASnapshot(unittest.TestCase):
    """QM/MM NAMD used to draw fresh Maxwell-Boltzmann velocities for every
    atom and rescale them to init_temp, whatever preceded it."""

    def _namd(self, tmp, name, md_extra):
        deck = Path(tmp) / f"{name}.oqp"
        deck.write_text(
            'mrsf(nstate=2)/bhhlyp/6-31g* namd(S1,scheme=Overlap)\n'
            f'md(nstep=1,dt=0.5,trajectory_interval=1{md_extra})\n'
            f'qmmm(pdb_file="{EXAMPLES / "formaldehyde_water.pdb"}",'
            f'forcefield_files="{EXAMPLES / "formaldehyde.xml"} {EXAMPLES / "tip3p.xml"}",'
            'qm_atoms="0-3",cutoff=NoCutoff)\n'
            f'geom="{ROOT / "examples" / "geometries" / "CH2O-2bc62dda4b8a.xyz"}"\n')
        subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(deck)], cwd=tmp,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        from oqp.library.namd import read_namd_trajectory
        _, records = read_namd_trajectory(str(Path(tmp) / f"{name}.namd.trj"))
        return (np.asarray(records["coordinates_bohr"][0]).reshape(-1, 3),
                np.asarray(records["velocities_au"][0]).reshape(-1, 3),
                float(records["e_kin_hartree"][0]),
                (Path(tmp) / f"{name}.log").read_text())

    def test_frame_zero_is_the_snapshot(self):
        snap = md_snapshot.read_snapshot(EXAMPLES / "formaldehyde_water.snapshot.npz")
        v_snap = snap["velocities_nm_ps"] / md_snapshot.AU_VELOCITY_TO_NM_PER_PS
        with tempfile.TemporaryDirectory() as tmp:
            r, v, ke, log = self._namd(
                tmp, "snap", f',snapshot="{EXAMPLES / "formaldehyde_water.snapshot.npz"}"')
            _, v_ctl, ke_ctl, log_ctl = self._namd(tmp, "ctl", ",seed=3")
        np.testing.assert_allclose(r, snap["positions_nm"] * md_snapshot.NM_TO_BOHR,
                                   atol=1e-7)
        # QM atoms carry no constraint: their velocities arrive unchanged
        np.testing.assert_allclose(v[:4], v_snap[:4], rtol=0, atol=1e-12)
        # rigid MM water: the snapshot is written constraint-consistent, so
        # the RATTLE projection at the start of the run has nothing to remove
        np.testing.assert_allclose(v[4:], v_snap[4:], rtol=0, atol=1e-7)
        self.assertIn("QM/MM NAMD initial conditions", log)
        self.assertIn("no velocities drawn or rescaled", log)
        # control: without a snapshot the run draws 300 K velocities, which is
        # neither the snapshot's kinetic energy nor its velocities
        self.assertNotIn("QM/MM NAMD initial conditions", log_ctl)
        self.assertFalse(np.allclose(v_ctl[:4], v_snap[:4], atol=1e-6))
        ndof = 3 * 19 - 15 - 3
        self.assertAlmostEqual(ke_ctl, 0.5 * ndof * 3.166811563e-6 * 300.0, places=6)
        self.assertNotAlmostEqual(ke, ke_ctl, places=4)

    def test_namd_does_not_checkpoint_over_its_input_snapshot(self):
        """restart_file=<the snapshot> would replace the starting point with a
        surface-hopping checkpoint at step 0."""
        import shutil
        with tempfile.TemporaryDirectory() as tmp:
            snapshot = Path(tmp) / "state.npz"
            shutil.copy(EXAMPLES / "formaldehyde_water.snapshot.npz", snapshot)
            before = snapshot.read_bytes()
            deck = Path(tmp) / "alias.oqp"
            deck.write_text(
                'mrsf(nstate=2)/bhhlyp/6-31g* namd(S1,scheme=Overlap)\n'
                f'md(nstep=1,dt=0.5,snapshot="{snapshot}",restart_file="{snapshot}")\n'
                f'qmmm(pdb_file="{EXAMPLES / "formaldehyde_water.pdb"}",'
                f'forcefield_files="{EXAMPLES / "formaldehyde.xml"} {EXAMPLES / "tip3p.xml"}",'
                'qm_atoms="0-3",cutoff=NoCutoff)\n'
                f'geom="{ROOT / "examples" / "geometries" / "CH2O-2bc62dda4b8a.xyz"}"\n')
            result = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(deck)],
                                    cwd=tmp, capture_output=True, text=True)
            after = snapshot.read_bytes()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("restart_file aliases simulation input snapshot",
                      result.stdout + result.stderr)
        self.assertEqual(after, before)

    def test_gas_phase_namd_rejects_a_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            deck = Path(tmp) / "gas.oqp"
            deck.write_text(
                'mrsf(nstate=2)/bhhlyp/6-31g* namd(S1,scheme=Overlap)\n'
                f'md(nstep=1,snapshot="{EXAMPLES / "formaldehyde_water.snapshot.npz"}")\n'
                f'geom="{ROOT / "examples" / "geometries" / "CH2O-2bc62dda4b8a.xyz"}"\n')
            result = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(deck)],
                                    cwd=tmp, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("QM/MM phase-space point", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
