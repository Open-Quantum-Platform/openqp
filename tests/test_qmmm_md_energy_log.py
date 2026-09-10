"""The ground-state QM/MM MD driver (runtype=md) must log the Hamiltonian.

The driver drives OpenMM with a CustomExternalForce whose energy expression is
LINEAR in the coordinates (``-F.r + E - ecorr``).  Read after the Verlet step,
OpenMM's potential is therefore only the first-order extrapolation
``E(r0) - F(r0).(r - r0)`` of the QM/MM energy, with a second-order
remainder that grows with the number of atoms.  The driver now samples E_pot
(the energy returned by the force backend) and the time-centred kinetic
energy at the positions the force was computed for, before the step.

The driver also has to honour ``[qmmm] rigidwater``: integrating the TIP3P
O-H stretch explicitly at 0.5 fs makes the total energy of an 18 500-atom
box fluctuate by ~90 kJ/mol (pure MM gives the same), against 0.8 kJ/mol with
the water constraints on.
"""
import io
import os
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "pyoqp" / "oqp" / "library" / "qmmm_md.py"
EXAMPLES = ROOT / "examples" / "QMMM"

try:
    import numpy as np
    import openmm  # noqa: F401
    _HAVE = True
except Exception:  # pragma: no cover
    _HAVE = False


class TestEnergyIsSampledBeforeTheStep(unittest.TestCase):
    """Order of operations pinned in the source: the energy row of a step is
    taken between the force update and the Verlet step, from the backend's
    energy, not from OpenMM's linearised potential."""

    def setUp(self):
        self.src = SRC.read_text()
        m = re.search(r"\n    def step\(self\):.*?\n    def _report_energies", self.src, re.S)
        self.assertIsNotNone(m)
        self.step_src = m.group(0)

    def test_sampling_precedes_the_verlet_step(self):
        i_update = self.step_src.index("self._update_qmmm_force(pos_pre)")
        i_energy = self.step_src.index("E_pot = self._qmmm_energy_kJ")
        i_kin = self.step_src.index("getKineticEnergy()")
        i_report = self.step_src.index("self._report_energies(")
        i_step = self.step_src.index("self.simulation_md.step(1)")
        self.assertLess(i_update, i_energy)
        self.assertLess(i_energy, i_step)
        self.assertLess(i_kin, i_step)
        self.assertLess(i_report, i_step)

    def test_no_openmm_potential_read_in_step(self):
        self.assertNotIn("getPotentialEnergy", self.step_src)
        self.assertNotIn("app.StateDataReporter(", self.src)

    def test_backend_energy_is_recorded_by_the_force_update(self):
        m = re.search(r"\n    def _update_qmmm_force\(self, positions\):.*?\n    def ", self.src, re.S)
        self.assertIn("self._qmmm_energy_kJ = _to_kJmol(qmmm_energy)", m.group(0))


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestRigidWaterIsWired(unittest.TestCase):
    def test_source(self):
        src = SRC.read_text()
        self.assertIn('qmmm_cfg.get("rigidwater", True)', src)
        self.assertIn("rigidWater=True", src)
        self.assertIn("self.system_md.addConstraint(p1, p2, dist)", src)
        self.assertIn("3 * self.system_md.getNumParticles() - self.n_constraints", src)


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestReportRow(unittest.TestCase):
    """The driver's own CSV writer: header, step-0 row, cadence, columns."""

    def _driver(self, ensemble="nve", interval=1):
        from oqp.library.qmmm_md import QMMM_MD
        d = QMMM_MD.__new__(QMMM_MD)
        d.ensemble = ensemble
        d.report_interval = interval
        d._log_handle = io.StringIO()
        d._wall_t0 = None
        return d

    def test_rows_and_cadence(self):
        from contextlib import redirect_stdout
        d = self._driver(interval=2)
        out = io.StringIO()
        with redirect_stdout(out):
            for s in range(5):
                d._report_energies(s, s * 0.0005, -100.0 + s, 10.0, -90.0 + s, 300.0, float("nan"))
        rows = d._log_handle.getvalue().splitlines()
        self.assertEqual(len(rows), 3)                       # steps 0, 2, 4
        cols = rows[0].split(",")
        self.assertEqual(len(cols), 6)                       # t, Epot, Ekin, Etot, T, speed
        self.assertAlmostEqual(float(cols[1]), -100.0)
        self.assertAlmostEqual(float(rows[2].split(",")[3]), -86.0)
        # stdout mirrors the file with a leading step index
        srows = out.getvalue().splitlines()
        self.assertEqual([r.split(",")[0] for r in srows], ["0", "2", "4"])

    def test_npt_adds_volume_column(self):
        d = self._driver(ensemble="npt")
        d._report_energies(0, 0.0, -1.0, 2.0, 1.0, 300.0, 27.5)
        cols = d._log_handle.getvalue().strip().split(",")
        self.assertEqual(len(cols), 7)
        self.assertAlmostEqual(float(cols[5]), 27.5)


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.library.qmmm_md import QMMM_MD  # noqa: F401
        import openmm  # noqa: F401
        return (EXAMPLES / "formaldehyde_water.pdb").exists()
    except Exception:
        return False


DECK = """[input]
system=
   6   0.000000   0.000000   0.000000
   8   0.000000   0.000000   1.203000
   1   0.000000   0.943000  -0.589000
   1   0.000000  -0.943000  -0.589000
charge=0
runtype=md
basis=sto-3g
method=hf
qmmm_flag=True
[guess]
type=huckel
[scf]
type=rhf
multiplicity=1
conv=1e-10
[qmmm]
pdb_file={pdb}
forcefield_files={ff} {tip}
qm_atoms=0-3
cutoff=NoCutoff
embedding=electrostatic
n_steps=6
timestep=0.5
ensemble=nve
report_interval=1
temperature=300.0
energy_file=e.npz
log_file=e.dat
trajectory_file=t.pdb
trajectory_format=pdb
"""


@unittest.skipUnless(_HAVE and _runtime_available(), "OpenMM or compiled OpenQP runtime unavailable")
class TestLoggedEnergyIsTheBackendEnergy(unittest.TestCase):
    """Six steps of formaldehyde in five waters: every logged E_pot equals the
    energy the force backend returns at the trajectory frame of that step,
    and the logged total is what step() returned."""

    def test_rows_match_backend_energy(self):
        import tempfile
        import openmm.unit as u
        from oqp.library.qmmm_md import QMMM_MD
        with tempfile.TemporaryDirectory() as tmp:
            deck = Path(tmp) / "md.inp"
            deck.write_text(DECK.format(pdb=EXAMPLES / "formaldehyde_water.pdb",
                                        ff=EXAMPLES / "formaldehyde.xml",
                                        tip=EXAMPLES / "tip3p.xml"))
            cwd = os.getcwd(); os.chdir(tmp)
            try:
                d = QMMM_MD(oqp_cfg=str(deck))
                d.setup()
                # independent reference at the start geometry, before any step
                e0, _ = d.oqp_driver.compute_force(d.pdb.positions, d.pdb.topology,
                                                   d.mm_systems, d.qm_atoms)
                e0 = float(e0.value_in_unit(u.kilojoule_per_mole)) if u.is_quantity(e0) else float(e0)
                # positions each step's force (and energy row) was computed for
                seen = []
                _upd = d._update_qmmm_force
                def upd(positions):
                    seen.append(positions); return _upd(positions)
                d._update_qmmm_force = upd
                returned = [d.step() for _ in range(6)]
                d._save_traj_data()
                z = np.load("e.npz")
                rows = np.loadtxt("e.dat", delimiter=",", skiprows=1)
                # row 3 re-evaluated independently at the positions of step 3
                # (the MM contexts read their positions from the caller)
                for sim in d.mm_systems.values():
                    if hasattr(sim, "context"):
                        sim.context.setPositions(seen[3])
                e3, _ = d.oqp_driver.compute_force(seen[3], d.pdb.topology, d.mm_systems, d.qm_atoms)
                e3 = float(e3.value_in_unit(u.kilojoule_per_mole)) if u.is_quantity(e3) else float(e3)
                # the OLD quantity (OpenMM's linear potential after the step) is
                # not what is logged: it differs from the exact energy at the
                # same positions by the second-order remainder
                for sim in d.mm_systems.values():
                    if hasattr(sim, "context"):
                        sim.context.setPositions(seen[-1])
                d.step()
                lin = d.simulation_md.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                d._log_handle.close()
            finally:
                os.chdir(cwd)
        self.assertEqual(d.system_md.getNumConstraints(), 15)      # 5 rigid TIP3P waters, QM H2CO free
        self.assertEqual(d.n_constraints, 15)
        self.assertEqual(list(z["step"]), [0, 1, 2, 3, 4, 5])
        self.assertAlmostEqual(z["E_pot"][0], e0, places=6)
        self.assertAlmostEqual(z["E_pot"][3], e3, places=5)
        self.assertEqual(len(seen), 7)
        self.assertNotAlmostEqual(lin, d._qmmm_energy_kJ, places=3)
        np.testing.assert_allclose(z["E_tot"], returned, rtol=0, atol=1e-9)
        np.testing.assert_allclose(rows[:, 3], z["E_tot"], rtol=0, atol=1e-7)
        np.testing.assert_allclose(z["E_pot"] + z["E_kin"], z["E_tot"], rtol=0, atol=1e-9)

    def test_rigidwater_false_leaves_water_flexible(self):
        import tempfile
        from oqp.library.qmmm_md import QMMM_MD
        with tempfile.TemporaryDirectory() as tmp:
            deck = Path(tmp) / "md.inp"
            deck.write_text(DECK.format(pdb=EXAMPLES / "formaldehyde_water.pdb",
                                        ff=EXAMPLES / "formaldehyde.xml",
                                        tip=EXAMPLES / "tip3p.xml").replace("ensemble=nve", "ensemble=nve\nrigidwater=false"))
            cwd = os.getcwd(); os.chdir(tmp)
            try:
                d = QMMM_MD(oqp_cfg=str(deck)); d.setup()
                self.assertEqual(d.system_md.getNumConstraints(), 0)
                self.assertEqual(d.n_constraints, 0)
                d._log_handle.close()
            finally:
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
