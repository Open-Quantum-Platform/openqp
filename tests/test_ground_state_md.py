import csv
from types import SimpleNamespace

import numpy as np
import pytest

try:
    from oqp.library import ground_state_md
except (ImportError, RuntimeError) as error:
    pytest.skip(f"installed OpenQP runtime required: {error}", allow_module_level=True)


class _Data(dict):
    pass


class _Mol:
    def __init__(self, tmp_path, velocity_file):
        self.log = str(tmp_path / "bomd.log")
        self.data = _Data(natom=2)
        self.config = {
            "md": {
                "nstep": 2,
                "dt": 1.0,
                "init_temp": 300.0,
                "velocity": str(velocity_file),
                "seed": 17,
                "rng_stream": 3,
                "thermostat": "off",
                "ensemble": "nve",
                "trajectory_file": str(tmp_path / "trajectory.xyz"),
                "energy_file": str(tmp_path / "energy.csv"),
                "mo_reuse": True,
            },
            "guess": {"type": "huckel"},
            "properties": {"grad": []},
        }
        self._coordinates = np.array([[-0.7, 0.0, 0.0], [0.7, 0.0, 0.0]])
        self.energies = [0.0]
        self.grads = [np.zeros((2, 3))]

    def get_mass(self):
        return np.array([1.0, 1.0])

    def get_atoms(self):
        return np.array([1, 1])

    def get_system(self):
        return self._coordinates.reshape(-1)

    def update_system(self, coordinates):
        self._coordinates = np.asarray(coordinates, dtype=float).reshape((2, 3))


class _ConstantForceMD(ground_state_md.GroundStateMD):
    def _electronic(self, *, continuation):
        return -1.0, np.zeros((self.natom, 3))


def test_ground_state_md_propagates_and_writes_paired_outputs(tmp_path, monkeypatch):
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0.0001 0 0\n-0.0001 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    monkeypatch.setattr(ground_state_md, "dump_log", lambda *args, **kwargs: None)

    driver = _ConstantForceMD(mol)
    initial = mol._coordinates.copy()
    driver.run()

    expected = initial + 2 * driver.dt * np.array(
        [[0.0001, 0.0, 0.0], [-0.0001, 0.0, 0.0]]
    )
    np.testing.assert_allclose(mol._coordinates, expected, atol=1.0e-14)
    np.testing.assert_allclose(mol.md_velocity, driver.velocity)

    xyz_lines = (tmp_path / "trajectory.xyz").read_text(encoding="utf-8").splitlines()
    assert len(xyz_lines) == 3 * (2 + 2)
    with (tmp_path / "energy.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert [int(row["step"]) for row in rows] == [0, 1, 2]
    assert all(float(row["potential_hartree"]) == -1.0 for row in rows)


def test_langevin_requires_positive_friction(tmp_path):
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.config["md"].update({"thermostat": "langevin", "friction": 0.0})
    mol.config["md"]["thermostat_friction"] = 0.0
    try:
        ground_state_md.GroundStateMD(mol)
    except ValueError as error:
        assert "friction must be positive" in str(error)
    else:
        raise AssertionError("zero-friction Langevin thermostat was accepted")


def test_ground_state_md_rejects_multiple_mpi_ranks(tmp_path):
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.mpi_manager = SimpleNamespace(use_mpi=1, size=2, rank=0)
    with pytest.raises(ValueError, match="requires a single MPI rank"):
        ground_state_md.GroundStateMD(mol)


def test_electronic_overrides_are_temporary_across_repeated_calls(
        tmp_path, monkeypatch):
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.config["properties"]["grad"] = [7]
    mol.config["input"] = {"basis": "6-31g"}
    mol.config["scf"] = {"init_scf": "rhf", "init_basis": "sto-3g"}
    seen = []

    class FakeSinglePoint:
        def __init__(self, current_mol):
            self.mol = current_mol

        def energy(self, do_init_scf=True):
            seen.append((
                self.mol.config["guess"]["type"],
                list(self.mol.config["properties"]["grad"]),
                do_init_scf,
            ))
            self.mol.energies = [-2.0]

    class FakeGradient:
        def __init__(self, current_mol):
            self.mol = current_mol

        def gradient(self):
            self.mol.grads = [np.zeros((2, 3))]

    class FakeLastStep:
        def __init__(self, current_mol):
            self.mol = current_mol

        def compute(self, current_mol, grad_list):
            assert current_mol is self.mol
            assert grad_list == [0]

    monkeypatch.setattr(ground_state_md, "SinglePoint", FakeSinglePoint)
    monkeypatch.setattr(ground_state_md, "Gradient", FakeGradient)
    monkeypatch.setattr(ground_state_md, "LastStep", FakeLastStep)
    driver = ground_state_md.GroundStateMD(mol)

    driver._electronic(continuation=True)
    assert mol.config["guess"]["type"] == "huckel"
    assert mol.config["properties"]["grad"] == [7]

    mol.config["md"]["mo_reuse"] = False
    driver._electronic(continuation=True)
    assert seen == [("previous", [0], False), ("huckel", [0], True)]
    assert mol.config["guess"]["type"] == "huckel"
    assert mol.config["properties"]["grad"] == [7]


def test_electronic_overrides_are_restored_after_failure(tmp_path, monkeypatch):
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.config["properties"]["grad"] = [9]

    class FailingSinglePoint:
        def __init__(self, current_mol):
            self.mol = current_mol

        def energy(self, do_init_scf=True):
            assert self.mol.config["guess"]["type"] == "previous"
            assert self.mol.config["properties"]["grad"] == [0]
            assert do_init_scf is False
            raise RuntimeError("electronic failure")

    monkeypatch.setattr(ground_state_md, "SinglePoint", FailingSinglePoint)
    driver = ground_state_md.GroundStateMD(mol)
    with pytest.raises(RuntimeError, match="electronic failure"):
        driver._electronic(continuation=True)
    assert mol.config["guess"]["type"] == "huckel"
    assert mol.config["properties"]["grad"] == [9]


@pytest.mark.parametrize(
    "override, message",
    (
        ({"restart": True}, "restart is not available"),
        ({"restart_file": "state.npz"}, "restart_file is not available"),
    ),
)
def test_ground_state_md_rejects_restart_requests(tmp_path, override, message):
    """restart/restart_file are public md(...) options this driver cannot honour.

    It writes no checkpoint and reads none, so an accepted-but-ignored request
    makes fresh velocities, starts at step zero and then opens the trajectory
    and energy files of the run it was asked to continue in write mode --
    losing exactly the outputs the user wanted extended.  Reject before the
    first electronic evaluation and before any file is opened.
    """
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.config["md"].update(override)

    # outputs of the "previous" run, which must survive the rejection
    trajectory = tmp_path / "trajectory.xyz"
    energy = tmp_path / "energy.csv"
    trajectory.write_text("previous\n", encoding="utf-8")
    energy.write_text("previous\n", encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        ground_state_md.GroundStateMD(mol)

    assert trajectory.read_text(encoding="utf-8") == "previous\n"
    assert energy.read_text(encoding="utf-8") == "previous\n"


def test_ground_state_md_accepts_the_materialised_restart_defaults(tmp_path):
    """The schema always supplies restart=False and restart_file="", so the
    rejection must not fire on an ordinary deck."""
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.config["md"].update({"restart": False, "restart_file": "",
                             "restart_interval": 10})
    driver = ground_state_md.GroundStateMD(mol)
    assert driver.nstep == 2


def _frames_and_rows(tmp_path):
    lines = (tmp_path / "trajectory.xyz").read_text(encoding="utf-8").splitlines()
    frames = len(lines) // 4                 # 2 atoms + count line + comment line
    with open(tmp_path / "energy.csv", newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))[1:]
    return frames, [int(row[0]) for row in rows]


def test_trajectory_interval_sets_the_xyz_cadence(tmp_path, monkeypatch):
    """md(trajectory_interval=N) was accepted for gas-phase MD and ignored:
    every step was written.  The xyz trajectory follows it now; the energy
    table keeps every step."""
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    monkeypatch.setattr(ground_state_md, "dump_log", lambda *args, **kwargs: None)
    mol = _Mol(tmp_path, velocity_file)
    mol.config["md"].update({"nstep": 5, "trajectory_interval": 2})
    _ConstantForceMD(mol).run()
    frames, steps = _frames_and_rows(tmp_path)
    assert frames == 3                       # steps 0, 2, 4
    assert steps == [0, 1, 2, 3, 4, 5]

    # the default writes every step, as before
    mol = _Mol(tmp_path, velocity_file)
    mol.config["md"]["nstep"] = 3
    _ConstantForceMD(mol).run()
    assert _frames_and_rows(tmp_path) == (4, [0, 1, 2, 3])


def test_trajectory_interval_zero_is_about_ten_femtoseconds(tmp_path):
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.config["md"].update({"dt": 0.5, "trajectory_interval": 0})
    assert ground_state_md.GroundStateMD(mol).trajectory_interval == 20
    mol.config["md"]["trajectory_interval"] = -1
    with pytest.raises(ValueError, match="trajectory_interval must be >= 0"):
        ground_state_md.GroundStateMD(mol)


@pytest.mark.parametrize("output", ("trajectory_file", "energy_file"))
@pytest.mark.parametrize("victim", ("velocity", "geometry", "deck", "log"))
def test_outputs_may_not_overwrite_an_input(tmp_path, output, victim):
    """run() opens both outputs with "w"; pointing one at the geometry, the
    velocity file, the input deck or the log would destroy that file."""
    velocity_file = tmp_path / "velocity.txt"
    velocity_file.write_text("0 0 0\n0 0 0\n", encoding="utf-8")
    geometry = tmp_path / "start.xyz"
    geometry.write_text("2\n\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    deck = tmp_path / "job.inp"
    deck.write_text("[input]\nruntype = md\n", encoding="utf-8")
    mol = _Mol(tmp_path, velocity_file)
    mol.input_file = str(deck)
    mol.config["input"] = {"system": str(geometry)}
    target = {"velocity": velocity_file, "geometry": geometry, "deck": deck,
              "log": tmp_path / "bomd.log"}[victim]
    before = target.read_text(encoding="utf-8") if target.exists() else None
    mol.config["md"][output] = str(target)
    with pytest.raises(ValueError, match="would overwrite"):
        ground_state_md.GroundStateMD(mol)
    if before is not None:
        assert target.read_text(encoding="utf-8") == before
