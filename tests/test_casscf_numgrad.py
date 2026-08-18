"""Source-only checks for CASSCF central-difference nuclear gradients."""

from contextlib import contextmanager
import importlib.util
from pathlib import Path
import sys
import types

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _load_input_checker(monkeypatch):
    oqp = sys.modules.get("oqp", types.ModuleType("oqp"))
    utils = sys.modules.get("oqp.utils", types.ModuleType("oqp.utils"))
    monkeypatch.setitem(sys.modules, "oqp", oqp)
    monkeypatch.setitem(sys.modules, "oqp.utils", utils)
    monkeypatch.setattr(oqp, "utils", utils, raising=False)
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
    mpi_utils.MPIManager = type("MPIManager", (), {"use_mpi": False, "size": 1})
    monkeypatch.setitem(sys.modules, "oqp.utils.mpi_utils", mpi_utils)
    name = "input_checker_casscf_numgrad_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/utils/input_checker.py")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


def _casscf_config(method="casscf", runtype="grad"):
    config = {
        "input": {"method": method, "runtype": runtype, "functional": ""},
        "scf": {"type": "rhf", "multiplicity": 1},
        "cas": {"active_electrons": 2, "active_orbitals": 2},
        "ci": {"nroot": 1},
        "casscf": {"root": 0},
        "properties": {"grad": [0]},
        "optimize": {"istate": 0},
    }
    if method in {"sa-casscf", "sacasscf"}:
        config["ci"]["nroot"] = 2
        config["state_average"] = {
            "enabled": True, "nstate": 2, "target_roots": [0, 1],
        }
    return config


def _check(config, monkeypatch):
    checker = _load_input_checker(monkeypatch)
    report = checker.CheckReport()
    checker._check_casci(config, report)
    return report


def test_checker_accepts_casscf_and_sa_casscf_gradient_driven_runtypes(
        monkeypatch):
    for method in ("casscf", "sa-casscf"):
        for runtype in ("grad", "optimize", "ts", "mep", "irc"):
            config = _casscf_config(method, runtype)
            if method == "sa-casscf" and runtype != "grad":
                config["optimize"]["istate"] = 1
            report = _check(config, monkeypatch)
            assert report.ok, report.to_text()


def test_checker_rejects_physical_root_as_state_specific_gradient_slot(
        monkeypatch):
    config = _casscf_config()
    config["ci"]["nroot"] = 2
    config["properties"]["grad"] = [1]
    report = _check(config, monkeypatch)
    assert not report.ok
    assert "public gradient slot 0" in report.to_text()


def test_checker_accepts_enabled_state_average_with_casscf_method(monkeypatch):
    config = _casscf_config()
    config["ci"]["nroot"] = 2
    config["state_average"] = {
        "enabled": True, "nstate": 2, "target_roots": [0, 1],
    }
    config["properties"]["grad"] = [1]
    report = _check(config, monkeypatch)
    assert report.ok, report.to_text()


def test_checker_rejects_json_orbitals_for_moving_casscf_workflow(monkeypatch):
    config = _casscf_config(runtype="optimize")
    config["cas"]["orbital_source"] = "json"
    config["cas"]["orbital_file"] = "central.json"
    config["cas"]["sort_orbitals"] = "none"
    report = _check(config, monkeypatch)
    assert not report.ok
    assert "cannot drive this moving casscf" in report.to_text()


def test_checker_rejects_unwired_casscf_workflows_and_bad_displacement(
        monkeypatch):
    for runtype in ("meci", "mecp", "neb", "hess"):
        report = _check(_casscf_config(runtype=runtype), monkeypatch)
        assert not report.ok
        assert "input.runtype" in report.to_text()

    config = _casscf_config(method="sa-casscf")
    config["casscf"]["grad_step"] = 0.0
    report = _check(config, monkeypatch)
    assert not report.ok
    assert "casscf.grad_step" in report.to_text()


class _Groups:
    ngroup = 1
    group = 0
    is_group_root = True

    def __init__(self, ntask):
        self.indices = range(ntask)

    @staticmethod
    def reduce_sum(value):
        return value

    @staticmethod
    def gather_list(value):
        return value


class _MPIManager:
    use_mpi = False

    @contextmanager
    def task_groups(self, ntask, **_kwargs):
        yield _Groups(ntask)

    @staticmethod
    def set_mpi_comm(_data):
        return None


def _load_wf_numgrad(monkeypatch):
    oqp = sys.modules.get("oqp", types.ModuleType("oqp"))
    utils = sys.modules.get("oqp.utils", types.ModuleType("oqp.utils"))
    library = sys.modules.get("oqp.library", types.ModuleType("oqp.library"))
    monkeypatch.setitem(sys.modules, "oqp", oqp)
    monkeypatch.setitem(sys.modules, "oqp.utils", utils)
    monkeypatch.setitem(sys.modules, "oqp.library", library)
    monkeypatch.setattr(oqp, "utils", utils, raising=False)
    monkeypatch.setattr(oqp, "library", library, raising=False)

    file_utils = types.ModuleType("oqp.utils.file_utils")
    file_utils.dump_log = lambda *_args, **_kwargs: None
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
    mpi_utils.MPIManager = _MPIManager
    single_point = types.ModuleType("oqp.library.single_point")
    single_point.SinglePoint = object
    monkeypatch.setitem(sys.modules, "oqp.utils.file_utils", file_utils)
    monkeypatch.setitem(sys.modules, "oqp.utils.mpi_utils", mpi_utils)
    monkeypatch.setitem(sys.modules, "oqp.library.single_point", single_point)

    name = "wf_numgrad_casscf_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/library/wf_numgrad.py")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


class _FakeMolecule:
    def __init__(self):
        self.config = {
            "input": {"method": "casscf"},
            "casscf": {"grad_step": 2.0e-4, "grad_guess": "cold"},
            "guess": {"type": "hcore"},
        }
        self.data = {"natom": 1}
        self.usempi = False
        self._x = np.array([0.3, -0.2, 0.1])
        self.energies = [self._energy()]

    def _energy(self):
        # The exact gradient is [1.2, -0.7, 0.4] at every geometry.
        return float(np.dot([1.2, -0.7, 0.4], self._x) - 5.0)

    def get_system(self):
        return self._x.copy()

    def update_system(self, coordinates):
        self._x = np.asarray(coordinates, dtype=float).copy()


class _FakeSinglePoint:
    def __init__(self, mol):
        self.mol = mol

    def energy(self, **_kwargs):
        self.mol.energies = [self.mol._energy()]
        return self.mol.energies


def test_explicit_casscf_fd_reference_restores_central_geometry(
        monkeypatch):
    numgrad = _load_wf_numgrad(monkeypatch)
    mol = _FakeMolecule()
    x0 = mol.get_system()
    gradients = numgrad.wavefunction_numerical_gradient(
        mol, [0], sp=_FakeSinglePoint(mol))

    np.testing.assert_allclose(gradients[0, 0], [1.2, -0.7, 0.4], atol=1.0e-11)
    np.testing.assert_array_equal(mol.get_system(), x0)
    assert mol.wf_numgrad_layout["ntask"] == 6
    assert not hasattr(mol, "pt2_numgrad_layout")
    assert "casscf" in numgrad.CASSCF_NUMGRAD_METHODS
    assert "sa-casscf" in numgrad.WF_NUMGRAD_METHODS
