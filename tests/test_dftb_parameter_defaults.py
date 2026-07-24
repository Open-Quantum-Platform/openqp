"""Source-only tests for the bundled OpenQP-DFTB parameter default."""

import importlib.util
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_input_checker(monkeypatch):
    oqp = types.ModuleType("oqp")
    utils = types.ModuleType("oqp.utils")
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        use_mpi = False
        size = 1

    mpi_utils.MPIManager = MPIManager
    monkeypatch.setitem(sys.modules, "oqp", oqp)
    monkeypatch.setitem(sys.modules, "oqp.utils", utils)
    monkeypatch.setitem(sys.modules, "oqp.utils.mpi_utils", mpi_utils)

    name = "_openqp_dftb_default_checker"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/utils/input_checker.py"
    )
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


def _dftb_config():
    return {
        "input": {"method": "dftb", "runtype": "energy"},
        "dftb": {"backend": "native", "type": "ground"},
        "properties": {"scf_prop": []},
        "pcm": {"enabled": False},
    }


def test_checker_accepts_omitted_parameter_path_when_bundle_exists(
    monkeypatch,
):
    checker = _load_input_checker(monkeypatch)
    package = types.ModuleType("openqp_dftb")
    package.default_parameter_path = lambda: "/bundle/OB2W0PT3"
    monkeypatch.setitem(sys.modules, "openqp_dftb", package)
    monkeypatch.delenv("OPENQP_DFTB_PARAMETER_PATH", raising=False)

    report = checker.CheckReport()
    checker._check_dftb(_dftb_config(), report)

    assert not any(item.path == "dftb.parameter_path" for item in report.errors)


def test_checker_reports_missing_bundle_and_override(monkeypatch):
    checker = _load_input_checker(monkeypatch)
    package = types.ModuleType("openqp_dftb")
    monkeypatch.setitem(sys.modules, "openqp_dftb", package)
    monkeypatch.delenv("OPENQP_DFTB_PARAMETER_PATH", raising=False)

    report = checker.CheckReport()
    checker._check_dftb(_dftb_config(), report)

    assert any(item.path == "dftb.parameter_path" for item in report.errors)
