"""[dftgrid] xi_* validation: c = p - alpha must stay in the validated kernel range.

The xi^alpha kernel (source/dftlib/xi_alpha_kernel.F90) is validated for
0 <= c <= 2; for very negative alpha the Gauss-Jacobi normalization overflows
while 2^(-c)/Gamma(c) underflows, giving NaN weights.  The input checker must
reject such orders before the native kernel runs.  Pure Python: executes the
production checker source without the compiled runtime.
"""
import importlib.util
from pathlib import Path
import sys
import types

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def checker(monkeypatch):
    mpi = types.ModuleType("oqp.utils.mpi_utils")
    mpi.MPIManager = lambda: types.SimpleNamespace(use_mpi=False, size=1)
    monkeypatch.setitem(sys.modules, "oqp.utils.mpi_utils", mpi)
    spec = importlib.util.spec_from_file_location(
        "_xi_checker", ROOT / "pyoqp/oqp/utils/input_checker.py")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, "_xi_checker", module)
    spec.loader.exec_module(module)
    return module


def _errors(checker, alpha, p):
    config = {
        "input": {"runtype": "energy", "method": "hf", "functional": "m06-2x"},
        "dftgrid": {"xi_mode": 1, "xi_alpha": alpha, "xi_p": p, "xi_cutoff": 0.0},
    }
    report = checker.CheckReport()
    checker._check_xi(config, report)
    return [d for d in report.diagnostics if d.severity == "ERROR"]


@pytest.mark.parametrize("alpha,p", [(-2.0, -1), (-2.0, 0), (-1.0, 1), (-0.5, -1), (0.5, 1), (1.0, -1)])
def test_validated_orders_accepted(checker, alpha, p):
    assert _errors(checker, alpha, p) == []


@pytest.mark.parametrize("alpha,p", [(-1100.0, -1), (-2.5, 0), (-2.0001, -1), (-1.5, 1)])
def test_orders_beyond_kernel_range_rejected(checker, alpha, p):
    errors = _errors(checker, alpha, p)
    assert any("<= 2" in d.message for d in errors), errors
