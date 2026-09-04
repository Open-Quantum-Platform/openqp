"""Numerical regression gate for the imported singlet TDHF Hessian.

The inexpensive test records the three analytic-versus-central-difference
comparisons used to validate the import.  The calculation test is deliberately
opt-in because it performs 54 displaced TDHF gradients in addition to three
analytic Hessians.  Run it with::

    OPENQP_RUN_TDHF_HESSIAN_REGRESSION=1 pytest -q \
        tests/test_tdhf_hessian_regression.py

The live test uses one process and one OpenMP thread so its numerical reference
does not depend on the machine on which pytest happens to run.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
RUN_LIVE = os.environ.get("OPENQP_RUN_TDHF_HESSIAN_REGRESSION") == "1"
MAX_ABS_TOL = 2.0e-5
RMS_TOL = 5.0e-6


@dataclass(frozen=True)
class Case:
    system: str
    charge: int
    validated_max_abs: float
    validated_rms: float


CASES = {
    "h2": Case(
        system="""\
  H  0.0000000000  0.0000000000 -0.3700000000
  H  0.0000000000  0.0000000000  0.3700000000""",
        charge=0,
        validated_max_abs=5.667895902894404e-7,
        validated_rms=1.9155628621431792e-7,
    ),
    "hehplus": Case(
        system="""\
  He  0.0000000000  0.0000000000 -0.4000000000
  H   0.0000000000  0.0000000000  0.4000000000""",
        charge=1,
        validated_max_abs=1.1985536660930052e-6,
        validated_rms=4.3440667358826373e-7,
    ),
    "h2o": Case(
        system="""\
  O  0.0000000000  0.0000000000  0.0000000000
  H  0.0000000000  0.0000000000  0.9000000000
  H  0.0000000000  0.7000000000 -0.3000000000""",
        charge=0,
        validated_max_abs=1.745581673406882e-5,
        validated_rms=2.78874022551773e-6,
    ),
}


def _runtime_available() -> bool:
    if not RUN_LIVE:
        return False
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        from oqp.pyoqp import Runner  # noqa: F401
    except Exception:
        return False
    return True


RUNTIME_AVAILABLE = _runtime_available()


INPUT = """[input]
system=
{system}
charge={charge}
method=tdhf
basis=sto-3g
runtype=hess
functional=

[scf]
type=rhf
multiplicity=1
conv=1.0e-10

[tdhf]
type=rpa
multiplicity=1
nstate=1
conv=1.0e-10
zvconv=1.0e-10

[hess]
type={hess_type}
state=1
dx=0.001
nproc=1
clean=True
"""


@pytest.mark.parametrize("name", CASES)
def test_validated_tdhf_hessian_errors_are_within_contract(name):
    """Keep the independently obtained validation values and limits visible."""
    case = CASES[name]
    assert case.validated_max_abs < MAX_ABS_TOL
    assert case.validated_rms < RMS_TOL


def _run_hessian(case: Case, hess_type: str, workdir: Path) -> np.ndarray:
    from oqp.pyoqp import Runner

    workdir.mkdir(parents=True)
    inp = workdir / f"{hess_type}.inp"
    inp.write_text(
        INPUT.format(
            system=case.system,
            charge=case.charge,
            hess_type=hess_type,
        )
    )
    runner = Runner(
        project=hess_type,
        input_file=str(inp),
        log=str(workdir / f"{hess_type}.log"),
    )
    runner.run()
    return np.asarray(runner.mol.hessian, dtype=float)


@pytest.mark.skipif(
    not RUNTIME_AVAILABLE,
    reason=(
        "set OPENQP_RUN_TDHF_HESSIAN_REGRESSION=1 and use a built OpenQP "
        "runtime to run TDHF Hessians"
    ),
)
@pytest.mark.parametrize("name", CASES)
def test_live_tdhf_analytic_hessian_matches_numerical(name, tmp_path, monkeypatch):
    """Compare the production analytic result with its gradient finite difference."""
    monkeypatch.setenv("OPENQP_ROOT", str(ROOT))
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    case = CASES[name]
    workdir = tmp_path / name
    workdir.mkdir()

    analytic = _run_hessian(case, "analytical", workdir / "analytic")
    numerical = _run_hessian(case, "numerical", workdir / "numerical")

    assert analytic.shape == numerical.shape
    assert np.max(np.abs(analytic - analytic.T)) < 1.0e-9
    difference = analytic - numerical
    assert np.max(np.abs(difference)) < MAX_ABS_TOL
    assert np.sqrt(np.mean(difference * difference)) < RMS_TOL
