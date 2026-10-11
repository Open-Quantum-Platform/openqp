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
import re
import subprocess
import sys
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
    # H2 and HeH+ in STO-3G have one occupied and one virtual orbital, so the
    # response space holds exactly one root; TDHF refuses nstate > nocc*nvir.
    nstate: int = 2
    # cc-pVDZ/cc-pVTZ carry spherical d (and f) shells, so nbf differs from
    # the Cartesian extent and an AO-indexing slip in the response arrays shows.
    basis: str = "sto-3g"


CASES = {
    "h2": Case(
        system="""\
  H  0.0000000000  0.0000000000 -0.3700000000
  H  0.0000000000  0.0000000000  0.3700000000""",
        charge=0,
        # Re-validated after the minimal-basis shortcuts were removed: the old
        # values compared two wrong Hessians (zz 0.2016 against 1.159 from
        # the state-1 energies).  chc r630, GCC 12.3 + MKL ILP64.
        validated_max_abs=8.281666865617154e-7,
        validated_rms=3.361459712471747e-7,
        nstate=1,
    ),
    "hehplus": Case(
        system="""\
  He  0.0000000000  0.0000000000 -0.4000000000
  H   0.0000000000  0.0000000000  0.4000000000""",
        charge=1,
        validated_max_abs=1.1985536660930052e-6,
        validated_rms=4.3440667358826373e-7,
        nstate=1,
    ),
    # Linear H3+: two response roots.  Its numerical Hessian was 1.4e-3 off
    # (zz 0.26306 vs 0.26447) because the displaced gradients met only the
    # default squared-residual Z-vector tolerance; see
    # _finite_difference_response_tolerances in single_point.py.
    "h3plus": Case(
        system="""\
  H  0.0000000000  0.0000000000 -0.8500000000
  H  0.0000000000  0.0000000000  0.0000000000
  H  0.0000000000  0.0000000000  0.8500000000""",
        charge=1,
        validated_max_abs=1.4069516913323143e-6,
        validated_rms=2.149666497298211e-7,
    ),
    # Spherical bases: 5d/7f against 6d/10f Cartesian extents.
    "h2o_ccpvdz": Case(
        system="""\
  O  0.0000000000  0.0000000000  0.0000000000
  H  0.0000000000  0.0000000000  0.9000000000
  H  0.0000000000  0.7000000000 -0.3000000000""",
        charge=0,
        validated_max_abs=1.0294054729798052e-5,
        validated_rms=2.4847766312638596e-6,
        basis="cc-pvdz",
    ),
    "h2o_ccpvtz": Case(
        system="""\
  O  0.0000000000  0.0000000000  0.0000000000
  H  0.0000000000  0.0000000000  0.9000000000
  H  0.0000000000  0.7000000000 -0.3000000000""",
        charge=0,
        validated_max_abs=4.798990391685631e-6,
        validated_rms=1.7211117762219943e-6,
        basis="cc-pvtz",
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
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    try:
        from oqp.pyoqp import Runner  # noqa: F401
    except Exception:
        # Setting OPENQP_RUN_TDHF_HESSIAN_REGRESSION=1 is an explicit request
        # to run the live comparison, and CI sets it. Swallowing an import
        # failure here would report every Hessian comparison as skipped and
        # leave the suite green with the gate never having run -- the exact
        # failure mode this gate exists to prevent. Fail loudly instead; the
        # same reasoning as OQP_REQUIRE_NATIVE_TESTS in .github/workflows/CI.yml.
        raise
    return True


RUNTIME_AVAILABLE = _runtime_available()


INPUT = """[input]
system=
{system}
charge={charge}
method=tdhf
basis={basis}
runtype=hess
functional=

[scf]
type=rhf
multiplicity=1
conv=1.0e-10

[tdhf]
type=rpa
multiplicity=1
nstate={nstate}
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
            basis=case.basis,
            system=case.system,
            charge=case.charge,
            nstate=case.nstate,
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
def test_live_dft_hessian_consumers_stay_inside_their_probed_blocks(tmp_path, monkeypatch):
    """Pin the block bookkeeping in explicit_channel_derivative_matrix.

    The XC polarization loop probes only the MO block each consumer declares it
    reads (occ-virt for the amplitude RHS, occ-occ|virt-virt for the Z RHS,
    occ-occ for the relaxed density). With OQP_TDHESS_POISON_UNPROBED=1 every
    unprobed element is NaN, so a consumer that reads outside its declared
    block cannot produce a finite Hessian. A pure-TDHF case would not exercise
    this: the XC loop only runs with a functional.
    """
    monkeypatch.setenv("OPENQP_ROOT", os.environ.get("OPENQP_ROOT", str(ROOT)))
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    monkeypatch.setenv("OQP_TDHESS_POISON_UNPROBED", "1")
    from oqp.pyoqp import Runner

    workdir = tmp_path / "poison"
    workdir.mkdir()
    inp = workdir / "svwn.inp"
    inp.write_text(
        INPUT.format(system=CASES["h2o"].system, charge=0, nstate=CASES["h2o"].nstate, basis="sto-3g",
                     hess_type="analytical")
        .replace("functional=\n", "functional=svwn\n")
    )
    assert "functional=svwn" in inp.read_text()
    runner = Runner(project="svwn", input_file=str(inp), log=str(workdir / "svwn.log"))
    runner.run()
    hessian = np.asarray(runner.mol.hessian, dtype=float)
    assert hessian.shape == (9, 9)
    assert np.all(np.isfinite(hessian)), "a consumer read an unprobed (NaN) block"


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
    # Respect an OPENQP_ROOT supplied by the environment: CI points it at the
    # *installed* package directory (see .github/workflows/CI.yml), which is
    # where liboqp lives there. Hard-overriding it to the source tree, which
    # carries no lib/, would break the run in exactly the job that can execute
    # it. Falling back to the repo root keeps source-tree runs working. This
    # matches tests/test_cam_hessian.py and tests/test_rhf_hessian_dfunc.py.
    monkeypatch.setenv("OPENQP_ROOT", os.environ.get("OPENQP_ROOT", str(ROOT)))
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


@pytest.mark.skipif(
    not RUNTIME_AVAILABLE,
    reason=(
        "set OPENQP_RUN_TDHF_HESSIAN_REGRESSION=1 and use a built OpenQP "
        "runtime to run TDHF Hessians"
    ),
)
def test_live_analytic_hessian_refuses_one_root_of_a_larger_space(tmp_path):
    """nstate=1 is enough only when the response space holds one root.

    H2O/STO-3G has 10 occupied-virtual pairs; with one computed root the
    isolation of the lowest root cannot be checked.  The input checker
    counts the pairs from the basis and refuses the run before any SCF work
    (tdhf_hessian_lowest_root_is_isolated is the runtime guard behind it).
    """
    inp = tmp_path / "h2o.inp"
    inp.write_text(INPUT.format(system=CASES["h2o"].system, charge=0, nstate=1, basis="sto-3g",
                                hess_type="analytical"))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.setdefault("OPENQP_ROOT", str(ROOT))
    proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)], cwd=tmp_path,
                          env=env, capture_output=True, text=True, timeout=600)
    log = (tmp_path / "h2o.log").read_text() if (tmp_path / "h2o.log").exists() else ""
    text = log + proc.stdout + proc.stderr
    assert proc.returncode != 0, text[-2000:]
    # Refused by the input check, with a message that reaches the user.
    assert "tdhf.nstate>=2" in text, text[-2000:]
    assert "Final RHF energy" not in text
    assert not (tmp_path / "fort.6").exists()


@pytest.mark.skipif(
    not RUNTIME_AVAILABLE,
    reason=(
        "set OPENQP_RUN_TDHF_HESSIAN_REGRESSION=1 and use a built OpenQP "
        "runtime to run TDHF Hessians"
    ),
)
def test_runtime_isolation_guard_message_reaches_stderr(tmp_path):
    """The native guard behind the input check must abort visibly.

    The input check is bypassed (as an API caller that skips it would) so the
    run reaches tdhf_hessian with one root of a 10-pair space.  show_message
    then writes while the log unit is closed; it used to land in a stray
    fort.6 file in the working directory, now it goes to stderr.
    """
    inp = tmp_path / "h2o.inp"
    inp.write_text(INPUT.format(system=CASES["h2o"].system, charge=0, nstate=1, basis="sto-3g",
                                hess_type="analytical"))
    script = (
        "import sys, oqp.pyoqp as P\n"
        "P.check_input_values = lambda *a, **k: None\n"
        "P.Runner(project='h2o', input_file=sys.argv[1], log=sys.argv[2]).run()\n"
    )
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.setdefault("OPENQP_ROOT", str(ROOT))
    proc = subprocess.run([sys.executable, "-c", script, str(inp), str(tmp_path / "h2o.log")],
                          cwd=tmp_path, env=env, capture_output=True, text=True, timeout=600)
    assert proc.returncode != 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "requires an isolated lowest root" in proc.stderr, proc.stderr[-2000:]
    assert not (tmp_path / "fort.6").exists()


ENERGY_INPUT = """[input]
system=
{system}
charge={charge}
method=tdhf
basis={basis}
runtype=energy
functional=

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
conv=1.0e-11

[tdhf]
type=rpa
multiplicity=1
nstate={nstate}
conv=1.0e-11
"""

STATE1 = re.compile(r"PyOQP state 1\s+(-?\d+\.\d+)")


def _state1_energy(case, system, workdir):
    workdir.mkdir(parents=True)
    inp = workdir / "e.inp"
    inp.write_text(ENERGY_INPUT.format(system=system, charge=case.charge, nstate=case.nstate,
                                         basis=case.basis))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)], cwd=workdir,
                          env=env, capture_output=True, text=True, timeout=600)
    log = (workdir / "e.log").read_text() if (workdir / "e.log").exists() else ""
    assert proc.returncode == 0, log[-2000:] + proc.stderr[-2000:]
    return float(STATE1.findall(log)[0])


@pytest.mark.skipif(
    not RUNTIME_AVAILABLE,
    reason=(
        "set OPENQP_RUN_TDHF_HESSIAN_REGRESSION=1 and use a built OpenQP "
        "runtime to run TDHF Hessians"
    ),
)
@pytest.mark.parametrize("name", CASES)
def test_live_analytic_hessian_matches_state_energy_curvature(name, tmp_path, monkeypatch):
    """Pin one diagonal element against independent state-1 energies.

    The analytic and numerical Hessians share the SCF/TD machinery, so they
    can agree on a wrong answer: for H2 both read stale data (0.20 Eh/bohr^2
    analytic, 0.34 numerical) while the state-1 energies, each from a fresh
    Huckel start, curve with 1.159 Eh/bohr^2 along the bond.  Displace the
    last atom along z by +-h and compare the second difference of the
    state-1 total energy with the analytic zz element of that atom.
    """
    monkeypatch.setenv("OPENQP_ROOT", os.environ.get("OPENQP_ROOT", str(ROOT)))
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    case = CASES[name]
    analytic = _run_hessian(case, "analytical", tmp_path / "analytic")
    lines = case.system.splitlines()
    sym, x, y, z = lines[-1].split()
    h_ang = 0.003
    energies = []
    for k, dz in enumerate((-h_ang, 0.0, h_ang)):
        moved = lines[:-1] + [f"  {sym}  {float(x):.10f}  {float(y):.10f}  {float(z) + dz:.10f}"]
        energies.append(_state1_energy(case, "\n".join(moved), tmp_path / f"e{k}"))
    h_bohr = h_ang / 0.529177210903
    curvature = (energies[0] - 2.0 * energies[1] + energies[2]) / h_bohr**2
    assert analytic[-1, -1] == pytest.approx(curvature, abs=2.0e-4)
