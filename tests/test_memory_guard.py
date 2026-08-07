"""The pre-flight memory guard refuses jobs that do not fit this machine.

The point of the guard is that the ceiling is the machine, not a constant
compiled into the source: the same binary runs on a laptop and on a 500 GB
node, and a fixed cap is wrong in both directions. These tests pin the two
halves of that -- it refuses when the budget is too small, and it runs the
same job when the budget is large -- by driving OQP_MEMORY_LIMIT_GB, which
overrides the probe.

Skipped unless a built OpenQP is importable, like the other end-to-end tests.
"""

import os
import subprocess
import sys

import pytest


def _have_native_oqp():
    try:
        import oqp
    except Exception:
        return False
    return hasattr(oqp, "ccsd_t_energy")


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(),
    reason="native OpenQP library with ccsd_t_energy not importable",
)

INPUT = """[input]
system=
 8   0.000000000   0.000000000  -0.041061554
 1  -0.533194329   0.533194329  -0.614469223
 1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g
method=ccsd(t)
functional=

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
conv=1.0e-9

[cc]
nfzc=1
"""


def _run(tmp_path, limit_gb):
    inp = tmp_path / "mem.inp"
    inp.write_text(INPUT)
    log = inp.with_suffix(".log")
    if log.exists():
        log.unlink()
    env = dict(os.environ, OQP_MEMORY_LIMIT_GB=str(limit_gb))
    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, cwd=str(tmp_path), env=env, timeout=1800,
    )
    return proc, (log.read_text() if log.exists() else "")


def test_refuses_when_the_budget_is_smaller_than_the_estimate(tmp_path):
    proc, log = _run(tmp_path, 0.0001)

    assert proc.returncode != 0, "a job that cannot fit must not run"
    assert "needs about" in log, log[-2000:]
    assert "is available on this machine" in log
    # The two numbers are the whole point; they must survive formatting
    # rather than both rendering as "0.00 GB".
    assert "0.00 GB but only 0.00 GB" not in log
    # And it must say what the user can do about it.
    assert "freeze more core orbitals" in log
    assert "OQP_MEMORY_LIMIT_GB" in log


def test_same_job_runs_when_the_budget_is_generous(tmp_path):
    proc, log = _run(tmp_path, 1000)

    assert "estimated peak memory" in log
    assert "needs about" not in log, "must not refuse a job that fits"
    assert "E(CCSD(T), total)" in log, log[-2000:]


def test_estimate_is_reported_even_on_a_run_that_proceeds(tmp_path):
    """A run that fits should still record what it expected to need."""
    _, log = _run(tmp_path, 1000)
    line = [ln for ln in log.splitlines() if "estimated peak memory" in ln]
    assert line, log[-2000:]
    assert "available" in line[0]


def test_budget_says_whether_resident_memory_still_counts(monkeypatch):
    """The two budget sources mean different things, and callers depend on it.

    The live probe reports what REMAINS -- MemAvailable, or a cgroup limit
    minus current usage -- so memory this process already holds has been
    subtracted from it.  A caller sizing a future peak that includes an
    already-allocated array (the open-shell CC guard charges the packed AO
    store) must not charge those bytes twice against such a budget.
    OQP_MEMORY_LIMIT_GB is the opposite: a flat number for the whole job with
    nothing subtracted, so there the full peak is the right comparison.

    Getting this backwards is silent in both directions -- refusing jobs that
    fit, or admitting ones that do not -- so pin which is which.
    """
    import oqp

    monkeypatch.delenv("OQP_MEMORY_LIMIT_GB", raising=False)
    assert oqp.lib.oqp_memory_budget_includes_resident() == 0, \
        "the live probe already excludes resident memory"

    monkeypatch.setenv("OQP_MEMORY_LIMIT_GB", "64")
    assert oqp.lib.oqp_memory_budget_includes_resident() == 1, \
        "an explicit budget still has to cover memory already allocated"

    # Surrounding whitespace is still a number.
    monkeypatch.setenv("OQP_MEMORY_LIMIT_GB", "  64  ")
    assert oqp.lib.oqp_memory_budget_includes_resident() == 1

    # Unparseable values fall through to the probe, so they follow its meaning.
    # "1O" is the one that matters: a mistyped 10 with a letter O. strtod stops
    # at the O and reports success for the "1", so a prefix-only parse silently
    # imposed a 1 GB budget and aborted jobs against a limit nobody set. The
    # whole string has to be a number.
    for bad in ("not-a-number", "-8", "0", "1O", "8GB", "64 GB", "1e400", ""):
        monkeypatch.setenv("OQP_MEMORY_LIMIT_GB", bad)
        assert oqp.lib.oqp_memory_budget_includes_resident() == 0, bad
