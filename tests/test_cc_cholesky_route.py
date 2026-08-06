"""The Cholesky and integral-direct CC routes are pinned by the ROUTE, not the energy.

The two decks under ``examples/CC`` that exist for these controls cannot, on
their own, tell anyone whether the controls work.  Cholesky at ``tol=1e-10`` is
exact to about 1e-11 by construction -- reproducing the explicit ladder result
is the entire point of the factorisation -- while the example harness compares
energies with ``np.round(diff, 4) > 0`` (``compare_data`` in
``pyoqp/oqp/molecule/molecule.py``), an absolute gate near 5e-5 Ha.  Seven
orders of magnitude separate the two, so a regression that silently dropped
back to the conventional route would reproduce every committed reference.

That is not hypothetical: ``cholesky_direct=true`` was twice found to be
silently ignored during review of this branch -- once when ``cholesky`` was
left at ``auto``, once for open-shell references -- and both would have passed
the examples.

So assert on what actually distinguishes the routes, the log the driver emits,
and keep the control arm in the test itself: each deck is run with its control
on AND off, and the marker must appear in exactly one of them.  A test that
only ran the "on" arm would pass just as happily if the flag did nothing.

The energies are compared too, but as the physical invariant that matters --
the routes must agree with each other far more tightly than the example gate
demands -- rather than as the thing that detects the route.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).parents[1] / "examples" / "CC"
STORED = EXAMPLES / "h2o_ccsd_t_cholesky_sto-3g.inp"
DIRECT = EXAMPLES / "h2o_ccsd_t_cholesky_direct_sto-3g.inp"

# Emitted once the factorisation has run, whichever way it was produced.
FACTORISED = "Cholesky vectors"
# Emitted only by the integral-direct decomposition.
INTEGRAL_DIRECT = "direct factorisation used"

ENERGY = re.compile(r"E\(CCSD\(T\), total\)\s*=\s*(-?\d+\.\d+)")


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    import oqp
    return hasattr(oqp, "ccsd_t_energy")


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(),
    reason="native OpenQP library with ccsd_t_energy not importable",
)


def _run(tmp_path, name, deck_text):
    """Run a deck and return (log text, total CCSD(T) energy)."""
    inp = tmp_path / f"{name}.inp"
    inp.write_text(deck_text)
    for xyz in EXAMPLES.glob("*.xyz"):
        (tmp_path / xyz.name).write_bytes(xyz.read_bytes())

    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, cwd=str(tmp_path), env=dict(os.environ), timeout=1800,
    )
    log_path = inp.with_suffix(".log")
    log = log_path.read_text() if log_path.exists() else ""
    assert proc.returncode == 0, f"{name} did not run:\n{log[-3000:]}"

    found = ENERGY.findall(log)
    assert found, f"{name} produced no CCSD(T) total energy:\n{log[-3000:]}"
    return log, float(found[-1])


def test_cholesky_control_actually_selects_the_factorised_route(tmp_path):
    """[cc] cholesky must change the route, and the energy must not notice."""
    deck = STORED.read_text()
    assert "cholesky=true" in deck, "the example stopped setting the control"

    on_log, on_e = _run(tmp_path, "chol_on", deck)
    off_log, off_e = _run(tmp_path, "chol_off",
                          deck.replace("cholesky=true", "cholesky=false"))

    assert FACTORISED in on_log, \
        "cholesky=true did not factorise -- the control is being ignored"
    assert FACTORISED not in off_log, \
        "cholesky=false still factorised -- the control is being ignored"

    # The routes are two ways to build the same integrals, so this is a real
    # equivalence, not a tolerance to be relaxed if it ever starts failing.
    assert abs(on_e - off_e) < 1e-8, (on_e, off_e)


def test_cholesky_direct_control_actually_skips_the_stored_integrals(tmp_path):
    """[cc] cholesky_direct must select the integral-direct decomposition.

    Distinct from the test above: the direct route also factorises, so
    ``Cholesky vectors`` alone cannot tell the two apart.  It is the
    integral-pass line that is unique to the direct decomposition.
    """
    deck = DIRECT.read_text()
    assert "cholesky_direct=true" in deck, "the example stopped setting the control"

    on_log, on_e = _run(tmp_path, "direct_on", deck)
    off_log, off_e = _run(
        tmp_path, "direct_off",
        deck.replace("cholesky_direct=true", "cholesky_direct=false"))

    assert INTEGRAL_DIRECT in on_log, \
        "cholesky_direct=true did not take the direct route"
    assert INTEGRAL_DIRECT not in off_log, \
        "cholesky_direct=false still took the direct route"

    assert abs(on_e - off_e) < 1e-8, (on_e, off_e)
