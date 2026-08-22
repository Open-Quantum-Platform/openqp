"""A gradient belongs to the run that produced it (PR #367 review, P2).

``data._data.grad`` is allocated once and never cleared, so "is there a
gradient?" is answered by a marker rather than by the buffer's contents.  The
marker was only ever set.  The legacy ``OPENQP`` API runs several calculations
through one ``Runner``/``Molecule`` (``op.run("grad")`` then
``op.run("energy")``), and there the marker left standing makes an energy-only
run republish the previous run's derivatives as its own.
"""
import pytest


def _molecule_class():
    try:
        from oqp.molecule.molecule import Molecule
    except Exception:                        # pragma: no cover - no runtime
        pytest.skip("liboqp is not available")
    return Molecule


def test_marker_round_trip():
    molecule = object.__new__(_molecule_class())
    assert molecule.has_grad() is False      # never written
    molecule.mark_grad_valid()
    assert molecule.has_grad() is True
    molecule.invalidate_grad()
    assert molecule.has_grad() is False


def test_runner_invalidates_the_marker_before_every_run():
    """The reset must happen where all three entry points meet.

    A stub Runner is driven through ``Runner.run`` far enough to dispatch; the
    dispatch raises, which both stops the calculation and records that the
    marker was already cleared by then.
    """
    try:
        from oqp.pyoqp import Runner
    except Exception:                        # pragma: no cover - no runtime
        pytest.skip("liboqp is not available")

    seen = {}

    class _Sentinel(Exception):
        pass

    class _Mol:
        def __init__(self):
            self.config = {"input": {"runtype": "energy"}, "tests": {}}
            self._grad_valid = True          # left over from a gradient run

        def invalidate_grad(self):
            self._grad_valid = False

    runner = object.__new__(Runner)
    runner.mol = _Mol()

    def _dispatch(mol):
        seen["grad_valid_at_dispatch"] = mol._grad_valid
        raise _Sentinel

    runner.run_func = {"energy": _dispatch}

    with pytest.raises(_Sentinel):
        runner.run()

    assert seen["grad_valid_at_dispatch"] is False, (
        "an energy run reached its driver still carrying the previous run's "
        "gradient marker"
    )
