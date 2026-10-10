"""Locate the openqp working tree that the source-contract tests read.

The tests next to this file assert on the *shape* of engine source -- that a
guard precedes a division, that a subroutine does not carry state between
calls, that a loop nest is ordered a particular way.  They are development
gates, not engine regressions: none of them runs a calculation, and one of them
failed on a pure hoist that changed no arithmetic.  They were therefore moved
out of openqp's own suite (openqp#405) and into this repository.

In the private superset the engine is the repository root. A standalone
historical devkit checkout can still point ``OPENQP_SOURCE_ROOT`` at an engine
checkout or keep the two checkouts side by side::

    parent/
      openqp/
      openqp-devkit/

Without either, the tests skip rather than fail: there is nothing to assert
against, which is not the same as an assertion being violated.
"""

import os
from pathlib import Path

import pytest

# Any openqp checkout has this file; a directory that does not is not one.
_MARKER = Path("source") / "modules" / "tdhf_mrsf_gradient.F90"


def _candidates(helper_path=None):
    helper_path = Path(helper_path or __file__).resolve()
    env = os.environ.get("OPENQP_SOURCE_ROOT")
    if env:
        yield Path(env).expanduser()
    # devkit/tests/openqp_checkout.py -> private OpenQP repository root.
    yield helper_path.parents[2]
    # Compatibility with a standalone historical openqp-devkit checkout.
    yield helper_path.parents[2] / "openqp"


def find_openqp_root(helper_path=None):
    """Return the openqp checkout to read, or None if there is not one."""
    for candidate in _candidates(helper_path):
        if (candidate / _MARKER).is_file():
            return candidate.resolve()
    return None


def require_openqp():
    """Return the openqp checkout, skipping the whole module if absent."""
    root = find_openqp_root()
    if root is None:
        pytest.skip(
            "these gates read openqp engine source; set OPENQP_SOURCE_ROOT to "
            "an openqp checkout, or use the integrated private superset",
            allow_module_level=True,
        )
    return root
