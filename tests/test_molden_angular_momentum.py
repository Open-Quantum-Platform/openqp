"""The Molden writer must not be handed shells it cannot express.

The Molden [GTO] section defines shell labels only up to g, and the
Cartesian/spherical reorderings in MoldenWriter are tabulated to l=4 for the
same reason. A basis that reaches h -- cc-pV5Z and larger -- used to index
past the end of SHELL_TYPES and abort the whole calculation with an opaque
IndexError, after the SCF had already converged.

Two properties are pinned here: the writer names the problem instead of
crashing, and the detection is a cheap up-front check so a caller can skip
the file rather than leaving a truncated one behind.

Pure Python -- the Molden writer has no native dependency, so these run in a
source-only checkout.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "pyoqp" / "oqp" / "molden" / "moldenwriter.py"
SPEC = importlib.util.spec_from_file_location("_openqp_moldenwriter", MODULE_PATH)
moldenwriter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = moldenwriter
SPEC.loader.exec_module(moldenwriter)

MoldenWriter = moldenwriter.MoldenWriter


def _basis(angs):
    """Minimal basis dict with one uncontracted shell per angular momentum."""
    n = len(angs)
    return {
        "centers": [0] * n,
        "angs": list(angs),
        "ncontr": [1] * n,
        "alpha": [1.0] * n,
        "coef": [1.0] * n,
        "nbf": sum(2 * a + 1 for a in angs),
    }


def test_supported_basis_reports_no_problem():
    # s through g is the whole of what Molden defines.
    assert MoldenWriter.unsupported_angular_momentum(_basis([0, 1, 2, 3, 4])) is None
    assert MoldenWriter.max_angular_momentum(_basis([0, 1, 2, 3, 4])) == 4


def test_h_shell_is_reported_not_crashed_on():
    basis = _basis([0, 1, 2, 3, 4, 5])
    assert MoldenWriter.unsupported_angular_momentum(basis) == 5


def test_i_shell_reports_the_actual_offender():
    assert MoldenWriter.unsupported_angular_momentum(_basis([0, 6])) == 6


def test_empty_basis_is_not_a_failure():
    assert MoldenWriter.unsupported_angular_momentum(_basis([])) is None


def test_write_basis_raises_a_named_error_rather_than_indexerror(tmp_path):
    """Defensive: a caller that skipped the check gets a message, not IndexError."""
    out = tmp_path / "x.molden"
    with open(out, "w") as fh:
        writer = MoldenWriter(fh)
        with pytest.raises(ValueError, match="l=5"):
            writer.write_basis(1, _basis([0, 5]))


def test_write_basis_still_works_for_g(tmp_path):
    out = tmp_path / "g.molden"
    with open(out, "w") as fh:
        MoldenWriter(fh).write_basis(1, _basis([0, 4]))
    text = out.read_text()
    assert "[GTO]" in text
    # s and g shells both labelled
    assert "\ns 1" in text and "\ng 1" in text
