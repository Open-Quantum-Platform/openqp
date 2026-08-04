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
import os
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


def test_molecule_write_molden_clears_stale_output(tmp_path):
    """End-to-end on the real method, with the native layers stubbed out."""
    import types
    import warnings as _warnings
    from oqp.molecule.molecule import Molecule

    stale = tmp_path / "orbitals.molden"
    stale.write_text("[Molden Format]\nstale\n")

    mol = Molecule.__new__(Molecule)
    mol.usempi = False
    mol.data = types.SimpleNamespace(get_basis=lambda: _basis([0, 5]))

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        # Call the undecorated function; @mpi_dump only gates on rank.
        Molecule.write_molden.__wrapped__(mol, str(stale)) \
            if hasattr(Molecule.write_molden, "__wrapped__") \
            else Molecule.write_molden(mol, str(stale))

    assert not stale.exists(), "stale Molden file must be removed when skipping"
    assert any("Skipping Molden output" in str(w.message) for w in caught)


def test_molecule_write_molden_skips_any_unportable_basis(tmp_path):
    """The gate is supports_portable_ordering, not just the l ceiling.

    A mixed cartesian/spherical basis has every l <= g, so the angular-momentum
    check alone would wave it through -- but there is no portable reordering
    for it and the MO coefficients would come out silently permuted. The
    warning has no single offending l to name, so it must fall back to the
    generic reason rather than formatting ``l=None``.
    """
    import types
    import warnings as _warnings
    from oqp.molecule.molecule import Molecule

    mixed = _basis([0, 1, 2, 3, 4])
    mixed["nbf"] += 1          # matches neither the cartesian nor spherical count
    assert MoldenWriter.unsupported_angular_momentum(mixed) is None
    assert not MoldenWriter.supports_portable_ordering(mixed)

    out = tmp_path / "mixed.molden"
    mol = Molecule.__new__(Molecule)
    mol.usempi = False
    mol.data = types.SimpleNamespace(get_basis=lambda: mixed)

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        Molecule.write_molden.__wrapped__(mol, str(out)) \
            if hasattr(Molecule.write_molden, "__wrapped__") \
            else Molecule.write_molden(mol, str(out))

    messages = [str(w.message) for w in caught]
    assert not out.exists(), "an unportable basis must not produce a Molden file"
    assert any("no portable Molden ordering" in m for m in messages), messages
    assert not any("l=None" in m for m in messages), messages


def test_molecule_write_molden_survives_unremovable_stale_output(tmp_path, monkeypatch):
    """An undeletable stale file must not take the calculation down with it.

    Skipping exists so a by-product cannot cost the user the energy; an
    os.remove that raises (viewer holding the file open, read-only directory)
    would have propagated straight out of the SCF driver.
    """
    import types
    import warnings as _warnings
    from oqp.molecule.molecule import Molecule

    stale = tmp_path / "orbitals.molden"
    stale.write_text("[Molden Format]\nstale\n")

    def refuse(path):
        raise OSError(13, "Permission denied")

    monkeypatch.setattr(os, "remove", refuse)

    mol = Molecule.__new__(Molecule)
    mol.usempi = False
    mol.data = types.SimpleNamespace(get_basis=lambda: _basis([0, 5]))

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        Molecule.write_molden.__wrapped__(mol, str(stale)) \
            if hasattr(Molecule.write_molden, "__wrapped__") \
            else Molecule.write_molden(mol, str(stale))

    messages = [str(w.message) for w in caught]
    assert any("Skipping Molden output" in m for m in messages)
    assert any("Could not remove the stale Molden file" in m for m in messages)
    assert stale.exists(), "the file we could not remove is still there"
