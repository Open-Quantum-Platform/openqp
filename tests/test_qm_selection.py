"""One definition of the QM region must be enough for a QM/MM job.

``[qmmm] qm_atoms`` (0-based, OpenMM topology) and the indices after the PDB
name in ``[input] system`` (1-based, the QM molecule) are the same atoms.
Requiring both invited a silent mismatch, and leaving the second out stopped
the run with "Atom list not defined!".
"""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "_qm_selection_direct", ROOT / "pyoqp" / "oqp" / "utils" / "qm_selection.py")
qm_selection = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = qm_selection
_SPEC.loader.exec_module(qm_selection)
reconcile = qm_selection.reconcile_qm_selection


def test_qm_atoms_alone_builds_the_qm_molecule_list():
    assert reconcile("ala.pdb", "8,9,16,17,18") == (
        "ala.pdb 9 10 17 18 19", "8,9,16,17,18")
    # ranges, spaces, lists and unsorted input all work; the result is sorted
    assert reconcile("/a/b/w.pdb", "0-2")[0] == "/a/b/w.pdb 1 2 3"
    assert reconcile("w.pdb", "0-10, 22 21")[0] == (
        "w.pdb 1 2 3 4 5 6 7 8 9 10 11 22 23")
    assert reconcile("w.pdb", [2, 0, 1])[0] == "w.pdb 1 2 3"


def test_array_and_scalar_selections_are_accepted():
    """A script may pass qm_atoms as a NumPy array (QMMM_MD accepted that
    before the reconciliation existed); stringifying it gave '[0 1]'."""
    import numpy as np
    assert reconcile("w.pdb", np.array([0, 1, 2]))[0] == "w.pdb 1 2 3"
    assert reconcile("w.pdb", (2, 0, 1))[0] == "w.pdb 1 2 3"
    assert reconcile("w.pdb", np.int64(4))[0] == "w.pdb 5"
    assert reconcile("w.pdb", 4)[0] == "w.pdb 5"
    array = np.array([0, 1])
    assert reconcile("w.pdb 1 2", array)[1] is array        # passed through
    with pytest.raises(ValueError, match="two lists differ"):
        reconcile("w.pdb 2 3", np.array([0, 1]))


def test_indices_after_the_pdb_alone_give_qm_atoms():
    assert reconcile("ala.pdb 9 10 17 18 19", "") == (
        "ala.pdb 9 10 17 18 19", "8 9 16 17 18")
    assert reconcile("ala.pdb 9-10 17 18-19", None)[1] == "8 9 16 17 18"
    with pytest.raises(ValueError, match="1-based"):
        reconcile("ala.pdb 0 1 2", "")


def test_both_given_must_name_the_same_atoms():
    # agreeing lists pass through untouched, in any order
    assert reconcile("ala.pdb 9 10 17 18 19", "8,9,16,17,18") == (
        "ala.pdb 9 10 17 18 19", "8,9,16,17,18")
    assert reconcile("ala.pdb 19 18 17 10 9", "8-9 16-18")[0] == "ala.pdb 19 18 17 10 9"
    # the classic slip: the same numbers in both places
    with pytest.raises(ValueError, match="defined twice and the two lists differ"):
        reconcile("ala.pdb 8 9 16 17 18", "8,9,16,17,18")


def test_no_system_at_all_falls_back_to_the_qmmm_pdb_file():
    # legacy decks that give only [qmmm] pdb_file and qm_atoms
    assert reconcile(None, "0-2", pdb_file="water_dimer.pdb") == (
        "water_dimer.pdb 1 2 3", "0-2")
    assert reconcile("", "0-2", pdb_file="water_dimer.pdb")[0] == "water_dimer.pdb 1 2 3"
    assert reconcile(None, "0-2") == (None, "0-2")          # nothing to build from


def test_other_geometries_are_left_alone():
    # an xyz file or an inline table carries its own atoms (QM/MM NAMD of a
    # whole-molecule QM region); neither list is touched
    assert reconcile("solute.xyz", "0-3") == ("solute.xyz", "0-3")
    inline = "6 0.0 0.0 0.0\n8 0.0 0.0 1.2"
    assert reconcile(inline, "0-1") == (inline, "0-1")
    assert reconcile("ala.pdb", "") == ("ala.pdb", "")        # neither given


def test_both_entry_points_reconcile_before_the_molecule_is_built():
    molecule = (ROOT / "pyoqp" / "oqp" / "molecule" / "molecule.py").read_text()
    load = molecule[molecule.index("    def load_config(self, input_source):"):]
    assert load.index("self._reconcile_qm_selection()") < load.index(
        "self.data.apply_config(self.config)")
    driver = (ROOT / "pyoqp" / "oqp" / "library" / "qmmm_md.py").read_text()
    assert driver.index("reconcile_qm_selection(") < driver.index(
        "qmmm_cfg, qm_cfg = _extract_qmmm_config(oqp_cfg=oqp_cfg, mol=mol)")


@pytest.mark.parametrize("deck", sorted(
    path.name for path in (ROOT / "examples" / "QMMM").glob("*.oqp")))
def test_shipped_decks_define_the_qm_region_consistently(deck):
    """Every QM/MM deck either gives qm_atoms alone or two lists that agree."""
    text = (ROOT / "examples" / "QMMM" / deck).read_text()
    import re
    geom = re.search(r'geom="([^"]*)"', text)
    atoms = re.search(r'qm_atoms="([^"]*)"', text)
    if geom is None or atoms is None:
        pytest.skip("no QM selection in this deck")
    reconcile(geom.group(1), atoms.group(1))      # raises on a mismatch
