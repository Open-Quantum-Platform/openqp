"""Tagged ``basis=library`` geometry parsing regressions."""

from types import SimpleNamespace

import pytest

from oqp.library.set_basis import BasisData


def _basis_data(system, library):
    data = BasisData.__new__(BasisData)
    data.mol = SimpleNamespace(
        config={
            "input": {
                "basis": "library",
                "system": system,
                "library": library,
            }
        }
    )
    data.basis_names = []
    return data


def test_nonempty_tagged_inline_first_atom_is_not_treated_as_filename():
    data = _basis_data(
        "C 0.0 0.0 0.0 c1\nH 0.0 0.0 1.0 h1",
        "c1 cc-pvdz\nh1 6-31g*",
    )

    assert data.get_basislist() == ["cc-pvdz", "6-31g*"]


def test_legacy_file_backed_tagged_xyz_still_works(tmp_path):
    xyz = tmp_path / "tagged.xyz"
    xyz.write_text(
        "2\ntagged basis geometry\n"
        "C 0.0 0.0 0.0 c1\n"
        "H 0.0 0.0 1.0 h1\n",
        encoding="utf-8",
    )
    data = _basis_data(
        str(xyz),
        "c1 cc-pvdz\nh1 6-31g*",
    )

    assert data.get_basislist() == ["cc-pvdz", "6-31g*"]


def test_inline_atom_without_basis_tag_gets_tag_diagnostic():
    data = _basis_data(
        "C 0.0 0.0 0.0\nH 0.0 0.0 1.0 h1",
        "c1 cc-pvdz\nh1 6-31g*",
    )

    with pytest.raises(FileNotFoundError, match="correctly add a tag"):
        data.get_basislist()
