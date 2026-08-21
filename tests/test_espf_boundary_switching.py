"""Boundary-aware ESPF_SWSCALE selection (issue #260).

The switching width lives in the native ESPF gradient, which is handed a grid
and a density and cannot tell whether the QM/MM partition cut a bond.  The
driver can -- it just built the link atoms -- so the choice is made there.

1.8 over-smooths at a covalent boundary and 1.5 is worse without one, so
neither is a safe global default; these tests pin which one is selected when.
"""

import os

import pytest

pytest.importorskip("openmm", reason="the QM/MM driver imports OpenMM at module level")

from oqp.library.qmmm_driver import OpenQpQMMM


class _Selector:
    """Just enough of the driver to exercise the selection in isolation."""

    _select_boundary_switching = OpenQpQMMM._select_boundary_switching
    SWSCALE_WHOLE_MOLECULE = OpenQpQMMM.SWSCALE_WHOLE_MOLECULE
    SWSCALE_COVALENT_BOUNDARY = OpenQpQMMM.SWSCALE_COVALENT_BOUNDARY

    def __init__(self, link_atoms):
        self.link_atoms = link_atoms


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    monkeypatch.delenv("ESPF_SWSCALE", raising=False)


def test_covalent_boundary_selects_the_tighter_switching():
    selector = _Selector(link_atoms=[object()])
    assert selector._select_boundary_switching() == pytest.approx(1.5)
    assert float(os.environ["ESPF_SWSCALE"]) == pytest.approx(1.5)


def test_whole_molecule_qm_region_keeps_the_shipped_default():
    selector = _Selector(link_atoms=[])
    assert selector._select_boundary_switching() == pytest.approx(1.8)
    assert float(os.environ["ESPF_SWSCALE"]) == pytest.approx(1.8)


def test_an_explicit_setting_is_never_overridden(monkeypatch):
    """Someone who set it is running a sweep or reproducing a number."""
    monkeypatch.setenv("ESPF_SWSCALE", "2.1")
    selector = _Selector(link_atoms=[object()])
    assert selector._select_boundary_switching() is None
    assert selector.espf_swscale is None
    assert os.environ["ESPF_SWSCALE"] == "2.1"


def test_a_blank_setting_does_not_count_as_explicit(monkeypatch):
    monkeypatch.setenv("ESPF_SWSCALE", "   ")
    selector = _Selector(link_atoms=[object()])
    assert selector._select_boundary_switching() == pytest.approx(1.5)


def test_the_two_widths_are_the_documented_ones():
    """docs/espf_qmmm_switching.md and the native default must agree with these."""
    assert OpenQpQMMM.SWSCALE_WHOLE_MOLECULE == pytest.approx(1.8)
    assert OpenQpQMMM.SWSCALE_COVALENT_BOUNDARY == pytest.approx(1.5)
