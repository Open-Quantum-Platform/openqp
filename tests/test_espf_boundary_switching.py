"""Boundary-aware ESPF_SWSCALE selection (issue #260).

The switching width lives in the native ESPF gradient, which is handed a grid
and a density and cannot tell whether the QM/MM partition cut a bond.  The
driver can -- it just built the link atoms -- so the choice is made there.

1.8 over-smooths at a covalent boundary and 1.5 is worse without one, so
neither is a safe global default; these tests pin which one is selected when.

The selection helper lives in ``oqp.library.qmmm_connectivity`` precisely so
these tests run without OpenMM, which is optional and absent in the usual
developer and CI environments.
"""

import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Loaded by path, exactly like test_qmmm_connectivity.py, so the selection is
# exercised without the compiled OpenQP backend or OpenMM.
_spec = importlib.util.spec_from_file_location(
    "espf_boundary_switching_under_test",
    ROOT / "pyoqp/oqp/library/qmmm_connectivity.py",
)
_connectivity = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _connectivity
_spec.loader.exec_module(_connectivity)

SWSCALE_AUTO_MARKER = _connectivity.SWSCALE_AUTO_MARKER
SWSCALE_COVALENT_BOUNDARY = _connectivity.SWSCALE_COVALENT_BOUNDARY
SWSCALE_WHOLE_MOLECULE = _connectivity.SWSCALE_WHOLE_MOLECULE
select_boundary_switching = _connectivity.select_boundary_switching


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    monkeypatch.delenv("ESPF_SWSCALE", raising=False)
    monkeypatch.delenv(SWSCALE_AUTO_MARKER, raising=False)


def test_covalent_boundary_selects_the_tighter_switching():
    assert select_boundary_switching(covalent=True) == pytest.approx(1.5)
    assert float(os.environ["ESPF_SWSCALE"]) == pytest.approx(1.5)


def test_whole_molecule_qm_region_keeps_the_shipped_default():
    assert select_boundary_switching(covalent=False) == pytest.approx(1.8)
    assert float(os.environ["ESPF_SWSCALE"]) == pytest.approx(1.8)


def test_an_explicit_setting_is_never_overridden(monkeypatch):
    """Someone who set it is running a sweep or reproducing a number."""
    monkeypatch.setenv("ESPF_SWSCALE", "2.1")
    assert select_boundary_switching(covalent=True) is None
    assert os.environ["ESPF_SWSCALE"] == "2.1"


def test_a_blank_setting_does_not_count_as_explicit(monkeypatch):
    monkeypatch.setenv("ESPF_SWSCALE", "   ")
    assert select_boundary_switching(covalent=True) == pytest.approx(1.5)


def test_a_prior_automatic_choice_is_not_mistaken_for_user_configuration():
    """Two systems in one process must each get their own selection.

    The first selection writes ESPF_SWSCALE into the process environment for
    the native code to read; a second driver in the same process must treat
    that as its predecessor's automatic choice, not as a user override, in
    both directions.
    """
    assert select_boundary_switching(covalent=False) == pytest.approx(1.8)
    assert select_boundary_switching(covalent=True) == pytest.approx(1.5)
    assert float(os.environ["ESPF_SWSCALE"]) == pytest.approx(1.5)
    assert select_boundary_switching(covalent=False) == pytest.approx(1.8)
    assert float(os.environ["ESPF_SWSCALE"]) == pytest.approx(1.8)


def test_a_user_override_installed_between_runs_still_wins(monkeypatch):
    assert select_boundary_switching(covalent=True) == pytest.approx(1.5)
    monkeypatch.setenv("ESPF_SWSCALE", "2.1")
    assert select_boundary_switching(covalent=False) is None
    assert os.environ["ESPF_SWSCALE"] == "2.1"


def test_the_two_widths_are_the_documented_ones():
    """docs/espf_qmmm_switching.md and the native default must agree with these."""
    assert SWSCALE_WHOLE_MOLECULE == pytest.approx(1.8)
    assert SWSCALE_COVALENT_BOUNDARY == pytest.approx(1.5)


def test_the_driver_reuses_the_shared_selection():
    """The OpenMM-dependent driver must delegate, not fork, this policy."""
    pytest.importorskip(
        "openmm", reason="the QM/MM driver imports OpenMM at module level")
    if str(ROOT / "pyoqp") not in sys.path:
        sys.path.insert(0, str(ROOT / "pyoqp"))
    try:
        from oqp.library.qmmm_driver import OpenQpQMMM
    except (ImportError, RuntimeError) as exc:
        pytest.skip(f"oqp package is not importable here: {exc}")

    assert OpenQpQMMM.SWSCALE_WHOLE_MOLECULE == pytest.approx(
        SWSCALE_WHOLE_MOLECULE)
    assert OpenQpQMMM.SWSCALE_COVALENT_BOUNDARY == pytest.approx(
        SWSCALE_COVALENT_BOUNDARY)
