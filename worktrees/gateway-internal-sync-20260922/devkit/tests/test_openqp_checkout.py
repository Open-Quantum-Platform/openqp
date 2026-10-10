"""Exercise integrated and historical standalone engine discovery."""

from pathlib import Path

from openqp_checkout import find_openqp_root


MARKER = Path("source") / "modules" / "tdhf_mrsf_gradient.F90"


def test_finds_integrated_engine(monkeypatch):
    monkeypatch.delenv("OPENQP_SOURCE_ROOT", raising=False)

    assert find_openqp_root() == Path(__file__).resolve().parents[2]


def test_finds_sibling_engine_from_standalone_layout(tmp_path, monkeypatch):
    monkeypatch.delenv("OPENQP_SOURCE_ROOT", raising=False)
    helper = tmp_path / "openqp-devkit" / "tests" / "openqp_checkout.py"
    helper.parent.mkdir(parents=True)
    helper.touch()
    marker = tmp_path / "openqp" / MARKER
    marker.parent.mkdir(parents=True)
    marker.touch()

    assert find_openqp_root(helper) == tmp_path / "openqp"
