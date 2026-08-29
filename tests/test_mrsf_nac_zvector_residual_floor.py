from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"


def test_rohf_nac_zvector_uses_documented_residual_floor():
    text = SOURCE.read_text()
    assert "tol=1.0e-20_dp" in text
    assert "residual(irhs) > 1.0e-16_dp" in text
    assert "residual norm <= 1e-8" in text
    assert "accepted after MINRES stagnation" in text
    assert "WITHOUT_ABORT" in text
