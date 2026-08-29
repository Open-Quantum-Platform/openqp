from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"


def test_rohf_nac_zvector_uses_scaled_machine_precision_floor():
    text = SOURCE.read_text()
    assert "tol=1.0e-20_dp" in text
    assert "10.0_dp*sqrt(epsilon(1.0_dp))" in text
    assert "max(1.0_dp, sum(rhs(:,irhs)*rhs(:,irhs)))" in text
    assert "scaled_residual_sq > fallback_rel_sq" in text
    assert "fallback relative norm=10*sqrt(epsilon)" in text
    assert "accepted after MINRES stagnation" in text
    assert "WITHOUT_ABORT" in text


def test_scaled_floor_accepts_roundoff_stagnation_but_rejects_large_residual():
    fallback_rel_sq = (10.0 * sys.float_info.epsilon**0.5) ** 2

    def accepted(residual_sq, rhs_norm_sq):
        return residual_sq / max(1.0, rhs_norm_sq) <= fallback_rel_sq

    assert accepted(1.36622964e-16, 1.0)
    assert accepted(1.0e-12, 100.0)
    assert not accepted(1.0e-12, 1.0)
