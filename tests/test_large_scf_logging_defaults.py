from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_default_guess_stays_huckel_for_legacy_large_hf_inputs():
    """Large legacy HF benchmarks depend on the old Hückel default unless overridden."""
    source = (ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py").read_text()

    assert "'guess': {" in source
    assert "'type': {'type': string, 'default': 'huckel'}" in source
    assert "'type': {'type': string, 'default': 'sap'}" not in source


def test_scf_iteration_log_uses_non_overflowing_scientific_fields():
    """Large anions can exceed fixed F-format widths during early SCF iterations."""
    source = (ROOT / "source" / "scf.F90").read_text()

    assert "es23.12" in source
    assert "es14.6" in source
    assert "i16" in source
    assert "f17.10" not in source
    assert "i13" not in source
