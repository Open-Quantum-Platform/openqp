"""SCF iteration table must print pathological values as numbers, not '*****'.

The table previously used bare f17.10/f14.8 edit descriptors; any value too
large for the column (e.g. the garbage energies produced by the dangling
ecp_zn buffer) printed as asterisks, hiding exactly the information needed to
debug.  The columns now route through fmt_real17/fmt_real14, which keep the
usual fixed-point form (and column width) while the value fits and fall back
to same-width scientific notation otherwise.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_iteration_table_routes_reals_through_width_preserving_formatters():
    source = (ROOT / "source" / "scf.F90").read_text()

    assert "function fmt_real17" in source
    assert "function fmt_real14" in source
    assert "fmt_real17(energy%etot)" in source
    assert "fmt_real14(diis_error)" in source


def test_formatters_fall_back_to_same_width_scientific_notation():
    source = (ROOT / "source" / "scf.F90").read_text()

    # fixed-point first choice, scientific fallback, identical column widths
    assert "(f17.10)" in source
    assert "(es17.8e3)" in source
    assert "(f14.8)" in source
    assert "(es14.5e3)" in source
