from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]
TOOL_DIR = ROOT / "tools" / "scf-converger-ml"


def test_distilled_scf_selector_matches_runtime_copy():
    training_model = TOOL_DIR / "model" / "scf_selector_model.py"
    runtime_model = ROOT / "pyoqp" / "oqp" / "library" / "scf_selector_model.py"

    assert training_model.read_text() == runtime_model.read_text()


def test_openqp_scf_selector_database_schema():
    db_path = TOOL_DIR / "data" / "db_openqp_postfix.csv"
    with db_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        fields = set(reader.fieldnames)

    required = {
        "system",
        "tier",
        "ref",
        "mult",
        "func",
        "converger",
        "natom",
        "charge",
        "converged",
        "niter",
    }
    assert required <= fields
