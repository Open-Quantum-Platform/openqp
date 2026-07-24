"""Clean optional-OpenMM failure and test-runner classification."""

from types import SimpleNamespace

import pytest

from oqp.utils import oqp_tester, qmmm


def test_legacy_qmmm_entrypoint_raises_clean_missing_openmm_error(monkeypatch):
    missing = ModuleNotFoundError("No module named 'openmm'")
    monkeypatch.setattr(qmmm, "_HAVE_OPENMM", False)
    monkeypatch.setattr(qmmm, "_OPENMM_IMPORT_ERROR", missing)
    monkeypatch.setattr(qmmm, "mm", None)
    monkeypatch.setattr(qmmm, "app", None)
    monkeypatch.setattr(qmmm, "unit", None)

    with pytest.raises(ModuleNotFoundError, match="openmm.*QM/MM"):
        qmmm.openmm_init([1])


def test_oqp_tester_classifies_missing_openmm_as_skipped(monkeypatch, tmp_path):
    class MissingOpenMMRunner:
        def __init__(self, **kwargs):
            pass

        def run(self, test_mod=False):
            raise ModuleNotFoundError(
                "openmm is required for QM/MM calculations"
            )

    monkeypatch.setattr(oqp_tester, "Runner", MissingOpenMMRunner)
    tester = oqp_tester.OQPTester.__new__(oqp_tester.OQPTester)
    tester.output_dir = str(tmp_path)
    tester.mpi_manager = SimpleNamespace(use_mpi=0, rank=0)

    result = tester.run_single_test("qmmm_optional.inp")

    assert result["status"] == "SKIPPED"
    assert "optional OpenMM backend" in result["message"]
