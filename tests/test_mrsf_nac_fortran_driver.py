"""Static contracts for the single-call resident MRSF NAC driver."""

from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
METRIC = ROOT / "source" / "modules" / "mrsf_nac_metric_data.F90"
HEADER = ROOT / "include" / "oqp.h"
PYTHON = ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py"


def test_single_production_entry_is_exported_and_python_is_thin():
    driver = DRIVER.read_text()
    header = HEADER.read_text()
    production = PYTHON.read_text()
    assert 'bind(C, name="mrsf_nac_lagrangian")' in driver
    assert "void mrsf_nac_lagrangian(struct oqp_handle_t *inf);" in header
    assert production.count("oqp.mrsf_nac_lagrangian(mol)") == 1
    for forbidden in (
        "oqp.mrsf_nac_metric_data(mol)",
        "oqp.mrsf_nac_wpair(mol",
        "oqp.mrsf_nac_amp_pair(mol",
        "oqp.mrsf_nac_response(mol)",
        "oqp.mrsf_nac_rohf_zvector(mol)",
        "for I in range(nstate)",
        "for J in range(nstate)",
        "NAC_ANALYTIC_FORWARD_GATE",
        "oqp.hf_hessian(mol)",
    ):
        assert forbidden not in production


def test_driver_uses_actual_resident_state_count_and_streamed_metric():
    driver = DRIVER.read_text()
    metric = METRIC.read_text()
    assert "nstate64 = infos%tddft%nstate" in driver
    assert "size(energies) /= nstate" in driver
    assert "size(bvec_mo,2) /= nstate" in driver
    assert "call mrsf_nac_metric_column(infos, jstate, gamma_column)" in driver
    assert "call mrsf_nac_metric_data" not in driver
    assert "raw_sij(noca,noca,nstate)" in metric
    assert "raw_sab(nvirb,nvirb,nstate)" in metric
    assert "raw_sia(noca,nvirb,nstate)" in metric
    assert "energies_saved = energies" in driver
    assert "gap = energies_saved(jstate)-energies_saved(istate)" in driver
    assert "default_int_limit64/state_pair_size64" in driver


def test_driver_pair_sequence_contains_one_adjoint_per_ordered_pair():
    driver = DRIVER.read_text()
    body = driver.split(
        "subroutine mrsf_nac_lagrangian(infos)", 1
    )[1].split("end subroutine mrsf_nac_lagrangian", 1)[0]
    expected = (
        "call mrsf_nac_wpair_impl(infos, istate, jstate)",
        "call mrsf_nac_amp(infos, istate, jstate)",
        "call mrsf_nac_esum(infos, istate, jstate)",
        "call mrsf_nac_response(infos)",
        "call mrsf_nac_rohf_pair_overlap(infos)",
        "call mrsf_nac_rohf_zvector(infos)",
        "call mrsf_nac_rohf_hf_adjoint(infos)",
        "call mrsf_nac_xc_adjoint(infos)",
        "call mrsf_nac_pair_accumulate(infos, istate, jstate)",
    )
    positions = [body.index(token) for token in expected]
    assert positions == sorted(positions)
    assert body.count("call mrsf_nac_rohf_zvector(infos)") == 1
    pair_loop = body.index("do istate = 1, nstate")
    pair_skip = body.index("if (istate == jstate) cycle", pair_loop)
    zvector_call = body.index("call mrsf_nac_rohf_zvector(infos)")
    pair_end = body.index("end do", zvector_call)
    assert pair_loop < pair_skip < zvector_call < pair_end
    assert "call cphf_solve_rohf" not in body
    assert "hf_hessian" not in body


def test_python_production_call_cannot_enable_forward_cphf(monkeypatch):
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "nac_analytic_no_forward_cphf", PYTHON
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    class Data(dict):
        def __init__(self):
            super().__init__({
                "OQP::td_bvec_mo": np.arange(6.0).reshape(2, 3),
                "OQP::td_energies": np.array([0.1, 0.2, 0.3]),
                "natom": 1,
            })
            self._data = SimpleNamespace(
                control=SimpleNamespace(int2e_cutoff=1.0e-12),
                tddft=SimpleNamespace(nstate=3),
            )

    mol = SimpleNamespace(
        config={
            "scf": {"conv": 1.0e-10},
            "tdhf": {"conv": 1.0e-10, "multiplicity": 1},
        },
        data=Data(),
    )
    calls = []

    def resident_driver(active_mol):
        calls.append("mrsf_nac_lagrangian")
        active_mol.data["OQP::nac_dcv"] = np.zeros(27)
        active_mol.data["OQP::nac_nacv"] = np.zeros(27)

    def forbidden_forward_driver(_):
        raise AssertionError("production analytic NAC called forward CPHF")

    fake_oqp = SimpleNamespace(
        mrsf_nac_lagrangian=resident_driver,
        hf_hessian=forbidden_forward_driver,
    )
    monkeypatch.setitem(sys.modules, "oqp", fake_oqp)
    monkeypatch.setenv("NAC_ANALYTIC_FORWARD_GATE", "1")
    nacv, dcv = module.analytic_nac(mol)

    assert calls == ["mrsf_nac_lagrangian"]
    assert nacv.shape == (3, 3, 1, 3)
    assert dcv.shape == (3, 3, 1, 3)


def test_driver_guards_scope_gaps_and_restores_mutated_state():
    driver = DRIVER.read_text()
    for required in (
        "infos%control%scftype /= 3",
        "infos%tddft%umrsf",
        "infos%tddft%mult /= 1",
        "infos%control%conv > 1.0e-8_dp",
        "infos%tddft%cnvtol > 1.0e-8_dp",
        ".not. infos%mol_energy%SCF_converged",
        ".not. infos%mol_energy%Davidson_converged",
        "abs(gap) <= gap_floor",
        "bvec_mo = bvec_saved",
        "infos%control%int2e_cutoff = cutoff_saved",
    ):
        assert required in driver
    assert "inquire(unit=iw, opened=log_was_open)" in driver
