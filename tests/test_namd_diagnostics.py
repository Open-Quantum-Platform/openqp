"""Small contract tests for user-facing NAMD diagnostic tools."""

import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_hop_tracer_requires_an_explicit_input_and_derives_job_name():
    source = (ROOT / "tools" / "diagnostics" / "trace_namd_hop.py").read_text()

    assert 'parser.add_argument("--input", type=Path, required=True)' in source
    assert 'default="thymine.inp"' not in source
    assert "project = input_path.stem" in source


def test_hop_tracer_is_state_count_agnostic_and_rejects_untraced_modes():
    source = (ROOT / "tools" / "diagnostics" / "trace_namd_hop.py").read_text()

    assert 'self.coef[2]' not in source
    assert 'cmhp[2,' not in source
    assert "for state in range(nstate)" in source
    assert "for left in range(nstate)" in source
    assert "get('runtype', '')" in source
    assert "hop tracing requires a NAMD input" in source
    assert '("NAMD_SOC", "_log_soc")' in source
    assert "original_failure = (error, error.__traceback__)" in source
    assert "raise error.with_traceback(traceback)" in source


def test_hop_tracer_records_the_row_before_preserving_an_nve_abort(tmp_path):
    import sys

    sys.path.insert(0, str(ROOT / "pyoqp"))
    try:
        from oqp.library import namd as namd_module
    finally:
        sys.path.pop(0)

    path = ROOT / "tools" / "diagnostics" / "trace_namd_hop.py"
    spec = importlib.util.spec_from_file_location("trace_namd_hop_test", path)
    assert spec is not None and spec.loader is not None
    tracer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tracer)

    original_log = namd_module.NAMD._log_step
    original_prepare = namd_module.NAMD._prepare_hop_step

    def failing_log(_self, _istep, _r, **_kwargs):
        raise RuntimeError("NVE gate abort")

    trace = tmp_path / "failure.csv"
    try:
        namd_module.NAMD._log_step = failing_log
        tracer.install_trace(trace, skip_first_hop=False)
        driver = namd_module.NAMD.__new__(namd_module.NAMD)
        driver.nstate = 2
        driver.dt_fs = 0.5
        driver._last_hop_random = 0.25
        driver.active = 1
        driver.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
        driver.mol = SimpleNamespace(data={
            "OQP::namd_params": np.array([0.0, 0.0, 0.0, 0.25]),
            "OQP::namd_results": np.array([0.0, 0.1, 0.2, 0.0]),
            "OQP::td_states_overlap": np.eye(2),
        })
        driver._compute_tdc = lambda overlap: np.asarray(overlap)

        with pytest.raises(RuntimeError, match="NVE gate abort"):
            driver._log_step(1, np.zeros((1, 3)), hopped=True)
    finally:
        namd_module.NAMD._log_step = original_log
        namd_module.NAMD._prepare_hop_step = original_prepare

    rows = trace.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 2
    assert rows[0].startswith("step,t_fs,kernel_called")
    assert rows[1].startswith("1,0.5,1,1,1,0.25")


def test_hop_tracer_protects_input_replay_and_log_paths(tmp_path):
    path = ROOT / "tools" / "diagnostics" / "trace_namd_hop.py"
    spec = importlib.util.spec_from_file_location("trace_namd_hop_paths", path)
    assert spec is not None and spec.loader is not None
    tracer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tracer)

    input_path = tmp_path / "job.inp"
    random_file = tmp_path / "random.txt"
    log_path = tmp_path / "job.log"
    for protected in (input_path, random_file, log_path):
        protected.write_text("protected", encoding="utf-8")

    with pytest.raises(ValueError, match="input deck"):
        tracer.validate_trace_output(
            input_path, input_path=input_path, random_file=random_file,
            log_path=log_path)

    random_alias = tmp_path / "random-alias.csv"
    random_alias.symlink_to(random_file)
    with pytest.raises(ValueError, match="random replay file"):
        tracer.validate_trace_output(
            random_alias, input_path=input_path, random_file=random_file,
            log_path=log_path)

    log_alias = tmp_path / "log-alias.csv"
    os.link(log_path, log_alias)
    with pytest.raises(ValueError, match="job log"):
        tracer.validate_trace_output(
            log_alias, input_path=input_path, random_file=random_file,
            log_path=log_path)
    assert all(
        protected.read_text(encoding="utf-8") == "protected"
        for protected in (input_path, random_file, log_path)
    )
