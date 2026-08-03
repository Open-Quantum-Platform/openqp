"""Small contract tests for user-facing NAMD diagnostic tools."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_hop_tracer_requires_an_explicit_input_and_derives_job_name():
    source = (ROOT / "tools" / "diagnostics" / "trace_namd_hop.py").read_text()

    assert 'parser.add_argument("--input", type=Path, required=True)' in source
    assert 'default="thymine.inp"' not in source
    assert "project = args.input.stem" in source


def test_hop_tracer_is_state_count_agnostic_and_rejects_untraced_modes():
    source = (ROOT / "tools" / "diagnostics" / "trace_namd_hop.py").read_text()

    assert 'self.coef[2]' not in source
    assert 'cmhp[2,' not in source
    assert "for index, coefficient in enumerate(self.coef):" in source
    assert "for left in range(nstate):" in source
    assert 'get("qmmm_flag", False)' in source
    assert 'get("soc", False)' in source
    assert "QM/MM is rejected" in source
    assert "SOC is rejected" in source
