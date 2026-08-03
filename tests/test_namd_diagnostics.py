"""Small contract tests for user-facing NAMD diagnostic tools."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_hop_tracer_requires_an_explicit_input_and_derives_job_name():
    source = (ROOT / "tools" / "diagnostics" / "trace_namd_hop.py").read_text()

    assert 'parser.add_argument("--input", type=Path, required=True)' in source
    assert 'default="thymine.inp"' not in source
    assert "project = args.input.stem" in source
