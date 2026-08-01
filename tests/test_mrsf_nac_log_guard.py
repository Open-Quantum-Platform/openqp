"""Contracts for avoiding process-global fort.6 during Python NAC calls."""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
INTERCHANGE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
ENERGY = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
GATE = ROOT / "tools" / "nac_lagrangian" / "fort6_guard_gate.py"


def _body(source: str, signature: str) -> str:
    end_name = signature.split("(", 1)[0]
    return source.split(f"subroutine {signature}", 1)[1].split(
        f"end subroutine {end_name}", 1
    )[0]


def _assert_log_guard(body: str, log_owner: str) -> None:
    assert "use io_constants, only: iw" in body
    assert "inquire(unit=iw, opened=log_was_open)" in body
    assert "if (.not. log_was_open)" in body
    assert f"file={log_owner}%log_filename" in body
    assert "position='append'" in body
    assert "if (.not. log_was_open) close(iw)" in body


def test_python_entry_points_preserve_the_fortran_log_state():
    _assert_log_guard(_body(ENERGY.read_text(), "mrsf_nac_response_C"), "inf")
    _assert_log_guard(_body(GRADIENT.read_text(), "mrsf_nac_esum_C"), "inf")
    _assert_log_guard(_body(INTERCHANGE.read_text(), "mrsf_nac_xc_adjoint_C"), "inf")
    _assert_log_guard(
        _body(INTERCHANGE.read_text(), "mrsf_nac_rohf_zvector_batch(infos"),
        "infos",
    )


def _run_gate(code: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(GATE), "--", sys.executable, "-c", code],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_fort6_gate_accepts_a_clean_command_in_a_fresh_workdir():
    result = _run_gate("from pathlib import Path; assert not Path('fort.6').exists()")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "PASS" in result.stdout


def test_fort6_gate_rejects_a_fallback_fortran_log():
    result = _run_gate("from pathlib import Path; Path('fort.6').write_text('bad')")
    assert result.returncode != 0
    assert "FAIL" in result.stderr
