"""Focused checks for deterministic signs of diagonalized molecular orbitals."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
GUESS = ROOT / "source" / "guess.F90"


def _production_phase_module() -> str:
    source = GUESS.read_text()
    tolerance = re.search(
        r"(?m)^\s*real\(dp\), parameter :: ORBITAL_PHASE_TIE_RTOL\s*=.*$",
        source,
    )
    routine = re.search(
        r"(?ms)^\s*pure subroutine canonicalize_orbital_phases\b.*?"
        r"^\s*end subroutine canonicalize_orbital_phases\s*$",
        source,
    )
    assert tolerance is not None
    assert routine is not None
    return (
        "module orbital_phase_test_m\n"
        "  use precision, only: dp\n"
        "  implicit none\n"
        f"  {tolerance.group(0).strip()}\n"
        "  public :: canonicalize_orbital_phases\n"
        "contains\n"
        f"{routine.group(0)}\n"
        "end module orbital_phase_test_m\n"
    )


def test_orbital_phase_canonicalization_selftest(tmp_path: Path) -> None:
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("GNU Fortran compiler is not available")

    module_source = tmp_path / "orbital_phase_test_m.F90"
    module_source.write_text(_production_phase_module())
    executable = tmp_path / "orbital_phase_selftest"
    subprocess.run(
        [
            compiler,
            "-std=f2018",
            "-Wall",
            "-Wextra",
            "-Werror=implicit-interface",
            "-J",
            str(tmp_path),
            str(ROOT / "source" / "precision.F90"),
            str(module_source),
            str(ROOT / "tests" / "fortran" / "test_guess_orbital_phase.F90"),
            "-o",
            str(executable),
        ],
        cwd=tmp_path,
        check=True,
    )
    output = subprocess.check_output([str(executable)], cwd=tmp_path, text=True)
    assert "orbital phase canonicalization selftest passed" in output


def test_orbitals_are_canonicalized_after_ao_back_transform() -> None:
    source = GUESS.read_text()
    routine = re.search(
        r"(?ms)^\s*subroutine get_ab_initio_orbital\b.*?"
        r"^\s*end subroutine get_ab_initio_orbital\s*$",
        source,
    )
    assert routine is not None
    body = routine.group(0)
    assert body.index("call dgemm") < body.index("call canonicalize_orbital_phases")

