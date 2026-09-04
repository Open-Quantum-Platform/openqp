"""Analytic MRSF Z-vector derivative primitives and matrix-free solution."""

from pathlib import Path
import platform
import re
import shutil
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_z_response.F90"
SF_SOURCE = ROOT / "source/tdhf_sf_lib.F90"
RESPONSE_SOURCE = ROOT / "source/modules/tdhf_hessian_response.F90"
FORTRAN_ORACLE = ROOT / "tests/fortran/test_mrsf_hessian_z_response.F90"


def _extract_subroutine(source: str, name: str) -> str:
    match = re.search(
        rf"^\s*subroutine\s+{name}\b.*?^\s*end\s+subroutine\s+{name}\s*$",
        source,
        re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )
    assert match is not None, f"could not extract {name} from tdhf_sf_lib.F90"
    return match.group(0)


def test_variablewise_polarization_is_the_exact_bilinear_derivative():
    rng = np.random.default_rng(149104101)
    fock = rng.normal(size=(5, 5))
    density = rng.normal(size=(5, 3))
    dfock = rng.normal(size=(5, 5))
    ddensity = rng.normal(size=(5, 3))

    expected = dfock @ density + fock @ ddensity
    separated = 0.5 * (
        (fock + dfock) @ density
        - (fock - dfock) @ density
        + fock @ (density + ddensity)
        - fock @ (density - ddensity)
    )
    np.testing.assert_allclose(separated, expected, atol=3.0e-14, rtol=0.0)

    combined_forward_change = (
        (fock + dfock) @ (density + ddensity) - fock @ density
    )
    assert np.max(np.abs(combined_forward_change - expected)) > 1.0e-3


def test_source_reuses_spin_adapted_rhs_lhs_and_matrix_free_solver():
    compact = "".join(SOURCE.read_text().lower().split())
    assert "usetdhf_sf_lib,only:sfrorhs,sfrolhs" in compact
    assert (
        "usetdhf_hessian_response_mod,only:"
        "solve_mrsf_z_response_matrix_free"
    ) in compact
    assert "callsfrorhs" in compact
    assert "callsfrolhs" in compact
    assert "callsolve_mrsf_z_response_matrix_free" in compact
    assert compact.count("callevaluate_mrsf_z_rhs") == 6
    assert compact.count("callevaluate_mrsf_z_operator_action") == 6


def test_source_has_separate_first_derivative_polarizations():
    source = SOURCE.read_text().lower()
    for phrase in (
        "direct response polarization",
        "fock polarization at fixed baseline difference densities",
        "difference-density polarization at fixed baseline fock matrices",
        "orbital-energy polarization at fixed z",
        "response-image polarization at fixed z",
    ):
        assert phrase in source
    assert "(df)(dt)" in source
    assert "never enter the first derivative" in source


def test_two_somo_multiplicity_and_method_boundaries_fail_closed():
    compact = "".join(SOURCE.read_text().lower().split())
    assert "nocca-noccb/=2" in compact
    assert "target_multiplicity/=1.and.target_multiplicity/=3" in compact
    assert "if(is_umrsf)then" in compact
    assert "mrsf_z_response_umrsf_unsupported" in compact
    assert "if(cam_flag)then" in compact
    assert "mrsf_z_response_cam_unsupported" in compact


def test_nakata_provenance_and_spin_adapted_topology_are_explicit():
    source = SOURCE.read_text()
    assert "Hiroya Nakata's TDHF/TDDFT Hessian" in source
    assert "CO(+), OV(+), CV(+), OO, CV(-), OV(-), CO(-)" in source


def test_new_fortran_sources_respect_132_columns():
    for path in (SOURCE, FORTRAN_ORACLE):
        violations = [
            (line_number, len(line))
            for line_number, line in enumerate(path.read_text().splitlines(), 1)
            if len(line) > 132
        ]
        assert violations == []


def test_fortran_directional_oracle_and_matrix_free_solution(tmp_path: Path):
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    assert compiler is not None, "GNU Fortran is required for the strict oracle"

    sf_source = SF_SOURCE.read_text()
    sf_stub = tmp_path / "tdhf_sf_lib_z_response_oracle.F90"
    sf_stub.write_text(
        "module tdhf_sf_lib\ncontains\n"
        + _extract_subroutine(sf_source, "sfrorhs")
        + "\n"
        + _extract_subroutine(sf_source, "sfrolhs")
        + "\nend module tdhf_sf_lib\n"
    )
    executable = tmp_path / "test_mrsf_hessian_z_response"
    command = [
        compiler,
        "-std=f2018",
        "-O0",
        "-Wall",
        "-Wextra",
        "-Werror=line-truncation",
        "-ffree-line-length-132",
        "-fcheck=all",
        "-fbacktrace",
        "-J",
        str(tmp_path),
        str(ROOT / "source/precision.F90"),
        str(sf_stub),
        str(RESPONSE_SOURCE),
        str(SOURCE),
        str(FORTRAN_ORACLE),
        "-o",
        str(executable),
    ]
    if platform.system() == "Darwin":
        command.extend(["-framework", "Accelerate"])
    else:
        command.append("-lblas")
    subprocess.run(command, cwd=tmp_path, check=True)
    subprocess.run([str(executable)], cwd=tmp_path, check=True)
