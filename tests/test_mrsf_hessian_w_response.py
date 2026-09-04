import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_w_response.F90"
FORTRAN_TEST = ROOT / "tests/fortran/test_tdhf_mrsf_hessian_w_response.F90"


def test_mrsf_w_response_strict_fortran(tmp_path: Path):
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    assert compiler is not None, "a Fortran compiler is required for this test"
    executable = tmp_path / "test_mrsf_hessian_w_response"
    subprocess.run(
        [
            compiler,
            "-std=f2018",
            "-Wall",
            "-Wextra",
            "-Werror",
            "-fcheck=all",
            str(ROOT / "source/precision.F90"),
            str(SOURCE),
            str(FORTRAN_TEST),
            "-o",
            str(executable),
        ],
        cwd=tmp_path,
        check=True,
    )
    subprocess.run([str(executable)], cwd=tmp_path, check=True)


def test_mrsf_w_response_scope_provenance_and_line_length():
    text = SOURCE.read_text()
    lowered = text.lower()

    assert "Hiroya Nakata's TDHF/TDDFT analytical Hessian formulation" in text
    assert "noca - nocb /= 2" in text
    assert "reference_multiplicity /= 3" in text
    assert "target_multiplicity /= 1 .and. target_multiplicity /= 3" in text
    assert "MRSF_W_UMRSF_UNSUPPORTED" in text
    assert "MRSF_W_CAM_UNSUPPORTED" in text
    assert "0.25_dp*(transformed + transpose(transformed))" in text

    forbidden = (
        "type(determinant",
        "slater_",
        "fock_space",
        "displaced_geometry",
        "geometry_finite_difference",
    )
    for token in forbidden:
        assert token not in lowered

    for path in (SOURCE, FORTRAN_TEST):
        long_lines = [
            (line_number, len(line))
            for line_number, line in enumerate(path.read_text().splitlines(), 1)
            if len(line) > 132
        ]
        assert not long_lines, f"{path.name} exceeds 132 columns: {long_lines}"


def test_mrsf_w_response_product_rule_has_no_first_derivative_products():
    text = SOURCE.read_text()

    required_terms = (
        "dfa(k,x)*scr(k,open_first:open_last)",
        "fa(k,x)*dscr(k,open_first:open_last)",
        "dmo_energy(1:nocb)*scr(1:nocb,x)",
        "mo_energy(1:nocb)*dscr(1:nocb,x)",
        "matmul(matmul(dmo(:,:,perturbation), wmo), transpose(mo))",
        "matmul(matmul(mo, dwmo), transpose(mo))",
        "matmul(matmul(mo, wmo), transpose(dmo(:,:,perturbation)))",
    )
    for term in required_terms:
        assert term in text

    forbidden_products = (
        "dfa(k,x)*dscr",
        "dfb(open_first,open_first)*dscr",
        "dmo_energy(1:nocb)*dscr",
        "matmul(matmul(dmo(:,:,perturbation), dwmo)",
    )
    for term in forbidden_products:
        assert term not in text
