"""Strict tests for spin-adapted MRSF Cartesian response rows."""

from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_rows.F90"
SELFTEST = ROOT / "tests/fortran/test_tdhf_mrsf_hessian_rows.F90"


def _executable_text(source: str) -> str:
    return "\n".join(line.split("!", 1)[0] for line in source.splitlines())


def _assert_required_production_seams(source: str) -> None:
    required = (
        "grd2_mrsf_compute_data_t",
        "call contraction%init()",
        "call contraction%build_cart(basis)",
        "petite=.false.",
        "dplus=reference_spin_density",
        "pplus=relaxed_spin_density+d_relaxed_spin_density",
        "splus=seven_density+d_seven_density",
        "dplus=d_reference_spin_density",
        "dminus=-d_reference_spin_density",
        "call check_mrsf_seven_density_alias",
        "infos%tddft%mult/=1 .and. infos%tddft%mult/=3",
        "MRSF_ROWS_XC_INCOMPLETE",
    )
    for token in required:
        assert token in source


def test_one_e_oracle_raw_rows_and_fail_closed_boundaries():
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        return
    with tempfile.TemporaryDirectory(prefix="oqp-mrsf-hessian-rows-") as tmp:
        executable = Path(tmp) / "test_tdhf_mrsf_hessian_rows"
        subprocess.run(
            [
                compiler,
                "-cpp",
                "-DOQP_MRSF_ROWS_PRIMITIVES_ONLY",
                str(ROOT / "source/precision.F90"),
                str(SOURCE),
                str(SELFTEST),
                "-std=f2018",
                "-Wall",
                "-Wextra",
                "-Werror",
                "-Wimplicit-interface",
                "-fcheck=all",
                "-fbacktrace",
                "-o",
                str(executable),
            ],
            cwd=tmp,
            check=True,
        )
        subprocess.run([str(executable)], cwd=tmp, check=True)


def test_production_uses_seven_density_gradient_contraction_and_raw_rows():
    source = SOURCE.read_text()
    _assert_required_production_seams(source)
    executable = _executable_text(source).lower()
    assert "grad_nn" not in executable
    assert "petite=.true." not in executable
    assert "rows=0.5" not in executable
    assert "transpose(rows)" not in executable


def test_reference_response_excludes_fixed_ground_state_hessian_term():
    source = SOURCE.read_text()
    mixed = source.split(
        "dplus=d_reference_spin_density(:,:,:,response)", 1
    )[1].split("two_reference_mixed(:,response)", 1)[0]
    assert "dminus=-d_reference_spin_density(:,:,:,response)" in mixed
    assert "pplus=relaxed_spin_density" in mixed
    assert "pminus=relaxed_spin_density" in mixed
    assert "dplus=reference_spin_density" not in mixed
    assert "dminus=reference_spin_density" not in mixed
    assert "hess1" not in source.lower()
    assert "grd2_hess" not in source.lower()


def test_xc_and_channel_7_boundaries_fail_closed_before_contractions():
    source = SOURCE.read_text()
    alias_gate = source.index("call check_mrsf_seven_density_alias(")
    first_one_e = source.index("call grad_ee_overlap(")
    first_two_e = source.index("call mrsf_two_e_gradient_row(")
    assert alias_gate < first_one_e < first_two_e
    assert source.index(".not.xc_complete") < first_one_e
    assert "infos%dft%cam_flag .or." in source
    assert "infos%functional%needTau .or. infos%functional%needLapl" in source


def test_negative_source_mutations_are_detected():
    source = SOURCE.read_text()
    _assert_required_production_seams(source)
    mutations = (
        source.replace("petite=.false.", "petite=.true.", 1),
        source.replace("dminus=-d_reference_spin_density", "dminus=0.0_dp", 1),
        source.replace(
            "call check_mrsf_seven_density_alias",
            "call skip_alias_check",
            1,
        ),
        source.replace("splus=seven_density+d_seven_density", "splus=seven_density", 1),
    )
    for mutated in mutations:
        try:
            _assert_required_production_seams(mutated)
        except AssertionError:
            continue
        raise AssertionError("a required row-assembly mutation escaped detection")


def test_no_geometry_displacement_or_state_expansion_enters_executable_code():
    executable = _executable_text(SOURCE.read_text()).lower()
    forbidden = (
        "geometry_plus",
        "geometry_minus",
        "nuclear_step",
        "determinant",
        "slater",
        "fock_space",
    )
    assert not any(token in executable for token in forbidden)


def test_fortran_source_respects_132_column_limit():
    for path in (SOURCE, SELFTEST):
        long_lines = [
            (number, len(line))
            for number, line in enumerate(path.read_text().splitlines(), 1)
            if len(line) > 132
        ]
        assert not long_lines, f"{path.name}: {long_lines}"
