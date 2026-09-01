"""Analytic spin-resolved semilocal XC Fock-derivative verification."""

from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
POINT = ROOT / "source/dftlib/dft_gridint_mrsf_xc_fock_deriv_point.F90"
DRIVER = ROOT / "source/modules/mrsf_xc_fock_total_derivative.F90"
SELFTEST = ROOT / "tests/fortran/test_mrsf_xc_fock_total_derivative.F90"


def test_exact_pointwise_lda_gga_oracles_and_translation():
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        return
    with tempfile.TemporaryDirectory(prefix="oqp-mrsf-xc-fock-") as tmp:
        executable = Path(tmp) / "test_mrsf_xc_fock_total_derivative"
        subprocess.run(
            [
                compiler,
                str(ROOT / "source/precision.F90"),
                str(POINT),
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


def test_production_driver_has_complete_analytic_decomposition():
    source = DRIVER.read_text()
    assert "partition_weight_nuclear_derivatives" in source
    assert "moving_ao_pair_derivative" in source
    assert "moving_density_derivative" in source
    assert "utddft_fxc" in source
    assert "derivative_a=derivative_a+kernel_a" in source
    assert "derivative_b=derivative_b+kernel_b" in source
    assert "nMtx=ncart" in source
    assert "include_weight_derivative" not in source


def test_meta_gga_and_cam_fail_before_grid_or_kernel_evaluation():
    source = DRIVER.read_text()
    meta_gate = source.index("infos%functional%needTau")
    cam_gate = source.index("infos%dft%cam_flag")
    grid_call = source.index("call run_xc")
    kernel_call = source.index("call utddft_fxc")
    assert meta_gate < grid_call < kernel_call
    assert cam_gate < grid_call
    assert "status=-2" in source[meta_gate:grid_call]
    assert "status=-3" in source[cam_gate:grid_call]


def test_density_and_full_orbital_response_interfaces_are_exposed():
    source = DRIVER.read_text()
    assert "public :: mrsf_xc_fock_total_derivative" in source
    assert "public :: mrsf_xc_fock_total_derivative_from_dmo" in source
    dmo_body = source.split(
        "subroutine mrsf_xc_fock_total_derivative_from_dmo", 1
    )[1].split("end subroutine mrsf_xc_fock_total_derivative_from_dmo", 1)[0]
    assert "occupation_a" in dmo_body
    assert "occupation_b" in dmo_body
    assert "dmo_a(i,orbital,k)*mo_a(j,orbital)" in dmo_body
    assert "mo_a(i,orbital)*dmo_a(j,orbital,k)" in dmo_body


def test_no_nuclear_finite_difference_or_state_representation_enters_driver():
    source = DRIVER.read_text().lower()
    forbidden = ("slater", "determinant", "geometry_plus", "geometry_minus")
    # The documentation explicitly states that determinant/Slater objects are
    # absent, so inspect executable text after removing comments.
    executable = "\n".join(
        line.split("!", 1)[0] for line in source.splitlines()
    )
    assert not any(token in executable for token in forbidden)
