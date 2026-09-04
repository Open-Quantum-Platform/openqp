"""Exact local chain-rule checks for the MRSF semilocal XC Hessian."""

from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
POINT = ROOT / "source/dftlib/dft_gridint_mrsf_xc_hessian_point.F90"
SELFTEST = ROOT / "tests/fortran/test_mrsf_xc_hessian_point.F90"


def test_fixed_hessian_and_response_row_against_complex_step_oracles():
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        return
    with tempfile.TemporaryDirectory(prefix="oqp-mrsf-xc-hessian-") as tmp:
        executable = Path(tmp) / "test_mrsf_xc_hessian_point"
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
        subprocess.run([str(executable)],cwd=tmp,check=True)


def test_driver_uses_physical_spin_densities_and_third_xc_derivatives():
    driver = (ROOT / "source/dftlib/dft_gridint_mrsf_xc_hessian.F90").read_text()
    assert "xc_der3_contr" in driver
    assert "density_d" in driver
    assert "density_p" in driver
    executable = "\n".join(
        line.split("!",1)[0] for line in driver.lower().splitlines()
    )
    assert "determinant" not in executable
    assert "slater" not in executable


def test_meta_gga_and_cam_remain_fail_closed():
    driver = (ROOT / "source/dftlib/dft_gridint_mrsf_xc_hessian.F90").read_text()
    gate = driver.index("infos%functional%needTau")
    grid_call = driver.index("call run_xc")
    assert "infos%functional%needLapl" in driver[gate:grid_call]
    assert "infos%dft%cam_flag" in driver[gate:grid_call]
    assert "status=-2" in driver[gate:grid_call]


def test_grid_point_workspace_is_reused_without_pointwise_allocation():
    driver = (ROOT / "source/dftlib/dft_gridint_mrsf_xc_hessian.F90").read_text()
    point = driver[driver.index("  subroutine accumulate_point"):]
    point = point[:point.index("  end subroutine accumulate_point")]
    executable = "\n".join(
        line.split("!", 1)[0] for line in point.lower().splitlines()
    )
    assert "workspace" in executable
    assert "allocate(" not in executable
    assert "deallocate(" not in executable


def test_all_coordinate_response_densities_are_contracted_in_one_ao_pass():
    driver = (ROOT / "source/dftlib/dft_gridint_mrsf_xc_hessian.F90").read_text()
    point = driver[driver.index("  subroutine accumulate_point"):]
    point = point[:point.index("  end subroutine accumulate_point")]
    assert "response_field_value_gradient_batch" in point
    assert point.count("call gga_density_nuclear_point_first_batch") == 1
    assert "response_combined" in point
    assert "matrix=k+(field-1)*ncart" in point
