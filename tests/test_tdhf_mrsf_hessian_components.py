import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_mrsf_hessian_components_strict_fortran():
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        return
    with tempfile.TemporaryDirectory(prefix="oqp-mrsf-hess-components-") as tmp:
        executable = Path(tmp) / "test_mrsf_hessian_components"
        subprocess.run(
            [
                compiler,
                "-std=f2018",
                "-Wall",
                "-Wextra",
                "-Werror",
                "-fcheck=all",
                str(ROOT / "source/precision.F90"),
                str(ROOT / "source/modules/tdhf_mrsf_hessian_components.F90"),
                str(ROOT / "tests/fortran/test_tdhf_mrsf_hessian_components.F90"),
                "-o",
                str(executable),
            ],
            cwd=tmp,
            check=True,
        )
        subprocess.run([str(executable)], cwd=tmp, check=True)


def test_components_record_initial_verified_scope_and_provenance_boundary():
    source = (
        ROOT / "source/modules/tdhf_mrsf_hessian_components.F90"
    ).read_text()
    for condition in (
        "scf_type==3",
        "reference_mult==3",
        "nocca-noccb==2",
        "target_mult==1 .or. target_mult==3",
        ".not.umrsf",
        ".not.needs_tau",
        ".not.needs_laplacian",
        ".not.range_separated",
    ):
        assert condition in source
    assert "row_asymmetry=maxval" in source
