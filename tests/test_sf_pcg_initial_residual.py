"""Exercise the real spin-flip PCG initializer at exact and nonzero residuals."""
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_exact_initial_solution_does_not_normalize_zero(tmp_path):
    compiler = shutil.which("gfortran") or shutil.which("gfortran-15")
    if not compiler:
        pytest.skip("GNU Fortran required for floating-point trap regression")
    text = (ROOT / "source/tdhf_sf_lib.F90").read_text()
    routine = re.search(r"  subroutine pcgrbpini\(.*?end subroutine pcgrbpini", text, re.S).group()
    fixture = '''module precision
  integer, parameter :: dp = kind(1.0d0)
end module
module real_initializer
contains
''' + routine + '''
end module
program check
  use precision
  use real_initializer
  use, intrinsic :: ieee_arithmetic
  implicit none
  real(dp) :: r(3), pk(3), error, d(3), ax(3), precond(3)
  integer :: repeat
  precond = [1.0_dp, 0.5_dp, 2.0_dp]
  do repeat = 1, 3
    d = 0.0_dp
    ax = 0.0_dp
    call pcgrbpini(r, pk, error, d, precond, ax)
    if (error /= 0.0_dp .or. any(pk /= 0.0_dp)) stop 1
    d = [1.0_dp, 2.0_dp, -1.0_dp]
    ax = d
    call pcgrbpini(r, pk, error, d, precond, ax)
    if (error /= 0.0_dp .or. any(pk /= 0.0_dp)) stop 2
    ax = 0.0_dp
    call pcgrbpini(r, pk, error, d, precond, ax)
    if (error /= 6.0_dp .or. any(.not. ieee_is_finite(pk))) stop 3
    if (maxval(abs(pk-[0.2_dp,0.2_dp,-0.4_dp])) > 1e-14_dp) stop 4
  end do
end program
'''
    src = tmp_path / "initial.f90"
    src.write_text(fixture)
    executable = tmp_path / "initial"
    subprocess.run([compiler, "-O0", "-ffpe-trap=invalid,zero,overflow",
                    str(src), "-o", str(executable)], check=True)
    subprocess.run([str(executable)], check=True)


def test_mrsf_cg_skips_iteration_for_converged_or_invalid_initial_state():
    text = (ROOT / "source/modules/tdhf_mrsf_z_vector.F90").read_text()
    body = text.split("subroutine run_mrsf_cg_zvector()", 1)[1].split("end subroutine run_mrsf_cg_zvector", 1)[0]
    guard = "if (mrsf_zvector_breakdown .or. error < cnvtol) return"
    assert body.index(guard) > body.index("non-finite initial PCG state")
    assert body.index(guard) < body.index("do iter = 1, infos%control%maxit_zv")
