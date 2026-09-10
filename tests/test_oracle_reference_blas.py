"""The oracle tests' reference BLAS must agree with NumPy."""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_BLAS = ROOT / "tests/fortran/oracle_reference_blas.F90"

DRIVER = """\
program check_reference_blas
  implicit none
  integer, parameter :: m = 4, n = 3, k = 5
  double precision :: a(m, k), at(k, m), b(k, n), bt(n, k), c(m, n), c0(m, n)
  double precision :: x(n), xm(m), y0(m), yn0(n), y(m), yn(n)
  read (*, *) a, b, c0, x, xm, y0, yn0
  at = transpose(a)
  bt = transpose(b)
  c = c0
  call dgemm('N', 'N', m, n, k, 0.7d0, a, m, b, k, -0.3d0, c, m)
  write (*, '(*(es26.17e3))') c
  c = c0
  call dgemm('T', 'N', m, n, k, 0.7d0, at, k, b, k, -0.3d0, c, m)
  write (*, '(*(es26.17e3))') c
  c = c0
  call dgemm('N', 'T', m, n, k, 0.7d0, a, m, bt, n, 1.0d0, c, m)
  write (*, '(*(es26.17e3))') c
  c = c0
  call dgemm('T', 'T', m, n, k, 0.7d0, at, k, bt, n, 0.0d0, c, m)
  write (*, '(*(es26.17e3))') c
  y = y0
  call dgemv('N', m, n, 0.4d0, c0, m, x, 1, 1.5d0, y, 1)
  write (*, '(*(es26.17e3))') y
  yn = yn0
  call dgemv('T', m, n, 0.4d0, c0, m, xm, 1, -0.5d0, yn, 1)
  write (*, '(*(es26.17e3))') yn
  y = y0
  call dgemv('N', m, n, 0.4d0, c0, m, x, -1, 0.0d0, y, 1)
  write (*, '(*(es26.17e3))') y
end program check_reference_blas
"""


def test_reference_blas_matches_numpy(tmp_path):
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("GNU Fortran compiler is required")
    rng = np.random.default_rng(20260910)
    m, n, k = 4, 3, 5
    a, b, c0 = rng.normal(size=(m, k)), rng.normal(size=(k, n)), rng.normal(size=(m, n))
    x, xm, y0, yn0 = rng.normal(size=n), rng.normal(size=m), rng.normal(size=m), rng.normal(size=n)
    driver = tmp_path / "check_reference_blas.F90"
    driver.write_text(DRIVER)
    exe = tmp_path / "check_reference_blas"
    subprocess.run([compiler, "-std=f2018", "-O0", "-Wall", "-Wextra", "-fcheck=all",
                    str(REFERENCE_BLAS), str(driver), "-o", str(exe)], cwd=tmp_path, check=True)
    stdin = " ".join(f"{v:.17e}" for arr in (a, b, c0, x, xm, y0, yn0)
                     for v in np.asarray(arr).ravel(order="F"))
    lines = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True,
                           check=True).stdout.split("\n")
    got = [np.array(line.split(), dtype=float) for line in lines if line.strip()]
    ab = a @ b
    expected = [
        (0.7 * ab - 0.3 * c0).ravel(order="F"),
        (0.7 * ab - 0.3 * c0).ravel(order="F"),
        (0.7 * ab + 1.0 * c0).ravel(order="F"),
        (0.7 * ab).ravel(order="F"),
        0.4 * c0 @ x + 1.5 * y0,
        0.4 * c0.T @ xm - 0.5 * yn0,
        0.4 * c0 @ x[::-1],
    ]
    assert len(got) == len(expected)
    for value, reference in zip(got, expected):
        np.testing.assert_allclose(value, reference, rtol=1.0e-13, atol=1.0e-14)
