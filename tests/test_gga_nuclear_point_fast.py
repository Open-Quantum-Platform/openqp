"""The atom-block GGA density nuclear derivatives equal the reference loops.

gga_density_nuclear_point (used per grid point by the TDDFT fixed-density XC
Hessian) contracts the density into per-AO vectors and atom-block sums instead
of looping over every atom pair.  This test compiles mod_dft_gga_nuclear_point
with -fcheck=all and compares it with gga_density_nuclear_point_reference (the
original quadruple loop) on random AO data, for symmetric and nonsymmetric
densities, several atom/AO counts, and AOs on every atom or only some atoms.
"""
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]

PROGRAM = """
program check
  use precision, only: fp
  use mod_dft_gga_nuclear_point
  implicit none
  integer :: case_id, nao, nat, i
  real(fp), allocatable :: d(:,:), aov(:), g1(:,:), g2(:,:), g3(:,:)
  integer, allocatable :: atom(:)
  real(fp), allocatable :: r1(:,:), rg1(:,:,:), r2(:,:,:,:), rg2(:,:,:,:,:)
  real(fp), allocatable :: f1(:,:), fg1(:,:,:), f2(:,:,:,:), fg2(:,:,:,:,:)
  real(fp) :: err, scale
  integer :: seed(64)
  seed = 4711
  call random_seed(put=seed(1:size(seed)))
  do case_id = 1, 6
    nat = 1 + mod(case_id, 4) + case_id/3
    nao = 3*nat + case_id
    allocate(d(nao,nao), aov(nao), g1(nao,3), g2(nao,6), g3(nao,10), atom(nao))
    call random_number(d); d = d - 0.5_fp
    if (mod(case_id,2) == 0) d = 0.5_fp*(d + transpose(d))
    call random_number(aov); call random_number(g1); call random_number(g2); call random_number(g3)
    g1 = g1 - 0.5_fp; g2 = g2 - 0.5_fp; g3 = g3 - 0.5_fp
    do i = 1, nao
      atom(i) = 1 + mod(i*7, nat)
      if (case_id == 5) atom(i) = 1 + mod(i, max(1,nat-1))   ! last atom without AOs
    end do
    allocate(r1(3,nat), rg1(3,3,nat), r2(3,3,nat,nat), rg2(3,3,3,nat,nat))
    allocate(f1(3,nat), fg1(3,3,nat), f2(3,3,nat,nat), fg2(3,3,3,nat,nat))
    call gga_density_nuclear_point_reference(d, atom, aov, g1, g2, g3, r1, rg1, r2, rg2)
    call gga_density_nuclear_point(d, atom, aov, g1, g2, g3, f1, fg1, f2, fg2)
    scale = max(1.0_fp, maxval(abs(rg2)))
    err = max(maxval(abs(f1-r1)), maxval(abs(fg1-rg1)), maxval(abs(f2-r2)), maxval(abs(fg2-rg2)))
    if (err > 1.0e-12_fp*scale) then
      print *, 'case', case_id, 'nat', nat, 'nao', nao, 'max error', err
      error stop 'fast and reference nuclear derivatives differ'
    end if
    deallocate(d, aov, g1, g2, g3, atom, r1, rg1, r2, rg2, f1, fg1, f2, fg2)
  end do
end program
"""


def test_fast_matches_reference(tmp_path):
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("GNU Fortran required")
    (tmp_path / "precision.f90").write_text(
        "module precision\n  integer, parameter :: fp = kind(1d0)\nend module\n")
    module = tmp_path / "gga.F90"
    module.write_text((ROOT / "source/dftlib/dft_gga_nuclear_point.F90").read_text())
    main = tmp_path / "check.f90"
    main.write_text(PROGRAM)
    exe = tmp_path / "check"
    subprocess.run([compiler, "-O0", "-g", "-fcheck=all", "-ffree-line-length-none",
                    str(tmp_path / "precision.f90"), str(module), str(main), "-o", str(exe)],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    result = subprocess.run([str(exe)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
