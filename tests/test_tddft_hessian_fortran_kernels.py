"""Compile-and-run checks for standalone TDDFT Hessian Fortran kernels.

The programs under ``tests/fortran`` intentionally remain outside ``liboqp``.
Each test builds its program in a temporary directory with GNU Fortran and runs
the resulting executable.  The AO test extracts ``compAOvggg`` verbatim from
``basis_tools.F90`` into a minimal test-only module, avoiding the unrelated
link dependencies of the full basis module while still testing production
code rather than a duplicate implementation.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FORTRAN_TESTS = ROOT / "tests" / "fortran"


@pytest.fixture(scope="module")
def gfortran() -> str:
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("GNU Fortran compiler is not available")
    return compiler


def _compile_and_run(compiler: str, tmp_path: Path, sources: list[Path]) -> str:
    executable = tmp_path / "selftest"
    command = [
        compiler,
        "-std=f2018",
        "-Wall",
        "-Wextra",
        "-Werror=implicit-interface",
        "-J",
        str(tmp_path),
        *map(str, sources),
        "-o",
        str(executable),
    ]
    subprocess.run(command, cwd=tmp_path, check=True)
    return subprocess.check_output([str(executable)], cwd=tmp_path, text=True)


def test_gga_density_nuclear_point_selftest(gfortran: str, tmp_path: Path) -> None:
    output = _compile_and_run(
        gfortran,
        tmp_path,
        [
            ROOT / "source" / "precision.F90",
            ROOT / "source" / "dftlib" / "dft_gga_nuclear_point.F90",
            FORTRAN_TESTS / "test_dft_gga_nuclear_point.F90",
        ],
    )
    assert "GGA point second-derivative max error" in output


def test_partition_weight_hessian_selftest(gfortran: str, tmp_path: Path) -> None:
    output = _compile_and_run(
        gfortran,
        tmp_path,
        [
            ROOT / "source" / "precision.F90",
            ROOT / "source" / "dftlib" / "dft_partfunc.F90",
            ROOT / "source" / "dftlib" / "dft_partition_hessian.F90",
            FORTRAN_TESTS / "test_dft_partition_hessian.F90",
        ],
    )
    assert "dft partition Hessian selftest passed" in output


def _aoval_test_module() -> str:
    source = (ROOT / "source" / "basis_tools.F90").read_text()
    match = re.search(
        r"(?ms)^  subroutine compAOvggg\b.*?^  end subroutine compAOvggg\s*$",
        source,
    )
    if match is None:
        raise AssertionError("compAOvggg was not found in basis_tools.F90")

    # Only the fields and setup methods exercised by the standalone program are
    # present here.  The numerical routine below is copied verbatim at test
    # collection time from the production module.
    preamble = r"""
module atomic_structure_m
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
  implicit none
  type, public :: atomic_structure
    real(c_double), allocatable :: zn(:), mass(:), grad(:,:), xyz(:,:)
  contains
    procedure :: init => atomic_structure_init
  end type atomic_structure
contains
  function atomic_structure_init(self, natoms) result(status)
    class(atomic_structure), intent(inout) :: self
    integer(c_int64_t), intent(in) :: natoms
    integer(c_int) :: status
    allocate(self%zn(natoms), self%mass(natoms), self%grad(3,natoms), &
             self%xyz(3,natoms), stat=status)
  end function atomic_structure_init
end module atomic_structure_m

module basis_tools
  use precision, only: fp
  use atomic_structure_m, only: atomic_structure
  implicit none
  integer, parameter :: BAS_MXANG=1
  integer, parameter :: NUM_CART_BF(0:BAS_MXANG)=[1,3]
  integer, parameter :: NUM_SPH_BF(0:BAS_MXANG)=[1,3]
  integer, parameter :: cart_x(3,0:BAS_MXANG)=reshape([0,0,0, 1,0,0],[3,2])
  integer, parameter :: cart_y(3,0:BAS_MXANG)=reshape([0,0,0, 0,1,0],[3,2])
  integer, parameter :: cart_z(3,0:BAS_MXANG)=reshape([0,0,0, 0,0,1],[3,2])
  logical, parameter :: HARMONIC_ACTIVE=.false.
  type :: basis_set
    real(fp), allocatable :: ex(:), cc(:), bfnrm(:)
    real(fp), allocatable :: prim_mx_dist2(:), shell_mx_dist2(:)
    integer, allocatable :: g_offset(:), origin(:), am(:), harmonic(:)
    integer, allocatable :: ncontr(:), ao_offset(:), naos(:)
    integer :: nshell=0, nprim=0, nbf=0, mxcontr=0, mxam=0
    type(atomic_structure), pointer :: atoms
  contains
    procedure :: reserve
    procedure :: set_screening
    procedure :: compAOvgg_test
    procedure :: compAOvggg
    generic :: aoval => compAOvgg_test, compAOvggg
  end type basis_set
contains
  subroutine reserve(basis, nshell, nprim, nbf)
    class(basis_set), intent(inout) :: basis
    integer, intent(in) :: nshell, nprim, nbf
    allocate(basis%ex(nprim), basis%cc(nprim), basis%bfnrm(nbf), &
      basis%prim_mx_dist2(nprim), basis%shell_mx_dist2(nshell), &
      basis%g_offset(nshell), basis%origin(nshell), basis%am(nshell), &
      basis%harmonic(nshell), basis%ncontr(nshell), &
      basis%ao_offset(nshell), basis%naos(nshell))
  end subroutine reserve

  subroutine set_screening(basis, threshold)
    class(basis_set), intent(inout) :: basis
    real(fp), intent(in) :: threshold
    basis%prim_mx_dist2 = huge(threshold)
    basis%shell_mx_dist2 = huge(threshold)
  end subroutine set_screening

  subroutine cart2sph_vec(cart, spherical, am)
    real(fp), intent(in) :: cart(:)
    real(fp), intent(out) :: spherical(:)
    integer, intent(in) :: am
    spherical = cart(:size(spherical)) + 0.0_fp*am
  end subroutine cart2sph_vec

  subroutine compAOvgg_test(basis, xyz, nnz, v, gx, gy, gz, &
      gxx, gyy, gzz, gxy, gyz, gxz)
    class(basis_set) :: basis
    real(fp), intent(in) :: xyz(3)
    integer, intent(out) :: nnz
    real(fp), contiguous, intent(out) :: v(:), gx(:), gy(:), gz(:)
    real(fp), contiguous, intent(out) :: gxx(:), gyy(:), gzz(:)
    real(fp), contiguous, intent(out) :: gxy(:), gyz(:), gxz(:)
    real(fp) :: g3(size(v),10)
    call basis%compAOvggg(xyz,nnz,v,gx,gy,gz,gxx,gyy,gzz,gxy,gyz,gxz, &
      g3(:,1),g3(:,2),g3(:,3),g3(:,4),g3(:,5),g3(:,6),g3(:,7), &
      g3(:,8),g3(:,9),g3(:,10))
  end subroutine compAOvgg_test
"""
    return preamble + "\n" + match.group(0) + "\nend module basis_tools\n"


def test_ao_third_derivative_selftest(gfortran: str, tmp_path: Path) -> None:
    module_source = tmp_path / "basis_tools_aoval_test.F90"
    module_source.write_text(_aoval_test_module())
    output = _compile_and_run(
        gfortran,
        tmp_path,
        [
            ROOT / "source" / "precision.F90",
            module_source,
            FORTRAN_TESTS / "test_aoval_third_derivative.F90",
        ],
    )
    assert "AO G3 max |analytic - FD(G2)|" in output
