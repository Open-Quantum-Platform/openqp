"""Integral buffers must allow first use, resizing and repeated cleanup."""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_pair_buffers_with_runtime_bounds_checks(tmp_path):
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("GNU Fortran required for runtime allocation/bounds checks")
    # Compile the complete production modules. Only their basis/constant
    # dependencies are reduced to the fields needed by these storage routines.
    fixture = tmp_path / "fixture.f90"
    fixture.write_text("""
module precision
 integer,parameter :: dp=kind(1d0)
end module
module constants
 integer,parameter :: NUM_CART_BF(4)=[1,3,6,10]
end module
module basis_tools
 use precision
 type atoms_t
  real(dp)::xyz(3,8)=0
 end type
 type basis_set
  integer::nshell=1,mxcontr=1
  integer::g_offset(8)=[1,2,3,4,5,6,7,8],ncontr(8)=1,am(8)=1
  integer::origin(8)=1,ao_offset(8)=1,naos(8)=1,harmonic(8)=0
  real(dp)::shell_centers(8,3)=0,cc(8)=1,ex(8)=1
  type(atoms_t)::atoms
 end type
end module
""")
    program = tmp_path / "check.f90"
    program.write_text("""
program check
 use basis_tools
 use mod_shell_tools
 use int2_pairs
 implicit none
 type(basis_set)::b,b2
 type(shpair_t)::sp,sp2
 type(int2_pair_storage)::pairs
 type(int2_cutoffs_t)::cut
 integer::n,k
 integer,parameter::sizes(5)=[1,4,2,8,1]
 cut%quartet_cutoff=1d-12
 cut%exponent_cutoff=30d0
 do k=1,size(sizes)
  n=sizes(k)
  b%nshell=n
  b%mxcontr=n
  b2%mxcontr=2
  call sp%alloc(b)
  call sp2%alloc2(b,b2)
  if (size(sp%p)<n*n .or. size(sp2%p)<2*n) error stop 'shell capacity'
  call pairs%alloc(b,cut)
  if (size(pairs%ppid,2)<n*(n+1)/2) error stop 'pair capacity'
  if (size(pairs%g)/=n*(n+1)/2) error stop 'primitive capacity'
 end do
 call pairs%clean()
 call pairs%clean()
 call pairs%alloc(b,cut)
 call pairs%clean()
end program
""")
    exe = tmp_path / "check"
    subprocess.run([compiler, "-O0", "-g", "-fcheck=all", "-fdefault-integer-8",
                    "-ffree-line-length-none", str(fixture),
                    str(ROOT / "source/integrals/mod_shell_tools.F90"),
                    str(ROOT / "source/integrals/int2_pairs.F90"),
                    str(program), "-o", str(exe)], cwd=tmp_path,
                   capture_output=True, text=True, check=True)
    subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True,
                   check=True, timeout=10)


def test_omitted_cam_option_does_not_reuse_the_previous_call(tmp_path):
    import re

    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("GNU Fortran required")
    source = (ROOT / "source/integrals/int2.F90").read_text()
    routine = re.search(r"  subroutine int2_run\(.*?  end subroutine int2_run",
                        source, re.S).group()
    program = tmp_path / "cam.f90"
    program.write_text("""
module fixture
 implicit none
 integer,parameter::dp=kind(1d0),ERR_CAM_PARAM=7,WITH_ABORT=1
 type int2_compute_data_t
 end type
 type int2_compute_t
  integer::generic_calls=0,cam_calls=0
 contains
  procedure::run=>int2_run
  procedure::run_cam
  procedure::run_generic
 end type
contains
 subroutine run_cam(this,data,a,b,m,ac,bc)
  class(int2_compute_t)::this
  class(int2_compute_data_t)::data
  real(dp)::a,b,m
  real(dp),optional::ac,bc
  this%cam_calls=this%cam_calls+1
 end subroutine
 subroutine run_generic(this,data)
  class(int2_compute_t)::this
  class(int2_compute_data_t)::data
  this%generic_calls=this%generic_calls+1
 end subroutine
 subroutine show_message(text,mode)
  character(*)::text
  integer::mode
  error stop 'unexpected missing CAM parameters'
 end subroutine
""" + routine + """
end module
program check
 use fixture
 implicit none
 type(int2_compute_t)::driver
 type(int2_compute_data_t)::data
 integer::status
 call driver%run(data,cam=.true.,alpha=0.2d0,beta=0.3d0,mu=0.4d0)
 call driver%run(data,stat=status)
 if (status/=0) error stop 'omitted CAM retained previous setting'
 call driver%run(data,cam=.false.)
 if (driver%cam_calls/=1 .or. driver%generic_calls/=2) error stop 'wrong integral mode'
end program
""")
    exe = tmp_path / "cam"
    subprocess.run([compiler, "-fcheck=all", "-ffree-line-length-none",
                    str(program), "-o", str(exe)], cwd=tmp_path, check=True,
                   capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True, timeout=10)
