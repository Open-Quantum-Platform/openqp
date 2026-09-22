"""Run production DFT buffer routines across functional and size changes."""
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_xc_buffers_and_empty_functionals(tmp_path):
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    if compiler is None:
        pytest.skip('GNU Fortran required for bounds checks')
    source = (ROOT / 'source/dftlib/dft_gridint_fxc.F90').read_text()
    routines = []
    for name in ('parallel_start', 'clean'):
        routines.append(re.search(r'  subroutine ' + name + r'\(.*?  end subroutine[^\n]*',
                                  source, re.S).group())
    fsrc = (ROOT / 'source/dftlib/functionals.F90').read_text()
    routines.append(re.search(r'  logical function can_calculate.*?  end function can_calculate',
                              fsrc, re.S).group())
    program = '''
module fixture
 implicit none
 integer,parameter::fp=kind(1d0),OQP_FUNTYP_MGGA=2
 type xc_engine_t
  integer::numAOs=2,maxPts=3,numTmpVec=2,funTyp=1
  logical::hasBeta=.false.
 end type
 type xc_consumer_tde_t
  integer::nMtx=2
  real(fp),allocatable::focks(:,:,:,:,:),mo(:,:,:,:,:),rRho(:,:,:,:),drRho(:,:,:,:,:),rTau(:,:,:,:)
  real(fp),allocatable::focks_(:,:),tmpMO_(:,:),tmpDensity_(:,:,:),tmp_(:,:),moG1_(:,:)
 contains
  procedure::clean
  procedure::parallel_start
 end type
 type functional_t
  integer,allocatable::functionals_list(:)
 contains
  procedure::can_calculate
 end type
contains
''' + '\n'.join(routines) + '''
end module
program check
 use fixture
 implicit none
 type(xc_engine_t)::e
 type(xc_consumer_tde_t)::c
 type(functional_t)::f
 integer::n,k
 integer,parameter::sizes(4)=[2,4,1,2]
 if(f%can_calculate())error stop 'empty functional'
 allocate(f%functionals_list(0))
 if(f%can_calculate())error stop 'zero functionals'
 deallocate(f%functionals_list)
 allocate(f%functionals_list(1))
 if(.not.f%can_calculate())error stop 'valid functional'
 do k=1,size(sizes)
  n=sizes(k);e%numAOs=n
  e%funTyp=1
  call c%parallel_start(e,1)
  c%focks=7
  call c%parallel_start(e,1)
  if(any(c%focks/=0))error stop 'stale focks'
  e%funTyp=OQP_FUNTYP_MGGA
  call c%parallel_start(e,1)
  if(.not.allocated(c%moG1_))error stop 'missing MGGA storage'
  call c%parallel_start(e,2)
  if(size(c%moG1_,2)/=2)error stop 'thread resize'
  e%funTyp=1
  call c%parallel_start(e,1)
  if(allocated(c%moG1_))error stop 'obsolete MGGA storage'
 end do
 call c%clean()
 call c%clean()
end program
'''
    path = tmp_path / 'check.f90'
    path.write_text(program)
    exe = tmp_path / 'check'
    subprocess.run([compiler, '-O0', '-fcheck=all', '-fdefault-integer-8',
                    '-ffree-line-length-none', str(path), '-o', str(exe)],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True, timeout=10)
