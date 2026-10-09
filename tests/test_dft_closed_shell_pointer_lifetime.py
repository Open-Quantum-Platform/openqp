"""Bounds/pointer regressions for closed-shell and fxc DFT scratch storage.

Each test compiles the production routine (extracted from the source tree)
into a small fixture with ``-fcheck=all`` and runs it, so a reintroduced
defect aborts with a Fortran runtime error:

* dft_gridint_energy ``resetOrbPointers``: fock_b must not be remapped onto
  the unallocated fb2 of a closed-shell run.
* dft_gridint_tdxc_grad ``resetGradPointers``: with do_fxc the ground-state
  and X+Y slabs of tmpV/tmpG1 must both stay addressable (1:nDeriv).
* dft.F90 ``dftder``: the closed-shell beta density passed to derexc_blk's
  explicit ``db`` dummy must be allocated.
"""
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
DFTLIB = ROOT / 'source/dftlib'


def _compiler():
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    if compiler is None:
        pytest.skip('GNU Fortran required for bounds checks')
    return compiler


def _routine(source, name):
    match = re.search(r'^ *subroutine ' + name + r'\(.*?^ *end subroutine[^\n]*',
                      source, re.S | re.M)
    assert match, f'{name} not found'
    return match.group()


def _run(tmp_path, program):
    path = tmp_path / 'check.f90'
    path.write_text(program)
    exe = tmp_path / 'check'
    subprocess.run([_compiler(), '-O0', '-g', '-fcheck=all', '-fdefault-integer-8',
                    '-ffree-line-length-none', str(path), '-o', str(exe)],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    result = subprocess.run([str(exe)], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def test_closed_shell_fock_b_not_remapped(tmp_path):
    source = (DFTLIB / 'dft_gridint_energy.F90').read_text()
    routines = '\n'.join(_routine(source, name)
                         for name in ('parallel_start', 'clean', 'resetOrbPointers'))
    program = '''
module fixture
 implicit none
 integer,parameter::fp=kind(1d0)
 type xc_engine_t
  integer::numAOs=4,numAOs_p=3,maxPts=5,numPts=4,numTmpVec=2
  logical::hasBeta=.false.
 end type
 type xc_consumer_ks_t
  real(fp),allocatable::fa2(:,:),fb2(:,:),focks_(:,:),tmp_(:,:)
 contains
  procedure::parallel_start
  procedure::clean
  procedure::resetOrbPointers
 end type
contains
''' + routines + '''
end module
program check
 use fixture
 implicit none
 type(xc_engine_t)::e
 type(xc_consumer_ks_t),target::c
 real(fp),pointer::focks(:,:,:),tmp(:,:,:),fa(:,:),fb(:,:)
 integer::t
 logical::beta
 do t=1,2
  beta=(t==2)
  e%hasBeta=beta
  call c%parallel_start(e,2)
  call c%resetOrbPointers(e,focks=focks,tmp=tmp,fock_a=fa,fock_b=fb,myThread=2)
  if(.not.associated(fa))error stop 'fock_a not associated'
  if(associated(fb).neqv.beta)error stop 'fock_b association does not follow hasBeta'
  fa=1
  if(beta)then
   fb=2
   if(any(c%fb2(1:e%numAOs**2,2)/=2))error stop 'fock_b not mapped onto fb2'
  end if
  if(size(focks,3)/=merge(2,1,beta))error stop 'focks spin extent'
 end do
 call c%clean()
end program
'''
    _run(tmp_path, program)


@pytest.mark.parametrize('do_fxc', [True, False])
@pytest.mark.parametrize('has_beta', [False, True])
def test_fxc_gradient_scratch_slabs(tmp_path, do_fxc, has_beta):
    source = (DFTLIB / 'dft_gridint_tdxc_grad.F90').read_text()
    routine = _routine(source, 'resetGradPointers')
    flags = (f'logical,parameter::DO_FXC=.{str(do_fxc).lower()}.,'
             f'HAS_BETA=.{str(has_beta).lower()}.')
    program = '''
module fixture
 implicit none
 integer,parameter::fp=kind(1d0)
 type xc_engine_t
  integer::numAOs=4,numAOs_p=3,maxPts=5,numPts=4
  logical::hasBeta=.false.
 end type
 type xc_consumer_tdg_t
  integer::nMtx=2
  logical::do_fxc=.true.
  real(fp),allocatable::tmpGrad_(:,:),tmpV_(:,:,:),tmpG1_(:,:,:)
 contains
  procedure::resetGradPointers
 end type
contains
''' + routine + '''
end module
program check
 use fixture
 implicit none
 ''' + flags + '''
 type(xc_engine_t)::e
 type(xc_consumer_tdg_t),target::c
 real(fp),pointer::g(:,:,:),v(:,:,:,:,:),g1(:,:,:,:,:,:)
 integer::nspin,nderiv,d,t
 e%hasBeta=HAS_BETA
 c%do_fxc=DO_FXC
 nspin=merge(2,1,HAS_BETA)
 nderiv=merge(2,1,DO_FXC)
 ! storage sized as in parallel_start: full AO/point extents, 2 threads
 allocate(c%tmpGrad_(e%numAOs*3*c%nMtx,2), &
          c%tmpV_(e%numAOs*e%maxPts*c%nMtx*nspin,nderiv,2), &
          c%tmpG1_(e%numAOs*e%maxPts*3*c%nMtx*nspin,nderiv,2),source=0.0_fp)
 call c%resetGradPointers(e,g,v,g1,2)
 if(lbound(v,5)/=1.or.ubound(v,5)/=nderiv)error stop 'tmpV derivative bounds'
 if(lbound(g1,6)/=1.or.ubound(g1,6)/=nderiv)error stop 'tmpG1 derivative bounds'
 ! every slab is writable and the slabs do not overlap
 do d=1,nderiv
  v(:,:,:,:,d)=d
  g1(:,:,:,:,:,d)=d
 end do
 do d=1,nderiv
  if(any(v(:,:,:,:,d)/=d))error stop 'tmpV slabs overlap'
  if(any(g1(:,:,:,:,:,d)/=d))error stop 'tmpG1 slabs overlap'
 end do
 ! the other thread's block is untouched
 if(any(c%tmpV_(:,:,1)/=0).or.any(c%tmpG1_(:,:,1)/=0))error stop 'thread overlap'
 t=0
end program
'''
    _run(tmp_path, program)


def test_closed_shell_beta_density_allocated(tmp_path):
    source = (DFTLIB / 'dft.F90').read_text()
    block = re.search(r'^    if \(\.not\.allocated\(tda\)\) then.*?'
                      r'^    if \(urohf\) then.*?^    end if\n',
                      source, re.S | re.M)
    assert block, 'dftder tda/tdb allocation block not found'
    program = '''
module fixture
 implicit none
 integer,parameter::dp=kind(1d0)
 logical,parameter::WITH_ABORT=.true.
contains
 subroutine show_message(msg,flag)
  character(*),intent(in)::msg
  logical,intent(in)::flag
  error stop msg
 end subroutine
 ! same dummy declaration as derexc_blk: db(nbf,*) read only for urohf
 subroutine derexc_stub(nbf,da,db,urohf)
  integer,intent(in)::nbf
  real(dp),intent(inout)::da(nbf,*),db(nbf,*)
  logical,intent(in)::urohf
  if(urohf)db(1,1)=db(1,1)+da(1,1)
 end subroutine
 subroutine dftder_alloc(urohf)
  logical,intent(in)::urohf
  real(kind=dp),allocatable::tda(:,:),tdb(:,:),dedft(:,:)
  integer::nbf,nat,iok
  nbf=3; nat=2
''' + block.group() + '''
  tda=1
  if(urohf)tdb=2
  call derexc_stub(nbf,tda,tdb,urohf)
 end subroutine
end module
program check
 use fixture
 implicit none
 call dftder_alloc(.false.)
 call dftder_alloc(.true.)
end program
'''
    _run(tmp_path, program)
