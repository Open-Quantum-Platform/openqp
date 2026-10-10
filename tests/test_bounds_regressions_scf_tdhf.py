"""Bounds-checked regressions for the RHF SCF, molgrid, TDDFT and DFT-gradient fixes.

Each fixture compiles the production routine (extracted from the source
tree) with ``-fcheck=all`` and runs it, so a reintroduced defect aborts with a
Fortran runtime error, as in test_dft_closed_shell_pointer_lifetime.py:

* scf.F90 ``build_scf_density``: an RHF density array has one column and no
  beta MOs; the beta column must never be referenced (was pdmat(:,2)).
* dft_molgrid.F90 ``is_vector_inside_shape``: GCC 14/15 -fcheck=bounds
  aborted on the vecs(:, ubound(vecs, 2)) section argument.
* tdhf_lib.F90 response-image slots: a one-copy apb/amb image read from a
  higher thread id must clamp to slot 1, with the slot taken into a local.
* dft_gridint_grad.F90 ``derexc_blk``: a closed-shell run must not associate
  xc_opts%wfBeta with the unallocated beta scratch density.

The source checks next to each fixture fail if the old construct returns,
since GCC 12/13 (and release builds) would not abort on it.
"""
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'source'


def _compiler():
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    if compiler is None:
        pytest.skip('GNU Fortran required for bounds checks')
    return compiler


def _routine(source, kind, name):
    match = re.search(r'^ *(?:pure +)?' + kind + r' ' + name + r'\(.*?^ *end ' + kind + r'[^\n]*',
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


def test_rhf_scf_density_has_no_beta_column(tmp_path):
    source = (SOURCE / 'scf.F90').read_text()
    # every SCF density rebuild goes through the column-count-aware helper
    assert 'get_ab_initio_density(pdmat(:,1)' not in source
    assert 'get_ab_initio_density(dens_prev(:,1)' not in source
    routine = _routine(source, 'subroutine', 'build_scf_density')
    program = '''
module precision
 implicit none
 integer,parameter::dp=kind(1d0)
end module
module types
 implicit none
 type information
  integer::dummy=0
 end type
end module
module basis_tools
 implicit none
 type basis_set
  integer::dummy=0
 end type
end module
module guess
 use precision, only: dp
 use types, only: information
 use basis_tools, only: basis_set
 implicit none
 integer::nbeta=0
contains
 subroutine get_ab_initio_density(alpha_density,alpha_orbital,beta_density,beta_orbital,infos,basis)
  type(information),intent(in)::infos
  type(basis_set),intent(in)::basis
  real(dp)::alpha_density(:),alpha_orbital(:,:)
  real(dp),optional::beta_density(:),beta_orbital(:,:)
  alpha_density=sum(alpha_orbital)
  if(present(beta_density))then
   beta_density=sum(beta_orbital)
   nbeta=nbeta+1
  end if
 end subroutine
end module
module fixture
 implicit none
contains
''' + routine + '''
end module
program check
 use precision, only: dp
 use types, only: information
 use basis_tools, only: basis_set
 use guess, only: nbeta
 use fixture
 implicit none
 type(information)::infos
 type(basis_set)::basis
 real(dp),allocatable::rhf(:,:),uhf(:,:),mo_a(:,:),mo_b(:,:)
 allocate(rhf(6,1),uhf(6,2),source=0.0_dp)
 allocate(mo_a(3,3),source=1.0_dp)
 allocate(mo_b(3,3),source=2.0_dp)
 ! RHF: one density column, beta MOs absent
 call build_scf_density(rhf,mo_a,infos=infos,basis=basis)
 if(nbeta/=0.or.any(rhf/=9))error stop 'RHF density'
 ! UHF/ROHF: both columns built
 call build_scf_density(uhf,mo_a,mo_b,infos,basis)
 if(nbeta/=1.or.any(uhf(:,1)/=9).or.any(uhf(:,2)/=18))error stop 'UHF density'
end program
'''
    _run(tmp_path, program)


def test_molgrid_inside_shape_bounds(tmp_path):
    source = (SOURCE / 'dftlib' / 'dft_molgrid.F90').read_text()
    routine = _routine(source, 'function', 'is_vector_inside_shape')
    # (code lines only; the explanatory comment names the old section)
    assert re.search(r'^[^!\n]*vecs\(:, *ubound\(vecs', routine, re.M) is None
    cross = _routine(source, 'function', 'cross_product')
    program = '''
module fixture
 implicit none
 integer,parameter::fp=kind(1d0)
contains
''' + cross + '\n' + routine + '''
end module
program check
 use fixture
 implicit none
 real(fp)::tri(3,3),test(3)
 ! spherical triangle around +z
 tri(:,1)=[1.0_fp,0.0_fp,1.0_fp]
 tri(:,2)=[-0.5_fp,0.8_fp,1.0_fp]
 tri(:,3)=[-0.5_fp,-0.8_fp,1.0_fp]
 test=[0.0_fp,0.0_fp,1.0_fp]
 if(.not.is_vector_inside_shape(test,tri))error stop 'inside'
 test=[5.0_fp,5.0_fp,1.0_fp]
 if(is_vector_inside_shape(test,tri))error stop 'outside'
end program
'''
    _run(tmp_path, program)


def test_tdhf_response_image_slot_clamped(tmp_path):
    source = (SOURCE / 'tdhf_lib.F90').read_text()
    # the image slot is never computed inside the associate selector
    assert 'this%apb(:,:,:,min(' not in source
    assert 'this%amb(:,:,:,min(' not in source
    slots = re.findall(r'^ *iapb = min\(mythread, size\(this%apb, 4\)\)\n'
                       r' *iamb = min\(mythread, size\(this%amb, 4\)\)\n', source, re.M)
    assert len(slots) == 3, 'response image slots must be clamped into locals'
    program = '''
module fixture
 implicit none
 integer,parameter::dp=kind(1d0)
 type images_t
  real(dp),allocatable::apb(:,:,:,:),amb(:,:,:,:)
 end type
contains
 subroutine update(this,mythread)
  type(images_t),target,intent(inout)::this
  integer,intent(in)::mythread
  integer::iapb,iamb
''' + slots[0] + '''
  associate(apb=>this%apb(:,:,:,iapb),amb=>this%amb(:,:,:,iamb))
   apb(1,1,1)=apb(1,1,1)+mythread
   amb(1,1,1)=amb(1,1,1)-mythread
  end associate
 end subroutine
end module
program check
 use fixture
 implicit none
 type(images_t),target::img
 ! active image per thread, inactive image with a single shared copy
 allocate(img%apb(2,2,1,4),img%amb(2,2,1,1),source=0.0_dp)
 call update(img,3)
 if(img%apb(1,1,1,3)/=3.or.img%amb(1,1,1,1)/=-3)error stop 'image slots'
end program
'''
    _run(tmp_path, program)


def test_closed_shell_gradient_beta_pointer_disassociated(tmp_path):
    source = (SOURCE / 'dftlib' / 'dft_gridint_grad.F90').read_text()
    routine = _routine(source, 'subroutine', 'derexc_blk')
    # the beta scratch density db2 is allocated only for urohf
    assert re.search(r'^ *xc_opts%wfBeta => db2', routine, re.M) is None, \
        'closed-shell derexc_blk must not associate wfBeta with unallocated db2'
    block = re.search(r'^ *xc_opts%wfAlpha => da2\n *if \(urohf\) xc_opts%wfBeta => db2\n',
                      routine, re.M)
    assert block, 'derexc_blk wfAlpha/wfBeta association block not found'
    program = '''
module fixture
 implicit none
 integer,parameter::fp=kind(1d0)
 type opts_t
  real(fp),pointer::wfAlpha(:,:)=>null(),wfBeta(:,:)=>null()
 end type
contains
 logical function beta_associated(urohf)
  logical,intent(in)::urohf
  real(fp),target,allocatable::da2(:,:),db2(:,:)
  type(opts_t)::xc_opts
  allocate(da2(3,3),source=1.0_fp)
  if(urohf)allocate(db2(3,3),source=2.0_fp)
''' + block.group() + '''
  beta_associated=associated(xc_opts%wfBeta)
  if(beta_associated)beta_associated=all(xc_opts%wfBeta==2.0_fp)
 end function
end module
program check
 use fixture
 implicit none
 if(beta_associated(.false.))error stop 'closed-shell wfBeta associated'
 if(.not.beta_associated(.true.))error stop 'open-shell wfBeta missing'
end program
'''
    _run(tmp_path, program)
