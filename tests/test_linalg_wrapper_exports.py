"""Compile the production wrapper modules with the supported ILP64 width."""
from pathlib import Path
import re
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def compiled_consumer(tmp_path):
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    assert compiler, 'GNU Fortran is required to check the ILP64 wrapper interface'
    # Only the diagnostic service is replaced. Integer kinds, wrappers and
    # the public module interface are compiled directly from production source.
    messages = tmp_path / 'messages.f90'
    messages.write_text('''module messages
 implicit none
 integer, parameter :: WITH_ABORT=1
contains
 subroutine show_message(text, level)
  character(*), intent(in) :: text
  integer, intent(in) :: level
  error stop 'unexpected diagnostic'
 end subroutine
end module
''')
    consumer = tmp_path / 'consumer.f90'
    source = (ROOT / 'source/modules/casscf_driver.F90').read_text()
    canonicalize = re.search(r'  subroutine canonicalize\(.*?  end subroutine canonicalize',
                             source, re.S).group()
    cbuf_decl = re.search(r'^    real\(dp\).*:: cbuf\(.*$', canonicalize, re.M).group()
    storage_decl = re.search(r'^    real\(dp\).*:: cstorage\(.*$', canonicalize, re.M).group()
    remap = re.search(r'^      cstorage\(.*=>.*$', canonicalize, re.M).group()
    rotation = re.search(r"      call dgemm\('T', 'N', m, n, m,.*?blk, m\)",
                         canonicalize, re.S).group()
    consumer.write_text('''subroutine consumer(a,b,c)
 use oqp_linalg, only: dgemm
 implicit none
 double precision :: a(2,3), b(3,4), c(2,4)
 call dgemm('N','N',2,4,3,1d0,a,2,b,3,0d0,c,2)
end subroutine
subroutine rotation(cbuf,sub,blk,n,m,lo,hi)
 use oqp_linalg, only: dgemm
 use precision, only: dp
 implicit none
 integer :: n,m,lo,hi
''' + cbuf_decl + '\n' + storage_decl + '''
 real(dp) :: sub(m,m),blk(m,n)
''' + remap + '\n' + rotation + '''
end subroutine
''')
    sources = [ROOT / 'source/precision.F90', messages,
               ROOT / 'source/mathlib/mathlib_types.F90',
               ROOT / 'source/mathlib/blas_wrap.F90',
               ROOT / 'source/mathlib/lapack_wrap.F90',
               ROOT / 'source/mathlib/oqp_linalg.F90']
    result = subprocess.run(
        [compiler, '-cpp', '-DOQP_BLAS_INT=8', '-fdefault-integer-8',
         '-fcheck=all', '-c', *map(str, sources)], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    # The packing contract belongs to the changed CASSCF operand, not
    # unchanged wrappers such as DGETRI's legacy pivot association.
    result = subprocess.run([compiler, '-fdefault-integer-8', '-fcheck=all',
                             '-Warray-temporaries', '-c', str(consumer)],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'array temporary' not in result.stderr.lower(), result.stderr
    return compiler, tmp_path / 'consumer.o'


def test_ilp64_dgemm_is_exported_by_oqp_linalg(compiled_consumer):
    _, consumer = compiled_consumer
    symbols = subprocess.check_output(['nm', str(consumer)], text=True)
    assert '__blas_wrap_MOD_oqp_dgemm_i64' in symbols



def test_native_rectangular_casscf_block_with_original_stride(compiled_consumer, tmp_path):
    from oqp.runtime import library_path, resolve_oqp_root
    root, suffix = resolve_oqp_root()
    library = library_path(root, suffix)
    assert library, 'installed native OpenQP required'
    compiler, consumer = compiled_consumer
    program = tmp_path / 'rectangular.f90'
    program.write_text("""program rectangular
 use precision, only: dp
 implicit none
 integer :: k,n,m,lo,hi,i,j
 integer, parameter :: ns(4)=[5,7,4,5], ms(4)=[2,3,1,2], ls(4)=[1,2,3,0]
 real(dp), allocatable :: cbuf(:,:), original(:,:), sub(:,:), blk(:,:), expected(:,:)
 interface
  subroutine rotation(cbuf,sub,blk,n,m,lo,hi)
   use precision, only: dp
   integer :: n,m,lo,hi
   real(dp), contiguous, target :: cbuf(0:,0:)
   real(dp) :: sub(m,m),blk(m,n)
  end subroutine
 end interface
 do k=1,4
  n=ns(k);m=ms(k);lo=ls(k);hi=lo+m-1
  allocate(cbuf(0:n-1,0:n-1), original(0:n-1,0:n-1), sub(m,m),blk(m,n),expected(m,n))
  do j=0,n-1
   do i=0,n-1
    cbuf(i,j)=sin(real(7*i+3*j+1,dp))+0.03_dp*i*j
   end do
  end do
  do j=1,m
   do i=1,m
    sub(i,j)=0.13_dp*i-0.07_dp*j+cos(real(i+2*j,dp))
   end do
  end do
  original=cbuf
  expected=matmul(transpose(sub),cbuf(lo:hi,:))
  call rotation(cbuf,sub,blk,n,m,lo,hi)
  if(maxval(abs(blk-expected))>1e-12_dp) error stop 'rectangular stride mismatch'
  if(any(cbuf/=original)) error stop 'borrowed input changed'
  deallocate(cbuf,original,sub,blk,expected)
 end do
end program
""")
    executable = tmp_path / 'rectangular'
    subprocess.run([compiler, '-fdefault-integer-8', '-fcheck=all', str(program),
                    str(consumer), str(library), '-Wl,-rpath,' + str(library.parent),
                    '-o', str(executable)], cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], cwd=tmp_path, check=True, capture_output=True, text=True)


def _production_routine(path, name):
    source = path.read_text()
    match = re.search(r'^\s*subroutine ' + name + r'\(.*?^\s*end subroutine(?: ' + name + r')?\s*$',
                      source, re.I | re.M | re.S)
    assert match, name
    return match.group()


@pytest.mark.parametrize('width', ['1', '4'])
def test_production_response_matrix_storage(compiled_consumer, tmp_path, width):
    """Execute complete changed routines, with strided parents and sentinels.

    Only their surrounding data types/diagnostic service are reduced. Small
    matrices pin array association; molecular validation remains separate.
    """
    import os
    from oqp.runtime import library_path, resolve_oqp_root
    root, suffix = resolve_oqp_root()
    library = library_path(root, suffix)
    assert library, 'installed native OpenQP required'
    compiler, _ = compiled_consumer
    support = tmp_path / 'support.f90'
    support.write_text('''module types
 type basis_t
  integer :: nbf
 end type
 type prop_t
  integer :: nelec_a,nelec_b
 end type
 type tddft_t
  integer :: mult
  logical :: debug_mode
 end type
 type information
  type(basis_t) :: basis
  type(prop_t) :: mol_prop
  type(tddft_t) :: tddft
 end type
end module
module io_constants
 integer,parameter :: iw=6
end module
module mathlib
 use precision,only: dp
 implicit none
contains
 subroutine orthogonal_transform(trans,n,u,a,b,wrk)
  character :: trans
  integer :: n
  real(dp) :: u(n,n),a(n,n),b(n,n),wrk(n,n)
  if(trans/='t') error stop 'unexpected transform in contract test'
  b=matmul(matmul(u,a),transpose(u))
 end subroutine
 subroutine pack_matrix(a,ap)
  real(dp) :: a(:,:),ap(:)
  integer :: i,j,k
  k=0
  do j=1,size(a,2)
   do i=1,j
    k=k+1;ap(k)=a(i,j)
   end do
  end do
 end subroutine
 subroutine triangular_to_full(a,n,uplo)
  integer :: n,i
  character :: uplo
  real(dp) :: a(n,n)
  if(uplo/='u') error stop 'unexpected triangle in contract test'
  do i=1,n
   a(i+1:,i)=a(i,i+1:)
  end do
 end subroutine
end module
''')
    mrsf = ROOT / 'source/tdhf_mrsf_lib.F90'
    routines = [_production_routine(mrsf, name) for name in
                ('mrsfcbc', 'umrsfcbc', 'mrsfmntoia', 'umrsfmntoia', 'mrsfesum')]
    routines += [_production_routine(ROOT / 'source/tdhf_lib.F90', 'tdhf_unrelaxed_density'),
                 _production_routine(ROOT / 'source/modules/tdhf_z_vector.F90', 'compute_w_mo')]
    # The same unchanged ILP64 DSYR2K wrapper pins the output offset too;
    # the package build separately executes production's registered raw alias.
    first = routines[0]
    ov = re.search(r'    tmp\(:,1\) = 0\.0_dp\n    if \(nbf > nocca\) then\n.*?    end if',
                   first, re.S).group()
    matrices = tmp_path / 'matrices.f90'
    matrices.write_text('''module matrix_contracts
 use precision,only: dp
 use oqp_linalg,only: dgemm
 use blas_wrap,only: dsyr2k=>oqp_dsyr2k_i64
 implicit none
contains
''' + '\n'.join(routines) + '''
 subroutine ov_contract(vb,bvec,tmp,nbf,nocca,lr2)
  integer :: nbf,nocca,lr2
  real(dp),target,intent(in) :: vb(nbf,nbf),bvec(nbf,nbf)
  real(dp) :: tmp(nbf,1)
  real(dp),contiguous,pointer :: vb_storage(:),bvec_storage(:)
  vb_storage(1:size(vb))=>vb
  bvec_storage(1:size(bvec))=>bvec
''' + ov + '''
 end subroutine
end module
''')
    executable = tmp_path / 'matrix_views'
    result = subprocess.run([compiler, '-fdefault-integer-8', '-fcheck=all',
                             '-ffree-line-length-none', str(support), str(matrices),
                             str(tmp_path / 'messages.o'),
                             str(ROOT / 'tests/native/linalg_matrix_views.F90'), str(library),
                             '-Wl,-rpath,' + str(library.parent), '-o', str(executable)],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=width, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    result = subprocess.run([str(executable)], env=env, cwd=tmp_path,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'MATRIX_STORAGE_CONTRACT_PASS' in result.stdout
