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
    rotation = re.search(r"      call dgemm\('T', 'N', m, n, m,.*?blk, m\)", source, re.S).group()
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
 real(dp), contiguous, target :: cbuf(0:,0:)
 real(dp), contiguous, pointer :: cstorage(:)
 real(dp) :: sub(m,m),blk(m,n)
 cstorage(1:n*n) => cbuf
''' + rotation + '''
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
