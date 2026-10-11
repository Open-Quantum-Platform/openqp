"""Compile the production wrapper modules with the supported ILP64 width."""
from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def test_ilp64_dgemm_is_exported_by_oqp_linalg(tmp_path):
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
 real(dp), contiguous :: cbuf(0:,0:)
 real(dp) :: sub(m,m),blk(m,n)
''' + rotation + '''
end subroutine
''')
    sources = [ROOT / 'source/precision.F90', messages,
               ROOT / 'source/mathlib/mathlib_types.F90',
               ROOT / 'source/mathlib/blas_wrap.F90',
               ROOT / 'source/mathlib/lapack_wrap.F90',
               ROOT / 'source/mathlib/oqp_linalg.F90', consumer]
    result = subprocess.run(
        [compiler, '-cpp', '-DOQP_BLAS_INT=8', '-fdefault-integer-8',
         '-c', *map(str, sources)], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    symbols = subprocess.check_output(['nm', str(tmp_path / 'consumer.o')], text=True)
    assert '__blas_wrap_MOD_oqp_dgemm_i64' in symbols
