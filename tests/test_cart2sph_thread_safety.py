"""The spherical-harmonic transform tables must be safe to use from threads.

source/integrals/cart2sph.F90 is called from OpenMP integral loops
(kin_ovl_ints for S and T, the ERI shell quartets), and the first call in a
process is normally made there by several threads at once. The l=5/6 tables
used to be built lazily into SAVE variables with an unsynchronised
"if (.not. allocated(t)) t = build()": a second thread could read the table
after its allocation but before its copy. aug-cc-pV5Z Li then gave wrong S and
T, NaN SCF energies and wrong MP2 energies intermittently at OMP_NUM_THREADS=4.

The module is compiled from production source. The first test pins the h/i
literal tables; the second fails deterministically on any static writable
storage in the module, which is the cause; the third reproduces the race
itself (61 of 500 fresh processes failed at 32 threads on Intel macOS before
the fix, but only about 0.5% on a 2-core Linux runner).
"""
from pathlib import Path
import os
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'source/precision.F90', ROOT / 'source/constants.F90',
           ROOT / 'source/integrals/cart2sph.F90']
THREADS = 32
TRIALS = 400

DRIVER = '''program c2s_first_use
  use precision, only: dp
  use cart2sph, only: c2s_selftest
  use omp_lib, only: omp_get_thread_num, omp_get_num_threads
  implicit none
  real(dp) :: err(0:255)
  integer :: t, nt
  err = huge(1.0_dp)
  nt = 0
!$omp parallel private(t)
  t = omp_get_thread_num()
!$omp single
  nt = omp_get_num_threads()
!$omp end single
!$omp barrier
  call c2s_selftest(err(t))
!$omp end parallel
  do t = 0, nt - 1
    if (.not. (err(t) < 1.0e-13_dp)) then
      print '(a,i0,a,es12.4)', 'thread ', t, ' c2s error ', err(t)
      error stop 1
    end if
  end do
  print '(a,i0)', 'threads ', nt
end program
'''


@pytest.fixture(scope='module')
def compiled(tmp_path_factory):
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    assert compiler, 'GNU Fortran is required to compile the production c2s tables'
    work = tmp_path_factory.mktemp('cart2sph')
    driver = work / 'driver.F90'
    driver.write_text(DRIVER)
    flags = ['-fdefault-integer-8', '-ffree-line-length-none', '-fopenmp', '-O2']
    result = subprocess.run([compiler, *flags, '-c', *map(str, SOURCES)],
                            cwd=work, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    result = subprocess.run([compiler, *flags, 'precision.o', 'constants.o', 'cart2sph.o',
                             str(driver), '-o', 'c2s_first_use'],
                            cwd=work, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return work


def _run(work, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='false')
    return subprocess.run([str(work / 'c2s_first_use')], cwd=work, env=env,
                          capture_output=True, text=True)


def test_h_and_i_tables_match_generator_and_are_orthonormal(compiled):
    run = _run(compiled, 1)
    assert run.returncode == 0, run.stdout + run.stderr
    assert 'threads 1' in run.stdout


def test_cart2sph_has_no_writable_static_storage(compiled):
    nm = shutil.which('nm')
    assert nm, 'nm is required to inspect the cart2sph object'
    if sys.platform == 'darwin':
        listing = subprocess.check_output([nm, '-m', 'cart2sph.o'], cwd=compiled, text=True)
        writable = [line for line in listing.splitlines() if '(__DATA,' in line]
    else:
        listing = subprocess.check_output([nm, 'cart2sph.o'], cwd=compiled, text=True)
        writable = [line for line in listing.splitlines()
                    if len(line.split()) >= 2 and line.split()[-2] in 'bBdDCgGsS']
    assert not writable, ('cart2sph is called from OpenMP integral loops and must not '
                          'keep static mutable state: ' + '; '.join(writable))


def test_concurrent_first_use_of_h_and_i_tables(compiled):
    failures = []
    for trial in range(TRIALS):
        run = _run(compiled, THREADS)
        if run.returncode != 0 or ('threads %d' % THREADS) not in run.stdout:
            failures.append((trial, (run.stdout + run.stderr).strip()))
    assert not failures, '%d of %d fresh processes read a corrupt table, first: %s' % (
        len(failures), TRIALS, failures[0])
