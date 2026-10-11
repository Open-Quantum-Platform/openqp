"""Compile the production rank-budget method with deterministic MPI services.

This standalone diagnostic is not an installed-package or real-MPI test. The
shared service deliberately lies or fails, reproducing an unsafe singleton
answer on a co-located communicator without relying on a particular MPI vendor.
"""
from pathlib import Path
import re
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def compile_budget(directory, source, enable_mpi):
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    assert compiler, 'GNU Fortran is required for the rank-budget diagnostic'
    method = re.search(r'  function par_env_t_local_size\(.*?  end function par_env_t_local_size',
                       source, re.S).group()
    services = '''module mpi
 use iso_fortran_env, only: int32
 implicit none
 integer, parameter :: MPI_INTEGER_KIND=int32
 integer(int32), parameter :: MPI_COMM_WORLD=10, MPI_COMM_TYPE_SHARED=20
 integer(int32), parameter :: MPI_INFO_NULL=0, MPI_SUCCESS=0
 integer(int32) :: ranks=1, size_error=0, split_error=0
contains
 subroutine MPI_Comm_size(comm,n,ierr)
  integer(int32), intent(in) :: comm
  integer(int32), intent(out) :: n,ierr
  n=1
  if(comm==MPI_COMM_WORLD) n=ranks
  ierr=size_error
 end subroutine
 subroutine MPI_Comm_split_type(comm,kind,key,info,node,ierr)
  integer(int32), intent(in) :: comm,kind,key,info
  integer(int32), intent(out) :: node,ierr
  node=30
  ierr=split_error
 end subroutine
 subroutine MPI_Comm_free(comm,ierr)
  integer(int32), intent(inout) :: comm
  integer(int32), intent(out) :: ierr
  ierr=0
 end subroutine
 subroutine MPI_Abort(comm,code,ierr)
  integer(int32), intent(in) :: comm,code
  integer(int32), intent(out) :: ierr
  if(comm/=MPI_COMM_WORLD.or.code/=1) error stop 'wrong abort target'
  print *, 'WORLD_ABORT'
  ierr=0 ! Return deliberately to exercise the terminal fallback too.
 end subroutine
end module
'''
    unit = directory / 'budget.F90'
    unit.write_text(services + '''module budget
 use iso_fortran_env, only: int32
 use mpi
 implicit none
 integer, parameter :: PARALLEL_INT=int32
 type par_env_t
  integer(int32) :: comm=MPI_COMM_WORLD
  logical :: use_mpi=.true.
 contains
  procedure :: local_size=>par_env_t_local_size
 end type
contains
''' + method + '''
end module
program check
 use budget
 implicit none
 type(par_env_t) :: pe
 character(32) :: arg
 integer :: usempi,n
 call get_command_argument(1,arg)
 read(arg,*) ranks
 call get_command_argument(2,arg)
 read(arg,*) size_error
 call get_command_argument(3,arg)
 read(arg,*) split_error
 call get_command_argument(4,arg)
 read(arg,*) usempi
 pe%use_mpi=usempi/=0
 n=pe%local_size()
 print *, 'COUNT',n
end program
''')
    executable = directory / 'budget'
    command = [compiler, '-cpp', '-fdefault-integer-8', '-fcheck=all']
    if enable_mpi:
        command.append('-DENABLE_MPI')
    result = subprocess.run(command + [str(unit), '-o', str(executable)],
                            cwd=directory, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return executable


@pytest.fixture(scope='module')
def budget_programs(tmp_path_factory):
    source = (ROOT / 'source/parallel.F90').read_text()
    return {mpi: compile_budget(tmp_path_factory.mktemp(f'budget-{mpi}'), source, mpi)
            for mpi in (False, True)}


@pytest.mark.parametrize('ranks,split_error', [(1, 0), (2, 0), (4, 0), (2, 7), (4, 7)])
def test_shared_singleton_or_failure_cannot_underestimate(budget_programs, ranks, split_error):
    result = subprocess.run([str(budget_programs[True]), str(ranks), '0', str(split_error), '1'],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(rf'COUNT\s+{ranks}\s*$', result.stdout)


@pytest.mark.parametrize('ranks,error', [(2, 7), (0, 0), (-1, 0)])
def test_unknown_rank_count_terminates_world(budget_programs, ranks, error):
    result = subprocess.run([str(budget_programs[True]), str(ranks), str(error), '0', '1'],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert 'WORLD_ABORT' in result.stdout
    assert 'cannot determine MPI rank count' in result.stderr
    assert 'COUNT' not in result.stdout


@pytest.mark.parametrize('mpi,usempi', [(False, 1), (True, 0)])
def test_non_mpi_budget_is_one(budget_programs, mpi, usempi):
    result = subprocess.run([str(budget_programs[mpi]), '4', '7', '7', str(usempi)],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r'COUNT\s+1\s*$', result.stdout)
