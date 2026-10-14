"""Compile the production rank-budget method with deterministic MPI services.

This standalone diagnostic is not an installed-package or real-MPI test. It
exercises processor-identity grouping across several node layouts and errors.
The shared service deliberately lies or fails, reproducing an unsafe singleton
answer on co-located ranks without relying on a particular MPI vendor.
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
 integer(int32), parameter :: MPI_MAX_PROCESSOR_NAME=64, MPI_CHARACTER=1
 integer(int32) :: ranks=1, size_error=0, split_error=0
 integer(int32) :: layout=0, caller=0, name_error=0, bad_bcast=-1
 integer(int32) :: name_length=-1, blank_root=-1, expected_comm=MPI_COMM_WORLD
 integer(int32) :: broadcast_calls=0, split_calls=0, free_calls=0
contains
 subroutine MPI_Comm_size(comm,n,ierr)
  integer(int32), intent(in) :: comm
  integer(int32), intent(out) :: n,ierr
  n=1
  if(comm==expected_comm) n=ranks
  ierr=size_error
 end subroutine
 subroutine MPI_Comm_split_type(comm,kind,key,info,node,ierr)
  integer(int32), intent(in) :: comm,kind,key,info
  integer(int32), intent(out) :: node,ierr
  node=30
  ierr=split_error
  split_calls=split_calls+1
 end subroutine
 subroutine MPI_Comm_free(comm,ierr)
  integer(int32), intent(inout) :: comm
  integer(int32), intent(out) :: ierr
  ierr=0
  free_calls=free_calls+1
 end subroutine
 function processor(rank) result(name)
  integer(int32), intent(in) :: rank
  character(MPI_MAX_PROCESSOR_NAME) :: name
  select case(layout)
  case(1)
   write(name,'(a,i0)') 'node-',rank/4
  case(2)
   write(name,'(a,i0)') 'node-',mod(rank,3_int32)
  case(3)
   name='node-AA'
   if(rank==0) name='node-A'
  case(4)
   name=repeat('A',int(MPI_MAX_PROCESSOR_NAME))
  case default
   name='node-A'
  end select
 end function
 subroutine MPI_Get_processor_name(name,length,ierr)
  character(*), intent(out) :: name
  integer(int32), intent(out) :: length,ierr
  character(MPI_MAX_PROCESSOR_NAME) :: actual
  actual=processor(caller)
  if(layout==5) actual='absent'
  length=int(len_trim(actual),int32)
  ! The MPI service need not initialize bytes beyond the returned length.
  name=repeat('#',len(name))
  name(1:length)=actual(1:length)
  if(name_length>=0) length=name_length
  if(name_length==-2) name=' '
  ierr=name_error
 end subroutine
 subroutine MPI_Bcast(name,count,datatype,root,comm,ierr)
  character(*), intent(inout) :: name
  integer(int32), intent(in) :: count,datatype,root,comm
  integer(int32), intent(out) :: ierr
  character(MPI_MAX_PROCESSOR_NAME) :: expected_local
  if(count/=MPI_MAX_PROCESSOR_NAME.or.datatype/=MPI_CHARACTER.or.comm/=expected_comm) &
   error stop 'wrong broadcast ABI/communicator'
  if(root/=broadcast_calls) error stop 'inconsistent collective order'
  expected_local=processor(caller)
  if(layout==5) expected_local='absent'
  if(name/=expected_local) error stop 'uninitialized or unpadded local name'
  broadcast_calls=broadcast_calls+1
  name=processor(root)
  if(root==blank_root) name=' '
  ierr=MPI_SUCCESS
  if(root==bad_bcast) ierr=7
 end subroutine
 subroutine MPI_Abort(comm,code,ierr)
  integer(int32), intent(in) :: comm,code
  integer(int32), intent(out) :: ierr
  if(comm/=MPI_COMM_WORLD.or.code/=1) error stop 'wrong abort target'
  print *, 'WORLD_ABORT'
  print *, 'BROADCASTS',broadcast_calls
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
 integer :: usempi,n,repeat_mode
 call get_command_argument(1,arg)
 read(arg,*) ranks
 call get_command_argument(2,arg)
 read(arg,*) size_error
 call get_command_argument(3,arg)
 read(arg,*) split_error
 call get_command_argument(4,arg)
 read(arg,*) usempi
 pe%use_mpi=usempi/=0
 if(command_argument_count()>=5) then
  call get_command_argument(5,arg)
  read(arg,*) layout
 endif
 if(command_argument_count()>=6) then
  call get_command_argument(6,arg)
  read(arg,*) caller
 endif
 if(command_argument_count()>=7) then
  call get_command_argument(7,arg)
  read(arg,*) name_error
 endif
 if(command_argument_count()>=8) then
  call get_command_argument(8,arg)
  read(arg,*) bad_bcast
 endif
 if(command_argument_count()>=9) then
  call get_command_argument(9,arg)
  read(arg,*) name_length
 endif
 if(command_argument_count()>=10) then
  call get_command_argument(10,arg)
  read(arg,*) blank_root
 endif
 if(command_argument_count()>=11) then
  call get_command_argument(11,arg)
  read(arg,*) expected_comm
  pe%comm=expected_comm
 endif
 n=pe%local_size()
 if(split_calls/=0.or.free_calls/=0) error stop 'unexpected communicator lifetime'
 print *, 'BROADCASTS',broadcast_calls
 print *, 'COUNT',n
 if(command_argument_count()>=12) then
  call get_command_argument(12,arg)
  read(arg,*) repeat_mode
  if(repeat_mode==1) then
   ! Reuse the same environment while node layouts grow and shrink.
   broadcast_calls=0
   ranks=16
   layout=1
   caller=15
   n=pe%local_size()
   if(n/=4.or.broadcast_calls/=16) error stop 'grown layout'
   broadcast_calls=0
   ranks=2
   layout=0
   caller=0
   n=pe%local_size()
   if(n/=2.or.broadcast_calls/=2) error stop 'shrunk layout'
   if(split_calls/=0.or.free_calls/=0) error stop 'repeated communicator lifetime'
   print *, 'REPEATED_LAYOUTS_OK'
  endif
 endif
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


@pytest.mark.parametrize('ranks,layout,caller,expected', [
    (16, 1, 0, 4), (16, 1, 15, 4),
    (7, 2, 0, 3), (7, 2, 1, 2), (7, 2, 6, 3),
    (4, 3, 0, 1), (4, 3, 1, 3), (4, 4, 2, 4),
])
def test_only_matching_processor_identities_count(budget_programs, ranks, layout, caller, expected):
    result = subprocess.run([str(budget_programs[True]), str(ranks), '0', '7', '1',
                             str(layout), str(caller)], capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(rf'COUNT\s+{expected}\s*$', result.stdout)
    assert re.search(rf'BROADCASTS\s+{ranks}\s', result.stdout)


@pytest.mark.parametrize('name_error,bad_bcast,length,blank_root', [
    (7, -1, -1, -1), (0, -1, 0, -1), (0, -1, 65, -1), (0, -1, -2, -1),
    (0, 0, -1, -1), (0, 3, -1, -1), (0, -1, -1, 2),
])
def test_processor_collective_errors_terminate_world(budget_programs, name_error, bad_bcast,
                                                     length, blank_root):
    result = subprocess.run([str(budget_programs[True]), '4', '0', '0', '1', '0', '0',
                             str(name_error), str(bad_bcast), str(length), str(blank_root)],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert 'WORLD_ABORT' in result.stdout
    assert 'cannot determine MPI rank count' in result.stderr
    assert 'COUNT' not in result.stdout


def test_grouping_uses_the_callers_communicator(budget_programs):
    result = subprocess.run([str(budget_programs[True]), '7', '0', '0', '1', '2', '1',
                             '0', '-1', '-1', '-1', '99'],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r'COUNT\s+2\s*$', result.stdout)


def test_missing_own_processor_identity_terminates_world(budget_programs):
    result = subprocess.run([str(budget_programs[True]), '4', '0', '0', '1', '5', '0'],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert 'WORLD_ABORT' in result.stdout
    assert re.search(r'BROADCASTS\s+4\s', result.stdout)
    assert 'cannot determine MPI rank count' in result.stderr
    assert 'COUNT' not in result.stdout


def test_repeated_calls_follow_changed_node_layouts(budget_programs):
    result = subprocess.run([str(budget_programs[True]), '7', '0', '7', '1', '2', '1',
                             '0', '-1', '-1', '-1', '99', '1'],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'REPEATED_LAYOUTS_OK' in result.stdout
