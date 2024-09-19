! buffer_types defines a list of data types used in parallel communication.
! Each entry specifies:
! (1) A human-readable name for the data type,
! (2) The corresponding MPI data type for communication,
! (3) The Fortran data type (with the appropriate kind),
! (4) The array rank or shape (e.g., scalar, 1D, 2D, etc.).
! These types are used for creating allreduce/bcast buffers in MPI operations.


#:set buffer_types = [&
   &("int32_scalar",  "MPI_INTEGER",      "integer(int32)",         ""),&
   &("int32_1d",      "MPI_INTEGER",      "integer(int32)",         "(:)"),&
   &("int64_scalar",  "MPI_INTEGER8",     "integer(int64)",         ""),&
   &("int64_1d",      "MPI_INTEGER8",     "integer(int64)",         "(:)"),&
   &("dp_scalar",     "MPI_DOUBLE_PRECISION", "real(kind=dp)",      ""),&
   &("dp_1d",         "MPI_DOUBLE_PRECISION", "real(kind=dp)",      "(:)"),&
   &("dp_2d",         "MPI_DOUBLE_PRECISION", "real(kind=dp)",      "(:,:)"),&
   &("dp_3d",         "MPI_DOUBLE_PRECISION", "real(kind=dp)",      "(:,:,:)"),&
   &("dp_4d",         "MPI_DOUBLE_PRECISION", "real(kind=dp)",      "(:,:,:,:)"),&
   &("byte",          "MPI_BYTE",         "character(kind=c_char, len=1)", "(*)"),&
   &("c_bool",        "MPI_C_BOOL",       "logical(c_bool)",        "")]

module parallel
  use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real64
  use precision, only: fp, dp
  use iso_c_binding, only: c_char, c_bool
#ifdef ENABLE_MPI
  use mpi
#endif
  implicit none

#ifndef ENABLE_MPI
  integer, parameter :: MPI_COMM_NULL = 0
  integer, parameter :: PARALLEL_INT = int32
#else
  integer, parameter :: PARALLEL_INT = MPI_INTEGER_KIND
#endif

  private
  public :: par_env_t
  public :: PARALLEL_INT
  public :: MPI_COMM_NULL


  type :: par_env_t

    integer(PARALLEL_INT) :: comm = MPI_COMM_NULL
    integer(PARALLEL_INT) :: rank = 0
    integer(PARALLEL_INT) :: size = 1
    integer(PARALLEL_INT) :: err = 0

    logical(c_bool) :: use_mpi = .false.

    contains

      procedure, pass(self) :: init => par_env_t_init
      procedure, pass(self) :: barrier => par_env_t_barrier
      procedure, pass(self) :: get_hostnames

#:for type_name, mpi_type, kind, rank in buffer_types
      procedure, pass(self) :: par_env_t_bcast_${type_name}$
      procedure, pass(self) :: par_env_t_allreduce_${type_name}$
#:endfor

      generic :: bcast => &
#:for type_name, mpi_type, kind, rank in buffer_types
              par_env_t_bcast_${type_name}$#{if type_name != buffer_types[-1][0]}#,&#{endif}#
#:endfor

      generic :: allreduce => &
#:for type_name, mpi_type, kind, rank in buffer_types
             par_env_t_allreduce_${type_name}$#{if type_name != buffer_types[-1][0]}#,&#{endif}#
#:endfor

end type par_env_t

contains

  subroutine par_env_t_init(self, comm, use_mpi)
    class(par_env_t), intent(inout) :: self
    integer(PARALLEL_INT), intent(in) :: comm
    logical(c_bool), intent(in) :: use_mpi

    self%use_mpi = use_mpi
#ifdef ENABLE_MPI
    self%comm = comm
    if (.not. use_mpi) return
    call MPI_Comm_rank(comm, self%rank, self%err)
    call MPI_Comm_size(comm, self%size, self%err)
#else
    self%use_mpi = .FALSE.
    self%rank = 0
    self%err = 0
#endif
  end subroutine par_env_t_init

  subroutine par_env_t_barrier(self)
    class(par_env_t) :: self
#ifndef ENABLE_MPI
    return
#else
    if (self%use_mpi) then
      call MPI_Barrier(self%comm, self%err)
    endif
#endif
  end subroutine par_env_t_barrier
!##################################
  subroutine get_hostnames(self, node_info_str)
    implicit none

    character(len=:), allocatable, intent(inout) :: node_info_str
    class(par_env_t), intent(inout) :: self
    integer(PARALLEL_INT), PARAMETER :: max_host_print = 4
#ifdef ENABLE_MPI
    integer(PARALLEL_INT) :: node_comm, local_rank, node_root_comm, ierr
    integer(PARALLEL_INT) :: num_nodes
    character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
    integer(PARALLEL_INT) :: node_info_len
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: gathered_node_names(:)
    integer :: i
    character(len=32) :: num_nodes_str
#endif
    integer(PARALLEL_INT) :: name_len
    name_len = 28
#ifndef ENABLE_MPI
    allocate(character(len=name_len) :: node_info_str)
    call hostnm(node_info_str)
#else
    if(self%use_mpi) then
      call MPI_Comm_split_type(self%comm, MPI_COMM_TYPE_SHARED, int(0, kind=PARALLEL_INT), MPI_INFO_NULL, node_comm, ierr)
      call MPI_Comm_rank(node_comm, local_rank, ierr)
      call MPI_Comm_split(self%comm, int(merge(0, 1, local_rank == 0), kind=PARALLEL_INT), self%rank, node_root_comm, ierr)
      if (local_rank == 0) then
        call MPI_Comm_size(node_root_comm, num_nodes, ierr)
      end if

      call MPI_Bcast(num_nodes, int(1, kind=PARALLEL_INT), MPI_INTEGER, int(0, kind=PARALLEL_INT), self%comm, ierr)
      call MPI_Get_processor_name(hostname, name_len, ierr)
      if (local_rank == 0 .and. num_nodes <= max_host_print) then
        allocate(gathered_node_names(num_nodes))
        call MPI_Gather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, gathered_node_names, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER,&
                        int(0, kind=PARALLEL_INT) , node_root_comm, ierr)
      end if

      if (self%rank == 0) then
        if (num_nodes > max_host_print) then
          write(num_nodes_str, '(I0)') num_nodes
          node_info_str = "Number of nodes " // trim(num_nodes_str)
        else
          node_info_str = ''
          do i = 1, num_nodes
            node_info_str = trim(node_info_str) // trim(gathered_node_names(i)) // ', '
          end do
          node_info_str = node_info_str(1:len_trim(node_info_str)-1)
        end if
      end if

      if (self%rank == 0) then
        node_info_len = len_trim(node_info_str)
      endif
      call MPI_Bcast(node_info_len, int(1, kind=PARALLEL_INT), MPI_INTEGER, int(0, kind=PARALLEL_INT), self%comm, ierr)
      if (self%rank /= 0) then
        allocate(character(len=node_info_len) :: node_info_str)
      end if
      call MPI_Bcast(node_info_str, node_info_len, MPI_CHARACTER, int(0, kind=PARALLEL_INT), self%comm, ierr)

      if (local_rank == 0 .and. num_nodes <= max_host_print) then
        deallocate(gathered_node_names)
        call MPI_Comm_free(node_root_comm, ierr)
      end if
      call MPI_Comm_free(node_comm, ierr)

    else
      allocate(character(len=name_len) :: node_info_str)
      call hostnm(node_info_str)
    endif
#endif
  end subroutine get_hostnames
!##################################
#:for type_name, mpi_type, kind, rank in buffer_types
  subroutine par_env_t_bcast_${type_name}$(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    ${kind}$, intent(inout) :: buffer${rank}$
    integer, intent(in) :: length
    integer(PARALLEL_INT), optional, intent(in) :: root
    integer(PARALLEL_INT) :: root_
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    root_ = 0
    if (present(root)) root_ = root
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), ${mpi_type}$, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_${type_name}$

  subroutine par_env_t_allreduce_${type_name}$(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    ${kind}$, intent(inout) :: buffer${rank}$
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), ${mpi_type}$, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_${type_name}$

#:endfor

end module parallel
