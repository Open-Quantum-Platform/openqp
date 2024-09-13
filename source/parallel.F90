! buffer_types defines a list of data types used in parallel communication.
! Each entry specifies:
! (1) A human-readable name for the data type,
! (2) The corresponding MPI data type for communication,
! (3) The Fortran data type (with the appropriate kind),
! (4) The array rank or shape (e.g., scalar, 1D, 2D, etc.).
! These types are used for creating allreduce/bcast buffers in MPI operations.



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

      procedure, pass(self) :: par_env_t_bcast_int32_scalar
      procedure, pass(self) :: par_env_t_allreduce_int32_scalar
      procedure, pass(self) :: par_env_t_bcast_int32_1d
      procedure, pass(self) :: par_env_t_allreduce_int32_1d
      procedure, pass(self) :: par_env_t_bcast_int64_scalar
      procedure, pass(self) :: par_env_t_allreduce_int64_scalar
      procedure, pass(self) :: par_env_t_bcast_int64_1d
      procedure, pass(self) :: par_env_t_allreduce_int64_1d
      procedure, pass(self) :: par_env_t_bcast_dp_scalar
      procedure, pass(self) :: par_env_t_allreduce_dp_scalar
      procedure, pass(self) :: par_env_t_bcast_dp_1d
      procedure, pass(self) :: par_env_t_allreduce_dp_1d
      procedure, pass(self) :: par_env_t_bcast_dp_2d
      procedure, pass(self) :: par_env_t_allreduce_dp_2d
      procedure, pass(self) :: par_env_t_bcast_dp_3d
      procedure, pass(self) :: par_env_t_allreduce_dp_3d
      procedure, pass(self) :: par_env_t_bcast_dp_4d
      procedure, pass(self) :: par_env_t_allreduce_dp_4d
      procedure, pass(self) :: par_env_t_bcast_byte
      procedure, pass(self) :: par_env_t_allreduce_byte
      procedure, pass(self) :: par_env_t_bcast_c_bool
      procedure, pass(self) :: par_env_t_allreduce_c_bool

      generic :: bcast => &
              par_env_t_bcast_int32_scalar,&
              par_env_t_bcast_int32_1d,&
              par_env_t_bcast_int64_scalar,&
              par_env_t_bcast_int64_1d,&
              par_env_t_bcast_dp_scalar,&
              par_env_t_bcast_dp_1d,&
              par_env_t_bcast_dp_2d,&
              par_env_t_bcast_dp_3d,&
              par_env_t_bcast_dp_4d,&
              par_env_t_bcast_byte,&
              par_env_t_bcast_c_bool

      generic :: allreduce => &
             par_env_t_allreduce_int32_scalar,&
             par_env_t_allreduce_int32_1d,&
             par_env_t_allreduce_int64_scalar,&
             par_env_t_allreduce_int64_1d,&
             par_env_t_allreduce_dp_scalar,&
             par_env_t_allreduce_dp_1d,&
             par_env_t_allreduce_dp_2d,&
             par_env_t_allreduce_dp_3d,&
             par_env_t_allreduce_dp_4d,&
             par_env_t_allreduce_byte,&
             par_env_t_allreduce_c_bool

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
  subroutine get_hostnames(self, unique_hostnames)
    implicit none
    character(len=:), allocatable, intent(inout) :: unique_hostnames
    class(par_env_t), intent(inout) :: self
#ifdef ENABLE_MPI
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: all_processor_names(:)
    integer(PARALLEL_INT) :: i, j, unique_count, ierr, name_len
    logical :: is_duplicate
#endif
    integer(PARALLEL_INT) :: total_length, root
    total_length = 28
    root = 0
#ifndef ENABLE_MPI
    allocate(character(len=total_length) :: unique_hostnames)
    call hostnm(unique_hostnames)
#else
    if(self%use_mpi) then
      call MPI_Get_processor_name(processor_name, name_len, ierr)

      if (self%rank == 0) then
        allocate(character(len=MPI_MAX_PROCESSOR_NAME) :: all_processor_names(self%size))
      end if

      call MPI_Gather(processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                      all_processor_names, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                      root, self%comm, ierr)

      if (self%rank == 0) then
        unique_count = 0
        total_length = 0
        do i = 1, self%size
          is_duplicate = .false.
          do j = 1, unique_count
            if (trim(all_processor_names(i)) == trim(all_processor_names(j))) then
              is_duplicate = .true.
              exit
            end if
          end do

          if (.not. is_duplicate) then
            unique_count = unique_count + 1
            all_processor_names(unique_count) = all_processor_names(i)
            total_length = total_length + len_trim(all_processor_names(i)) + 1
          end if
        end do

        allocate(character(len=total_length) :: unique_hostnames)
        unique_hostnames = ''

        do i = 1, unique_count
          if (len_trim(unique_hostnames) > 0) then
            unique_hostnames = trim(unique_hostnames) // ',' // trim(adjustl(all_processor_names(i)))
          else
            unique_hostnames = trim(adjustl(all_processor_names(i)))
          end if
        end do

        deallocate(all_processor_names)
      end if

      call MPI_Bcast(total_length,int(1, kind=PARALLEL_INT), MPI_INTEGER, root, self%comm, ierr)

      if (self%rank /= 0) then
        allocate(character(len=total_length) :: unique_hostnames)
      end if
      call MPI_Bcast(unique_hostnames, total_length, MPI_CHARACTER, root, self%comm, ierr)

    else
      allocate(character(len=total_length) :: unique_hostnames)
      call hostnm(unique_hostnames)
    endif
#endif
  end subroutine get_hostnames

!##################################
  subroutine par_env_t_bcast_int32_scalar(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    integer(int32), intent(inout) :: buffer
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_int32_scalar

  subroutine par_env_t_allreduce_int32_scalar(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    integer(int32), intent(inout) :: buffer
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_int32_scalar

  subroutine par_env_t_bcast_int32_1d(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    integer(int32), intent(inout) :: buffer(:)
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_int32_1d

  subroutine par_env_t_allreduce_int32_1d(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    integer(int32), intent(inout) :: buffer(:)
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_int32_1d

  subroutine par_env_t_bcast_int64_scalar(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    integer(int64), intent(inout) :: buffer
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER8, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_int64_scalar

  subroutine par_env_t_allreduce_int64_scalar(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    integer(int64), intent(inout) :: buffer
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER8, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_int64_scalar

  subroutine par_env_t_bcast_int64_1d(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    integer(int64), intent(inout) :: buffer(:)
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER8, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_int64_1d

  subroutine par_env_t_allreduce_int64_1d(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    integer(int64), intent(inout) :: buffer(:)
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_INTEGER8, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_int64_1d

  subroutine par_env_t_bcast_dp_scalar(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_dp_scalar

  subroutine par_env_t_allreduce_dp_scalar(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_dp_scalar

  subroutine par_env_t_bcast_dp_1d(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:)
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_dp_1d

  subroutine par_env_t_allreduce_dp_1d(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:)
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_dp_1d

  subroutine par_env_t_bcast_dp_2d(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:,:)
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_dp_2d

  subroutine par_env_t_allreduce_dp_2d(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:,:)
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_dp_2d

  subroutine par_env_t_bcast_dp_3d(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:,:,:)
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_dp_3d

  subroutine par_env_t_allreduce_dp_3d(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:,:,:)
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_dp_3d

  subroutine par_env_t_bcast_dp_4d(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:,:,:,:)
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_dp_4d

  subroutine par_env_t_allreduce_dp_4d(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    real(kind=dp), intent(inout) :: buffer(:,:,:,:)
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_dp_4d

  subroutine par_env_t_bcast_byte(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    character(kind=c_char, len=1), intent(inout) :: buffer(*)
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_BYTE, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_byte

  subroutine par_env_t_allreduce_byte(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    character(kind=c_char, len=1), intent(inout) :: buffer(*)
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_BYTE, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_byte

  subroutine par_env_t_bcast_c_bool(self, buffer, length, root)
    class(par_env_t), intent(inout) :: self
    logical(c_bool), intent(inout) :: buffer
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
    call MPI_Bcast(buffer, int(length, kind=PARALLEL_INT), MPI_C_BOOL, root_, self%comm, self%err)
#endif
  end subroutine par_env_t_bcast_c_bool

  subroutine par_env_t_allreduce_c_bool(self, buffer, length)
    class(par_env_t), intent(inout) :: self
    logical(c_bool), intent(inout) :: buffer
    integer, intent(in) :: length
#ifndef ENABLE_MPI
    self%err = 0
    return
#else
    if (.not. self%use_mpi) return
    call MPI_Allreduce(MPI_IN_PLACE, buffer, int(length, kind=PARALLEL_INT), MPI_C_BOOL, MPI_SUM, self%comm, self%err)
#endif
  end subroutine par_env_t_allreduce_c_bool


end module parallel
