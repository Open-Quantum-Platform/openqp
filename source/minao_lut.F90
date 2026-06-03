!> @brief MINAO atomic-density lookup table.
!>
!> Loads spherically-averaged neutral-atom density matrices in a fixed minimal
!> reference basis (STO-3G, Cartesian) from the OpenQP data file
!> (basis_sets/minao_sto3g.dat, generated offline by
!> tools/minao/generate_minao_data.py from PySCF's atomic HF densities). These
!> are superposed block-diagonally and projected onto the target basis to form
!> the MINAO initial guess.
module minao_lut

  use precision, only: dp

  implicit none

  private
  public :: minao_table_t

  type :: minao_elem_t
    integer :: nao = 0
    real(kind=dp) :: nelec = 0.0_dp
    real(kind=dp), allocatable :: dm(:,:)   !< (nao,nao) atomic density, AO basis
  end type

  type :: minao_table_t
    integer :: zmax = 0
    type(minao_elem_t), allocatable :: elem(:)
  contains
    procedure :: load => minao_load
    procedure :: clean => minao_clean
  end type

contains

  subroutine minao_load(self, filename, err)
    class(minao_table_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    logical, intent(out) :: err

    integer :: u, ios, z, zz, nao, i, j
    character(len=256) :: line
    real(kind=dp), allocatable :: flat(:)

    err = .false.
    call self%clean()

    open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      err = .true.
      return
    end if

    ! first non-comment line: zmax
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) then
        err = .true.; close(u); return
      end if
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      exit
    end do
    read(line, *, iostat=ios) self%zmax
    if (ios /= 0 .or. self%zmax <= 0) then
      err = .true.; close(u); return
    end if

    allocate(self%elem(self%zmax))

    do z = 1, self%zmax
      read(u, *, iostat=ios) zz, nao, self%elem(z)%nelec
      if (ios /= 0 .or. zz /= z) then
        err = .true.; close(u); return
      end if
      self%elem(z)%nao = nao
      allocate(self%elem(z)%dm(nao, nao), flat(nao*nao))
      read(u, *, iostat=ios) flat
      if (ios /= 0) then
        err = .true.; close(u); return
      end if
      ! data written row-major
      do i = 1, nao
        do j = 1, nao
          self%elem(z)%dm(i, j) = flat((i-1)*nao + j)
        end do
      end do
      deallocate(flat)
    end do
    close(u)
  end subroutine minao_load

  subroutine minao_clean(self)
    class(minao_table_t), intent(inout) :: self
    integer :: z
    if (allocated(self%elem)) then
      do z = 1, size(self%elem)
        if (allocated(self%elem(z)%dm)) deallocate(self%elem(z)%dm)
      end do
      deallocate(self%elem)
    end if
    self%zmax = 0
  end subroutine minao_clean

end module minao_lut
