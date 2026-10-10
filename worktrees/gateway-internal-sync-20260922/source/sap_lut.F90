!> @brief Radial superposition-of-atomic-potentials (SAP) lookup table.
!>
!> Loads Susi Lehtola's tabulated effective atomic charges Z_eff(r) from the
!> OpenQP data file (basis_sets/sap_grasp.dat, generated offline by
!> tools/sap/generate_sap_data.py from pyscf.dft.sap_data) and evaluates the
!> SAP potential of a neutral atom,
!>     V_A(r) = -Z_eff(r) / r,
!> by linear interpolation, matching the reference implementation
!> (S. Lehtola, J. Chem. Theory Comput. 15, 1593 (2019)).
module sap_lut

  use precision, only: dp

  implicit none

  private
  public :: sap_table_t

  type :: sap_table_t
    integer :: nr = 0          !< number of radial points
    integer :: zmax = 0        !< highest supported atomic number
    real(kind=dp), allocatable :: r(:)        !< radial grid (nr), bohr
    real(kind=dp), allocatable :: zeff(:,:)   !< Z_eff(r) (nr, zmax)
    real(kind=dp) :: rmax = 0.0_dp            !< largest tabulated r
  contains
    procedure :: load => sap_load
    procedure :: potential => sap_potential
    procedure :: clean => sap_clean
  end type

contains

  !> @brief Load the radial SAP table from a data file.
  subroutine sap_load(self, filename, err)
    use io_constants, only: IW
    class(sap_table_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    logical, intent(out) :: err

    integer :: u, ios, z
    character(len=256) :: line

    err = .false.
    call self%clean()

    open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      err = .true.
      return
    end if

    ! Skip comment lines beginning with '#'
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) then
        err = .true.
        close(u)
        return
      end if
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      exit
    end do

    ! 'line' now holds the header: zmax nr
    read(line, *, iostat=ios) self%zmax, self%nr
    if (ios /= 0 .or. self%zmax <= 0 .or. self%nr <= 0) then
      err = .true.
      close(u)
      return
    end if

    allocate(self%r(self%nr), self%zeff(self%nr, self%zmax))

    read(u, *, iostat=ios) self%r
    if (ios /= 0) then
      err = .true.
      close(u)
      return
    end if
    do z = 1, self%zmax
      read(u, *, iostat=ios) self%zeff(:, z)
      if (ios /= 0) then
        err = .true.
        close(u)
        return
      end if
    end do
    close(u)

    self%rmax = self%r(self%nr)
  end subroutine sap_load

  !> @brief SAP potential V(r) = -Z_eff(r)/r for an atom of nuclear charge z.
  !>
  !> Uses linear interpolation on the tabulated grid. Returns 0 beyond the
  !> tabulated range (Z_eff has decayed to 0 there) and the bare nuclear
  !> limit -z/r is approached as r -> 0 because Z_eff(0) = z.
  pure function sap_potential(self, z, dist) result(v)
    class(sap_table_t), intent(in) :: self
    integer, intent(in) :: z
    real(kind=dp), intent(in) :: dist
    real(kind=dp) :: v

    real(kind=dp) :: zeff_r
    integer :: lo, hi, mid
    real(kind=dp) :: t

    v = 0.0_dp
    if (z < 1 .or. z > self%zmax) return
    if (dist <= 0.0_dp) return
    if (dist >= self%rmax) return

    ! Binary search for the bracketing interval [r(lo), r(lo+1)]
    lo = 1
    hi = self%nr
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (self%r(mid) <= dist) then
        lo = mid
      else
        hi = mid
      end if
    end do

    t = (dist - self%r(lo)) / (self%r(lo+1) - self%r(lo))
    zeff_r = self%zeff(lo, z) + t * (self%zeff(lo+1, z) - self%zeff(lo, z))

    v = -zeff_r / dist
  end function sap_potential

  subroutine sap_clean(self)
    class(sap_table_t), intent(inout) :: self
    if (allocated(self%r)) deallocate(self%r)
    if (allocated(self%zeff)) deallocate(self%zeff)
    self%nr = 0
    self%zmax = 0
    self%rmax = 0.0_dp
  end subroutine sap_clean

end module sap_lut
