module physical_constants

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private

  real(real64), public, parameter :: ELECTRON_CHARGE = 1.602176634d-19

!   Recent (2018) NIST CODATA values (SI)
  real(real64), public, parameter :: BOHR_RADIUS = 5.29177210903d-11
  real(real64), public, parameter :: ANGSTROM = 1.00000000000d-10
  real(real64), public, parameter :: NANOMETER = 1.00000000000d-09

  real(real64), public, parameter :: HARTREE = 4.3597447222071d-18
  real(real64), public, parameter :: ELECTRONVOLT = 1.602176634d-19
  real(real64), public, parameter :: JOULE = 1.0d+00
  real(real64), public, parameter :: CALORIE = 4.184d+00


  real(real64), public, parameter :: STATC = 1.0d0/2997924580.0d0
  real(real64), public, parameter :: CENTIMETER = 1.0d-02

  real(real64), public, parameter :: DEBYE = 1.0d-10*STATC*ANGSTROM
  real(real64), public, parameter :: BUCKINGHAM = DEBYE*ANGSTROM
  real(real64), public, parameter :: CGS_OCT = BUCKINGHAM*ANGSTROM

  real(real64), public, parameter :: UNITS_DIPOLE = ELECTRON_CHARGE*BOHR_RADIUS
  real(real64), public, parameter :: UNITS_QUADRUPOLE = UNITS_DIPOLE*BOHR_RADIUS
  real(real64), public, parameter :: UNITS_OCTOPOLE = UNITS_QUADRUPOLE*BOHR_RADIUS

  real(real64), public, parameter :: K_BOLTZMANN = 1.380649d-23
  real(real64), public, parameter :: N_AVOGADRO = 6.02214076d+23

!   Internal default units are atomic units (length: Bohr, energy: Hartree)
!   Units of length
  real(real64), public, parameter :: UNITS_BOHR = 1.0d0
  real(real64), public, parameter :: UNITS_ANGSTROM = ANGSTROM/BOHR_RADIUS
  real(real64), public, parameter :: UNITS_NM = NANOMETER/BOHR_RADIUS

!   Units of energy
  real(real64), public, parameter :: UNITS_HARTREE = 1.0d0
  real(real64), public, parameter :: UNITS_EV = ELECTRONVOLT/HARTREE
  real(real64), public, parameter :: UNITS_KCALMOL = (1000*CALORIE/N_AVOGADRO)/HARTREE
  real(real64), public, parameter :: UNITS_KJMOL = (1000*JOULE/N_AVOGADRO)/HARTREE

!   Common conversion constants
  real(real64), public, parameter :: BOHR_TO_ANGSTROM = UNITS_BOHR/UNITS_ANGSTROM
  real(real64), public, parameter :: ANGSTROM_TO_BOHR = UNITS_ANGSTROM/UNITS_BOHR

!  Electric moments
  real(real64), public, parameter :: AU_TO_DEBYE = UNITS_DIPOLE/DEBYE
  real(real64), public, parameter :: AU_TO_BUCK = UNITS_QUADRUPOLE/BUCKINGHAM
  real(real64), public, parameter :: AU_TO_OCT = UNITS_OCTOPOLE/CGS_OCT

!   COMPATIBILITY:
  real(real64), public, parameter :: AU2ANG = BOHR_TO_ANGSTROM
  real(real64), public, parameter :: EV2HTREE = HARTREE/ELECTRONVOLT

end module
