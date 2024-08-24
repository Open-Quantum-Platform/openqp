module elements

  use strings, only: to_upper
  use iso_fortran_env, only: real64
  use physical_constants, only: ANGSTROM_TO_BOHR

  implicit none
  private
  public get_element_id

  integer, parameter, public :: MAX_ELEMENT_Z = 110

!   Short elements name, Camel Case
  character(len=2), parameter, public :: ELEMENTS_ATOMNAME(MAX_ELEMENT_Z) = [ &
                                 'H ', 'He', 'Li', 'Be', 'B ', 'C ', &
                                 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', &
                                 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
                                 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
                                 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', &
                                 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', &
                                 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', &
                                 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', &
                                 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
                                 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', &
                                 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', &
                                 'Bh', 'Hs', 'Mt', 'Ds']

!   Short elements name, UPPER CASE
  character(len=4), parameter, public :: ELEMENTS_SHORT_NAME(MAX_ELEMENT_Z) = [ &
                                 "H   ", "HE  ", "LI  ", "BE  ", "B   ", "C   ", "N   ", "O   ", "F   ", "NE  ", &
                                 "NA  ", "MG  ", "AL  ", "SI  ", "P   ", "S   ", "CL  ", "AR  ", &
                                 "K   ", "CA  ", "SC  ", "TI  ", "V   ", "CR  ", "MN  ", "FE  ", "CO  ", "NI  ", &
                                 "CU  ", "ZN  ", "GA  ", "GE  ", "AS  ", "SE  ", "BR  ", "KR  ", &
                                 "RB  ", "SR  ", "Y   ", "ZR  ", "NB  ", "MO  ", "TC  ", "RU  ", "RH  ", "PD  ", &
                                 "AG  ", "CD  ", "IN  ", "SN  ", "SB  ", "TE  ", "I   ", "XE  ", &
                                 "CS  ", "BA  ", "LA  ", "CE  ", "PR  ", "ND  ", "PM  ", "SM  ", "EU  ", "GD  ", &
                                 "TB  ", "DY  ", "HO  ", "ER  ", "TM  ", "YB  ", "LU  ", &
                                 "HF  ", "TA  ", "W   ", "RE  ", "OS  ", "IR  ", "PT  ", &
                                 "AU  ", "HG  ", "TL  ", "PB  ", "BI  ", "PO  ", "AT  ", "RN  ", &
                                 "FR  ", "RA  ", "AC  ", "TH  ", "PA  ", "U   ", "NP  ", "PU  ", "AM  ", "CM  ", &
                                 "BK  ", "CF  ", "ES  ", "FM  ", "MD  ", "NO  ", "LR  ", &
                                 "RF  ", "DB  ", "SG  ", "BH  ", "HS  ", "MT  ", "DS  "]

!   Long elements name, UPPER CASE
  character(len=16), parameter, public :: ELEMENTS_LONG_NAME(MAX_ELEMENT_Z) = [ &
                                  "HYDROGEN        ", "HELIUM          ", "LITHIUM         ", "BERYLLIUM       ", &
                                  "BORON           ", "CARBON          ", "NITROGEN        ", "OXYGEN          ", &
                                  "FLUORINE        ", "NEON            ", "SODIUM          ", "MAGNESIUM       ", &
                                  "ALUMINIUM       ", "SILICON         ", "PHOSPHORUS      ", "SULFUR          ", &
                                  "CHLORINE        ", "ARGON           ", "POTASSIUM       ", "CALCIUM         ", &
                                  "SCANDIUM        ", "TITANIUM        ", "VANADIUM        ", "CHROMIUM        ", &
                                  "MANGANESE       ", "IRON            ", "COBALT          ", "NICKEL          ", &
                                  "COPPER          ", "ZINC            ", "GALLIUM         ", "GERMANIUM       ", &
                                  "ARSENIC         ", "SELENIUM        ", "BROMINE         ", "KRYPTON         ", &
                                  "RUBIDIUM        ", "STRONTIUM       ", "YTTRIUM         ", "ZIRCONIUM       ", &
                                  "NIOBIUM         ", "MOLYBDENUM      ", "TECHNETIUM      ", "RUTHENIUM       ", &
                                  "RHODIUM         ", "PALLADIUM       ", "SILVER          ", "CADMIUM         ", &
                                  "INDIUM          ", "TIN             ", "ANTIMONY        ", "TELLURIUM       ", &
                                  "IODINE          ", "XENON           ", "CAESIUM         ", "BARIUM          ", &
                                  "LANTHANUM       ", "CERIUM          ", "PRASEODYMIUM    ", "NEODYMIUM       ", &
                                  "PROMETHIUM      ", "SAMARIUM        ", "EUROPIUM        ", "GADOLINIUM      ", &
                                  "TERBIUM         ", "DYSPROSIUM      ", "HOLMIUM         ", "ERBIUM          ", &
                                  "THULIUM         ", "YTTERBIUM       ", "LUTETIUM        ", "HAFNIUM         ", &
                                  "TANTALUM        ", "TUNGSTEN        ", "RHENIUM         ", "OSMIUM          ", &
                                  "IRIDIUM         ", "PLATINUM        ", "GOLD            ", "MERCURY         ", &
                                  "THALLIUM        ", "LEAD            ", "BISMUTH         ", "POLONIUM        ", &
                                  "ASTATINE        ", "RADON           ", "FRANCIUM        ", "RADIUM          ", &
                                  "ACTINIUM        ", "THORIUM         ", "PROTACTINIUM    ", "URANIUM         ", &
                                  "NEPTUNIUM       ", "PLUTONIUM       ", "AMERICIUM       ", "CURIUM          ", &
                                  "BERKELIUM       ", "CALIFORNIUM     ", "EINSTEINIUM     ", "FERMIUM         ", &
                                  "MENDELEVIUM     ", "NOBELIUM        ", "LAWRENCIUM      ", "RUTHERFORDIUM   ", &
                                  "DUBNIUM         ", "SEABORGIUM      ", "BOHRIUM         ", "HASSIUM         ", &
                                  "MEITNERIUM      ", "DARMSTADTIUM    "]

! M. Mantina et. al., J. Phys. Chem. A, Vol. 113, No. 19, 2009: H-Ca, Ga-Sr, In-Ba, Tl-Ra
! S. Batsanov, Inorganic Materials, Vol. 37, No. 9, 2001, pp. 871–885: Sc-Zn, Y-Cd, Hf-Hg
! S.-Z. Hu et. al., Z. Kristallogr.224(2009) 375–383: La-Lu, Ac-Am
! Cm-Ds : 2.4
  real(real64), parameter, public :: ELEMENTS_VDW_RADII(MAX_ELEMENT_Z) = ANGSTROM_TO_BOHR*[ &
                         1.10, 1.40, 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27, &
                         1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, 2.30, 2.15, &
                         2.05, 2.05, 2.05, 2.05, 2.00, 2.00, 2.00, 2.10, 1.87, 2.11, 1.85, &
                         1.90, 1.83, 2.02, 3.03, 2.49, 2.40, 2.30, 2.15, 2.10, 2.05, 2.05, &
                         2.00, 2.05, 2.10, 2.20, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16, 3.43, &
                         2.68, 2.43, 2.42, 2.40, 2.39, 2.38, 2.36, 2.35, 2.34, 2.33, 2.31, &
                         2.30, 2.29, 2.27, 2.26, 2.24, 2.25, 2.20, 2.10, 2.05, 2.00, 2.00, &
                         2.05, 2.10, 2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, &
                         2.47, 2.45, 2.43, 2.41, 2.39, 2.37, 2.35, &
                         2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4 &
          ]

contains

  integer function get_element_id(element_name) result(id)
    character(len=*), intent(in) :: element_name
    logical :: found
    character(len=16) :: name_upcase
    found = .false.
    name_upcase(:) = to_upper(element_name)
    do id = 1, MAX_ELEMENT_Z
      if (name_upcase == ELEMENTS_LONG_NAME(id) &
          .or. name_upcase == ELEMENTS_SHORT_NAME(id)) then
        found = .true.
        exit
      end if
    end do
    if (.not. found) id = -1
  end function

end module
